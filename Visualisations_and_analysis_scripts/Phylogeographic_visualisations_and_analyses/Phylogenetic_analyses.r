library(diagram)
library(lubridate)
library(seraphim)
library(treeio)
library(viridis)

analysis_date = "270221"
writingFiles = FALSE; savingPlots = FALSE

postcodes = shapefile("Shapefiles_study_area/England_postcode_districts.shp"); crs(postcodes) = CRS("+init=epsg:27700")
UTLAs = crop(spTransform(shapefile("Shapefiles_study_area/UTLA_administrative_areas.shp"), crs(postcodes)), postcodes)
regions = shapefile("Shapefiles_study_area/England_defined_regions_2.shp")
pop = projectRaster(raster("WorldPop_pop_raster.tif"), crs=crs(postcodes)); pop[] = log(pop[])

# 1. Preparing the continuous phylogeographic analysis (RRW, Cauchy model)

tree = read.nexus(paste0("TreeTime_",analysis_date,".tre"))
metadata = read.csv(paste0("Metadata_England.csv"))
all_collection_dates = dmy(gsub("\\/","-",metadata[,"sample_date"]))
print(c(min(all_collection_dates),max(all_collection_dates)))
if (!file.exists(paste0("Sampling_",analysis_date,".csv")))
	{
		samplingData = matrix(nrow=length(tree$tip.label), ncol=5)
		colnames(samplingData) = c("sequence_ID","collection_date","postcode","longitude","latitude")
		samplingData[,"sequence_ID"] = gsub("'","",tree$tip.label); i = 1
		for (i in i:dim(samplingData)[1])
			{
				index = which(metadata[,"name"]==samplingData[i,"sequence_ID"])
				if (length(index) == 1)
					{
						samplingData[i,"collection_date"] = decimal_date(dmy(gsub("\\/","-",metadata[index,"sample_date"])))
						samplingData[i,"postcode"] = metadata[index,"postcode"]
						indices = which(postcodes@data[,"PostDist"]==samplingData[i,"postcode"])
						if (length(indices) > 0)
							{
								maxArea = 0; polIndex1 = 0; polIndex2 = 0
								for (j in 1:length(indices))
									{
										for (k in 1:length(postcodes@polygons[[indices[j]]]@Polygons))
											{
												if (maxArea < postcodes@polygons[[indices[j]]]@Polygons[[k]]@area)
													{
														maxArea = postcodes@polygons[[indices[j]]]@Polygons[[k]]@area; polIndex1 = indices[j]; polIndex2 = k
													}
											}
									}
								pol = postcodes@polygons[[polIndex1]]@Polygons[[polIndex2]]
								p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
								pol = sps; # proj4string(pol) = postcodes@proj4string
								# samplingData[i,c("longitude","latitude")] = coordinates(pol) # to avoid a jitter:
								samplingData[i,c("longitude","latitude")] = spsample(pol, 1, type="random")@coords
							}
					}
			}
		write.csv(samplingData, paste0("Sampling_",analysis_date,".csv"), quote=F, row.names=F); print(c(i,dim(samplingData)[1]))
	}	
samplingData = read.csv(paste0("Sampling_",analysis_date,".csv"), head=T); clusters2 = list(); clusters2[[1]] = samplingData
template = scan("RRW_template_file.xml", what="", sep="\n", quiet=T, blank.lines.skip=F)
template = gsub("All_clades","TreeTime_090221",template); xml = c()
sink(file=paste0("TreeTime_",analysis_date,".xml"))
for (i in 1:length(template))
	{
		cat(template[i],"\n")
		if (grepl("Insert taxa blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<taxa id=\"taxa_",j,"\">","\n"))
								for (k in 1:dim(clusters2[[j]])[1])
									{
										if (!is.na(clusters2[[j]][k,"longitude"]))
											{
												cat(paste0("\t\t<taxon id=\"",clusters2[[j]][k,"sequence_ID"],"\">","\n"))
												cat(paste0("\t\t\t<date value=\"",clusters2[[j]][k,"collection_date"],"\" direction=\"forwards\" units=\"years\"/>","\n"))
												cat("\t\t\t<attr name=\"latitude\">\n")
												cat(paste0("\t\t\t\t",clusters2[[j]][k,"latitude"],"\n"))
												cat("\t\t\t</attr>\n")
												cat("\t\t\t<attr name=\"longitude\">\n")
												cat(paste0("\t\t\t\t",clusters2[[j]][k,"longitude"],"\n"))
												cat("\t\t\t</attr>\n")
												cat("\t\t\t<attr name=\"coordinates\">\n")
												cat(paste0("\t\t\t\t",clusters2[[j]][k,"latitude"]," ",clusters2[[j]][k,"longitude"],"\n"))
												cat("\t\t\t</attr>\n")
												cat("\t\t</taxon>\n")
											}
									}
								cat("\t</taxa>","\n")
							}
					}
			}
		if (grepl("Insert alignment blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<alignment id=\"alignment_",j,"\" dataType=\"nucleotide\">","\n"))
								for (k in 1:dim(clusters2[[j]])[1])
									{
										if (!is.na(clusters2[[j]][k,"longitude"]))
											{
												cat("\t\t<sequence>\n")
												cat(paste0("\t\t\t<taxon idref=\"",clusters2[[j]][k,"sequence_ID"],"\"/>","\n"))
												cat("\t\t\tNNNN\n")
												cat("\t\t</sequence>\n")
											}
									}
								cat("\t</alignment>","\n")
							}
					}
			}
		if (grepl("Insert pattern blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<patterns id=\"patterns_",j,"\" from=\"1\" strip=\"false\">","\n"))
								cat(paste0("\t\t<alignment idref=\"alignment_",j,"\"/>","\n"))
								cat("\t</patterns>","\n")
							}
					}
			}
		if (grepl("Insert starting tree blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{	
								tips = clusters2[[j]][,"sequence_ID"]; tips = tips[which(!is.na(clusters2[[j]][,"longitude"]))]
								tre = tree; tips_to_drop = tre$tip.label[which(!gsub("'","",tre$tip.label)%in%tips)]
								if (length(tips_to_drop) > 0) tre = ape::drop.tip(tre, tips_to_drop)
								tre = multi2di(tre); tree$edge.length[which(tree$edge.length==0)] = 0.5
								write.tree(tre, paste0("Empirical_tree_RRW.tre"))
								tre = scan(paste0("Empirical_tree_RRW.tre"), what="", sep="\n", quiet=T)
								txt = c("#NEXUS","begin trees;",paste0("\ttree tree_1 = [&R] ",tre),"end;")
								write(txt, paste0("Empirical_tree_RRW.tre"))
								cat(paste0("\t<empiricalTreeDistributionModel id=\"treeModel_",j,"\" fileName=\"Empirical_tree_RRW.tre\">","\n"))
								cat(paste0("\t\t<taxa idref=\"taxa_",j,"\"/>","\n"))
								cat("\t</empiricalTreeDistributionModel>","\n")
							}
					}
			}
		if (grepl("Insert tree model blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<treeModel id=\"treeModel_",j,"\">","\n"))
								cat(paste0("\t\t<coalescentTree idref=\"startingTree_",j,"\"/>","\n"))
								cat("\t\t<rootHeight>","\n")
								cat(paste0("\t\t\t<parameter id=\"treeModel.rootHeight_",j,"\"/>","\n"))
								cat("\t\t</rootHeight>","\n")
								cat("\t\t<nodeHeights internalNodes=\"true\">","\n")
								cat(paste0("\t\t\t<parameter id=\"treeModel.internalNodeHeights_",j,"\"/>","\n"))
								cat("\t\t</nodeHeights>","\n")
								cat("\t\t<nodeHeights internalNodes=\"true\" rootNode=\"true\">","\n")
								cat(paste0("\t\t\t<parameter id=\"treeModel.allInternalNodeHeights_",j,"\"/>","\n"))
								cat("\t\t</nodeHeights>","\n")
								cat("\t</treeModel>","\n")
							}
					}
			}
		if (grepl("Insert arbitraryBranchRates blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<arbitraryBranchRates id=\"coordinates.diffusion.branchRates",j,"\">","\n"))
								cat(paste0("\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
								cat("\t\t<rates>","\n")
								cat(paste0("\t\t\t<parameter id=\"coordinates.diffusion.rates",j,"\" lower=\"0.0\"/>","\n"))
								cat("\t\t</rates>","\n")
								cat("\t</arbitraryBranchRates>","\n")
							}
					}
			}
		if (grepl("Insert distributionLikelihood blocks 1",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<distributionLikelihood id=\"coordinates.diffusion.prior",j,"\">","\n"))
								cat("\t\t<data>","\n")
								cat(paste0("\t\t\t<parameter idref=\"coordinates.diffusion.rates",j,"\"/>","\n"))
								cat("\t\t</data>","\n")
								cat("\t\t<distribution>","\n")
								cat(paste0("\t\t\t<onePGammaDistributionModel>","\n"))
								cat("\t\t\t\t<shape>","\n")
								cat("\t\t\t\t\t<parameter value=\"0.5\"/>","\n")
								cat("\t\t\t\t</shape>","\n")
								cat("\t\t\t</onePGammaDistributionModel>","\n")
								cat("\t\t</distribution>","\n")
								cat("\t</distributionLikelihood>","\n")
							}
					}
			}
		if (grepl("Insert coordinates.traitLikelihood blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<multivariateTraitLikelihood id=\"coordinates.traitLikelihood",j,"\" traitName=\"coordinates\" useTreeLength=\"true\" scaleByTime=\"true\" reportAsMultivariate=\"true\" reciprocalRates=\"true\" integrateInternalTraits=\"true\">","\n"))
								cat("\t\t<multivariateDiffusionModel idref=\"coordinates.diffusionModel\"/>","\n")
								cat(paste0("\t\t<treeModel idref=\"treeModel_",j,"\"/>"))
								cat("\t\t<traitParameter>","\n")
								cat(paste0("\t\t\t<parameter id=\"leaf.coordinates",j,"\"/>","\n"))
								cat("\t\t</traitParameter>","\n")
								cat("\t\t<conjugateRootPrior>","\n")
								cat("\t\t\t<meanParameter>","\n")
								cat("\t\t\t\t<parameter value=\"0.0 0.0\"/>","\n")
								cat("\t\t\t</meanParameter>","\n")
								cat("\t\t\t<priorSampleSize>","\n")
								cat("\t\t\t\t<parameter value=\"0.000001\"/>","\n")
								cat("\t\t\t</priorSampleSize>","\n")
								cat("\t\t</conjugateRootPrior>","\n")
								cat(paste0("\t\t<arbitraryBranchRates idref=\"coordinates.diffusion.branchRates",j,"\"/>","\n"))
								cat("\t</multivariateTraitLikelihood>","\n")
							}
					}
			}
		if (grepl("Insert continuousDiffusionStatistic blocks 1",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t<continuousDiffusionStatistic id=\"coordinates.diffusionRate",j,"\" greatCircleDistance=\"true\">","\n"))
								cat(paste0("\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
								cat("\t</continuousDiffusionStatistic>","\n")
							}
					}
			}
		if (grepl("Insert scaleOperator blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t<scaleOperator scaleFactor=\"0.75\" weight=\"30\">","\n"))
								cat(paste0("\t\t\t<parameter idref=\"coordinates.diffusion.rates",j,"\"/>","\n"))
								cat("\t\t</scaleOperator>","\n")
							}
					}
			}
		if (grepl("Insert precisionGibbsOperator blocks",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t<precisionGibbsOperator weight=\"2\">","\n"))
								cat(paste0("\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
								cat("\t\t\t<multivariateWishartPrior idref=\"coordinates.precisionPrior\"/>","\n")
								cat("\t\t</precisionGibbsOperator>","\n")
							}
					}
			}
		if (grepl("Insert distributionLikelihood blocks 2",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t\t\t<distributionLikelihood idref=\"coordinates.diffusion.prior",j,"\"/>","\n"))
							}
					}
			}
		if (grepl("Insert multivariateTraitLikelihood blocks 1",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
							}
					}
			}
		if (grepl("Insert continuousDiffusionStatistic blocks 2",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t\t\t<continuousDiffusionStatistic idref=\"coordinates.diffusionRate",j,"\"/>","\n"))
							}
					}
			}
		if (grepl("Insert multivariateTraitLikelihood blocks 2",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
							}
					}
			}
		if (grepl("<!-- Insert logTree blocks -->",template[i]))
			{
				for (j in 1:length(clusters2))
					{
						if ((dim(clusters2[[j]])[1] >= 3)&(sum(!is.na(clusters2[[j]][,"longitude"])) >= 3))
							{
								cat(paste0("\t\t<logTree id=\"treeFileLog",j,"\" logEvery=\"100000\" nexusFormat=\"true\" fileName=\"TreeTime_",analysis_date,".trees\" sortTranslationTable=\"true\">","\n"))
								cat(paste0("\t\t\t<treeModel idref=\"treeModel_",j,"\"/>","\n"))
								cat("\t\t\t<joint idref=\"joint\"/>","\n")
								cat("\t\t\t<trait name=\"coordinates\" tag=\"coordinates\">","\n")
								cat(paste0("\t\t\t\t<multivariateTraitLikelihood idref=\"coordinates.traitLikelihood",j,"\"/>","\n"))
								cat("\t\t\t</trait>","\n")
								cat("\t\t\t<multivariateDiffusionModel idref=\"coordinates.diffusionModel\"/>","\n")
								cat("\t\t\t<trait name=\"rate\" tag=\"coordinates.rate\">","\n")
								cat(paste0("\t\t\t\t<arbitraryBranchRates idref=\"coordinates.diffusion.branchRates",j,"\"/>","\n"))
								cat("\t\t\t</trait>","\n")
								cat("\t\t</logTree>","\n")
							}
					}
			}
	}
sink(NULL)

# 2. Extracting spatio-temporal information embedded in MCC and posterior trees

if (!file.exists(paste0("TreeTime_",analysis_date,".csv")))
	{
		source("MCC_tree_extraction.r")
		samplingData = read.csv(paste0("Sampling_",analysis_date,".csv"), head=T)
		samplingData = samplingData[which(!is.na(samplingData[,"longitude"])),]
		mostRecentSamplingDatum = max(samplingData[,"collection_date"])
		mcc_tre = readAnnotatedNexus(paste0("TreeTime_",analysis_date,".tree"))
		mcc = MCC_tree_extraction(mcc_tre, mostRecentSamplingDatum)
		write.csv(mcc, paste0("TreeTime_",analysis_date,".csv"), quote=F, row.names=F)
			# before continuing, longitude/latitude columns have to be inverted in the CSV file
		midYears = matrix(nrow=dim(mcc)[1], ncol=1)
		midYears[,1] = (mcc[,"startYear"]+mcc[,"endYear"])/2
		colnames(midYears) = c("midYear"); mcc = cbind(mcc, midYears)
		start_and_end_UTLAs = matrix(nrow=dim(mcc)[1], ncol=2)
		colnames(start_and_end_UTLAs) = c("startUTLA","endUTLA")
		for (i in 1:dim(UTLAs@data)[1])
			{
				maxArea = 0; polIndex = 0; # print(i)
				for (j in 1:length(UTLAs@polygons[[i]]@Polygons))
					{
						if (maxArea < UTLAs@polygons[[i]]@Polygons[[j]]@area)
							{
								maxArea = UTLAs@polygons[[i]]@Polygons[[j]]@area; polIndex = j
							}
					}
				pol = UTLAs@polygons[[i]]@Polygons[[polIndex]]
				indices = which(point.in.polygon(mcc[,"startLon"],mcc[,"startLat"],pol@coords[,1],pol@coords[,2])==1)
				start_and_end_UTLAs[indices,"startUTLA"] = UTLAs@data[i,"NUMBER"]
				indices = which(point.in.polygon(mcc[,"endLon"],mcc[,"endLat"],pol@coords[,1],pol@coords[,2])==1)
				start_and_end_UTLAs[indices,"endUTLA"] = UTLAs@data[i,"NUMBER"]
			}
		start_and_end_regions = matrix(nrow=dim(mcc)[1], ncol=2)
		colnames(start_and_end_regions) = c("startRegion","endRegion")
		for (i in 1:(dim(regions@data)[1]-1))
			{
				maxArea = 0; polIndex = 0; # print(i)
				for (j in 1:length(regions@polygons[[i]]@Polygons))
					{
						if (maxArea < regions@polygons[[i]]@Polygons[[j]]@area)
							{
								maxArea = regions@polygons[[i]]@Polygons[[j]]@area; polIndex = j
							}
					}
				pol = regions@polygons[[i]]@Polygons[[polIndex]]
				indices = which(point.in.polygon(mcc[,"startLon"],mcc[,"startLat"],pol@coords[,1],pol@coords[,2])==1)
				if (length(indices) > 0)
					{
						if (sum(!is.na(start_and_end_regions[indices,"startRegion"])) != 0) print(c(i,1))
						start_and_end_regions[indices,"startRegion"] = regions@data[i,"Merged_loc"]
					}
				indices = which(point.in.polygon(mcc[,"endLon"],mcc[,"endLat"],pol@coords[,1],pol@coords[,2])==1)
				if (length(indices) > 0)
					{
						if (sum(!is.na(start_and_end_regions[indices,"endRegion"])) != 0) print(c(i,2))
						start_and_end_regions[indices,"endRegion"] = regions@data[i,"Merged_loc"]
					}
			}
		mcc = cbind(mcc, start_and_end_UTLAs, start_and_end_regions)
		startCoords = mcc[,c("startLon","startLat")]; coordinates(startCoords) = c("startLon","startLat")
		proj4string(startCoords) = CRS("+init=epsg:27700"); startCoords = spTransform(startCoords, CRS("+init=epsg:4326"))
		endCoords = mcc[,c("endLon","endLat")]; coordinates(endCoords) = c("endLon","endLat")
		proj4string(endCoords) = CRS("+init=epsg:27700"); endCoords = spTransform(endCoords, CRS("+init=epsg:4326"))
		geographicDistances = matrix(nrow=dim(startCoords@coords)[1], ncol=1); colnames(geographicDistances) = c("geoDist_km")
		for (i in 1:dim(geographicDistances)[1])
			{
				geographicDistances[i,1] = rdist.earth(cbind(startCoords@coords[i,1],startCoords@coords[i,2]), cbind(endCoords@coords[i,1],endCoords@coords[i,2]), miles=F)
			}
		mcc = cbind(mcc, geographicDistances)
		write.csv(mcc, paste0("TreeTime_",analysis_date,".csv"), quote=F, row.names=F)
	}	else		{
		mcc = read.csv(paste0("TreeTime_",analysis_date,".csv"), head=T)
	}
mcc = mcc[-which((is.na(mcc[,"startUTLA"]))|(is.na(mcc[,"endUTLA"]))),]
mcc = mcc[order(mcc[,"endYear"]),]

# 3. Generating the dispersal history graphs (mapped MCC trees, 80% HPD polygons)

if (savingPlots)
	{
		colourScale = rev(colorRampPalette(brewer.pal(11,"BrBG"))(141)[16:116])
		colourScale = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(141)[16:116])
		colourScale = gsub("FF","",viridis::viridis(101)[1:101])
		minYear = min(mcc[,"startYear"]); maxYear = max(mcc[,"endYear"])
		startYears_indices = (((mcc[,"startYear"]-minYear)/(maxYear-minYear))*100)+1
		endYears_indices = (((mcc[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
		startYears_colours = colourScale[startYears_indices]
		endYears_colours = colourScale[endYears_indices]
		plottingInternalNodes = TRUE; plotLegend = TRUE
		UTLA_sampled_sequences = matrix(nrow=dim(UTLAs@data)[1], ncol=1)
		for (i in 1:dim(UTLA_sampled_sequences)[1])
			{
				UTLA_sampled_sequences[i,1] = length(which((!mcc[,"node2"]%in%mcc[,"node1"])&(mcc[,"endUTLA"]==UTLAs@data[i,"NUMBER"])))
			}
		UTLA_all_introductions = matrix(nrow=dim(UTLAs@data)[1], ncol=1)
		for (i in 1:dim(UTLA_all_introductions)[1])
			{
				UTLA_all_introductions[i,1] = length(which((mcc[,"startUTLA"]!=UTLAs@data[i,"NUMBER"])&(mcc[,"endUTLA"]==UTLAs@data[i,"NUMBER"])))
			}
		pdf("Figure_A4_OLD.pdf", width=3, height=3.5)
		par(oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
		minVal = 0; maxVal = max(UTLA_sampled_sequences)
		indices = (((UTLA_sampled_sequences[,1]-minVal)/(maxVal-minVal))*100)+1
		cols = rev(colourScale)[indices]; cols[which(UTLA_sampled_sequences[,1]==0)] = NA
		plot(postcodes, border=NA, col="gray90", lwd=0.01)
		plot(UTLAs, border="white", col=cols, add=T, lwd=0.2)
		if (plotLegend == TRUE)
			{
				rast = raster(matrix(nrow=1, ncol=2)); rast[1] = 1; rast[2] = max(UTLA_sampled_sequences)
				plot(rast, legend.only=T, add=T, col=rev(colourScale), legend.width=0.5, legend.shrink=0.3, smallplot=c(0.45,0.94,0.114,0.125),
					 legend.args=list(text="", cex=0.6, line=0.3, col="gray30"), horizontal=T,
				     axis.args=list(cex.axis=0.5, lwd=0, lwd.tick=0.2, tck=-0.9, col.axis="gray30", line=0, mgp=c(0,-0.20,0)))
			}
		dev.off()
		pdf("Figure_A4_NEW.pdf", width=3, height=3.5)
		par(oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
		minVal = 0; maxVal = max(UTLA_all_introductions)
		indices = (((UTLA_all_introductions[,1]-minVal)/(maxVal-minVal))*100)+1
		cols = rev(colourScale)[indices]; cols[which(UTLA_all_introductions[,1]==0)] = NA
		plot(postcodes, border=NA, col="gray90", lwd=0.01)
		plot(UTLAs, border="white", col=cols, add=T, lwd=0.2)
		if (plotLegend == TRUE)
			{
				rast = raster(matrix(nrow=1, ncol=2)); rast[1] = 1; rast[2] = max(UTLA_all_introductions)
				plot(rast, legend.only=T, add=T, col=rev(colourScale), legend.width=0.5, legend.shrink=0.3, smallplot=c(0.45,0.94,0.114,0.125),
					 legend.args=list(text="", cex=0.6, line=0.3, col="gray30"), horizontal=T,
				     axis.args=list(cex.axis=0.5, lwd=0, lwd.tick=0.2, tck=-0.9, col.axis="gray30", line=0, mgp=c(0,-0.20,0)))
			}
		dev.off()	
		pdf("Figure_A2_NEW.pdf", width=6, height=7)
		par(oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
		plot(postcodes, border=NA, col="gray90", lwd=0.01)
		plot(UTLAs, border="white", col=NA, add=T, lwd=0.4)
		selectedBranches = 1:dim(mcc)[1]; cexNode = 0.3
		for (i in selectedBranches)
			{
				curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
						    arr.width=0, lwd=0.1, lty=1, lcol="gray10", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
			}
		for (i in rev(selectedBranches))
			{
				if (plottingInternalNodes == TRUE)
					{
						if (mcc[i,"node2"]%in%mcc[,"node1"])
							{
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=cexNode)
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray10", cex=cexNode, lwd=0.1)
							}
					}
				if (!mcc[i,"node2"]%in%mcc[,"node1"])
					{
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=cexNode)
						points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray10", cex=cexNode, lwd=0.1)
					}
			}
		rootIndex = which(!mcc[,"node1"]%in%mcc[,"node2"])[1]
		points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=startYears_colours[rootIndex], cex=cexNode)
		points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray10", cex=cexNode, lwd=0.1)
		if (plotLegend == TRUE)
			{
				selectedDates = decimal_date(ymd(c("2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01")))
				selectedLabels = c("01-09-20","01-10-20","01-11-20","01-12-20","01-01-21")
				rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(mcc[,"startYear"]); rast[2] = max(mcc[,"endYear"])
				plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.46,0.92,0.106,0.112),
					 legend.args=list(text="", cex=0.55, line=0.3, col="gray30"), horizontal=T,
				     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-1, col.axis="gray30", line=0, mgp=c(0,-0.1,0),
				     at=selectedDates, labels=selectedLabels))
			}
		dev.off()
		pdf("Figure_A1_NEW.pdf", width=8, height=9)
		par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30"); cexNode = 0.35
		cutOffs = c(decimal_date(dmy(c("05-11-2020","01-12-2020","20-12-2020","12-01-2021"))))
		dates = c("05-11-2020","01-12-2020","20-12-2020","12-01-2021")
		for (h in 1:length(cutOffs))
			{
				plot(postcodes, border=NA, col="gray95", lwd=0.01)
				plot(UTLAs, border="white", col=NA, add=T, lwd=0.4)
				selectedBranches = which(mcc[,"endYear"]<cutOffs[h])
				for (i in selectedBranches)
					{
						curvedarrow(cbind(mcc[i,"startLon"],mcc[i,"startLat"]), cbind(mcc[i,"endLon"],mcc[i,"endLat"]), arr.length=0,
								    arr.width=0, lwd=0.1, lty=1, lcol="gray10", arr.col=NA, arr.pos=F, curve=0.1, dr=NA, endhead=F)
					}
				for (i in rev(selectedBranches))
					{
						if (plottingInternalNodes == TRUE)
							{
								if (mcc[i,"node2"]%in%mcc[,"node1"])
									{
										points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=cexNode)
										points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray10", cex=cexNode, lwd=0.1)
									}
							}
						if (!mcc[i,"node2"]%in%mcc[,"node1"])
							{
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=16, col=endYears_colours[i], cex=cexNode)
								points(mcc[i,"endLon"], mcc[i,"endLat"], pch=1, col="gray10", cex=cexNode, lwd=0.1)
							}
					}
				rootIndex = which(!mcc[,"node1"]%in%mcc[,"node2"])[1]
				points(mcc[i,"startLon"], mcc[i,"startLat"], pch=16, col=startYears_colours[rootIndex], cex=cexNode)
				points(mcc[i,"startLon"], mcc[i,"startLat"], pch=1, col="gray10", cex=cexNode, lwd=0.1)
				mtext(dates[h], line=-18.5, at=225000, cex=0.7, col="gray30")
			}
		dev.off()
	}

# 4. Analysing of B.1.1.7 dispersal events on a weekly basis

tab1 = read.csv("London_stats_MK.csv", head=T)
dates = decimal_date(dmy(gsub("\\/","-",tab1[,"week_start"])))
tab1 = tab1[order(dates),]; dates = dates[order(dates)]
tab2 = matrix(0, nrow=dim(tab1)[1], ncol=5)
colnames(tab2) = c("estimated_exports_London","all_lineage_dispersal","lineage_exports_London","median_distance_London","median_distance_notLondon")
tab2[,"estimated_exports_London"] = tab1[,"B.1.1.7_exports"]
samplingData = read.csv(paste0("Sampling_",analysis_date,".csv"), head=T)
samplingData = samplingData[which(!is.na(samplingData[,"longitude"])),]
mostRecentSamplingDatum = max(samplingData[,"collection_date"])
usingMidBranchTimes = FALSE
for (i in 1:dim(tab1)[1])
	{
		if (i != dim(tab1)[1])
			{
				ending_day = dates[i+1]
			}	else	{
				ending_day = mostRecentSamplingDatum+1/365
			}
		if (usingMidBranchTimes == TRUE) sub1 = mcc[which((mcc[,"midYear"]>=dates[i])&(mcc[,"midYear"]<ending_day)),]
		if (usingMidBranchTimes == FALSE) sub1 = mcc[which((mcc[,"endYear"]>=dates[i])&(mcc[,"endYear"]<ending_day)),]
		if (dim(sub1)[1] > 0)
			{
				tab2[i,"all_lineage_dispersal"] = dim(sub1)[1]
				sub2 = sub1[which(sub1[,"startRegion"]=="Greater_London"),]
				if (dim(sub2)[1] > 0)
					{
						tab2[i,"lineage_exports_London"] = dim(sub2[which(sub2[,"endRegion"]!="Greater_London"),])[1]
						sub3 = sub2[which(sub2[,"geoDist_km"]>=5),]
						if (dim(sub3)[1] > 0)
							{
								tab2[i,"median_distance_London"] = median(sub3[,"geoDist_km"])
							}
					}
				sub2 = sub1[which(sub1[,"startRegion"]!="Greater_London"),]
				if (dim(sub2)[1] > 0)
					{
tab2[i,"lineage_exports_London"] = dim(sub2[which(sub2[,"endRegion"]!="Greater_London"),])[1]
						sub3 = sub2[which(sub2[,"geoDist_km"]>=5),]
						if (dim(sub3)[1] > 0)
							{
								tab2[i,"median_distance_notLondon"] = median(sub3[,"geoDist_km"])
							}
					}
			}
	}
if (writingFiles)
	{
		tmp = matrix(nrow=dim(tab2)[1], ncol=1); tmp[,1] = tab2[,"lineage_exports_London"]
		tmp[,1] = round(tab2[,"lineage_exports_London"]/tab2[,"all_lineage_dispersal"], 3)
		allLabels = gsub("\\/","-",gsub("\\/20","\\/",tab1[,"week_start"]))
		row.names(tmp) = allLabels; colnames(tmp) = c("ratio_exports_from_London_vs_all")
		write.csv(tmp, "Ratios_between_exports_from_London_and_all_dispersal_events.csv", quote=F)
	}
tab3 = matrix(nrow=dim(tab1)[1], ncol=9)
colnames(tab3) = regions@data[1:9,"Merged_loc"]
for (i in 1:dim(tab1)[1])
	{
		if (i != dim(tab1)[1])
			{
				ending_day = dates[i+1]
			}	else	{
				ending_day = mostRecentSamplingDatum+1/365
			}
		sub1 = mcc[which((mcc[,"endYear"]>=dates[i])&(mcc[,"endYear"]<ending_day)),]
		if (dim(sub1)[1] > 0)
			{
				for (j in 1:dim(tab3)[2])
					{
						sub2 = sub1[which((sub1[,"startRegion"]!=colnames(tab3)[j])&(sub1[,"endRegion"]==colnames(tab3)[j])),]
						if (dim(sub2)[1] > 0)
							{
								sub3 = mcc[which((mcc[,"startRegion"]==colnames(tab3)[j])&(mcc[,"endRegion"]==colnames(tab3)[j])),]
								nberOfSamplesConnectedToEachIntroduction = rep(NA, dim(sub2)[1])
								for (k in 1:dim(sub2)[1])
									{
										sub4 = sub2[k,]; allBranchesScreened = FALSE
										while (allBranchesScreened == FALSE)
											{
												if (sum((sub3[,"node1"]%in%sub4[,"node2"])&(!sub3[,"node1"]%in%sub4[,"node1"])) > 0)
													{
														sub4 = rbind(sub4, sub3[which((sub3[,"node1"]%in%sub4[,"node2"])&(!sub3[,"node1"]%in%sub4[,"node1"])),])
													}	else	{
														allBranchesScreened = TRUE
													}
											}
										sub5 = sub4[which(!sub4[,"node2"]%in%sub4[,"node1"]),]
										nberOfSamplesConnectedToEachIntroduction[k] = dim(sub5)[1]
									}
								tab3[i,j] = round(mean(nberOfSamplesConnectedToEachIntroduction),2)
							}
					}
			}
	}
tab3 = tab3[,c("Kent","Greater_London","Norfolk","Birmingham","Greater_Manchester","Cumbria")]
allLabels = gsub("\\/","-",gsub("\\/20","\\/",tab1[,"week_start"])); row.names(tab3) = allLabels
if (writingFiles) write.csv(tab3, "Mean_number_of_samples_resulting_from_introduction_events.csv", quote=F)
if (savingPlots)
	{
		pdf("Figure_A3_NEW.pdf", width=3.75, height=2.3)
		par(oma=c(0,0,0,0), mar=c(3,2,0,2), lwd=0.2, col="gray30")
		plot(1:dim(tab2)[1], tab2[,"lineage_exports_London"], type="l", lwd=0.3, col="red", axes=F, ann=F)
		points(1:dim(tab2)[1], tab2[,"lineage_exports_London"], cex=0.6, pch=16, col="red")
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.45,0), lwd=0.2, tck=-0.03, col.tick="gray30", col.axis="gray30", col="gray30", las=2, at=1:dim(tab2)[1], labels=allLabels)
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.2, tck=-0.03, col.tick="gray30", col.axis="red", col="gray30")
		title(ylab="B.1.1.7 lineage exports", cex.lab=0.7, mgp=c(1.1,0,0), col.lab="red")
		par(new=T)
		plot(1:dim(tab2)[1], tab2[,1], type="l", lwd=0.3, col="gray30", axes=F, ann=F); points(1:dim(tab2)[1], tab2[,1], cex=0.6, pch=16, col="gray30")
		axis(side=4, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.10,0), lwd=0.2, tck=-0.03, col.tick="gray30", col.axis="gray30", col="gray30")
		dev.off()
		pdf("Figure_B3_NEW.pdf", width=7, height=2.0)
		par(oma=c(0,0,0,0), mar=c(3,2,0,0), lwd=0.2, col="gray30")
		plot(1:(dim(tab2)[1]-1), tab2[1:(dim(tab2)[1]-1),"median_distance_notLondon"], ylim=c(0,30), type="l", lwd=0.3, col="gray30", axes=F, ann=F)
		points(1:(dim(tab2)[1]-1), tab2[1:(dim(tab2)[1]-1),"median_distance_notLondon"], cex=0.6, pch=16, col="gray30")
		lines(1:(dim(tab2)[1]-1), tab2[1:(dim(tab2)[1]-1),"median_distance_London"], lwd=0.3, col="red")
		points(1:(dim(tab2)[1]-1), tab2[1:(dim(tab2)[1]-1),"median_distance_London"], cex=0.6, pch=16, col="red")
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.45,0), lwd=0.2, tck=-0.03, col.tick="gray30", col.axis="gray30", col="gray30", las=2, at=1:dim(tab2)[1], labels=allLabels)
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.2, tck=-0.03, col.tick="gray30", col.axis="gray30", col="gray30")
		title(ylab="median distance (km)", cex.lab=0.7, mgp=c(1.1,0,0), col.lab="gray30")
		dev.off()
		tmp = matrix(nrow=dim(tab2)[1], ncol=1); tmp[,1] = tab2[,"lineage_exports_London"]
		tmp[,1] = round(tab2[,"lineage_exports_London"]/tab2[,"all_lineage_dispersal"], 3)
		pdf("Figure_B4_NEW.pdf", width=7, height=2.0)
		par(oma=c(0,0,0,0), mar=c(3,2,0,0), lwd=0.2, col="gray30")
		plot(1:(dim(tab2)[1]-1), tmp[1:(dim(tab2)[1]-1),1], ylim=c(0,0.27), type="l", lwd=0.3, col="gray30", axes=F, ann=F)
		points(1:(dim(tab2)[1]-1), tmp[1:(dim(tab2)[1]-1),1], cex=0.6, pch=16, col="gray30")
		axis(side=1, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.45,0), lwd=0.2, tck=-0.03, col.tick="gray30", col.axis="gray30", col="gray30", las=2, at=1:dim(tab2)[1], labels=allLabels)
		axis(side=2, lwd.tick=0.2, cex.axis=0.6, mgp=c(0,0.25,0), lwd=0.2, tck=-0.03, col.tick="gray30", col.axis="gray30", col="gray30")
		title(ylab="ratio London/all exports", cex.lab=0.7, mgp=c(1.1,0,0), col.lab="gray30")
		dev.off()
	}

# 5. Analysing the introduction events in the different UTLA or region polygons

months = c("09","10","11","12","01")
years = c(2020,2020,2020,2020,2021)
maxDays = c(30,31,30,31,13)
days1 = c(); days2 = c()
for (i in 1:length(months))
	{
		if (years[i] == 2020) day_time = 1/366
		if (years[i] == 2021) day_time = 1/365
		days1 = c(days1, paste(c(paste0("0",1:9),10:maxDays[i]),months[i],years[i],sep="-"))
		days2 = c(days2, decimal_date(dmy(paste(c(paste0("0",1:9),10:maxDays[i]),months[i],years[i],sep="-")))+day_time)
	}

	# 5.1. Analysing the daily numbers of introductions per UTLA

UTLA_introductions = matrix(nrow=dim(UTLAs@data)[1], ncol=length(days2)-1)
row.names(UTLA_introductions) = UTLAs@data[,"NUMBER"]
colnames(UTLA_introductions) = days1[1:(length(days1)-1)]
for (i in 1:dim(UTLA_introductions)[1])
	{
		sub = mcc[which((mcc[,"startUTLA"]!=UTLAs@data[i,"NUMBER"])&(mcc[,"endUTLA"]==UTLAs@data[i,"NUMBER"])),]
		for (d in 2:length(days2))
			{
				UTLA_introductions[i,d-1] = length(which((sub["endYear"]>days2[d-1])&(sub["endYear"]<=days2[d])))
			}
	}
write.csv(UTLA_introductions, "UTLA_introductions.csv", quote=F)
if (savingPlots)
	{
		days = decimal_date(dmy(colnames(UTLA_introductions)))
		UTLA_cumulative_introductions = matrix(nrow=dim(UTLA_introductions)[1], ncol=dim(UTLA_introductions)[2])
		UTLA_cumulative_introductions[,1] = UTLA_introductions[,2]
		for (i in 2:dim(UTLA_introductions)[2]) UTLA_cumulative_introductions[,i] = UTLA_cumulative_introductions[,i-1]+UTLA_introductions[,i]
		UTLA_all_introductions = matrix(nrow=dim(UTLA_introductions)[1], ncol=1)
		for (i in 1:dim(UTLA_introductions)[1]) UTLA_all_introductions[i,1] = sum(UTLA_introductions[i,])
		minVal = 0; maxVal = max(UTLA_cumulative_introductions)
		colourScale = rev(gsub("FF","",viridis(101))); plotLegend = TRUE
		cutOffs = c(decimal_date(dmy(c("05-11-2020","01-12-2020","20-12-2020","12-01-2021"))))
		dates = c("05-11-2020","01-12-2020","20-12-2020","12-01-2021")
		pdf("Figure_A1_alternative.pdf", width=3, height=3.5)
		par(oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
		index = which(days<=cutOffs[1])[length(which(days<=cutOffs[1]))]
		indices = (((UTLA_cumulative_introductions[,index]-minVal)/(maxVal-minVal))*100)+1
		cols = colourScale[indices]; cols[which(UTLA_cumulative_introductions[,index]==0)] = NA
		plot(postcodes, border=NA, col="gray90", lwd=0.01)
		plot(UTLAs, border="white", col=cols, add=T, lwd=0.2)
		mtext(dates[1], line=-12, at=225000, cex=0.65, col="gray30")
		dev.off()
		pdf("Figure_C4_alternative.pdf", width=3, height=3.5)
		par(oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
		index = which(days<=cutOffs[4])[length(which(days<=cutOffs[4]))]
		indices = (((UTLA_cumulative_introductions[,index]-minVal)/(maxVal-minVal))*100)+1
		cols = colourScale[indices]; cols[which(UTLA_cumulative_introductions[,index]==0)] = NA
		plot(postcodes, border=NA, col="gray90", lwd=0.01)
		plot(UTLAs, border="white", col=cols, add=T, lwd=0.2)
		if (plotLegend == TRUE)
			{
				rast = raster(matrix(nrow=1, ncol=2)); rast[1] = 1; rast[2] = maxVal
				plot(rast, legend.only=T, add=T, col=rev(colourScale), legend.width=0.5, legend.shrink=0.3, smallplot=c(0.45,0.94,0.114,0.125),
					 legend.args=list(text="", cex=0.6, line=0.3, col="gray30"), horizontal=T,
				     axis.args=list(cex.axis=0.5, lwd=0, lwd.tick=0.2, tck=-0.9, col.axis="gray30", line=0, mgp=c(0,-0.20,0)))
			}
		mtext(dates[4], line=-12, at=225000, cex=0.65, col="gray30")
		dev.off()		
		pdf("Figure_D_alternative.pdf", width=8, height=9)
		par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30"); cexNode = 0.35
		for (h in 1:length(cutOffs))
			{
				index = which(days<=cutOffs[h])[length(which(days<=cutOffs[h]))]
				indices = (((UTLA_cumulative_introductions[,index]-minVal)/(maxVal-minVal))*100)+1
				cols = colourScale[indices]; cols[which(UTLA_cumulative_introductions[,index]==0)] = NA
				plot(postcodes, border=NA, col="gray90", lwd=0.01)
				plot(UTLAs, border="white", col=cols, add=T, lwd=0.1)
				if ((plotLegend == TRUE)&(h==length(cutOffs)))
					{
						rast = raster(matrix(nrow=1, ncol=2)); rast[1] = minVal; rast[2] = maxVal
						plot(rast, legend.only=T, add=T, col=colourScale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.45,0.88,0.098,0.110),
							 legend.args=list(text="", cex=0.8, line=0.3, col="gray30"), horizontal=T,
				  		   	 axis.args=list(cex.axis=0.8, lwd=0, lwd.tick=0.2, tck=-0.9, col.axis="gray30", line=0, mgp=c(0,0.15,0), at=seq(0,600,200)))
					}
				mtext(dates[h], line=-18.5, at=225000, cex=0.7, col="gray30")
			}
		dev.off()
	}

	# 5.2. Analysing the daily numbers of lineage transitions between UTLAs or regions

if (!file.exists("UTLA_allTransitions.rds"))
	{
		UTLA_transition_matrices = list(); dates = days1[1:(length(days1)-1)]
		for (d in 2:length(days2))
			{
				sub = mcc[which((mcc["endYear"]>days2[d-1])&(mcc["endYear"]<=days2[d])),]
				UTLA_transition_matrix = matrix(nrow=dim(UTLAs)[1], ncol=dim(UTLAs)[1])
				row.names(UTLA_transition_matrix) = UTLAs@data[,"NUMBER"]
				colnames(UTLA_transition_matrix) = UTLAs@data[,"NUMBER"]
				for (i in 1:dim(UTLAs)[1])
					{
						for (j in 1:dim(UTLAs)[1])
							{
								UTLA_transition_matrix[i,j] = length(which((sub[,"startUTLA"]==UTLAs@data[i,"NUMBER"])&(sub[,"endUTLA"]==UTLAs@data[j,"NUMBER"])))
							}
					}
				UTLA_transition_matrices[[d-1]] = UTLA_transition_matrix
			}
		rdsObject = list(); rdsObject[[1]] = dates; rdsObject[[2]] = UTLA_transition_matrices
		saveRDS(rdsObject, "UTLA_allTransitions.rds")
	}	else		{
		rdsObject = readRDS("UTLA_allTransitions.rds")
		UTLA_transition_matrices = rdsObject[[2]]
	}
if (!file.exists("Region_allTransitions.rds"))
	{
		region_transition_matrices = list(); dates = days1[1:(length(days1)-1)]
		for (d in 2:length(days2))
			{
				sub = mcc[which((mcc["endYear"]>days2[d-1])&(mcc["endYear"]<=days2[d])),]
				region_transition_matrix = matrix(nrow=dim(regions@data)[1], ncol=dim(regions@data)[1])
				row.names(region_transition_matrix) = regions@data[,"Merged_loc"]
				colnames(region_transition_matrix) = regions@data[,"Merged_loc"]
				for (i in 1:dim(regions@data)[1])
					{
						for (j in 1:dim(regions@data)[1])
							{
								region_transition_matrix[i,j] = length(which((sub[,"startRegion"]==regions@data[i,"Merged_loc"])&(sub[,"endRegion"]==regions@data[j,"Merged_loc"])))
							}
					}
				region_transition_matrices[[d-1]] = region_transition_matrix
			}
		rdsObject = list(); rdsObject[[1]] = dates; rdsObject[[2]] = region_transition_matrices
		saveRDS(rdsObject, "Region_allTransitions.rds")
	}	else		{
		rdsObject = readRDS("Region_allTransitions.rds")
		region_transition_matrices = rdsObject[[2]]
	}
centroidRegions = matrix(nrow=9, ncol=2)
for (i in 1:dim(centroidRegions)[1])
	{
		maxArea = 0; polIndex = 0
		for (j in 1:length(regions@polygons[[i]]@Polygons))
			{
				if (maxArea < regions@polygons[[i]]@Polygons[[j]]@area)
					{
						maxArea = regions@polygons[[i]]@Polygons[[j]]@area; polIndex = j
					}
			}
		pol = regions@polygons[[i]]@Polygons[[polIndex]]
		p = Polygon(pol@coords); ps = Polygons(list(p),1); sps = SpatialPolygons(list(ps))
		centroidRegions[i,] = coordinates(sps)
	}
if (savingPlots)
	{
		dates1 = c("01-09-20","05-11-20","01-12-20","20-12-20"); 
		dates2 = c("05-11-20","01-12-20","20-12-20","12-01-21")
		regions_mod = gBuffer(regions, byid=T, width=0); regions_mod = gSimplify(regions_mod, 10)
		cumulated_transitions = matrix(0, nrow=dim(regions@data)[1], ncol=dim(regions@data)[1])
		for (i in 1:(length(days2)-1)) cumulated_transitions[] = cumulated_transitions[]+region_transition_matrices[[i]][]
		maxNberOfCumulatedLocalTransitions = max(cumulated_transitions)
		diag(cumulated_transitions) = 0; maxNberOfCumulatedTransitions = max(cumulated_transitions)
		cumulated_transitions[cumulated_transitions[]==0] = NA; minNberOfCumulatedTransitions = min(cumulated_transitions, na.rm=T)
		pdf("Figure_B1a_NEW.pdf", width=8, height=2.5)
		par(mfrow=c(1,4), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
		cutOffs = c(decimal_date(dmy(c("05-11-2020","01-12-2020","20-12-2020","12-01-2021"))))
		vMax1 = 0; vMax2 = 0; ratios = c(); cexNode = 0.35
		for (h in 1:length(cutOffs))
			{ 
				if (h > 1)
					{
						starting_day = cutOffs[h-1]
					}	else	{
						starting_day = decimal_date(dmy("01_09-2020"))
					}
				cumulated_transitions = matrix(0, nrow=dim(regions@data)[1], ncol=dim(regions@data)[1])
				index1 = which(days2>starting_day)[1]; index2 = which(days2<=cutOffs[h])[length(which(days2<=cutOffs[h]))]
				for (i in index1:index2) cumulated_transitions[] = cumulated_transitions[]+region_transition_matrices[[i]][]
				cumulated_transitions = cumulated_transitions[1:9,1:9]
				plot(regions_mod, border="white", col=c(rep("gray85",9),"gray95"), lwd=0.4)
				for (i in 1:dim(cumulated_transitions)[1])
					{
						if (cumulated_transitions[i,i] != 0)
							{
								indices = c(1:dim(cumulated_transitions)[1]); indices = indices[which(indices!=i)]
								cexValue = (cumulated_transitions[i,i]/maxNberOfCumulatedLocalTransitions)*40
								if (vMax1 < (cumulated_transitions[i,i]/sum(cumulated_transitions[indices,i])))
									{
										if (sum(cumulated_transitions[indices,i]) != 0)
											{
												vMax1 = (cumulated_transitions[i,i]/sum(cumulated_transitions[indices,i]))
											}
									}
								cexValue = (cumulated_transitions[i,i]/sum(cumulated_transitions[indices,i]))/5.3
								ratios = c(ratios, (cumulated_transitions[i,i]/sum(cumulated_transitions[indices,i])))
								points(centroidRegions[i,1], centroidRegions[i,2], cex=cexValue, pch=16, lwd=0.0, col=rgb(1,0,0,0.3,1))
								points(centroidRegions[i,1], centroidRegions[i,2], cex=cexValue, pch=1, lwd=0.3, col=rgb(1,0,0,1,1))
							}
					}
				for (i in 1:dim(cumulated_transitions)[1])
					{
						for (j in 1:dim(cumulated_transitions)[1])
							{
								if ((i != j)&(cumulated_transitions[i,j] >= 1))
									{
										arrCol = "gray30"; arrPos = 0.45; arrType = "triangle"
										if (vMax2 < cumulated_transitions[i,j]) vMax2 = cumulated_transitions[i,j]
										LWD = (((log(cumulated_transitions[i,j])-4)/(log(maxNberOfCumulatedTransitions)-4))*6)+0.1
										arrow = (0.15*(log(cumulated_transitions[i,j])/log(maxNberOfCumulatedTransitions)))+0.002
										if (cumulated_transitions[i,j] < 3) { arrow = 0; arrCol = NA; arrPos = F; arrType = "none" }
										curvedarrow(cbind(centroidRegions[i,1],centroidRegions[i,2]), cbind(centroidRegions[j,1],centroidRegions[j,2]), arr.length=0.1,
								   					arr.width=arrow, lwd=LWD, lty=1, lcol="gray30", arr.col=arrCol, arr.pos=arrPos, curve=0.13, dr=NA, endhead=F, arr.type=arrType)
									}
							}
					}
				mtext(dates1[h], line=-12.8, at=225000, cex=0.55, col="gray30")
				mtext(dates2[h], line=-13.7, at=225000, cex=0.55, col="gray30")
			}
		dev.off()
		pdf("Figure_B1b_NEW.pdf", width=8, height=2.5)
		par(mfrow=c(1,4), oma=c(0,0,0,0), mar=c(0,0,0,0), lwd=0.2, col="gray30")
		plot(regions_mod, border="white", col=c(rep("gray85",9),"gray95"), lwd=0.4)
		v = 10; LWD = (((log(v)-4)/(log(maxNberOfCumulatedTransitions)-4))*6)+0.1; arrow = (0.15*(log(v)/log(maxNberOfCumulatedTransitions)))+0.002
		curvedarrow(cbind(300000,600000),cbind(400000,600000), arr.length=0.1,arr.width=arrow, lwd=LWD, lty=1, lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=0, dr=NA, endhead=F, arr.type="triangle")
		v = 50; LWD = (((log(v)-4)/(log(maxNberOfCumulatedTransitions)-4))*6)+0.1; arrow = (0.15*(log(v)/log(maxNberOfCumulatedTransitions)))+0.002
		curvedarrow(cbind(300000,500000),cbind(400000,500000), arr.length=0.1,arr.width=arrow, lwd=LWD, lty=1, lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=0, dr=NA, endhead=F, arr.type="triangle")
		v = 100; LWD = (((log(v)-4)/(log(maxNberOfCumulatedTransitions)-4))*6)+0.1; arrow = (0.15*(log(v)/log(maxNberOfCumulatedTransitions)))+0.002
		curvedarrow(cbind(300000,400000),cbind(400000,400000), arr.length=0.1,arr.width=arrow, lwd=LWD, lty=1, lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=0, dr=NA, endhead=F, arr.type="triangle")
		v = 200; LWD = (((log(v)-4)/(log(maxNberOfCumulatedTransitions)-4))*6)+0.1; arrow = (0.15*(log(v)/log(maxNberOfCumulatedTransitions)))+0.002
		curvedarrow(cbind(300000,300000),cbind(400000,300000), arr.length=0.1,arr.width=arrow, lwd=LWD, lty=1, lcol="gray30", arr.col="gray30", arr.pos=0.5, curve=0, dr=NA, endhead=F, arr.type="triangle")
		plot(regions_mod, border="white", col=c(rep("gray85",9),"gray95"), lwd=0.4)
		v = 10; cexValue = v/5.3; points(300000, 600000, cex=cexValue, pch=16, lwd=0.0, col=rgb(1,0,0,0.3,1)); points(300000, 600000, cex=cexValue, pch=1, lwd=0.3, col=rgb(1,0,0,1,1))
		v = 20; cexValue = v/5.3; points(300000, 500000, cex=cexValue, pch=16, lwd=0.0, col=rgb(1,0,0,0.3,1)); points(300000, 500000, cex=cexValue, pch=1, lwd=0.3, col=rgb(1,0,0,1,1))
		v = 50; cexValue = v/5.3; points(300000, 400000, cex=cexValue, pch=16, lwd=0.0, col=rgb(1,0,0,0.3,1)); points(300000, 400000, cex=cexValue, pch=1, lwd=0.3, col=rgb(1,0,0,1,1))
		dev.off()
		pdf("Figure_B2_NEW.pdf", width=8, height=2.0); xMin = 50; xMax = 550
		par(mfrow=c(1,4), oma=c(0,0,0,0), mar=c(2.2,1.5,1,1), lwd=0.2, col="gray30"); cexNode = 0.35
		cutOffs = c(decimal_date(dmy(c("05-11-2020","01-12-2020","20-12-2020","12-01-2021")))); cols1 = list(); cols2 = list()
		cols1[[1]] = rgb(222,67,39,255,maxColorValue=255); cols2[[1]] = rgb(222,67,39,150,maxColorValue=255)
		cols1[[2]] = rgb(150,150,150,255,maxColorValue=255); cols2[[2]] = rgb(150,150,150,150,maxColorValue=255)
		for (h in 1:length(cutOffs))
			{
				if (h > 1)
					{
						starting_day = cutOffs[h-1]
					}	else	{
						starting_day = 0
					}
				sub = mcc[which((mcc[,"midYear"]>starting_day)&(mcc[,"midYear"]<=cutOffs[h])),]
				distances1 = sub[which(sub[,"startRegion"]=="Greater_London"),"geoDist_km"]
				distances2 = sub[which(sub[,"startRegion"]!="Greater_London"),"geoDist_km"]
				distances1 = distances1[distances1[]>=xMin]; distances2 = distances2[distances2[]>=xMin]
				hist(distances1, breaks=seq(xMin,xMax,by=(xMax-xMin)/50), col=cols2[[1]], border="gray30", lwd=0.1, xlim=c(xMin,xMax), axes=F, ann=F)
				hist(distances2, breaks=seq(xMin,xMax,by=(xMax-xMin)/50), col=cols2[[2]], border="gray30", lwd=0.1, add=T)
				axis(side=1, lwd.tick=0.2, cex.axis=0.7, mgp=c(0,0.10,0), lwd=0.2, tck=-0.02, col.tick="gray30", col.axis="gray30", col="gray30", at=seq(0,600,100))
				axis(side=2, lwd.tick=0.2, cex.axis=0.7, mgp=c(0,0.25,0), lwd=0.2, tck=-0.02, col.tick="gray30", col.axis="gray30", col="gray30")
				title(xlab="distance (km)", cex.lab=0.9, mgp=c(1.1,0,0), col.lab="gray30")
				if (h == 1) title(ylab="frequency", cex.lab=0.9, mgp=c(1.5,0,0), col.lab="gray30")
				mtext(dates1[h], line=-1.1, at=300, cex=0.6, col="gray30")
				mtext(dates2[h], line=-2.0, at=300, cex=0.6, col="gray30")
			}
		dev.off()
	}

# 6. Additional analysis of B.1.1.7 introduction events (on a weekly basis)

fourRegions = shapefile("Shapefiles_study_area/England_4selected_regions.shp")
four_selected_regions = gsub(" ","_",fourRegions@data[,"LAD20NM"])
start_end_selected_regions = matrix(nrow=dim(mcc)[1], ncol=2)
colnames(start_end_selected_regions) = c("startSelectedRegion","endSelectedRegion")
for (i in 1:dim(fourRegions@data)[1])
	{
		maxArea = 0; polIndex = 0; # print(i)
		for (j in 1:length(fourRegions@polygons[[i]]@Polygons))
			{
				if (maxArea < fourRegions@polygons[[i]]@Polygons[[j]]@area)
					{
						maxArea = fourRegions@polygons[[i]]@Polygons[[j]]@area; polIndex = j
					}
			}
		pol = fourRegions@polygons[[i]]@Polygons[[polIndex]]
		indices = which(point.in.polygon(mcc[,"startLon"],mcc[,"startLat"],pol@coords[,1],pol@coords[,2])==1)
		if (length(indices) > 0)
			{
				if (sum(!is.na(start_end_selected_regions[indices,"startSelectedRegion"])) != 0) print(c(i,2))
				start_end_selected_regions[indices,"startSelectedRegion"] = gsub(" ","_",fourRegions@data[i,"LAD20NM"])
			}
		indices = which(point.in.polygon(mcc[,"endLon"],mcc[,"endLat"],pol@coords[,1],pol@coords[,2])==1)
		if (length(indices) > 0)
			{
				if (sum(!is.na(start_end_selected_regions[indices,"endSelectedRegion"])) != 0) print(c(i,2))
				start_end_selected_regions[indices,"endSelectedRegion"] = gsub(" ","_",fourRegions@data[i,"LAD20NM"])
			}
	}
start_end_selected_regions[which(is.na(start_end_selected_regions[,1])),1] = "other"
start_end_selected_regions[which(is.na(start_end_selected_regions[,2])),2] = "other"
mcc = mcc[,1:17]; mcc = cbind(mcc, start_end_selected_regions)
tab1 = read.csv("London_stats_MK.csv", head=T)
dates = decimal_date(dmy(gsub("\\/","-",tab1[,"week_start"])))
tab1 = tab1[order(dates),]; dates = dates[order(dates)]
tab2 = matrix(0, nrow=dim(tab1)[1], ncol=length(four_selected_regions))
colnames(tab2) = four_selected_regions
samplingData = read.csv(paste0("Sampling_",analysis_date,".csv"), head=T)
samplingData = samplingData[which(!is.na(samplingData[,"longitude"])),]
mostRecentSamplingDatum = max(samplingData[,"collection_date"])
usingMidBranchTimes = FALSE
for (i in 1:dim(tab1)[1])
	{
		if (i != dim(tab1)[1])
			{
				ending_day = dates[i+1]
			}	else	{
				ending_day = mostRecentSamplingDatum+1/365
			}
		if (usingMidBranchTimes == TRUE) sub1 = mcc[which((mcc[,"midYear"]>=dates[i])&(mcc[,"midYear"]<ending_day)),]
		if (usingMidBranchTimes == FALSE) sub1 = mcc[which((mcc[,"endYear"]>=dates[i])&(mcc[,"endYear"]<ending_day)),]
		if (dim(sub1)[1] > 0)
			{
				for (j in 1:length(four_selected_regions))
					{
						sub2 = sub1[which(sub1[,"endSelectedRegion"]==four_selected_regions[j]),]
						if (dim(sub2)[1] > 0)
							{
								tab2[i,four_selected_regions[j]] = dim(sub2[which(sub2[,"startSelectedRegion"]!=four_selected_regions[j]),])[1]
							}
					}
			}
	}
if (writingFiles)
	{
		allLabels = gsub("\\/","-",gsub("\\/20","\\/",tab1[,"week_start"])); row.names(tab2) = allLabels
		write.csv(tab2, "B117_lineage_introduction_events_in_4_specific_regions.csv", quote=F)
	}

tab1 = read.csv("London_stats_MK.csv", head=T)
dates = decimal_date(dmy(gsub("\\/","-",tab1[,"week_start"])))
tab1 = tab1[order(dates),]; dates = dates[order(dates)]
tab2 = matrix(0, nrow=dim(tab1)[1], ncol=dim(UTLAs@data)[1])
colnames(tab2) = UTLAs@data[,"CODE"]
samplingData = read.csv(paste0("Sampling_",analysis_date,".csv"), head=T)
samplingData = samplingData[which(!is.na(samplingData[,"longitude"])),]
mostRecentSamplingDatum = max(samplingData[,"collection_date"])
usingMidBranchTimes = FALSE
for (i in 1:dim(tab1)[1])
	{
		if (i != dim(tab1)[1])
			{
				ending_day = dates[i+1]
			}	else	{
				ending_day = mostRecentSamplingDatum+1/365
			}
		if (usingMidBranchTimes == TRUE) sub1 = mcc[which((mcc[,"midYear"]>=dates[i])&(mcc[,"midYear"]<ending_day)),]
		if (usingMidBranchTimes == FALSE) sub1 = mcc[which((mcc[,"endYear"]>=dates[i])&(mcc[,"endYear"]<ending_day)),]
		if (dim(sub1)[1] > 0)
			{
				for (j in 1:dim(UTLAs@data)[1])
					{
						sub2 = sub1[which(sub1[,"endUTLA"]==UTLAs@data[j,"NUMBER"]),]
						if (dim(sub2)[1] > 0)
							{
								tab2[i,j] = dim(sub2[which(sub2[,"startUTLA"]!=UTLAs@data[j,"NUMBER"]),])[1]
							}
					}
			}
	}
if (writingFiles)
	{
		allLabels = gsub("\\/","-",gsub("\\/20","\\/",tab1[,"week_start"])); row.names(tab2) = allLabels
		write.csv(tab2, "B117_lineage_introduction_events_in_UTLAs.csv", quote=F)
	}

