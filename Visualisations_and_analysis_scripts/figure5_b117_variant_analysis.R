#MOUGK and SV Scarpino
#Dec. 2020
#UK Variant analysis
#Regression analysis and plots

rm(list = ls())

###########
#libraries#
###########
library(ggplot2)
library(glmulti)
library(zoo)
library(binom)
library(RColorBrewer)
library(forcats)

#########
#Globals#
#########
do_plots <- FALSE
down_weight <- 0.5 #how much to weight the posterior from the previous sampling point
lockdown_start <- as.POSIXct(strptime("2020-11-18", format = "%Y-%m-%d"))
lockdown_stop <- as.POSIXct(strptime("2020-12-02", format = "%Y-%m-%d"))
stop_date <- as.POSIXct(strptime("2020-12-16", format = "%Y-%m-%d")) #this is from the pcr testing and growth rates used in the analysis
new_date <- as.POSIXct(strptime("2020-12-30", format = "%Y-%m-%d"))

######
#Data#
######
start_date_rates <- as.POSIXct(strptime("2020-12-10", format = "%Y-%m-%d"))
stop_date_rates <- start_date_rates + 60*60*24*7
mobility_focal <- "London"

#source("build_reg_data.R")

var_freq <- read.csv('b117_variant_data_F5-1.csv')

#mt_cases <- match(var_freq$utla, names(cases_by_utla_by))
#mt_mob <- match(var_freq$utla, names(mob_by))
#mt_pops <- match(var_freq$utla, pops$Code)
#mt_variant <- match(var_freq$utla, dat.reg$UTLA)

#var_freq$arrival_days_later <- dat.reg$days_since_variant[mt_variant]
#var_freq$arrival <- dat.reg$variant_date[mt_variant]
#var_freq$growth <- dat.reg$growth[mt_variant]
#var_freq$mob <- as.numeric(mob_by)[mt_mob]
#var_freq$cases_pre_nov <- as.numeric(cases_by_utla_by)[mt_cases]  
#var_freq$cases_pre_jan <- as.numeric(cases_by_utla_jan_by)[mt_cases]  
#var_freq$pop <- pops$All.ages[mt_pops]
#var_freq$S_in_nov <- var_freq$pop - var_freq$cases_pre_nov
#var_freq$S_in_jan <- var_freq$pop - var_freq$cases_pre_jan
#var_freq$inc_pre_nov <- var_freq$cases_pre_nov/var_freq$pop
#var_freq$inc_pre_jan <- var_freq$cases_pre_jan/var_freq$pop
#var_freq$propu20 <- dat.reg$Pop19_prop_u20[mt_variant]  
#var_freq$u20 <- dat.reg$Pop19[mt_variant] 

#add_to <- "var_freq"
#source("process_mobility_data.R")

##########
#Analysis#
##########
#variant samps
ps1 <- 0.1 #biased low in first time point
ps2 <- 1

vals <- seq(0.001, 0.999, length.out = 1000)
dens <- c()
locs <- c()
times <- c()
values <- c()
map_change_early_late <- rep(NA, nrow(var_freq))
map_change_late_later <- rep(NA, nrow(var_freq))
map_change_early_later <- rep(NA, nrow(var_freq))
map_change_later_latest <- rep(NA, nrow(var_freq))
map_change_latest_jan21 <- rep(NA, nrow(var_freq))
map_early <- rep(NA, nrow(var_freq))
map_late <- rep(NA, nrow(var_freq))
map_later <- rep(NA, nrow(var_freq))
map_latest <- rep(NA, nrow(var_freq))
map_jan21 <- rep(NA, nrow(var_freq))
for(i in 1:nrow(var_freq)){
  early.posts1.i <- var_freq$Early_count_elephants[i] + ps1
  early.posts2.i <- var_freq$Early_count_all[i] - var_freq$Early_count_elephants[i] + ps2
  early.den.i <- dbeta(vals, early.posts1.i, early.posts2.i)/sum(dbeta(vals,early.posts1.i, early.posts2.i))
  
  late.posts1.i <- var_freq$Late_count_elephants[i] + early.posts1.i*down_weight
  late.posts2.i <- var_freq$Late_count_all[i] - var_freq$Late_count_elephants[i] + early.posts2.i*down_weight
  late.den.i <- dbeta(vals, late.posts1.i, late.posts2.i)/sum(dbeta(vals, late.posts1.i, late.posts2.i))
  
  later.posts1.i <- var_freq$Later_count_elephants[i] + late.posts1.i*down_weight
  later.posts2.i <- var_freq$Later_count_all[i] - var_freq$Later_count_elephants[i] + late.posts2.i*down_weight
  later.den.i <- dbeta(vals, later.posts1.i, later.posts2.i)/sum(dbeta(vals, later.posts1.i, later.posts2.i))
  
  latest.posts1.i <- var_freq$Latest_count_elephants[i] + later.posts1.i*down_weight
  latest.posts2.i <- var_freq$Latest_count_all[i] - var_freq$Latest_count_elephants[i] + later.posts2.i*down_weight
  latest.den.i <- dbeta(vals, latest.posts1.i, latest.posts2.i)/sum(dbeta(vals, latest.posts1.i, latest.posts2.i))
  
  jan.posts1.i <- var_freq$Jan2021_count_elephants[i] + latest.posts1.i*down_weight
  jan.posts2.i <- var_freq$Jan2021_count_all[i] - var_freq$Jan2021_count_elephants[i] + latest.posts2.i*down_weight
  jan.den.i <- dbeta(vals, jan.posts1.i, jan.posts2.i)/sum(dbeta(vals, jan.posts1.i, jan.posts2.i))
  
  if(is.na(early.posts2.i) == TRUE){
    map.early.i <- NA
  }else{
    #map.early.i <- vals[which.max(early.den.i)]
    map.early.i <- weighted.mean(x = vals, w = early.den.i/sum(early.den.i))
  }
  
  if(is.na(late.posts2.i) == TRUE){
    map.late.i <- NA
  }else{
    #map.late.i <- vals[which.max(late.den.i)]
    map.late.i <- weighted.mean(x = vals, w = late.den.i/sum(late.den.i))
  }
  
  if(is.na(later.posts2.i) == TRUE){
    map.later.i <- NA
  }else{
    #map.later.i <- vals[which.max(later.den.i)]
    map.later.i <- weighted.mean(x = vals, w = later.den.i/sum(later.den.i))
  }
  
  if(is.na(latest.posts2.i) == TRUE){
    map.latest.i <- NA
  }else{
    #map.later.i <- vals[which.max(later.den.i)]
    map.latest.i <- weighted.mean(x = vals, w = latest.den.i/sum(latest.den.i))
  }
  
  if(is.na(jan.posts2.i) == TRUE){
    map.jan.i <- NA
  }else{
    #map.later.i <- vals[which.max(later.den.i)]
    map.jan.i <- weighted.mean(x = vals, w = jan.den.i/sum(jan.den.i))
  }
  
  map_change_early_late[i] <-  map.late.i-map.early.i
  map_change_late_later[i] <-  map.later.i-map.late.i
  map_change_early_later[i] <- map.later.i-map.early.i
  map_change_later_latest[i] <- map.latest.i-map.later.i
  map_change_latest_jan21[i] <- map.jan.i-map.latest.i
  map_early[i] <- map.early.i
  map_late[i] <- map.late.i
  map_later[i] <- map.later.i
  map_latest[i] <- map.latest.i
  map_jan21[i] <- map.jan.i
  
  dens <- c(dens, early.den.i, late.den.i, later.den.i, latest.den.i, jan.den.i)
  locs <- c(locs, rep(var_freq$utla[i], length(vals)*5))
  times <- c(times, rep("early", length(vals)), rep("late", length(vals)), rep("later", length(vals)), rep("latest", length(vals)), rep("Jan2021", length(vals)))
  values <- c(values, vals, vals, vals, vals, vals)
}

var_freq$map_change_early_late <- map_change_early_late
var_freq$map_change_late_later <- map_change_late_later
var_freq$map_change_early_later <- map_change_early_later
var_freq$map_change_later_latest <- map_change_later_latest
var_freq$map_change_latest_jan2021 <- map_change_latest_jan21
var_freq$map_early <- map_early 
var_freq$map_late <- map_late
var_freq$map_later <- map_later 
var_freq$map_latest <- map_latest 
var_freq$map_jan2021 <- map_jan21 

pal_col <- c("#4575b4", "#ffffbf", "#d73027")
pal_fill <- paste0(pal_col, "75")

dat.plot <- data.frame(locs, times, dens, values)

if(do_plots == TRUE){
  use_dens <- which(dat.plot$locs %in% sample(unique(dat.plot$locs), 4, replace = FALSE))
  quartz(width = 10, height = 8)
  dat.plot$times <- fct_relevel(dat.plot$times, "Jan2021", "latest", "later", "late", "early")
  ggplot(dat.plot[use_dens,], aes(x = values, y = dens, fill = times)) + geom_density(stat = "identity") + facet_wrap(~locs) + scale_color_brewer(palette = "RdBu") + scale_fill_brewer(palette = "RdBu") + xlab("Variant strain proportion") + ylab("Density") + theme(legend.position = c(0.9,0.8), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 10), axis.title = element_text(colour = "black", size = 16), panel.grid.minor = element_line(colour = "#00000000",linetype = 3), panel.grid.major = element_line(colour = "#00000000", linetype = 3)) + labs(fill = "Sampling point") + scale_y_continuous(limits = c(0, 0.07))
  
  quartz(width = 6, height = 6)
  ggplot(var_freq, aes(x = log(mob), y = log(map_change_late_later))) + geom_point(color = "#4d4d4d", size = 2) + xlab(paste0("Mobility from ", mobility_focal, " (log-scale)")) + ylab("Increase in variant frequency (log-scale)") + theme(legend.position = "none", legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 10), axis.title = element_text(colour = "black", size = 16), panel.grid.minor = element_line(colour = "#00000000",linetype = 3), panel.grid.major = element_line(colour = "#00000000", linetype = 3)) + geom_smooth(method = 'lm', color = "#d6604d")
  
  summary(lm(log(map_change_later_latest) ~ log(mob), data = var_freq))
  
  quartz(width = 6, height = 6)
  ggplot(var_freq, aes(x = log(map_later), y = log(map_latest), color = )) + geom_point() + ylab("Variant frequency (Dec 17th - Dec 30th) log-scale") + xlab("Variant frequency (Dec 2nd - Dec 16th) log-scale") + theme(legend.position = c(0.2,0.8), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 10), axis.title = element_text(colour = "black", size = 16), panel.grid.minor = element_line(colour = "#00000000",linetype = 3), panel.grid.major = element_line(colour = "#00000000", linetype = 3)) + scale_x_continuous() + geom_smooth(method = "lm", color = "#4d4d4d80", lty = 2, alpha = 0.2)+ labs(color = "Mobility from London (log)")
  
  quartz(width = 6, height = 6)
  ggplot(var_freq, aes(x = arrival_days_later, y = log(map_latest))) + geom_point() + ylab("Variant frequency (Dec 17th - Dec 30th) log-scale") + xlab("Estimated time since B.1.1.7 arrival (days)") + theme(legend.position = c(0.2,0.8), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 10), axis.title = element_text(colour = "black", size = 16), panel.grid.minor = element_line(colour = "#00000000",linetype = 3), panel.grid.major = element_line(colour = "#00000000", linetype = 3)) + scale_x_continuous() + geom_smooth(method = "lm", color = "#4d4d4d80", lty = 2, alpha = 0.2)
  
  quartz(width = 6, height = 6)
  ggplot(var_freq, aes(x = log(arrival_days_later), y = log(map_latest))) + geom_point() + ylab("Variant frequency (Dec 17th - Dec 30th) log-scale") + xlab("Estimated time since B.1.1.7 arrival (days) log-scale") + theme(legend.position = c(0.2,0.8), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 10), axis.title = element_text(colour = "black", size = 16), panel.grid.minor = element_line(colour = "#00000000",linetype = 3), panel.grid.major = element_line(colour = "#00000000", linetype = 3)) + scale_x_continuous() + geom_smooth(method = "lm", color = "#4d4d4d80", lty = 2, alpha = 0.2)
  
  quartz(width = 6, height = 6)
  ggplot(var_freq, aes(x = inc_pre_nov, y = log(map_latest), color = log(later_mob))) + geom_point() + xlab("Attack rate pre Nov.") + ylab("Variant frequency (recent)") + theme(legend.position = c(0.8,0.8), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 10), axis.title = element_text(colour = "black", size = 16), panel.grid.minor = element_line(colour = "#00000000",linetype = 3), panel.grid.major = element_line(colour = "#00000000", linetype = 3)) + geom_smooth(method = "lm", col = "#4d4d4d", lty = 2)+ labs(color = "Mobility from London (log)") + scale_y_continuous(limits = c(-6, 0))
  
  quartz(width = 6, height = 6)
  ggplot(var_freq, aes(x = log(later_mob), y = map_later, color = log(inc_pre_nov))) + geom_point() + xlab("Mobility from London") + ylab("Variant frequency (Dec 15th)") + theme(legend.position = c(0.2,0.8), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 10), axis.title = element_text(colour = "black", size = 16), panel.grid.minor = element_line(colour = "#00000000",linetype = 3), panel.grid.major = element_line(colour = "#00000000", linetype = 3)) + geom_smooth(method = "lm", col = "#4d4d4d", lty = 2)+ labs(color = "Mobility from London (log)")
  
  plot_freq <- data.frame(c(var_freq$utla, var_freq$utla, var_freq$utla, var_freq$utla, var_freq$utla), c(var_freq$map_early, var_freq$map_late, var_freq$map_later, var_freq$map_latest, var_freq$map_jan2021), c(var_freq$map_late - var_freq$map_early, var_freq$map_later - var_freq$map_early, var_freq$map_later - var_freq$map_late, var_freq$map_latest - var_freq$map_later, var_freq$map_jan2021 - var_freq$map_latest), c(rep("1-Pre-lockdown", nrow(var_freq)), rep("2-Lockdown", nrow(var_freq)), rep("3-Post-lockdown", nrow(var_freq)),  rep("4-Late-Dec", nrow(var_freq)),rep("5-Early-Jan 21", nrow(var_freq))))
  colnames(plot_freq) <- c("Loc", "Map", "Diff", "Samp")
  
  quartz()
  ggplot(plot_freq, aes(x = Samp, y = Map,  group = Loc, color = Diff)) + geom_point(color = "#00000075") +  geom_line() + scale_color_gradientn(colors = rev(brewer.pal(n = 10, name = "RdBu")), limits = c(-1,1)) + xlab("Sampling point") + ylab("Variant frequency") + theme(legend.position = c(0.2,0.8), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 10), axis.title = element_text(colour = "black", size = 16), panel.grid.minor = element_line(colour = "#00000000",linetype = 3), panel.grid.major = element_line(colour = "#00000000", linetype = 3)) + scale_x_discrete(expand = c(0.05,0.05)) + scale_y_continuous(limits = c(0,1)) + labs(color = "Prop. change") 
  
  glmulti.lm.b117 <-
    glmulti(log(map_change_later_latest) ~ log(later_mob+1e-5) + log(propu20+0.001)+ log(inc_pre_nov), data = var_freq,
            level = 2,               # 2 = pairwise interaction considered
            method = "h",            # h = Exhaustive approach
            crit = "bic",            # bic = BIC as criteria
            confsetsize = 10,        # Keep # best models
            plotty = F, report = F,  # Plot or interim reports
            fitfunction = lm)
  best.glmulti.lm.b117 <- slot(object = glmulti.lm.b117, name = "objects")[[1]]
  summary(best.glmulti.lm.b117)
  
  glmulti.lm.growth <-
    glmulti(growth ~ log(later_mob+1e-5) + log(u20+0.001)+ log(inc_pre_nov), data = var_freq,
            level = 2,               # 2 = pairwise interaction considered
            method = "h",            # h = Exhaustive approach
            crit = "bic",            # bic = BIC as criteria
            confsetsize = 10,        # Keep # best models
            plotty = F, report = F,  # Plot or interim reports
            fitfunction = lm)
  best.glmulti.lm.growth <- slot(object = glmulti.lm.growth, name = "objects")[[1]]
  summary(best.glmulti.lm.growth)
  
  glmulti.lm.b117.later <-
    glmulti(log(map_later) ~ log(mob) + log(map_late)+ log(arrival_days_later+1) + log(inc_pre_nov), data = var_freq,
            level = 2,               # 2 = pairwise interaction considered
            method = "h",            # h = Exhaustive approach
            crit = "bic",            # bic = BIC as criteria
            confsetsize = 10,        # Keep # best models
            plotty = F, report = F,  # Plot or interim reports
            fitfunction = lm)
  best.glmulti.lm.b117.later <- slot(object = glmulti.lm.b117.later, name = "objects")[[1]]
  summary(best.glmulti.lm.b117.later)
  
}


