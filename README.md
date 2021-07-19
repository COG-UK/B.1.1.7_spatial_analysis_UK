# B.1.1.7_spatial_analysis_UK
[![DOI](https://zenodo.org/badge/384450241.svg)](https://zenodo.org/badge/latestdoi/384450241)

## Commonly used acronyms

UTLA: Upper tier local authority 

LTLA: Lower tier local authority 

LAD: Local administrative district


## Files in Data
- aggregate_Nsgtf.csv: Numbers of non-SGTF cases by UTLA and date
- aggregate_sgtf.csv: Numbers of SGTF cases by UTLA and date
- b117_variant_data_5.csv: Counts of B.1.1.7 by UTLAs, used for generating figure 5
- cases_ltla_2021-05-10.csv: Number of cases by LTLA by date up until 2021-05-10
- earliest_B.1.1.7_sequences.csv : Contains the earliest B.1.1.7 sequence date for each UTLA, produced by calculate_B.1.1.7_earliest_dates.py
- estimated_growth_rates.csv: Estimated growth rate of SGTF and non-SGTF cases by LTLA
- merged_data_4c.csv: data containing growth rates of SGTF/Non-SGTF cases by LTLA, number of movements and introductions from the continuous phylogeographic analysis
- SGTF_growths_4b.csv: Contains growth rates by SGTF status and number of exports from the continuous phylogeography
- weekly_B.1.1.7_sequences.csv: Contains the number of B.1.1.7 and non-B.1.1.7 sequences within each area in each week, produced by extract_B.1.1.7_sequence_trajectories.py

In geographic_files:
- gadm36_GBR_2.json and gadm36_GBR_3.json : map files for adm2 and adm3 in the UK, taken from the Global Adminstrative Database (https://gadm.org/)
- gadm36_IRL_1.json : map for the Republic of Ireland, taken from the Global Administrative Database (https://gadm.org/)
- NI_counties.geojson : geographical data for Northern Irish counties (from Open Data NI https://www.opendatani.gov.uk/dataset?tags=Counties)
- UTLA_administrative_areas.* : Shapefiles for UTLAs (from data.gov.uk)
- England_postcoode_districts.* : Shapefiles for postcode districts in England (Addy 2017, https://datashare.ed.ac.uk/handle/10283/2597)
- England_defined_regions1.* : Shapefiles where all but a few specific of the most urban (Birmingham, Greater Manchester) and most rural (Cumbria, Norfolk and Cornwall), as well as Greater London and Kent are merged together.
- England_defined_regions2.* : The same as above, but also including Northumberland and Bristol
- England_4selected_regions.* : Shapefiles only containing four regions of England (Newcastle upon Tyne, Birmingham, Manchester, Liverpool)


In lookup_tables:

- adm2_aggregation.csv : details of aggregating official adm2s to more commonly used ones to avoid false missing data
- LAD_UTLA_adm2.csv : lookup table for comparing different levels of UK geography, namely Local Adminstrative Areas, Upper Tier Local Authorities and Administrative level 2 areas.
- nuts_to_adm2.tsv : lookup table for NUTS1 regions to adminstrative level 2 areas
- postcode_to_adm2.tsv : lookup table for postcode region (ie first half of a UK postcode) to adm2


## Scripts in Visualisations_and_analysis_scripts

- analysisgrowthrate.R: Used to get growth rates per LTLA per date
- analysis4c.R: Used to generate figure 4c
- calculate_B.1.1.7_earliest_dates.py : Used to calculate the earliest B.1.1.7 sequence from each UTLA. Run on 20210122_Genomes_UTLA_corrected_clustered_adm2.csv
- Case corelations.ipynb : used to analyse and create figures of correlations between sequence numbers, case counts and introductions
- extract_B.1.1.7_sequence_trajectories.py : Used to calculate the number of B.1.1.7 sequences in each area through time. Run on 20210122_Genomes_UTLA_corrected_clustered_adm2.csv
- figure5_b117_variant_analysis.R: Used to plot Figure 5. Run onE arly_and_late_and_later_frequencies_updated.csv
- Figure_4b.R: runs on SGTF_growths_F4b.csv
- Figures and analysis.ipynb : used for some descriptive figures
- mobility-uk-b117.ipynb : Used to plot figure 1d and 1e
- plot_B.1.1.7_earliest_dates.R : Used to plot Figure 1B and calculate the number of new areas with B.1.1.7 detected each week. Run on earliest_B.1.1.7_sequences.csv
- plot_B.1.1.7_sequences_trajectories.R : Used to plot Figure S10. Run on weekly_B.1.1.7_sequences.csv
- Tree_visualisation.ipynb : Used to generate Figure S3

Also contains *Phylogeographic_visualisations_and_analyses* which contains analysis files for the continuous phylogeographic analysis.

## Data availability and access

Aggregated epidemiological data used in this study are available from https://coronavirus.data.gov.uk/details/download. SARS-CoV-2 infection survey data are available via the Office of National Statistics (ONS) and available from https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19infectionsurveydata. Raw epidemiological SARS-CoV-2 line list data are available via Public Health England (PHE) and aggregated statistics are available via Github. All genomes, phylogenetic trees, basic metadata are available from the COG-UK consortium website (https://www.cogconsortium.uk/data). The O2 aggregated, anonymised mobile data insights dataset is not publicly available owing to stringent licensing agreements. Information on the process of requesting access to the O2 aggregated mobile data insights dataset is available o2@businesso2.co.uk. The Google COVID-19 Aggregated Mobility Research Dataset is not publicly available owing to stringent licensing agreements. Information on the process of requesting access to the Google mobility data is available from sadilekadam@google.com.

Genome sequences used were generated by the COG-UK consortium (https://www.cogconsortium.uk/). Data linking COG-IDs to location have been removed to protect privacy, however if you require this data please visit https://www.cogconsortium.uk/contact/ for information on accessing consortium-only data.
