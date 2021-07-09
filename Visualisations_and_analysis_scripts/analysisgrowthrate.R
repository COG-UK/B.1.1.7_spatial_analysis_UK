
################################################################################
# READ IN AGGREGATED SGTF AND NON-SGTF CASE COUNT DATA
# READ IN TOTAL OBSERVED CASES BY LTLA AND DATE (NOT SPLIT INTO SGTF AND NON_SGTF)
# SCALEUP THE SGTF AND NON-SGTF CASES
# CALCULATE EXCESS GROWTH RATE (FROM SCALED UP SGTF AND NON-SGTF COUNTS)
################################################################################

rm(list = ls())

library("lubridate")
library(plyr)
library("gridExtra")
library(ggplot2)
library(scales)
library(ggmuller)
library(zoo)
library(tidyverse)
library(lme4)
library(data.table) # to get dates between two dates
library(reshape2) # to convert wide to long format (B117 intro dataset)

## SET YOUR WORKING DIRECTORY WHERE YOU HAVE THE DATA STORED

## GROWTH RATE FUNCTION
growthrate <- function(aggregate_dat){
  diff.df <- aggregate_dat %>% rename(Date = specimen_date.x,
                                      County = LTLA_code,
                                      New = count) %>%
    select(Date, County, New) %>%
    arrange(County)
  
  dates <- unique(diff.df$Date)
  dates <- dates[order(dates)]
  
  # why we need random intercept and slope (different starting points and slope)
  ggplot(data = diff.df[1:1000,], aes(x = Date, y = New, fill = County))+
    geom_line(aes(color = County))
  
  doubling_prov <- list()
  doubling_fixed <- list()
  for(i in 6:(length(dates))){
    use.i <- which(diff.df$Date >= dates[i-5] & diff.df$Date <= dates[i])
    
    data.i <- diff.df[use.i, ]
    
    data.i$DATE <- as.numeric(data.i$Date - dates[i-5], unit = "days")
    # model for random intercept and slope
    mod3.i <- try(lmer(data = data.i, log(New + 1) ~ DATE + (DATE |County)), silent = TRUE)
    
    if(is(mod3.i)[1] == "try-error"){
      fixed.i <- NA
      doubling.i <- NA
      mob.i <- NA
      names.i <- NA
    }else{
      fixed.i <-   fixef(mod3.i)["DATE"]
      doubling.i <- ranef(mod3.i)$County$DATE+fixed.i
      names.i <- rownames(ranef(mod3.i)$County)
    }
    
    doubling_prov[[i-5]] <- doubling.i
    names(doubling_prov[[i-5]]) <- names.i
    doubling_fixed[[i-5]] <- fixed.i
  }
  
  rates <- unlist(lapply(doubling_prov, function(x) x))
  ltla_code <- unlist(lapply(doubling_prov, function(x) names(x)))
  date <- rep(dates[6:(length(dates))], times = unlist(lapply(doubling_prov, function(x) length(x))))
  growthrate_df <- data.frame(rates, date, ltla_code)
  head(growthrate_df)
  
  return(growthrate_df)
}


## SET YOUR WORKING DIRECTORY TO WHERE THE DATASETS ARE STORED
## READ IN OBSERVED AGGREGATED SGTF AND NON-SGTF COUNT DATA BY LTLA AND DATE
aggregate_sgtf <- read_csv("aggregate_sgtf.csv")
aggregate_Nsgtf <- read_csv("aggregate_Nsgtf.csv")

aggregate <- merge(aggregate_sgtf, aggregate_Nsgtf, by = c("specimen_date.x", "LTLA_code"), all = TRUE) %>%
             select("specimen_date.x", "LTLA_code", "count.x", "count.y") %>%
             rename(date = specimen_date.x, ltla_code = LTLA_code,
                    sgtf_count = count.x, Nsgtf_count = count.y) %>%
             filter(date>='2020-10-5' & date<='2021-01-16')
aggregate$sample_size <- rowSums(aggregate[, c("sgtf_count", "Nsgtf_count")], na.rm = TRUE)

## READ IN DAILY CASES RECORDED (NOT DIVIDED INTO SGTF AND NON-SGTF)
df_cases <- read_csv("cases_ltla_2021-05-10.csv") %>%
  mutate(date = as.Date(date)) %>%
  filter(date>'2020-09-19') %>%
  rename(ltla_code = areaCode)

df <- merge(aggregate, df_cases, by = c("date", "ltla_code"))

## CREATE "NEW" SCALED UP COUNTS OF SGTF AND NON-SGTF CASES
df_new <- df %>%
            mutate(prop_sampled = sample_size/newCasesBySpecimenDate,
                   count_new_sgtf = ifelse(prop_sampled<1, round(sgtf_count/prop_sampled), sgtf_count),
                   count_new_Nsgtf = ifelse(prop_sampled<1, round(Nsgtf_count/prop_sampled), Nsgtf_count))

# RENAME THE "NEW" COUNTS OF SGTF AND NON_SGTF CASES TO CALCULATE THE GROWTH RATE
aggregate_sgtf_new <- df_new %>% rename(count = count_new_sgtf, specimen_date.x = date, LTLA_code = ltla_code)
aggregate_Nsgtf_new <- df_new %>% rename(count = count_new_Nsgtf, specimen_date.x = date, LTLA_code = ltla_code)

# FOR NEW SCALED UP SGTF AND NON-SGTF GROWTH RATES
growthrate_sgft_df_new <- growthrate(aggregate_sgtf_new) %>% rename(rates_sgtf = rates)
growthrate_Nsgft_df_new <- growthrate(aggregate_Nsgtf_new) %>% rename(rates_Nsgtf = rates)

df_final <- merge(growthrate_sgft_df_new, growthrate_Nsgft_df_new, by = c("date", "ltla_code")) %>%
  mutate(dateltlacode = paste0(date, ltla_code),
         excess_growthrate = rates_sgtf - rates_Nsgtf)

write.csv(df_final, "estimated_growth_rates.csv")

