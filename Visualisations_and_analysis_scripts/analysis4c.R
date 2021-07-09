
################################################################################
# READ IN MERGED FINAL DATASET
# SMOOTH THE VARIABLES USED IN REGRESSION
# RUN REGRESSION ON 6 WINDOW ROLLING PERIODS
# PLOT THE CONDITIONAL R SQUARED
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
library(wesanderson)
library(MuMIn) # for getting R2 and p value of lmer

## SET YOUR WORKING DIRECTORY WHERE YOU HAVE THE DATA STORED

## READ IN THE COMBINED DATASET
df <- read_csv("merged_data_4c.csv")

## FUNCTION FOR CALCULATING AND PLOTTING R2 FOR 6 DAYS (ROLLING WINDOWS)- FIG 4C
mod_fun <- function(data_input){
  df <- data_input %>% filter(date>='2020-11-15' & date<= '2020-12-20')
  
  dates <- unique(df$date)
  dates <- dates[order(dates)]
  cR2_full <- NA
  cR2_wo_intro <- NA
  cR2_diff <- NA
  coef_mob <- NA
  
  for(i in 6:(length(dates))){
    use.i <- which(df$date >= dates[i-5] & df$date <= dates[i])
    
    data.i <- df[use.i, ]
    
    data.i$DATE <- as.numeric(data.i$date - dates[i-5], unit = "days")
    # Two models for random intercept and slope
    # Full model with introductions, and other without
    mod_full.i <- try(lmer(data = data.i, excess_growthrate ~ inf_mov_log + 
                             b117_intro_count_daily + (1|ltla_code)), silent = TRUE)
    r2_mod_full.i <- r.squaredGLMM(mod_full.i)
    mod_wo_intro.i <- try(lmer(data = data.i, excess_growthrate ~inf_mov_log + 
                                 (1|ltla_code)), silent = TRUE)
    r2_mod_wo_intro.i <- r.squaredGLMM(mod_wo_intro.i)  
    
    if(is(mod_full.i)[1] == "try-error"){
      cR2_full.i <- NA
      cR2_wo_intro.i <- NA
      cR2_diff.i <- NA
      coef_mob.i <- NA
      
    }else{
      cR2_full.i <- r2_mod_full.i[2]
      cR2_wo_intro.i <- r2_mod_wo_intro.i[2]
      cR2_diff.i <- r2_mod_full.i[2] - r2_mod_wo_intro.i[2]
      coef_mob.i <- fixef(mod_full.i)["inf_mov_log"]
      
    }
    
    cR2_full[i-5] <- cR2_full.i
    cR2_wo_intro[i-5] <- cR2_wo_intro.i
    cR2_diff[i-5] <- cR2_diff.i
    coef_mob[i-5] <- coef_mob.i
  }
  
  date <- rep(dates[6:(length(dates))], times = unlist(lapply(cR2_full, function(x) length(x))))
  rollingwindow_df <- data.frame(date, cR2_full, cR2_wo_intro, cR2_diff, coef_mob)
  
  
  rollingwindow_df_melt <- melt(rollingwindow_df[,c("date", "cR2_wo_intro")], id.vars = "date")
  
  a <-ggplot(data = rollingwindow_df_melt) +
    geom_point(aes(x = date, y = value, color = variable), show.legend = FALSE, alpha = 0.6) +
    scale_colour_grey() +
    geom_rect(data = data.frame(xmin = as.Date("2020-11-15"),
                                xmax = as.Date("2020-12-01"),
                                ymin = -Inf,
                                ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "grey", alpha = 0.5) + 
    ylim(0, 1) + 
    scale_x_date(limits = c(as.Date('2020-11-15'), as.Date('2020-12-20'))) + 
    theme_bw() +
    labs(x = "Date", y = expression(paste("Conditional ", R^2))) +
    theme(text = element_text(size=15))
  
  plot(a)
  
}

## SMOOTH VALUES 
df_smooth <- df %>%
  dplyr::arrange(desc(ltla_code)) %>%
  dplyr::group_by(ltla_code) %>%
  dplyr::arrange(date) %>%
  dplyr::mutate(rates_sgtf6 = zoo::rollmean(rates_sgtf, k = 6, fill = NA),
                excess_growthrate6 = zoo::rollmean(excess_growthrate, k = 6, fill = NA),
                rates_Nsgtf6 = zoo::rollmean(rates_Nsgtf, k = 6, fill = NA),
                inf_mov_log6 = zoo::rollmean(inf_mov_log, k = 6, fill = NA),
                b117_intro_count_daily6 = zoo::rollmean(b117_intro_count_daily, k = 6, fill = NA)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(ltla_code, date)  


## USE DF_SMOOTH HERE TO USE THE PLOTTING FUNCTION (RENAME SMOOTH VARIBLES TO USE IN THE FUNCTION)
Fig4c <- mod_fun(data_input = df_smooth %>% 
                              select(date, ltla_code, rates_sgtf6, rates_Nsgtf6,
                                     excess_growthrate6, inf_mov_log6, b117_intro_count_daily6,
                                     LTLA_name) %>%
                              filter(date>='2020-11-15' & date<= '2020-12-20') %>%
                              rename(excess_growthrate = excess_growthrate6,
                                     rates_sgtf = rates_sgtf6,
                                     rates_Nsgtf = rates_Nsgtf6,
                                     inf_mov_log = inf_mov_log6,
                                     b117_intro_count_daily = b117_intro_count_daily6))


## SAVE THE PLOT
pdf(file = "Fig4c.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

grid.arrange(Fig4c, ncol=2, nrow = 4)

dev.off()


