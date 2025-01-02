#Project: NAMASTE National Univariate Exploratory Analysis
#Script: Script 6 - Mixed Linear w/out Time Models to Investigate Changes in the Vegetation Community by Zone

#Author: Grant McKown, Jackson Estuarine Laboratory, University of New Hampshire
#Contact: james.mckown@unh.edu

#Project Description: The NERR system collects vegetation data along permanent transects and plots
#   in salt marshes throughout the United States. Vegetation data was synthesized, formatted, and 
#   filtered (based on field monitoring requirements). Reserve, site, transect, and plot-level characteristics
#   were also collected including hydrology, landscape geomorphology, elevation, and climate. The 
#   series of scripts conducts an exploratory analysis of different vegetation response factors to 
#   a series of explanatory and grouping factors. 

#Script Description: Script 5 formats the site-level vegetation data, calculates the percent change across plot, 
#                   summarizes percent change across site for each vegetation zone. Lastly, it calculates 
#                   linear regression models by Vegetation Zone and other Site Characteristics

#Note - Script 5 is very similar to Script 3, since it accomplishes very similar analyses and research goals. 

#-----------------------------------
#Chapter 1: Import package library
#-----------------------------------


#Data Analysis Packages

library(dplyr)
library(tidyr)
library(lme4)
library(afex)
library(rstatix)
library(broom)
library(broom.mixed)
library(MuMIn)
library(purrr)
library(ggResidpanel)
library(ggeffects)

#Graphing and Data Visualization Packages

library(ggplot2)
library(patchwork)

#---------------------------------------------------------------------------------------------
#Chapter 2: Import the Formatted Vegetation Plot Dataset
#---------------------------------------------------------------------------------------------

veg <- read.csv("Formatted Datasets\\Veg Plot Dataframe Formatted.csv") %>%
  select(-X)

glimpse(veg)


#-------------------------------------------------------------------------------------------
#Chapter 3: Calculate the Annualized Change Based on the Last and First Year
#-----------------------------------------------------------------------------------------

#For this analysis, we are not concerned about the inter-annual changes in vegetation metrics
# over time. We are just concerned about the change between the first year in the dataset
# and last year monitoring. Percent change provided very wacky results, so changed to 
# calculating the annualized change over time. 

#Task 1: Calculate Percent Change between first and last year of each plot

veg_change <- veg %>%
  group_by(Reserve, SiteID, TransectID, PlotID, 
           Vegetation_Zone, Region, Geomorphology, tidal_range, salinity) %>%
  arrange(Year) %>%
  summarise(across(abiotic_cover:salt_ratio,
                   list(
                    change = ~(last(.) - first(.))/(max(Year) - min(Year))))) %>%
  ungroup() %>%
  mutate(across(abiotic_cover_change:salt_ratio_change,
                ~round(., 2)))

write.csv(veg_change,
          "Formatted Datasets\\Annualized Change of Vegetation Metrics at Plot Level.csv")


#Task 2: Calculate the mean change across vegetation metrics for each Site in the dataset

# This is similar to how the original vegetation datasets were used for the Linear Mixed Time Models

veg_change_site <- veg_change %>%
  group_by(Reserve, SiteID, 
           Vegetation_Zone, Region, Geomorphology, tidal_range, salinity) %>%
  summarise(across(abiotic_cover_change:salt_ratio_change,
            list(~mean(., na.rm = TRUE)))) %>%
  ungroup() %>%
  mutate(across(abiotic_cover_change_1:salt_ratio_change_1,
                ~round(., 3))) %>%
  rename(
    abiotic_cover = abiotic_cover_change_1,
    live_cover = live_cover_change_1,
    halophyte = halophyte_change_1,
    freshwater = freshwater_change_1,
    invasive_cover = invasive_cover_change_1,
    salt_ratio = salt_ratio_change_1,
    richness = Richness_change_1,
    SWdiv = SWdiv_change_1,
    EIR = EIR_change_1)

write.csv(veg_change,
          "Formatted Datasets\\Annualized Change of Vegetation Metrics at Site Level.csv")            


#-------------------------------------------------------------------------
#Chapter 3: Calculate Mixed Linear Models without Time as a Fixed Factor
#-------------------------------------------------------------------------

#First, the veg_change_site dataframe needs to be gathered twice - by Vegetation Metric and Site Characteristic

#Second, mixed linear models will be conducted with a mix of purr and broom as in Scripts 3 & 4

# The mixed linear models will test whether vegetation zone and site characteristic
# had an effect on the rate of change of each vegetation metric

#Task 1: Format the dataset from Wide --> Long for Vegetation Metrics and Site Characteristics

veg_change_format <- veg_change_site %>%
  gather(Region:salinity, 
         key = "Site_variable", 
         value = "Category") %>%
  gather(abiotic_cover:salt_ratio,
         key = "Metric",
         value = "Value")


#Task 2: Linear mixed models using purr and broom.mixed

# The overall model is: Vegetation Metric ~ Vegetation Zone * Site Characteristic + (1|SiteID\TransectID)

regression_veg_change <- veg_change_format %>%
  group_by(Metric, Site_variable) %>%
  mutate(SampleSize = length(Reserve)) %>%
  ungroup() %>%
  nest(data = c(Reserve:Vegetation_Zone, Value, Category)) %>%
  mutate(regression = map(.x = data,
                          ~lmer(Value ~ Vegetation_Zone * Category + (1|Reserve/SiteID),
                                data = .x) %>%
                            anova() %>%
                            select(-NumDF, -DenDF) %>%
                            tidy())) %>%
  
  unnest(regression) %>%
  select(-data) %>%
  mutate(across(sumsq:p.value, ~round(., 3))) %>%
  ungroup() 

glimpse(regression_veg_change)