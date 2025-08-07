#Project: NAMASTE National Univariate Exploratory Analysis
#Script: Script 11 - Extract Effect Size of Site Regressions

#Author: Grant McKown, Jackson Estuarine Laboratory, University of New Hampshire
#Contact: james.mckown@unh.edu

#Project Description: The NERR system collects vegetation data along permanent transects and plots
#   in salt marshes throughout the United States. Vegetation data was synthesized, formatted, and 
#   filtered (based on field monitoring requirements). Reserve, site, transect, and plot-level characteristics
#   were also collected including hydrology, landscape geomorphology, elevation, and climate. The 
#   series of scripts conducts an exploratory analysis of different vegetation response factors to 
#   a series of explanatory and grouping factors. 

#Script Description: Script 11 rehashes the code from Script 3 but extracts the effect
# size and standard error of effect size for the site-level and site - zone level
# regressions for the meta-analysis. Analysis only focuses on EMI. 


#Chapter 1: Import package library ---------------------------------------------

#Data Analysis Packages

library(dplyr)
library(tidyr)
library(lme4)
library(afex)
library(rstatix)
library(broom)
library(broom.mixed)
library(purrr)
library(ggResidpanel)

#Graphing and Data Visualization Packages

library(ggplot2)
library(patchwork)


#Chapter 2: Format the vegetation dataset for broom - map functions to create numerous models ------------------------

#Task 1: Import the national plot dataframe
veg <- read.csv("Formatted Datasets\\Veg Plot Dataframe Formatted.csv") %>%
  select(-X, -invasive_cover) %>%
  select(Reserve:Year, Region:SLR_Rate_19yrs, abiotic_cover:salt_ratio)

glimpse(veg)


#Task 2: Format the national plot dataframe by renaming vegetation metrics and changing
#        the format of the dataset from Wide --> Long
veg_format <- veg %>%
  gather(abiotic_cover:salt_ratio,
         key = "Metric",
         value = "Value") %>%
  filter(!is.na(Value)) %>%
  mutate(Metric = ifelse(Metric == "abiotic_cover", "Abiotic Cover", 
                         ifelse(Metric == "live_cover", "Live Cover",
                                ifelse(Metric == "halophyte", "Halophyte Cover",
                                       ifelse(Metric == "SWdiv", "Shannon-Weiner Diversity", 
                                              ifelse(Metric == "freshwater", "Freshwater Cover",
                                                     ifelse(Metric == "salt_ratio", "Salt Ratio", 
                                                            ifelse(Metric == "invasive_cover", "Invasive Cover",
                                                                   Metric)))))))) %>%
  mutate(Metric = factor(Metric, levels = c("Abiotic Cover", "Live Cover", "Halophyte Cover",
                                            "Freshwater Cover", "Invasive Cover", "EMI", 
                                            "Salt Ratio", "Shannon-Weiner Diversity", "Richness"))) %>%
  filter(Metric == "EMI")

glimpse(veg_format)



#Chapter 3: Effect Size of Zone - Site Regressions----------------------------

#Task 1: Effect Size for sites with multiple vegetation zones

# This is accomplished by nesting the data essentially by Site, Vegetation Zone, and Metric
effect_multi_zones <- veg_format %>%
  group_by(Reserve, SiteID) %>%
  mutate(Veg_Zone_Count = length(unique(Vegetation_Zone)),
         SampleSize = n()) %>%
  ungroup() %>%
  filter(Veg_Zone_Count > 1) %>%
  group_by(Reserve, SiteID, Metric, SampleSize, Region, salinity, 
           tidal_range, Geomorphology, SLR_Rate_19yrs) %>%
  nest() %>%
  mutate(regression = map(.x = data,
                          ~lm(Value ~ Year * Vegetation_Zone,
                                    data = .x) %>%
                            tidy())) %>%
  unnest(regression) %>%
  select(-data, -statistic, -p.value) %>%
  mutate(
    term = ifelse(term == "(Intercept)",
                  "VegetationZone_Low", term),
    Regression_Type = "ANCOVA")


glimpse(effect_multi_zones)

write.csv(effect_multi_zones,
          "Output Stats\\Site - Zone Compilation - Veg Zone 2 or more - Effect Size.csv")


#Task 2: Calculate the p-values of the regressions of sites with only 1 site using simple linear 
# regressions

# This is accomplished by nesting the data essentially by Site, Vegetation Zone, and Metric
regression_single_zone <- veg_format %>%
  group_by(Reserve, SiteID) %>%
  mutate(Veg_Zone_Count = length(unique(Vegetation_Zone)),
         SampleSize = n()) %>%
  ungroup() %>%
  filter(Veg_Zone_Count == 1) %>%
  group_by(Reserve, SiteID, Metric, SampleSize, Region, salinity, 
           tidal_range, Geomorphology, SLR_Rate_19yrs, 
           Vegetation_Zone) %>%
  nest() %>%
  mutate(regression = map(.x = data,
                          ~lm(Value ~ Year,
                              data = .x) %>%
                            tidy())) %>%
  unnest(regression) %>%
  select(-data, -statistic, -p.value) %>%
  filter(term != "(Intercept)") %>%
  mutate(Regression_Type = "1 Veg Zone - Linear Regression")


glimpse(regression_single_zone)

write.csv(regression_single_zone,
          "Output Stats\\Site - Zone Compilation - Veg Zone 1 - Effect Size.csv")


#Chapter 4: Size-level regressions -----------------------------------------


#Task 1: Formatting and organizing the dataset for regression slopes including:

# This is accomplished by nesting the data essentially by Site and Metric
effect_site <- veg_format %>%
  group_by(Reserve, SiteID) %>%
  mutate(SampleSize = n()) %>%
  ungroup() %>%
  group_by(Reserve, SiteID, Metric, SampleSize, Region, 
           salinity, tidal_range, Geomorphology, SLR_Rate_19yrs) %>%
  nest() %>%
  mutate(regression = map(.x = data,
                          ~lm(Value ~ Year,
                              data = .x) %>%
                            tidy())) %>%
  unnest(regression) %>%
  select(-data, -statistic, -p.value) %>%
  filter(term != "(Intercept)") %>%
  mutate(Regression_Type = "Site - Level Linear Regression")


glimpse(effect_site)

write.csv(effect_site,
          "Output Stats\\Site Regressions - Effect Size.csv")

