#Project: NAMASTE National Univariate Exploratory Analysis
#Script: Script 1 - Data Import and Formatting

#Author: Grant McKown, Jackson Estuarine Laboratory, University of New Hampshire
#Contact: james.mckown@unh.edu

#Project Description: The NERR system collects vegetation data along permanent transects and plots
#   in salt marshes throughout the United States. Vegetation data was synthesized, formatted, and 
#   filtered (based on field monitoring requirements). Reserve, site, transect, and plot-level characteristics
#   were also collected including hydrology, landscape geomorphology, elevation, and climate. The 
#   series of scripts conducts an exploratory analysis of different vegetation response factors to 
#   a series of explanatory and grouping factors. 

#Script Description: Script 1 imports the vegetation datasets compiled across all of the NERRS. Any
#   data formatting required for the regression analysis will be conducted here. Formatted datasets
#   are exported. 


#Chapter 1: Import package library

#Data Analysis

library(dplyr)
library(tidyr)
library(lme4)
library(afex)
library(rstatix)

#Chapter 1: Import and Explore the Vegetation Datasets

#Page 1: Import and format the 'veg_and_expl' dataset
 
  # The 'veg_and_expl' dataset contains the calculated vegetation characteristics and the site characteristics

veg <- read.csv("Input Data\\veg_and_expl.csv")

glimpse(veg)


veg_format <- veg %>%
  select(Reserve, SiteID, TransectID, PlotID, Vegetation_Zone, Year,
         Unique_ID, Years_sinceStart, Total.unvegetated, Total.live.veg,
         F.Freshwater, EIR, Invasive_Cover, Richness, SWdiv, Salt_to_Total, 
         Geomorphology, Tidal.Range, Salinity.category) %>%
  rename(
    abiotic_cover = Total.unvegetated,
    live_cover = Total.live.veg,
    freshwater = F.Freshwater,
    invasive_cover = Invasive_Cover,
    salt_ratio = Salt_to_Total,
    tidal_range = Tidal.Range,
    salinity = Salinity.category)

glimpse(veg_format)

write.csv(veg_format,
          "Formatted Datasets\\Veg Plot Dataframe Formatted.csv")
