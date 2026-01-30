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

#Script Description: Script 1 imports the several vegetation datasets compiled across all of the NERRS. 
#   Any data formatting required for the regression analysis will be conducted here. Formatted datasets
#   are exported as CSV. First, the national plot dataset is formatted and then averaged by (1) Vegetation 
#   Zone within each Site and (2) within each Site. Second, vegetation slopes by plot datasets were 
#   averaged by (1) Vegetation Zone within each Site and (2) Vegetation Zone by other NERR staff. In 
#   this script, they are then formatted for future analyses. Skimr data summary are conducted at the end of script. 



# Chapter 1: Import package library --------------------------------------------


#Data Analysis

library(dplyr)
library(tidyverse)
library(skimr)


# Chapter 2: Import and Format the National Plot Dataframe --------------------------

#Task 1: Import and format the 'veg_and_expl' dataset (aka the National Plot Dataset)
 
  # The 'veg_and_expl' dataset contains the calculated vegetation characteristics 
  # and the site characteristics

veg <- read.csv("Input Data\\Finalized National Plot Dataframe July 2025.csv")

glimpse(veg)

#Task 2: Format the National Plot Data Frame

veg_format <- veg %>%
  #Select the relevant site identifiers, site characteristics, and vegetation metrics for future analyses
  select(Reserve, SiteID, TransectID, PlotID, Vegetation_Zone, Year,
         Longitude, Latitude, 
         Total.unvegetated, Total.live.veg, H.Halophyte,
         F.Freshwater, EMI, Invasive_Cover, Richness, SWdiv, Salt_to_Total, 
         NERR.Region, Geomorphology, Tidal.Range, Salinity.category,
         Local.linear.water.level.change.rate...19.yr.rate) %>%
  #Rename columns for ease of coding for the rest of the scripts
  rename(
    Region = NERR.Region,
    abiotic_cover = Total.unvegetated,
    live_cover = Total.live.veg,
    halophyte = H.Halophyte,
    freshwater = F.Freshwater,
    invasive_cover = Invasive_Cover,
    salt_ratio = Salt_to_Total,
    tidal_range = Tidal.Range,
    salinity = Salinity.category,
    SLR_Rate_19yrs = Local.linear.water.level.change.rate...19.yr.rate) %>%
#Change all NAs of visual cover types to zero
mutate(across(abiotic_cover:salt_ratio,
              ~ifelse(is.na(.) == TRUE, 0,.)))

glimpse(veg_format)

write.csv(veg_format,
          "Formatted Datasets\\Veg Plot Dataframe Formatted.csv")

#Task 3: Summarize the national plot data frame by vegetation zone and site

veg_zone <- veg_format %>%
  select(-TransectID, -PlotID) %>%
  group_by(Reserve, SiteID, Year, Vegetation_Zone, Region, Geomorphology, 
           tidal_range, salinity, SLR_Rate_19yrs) %>%
  summarise(across(abiotic_cover:salt_ratio, 
                   ~mean(., na.rm = TRUE))) %>%
  mutate(across(abiotic_cover:salt_ratio, 
                ~round(., 2))) %>%
  ungroup()

glimpse(veg_zone)

write.csv(veg_zone, "Formatted Datasets\\Veg Dataframe Summarised By Site and Zone.csv")


veg_zone_info <- veg_zone %>%
  group_by(Vegetation_Zone) %>%
  summarise(
    Site_total = length(unique(SiteID))
  )

glimpse(veg_zone_info)

#Task 4: Summarize the national plot data frame by site

veg_site <- veg_format %>%
  select(-TransectID, -PlotID) %>%
  group_by(Reserve, SiteID, Year, Region, Geomorphology, 
           tidal_range, salinity, SLR_Rate_19yrs) %>%
  summarise(across(abiotic_cover:salt_ratio, 
                   ~mean(., na.rm = TRUE))) %>%
  mutate(across(abiotic_cover:salt_ratio, 
                ~round(., 2))) %>%
  ungroup()

glimpse(veg_site)

write.csv(veg_site, "Formatted Datasets\\Veg Dataframe Summarised By Site.csv")


#Chapter 3: Format the Slope - Site Data Frame ---------------------------------

#Task 1: Import the vegetation metric slope data frame summarized by site 

slope_site <- read.csv("Input Data\\Veg Metric Slopes by Site.csv")

#Task 2: Format the vegetation metric slope data frame

slope_format <- slope_site %>%
  # Select the relevant site identifiers, characteristics, explanatory factors, and vegetation metrics
  select(Reserve, SiteID,
         Total.unvegetated_slope, Total.live.veg_slope, H.Halophyte_slope,
         F.Freshwater_slope, EMI_slope, Invasive_Cover_slope, Richness_slope, 
         SWdiv_slope, Salt_to_Total_slope, 
         NERR_Region, Geomorphology, Tidal_Range, Salinity_category,
         SLR_last19yrs, NERRs_Landscape_resiliency_condition_sum_quantile) %>%
  # Rename columns for ease of coding in future scripts
  rename(
    Region = NERR_Region,
    abiotic_cover = Total.unvegetated_slope,
    live_cover = Total.live.veg_slope,
    halophyte = H.Halophyte_slope,
    freshwater = F.Freshwater_slope,
    invasive_cover = Invasive_Cover_slope,
    salt_ratio = Salt_to_Total_slope,
    tidal_range = Tidal_Range,
    salinity = Salinity_category)

glimpse(slope_format)

write.csv(slope_format,
          "Formatted Datasets\\Veg Slope by Site Formatted.csv")



# Chapter 4: Format the Slope - Site - Veg Zone Data Frame ---------------------------

#Task 1: Import the vegetation metric slope data frame summarized by site and veg zone

slope_zone <- read.csv("Input Data\\Veg Metric Slopes by Site and Zone.csv")

#Task 2: Format the vegetation metric slope data frame

slope_format <- slope_zone %>%
  # Select the relevant site identifiers, characteristics, explanatory factors, and vegetation metrics
  select(Reserve, SiteID, Vegetation_Zone,
         Total.unvegetated_slope, Total.live.veg_slope, H.Halophyte_slope,
         F.Freshwater_slope, EMI_slope, Invasive_Cover_slope, Richness_slope, 
         SWdiv_slope, Salt_to_Total_slope, 
         NERR_Region, Geomorphology, Tidal_Range, Salinity_category,
         SLR_last19yrs, NERRs_Landscape_resiliency_condition_sum_quantile) %>%
  # Rename columns for ease of coding in future scripts
  rename(
    Region = NERR_Region,
    abiotic_cover = Total.unvegetated_slope,
    live_cover = Total.live.veg_slope,
    halophyte = H.Halophyte_slope,
    freshwater = F.Freshwater_slope,
    invasive_cover = Invasive_Cover_slope,
    salt_ratio = Salt_to_Total_slope,
    tidal_range = Tidal_Range,
    salinity = Salinity_category)

glimpse(slope_format)

write.csv(slope_format,
          "Formatted Datasets\\Veg Slope by Site and Zone Formatted.csv")


#Chapter 5: Summary Statistics for the Plot & Site Datasets -----------------------------


#Create a data summary table for the entire plot dataframe with the skim() function
# from the skimr package

#Task 1: Create and export a summary table for the Plot Dataframe

#Note - Used the original, non-formatted dataset, to conduct skimr data summary
# We wanted to analyze all of the vegetation cover types and retain the NAs in the original dataset
veg_plot_skim <- veg %>%
  select(Total.unvegetated:Salicornia.pacifica) %>%
  #Convert all zeroes into NAs to analyze the presence/absence of all vegetation metrics
  mutate(across(Total.unvegetated:Salicornia.pacifica,
                ~ifelse(. == 0, NA, .))) %>%
  skim_without_charts(.) %>%
  mutate(across(where(is.numeric),
                ~round(., 2)))
       

veg_plot_skim

write.csv(veg_plot_skim,
          "Output Stats\\Dataset Summary - Skim - Plot Dataframe.csv")


# Task 2: Summarise the distribution of plots and sites for each of the categorical variables

namaste_info <- veg_format %>%
  mutate(
    SiteID = as.factor(SiteID),
    Site_Year = paste(SiteID, Year, " - ")) %>%
  gather(c(Region:salinity, Vegetation_Zone),
         key = "Explanatory_Metric",
         value = "Factor") %>%
  group_by(Explanatory_Metric, Factor) %>%
  summarise(
    Plot_Count = n(),
    Site_Count = length(unique(SiteID)),
    Site_Year_Count = length(unique(Site_Year))) %>%
  ungroup()

glimpse(namaste_info)

write.csv(namaste_info,
          "Output Stats\\Dataset Summary - Summary of Site Counts Dataframe.csv")











