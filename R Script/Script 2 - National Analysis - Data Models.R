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

#Script Description: Script 2 imports the formatted vegetation dataset and conducts the models for
#   for Abiotic Cover

#Chapter 1: Import package library

#Data Analysis

library(dplyr)
library(tidyr)
library(lme4)
library(afex)
library(rstatix)
library(broom)
library(broom.mixed)
library(purrr)

#Chapter 2: Format the vegetation dataset for broom - map functions to create numerous models

veg <- read.csv("Formatted Datasets\\Veg Plot Dataframe Formatted.csv")

#Chapter 1: Time Model

time_mod <- lmer(abiotic_cover ~ Year + Vegetation_Zone + (1|SiteID) + (1|Reserve),
                  data = veg)

time_anova <- tidy(anova(time_mod))



  