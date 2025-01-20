#Project: NAMASTE National Univariate Exploratory Analysis
#Script: Script 10  - Overview Figure of Univariate Analysis


#Author: Grant McKown, Jackson Estuarine Laboratory, University of New Hampshire
#Contact: james.mckown@unh.edu

#Project Description: The NERR system collects vegetation data along permanent transects and plots
#   in salt marshes throughout the United States. Vegetation data was synthesized, formatted, and 
#   filtered (based on field monitoring requirements). Reserve, site, transect, and plot-level characteristics
#   were also collected including hydrology, landscape geomorphology, elevation, and climate. The 
#   series of scripts conducts an exploratory analysis of different vegetation response factors to 
#   a series of explanatory and grouping factors. 

#Script Description: Script 10 creates an overview figure of the Univariate Analysis with 
# halophyte cover as the example vegetation metric.


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
library(purrr)
library(ggResidpanel)

#Graphing and Data Visualization Packages

library(ggplot2)
library(patchwork)



#-----------------------------------------------------------------------------
# Chapter 2: Create compilation of site-level regressions
#-----------------------------------------------------------------------------

#Task 1: Import the Formatted Plot Data frame

plots <- read.csv("Formatted Datasets\\Veg Plot Dataframe Formatted.csv") %>%
  select(-X) %>%
  mutate(Reserve_Site = paste(Reserve, SiteID, delimiter = " - "))

glimpse(plots)


#Task 2: Graph simple linear regressions of halophyte cover

site_regression <- ggplot(
  data = plots,
  aes(x = Year,
      y = EIR)) +
  geom_point(
    colour = "darkblue", alpha = 0.50) +
  geom_smooth(
    method = "lm",
    linewidth = 1.25, colour = "orange", 
    fill = "gray") +
  scale_y_continuous(limits = c(-0.05, 1.05),
                     breaks = seq(0, 1, 0.20),
                     expand = c(0,0)) +
  scale_x_continuous(limits = c(2005, 2023),
                     breaks = seq(2006, 2022, 2),
                     expand = c(0,0)) +
  labs(y = "",
       x = "") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 10)) +
  facet_wrap(~Reserve_Site,
             ncol = 8,
             nrow = 7)

site_regression  

ggsave(site_regression,
       filename = "Output Figures\\Overview Figure - Site Regression - EMI.jpg",
       units = "in",
       height = 8, width = 12, dpi = 300, limitsize = FALSE)

#-----------------------------------------------------------------------------
# Chapter 2: Create compilation of region-level regressions
#-----------------------------------------------------------------------------

#Task 1: Import the Formatted Plot Dataframe

sites <- read.csv("Formatted Datasets\\Veg Dataframe Summarised by Site.csv") %>%
  select(-X) %>%
  mutate(Region = factor(Region, levels = c("Northeast", "Mid-Atlantic", "Southeast",
                                            "Gulf Coast", "West Coast")))

glimpse(sites)


#Task 2: Graph simple linear regressions of halophyte cover

region_regression <- ggplot(
  data = sites,
  aes(x = Year,
      y = EIR)) +
  geom_point(
    colour = "darkblue", alpha = 0.50, size = 3.5) +
  geom_smooth(
    method = "lm",
    linewidth = 1.25, colour = "orange", 
    fill = "gray") +
  scale_y_continuous(limits = c(-0.1, 1.05),
                     breaks = seq(0, 1, 0.20),
                     expand = c(0,0)) +
  scale_x_continuous(limits = c(2005, 2023),
                     breaks = seq(2006, 2022, 2),
                     expand = c(0,0)) +
  labs(y = "",
       x = "") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, colour = "black")) +
  facet_wrap(~Region,
             ncol = 3,
             nrow = 2)

region_regression  

ggsave(region_regression,
       filename = "Output Figures\\Overview Figure - Region Regression - EMI.jpg",
       units = "in",
       height = 8, width = 12, dpi = 300, limitsize = FALSE)

  
  
  
  
  
  
  









