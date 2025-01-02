#Project: NAMASTE National Univariate Exploratory Analysis
#Script: Script 4 - National Mixed Linear Regressions by Vegetation Zone

#Author: Grant McKown, Jackson Estuarine Laboratory, University of New Hampshire
#Contact: james.mckown@unh.edu

#Project Description: The NERR system collects vegetation data along permanent transects and plots
#   in salt marshes throughout the United States. Vegetation data was synthesized, formatted, and 
#   filtered (based on field monitoring requirements). Reserve, site, transect, and plot-level characteristics
#   were also collected including hydrology, landscape geomorphology, elevation, and climate. The 
#   series of scripts conducts an exploratory analysis of different vegetation response factors to 
#   a series of explanatory and grouping factors. 

#Script Description: Script 4 formats the site-level vegetation data, conducts linear mixed models at the 
#                   the vegetation zone level across the nation. The slopes of significant trends are
#                   calculated. Model summaries and regression slopes are exported as CSV. Graphs are
#                   created at the end of the script. 


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
#Chapter 2: Format the vegetation dataset for broom - map functions to create numerous models
#---------------------------------------------------------------------------------------------

veg <- read.csv("Formatted Datasets\\Veg Dataframe Summarised By Site and Zone.csv") %>%
  select(-X) %>%
  dplyr::rename(EMI = EIR)

glimpse(veg)



#-----------------------------------------------------------------------------
#Chapter 3: Calculate live cover slope of each individual site, vegetation zone
#-----------------------------------------------------------------------------

# Formatting and organizing the dataset for regression slopes including:
# Gathering the Site Variables and gathering the vegetation metrics, thus transforming the dataset
# into a very LONG dataset. Essentially, breaking down each row into ~30 rows. This will allow
# for simultaneous linear regressions for each vegetation metric for each site variable!

veg_format <- veg %>%
  #Transform the summarized data frame from Wide --> Long for vegetation metrics
  gather(abiotic_cover:salt_ratio,
         key = "Metric",
         value = "Value") %>%
  #Remove pesky NA values from the dataset -- mostly in invasive and freshwater cover
  filter(!is.na(Value)) %>%
  #Rename the vegetation metrics for ease of coding later in script
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
                                            "Salt Ratio", "Shannon-Weiner Diversity", "Richness")),
         Vegetation_Zone = factor(Vegetation_Zone,
                                  levels = c("Low", "Mid", "Up")))

glimpse(veg_format)


#------------------------------------------------------------------------------------------
#Chapter 4: National mixed linear regressions of each vegetation metric by vegetation zone
#------------------------------------------------------------------------------------------

# This is accomplished by nesting the data essentially by Vegetation Zone and Vegetation Metric (long dataset)

# Mixed linear modeling is a simple time model with random factors of Reserve and nested Site 

regression_models_veg <- veg_format %>%
  group_by(Metric) %>%
  #Calculate the sample size for each national regression
  mutate(SampleSize = length(Reserve),
         Vegetation_Zone = as.character(Vegetation_Zone)) %>%
  ungroup() %>%
  #Group and nest the data set by Vegetation Metric
  group_by(Metric, SampleSize) %>%
  nest() %>%
  #Calculate the regressions using map() and broom.mixed package 
  mutate(regression = map(.x = data,
                          ~lmer(Value ~ Year * Vegetation_Zone + (1|Reserve/SiteID),
                                      data = .x) %>%
                            anova() %>%
                            select(-NumDF, -DenDF) %>%
                            tidy())) %>%

  unnest(regression) %>%
  select(-data) %>%
    mutate(across(sumsq:p.value, ~round(., 3))) %>%
  ungroup() 

glimpse(regression_models_veg)


write.csv(regression_models_veg,
          "Output Stats\\Mixed Linear Time Model - Regression Model Outputs by Vegetation Zone.csv")



#---------------------------------------------------------------------------------------
#Chapter 5: Calculate the Slopes of each Vegetation Metric Across All Vegetation Zones 
#---------------------------------------------------------------------------------------

# Run a dplyr loop to predict the values of the regressions across time itself (not considering zone)

regression_predicted <- veg_format %>%
  group_by(Metric) %>%
  nest() %>%
  #Similar code format so far with map() and broom.mixed package as in Chapter 4 of Script
  mutate(regression = map(.x = data,
                          ~lmer(Value ~ Year * Vegetation_Zone + (1|Reserve/SiteID),
                                data = .x) %>%
                            #Using ggpredict() to calculate national trends, irrespective of zone 
                            ggpredict(., 
                                      terms = c("Year [all]"),
                                      type = "fixed", interval = "confidence"))) %>%
  select(-data) %>%
  unnest(regression) %>%
  #Rename the default column names produced by ggpredict()
  rename(Year = x,
         Value_pred = predicted) %>%
  mutate(across(Value_pred:conf.high, ~round(., 3)))

glimpse(regression_predicted)


#Calculate the slopes of each vegetation metric across Vegetation Zone

regression_slopes <- regression_predicted %>%
  group_by(Metric) %>%
  summarise(slope = (Value_pred[which(Year == max(Year))] - Value_pred[which(Year == min(Year))]) / ((max(Year) - min(Year)))) %>%
  ungroup() %>%
  mutate(slope = round(slope, 3))

write.csv(regression_slopes,
          "Output Stats\\Mixed Linear Time Model - Slopes Across All Zones.csv")



#---------------------------------------------------------------------------------------
#Chapter 6: Calculate the Slopes of each Vegetation Metric For Each Vegetation Zone 
#---------------------------------------------------------------------------------------

#Re-run the combination of map, ggpredict chunk of code, except with the addition of the "Vegetation_Zone [all]"
# to tell the software to also include the Vegetation Zone factor in predicted values

regression_predicted_veg <- veg_format %>%
  group_by(Metric) %>%
  nest() %>%
  mutate(regression = map(.x = data,
                          ~lmer(Value ~ Year * Vegetation_Zone + (1|Reserve/SiteID),
                                data = .x) %>%
        #Note the inclusion of Vegetation Zone in the ggpredict() function to predict
          # the trends for each vegetation zone for each vegetation metric
        ggpredict(., 
                  terms = c("Year [all]", "Vegetation_Zone [all]"),
                  type = "fixed", interval = "confidence"))) %>%
  select(-data) %>%
  unnest(regression) %>%
  rename(Year = x,
         Value_pred = predicted,
         Vegetation_Zone = group) %>%
  mutate(across(Value_pred:conf.high, ~round(., 3)))

glimpse(regression_predicted_veg)


#Calculate the slopes of each vegetation metric across Vegetation Zone

regression_slopes_veg <- regression_predicted_veg %>%
  group_by(Vegetation_Zone, Metric) %>%
  summarise(slope = (Value_pred[which(Year == max(Year))] - Value_pred[which(Year == min(Year))]) / ((max(Year) - min(Year)))) %>%
ungroup() %>%
  mutate(slope = round(slope, 3))

write.csv(regression_slopes_veg,
          "Output Stats\\Mixed Linear Time Model - Slopes by Vegetation Zone.csv")



#-----------------------------------------------------------------------------
#Chapter 7: Graph the General Trends Across Time and Vegetation Zone
#----------------------------------------------------------------------------

#Data visualization is broken down into 4 graphs:
  # (1) Vegetation Cover - Live, Abiotic, Freshwater, and Halophyte
  # (2) Species - Richness, Diversity
  # (3) Ratios - Salt Ratio, EMI
  # (4) Combine all three graphs into one comprehensive, national time model graph



#Task 1: Graph the results of the predicted values for the Abiotic, Live Cover, Halophyte Cover, and Freshwater Cover

regression_predicted_cover <- regression_predicted %>%
  filter(Metric == "Abiotic Cover" |
          Metric == "Live Cover" |
          Metric == "Halophyte Cover" |
          Metric == "Freshwater Cover") %>%
  mutate(conf.high = ifelse(conf.high > 100, 100, conf.high),
         conf.low = ifelse(conf.low < 0, 0, conf.low))

national_cover_graph <- ggplot(data = regression_predicted_cover,
                                 aes(x = Year,
                                     y = Value_pred,
                                     group = Metric)) + 
  geom_line(aes(colour = Metric),
            linewidth = 1.25, linetype = 'dashed') + 
  geom_ribbon(aes(x = Year, ymin = conf.low, ymax = conf.high,
              fill = Metric),
              alpha = 0.5) +   
  scale_y_continuous(limits = c(0, 102),
                      breaks = seq(0, 100, 20),
                      expand = c(0,0)) +
  scale_x_continuous(limits = c(2005, 2023),
                     breaks = seq(2006, 2022, 2),
                     expand = c(0,0)) +
  labs(y = "Visual Cover (%)",
       x = "") +
  theme_bw() +
  theme(
    legend.position = c(0.15, 0.875),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 15, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18))


national_cover_graph


ggsave(national_cover_graph,
       filename = "Output Figures\\National Time Mixed Model - Veg Covers - Veg Zone.jpg",
       units = "in",
       height = 8, width = 12, dpi = 300, limitsize = FALSE)



#Task 2 - Graph the results of the predicted values for the Shannon-Diversity & Richness

regression_predicted_richness <- regression_predicted %>%
  filter(Metric == "Richness" |
           Metric == "Shannon-Weiner Diversity")

national_richness_graph <- ggplot(data = regression_predicted_richness,
                               aes(x = Year,
                                   y = Value_pred,
                                   group = Metric)) + 
  geom_line(aes(colour = Metric),
            linewidth = 1.25, linetype = "dashed") + 
  geom_ribbon(aes(x = Year, ymin = conf.low, ymax = conf.high,
                  fill = Metric),
              alpha = 0.5) +   
  scale_y_continuous(limits = c(0, 3.5),
                     breaks = seq(0, 3.5, 0.5),
                     expand = c(0,0)) +
  scale_x_continuous(limits = c(2005, 2024),
                     breaks = seq(2006, 2022, 2),
                     expand = c(0,0)) +
  labs(y = "",
       x = "") +
  theme_bw() +
  theme(
    legend.position = c(0.15, 0.875),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 15, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18))


national_richness_graph


ggsave(national_richness_graph,
       filename = "Output Figures\\National Time Mixed Model - Richness & Diversity.jpg",
       units = "in",
       height = 8, width = 12, dpi = 300, limitsize = FALSE)



#Task 3 - Graph the results of the predicted values for the Salt Ratio and EMI

regression_predicted_ratio <- regression_predicted_veg %>%
  filter(Metric == "EMI") %>%
  mutate(conf.high = ifelse(conf.high > 1, 1, conf.high))
  

national_ratio_graph <- ggplot(data = regression_predicted_ratio,
                                  aes(x = Year,
                                      y = Value_pred,
                                      group = Vegetation_Zone)) + 
  geom_line(aes(colour = Vegetation_Zone),
            linewidth = 1.25, linetype = "dashed") + 
  geom_ribbon(aes(x = Year, ymin = conf.low, ymax = conf.high,
                  fill = Vegetation_Zone),
              alpha = 0.5) +   
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.20),
                     expand = c(0,0)) +
  scale_x_continuous(limits = c(2005, 2023),
                     breaks = seq(2006, 2022, 2),
                     expand = c(0,0)) +
  labs(y = "EMI",
       x = "") +
  theme_bw() +
  theme(
    legend.position = c(0.10, 0.875),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 15, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18))


national_ratio_graph


ggsave(national_ratio_graph,
       filename = "Output Figures\\National Time Mixed Model - EMI.jpg",
       units = "in",
       height = 8, width = 12, dpi = 300, limitsize = FALSE)




#Task 4 - Combine the three graphs and export it for the ultimate Natinal Time Mixed Model Figure


national_model_graph <- national_cover_graph / national_richness_graph / national_ratio_graph

national_model_graph


ggsave(national_model_graph,
       filename = "Output Figures\\National Time Mixed Model - Combined Graph.jpg",
       units = "in",
       height = 14, width = 20, dpi = 300, limitsize = FALSE)
