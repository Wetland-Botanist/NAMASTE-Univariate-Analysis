#Project: NAMASTE National Univariate Exploratory Analysis
#Script: Script 6 - National Mixed Linear Regressions Across Study

#Author: Grant McKown, Jackson Estuarine Laboratory, University of New Hampshire
#Contact: james.mckown@unh.edu

#Project Description: The NERR system collects vegetation data along permanent transects and plots
#   in salt marshes throughout the United States. Vegetation data was synthesized, formatted, and 
#   filtered (based on field monitoring requirements). Reserve, site, transect, and plot-level characteristics
#   were also collected including hydrology, landscape geomorphology, elevation, and climate. The 
#   series of scripts conducts an exploratory analysis of different vegetation response factors to 
#   a series of explanatory and grouping factors. 

#Script Description: Script 6 conducts national regression across the entire study region. Script 6
# is very similar to scripts 4 & 5. 



#Chapter 1: Import package library --------------------------------------------


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


#Chapter 2: Format the vegetation dataset for broom - map functions ------------------


veg <- read.csv("Formatted Datasets\\Veg Dataframe Summarised By Site.csv") %>%
  select(-X)

glimpse(veg)

# Formatting and organizing the dataset for regression slopes including:
# Gathering the Site Variables and gathering the vegetation metrics, thus transforming the dataset
# into a very LONG dataset. Essentially, breaking down each row into ~30 rows. This will allow
# for simultaneous linear regressions for each vegetation metric for each site variable!

veg_format <- veg %>%
  # Format the data table from wide to long
  gather(abiotic_cover:salt_ratio,
         key = "Metric",
         value = "Value") %>%
  # Remove any NA values from the dataset (QAQCing the dataset before analysis)
  filter(!is.na(Value)) %>%
  # Long winded code to rename the variables 
  mutate(Metric = ifelse(Metric == "abiotic_cover", "Abiotic Cover", 
                         ifelse(Metric == "live_cover", "Live Cover",
                                ifelse(Metric == "halophyte", "Halophyte Cover",
                                       ifelse(Metric == "SWdiv", "Shannon-Weiner Diversity", 
                                              ifelse(Metric == "freshwater", "Freshwater Cover",
                                                     ifelse(Metric == "salt_ratio", "Salt Ratio", 
                                                            ifelse(Metric == "invasive_cover", "Invasive Cover",
                                                                   Metric)))))))) %>%
  # Assign factors to each vegetation metric for review and graphing purposes
  mutate(Metric = factor(Metric, levels = c("Abiotic Cover", "Live Cover", "Halophyte Cover",
                                            "Freshwater Cover", "Invasive Cover", "EMI", 
                                            "Salt Ratio", "Shannon-Weiner Diversity", "Richness")))

glimpse(veg_format)



#Chapter 3: Mixed linear regressions of each vegetation metric -------------------------------


# This is accomplished by nesting the data essentially by Vegetation Zone 
# and Vegetation Metric (long dataset)

# Mixed linear modeling is a simple time model with random factors of Reserve and nested Site 

regression_models_site <- veg_format %>%
  group_by(Metric) %>%
  # Calculate the sample size for reporting in the ANOVA table
  mutate(SampleSize = n()) %>%
  ungroup() %>%
  # Group the metrics to create mixed linear models over time
  group_by(Metric, SampleSize) %>%
  nest() %>%
  # Mixed linear models for each vegetation metric
  mutate(regression = map(.x = data,
                          ~lmer(Value ~ Year + (1|Reserve/SiteID),
                                data = .x) %>%
                            # Create the ANOVA Table for each model
                            anova() %>%
                            # select(-NumDF, -DenDF) %>%
                            # Format the ANOVA table into a dataframe
                            tidy())) %>%
  # Unnest and keep the ANOVA tables of the mixed linear models
  unnest(regression) %>%
  select(-data) %>%
  # Format the ANOVA table by rounding all values to 3 decimal places
  mutate(across(sumsq:p.value, ~round(., 3))) %>%
  ungroup() 

glimpse(regression_models_site)

# Export and save the ANOVA tables of the time series mixed linear models
write.csv(regression_models_site,
          "Output Stats\\Mixed Linear Time Model - Regression Model Outputs Across Study.csv")




# Chapter 4: Calculate the Slopes of each Vegetation Metric -------------------------

# Run a dplyr loop to predict the values of the national vegetation regressions

regression_predicted <- veg_format %>%
  group_by(Metric) %>%
  nest() %>%
  # Again, create the mixed models and then use the models to predict the slopes
  mutate(regression = map(.x = data,
                          # Create the mixed linear models
                          ~lmer(Value ~ Year + (1|Reserve/SiteID),
                                data = .x) %>%
                            # Use ggpredict to predict values of each model over time
                            ggpredict(., 
                                      terms = c("Year [all]"),
                                      type = "fixed", interval = "confidence"))) %>%
  select(-data) %>%
  # Pull out the table of predicted values
  unnest(regression) %>%
  rename(Year = x,
         Value_pred = predicted) %>%
  # Calculate the standard error low and high for graphing purposes
  mutate(stderr.low = Value_pred - std.error,
         stderr.high = Value_pred + std.error) %>%
  select(Metric, Year, Value_pred, std.error, stderr.low, stderr.high, conf.low, conf.high) %>%
  mutate(across(Value_pred:stderr.high, ~round(., 3)))

glimpse(regression_predicted)


#Calculate the slopes of each vegetation metric

regression_slopes <- regression_predicted %>%
  group_by(Metric) %>%
  # Use the predicted values to calculate the slope based on first and last 
  # predicted value
  summarise(slope = (Value_pred[which(Year == max(Year))] - Value_pred[which(Year == min(Year))]) / ((max(Year) - min(Year)))) %>%
  ungroup() %>%
  mutate(slope = round(slope, 3))

write.csv(regression_slopes,
          "Output Stats\\Mixed Linear Time Model - Slopes Across Across Study.csv")




#Chapter 5: Graph the General Trends Across Time -------------------------------

#Data visualization is broken down into 4 graphs:
# (1) Vegetation Cover - Live, Abiotic, Freshwater, and Halophyte
# (2) Species - Richness, Diversity
# (3) Ratios - Salt Ratio, EMI
# (4) Combine all three graphs into one comprehensive, national time model graph



#Task 1: Graph the results of the predicted values for the Abiotic, Live Cover, Halophyte Cover, 
#and Freshwater Cover

regression_predicted_cover <- regression_predicted %>%
  filter(Metric == "Halophyte Cover") %>%
  mutate(conf.high = ifelse(conf.high > 100, 100, conf.high), 
         conf.low = ifelse(conf.low < 0, 0, conf.low))
  
veg_graph <- veg_format %>%
  filter(Metric == "Halophyte Cover")
         
         
national_cover_graph <- ggplot(data = regression_predicted_cover,
                                        aes(x = Year,
                                            y = Value_pred)) + 
          geom_point(data = veg_graph,
                     aes(x = Year,
                         y = Value),
                     colour = "darkblue",
                     size = 4.5, alpha = 0.35) +
  geom_ribbon(aes(x = Year, 
                  ymin = conf.low, ymax = conf.high),
              alpha = 0.75, fill = "gray") + 
           geom_line(linewidth = 1.5, 
                     colour = "black", linetype = "dashed") + 
           scale_y_continuous(limits = c(0, 105),
                              breaks = seq(0, 100, 20),
                              expand = c(0,0)) +
           scale_x_continuous(limits = c(2005, 2022.5),
                              breaks = seq(2006, 2022, 2),
                              expand = c(0,0)) +
           labs(y = "Halophoyte Cover (%)",
                x = "") +
           theme_bw() +
           theme(
             legend.position = c(0.10, 0.60),
             legend.title = element_blank(),
             legend.text = element_text(size = 18, colour = "black"),
             panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank(),
             axis.title = element_text(size = 24, colour = "black"),
             axis.text = element_text(size = 22, colour = "black"),
             strip.background = element_blank(),
             strip.text = element_blank()) +
  facet_wrap(~Metric)
         
  national_cover_graph
         

ggsave(national_cover_graph,
          filename = "Output Figures\\National Time Mixed Model - Halophyte Color - Whole Study.jpg",
          units = "in",
          height = 8, width = 12, dpi = 300, limitsize = FALSE)
         
         
#Task 2 - Graph the results of the predicted values for the Shannon-Diversity & Richness
         
         
 # Note - due to no significant trend found for either metric, the regression_predicted, not broken down by 
 # different regions, will be shown. 
         
regression_predicted_species <- regression_predicted %>%
           filter(Metric == "EMI") %>%
           mutate(conf.low = ifelse(conf.low < 0, 0, conf.low))

veg_graph <- veg_format %>%
           filter(Metric == "EMI")
         
         national_cover_graph <- ggplot(data = regression_predicted_species,
                                        aes(x = Year,
                                            y = Value_pred)) + 
           geom_point(data = veg_graph,
                      aes(x = Year,
                          y = Value),
                      size = 4.5, alpha = 0.35, colour = "darkblue") +
           geom_ribbon(aes(x = Year, ymin = conf.low, ymax = conf.high),
                       alpha = 0.75, fill = "gray") + 
           geom_line(linewidth = 1.5,
                     colour = "black", linetype = "dashed") + 
           scale_y_continuous(limits = c(0, 1.1),
                              breaks = seq(0, 1, 0.25),
                              expand = c(0,0)) +
           scale_x_continuous(limits = c(2005, 2022.5),
                              breaks = seq(2006, 2022, 2),
                              expand = c(0,0)) +
           labs(y = "Ecotone Migration Index",
                x = "") +
           theme_bw() +
           theme(
             legend.position = c(0.10, 0.60),
             legend.title = element_blank(),
             legend.text = element_text(size = 18, colour = "black"),
             panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank(),
             axis.title = element_text(size = 18, colour = "black"),
             axis.text = element_text(size = 18, colour = "black"),
             strip.background = element_blank(),
             strip.text = element_blank()) +
           facet_wrap(~Metric)
         
         national_cover_graph
         
         
         ggsave(national_cover_graph,
                filename = "Output Figures\\National Time Mixed Model - EMI - Whole Study.jpg",
                units = "in",
                height = 8, width = 12, dpi = 300, limitsize = FALSE)
         
         