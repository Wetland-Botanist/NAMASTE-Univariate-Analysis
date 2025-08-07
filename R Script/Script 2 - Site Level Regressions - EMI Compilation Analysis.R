#Project: NAMASTE National Univariate Exploratory Analysis
#Script: Script 2 - Site Level Regressions and Compilation for Ecotone Migration Index

#Author: Grant McKown, Jackson Estuarine Laboratory, University of New Hampshire
#Contact: james.mckown@unh.edu

#Project Description: The NERR system collects vegetation data along permanent transects and plots
#   in salt marshes throughout the United States. Vegetation data was synthesized, formatted, and 
#   filtered (based on field monitoring requirements). Reserve, site, transect, and plot-level characteristics
#   were also collected including hydrology, landscape geomorphology, elevation, and climate. The 
#   series of scripts conducts an exploratory analysis of different vegetation response factors to 
#   a series of explanatory and grouping factors. 

#Script Description: Script 2 makes small formatting adjustments to the 
# 'Veg Plot Formatted Dataframe' from Script 1. The main objectives of the 
#  script are to run simple linear regressions over time of vegetation metrics
#  within each individual site or vegetation zone of each site. Next, a compilation analysis
#  is conducted to weigh the directionality and significance of the trends of EMI. Compilation analysis
#  is then exported in CSV and graphed accordingly. 


#Chapter 1: Import package library ----------------------------------------------

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


#Chapter 2: Format the vegetation dataset ---------------------------------------

veg <- read.csv("Formatted Datasets\\Veg Plot Dataframe Formatted.csv") %>%
  select(-X) %>%
  select(Reserve:Year, EMI, Region:salinity)

glimpse(veg)



#Chapter 3: EMI Compilation Analysis by Vegetation Zone -------------------------

#Task 1: Formatting and organizing the plot dataset for regression slopes

veg_format <- veg %>%
  #For data visualization down the road, assign factors to the Vegetation Zone
  mutate(Vegetation_Zone = factor(Vegetation_Zone,
                                  levels = c("Low", "Mid", "Up")))

glimpse(veg_format)


#Task 2: Reduce the formatted vegetation dataset to only focus on vegetation zone

veg_format_zone <- veg_format %>%
  select(-c(Region, Geomorphology, tidal_range, salinity))

#Task 3: Calculate the slopes and p-values of individual linear regressions for each site across all 
# vegetation metrics and vegetation zones.

# This is accomplished by nesting the data essentially by Site, Vegetation Zone, and Metric
regression_slopes_zone <- veg_format_zone %>%
  group_by(Reserve, SiteID, Vegetation_Zone) %>%
  mutate(Sample_Size = n()) %>%
  ungroup() %>%
  group_by(Reserve, SiteID, Vegetation_Zone, Sample_Size) %>%
  nest() %>%
  # Using map() function and broom package, a regression is run for each nested dataset.
  # The model summary is output, which provides slope and p-value of regression
  mutate(regression = map(.x = data,
                                    ~summary(lm(EMI ~ Year,
                                                  data = .x)) %>%
                                      tidy())) %>%
  unnest(regression) %>%
  #Format the summary table by removing unwanted terms, keeping the "Year" term, and 
  # formatting the number of decimal points
  filter(term == 'Year') %>%
  select(-data) %>% 
  mutate(across(estimate:p.value, ~round(., 4))) %>%
  ungroup() %>%
  #There were a handful of sites - vegetation zones with "perfect" EMI regressions (aka EMI value never
  # changed and was typically 1 or 0), so manually assigning a p-value of 1, since they are perfect
  # horizontal lines by definition of a statistical regression. 
  mutate(p.value = ifelse(is.na(p.value), 1, p.value)) %>%
  rename(slope = estimate,
         pvalue = p.value)


glimpse(regression_slopes_zone)



#Task 4: Calculate the compilation analysis scores based on the p-scores and the slopes direction

# The compilation analysis score is based on two criteria: p-value and direction of slope
  # The magnitude of the score is based on the p-value 
  # The score is then multiplied by the direction of the score (i.e., negative slope = negative score)
regression_scores_zone <- regression_slopes_zone %>%
         mutate(score = ifelse(pvalue <= 0.05, 2, 
                        ifelse(pvalue > 0.20, 0, 1)),
         score = ifelse(slope < 0, score * -1, score))

  
glimpse(regression_scores_zone)

write.csv(regression_scores_zone,
          "Output Stats\\Compilation Analysis - EMI Scores & Slopes - Veg Zone.csv")



#Task 5: Calculate descriptive statistics for the slope and score by vegetation zone

score_stats_zone <- regression_scores_zone %>%
  group_by(Vegetation_Zone) %>%
  #Average and standard error of the slopes and compilation analysis scores
  summarise(slope.m = mean(slope, na.rm = TRUE),
            slope.se = sd(slope, na.rm = TRUE)/sqrt(n()),
            
            score.m = mean(score, na.rm = TRUE),
            score.se = sd(score, na.rm = TRUE)/sqrt(n())) %>%
  mutate(across(slope.m:score.se, ~round(., 3)))

glimpse(score_stats_zone)


write.csv(score_stats_zone,
          "Output Stats\\Compilation Analysis - EMI Score Stats - Veg Zone.csv")


#Task 6: Graph the EMI compilation analysis by vegetation zone

EMI_veg_zone_graph <- ggplot(data = score_stats_zone,
                             aes(x = Vegetation_Zone, 
                                 y = score.m)) +
  geom_bar(aes(fill = Vegetation_Zone),
           stat = 'identity', position = position_dodge(0.9),
           linewidth = 1.5, colour = "black") + 
  geom_errorbar(aes(x = Vegetation_Zone,
                    ymax = score.m + score.se, 
                    ymin = score.m - score.se,
                    group = Vegetation_Zone),
                position = position_dodge(0.9),
                width = 0.5, linewidth = 1.25) + 
  scale_y_continuous(limits = c(-1.1, 1.6),
                     breaks = seq(-1.0, 1.5, 0.5),
                     expand = c(0,0)) +
  labs(y = "EMI Score",
       x = "Vegetation Zone") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 20, colour = "black"),
    axis.text = element_text(size = 20, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 20))

EMI_veg_zone_graph


ggsave(EMI_veg_zone_graph,
       filename = "Output Figures\\EMI Compilation Analysis by Vegetation Zone.jpg",
       units = "in",
       height = 8, width = 12, dpi = 300, limitsize = FALSE)





#Chapter 4: EMI Compilation Analysis by Site Characteristics -----------------------------

# Note the code in Chapter 3 is very similar to Chapter 2 with minor changes to accommodate
# the change in the geographic scope of the analysis: from Vegetation Zone --> Site Level.  


#Task 1: Formatting and organizing the dataset for regression slopes including:

veg_format_site <- veg %>%
  #Format the dataset from WIDE --> LONG by consolidating the site characteristics
  gather(Region:salinity,
         key = "site_variable",
         value = "Category") %>%
  #Filter out all sites with less than 3 years of data, must span at least 5 years
  select(-Vegetation_Zone)

glimpse(veg_format_site)

#Task 3: Calculate the slopes and p-values of individual linear regressions for each site across all 
# vegetation metrics and site characteristics.

# This is accomplished by nesting the data essentially by Site, Vegetation Zone, and Metric
regression_slopes_site <- veg_format_site %>%
  group_by(Reserve, SiteID, Category, site_variable) %>%
  mutate(SampleSize = n()) %>%
  ungroup() %>%
  group_by(Reserve, SiteID, Category, site_variable, SampleSize) %>%
  nest() %>%
  # Using map() function and broom package, a regression is run for each nested dataset.
  # The model summary is output, which provides slope and p-value of regression
  mutate(regression = map(.x = data,
                          ~summary(lm(EMI ~ Year,
                                      data = .x)) %>%
                            tidy())) %>%
  unnest(regression) %>%
  #Format the summary table by removing unwanted terms, keeping the "Year" term, and 
  # formatting the number of decimal points
  filter(term == 'Year') %>%
  select(-data) %>% 
  mutate(across(estimate:p.value, ~round(., 4))) %>%
  ungroup() %>%
  #There were a handful of sites - vegetation zones with "perfect" EMI regressions (aka EMI value never
  # changed and was typically 1 or 0), so manually assigning a p-value of 1, since they are perfect
  # horizontal lines by definition of a statistical regdression. 
  mutate(p.value = ifelse(is.na(p.value), 1, p.value)) %>%
  rename(slope = estimate,
         pvalue = p.value)


glimpse(regression_slopes_site)



#Task 4: Calculate the compilation analysis scores based on the p-scores and the slopes direction

# The compilation analysis score is based on two criteria: p-value and direction of slope
# The magnitude of the score is based on the p-value 
# The score is then multiplied by the direction of the score (i.e., negative slope = negative score)

regression_scores_site <- regression_slopes_site %>%
  mutate(score = ifelse(pvalue <= 0.05, 2, 
                        ifelse(pvalue > 0.20, 0, 1)),
         score = ifelse(slope < 0, score * -1, score))


glimpse(regression_scores_site)

write.csv(regression_scores_site,
          "Output Stats\\Compilation Analysis - EMI Scores & Slopes - Site-level.csv")

#Task 5: Re-format the compilation analysis dataset from LONG --> WIDE for easier viewing in Excel

regression_scores_format <- regression_scores_site %>%
  spread(key = site_variable,
         value = Category)

write.csv(regression_scores_site,
          "Output Stats\\Compilation Analysis - EMI Scores & Slopes - Site-level FORMATTED.csv")

#Task 6: Calculate the descriptive statistics of the EMI compilation analysis by site characteristic
score_stats_site <- regression_scores_site %>%
  group_by(site_variable, Category) %>%
  #Mean and standard error of the compilation analysis scores and slopes
  summarise(slope.m = mean(slope, na.rm = TRUE),
            slope.se = sd(slope, na.rm = TRUE)/sqrt(n()),
            
            score.m = mean(score, na.rm = TRUE),
            score.se = sd(score, na.rm = TRUE)/sqrt(n())) %>%
  mutate(across(slope.m:score.se, ~round(., 3)))

glimpse(score_stats_site)

write.csv(score_stats_site,
          "Output Stats\\Compilation Analysis - EMI Score Stats - Site-level.csv")



#Chapter 5: Graphing the Compilation Analysis by Site Characteristics ----------------------------

#Task 1: Graph the compilation analysis score by Region

#Formatting the EMI scores for visualization by Region
score_stats_region <- score_stats_site %>%
  filter(site_variable == "Region") %>%
  #Ordering the regions by East --> West
  mutate(Category = factor(Category,
                           levels = c("Northeast", "Mid-Atlantic", "Southeast", 
                                      "Gulf Coast", "West Coast")))
  
glimpse(score_stats_region)

#Graphing the EMI scores by Region

meta_region_graph <- ggplot(data = score_stats_region,
                                aes(x = Category, 
                                    y = score.m)) +
  geom_bar(aes(fill = Category),
           stat = 'identity', position = position_dodge(0.9),
           linewidth = 1.5, colour = "black") + 
  geom_errorbar(aes(x = Category,
                    ymax = score.m + score.se, 
                    ymin = score.m - score.se,
                    group = Category),
                position = position_dodge(0.9),
                width = 0.5, linewidth = 1.25) + 
  scale_y_continuous(limits = c(-1.1, 2.05),
                     breaks = seq(-1.0, 2.0, 0.5),
                     expand = c(0,0)) +
  labs(y = "EMI Score",
       x = "") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 20, colour = "black"),
    axis.text = element_text(size = 20, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 20))

meta_region_graph


ggsave(meta_region_graph,
      filename = "Output Figures\\Compilation Analysis Graph by Region - EMI.jpg",
       units = "in",
       height = 8, width = 12, dpi = 300, limitsize = FALSE)




#Task 2: Graph the compilation analysis score by Tidal Range

#Formatting the EMI scores for visualization by Tidal Range
score_stats_tides <- score_stats_site %>%
  filter(site_variable == "tidal_range") %>%
  #Order tidal regime by increasing tidal amplitude
  mutate(Category = factor(Category, 
                           levels = c("Microtidal", "Mesotidal", "Macrotidal")))

glimpse(score_stats_tides)


meta_tides_graph <- ggplot(data = score_stats_tides,
                            aes(x = Category, 
                                y = score.m)) +
  geom_bar(aes(fill = Category),
           stat = 'identity', position = position_dodge(0.9),
           linewidth = 1.5, colour = "black") + 
  geom_errorbar(aes(x = Category,
                    ymax = score.m + score.se, 
                    ymin = score.m - score.se,
                    group = Category),
                position = position_dodge(0.9),
                width = 0.5, linewidth = 1.25) + 
  scale_y_continuous(limits = c(-1.1, 1.6),
                     breaks = seq(-1.0, 1.5, 0.5),
                     expand = c(0,0)) +
  labs(y = "EMI Score",
       x = "") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 20, colour = "black"),
    axis.text = element_text(size = 20, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 20))

meta_tides_graph


ggsave(meta_tides_graph,
       filename = "Output Figures\\Compilation Analysis Graph by Tide Range - EMI.jpg",
       units = "in",
       height = 8, width = 12, dpi = 300, limitsize = FALSE)




#Task 3: Graph the compilation analysis score by Salinity Regime

#Formatting the EMI scores for visualization by Salinity Regime
score_stats_salinity <- score_stats_site %>%
  filter(site_variable == "salinity") %>%
  #Ordering by increasing salinity concentration
  mutate(Category = factor(Category,
                           levels = c("Oligohaline", "Mesohaline", "Polyhaline")))

glimpse(score_stats_salinity)


meta_salinity_graph <- ggplot(data = score_stats_salinity,
                           aes(x = Category, 
                               y = score.m)) +
  geom_bar(aes(fill = Category),
           stat = 'identity', position = position_dodge(0.9),
           linewidth = 1.5, colour = "black") + 
  geom_errorbar(aes(x = Category,
                    ymax = score.m + score.se, 
                    ymin = score.m - score.se,
                    group = Category),
                position = position_dodge(0.9),
                width = 0.5, linewidth = 1.25) + 
  scale_y_continuous(limits = c(-1.1, 1.6),
                     breaks = seq(-1.0, 1.5, 0.5),
                     expand = c(0,0)) +
  labs(y = "EMI Score",
       x = "") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 20, colour = "black"),
    axis.text = element_text(size = 20, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 20))
meta_salinity_graph


ggsave(meta_salinity_graph,
       filename = "Output Figures\\compilation analysis Graph by Salinity Regime - EMI.jpg",
       units = "in",
       height = 8, width = 12, dpi = 300, limitsize = FALSE)




#Task 4: Graph the compilation analysis score by Geomorphology

#Formatting the EMI scores for visualization by Geomorphology
score_stats_geomorphology <- score_stats_site %>%
  filter(site_variable == "Geomorphology") %>%
  #Ordering by alphabetical order
  mutate(Category == factor(Category, 
                            levels = c("Back Barrier", "Bay Front", "Riverine")))

glimpse(score_stats_geomorphology)


meta_geomorphology_graph <- ggplot(data = score_stats_geomorphology,
                              aes(x = Category, 
                                  y = score.m)) +
  geom_bar(aes(fill = Category),
           stat = 'identity', position = position_dodge(0.9),
           linewidth = 1.5, colour = "black") + 
  geom_errorbar(aes(x = Category,
                    ymax = score.m + score.se, 
                    ymin = score.m - score.se,
                    group = Category),
                position = position_dodge(0.9),
                width = 0.5, linewidth = 1.25) + 
  scale_y_continuous(limits = c(-1.1, 1.6),
                     breaks = seq(-1.0, 1.5, 0.5),
                     expand = c(0,0)) +
  labs(y = "EMI Score",
       x = "") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 20, colour = "black"),
    axis.text = element_text(size = 20, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 20))


meta_geomorphology_graph


ggsave(meta_geomorphology_graph,
       filename = "Output Figures\\compilation analysis Graph by Geomorphology - EMI.jpg",
       units = "in",
       height = 8, width = 12, dpi = 300, limitsize = FALSE)




#Task 5: Combine the Bar Graphs into One Graphic

meta_bargraph_compiled <- (meta_region_graph + meta_geomorphology_graph) / (meta_salinity_graph + meta_tides_graph)

meta_bargraph_compiled

ggsave(meta_bargraph_compiled,
       filename = "Output Figures\\compilation analysis Graph Compiled by Metrics - EMI.jpg",
       units = "in",
       height = 12, width = 18, dpi = 300, limitsize = FALSE)

