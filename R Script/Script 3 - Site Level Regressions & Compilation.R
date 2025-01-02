#Project: NAMASTE National Univariate Exploratory Analysis
#Script: Script 3 - Site & Vegetation Zone Level Regressions and Compilation Analysis for Vegetation Metrics

#Author: Grant McKown, Jackson Estuarine Laboratory, University of New Hampshire
#Contact: james.mckown@unh.edu

#Project Description: The NERR system collects vegetation data along permanent transects and plots
#   in salt marshes throughout the United States. Vegetation data was synthesized, formatted, and 
#   filtered (based on field monitoring requirements). Reserve, site, transect, and plot-level characteristics
#   were also collected including hydrology, landscape geomorphology, elevation, and climate. The 
#   series of scripts conducts an exploratory analysis of different vegetation response factors to 
#   a series of explanatory and grouping factors. 

#Script Description: Script 3 formats the plot-level vegetation dataset, conducts simple linear regressions 
# of vegetation metrics for each site and vegetation zone within each site. The p-values and slopes
# of the regressions are then averaged +/- standard error by vegetation zone and site characteristic.
# The descriptive statistics and regression model summaries are exported as CSV files and graphed. 


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

#---------------------------------------------------------------------------------------------
#Chapter 2: Format the vegetation dataset for broom - map functions to create numerous models
#---------------------------------------------------------------------------------------------

#Task 1: Import the national plot dataframe
veg <- read.csv("Formatted Datasets\\Veg Plot Dataframe Formatted.csv") %>%
  select(-X, -invasive_cover) %>%
  select(Reserve:Year, Region:salinity, abiotic_cover:salt_ratio) %>%
  rename(EMI = EIR) %>%
  mutate(across(abiotic_cover:salt_ratio,
                ~ifelse(is.na(.), 0, .)))

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
                                            "Salt Ratio", "Shannon-Weiner Diversity", "Richness")))

glimpse(veg_format)


#-----------------------------------------------------------------------------------
#Chapter 3: Compilation Analysis by Vegetation Zone
#-----------------------------------------------------------------------------------

#Task 1: Formatting and organizing the dataset for regression slopes including:

veg_zone <- veg_format %>%
# For data visualization down the road, assign factors to the Vegetation Zone
  mutate(Vegetation_Zone = factor(Vegetation_Zone,
                                  levels = c("Low", "Mid", "Up"))) %>%
# Reduce the formatted vegetation dataset to only focus on vegetation zone
  select(-c(Region, Geomorphology, tidal_range, salinity))

glimpse(veg_zone)


#Task 2: Calculate the slopes and p-values of individual linear regressions for each vegetation
# zone within each site across all vegetation metrics 

# This is accomplished by nesting the data essentially by Site, Vegetation Zone, and Metric
regression_slopes_zone <- veg_zone %>%
  group_by(Reserve, SiteID, Vegetation_Zone, Metric) %>%
  mutate(SampleSize = n()) %>%
  ungroup() %>%
  filter(SampleSize > 2) %>%
  group_by(Reserve, SiteID, Vegetation_Zone, Metric, SampleSize) %>%
  nest() %>%
  mutate(regression = map(.x = data,
                          ~summary(lm(Value ~ Year,
                                      data = .x)) %>%
                            tidy())) %>%
  unnest(regression) %>%
  #Format the summary table by removing unwanted terms, keeping the "Year" term, and 
  # formatting the number of decimal points
  filter(term == 'Year') %>%
  select(-data) %>% 
  mutate(across(estimate:p.value, ~round(., 4))) %>%
  ungroup() %>%
  #There were a handful of sites - vegetation zones with "perfect" regressions (aka values never
  # changed and was typically 1 or 0), so manually assigning a p-value of 1, since they are perfect
  # horizontal lines by definition of a statistical regression. 
  mutate(p.value = ifelse(is.na(p.value), 1, p.value)) %>%
  rename(slope = estimate,
         pvalue = p.value)


glimpse(regression_slopes_zone)

write.csv(regression_slopes_zone,
          "Output Stats\\Compilation Analysis - Veg Metrics Trends - Vegetation Zone.csv")

#Task 3: Calculate descriptive statistics for the slope and score by vegetation zone

# The compilation analysis for vegetation metrics, not EMI, are calculated differently:
# (1) slopes of the regressions are averaged +/- standard error for each vegetation zone
# (2) p-values are averaged +/- standard for each vegetation zone

regression_stats_zone <- regression_slopes_zone %>%
  group_by(Vegetation_Zone, Metric) %>%
  summarise(slope.m = mean(slope, na.rm = TRUE),
            slope.se = sd(slope, na.rm = TRUE)/sqrt(n()),
            #Calculates the % of sites that have a significant trend
            pvalue.score = (length(pvalue[pvalue <= 0.05]) / n()) * 100 ) %>%
  ungroup() %>%
  mutate(across(slope.m:pvalue.score, ~round(., 3)))

glimpse(regression_stats_zone)


write.csv(regression_stats_zone,
          "Output Stats\\Compilation Analysis - Veg Metrics Trend Stats - Veg Zone.csv")


#Task 4: Graph the p-value and slopes of the metrics besides EMI

#For the following graphs, the code is written to be a 'plug n play' by editing 
# the metrics filtered out in the first chunk of code. The user can then edit the 
# larger graph to their own purposes / aesthetics. 


# Filter out the vegetation metrics not to be contained in the pannelled graph

pvalue_vegzone <- regression_stats_zone %>%
  filter(Metric != "EMI",
           Metric != "Abiotic Cover")


#Graph of the percent of significant regressions based on p-value < 0.05

pvalue_vegzone_graph <- ggplot(data = pvalue_vegzone,
                               aes(x = Vegetation_Zone,
                                   y = pvalue.score,
                                   group = Vegetation_Zone)) + 
  geom_bar(aes(fill = Vegetation_Zone),
                    stat = 'identity', position = position_dodge(0.9),
                    linewidth = 1.5, colour = "black") +  
  labs(y = "Percent Regressions Significant (p < 0.05)",
       x = "Vegetation Zone") +
  #scale_y_continuous(limits = c(0, 100),
  #                   breaks = seq(0, 100, 20)) + 
  theme_bw() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 18, colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 18, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18)) +
facet_wrap(~Metric,
           nrow = 2, ncol = 3)

pvalue_vegzone_graph

ggsave(pvalue_vegzone_graph,
       filename = "Output Figures\\Site Level Regressions - significant pvalues - Veg Zone 40 x-axis.jpg",
       units = "in",
       height = 10, width = 14, dpi = 300, limitsize = FALSE)



#Graph of the descriptive statistics of the slopes

#For the following graphs, the code is written to be a 'plug n play' by editing 
# the metrics filtered out in the first chunk of code. The user can then edit the 
# larger graph to their own purposes / aesthetics. 

# Filter out the vegetation metrics not to be contained in the visual cover 
# pannelled graph
pvalue_vegzone_cover <- regression_stats_zone %>%
  filter(Metric == "Live Cover" |
           Metric == "Halophyte Cover" |
           Metric == "Freshwater Cover")


slope_vegzone_cover_graph <- ggplot(data = pvalue_vegzone_cover,
                               aes(x = Vegetation_Zone,
                                   y = slope.m,
                                   group = Vegetation_Zone)) + 
  geom_bar(aes(fill = Vegetation_Zone),
           stat = 'identity', position = position_dodge(0.9),
           linewidth = 1.5, colour = "black") +  
  geom_errorbar(aes(x = Vegetation_Zone,
                    ymax = slope.m + slope.se, 
                    ymin = slope.m - slope.se,
                    group = Vegetation_Zone),
                position = position_dodge(0.9),
                width = 0.5, linewidth = 1.25) + 
  labs(y = "Mean Slopes",
       x = "Vegetation Zone") +
  scale_y_continuous(limits = c(-1.5, 1.0),
                     breaks = seq(-1.5, 1.0, 0.5)) + 
  theme_bw() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 18, colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 18, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18)) +
  facet_wrap(~Metric,
             nrow = 2, ncol = 3)

slope_vegzone_cover_graph

ggsave(slope_vegzone_cover_graph,
       filename = "Output Figures\\Site Level Regressions - Cover Slopes - Veg Zone.jpg",
       units = "in",
       height = 10, width = 14, dpi = 300, limitsize = FALSE)


# Filter out the vegetation metrics not to be contained in the misc. metrics
# pannelled graph

pvalue_vegzone_misc <- regression_stats_zone %>%
  filter(Metric == "Richness" |
           Metric == "Salt Ratio" |
           Metric == "Shannon-Weiner Diversity")


slope_vegzone_misc_graph <- ggplot(data = pvalue_vegzone_misc,
                              aes(x = Vegetation_Zone,
                                  y = slope.m,
                                  group = Vegetation_Zone)) + 
  geom_bar(aes(fill = Vegetation_Zone),
           stat = 'identity', position = position_dodge(0.9),
           linewidth = 1.5, colour = "black") +  
  geom_errorbar(aes(x = Vegetation_Zone,
                    ymax = slope.m + slope.se, 
                    ymin = slope.m - slope.se,
                    group = Vegetation_Zone),
                position = position_dodge(0.9),
                width = 0.5, linewidth = 1.25) + 
  labs(y = "Mean Slopes",
       x = "Vegetation Zone") +
  scale_y_continuous(limits = c(-0.04, 0.10),
                     breaks = seq(-0.04, 0.10, 0.02)) + 
  theme_bw() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 18, colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 18, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18)) +
  facet_wrap(~Metric,
             nrow = 2, ncol = 3)

slope_vegzone_misc_graph

ggsave(slope_vegzone_misc_graph,
       filename = "Output Figures\\Site Level Regressions - Misc Slopes - Veg Zone.jpg",
       units = "in",
       height = 10, width = 14, dpi = 300, limitsize = FALSE)

slope_vegzone <- slope_vegzone_cover_graph / slope_vegzone_misc_graph

slope_vegzone


ggsave(slope_vegzone,
       filename = "Output Figures\\Site Level Regressions - Pannel Slopes - Veg Zone.jpg",
       units = "in",
       height = 10, width = 14, dpi = 300, limitsize = FALSE)

#-----------------------------------------------------------------------------------
#Chapter 4: Compilation Analysis by Site
#-----------------------------------------------------------------------------------

#Task 1: Formatting and organizing the dataset for regression slopes including:

veg_site <- veg_format %>%
  #Format the dataset from WIDE --> LONG by consolidating the site characteristics
  gather(Region:salinity,
         key = "site_variable",
         value = "Category") %>%
  #with less than 3 years of data, must span at least 5 years
  select(-c(Vegetation_Zone))

glimpse(veg_site)


#Task 4: Calculate the slopes and p-values of individual linear regressions for each site across all 
# vegetation metrics

# This is accomplished by nesting the data essentially by Site and Metric
regression_slopes_site <- veg_site %>%
  group_by(Reserve, SiteID, Metric) %>%
  mutate(SampleSize = n()) %>%
  ungroup() %>%
  filter(SampleSize > 2) %>%
  group_by(Reserve, SiteID, site_variable, Category, Metric, SampleSize) %>%
  nest() %>%
  mutate(regression = map(.x = data,
                          ~summary(lm(Value ~ Year,
                                      data = .x)) %>%
                            tidy())) %>%
  unnest(regression) %>%
  #Format the summary table by removing unwanted terms, keeping the "Year" term, and 
  # formatting the number of decimal points
  filter(term == 'Year') %>%
  select(-data) %>% 
  mutate(across(estimate:p.value, ~round(., 4))) %>%
  ungroup() %>%
  #There were a handful of sites - vegetation zones with "perfect" regressions (aka values never
  # changed and was typically 1 or 0), so manually assigning a p-value of 1, since they are perfect
  # horizontal lines by definition of a statistical regression. 
  mutate(p.value = ifelse(is.na(p.value), 1, p.value)) %>%
  rename(slope = estimate,
         pvalue = p.value)


glimpse(regression_slopes_site)

write.csv(regression_slopes_site,
          "Output Stats\\Compilation Analysis - Veg Metrics Trends - Site Level.csv")

#Task 5: Calculate descriptive statistics for the slope and score by site characteristic

# The compilation analysis for vegetation metrics, not EMI, are calculated differently:
# (1) slopes of the regressions are averaged +/- standard error for each vegetation zone
# (2) p-values are averaged +/- standard for each vegetation zone

regression_stats_site <- regression_slopes_site %>%
  group_by(site_variable, Category, Metric) %>%
  summarise(slope.m = mean(slope, na.rm = TRUE),
            slope.se = sd(slope, na.rm = TRUE)/sqrt(n()),
            #Calculates the % of sites that have a significant trend
            pvalue.score = (length(pvalue[pvalue <= 0.05]) / n()) * 100 ) %>%
  ungroup() %>%
  mutate(across(slope.m:pvalue.score, ~round(., 3)))

glimpse(regression_stats_site)


write.csv(regression_stats_site,
          "Output Stats\\Compilation Analysis - Veg Metrics Trend Stats - Site Level.csv")




#Task 6: Graph the p-value and slopes of vegetation metrics

#For the following graphs, the code is written to be a 'plug n play' by editing 
# the metrics filtered out in the first chunk of code. The user can then edit the 
# larger graph to their own purposes / aesthetics. 


# Filter out the vegetation metrics not to be contained in the pannelled graph

#Graph of the percent of significant regressions based on p-value < 0.05

# Filter out the site characteristics not to be contained in the pannelled graph

pvalue_site <- regression_stats_site %>%
  filter(Metric != "EMI",
         Metric != "Abiotic Cover",
         site_variable == "Geomorphology") %>%
  mutate(Category = factor(Category,
                           levels = c("Back Barrier", "Bay Front", "Riverine")))


pvalue_site_graph <- ggplot(data = pvalue_site,
                               aes(x = Category,
                                   y = pvalue.score,
                                   group = Category)) + 
  geom_bar(aes(fill = Category),
           stat = 'identity', position = position_dodge(0.9),
           linewidth = 1.5, colour = "black") +  
  labs(y = "Percent Regressions Significant (p < 0.05)",
       x = "") +
  scale_y_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, 20)) + 
  theme_bw() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 18, colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 18, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18)) +
  facet_wrap(~Metric,
             nrow = 2, ncol = 3)

pvalue_site_graph

ggsave(pvalue_site_graph,
       filename = "Output Figures\\Site Level Regressions - significant pvalues - Geomorph - 100 xaxis.jpg",
       units = "in",
       height = 10, width = 14, dpi = 300, limitsize = FALSE)



#Graph of the descriptive statistics of the slopes

#For the following graphs, the code is written to be a 'plug n play' by editing 
# the metrics filtered out in the first chunk of code. The user can then edit the 
# larger graph to their own purposes / aesthetics. 

# Filter out the vegetation metrics not to be contained in the visual cover 
# pannelled graph
pvalue_site_cover <- pvalue_site %>%
  filter(Metric == "Live Cover" |
           Metric == "Halophyte Cover" |
           Metric == "Freshwater Cover")


slope_site_cover_graph <- ggplot(data = pvalue_site_cover,
                                    aes(x = Category,
                                        y = slope.m,
                                        group = Category)) + 
  geom_bar(aes(fill = Category),
           stat = 'identity', position = position_dodge(0.9),
           linewidth = 1.5, colour = "black") +  
  geom_errorbar(aes(x = Category,
                    ymax = slope.m + slope.se, 
                    ymin = slope.m - slope.se,
                    group = Category),
                position = position_dodge(0.9),
                width = 0.5, linewidth = 1.25) + 
  labs(y = "Mean Slopes",
       x = "") +
  scale_y_continuous(limits = c(-3.5, 1.0),
                     breaks = seq(-4.0, 1.0, 1.0)) + 
  theme_bw() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 18, colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 18, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18)) +
  facet_wrap(~Metric,
             nrow = 2, ncol = 3)

slope_site_cover_graph

ggsave(slope_site_cover_graph,
       filename = "Output Figures\\Site Level Regressions - Cover Slopes - Geomoprh.jpg",
       units = "in",
       height = 10, width = 14, dpi = 300, limitsize = FALSE)




pvalue_site_misc <- pvalue_site %>%
  filter(Metric == "Salt Ratio" |
           Metric == "Richness" |
           Metric == "Shannon-Weiner Diversity")


slope_site_misc_graph <- ggplot(data = pvalue_site_misc,
                                 aes(x = Category,
                                     y = slope.m,
                                     group = Category)) + 
  geom_bar(aes(fill = Category),
           stat = 'identity', position = position_dodge(0.9),
           linewidth = 1.5, colour = "black") +  
  geom_errorbar(aes(x = Category,
                    ymax = slope.m + slope.se, 
                    ymin = slope.m - slope.se,
                    group = Category),
                position = position_dodge(0.9),
                width = 0.5, linewidth = 1.25) + 
  labs(y = "Mean Slopes",
       x = "") +
  scale_y_continuous(limits = c(-0.10, 0.15),
                     breaks = seq(-0.10, 0.15, 0.05)) + 
  theme_bw() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 18, colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 18, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18)) +
  facet_wrap(~Metric,
             nrow = 2, ncol = 3)

slope_site_misc_graph

ggsave(slope_site_misc_graph,
       filename = "Output Figures\\Site Level Regressions - Misc Slopes - Geomorph.jpg",
       units = "in",
       height = 10, width = 14, dpi = 300, limitsize = FALSE)

slope_site <- slope_site_cover_graph / slope_site_misc_graph

slope_site


ggsave(slope_site,
       filename = "Output Figures\\Site Level Regressions - Pannel Slopes - Geomorph.jpg",
       units = "in",
       height = 10, width = 14, dpi = 300, limitsize = FALSE)

