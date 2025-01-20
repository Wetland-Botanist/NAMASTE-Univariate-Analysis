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
  rename(EMI = EIR)

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
#Chapter 3: Compilation Analysis with Vegetation Zone as an Interaction (for pvalues)
#-----------------------------------------------------------------------------------

#To calculate the p-value of the compilation analysis, regressions will be run at the site-level
# with models of: Vegetation Metric ~ Year + Veg_Zone * Veg_Zone

#However, there are a handful of sites with only 1 vegetation zone, which is not compatible with
# an ANCOVA analysis. Therefore, those sites will be run with a model of: Metric ~ Year. A
# simple linear regression. 

#Then, the p-values will be counted for each site. In the interaction model, a site will be
# assigned a score of 1 if Year or Interaction terms are significant. In the linear regression models,
# a site will be assigned a score of 1 if the Year term is significant.

#After all the p-value scores are assigned, the two datasets will be row combined and then 
# the percent of sites with significant p-values will be calculated and visualized



#Task 1: Calculate the p-values of the regressions of sites with 2 or more vegetation zones
# The model is Year * Vegetation Zone, which will provide three instances for 

# This is accomplished by nesting the data essentially by Site, Vegetation Zone, and Metric
regression_multi_zones <- veg_format %>%
  group_by(Reserve, SiteID, Metric) %>%
  mutate(Veg_Zone_Count = length(unique(Vegetation_Zone)),
         SampleSize = n()) %>%
  ungroup() %>%
  filter(Veg_Zone_Count > 1) %>%
  group_by(Reserve, SiteID, Metric, SampleSize, Region, salinity, tidal_range, Geomorphology) %>%
  nest() %>%
  mutate(regression = map(.x = data,
                          ~anova(lm(Value ~ Year * Vegetation_Zone,
                                      data = .x)) %>%
                            tidy())) %>%
  unnest(regression) %>%
  #Format the summary table by removing unwanted terms, keeping the "Year" term, and 
  # formatting the number of decimal points
  select(-data) %>%
  mutate(across(sumsq:p.value, ~round(., 4))) %>%
  ungroup() %>%
  #There were a handful of sites - vegetation zones with "perfect" regressions (aka values never
  # changed and was typically 1 or 0), so manually assigning a p-value of 1, since they are perfect
  # horizontal lines by definition of a statistical regression. 
  mutate(p.value = ifelse(is.na(p.value), 1, p.value)) %>%
  rename(pvalue = p.value)


glimpse(regression_multi_zones)

write.csv(regression_multi_zones,
          "Output Stats\\Site - Zone Compilation - Veg Zone 2 or more - ANOVA Tables.csv")


#Determine the number of significant p-values for each site based on the vegetation zone interaction
# regressions

pvalue_multi_zones <- regression_multi_zones %>%
  filter(term != "Vegetation_Zone",
         term != "Residuals") %>%
  group_by(Reserve, SiteID, Metric, Region, salinity, tidal_range, Geomorphology) %>%
  summarise(pvalue_count = sum(pvalue <= 0.05)) %>%
  ungroup() %>%
  mutate(pvalue_check = ifelse(pvalue_count > 0, 1, 0))


#Task 2: Calculate the p-values of the regressions of sites with only 1 site using simple linear 
# regressions

# This is accomplished by nesting the data essentially by Site, Vegetation Zone, and Metric
regression_single_zone <- veg_format %>%
  group_by(Reserve, SiteID, Metric) %>%
  mutate(Veg_Zone_Count = length(unique(Vegetation_Zone)),
         SampleSize = n()) %>%
  ungroup() %>%
  filter(Veg_Zone_Count == 1) %>%
  group_by(Reserve, SiteID, Metric, SampleSize, Region, salinity, tidal_range, Geomorphology) %>%
  nest() %>%
  mutate(regression = map(.x = data,
                          ~anova(lm(Value ~ Year,
                                    data = .x)) %>%
                            tidy())) %>%
  unnest(regression) %>%
  #Format the summary table by removing unwanted terms, keeping the "Year" term, and 
  # formatting the number of decimal points
  select(-data) %>%
mutate(across(sumsq:p.value, ~round(., 4))) %>%
  ungroup() %>%
  #There were a handful of sites - vegetation zones with "perfect" regressions (aka values never
  # changed and was typically 1 or 0), so manually assigning a p-value of 1, since they are perfect
  # horizontal lines by definition of a statistical regression. 
  mutate(p.value = ifelse(is.na(p.value), 1, p.value)) %>%
  rename(pvalue = p.value)


glimpse(regression_single_zone)

write.csv(regression_single_zone,
          "Output Stats\\Site - Zone Compilation - Veg Zone 1 - ANOVA Tables.csv")


#Determine the number of significant p-values for each site based on the vegetation zone interaction
# regressions

pvalue_single_zone <- regression_single_zone %>%
  filter(term != "Residuals") %>%
  group_by(Reserve, SiteID, Metric, Region, salinity, tidal_range, Geomorphology) %>%
  summarise(pvalue_count = sum(pvalue <= 0.05)) %>%
  ungroup() %>%
  mutate(pvalue_check = ifelse(pvalue_count > 0, 1, 0))


#Task 3: Combine the two pvalue score datasets of the single and multi vegetation zone models
# into one dataset using row combine

pvalue_scores <- rbind(pvalue_multi_zones, pvalue_single_zone)

glimpse(pvalue_scores)


write.csv(pvalue_scores,
          "Output Stats\\Site - Zone Compilation - pvalue count.csv")


#Task 3: Calculate the percent of significant p-values for sites across metrics and site categories

# The compilation analysis for vegetation metrics, not EMI, are calculated differently:
# (2) p-values are averaged +/- standard for each vegetation zone

regression_stats_zone <- pvalue_scores %>%
  gather(Region:Geomorphology,
         key = site_variable,
         value = Category) %>%
  group_by(Metric, site_variable, value = Category) %>%
  summarise(pvalue.percent = (length(pvalue_check[pvalue_check == 1]) / n()) * 100,
            SampleSize = n()) %>%
  ungroup() %>%
  mutate(pvalue.percent = round(pvalue.percent, 3))

glimpse(regression_stats_zone)


write.csv(regression_stats_zone,
          "Output Stats\\Site - Zone Compilation - pvalue summary statistics.csv")


#Task 4: Graph the p-value and slopes of the metrics besides EMI

#For the following graphs, the code is written to be a 'plug n play' by editing 
# the metrics filtered out in the first chunk of code. The user can then edit the 
# larger graph to their own purposes / aesthetics. 


# Filter out the vegetation metrics not to be contained in the pannelled graph

pvalue_vegzone <- regression_stats_zone %>%
  filter(Metric != "Abiotic Cover",
         site_variable == "Geomorphology") %>%
  mutate(value = factor(value, levels = c("Back Barrier", "Bay Front", "Riverine")))


#Graph of the percent of significant regressions based on p-value < 0.05

pvalue_vegzone_graph <- ggplot(data = pvalue_vegzone,
                               aes(x = value,
                                   y = pvalue.percent)) + 
  geom_bar(aes(fill = value),
                    stat = 'identity', position = position_dodge(0.9),
                    linewidth = 1.5, colour = "black") + 
  geom_text(aes(label = SampleSize),
            stat = 'identity', vjust = -0.75,
            size = 6, fontface = "bold") +
  labs(y = "Percent Regressions Significant (p < 0.05)",
       x = "") +
  scale_y_continuous(limits = c(0, 102),
                     breaks = seq(0, 100, 20),
                     expand = c(0,0)) + 
  theme_bw() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 16, colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 16, colour = "black"),
    axis.text =  element_text(size = 16, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 16)) +
facet_wrap(~Metric,
           nrow = 2, ncol = 4)

pvalue_vegzone_graph

ggsave(pvalue_vegzone_graph,
       filename = "Output Figures\\Site Level Regressions - geomorphology - tides - label.jpg",
       units = "in",
       height = 10, width = 18, dpi = 300, limitsize = FALSE)



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
          "Output Stats\\Site Compilation - ANOVA Tables.csv")

slopes_site_stats <- regression_slopes_site %>%
  filter(site_variable == "Region") %>%
  group_by(Metric) %>%
  summarise(slope.mean = mean(slope, na.rm = TRUE),
            slope.se = sd(slope, na.rm = TRUE)/sqrt(n())) %>%
  mutate(across(slope.mean:slope.se, ~round(., 3)))
  

write.csv(slopes_site_stats,
          "Output Stats\\Site Compilation - National Metric - Slope Summary Statistics.csv")

#Task 5: Calculate descriptive statistics for the slope and score by site characteristic

# The compilation analysis for vegetation metrics, not EMI, are calculated differently:
# (1) slopes of the regressions are averaged +/- standard error for each vegetation zone
# (2) p-values are averaged +/- standard for each vegetation zone

regression_stats_site <- regression_slopes_site %>%
  group_by(site_variable, Category, Metric) %>%
  summarise(slope.m = mean(slope, na.rm = TRUE),
            slope.se = sd(slope, na.rm = TRUE)/sqrt(n())) %>%
  ungroup() %>%
  mutate(across(slope.m:pvalue.score, ~round(., 3)))

glimpse(regression_stats_site)


write.csv(regression_stats_site,
          "Output Stats\\Site Compilation - Site Characteristic - Slope Summary Statistics.csv")




#Task 6: Graph the slopes of vegetation metrics

#Graph of the descriptive statistics of the slopes

#For the following graphs, the code is written to be a 'plug n play' by editing 
# the metrics filtered out in the first chunk of code. The user can then edit the 
# larger graph to their own purposes / aesthetics. 

# Filter out the vegetation metrics not to be contained in the visual cover 
# pannelled graph

pvalue_site <- regression_stats_site %>%
  filter(site_variable == "tidal_range") %>%
  mutate(Category = factor(Category, levels = c("Microtidal", "Mesotidal", "Macrotidal")))


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
  labs(y = "Slopes",
       x = "") +
  scale_y_continuous(limits = c(-3, 2),
                     breaks = seq(-3.0, 2, 1.0)) + 
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
       filename = "Output Figures\\Site Level Regressions - Cover Slopes - tides.jpg",
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
  labs(y = "Slopes",
       x = "") +
  scale_y_continuous(limits = c(-0.10, 0.1),
                     breaks = seq(-0.10, 0.1, 0.05)) + 
  theme_bw() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 18, colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text.x = element_text(size = 18, colour = "black", angle = 25, hjust = 1.0),
    axis.text.y = element_text(size = 18, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18)) +
  facet_wrap(~Metric,
             nrow = 2, ncol = 3)

slope_site_misc_graph

ggsave(slope_site_misc_graph,
       filename = "Output Figures\\Site Level Regressions - Misc Slopes - tides.jpg",
       units = "in",
       height = 10, width = 14, dpi = 300, limitsize = FALSE)

slope_site <- slope_site_cover_graph / slope_site_misc_graph

slope_site


ggsave(slope_site,
       filename = "Output Figures\\Site Level Regressions - Pannel Slopes - tides.jpg",
       units = "in",
       height = 10, width = 14, dpi = 300, limitsize = FALSE)

