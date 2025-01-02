#Project: NAMASTE National Univariate Exploratory Analysis
#Script: Script 2 - Meta-analysis of Site-level Linear Regressions of Plot Data

#Author: Grant McKown, Jackson Estuarine Laboratory, University of New Hampshire
#Contact: james.mckown@unh.edu

#Project Description: The NERR system collects vegetation data along permanent transects and plots
#   in salt marshes throughout the United States. Vegetation data was synthesized, formatted, and 
#   filtered (based on field monitoring requirements). Reserve, site, transect, and plot-level characteristics
#   were also collected including hydrology, landscape geomorphology, elevation, and climate. The 
#   series of scripts conducts an exploratory analysis of different vegetation response factors to 
#   a series of explanatory and grouping factors. 

#Script Description: Script 2 formats the plot-level vegetation dataset, conducts linear regressions of 
#                    vegetation metrics for each site, and then assigning meta-analysis scores based on 
#                    slopes and regression p-values.

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

#---------------------------------------------------------------------------------------------
#Chapter 2: Format the vegetation dataset for broom - map functions to create numerous models
#---------------------------------------------------------------------------------------------

veg <- read.csv("Formatted Datasets\\Veg Plot Dataframe Formatted.csv") %>%
  select(-X) %>%
  select(Reserve:Year, EIR, Region:salinity)

glimpse(veg)


#-----------------------------------------------------------------------------
#Chapter 3: Calculate live cover slope of each individual site, vegetation zone
#-----------------------------------------------------------------------------

# Formatting and organizing the dataset for regression slopes including:
  # Gathering the Site Variables and gathering the vegetation metrics, thus transforming the dataset
  # into a very LONG dataset. Essentially, breaking down each row into ~30 rows. This will allow
  # for simultaneous linear regressions for each vegetation metric for each site variable!

veg_format <- veg %>%
  gather(Region:salinity, 
         key = "Site_variable", 
         value = "Category") %>%
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
                                       "Freshwater Cover", "Invasive Cover", "EIR", 
                                       "Salt Ratio", "Shannon-Weiner Diversity", "Richness")),
         Vegetation_Zone = factor(Vegetation_Zone,
                                  levels = c("Low", "Mid", "Up")))

glimpse(veg_format)

#Calculate the slopes and p-values of individual linear regressions for each site across all 
# vegetation metrics and vegetation zones.

# This is accomplished by nesting the data essentially by Site, Vegetation Zone, and Metric

regression_slopes <- veg_format %>%
  group_by(Reserve, SiteID, Vegetation_Zone, Site_variable, Category, Metric) %>%
#Filter out all sites with less than 3 years of data, must span at least 5 years
  filter(length(Year) > 2) %>%
  filter(Year[length(Year)] - Year[1] >= 5) %>%
  nest() %>%
  mutate(regression = map(.x = data,
                                    ~summary(lm(Value ~ Year,
                                                  data = .x)) %>%
                                      tidy())) %>%
  unnest(regression) %>%
  filter(term == 'Year') %>%
  select(-std.error, -statistic, -term, -data) %>%
  mutate(estimate = round(estimate, 3),
         p.value = round(p.value, 3)) %>%
  ungroup() 

glimpse(regression_slopes)

#Third, calculate the meta-analysis scores based on the p-scores and the slopes direction

regression_scores <- regression_slopes %>%
         mutate(score = ifelse(p.value <= 0.05, 2, 
                        ifelse(p.value > 0.20, 0, 1)),
         score = ifelse(estimate < 0, score * -1, score)) %>%
  rename(slope = estimate,
         pvalue = p.value) %>%
  filter(!is.na(score)) %>%
  filter(!is.na(slope))
  
glimpse(regression_scores)

write.csv(regression_scores,
          "Output Stats\\Site Linear Regressions & Meta-analysis Scores.csv")


regression_scores_readable <- regression_scores %>%
  spread(key = Site_variable, value = Category)

glimpse(regression_scores_readable)


write.csv(regression_scores_readable,
          "Output Stats\\Site Linear Regressions & Meta-analysis Scores - FORMATTEd.csv")

#-----------------------------------------------------------------------------------------------
#Chapter 4: Calculate the mean and standard error of the slope and scores across vegetation zone
#-----------------------------------------------------------------------------------------------

score_stats_veg <- regression_scores %>%
  group_by(Vegetation_Zone, Metric) %>%
  summarise(slope.m = mean(slope, na.rm = TRUE),
            slope.se = sd(slope, na.rm = TRUE)/sqrt(n()),
            
            score.m = mean(score, na.rm = TRUE),
            score.se = sd(score, na.rm = TRUE)/sqrt(n())) %>%
  mutate(across(slope.m:score.se, ~round(., 3))) %>%
  mutate(Vegetation_Zone = factor(Vegetation_Zone,
                                  levels = c("Low", "Mid", "Up")))

glimpse(score_stats_veg)


write.csv(score_stats_veg,
          "Output Stats\\Meta-analysis Score Stats by Veg Zone.csv")





#-------------------------------------------------------------
#Chapter 5: Graph the Meta-analysis Scores by Vegetation Zone
#--------------------------------------------------------------


# Nine-pannelled graph of the meta-analysis score descriptive statistics with each panel a vegetaiton metric

score_stats_veg_mod <- score_stats_veg %>%
  filter(Metric == "EIR")
  


meta_veg_graph <- ggplot(data = score_stats_veg_mod,
                             aes(x = Vegetation_Zone, 
                                 y = score.m)) +
  geom_bar(aes(fill = Vegetation_Zone),
           stat = 'identity', position = position_dodge(0.9),
           linewidth = 1.25, colour = "black") + 
  geom_errorbar(aes(x = Vegetation_Zone,
                    ymax = score.m + score.se, 
                    ymin = score.m - score.se,
                    group = Vegetation_Zone),
                position = position_dodge(0.9),
                width = 0.25, linewidth = 1.25) + 
  scale_y_continuous(limits = c(-2.1, 2.1),
                     breaks = seq(-2.0, 2.0, 1.0),
                     expand = c(0,0)) +
  labs(y = "Score",
       x = "Vegetation Zone") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 18, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18))
 # facet_wrap(~Metric,
 #             ncol = 3, nrow = 3)

meta_veg_graph


ggsave(meta_veg_graph,
       filename = "Output Figures\\Meta-analysis Graph by Vegetation Zone - EIR X-Axis 2.jpg",
       units = "in",
       height = 10, width = 16, dpi = 300, limitsize = FALSE)






#-----------------------------------------------------------------------------
#Chapter 6: Calculate the Meta-analysis Score Descriptive Stats Across Site Characteristics
#-----------------------------------------------------------------------------

score_stats_site <- regression_scores %>%
  group_by(Site_variable, Category, Metric) %>%
  summarise(slope.m = mean(slope, na.rm = TRUE),
            slope.se = sd(slope, na.rm = TRUE)/sqrt(n()),
            
            score.m = mean(score, na.rm = TRUE),
            score.se = sd(score, na.rm = TRUE)/sqrt(n())) %>%
  mutate(across(slope.m:score.se, ~round(., 3)))
  
glimpse(score_stats_site)

write.csv(score_stats_site,
          "Output Stats\\Meta-analysis Score Stats by Site Characteristics.csv")





#-----------------------------------------------------------------------------
#Chapter 7: Graph the Meta-analysis by individual site characteristics by Bar Graph
#-----------------------------------------------------------------------------
#------------------------------------------------
#Task 1: Graph the meta-analysis score by Region
#------------------------------------------------


# Nine-pannelled graph of the meta-analysis score descriptive statistics with each panel a vegetation metric

score_stats_region <- score_stats_site %>%
  filter(Site_variable == "Region") %>%
  filter(Metric == "EIR") %>%
  mutate(Category = factor(Category,
                           levels = c("Northeast", "Mid-Atlantic", "Southeast", "Gulf Coast", "West Coast")))
  
glimpse(score_stats_region)


meta_region_graph <- ggplot(data = score_stats_region,
                                aes(x = Category, 
                                    y = score.m)) +
  geom_bar(aes(fill = Category),
           stat = 'identity', position = position_dodge(0.9),
           linewidth = 1.25, colour = "black") + 
  geom_errorbar(aes(x = Category,
                    ymax = score.m + score.se, 
                    ymin = score.m - score.se,
                    group = Category),
                position = position_dodge(0.9),
                width = 0.25, linewidth = 1.25) + 
  scale_y_continuous(limits = c(-2.1, 2.1),
                     breaks = seq(-2.0, 2.0, 1.0),
                     expand = c(0,0)) +
  labs(y = "Score",
       x = "Region") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 17, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18)) +
  facet_wrap(~Metric,
             ncol = 3, nrow = 3)

meta_region_graph


ggsave(meta_region_graph,
      filename = "Output Figures\\Meta-analysis Graph by Region - EIR & X-axis 2.jpg",
       units = "in",
       height = 10, width = 16, dpi = 300, limitsize = FALSE)



#------------------------------------------------------
#Task 2: Graph the meta-analysis score by Tidal Range
#------------------------------------------------------


# Nine-pannelled graph of the meta-analysis score descriptive statistics with each panel a vegetation metric

score_stats_tides <- score_stats_site %>%
  filter(Site_variable == "tidal_range") %>%
  filter(Metric == "EIR") %>%
  mutate(Category = factor(Category, 
                           levels = c("Microtidal", "Mesotidal", "Macrotidal")))

glimpse(score_stats_tides)


meta_tides_graph <- ggplot(data = score_stats_tides,
                            aes(x = Category, 
                                y = score.m)) +
  geom_bar(aes(fill = Category),
           stat = 'identity', position = position_dodge(0.9),
           linewidth = 1.25, colour = "black") + 
  geom_errorbar(aes(x = Category,
                    ymax = score.m + score.se, 
                    ymin = score.m - score.se,
                    group = Category),
                position = position_dodge(0.9),
                width = 0.25, linewidth = 1.25) + 
  scale_y_continuous(limits = c(-1.1, 1.1),
                     breaks = seq(-2.0, 2.0, 1.0),
                     expand = c(0,0)) +
  labs(y = "Score",
       x = "Tide Range") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 18, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18))
 # facet_wrap(~Metric,
#             ncol = 3, nrow = 3)

meta_tides_graph


ggsave(meta_tides_graph,
       filename = "Output Figures\\Meta-analysis Graph by Tide Range - EIR & X-axis 1.jpg",
       units = "in",
       height = 10, width = 16, dpi = 300, limitsize = FALSE)





#---------------------------------------------------------
#Task 3: Graph the meta-analysis score by Salinity Regime
#---------------------------------------------------------


# Nine-pannelled graph of the meta-analysis score descriptive statistics with each panel a vegetation metric

score_stats_salinity <- score_stats_site %>%
  filter(Site_variable == "salinity") %>%
  filter(Metric == "EIR") %>%
  mutate(Category = factor(Category,
                           levels = c("Oligohaline", "Mesohaline", "Polyhaline")))

glimpse(score_stats_salinity)


meta_salinity_graph <- ggplot(data = score_stats_salinity,
                           aes(x = Category, 
                               y = score.m)) +
  geom_bar(aes(fill = Category),
           stat = 'identity', position = position_dodge(0.9),
           linewidth = 1.25, colour = "black") + 
  geom_errorbar(aes(x = Category,
                    ymax = score.m + score.se, 
                    ymin = score.m - score.se,
                    group = Category),
                position = position_dodge(0.9),
                width = 0.25, linewidth = 1.25) + 
  scale_y_continuous(limits = c(-1.1, 1.1),
                     breaks = seq(-2.0, 2.0, 1.0),
                     expand = c(0,0)) +
  labs(y = "Score",
       x = "Salinity Regime") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 18, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18))
  #facet_wrap(~Metric,
   #          ncol = 3, nrow = 3)

meta_salinity_graph


ggsave(meta_salinity_graph,
       filename = "Output Figures\\Meta-analysis Graph by Salinity Regime - EIR & X-axis 1.jpg",
       units = "in",
       height = 12, width = 18, dpi = 300, limitsize = FALSE)





#---------------------------------------------------------
#Task 4: Graph the meta-analysis score by Geomorphology
#---------------------------------------------------------


# Nine-pannelled graph of the meta-analysis score descriptive statistics with each panel a vegetation metric

score_stats_geomorphology <- score_stats_site %>%
  filter(Site_variable == "Geomorphology") %>%
  filter(Metric == "EIR") %>%
  mutate(Category == factor(Category, 
                            levels = c("Back Barrier", "Bay Front", "Riverine")))

glimpse(score_stats_geomorphology)


meta_geomorphology_graph <- ggplot(data = score_stats_geomorphology,
                              aes(x = Category, 
                                  y = score.m)) +
  geom_bar(aes(fill = Category),
           stat = 'identity', position = position_dodge(0.9),
           linewidth = 1.25, colour = "black") + 
  geom_errorbar(aes(x = Category,
                    ymax = score.m + score.se, 
                    ymin = score.m - score.se,
                    group = Category),
                position = position_dodge(0.9),
                width = 0.25, linewidth = 1.25) + 
  scale_y_continuous(limits = c(-1.1, 1.1),
                     breaks = seq(-2.0, 2.0, 1.0),
                     expand = c(0,0)) +
  labs(y = "Score",
       x = "Geomorphology") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 18, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18))
  #facet_wrap(~Metric,
  #           ncol = 3, nrow = 3)

meta_geomorphology_graph


ggsave(meta_geomorphology_graph,
       filename = "Output Figures\\Meta-analysis Graph by Geomorphology - EIR & X-axis 1.jpg",
       units = "in",
       height = 10, width = 16, dpi = 300, limitsize = FALSE)




#-------------------------------------------------------------------------------
#Chapter 8: Graph the Meta-analysis by individual site characteristics by boxplot
#-----------------------------------------------------------------------------

#For Box Plots, we need to use the non-summarized data, Regression_Scores


#------------------------------------------------
#Task 1: Graph the meta-analysis score by Region
#------------------------------------------------


# Nine-pannelled graph of the meta-analysis score descriptive statistics with each panel a vegetation metric

regression_scores_region <- regression_scores %>%
  filter(Site_variable == "Region") %>%
  filter(Metric == "EIR") %>%
  mutate(Category = factor(Category,
                           levels = c("Northeast", "Mid-Atlantic", "Southeast", "Gulf Coast", "West Coast")))


glimpse(regression_scores_region)


meta_region_boxplot <- ggplot(data = regression_scores_region,
                            aes(x = Category, 
                                y = score,
                                group = Category)) +
  geom_boxplot(aes(fill = Category),
               linewidth = 1.10) +
  scale_y_continuous(limits = c(-2.1, 2.1),
                     breaks = seq(-2.0, 2.0, 1.0),
                     expand = c(0,0)) +
  labs(y = "Score",
       x = "Region") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 18, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18))
#  facet_wrap(~Metric,
#             ncol = 3, nrow = 3)

meta_region_boxplot


ggsave(meta_region_boxplot,
       filename = "Output Figures\\Meta-analysis Graph by Region Boxplot - EIR & X-axis 2.jpg",
       units = "in",
       height = 10, width = 16, dpi = 300, limitsize = FALSE)



#------------------------------------------------------
#Task 2: Graph the meta-analysis score by Tidal Range
#------------------------------------------------------


# Nine-pannelled graph of the meta-analysis score descriptive statistics with each panel a vegetation metric

regression_score_tides <- regression_scores %>%
  filter(Site_variable == "tidal_range") %>%
  filter(Metric == "EIR") %>%
  mutate(Category = factor(Category, 
                           levels = c("Microtidal", "Mesotidal", "Macrotidal")))

glimpse(regression_score_tides)


meta_tides_boxplot <- ggplot(data = regression_score_tides,
                              aes(x = Category, 
                                  y = score,
                                  group = Category)) +
  geom_boxplot(aes(fill = Category),
               linewidth = 1.10) +
  scale_y_continuous(limits = c(-2.1, 2.1),
                     breaks = seq(-2.0, 2.0, 1.0),
                     expand = c(0,0)) +
  labs(y = "Score",
       x = "Tide Range") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 18, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18))
  #facet_wrap(~Metric,
  #           ncol = 3, nrow = 3)


meta_tides_boxplot


ggsave(meta_tides_boxplot,
       filename = "Output Figures\\Meta-analysis Graph by Tide Range Boxplot - EIR & X-axis 2.jpg",
       units = "in",
       height = 12, width = 18, dpi = 300, limitsize = FALSE)





#---------------------------------------------------------
#Task 3: Graph the meta-analysis score by Salinity Regime
#---------------------------------------------------------


# Nine-pannelled graph of the meta-analysis score descriptive statistics with each panel a vegetation metric

regression_score_salinity <- regression_scores %>%
  filter(Site_variable == "salinity") %>%
  filter(Metric == "EIR") %>%
  mutate(Category = factor(Category,
                           levels = c("Oligohaline", "Mesohaline", "Polyhaline")))

glimpse(regression_score_salinity)


meta_salinity_boxplot <- ggplot(data = regression_score_salinity,
                             aes(x = Category, 
                                 y = score,
                                 group = Category)) +
  geom_boxplot(aes(fill = Category),
               linewidth = 1.10) +
  scale_y_continuous(limits = c(-2.1, 2.1),
                     breaks = seq(-2.0, 2.0, 1.0),
                     expand = c(0,0)) +
  labs(y = "Score",
       x = "Salinity Regime") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 18, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18))
#  facet_wrap(~Metric,
#             ncol = 3, nrow = 3)


meta_salinity_boxplot



ggsave(meta_salinity_boxplot,
       filename = "Output Figures\\Meta-analysis Graph by Salinity Regime by Boxplot - EIR & X-axis 2.jpg",
       units = "in",
       height = 10, width = 16, dpi = 300, limitsize = FALSE)





#---------------------------------------------------------
#Task 4: Graph the meta-analysis score by Geomorphology
#---------------------------------------------------------


# Nine-pannelled graph of the meta-analysis score descriptive statistics with each panel a vegetation metric

regression_score_geomorph <- regression_scores %>%
  filter(Site_variable == "Geomorphology") %>%
  filter(Metric == "EIR") %>%
  mutate(Category == factor(Category, 
                            levels = c("Back Barrier", "Bay Front", "Riverine")))


glimpse(regression_score_geomorph)


meta_geomorphology_boxplot <- ggplot(data = regression_score_geomorph,
                                   aes(x = Category, 
                                       y = score,
                                       group = Category)) +
  geom_boxplot(aes(fill = Category),
               linewidth = 1.10) +
  scale_y_continuous(limits = c(-2.1, 2.1),
                     breaks = seq(-2.0, 2.0, 1.0),
                     expand = c(0,0)) +
  labs(y = "Score",
       x = "Geomorphology") +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 18, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18))
  #facet_wrap(~Metric,
   #          ncol = 3, nrow = 3)

meta_geomorphology_boxplot


ggsave(meta_geomorphology_boxplot,
       filename = "Output Figures\\Meta-analysis Graph by Geomorphology by Boxplot - EIR & X-axis 2.jpg",
       units = "in",
       height = 10, width = 16, dpi = 300, limitsize = FALSE)





#-------------------------------------------------------
#Chapter 8: Super Cool Meta-analysis Graph Across Sites
#-------------------------------------------------------

#Task 1: Prep the regression_slopes dataset

regression_formatted <- regression_scores %>%
  mutate(slope_sig = ifelse(pvalue <= 0.05, "yes", "no"))

glimpse(regression_formatted)

#Task 2: Graph the Super cool Meta-analysis Graph for all nine vegetation metrics

# Create a 9 pannelled graph

ggplot(aes(x = SiteID, 
           y = slope)) +  # use forcats to turn that thing that's nicely in order into a factor
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high,
                      col = sig_trend)) +  # geom_pointrange for the point + CIs; color by trend category (this is where Namaste could have 3)
  geom_hline(yintercept = 0) +
  labs(title = "Slopes for all stations",
       subtitle = "slopes are in [parameter's measurement units] / year",
       x = "Station",
       col = "Significant slope?") +
  theme(axis.text.x = element_blank(), 
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.position = "bottom",   # move the legend around
        legend.justification = "left",
        legend.direction = "horizontal") +
facet_wrap(~parameter, 
           scales = "free") +  # wrap it by parameter (or region, whatever you want)


