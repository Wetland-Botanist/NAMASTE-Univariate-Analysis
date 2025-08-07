#Project: NAMASTE National Univariate Exploratory Analysis
#Script: Script 7 - National Regressions of Explanatory Factors of Site-level Slopes 

#Author: Grant McKown, Jackson Estuarine Laboratory, University of New Hampshire
#Contact: james.mckown@unh.edu

#Project Description: The NERR system collects vegetation data along permanent transects and plots
#   in salt marshes throughout the United States. Vegetation data was synthesized, formatted, and 
#   filtered (based on field monitoring requirements). Reserve, site, transect, and plot-level characteristics
#   were also collected including hydrology, landscape geomorphology, elevation, and climate. The 
#   series of scripts conducts an exploratory analysis of different vegetation response factors to 
#   a series of explanatory and grouping factors. 

# Script Description: Script 7 formats the calculated slope by site dataset, conducts 
#   national simple linear regressions of the slopes for each vegetation metric against
#   sea level rise in the past 19 years and NERR Landscape Condition. The 
#   slopes of the regressions are extracted from the model summaries, exported by CSV, 
#   and results are graphed. 



#Chapter 1: Import package library --------------------------------------------


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


#Chapter 2: Import and Format the vegetation dataset -----------------------------------


#Task 1 - Import the dataset and remove extraneous columns
veg <- read.csv("Formatted Datasets\\Veg Slope by Site Formatted.csv") %>%
  select(-X) %>%
  select(-invasive_cover)

glimpse(veg)


#Task 2 - Gather the vegetation metrics and rename the values

veg_format <- veg %>%
  gather(abiotic_cover:salt_ratio,
         key = "Metric",
         value = "Value") %>%
  mutate(Metric = ifelse(Metric == "abiotic_cover", "Abiotic Cover", 
                         ifelse(Metric == "live_cover", "Live Cover",
                                ifelse(Metric == "halophyte", "Halophyte Cover",
                                       ifelse(Metric == "SWdiv_slope", "Shannon-Weiner Diversity", 
                                              ifelse(Metric == "freshwater", "Freshwater Cover",
                                                     ifelse(Metric == "salt_ratio", "Salt Ratio", 
                                                            ifelse(Metric == "EMI_slope", "EMI",
                                                                   ifelse(Metric == "Richness_slope", "Richness",
                                                                          Metric))))))))) %>%
  mutate(Metric = factor(Metric, levels = c("Abiotic Cover", "Live Cover", "Halophyte Cover",
                                            "Freshwater Cover", "EMI", 
                                            "Salt Ratio", "Shannon-Weiner Diversity", "Richness"))) %>%
  rename(SLR = SLR_last19yrs,
         Landscape = NERRs_Landscape_resiliency_condition_sum_quantile)



#Chapter 3: Sea Level Rise Regressions for the slope dataset ------------------------------


#Task 1 - Run the national regressions for Sea Level Rise to vegetation metrics

slr_regression <- veg_format %>%
  filter(!is.na(Value)) %>%
  group_by(Metric) %>%
  mutate(SampleSize = length(Value)) %>%
  ungroup() %>%
  group_by(Metric, SampleSize) %>%
  nest() %>%
  mutate(regression = map(.x = data,
                          ~lm(Value ~ SLR,
                             data = .x) %>%
                            tidy()),
         rsquared = map(.x = data,
                            ~lm(Value ~ SLR,
                                data = .x) %>%
                              glance() %>%
                              select(r.squared, adj.r.squared))) %>%
  unnest(c(regression, rsquared)) %>%
  select(-data) %>%
  filter(term == "SLR") %>%
  mutate(across(estimate:adj.r.squared, ~round(., 3)),
         p.value = round(p.value, 4))

glimpse(slr_regression)


write.csv(slr_regression,
          "Output Stats\\National Regression of Slopes & SLR - Model Summary.csv")

#Task 2 - Graph the statistically appropriate regressions

#Instead of predicting the slopes and the values with ggpredict, will use the geom_stat function 
# in ggplot2, since the regressions are simple linear regressions

veg_slr <- veg_format %>%
  filter(Metric == "EMI")

slr_graph <- ggplot(aes(x = SLR,
                        y = Value),
                    data = veg_slr) +
  geom_point(aes(fill = Region),
             alpha = 0.75, size = 7, shape = 21) +
  geom_smooth(method = "lm",
              alpha = 0.5,
              colour = "black",
              linewidth = 1.5,
              linetype = "dashed") +
  scale_x_continuous(limits = c(-8, 12),
                     breaks = seq(-8, 12, 2)) + 
  labs(x = "Sea Level Rise Last 19 Years (mm / yr)",
       y = "EMI (per yr)") +
  theme_bw() +
  theme(
    legend.position = c(0.125, 0.875),
    legend.title = element_blank(),
    legend.text = element_text(size = 16, colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 18, colour = "black"))

slr_graph

ggsave(slr_graph,
       filename = "Output Figures\\National Regression - SLR & EMI.jpg",
       units = "in",
       height = 8, width = 12, dpi = 300, limitsize = FALSE)




#Chapter 4: Landscape Condition Regressions for the slope dataset --------------------------


#Task 1 - Run the national regressions for Sea Level Rise to vegetation metrics

landscape_regression <- veg_format %>%
  group_by(Metric) %>%
  nest() %>%
mutate(regression = map(.x = data,
                        ~lm(Value ~ Landscape,
                            data = .x) %>%
                          summary() %>%
                          tidy()),
       rsquared = map(.x = data,
                      ~lm(Value ~ Landscape,
                          data = .x) %>%
                        glance() %>%
                        select(r.squared, adj.r.squared))) %>%
  unnest(c(regression, rsquared)) %>%
  select(-data) %>%
  filter(term == "Landscape") %>%
  mutate(across(estimate:adj.r.squared, ~round(., 3)),
         p.value = round(p.value, 4))

glimpse(landscape_regression)

write.csv(landscape_regression,
          "Output Stats\\National Regression of Slopes & Landscape - Model Summary.csv")

#Task 2 - Graph the statistically appropriate regressions

#Instead of predicting the slopes and the values with ggpredict, will use the geom_stat function 
# in ggplot2, since the regressions are simple linear regressions

veg_landscape <- veg_format %>%
  filter(Metric == "Shannon-Weiner Diversity")

landscape_graph <- ggplot(aes(x = Landscape,
                        y = Value),
                    data = veg_landscape) +
  geom_point(aes(fill = Region),
             alpha = 0.75, size = 7, shape = 21) +
  geom_smooth(method = "lm",
              alpha = 0.5,
              colour = "black",
              linewidth = 1.5,
              linetype = "dashed") +
  scale_x_continuous(limits = c(1.5, 10.5),
                     breaks = seq(2, 10, 2)) + 
  labs(x = "NERR Landscape Condition",
       y = "Shannon-Weiner Diversity (per yr)") +
  theme_bw() +
  theme(
    legend.position = c(0.125, 0.875),
    legend.title = element_blank(),
    legend.text = element_text(size = 16, colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 18, colour = "black"))

landscape_graph

ggsave(landscape_graph,
       filename = "Output Figures\\National Regression - Landscape & Diversity.jpg",
       units = "in",
       height = 8, width = 12, dpi = 300, limitsize = FALSE)
