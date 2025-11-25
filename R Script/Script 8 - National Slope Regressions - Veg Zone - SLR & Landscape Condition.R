#Project: NAMASTE National Univariate Exploratory Analysis
#Script: Script 8 - National Regressions of Explanatory Factors of Vegetation Zone Slopes

#Author: Grant McKown, Jackson Estuarine Laboratory, University of New Hampshire
#Contact: james.mckown@unh.edu

#Project Description: The NERR system collects vegetation data along permanent transects and plots
#   in salt marshes throughout the United States. Vegetation data was synthesized, formatted, and 
#   filtered (based on field monitoring requirements). Reserve, site, transect, and plot-level characteristics
#   were also collected including hydrology, landscape geomorphology, elevation, and climate. The 
#   series of scripts conducts an exploratory analysis of different vegetation response factors to 
#   a series of explanatory and grouping factors. 

#Script Description: Script 8 formats the calculated slope by site dataset, conducts 
#   national simple linear regressions of the slopes for each vegetation metric against
#   sea level rise in the past 19 years and NERR Landscape Condition. The 
#   slopes of the regressions are extracted from the model summaries, exported by CSV, 
#   and results are graphed. 



#Chapter 1: Import package library --------------------------------------------


#Data Analysis Packages

library(dplyr)
library(tidyr)
library(rstatix)
library(broom)
library(purrr)
library(ggResidpanel)
library(ggeffects)

#Graphing and Data Visualization Packages

library(ggplot2)
library(patchwork)


#Chapter 2: Import and Format the vegetation dataset -----------------------------------


#Task 1 - Import the dataset and remove extraneous columns
veg <- read.csv("Formatted Datasets\\Veg Slope by Site and Zone Formatted.csv") %>%
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

glimpse(veg_format)


#Chapter 3: Sea Level Rise Regressions for the slope dataset ---------------------------------


#Task 1 - Run the national regressions for Sea Level Rise to vegetation metrics

slr_regression <- veg_format %>%
  filter(!is.na(Value)) %>%
  group_by(Metric) %>%
  mutate(SampleSize = length(Value)) %>%
  ungroup() %>%
  group_by(Metric, SampleSize) %>%
  nest() %>%
  mutate(regression = map(.x = data,
                          ~lm(Value ~ SLR * Vegetation_Zone,
                             data = .x) %>%
                            anova() %>%
                            tidy()),
    rsquared = map(.x = data,
                   ~lm(Value ~ SLR,
                       data = .x) %>%
                     glance() %>%
                     select(r.squared, adj.r.squared))) %>%
  unnest(c(regression, rsquared)) %>%
  select(-data) %>%
  mutate(across(sumsq:adj.r.squared, ~round(., 3)),
         p.value = round(p.value, 4))

glimpse(slr_regression)


write.csv(slr_regression,
          "Output Stats\\National Regression of Slopes & SLR - Veg Zone - ANOVA Table.csv")


#Task 2 - Predict the linear regressions of the statistically significant interactions

# Run a dplyr loop to predict the values of the regressions

regression_predicted <- veg_format %>%
  mutate(Vegetation_Zone = factor(Vegetation_Zone)) %>%
  group_by(Metric) %>%
  nest() %>%
  mutate(regression = map(.x = data,
                          ~lm(Value ~ SLR * Vegetation_Zone,
                                data = .x) %>%
                            ggpredict(., 
                                      terms = c("SLR [all]", "Vegetation_Zone [all]"),
                                      type = "fixed", interval = "confidence"))) %>%
  select(-data) %>%
  unnest(regression) %>%
  rename(SLR = x,
         Value_pred = predicted,
         Vegetation_Zone = group) %>%
  filter(Metric == "EMI" |
           Metric == "Freshwater Cover" |
           Metric == "Shannon-Weiner Diversity")

glimpse(regression_predicted)


#Calculate the slopes of each vegetation metric across Vegetation Zone

regression_slopes <- regression_predicted %>%
  group_by(Metric, Vegetation_Zone) %>%
  summarise(slope = (Value_pred[which(SLR == max(SLR))] - Value_pred[which(SLR == min(SLR))]) / ((max(SLR) - min(SLR)))) %>%
  ungroup() %>%
  mutate(slope = round(slope, 3))

write.csv(regression_slopes,
          "Output Stats\\National Slope Models - SLR Slopes Across All Zones - Veg Zones.csv")


#Task 2 - Graph the statistically appropriate regressions

#Instead of predicting the slopes and the values with ggpredict, will use the geom_stat function 
# in ggplot2, since the regressions are simple linear regressions

veg_slr <- veg_format %>%
  filter(Metric == "Shannon-Weiner Diversity")

regression_slr <- regression_predicted %>%
  filter(Metric == "Shannon-Weiner Diversity")

slr_graph <- ggplot() +
  geom_point(data = veg_slr,
             aes(x = SLR,
                 y = Value,
                 fill = Vegetation_Zone),
             alpha = 0.75, size = 7, shape = 21,
             position = "jitter") +
  geom_ribbon(data = regression_slr,
              aes(x = SLR,
                  ymax = Value_pred + std.error,
                  ymin = Value_pred - std.error,
                  fill = Vegetation_Zone),
              alpha = 0.35) + 
  geom_line(data = regression_slr,
            aes(x = SLR, y = Value_pred,
                colour = Vegetation_Zone),
            linetype = "dashed", size = 1.25) +
  scale_x_continuous(limits = c(-8, 12),
                     breaks = seq(-8, 12, 2)) + 
  labs(x = "Sea Level Rise Last 19 Years (mm / yr)",
       y = "EMI(per yr)") +
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
       filename = "Output Figures\\National Slope Regression - Veg zone - SLR & Diversity.jpg",
       units = "in",
       height = 8, width = 12, dpi = 300, limitsize = FALSE)



#Chapter 4: Landscape Conditions Regressions for the slope dataset --------------------------


#Task 1 - Run the national regressions for Landscape Condition to vegetation metrics

landscape_regression <- veg_format %>%
  group_by(Metric) %>%
  nest() %>%
  mutate(regression = map(.x = data,
                          ~lm(Value ~ Landscape * Vegetation_Zone,
                              data = .x) %>%
                            anova() %>%
                            tidy()), 
         rsquared = map(.x = data,
                        ~lm(Value ~ Landscape * Vegetation_Zone,
                            data = .x) %>%
                          glance() %>%
                          select(r.squared, adj.r.squared))) %>%
  unnest(c(regression, rsquared)) %>%
  select(-data) %>%
  mutate(across(sumsq:adj.r.squared, ~round(., 3)),
         p.value = round(p.value, 4))

glimpse(landscape_regression)


write.csv(landscape_regression,
          "Output Stats\\National Regression of Slopes & Landscape - Veg Zone - ANOVA Table.csv")


#Task 2 - Predict the linear regressions of the statistically significant interactions

# Run a dplyr loop to predict the values of the regressions

regression_predicted <- veg_format %>%
  mutate(Vegetation_Zone = factor(Vegetation_Zone)) %>%
  group_by(Metric) %>%
  nest() %>%
  mutate(regression = map(.x = data,
                          ~lm(Value ~ Landscape * Vegetation_Zone,
                              data = .x) %>%
                            ggpredict(., 
                                      terms = c("Landscape [all]", "Vegetation_Zone [all]"),
                                      type = "fixed", interval = "confidence"))) %>%
  select(-data) %>%
  unnest(regression) %>%
  rename(Landscape = x,
         Value_pred = predicted,
         Vegetation_Zone = group) %>%
  filter(Metric == "EMI" |
           Metric == "Freshwater Cover" |
           Metric == "Shannon-Weiner Diversity" |
           Metric == "Richness")

glimpse(regression_predicted)


#Calculate the slopes of each vegetation metric across Vegetation Zone

regression_slopes <- regression_predicted %>%
  group_by(Metric, Vegetation_Zone) %>%
  summarise(slope = (Value_pred[which(Landscape == max(Landscape))] - Value_pred[which(Landscape == min(Landscape))]) / ((max(Landscape) - min(Landscape)))) %>%
  ungroup() %>%
  mutate(slope = round(slope, 3))

write.csv(regression_slopes,
          "Output Stats\\National Slope Models - Landscape Slopes Across All Zones - Veg Zones.csv")


#Task 2 - Graph the statistically appropriate regressions

#Instead of predicting the slopes and the values with ggpredict, will use the geom_stat function 
# in ggplot2, since the regressions are simple linear regressions

veg_landscape <- veg_format %>%
  filter(Metric == "Shannon-Weiner Diversity")

regression_landscape <- regression_predicted %>%
  filter(Metric == "Shannon-Weiner Diversity")

landscape_graph <- ggplot() +
  geom_point(data = veg_landscape,
             aes(x = Landscape,
                 y = Value,
                 fill = Vegetation_Zone),
             alpha = 0.75, size = 7, shape = 21,
             position = "jitter") +
  geom_ribbon(data = regression_landscape,
              aes(x = Landscape,
                  ymax = Value_pred + std.error,
                  ymin = Value_pred - std.error,
                  fill = Vegetation_Zone),
              alpha = 0.35) + 
  geom_line(data = regression_landscape,
            aes(x = Landscape, y = Value_pred,
                colour = Vegetation_Zone),
            linetype = "dashed", size = 1.25) +
  scale_x_continuous(limits = c(1.75, 10.25),
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
       filename = "Output Figures\\National Slope Regression - Veg zone - Landscape & Diversity.jpg",
       units = "in",
       height = 8, width = 12, dpi = 300, limitsize = FALSE)

