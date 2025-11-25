#Project: NAMASTE National Univariate Exploratory Analysis
#Script: Script 5 - National Mixed Linear Regressions by Region

#Author: Grant McKown, Jackson Estuarine Laboratory, University of New Hampshire
#Contact: james.mckown@unh.edu

#Project Description: The NERR system collects vegetation data along permanent transects and plots
#   in salt marshes throughout the United States. Vegetation data was synthesized, formatted, and 
#   filtered (based on field monitoring requirements). Reserve, site, transect, and plot-level characteristics
#   were also collected including hydrology, landscape geomorphology, elevation, and climate. The 
#   series of scripts conducts an exploratory analysis of different vegetation response factors to 
#   a series of explanatory and grouping factors. 

#Script Description: Script 5 formats the site-level vegetation data, conducts linear mixed models at the 
#                   the vegetation zone level across the nation. The slopes of significant trends are
#                   calculated. Model summaries and regression slopes are exported as CSV. Graphs are
#                   created at the end of the script. 

#Note - Script 5 is very similar to Script 3, since it accomplishes very similar analyses and research goals. 


#Chapter 1: Import package library -----------------------------------------------


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



#Chapter 2: Format the vegetation dataset ---------------------------------------


veg <- read.csv("Formatted Datasets\\Veg Dataframe Summarised By Site.csv") %>%
  select(-X)

glimpse(veg)




#Chapter 3: Calculate live cover slope of each individual site, vegetation zone -------------------------

# Formatting and organizing the dataset for regression slopes including:
# Gathering the Site Variables and gathering the vegetation metrics, thus transforming the dataset
# into a very LONG dataset. Essentially, breaking down each row into ~30 rows. This will allow
# for simultaneous linear regressions for each vegetation metric for each site variable!

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





#Chapter 4: #Conduct mixed linear regressions of each vegetation metric by vegetation zone -------------------------


# This is accomplished by nesting the data essentially by Vegetation Zone and Vegetation Metric (long dataset)

# Mixed linear modeling is a simple time model with random factors of Reserve and nested Site 

regression_models_site <- veg_format %>%
  gather(Region:salinity, 
         key = "Site_variable", 
         value = "Category") %>%
  group_by(Metric, Site_variable) %>%
  mutate(SampleSize = n()) %>%
  ungroup() %>%
  group_by(Metric, Site_variable, SampleSize) %>%
  nest() %>%
  mutate(regression = map(.x = data,
                          ~lmer(Value ~ Year * Category + (1|Reserve/SiteID),
                                      data = .x) %>%
                            anova() %>%
                            tidy())) %>%
        
  unnest(regression) %>%
  select(-data) %>%
    mutate(across(sumsq:p.value, ~round(., 3))) %>%
  ungroup() 

glimpse(regression_models_site)


write.csv(regression_models_site,
          "Output Stats\\Mixed Linear Time Model - Regression Model Outputs by Site Characteristic.csv")



#Chapter 5: Calculate the Slopes of each Vegetation Metric Across All Site Characteristics -------------------------

# Run a dplyr loop to predict the values of the national vegetation regressions

regression_predicted <- veg_format %>%
  gather(Region:salinity, 
         key = "Site_variable", 
         value = "Category") %>%
  group_by(Metric, Site_variable) %>%
  nest() %>%
  mutate(regression = map(.x = data,
                          ~lmer(Value ~ Year * Category + (1|Reserve/SiteID),
                                data = .x) %>%
                            ggpredict(., 
                                      terms = c("Year [all]"),
                                      type = "fixed", interval = "confidence"))) %>%
  select(-data) %>%
  unnest(regression) %>%
  rename(Year = x,
         Value_pred = predicted) %>%
  select(Metric, Site_variable, Year, Value_pred, std.error, conf.low, conf.high) %>%
  mutate(across(Value_pred:conf.high, ~round(., 3)))

glimpse(regression_predicted)


#Calculate the slopes of each vegetation metric across site characteristics

regression_slopes <- regression_predicted %>%
  group_by(Site_variable, Metric) %>%
  summarise(slope = (Value_pred[which(Year == max(Year))] - Value_pred[which(Year == min(Year))]) / ((max(Year) - min(Year)))) %>%
  ungroup() %>%
  mutate(slope = round(slope, 3))

write.csv(regression_slopes,
          "Output Stats\\Mixed Linear Time Model - Slopes Across All Site Characteristics.csv")




#Chapter 6: Calculate the Slopes of each Vegetation Metric For Each Site Characteristic ------------------------------


#Re-run the combination of map, ggpredict chunk of code, except with the addition of the 
# Category [all] to tell the software to also include individual site categories for each
# site characteristic (Region, salinity, etc.)

regression_predicted_site <- veg_format %>%
  gather(Region:salinity, 
         key = "Site_variable", 
         value = "Category") %>%
  group_by(Metric, Site_variable) %>%
  nest() %>%
  mutate(regression = map(.x = data,
                          ~lmer(Value ~ Year * Category + (1|Reserve/SiteID),
                                data = .x) %>%
        ggpredict(., 
                  terms = c("Year [all]", "Category [all]"),
                  type = "fixed", interval = "confidence"))) %>%
  select(-data) %>%
  unnest(regression) %>%
  rename(Year = x,
         Value_pred = predicted,
         Category = group) %>%
  mutate(stderr.high = Value_pred + std.error,
         stderr.low = Value_pred - std.error) %>%
  select(Metric, Site_variable, Category, Year, Value_pred, std.error, conf.low, conf.high, stderr.low, stderr.high) %>%
  mutate(across(Value_pred:stderr.high, ~round(., 3)))

glimpse(regression_predicted_site)


#Calculate the slopes of each vegetation metric across Site Characteristics

regression_slopes_site <- regression_predicted_site %>%
  group_by(Site_variable, Category, Metric) %>%
  summarise(slope = (Value_pred[which(Year == max(Year))] - Value_pred[which(Year == min(Year))]) / ((max(Year) - min(Year)))) %>%
ungroup() %>%
  mutate(slope = round(slope, 3))

glimpse(regression_slopes_site)

write.csv(regression_slopes_site,
          "Output Stats\\Mixed Linear Time Model - Slopes by Site Characteristic.csv")




#Chapter 7: Graph the General Trends Across Time and Site Characteristic ----------------------

#Data visualization is broken down into 4 graphs:
  # (1) Vegetation Cover - Live, Abiotic, Freshwater, and Halophyte
  # (2) Species - Richness, Diversity
  # (3) Ratios - Salt Ratio, EMI
  # (4) Combine all three graphs into one comprehensive, national time model graph

#Task 1: Graph the results of the predicted values for the Abiotic, Live Cover, Halophyte Cover, 
          #and Freshwater Cover

regression_predicted_cover <- regression_predicted_site %>%
  filter(Site_variable == "tidal_range",
         Metric == "Halophyte Cover", ) %>%
  mutate(stderr.high = ifelse(stderr.high > 100, 100, stderr.high),
         stderr.low = ifelse(stderr.low < 0, 0, stderr.low),
         Value_pred = ifelse(Value_pred < 0, 0, Value_pred),
         Category = factor(Category,
                           levels = c("Microtidal", "Mesotidal", "Macrotidal")))


national_cover_graph <- ggplot(data = regression_predicted_cover,
                                 aes(x = Year,
                                     y = Value_pred)) + 
  geom_point(data = veg,
             aes(x = Year,
                 y = halophyte,
                 colour = Region),
             size = 4.5, alpha = 0.35, colour = "darkblue") +
  geom_ribbon(aes(x = Year, ymin = stderr.low, 
                  ymax = stderr.high),
              alpha = 0.85, fill = "grey") + 
  geom_line(colour = "black", 
            linewidth = 1.25, linetype = "dashed") + 
  scale_y_continuous(limits = c(0, 102),
                     breaks = seq(0, 100, 20),
                     expand = c(0,0)) +
  scale_x_continuous(limits = c(2005, 2022.5),
                     breaks = seq(2006, 2022, 2),
                     expand = c(0,0)) +
  labs(y = "Halophyte Cover (%)",
       x = "") +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 18, colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    strip.text = element_text(size = 18, colour = "black"),
    strip.background = element_blank()) + 
  facet_wrap(~Category, 
             ncol = 3)


national_cover_graph


ggsave(national_cover_graph,
       filename = "Output Figures\\National Time Mixed Model - Tidal Range - Halophyte Cover.jpg",
       units = "in",
       height = 6, width = 16, dpi = 300, limitsize = FALSE)



#Task 2 - Graph the results of the predicted values for the Shannon-Diversity & Richness


# Note - due to no significant trend found for either metric, the regression_predicted, not broken down by 
# different regions, will be shown. 


regression_predicted_richness <- regression_predicted_site %>%
  filter(Metric == "Shannon-Weiner Diversity", 
         Site_variable == "Region") %>%
    mutate(stderr.low = ifelse(stderr.low < 0, 0, stderr.low), 
           Category = factor(Category, levels = 
                               c("Northeast", "Mid-Atlantic", "Southeast", "Gulf Coast",
                                          "West Coast")))

national_richness_graph <- ggplot(data = regression_predicted_richness,
                               aes(x = Year,
                                   y = Value_pred)) + 
  geom_point(data = veg,
             aes(x = Year,
                 y = EMI,
                 colour = Region),
             size = 4.5, alpha = 0.35, colour = "darkblue") +
  geom_ribbon(aes(x = Year, ymin = stderr.low, 
                  ymax = stderr.high),
              alpha = 0.85, fill = "grey") + 
  geom_line(colour = "black", 
            linewidth = 1.25, linetype = "dashed") + 
  scale_y_continuous(limits = c(0, 1.05),
                     breaks = seq(0, 1, 0.20),
                     expand = c(0,0)) +
  scale_x_continuous(limits = c(2005, 2022.5),
                     breaks = seq(2006, 2022, 2),
                     expand = c(0,0)) +
  labs(y = "Shannon-Weiner Diversity",
       x = "") +
  theme_bw() +
  theme(
    legend.position = c(0.15, 0.875),
    legend.title = element_blank(),
    legend.text = element_text(size = 18, colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18)) + 
  facet_wrap(~Category,
             nrow = 2, ncol = 3)

national_richness_graph


ggsave(national_richness_graph,
       filename = "Output Figures\\National Time Mixed Model - Region - Diversity.jpg",
       units = "in",
       height = 8, width = 16, dpi = 300, limitsize = FALSE)



#Task 3 - Graph the results of the predicted values for the Salt Ratio and EMI

#A significant trend was found for EMI in region, so the regressions will be shown in panels brokn down
# by regions

regression_predicted_ratio <- regression_predicted_site %>%
  filter(Metric == "EMI") %>%
  filter(Site_variable == "Region") %>%
  mutate(stderr.high = ifelse(stderr.high > 1, 1, stderr.high),
         stderr.low = ifelse(stderr.low < 0, 0, stderr.low),
         Category = factor(Category, levels = 
                             c("Northeast", "Mid-Atlantic", "Southeast",
                               "Gulf Coast", "West Coast")))
  

national_ratio_graph <- ggplot(data = regression_predicted_ratio,
                                  aes(x = Year,
                                      y = Value_pred)) + 
  geom_point(data = veg,
             aes(x = Year,
                 y = EMI,
                 colour = Region),
             size = 4.5, alpha = 0.35, colour = "darkblue") +
  geom_ribbon(aes(x = Year, ymin = stderr.low, 
                  ymax = stderr.high),
              alpha = 0.85, fill = "grey") + 
  geom_line(colour = "black", 
            linewidth = 1.25, linetype = "dashed") + 
  scale_y_continuous(limits = c(0, 1.05),
                     breaks = seq(0, 1, 0.20),
                     expand = c(0,0)) +
  scale_x_continuous(limits = c(2005, 2023),
                     breaks = seq(2006, 2022, 2),
                     expand = c(0,0)) +
  labs(y = "Ecotone Migration Index",
       x = "") +
  theme_bw() +
  theme(
    legend.position = c(0.80, 0.125),
    legend.title = element_blank(),
    legend.text = element_text(size = 18, colour = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 18, colour = "black"),
    axis.text = element_text(size = 14, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 18)) +
  facet_wrap(~Category,
            nrow = 2,
            ncol = 3)


national_ratio_graph


ggsave(national_ratio_graph,
       filename = "Output Figures\\National Time Mixed Model - Region - EMI.jpg",
       units = "in",
       height = 8, width = 16, dpi = 300, limitsize = FALSE)



