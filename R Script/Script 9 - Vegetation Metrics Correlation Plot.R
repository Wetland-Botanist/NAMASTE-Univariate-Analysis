#Project: NAMASTE National Univariate Exploratory Analysis
#Script: Script 9 - Correlation Analysis between Vegetation Metrics

#Author: Grant McKown, Jackson Estuarine Laboratory, University of New Hampshire
#Contact: james.mckown@unh.edu

#Project Description: The NERR system collects vegetation data along permanent transects and plots
#   in salt marshes throughout the United States. Vegetation data was synthesized, formatted, and 
#   filtered (based on field monitoring requirements). Reserve, site, transect, and plot-level characteristics
#   were also collected including hydrology, landscape geomorphology, elevation, and climate. The 
#   series of scripts conducts an exploratory analysis of different vegetation response factors to 
#   a series of explanatory and grouping factors. 

#Script Description: Script 9 conducts a correlation analysis between vegetation metrics at the plot 
# level to better understand how vegeation metrics respond. 


#-----------------------------------
#Chapter 1: Import package library
#-----------------------------------

#Data Analysis Packages

library(dplyr)
library(tidyr)

#Graphing and Data Visualization Packages

library(corrplot)

#---------------------------------------------------------------------------------------------
#Chapter 2: Format the vegetation dataset
#---------------------------------------------------------------------------------------------

#Task 1: Import the national plot dataframe
veg <- read.csv("Formatted Datasets\\Veg Plot Dataframe Formatted.csv") %>%
  select(-X, -invasive_cover) %>%
  select(Reserve:Year, Region:salinity, abiotic_cover:salt_ratio) %>%
  rename(EMI = EIR) %>%
  mutate(across(abiotic_cover:salt_ratio,
                ~ifelse(is.na(.), 0, .))) %>%
  rename(Live_Cover = live_cover,
         Abiotic_Cover = abiotic_cover,
         Halophyte_Cover = halophyte,
         Freshwater_Cover = freshwater,
         Ecotone_Migraiton_Index = EMI,
         Salt_Ratio = salt_ratio,
         Shannon_Weiner = SWdiv,
         Species_Richness = Richness)

glimpse(veg)


#-------------------------------------------------------------
# Chapter 3: Create Correlation Matrix, Create Plot
#------------------------------------------------------------

#subset out the columns you want to correlate and run that matrix through “cor” function

selected_h3_data <- veg %>%
  select(Abiotic_Cover:Salt_Ratio)

correlation_matrix <- cor(selected_h3_data,  use = "pairwise.complete.obs") 


# Create a correlation heatmap

corrplot(correlation_matrix,
         
         method='ellipse',
         
         type='upper',  #prevents repeated variables in the plot (just upper half of the correlation)
         
         order='hclust', #help function for corrplot gives all the order choices
         
         tl.col="black",
         
         tl.srt=45,
         
         diag=FALSE,
         
         sig.level=0.05,
         
         insig='label_sig')

#Then re-run to get the asterisks for significance placed ontop

numeric_columns <- sapply(selected_h3_data, is.numeric)

# Subset data to include only numeric columns

numeric_data <- selected_h3_data[, numeric_columns]

# Calculate the correlation matrix and obtain p-values

correlation_results <- Hmisc::rcorr((as.matrix(numeric_data)))

# Extract correlation matrix

correlation_matrix <- correlation_results$r

print(correlation_matrix)

# Extract p-values matrix

p_values_matrix <- correlation_results$P

# Display the p-values matrix

print(p_values_matrix)

significance_level <- 0.05

# Create a matrix indicating whether p-values are below the significance level

significant_matrix <- p_values_matrix < significance_level

# Display the matrix of significant correlations

print(significant_matrix)

png(filename = "mycorrplot.png", width = 1200, height = 800)

corrplot(
  
  correlation_matrix,
  
  method='ellipse',
  
  type='upper',
  
  order='hclust',
  
  tl.col="black",
  
  tl.srt=45,
  
  tl.cex = 1.2,  #
  
  diag=FALSE,
  
  sig.level=0.05,  # You can adjust the significance level as needed
  
  insig='label_sig',  # Display non-significant correlations as blank
  
  p.mat = p_values_matrix,
  
  cl.cex = 1.15# Provide the matrix of p-values
  
)


dev.off()
