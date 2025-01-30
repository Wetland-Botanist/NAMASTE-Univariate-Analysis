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

#Script Description: Script 9 conducts a Spearman Rank correlation analysis between vegetation metrics at the plot 
# level to better understand how vegetation metrics respond. Several types of correlation graphs are created. 


#-----------------------------------
#Chapter 1: Import package library
#-----------------------------------

#Data Analysis Packages

library(tidyverse)
library(Hmisc)

#Graphing and Data Visualization Packages

library(corrplot)
library(GGally)

#--------------------------------------------
#Chapter 2: Format the vegetation dataset
#------------------------------------------

veg <- read.csv("Formatted Datasets\\Veg Dataframe Summarised By Site.csv") %>%
  select(-X)
  
glimpse(veg)

veg_format <- veg %>%
  select(abiotic_cover, live_cover, halophyte, freshwater, EIR, salt_ratio, Richness, SWdiv) %>%
  rename(
    Abiotic_Cover = abiotic_cover,
    Live_Cover = live_cover,
    Halophyte_Cover = halophyte,
    Freshwater_Cover = freshwater,
    EMI = EIR,
    Salt_Ratio = salt_ratio,
    Species_Richness = Richness,
    Shannon_Weiner = SWdiv)


#-------------------------------------------
#Chapter 3: Calculate Correlation Matrix
#-------------------------------------------

spearman_corr <- rcorr(x = as.matrix(veg_format),
                     type = "spearman")

glimpse(spearman_corr)

spearman_r <- tidy(spearman_corr$r)

spearman_p <- tidy(spearman_corr$P)

write.csv(spearman_r, 
          "Output Stats\\Vegetation Spearman Correlation - r values.csv")

write.csv(spearman_p,
          "Output Stats\\Vegetation Spearman Correlation - p values.csv")


#--------------------------------------------------------------------------
#Chapter 4: Create Correlation Matrix, Create Plot - Cressman Method
#------------------------------------------------------------------------

# function to format the scatterplots
# could change method to lm, gam, whatever other things you
# can use in geom_smooth. can also change colors, transparency, etc.
lowerFn <- function(data, mapping, method = "lm", ...) {
  p <- ggplot(data = data, mapping = mapping) +
    geom_point(colour = "navy", alpha = 0.8) +
    geom_smooth(method = method, color = "darkorange", 
                se = FALSE, linewidth = 0.9, ...)
  p
}

# build the plot and correlation matrix
# with 'method' can use pearson, spearman, kendall
# make sure your data frame is either all numeric or
# that you subset the numeric columns
veg_spearman_corr <- veg_format |> 
  ggpairs(upper = list(continuous = wrap("cor", 
                                         method = "spearman", 
                                         use = "pairwise.complete.obs",
                                         colour = "black")),
          lower = list(continuous = wrap(lowerFn)))

veg_spearman_graph <- veg_spearman_corr + 
  theme(
    axis.text = element_text(size = 10, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(size = 10, colour = "black")
  )

veg_spearman_graph

ggsave(veg_spearman_graph,
       height = 10, width = 14, dpi = 300,
       filename = "Output Figures\\Veg Spearman Correlation - Mini Graph Compilation.jpg")

#-------------------------------------------------------------
# Chapter 5: Create Correlation Figure - Oval Type
#------------------------------------------------------------

#subset out the columns you want to correlate and run that matrix through “cor” function

correlation_matrix <- cor(veg_format,  use = "pairwise.complete.obs",
                          method = "spearman") 


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

numeric_columns <- sapply(veg_format, is.numeric)

# Subset data to include only numeric columns

numeric_data <- veg_format[, numeric_columns]

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


