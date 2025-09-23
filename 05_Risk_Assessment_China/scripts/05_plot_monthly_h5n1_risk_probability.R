# ==============================================================================
# Description: This script visualizes the monthly predicted H5N1 risk probability 
#              distribution across China for a specific scenario and year. 
#              It reads monthly MaxEnt prediction results (January to December), 
#              processes the data, and creates a faceted map showing the spatial 
#              distribution of predicted risk probabilities for each month.
#              This example uses SSP126 scenario for the year 2025.
# ==============================================================================

# Load required libraries
library(terra)
library(sf)
library(tidyverse)

# Set working directory
setwd("model_projection_2025/permonth_projection_output")

# Create vector containing English month names (January to December)
months_vector <- c("January", "February", "March", "April", "May", "June", 
                   "July", "August", "September", "October", "November", "December")

# Read China boundary shapefile
china_map <- read_sf('map_of_China.json')

# Initialize data frame to store prediction results
prediction_df <- data.frame()

# Process each monthly prediction file (using SSP126 scenario 2025 as example)
for (i in months_vector) {
  # Read monthly prediction raster file
  month_i_rast <- rast(paste('Predicted_distribution_of_H5N1_in_China_', i, '.asc', sep = ''))
  
  # Convert raster to data frame with cell coordinates
  month_i_df <- as.data.frame(month_i_rast, cell = T, xy = T)
  
  # Rename probability column and add month information
  names(month_i_df)[4] <- 'Predicted_probability'
  month_i_df$month <- i
  prediction_df <- bind_rows(prediction_df, month_i_df)
}

# Convert month to factor with correct order
prediction_df$month <- factor(prediction_df$month, levels = months_vector)

# Categorize predicted probability into bins (0-1 with 0.1 intervals)
prediction_df <- prediction_df |> 
  mutate(Predicted_category = factor(cut(Predicted_probability,
                                         breaks = seq(0, 1, by = 0.1),
                                         include.lowest = TRUE),
                                     levels = rev(levels(cut(Predicted_probability,
                                                             breaks = seq(0, 1, by = 0.1),
                                                             include.lowest = TRUE)))))

# Create faceted plot showing monthly risk probability distribution
prediction_results <- 
  ggplot() +
  # Plot probability as colored tiles
  geom_tile(data = prediction_df, aes(x = x, y = y, fill = Predicted_probability)) +
  # Create facets for each month (2 rows, 6 columns)
  facet_wrap(~ month, nrow = 2, ncol = 6) + 
  # Set color gradient from white to orange
  scale_fill_gradient(
    low = "white",           # Start with white
    high = "#fe9929",        # Medium orange (contrasts with blue)
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    name = "Probability"     # Legend title
  ) +
  # Add China boundary overlay
  geom_sf(data = china_map, fill = NA, color = "darkgray", size = 0.1) +
  theme_classic() +
  theme(strip.background = element_blank(),  
        strip.placement = "outside",  
        axis.line = element_blank(),         
        axis.ticks = element_blank(),      
        axis.text.x = element_blank(),   
        axis.text.y = element_blank(),       
        axis.title.x = element_blank(),       
        axis.title.y = element_blank(),
        strip.text = element_text(colour = 'black', size = 13),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.key.spacing.y = unit(2, 'mm'))   

# Save the plot as PDF
# Save the plot as PDF
ggsave(plot = prediction_results,
       filename = 'H5N1_risk_probability_SSP126_2025_January_to_December.pdf',
       width = 10, height = 8)


