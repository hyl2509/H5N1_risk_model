# Load required libraries
library(tidyverse)
library(lubridate)
library(scales)

# ---------------------------
# DATA IMPORT AND PREPARATION
# ---------------------------

# Read Markov jump statistics data (2021-2023)
jump_counts <- read_tsv('europe_markov_jumps_2021_2023.tsv')

# Set month factor levels (July-December first, then January-June)
jump_counts$month <- factor(jump_counts$month, levels = c(7:12, 1:6))

# ---------------------------
# DATA VISUALIZATION
# ---------------------------

# Create time series plot of jumps originating from North-Central Europe in 2022
jump_counts |> 
  # Filter for 2022 data from North-Central region
  filter(year == 2022) |> 
  filter(from == 'North-Central') |> 
  # Create date objects and transmission path labels
  mutate(
    date = as.Date(paste(year, month, 1, sep = '-')),
    path = paste(from, 'to', to, sep = ' ')
  ) |> 
  # Generate line plot
  ggplot() + 
  geom_line(
    aes(x = month, y = n_jump, color = path, group = path), 
    linewidth = 0.7
  ) +
  # Set custom color scheme for transmission paths
  scale_color_manual(
    values = c(
      "North-Central to North-Eastern" = "#4B65AF",  # Blue
      "North-Central to South-Western" = "#EF673F"   # Orange
    )
  ) +
  # Configure axis and legend labels
  labs(
    x = 'Month', 
    y = 'Jump counts', 
    color = 'Transmission path'
  ) +
  # Apply black-and-white theme with customizations
  theme_bw() +
  theme(
    axis.text = element_text(color = 'black', size = 14),
    axis.title = element_text(color = 'black', size = 13),
    legend.title = element_text(color = 'black', size = 13),
    legend.text = element_text(color = 'black', size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = 'darkgray'),
    legend.position = 'none'  # Remove legend
  )

# ---------------------------
# OUTPUT SAVING
# ---------------------------

# Save the visualization as PDF
ggsave(
  path = '.',
  filename = 'europe_north_central_jump.pdf',
   dpi = 300,
   width = 5, height = 4)
  
  