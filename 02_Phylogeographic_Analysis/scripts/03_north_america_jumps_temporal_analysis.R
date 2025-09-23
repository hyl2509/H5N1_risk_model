# Load required libraries
library(tidyverse)
library(lubridate)
library(scales)

# ---------------------------
# MARKOV JUMP TEMPORAL ANALYSIS
# ---------------------------

# Read Markov jump statistics (2022-2023)
jump_counts <- read_tsv('north_america_markov_jumps_2022_2023.tsv')

# Set month factor levels (July-December first, then January-June)
jump_counts$month <- factor(jump_counts$month, levels = c(7:12, 1:6))

# ---------------------------
# CENTRAL REGION JUMP DYNAMICS (2022)
# ---------------------------

# Plot transmission paths originating from Central region in 2022
jump_counts |> 
  filter(year == 2022) |> 
  filter(from == 'Central') |> 
  mutate(
    date = as.Date(paste(year, month, 1, sep = '-')),
    path = paste(from, 'to', to, sep = ' ')
  ) |> 
  ggplot() + 
  geom_line(
    aes(x = month, y = n_jump, color = path, group = path), 
    linewidth = 0.7
  ) +
  # Color scheme for transmission paths
  scale_color_manual(
    values = c(
      "Central to North" = "#4B65AF",  # Blue
      "Central to South" = "#EF673F"   # Orange
    )
  ) +
  labs(
    x = 'Month', 
    y = 'Jump counts', 
    color = 'Transmission path'
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(color = 'black', size = 14),
    axis.title = element_text(color = 'black', size = 13),
    legend.title = element_text(color = 'black', size = 13),
    legend.text = element_text(color = 'black', size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = 'darkgray'),
    legend.position = 'none'
  )

# Save the plot
ggsave(
  filename = 'north_america_central_jump.pdf',
  width = 5, 
  height = 4,
  dpi = 300
)
