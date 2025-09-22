# ==============================================================================
# Description: Processes and visualizes monthly migratory bird abundance data
#              for Europe and North America regions
#
# Input: Monthly TIFF files in /globally_monthly_all_species_sum_of_abundance/
# Output: PDF maps of monthly abundance distributions in /result/
#
# Dependencies: See library imports below
# ==============================================================================

# --------------------------
# 1. Load Required Packages
# --------------------------
library(tidyverse)    
library(sf)          
library(terra)       
library(rnaturalearth) 
library(rnaturalearthhires) 
library(raster)      
library(lubridate)   
library(scales)       
library(RColorBrewer) 
library(patchwork)    

# --------------------------
# 2. Data Loading
# --------------------------

# Load global monthly bird abundance data from TIFF files
bird_abundance_files <- list.files(
  path = 'globally_monthly_all_species_sum_of_abundance/',
  full.names = TRUE, 
  recursive = TRUE
)

# Read all raster files and store in list
bird_abundance <- lapply(bird_abundance_files, rast)

# Name list elements by month (removing filename prefixes/suffixes)
names(bird_abundance) <- basename(bird_abundance_files) |> 
  str_remove_all('Sum_of_bird_abundance_in_|.tif')

# --------------------------
# 3. Europe Analysis
# --------------------------

# Download European country boundaries
ne_countries <- ne_download(
  scale = 110, 
  category = "cultural",
  type = "countries",
  returnclass = "sf"
)

# Filter for European continent
europea <- ne_countries |> 
  filter(CONTINENT == 'Europe')

# Initialize empty dataframe for European results
europea_abundance_df <- data.frame()

# Process each month's data
for (i in month.name) {
  print(paste('Processing ', i, sep = ''))
  
  # Project raster to match European boundary CRS
  month_i_abundance <- bird_abundance[[i]]
  month_i_abundance_project <- project(month_i_abundance, crs(europea))
  
  # Mask and crop to European boundaries
  europea_month_i_abundance <- mask(month_i_abundance_project, europea)
  europea_month_i_abundance <- crop(europea_month_i_abundance, europea)
  
  # Convert to dataframe with coordinates
  raster_df <- as.data.frame(europea_month_i_abundance, xy = TRUE)
  raster_df$month <- i
  europea_abundance_df <- bind_rows(europea_abundance_df, raster_df)
}

# Apply log2 transformation for visualization
europea_abundance_df <- europea_abundance_df %>% 
  mutate(log_sum = log2(sum + 1))

# Ensure months are ordered correctly
europea_abundance_df$month <- factor(europea_abundance_df$month, 
                                     levels = month.name)

# Define visualization parameters
original_breaks <- c(0, 3, 15, 63, 255, 1023)
log_breaks <- log2(original_breaks + 1)

# Create European abundance visualization
ggplot() +
  geom_tile_rast(data = europea_abundance_df, 
                 aes(x = x, y = y, fill = log_sum)) +
  facet_wrap(~ month, nrow = 2, ncol = 6) +
  scale_fill_gradientn(
    colours = rev(viridis::plasma(n = 100))[c(10, 19, 28, 46, 55, 64, 73, 82, 86)],
    values = scales::rescale(c(0, 1, 2, 3, 4, 5, 6, 7, 8)),  
    breaks = log_breaks,
    labels = original_breaks
  ) + 
  geom_sf(data = europea, fill = NA, lwd = 0.1, color = 'darkgray') +
  labs(fill = "Abundance", x = NULL, y = NULL) +
  theme_light() +
  coord_sf(xlim = c(-25, 60), ylim = c(35, 75)) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    strip.text = element_text(colour = 'black', size = 9),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10)
  )

# Save European visualization
ggsave(
  path = 'result',
  filename = 'Monthly_abundance_of_migratory_birds_in_Europe.pdf',
  width = 8, height = 4,
  dpi = 300
)

# --------------------------
# 4. North America Analysis
# --------------------------

# Download North American country boundaries
na_countries <- rnaturalearth::ne_countries(
  scale = 110, 
  type = 'countries',
  continent = 'North America'
)

# Exclude USA and Canada (will process states/provinces separately)
na_countries_exclude_usa_can <- na_countries |> 
  filter(!sovereignt %in% c('United States of America', 'Canada')) |> 
  select(country = sovereignt)

# Download US and Canada state-level boundaries
usa_can_state <- ne_states(country = c('United States of America', 'Canada')) |>
  select(country = admin, state = name)

# Combine country and state boundaries
na_sf <- bind_rows(na_countries_exclude_usa_can, usa_can_state)

# Classify regions (North, Central, South)
na_sf_region <- na_sf |> 
  mutate(region = case_when(
    country == 'United States of America' & state == 'Alaska' ~ 'North', 
    country == 'United States of America' & 
      state %in% c('California', 'Arizona', 'New Mexico', 'Texas', 'Louisiana', 
                   'Oklahoma', 'Arkansas', 'Mississippi', 'Alabama', 'Tennessee', 
                   'North Carolina', 'South Carolina', 'Georgia', 'Florida', 
                   'Hawaii') ~ 'South',
    country == 'United States of America' ~ 'Central', 
    country %in% c('Canada', 'Denmark') ~ 'North',
    .default = 'South'
  ))

# Define crop box for North America
crop_box <- st_as_sfc(st_bbox(c(xmin = -0.1, xmax = -180, ymin = -90, ymax = 90), 
                              crs = st_crs(na_sf_region)))

# Crop data to bounding box
na_croped_df <- st_intersection(na_sf_region, crop_box)
na_croped_df$region <- factor(na_croped_df$region, 
                              levels = c('North', 'Central', 'South'))

# Union geometries by region
na_croped_df_union <- na_croped_df |> 
  group_by(region) |> 
  summarize(
    geometry = st_union(geometry, is_coverage = TRUE)
  ) |> 
  st_make_valid()

# Initialize empty dataframe for North American results
north_america_abundance_df <- data.frame()

# Process each month's data
for (i in month.name) {
  print(paste('Processing ', i, sep = ''))
  
  # Project raster to match North American boundary CRS
  month_i_abundance <- bird_abundance[[i]]
  month_i_abundance_project <- project(month_i_abundance, crs(na_croped_df))
  
  # Mask and crop to North American boundaries
  north_america_month_i_abundance <- mask(month_i_abundance_project, na_croped_df)
  north_america_month_i_abundance <- crop(north_america_month_i_abundance, na_croped_df)
  
  # Convert to dataframe with coordinates
  raster_df <- as.data.frame(north_america_month_i_abundance, xy = TRUE)
  raster_df$month <- i
  north_america_abundance_df <- bind_rows(north_america_abundance_df, raster_df)
}

# Apply log2 transformation for visualization
north_america_abundance_df <- north_america_abundance_df |> 
  mutate(log_sum = log2(sum + 1))

# Ensure months are ordered correctly
north_america_abundance_df$month <- factor(north_america_abundance_df$month, 
                                           levels = month.name)

# Define visualization parameters
original_breaks <- c(0, 3, 15, 63, 255, 1023, 4095)
log_breaks <- log2(original_breaks + 1)

# Create North American abundance visualization
ggplot() +
  geom_tile_rast(data = north_america_abundance_df, 
                 aes(x = x, y = y, fill = log_sum)) +
  facet_wrap(~ month, nrow = 2, ncol = 6) +
  scale_fill_gradientn(
    colours = rev(viridis::plasma(n = 100))[c(10, 19, 28, 46, 55, 64, 73, 82)],
    values = scales::rescale(seq(0, 10, length.out = 9)),  
    breaks = log_breaks,
    labels = original_breaks
  ) +
  geom_sf(data = na_countries, fill = NA, lwd = 0.1, color = 'darkgray') +
  labs(fill = "Abundance", x = NULL, y = NULL) +
  theme_light() +
  coord_sf(expand = FALSE) +
  theme(
    strip.background = element_blank(),  
    strip.placement = "outside",  
    axis.line = element_blank(),         
    axis.ticks = element_blank(),      
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),       
    axis.title.y = element_blank(),
    strip.text = element_text(colour = 'black', size = 9),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10),
    panel.grid = element_blank()
  )

# Save North American visualization
ggsave(
  path = 'result',
  filename = 'Monthly_abundance_of_migratory_birds_in_North_America.pdf',
  dpi = 300,
  width = 8, 
  height = 4

)

