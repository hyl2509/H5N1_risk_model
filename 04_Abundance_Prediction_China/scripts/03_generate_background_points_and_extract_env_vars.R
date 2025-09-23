# ==============================================================================
# Description: This script processes bird observation data for species distribution 
#              modeling, using Common Merganser as an example. It generates pseudo-absence 
#              (background) points that are spatiotemporally separated from actual 
#              detection points, then extracts multiple environmental variables 
#              for both detection and background points including land cover, 
#              elevation, NDVI, precipitation, temperature, and relative humidity.
#
# Key steps:
# 1. Load and preprocess bird observation data
# 2. Generate spatiotemporally constrained background points
# 3. Extract environmental variables for all points (detection + background)
# 4. Combine all data and export for species distribution modeling
#
# Inputs: Bird observation checklists, China boundary, environmental raster data
# Outputs: CSV file with observation data and extracted environmental variables
# ==============================================================================

# Load required libraries
library(tidyverse)
library(sf)
library(terra)
library(units)
library(readxl)
library(exactextractr)
library(landscapemetrics)
library(raster)
library(ncdf4)

# Set seed for reproducibility and define target bird species
set.seed(2024)
bird_name <- 'Common_Merganser'

# Set working directory to species-specific folder
setwd(paste("42_birds_observation_data/", bird_name, 
            '/data_after_raw_data_preprocess', sep = ''))

# Load China boundary map and bird observation data
china_map <- read_sf("map_of_China.json")
detection <- read_csv(paste('checklists_', bird_name, '.csv', sep = ''))

# Filter observations with valid effort hours
detection <- detection |> 
  filter(effort_hours >= 0)

# ==============================================================================
# SECTION 1: Process detection data and create spatial buffers
# ==============================================================================

# Convert detection data to spatial format
detection_sf <- detection |> 
  dplyr::select(observation_date, latitude, longitude) |> 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Visualize detection points on China map
ggplot(china_map) + 
  geom_sf() +
  geom_sf(data = detection_sf, color = 'red', size = 0.5)

# Ensure detection points are within China boundaries
sf_use_s2(F)
detection_sf_within_china <- st_intersection(detection_sf, china_map)
sf_use_s2(T)

# Create 6km buffers around detection points
detection_sf_buffers <- st_buffer(detection_sf_within_china, dist = set_units(6, "km"))

# ==============================================================================
# SECTION 2: Generate background (pseudo-absence) points
# ==============================================================================

# Generate random dates within the detection period
start_date <- as.Date(min(detection_sf_buffers$observation_date))
end_date <- as.Date(max(detection_sf_buffers$observation_date))

random_dates <- sample(seq.Date(start_date, end_date, by = "day"), 
                       size = nrow(detection_sf_within_china), replace = TRUE)

nondetection <- data.frame(observation_date = random_dates)

# Count occurrences of each date
nondetection_unique_date <- nondetection |> 
  count(observation_date)

# Generate background points outside spatiotemporal buffers of detections
nondetection_random_point <- data.frame()

for (i in 1:nrow(nondetection_unique_date)) {
  print(paste('process row', i, sep = ' '))
  nondetection_i <- nondetection_unique_date[i, ]
  
  # Find detections within ±15 days of current background date
  detection_sf_buffers_i <- detection_sf_buffers |> 
    filter(observation_date >= nondetection_i$observation_date - 15, 
           observation_date <= nondetection_i$observation_date + 15)
  
  # Create union of detection buffers and ensure valid geometry
  combined_buffers <- st_union(detection_sf_buffers_i)
  combined_buffers <- st_make_valid(combined_buffers)
  china_map <- st_make_valid(china_map)
  
  # Project to planar coordinate system for geometric operations
  china_map_projected <- st_transform(china_map, crs = 3857)
  combined_buffers_projected <- st_transform(combined_buffers, crs = 3857)
  china_map_projected <- st_make_valid(china_map_projected)
  combined_buffers_projected <- st_make_valid(combined_buffers_projected)
  
  # Calculate area outside detection buffers
  non_buffer_area <- st_difference(china_map_projected, combined_buffers_projected)
  
  # Generate random points in non-detection area
  random_points <- st_sample(non_buffer_area, nondetection_i$n)
  
  # Convert to proper spatial format and transform to geographic coordinates
  random_points_sf <- st_as_sf(random_points)
  random_points_sf <- st_transform(random_points_sf, crs = 4326)
  
  # Extract coordinates and add date information
  coordinates <- st_coordinates(random_points_sf)
  coordinates_df <- as.data.frame(coordinates)
  colnames(coordinates_df) <- c("longitude", "latitude")
  coordinates_df_add_date <- coordinates_df |> 
    mutate(observation_date = nondetection_i$observation_date)
  nondetection_random_point <- bind_rows(nondetection_random_point, coordinates_df_add_date)
}

# Add required variables for background points
nondetection_random_point_add_variable <- nondetection_random_point |> 
  mutate(checklist_id = str_c('F', row.names(nondetection_random_point), sep = ''),
         observation_count = 0,
         species_observed = FALSE)

# Add random sampling effort and temporal variables
nondetection_random_point_add_variable$hours_of_day <- runif(nrow(nondetection_random_point_add_variable), min = 0, max = 24)
nondetection_random_point_add_variable$effort_hours <- runif(nrow(nondetection_random_point_add_variable), min = 0, max = 10)
nondetection_random_point_add_variable$year <- year(nondetection_random_point_add_variable$observation_date)
nondetection_random_point_add_variable$day_of_year <- yday(nondetection_random_point_add_variable$observation_date)
nondetection_random_point_add_variable$month <- month(nondetection_random_point_add_variable$observation_date)
nondetection_random_point_add_variable$month <- sprintf("%02d", nondetection_random_point_add_variable$month)

nondetection <- nondetection_random_point_add_variable

# Combine detection and background data
observation_data <- bind_rows(detection, nondetection)

# Ensure month format is consistent
observation_data$month <- month(observation_data$observation_date)
observation_data$month <- sprintf("%02d", observation_data$month)

# ==============================================================================
# SECTION 3: Create monthly distribution maps
# ==============================================================================

# Create vector of English month names for file naming
months_vector <- c("January", "February", "March", "April", "May", "June", 
                   "July", "August", "September", "October", "November", "December")
names(months_vector) <- 1:12

# Convert observation data to spatial format for mapping
checklists_sf <- observation_data |> 
  dplyr::select(observation_date, species_observed, longitude, latitude) |> 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Generate monthly distribution maps
for (i in 1:12) {
  month_i <- checklists_sf |> 
    filter(month(observation_date) == i)
  
  distribution_plot <- ggplot(china_map) + 
    geom_sf() + 
    geom_sf(data = month_i, aes(color = species_observed), size = 1) +
    labs(title = months_vector[[i]]) +
    theme(plot.title = element_text(hjust = 0.5))
  
  plot_name <- paste("42_birds_observation_data/", bird_name,
                     '/permonth_observation_point_distribution_in_map/', months_vector[[i]], '.tif',
                     sep = '')
  ggsave(plot = distribution_plot, filename = plot_name, width = 5, height = 4)
}

# ==============================================================================
# SECTION 4: Extract environmental variables for all observation points
# ==============================================================================

# Create spatial points and 3km buffers for environmental variable extraction
observation_data_sf <- observation_data |> 
  dplyr::select(checklist_id, year, latitude, longitude) |> 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

buffers <- st_buffer(observation_data_sf, dist = set_units(3, "km"))

# ------------------------------------------------------------------------------
# 4.1 Extract land cover data and calculate landscape metrics
# ------------------------------------------------------------------------------

# Prepare buffers for land cover extraction with year constraints
buffers_for_landcover <- buffers |> 
  mutate(year_lc = case_when(year < 2000 ~ '2000',
                             year > 2023 ~ '2023',
                             .default = as.character(year)))

# Load land cover class definitions and data
lc_classes <- read_xlsx("F:/实验数据/计算中国H5N1主要候鸟宿主丰度/land_cover_for_china_2000-2023/landcover_classes.xlsx")
landcover_files <- list.files(path = "F:/实验数据/计算中国H5N1主要候鸟宿主丰度/land_cover_for_china_2000-2023",
                              pattern = '.tif', full.names = T)
landcover <- rast(landcover_files)

# Transform buffers to match land cover projection
buffers_for_landcover <- st_transform(buffers_for_landcover, crs = crs(landcover))

# Calculate landscape metrics for each buffer
lsm <- NULL 

for (i in 1:nrow(buffers_for_landcover)) {
  print(paste('process landcover row ', i, sep = ''))
  landcover_buffer_i <- buffers_for_landcover[i, ]
  year <- as.character(landcover_buffer_i$year_lc)
  layer_name <- paste('CLCD_v01_', year, '_albert', sep = '')
  
  # Use tryCatch to handle potential errors
  tryCatch(
    {
      lsm[[i]] <- crop(landcover[[layer_name]], landcover_buffer_i) |> 
        mask(landcover_buffer_i) |> 
        calculate_lsm(level = "class", metric = c("pland", "ed")) |> 
        mutate(checklist_id = landcover_buffer_i$checklist_id, 
               year = landcover_buffer_i$year) |> 
        dplyr::select(checklist_id, year, class, metric, value)
    },
    error = function(e) {
      message(paste("Error processing row", i, ":", e$message))
    }
  )
}

lsm <- bind_rows(lsm)

# Transform landscape metrics to wide format
lsm_wide <- lsm |> 
  complete(nesting(checklist_id, year),
           class = lc_classes$class,
           metric = c("ed", "pland"),
           fill = list(value = 0)) |> 
  inner_join(dplyr::select(lc_classes, class, label), by = "class") |> 
  pivot_wider(values_from = value,
              names_from = c(class, label, metric),
              names_glue = "{metric}_c{str_pad(class, 1)}_{label}",
              names_sort = TRUE) |> 
  arrange(checklist_id, year)

# ------------------------------------------------------------------------------
# 4.2 Extract elevation data
# ------------------------------------------------------------------------------

# Load elevation raster and transform buffers
elevation <- rast("F:/实验数据/计算中国H5N1主要候鸟宿主丰度/1km_global_elevation/1km_global_elevation.tif")
buffers_for_elevation <- st_transform(buffers, crs = crs(elevation))

# Extract mean and standard deviation of elevation within each buffer
elev_buffer <- exact_extract(elevation, buffers_for_elevation, fun = c("mean", "stdev"),
                             progress = FALSE) |> 
  mutate(checklist_id = buffers_for_elevation$checklist_id, year = buffers_for_elevation$year) |> 
  dplyr::select(checklist_id, year, 
                elevation_mean = mean,
                elevation_sd = stdev)

# ------------------------------------------------------------------------------
# 4.3 Extract monthly NDVI (Normalized Difference Vegetation Index)
# ------------------------------------------------------------------------------

# Prepare buffers for monthly variable extraction
checklists_selected_sf <- observation_data |> 
  distinct(checklist_id, year, month, latitude, longitude) |> 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

checklists_selected_sf_buffers <- st_buffer(checklists_selected_sf, dist = set_units(3, "km"))

# Apply year constraints for NDVI data availability
checklists_selected_sf_buffers_ndvi <- checklists_selected_sf_buffers |>
  mutate(year_ndvi = case_when(year < 2000 ~ '2000',
                               year > 2023 ~ '2023',
                               .default = as.character(year)))

ndvi_buffer <- data.frame()

# Extract NDVI values for each buffer
for (i in seq_len(nrow(checklists_selected_sf_buffers_ndvi))) {
  print(paste('process ndvi row ', i, sep = ''))
  checklists_selected_sf_buffers_i <- checklists_selected_sf_buffers_ndvi[i, ]
  ndvi_month_corr_i <- rast(
    list.files(
      path = paste('China_regional_250m_fractional_vegetation_cover_data_set_2000-2023/', checklists_selected_sf_buffers_i$year_ndvi, sep = ''),
      pattern = paste('250m_', checklists_selected_sf_buffers_i$year_ndvi, checklists_selected_sf_buffers_i$month, '_', sep = ''), 
      full.names = TRUE))
  
  checklists_selected_sf_buffers_i <- st_transform(checklists_selected_sf_buffers_i, crs = crs(ndvi_month_corr_i))
  
  ndvi_i <- exact_extract(ndvi_month_corr_i, checklists_selected_sf_buffers_i, fun = c("mean", "stdev"),
                          progress = FALSE) |> 
    mutate(checklist_id = checklists_selected_sf_buffers_i$checklist_id, 
           year = checklists_selected_sf_buffers_i$year,
           month = checklists_selected_sf_buffers_i$month) |> 
    dplyr::select(checklist_id, year, month,
                  ndvi_mean = mean,
                  ndvi_sd = stdev) 
  
  ndvi_buffer <- bind_rows(ndvi_buffer, ndvi_i)
}

# ------------------------------------------------------------------------------
# 4.4 Extract monthly precipitation data
# ------------------------------------------------------------------------------

checklists_selected_sf_buffers_pre <- checklists_selected_sf_buffers |>
  mutate(year_pre = case_when(year < 2000 ~ '2000',
                              year > 2023 ~ '2023',
                              .default = as.character(year)))

precipitation_buffer <- data.frame()

for (i in seq_len(nrow(checklists_selected_sf_buffers_pre))) {
  print(paste('process pre row ', i, sep = ''))
  checklists_selected_sf_buffers_i <- checklists_selected_sf_buffers_pre[i, ]
  file_path <- list.files(
    path = paste('1km_monthly_precipitation_dataset_for_China_2000-2023/', 'pre_', 
                 checklists_selected_sf_buffers_i$year_pre, sep = ''), full.names = T)
  raster_date <- brick(file_path)
  raster_date <- rast(raster_date)
  precipitation_month_corr_i <- raster_date[[paste('X', as.integer(checklists_selected_sf_buffers_i$month),
                                                   sep = '')]]
  
  checklists_selected_sf_buffers_i <- st_transform(checklists_selected_sf_buffers_i, 
                                                   crs = crs(precipitation_month_corr_i))
  
  precipitation_i <- exact_extract(precipitation_month_corr_i, checklists_selected_sf_buffers_i, 
                                   fun = c("mean", "stdev"), progress = FALSE) |> 
    mutate(checklist_id = checklists_selected_sf_buffers_i$checklist_id, 
           year = checklists_selected_sf_buffers_i$year,
           month = checklists_selected_sf_buffers_i$month) |> 
    dplyr::select(checklist_id, year, month,
                  precipitation_mean = mean,
                  precipitation_sd = stdev) 
  
  precipitation_buffer <- bind_rows(precipitation_buffer, precipitation_i)
}

# ------------------------------------------------------------------------------
# 4.5 Extract monthly temperature data
# ------------------------------------------------------------------------------

checklists_selected_sf_buffers_temp <- checklists_selected_sf_buffers |>
  mutate(year_temp = case_when(year < 2000 ~ '2000',
                               year > 2023 ~ '2023',
                               .default = as.character(year)))

temperature_buffer <- data.frame()

for (i in seq_len(nrow(checklists_selected_sf_buffers_temp))) {
  print(paste('process temp row ', i, sep = ''))
  checklists_selected_sf_buffers_i <- checklists_selected_sf_buffers_temp[i, ]
  file_path <- list.files(
    path = paste('1km_monthly_mean_temperature_dataset_for_china_2000-2023/', 'tmp_', 
                 checklists_selected_sf_buffers_i$year_temp, sep = ''), full.names = T)
  raster_date <- brick(file_path)
  raster_date <- rast(raster_date)
  temperature_month_corr_i <- raster_date[[paste('X', as.integer(checklists_selected_sf_buffers_i$month),
                                                 sep = '')]]
  
  checklists_selected_sf_buffers_i <- st_transform(checklists_selected_sf_buffers_i, 
                                                   crs = crs(temperature_month_corr_i))
  
  temperature_i <- exact_extract(temperature_month_corr_i, checklists_selected_sf_buffers_i, 
                                 fun = c("mean", "stdev"), progress = FALSE) |> 
    mutate(checklist_id = checklists_selected_sf_buffers_i$checklist_id, 
           year = checklists_selected_sf_buffers_i$year,
           month = checklists_selected_sf_buffers_i$month) |> 
    dplyr::select(checklist_id, year, month,
                  temperature_mean = mean,
                  temperature_sd = stdev) 
  
  temperature_buffer <- bind_rows(temperature_buffer, temperature_i)
}

# ------------------------------------------------------------------------------
# 4.6 Extract monthly relative humidity data (2003-2020)
# ------------------------------------------------------------------------------

# Apply year constraints for relative humidity data availability (2003-2020)
checklists_selected_sf_buffers_rh <- checklists_selected_sf_buffers |>
  mutate(year_rh = case_when(year < 2003 ~ '2003',
                             year > 2020 ~ '2020',
                             .default = as.character(year)))

rh_buffer <- data.frame() 

for (i in 1:nrow(checklists_selected_sf_buffers_rh)) {
  print(paste('process rh row ', i, sep = ''))
  checklists_selected_sf_buffers_i <- checklists_selected_sf_buffers_rh[i, ]
  file_path <- paste('1km_atmospheric_moisture_index_collection_over_China_2003–2020/', 'HiMIC_Monthly_China_RH_', 
                     checklists_selected_sf_buffers_i$year_rh, '01_', checklists_selected_sf_buffers_i$year_rh, '12_GeoTIFF/',
                     'HiMIC_Monthly_China_RH_', checklists_selected_sf_buffers_i$year_rh, '-', checklists_selected_sf_buffers_i$month,
                     '-01.tif', sep = '')
  
  # Load relative humidity data and convert from percentage to decimal
  rh_month_corr_i <- rast(file_path) / 100
  
  checklists_selected_sf_buffers_i <- st_transform(checklists_selected_sf_buffers_i, 
                                                   crs = crs(rh_month_corr_i))
  
  rh_i <- exact_extract(rh_month_corr_i, checklists_selected_sf_buffers_i, 
                        fun = c("mean", "stdev"), progress = FALSE) |> 
    mutate(checklist_id = checklists_selected_sf_buffers_i$checklist_id, 
           year = checklists_selected_sf_buffers_i$year,
           month = checklists_selected_sf_buffers_i$month) |> 
    dplyr::select(checklist_id, year, month,
                  rh_mean = mean,
                  rh_sd = stdev) 
  
  rh_buffer <- bind_rows(rh_buffer, rh_i)
}

# ==============================================================================
# SECTION 5: Combine all environmental variables and export final dataset
# ==============================================================================

# Merge all environmental variables
env_variables <- inner_join(lsm_wide, elev_buffer, by = c("checklist_id", "year"))
env_variables <- left_join(rh_buffer, env_variables, by = c("checklist_id", "year"))
env_variables <- left_join(env_variables, ndvi_buffer)
env_variables <- left_join(env_variables, temperature_buffer)
env_variables <- left_join(env_variables, precipitation_buffer)

# Join environmental variables with observation data
env_variables <- left_join(observation_data, env_variables, by = c("checklist_id", "year", 'month'))

# Export final dataset
file_path_to_save <- paste("42_birds_observation_data/", bird_name,
                           '/generate_random_false_point_and_extract_all_point_env_envirable/observation_data_env_variables.csv',
                           sep = '')

write_csv(env_variables, file_path_to_save)

