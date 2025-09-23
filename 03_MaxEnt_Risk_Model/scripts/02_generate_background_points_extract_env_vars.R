# H5N1 Background Point Generation and Environmental Variable Extraction
#
# Purpose: This script generates pseudo-absence (background) points for H5N1 modeling by:
# 1. Calculating average US county area to determine buffer radius (30km)
# 2. Generating random points across the US outside H5N1 detection areas
# 3. Extracting environmental variables (land cover, elevation, GPP, temperature,
#    precipitation, humidity, and bird abundance) for each background point
# 
# Inputs:
# - US county shapefile
# - USA H5N1 surveillance data (GPKG format)
# - Global land cover data (2020)
# - Elevation data (WorldClim)
# - Monthly GPP, temperature, precipitation, humidity rasters
# - Monthly bird abundance rasters
#
# Outputs:
# - CSV files with background points and extracted environmental variables

# Load required libraries
library(terra)
library(sf)
library(ncdf4)
library(tidyverse)
library(units)
library(rnaturalearth)
library(landscapemetrics)
library(exactextractr)
library(doParallel)
library(foreach)

# Set seed for reproducibility
set.seed(2025)

# ==============================================================================
# SECTION 1: Calculate average US county area to determine buffer radius
# ==============================================================================

# Load US county boundaries
usa_county <- read_sf('map_data/usa_county_map/cb_2022_us_county_500k.shp')

# Remove non-county entities (cities that are independent of counties)
usa_county <- usa_county |> 
  filter(!NAMELSAD %in% c('St. Louis city', 'Baltimore city', 'Roanoke city', 'Richmond city'))

# Transform to standard coordinate reference system
usa_county <- st_transform(usa_county, 4326)

# Calculate area for each county in square kilometers
usa_county$area_km2 <- set_units(st_area(usa_county), km^2)

# Calculate mean county area to determine appropriate buffer radius
mean_area <- mean(usa_county$area_km2)
# Based on mean area, a 30km buffer radius was determined to be appropriate

# ==============================================================================
# SECTION 2: Generate random background points
# ==============================================================================

# Load base map and H5N1 occurrence data
usa_map <- ne_countries(scale = 110, country = 'United States of America')
h5n1_info_sf <- read_sf('usa_h5n1_avian_surveillance_processed.gpkg')
h5n1_info_sf$Collection_date <- as.Date(h5n1_info_sf$Collection_date, '%m/%d/%Y')

# Generate random dates within the range of H5N1 detection dates
start_date <- as.Date(min(h5n1_info_sf$Collection_date))
end_date <- as.Date(max(h5n1_info_sf$Collection_date))
random_dates <- sample(seq.Date(start_date, end_date, by = "day"), size = 6000, replace = TRUE)

# Create data frame with random dates and their frequencies
background_point_date <- data.frame(Collection_date = random_dates)
background_point_unique_date <- background_point_date |> 
  count(Collection_date)

# Generate random points using parallel processing
# Points are generated outside spatiotemporal buffers around H5N1 detections
nondetection_random_point <- data.frame()

# Set up parallel processing
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Generate points in parallel
nondetection_random_point <- foreach(
  i = 1:nrow(background_point_unique_date),
  .combine = bind_rows,
  .packages = c("sf", "dplyr"),
  .export = c("background_point_unique_date", "h5n1_info_sf", "usa_map")
) %dopar% {
  # Get current date and number of points to generate
  current_row <- background_point_unique_date[i, ]
  
  # Find H5N1 detections within ±30 days of current date
  h5n1_subset <- h5n1_info_sf |> 
    filter(Collection_date >= current_row$Collection_date - 30, 
           Collection_date <= current_row$Collection_date + 30)
  
  # Create non-detection area by excluding areas around H5N1 detections
  if (nrow(h5n1_subset) > 0) {
    combined_area <- st_union(h5n1_subset)
    non_detection_area <- st_difference(usa_map, combined_area)
  } else {
    non_detection_area <- usa_map
  }
  
  # Generate random points in non-detection area
  random_points <- st_sample(non_detection_area, size = current_row$n)
  random_points_sf <- st_as_sf(random_points) |> 
    st_transform(4326)
  
  # Extract coordinates and add date information
  coordinates <- st_coordinates(random_points_sf) |> 
    as.data.frame() |> 
    setNames(c("longitude", "latitude")) |> 
    mutate(Collection_date = current_row$Collection_date)
  
  return(coordinates)
}

# Close parallel cluster
stopCluster(cl)

# Convert to spatial object
nondetection_random_point_sf <- st_as_sf(nondetection_random_point,
                                         coords = c('longitude', 'latitude'), 
                                         crs = 4326, remove = F)

# Create 30km buffers around each point for environmental variable extraction
nondetection_random_point_sf_buffers <- st_buffer(nondetection_random_point_sf, dist = set_units(30, "km"))

# ==============================================================================
# SECTION 3: Extract environmental variables for background points
# ==============================================================================

# ------------------------------------------------------------------------------
# 3.1 Extract land cover data and calculate landscape metrics
# ------------------------------------------------------------------------------

landcover <- rast('result20.tif')
lc_classes <- read_csv('landcover_classes.csv')

lsm <- NULL
for (i in seq_len(nrow(nondetection_random_point_sf_buffers))) {
  print(paste('process row', i, sep = ''))
  geom_i <- st_transform(nondetection_random_point_sf_buffers[i, ], crs = crs(landcover))
  
  # Crop and mask landcover raster to buffer area
  lsm[[i]] <- crop(landcover, geom_i) |> 
    mask(geom_i) |> 
    # Calculate landscape metrics (percentage of landscape and edge density)
    calculate_lsm(level = "class", metric = c("pland", "ed")) |> 
    # Add identifying variables
    mutate(longitude = geom_i$longitude, 
           latitude = geom_i$latitude,
           Collection_date = geom_i$Collection_date) |> 
    select(longitude, latitude, Collection_date, class, metric, value)
}

lsm <- bind_rows(lsm)

# Transform landscape metrics to wide format
lsm_wide <- lsm |> 
  # Fill missing classes with zeros
  complete(nesting(longitude, latitude, Collection_date),
           class = lc_classes$class,
           metric = c("ed", "pland"),
           fill = list(value = 0)) |> 
  # Join with descriptive class names
  inner_join(select(lc_classes, class, name), by = "class") |> 
  # Convert from long to wide format
  pivot_wider(values_from = value,
              names_from = c(class, name, metric),
              names_glue = "{metric}_c{str_pad(class, 2, pad = '0')}_{name}",
              names_sort = TRUE) 

# ------------------------------------------------------------------------------
# 3.2 Extract elevation data
# ------------------------------------------------------------------------------

elevation <- rast('raw_elevation_data_from_worldclim/wc2.1_30s_elev.tif')

geom_elev <- exact_extract(elevation, nondetection_random_point_sf_buffers, 
                           fun = c("mean", "stdev"), progress = FALSE) |> 
  # Add identifying variables
  mutate(longitude = nondetection_random_point_sf_buffers$longitude, 
         latitude = nondetection_random_point_sf_buffers$latitude,
         Collection_date = nondetection_random_point_sf_buffers$Collection_date) |> 
  select(longitude, latitude, Collection_date,
         elevation_mean = mean,
         elevation_sd = stdev)

# ------------------------------------------------------------------------------
# 3.3 Extract Gross Primary Productivity (GPP) data
# ------------------------------------------------------------------------------

geom_time_gpp <- NULL

for (i in seq(nrow(nondetection_random_point_sf_buffers))) {
  print(paste('process gpp row ', i, sep = ''))
  geom_time_i <- nondetection_random_point_sf_buffers[i, ]
  
  # Construct GPP file name based on date
  gpp_name <- paste('monthly_gpp_2021-2024/processed_gpp/GOSIF_GPP_',
                    if_else(year(geom_time_i$Collection_date) == 2025, 2024, year(geom_time_i$Collection_date)), 
                    '.M', str_pad(month(geom_time_i$Collection_date), 2, pad = '0'), '_Mean.tif', sep = '')
  
  gpp <- rast(gpp_name)
  
  geom_time_i_gpp <- exact_extract(gpp, geom_time_i, fun = c("mean", "stdev"), progress = FALSE) |> 
    # Add identifying variables
    mutate(longitude = geom_time_i$longitude, 
           latitude = geom_time_i$latitude,
           Collection_date = geom_time_i$Collection_date) |> 
    select(longitude, latitude, Collection_date,
           gpp_mean = mean,
           gpp_sd = stdev)
  geom_time_gpp <- bind_rows(geom_time_gpp, geom_time_i_gpp)
}

# ------------------------------------------------------------------------------
# 3.4 Extract temperature data (parallel processing)
# ------------------------------------------------------------------------------

n_cores <- detectCores() - 1 
cl <- makeCluster(n_cores)
registerDoParallel(cl)

geom_time_temp_list <- vector("list", nrow(nondetection_random_point_sf_buffers))

geom_time_temp_list <- foreach(
  i = seq_len(nrow(nondetection_random_point_sf_buffers)),
  .packages = c("sf", "terra", "exactextractr", "dplyr", "lubridate"),
  .export = c("nondetection_random_point_sf_buffers"),
  .combine = "c",  # Combine sublists into one list
  .inorder = TRUE  # Maintain original order
) %dopar% {
  # Process single geometry
  geom_i <- nondetection_random_point_sf_buffers[i, ]
  
  # Temperature data processing
  result <- tryCatch({
    temp_name <- paste0('monthly_temperature_2021-2025/temp_',
                        year(geom_i$Collection_date), '.tif')
    temp_month <- month(geom_i$Collection_date)
    temp <- rast(temp_name)[[temp_month]]
    temp_celsius <- temp - 273.15  # Convert from Kelvin to Celsius
    
    # Coordinate transformation and value extraction
    geom_transformed <- st_transform(geom_i, crs = crs(temp_celsius))
    exact_extract(temp_celsius, geom_transformed, 
                  fun = c("mean", "stdev"),
                  progress = FALSE) |> 
      mutate(longitude = geom_i$longitude, 
             latitude = geom_i$latitude,
             Collection_date = geom_i$Collection_date) |> 
      select(longitude, latitude, Collection_date,
             temp_mean = mean,
             temp_sd = stdev)
  }, error = function(e) {
    message("Error in row ", i, ": ", e$message)
    return(NULL)
  })
  
  # Return single-element list
  list(result)
}

# Close parallel cluster
stopCluster(cl)

geom_time_temp <- bind_rows(geom_time_temp_list)

# ------------------------------------------------------------------------------
# 3.5 Extract precipitation data
# ------------------------------------------------------------------------------

geom_time_pre <- NULL

for (i in seq(nrow(nondetection_random_point_sf_buffers))) {
  print(paste('process row', i, sep = ''))
  geom_time_for_pre_i <- nondetection_random_point_sf_buffers[i, ]
  
  # Construct precipitation file name
  pre_name <- paste('monthly_precipitation_2021-2025/precipitation_',
                    year(geom_time_for_pre_i$Collection_date), '.tif', sep = '')
  pre_month <- month(geom_time_for_pre_i$Collection_date)
  pre <- rast(pre_name)[[pre_month]]
  
  geom_time_for_pre_i <- st_transform(geom_time_for_pre_i, crs = crs(pre))
  i_pre <- exact_extract(pre, geom_time_for_pre_i, fun = c("mean", "stdev"), progress = FALSE) |> 
    mutate(longitude = geom_time_for_pre_i$longitude, 
           latitude = geom_time_for_pre_i$latitude,
           Collection_date = geom_time_for_pre_i$Collection_date) |> 
    select(longitude, latitude, Collection_date,
           pre_mean = mean,
           pre_sd = stdev)
  geom_time_pre <- bind_rows(geom_time_pre, i_pre)
}

# ------------------------------------------------------------------------------
# 3.6 Extract relative humidity data
# ------------------------------------------------------------------------------

geom_time_rh <- NULL

for (i in seq(nrow(nondetection_random_point_sf_buffers))) {
  print(paste('process row rh ', i, sep = ''))
  geom_time_for_rh_i <- nondetection_random_point_sf_buffers[i, ]
  
  # Construct relative humidity file name
  rh_name <- paste('calculated_monthly_relative_humidity_2021-2025/relative_humidity_',
                   year(geom_time_for_rh_i$Collection_date), '.tif', sep = '')
  rh_month <- month(geom_time_for_rh_i$Collection_date)
  rh <- rast(rh_name)[[rh_month]]
  
  geom_time_for_rh_i <- st_transform(geom_time_for_rh_i, crs = crs(rh))
  i_rh <- exact_extract(rh, geom_time_for_rh_i, fun = c("mean", "stdev"), progress = FALSE) |> 
    mutate(longitude = geom_time_for_rh_i$longitude, 
           latitude = geom_time_for_rh_i$latitude,
           Collection_date = geom_time_for_rh_i$Collection_date) |> 
    select(longitude, latitude, Collection_date,
           rh_mean = mean,
           rh_sd = stdev)
  geom_time_rh <- bind_rows(geom_time_rh, i_rh)
}

# ------------------------------------------------------------------------------
# 3.7 Extract migratory bird abundance data
# ------------------------------------------------------------------------------

# Create month name mapping vector
months_vector <- c("January", "February", "March", "April", "May", "June", 
                   "July", "August", "September", "October", "November", "December")
names(months_vector) <- 1:12

geom_time_bird_abundance <- NULL

for (i in seq(nrow(nondetection_random_point_sf_buffers))) {
  print(paste('process row', i, sep = ''))
  geom_time_for_abundacne_i <- nondetection_random_point_sf_buffers[i, ]
  
  # Construct bird abundance file name based on month
  abundacne_name <- paste('The_sum_of_all_species_abundances_per_month/Sum_of_bird_abundance_in_',
                          months_vector[[month(geom_time_for_abundacne_i$Collection_date)]], '.tif', sep = '')
  
  abundacne <- rast(abundacne_name)
  
  geom_time_for_abundacne_i <- st_transform(geom_time_for_abundacne_i, crs = crs(abundacne))
  i_abundacne <- exact_extract(abundacne, geom_time_for_abundacne_i, fun = c("mean", "stdev"), progress = FALSE) |> 
    mutate(longitude = geom_time_for_abundacne_i$longitude, 
           latitude = geom_time_for_abundacne_i$latitude,
           Collection_date = geom_time_for_abundacne_i$Collection_date) |> 
    select(longitude, latitude, Collection_date,
           bird_abundance_mean = mean,
           bird_abundance_sd = stdev)
  geom_time_bird_abundance <- bind_rows(geom_time_bird_abundance, i_abundacne)
}

# ==============================================================================
# SECTION 4: Combine all environmental variables and export results
# ==============================================================================

# Merge all environmental variables with the background points data
nondetection_random_point <- left_join(nondetection_random_point, lsm_wide)
nondetection_random_point <- left_join(nondetection_random_point, geom_elev)
nondetection_random_point <- left_join(nondetection_random_point, geom_time_gpp)
nondetection_random_point <- left_join(nondetection_random_point, geom_time_temp)
nondetection_random_point <- left_join(nondetection_random_point, geom_time_pre)
nondetection_random_point <- left_join(nondetection_random_point, geom_time_rh)
nondetection_random_point <- left_join(nondetection_random_point, geom_time_bird_abundance)

# Add species identifier for background points
nondetection_random_point$species <- 'background'
nondetection_random_point <- nondetection_random_point |> 
  relocate(species)

# Remove geometry column and export final training dataset
nondetection_random_point <- nondetection_random_point[, -4]
write_csv(drop_na(nondetection_random_point),
          'random_background_points_env_vars_for_training.csv')





