# ==============================================================================
# Description: This script extracts monthly environmental variables for China 
#              under SSP245 scenario. It processes land cover, elevation, GPP, 
#              temperature, precipitation, relative humidity, and bird abundance
#              data for each grid cell buffer zone. The script calculates mean 
#              and standard deviation for continuous variables, and landscape
#              metrics (pland and edge density) for land cover classes.
# ==============================================================================

# Load required libraries
library(tidyverse)
library(terra)
library(sf)
library(landscapemetrics)
library(exactextractr)

# ==============================================================================
# FUNCTION: extract_environmental_variables
# Description: Extracts environmental variables for a specific year and month under SSP245 scenario
# Parameters:
#   - year: Integer representing the target year
#   - month: Integer (1-12) representing the target month
#   - pg_buffer: Spatial polygons data frame containing grid cell buffers
# Returns: Data frame with all extracted environmental variables
# ==============================================================================
extract_environmental_variables <- function(year, month, pg_buffer) {
  
  # Create month name vector for file naming
  months_vector <- c("January", "February", "March", "April", "May", "June", 
                     "July", "August", "September", "October", "November", "December")
  names(months_vector) <- 1:12
  
  # Extract land cover metrics (pland and edge density)
  landcover <- rast('../predictor/result20.tif')
  lc_classes <- read_csv('../predictor/landcover_classes.csv')
  
  lsm <- NULL
  for (i in seq_len(nrow(pg_buffer))) {
    print(paste('Processing landcover row ', i, sep = ''))
    geom_i <- st_transform(pg_buffer[i, ], crs = crs(landcover))
    
    # Crop and mask landcover raster so all values outside buffer are missing
    lsm[[i]] <- crop(landcover, geom_i) |>
      mask(geom_i) |>
      # Calculate landscape metrics
      calculate_lsm(level = "class", metric = c("pland", "ed")) |>
      # Add variables to uniquely identify each point
      mutate(cell_id = geom_i$cell_id) |>
      select(cell_id, class, metric, value)
  }
  
  lsm <- bind_rows(lsm)
  lsm_wide <- lsm |>
    # Fill missing classes with zeros
    complete(nesting(cell_id),
             class = lc_classes$class,
             metric = c("ed", "pland"),
             fill = list(value = 0)) |>
    # Bring in more descriptive names
    inner_join(select(lc_classes, class, name), by = "class") |>
    # Transform from long to wide format
    pivot_wider(values_from = value,
                names_from = c(class, name, metric),
                names_glue = "{metric}_c{str_pad(class, 2, pad = '0')}_{name}",
                names_sort = TRUE)
  
  # Extract elevation statistics
  elevation <- rast('../predictor/wc2.1_30s_elev.tif')
  geom_elev <- exact_extract(elevation, pg_buffer, fun = c("mean", "stdev"),
                             progress = FALSE) |>
    # Add variables to uniquely identify each point
    mutate(cell_id = pg_buffer$cell_id) |>
    select(cell_id,
           elevation_mean = mean,
           elevation_sd = stdev)
  
  # Extract GPP statistics for specified year and month
  geom_time_gpp <- NULL
  gpp_file <- paste0('../predictor/GPP_mon_ssp245_', year, '.tif')
  gpp <- rast(gpp_file)
  gpp <- gpp[[month]]
  
  for (i in seq(nrow(pg_buffer))) {
    print(paste('Processing gpp row ', i, sep = ' '))
    geom_time_i <- pg_buffer[i, ]
    geom_time_i_gpp <- exact_extract(gpp, geom_time_i, fun = c("mean", "stdev"),
                                     progress = FALSE) |>
      # Add variables to uniquely identify each point
      mutate(cell_id = geom_time_i$cell_id) |>
      select(cell_id,
             gpp_mean = mean,
             gpp_sd = stdev)
    geom_time_gpp <- bind_rows(geom_time_gpp, geom_time_i_gpp)
  }
  
  # Extract temperature statistics for specified year and month
  geom_time_temp <- NULL
  temp_file <- paste0('../predictor/EC-Earth3_ssp245_tmp-30s-', year, '.tif')
  temp <- rast(temp_file)
  temp <- temp[[month]]
  
  for (i in seq(nrow(pg_buffer))) {
    print(paste('Processing temp row ', i, sep = ''))
    geom_time_for_temp_i <- pg_buffer[i, ]
    
    i_temp <- exact_extract(temp, geom_time_for_temp_i, fun = c("mean", "stdev"),
                            progress = FALSE) |>
      mutate(cell_id = geom_time_for_temp_i$cell_id) |>
      select(cell_id,
             temp_mean = mean,
             temp_sd = stdev)
    geom_time_temp <- bind_rows(geom_time_temp, i_temp)
  }
  
  # Extract precipitation statistics for specified year and month
  geom_time_pre <- NULL
  pre_file <- paste0('../predictor/EC-Earth3_ssp245_pre-30s-', year, '.tif')
  pre <- rast(pre_file)
  pre <- pre[[month]]
  
  for (i in seq(nrow(pg_buffer))) {
    print(paste('Processing pre row ', i, sep = ''))
    geom_time_for_pre_i <- pg_buffer[i, ]
    i_pre <- exact_extract(pre, geom_time_for_pre_i, fun = c("mean", "stdev"),
                           progress = FALSE) |>
      mutate(cell_id = geom_time_for_pre_i$cell_id) |>
      select(cell_id,
             pre_mean = mean,
             pre_sd = stdev)
    geom_time_pre <- bind_rows(geom_time_pre, i_pre)
  }
  
  # Extract relative humidity statistics for specified year and month
  geom_time_rh <- NULL
  rh_file <- paste0('../predictor/hurs_Amon_EC-Earth3-HR_ssp245_r1i1p1f1_gr_', 
                    year, '01-', year, '12.tif')
  rh <- rast(rh_file)
  rh <- rh[[month]]
  
  for (i in seq(nrow(pg_buffer))) {
    print(paste('Processing rh row ', i, sep = ''))
    geom_time_for_rh_i <- pg_buffer[i, ]
    
    i_rh <- exact_extract(rh, geom_time_for_rh_i, fun = c("mean", "stdev"),
                          progress = FALSE) |>
      mutate(cell_id = geom_time_for_rh_i$cell_id) |>
      select(cell_id,
             rh_mean = mean,
             rh_sd = stdev)
    geom_time_rh <- bind_rows(geom_time_rh, i_rh)
  }
  
  # Extract bird abundance statistics for specified month
  geom_time_bird_abundance <- NULL
  abundance <- rast(paste('../predictor/Sum_of_bird_abundance_in_',
                          months_vector[[month]], '.tif', sep = ''))
  
  for (i in seq(nrow(pg_buffer))) {
    print(paste('Processing bird row ', i, sep = ''))
    geom_time_for_abundance_i <- pg_buffer[i, ]
    
    geom_time_for_abundance_i <- st_transform(geom_time_for_abundance_i, crs = crs(abundance))
    i_abundance <- exact_extract(abundance, geom_time_for_abundance_i, fun = c("mean", "stdev"),
                                 progress = FALSE) |>
      mutate(cell_id = geom_time_for_abundance_i$cell_id) |>
      select(cell_id,
             bird_abundance_mean = mean,
             bird_abundance_sd = stdev)
    geom_time_bird_abundance <- bind_rows(geom_time_bird_abundance, i_abundance)
  }
  
  # Combine all environmental variables into final data frame
  pg_buffer_env_variables <- left_join(lsm_wide, geom_elev)
  pg_buffer_env_variables <- left_join(pg_buffer_env_variables, geom_time_gpp)
  pg_buffer_env_variables <- left_join(pg_buffer_env_variables, geom_time_temp)
  pg_buffer_env_variables <- left_join(pg_buffer_env_variables, geom_time_pre)
  pg_buffer_env_variables <- left_join(pg_buffer_env_variables, geom_time_rh)
  pg_buffer_env_variables <- left_join(pg_buffer_env_variables, geom_time_bird_abundance)
  
  return(pg_buffer_env_variables)
}

# ==============================================================================
# MAIN EXECUTION SECTION
# ==============================================================================

# Read buffer polygons
pg_buffer <- read_sf('../pg_buffer/buffers_pg.gpkg')

# Define target years and months
target_years <- c(2025, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100)
months_vector <- c("January", "February", "March", "April", "May", "June", 
                   "July", "August", "September", "October", "November", "December")

# Process all years and months
for (year in target_years) {
  for (month in 1:12) {
    print(paste("Processing year:", year, "month:", months_vector[month]))
    
    # Extract environmental variables for current year and month
    monthly_data <- extract_environmental_variables(year, month, pg_buffer)
    
    # Save results to CSV file
    output_file <- paste0('ssp245_env_variables_', year, '_', 
                          sprintf("%02d", month), '_', months_vector[month], '.csv')
    write_csv(monthly_data, output_file)
    
    print(paste("Completed year:", year, "month:", months_vector[month]))
  }
}

