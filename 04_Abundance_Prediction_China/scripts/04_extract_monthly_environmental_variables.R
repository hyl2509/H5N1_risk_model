# ==============================================================================
# Description: This script extracts monthly environmental variables for each 
#              grid cell across China. It processes 12 months of data including
#              elevation, land cover, relative humidity, NDVI, precipitation, 
#              and temperature. For each month, it calculates mean and standard 
#              deviation within buffer zones around grid cells, and computes 
#              landscape metrics for land cover classes.
# ==============================================================================

# Load required libraries
library(landscapemetrics)
library(sf)
library(terra)
library(tidyverse)
library(ncdf4)
library(raster)
library(exactextractr)
library(readxl)

# ==============================================================================
# FUNCTION: extract_environmental_variables
# Description: Extracts environmental variables for a specific month
# Parameters:
#   - month: Integer (1-12) representing the month to process
#   - buffer_pg: Spatial polygons data frame containing grid cell buffers
# Returns: Data frame with all extracted environmental variables for the month
# ==============================================================================
extract_environmental_variables <- function(month, buffer_pg) {
  
  # Define file paths for environmental variables
  elevation <- rast('environmental_variable_for_january_extract/1km_global_elevation.tif')
  landcover <- rast("environmental_variable_for_january_extract/CLCD_v01_2023_albert.tif")
  rh <- rast("environmental_variable_for_january_extract/HiMIC_Monthly_China_RH_2020-01-01.tif")
  ndvi <- rast("environmental_variable_for_january_extract/HXPT_FVC_MONTH_MAX_250m_202301.tif")
  
  # Read precipitation and temperature data from NetCDF files
  pre <- brick("environmental_variable_for_january_extract/pre_2023.nc") |> rast()
  pre <- pre[[month]]
  
  tmp <- brick("environmental_variable_for_january_extract/tmp_2023.nc") |> rast()
  tmp <- tmp[[month]]
  
  # Transform buffer polygons to match CRS of each environmental variable
  buffer_pg_elevation <- st_transform(buffer_pg, crs = crs(elevation))
  buffer_pg_landcover <- st_transform(buffer_pg, crs = crs(landcover))
  buffer_pg_rh <- st_transform(buffer_pg, crs = crs(rh))
  buffer_pg_ndvi <- st_transform(buffer_pg, crs = crs(ndvi))
  buffer_pg_pre <- st_transform(buffer_pg, crs = crs(pre))
  buffer_pg_tmp <- st_transform(buffer_pg, crs = crs(tmp))
  
  # Initialize lists to store results for each buffer
  elev_buffer <- NULL
  landcover_buffer <- NULL
  rh_buffer <- NULL
  ndvi_buffer <- NULL
  pre_buffer <- NULL
  tmp_buffer <- NULL
  
  # Process each grid cell buffer
  for (i in 1:nrow(buffer_pg)) {
    print(paste('Processing row ', i, sep = ''))
    
    # Extract elevation statistics (mean and standard deviation)
    buffer_pg_elevation_i <- buffer_pg_elevation[i, ]
    elev_buffer[[i]] <- exact_extract(elevation, buffer_pg_elevation_i, 
                                      fun = c("mean", "stdev"), progress = FALSE) |> 
      mutate(cell_id = buffer_pg_elevation_i$cell_id) |> 
      dplyr::select(cell_id,
                    elevation_mean = mean,
                    elevation_sd = stdev)
    
    # Extract landcover metrics (pland and edge density) using landscapemetrics
    buffer_pg_landcover_i <- buffer_pg_landcover[i, ]
    landcover_buffer[[i]] <- crop(landcover, buffer_pg_landcover_i) |> 
      mask(buffer_pg_landcover_i) |> 
      calculate_lsm(level = "class", metric = c("pland", "ed")) |> 
      mutate(cell_id = buffer_pg_landcover_i$cell_id) |> 
      dplyr::select(cell_id, class, metric, value)
    
    # Extract relative humidity statistics
    buffer_pg_rh_i <- buffer_pg_rh[i, ]
    rh_buffer[[i]] <- exact_extract(rh, buffer_pg_rh_i, 
                                    fun = c("mean", "stdev"), progress = FALSE) |> 
      mutate(cell_id = buffer_pg_rh_i$cell_id) |> 
      dplyr::select(cell_id,
                    rh_mean = mean,
                    rh_sd = stdev)
    
    # Extract NDVI statistics
    buffer_pg_ndvi_i <- buffer_pg_ndvi[i, ]
    ndvi_buffer[[i]] <- exact_extract(ndvi, buffer_pg_ndvi_i, 
                                      fun = c("mean", "stdev"), progress = FALSE) |> 
      mutate(cell_id = buffer_pg_ndvi_i$cell_id) |> 
      dplyr::select(cell_id,
                    ndvi_mean = mean,
                    ndvi_sd = stdev) 
    
    # Extract precipitation statistics
    buffer_pg_pre_i <- buffer_pg_pre[i, ]
    pre_buffer[[i]] <- exact_extract(pre, buffer_pg_pre_i, 
                                     fun = c("mean", "stdev"), progress = FALSE) |> 
      mutate(cell_id = buffer_pg_pre_i$cell_id) |> 
      dplyr::select(cell_id,
                    precipitation_mean = mean,
                    precipitation_sd = stdev) 
    
    # Extract temperature statistics
    buffer_pg_tmp_i <- buffer_pg_tmp[i, ]
    tmp_buffer[[i]] <- exact_extract(tmp, buffer_pg_tmp_i, 
                                     fun = c("mean", "stdev"), progress = FALSE) |> 
      mutate(cell_id = buffer_pg_tmp_i$cell_id) |> 
      dplyr::select(cell_id,
                    temperature_mean = mean,
                    temperature_sd = stdev)  
  }
  
  # Combine results from all buffers
  elev_buffer <- bind_rows(elev_buffer)
  landcover_buffer <- bind_rows(landcover_buffer)
  rh_buffer <- bind_rows(rh_buffer)
  ndvi_buffer <- bind_rows(ndvi_buffer)
  pre_buffer <- bind_rows(pre_buffer)
  tmp_buffer <- bind_rows(tmp_buffer)
  
  # Process landcover data: complete missing classes and transform to wide format
  lc_classes <- read_xlsx("environmental_variable_for_january_extract/landcover_classes.xlsx")
  landcover_buffer_wide <- landcover_buffer |> 
    complete(nesting(cell_id),
             class = lc_classes$class,
             metric = c("ed", "pland"),
             fill = list(value = 0)) |> 
    inner_join(dplyr::select(lc_classes, class, label), by = "class") |> 
    pivot_wider(values_from = value,
                names_from = c(class, label, metric),
                names_glue = "{metric}_c{str_pad(class, 1)}_{label}",
                names_sort = TRUE) 
  
  # Combine all environmental variables into final data frame
  buffer_pg_env_variable <- left_join(elev_buffer, landcover_buffer_wide) |> 
    left_join(rh_buffer) |> 
    left_join(ndvi_buffer) |> 
    left_join(pre_buffer) |> 
    left_join(tmp_buffer)
  
  return(buffer_pg_env_variable)
}

# ==============================================================================
# MAIN EXECUTION SECTION
# ==============================================================================

# Read buffer polygons from GeoPackage file
file_path <- "buffers_pg.gpkg"
buffer_pg <- read_sf(file_path)

# Process all 12 months
for (month in 1:12) {
  print(paste("Processing month:", month))
  
  # Extract environmental variables for current month
  monthly_data <- extract_environmental_variables(month, buffer_pg)
  
  # Save results to CSV file
  output_file <- paste0('result_of_prediction_grid_env_variable/buffer_pg_env_variable_month_', 
                        sprintf("%02d", month), '.csv')
  write_csv(monthly_data, output_file)
}

