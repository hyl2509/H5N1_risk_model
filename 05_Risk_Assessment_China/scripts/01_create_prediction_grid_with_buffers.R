# ==============================================================================
# Description: This script creates a 30km resolution prediction grid covering China
#              and generates 30km buffer zones around each grid cell center.
#              The grid and buffers are used for extracting environmental variables
#              in H5N1 risk prediction modeling.
# Steps:
#   1. Read China boundary shapefile
#   2. Create 30km resolution raster template
#   3. Rasterize China boundary to create prediction grid
#   4. Generate 30km buffer zones around each grid cell center
#   5. Save results as GeoTIFF and GeoPackage files
# ==============================================================================

# Load required libraries
library(sf)
library(terra)
library(tidyverse)
library(units)

# Set working directory
setwd("model_projection/generate_China_30km_prediction_grid")

# Read China boundary shapefile
china_sf <- read_sf('map_of_China.json')

# Create a raster template with 30km resolution (approximately 0.26949462 degrees)
raster_template <- rast(ext(china_sf), resolution = 0.26949462,
                        crs = crs(china_sf))

# Rasterize the China boundary to create prediction grid
raster_output <- terra::rasterize(china_sf, raster_template, touches = TRUE) |> 
  setNames("projection_region")

# Save the prediction grid as GeoTIFF
writeRaster(raster_output, '30km_prediction_grid.tif')

# Generate 30km buffer zones around each grid cell center for environmental variable extraction
r_df <- as.data.frame(raster_output, cells = TRUE, xy = TRUE)

buffers_pg <- NULL 
for (i in 1:nrow(r_df)) {
  print(paste('Processing row ', i, sep = ''))
  buffers_pg_i <- r_df[i, ] |> 
    select(cell_id = cell, x, y) |> 
    st_as_sf(coords = c("x", "y"), crs = 4326, remove = FALSE) |> 
    st_buffer(set_units(30, "km"))
  buffers_pg <- bind_rows(buffers_pg, buffers_pg_i)
}

# Save buffer zones as GeoPackage file
st_write(buffers_pg, 'prediction_grid_point_buffer/buffers_pg.gpkg')


