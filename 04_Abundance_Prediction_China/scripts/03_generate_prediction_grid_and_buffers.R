# ==============================================================================
# Circular Buffer Generation for Environmental Variable Extraction
# 
# Description:
# This script generates 3km circular buffers around grid cell centers for 
# extracting environmental variables used in predicting relative abundance of 
# migratory birds. The workflow includes:
# 1. Processing bird observation checklists and assigning unique location IDs
# 2. Creating a prediction grid covering China with 3km resolution
# 3. Generating circular buffers around each grid cell center
# 4. Saving buffers for subsequent environmental variable extraction
#
# Key features:
# - Uses Lambert Azimuthal Equal Area projection optimized for China
# - Handles temporal data range (2000-2023)
# - Creates consistent 3km buffers for spatial analysis
# - Outputs spatial data in standardized formats for further processing
#
# ==============================================================================

# Load required libraries
library(tidyverse)
library(sf)
library(terra)
library(landscapemetrics)
library(units)
library(readxl)
library(exactextractr)

# Read bird observation checklists
checklists <- read_csv("checklists_Lesser_Whistling_Duck.csv")

# Standardize year values to available land cover data range (2000-2023)
checklists <- checklists |> 
  mutate(year_lc = case_when(
    year < 2000 ~ '2000',
    year > 2023 ~ '2023',
    .default = as.character(year)
  ))

# Identify unique geographic locations across all checklists
checklists_distinct_location <- checklists |> 
  distinct(latitude, longitude)

# Assign unique ID to each distinct location
checklists_distinct_location <- checklists_distinct_location |> 
  mutate(locality_id = row.names(checklists_distinct_location))

# Join location IDs back to original checklists
checklists <- left_join(checklists, checklists_distinct_location)

# Generate circular neighborhoods for all checklists
checklists_sf <- checklists |> 
  # Identify unique location/year combinations
  distinct(locality_id, year_lc, latitude, longitude) |> 
  # Convert to spatial points
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Calculate centroid of all observation points for projection optimization
centroid_coords <- checklists_sf |> 
  st_union() |> 
  st_centroid() |> 
  st_coordinates() |> 
  round(1)

# Define Lambert's Azimuthal Equal Area projection optimized for China
# Centered at approximate centroid of observation data (111.4°E, 32°N)
laea_crs <- st_crs("+proj=laea +lat_0=32 +lon_0=111.4")

# Load and project China boundary map
study_region <- read_sf("map_of_China.json") |> 
  st_transform(crs = laea_crs)

# Create raster template covering China with 3000m resolution
r <- rast(study_region, res = c(3000, 3000))

# Rasterize study region, assigning value 1 to cells inside China
r <- rasterize(study_region, r, values = 1) |> 
  setNames("study_region")

# Save prediction grid for later use
writeRaster(r, "3km_prediction_grid.tif", overwrite = TRUE)

# Convert raster to dataframe with cell coordinates and values
r_df <- as.data.frame(r, cells = TRUE, xy = TRUE)

# Generate 3km circular buffers for each prediction grid cell center
buffers_pg <- NULL
for (i in 1:nrow(r_df)) {
  # Print progress indicator
  print(paste('Processing row ', i, sep = ''))
  
  # Process each grid cell
  buffers_pg_i <- r_df[i, ] |> 
    select(cell_id = cell, x, y) |> 
    # Convert to spatial points in LAEA projection
    st_as_sf(coords = c("x", "y"), crs = laea_crs, remove = FALSE) |> 
    # Transform to geographic coordinates for buffer creation
    st_transform(crs = 4326) |> 
    # Create 3km circular buffer around each point
    st_buffer(set_units(3, "km"))
  
  # Combine buffers into single dataset
  buffers_pg <- bind_rows(buffers_pg, buffers_pg_i)
}

# Save buffers to GeoPackage file for environmental variable extraction
st_write(buffers_pg, "buffers_pg.gpkg")


