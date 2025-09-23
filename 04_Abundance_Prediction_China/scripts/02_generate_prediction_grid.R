# ==============================================================================
# Description: This script creates a standardized prediction grid for China
#              to be used in species distribution modeling for H5N1 host birds.
#              It processes bird observation checklists and generates a 5km
#              resolution grid with circular buffers for environmental variable
#              extraction.
#
# Key steps:
# 1. Load and preprocess bird observation checklist data
# 2. Create unique location identifiers and convert to spatial data
# 3. Define equal-area projection centered on China
# 4. Generate 5km resolution prediction grid covering China
# 5. Create circular buffers around grid cell centers for analysis
#
# Inputs: Bird observation checklists, China boundary map
# Outputs: Prediction grid raster and buffer shapefile
# ==============================================================================

# Load required libraries
library(tidyverse)
library(sf)
library(terra)
library(landscapemetrics)
library(units)
library(readxl)
library(exactextractr)

# ==============================================================================
# SECTION 1: Load and preprocess bird observation data
# ==============================================================================

# Read bird observation checklists for Lesser Whistling Duck
checklists <- read_csv(
  "checklists_Lesser_Whistling_Duck.csv")

# Standardize year values for land cover data compatibility
# Constrain years between 2000 and 2023 to match available land cover data
checklists <- checklists |> 
  mutate(year_lc = case_when(
    year < 2000 ~ '2000',
    year > 2023 ~ '2023',
    .default = as.character(year)
  ))

# Identify unique geographic locations from the checklists
checklists_distinct_location <- checklists |> 
  distinct(latitude, longitude)

# Assign unique identifier to each location
checklists_distinct_location <- checklists_distinct_location |> 
  mutate(locality_id = row.names(checklists_distinct_location))

# Join location IDs back to the main checklist data
checklists <- left_join(checklists, checklists_distinct_location)

# ==============================================================================
# SECTION 2: Convert checklist data to spatial format
# ==============================================================================

# Create spatial points from checklist locations
checklists_sf <- checklists |> 
  # Identify unique location/year combinations
  distinct(locality_id, year_lc, latitude, longitude) |> 
  # Convert to spatial points with WGS84 coordinate system
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Calculate centroid of all observation points for projection definition
checklists_sf |> 
  st_union() |> 
  st_centroid() |> 
  st_coordinates() |> 
  round(1)

# ==============================================================================
# SECTION 3: Define projection and study region
# ==============================================================================

# Lambert's Azimuthal Equal Area projection centered on China
# This projection preserves area measurements, important for spatial analysis
laea_crs <- st_crs("+proj=laea +lat_0=32 +lon_0=111.4")

# Load China boundary and transform to equal-area projection
study_region <- read_sf("map_of_China.json") |> 
  st_transform(crs = laea_crs)

# ==============================================================================
# SECTION 4: Create prediction grid
# ==============================================================================

# Create a raster template covering China with 5km resolution
r <- rast(study_region, res = c(5000, 5000))

# Rasterize study region - assign value 1 to cells inside China
r <- rasterize(study_region, r, values = 1) |> 
  setNames("study_region")

# Save prediction grid for later use
r <- writeRaster(r, "5km_prediction_grid.tif",
                 overwrite = TRUE)

# ==============================================================================
# SECTION 5: Create analysis buffers around grid cells
# ==============================================================================

# Generate circular buffers around prediction grid cell centers
# 5km radius buffers will be used for environmental variable extraction
buffers_pg <- as.data.frame(r, cells = TRUE, xy = TRUE) |> 
  select(cell_id = cell, x, y) |> 
  # Convert to spatial points in LAEA projection
  st_as_sf(coords = c("x", "y"), crs = laea_crs, remove = FALSE) |> 
  # Transform to geographic coordinates for compatibility
  st_transform(crs = 4326) |> 
  # Create 5km circular buffers around each grid point
  st_buffer(set_units(5, "km"))

# Save buffers as shapefile for future analysis
write_sf(buffers_pg, "buffers_prediction_grid.shp")

