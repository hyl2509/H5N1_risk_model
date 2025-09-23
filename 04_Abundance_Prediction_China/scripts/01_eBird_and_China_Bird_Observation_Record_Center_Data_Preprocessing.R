# ==============================================================================
# Description: 
#   This script preprocesses bird observation data from two sources:
#   1. eBird (via EBird Basic Dataset - EBD and Sampling Event Data - SED)
#   2. China Bird Observation Record Center (Excel files)
#   The processed data is filtered, cleaned, and combined for downstream analysis.
#   Example species: Lesser_Whistling_Duck (Dendrocygna javanica)
#
# Inputs:
#   - eBird SED files (*_sampling.txt)
#   - eBird EBD files (*-2024.txt)
#   - China Bird Observation Record Center Excel files
#   - China boundary shapefile
#
# Outputs:
#   - Combined checklist CSV file per species
#
# Dependencies:
#   See required libraries below
# ==============================================================================

# Load required libraries
library(auk)      
library(dplyr)    
library(ggplot2) 
library(gridExtra) 
library(lubridate) 
library(readr)     
library(sf)       
library(readxl)    
library(tidyverse)  
library(tidygeocoder) 

# ==============================================================================
# 0. Initial Setup
# ==============================================================================

# Target species and working directory setup
bird_name <- 'Lesser_Whistling_Duck'
work_dir <- paste('42_birds_observation_data/', 
                  bird_name, sep = '')
setwd(work_dir)

# ==============================================================================
# 1. Read Sampling Event Data (SED)
# ==============================================================================

# Locate and read sampling event files
f_sed <- list.files(path = 'raw_data', recursive = TRUE, 
                    full.names = TRUE, pattern = '_sampling.txt')
checklists <- read_sampling(f_sed)

# ==============================================================================
# 2. Read EBD (eBird Basic Dataset)
# ==============================================================================

# Locate and read observation files
f_ebd <- list.files(path = 'raw_data', recursive = TRUE, 
                    full.names = TRUE, pattern = '-2024.txt')
observations <- read_ebd(f_ebd)

# ==============================================================================
# 3. Temporal Filtering (1998-2024)
# ==============================================================================

# Filter checklists to complete reports within date range
checklists <- checklists |> 
  filter(all_species_reported,
         between(year(observation_date), 1998, 2024))

# Filter observations to complete reports within date range
observations <- observations |> 
  filter(all_species_reported,
         between(year(observation_date), 1998, 2024))

# ==============================================================================
# 4. Spatial Filtering (Remove Offshore Checklists)
# ==============================================================================

# Convert checklist locations to spatial points
checklists_sf <- checklists |> 
  select(checklist_id, latitude, longitude) |> 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Load and prepare China boundary with 1km buffer
sf_use_s2(FALSE)
study_region <- read_sf("map_of_China.json") |> 
  st_make_valid() |> 
  st_union() |> 
  st_transform(crs = st_crs(checklists_sf))
sf_use_s2(TRUE)
study_region_buffered <- st_buffer(study_region, dist = 1000)

# Spatially subset checklists within China
in_region <- checklists_sf[study_region_buffered, ]

# Remove checklists outside study region
checklists <- semi_join(checklists, in_region, by = "checklist_id")
observations <- semi_join(observations, in_region, by = "checklist_id")

# Remove observations without matching checklists
observations <- semi_join(observations, checklists, by = "checklist_id")

# ==============================================================================
# 5. Zero-Filling Processing
# ==============================================================================

# Perform zero-filling (presence/absence records)
zf <- auk_zerofill(observations, checklists, collapse = TRUE)

# Helper function: Convert time to decimal hours
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

# Clean and transform variables
zf <- zf |> 
  mutate(
    observation_count = as.integer(observation_count),
    effort_distance_km = if_else(protocol_type == "Stationary", 
                                 0, effort_distance_km),
    effort_hours = duration_minutes / 60,
    effort_speed_kmph = effort_distance_km / effort_hours,
    hours_of_day = time_to_decimal(time_observations_started),
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )

# Apply additional quality filters
zf_filtered <- zf |> 
  filter(protocol_type %in% c("Stationary", "Traveling"),
         effort_hours <= 6,
         effort_distance_km <= 10,
         effort_speed_kmph <= 100,
         number_observers <= 10)

# Prepare final checklist data
checklists <- zf_filtered |> 
  select(checklist_id, 
         observation_count, species_observed, 
         latitude, longitude,
         observation_date, year, day_of_year,
         hours_of_day, 
         effort_hours)

# Keep only checklists with positive observations
checklists <- checklists |> 
  filter(species_observed == 'TRUE')

# ==============================================================================
# 6. Process China Bird Observation Record Center
# ==============================================================================

# Locate and read China Bird Observation Record Center data
birdreport_file <- list.files(path = 'raw_data/observation_from_China_Bird_report', 
                              pattern = '.xlsx', full.names = TRUE)
observations_china <- read_xlsx(birdreport_file)

# Standardize column names and types
observations_china <- observations_china |> 
  rename(observation_count = observation_number) |> 
  select(-c(1:2))

# Convert and calculate time-related variables
observations_china <- observations_china |> 
  mutate(datetime_observations_started = as_datetime(observation_start_time),
         datetime_observations_ended = as_datetime(observation_end_time))

# Create eBird-compatible fields
observations_china <- observations_china |> 
  mutate(
    effort_hours = as.numeric(datetime_observations_ended - datetime_observations_started, 
                              units = "hours"),
    observation_date = date(datetime_observations_started),
    year = year(observation_date),
    time_observations_started = format(datetime_observations_started, "%H:%M:%S"),
    day_of_year = yday(observation_date),
    checklist_id = paste('B', rownames(observations_china), sep = ''),
    species_observed = as.logical('TRUE'),
    hours_of_day = time_to_decimal(time_observations_started),
    observation_count = as.integer(observation_count)
  ) |> 
  filter(effort_hours <= 10)

# Select matching columns with eBird data
observations_china <- observations_china |> 
  select(checklist_id, 
         observation_count, species_observed, 
         latitude, longitude,
         observation_date, year, day_of_year,
         hours_of_day, 
         effort_hours)

# ==============================================================================
# 7. Combine Data Sources and Export
# ==============================================================================

# Merge eBird and China Bird Report data
combined <- bind_rows(checklists, observations_china)

# Write processed data to CSV
file_name <- paste('data_after_raw_data_preprocess/checklists_', bird_name, '.csv', sep = '')
write_csv(drop_na(combined), file_name)



