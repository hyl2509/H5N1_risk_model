# H5N1 Environmental Variable Extraction Script
# 
# Purpose: This script processes H5N1 surveillance data from wild birds in the US,
# extracting various environmental variables (land cover, elevation, GPP, temperature,
# precipitation, humidity, and bird abundance) for each sampling location.
# The output is used for MaxEnt modeling to predict H5N1 high-risk areas in China.
#
# Inputs:
# - USA H5N1 surveillance data (GPKG format)
# - Global land cover data (2020)
# - Elevation data (WorldClim)
# - Monthly GPP, temperature, precipitation, humidity rasters
# - Monthly bird abundance rasters
#
# Outputs:
# - CSV file for model training


# Load required packages
library(sf)           
library(terra)       
library(tidyverse)    
library(landscapemetrics) 
library(exactextractr) 


# Load and preprocess H5N1 surveillance data
h5n1_info_sf <- read_sf('usa_h5n1_avian_surveillance_processed.gpkg')
h5n1_info_sf$Collection_date <- as.Date(h5n1_info_sf$Collection_date, "%m/%d/%Y")
h5n1_info_sf$year <- year(h5n1_info_sf$Collection_date)
h5n1_info_sf$month <- str_pad(month(h5n1_info_sf$Collection_date), 2, pad = '0')

# 1. Extract unique geographical locations
# For time-independent variables
unique_geom <- h5n1_info_sf |> 
  distinct(State, County, geom)

# For time-dependent variables
unique_geom_time <- h5n1_info_sf |> 
  distinct(State, County, year, month, geom)

# 2. Extract environmental variables for each location

## Land cover data extraction
landcover <- rast('result20.tif')
lc_classes <- read_csv('landcover_classes.csv')

lsm <- NULL
for (i in seq_len(nrow(unique_geom))) {
  print(paste('process row', i, sep = ''))
  geom_i <- st_transform(unique_geom[i, ], crs = crs(landcover))
  
  # Process landcover data for current location
  lsm[[i]] <- crop(landcover, geom_i) |> 
    mask(geom_i) |> 
    calculate_lsm(level = "class", metric = c("pland", "ed")) |> 
    mutate(State = geom_i$State, 
           County = geom_i$County) |> 
    select(State, County, class, metric, value)
}

lsm <- bind_rows(lsm)

# Transform landscape metrics to wide format
lsm_wide <- lsm |> 
  complete(nesting(State, County),
           class = lc_classes$class,
           metric = c("ed", "pland"),
           fill = list(value = 0)) |> 
  inner_join(select(lc_classes, class, name), by = "class") |> 
  pivot_wider(values_from = value,
              names_from = c(class, name, metric),
              names_glue = "{metric}_c{str_pad(class, 2, pad = '0')}_{name}",
              names_sort = TRUE) |> 
  arrange(State, County)

## Elevation data extraction
elevation <- rast('wc2.1_30s_elev.tif')

unique_geom_elev <- exact_extract(elevation, unique_geom, fun = c("mean", "stdev"),
                                  progress = FALSE) |> 
  mutate(State = unique_geom$State, 
         County = unique_geom$County) |> 
  select(State, County, 
         elevation_mean = mean,
         elevation_sd = stdev)

## Gross Primary Productivity (GPP) extraction
unique_geom_time_gpp <- NULL

for (i in seq(nrow(unique_geom_time))) {
  print(paste('process gpp row ', i, sep = ''))
  unique_geom_time_i <- unique_geom_time[i, ]
  gpp_name <- paste('monthly_gpp_2021-2024/processed_gpp/GOSIF_GPP_',
                    if_else(unique_geom_time_i$year == 2025, 2024, unique_geom_time_i$year), '.M', unique_geom_time_i$month, '_Mean.tif', sep = '')
  gpp <- rast(gpp_name)
  
  unique_geom_time_i_gpp <- exact_extract(gpp, unique_geom_time_i, fun = c("mean", "stdev"),
                                          progress = FALSE) |> 
    mutate(State = unique_geom_time_i$State, 
           County = unique_geom_time_i$County,
           year = unique_geom_time_i$year,
           month = unique_geom_time_i$month) |> 
    select(State, County, year, month,
           gpp_mean = mean,
           gpp_sd = stdev)
  unique_geom_time_gpp <- bind_rows(unique_geom_time_gpp, unique_geom_time_i_gpp)
}

## Temperature data extraction (converted from Kelvin to Celsius)
unique_geom_time_temp <- NULL
unique_geom_time_for_temp <- unique_geom_time |> 
  mutate(month = as.numeric(month))

for (i in seq(nrow(unique_geom_time_for_temp))) {
  print(paste('process row', i, sep = ''))
  unique_geom_time_for_temp_i <- unique_geom_time_for_temp[i, ]
  temp_name <- paste('monthly_temperature_2021-2025/temp_',
                     unique_geom_time_for_temp_i$year, '.tif', sep = '')
  temp_month <- unique_geom_time_for_temp_i$month
  temp <- rast(temp_name)[[temp_month]]
  temp_celsius <- temp - 273.15  # Convert from Kelvin to Celsius
  unique_geom_time_for_temp_i <- st_transform(unique_geom_time_for_temp_i, crs = crs(temp_celsius))
  i_temp <- exact_extract(temp_celsius, unique_geom_time_for_temp_i, fun = c("mean", "stdev"),
                          progress = FALSE) |> 
    mutate(State = unique_geom_time_for_temp_i$State, 
           County = unique_geom_time_for_temp_i$County,
           year = unique_geom_time_for_temp_i$year,
           month = unique_geom_time_for_temp_i$month) |> 
    select(State, County, year, month,
           temp_mean = mean,
           temp_sd = stdev)
  unique_geom_time_temp <- bind_rows(unique_geom_time_temp, i_temp)
}

## Precipitation data extraction
unique_geom_time_pre <- NULL
unique_geom_time_for_pre <- unique_geom_time |> 
  mutate(month = as.numeric(month))

for (i in seq(nrow(unique_geom_time_for_pre))) {
  print(paste('process row', i, sep = ''))
  unique_geom_time_for_pre_i <- unique_geom_time_for_pre[i, ]
  pre_name <- paste('monthly_precipitation_2021-2025/precipitation_',
                    unique_geom_time_for_pre_i$year, '.tif', sep = '')
  pre_month <- unique_geom_time_for_pre_i$month
  pre <- rast(pre_name)[[pre_month]]
  
  unique_geom_time_for_pre_i <- st_transform(unique_geom_time_for_pre_i, crs = crs(pre))
  i_pre <- exact_extract(pre, unique_geom_time_for_pre_i, fun = c("mean", "stdev"),
                         progress = FALSE) |> 
    mutate(State = unique_geom_time_for_pre_i$State, 
           County = unique_geom_time_for_pre_i$County,
           year = unique_geom_time_for_pre_i$year,
           month = unique_geom_time_for_pre_i$month) |> 
    select(State, County, year, month,
           pre_mean = mean,
           pre_sd = stdev)
  unique_geom_time_pre <- bind_rows(unique_geom_time_pre, i_pre)
}

## Relative humidity data extraction (continued)
rh <- rast(rh_name)[[rh_month]]

unique_geom_time_for_rh_i <- st_transform(unique_geom_time_for_rh_i, crs = crs(rh))
i_rh <- exact_extract(rh, unique_geom_time_for_rh_i, fun = c("mean", "stdev"),
                      progress = FALSE) |> 
  mutate(State = unique_geom_time_for_rh_i$State, 
         County = unique_geom_time_for_rh_i$County,
         year = unique_geom_time_for_rh_i$year,
         month = unique_geom_time_for_rh_i$month) |> 
  select(State, County, year, month,
         rh_mean = mean,
         rh_sd = stdev)
unique_geom_time_rh <- bind_rows(unique_geom_time_rh, i_rh)
}

## Migratory bird abundance data extraction
# Create month name vector for file naming
months_vector <- c("January", "February", "March", "April", "May", "June", 
                   "July", "August", "September", "October", "November", "December")
names(months_vector) <- 1:12

unique_geom_time_bird_abundance <- NULL

for (i in seq(nrow(unique_geom_time_for_pre))) {
  print(paste('process row', i, sep = ''))
  unique_geom_time_for_pre_i <- unique_geom_time_for_pre[i, ]
  pre_name <- paste('The_sum_of_all_species_abundances_per_month/Sum_of_bird_abundance_in_',
                    months_vector[[unique_geom_time_for_pre_i$month]], '.tif', sep = '')
  
  pre <- rast(pre_name)
  
  unique_geom_time_for_pre_i <- st_transform(unique_geom_time_for_pre_i, crs = crs(pre))
  i_pre <- exact_extract(pre, unique_geom_time_for_pre_i, fun = c("mean", "stdev"),
                         progress = FALSE) |> 
    mutate(State = unique_geom_time_for_pre_i$State, 
           County = unique_geom_time_for_pre_i$County,
           year = unique_geom_time_for_pre_i$year,
           month = unique_geom_time_for_pre_i$month) |> 
    select(State, County, year, month,
           bird_abundance_mean = mean,
           bird_abundance_sd = stdev)
  unique_geom_time_bird_abundance <- bind_rows(unique_geom_time_bird_abundance, i_pre)
}

# 3. Merge all environmental variables with original H5N1 data
h5n1_info_sf <- left_join(h5n1_info_sf, lsm_wide)
h5n1_info_sf <- left_join(h5n1_info_sf, unique_geom_time_gpp)
h5n1_info_sf <- left_join(h5n1_info_sf, unique_geom_elev)
h5n1_info_sf$month <- as.numeric(h5n1_info_sf$month)
h5n1_info_sf <- left_join(h5n1_info_sf, unique_geom_time_temp)
h5n1_info_sf <- left_join(h5n1_info_sf, unique_geom_time_pre)
h5n1_info_sf <- left_join(h5n1_info_sf, unique_geom_time_rh)
h5n1_info_sf <- left_join(h5n1_info_sf, unique_geom_time_bbird_abundance)
h5n1_info_sf$species <- 'H5N1'
h5n1_info_sf <- as.data.frame(h5n1_info_sf)


# Prepare and save simplified CSV for model training
h5n1_info_sf <- h5n1_info_sf |> 
  select(-c(State, County, Collection_date, geom, year, month)) |> 
  relocate(species)

write_csv(drop_na(h5n1_info_sf),
          'usa_h5n1_surveillance_env_vars_for_training.csv')



