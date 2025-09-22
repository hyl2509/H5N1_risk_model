# Processing Migratory Bird Data from eBird
# Description: This script processes abundance and migration data for Anseriformes and Charadriiformes birds from eBird.

# Load required packages
library(tidyverse)
library(ebirdst)
library(fields)
library(lubridate)
library(rnaturalearth)
library(sf)
library(terra)
library(readxl)

### 1. Get bird taxonomy list from eBird and filter target orders
ebird_taxonomy <- read.csv('data/ebird-taxonomy.csv', header = TRUE)

# Select Anseriformes and Charadriiformes birds
ebird_order <- ebird_taxonomy |>
  select(species_code, order) |> 
  filter(order == 'Charadriiformes' | order == 'Anseriformes') |> 
  left_join(ebirdst_runs)

### 2. Identify high-quality migratory species and download abundance data (27km resolution)

# Filter for migratory species with good data quality
migration_birds <- ebird_order |>
  filter(is_resident == FALSE) |> 
  filter(breeding_quality > 1 &
           nonbreeding_quality > 1 &
           postbreeding_migration_quality > 1 &
           prebreeding_migration_quality > 1)

species_target <- migration_birds$species_code

# Save species information
write_tsv(migration_birds, 
          'data/Information_on_151_species_of_migratory_birds.tsv')

# Download abundance data with error handling
error_species <- c()  # Vector to store failed downloads

for (i in species_target) {
  tryCatch({
    ebirdst_download_status(species = i, 
                          pattern = "abundance_median_27km",
                          path = "downloads/abundance_data")
  }, error = function(e) {
    message(paste("Error downloading species:", i, " - ", e$message))
    error_species <<- c(error_species, i)
  })
}

### 3. Calculate weekly total abundance across all species

# Load raster files
raster_file_median <- list.files(path = 'downloads/abundance_data',
                               pattern = '_abundance_median_27km_2023',
                               recursive = TRUE, full.names = TRUE)

bird_weekly_median <- lapply(raster_file_median, rast)

# Sum abundance across species for each week
for (i in names(bird_weekly_median[[1]])) {
  print(paste('Processing week:', i))
  empty_raster <- rast(ext = ext(bird_weekly_median[[1]]), 
                      resolution = res(bird_weekly_median[[1]]), 
                      crs = crs(bird_weekly_median[[1]]))
  values(empty_raster) <- 0
  names(empty_raster) <- i
  
  for (j in 1:length(bird_weekly_median)) {
    rast_j <- bird_weekly_median[[j]][i]
    empty_raster <- sum(empty_raster, rast_j, na.rm = TRUE)
  }
  
  writeRaster(empty_raster, 
             paste('globally_weekly_all_species_sum_of_abundance/', i, '.tif', sep = ''))
}

### 4. Calculate monthly abundance for each species (Jan-Dec)

# Create date reference dataframe for 2023
start_date <- as.Date("2023-01-01")
end_date <- as.Date("2023-12-31")
date_sequence <- seq(start_date, end_date, by = "day")
date_df <- data.frame(Date = date_sequence) |> 
  mutate(month = month(Date))

# Map dates to eBird weeks
bird_week_date <- as.Date(names(bird_weekly_median[[1]]))
names(bird_week_date) <- 1:52

date_df <- date_df %>%
  mutate(week = findInterval(Date, bird_week_date)) |> 
  mutate(week = if_else(week == 0, 52, week))

# English month names for file naming
months_vector <- c("January", "February", "March", "April", "May", "June", 
                  "July", "August", "September", "October", "November", "December")
names(months_vector) <- 1:12

# Function to calculate weighted monthly abundance
monthly_abundance <- function(month_num) {
  month_data <- date_df[date_df$month == month_num, ]
  count_month <- month_data |> 
    group_by(week) |> 
    count() |> 
    ungroup()
  
  for (j in 1:length(bird_weekly_median)) {
    j_species_data <- bird_weekly_median[[j]]
    values <- c()
    
    for (i in 1:nrow(count_month)) {
      row_i <- count_month[i, ]
      weekly_abundance <- j_species_data[[row_i$week]]
      values <- c(values, weekly_abundance)
    }
    
    values <- rast(values)
    weighted_mean_abundance <- terra::weighted.mean(x = values, w = count_month$n, na.rm = TRUE)
    
    # Generate output filename
    file_name <- sources(j_species_data) |> basename()
    species_name <- str_split(file_name, '_')[[1]][1]
    output_file <- paste(species_name, '_mean_abundance_in_', months_vector[month_num], sep = '')
    output_path <- paste("globally_every_bird_January_to_december_abundance/",
                        months_vector[month_num], '/', output_file, '.tif', sep = '')
    
    writeRaster(weighted_mean_abundance, filename = output_path, overwrite = TRUE)
  }
}

# Process all months
for (k in 1:12) {
  print(paste('Processing month:', k))
  monthly_abundance(k)
}

### 5. Calculate monthly total abundance across all species

for (i in 1:12) {
  bird_files <- list.files(
    path = paste0("globally_every_bird_January_to_december_abundance/",
                 months_vector[[i]]),
    full.names = TRUE)
  
  month_data <- rast(bird_files)
  month_sum_abundance <- sum(month_data, na.rm = TRUE)
  
  output_file <- paste("globally_monthly_all_species_sum_of_abundance/",
                      'Sum_of_bird_abundance_in_', months_vector[[i]], '.tif', sep = '')
  writeRaster(month_sum_abundance, filename = output_file, overwrite = TRUE)
}
