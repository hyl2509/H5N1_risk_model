# Processing Migratory Bird Data from eBird

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
taxonomy_ref <- read.csv('data/ebird-taxonomy.csv', header = TRUE)

# Select Anseriformes and Charadriiformes birds
target_orders <- taxonomy_ref |>
  select(species_code, order) |> 
  filter(order %in% c('Charadriiformes', 'Anseriformes')) |> 
  left_join(ebirdst_runs)

### 2. Identify high-quality migratory species and download abundance data
migratory_species <- target_orders |>
  filter(is_resident == FALSE) |> 
  filter(breeding_quality > 1 &
           nonbreeding_quality > 1 &
           postbreeding_migration_quality > 1 &
           prebreeding_migration_quality > 1)

species_list <- migratory_species$species_code

# Save filtered species information
write_tsv(migratory_species, 
          'data/Information_on_151_species_of_migratory_birds.tsv')

# Download abundance data with error handling
failed_downloads <- c()

for (species in species_list) {
  tryCatch({
    ebirdst_download_status(
      species = species, 
      pattern = "abundance_median_27km",
      path = "downloads/abundance_data"
    )
  }, error = function(e) {
    message(paste("Failed download:", species, "-", e$message))
    failed_downloads <<- c(failed_downloads, species)
  })
}

### 3. Calculate weekly total abundance across all species
abundance_files <- list.files(
  path = 'downloads/abundance_data',
  pattern = '_abundance_27km_2023\\.tif$',
  recursive = TRUE, 
  full.names = TRUE
)

weekly_abundance <- lapply(abundance_files, rast)

# [Rest of your processing code remains the same...]
# [Update output filenames accordingly below]

