# Load required libraries
library(rnaturalearth)
library(sf)
library(tidyverse)
library(jsonlite)

### 1. Europe Regional Partitioning ###

# Download and filter European countries
ne_countries <- ne_download(scale = 50, 
                            category = "cultural",
                            type = "countries",
                            returnclass = "sf")

europea <- ne_countries |> filter(CONTINENT == 'Europe')

# Define regional groupings
South_Western <- c('France', 'Portugal', 'Spain', 'Monaco', 'Italy', 
                   'Switzerland', 'Andorra', 'San Marino', 'Vatican')

North_Central <- c("Czechia", "Germany", "Belgium", "Ireland", 
                   "Luxembourg", "Netherlands", "United Kingdom", 
                   "Iceland", "Isle of Man", "Jersey", 
                   "Albania", "Austria", "Croatia", "Greece", 
                   "Hungary", "Kosovo", "Slovakia", "Slovenia", 
                   "Bosnia and Herzegovina", "Malta", "Montenegro", 
                   "Macedonia", "Serbia", "Bulgaria", "Moldova", 
                   "Romania", "Cyprus", 'North Macedonia', 
                   'Republic of Serbia', 'Liechtenstein')

North_Eastern <- c("Estonia", "Finland", "Latvia", "Lithuania", 
                   "Russia", "Belarus", "Norway", "Sweden", 
                   "Ukraine", "Poland", "Denmark")

# Assign regions to countries
eu_region <- europea |> 
  mutate(region = case_when(
    SOVEREIGNT %in% South_Western ~ 'South-Western',
    SOVEREIGNT %in% North_Central ~ 'North-Central',
    SOVEREIGNT %in% North_Eastern ~ 'North-Eastern',
    .default = SOVEREIGNT
  )) |> 
  relocate(region)

# Set factor levels for regions
eu_region$region <- factor(eu_region$region, 
                           levels = c('South-Western', 'North-Central', 'North-Eastern'))

## Calculate central coordinates for each region ##

# Disable S2 geometry processing for compatibility
sf_use_s2(FALSE)

# Define bounding box for Europe
bbox <- st_bbox(c(xmin = -25, xmax = 60, ymax = 75, ymin = 35), 
                crs = st_crs(eu_region))

# Crop and merge geometries by region
eu_region_cropped <- st_crop(eu_region, bbox)

eu_region_cropped_merged <- eu_region_cropped |> 
  group_by(region) |> 
  summarise(geometry = st_union(geometry))

# Save merged regions as GeoJSON
st_write(eu_region_cropped_merged, "./data/geographic/europe/europe_regions.geojson")

# Calculate centroids for each region
eu_region_cropped_central_coordinate <- eu_region_cropped_merged |> 
  st_centroid()

coordinates <- st_coordinates(eu_region_cropped_central_coordinate$geometry)

# Create and save coordinate data frame
coords_df <- data.frame(
  region = eu_region_cropped_central_coordinate$region,
  latitude = coordinates[, 2],
  longitude = coordinates[, 1]
)

write_tsv(coords_df, './data/geographic/europe/europe_region_centroids.tsv')

## Plot regional map ##
eu_region_cropped$region <- factor(eu_region_cropped$region, 
                                   levels = c('North-Eastern', 'North-Central', 'South-Western'))

ggplot(eu_region_cropped) +
  geom_sf(aes(fill = region)) +
  scale_fill_manual(values = c(
    'North-Eastern' = '#FED297', 
    'North-Central' = '#CBE9B1',
    'South-Western' = '#ABE0E4'
  )) +
  labs(fill = NULL) +
  theme_light() +
  theme(
    legend.key.spacing.y = unit(0.2, 'cm'),
    plot.margin = margin(t = 2, r = 0, b = 2, l = 0, unit = 'mm'),
    legend.text = element_text(size = 13),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) 

# Save plot
ggsave("./figures/europe_regions_map.pdf", width=6, height=4)

### 2. North America Regional Partitioning ###

# Get North American countries
na_countries <- rnaturalearth::ne_countries(
  scale = 50, 
  type = 'countries',
  continent = 'North America'
)

# Exclude USA and Canada (will process states separately)
na_countries_exclude_usa_can <- na_countries |> 
  filter(!sovereignt %in% c('United States of America', 'Canada')) |> 
  select(country = sovereignt)

# Get US and Canada states
usa_can_state <- ne_states(country = c('United States of America', 'Canada'))

usa_can_state <- usa_can_state |> 
  select(country = admin, state = name)

# Combine countries and states
na_sf <- bind_rows(na_countries_exclude_usa_can, usa_can_state)

## Define North American regions ##
na_sf_region <- na_sf |> 
  mutate(region = case_when(
    country == 'United States of America' & state == 'Alaska' ~ 'North',
    country == 'United States of America' & 
      state %in% c('California', 'Arizona', 'New Mexico', 'Texas', 
                   'Louisiana', 'Oklahoma', 'Arkansas', 'Mississippi', 
                   'Alabama', 'Tennessee', 'North Carolina',
                   'South Carolina', 'Georgia', 'Florida', 'Hawaii') ~ 'South',
    country == 'United States of America' ~ 'Central',
    country %in% c('Canada', 'Denmark') ~ 'North',
    .default = 'South'
  ))

# Define crop box for North America
crop_box <- st_as_sfc(st_bbox(c(xmin = -0.1, xmax = -180, 
                                ymin = -90, ymax = 90), 
                              crs = st_crs(na_sf_region)))

na_croped_df <- st_intersection(na_sf_region, crop_box)

# Set factor levels for regions
na_croped_df$region <- factor(na_croped_df$region, 
                              levels = c('North', 'Central', 'South'))

## Calculate central coordinates for each region ##

# Merge geometries by region
na_region_merged <- na_croped_df |> 
  group_by(region) |> 
  summarise(geometry = st_union(geometry))

# Save merged regions as GeoJSON
st_write(na_region_merged, "./data/geographic/north_america/north_america_regions.geojson")

# Calculate centroids for each region
na_region_merged_central_coordinate <- na_region_merged |> 
  st_centroid()

coordinates <- st_coordinates(na_region_merged_central_coordinate$geometry)

# Create and save coordinate data frame
coords_df <- data.frame(
  region = na_region_merged_central_coordinate$region,
  latitude = coordinates[, 2],
  longitude = coordinates[, 1]
)

write_tsv(coords_df, './data/geographic/north_america/north_america_region_centroids.tsv')

## Plot regional map ##
ggplot(na_croped_df) +
  geom_sf(aes(fill = region)) +
  scale_fill_manual(values = c(
    'North' = '#FED297', 
    'Central' = '#CBE9B1',
    'South' = '#ABE0E4'
  )) +
  labs(fill = NULL) +
  theme_light() +
  coord_sf(expand = FALSE) +
  theme(
    legend.key.spacing.y = unit(0.2, 'cm'),
    legend.text = element_text(size = 13),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    aspect.ratio = 0.7
  )

# Save plot
ggsave("./figures/north_america_regions_map.pdf", width=6, height=4)
