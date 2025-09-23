# ==============================================================================
# 
# Description: 
# This script demonstrates the complete workflow for modeling and predicting 
# monthly relative abundance of bird species across China. Using Common Merganser 
# as an example, it includes:
# 1. Data preparation and spatial subsampling
# 2. Encounter rate modeling with random forest and calibration
# 3. Expected count modeling 
# 4. Model evaluation and validation
# 5. Monthly relative abundance prediction across China
#
# Key steps:
# - Environmental variable extraction
# - Spatiotemporal stratified sampling
# - Machine learning model training (Random Forest)
# - Model calibration and threshold optimization
# - Partial dependence analysis
# - Monthly prediction raster generation
#

# ==============================================================================

# Load required libraries
library(tidyverse)
library(sf)
library(ebirdst)
library(caret)
library(ROCR)
library(scam)
library(mccf1)
library(gridExtra)
library(terra)
library(fields)
library(ranger)

# Set seed for reproducibility
set.seed(2024)

# Define target bird species
bird_name <- 'Common_Merganser'

# Set working directory
setwd(paste("42_birds_observation_data/", bird_name,
            '/generate_random_false_point_and_extract_all_point_env_envirable', sep = ''))

# Load China boundary map
china_map <- read_sf("map_of_China.json")

# Load environmental variables for observation data
env_variable <- read_csv('observation_data_env_variables.csv')
env_variable <- drop_na(env_variable)

# Split data into training (80%) and testing (20%) sets
env_variable$type <- if_else(runif(nrow(env_variable)) <= 0.80, "train", "test")

# ==============================================================================
# ENCOUNTER RATE MODEL TRAINING
# ==============================================================================

# 1. Spatiotemporal subsampling of checklists from January to December
env_variable_subsample <- grid_sample_stratified(env_variable,
                                                 obs_column = "species_observed",
                                                 sample_by = "type")

# Convert to spatial object for visualization
env_variable_subsample_sf <- env_variable_subsample |>
  select(species_observed, longitude, latitude) |>
  st_as_sf(coords = c('longitude', 'latitude'), crs = 4326)

# Plot sampling distribution
ggplot(china_map) +
  geom_sf() +
  geom_sf(data = env_variable_subsample_sf, size = 0.2, aes(color = species_observed))

# 2. Random Forest model training for encounter rate
checklists_train <- env_variable_subsample |>
  filter(type == 'train') |>
  # Select only columns used in the model
  select(species_observed,
         year, day_of_year, hours_of_day,
         effort_hours,
         starts_with("pland_"),
         starts_with("ed_"),
         starts_with("elevation_"),
         starts_with("rh_"),
         starts_with("ndvi_"),
         starts_with("temperature_"),
         starts_with("precipitation_"))

# Fix column names with special characters
checklists_train <- checklists_train |>
  rename(pland_c6_sonw_or_Ice = `pland_c6_sonw/Ice`,
         ed_c6_sonw_or_Ice = `ed_c6_sonw/Ice`)

# Calculate detection frequency
detection_freq <- mean(checklists_train$species_observed)

# Train Random Forest model (ranger requires factor response for classification)
er_model <- ranger(formula = as.factor(species_observed) ~ .,
                   data = checklists_train,
                   importance = "impurity",
                   probability = TRUE,
                   replace = TRUE,
                   sample.fraction = c(detection_freq, detection_freq))

# 3. Model calibration
# Predicted encounter rate based on out-of-bag samples
er_pred <- er_model$predictions[, 2]

# Observed detection (converted back from factor)
det_obs <- as.integer(checklists_train$species_observed)

# Construct data frame for calibration model training
obs_pred <- data.frame(obs = det_obs, pred = er_pred)

# Train calibration model using SCAM (Shape Constrained Additive Models)
calibration_model <- scam(obs ~ s(pred, k = 6, bs = "mpi"),
                          gamma = 2,
                          data = obs_pred)

# Group predicted encounter rate into bins of width 0.02
# Calculate mean observed encounter rates in each bin
er_breaks <- seq(0, 1, by = 0.02)
mean_er <- obs_pred |>
  mutate(er_bin = cut(pred, breaks = er_breaks, include.lowest = TRUE)) |>
  group_by(er_bin) |>
  summarise(n_checklists = n(),
            pred = mean(pred),
            obs = mean(obs),
            .groups = "drop")

# Make predictions from calibration model
calibration_curve <- data.frame(pred = er_breaks)
cal_pred <- predict(calibration_model, calibration_curve, type = "response")
calibration_curve$calibrated <- cal_pred

# Compare binned mean encounter rates to calibration model
ggplot(calibration_curve) +
  aes(x = pred, y = calibrated) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_line(color = "blue") +
  geom_point(data = mean_er,
             aes(x = pred, y = obs),
             size = 2, alpha = 0.6,
             show.legend = FALSE) +
  labs(x = "Estimated encounter rate",
       y = "Observed encounter rate",
       title = "Calibration model") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  theme(plot.title = element_text(hjust = 0.5))

# Save calibration plot
ggsave(path = paste("42_birds_observation_data/", bird_name,
                    '/prediction_resulits_of_relative_abundance/model_estimation_results', sep = ''),
       filename = 'Calibration_model.tif',
       height = 4, width = 7)

# 4. Threshold optimization using MCC and F1 scores
mcc_f1 <- mccf1(
  # Observed detection/non-detection
  response = obs_pred$obs,
  # Predicted encounter rate from random forest
  predictor = obs_pred$pred)

# Identify best threshold
mcc_f1_summary <- summary(mcc_f1)
threshold <- mcc_f1_summary$best_threshold[1]

# 5. Model evaluation on test set
checklists_test <- filter(env_variable_subsample, type == "test") |>
  select(species_observed,
         year, day_of_year, hours_of_day,
         effort_hours,
         starts_with("pland_"),
         starts_with("ed_"),
         starts_with("elevation_"),
         starts_with("rh_"),
         starts_with("ndvi_"),
         starts_with("temperature_"),
         starts_with("precipitation_")) |>
  rename(pland_c6_sonw_or_Ice = `pland_c6_sonw/Ice`,
         ed_c6_sonw_or_Ice = `ed_c6_sonw/Ice`) |>
  mutate(species_observed = as.integer(species_observed))

# Predict on test data using random forest model
pred_er <- predict(er_model, data = checklists_test, type = "response")
pred_er <- pred_er$predictions[, 2]  # Extract probability of detection

# Convert predictions to binary (presence/absence) using threshold
pred_binary <- as.integer(pred_er > threshold)

# Apply calibration
pred_calibrated <- predict(calibration_model,
                           newdata = data.frame(pred = pred_er),
                           type = "response") |>
  as.numeric()

# Constrain probabilities to 0-1 range
pred_calibrated[pred_calibrated < 0] <- 0
pred_calibrated[pred_calibrated > 1] <- 1

# Combine observations and predictions
obs_pred_test <- data.frame(id = seq_along(pred_calibrated),
                            # Actual detection/non-detection
                            obs = as.integer(checklists_test$species_observed),
                            # Binary detection/non-detection prediction
                            pred_binary = pred_binary,
                            # Calibrated encounter rate
                            pred_calibrated = pred_calibrated)

# Calculate evaluation metrics
# Mean squared error (MSE)
mse <- mean((obs_pred_test$obs - obs_pred_test$pred_calibrated)^2, na.rm = TRUE)

# Precision-recall AUC
em <- precrec::evalmod(scores = obs_pred_test$pred_binary,
                       labels = obs_pred_test$obs)
pr_auc <- precrec::auc(em) |>
  filter(curvetypes == "PRC") |>
  pull(aucs)

# Calculate binary prediction metrics: sensitivity, specificity
pa_metrics <- obs_pred_test |>
  select(id, obs, pred_binary) |>
  PresenceAbsence::presence.absence.accuracy(na.rm = TRUE, st.dev = FALSE)

# MCC and F1 scores
mcc_f1 <- calculate_mcc_f1(obs_pred_test$obs, obs_pred_test$pred_binary)

# Combine all performance metrics
ppms <- data.frame(
  mse = mse,
  sensitivity = pa_metrics$sensitivity,
  specificity = pa_metrics$specificity,
  pr_auc = pr_auc,
  mcc = mcc_f1$mcc,
  f1 = mcc_f1$f1
)

# Display metrics
knitr::kable(pivot_longer(ppms, everything()), digits = 3)

# Save evaluation metrics
ppms_longer <- ppms |>
  pivot_longer(cols = 1:6,
               names_to = 'name',
               values_to = 'value')

write_csv(ppms_longer,
          paste("42_birds_observation_data/", bird_name,
                '/prediction_resulits_of_relative_abundance/model_estimation_results/encounter_rate_model_evaluation_indicators.csv', sep = ''))

# 6. Predictor importance analysis
# Extract predictor importance from random forest model
pred_imp <- er_model$variable.importance

pred_imp <- data.frame(predictor = names(pred_imp),
                       importance = pred_imp) |>
  arrange(desc(importance))

# Plot importance of top 35 predictors
ggplot(head(pred_imp, 35)) +
  aes(x = reorder(predictor, importance), y = importance) +
  geom_col() +
  geom_hline(yintercept = 0, linewidth = 2, colour = "#555555") +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  labs(x = NULL,
       y = "Predictor Importance (Gini Index)",
       title = 'Top 35 predictors of importance') +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(colour = "#cccccc", linewidth = 0.5),
        plot.title = element_text(hjust = 0.5))

# Save importance plot
ggsave(path = paste("42_birds_observation_data/", bird_name,
                    '/prediction_resulits_of_relative_abundance/model_estimation_results', sep = ''),
       filename = 'encounter_rate_model_predictors_of_importance.tif',
       height = 5, width = 8)

# 7. Partial dependence analysis
# Function to calculate partial dependence for a single predictor
calculate_pd <- function(predictor, er_model, calibration_model,
                         data, x_res = 25, n = 1000) {
  # Create prediction grid using quantiles
  x_grid <- quantile(data[[predictor]],
                     probs = seq(from = 0, to = 1, length = x_res),
                     na.rm = TRUE)
  # Remove duplicates
  x_grid <- x_grid[!duplicated(signif(x_grid, 8))]
  x_grid <- unname(unique(x_grid))
  grid <- data.frame(predictor = predictor, x = x_grid)
  names(grid) <- c("predictor", predictor)
  
  # Subsample training data
  n <- min(n, nrow(data))
  data <- data[sample(seq.int(nrow(data)), size = n, replace = FALSE), ]
  
  # Drop focal predictor from data
  data <- data[names(data) != predictor]
  grid <- merge(grid, data, all = TRUE)
  
  # Predict encounter rate
  p <- predict(er_model, data = grid)
  
  # Summarize results
  pd <- grid[, c("predictor", predictor)]
  names(pd) <- c("predictor", "x")
  pd$encounter_rate <- p$predictions[, 2]
  pd <- dplyr::group_by(pd, predictor, x)
  pd <- dplyr::summarise(pd,
                         encounter_rate = mean(encounter_rate, na.rm = TRUE),
                         .groups = "drop")
  
  # Apply calibration
  pd$encounter_rate <- predict(calibration_model,
                               newdata = data.frame(pred = pd$encounter_rate),
                               type = "response")
  pd$encounter_rate <- as.numeric(pd$encounter_rate)
  
  # Constrain to 0-1 range
  pd$encounter_rate[pd$encounter_rate < 0] <- 0
  pd$encounter_rate[pd$encounter_rate > 1] <- 1
  
  return(pd)
}

# Calculate partial dependence for each of the top 6 predictors
pd <- NULL
for (predictor in head(pred_imp$predictor)) {
  pd <- calculate_pd(predictor,
                     er_model = er_model,
                     calibration_model = calibration_model,
                     data = checklists_train) |>
    bind_rows(pd)
}

# Plot partial dependence
ggplot(pd) +
  aes(x = x, y = encounter_rate) +
  geom_line() +
  geom_point() +
  facet_wrap(~ factor(predictor, levels = rev(unique(predictor))),
             ncol = 2, scales = "free") +
  labs(x = NULL, y = "Encounter Rate",
       title = 'Partial dependence for each of the top 6 predictors') +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "grey60"),
        axis.ticks  = element_line(color = "grey60"),
        plot.title = element_text(hjust = 0.5))

# Save partial dependence plot
ggsave(path = paste("42_birds_observation_data/", bird_name,
                    '/prediction_resulits_of_relative_abundance/model_estimation_results', sep = ''),
       filename = 'encounter_rate_model_partial_dependence_for_each_of_the_top_6_predictors.tif',
       height = 5, width = 8)

# 8. Standardize effort variables
# Find peak time of day from partial dependence
pd_time <- calculate_pd("hours_of_day",
                        er_model = er_model,
                        calibration_model = calibration_model,
                        data = checklists_train) |>
  select(hours_of_day = x, encounter_rate)

# Create histogram of observation times
g_hist <- ggplot(checklists_train) +
  aes(x = hours_of_day) +
  geom_histogram(binwidth = 1, center = 0.5, color = "grey30",
                 fill = "grey50") +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "Hours since midnight",
       y = "# checklists",
       title = "Distribution of observation start times")

# Create partial dependence plot for time
g_pd <- ggplot(pd_time) +
  aes(x = hours_of_day, y = encounter_rate) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  labs(x = "Hours since midnight",
       y = "Encounter rate",
       title = "Observation start time partial dependence")

# Combine plots
grid.arrange(g_hist, g_pd)

# Trim ends of partial dependence curve
pd_time_trimmed <- pd_time[c(-1, -nrow(pd_time)), ]

# Plot trimmed partial dependence
ggplot(pd_time_trimmed) +
  aes(x = hours_of_day, y = encounter_rate) +
  geom_line() +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 24, by = 3)) +
  labs(x = "Hours since midnight",
       y = "Occurrence probability") +
  theme(panel.grid.minor = element_blank())

# Save time dependence plot
ggsave(filename = 'The probability varies with the start time of observation.pdf',
       path = paste("42_birds_observation_data/", bird_name,
                    '/prediction_resulits_of_relative_abundance/model_estimation_results', sep = ''),
       width = 4, height = 3)

# Identify time maximizing encounter rate
pd_time_trimmed <- arrange(pd_time_trimmed, desc(encounter_rate))
t_peak <- pd_time_trimmed$hours_of_day[1]
print(t_peak)

# Find peak effort hours from partial dependence
pd_effort_hours <- calculate_pd("effort_hours",
                                er_model = er_model,
                                calibration_model = calibration_model,
                                data = checklists_train) |>
  select(effort_hours = x, encounter_rate)

# Create histogram of effort hours
g_effort_hist <- ggplot(checklists_train) +
  aes(x = effort_hours) +
  geom_histogram(binwidth = 0.5, boundary = 0, color = "grey30",
                 fill = "grey50") +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "effort hours",
       y = "# checklists",
       title = "Distribution of effort hours")

# Trim ends of partial dependence curve for effort
pd_effort_trimmed <- pd_effort_hours[c(-1, -nrow(pd_effort_hours)), ]

# Plot effort dependence
ggplot(pd_effort_trimmed) +
  aes(x = effort_hours, y = encounter_rate) +
  geom_line() +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 10, by = 1)) +
  labs(x = "Observation duration",
       y = "Occurrence probability") +
  theme(panel.grid.minor = element_blank())

# Save effort dependence plot
ggsave(filename = 'The probability varies with the duration of observation.pdf',
       path = paste("42_birds_observation_data/", bird_name,
                    '/prediction_resulits_of_relative_abundance/model_estimation_results', sep = ''),
       width = 4, height = 3)

# Identify effort duration maximizing encounter rate
pd_effort_trimmed <- arrange(pd_effort_trimmed, desc(encounter_rate))
effort_peak <- pd_effort_trimmed$effort_hours[1]
print(effort_peak)

# ==============================================================================
# EXPECTED COUNT MODEL TRAINING
# ==============================================================================

# 1. Train expected count model
train_count <- env_variable_subsample |>
  filter(type == 'train') |>
  # Select columns for count model
  select(species_observed,
         observation_count,
         year, day_of_year, hours_of_day,
         effort_hours,
         starts_with("pland_"),
         starts_with("ed_"),
         starts_with("elevation_"),
         starts_with("rh_"),
         starts_with("ndvi_"),
         starts_with("temperature_"),
         starts_with("precipitation_")) |>
  rename(pland_c6_sonw_or_Ice = `pland_c6_sonw/Ice`,
         ed_c6_sonw_or_Ice = `ed_c6_sonw/Ice`)

# Add predicted encounter rate
train_count$pred_er <- er_model$predictions[, 2]

# Subset to only observed or predicted detections
train_count <- train_count |>
  filter(!is.na(observation_count),
         observation_count > 0 | pred_er > threshold) |>
  select(-species_observed, -pred_er)

# Predict encounter rate for count model training
predicted_er <- predict(er_model, data = train_count, type = "response")
predicted_er <- predicted_er$predictions[, 2]
train_count$predicted_er <- predicted_er

# Log-transform count data
train_count$observation_count <- log10(train_count$observation_count + 1)

# Train count model
count_model <- ranger(formula = observation_count ~ .,
                      data = train_count,
                      importance = "impurity",
                      replace = TRUE)

# 2. Model evaluation
checklists_test <- filter(env_variable_subsample, type == "test") |>
  select(species_observed,
         observation_count,
         year, day_of_year, hours_of_day,
         effort_hours,
         starts_with("pland_"),
         starts_with("ed_"),
         starts_with("elevation_"),
         starts_with("rh_"),
         starts_with("ndvi_"),
         starts_with("temperature_"),
         starts_with("precipitation_")) |>
  rename(pland_c6_sonw_or_Ice = `pland_c6_sonw/Ice`,
         ed_c6_sonw_or_Ice = `ed_c6_sonw/Ice`) |>
  mutate(species_observed = as.integer(species_observed)) |>
  filter(!is.na(observation_count))

# Estimate encounter rate for test data
pred_er <- predict(er_model, data = checklists_test, type = "response")
pred_er <- pred_er$predictions[, 2]

# Convert to binary using threshold
pred_binary <- as.integer(pred_er > threshold)

# Apply calibration
pred_calibrated <- predict(calibration_model,
                           newdata = data.frame(pred = pred_er),
                           type = "response") |>
  as.numeric()

# Constrain probabilities to 0-1 range
pred_calibrated[pred_calibrated < 0] <- 0
pred_calibrated[pred_calibrated > 1] <- 1

# Add predicted encounter rate for count estimates
checklists_test$predicted_er <- pred_er

# Estimate count
pred_count <- predict(count_model, data = checklists_test, type = "response")
pred_count <- pred_count$predictions
pred_count <- (10^pred_count - 1)  # Back-transform from log scale

# Calculate relative abundance (product of encounter rate and count)
pred_abundance <- pred_calibrated * pred_count

# Combine observations and predictions
obs_pred_test <- data.frame(
  id = seq_along(pred_abundance),
  # Actual detection/non-detection
  obs_detected = as.integer(checklists_test$species_observed),
  obs_count = checklists_test$observation_count,
  # Model estimates
  pred_binary = pred_binary,
  pred_er = pred_calibrated,
  pred_count = pred_count,
  pred_abundance = pred_abundance)

# Subset to only checklists where detection occurred
detections_test <- filter(obs_pred_test, obs_detected > 0)

# Calculate count metrics (only on checklists with detections)
count_spearman <- cor(detections_test$pred_count,
                      detections_test$obs_count,
                      method = "spearman")

log_count_pearson <- cor(log(detections_test$pred_count + 1),
                         log(detections_test$obs_count + 1),
                         method = "pearson")

# Calculate abundance metrics
abundance_spearman <- cor(detections_test$pred_abundance,
                          detections_test$obs_count,
                          method = "spearman")

log_abundance_pearson <- cor(log(detections_test$pred_abundance + 1),
                             log(detections_test$obs_count + 1),
                             method = "pearson")

# Combine performance metrics
ppms <- data.frame(
  count_spearman = count_spearman,
  log_count_pearson = log_count_pearson,
  abundance_spearman = abundance_spearman,
  log_abundance_pearson = log_abundance_pearson)

# Display metrics
knitr::kable(pivot_longer(ppms, everything()), digits = 3)

# Save evaluation metrics
ppms_longer <- ppms |>
  pivot_longer(cols = 1:4,
               names_to = 'name',
               values_to = 'value')

write_csv(ppms_longer,
          paste("42_birds_observation_data/", bird_name,
                '/prediction_resulits_of_relative_abundance/model_estimation_results/expected_observation_number_model_evaluation_indicators.csv', sep = ''))

# ==============================================================================
# MONTHLY RELATIVE ABUNDANCE PREDICTION ACROSS CHINA
# ==============================================================================

# Create vector of English month names
months_vector <- c("January", "February", "March", "April", "May", "June",
                   "July", "August", "September", "October", "November", "December")

# Assign names to each element
names(months_vector) <- 1:12

# Loop through each month for prediction
for (i in 1:12) {
  # Load prediction grid environmental variables for current month
  pred_grid_file <- paste('extract_prediction_grid_environmental_variables/',
                          months_vector[[i]], '/result_of_prediction_grid_env_variable/', 
                          tolower(months_vector[[i]]), '_prediction_grid_env.csv', sep = '')
  pred_grid <- read_csv(pred_grid_file)
  
  # Prepare prediction grid with standardized effort variables
  pred_grid_eff <- pred_grid |>
    mutate(observation_date = ymd(paste('2023-', i, '-15', sep = '')),
           year = year(observation_date),
           day_of_year = yday(observation_date),
           # Use optimal time for detection identified earlier
           hours_of_day = t_peak,
           effort_hours = effort_peak) |>
    rename(pland_c6_sonw_or_Ice = `pland_c6_sonw/Ice`,
           ed_c6_sonw_or_Ice = `ed_c6_sonw/Ice`)
  
  # Estimate encounter rate
  pred_er <- predict(er_model, data = pred_grid_eff, type = "response")
  pred_er <- pred_er$predictions[, 2]
  
  # Define range boundary using threshold
  pred_binary <- as.integer(pred_er > threshold)
  
  # Apply calibration
  pred_er_cal <- predict(calibration_model,
                         data.frame(pred = pred_er),
                         type = "response") |>
    as.numeric()
  
  # Constrain to 0-1 range
  pred_er_cal[pred_er_cal < 0] <- 0
  pred_er_cal[pred_er_cal > 1] <- 1
  
  # Add predicted encounter rate for count estimates
  pred_grid_eff$predicted_er <- pred_er
  
  # Estimate count
  pred_count <- predict(count_model, data = pred_grid_eff, type = "response")
  pred_count <- pred_count$predictions
  pred_count <- (10^pred_count - 1)  # Back-transform from log scale
  
  # Combine predictions with coordinates
  predictions <- data.frame(cell_id = pred_grid_eff$cell_id,
                            x = pred_grid_eff$x,
                            y = pred_grid_eff$y,
                            in_range = pred_binary,
                            encounter_rate = pred_er_cal,
                            count = pred_count)
  
  # Calculate relative abundance (encounter rate × count)
  predictions$abundance <- predictions$encounter_rate * predictions$count
  
  # Rasterize predictions
  layers <- c("in_range", "encounter_rate", "count", "abundance")
  r <- rast("5km_prediction_grid.tif")
  crs <- st_crs(r)
  
  r_pred <- predictions |>
    # Convert to spatial features
    st_as_sf(coords = c("x", "y"), crs = crs) |>
    select(all_of(layers)) |>
    # Rasterize
    rasterize(r, field = layers)
  
  # Save prediction raster
  output_file <- paste("42_birds_observation_data/", bird_name,
                       '/prediction_resulits_of_relative_abundance/prediction_relative_abundance_raster_file/',
                       months_vector[[i]], '_relative_abundance.tif',
                       sep = '')
  
  terra::writeRaster(r_pred, filename = output_file, overwrite = TRUE)
  
  # Create visualization
  study_region <- st_transform(china_map, crs = crs)
  
  # Define quantile breaks, excluding zeros
  brks <- ifel(r_pred[["abundance"]] > 0, r_pred[["abundance"]], NA) |>
    global(fun = quantile,
           probs = seq(0, 1, 0.1), na.rm = TRUE) |>
    as.numeric() |>
    unique()
  
  # Label bottom, middle, and top values
  lbls <- round(c(min(brks), median(brks), max(brks)), 2)
  
  # Use eBird Status and Trends color palette
  pal <- ebirdst_palettes(length(brks) - 1)
  
  # Set output file path for plot
  output_tiff <- paste("42_birds_observation_data/", bird_name,
                       '/prediction_resulits_of_relative_abundance/prediction_relative_abundance_plot/',
                       months_vector[[i]], '_relative_abundance.tif',
                       sep = '')
  
  # Start TIFF device for high-resolution plot
  tiff(filename = output_tiff, width = 600, height = 550, units = 'mm', res = 300)
  
  # Create abundance map
  plot(r_pred[["abundance"]], col = c("#e6e6e6", pal),
       breaks = c(0, brks), legend = FALSE,
       axes = FALSE, bty = "n")
  
  # Add China boundary
  plot(study_region, border = "#000000", col = NA, lwd = 1, add = TRUE)
  
  # Add legend
  par(new = TRUE, mar = c(0, 0, 0, 0))
  title <- paste(str_replace_all(bird_name, '_', ' '), ' Relative Abundance (', months_vector[[i]], ' 2023)', sep = '')
  
  image.plot(zlim = c(0, 1), legend.only = TRUE,
             col = pal, breaks = seq(0, 1, length.out = length(brks)),
             smallplot = c(0.25, 0.75, 0.05, 0.1),  # Adjust legend position
             horizontal = TRUE,
             axis.args = list(at = c(0, 0.5, 1), labels = lbls,
                              fg = "black", col.axis = "black",
                              cex.axis = 3, lwd.ticks = 0.5,
                              padj = 0.1),  # Move text upward
             legend.args = list(text = title,
                                side = 3, col = "black",
                                cex = 3, line = -0.1, padj = -0.8))  # Adjust text position
  
  # Close TIFF device
  dev.off()
}