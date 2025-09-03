# ==========================================================
# Gaussian Process Gravity Model Analysis: Dubai Training & NYC Prediction
# ==========================================================

# Load required libraries
library(tidyverse)
library(kernlab)  # For Gaussian Process regression
library(scales)

# Constants
epsilon <- 1e-6  # small constant to avoid log(0)

# ==========================================================
# UTILITY FUNCTIONS
# ==========================================================

check_required_cols <- function(df, required_cols, df_name = "dataframe") {
  missing_cols <- setdiff(required_cols, colnames(df))
  if(length(missing_cols) > 0) {
    stop(paste0("Missing columns in ", df_name, ": ", paste(missing_cols, collapse = ", ")))
  }
}

prepare_flow_data <- function(flow_df, census_df, scale_vars = TRUE) {
  # Merge flow data with census variables and create log-transformed predictors
  
  # Check required columns
  check_required_cols(flow_df, c("cost", "origin", "destination", "flow"), "flow_df")
  check_required_cols(census_df, c("h3_index", "num_vendors", "pop_share", "mean_rent"), "census_df")
  
  # Merge census data for origins and destinations
  merged_df <- flow_df %>%
    left_join(census_df %>% 
                select(h3_index, num_vendors) %>% 
                rename(origin = h3_index, num_vendors_origin = num_vendors),
              by = "origin") %>%
    left_join(census_df %>% 
                select(h3_index, pop_share, mean_rent) %>% 
                rename(destination = h3_index, 
                       pop_share_dest = pop_share,
                       mean_rent_dest = mean_rent),
              by = "destination") %>%
    # Create log-transformed variables
    mutate(
      log_flow = log(flow + 1),
      log_num_vendors_origin = log(num_vendors_origin + epsilon),
      log_pop_share_dest = log(pop_share_dest),
      log_mean_rent_dest = log(mean_rent_dest + epsilon)
    ) %>%
    filter(!is.na(log_flow), !is.na(cost), 
           !is.na(log_num_vendors_origin), 
           !is.na(log_pop_share_dest), 
           !is.na(log_mean_rent_dest))
  
  # Optionally scale variables for GP
  if(scale_vars) {
    merged_df <- merged_df %>%
      mutate(
        cost_scaled = as.numeric(scale(cost)),
        log_num_vendors_origin_scaled = as.numeric(scale(log_num_vendors_origin)),
        log_pop_share_dest_scaled = as.numeric(scale(log_pop_share_dest)),
        log_mean_rent_dest_scaled = as.numeric(scale(log_mean_rent_dest))
      )
    
    # Store scaling parameters
    scaling_params <- list(
      cost_mean = mean(merged_df$cost, na.rm = TRUE),
      cost_sd = sd(merged_df$cost, na.rm = TRUE),
      log_num_vendors_mean = mean(merged_df$log_num_vendors_origin, na.rm = TRUE),
      log_num_vendors_sd = sd(merged_df$log_num_vendors_origin, na.rm = TRUE),
      log_pop_share_mean = mean(merged_df$log_pop_share_dest, na.rm = TRUE),
      log_pop_share_sd = sd(merged_df$log_pop_share_dest, na.rm = TRUE),
      log_mean_rent_mean = mean(merged_df$log_mean_rent_dest, na.rm = TRUE),
      log_mean_rent_sd = sd(merged_df$log_mean_rent_dest, na.rm = TRUE)
    )
    
    attr(merged_df, "scaling_params") <- scaling_params
  }
  
  return(merged_df)
}

# ==========================================================
# GAUSSIAN PROCESS MODEL FUNCTIONS
# ==========================================================

fit_gaussian_process <- function(flow_df, census_df, dist_max = 30, kernel_type = "rbfdot", 
                                 sample_size = 5000, use_scaled = TRUE) {
  # Fit Gaussian Process regression model for flow prediction
  
  cat("Preparing data for GP model...\n")
  
  # Prepare data
  df <- prepare_flow_data(flow_df, census_df, scale_vars = TRUE) %>%
    filter(cost <= dist_max)
  
  # Sample data if too large for GP
  if(nrow(df) > sample_size) {
    cat("Sampling", sample_size, "observations from", nrow(df), "total observations\n")
    set.seed(123)
    df <- df %>% sample_n(sample_size)
  }
  
  # Prepare feature matrix
  if(use_scaled) {
    X <- df %>% 
      select(cost_scaled, log_num_vendors_origin_scaled, 
             log_pop_share_dest_scaled, log_mean_rent_dest_scaled) %>%
      as.matrix()
    feature_names <- c("cost_scaled", "log_num_vendors_origin_scaled", 
                       "log_pop_share_dest_scaled", "log_mean_rent_dest_scaled")
  } else {
    X <- df %>% 
      select(cost, log_num_vendors_origin, 
             log_pop_share_dest, log_mean_rent_dest) %>%
      as.matrix()
    feature_names <- c("cost", "log_num_vendors_origin", 
                       "log_pop_share_dest", "log_mean_rent_dest")
  }
  
  y <- df$log_flow
  
  cat("Feature matrix dimensions:", dim(X), "\n")
  cat("Target variable range:", range(y, na.rm = TRUE), "\n")
  
  # Remove any remaining NA values
  complete_cases <- complete.cases(X, y)
  X <- X[complete_cases, , drop = FALSE]
  y <- y[complete_cases]
  
  cat("Complete cases:", length(y), "\n")
  
  # Define kernel
  if(kernel_type == "rbfdot") {
    kernel <- rbfdot(sigma = 0.1)  # RBF kernel
  } else if(kernel_type == "polydot") {
    kernel <- polydot(degree = 2)  # Polynomial kernel
  } else {
    kernel <- rbfdot(sigma = 0.1)  # Default to RBF
  }
  
  cat("Fitting GP model with", kernel_type, "kernel...\n")
  
  # Fit Gaussian Process
  gp_model <- gausspr(X, y, kernel = kernel, variance.model = TRUE)
  
  # Calculate model performance on training data
  y_pred <- predict(gp_model, X)
  mse <- mean((y - y_pred)^2, na.rm = TRUE)
  rmse <- sqrt(mse)
  r_squared <- 1 - sum((y - y_pred)^2, na.rm = TRUE) / sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  
  cat("Training performance:\n")
  cat("  RMSE:", round(rmse, 4), "\n")
  cat("  R²:", round(r_squared, 4), "\n")
  
  # Return model and metadata
  return(list(
    model = gp_model,
    training_data = df,
    feature_names = feature_names,
    use_scaled = use_scaled,
    scaling_params = attr(df, "scaling_params"),
    performance = list(rmse = rmse, r_squared = r_squared),
    kernel_type = kernel_type
  ))
}

predict_gp_flows <- function(gp_result, new_data, return_variance = TRUE) {
  # Make predictions using fitted GP model
  
  # Prepare feature matrix for new data
  if(gp_result$use_scaled) {
    # Apply same scaling as training data
    scaling <- gp_result$scaling_params
    X_new <- new_data %>%
      mutate(
        cost_scaled = (cost - scaling$cost_mean) / scaling$cost_sd,
        log_num_vendors_origin_scaled = (log_num_vendors_origin - scaling$log_num_vendors_mean) / scaling$log_num_vendors_sd,
        log_pop_share_dest_scaled = (log_pop_share_dest - scaling$log_pop_share_mean) / scaling$log_pop_share_sd,
        log_mean_rent_dest_scaled = (log_mean_rent_dest - scaling$log_mean_rent_mean) / scaling$log_mean_rent_sd
      ) %>%
      select(all_of(gp_result$feature_names)) %>%
      as.matrix()
  } else {
    X_new <- new_data %>%
      select(all_of(gp_result$feature_names)) %>%
      as.matrix()
  }
  
  # Make predictions
  if(return_variance) {
    pred_result <- predict(gp_result$model, X_new, type = "response")
    # Note: kernlab doesn't directly return variance, but we can estimate uncertainty
    pred_mean <- as.numeric(pred_result)
    
    # Simple uncertainty estimate based on distance from training data
    training_X <- gp_result$model@xmatrix
    min_dist <- apply(X_new, 1, function(x) {
      min(sqrt(rowSums((training_X - matrix(x, nrow = nrow(training_X), ncol = length(x), byrow = TRUE))^2)))
    })
    
    # Heuristic: higher uncertainty for points far from training data
    pred_var <- pmax(0.1, min_dist * 0.5)
    
    return(list(
      mean = pred_mean,
      variance = pred_var,
      log_flow_pred = pred_mean,
      flow_pred = pmax(exp(pred_mean) - 1, 0)
    ))
  } else {
    pred_mean <- as.numeric(predict(gp_result$model, X_new))
    return(list(
      mean = pred_mean,
      log_flow_pred = pred_mean,
      flow_pred = pmax(exp(pred_mean) - 1, 0)
    ))
  }
}

# ==========================================================
# VISUALIZATION FUNCTIONS
# ==========================================================

plot_gp_distance_decay <- function(gp_result, n_points = 200, dist_range = c(0, 30)) {
  # Plot distance decay curve from GP model
  
  # Create prediction grid varying only distance
  training_data <- gp_result$training_data
  
  # Use median values for other variables
  median_vendors <- median(training_data$log_num_vendors_origin, na.rm = TRUE)
  median_pop <- median(training_data$log_pop_share_dest, na.rm = TRUE)
  median_rent <- median(training_data$log_mean_rent_dest, na.rm = TRUE)
  
  pred_grid <- data.frame(
    cost = seq(dist_range[1], dist_range[2], length.out = n_points),
    log_num_vendors_origin = median_vendors,
    log_pop_share_dest = median_pop,
    log_mean_rent_dest = median_rent
  )
  
  # Make predictions
  preds <- predict_gp_flows(gp_result, pred_grid, return_variance = TRUE)
  
  pred_grid$pred_log_flow <- preds$mean
  pred_grid$pred_flow <- preds$flow_pred
  pred_grid$pred_var <- preds$variance
  pred_grid$pred_lower <- preds$flow_pred - 1.96 * sqrt(preds$variance)
  pred_grid$pred_upper <- preds$flow_pred + 1.96 * sqrt(preds$variance)
  
  # Create plot
  p <- ggplot(pred_grid, aes(x = cost)) +
    geom_ribbon(aes(ymin = pmax(pred_lower, 0), ymax = pred_upper), 
                fill = "lightblue", alpha = 0.3) +
    geom_line(aes(y = pred_flow), color = "red", size = 1.2) +
    labs(title = "GP Model: Predicted Flow vs Distance",
         subtitle = "Other variables held at median values",
         x = "Distance (cost)", 
         y = "Predicted Flow") +
    theme_minimal()
  
  return(list(plot = p, predictions = pred_grid))
}

plot_gp_feature_effects <- function(gp_result, feature_name, n_points = 100) {
  # Plot individual feature effects from GP model
  
  training_data <- gp_result$training_data
  
  # Get feature range
  if(feature_name == "cost") {
    feature_range <- range(training_data$cost, na.rm = TRUE)
    feature_values <- seq(feature_range[1], feature_range[2], length.out = n_points)
  } else if(feature_name == "log_num_vendors_origin") {
    feature_range <- range(training_data$log_num_vendors_origin, na.rm = TRUE)
    feature_values <- seq(feature_range[1], feature_range[2], length.out = n_points)
  } else if(feature_name == "log_pop_share_dest") {
    feature_range <- range(training_data$log_pop_share_dest, na.rm = TRUE)
    feature_values <- seq(feature_range[1], feature_range[2], length.out = n_points)
  } else if(feature_name == "log_mean_rent_dest") {
    feature_range <- range(training_data$log_mean_rent_dest, na.rm = TRUE)
    feature_values <- seq(feature_range[1], feature_range[2], length.out = n_points)
  } else {
    stop("Invalid feature name")
  }
  
  # Create prediction grid with other variables at median
  pred_grid <- data.frame(
    cost = rep(median(training_data$cost, na.rm = TRUE), n_points),
    log_num_vendors_origin = rep(median(training_data$log_num_vendors_origin, na.rm = TRUE), n_points),
    log_pop_share_dest = rep(median(training_data$log_pop_share_dest, na.rm = TRUE), n_points),
    log_mean_rent_dest = rep(median(training_data$log_mean_rent_dest, na.rm = TRUE), n_points)
  )
  
  # Vary the target feature
  pred_grid[[feature_name]] <- feature_values
  
  # Make predictions
  preds <- predict_gp_flows(gp_result, pred_grid, return_variance = TRUE)
  
  pred_grid$pred_flow <- preds$flow_pred
  pred_grid$pred_var <- preds$variance
  pred_grid$pred_lower <- preds$flow_pred - 1.96 * sqrt(preds$variance)
  pred_grid$pred_upper <- preds$flow_pred + 1.96 * sqrt(preds$variance)
  
  # Create plot
  p <- ggplot(pred_grid, aes_string(x = feature_name)) +
    geom_ribbon(aes(ymin = pmax(pred_lower, 0), ymax = pred_upper), 
                fill = "lightblue", alpha = 0.3) +
    geom_line(aes(y = pred_flow), color = "red", size = 1.2) +
    labs(title = paste("GP Model: Effect of", feature_name),
         subtitle = "Other variables held at median values",
         x = feature_name, 
         y = "Predicted Flow") +
    theme_minimal()
  
  return(p)
}

# ==========================================================
# MAIN ANALYSIS WORKFLOW
# ==========================================================

run_dubai_gp_analysis <- function(flow_filepath, census_filepath, city_name = "Dubai", 
                                  dist_max = 30, kernel_type = "rbfdot", sample_size = 5000) {
  # Run complete GP analysis for Dubai data
  
  cat("=== Loading Dubai Data ===\n")
  flow_df <- read_csv(path.expand(flow_filepath), show_col_types = FALSE)
  census_df <- read_csv(path.expand(census_filepath), show_col_types = FALSE)
  
  cat("Flow data:", nrow(flow_df), "rows\n")
  cat("Census data:", nrow(census_df), "rows\n")
  
  cat("=== Fitting Gaussian Process Model ===\n")
  gp_result <- fit_gaussian_process(flow_df, census_df, dist_max = dist_max, 
                                    kernel_type = kernel_type, sample_size = sample_size)
  
  # Print model summary
  cat("\n", rep("=", 50), "\n", sep = "")
  cat("GP Model Summary for", city_name, "\n")
  cat(rep("=", 50), "\n", sep = "")
  cat("Kernel:", gp_result$kernel_type, "\n")
  cat("Training samples:", length(gp_result$model@ymatrix), "\n")
  cat("Features:", paste(gp_result$feature_names, collapse = ", "), "\n")
  cat("RMSE:", round(gp_result$performance$rmse, 4), "\n")
  cat("R²:", round(gp_result$performance$r_squared, 4), "\n")
  cat("Scaled features:", gp_result$use_scaled, "\n")
  
  # Create visualizations
  cat("\n=== Creating Visualizations ===\n")
  
  # Distance decay plot
  decay_result <- plot_gp_distance_decay(gp_result)
  print(decay_result$plot)
  
  # Individual feature effect plots
  feature_plots <- list()
  for(feature in c("cost", "log_num_vendors_origin", "log_pop_share_dest", "log_mean_rent_dest")) {
    feature_plots[[feature]] <- plot_gp_feature_effects(gp_result, feature)
    print(feature_plots[[feature]])
  }
  
  return(list(
    gp_result = gp_result,
    decay_predictions = decay_result$predictions,
    feature_plots = feature_plots
  ))
}

# ==========================================================
# NYC PREDICTION USING DUBAI GP MODEL
# ==========================================================

predict_nyc_flows_gp <- function(dubai_gp_result, empty_od_file, census_file, 
                                 T_total = 511000000, n_sims = 1000) {
  # Predict NYC flows using Dubai GP model with Monte Carlo sampling
  
  cat("=== NYC Flow Prediction using Dubai GP Model ===\n")
  
  # Load NYC data
  cat("Loading NYC data...\n")
  od_df <- read_csv(empty_od_file, show_col_types = FALSE) %>%
    rename(origin_h3 = origin, dest_h3 = destination)
  
  census_df <- read_csv(census_file, show_col_types = FALSE) %>%
    mutate(
      log_num_vendors = log(num_express_vendors + epsilon),
      log_pop_share = log(pop_share + epsilon),
      log_med_rent = log(med_gross_rent + epsilon)
    ) %>%
    select(h3_index, log_num_vendors, log_pop_share, log_med_rent)
  
  # Join census data to OD pairs (matching Dubai variable names)
  od_df <- od_df %>%
    left_join(
      census_df %>%
        rename(origin_h3 = h3_index,
               log_num_vendors_origin = log_num_vendors),
      by = "origin_h3"
    ) %>%
    left_join(
      census_df %>%
        rename(dest_h3 = h3_index,
               log_pop_share_dest = log_pop_share,
               log_mean_rent_dest = log_med_rent),
      by = "dest_h3"
    ) %>%
    # Rename distance to match Dubai model
    rename(cost = distance_km) %>%
    filter(!is.na(log_num_vendors_origin), 
           !is.na(log_pop_share_dest), 
           !is.na(log_mean_rent_dest),
           !is.na(cost))
  
  cat("NYC data prepared:\n")
  cat("  OD pairs:", nrow(od_df), "\n")
  cat("  Origins:", length(unique(od_df$origin_h3)), "\n")
  cat("  Destinations:", length(unique(od_df$dest_h3)), "\n")
  
  # Make base predictions using GP model
  cat("Making GP predictions...\n")
  gp_preds <- predict_gp_flows(dubai_gp_result$gp_result, od_df, return_variance = TRUE)
  
  od_df$pred_log_flow_mean <- gp_preds$mean
  od_df$pred_log_flow_var <- gp_preds$variance
  
  # Monte Carlo simulation incorporating GP uncertainty
  cat("Running Monte Carlo simulation with GP uncertainty...\n")
  set.seed(1234)
  
  # Storage for simulations
  flow_sims <- matrix(NA, nrow = nrow(od_df), ncol = n_sims)
  
  for(s in 1:n_sims) {
    # Sample from GP posterior (approximate)
    log_flow_sample <- rnorm(nrow(od_df), 
                             mean = od_df$pred_log_flow_mean,
                             sd = sqrt(od_df$pred_log_flow_var))
    
    # Transform to flow scale
    flow_sample <- pmax(exp(log_flow_sample) - 1, 0)
    
    # Normalize to total flow
    total_flow <- sum(flow_sample, na.rm = TRUE)
    if(total_flow > 0) {
      flow_sample <- flow_sample * (T_total / total_flow)
    }
    
    flow_sims[, s] <- flow_sample
    
    if(s %% 100 == 0) cat("Completed simulation", s, "/", n_sims, "\n")
  }
  
  # Compute summary statistics
  od_df$mean_flow <- rowMeans(flow_sims, na.rm = TRUE)
  od_df$sd_flow <- apply(flow_sims, 1, sd, na.rm = TRUE)
  od_df$q025 <- apply(flow_sims, 1, quantile, probs = 0.025, na.rm = TRUE)
  od_df$q975 <- apply(flow_sims, 1, quantile, probs = 0.975, na.rm = TRUE)
  
  # Print diagnostics
  cat("Prediction Summary:\n")
  cat("  Total predicted flow:", sum(od_df$mean_flow, na.rm = TRUE), "(target:", T_total, ")\n")
  cat("  Mean flow range:", range(od_df$mean_flow, na.rm = TRUE), "\n")
  cat("  Top 10 flows:\n")
  
  top_flows <- od_df %>%
    arrange(desc(mean_flow)) %>%
    head(10) %>%
    select(origin_h3, dest_h3, cost, mean_flow, sd_flow, q025, q975)
  print(top_flows)
  
  # Save results
  output_file <- "~/imperial/data/nyc_od_predicted_flows_gp.csv"
  write_csv(od_df, output_file)
  cat("Results saved to:", output_file, "\n")
  
  return(od_df)
}

# ==========================================================
# MAIN EXECUTION
# ==========================================================

main_gp <- function() {
  # Main execution function for GP analysis
  
  # File paths
  dubai_flow_file <- "~/imperial/delivery-tunnels/flow_df_for_r_dubai.csv"
  dubai_census_file <- "~/imperial/data/census_hex_dubai.csv"
  nyc_empty_od_file <- "~/imperial/data/empty_od_nyc.csv"
  nyc_census_file <- "~/imperial/data/census_hex_nyc.csv"
  
  # Run Dubai GP analysis
  cat("Starting Dubai GP model training...\n")
  dubai_gp_results <- run_dubai_gp_analysis(
    flow_filepath = dubai_flow_file,
    census_filepath = dubai_census_file,
    city_name = "Dubai",
    dist_max = 30,
    kernel_type = "rbfdot",  # Can also try "polydot"
    sample_size = 5000  # Adjust based on computational resources
  )
  
  # Use Dubai GP model to predict NYC flows
  cat("\nStarting NYC flow prediction using Dubai GP model...\n")
  nyc_results <- predict_nyc_flows_gp(
    dubai_gp_result = dubai_gp_results,
    empty_od_file = nyc_empty_od_file,
    census_file = nyc_census_file,
    T_total = 511000000,
    n_sims = 1000
  )
  
  cat("\nGP Analysis complete!\n")
  
  return(list(
    dubai = dubai_gp_results,
    nyc = nyc_results
  ))
}

# Run the GP analysis
gp_results <- main_gp()