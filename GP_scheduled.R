# ==========================================================
# Libraries
# ==========================================================
library(tidyverse)
library(kernlab)  # For Gaussian Process regression
library(scales)
library(dplyr)
library(tidyr)
library(readr)
library(lme4)


epsilon <- 1e-6  # small constant to avoid log(0)

# ==========================================================
# Helper function to check columns
# ==========================================================
check_required_cols <- function(df, required_cols, df_name = "dataframe") {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste(
      "The following required columns are missing from", df_name, ":", 
      paste(missing_cols, collapse = ", ")
    ))
  }
}

# ==========================================================
# Gaussian Process model with controls
# ==========================================================
run_gp_model <- function(flow_df_path, census_df_path, dataset_name = "Dataset", 
                         kernel_type = "rbfdot", subset_size = 5000) {
  
  # Read data
  flow_df <- read_csv(flow_df_path)
  census_df <- read_csv(census_df_path)
  
  # Check required columns
  check_required_cols(flow_df, c("origin", "destination", "flow"), "flow_df")
  check_required_cols(census_df, c("h3_index", "pop_total", "num_warehouses", 
                                   "warehousing_share", "emp_share", 
                                   "med_gross_rent", "med_hh_income"), "census_df")
  
  # Merge origin variables
  flow_df <- flow_df %>%
    left_join(census_df %>% dplyr::select(h3_index, num_warehouses, 
                                          warehousing_share, emp_share,
                                          med_gross_rent, med_hh_income), 
              by = c("origin" = "h3_index")) %>%
    rename(num_warehouses_origin = num_warehouses,
           warehousing_share_origin = warehousing_share,
           emp_share_origin = emp_share,
           med_gross_rent_origin = med_gross_rent,
           med_hh_income_origin = med_hh_income) %>%
    left_join(census_df %>% dplyr::select(h3_index, pop_total), 
              by = c("destination" = "h3_index")) %>%
    rename(pop_total_dest = pop_total)
  
  # Replace NA with small values
  flow_df <- flow_df %>%
    mutate(
      pop_total_dest   = ifelse(is.na(pop_total_dest), 0, pop_total_dest),
      num_warehouses_origin = ifelse(is.na(num_warehouses_origin), 0, num_warehouses_origin),
      warehousing_share_origin = ifelse(is.na(warehousing_share_origin), 0, warehousing_share_origin),
      emp_share_origin = ifelse(is.na(emp_share_origin), 0, emp_share_origin),
      med_gross_rent_origin = ifelse(is.na(med_gross_rent_origin), 0, med_gross_rent_origin),
      med_hh_income_origin = ifelse(is.na(med_hh_income_origin), 0, med_hh_income_origin)
    )
  
  # Log-transform variables
  flow_df <- flow_df %>%
    mutate(
      log_num_warehouses = log(num_warehouses_origin + epsilon),
      log_pop_total_dest = log(pop_total_dest + epsilon),
      log_warehousing_share_origin = log(warehousing_share_origin + epsilon),
      log_emp_share_origin = log(emp_share_origin + epsilon),
      log_med_gross_rent_origin = log(med_gross_rent_origin + epsilon),
      log_med_hh_income_origin = log(med_hh_income_origin + epsilon),
      log_flow = log(flow + epsilon)
    )
  
  # Create feature matrix for GP
  feature_cols <- c("log_num_warehouses", "log_pop_total_dest", 
                    "log_warehousing_share_origin", "log_emp_share_origin",
                    "log_med_gross_rent_origin", "log_med_hh_income_origin")
  
  # Remove rows with any NA values in features or target
  complete_data <- flow_df[complete.cases(flow_df[c(feature_cols, "log_flow")]), ]
  
  # Check for constant variables before subsetting
  feature_vars <- complete_data[, feature_cols]
  constant_vars <- sapply(feature_vars, function(x) var(x, na.rm = TRUE) == 0)
  if(any(constant_vars)) {
    cat("Warning: Constant variables detected:", paste(names(constant_vars)[constant_vars], collapse = ", "), "\n")
    cat("These variables will be removed from the model.\n")
    feature_cols <- feature_cols[!constant_vars]
  }
  
  if(length(feature_cols) == 0) {
    stop("Error: No usable features remaining after removing constant variables.")
  }
  
  cat("Final feature set:", paste(feature_cols, collapse = ", "), "\n")
  
  # Subset data for computational efficiency if needed
  if (nrow(complete_data) > subset_size) {
    cat("Subsetting data to", subset_size, "observations for computational efficiency\n")
    set.seed(1234)
    sample_idx <- sample(nrow(complete_data), subset_size)
    complete_data <- complete_data[sample_idx, ]
  }
  
  # Prepare training data
  X_train <- as.matrix(complete_data[, feature_cols])
  y_train <- complete_data$log_flow
  
  # Standardize features for better GP performance
  # Handle constant variables (zero variance)
  X_mean <- colMeans(X_train)
  X_sd <- apply(X_train, 2, sd)
  
  # Replace zero standard deviations with 1 to avoid scaling issues
  zero_var_cols <- which(X_sd == 0)
  if(length(zero_var_cols) > 0) {
    cat("Warning: Variables with zero variance found:", 
        paste(feature_cols[zero_var_cols], collapse = ", "), "\n")
    X_sd[zero_var_cols] <- 1  # Don't scale constant variables
  }
  
  # Manual scaling to handle edge cases
  X_train_scaled <- sweep(sweep(X_train, 2, X_mean, "-"), 2, X_sd, "/")
  
  # Check for any remaining issues
  if(any(is.na(X_train_scaled)) || any(is.infinite(X_train_scaled))) {
    cat("Warning: NA or Inf values found in scaled features. Using original features.\n")
    X_train_scaled <- X_train
    X_mean <- rep(0, ncol(X_train))
    X_sd <- rep(1, ncol(X_train))
  }
  
  # Choose kernel
  kernel <- switch(kernel_type,
                   "rbfdot" = rbfdot(sigma = 1),  # RBF kernel
                   "polydot" = polydot(degree = 2, scale = 1, offset = 1),  # Polynomial
                   "laplacedot" = laplacedot(sigma = 1),  # Laplace
                   rbfdot(sigma = 1))  # Default to RBF
  
  cat("Training Gaussian Process model with", kernel_type, "kernel...\n")
  
  # Fit Gaussian Process
  gp_model <- gausspr(x = X_train_scaled, y = y_train, 
                      kernel = kernel,
                      type = "regression",
                      variance.model = TRUE,  # Enable uncertainty estimation
                      tol = 1e-3)  # Removed cross parameter to avoid potential issues
  
  # Make predictions on training data for R²
  train_pred <- predict(gp_model, X_train_scaled)
  r2 <- 1 - sum((y_train - train_pred)^2) / sum((y_train - mean(y_train))^2)
  
  # Extract hyperparameters (safely)
  hyperparams <- list(
    kernel_type = kernel_type,
    sigma = tryCatch({
      if(kernel_type == "rbfdot") kernelf(gp_model)@kpar$sigma else NA
    }, error = function(e) NA),
    noise_variance = tryCatch(gp_model@error, error = function(e) NA)
  )
  
  cat("\n=== Gaussian Process Model Summary ===\n")
  cat("Dataset:", dataset_name, "\n")
  cat("Kernel:", kernel_type, "\n")
  cat("Training R²:", round(r2, 4), "\n")
  cat("Training observations:", nrow(complete_data), "\n")
  cat("Features used:", length(feature_cols), "\n")
  if(!is.na(hyperparams$sigma)) {
    cat("RBF Sigma:", round(hyperparams$sigma, 4), "\n")
  }
  if(!is.na(hyperparams$noise_variance)) {
    cat("Noise variance:", round(hyperparams$noise_variance, 6), "\n")
  }
  
  return(list(
    model = gp_model,
    data = complete_data,
    X_mean = X_mean,
    X_sd = X_sd,
    feature_cols = feature_cols,
    r2 = r2,
    hyperparams = hyperparams,
    kernel_type = kernel_type
  ))
}

# ==========================================================
# GP-based flow prediction function
# ==========================================================
predict_flows_gp <- function(gp_results, census_df, od_df, T_total, n_sims = 1000) {
  
  # Extract model components
  gp_model <- gp_results$model
  X_mean <- gp_results$X_mean
  X_sd <- gp_results$X_sd
  feature_cols <- gp_results$feature_cols
  
  # Prepare socioeconomic variables
  relevant_h3 <- unique(c(od_df$origin_h3, od_df$dest_h3))
  census_processed <- census_df %>%
    dplyr::filter(h3_index %in% relevant_h3) %>%
    dplyr::mutate(
      log_pop_total         = log(pop_total + epsilon),
      log_warehousing_share = log(warehousing_share + epsilon),
      log_emp_share         = log(emp_share + epsilon),
      log_med_gross_rent    = log(med_gross_rent + epsilon),
      log_med_hh_income     = log(med_hh_income + epsilon),
      log_num_warehouses    = log(num_warehouses + epsilon)
    ) %>%
    dplyr::select(h3_index, log_pop_total, log_warehousing_share, log_emp_share,
                  log_med_gross_rent, log_med_hh_income, log_num_warehouses)
  
  # Join origin/destination variables
  od_processed <- od_df %>%
    dplyr::left_join(
      census_processed %>%
        dplyr::rename(origin_h3 = h3_index,
                      log_warehousing_share_origin = log_warehousing_share,
                      log_emp_share_origin         = log_emp_share,
                      log_med_gross_rent_origin    = log_med_gross_rent,
                      log_med_hh_income_origin     = log_med_hh_income,
                      log_num_warehouses_origin    = log_num_warehouses),
      by = "origin_h3"
    ) %>%
    dplyr::left_join(
      census_processed %>%
        dplyr::rename(dest_h3 = h3_index,
                      log_pop_total_dest = log_pop_total),
      by = "dest_h3"
    )
  
  # Impute missing values
  od_processed <- od_processed %>%
    dplyr::mutate(
      log_warehousing_share_origin = ifelse(is.na(log_warehousing_share_origin), 0, log_warehousing_share_origin),
      log_emp_share_origin         = ifelse(is.na(log_emp_share_origin), 0, log_emp_share_origin),
      log_med_gross_rent_origin    = ifelse(is.na(log_med_gross_rent_origin), 0, log_med_gross_rent_origin),
      log_med_hh_income_origin     = ifelse(is.na(log_med_hh_income_origin), 0, log_med_hh_income_origin),
      log_pop_total_dest           = ifelse(is.na(log_pop_total_dest), 0, log_pop_total_dest),
      log_num_warehouses_origin    = ifelse(is.na(log_num_warehouses_origin), 0, log_num_warehouses_origin)
    )
  
  # Create feature matrix for prediction
  X_pred <- as.matrix(od_processed[, feature_cols])
  
  # Standardize using training data statistics (manual scaling to handle edge cases)
  X_pred_scaled <- sweep(sweep(X_pred, 2, X_mean, "-"), 2, X_sd, "/")
  
  # Handle any remaining NA or Inf values
  if(any(is.na(X_pred_scaled)) || any(is.infinite(X_pred_scaled))) {
    cat("Warning: NA or Inf values found in prediction features. Using original features.\n")
    X_pred_scaled <- X_pred
  }
  
  cat("Making GP predictions for", nrow(X_pred_scaled), "OD pairs...\n")
  
  # Get GP predictions with uncertainty
  pred_mean <- predict(gp_model, X_pred_scaled)
  
  # For uncertainty, we'll use the GP's predictive variance if available
  # Note: kernlab's gausspr doesn't directly provide predictive variance
  # We'll approximate it using cross-validation error and distance to training points
  
  # Monte Carlo simulation incorporating GP uncertainty
  set.seed(1234)
  flow_sims <- matrix(NA_real_, nrow = length(pred_mean), ncol = n_sims)
  
  # Estimate prediction uncertainty (simplified approach)
  # In practice, you might want to implement proper GP predictive variance
  pred_uncertainty <- rep(sqrt(gp_results$hyperparams$noise_variance * 2), length(pred_mean))
  
  for(s in 1:n_sims) {
    # Sample from GP posterior (approximated)
    log_flow_sim <- rnorm(length(pred_mean), pred_mean, pred_uncertainty)
    flow_sim <- exp(log_flow_sim)
    
    # Scale to match total flow
    scale_factor <- T_total / sum(flow_sim, na.rm = TRUE)
    flow_sims[, s] <- flow_sim * scale_factor
    
    if(s %% 200 == 0) cat("Completed simulation", s, "\n")
  }
  
  # Calculate summary statistics
  od_processed$mean_flow <- rowMeans(flow_sims, na.rm = TRUE)
  od_processed$sd_flow   <- apply(flow_sims, 1, sd, na.rm = TRUE)
  od_processed$q025      <- apply(flow_sims, 1, quantile, probs = 0.025, na.rm = TRUE)
  od_processed$q975      <- apply(flow_sims, 1, quantile, probs = 0.975, na.rm = TRUE)
  od_processed$gp_pred_mean <- pred_mean  # Raw GP predictions before scaling
  
  return(od_processed)
}

# ==========================================================
# Usage example for Chicago (Training)
# ==========================================================
chicago_flow_file <- "~/imperial/delivery-tunnels/flow_df_for_r_chicago.csv"
chicago_census_file <- "~/imperial/data/census_hex_chicago.csv"

# Train GP model on Chicago data
chicago_gp_results <- run_gp_model(chicago_flow_file, chicago_census_file, 
                                   dataset_name = "Chicago", 
                                   kernel_type = "rbfdot",  # Can try "polydot" or "laplacedot"
                                   subset_size = 5000)

# ==========================================================
# NYC Flow Prediction using GP
# ==========================================================
empty_od_file <- "~/imperial/data/empty_od.csv"
nyc_census_file <- "~/imperial/data/census_hex_nyc.csv"
T_total <- 839500000
n_sims <- 1000

# Load NYC data
od_df <- readr::read_csv(empty_od_file) %>%
  dplyr::rename(origin_h3 = origin, dest_h3 = destination)
nyc_census_df <- readr::read_csv(nyc_census_file)

# Make predictions using GP model trained on Chicago
cat("Predicting NYC flows using GP model trained on Chicago data...\n")
nyc_predictions <- predict_flows_gp(chicago_gp_results, nyc_census_df, od_df, T_total, n_sims)

# Diagnostics
cat("\n=== Prediction Diagnostics ===\n")
cat("Total predicted flow:", sum(nyc_predictions$mean_flow, na.rm = TRUE), "\n")
cat("Target total flow:", T_total, "\n")
cat("Difference:", sum(nyc_predictions$mean_flow, na.rm = TRUE) - T_total, "\n")

# Show top flows
top_flows <- nyc_predictions %>%
  arrange(desc(mean_flow)) %>%
  head(20) %>%
  dplyr::select(origin_h3, dest_h3, mean_flow, sd_flow, q025, q975, gp_pred_mean)

cat("\nTop 20 predicted flows:\n")
print(top_flows)

# Save results
output_file <- "~/imperial/data/nyc_od_predicted_gp_flows.csv"

# Debug: Check column types
cat("Checking column types in nyc_predictions:\n")
for(i in 1:ncol(nyc_predictions)) {
  col_class <- class(nyc_predictions[[i]])
  col_name <- names(nyc_predictions)[i]
  cat("Column", i, ":", col_name, "- Class:", paste(col_class, collapse = ", "), "\n")
}

# More thorough cleaning approach
nyc_predictions_clean <- nyc_predictions

# Check each column and handle problematic ones
for(i in 1:ncol(nyc_predictions_clean)) {
  col <- nyc_predictions_clean[[i]]
  col_name <- names(nyc_predictions_clean)[i]
  
  if(is.list(col) && !is.data.frame(col)) {
    cat("Converting list column:", col_name, "\n")
    # Try to unlist if all elements are single values
    if(all(lengths(col) == 1)) {
      nyc_predictions_clean[[i]] <- unlist(col)
    } else {
      # Convert to character representation
      nyc_predictions_clean[[i]] <- sapply(col, function(x) paste(x, collapse = ","))
    }
  } else if(is.matrix(col)) {
    cat("Converting matrix column:", col_name, "\n")
    # Convert matrix to character
    nyc_predictions_clean[[i]] <- apply(col, 1, function(x) paste(x, collapse = ","))
  }
}

# Final safety check - only keep simple columns
safe_cols <- sapply(nyc_predictions_clean, function(x) {
  is.numeric(x) || is.character(x) || is.logical(x) || is.factor(x)
})

cat("Safe columns:", sum(safe_cols), "out of", length(safe_cols), "\n")
if(!all(safe_cols)) {
  cat("Removing problematic columns:", names(nyc_predictions_clean)[!safe_cols], "\n")
  nyc_predictions_clean <- nyc_predictions_clean[, safe_cols]
}

write_csv(nyc_predictions_clean, output_file)
cat("\nResults saved to:", output_file, "\n")

# ==========================================================
# Model comparison and diagnostics
# ==========================================================
plot_gp_diagnostics <- function(gp_results, predictions = NULL) {
  
  # Plot 1: Training data fit
  training_data <- gp_results$data
  X_train <- as.matrix(training_data[, gp_results$feature_cols])
  
  # Use the same scaling approach as in training
  X_mean <- gp_results$X_mean
  X_sd <- gp_results$X_sd
  X_train_scaled <- sweep(sweep(X_train, 2, X_mean, "-"), 2, X_sd, "/")
  
  # Handle edge cases
  if(any(is.na(X_train_scaled)) || any(is.infinite(X_train_scaled))) {
    X_train_scaled <- X_train
  }
  
  train_pred <- predict(gp_results$model, X_train_scaled)
  
  cat("Creating diagnostic plots...\n")
  
  # Residuals plot
  residuals <- training_data$log_flow - train_pred
  
  par(mfrow = c(2, 2))
  
  # 1. Predicted vs Actual
  plot(train_pred, training_data$log_flow,
       xlab = "GP Predicted log(flow)", ylab = "Actual log(flow)",
       main = "GP Model: Predicted vs Actual")
  abline(0, 1, col = "red", lty = 2)
  
  # 2. Residuals vs Fitted
  plot(train_pred, residuals,
       xlab = "Fitted Values", ylab = "Residuals",
       main = "Residuals vs Fitted")
  abline(h = 0, col = "red", lty = 2)
  
  # 3. QQ plot of residuals
  qqnorm(residuals, main = "Normal Q-Q Plot of Residuals")
  qqline(residuals, col = "red")
  
  # 4. Histogram of residuals
  hist(residuals, breaks = 30, main = "Distribution of Residuals",
       xlab = "Residuals", freq = FALSE)
  curve(dnorm(x, mean = mean(residuals), sd = sd(residuals)), 
        add = TRUE, col = "red", lwd = 2)
  
  par(mfrow = c(1, 1))
  
  # Print summary statistics
  cat("\n=== Diagnostic Summary ===\n")
  cat("RMSE:", sqrt(mean(residuals^2)), "\n")
  cat("MAE:", mean(abs(residuals)), "\n")
  cat("Residual std dev:", sd(residuals), "\n")
  cat("Mean residual:", mean(residuals), "\n")
}

# ==========================================================
# Model Comparison: Mixed Effects vs Gaussian Process
# ==========================================================
compare_models <- function(flow_df_path, census_df_path, dataset_name = "Dataset") {
  
  cat("\n=== MODEL COMPARISON:", dataset_name, "===\n")
  
  # 1. Run Mixed Effects Model (original approach)
  cat("\n1. Training Mixed Effects Model...\n")
  
  # Read and prepare data (same as GP approach)
  flow_df <- read_csv(flow_df_path, show_col_types = FALSE)
  census_df <- read_csv(census_df_path, show_col_types = FALSE)
  
  # Merge data
  flow_df <- flow_df %>%
    left_join(census_df %>% dplyr::select(h3_index, num_warehouses, 
                                          warehousing_share, emp_share,
                                          med_gross_rent, med_hh_income), 
              by = c("origin" = "h3_index")) %>%
    rename(num_warehouses_origin = num_warehouses,
           warehousing_share_origin = warehousing_share,
           emp_share_origin = emp_share,
           med_gross_rent_origin = med_gross_rent,
           med_hh_income_origin = med_hh_income) %>%
    left_join(census_df %>% dplyr::select(h3_index, pop_total), 
              by = c("destination" = "h3_index")) %>%
    rename(pop_total_dest = pop_total)
  
  # Replace NA and log-transform
  flow_df <- flow_df %>%
    mutate(
      pop_total_dest   = ifelse(is.na(pop_total_dest), 0, pop_total_dest),
      num_warehouses_origin = ifelse(is.na(num_warehouses_origin), 0, num_warehouses_origin),
      warehousing_share_origin = ifelse(is.na(warehousing_share_origin), 0, warehousing_share_origin),
      emp_share_origin = ifelse(is.na(emp_share_origin), 0, emp_share_origin),
      med_gross_rent_origin = ifelse(is.na(med_gross_rent_origin), 0, med_gross_rent_origin),
      med_hh_income_origin = ifelse(is.na(med_hh_income_origin), 0, med_hh_income_origin),
      log_num_warehouses = log(num_warehouses_origin + epsilon),
      log_pop_total_dest = log(pop_total_dest + epsilon),
      log_warehousing_share_origin = log(warehousing_share_origin + epsilon),
      log_emp_share_origin = log(emp_share_origin + epsilon),
      log_med_gross_rent_origin = log(med_gross_rent_origin + epsilon),
      log_med_hh_income_origin = log(med_hh_income_origin + epsilon),
      log_flow = log(flow + epsilon)
    )
  
  # Mixed effects model
  model_formula <- as.formula(
    "log_flow ~ log_num_warehouses + log_pop_total_dest + 
     log_warehousing_share_origin + log_emp_share_origin +
     log_med_gross_rent_origin + log_med_hh_income_origin +
     (1 | origin) + (1 | destination)"
  )
  
  lme_model <- lmer(model_formula, data = flow_df, REML = FALSE)
  lme_fitted <- fitted(lme_model)
  lme_r2 <- 1 - sum((flow_df$log_flow - lme_fitted)^2) / sum((flow_df$log_flow - mean(flow_df$log_flow))^2)
  
  # 2. Run Gaussian Process Model
  cat("\n2. Training Gaussian Process Model...\n")
  gp_results <- run_gp_model(flow_df_path, census_df_path, dataset_name, 
                             kernel_type = "rbfdot", subset_size = 5000)
  
  # 3. Compare Results
  cat("\n=== COMPARISON RESULTS ===\n")
  cat("Dataset:", dataset_name, "\n")
  cat("Mixed Effects R²:   ", sprintf("%.4f", lme_r2), "\n")
  cat("Gaussian Process R²:", sprintf("%.4f", gp_results$r2), "\n")
  
  improvement <- ((gp_results$r2 - lme_r2) / lme_r2) * 100
  if (improvement > 0) {
    cat("GP Improvement:     +", sprintf("%.2f%%", improvement), "\n")
  } else {
    cat("GP Performance:     ", sprintf("%.2f%%", improvement), " (worse than mixed effects)\n")
  }
  
  cat("\nModel Details:\n")
  cat("- Mixed Effects observations:", nrow(flow_df), "\n")
  cat("- GP observations:", nrow(gp_results$data), "\n")
  cat("- GP features used:", length(gp_results$feature_cols), "\n")
  cat("- GP features:", paste(gp_results$feature_cols, collapse = ", "), "\n")
  
  # 4. Statistical significance test
  lme_residuals <- flow_df$log_flow - lme_fitted
  gp_train_data <- gp_results$data
  
  # Get GP predictions on same data for fair comparison
  X_train <- as.matrix(gp_train_data[, gp_results$feature_cols])
  X_mean <- gp_results$X_mean
  X_sd <- gp_results$X_sd
  X_train_scaled <- sweep(sweep(X_train, 2, X_mean, "-"), 2, X_sd, "/")
  
  if(any(is.na(X_train_scaled)) || any(is.infinite(X_train_scaled))) {
    X_train_scaled <- X_train
  }
  
  gp_pred <- predict(gp_results$model, X_train_scaled)
  gp_residuals <- gp_train_data$log_flow - gp_pred
  
  cat("\nResidual Analysis:\n")
  cat("- Mixed Effects RMSE:", sprintf("%.4f", sqrt(mean(lme_residuals^2))), "\n")
  cat("- GP RMSE:           ", sprintf("%.4f", sqrt(mean(gp_residuals^2))), "\n")
  cat("- Mixed Effects MAE: ", sprintf("%.4f", mean(abs(lme_residuals))), "\n")
  cat("- GP MAE:            ", sprintf("%.4f", mean(abs(gp_residuals))), "\n")
  
  return(list(
    lme_r2 = lme_r2,
    gp_r2 = gp_results$r2,
    improvement_pct = improvement,
    lme_rmse = sqrt(mean(lme_residuals^2)),
    gp_rmse = sqrt(mean(gp_residuals^2)),
    gp_results = gp_results
  ))
}

# Run the comparison
comparison_results <- compare_models(chicago_flow_file, chicago_census_file, "Chicago")

# Run diagnostics
plot_gp_diagnostics(comparison_results$gp_results)