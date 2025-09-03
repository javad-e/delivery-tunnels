library(tidyverse)
library(lme4)
library(splines)
library(scales)

epsilon <- 1e-6

# ==========================================================
# HELPER FUNCTIONS
# ==========================================================

check_required_cols <- function(df, required_cols, df_name = "dataframe") {
  missing_cols <- setdiff(required_cols, colnames(df))
  if(length(missing_cols) > 0) {
    stop(paste0("Missing columns in ", df_name, ": ", paste(missing_cols, collapse = ", ")))
  }
}

prepare_flow_data <- function(flow_df, census_df) {
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
      origin = factor(origin),
      destination = factor(destination),
      log_flow = log(flow + 1),
      log_num_vendors_origin = log(num_vendors_origin + epsilon),
      log_pop_share_dest = log(pop_share_dest + epsilon),
      log_mean_rent_dest = log(mean_rent_dest + epsilon)
    ) %>%
    droplevels()
  
  return(merged_df)
}

extract_model_params <- function(model) {
  # Extract key parameters from fitted mixed effects model
  
  # Fixed effects
  fixed_effects <- fixef(model)
  
  # Variance components
  varcomp <- as.data.frame(VarCorr(model))
  origin_var <- varcomp %>% filter(grp == "origin" & var1 == "(Intercept)") %>% pull(vcov)
  dest_var <- varcomp %>% filter(grp == "destination" & var1 == "(Intercept)") %>% pull(vcov)
  
  # Random effects means
  re <- ranef(model)
  origin_mean <- mean(re$origin[, "(Intercept)"])
  dest_mean <- mean(re$destination[, "(Intercept)"])
  
  # Model fit
  df_used <- model.frame(model)
  pred_log_flow <- predict(model, newdata = df_used)
  ss_total <- sum((df_used$log_flow - mean(df_used$log_flow))^2)
  ss_res <- sum((df_used$log_flow - pred_log_flow)^2)
  R2 <- 1 - ss_res / ss_total
  
  return(list(
    fixed_effects = fixed_effects,
    origin_var = origin_var,
    dest_var = dest_var,
    origin_mean = origin_mean,
    dest_mean = dest_mean,
    R2 = R2,
    model = model
  ))
}

# ==========================================================
# SPLINE MODEL FUNCTIONS
# ==========================================================

run_spline_model <- function(flow_df, census_df, dist_min = 0, dist_max = 25, 
                             spline_step = 5, pred_points = 200) {
  # Fit spline-based gravity model with mixed effects
  
  # Prepare data
  df <- prepare_flow_data(flow_df, census_df) %>%
    filter(cost >= dist_min, cost <= dist_max)
  
  # Create spline basis
  knots <- seq(dist_min + spline_step, dist_max - spline_step, by = spline_step)
  spline_basis <- splines::bs(df$cost, knots = knots, degree = 3, intercept = FALSE,
                              Boundary.knots = c(dist_min, dist_max))
  colnames(spline_basis) <- paste0("spline_", seq_len(ncol(spline_basis)) - 1)
  df <- bind_cols(df, as.data.frame(spline_basis))
  
  # Build model formula
  spline_terms <- paste0("spline_", seq_len(ncol(spline_basis)) - 1, collapse = " + ")
  formula_text <- paste0(
    "log_flow ~ 1 + ", spline_terms,
    " + log_num_vendors_origin + log_pop_share_dest + log_mean_rent_dest",
    " + (1 | origin) + (1 | destination)"
  )
  
  # Fit model
  model <- lmer(as.formula(formula_text), data = df)
  
  # Extract parameters
  params <- extract_model_params(model)
  
  # Generate predictions for plotting
  cost_grid <- create_prediction_grid(df, knots, dist_min, dist_max, pred_points)
  cost_grid$pred_log_flow <- predict(model, newdata = cost_grid, re.form = NA)
  cost_grid$pred_flow <- pmax(exp(cost_grid$pred_log_flow) - 1, 0)
  
  result <- cost_grid %>%
    select(cost, pred_flow) %>%
    mutate(r2 = params$R2)
  
  attr(result, "params") <- params
  return(result)
}

create_prediction_grid <- function(df, knots, dist_min, dist_max, pred_points) {
  # Create prediction grid for spline model
  
  cost_grid <- data.frame(cost = seq(dist_min, dist_max, length.out = pred_points))
  
  # Recreate spline basis
  spline_basis_grid <- splines::bs(cost_grid$cost, knots = knots, degree = 3, 
                                   intercept = FALSE, Boundary.knots = c(dist_min, dist_max))
  colnames(spline_basis_grid) <- paste0("spline_", seq_len(ncol(spline_basis_grid)) - 1)
  cost_grid <- bind_cols(cost_grid, as.data.frame(spline_basis_grid))
  
  # Set control variables to means
  cost_grid$log_num_vendors_origin <- mean(df$log_num_vendors_origin, na.rm = TRUE)
  cost_grid$log_pop_share_dest <- mean(df$log_pop_share_dest, na.rm = TRUE)
  cost_grid$log_mean_rent_dest <- mean(df$log_mean_rent_dest, na.rm = TRUE)
  
  # Set factor levels
  cost_grid$origin <- factor(rep(levels(df$origin)[1], nrow(cost_grid)), levels = levels(df$origin))
  cost_grid$destination <- factor(rep(levels(df$destination)[1], nrow(cost_grid)), levels = levels(df$destination))
  
  return(cost_grid)
}

# ==========================================================
# EXPONENTIAL MODEL FUNCTIONS
# ==========================================================

run_exponential_model <- function(flow_df, census_df) {
  # Fit exponential gravity model with mixed effects
  
  # Prepare data
  df <- prepare_flow_data(flow_df, census_df)
  
  # Fit linear mixed model (exponential in original scale)
  model <- lmer(log_flow ~ cost + log_num_vendors_origin + log_pop_share_dest + 
                  log_mean_rent_dest + (1 | origin) + (1 | destination),
                data = df)
  
  # Extract parameters
  params <- extract_model_params(model)
  
  return(params)
}

# ==========================================================
# VISUALIZATION FUNCTIONS
# ==========================================================
plot_random_effects <- function(model, city_name) {
  re <- ranef(model)
  df_re <- bind_rows(
    tibble(group = "origin", effect = re$origin[, "(Intercept)"]),
    tibble(group = "destination", effect = re$destination[, "(Intercept)"])
  )
  
  ggplot(df_re, aes(x = effect)) +
    geom_density(fill = "lightblue", alpha = 0.6) +
    facet_wrap(~group, scales = "free") +
    labs(title = paste("Random Effects Distribution for", city_name),
         x = "Random Effect", y = "Density") +
    theme_minimal()
}

plot_distance_decay <- function(spline_results, city_name) {
  # Scale predictions for plotting
  spline_results$pred_flow_scaled <- scales::rescale(spline_results$pred_flow)
  
  # Fit exponential curve
  exp_fit <- nls(pred_flow_scaled ~ exp(-b * cost), data = spline_results, start = list(b = 0.1))
  b_est <- coef(exp_fit)["b"]
  spline_results$pred_exp_scaled <- exp(-b_est * spline_results$cost)
  
  plot <- ggplot(spline_results, aes(x = cost)) +
    geom_line(aes(y = pred_flow_scaled, color = "Spline"), size = 1.2) +
    geom_line(aes(y = pred_exp_scaled, color = "Exp(-b*cost)"), linetype = "dashed", size = 1.2) +
    scale_color_manual("", values = c("Spline" = "red", "Exp(-b*cost)" = "blue")) +
    labs(title = paste("Predicted Flow vs Cost for", city_name),
         x = "Cost", y = "Predicted Flow (rescaled for plotting)") +
    theme_minimal()
  
  return(list(plot = plot, b_estimate = b_est))
}

# ==========================================================
# MAIN ANALYSIS WORKFLOW
# ==========================================================

run_dubai_analysis <- function(flow_filepath, census_filepath, city_name = "Dubai", 
                               dist_max = 30, spline_step = 5) {
  # Load data
  flow_df <- read_csv(path.expand(flow_filepath), show_col_types = FALSE)
  census_df <- read_csv(path.expand(census_filepath), show_col_types = FALSE)
  
  # Fit models
  cat("=== Fitting Spline Model ===\n")
  spline_results <- run_spline_model(flow_df, census_df, dist_max = dist_max, spline_step = spline_step)
  spline_params <- attr(spline_results, "params")
  
  cat("=== Fitting Exponential Model ===\n")
  exp_params <- run_exponential_model(flow_df, census_df)
  
  # Print results
  print_model_summary(spline_params, "Spline Model", city_name)
  print_model_summary(exp_params, "Exponential Model", city_name)
  
  # Create plots
  p1 <- plot_random_effects(exp_params$model, city_name)
  decay_plot <- plot_distance_decay(spline_results, city_name)
  
  print(p1)
  print(decay_plot$plot)
  cat(sprintf("Estimated b for y = exp(-b*cost): %.4f\n", decay_plot$b_estimate))
  
  return(list(
    spline_params = spline_params,
    exp_params = exp_params,
    spline_results = spline_results
  ))
}

print_model_summary <- function(params, model_name, city_name) {
  cat("\n", rep("=", 50), "\n", sep = "")
  cat(model_name, "Summary for", city_name, "\n")
  cat(rep("=", 50), "\n", sep = "")
  
  cat(sprintf("R²: %.4f\n", params$R2))
  cat(sprintf("Origin RE variance: %.6f\n", params$origin_var))
  cat(sprintf("Destination RE variance: %.6f\n", params$dest_var))
  cat(sprintf("Origin RE mean: %.6f\n", params$origin_mean))
  cat(sprintf("Destination RE mean: %.6f\n", params$dest_mean))
  
  cat("\nFixed Effects Coefficients:\n")
  print(params$fixed_effects)
}

# ==========================================================
# NYC PREDICTION USING DUBAI PARAMETERS
# ==========================================================

predict_nyc_flows <- function(dubai_params, empty_od_file, census_file, T_total = 511000000, n_sims = 1000) {
  # Extract parameters from Dubai exponential model
  fixed_effects <- dubai_params$fixed_effects
  sigma_u <- sqrt(dubai_params$origin_var)    # origin RE SD
  sigma_v <- sqrt(dubai_params$dest_var)      # destination RE SD
  
  cat("Using Dubai model parameters:\n")
  cat(sprintf("Intercept: %.6f\n", fixed_effects["(Intercept)"]))
  cat(sprintf("Distance coefficient: %.6f\n", fixed_effects["cost"]))
  cat(sprintf("Origin vendors coefficient: %.6f\n", fixed_effects["log_num_vendors_origin"]))
  cat(sprintf("Dest population coefficient: %.6f\n", fixed_effects["log_pop_share_dest"]))
  cat(sprintf("Dest rent coefficient: %.6f\n", fixed_effects["log_mean_rent_dest"]))
  cat(sprintf("Origin RE SD: %.6f\n", sigma_u))
  cat(sprintf("Destination RE SD: %.6f\n", sigma_v))
  
  # Load NYC data
  cat("\nLoading NYC data...\n")
  od_df <- read_csv(empty_od_file, show_col_types = FALSE) %>%
    rename(origin_h3 = origin, dest_h3 = destination)
  
  census_df <- read_csv(census_file, show_col_types = FALSE) %>%
    mutate(
      log_num_vendors = log(num_express_vendors + epsilon),
      log_pop_share = log(pop_share + epsilon),
      log_med_rent = log(med_gross_rent + epsilon)
    ) %>%
    select(h3_index, log_num_vendors, log_pop_share, log_med_rent)
  
  # Join census data to OD pairs
  od_df <- od_df %>%
    left_join(
      census_df %>%
        rename(origin_h3 = h3_index,
               log_num_vendors_origin = log_num_vendors,
               log_pop_share_origin = log_pop_share),
      by = "origin_h3"
    ) %>%
    left_join(
      census_df %>%
        rename(dest_h3 = h3_index,
               log_num_vendors_dest = log_num_vendors,
               log_pop_share_dest = log_pop_share,
               log_med_rent_dest = log_med_rent),
      by = "dest_h3"
    )
  
  # Compute base linear predictor
  od_df <- od_df %>%
    mutate(LP_base = fixed_effects["(Intercept)"] + 
             fixed_effects["cost"] * distance_km +
             fixed_effects["log_num_vendors_origin"] * log_num_vendors_origin +
             fixed_effects["log_pop_share_dest"] * log_pop_share_dest +
             fixed_effects["log_mean_rent_dest"] * log_med_rent_dest)
  
  cat("Data preparation complete.\n")
  cat("Origins:", length(unique(od_df$origin_h3)), 
      "Destinations:", length(unique(od_df$dest_h3)), "\n")
  cat("OD pairs:", nrow(od_df), "\n")
  cat("Missing values in LP_base:", sum(is.na(od_df$LP_base)), "\n")
  
  cat("Running Monte Carlo simulation...\n")
  results <- run_monte_carlo_prediction(od_df, sigma_u, sigma_v, T_total, n_sims)
  
  output_file <- "~/imperial/data/nyc_od_predicted_express_flows.csv"
  write_csv(results, output_file)
  cat("Results saved to:", output_file, "\n")
  
  return(results)
}

run_monte_carlo_prediction <- function(od_df, sigma_u, sigma_v, T_total, n_sims) {
  # Create index mappings
  origins <- unique(od_df$origin_h3)
  dests <- unique(od_df$dest_h3)
  n_orig <- length(origins)
  n_dest <- length(dests)
  
  origin_map <- setNames(seq_len(n_orig), origins)
  dest_map <- setNames(seq_len(n_dest), dests)
  
  # Reshape LP_base to matrix
  LP_mat <- matrix(NA_real_, nrow = n_orig, ncol = n_dest)
  for(i in seq_len(nrow(od_df))) {
    ori_idx <- origin_map[[od_df$origin_h3[i]]]
    dest_idx <- dest_map[[od_df$dest_h3[i]]]
    LP_mat[ori_idx, dest_idx] <- od_df$LP_base[i]
  }
  
  # Monte Carlo simulations
  set.seed(1234)
  preds_array <- array(NA_real_, dim = c(n_orig, n_dest, n_sims))
  
  for(s in seq_len(n_sims)) {
    # Sample random effects
    u <- rnorm(n_orig, mean = 0, sd = sigma_u)
    v <- rnorm(n_dest, mean = 0, sd = sigma_v)
    
    # Add random effects to linear predictor
    LP_sim <- LP_mat + outer(u, rep(1, n_dest)) + outer(rep(1, n_orig), v)
    
    # Transform to flow scale and normalize
    pred_sim <- exp(LP_sim)
    total_pred <- sum(pred_sim, na.rm = TRUE)
    
    if(total_pred > 0 && !is.na(total_pred)) {
      scale_factor <- T_total / total_pred
      pred_sim <- pred_sim * scale_factor
      preds_array[,,s] <- pred_sim
    }
    
    if(s %% 100 == 0) cat("Completed simulation", s, "/", n_sims, "\n")
  }
  
  # Compute summary statistics
  mean_pred <- apply(preds_array, c(1,2), function(x) mean(x, na.rm = TRUE))
  sd_pred <- apply(preds_array, c(1,2), function(x) sd(x, na.rm = TRUE))
  q_low <- apply(preds_array, c(1,2), function(x) quantile(x, probs = 0.025, na.rm = TRUE))
  q_high <- apply(preds_array, c(1,2), function(x) quantile(x, probs = 0.975, na.rm = TRUE))
  
  # Map back to original dataframe
  od_df$mean_flow <- NA
  od_df$sd_flow <- NA
  od_df$q025 <- NA
  od_df$q975 <- NA
  
  for(i in seq_len(nrow(od_df))) {
    ori_idx <- origin_map[[od_df$origin_h3[i]]]
    dest_idx <- dest_map[[od_df$dest_h3[i]]]
    od_df$mean_flow[i] <- mean_pred[ori_idx, dest_idx]
    od_df$sd_flow[i] <- sd_pred[ori_idx, dest_idx]
    od_df$q025[i] <- q_low[ori_idx, dest_idx]
    od_df$q975[i] <- q_high[ori_idx, dest_idx]
  }
  
  # Print diagnostics
  cat("Total predicted flow:", sum(od_df$mean_flow, na.rm = TRUE), "(target:", T_total, ")\n")
  cat("Top 10 flows:\n")
  print(od_df %>% 
          arrange(desc(mean_flow)) %>% 
          head(10) %>%
          select(origin_h3, dest_h3, distance_km, mean_flow, sd_flow))
  
  return(od_df)
}

# ==========================================================
# MAIN EXECUTION
# ==========================================================

main <- function() {
  # File paths
  dubai_flow_file <- "~/imperial/delivery-tunnels/flow_df_for_r_dubai.csv"
  dubai_census_file <- "~/imperial/data/census_hex_dubai.csv"
  nyc_empty_od_file <- "~/imperial/data/empty_od_nyc.csv"
  nyc_census_file <- "~/imperial/data/census_hex_nyc.csv"
  
  # Run Dubai analysis
  cat("Starting Dubai model training...\n")
  dubai_results <- run_dubai_analysis(dubai_flow_file, dubai_census_file, "Dubai")
  
  # Use Dubai parameters to predict NYC flows
  cat("\nStarting NYC flow prediction...\n")
  nyc_results <- predict_nyc_flows(
    dubai_params = dubai_results$exp_params,
    empty_od_file = nyc_empty_od_file,
    census_file = nyc_census_file,
    T_total = 511000000,
    n_sims = 1000
  )
  
  cat("\nAnalysis complete!\n")
  
  return(list(dubai = dubai_results, nyc = nyc_results))
}

# Run the analysis
results <- main()
