library(tidyverse)
library(lme4)
library(splines)
library(merTools)
library(scales)

epsilon <- 1e-6  # small constant to avoid log(0)

check_required_cols <- function(df, required_cols, df_name = "dataframe") {
  missing_cols <- setdiff(required_cols, colnames(df))
  if(length(missing_cols) > 0) {
    stop(paste0("Missing columns in ", df_name, ": ", paste(missing_cols, collapse = ", ")))
  }
}

run_spline_model_baseline <- function(flow_df, dist_min = 0, dist_max = 25, spline_step = 5, pred_points = 200) {
  
  # Check columns
  check_required_cols(flow_df, c("cost", "origin", "destination", "flow"), "flow_df")
  
  df <- flow_df %>%
    dplyr::filter(cost >= dist_min, cost <= dist_max) %>%
    dplyr::mutate(
      origin       = factor(origin),
      destination  = factor(destination),
      log_flow     = log(flow + 1)
    ) %>%
    droplevels()
  
  knots <- seq(dist_min + spline_step, dist_max - spline_step, by = spline_step)
  spline_basis <- splines::bs(df$cost,
                              knots = knots,
                              degree = 3,
                              intercept = FALSE,
                              Boundary.knots = c(dist_min, dist_max))
  colnames(spline_basis) <- paste0("spline_", seq_len(ncol(spline_basis)) - 1)
  df <- dplyr::bind_cols(df, as.data.frame(spline_basis))
  
  spline_terms <- paste0("spline_", seq_len(ncol(spline_basis)) - 1, collapse = " + ")
  formula_text <- paste0("log_flow ~ 1 + ", spline_terms, " + (1 | origin) + (1 | destination)")
  
  model <- lme4::lmer(as.formula(formula_text), data = df)
  
  cost_grid <- data.frame(cost = seq(dist_min, dist_max, length.out = pred_points))
  spline_basis_grid <- splines::bs(cost_grid$cost,
                                   knots = knots,
                                   degree = 3,
                                   intercept = FALSE,
                                   Boundary.knots = c(dist_min, dist_max))
  colnames(spline_basis_grid) <- paste0("spline_", seq_len(ncol(spline_basis_grid)) - 1)
  cost_grid <- dplyr::bind_cols(cost_grid, as.data.frame(spline_basis_grid))
  
  # Use the most common origin/destination levels for prediction grid
  most_common_origin <- names(which.max(table(df$origin)))
  most_common_dest <- names(which.max(table(df$destination)))
  
  cost_grid$origin <- factor(rep(most_common_origin, nrow(cost_grid)), levels = levels(df$origin))
  cost_grid$destination <- factor(rep(most_common_dest, nrow(cost_grid)), levels = levels(df$destination))
  
  cost_grid$pred_log_flow <- predict(model, newdata = cost_grid, re.form = NA)
  cost_grid$pred_flow <- pmax(exp(cost_grid$pred_log_flow) - 1, 0)
  
  df_fit <- model.frame(model)
  preds <- predict(model, newdata = df_fit, re.form = NULL)
  df_fit$pred_log_flow_full <- preds
  
  ss_total <- sum((df_fit$log_flow - mean(df_fit$log_flow))^2)
  ss_res <- sum((df_fit$log_flow - df_fit$pred_log_flow_full)^2)
  R2_full <- 1 - ss_res / ss_total
  
  res <- dplyr::select(as_tibble(cost_grid), cost, pred_flow) %>%
    dplyr::mutate(pred_flow = scales::rescale(pred_flow), r2 = R2_full)
  attr(res, "model") <- model
  attr(res, "df") <- df
  return(res)
}

run_spline_model_with_controls <- function(flow_df, census_df, dist_min = 0, dist_max = 25, spline_step = 5, pred_points = 200) {
  
  # Check columns
  check_required_cols(flow_df, c("cost", "origin", "destination", "flow"), "flow_df")
  check_required_cols(census_df, c("h3_index", "warehousing_emp_tot", "pop_total", "med_hh_income"), "census_df")
  
  flow_df <- flow_df %>%
    dplyr::left_join(
      census_df %>%
        dplyr::select(h3_index, warehousing_emp_tot) %>%
        dplyr::rename(origin = h3_index,
                      warehousing_emp_tot_origin = warehousing_emp_tot),
      by = "origin"
    ) %>%
    dplyr::left_join(
      census_df %>%
        dplyr::select(h3_index, pop_total, med_hh_income) %>%
        dplyr::rename(destination = h3_index,
                      pop_total_dest = pop_total,
                      med_hh_income_dest = med_hh_income),
      by = "destination"
    )
  
  # Debug prints (optional)
  cat("After joins:\n")
  cat("Total rows:", nrow(flow_df), "\n")
  cat("Missing warehousing_emp_tot_origin:", sum(is.na(flow_df$warehousing_emp_tot_origin)), "\n")
  cat("Missing pop_total_dest:", sum(is.na(flow_df$pop_total_dest)), "\n")
  cat("Missing med_hh_income_dest:", sum(is.na(flow_df$med_hh_income_dest)), "\n")
  cat("Unique origins:", length(unique(flow_df$origin)), "\n")
  cat("Unique destinations:", length(unique(flow_df$destination)), "\n")
  
  df <- flow_df %>%
    dplyr::filter(cost >= dist_min, cost <= dist_max) %>%
    dplyr::mutate(
      origin       = factor(origin),
      destination  = factor(destination),
      log_flow     = log(flow + 1),
      log_pop_total_dest = log(pop_total_dest + epsilon),
      log_med_hh_income_dest = log(med_hh_income_dest + epsilon),
      log_warehousing_emp_origin = log(warehousing_emp_tot_origin + epsilon)
    ) %>%
    # Remove rows with NA or infinite values after log transform
    filter(
      !is.na(log_pop_total_dest), !is.infinite(log_pop_total_dest),
      !is.na(log_med_hh_income_dest), !is.infinite(log_med_hh_income_dest),
      !is.na(log_warehousing_emp_origin), !is.infinite(log_warehousing_emp_origin)
    ) %>%
    droplevels()
  
  knots <- seq(dist_min + spline_step, dist_max - spline_step, by = spline_step)
  spline_basis <- splines::bs(df$cost,
                              knots = knots,
                              degree = 3,
                              intercept = FALSE,
                              Boundary.knots = c(dist_min, dist_max))
  colnames(spline_basis) <- paste0("spline_", seq_len(ncol(spline_basis)) - 1)
  df <- dplyr::bind_cols(df, as.data.frame(spline_basis))
  
  spline_terms <- paste0("spline_", seq_len(ncol(spline_basis)) - 1, collapse = " + ")
  
  formula_text <- paste0(
    "log_flow ~ 1 + ", spline_terms,
    " + log_pop_total_dest + log_med_hh_income_dest + log_warehousing_emp_origin + (1 | origin) + (1 | destination)"
  )
  
  model <- lme4::lmer(as.formula(formula_text), data = df)
  
  cost_grid <- data.frame(cost = seq(dist_min, dist_max, length.out = pred_points))
  spline_basis_grid <- splines::bs(cost_grid$cost,
                                   knots = knots,
                                   degree = 3,
                                   intercept = FALSE,
                                   Boundary.knots = c(dist_min, dist_max))
  colnames(spline_basis_grid) <- paste0("spline_", seq_len(ncol(spline_basis_grid)) - 1)
  cost_grid <- dplyr::bind_cols(cost_grid, as.data.frame(spline_basis_grid))
  
  mean_log_pop <- mean(df$log_pop_total_dest, na.rm = TRUE)
  mean_log_income <- mean(df$log_med_hh_income_dest, na.rm = TRUE)
  mean_log_warehousing <- mean(df$log_warehousing_emp_origin, na.rm = TRUE)
  
  cost_grid$log_pop_total_dest <- mean_log_pop
  cost_grid$log_med_hh_income_dest <- mean_log_income
  cost_grid$log_warehousing_emp_origin <- mean_log_warehousing
  
  most_common_origin <- names(which.max(table(df$origin)))
  most_common_dest <- names(which.max(table(df$destination)))
  
  cost_grid$origin <- factor(rep(most_common_origin, nrow(cost_grid)), levels = levels(df$origin))
  cost_grid$destination <- factor(rep(most_common_dest, nrow(cost_grid)), levels = levels(df$destination))
  
  cost_grid$pred_log_flow <- predict(model, newdata = cost_grid, re.form = NA)
  cost_grid$pred_flow <- pmax(exp(cost_grid$pred_log_flow) - 1, 0)
  
  df_fit <- model.frame(model)
  preds <- predict(model, newdata = df_fit, re.form = NULL)
  df_fit$pred_log_flow_full <- preds
  
  ss_total <- sum((df_fit$log_flow - mean(df_fit$log_flow))^2)
  ss_res <- sum((df_fit$log_flow - df_fit$pred_log_flow_full)^2)
  R2_full <- 1 - ss_res / ss_total
  
  res <- dplyr::select(as_tibble(cost_grid), cost, pred_flow) %>%
    dplyr::mutate(pred_flow = scales::rescale(pred_flow), r2 = R2_full)
  
  attr(res, "model") <- model
  attr(res, "df") <- df
  return(res)
}


compare_models <- function(flow_filepath, census_filepath, city_name, dist_max = 30, spline_step = 5) {
  flow_df <- read_csv(path.expand(flow_filepath))
  census_df <- read_csv(path.expand(census_filepath))
  
  baseline_res <- run_spline_model_baseline(flow_df, dist_max = dist_max, spline_step = spline_step)
  controls_res <- run_spline_model_with_controls(flow_df, census_df, dist_max = dist_max, spline_step = spline_step)
  
  baseline_model <- attr(baseline_res, "model")
  controls_model <- attr(controls_res, "model")
  
  varcomp_baseline <- as.data.frame(VarCorr(baseline_model))
  varcomp_controls <- as.data.frame(VarCorr(controls_model))
  
  origin_var_baseline <- varcomp_baseline %>% filter(grp == "origin" & var1 == "(Intercept)") %>% pull(vcov)
  dest_var_baseline <- varcomp_baseline %>% filter(grp == "destination" & var1 == "(Intercept)") %>% pull(vcov)
  
  origin_var_controls <- varcomp_controls %>% filter(grp == "origin" & var1 == "(Intercept)") %>% pull(vcov)
  dest_var_controls <- varcomp_controls %>% filter(grp == "destination" & var1 == "(Intercept)") %>% pull(vcov)
  
  origin_var_red <- 100 * (origin_var_baseline - origin_var_controls) / origin_var_baseline
  dest_var_red <- 100 * (dest_var_baseline - dest_var_controls) / dest_var_baseline
  
  r2_baseline <- unique(baseline_res$r2)
  r2_controls <- unique(controls_res$r2)
  r2_improve <- r2_controls - r2_baseline
  
  cat("Model Comparison for", city_name, "\n")
  cat(sprintf("R² Baseline: %.4f | R² With Controls: %.4f | Improvement: %.4f\n", r2_baseline, r2_controls, r2_improve))
  cat(sprintf("Origin RE variance reduction: %.2f%%\n", origin_var_red))
  cat(sprintf("Destination RE variance reduction: %.2f%%\n", dest_var_red))
  
  # Extract random effects (BLUPs)
  re_baseline <- ranef(baseline_model)
  re_controls <- ranef(controls_model)
  
  # Prepare data for plotting distributions
  df_re <- bind_rows(
    tibble(
      group = "origin",
      model = "Baseline",
      effect = re_baseline$origin[, "(Intercept)"]
    ),
    tibble(
      group = "origin",
      model = "With Controls",
      effect = re_controls$origin[, "(Intercept)"]
    ),
    tibble(
      group = "destination",
      model = "Baseline",
      effect = re_baseline$destination[, "(Intercept)"]
    ),
    tibble(
      group = "destination",
      model = "With Controls",
      effect = re_controls$destination[, "(Intercept)"]
    )
  )
  
  # Plot distribution of random effects
  re_plot <- ggplot(df_re, aes(x = effect, fill = model)) +
    geom_density(alpha = 0.4) +
    facet_grid(group ~ ., scales = "free") +
    labs(
      title = paste("Distribution of Random Effects (Intercepts) for", city_name),
      x = "Random Effect Value",
      y = "Density",
      fill = "Model"
    ) +
    theme_minimal()
  
  print(re_plot)
  
  baseline_res <- baseline_res %>% mutate(city = city_name, model = "Baseline")
  controls_res <- controls_res %>% mutate(city = city_name, model = "With Controls")
  
  return(bind_rows(baseline_res, controls_res))  # Explicit return added here
}

# --- Usage example ---

scheduled_files <- c(
  "Chicago" = "~/imperial/delivery-tunnels/flow_df_for_r_illinois.csv"
)
census_file <- "~/imperial/data/census_hex_chicago.csv"

comparison_curves <- imap_dfr(scheduled_files, ~ compare_models(.x, census_file, .y, spline_step = 5))

ggplot(comparison_curves, aes(x = cost, y = pred_flow, color = model)) +
  geom_line(size = 1) +
  facet_wrap(~city) +
  labs(
    title = "Scaled Predicted Flow vs Cost: Baseline vs Controls (Chicago)",
    x = "Distance (km)",
    y = "Predicted Flow (scaled 0–1)",
    color = "Model"
  ) +
  theme_minimal()





library(lme4)

flow_df <- read_csv("~/imperial/delivery-tunnels/flow_df_for_r_illinois.csv")
census_df <- read_csv("~/imperial/data/census_hex_chicago.csv")

library(dplyr)
library(lme4)
library(splines)

evaluate_controls_lrt_summary <- function(flow_df, census_df, dist_min = 0, dist_max = 25, spline_step = 5) {
  # Join controls including med_home_value
  flow_df <- flow_df %>%
    dplyr::left_join(
      census_df %>%
        dplyr::select(h3_index, warehousing_emp_tot) %>%
        dplyr::rename(origin = h3_index,
                      warehousing_emp_tot_origin = warehousing_emp_tot),
      by = "origin"
    ) %>%
    dplyr::left_join(
      census_df %>%
        dplyr::select(h3_index, pop_total, med_hh_income, med_home_value) %>%
        dplyr::rename(destination = h3_index,
                      pop_total_dest = pop_total,
                      med_hh_income_dest = med_hh_income,
                      med_home_value_dest = med_home_value),
      by = "destination"
    )
  
  # Filter and transform (exclude zeros and NAs before log)
  df <- flow_df %>%
    dplyr::filter(
      cost >= dist_min, cost <= dist_max,
      !is.na(pop_total_dest), pop_total_dest > 0,
      !is.na(med_hh_income_dest), med_hh_income_dest > 0,
      !is.na(med_home_value_dest), med_home_value_dest > 0,
      !is.na(warehousing_emp_tot_origin), warehousing_emp_tot_origin > 0,
      !is.na(flow)
    ) %>%
    dplyr::mutate(
      origin = factor(origin),
      destination = factor(destination),
      log_flow = log(flow + 1),
      log_pop_total_dest = log(pop_total_dest),
      log_med_hh_income_dest = log(med_hh_income_dest),
      log_med_home_value_dest = log(med_home_value_dest),
      log_warehousing_emp_origin = log(warehousing_emp_tot_origin)
    ) %>%
    droplevels()
  
  # Create spline basis
  knots <- seq(dist_min + spline_step, dist_max - spline_step, by = spline_step)
  spline_basis <- splines::bs(df$cost,
                              knots = knots,
                              degree = 3,
                              intercept = FALSE,
                              Boundary.knots = c(dist_min, dist_max))
  colnames(spline_basis) <- paste0("spline_", seq_len(ncol(spline_basis)) - 1)
  df <- dplyr::bind_cols(df, as.data.frame(spline_basis))
  
  spline_terms <- paste0("spline_", seq_len(ncol(spline_basis)) - 1, collapse = " + ")
  
  # Baseline model (no controls)
  formula_base <- paste0("log_flow ~ 1 + ", spline_terms, " + (1 | origin) + (1 | destination)")
  model_base <- lme4::lmer(as.formula(formula_base), data = df, REML = FALSE)
  
  # Controls to test including med_home_value
  controls <- list(
    pop = "log_pop_total_dest",
    income = "log_med_hh_income_dest",
    home_value = "log_med_home_value_dest",
    warehousing = "log_warehousing_emp_origin"
  )
  
  results_list <- list()
  for (ctrl_name in names(controls)) {
    ctrl_var <- controls[[ctrl_name]]
    formula_ctrl <- paste0("log_flow ~ 1 + ", spline_terms, " + ", ctrl_var, " + (1 | origin) + (1 | destination)")
    model_ctrl <- lme4::lmer(as.formula(formula_ctrl), data = df, REML = FALSE)
    
    lrt <- anova(model_base, model_ctrl)
    
    results_list[[ctrl_name]] <- tibble(
      control_variable = ctrl_name,
      chisq = lrt$Chisq[2],
      df = lrt$Df[2],
      p_value = lrt$`Pr(>Chisq)`[2]
    )
  }
  
  results_tbl <- dplyr::bind_rows(results_list)
  return(results_tbl)
}


results_table <- evaluate_controls_lrt_summary(flow_df, census_df, dist_min = 0, dist_max = 25, spline_step = 5)
print(results_table)



