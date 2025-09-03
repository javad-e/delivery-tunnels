library(tidyverse)
library(lme4)
library(splines)
library(scales)

epsilon <- 1e-6  # small constant to avoid log(0)

check_required_cols <- function(df, required_cols, df_name = "dataframe") {
  missing_cols <- setdiff(required_cols, colnames(df))
  if(length(missing_cols) > 0) {
    stop(paste0("Missing columns in ", df_name, ": ", paste(missing_cols, collapse = ", ")))
  }
}

run_spline_model_dubai <- function(flow_df, census_df = NULL, use_controls = FALSE, dist_min = 0, dist_max = 25, spline_step = 5, pred_points = 200) {
  
  # Check columns
  check_required_cols(flow_df, c("cost", "origin", "destination", "flow"), "flow_df")
  
  if(use_controls) {
    check_required_cols(census_df, c("h3_index", "num_vendors", "pop_share"), "census_df")
    # Join controls
    flow_df <- flow_df %>%
      left_join(
        census_df %>%
          dplyr::select(h3_index, num_vendors) %>%
          rename(origin = h3_index, num_vendors_origin = num_vendors),
        by = "origin"
      ) %>%
      left_join(
        census_df %>%
          dplyr::select(h3_index, pop_share) %>%
          rename(destination = h3_index, pop_share_dest = pop_share),
        by = "destination"
      )
  }
  
  df <- flow_df %>%
    filter(cost >= dist_min, cost <= dist_max) %>%
    mutate(
      origin = factor(origin),
      destination = factor(destination),
      log_flow = log(flow + 1)
    )
  
  if(use_controls) {
    df <- df %>%
      mutate(
        log_num_vendors_origin = log(num_vendors_origin + epsilon),
        log_pop_share_dest = log(pop_share_dest + epsilon)
      ) %>%
      filter(
        !is.na(log_num_vendors_origin), !is.infinite(log_num_vendors_origin),
        !is.na(log_pop_share_dest), !is.infinite(log_pop_share_dest)
      ) %>%
      droplevels()
  } else {
    df <- df %>% droplevels()
  }
  
  # Create spline basis
  knots <- seq(dist_min + spline_step, dist_max - spline_step, by = spline_step)
  spline_basis <- splines::bs(df$cost,
                              knots = knots,
                              degree = 3,
                              intercept = FALSE,
                              Boundary.knots = c(dist_min, dist_max))
  colnames(spline_basis) <- paste0("spline_", seq_len(ncol(spline_basis)) - 1)
  df <- bind_cols(df, as.data.frame(spline_basis))
  
  spline_terms <- paste0("spline_", seq_len(ncol(spline_basis)) - 1, collapse = " + ")
  
  # Build formula depending on controls
  if(use_controls) {
    formula_text <- paste0(
      "log_flow ~ 1 + ", spline_terms,
      " + log_num_vendors_origin + log_pop_share_dest + (1 | origin) + (1 | destination)"
    )
  } else {
    formula_text <- paste0(
      "log_flow ~ 1 + ", spline_terms,
      " + (1 | origin) + (1 | destination)"
    )
  }
  
  model <- lme4::lmer(as.formula(formula_text), data = df)
  
  # Prediction grid
  cost_grid <- data.frame(cost = seq(dist_min, dist_max, length.out = pred_points))
  spline_basis_grid <- splines::bs(cost_grid$cost,
                                   knots = knots,
                                   degree = 3,
                                   intercept = FALSE,
                                   Boundary.knots = c(dist_min, dist_max))
  colnames(spline_basis_grid) <- paste0("spline_", seq_len(ncol(spline_basis_grid)) - 1)
  cost_grid <- bind_cols(cost_grid, as.data.frame(spline_basis_grid))
  
  if(use_controls) {
    cost_grid$log_num_vendors_origin <- mean(df$log_num_vendors_origin, na.rm = TRUE)
    cost_grid$log_pop_share_dest <- mean(df$log_pop_share_dest, na.rm = TRUE)
  }
  
  most_common_origin <- levels(df$origin)[1]
  most_common_dest <- levels(df$destination)[1]
  
  cost_grid$origin <- factor(rep(most_common_origin, nrow(cost_grid)), levels = levels(df$origin))
  cost_grid$destination <- factor(rep(most_common_dest, nrow(cost_grid)), levels = levels(df$destination))
  
  cost_grid$pred_log_flow <- predict(model, newdata = cost_grid, re.form = NA)
  cost_grid$pred_flow <- pmax(exp(cost_grid$pred_log_flow) - 1, 0)
  
  # Calculate R² on training data
  df$pred_log_flow_full <- predict(model, re.form = NULL)
  ss_total <- sum((df$log_flow - mean(df$log_flow))^2)
  ss_res <- sum((df$log_flow - df$pred_log_flow_full)^2)
  R2_full <- 1 - ss_res / ss_total
  
  res <- cost_grid %>%
    dplyr::select(cost, pred_flow) %>%
    mutate(pred_flow = scales::rescale(pred_flow), r2 = R2_full)
  
  attr(res, "model") <- model
  attr(res, "df") <- df
  return(res)
}

compare_models_dubai <- function(flow_filepath, census_filepath, city_name, dist_max = 30, spline_step = 5) {
  flow_df <- read_csv(path.expand(flow_filepath), show_col_types = FALSE)
  census_df <- read_csv(path.expand(census_filepath), show_col_types = FALSE)
  
  # Model WITHOUT controls
  res_no_ctrl <- run_spline_model_dubai(flow_df, NULL, use_controls = FALSE, dist_max = dist_max, spline_step = spline_step)
  model_no_ctrl <- attr(res_no_ctrl, "model")
  varcomp_no_ctrl <- as.data.frame(VarCorr(model_no_ctrl))
  origin_var_no_ctrl <- varcomp_no_ctrl %>% filter(grp == "origin" & var1 == "(Intercept)") %>% pull(vcov)
  dest_var_no_ctrl <- varcomp_no_ctrl %>% filter(grp == "destination" & var1 == "(Intercept)") %>% pull(vcov)
  
  # Model WITH controls
  res_ctrl <- run_spline_model_dubai(flow_df, census_df, use_controls = TRUE, dist_max = dist_max, spline_step = spline_step)
  model_ctrl <- attr(res_ctrl, "model")
  varcomp_ctrl <- as.data.frame(VarCorr(model_ctrl))
  origin_var_ctrl <- varcomp_ctrl %>% filter(grp == "origin" & var1 == "(Intercept)") %>% pull(vcov)
  dest_var_ctrl <- varcomp_ctrl %>% filter(grp == "destination" & var1 == "(Intercept)") %>% pull(vcov)
  
  cat("Model Comparison for", city_name, "\n")
  cat(sprintf("R² without controls: %.4f\n", unique(res_no_ctrl$r2)))
  cat(sprintf("Origin RE variance without controls: %.6f\n", origin_var_no_ctrl))
  cat(sprintf("Destination RE variance without controls: %.6f\n\n", dest_var_no_ctrl))
  
  cat(sprintf("R² with controls: %.4f\n", unique(res_ctrl$r2)))
  cat(sprintf("Origin RE variance with controls: %.6f\n", origin_var_ctrl))
  cat(sprintf("Destination RE variance with controls: %.6f\n\n", dest_var_ctrl))
  
  # Plot random effects for model with controls
  re <- ranef(model_ctrl)
  df_re <- bind_rows(
    tibble(group = "origin", effect = re$origin[, "(Intercept)"]),
    tibble(group = "destination", effect = re$destination[, "(Intercept)"])
  )
  
  library(ggplot2)
  p <- ggplot(df_re, aes(x = effect)) +
    geom_density(fill = "lightblue", alpha = 0.6) +
    facet_wrap(~group, scales = "free") +
    labs(title = paste("Random Effects Distribution for", city_name, "(With Controls)"),
         x = "Random Effect", y = "Density") +
    theme_minimal()
  print(p)
  
  list(
    without_controls = res_no_ctrl %>% mutate(city = city_name),
    with_controls = res_ctrl %>% mutate(city = city_name)
  )
}

# --- Usage example ---

dubai_flow_file <- "~/imperial/delivery-tunnels/flow_df_for_r_dubai.csv"
dubai_census_file <- "~/imperial/data/dubai_pop_merged.csv"

dubai_results <- compare_models_dubai(dubai_flow_file, dubai_census_file, "Dubai")
