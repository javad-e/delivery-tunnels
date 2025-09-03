# ==========================================================
# Libraries
# ==========================================================
library(tidyverse)
library(lme4)
library(scales)
library(dplyr)
library(tidyr)
library(readr)

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
# Mixed-effects model with controls
# ==========================================================
run_controls_only_model <- function(flow_df_path, census_df_path, dataset_name = "Dataset") {
  
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
      log_pop_total_dest = log(pop_total_dest + epsilon),   # destination population
      log_warehousing_share_origin = log(warehousing_share_origin + epsilon),
      log_emp_share_origin = log(emp_share_origin + epsilon),
      log_med_gross_rent_origin = log(med_gross_rent_origin + epsilon),
      log_med_hh_income_origin = log(med_hh_income_origin + epsilon),
      log_flow = log(flow + epsilon)
    )
  
  # Model formula including new variables
  model_formula <- as.formula(
    "log_flow ~ log_num_warehouses + log_pop_total_dest + 
     log_warehousing_share_origin + log_emp_share_origin +
     log_med_gross_rent_origin + log_med_hh_income_origin +
     (1 | origin) + (1 | destination)"
  )
  
  model <- lmer(model_formula, data = flow_df, REML = FALSE)
  
  # Calculate R² manually
  fitted_vals <- fitted(model)
  r2 <- 1 - sum((flow_df$log_flow - fitted_vals)^2) / sum((flow_df$log_flow - mean(flow_df$log_flow))^2)
  
  # Summarise coefficients
  summary_df <- summary(model)$coefficients %>%
    as.data.frame() %>%
    rownames_to_column("term")
  
  # Extract random effects std. dev.
  rand_eff <- as.data.frame(VarCorr(model))
  sigma_u <- rand_eff$sdcor[rand_eff$grp == "origin" & rand_eff$var1 == "(Intercept)"]
  sigma_v <- rand_eff$sdcor[rand_eff$grp == "destination" & rand_eff$var1 == "(Intercept)"]
  
  cat("\n=== Controls-Only Model Summary ===\n")
  cat("Dataset:", dataset_name, "\n")
  cat("R²:", round(r2, 4), "\n")
  print(summary_df)
  
  return(list(
    model = model,
    data = flow_df,
    summary = summary_df,
    r2 = r2,
    sigma_u = sigma_u,
    sigma_v = sigma_v
  ))
}

# ==========================================================
# Usage example for Chicago
# ==========================================================
chicago_flow_file <- "~/imperial/delivery-tunnels/flow_df_for_r_chicago.csv"
chicago_census_file <- "~/imperial/data/census_hex_chicago.csv"
chicago_results <- run_controls_only_model(chicago_flow_file, chicago_census_file, "Chicago")

# ==========================================================
# NYC OD Flow Prediction
# ==========================================================
empty_od_file <- "~/imperial/data/empty_od.csv"
census_file   <- "~/imperial/data/census_hex_nyc.csv"
T_total       <- 839500000
n_sims        <- 1000

# Estimated model parameters
beta0 <- 1.9864 
beta_vendors <- 0
beta_pop <- 0.3564
beta_warehousing_share <- 1.0409
beta_emp <- -1.6345
sigma_u <- sqrt(0)
sigma_v <- sqrt(0.3564)

# Load OD and Census data
od_df <- readr::read_csv(empty_od_file) %>%
  dplyr::rename(origin_h3 = origin, dest_h3 = destination)
census_df <- readr::read_csv(census_file)

# Prepare socioeconomic variables
relevant_h3 <- unique(c(od_df$origin_h3, od_df$dest_h3))
census_df <- census_df %>%
  dplyr::filter(h3_index %in% relevant_h3) %>%
  dplyr::mutate(
    log_pop_total         = log(pop_total + 1e-9),          # destination population
    log_warehousing_share = log(warehousing_share + 1e-9),
    log_emp_share         = log(emp_share + 1e-9),
    log_med_gross_rent    = log(med_gross_rent + 1e-9),
    log_med_hh_income     = log(med_hh_income + 1e-9)
  ) %>%
  dplyr::select(h3_index, log_pop_total, log_warehousing_share, log_emp_share,
                log_med_gross_rent, log_med_hh_income)

# Join origin/destination variables
od_df <- od_df %>%
  dplyr::left_join(
    census_df %>%
      dplyr::rename(origin_h3 = h3_index,
                    log_warehousing_share_origin = log_warehousing_share,
                    log_emp_share_origin         = log_emp_share,
                    log_med_gross_rent_origin    = log_med_gross_rent,
                    log_med_hh_income_origin     = log_med_hh_income),
    by = "origin_h3"
  ) %>%
  dplyr::left_join(
    census_df %>%
      dplyr::rename(dest_h3 = h3_index,
                    log_pop_total_dest = log_pop_total),  # destination population
    by = "dest_h3"
  )

# Impute missing and construct LP_base
od_df <- od_df %>%
  dplyr::mutate(
    log_warehousing_share_origin = ifelse(is.na(log_warehousing_share_origin), 0, log_warehousing_share_origin),
    log_emp_share_origin         = ifelse(is.na(log_emp_share_origin), 0, log_emp_share_origin),
    log_med_gross_rent_origin    = ifelse(is.na(log_med_gross_rent_origin), 0, log_med_gross_rent_origin),
    log_med_hh_income_origin     = ifelse(is.na(log_med_hh_income_origin), 0, log_med_hh_income_origin),
    log_pop_total_dest           = ifelse(is.na(log_pop_total_dest), 0, log_pop_total_dest),
    LP_base = beta0 +
      beta_pop * log_pop_total_dest +
      beta_warehousing_share * log_warehousing_share_origin +
      beta_emp * log_emp_share_origin +
      beta_vendors * 0 +
      0.1 * log_med_gross_rent_origin +
      0.1 * log_med_hh_income_origin
  )

# Prepare indices
origins <- sort(unique(od_df$origin_h3))
dests   <- sort(unique(od_df$dest_h3))
n_orig <- length(origins)
n_dest <- length(dests)

origin_map <- setNames(seq_along(origins), origins)
dest_map   <- setNames(seq_along(dests), dests)

LP_mat <- matrix(NA_real_, nrow = n_orig, ncol = n_dest)
for(i in seq_len(nrow(od_df))) {
  ori_idx  <- origin_map[[od_df$origin_h3[i]]]
  dest_idx <- dest_map[[od_df$dest_h3[i]]]
  LP_mat[ori_idx, dest_idx] <- od_df$LP_base[i]
}

# Monte Carlo simulations with REs
set.seed(1234)
preds_array <- array(NA_real_, dim = c(n_orig, n_dest, n_sims))
for(s in seq_len(n_sims)) {
  u <- rnorm(n_orig, mean = 0, sd = sigma_u)
  v <- rnorm(n_dest, mean = 0, sd = sigma_v)
  
  LP_sim <- LP_mat + outer(u, rep(1, n_dest)) + outer(rep(1, n_orig), v)
  pred_sim <- exp(LP_sim)
  
  scale_factor <- T_total / sum(pred_sim, na.rm = TRUE)
  pred_sim <- pred_sim * scale_factor
  preds_array[,,s] <- pred_sim
  
  if(s %% 100 == 0) cat("Completed simulation", s, "\n")
}

# Compute mean, SD, and quantiles per OD pair
mean_pred <- apply(preds_array, c(1,2), mean, na.rm = TRUE)
sd_pred   <- apply(preds_array, c(1,2), sd, na.rm = TRUE)
q_low     <- apply(preds_array, c(1,2), quantile, probs = 0.025, na.rm = TRUE)
q_high    <- apply(preds_array, c(1,2), quantile, probs = 0.975, na.rm = TRUE)

# Map results back
od_df$mean_flow <- NA
od_df$sd_flow   <- NA
od_df$q025      <- NA
od_df$q975      <- NA

for(i in seq_len(nrow(od_df))) {
  ori_idx  <- origin_map[[od_df$origin_h3[i]]]
  dest_idx <- dest_map[[od_df$dest_h3[i]]]
  od_df$mean_flow[i] <- mean_pred[ori_idx, dest_idx]
  od_df$sd_flow[i]   <- sd_pred[ori_idx, dest_idx]
  od_df$q025[i]      <- q_low[ori_idx, dest_idx]
  od_df$q975[i]      <- q_high[ori_idx, dest_idx]
}

# Diagnostics
cat("Check total flow (should equal T_total):", sum(od_df$mean_flow, na.rm = TRUE), "\n")
top_flows <- od_df %>%
  arrange(desc(mean_flow)) %>%
  head(50) %>%
  dplyr::select(origin_h3, dest_h3, mean_flow, sd_flow, q025, q975)
print(top_flows)

# Save results
write_csv(od_df, "~/imperial/data/nyc_od_predicted_scheduled_flows.csv")
