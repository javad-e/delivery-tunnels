library(tidyverse)
library(lme4)
library(splines)
library(merTools)
library(scales) 

run_spline_model_return_curve <- function(df, dist_min = 0, dist_max = 25, spline_step = 5, pred_points = 200) {
  df <- df %>%
    dplyr::filter(cost >= dist_min, cost <= dist_max) %>%
    dplyr::mutate(
      origin       = factor(origin),
      destination  = factor(destination),
      log_flow     = log(flow + 1)
    )
  
  knots <- seq(dist_min + spline_step, dist_max - spline_step, by = spline_step)
  spline_basis <- splines::bs(df$cost,
                              knots = knots,
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
  cost_grid$origin <- factor(rep(levels(df$origin)[1], nrow(cost_grid)), levels = levels(df$origin))
  cost_grid$destination <- factor(rep(levels(df$destination)[1], nrow(cost_grid)), levels = levels(df$destination))
  
  cost_grid$pred_log_flow <- predict(model, newdata = cost_grid, re.form = NA)
  cost_grid$pred_flow     <- pmax(exp(cost_grid$pred_log_flow) - 1, 0)
  
  # R²
  df$pred_log_flow_full <- predict(model, re.form = NULL)
  ss_total <- sum((df$log_flow - mean(df$log_flow))^2)
  ss_res   <- sum((df$log_flow - df$pred_log_flow_full)^2)
  R2_full  <- 1 - ss_res / ss_total
  
  return(dplyr::select(as_tibble(cost_grid), cost, pred_flow) %>%
           dplyr::mutate(pred_flow = scales::rescale(pred_flow), r2 = R2_full))
}




process_dataset <- function(filepath, city_name, dist_max = 30, spline_step = 5) {
  df <- read_csv(filepath)
  run_spline_model_return_curve(df, dist_max = dist_max, spline_step = spline_step) %>%
    mutate(city = city_name)
}

# --- Express Orders ---
express_files <- c(
  "Dubai"     = "~/imperial/delivery-tunnels/flow_df_for_r_dubai.csv",
  "Shanghai"  = "~/imperial/delivery-tunnels/flow_df_for_r_lade_sh.csv",
  "Chongqing" = "~/imperial/delivery-tunnels/flow_df_for_r_lade_cq.csv",
  "Hangzhou"  = "~/imperial/delivery-tunnels/flow_df_for_r_lade_hz.csv"
)

express_curves <- imap_dfr(express_files, process_dataset)

# Plot Express Orders
ggplot(express_curves, aes(x = cost, y = pred_flow, color = city)) +
  geom_line(size = 1) +
  labs(
    title = "Scaled Predicted Flow vs Cost (Express Orders)",
    x = "Distance (km)",
    y = "Predicted Flow (scaled 0–1)",
    color = "City"
  ) +
  theme_minimal()

# Print R2 for Express Orders
express_r2 <- express_curves %>%
  group_by(city) %>%
  summarise(R2 = unique(r2))
express_r2_string <- paste0(express_r2$city, ": ", round(express_r2$R2, 3), collapse = " | ")
cat("Model R² for Express Orders:\n", express_r2_string, "\n")

# Calculate combined mean curve + CI ribbon across cities
combined_ci_df <- express_curves %>%
  group_by(cost) %>%
  summarise(
    mean_flow = mean(pred_flow),
    sd_flow = sd(pred_flow),
    n = n(),
    se_flow = sd_flow / sqrt(n),
    lower = mean_flow - 1.96 * se_flow,
    upper = mean_flow + 1.96 * se_flow
  )

# Plot individual city curves + combined CI ribbon
ggplot() +
  geom_line(data = express_curves, aes(x = cost, y = pred_flow, color = city), size = 1) +
  geom_ribbon(data = combined_ci_df, aes(x = cost, ymin = lower, ymax = upper), fill = "grey70", alpha = 0.3) +
  labs(
    title = "Predicted Flow Curves with Combined 95% CI Ribbon",
    x = "Distance (km)",
    y = "Predicted Flow (scaled 0–1)",
    color = "City"
  ) +
  theme_minimal()



# --- Scheduled Orders ---
scheduled_files <- c(
  "Chicago"     = "~/imperial/delivery-tunnels/flow_df_for_r_illinois.csv",
  "Boston"      = "~/imperial/delivery-tunnels/flow_df_for_r_massachusetts.csv",
  "Los Angeles" = "~/imperial/delivery-tunnels/flow_df_for_r_california.csv"
)

scheduled_curves <- imap_dfr(scheduled_files, process_dataset, spline_step = 5)

# Plot Scheduled Orders
ggplot(scheduled_curves, aes(x = cost, y = pred_flow, color = city)) +
  geom_line(size = 1) +
  labs(
    title = "Scaled Predicted Flow vs Cost (Scheduled Orders)",
    x = "Distance (km)",
    y = "Predicted Flow (scaled 0–1)",
    color = "City"
  ) +
  theme_minimal()

# Print R2 for Scheduled Orders
scheduled_r2 <- scheduled_curves %>%
  group_by(city) %>%
  summarise(R2 = unique(r2))
scheduled_r2_string <- paste0(scheduled_r2$city, ": ", round(scheduled_r2$R2, 3), collapse = " | ")
cat("Model R² for Scheduled Orders:\n", scheduled_r2_string, "\n")
