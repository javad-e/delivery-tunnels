##########################################
# 1️⃣ Load libraries
##########################################
library(lme4)
library(readr)
library(splines)
library(ggplot2)
library(dplyr)

##########################################
# 2️⃣ Read talabat_sample.csv and filter to cost ≤ 15 km
##########################################
df <- readr::read_csv("~/imperial/delivery-tunnels/talabat_sample.csv") %>%
  dplyr::filter(cost <= 15) %>%
  dplyr::mutate(
    origin = paste0(vendor_lat, "_", vendor_lng),       # unique ID for origin
    destination = paste0(customer_lat, "_", customer_lng), # unique ID for destination
    origin = as.factor(origin),
    destination = as.factor(destination),
    log_flow = log(flow + 1)   # assuming `flow` column exists
  )

##########################################
# 3️⃣ Define spline knots every 2 km up to 15 km
##########################################
knots <- seq(2, 14, by = 2)  # knots inside boundary

spline_basis <- splines::bs(df$cost,
                            knots = knots,
                            degree = 3,
                            intercept = FALSE,
                            Boundary.knots = c(min(df$cost), 15))

colnames(spline_basis) <- paste0("spline_", 0:(ncol(spline_basis)-1))

df <- df %>%
  dplyr::select(-dplyr::starts_with("spline_")) %>%
  dplyr::bind_cols(as.data.frame(spline_basis))

##########################################
# 4️⃣ Fit mixed-effects model (NO extra controls)
##########################################
num_splines <- ncol(spline_basis)
spline_terms <- paste0("spline_", 0:(num_splines-1), collapse = " + ")

formula_text <- paste0(
  "log_flow ~ ", spline_terms, " + (1 | origin) + (1 | destination)"
)

model <- lme4::lmer(as.formula(formula_text), data = df)

summary(model)

##########################################
# 5️⃣ Create prediction grid restricted to ≤ 15 km
##########################################
cost_grid <- data.frame(cost = seq(min(df$cost), 15, length.out = 200))

spline_basis_grid <- splines::bs(cost_grid$cost,
                                 knots = knots,
                                 degree = 3,
                                 intercept = FALSE,
                                 Boundary.knots = c(min(df$cost), 15))

colnames(spline_basis_grid) <- paste0("spline_", 0:(ncol(spline_basis_grid)-1))

cost_grid <- dplyr::bind_cols(cost_grid, as.data.frame(spline_basis_grid))

##########################################
# 6️⃣ Predict fixed effects only
##########################################
cost_grid$pred_log_flow <- predict(model, newdata = cost_grid, re.form = NA)

# Back-transform to flow scale (subtract 1 to undo +1 offset)
cost_grid$pred_flow <- exp(cost_grid$pred_log_flow) - 1
cost_grid$pred_flow[cost_grid$pred_flow < 0] <- 0  # clamp negatives to 0 if any

##########################################
# 7️⃣ Plot results
##########################################
ggplot2::ggplot(df, ggplot2::aes(x = cost, y = log_flow)) +
  ggplot2::geom_point(alpha = 0.3) +
  ggplot2::geom_line(data = cost_grid, ggplot2::aes(x = cost, y = pred_log_flow),
                     color = "red", size = 1.2) +
  ggplot2::geom_vline(xintercept = knots, linetype = "dotted", color = "blue") +
  ggplot2::labs(title = "Spline Fit (cost ≤ 15 km)",
                x = "Cost (km)",
                y = "log(flow)") +
  ggplot2::theme_minimal()

ggplot2::ggplot(cost_grid, ggplot2::aes(x = cost, y = pred_flow)) +
  ggplot2::geom_line(color = "darkgreen", size = 1.2) +
  ggplot2::geom_vline(xintercept = knots, linetype = "dotted", color = "blue") +
  ggplot2::labs(title = "Predicted Flow vs Cost (cost ≤ 15 km, original flow scale)",
                x = "Cost (km)",
                y = "Predicted flow") +
  ggplot2::theme_minimal()
