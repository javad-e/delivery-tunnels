##########################################
# 1️⃣ Load libraries
##########################################
library(lme4)
library(readr)
library(splines)
library(dplyr)

##########################################
# 2️⃣ Read data and filter cost ≤ 10 km
##########################################
df <- readr::read_csv("~/imperial/delivery-tunnels/df_flow_for_r.csv") %>%
  dplyr::filter(cost <= 10) %>%
  dplyr::mutate(
    origin = as.factor(origin),
    destination = as.factor(destination),
    log_flow = log(flow + 1)  # add 1 to flow before log transform
  )

##########################################
# 3️⃣ Define spline knots every 2 km up to 10km
##########################################
knots <- seq(2, 14, by = 2)

spline_basis <- splines::bs(df$cost,
                            knots = knots,
                            degree = 3,
                            intercept = FALSE,
                            Boundary.knots = c(min(df$cost), 10))

colnames(spline_basis) <- paste0("spline_", 0:(ncol(spline_basis)-1))

df <- df %>%
  dplyr::select(-dplyr::starts_with("spline_")) %>%
  dplyr::bind_cols(as.data.frame(spline_basis))

num_splines <- ncol(spline_basis)
spline_terms <- paste0("spline_", 0:(num_splines-1), collapse = " + ")

##########################################
# 4️⃣ Read population and vendor data and merge (swapped controls)
##########################################
pop_data <- readr::read_csv("~/imperial/delivery-tunnels/dubai_pop_merged.csv")

# Vendors for origin (control for origin)
df <- df %>%
  dplyr::left_join(pop_data %>% dplyr::select(h3_index, num_vendors),
                   by = c("origin" = "h3_index")) %>%
  dplyr::rename(origin_vendors = num_vendors)

# Population for destination (control for destination)
df <- df %>%
  dplyr::left_join(pop_data %>% dplyr::select(h3_index, pop_share),
                   by = c("destination" = "h3_index")) %>%
  dplyr::rename(dest_pop = pop_share)

# Create logged controls
df <- df %>%
  dplyr::mutate(
    log_origin_vendors = log(origin_vendors + 1),
    log_dest_pop = log(dest_pop + 1)
  )

##########################################
# 5️⃣ Fit baseline model (no controls)
##########################################
model_base <- lme4::lmer(
  as.formula(paste0("log_flow ~ ", spline_terms, " + (1 | origin) + (1 | destination)")),
  data = df
)

##########################################
# 6️⃣ Fit full model (with controls)
##########################################
model_full <- lme4::lmer(
  as.formula(paste0(
    "log_flow ~ ", spline_terms, " + log_origin_vendors + log_dest_pop + (1 | origin) + (1 | destination)"
  )),
  data = df
)

##########################################
# 7️⃣ Extract and compare random intercept variances
##########################################
varcomp_base <- lme4::VarCorr(model_base)
varcomp_full <- lme4::VarCorr(model_full)

var_origin_base <- attr(varcomp_base$origin, "stddev")^2
var_origin_full <- attr(varcomp_full$origin, "stddev")^2

var_dest_base <- attr(varcomp_base$destination, "stddev")^2
var_dest_full <- attr(varcomp_full$destination, "stddev")^2

origin_var_reduction <- (var_origin_base - var_origin_full) / var_origin_base
dest_var_reduction <- (var_dest_base - var_dest_full) / var_dest_base

cat("Origin random intercept variance reduction:", round(origin_var_reduction * 100, 2), "%\n")
cat("Destination random intercept variance reduction:", round(dest_var_reduction * 100, 2), "%\n")

##########################################
# 8️⃣ Likelihood ratio test for model comparison
##########################################
anova_result <- anova(model_base, model_full)
print(anova_result)


library(ggplot2)
library(tibble)
library(dplyr)

# Extract random intercepts for origins
re_origin_base <- ranef(model_base)$origin %>%
  tibble::rownames_to_column("origin") %>%
  rename(random_intercept = `(Intercept)`) %>%
  mutate(model = "Baseline")

re_origin_full <- ranef(model_full)$origin %>%
  tibble::rownames_to_column("origin") %>%
  rename(random_intercept = `(Intercept)`) %>%
  mutate(model = "Full")

re_origin_all <- bind_rows(re_origin_base, re_origin_full)

# Plot density for origin random intercepts
ggplot(re_origin_all, aes(x = random_intercept, fill = model)) +
  geom_density(alpha = 0.4) +
  labs(title = "Density of Origin Random Intercepts",
       x = "Random Intercept",
       y = "Density",
       fill = "Model") +
  theme_minimal()

# Extract random intercepts for destinations
re_dest_base <- ranef(model_base)$destination %>%
  tibble::rownames_to_column("destination") %>%
  rename(random_intercept = `(Intercept)`) %>%
  mutate(model = "Baseline")

re_dest_full <- ranef(model_full)$destination %>%
  tibble::rownames_to_column("destination") %>%
  rename(random_intercept = `(Intercept)`) %>%
  mutate(model = "Full")

re_dest_all <- bind_rows(re_dest_base, re_dest_full)

# Plot density for destination random intercepts
ggplot(re_dest_all, aes(x = random_intercept, fill = model)) +
  geom_density(alpha = 0.4) +
  labs(title = "Density of Destination Random Intercepts",
       x = "Random Intercept",
       y = "Density",
       fill = "Model") +
  theme_minimal()

