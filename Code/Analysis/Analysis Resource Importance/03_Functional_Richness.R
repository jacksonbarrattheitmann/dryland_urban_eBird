library(tidyverse)
library(glmmTMB)

# Settings

min_checklists <- 10
min_localities <- 5
water_min_range <- 1
refit_models <- FALSE
habitat_step <- 0.1
n_aridity_points <- 200

input_file <- fs::path(
  "Outputs", "Analysis_Resource_Importance", "01_Data_Setup", "Data",
  "Resource_Importance_Analysis_Base.rds"
)

output_root <- fs::path(
  "Outputs", "Analysis_Resource_Importance", "03_Functional_Richness"
)

data_dir <- fs::path(output_root, "Data")
model_dir <- fs::path(output_root, "Models")
table_dir <- fs::path(output_root, "Table")
figure_dir <- fs::path(output_root, "Figure")

walk(
  c(data_dir, model_dir, table_dir, figure_dir),
  fs::dir_create
)

# Functional richness sample

analysis_base <- read_rds(input_file)

model_vars <- c(
  "checklist_id", "FRic", "sr", "avg_aridity_index",
  "water_pp", "shrub_pp", "built_pp",
  "duration_10", "effort_distance_km", "protocol_name", "NAME10",
  "locality_id", "state_name"
)

missing_vars <- setdiff(model_vars, names(analysis_base))
if (length(missing_vars) > 0) {
  stop("Missing variables: ", paste(missing_vars, collapse = ", "))
}

if (anyDuplicated(analysis_base$checklist_id) > 0) {
  stop("Duplicate checklist IDs remain after data setup.")
}

fric_pool <- analysis_base %>%
  filter(if_all(all_of(model_vars), ~ !is.na(.x)))

rm(analysis_base)
gc()

numeric_vars <- c(
  "FRic", "sr", "avg_aridity_index", "water_pp", "shrub_pp", "built_pp",
  "duration_10", "effort_distance_km"
)

finite_check <- fric_pool %>%
  summarise(across(all_of(numeric_vars), ~ all(is.finite(.x)))) %>%
  unlist(use.names = FALSE)

if (!all(finite_check)) {
  stop("Non-finite values remain in the functional-richness analysis pool.")
}

eligible_localities <- fric_pool %>%
  count(NAME10, locality_id, name = "n_checklists") %>%
  filter(n_checklists >= min_checklists)

eligible_cities <- eligible_localities %>%
  count(NAME10, name = "n_localities") %>%
  filter(n_localities >= min_localities)

fric_data <- fric_pool %>%
  semi_join(eligible_localities, by = c("NAME10", "locality_id")) %>%
  semi_join(eligible_cities, by = "NAME10")

rm(fric_pool, eligible_localities, eligible_cities)
gc()

if (any(fric_data$FRic <= 0 | fric_data$FRic >= 1)) {
  stop("Final FRic sample must be strictly inside (0,1) for beta regression.")
}

# Display baseline and final habitat centering

city_aridity <- fric_data %>%
  distinct(NAME10, avg_aridity_index)

if (nrow(city_aridity) != n_distinct(fric_data$NAME10)) {
  stop("At least one retained city has more than one aridity value.")
}

mean_city_aridity <- mean(city_aridity$avg_aridity_index)

city_baseline <- fric_data %>%
  group_by(NAME10) %>%
  summarise(
    city_mean_FRic_checklist = mean(FRic),
    .groups = "drop"
  )

if (any(city_baseline$city_mean_FRic_checklist <= 0 |
        city_baseline$city_mean_FRic_checklist >= 1)) {
  stop("At least one retained city has an invalid FRic baseline.")
}

typical_city_FRic_baseline <- median(city_baseline$city_mean_FRic_checklist)

locality_habitat <- fric_data %>%
  distinct(NAME10, locality_id, water_pp, shrub_pp, built_pp)

city_habitat <- locality_habitat %>%
  group_by(NAME10) %>%
  summarise(
    city_mean_water_pp = mean(water_pp),
    city_mean_shrub_pp = mean(shrub_pp),
    city_mean_built_pp = mean(built_pp),
    city_range_water_pp = max(water_pp) - min(water_pp),
    city_range_shrub_pp = max(shrub_pp) - min(shrub_pp),
    city_range_built_pp = max(built_pp) - min(built_pp),
    .groups = "drop"
  )

mean_city_water <- mean(city_habitat$city_mean_water_pp)
mean_city_shrub <- mean(city_habitat$city_mean_shrub_pp)
mean_city_built <- mean(city_habitat$city_mean_built_pp)

fric_analysis_data <- fric_data %>%
  select(-any_of(c(
    "city_mean_water_pp", "city_mean_shrub_pp", "city_mean_built_pp"
  ))) %>%
  left_join(
    city_habitat %>%
      select(NAME10, city_mean_water_pp, city_mean_shrub_pp, city_mean_built_pp),
    by = "NAME10"
  ) %>%
  mutate(
    aridity_centered = avg_aridity_index - mean_city_aridity,
    water_within_model_10 = (water_pp - city_mean_water_pp) / 10,
    shrub_within_model_10 = (shrub_pp - city_mean_shrub_pp) / 10,
    built_within_model_10 = (built_pp - city_mean_built_pp) / 10,
    city_context_water_10 = (city_mean_water_pp - mean_city_water) / 10,
    city_context_shrub_10 = (city_mean_shrub_pp - mean_city_shrub) / 10,
    city_context_built_10 = (city_mean_built_pp - mean_city_built) / 10,
    NAME10 = factor(NAME10),
    locality_id = factor(locality_id),
    protocol_name = factor(protocol_name, levels = c("Stationary", "Traveling"))
  )

rm(fric_data, locality_habitat, city_baseline)
gc()

centering_qc <- fric_analysis_data %>%
  distinct(
    NAME10, locality_id,
    water_within_model_10,
    shrub_within_model_10,
    built_within_model_10
  ) %>%
  group_by(NAME10) %>%
  summarise(
    across(
      c(water_within_model_10, shrub_within_model_10, built_within_model_10),
      mean
    ),
    .groups = "drop"
  )

if (any(abs(as.matrix(centering_qc %>% select(-NAME10))) > 1e-10)) {
  stop("Final within-city habitat centering failed.")
}

write_rds(
  fric_analysis_data,
  fs::path(data_dir, "Functional_Richness_Resource_Importance_Data.rds")
)

sample_summary <- tibble(
  n_checklists = nrow(fric_analysis_data),
  n_localities = n_distinct(fric_analysis_data$locality_id),
  n_cities = n_distinct(fric_analysis_data$NAME10),
  typical_city_FRic_baseline = typical_city_FRic_baseline,
  mean_city_aridity = mean_city_aridity,
  min_aridity = min(city_aridity$avg_aridity_index),
  max_aridity = max(city_aridity$avg_aridity_index)
)

# Model setup

habitat_info <- tribble(
  ~habitat,          ~stub,              ~within_term,               ~context_term,
  "Water",           "Water",           "water_within_model_10",   "city_context_water_10",
  "Shrub and scrub", "Shrub_and_scrub", "shrub_within_model_10",   "city_context_shrub_10",
  "Built",           "Built",           "built_within_model_10",   "city_context_built_10"
)

make_formula <- function(within_term, context_term) {
  as.formula(
    paste0(
      "FRic ~ ", within_term, " * aridity_centered + ", context_term,
      " + duration_10 + effort_distance_km + protocol_name + ",
      "(1 + ", within_term, " || NAME10) + (1 | locality_id)"
    )
  )
}

model_ok <- function(fit) {
  !inherits(fit, "error") &&
    isTRUE(fit$fit$convergence == 0) &&
    isTRUE(fit$sdr$pdHess)
}

model_matches <- function(fit, model_formula, n) {
  if (!model_ok(fit) || nobs(fit) != n) return(FALSE)
  
  saved_formula <- paste(deparse(stats::formula(fit)), collapse = "")
  expected_formula <- paste(deparse(model_formula), collapse = "")
  
  identical(gsub("\\s+", "", saved_formula), gsub("\\s+", "", expected_formula))
}

fit_model <- function(data, within_term, context_term, model_path) {
  model_formula <- make_formula(within_term, context_term)
  
  if (!refit_models && fs::file_exists(model_path)) {
    fit <- tryCatch(read_rds(model_path), error = function(e) e)
    if (model_matches(fit, model_formula, nrow(data))) return(fit)
    rm(fit)
    gc()
  }
  
  fit <- tryCatch(
    glmmTMB(
      model_formula,
      data = data,
      family = beta_family(link = "logit"),
      control = glmmTMBControl(
        optCtrl = list(iter.max = 1000, eval.max = 1000)
      )
    ),
    error = function(e) e
  )
  
  if (!model_ok(fit)) {
    rm(fit)
    gc()
    
    fit <- tryCatch(
      glmmTMB(
        model_formula,
        data = data,
        family = beta_family(link = "logit"),
        control = glmmTMBControl(
          optimizer = optim,
          optArgs = list(method = "BFGS")
        )
      ),
      error = function(e) e
    )
  }
  
  if (!model_ok(fit)) {
    stop("Model failed to converge: ", fs::path_file(model_path))
  }
  
  write_rds(fit, model_path)
  fit
}

interaction_term <- function(fit, within_term) {
  candidates <- c(
    paste0(within_term, ":aridity_centered"),
    paste0("aridity_centered:", within_term)
  )
  
  term <- intersect(candidates, names(fixef(fit)$cond))
  if (length(term) != 1) stop("Could not identify Habitat x Aridity term.")
  term
}

extract_interaction <- function(fit, habitat, within_term, model = "Primary Beta") {
  coef_table <- summary(fit)$coefficients$cond
  term <- interaction_term(fit, within_term)
  x <- coef_table[term, , drop = FALSE]
  
  tibble(
    habitat = habitat,
    model = model,
    Estimate = x[1, "Estimate"],
    SE = x[1, "Std. Error"],
    z = x[1, "z value"],
    p = x[1, "Pr(>|z|)"],
    CI_lower = Estimate - 1.96 * SE,
    CI_upper = Estimate + 1.96 * SE,
    Estimate_1pp = Estimate / 10,
    SE_1pp = SE / 10,
    CI_lower_1pp = CI_lower / 10,
    CI_upper_1pp = CI_upper / 10
  )
}

extract_qc <- function(fit, habitat, within_term) {
  beta <- fixef(fit)$cond
  pearson <- residuals(fit, type = "pearson")
  pearson_ratio <- sum(pearson^2, na.rm = TRUE) / max(1, nobs(fit) - length(beta))
  
  random_sd <- attr(VarCorr(fit)$cond$NAME10, "stddev")
  slope_sd <- unname(random_sd[within_term])
  
  vc <- vcov(fit)$cond
  cor_mat <- cov2cor(vc)
  diag(cor_mat) <- NA_real_
  
  tibble(
    habitat = habitat,
    n_checklists = nobs(fit),
    convergence = fit$fit$convergence,
    pdHess = isTRUE(fit$sdr$pdHess),
    beta_precision = as.numeric(sigma(fit)),
    pearson_dispersion = pearson_ratio,
    city_random_slope_SD = slope_sd,
    max_fixed_effect_correlation = max(abs(cor_mat), na.rm = TRUE)
  )
}

make_curve <- function(fit, habitat, within_term, aridity_values) {
  beta <- fixef(fit)$cond
  vc <- vcov(fit)$cond
  term <- interaction_term(fit, within_term)
  centered <- aridity_values - mean_city_aridity
  
  delta <- habitat_step * (beta[within_term] + beta[term] * centered)
  delta_var <- habitat_step^2 * (
    vc[within_term, within_term] +
      centered^2 * vc[term, term] +
      2 * centered * vc[within_term, term]
  )
  delta_se <- sqrt(pmax(delta_var, 0))
  
  baseline_logit <- qlogis(typical_city_FRic_baseline)
  fric_after <- plogis(baseline_logit + delta)
  fric_lower <- plogis(baseline_logit + delta - 1.96 * delta_se)
  fric_upper <- plogis(baseline_logit + delta + 1.96 * delta_se)
  
  tibble(
    habitat = habitat,
    avg_aridity_index = aridity_values,
    aridity_centered = centered,
    typical_city_FRic_baseline = typical_city_FRic_baseline,
    delta_logit_FRic = delta,
    delta_logit_FRic_SE = delta_se,
    odds_ratio = exp(delta),
    odds_ratio_CI_lower = exp(delta - 1.96 * delta_se),
    odds_ratio_CI_upper = exp(delta + 1.96 * delta_se),
    relative_change_percent_at_typical_baseline =
      100 * (fric_after / typical_city_FRic_baseline - 1),
    relative_CI_lower_percent_at_typical_baseline =
      100 * (fric_lower / typical_city_FRic_baseline - 1),
    relative_CI_upper_percent_at_typical_baseline =
      100 * (fric_upper / typical_city_FRic_baseline - 1)
  )
}

make_city_effects <- function(fit, habitat, within_term) {
  beta <- fixef(fit)$cond
  term <- interaction_term(fit, within_term)
  city_re <- ranef(fit)$cond$NAME10
  
  if (is.null(city_re) || !within_term %in% colnames(city_re)) {
    stop("Could not recover city random slopes for ", habitat, ".")
  }
  
  tibble(
    NAME10 = rownames(city_re),
    random_slope = as.numeric(city_re[, within_term])
  ) %>%
    left_join(
      city_aridity %>% mutate(NAME10 = as.character(NAME10)),
      by = "NAME10"
    ) %>%
    mutate(
      habitat = habitat,
      aridity_centered = avg_aridity_index - mean_city_aridity,
      fixed_slope = beta[within_term] + beta[term] * aridity_centered,
      city_slope = fixed_slope + random_slope,
      delta_logit = habitat_step * city_slope,
      relative_change_percent = 100 * (
        plogis(qlogis(typical_city_FRic_baseline) + delta_logit) /
          typical_city_FRic_baseline - 1
      )
    )
}

# Primary models

model_data <- fric_analysis_data %>%
  select(
    FRic,
    water_within_model_10, shrub_within_model_10, built_within_model_10,
    aridity_centered,
    city_context_water_10, city_context_shrub_10, city_context_built_10,
    duration_10, effort_distance_km, protocol_name,
    NAME10, locality_id
  )

rm(fric_analysis_data)
gc()

aridity_values <- seq(
  min(city_aridity$avg_aridity_index),
  max(city_aridity$avg_aridity_index),
  length.out = n_aridity_points
)

process_habitat <- function(habitat, stub, within_term, context_term) {
  model_path <- fs::path(model_dir, paste0(stub, "_model.rds"))
  fit <- fit_model(model_data, within_term, context_term, model_path)
  
  out <- list(
    model_path = model_path,
    interaction = extract_interaction(fit, habitat, within_term),
    qc = extract_qc(fit, habitat, within_term),
    curve = make_curve(fit, habitat, within_term, aridity_values),
    city_effects = make_city_effects(fit, habitat, within_term)
  )
  
  rm(fit)
  gc()
  
  out
}

primary_outputs <- pmap(habitat_info, process_habitat)

interaction_results <- map_dfr(primary_outputs, "interaction")
model_qc <- map_dfr(primary_outputs, "qc")
resource_curves <- map_dfr(primary_outputs, "curve")
city_effects <- map_dfr(primary_outputs, "city_effects")

model_registry <- habitat_info %>%
  mutate(model_path = map_chr(primary_outputs, "model_path"))

# Water support sensitivity

water_support <- city_habitat %>%
  transmute(
    NAME10,
    water_range_pp = city_range_water_pp,
    retained = water_range_pp >= water_min_range
  )

water_cities <- water_support %>%
  filter(retained) %>%
  pull(NAME10)

water_data <- model_data %>%
  filter(as.character(NAME10) %in% water_cities) %>%
  mutate(
    NAME10 = fct_drop(NAME10),
    locality_id = fct_drop(locality_id),
    protocol_name = fct_drop(protocol_name)
  )

rm(model_data)
gc()

water_info <- habitat_info %>%
  filter(habitat == "Water")

water_support_path <- fs::path(
  model_dir,
  paste0("Water_min_", water_min_range, "pp_model.rds")
)

water_support_fit <- fit_model(
  water_data,
  water_info$within_term,
  water_info$context_term,
  water_support_path
)

water_sensitivity <- bind_rows(
  interaction_results %>%
    filter(habitat == "Water") %>%
    mutate(
      model = "Primary Beta: all cities",
      n_cities = n_distinct(city_aridity$NAME10),
      n_checklists = sample_summary$n_checklists,
      .before = Estimate
    ),
  extract_interaction(
    water_support_fit,
    "Water",
    water_info$within_term,
    paste0("Beta: Water range >= ", water_min_range, " pp")
  ) %>%
    mutate(
      n_cities = n_distinct(water_data$NAME10),
      n_checklists = nrow(water_data),
      .before = Estimate
    )
)

water_support_qc <- extract_qc(
  water_support_fit,
  "Water",
  water_info$within_term
) %>%
  mutate(
    n_cities = n_distinct(water_data$NAME10),
    n_localities = n_distinct(water_data$locality_id),
    .after = n_checklists
  )

rm(water_support_fit, water_data)
gc()

# Figures

habitat_levels <- c("Water", "Shrub and scrub", "Built")

plot_curves <- resource_curves %>%
  mutate(habitat = factor(habitat, levels = habitat_levels))

plot_city_effects <- city_effects %>%
  mutate(habitat = factor(habitat, levels = habitat_levels))

resource_figure <- ggplot(
  plot_curves,
  aes(x = avg_aridity_index, y = relative_change_percent_at_typical_baseline)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_ribbon(
    aes(
      ymin = relative_CI_lower_percent_at_typical_baseline,
      ymax = relative_CI_upper_percent_at_typical_baseline
    ),
    alpha = 0.2
  ) +
  geom_point(
    data = plot_city_effects,
    aes(x = avg_aridity_index, y = relative_change_percent),
    inherit.aes = FALSE,
    alpha = 0.45,
    size = 1.4
  ) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ habitat, nrow = 1) +
  scale_x_reverse() +
  labs(
    x = "De Martonne aridity index",
    y = "Expected change in response (%)",
    caption = "Less arid  <-                                   ->  More arid"
  ) +
  theme_classic(base_size = 12) +
  theme(plot.caption = element_text(hjust = 0.5))

interaction_figure <- interaction_results %>%
  mutate(habitat = factor(habitat, levels = rev(habitat_levels))) %>%
  ggplot(aes(x = Estimate_1pp, y = habitat)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(
    aes(xmin = CI_lower_1pp, xmax = CI_upper_1pp),
    orientation = "y",
    width = 0.12
  ) +
  geom_point(size = 2.2) +
  labs(
    x = "Habitat x aridity interaction estimate per +1 percentage point (95% CI)",
    y = NULL
  ) +
  theme_classic(base_size = 12)

# Save outputs

write_csv(
  interaction_results,
  fs::path(table_dir, "Functional_Richness_Habitat_x_Aridity.csv")
)
write_csv(
  model_qc,
  fs::path(table_dir, "Functional_Richness_Model_QC.csv")
)
write_csv(
  resource_curves,
  fs::path(table_dir, "Functional_Richness_Effect_Curves.csv")
)
write_csv(
  city_effects,
  fs::path(table_dir, "Functional_Richness_City_Conditional_Effects.csv")
)
write_csv(
  water_sensitivity,
  fs::path(table_dir, "Functional_Richness_Water_Support_Sensitivity.csv")
)
write_csv(
  water_support_qc,
  fs::path(table_dir, "Functional_Richness_Water_Support_Model_QC.csv")
)
write_csv(
  water_support,
  fs::path(table_dir, "Functional_Richness_Water_Support_By_City.csv")
)
write_csv(
  sample_summary,
  fs::path(table_dir, "Functional_Richness_Sample_Summary.csv")
)

model_bundle <- list(
  response = "Functional richness",
  model_family = "Beta",
  link = "logit",
  primary_model_paths = model_registry$model_path,
  model_registry = model_registry,
  typical_city_FRic_baseline_for_display = typical_city_FRic_baseline,
  mean_city_aridity = mean_city_aridity,
  primary_resource_curves = resource_curves
)

write_rds(
  model_bundle,
  fs::path(data_dir, "Functional_Richness_Resource_Importance_Models.rds")
)

write_rds(
  list(
    primary_resource_curves = resource_curves,
    model_diagnostics = model_qc,
    water_support_sensitivity = water_sensitivity,
    water_support_model_qc = water_support_qc,
    water_support_model_path = water_support_path
  ),
  fs::path(data_dir, "Functional_Richness_Diagnostic_Bundle.rds")
)

ggsave(
  fs::path(figure_dir, "Figure_01_Functional_Richness_Resource_Importance.png"),
  resource_figure,
  width = 12,
  height = 4.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  fs::path(figure_dir, "Figure_01_Functional_Richness_Resource_Importance.pdf"),
  resource_figure,
  width = 12,
  height = 4.5,
  units = "in",
  bg = "white"
)

ggsave(
  fs::path(figure_dir, "Figure_02_Functional_Richness_Habitat_x_Aridity.png"),
  interaction_figure,
  width = 7,
  height = 4,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  fs::path(figure_dir, "Figure_02_Functional_Richness_Habitat_x_Aridity.pdf"),
  interaction_figure,
  width = 7,
  height = 4,
  units = "in",
  bg = "white"
)
