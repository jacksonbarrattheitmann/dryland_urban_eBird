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
  "Outputs", "Analysis_Resource_Importance", "06_CWM_Body_Mass"
)

data_dir <- fs::path(output_root, "Data")
model_dir <- fs::path(output_root, "Models")
table_dir <- fs::path(output_root, "Table")
figure_dir <- fs::path(output_root, "Figure")

walk(
  c(data_dir, model_dir, table_dir, figure_dir),
  fs::dir_create
)

# CWM body-mass sample

analysis_base <- read_rds(input_file)

model_vars <- c(
  "checklist_id", "CWM.Mass", "sr", "avg_aridity_index",
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

mass_pool <- analysis_base %>%
  filter(if_all(all_of(model_vars), ~ !is.na(.x)))

rm(analysis_base)
gc()

numeric_vars <- c(
  "CWM.Mass", "sr", "avg_aridity_index",
  "water_pp", "shrub_pp", "built_pp",
  "duration_10", "effort_distance_km"
)

finite_check <- mass_pool %>%
  summarise(across(all_of(numeric_vars), ~ all(is.finite(.x)))) %>%
  unlist(use.names = FALSE)

if (!all(finite_check)) {
  stop("Non-finite values remain in the CWM body-mass analysis pool.")
}

if (any(mass_pool$CWM.Mass <= 0)) {
  stop("CWM body mass must be strictly positive for Gamma regression.")
}

eligible_localities <- mass_pool %>%
  count(NAME10, locality_id, name = "n_checklists") %>%
  filter(n_checklists >= min_checklists)

eligible_cities <- eligible_localities %>%
  count(NAME10, name = "n_localities") %>%
  filter(n_localities >= min_localities)

mass_data <- mass_pool %>%
  semi_join(eligible_localities, by = c("NAME10", "locality_id")) %>%
  semi_join(eligible_cities, by = "NAME10")

rm(mass_pool, eligible_localities, eligible_cities)
gc()

# City baseline and final habitat centering

city_aridity <- mass_data %>%
  distinct(NAME10, avg_aridity_index)

if (nrow(city_aridity) != n_distinct(mass_data$NAME10)) {
  stop("At least one retained city has more than one aridity value.")
}

mean_city_aridity <- mean(city_aridity$avg_aridity_index)

city_baseline <- mass_data %>%
  group_by(NAME10) %>%
  summarise(
    city_mean_CWM_Mass_checklist = mean(CWM.Mass),
    .groups = "drop"
  )

if (any(!is.finite(city_baseline$city_mean_CWM_Mass_checklist) |
        city_baseline$city_mean_CWM_Mass_checklist <= 0)) {
  stop("At least one retained city has an invalid CWM body-mass baseline.")
}

locality_habitat <- mass_data %>%
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

mass_analysis_data <- mass_data %>%
  select(-any_of(c(
    "city_mean_CWM_Mass_checklist",
    "city_mean_water_pp", "city_mean_shrub_pp", "city_mean_built_pp"
  ))) %>%
  left_join(city_baseline, by = "NAME10") %>%
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
    city_CWM_Mass_log_offset = log(city_mean_CWM_Mass_checklist),
    NAME10 = factor(NAME10),
    locality_id = factor(locality_id),
    protocol_name = factor(protocol_name, levels = c("Stationary", "Traveling"))
  )

rm(mass_data, locality_habitat)
gc()

centering_qc <- mass_analysis_data %>%
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
  mass_analysis_data,
  fs::path(data_dir, "CWM_Body_Mass_Resource_Importance_Data.rds")
)

sample_summary <- tibble(
  n_checklists = nrow(mass_analysis_data),
  n_localities = n_distinct(mass_analysis_data$locality_id),
  n_cities = n_distinct(mass_analysis_data$NAME10),
  mean_CWM_Mass = mean(mass_analysis_data$CWM.Mass),
  median_CWM_Mass = median(mass_analysis_data$CWM.Mass),
  mean_city_aridity = mean_city_aridity,
  min_aridity = min(city_aridity$avg_aridity_index),
  max_aridity = max(city_aridity$avg_aridity_index)
)

# Model setup

habitat_info <- tribble(
  ~habitat,          ~stub,              ~within_term,              ~context_term,
  "Water",           "Water",           "water_within_model_10", "city_context_water_10",
  "Shrub and scrub", "Shrub_and_scrub", "shrub_within_model_10", "city_context_shrub_10",
  "Built",           "Built",           "built_within_model_10", "city_context_built_10"
)

make_formula <- function(within_term, context_term) {
  as.formula(
    paste0(
      "CWM.Mass ~ ", within_term, " * aridity_centered + ", context_term,
      " + duration_10 + effort_distance_km + protocol_name + ",
      "offset(city_CWM_Mass_log_offset) + ",
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
  family_info <- family(fit)
  
  identical(gsub("\\s+", "", saved_formula), gsub("\\s+", "", expected_formula)) &&
    identical(as.character(family_info$family), "Gamma") &&
    identical(as.character(family_info$link), "log")
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
      family = Gamma(link = "log"),
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
        family = Gamma(link = "log"),
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

extract_interaction <- function(fit, habitat, within_term, model = "Primary Gamma") {
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
  pearson_ratio <- sum(pearson^2, na.rm = TRUE) /
    max(1, nobs(fit) - length(beta))
  
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
    gamma_dispersion = as.numeric(sigma(fit)),
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
  
  tibble(
    habitat = habitat,
    avg_aridity_index = aridity_values,
    aridity_centered = centered,
    habitat_increment_percent_cover = habitat_step * 10,
    delta_eta = delta,
    delta_eta_SE = delta_se,
    response_ratio = exp(delta),
    response_ratio_CI_lower = exp(delta - 1.96 * delta_se),
    response_ratio_CI_upper = exp(delta + 1.96 * delta_se),
    relative_change_percent = 100 * (exp(delta) - 1),
    relative_CI_lower_percent = 100 * (exp(delta - 1.96 * delta_se) - 1),
    relative_CI_upper_percent = 100 * (exp(delta + 1.96 * delta_se) - 1)
  )
}

# Primary models

model_data <- mass_analysis_data %>%
  select(
    CWM.Mass,
    water_within_model_10, shrub_within_model_10, built_within_model_10,
    aridity_centered,
    city_context_water_10, city_context_shrub_10, city_context_built_10,
    duration_10, effort_distance_km, protocol_name,
    city_CWM_Mass_log_offset, NAME10, locality_id
  )

rm(mass_analysis_data)
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
    curve = make_curve(fit, habitat, within_term, aridity_values)
  )
  
  rm(fit)
  gc()
  
  out
}

primary_outputs <- pmap(habitat_info, process_habitat)

interaction_results <- map_dfr(primary_outputs, "interaction")
model_qc <- map_dfr(primary_outputs, "qc")
resource_curves <- map_dfr(primary_outputs, "curve")

model_registry <- habitat_info %>%
  mutate(
    response = "CWM.Mass",
    family = "Gamma",
    link = "log",
    model_path = map_chr(primary_outputs, "model_path"),
    .before = within_term
  )

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
      model = "Primary Gamma: all cities",
      n_cities = n_distinct(city_aridity$NAME10),
      n_localities = sample_summary$n_localities,
      n_checklists = sample_summary$n_checklists,
      .before = Estimate
    ),
  extract_interaction(
    water_support_fit,
    "Water",
    water_info$within_term,
    paste0("Gamma: Water range >= ", water_min_range, " pp")
  ) %>%
    mutate(
      n_cities = n_distinct(water_data$NAME10),
      n_localities = n_distinct(water_data$locality_id),
      n_checklists = nrow(water_data),
      .before = Estimate
    )
) %>%
  mutate(
    same_interaction_sign = sign(Estimate) == sign(first(Estimate)),
    significant_at_0_05 = p < 0.05
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

resource_figure <- ggplot(
  plot_curves,
  aes(x = avg_aridity_index, y = relative_change_percent)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_ribbon(
    aes(
      ymin = relative_CI_lower_percent,
      ymax = relative_CI_upper_percent
    ),
    alpha = 0.2
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
  fs::path(table_dir, "CWM_Body_Mass_Habitat_x_Aridity.csv")
)
write_csv(
  model_qc,
  fs::path(table_dir, "CWM_Body_Mass_Model_QC.csv")
)
write_csv(
  resource_curves,
  fs::path(table_dir, "CWM_Body_Mass_Effect_Curves.csv")
)
write_csv(
  water_sensitivity,
  fs::path(table_dir, "CWM_Body_Mass_Water_Support_Sensitivity.csv")
)
write_csv(
  water_support_qc,
  fs::path(table_dir, "CWM_Body_Mass_Water_Support_Model_QC.csv")
)
write_csv(
  water_support,
  fs::path(table_dir, "CWM_Body_Mass_Water_Support_By_City.csv")
)
write_csv(
  sample_summary,
  fs::path(table_dir, "CWM_Body_Mass_Sample_Summary.csv")
)

analysis_bundle <- list(
  response = "Community-weighted mean body mass",
  response_variable = "CWM.Mass",
  model_family = "Gamma",
  link = "log",
  city_relative_estimand = paste0(
    "Expected checklist CWM body mass multiplicatively relative to its ",
    "city's observed checklist-weighted CWM body-mass baseline"
  ),
  primary_city_baseline =
    "Checklist-weighted mean CWM.Mass within each final eligible city",
  primary_offset_variable = "city_CWM_Mass_log_offset",
  primary_offset_definition = "log(city_mean_CWM_Mass_checklist)",
  habitat_increment_percent_cover = habitat_step * 10,
  mean_city_aridity = mean_city_aridity,
  primary_model_paths = model_registry$model_path,
  model_registry = model_registry,
  primary_interactions = interaction_results,
  primary_model_qc = model_qc,
  primary_resource_curves = resource_curves,
  water_support_sensitivity = water_sensitivity,
  water_support_model_qc = water_support_qc,
  water_support_model_path = water_support_path
)

write_rds(
  analysis_bundle,
  fs::path(data_dir, "CWM_Body_Mass_Resource_Importance_Analysis.rds")
)

ggsave(
  fs::path(figure_dir, "Figure_01_CWM_Body_Mass_Resource_Importance.png"),
  resource_figure,
  width = 12,
  height = 4.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  fs::path(figure_dir, "Figure_01_CWM_Body_Mass_Resource_Importance.pdf"),
  resource_figure,
  width = 12,
  height = 4.5,
  units = "in",
  bg = "white"
)

ggsave(
  fs::path(figure_dir, "Figure_02_CWM_Body_Mass_Habitat_x_Aridity.png"),
  interaction_figure,
  width = 7,
  height = 4,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  fs::path(figure_dir, "Figure_02_CWM_Body_Mass_Habitat_x_Aridity.pdf"),
  interaction_figure,
  width = 7,
  height = 4,
  units = "in",
  bg = "white"
)
