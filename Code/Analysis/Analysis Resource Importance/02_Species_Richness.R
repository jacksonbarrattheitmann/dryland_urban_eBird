library(dplyr)
library(ggplot2)
library(glmmTMB)
library(purrr)
library(readr)
library(tibble)

# Settings

min_checklists <- 10
min_localities <- 5
water_min_range <- 1
refit_models <- FALSE
run_dharma <- TRUE
dharma_models <- "all"
dharma_n <- 250
dharma_sample_n <- 10000
dharma_outlier_boot <- 100
dharma_seed <- 123
overwrite_dharma <- FALSE

if (run_dharma && !requireNamespace("DHARMa", quietly = TRUE)) {
  stop("Install the DHARMa package to run model diagnostics.")
}

input_file <- fs::path(
  "Outputs", "Analysis_Resource_Importance", "01_Data_Setup", "Data",
  "Resource_Importance_Analysis_Base.rds"
)

output_root <- fs::path(
  "Outputs", "Analysis_Resource_Importance", "02_Species_Richness"
)

data_dir <- fs::path(output_root, "Data")
table_dir <- fs::path(output_root, "Table")
model_dir <- fs::path(output_root, "Models")
figure_dir <- fs::path(output_root, "Figure")

walk(c(data_dir, table_dir, model_dir, figure_dir), fs::dir_create)

if (!fs::file_exists(input_file)) {
  stop("Resource_Importance_Analysis_Base.rds was not found.")
}

if (!refit_models) {
  expected_models <- c(
    fs::path(model_dir, c(
      "Water_model.rds",
      "Shrub_and_scrub_model.rds",
      "Built_model.rds"
    )),
    fs::path(model_dir, paste0("Water_min_", water_min_range, "pp_model.rds"))
  )
  missing_models <- expected_models[!fs::file_exists(expected_models)]
  
  if (length(missing_models) > 0) {
    stop("Saved model(s) missing. No models were refitted: ",
         paste(missing_models, collapse = ", "))
  }
}

# Species richness sample

model_vars <- c(
  "sr", "avg_aridity_index", "water_pp", "shrub_pp", "built_pp",
  "duration_10", "effort_distance_km", "protocol_name", "NAME10",
  "locality_id"
)

required_vars <- c("checklist_id", model_vars)
analysis_base <- read_rds(input_file)

missing_vars <- setdiff(required_vars, names(analysis_base))
if (length(missing_vars) > 0) {
  stop("Missing variables: ", paste(missing_vars, collapse = ", "))
}

if (anyDuplicated(analysis_base$checklist_id) > 0) {
  stop("Duplicate checklist IDs remain after data setup.")
}

analysis_base <- analysis_base %>%
  select(all_of(required_vars))

sr_pool <- analysis_base %>%
  filter(if_all(all_of(model_vars), ~ !is.na(.x)))

rm(analysis_base)
gc()

numeric_vars <- c(
  "sr", "avg_aridity_index", "water_pp", "shrub_pp", "built_pp",
  "duration_10", "effort_distance_km"
)

finite_check <- sr_pool %>%
  summarise(across(all_of(numeric_vars), ~ all(is.finite(.x)))) %>%
  unlist(use.names = FALSE)

if (!all(finite_check) || any(sr_pool$sr <= 0)) {
  stop("Invalid values remain in the species-richness analysis pool.")
}

eligible_localities <- sr_pool %>%
  count(NAME10, locality_id, name = "n_checklists") %>%
  filter(n_checklists >= min_checklists)

eligible_cities <- eligible_localities %>%
  count(NAME10, name = "n_localities") %>%
  filter(n_localities >= min_localities)

sr_data <- sr_pool %>%
  semi_join(eligible_localities, by = c("NAME10", "locality_id")) %>%
  semi_join(eligible_cities, by = "NAME10")

rm(sr_pool, eligible_localities, eligible_cities, finite_check)
gc()

# Final aridity and habitat centering

city_aridity <- sr_data %>%
  distinct(NAME10, avg_aridity_index)

if (nrow(city_aridity) != n_distinct(sr_data$NAME10)) {
  stop("At least one retained city has more than one aridity value.")
}

mean_city_aridity <- mean(city_aridity$avg_aridity_index)

locality_habitat <- sr_data %>%
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

sr_model_data <- sr_data %>%
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

centering_qc <- sr_model_data %>%
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

sample_summary <- tibble(
  n_checklists = nrow(sr_model_data),
  n_localities = n_distinct(sr_model_data$locality_id),
  n_cities = n_distinct(sr_model_data$NAME10),
  mean_city_aridity = mean_city_aridity,
  min_aridity = min(city_aridity$avg_aridity_index),
  max_aridity = max(city_aridity$avg_aridity_index)
)

write_rds(
  sr_model_data,
  fs::path(data_dir, "Species_Richness_Resource_Importance_Data.rds")
)

rm(sr_data, locality_habitat, centering_qc)
gc()

# Model data
sr_model_data <- sr_model_data %>%
  select(
    sr,
    aridity_centered,
    water_within_model_10,
    shrub_within_model_10,
    built_within_model_10,
    city_context_water_10,
    city_context_shrub_10,
    city_context_built_10,
    duration_10,
    effort_distance_km,
    protocol_name,
    NAME10,
    locality_id
  )

gc()

# Model setup

habitat_info <- tribble(
  ~habitat,          ~stub,              ~within_term,                ~context_term,
  "Water",           "Water",           "water_within_model_10",    "city_context_water_10",
  "Shrub and scrub", "Shrub_and_scrub", "shrub_within_model_10",    "city_context_shrub_10",
  "Built",           "Built",           "built_within_model_10",    "city_context_built_10"
)

make_formula <- function(within_term, context_term) {
  as.formula(
    paste0(
      "sr ~ ", within_term, " * aridity_centered + ", context_term,
      " + duration_10 + effort_distance_km + protocol_name + ",
      "(1 + ", within_term, " || NAME10) + (1 | locality_id)"
    )
  )
}

model_ok <- function(fit) {
  inherits(fit, "glmmTMB") &&
    isTRUE(fit$fit$convergence == 0) &&
    isTRUE(fit$sdr$pdHess)
}

model_matches <- function(fit, expected_formula, n) {
  if (!model_ok(fit)) return(FALSE)
  
  tryCatch({
    saved_formula <- paste(deparse(formula(fit)), collapse = "")
    expected_formula <- paste(deparse(expected_formula), collapse = "")
    family_info <- family(fit)
    
    identical(gsub("\\s+", "", saved_formula),
              gsub("\\s+", "", expected_formula)) &&
      nobs(fit) == n &&
      identical(family_info$family, "nbinom2") &&
      identical(family_info$link, "log")
  }, error = function(e) FALSE)
}

fit_model <- function(data, within_term, context_term, model_path) {
  model_formula <- make_formula(within_term, context_term)
  
  if (!refit_models) {
    fit <- tryCatch(read_rds(model_path), error = function(e) e)
    
    if (!model_matches(fit, model_formula, nrow(data))) {
      stop("Saved model could not be validated. No refitting performed: ",
           fs::path_file(model_path))
    }
    
    message("Loaded model: ", fs::path_file(model_path))
    return(fit)
  }
  
  fit <- tryCatch(
    glmmTMB(
      model_formula,
      data = data,
      family = nbinom2(link = "log"),
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
        family = nbinom2(link = "log"),
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

extract_coefficients <- function(fit, habitat, model) {
  x <- summary(fit)$coefficients$cond
  
  tibble(
    habitat = habitat,
    model = model,
    term = rownames(x),
    Estimate = unname(x[, "Estimate"]),
    SE = unname(x[, "Std. Error"]),
    z = unname(x[, "z value"]),
    p = unname(x[, "Pr(>|z|)"]),
    CI_lower = Estimate - 1.96 * SE,
    CI_upper = Estimate + 1.96 * SE
  )
}

extract_model_summary <- function(fit, habitat, model) {
  tibble(
    habitat = habitat,
    model = model,
    family = family(fit)$family,
    link = family(fit)$link,
    n_checklists = nobs(fit),
    AIC = AIC(fit),
    BIC = BIC(fit),
    logLik = as.numeric(logLik(fit)),
    nb2_theta = as.numeric(sigma(fit)),
    convergence = fit$fit$convergence,
    pdHess = isTRUE(fit$sdr$pdHess)
  )
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

extract_interaction <- function(fit, habitat, within_term, model = "Primary NB2") {
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
    CI_upper = Estimate + 1.96 * SE
  )
}

extract_qc <- function(fit, habitat, within_term) {
  pearson <- residuals(fit, type = "pearson")
  residual_df <- df.residual(fit)
  pearson_ratio <- sum(pearson^2, na.rm = TRUE) / residual_df
  
  random_sd <- attr(VarCorr(fit)$cond$NAME10, "stddev")
  slope_sd <- unname(random_sd[within_term])
  
  vc <- vcov(fit)$cond
  cor_mat <- cov2cor(vc)
  diag(cor_mat) <- NA_real_
  
  out <- tibble(
    habitat = habitat,
    n_checklists = nobs(fit),
    convergence = fit$fit$convergence,
    pdHess = isTRUE(fit$sdr$pdHess),
    pearson_dispersion = pearson_ratio,
    city_random_slope_SD = slope_sd,
    max_fixed_effect_correlation = max(abs(cor_mat), na.rm = TRUE)
  )
  
  rm(pearson)
  gc()
  
  out
}

make_city_effects <- function(fit, habitat, within_term) {
  beta <- fixef(fit)$cond
  int_term <- interaction_term(fit, within_term)
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
      fixed_slope = beta[within_term] + beta[int_term] * aridity_centered,
      city_slope = fixed_slope + random_slope,
      relative_change_percent = 100 * (exp(0.1 * city_slope) - 1)
    )
}

make_curve <- function(fit, habitat, within_term, aridity_values) {
  beta <- fixef(fit)$cond
  vc <- vcov(fit)$cond
  int_term <- interaction_term(fit, within_term)
  centered <- aridity_values - mean_city_aridity
  
  b_habitat <- beta[within_term]
  b_interaction <- beta[int_term]
  
  effect <- 0.1 * (b_habitat + b_interaction * centered)
  effect_var <- 0.1^2 * (
    vc[within_term, within_term] +
      centered^2 * vc[int_term, int_term] +
      2 * centered * vc[within_term, int_term]
  )
  effect_se <- sqrt(pmax(effect_var, 0))
  
  tibble(
    habitat = habitat,
    avg_aridity_index = aridity_values,
    relative_change_percent = 100 * (exp(effect) - 1),
    relative_CI_lower_percent = 100 * (exp(effect - 1.96 * effect_se) - 1),
    relative_CI_upper_percent = 100 * (exp(effect + 1.96 * effect_se) - 1)
  )
}

# Primary models

aridity_values <- seq(
  min(city_aridity$avg_aridity_index),
  max(city_aridity$avg_aridity_index),
  length.out = 200
)

model_registry <- habitat_info %>%
  mutate(
    model_set = "Primary no response offset",
    model_path = fs::path(model_dir, paste0(stub, "_model.rds"))
  )

run_primary_model <- function(habitat, within_term, context_term, model_path) {
  fit <- fit_model(sr_model_data, within_term, context_term, model_path)
  
  out <- list(
    interaction = extract_interaction(fit, habitat, within_term),
    qc = extract_qc(fit, habitat, within_term),
    curve = make_curve(fit, habitat, within_term, aridity_values),
    city_effects = make_city_effects(fit, habitat, within_term),
    coefficients = extract_coefficients(fit, habitat, "Primary NB2"),
    summary = extract_model_summary(fit, habitat, "Primary NB2")
  )
  
  rm(fit)
  gc()
  
  out
}

primary_results <- model_registry %>%
  select(habitat, within_term, context_term, model_path) %>%
  pmap(run_primary_model)

interaction_results <- map_dfr(primary_results, "interaction")
model_qc <- map_dfr(primary_results, "qc")
resource_curves <- map_dfr(primary_results, "curve")
city_effects <- map_dfr(primary_results, "city_effects")
model_coefficients <- map_dfr(primary_results, "coefficients")
model_summaries <- map_dfr(primary_results, "summary")

rm(primary_results)
gc()

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

water_data <- sr_model_data %>%
  filter(as.character(NAME10) %in% water_cities) %>%
  select(
    sr,
    water_within_model_10,
    aridity_centered,
    city_context_water_10,
    duration_10,
    effort_distance_km,
    protocol_name,
    NAME10,
    locality_id
  ) %>%
  droplevels()

rm(sr_model_data)
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

water_support_qc <- extract_qc(
  water_support_fit,
  "Water support",
  water_info$within_term
)

model_coefficients <- bind_rows(
  model_coefficients,
  extract_coefficients(water_support_fit, "Water support", "Water support NB2")
)

model_summaries <- bind_rows(
  model_summaries,
  extract_model_summary(water_support_fit, "Water support", "Water support NB2")
)

water_sensitivity <- bind_rows(
  interaction_results %>%
    filter(habitat == "Water") %>%
    mutate(
      model = "Primary NB2: all cities",
      n_cities = sample_summary$n_cities,
      n_checklists = sample_summary$n_checklists,
      .before = Estimate
    ),
  extract_interaction(
    water_support_fit,
    "Water",
    water_info$within_term,
    paste0("NB2: Water range >= ", water_min_range, " pp")
  ) %>%
    mutate(
      n_cities = n_distinct(water_data$NAME10),
      n_checklists = nrow(water_data),
      .before = Estimate
    )
)

rm(water_support_fit, water_data)
gc()

# Figures

resource_figure <- ggplot(
  resource_curves,
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
  geom_point(
    data = city_effects,
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
  mutate(
    habitat = factor(habitat, levels = rev(habitat_info$habitat)),
    estimate_1pp = Estimate / 10,
    ci_lower_1pp = CI_lower / 10,
    ci_upper_1pp = CI_upper / 10
  ) %>%
  ggplot(aes(x = estimate_1pp, y = habitat)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(
    aes(xmin = ci_lower_1pp, xmax = ci_upper_1pp),
    orientation = "y",
    width = 0.12
  ) +
  geom_point(size = 2.2) +
  labs(
    x = "Habitat x aridity interaction estimate (95% CI)",
    y = NULL
  ) +
  theme_classic(base_size = 12)

# Save outputs

write_csv(
  interaction_results,
  fs::path(table_dir, "Species_Richness_Habitat_x_Aridity.csv")
)
write_csv(
  model_qc,
  fs::path(table_dir, "Species_Richness_Model_QC.csv")
)
write_csv(
  water_support_qc,
  fs::path(table_dir, "Species_Richness_Water_Support_Model_QC.csv")
)
write_csv(
  model_coefficients,
  fs::path(table_dir, "Species_Richness_All_Model_Coefficients.csv")
)
write_csv(
  model_summaries,
  fs::path(table_dir, "Species_Richness_All_Model_Summaries.csv")
)
write_csv(
  resource_curves,
  fs::path(table_dir, "Species_Richness_Effect_Curves.csv")
)
write_csv(
  city_effects,
  fs::path(table_dir, "Species_Richness_City_Conditional_Effects.csv")
)
write_csv(
  water_sensitivity,
  fs::path(table_dir, "Species_Richness_Water_Support_Sensitivity.csv")
)
write_csv(
  water_support,
  fs::path(table_dir, "Species_Richness_Water_Support_By_City.csv")
)
write_csv(
  sample_summary,
  fs::path(table_dir, "Species_Richness_Sample_Summary.csv")
)

ggsave(
  fs::path(figure_dir, "Figure_01_Species_Richness_Resource_Importance.png"),
  resource_figure,
  width = 12,
  height = 4.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  fs::path(figure_dir, "Figure_01_Species_Richness_Resource_Importance.pdf"),
  resource_figure,
  width = 12,
  height = 4.5,
  units = "in",
  bg = "white"
)

ggsave(
  fs::path(figure_dir, "Figure_02_Species_Richness_Habitat_x_Aridity.png"),
  interaction_figure,
  width = 7,
  height = 4,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  fs::path(figure_dir, "Figure_02_Species_Richness_Habitat_x_Aridity.pdf"),
  interaction_figure,
  width = 7,
  height = 4,
  units = "in",
  bg = "white"
)

# Model index

write_rds(
  list(
    primary_model_paths = model_registry$model_path,
    model_registry = model_registry
  ),
  fs::path(data_dir, "Species_Richness_Resource_Importance_Models.rds")
)


# DHARMa diagnostics

if (run_dharma) {
  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    stop("Install the DHARMa package to run model diagnostics.")
  }
  
  diagnostic_registry <- bind_rows(
    model_registry %>%
      select(habitat, stub, model_path),
    tibble(
      habitat = "Water support",
      stub = paste0("Water_min_", water_min_range, "pp"),
      model_path = water_support_path
    )
  )
  
  if (identical(dharma_models, "all")) {
    selected_registry <- diagnostic_registry
  } else {
    invalid_names <- setdiff(dharma_models, diagnostic_registry$stub)
    if (length(invalid_names) > 0) {
      stop("Unknown DHARMa model name(s): ",
           paste(invalid_names, collapse = ", "))
    }
    selected_registry <- diagnostic_registry %>%
      filter(stub %in% dharma_models)
  }
  
  if (length(dharma_sample_n) != 1 || !is.finite(dharma_sample_n) ||
      dharma_sample_n < 100 || dharma_sample_n != as.integer(dharma_sample_n)) {
    stop("dharma_sample_n must be an integer of at least 100.")
  }
  
  if (length(dharma_n) != 1 || !is.finite(dharma_n) ||
      dharma_n < 20 || dharma_n != as.integer(dharma_n)) {
    stop("dharma_n must be an integer of at least 20.")
  }
  
  if (length(dharma_outlier_boot) != 1 || !is.finite(dharma_outlier_boot) ||
      dharma_outlier_boot < 20 ||
      dharma_outlier_boot != as.integer(dharma_outlier_boot) ||
      dharma_outlier_boot > dharma_n) {
    stop("dharma_outlier_boot must be an integer between 20 and dharma_n.")
  }
  
  run_dharma_model <- function(habitat, stub, model_path) {
    test_path <- fs::path(
      table_dir, paste0("Species_Richness_DHARMa_", stub, "_Tests.csv")
    )
    residual_path <- fs::path(
      data_dir, paste0("Species_Richness_DHARMa_", stub, "_Residuals.csv")
    )
    figure_path <- fs::path(
      figure_dir, paste0("Species_Richness_DHARMa_", stub, ".png")
    )
    predictor_figure_path <- fs::path(
      figure_dir, paste0("Species_Richness_DHARMa_", stub, "_Predictors.png")
    )
    object_path <- fs::path(
      data_dir, paste0("Species_Richness_DHARMa_", stub, ".rds")
    )
    
    diagnostic_files <- c(
      test_path, residual_path, figure_path, predictor_figure_path, object_path
    )
    dharma_version <- as.character(utils::packageVersion("DHARMa"))
    
    if (!overwrite_dharma && all(fs::file_exists(diagnostic_files))) {
      old_tests <- tryCatch(
        read_csv(test_path, show_col_types = FALSE),
        error = function(e) NULL
      )
      
      if (!is.null(old_tests) && nrow(old_tests) == 4 &&
          all(c(
            "n_simulations", "n_sampled", "seed", "simulation",
            "DHARMa_version", "n_checklists", "outlier_method", "n_boot"
          ) %in% names(old_tests)) &&
          all(old_tests$n_simulations == dharma_n) &&
          all(old_tests$outlier_method == "bootstrap") &&
          all(old_tests$n_boot == dharma_outlier_boot) &&
          all(old_tests$n_sampled == pmin(dharma_sample_n, old_tests$n_checklists)) &&
          all(old_tests$seed == dharma_seed) &&
          all(old_tests$simulation == "Conditional NB2") &&
          all(old_tests$DHARMa_version == dharma_version) &&
          file.info(test_path)$mtime >= file.info(model_path)$mtime) {
        message("Loaded DHARMa results: ", stub)
        return(old_tests)
      }
    }
    
    fit <- read_rds(model_path)
    frame <- fit$frame
    if (is.null(frame) || nrow(frame) != nobs(fit)) {
      stop("The saved model frame is unavailable: ", stub)
    }
    
    if (!identical(family(fit)$family, "nbinom2") ||
        !identical(family(fit)$link, "log") ||
        length(fixef(fit)$zi) > 0) {
      stop("Conditional NB2 diagnostics require an NB2 log model without zero inflation.")
    }
    
    theta <- as.numeric(sigma(fit))
    if (length(theta) != 1 || !is.finite(theta) || theta <= 0) {
      stop("The model has no valid constant NB2 dispersion parameter: ", stub)
    }
    
    habitat_term <- if (stub == "Shrub_and_scrub") {
      "shrub_within_model_10"
    } else if (stub == "Built") {
      "built_within_model_10"
    } else {
      "water_within_model_10"
    }
    
    if (!all(c("sr", "aridity_centered", habitat_term, "NAME10") %in%
             names(frame))) {
      stop("Required diagnostic columns are missing: ", stub)
    }
    
    set.seed(dharma_seed)
    sampled_rows <- sort(sample.int(
      nobs(fit), size = min(dharma_sample_n, nobs(fit))
    ))
    sampled_frame <- frame[sampled_rows, , drop = FALSE]
    
    predicted_response <- as.numeric(predict(
      fit,
      newdata = sampled_frame,
      type = "response",
      re.form = NULL,
      allow.new.levels = FALSE
    ))
    observed_response <- as.numeric(sampled_frame$sr)
    
    if (length(predicted_response) != length(observed_response) ||
        any(!is.finite(predicted_response)) ||
        any(predicted_response <= 0) ||
        any(!is.finite(observed_response))) {
      stop("Invalid sampled predictions or responses: ", stub)
    }
    
    message(
      "DHARMa: ", stub, " - ", length(sampled_rows), " of ", nobs(fit),
      " checklists, ", dharma_n, " conditional NB2 simulations"
    )
    
    simulated_response <- replicate(
      dharma_n,
      stats::rnbinom(
        n = length(predicted_response),
        mu = predicted_response,
        size = theta
      )
    )
    
    sim <- DHARMa::createDHARMa(
      simulatedResponse = simulated_response,
      observedResponse = observed_response,
      fittedPredictedResponse = predicted_response,
      integerResponse = TRUE,
      seed = dharma_seed
    )
    
    tests <- list(
      Uniformity = DHARMa::testUniformity(sim, plot = FALSE),
      Dispersion = DHARMa::testDispersion(sim, plot = FALSE),
      Zero_inflation = DHARMa::testZeroInflation(sim, plot = FALSE),
      Outliers = DHARMa::testOutliers(
        sim, type = "bootstrap", nBoot = dharma_outlier_boot, plot = FALSE
      )
    )
    
    test_results <- imap_dfr(tests, function(result, test_name) {
      if (!is.list(result) || is.null(result$statistic) ||
          is.null(result$p.value)) {
        stop("Unexpected DHARMa test result for: ", test_name)
      }
      
      tibble(
        habitat = habitat,
        stub = stub,
        test = test_name,
        statistic = as.numeric(result$statistic)[1],
        p = as.numeric(result$p.value)[1],
        n_checklists = nobs(fit),
        n_sampled = length(sampled_rows),
        n_simulations = dharma_n,
        simulation = "Conditional NB2",
        outlier_method = "bootstrap",
        n_boot = dharma_outlier_boot,
        seed = dharma_seed,
        DHARMa_version = dharma_version
      )
    })
    
    diagnostic_data <- tibble(
      model_row = sampled_rows,
      NAME10 = as.character(sampled_frame$NAME10),
      aridity_centered = sampled_frame$aridity_centered,
      habitat_within_model_10 = sampled_frame[[habitat_term]],
      scaled_residual = sim$scaledResiduals,
      observed_response = observed_response,
      predicted_response = predicted_response
    )
    
    write_csv(test_results, test_path)
    write_csv(diagnostic_data, residual_path)
    write_rds(sim, object_path)
    
    grDevices::png(figure_path, width = 1600, height = 800, res = 150)
    tryCatch({
      graphics::par(mfrow = c(1, 2))
      DHARMa::plotQQunif(sim, testOutliers = FALSE)
      DHARMa::plotResiduals(sim, quantreg = FALSE)
    }, finally = grDevices::dev.off())
    
    grDevices::png(
      predictor_figure_path, width = 1600, height = 800, res = 150
    )
    tryCatch({
      graphics::par(mfrow = c(1, 2))
      DHARMa::plotResiduals(
        sim, form = diagnostic_data$aridity_centered,
        quantreg = FALSE, rank = FALSE, xlab = "Centered aridity"
      )
      DHARMa::plotResiduals(
        sim, form = diagnostic_data$habitat_within_model_10,
        quantreg = FALSE, rank = FALSE, xlab = "Within-city habitat (per 10 pp)"
      )
    }, finally = grDevices::dev.off())
    
    rm(fit, frame, sampled_frame, sim, simulated_response, diagnostic_data)
    gc()
    
    test_results
  }
  
  dharma_results <- pmap_dfr(selected_registry, run_dharma_model)
  write_csv(
    dharma_results,
    fs::path(table_dir, "Species_Richness_DHARMa_Tests.csv")
  )
}
