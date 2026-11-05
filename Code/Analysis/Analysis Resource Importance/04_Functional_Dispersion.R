library(dplyr)
library(ggplot2)
library(glmmTMB)
library(purrr)
library(readr)
library(tibble)
library(tidyr)

# Settings

min_checklists <- 10
min_localities <- 5
water_min_range <- 1
habitat_step <- 0.1
n_aridity_points <- 200
run_dharma <- TRUE
dharma_models <- "all"
dharma_n <- 250
dharma_sample_n <- 10000
dharma_seed <- 123

if (run_dharma && !requireNamespace("DHARMa", quietly = TRUE)) {
  stop("Install the DHARMa package to run model diagnostics.")
}

input_file <- fs::path(
  "Outputs", "Analysis_Resource_Importance", "01_Data_Setup", "Data",
  "Resource_Importance_Analysis_Base.rds"
)

output_root <- fs::path(
  "Outputs", "Analysis_Resource_Importance", "04_Functional_Dispersion"
)

data_dir <- fs::path(output_root, "Data")
model_dir <- fs::path(output_root, "Models")
table_dir <- fs::path(output_root, "Table")
figure_dir <- fs::path(output_root, "Figure")

walk(
  c(data_dir, model_dir, table_dir, figure_dir),
  fs::dir_create
)

if (!fs::file_exists(input_file)) {
  stop("Resource_Importance_Analysis_Base.rds was not found.")
}

legacy_model_files <- fs::dir_ls(
  model_dir,
  recurse = TRUE,
  type = "file",
  regexp = "\\.rds$",
  fail = FALSE
)

if (length(legacy_model_files) > 0) {
  warning(
    "Legacy fitted-model RDS files remain in the FDis Models folder. ",
    "This script will not read, overwrite, or create fitted-model RDS files."
  )
}

# Functional dispersion sample

analysis_base <- read_rds(input_file)

model_vars <- c(
  "checklist_id", "FDis", "sr", "avg_aridity_index",
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

fdis_pool <- analysis_base %>%
  filter(sr >= 2) %>%
  filter(if_all(all_of(model_vars), ~ !is.na(.x)))

rm(analysis_base)
gc()

numeric_vars <- c(
  "FDis", "sr", "avg_aridity_index", "water_pp", "shrub_pp", "built_pp",
  "duration_10", "effort_distance_km"
)

finite_check <- fdis_pool %>%
  summarise(across(all_of(numeric_vars), ~ all(is.finite(.x)))) %>%
  unlist(use.names = FALSE)

if (!all(finite_check)) {
  stop("Non-finite values remain in the functional-dispersion analysis pool.")
}

if (any(fdis_pool$FDis < 0)) {
  stop("Negative FDis values remain in the functional-dispersion analysis pool.")
}

eligible_localities <- fdis_pool %>%
  count(NAME10, locality_id, name = "n_checklists") %>%
  filter(n_checklists >= min_checklists)

eligible_cities <- eligible_localities %>%
  count(NAME10, name = "n_localities") %>%
  filter(n_localities >= min_localities)

fdis_data <- fdis_pool %>%
  semi_join(eligible_localities, by = c("NAME10", "locality_id")) %>%
  semi_join(eligible_cities, by = "NAME10")

rm(fdis_pool, eligible_localities, eligible_cities)
gc()

# City aridity and final habitat centering

city_aridity <- fdis_data %>%
  distinct(NAME10, avg_aridity_index)

if (nrow(city_aridity) != n_distinct(fdis_data$NAME10)) {
  stop("At least one retained city has more than one aridity value.")
}

mean_city_aridity <- mean(city_aridity$avg_aridity_index)

locality_habitat <- fdis_data %>%
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

fdis_analysis_data <- fdis_data %>%
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

rm(fdis_data, locality_habitat)
gc()

centering_qc <- fdis_analysis_data %>%
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
  fdis_analysis_data,
  fs::path(data_dir, "Functional_Dispersion_Resource_Importance_Data.rds"),
  compress = "gz"
)

sample_summary <- tibble(
  n_checklists = nrow(fdis_analysis_data),
  n_localities = n_distinct(fdis_analysis_data$locality_id),
  n_cities = n_distinct(fdis_analysis_data$NAME10),
  mean_FDis = mean(fdis_analysis_data$FDis),
  median_FDis = median(fdis_analysis_data$FDis),
  min_FDis = min(fdis_analysis_data$FDis),
  max_FDis = max(fdis_analysis_data$FDis),
  proportion_FDis_zero = mean(fdis_analysis_data$FDis == 0),
  proportion_FDis_one = mean(fdis_analysis_data$FDis == 1),
  mean_city_aridity = mean_city_aridity,
  min_aridity = min(city_aridity$avg_aridity_index),
  max_aridity = max(city_aridity$avg_aridity_index)
)

gamma_eligible <- sample_summary$min_FDis > 0

# Model setup

habitat_info <- tribble(
  ~habitat,          ~stub,              ~within_term,              ~context_term,
  "Water",           "Water",           "water_within_model_10", "city_context_water_10",
  "Shrub and scrub", "Shrub_and_scrub", "shrub_within_model_10", "city_context_shrub_10",
  "Built",           "Built",           "built_within_model_10", "city_context_built_10"
)

family_info <- tribble(
  ~family_key, ~family_label, ~family_name, ~link,      ~eligible,
  "Gaussian", "Gaussian",    "gaussian",   "identity", TRUE,
  "Gamma",    "Gamma",       "Gamma",      "log",      gamma_eligible
)

make_family <- function(family_key) {
  if (family_key == "Gaussian") {
    return(gaussian(link = "identity"))
  }
  
  if (family_key == "Gamma") {
    return(Gamma(link = "log"))
  }
  
  stop("Unknown family key: ", family_key)
}

make_formula <- function(within_term, context_term) {
  as.formula(
    paste0(
      "FDis ~ ", within_term, " * aridity_centered + ", context_term,
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

fit_model <- function(data, within_term, context_term, family_key) {
  model_formula <- make_formula(within_term, context_term)
  
  fit <- tryCatch(
    glmmTMB(
      model_formula,
      data = data,
      family = make_family(family_key),
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
        family = make_family(family_key),
        control = glmmTMBControl(
          optimizer = optim,
          optArgs = list(method = "BFGS")
        )
      ),
      error = function(e) e
    )
  }
  
  if (!model_ok(fit)) {
    return(NULL)
  }
  
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

extract_candidate_summary <- function(fit, family_key, family_label, habitat) {
  if (is.null(fit)) {
    return(tibble(
      family_key = family_key,
      family_label = family_label,
      habitat = habitat,
      converged = FALSE,
      pdHess = FALSE,
      n_checklists = NA_integer_,
      AIC = NA_real_,
      BIC = NA_real_,
      logLik = NA_real_
    ))
  }
  
  tibble(
    family_key = family_key,
    family_label = family_label,
    habitat = habitat,
    converged = fit$fit$convergence == 0,
    pdHess = isTRUE(fit$sdr$pdHess),
    n_checklists = nobs(fit),
    AIC = AIC(fit),
    BIC = BIC(fit),
    logLik = as.numeric(logLik(fit))
  )
}

extract_coefficients <- function(fit, habitat, model) {
  x <- summary(fit)$coefficients$cond
  
  tibble(
    habitat = habitat,
    model = model,
    family = family(fit)$family,
    link = family(fit)$link,
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
  fam <- family(fit)
  sigma_value <- as.numeric(sigma(fit))
  is_gamma <- identical(fam$family, "Gamma")
  
  tibble(
    habitat = habitat,
    model = model,
    family = fam$family,
    link = fam$link,
    n_checklists = nobs(fit),
    AIC = AIC(fit),
    BIC = BIC(fit),
    logLik = as.numeric(logLik(fit)),
    dispersion = sigma_value,
    gamma_shape = if (is_gamma) 1 / sigma_value^2 else NA_real_,
    convergence = fit$fit$convergence,
    pdHess = isTRUE(fit$sdr$pdHess)
  )
}

extract_interaction <- function(fit, habitat, within_term, model) {
  coef_table <- summary(fit)$coefficients$cond
  term <- interaction_term(fit, within_term)
  x <- coef_table[term, , drop = FALSE]
  
  tibble(
    habitat = habitat,
    model = model,
    family = family(fit)$family,
    link = family(fit)$link,
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
  
  pearson_ratio <- tryCatch({
    pearson <- residuals(fit, type = "pearson")
    sum(pearson^2, na.rm = TRUE) / max(1, nobs(fit) - length(beta))
  }, error = function(e) NA_real_)
  
  random_sd <- attr(VarCorr(fit)$cond$NAME10, "stddev")
  slope_sd <- unname(random_sd[within_term])
  
  vc <- vcov(fit)$cond
  cor_mat <- cov2cor(vc)
  diag(cor_mat) <- NA_real_
  
  tibble(
    habitat = habitat,
    family = family(fit)$family,
    link = family(fit)$link,
    n_checklists = nobs(fit),
    convergence = fit$fit$convergence,
    pdHess = isTRUE(fit$sdr$pdHess),
    dispersion = as.numeric(sigma(fit)),
    pearson_dispersion = pearson_ratio,
    city_random_slope_SD = slope_sd,
    max_fixed_effect_correlation = max(abs(cor_mat), na.rm = TRUE)
  )
}

expected_response <- function(fit, eta) {
  family(fit)$linkinv(eta)
}

make_curve <- function(fit, habitat, within_term, context_term, aridity_values) {
  beta <- fixef(fit)$cond
  vc <- vcov(fit)$cond
  term <- interaction_term(fit, within_term)
  centered <- aridity_values - mean_city_aridity
  frame <- fit$frame
  
  reference_protocol <- frame %>%
    count(protocol_name, sort = TRUE) %>%
    slice(1) %>%
    pull(protocol_name) %>%
    as.character()
  
  template <- frame[rep(1, length(aridity_values)), , drop = FALSE]
  template[[within_term]] <- 0
  template[[context_term]] <- 0
  template$aridity_centered <- centered
  template$duration_10 <- median(frame$duration_10, na.rm = TRUE)
  template$effort_distance_km <- median(frame$effort_distance_km, na.rm = TRUE)
  template$protocol_name <- factor(
    reference_protocol,
    levels = levels(frame$protocol_name)
  )
  
  baseline_link <- as.numeric(predict(
    fit,
    newdata = template,
    type = "link",
    re.form = NA,
    allow.new.levels = FALSE
  ))
  
  delta <- habitat_step * (beta[within_term] + beta[term] * centered)
  delta_var <- habitat_step^2 * (
    vc[within_term, within_term] +
      centered^2 * vc[term, term] +
      2 * centered * vc[within_term, term]
  )
  delta_se <- sqrt(pmax(delta_var, 0))
  
  baseline_response <- expected_response(fit, baseline_link)
  after_response <- expected_response(fit, baseline_link + delta)
  lower_response <- expected_response(fit, baseline_link + delta - 1.96 * delta_se)
  upper_response <- expected_response(fit, baseline_link + delta + 1.96 * delta_se)
  
  if (any(!is.finite(baseline_response))) {
    stop("Non-finite population FDis baseline encountered for ", habitat, ".")
  }
  
  tibble(
    habitat = habitat,
    family = family(fit)$family,
    link = family(fit)$link,
    avg_aridity_index = aridity_values,
    aridity_centered = centered,
    habitat_increment_percent_cover = habitat_step * 10,
    baseline_response = baseline_response,
    baseline_positive = baseline_response > 0,
    delta_eta = delta,
    delta_eta_SE = delta_se,
    relative_change_percent = if_else(
      baseline_response > 0,
      100 * (after_response / baseline_response - 1),
      NA_real_
    ),
    relative_CI_lower_percent = if_else(
      baseline_response > 0,
      100 * (lower_response / baseline_response - 1),
      NA_real_
    ),
    relative_CI_upper_percent = if_else(
      baseline_response > 0,
      100 * (upper_response / baseline_response - 1),
      NA_real_
    )
  )
}

city_context_reference <- city_aridity %>%
  left_join(city_habitat, by = "NAME10") %>%
  mutate(
    aridity_centered = avg_aridity_index - mean_city_aridity,
    city_context_water_10 = (city_mean_water_pp - mean_city_water) / 10,
    city_context_shrub_10 = (city_mean_shrub_pp - mean_city_shrub) / 10,
    city_context_built_10 = (city_mean_built_pp - mean_city_built) / 10
  )

make_city_effects <- function(fit, habitat, within_term, context_term) {
  beta <- fixef(fit)$cond
  term <- interaction_term(fit, within_term)
  frame <- fit$frame
  city_re <- ranef(fit)$cond$NAME10
  
  if (is.null(city_re) ||
      !all(c("(Intercept)", within_term) %in% colnames(city_re))) {
    stop("Could not recover city random intercept/slope for ", habitat, ".")
  }
  
  reference_protocol <- frame %>%
    count(protocol_name, sort = TRUE) %>%
    slice(1) %>%
    pull(protocol_name) %>%
    as.character()
  
  city_data <- city_context_reference %>%
    transmute(
      NAME10 = as.character(NAME10),
      avg_aridity_index,
      aridity_centered,
      context_value = .data[[context_term]]
    )
  
  template <- frame[rep(1, nrow(city_data)), , drop = FALSE]
  template[[within_term]] <- 0
  template[[context_term]] <- city_data$context_value
  template$aridity_centered <- city_data$aridity_centered
  template$duration_10 <- median(frame$duration_10, na.rm = TRUE)
  template$effort_distance_km <- median(frame$effort_distance_km, na.rm = TRUE)
  template$protocol_name <- factor(
    reference_protocol,
    levels = levels(frame$protocol_name)
  )
  template$NAME10 <- factor(
    city_data$NAME10,
    levels = levels(frame$NAME10)
  )
  template$locality_id <- factor(
    levels(frame$locality_id)[1],
    levels = levels(frame$locality_id)
  )
  
  fixed_baseline_link <- as.numeric(predict(
    fit,
    newdata = template,
    type = "link",
    re.form = NA,
    allow.new.levels = FALSE
  ))
  
  re_table <- tibble(
    NAME10 = rownames(city_re),
    random_intercept = as.numeric(city_re[, "(Intercept)"]),
    random_slope = as.numeric(city_re[, within_term])
  )
  
  city_data %>%
    mutate(fixed_baseline_link = fixed_baseline_link) %>%
    left_join(re_table, by = "NAME10") %>%
    mutate(
      habitat = habitat,
      family = family(fit)$family,
      link = family(fit)$link,
      fixed_slope = beta[within_term] + beta[term] * aridity_centered,
      city_slope = fixed_slope + random_slope,
      baseline_link = fixed_baseline_link + random_intercept,
      after_link = baseline_link + habitat_step * city_slope,
      baseline_response = expected_response(fit, baseline_link),
      after_response = expected_response(fit, after_link),
      relative_change_percent = if_else(
        baseline_response > 0,
        100 * (after_response / baseline_response - 1),
        NA_real_
      )
    ) %>%
    select(
      NAME10, habitat, family, link, avg_aridity_index, aridity_centered,
      random_intercept, random_slope, fixed_slope, city_slope,
      baseline_response, relative_change_percent
    )
}

# Model data

model_data <- fdis_analysis_data %>%
  select(
    FDis,
    water_within_model_10, shrub_within_model_10, built_within_model_10,
    aridity_centered,
    city_context_water_10, city_context_shrub_10, city_context_built_10,
    duration_10, effort_distance_km, protocol_name,
    NAME10, locality_id
  )

rm(fdis_analysis_data)
gc()

aridity_values <- seq(
  min(city_aridity$avg_aridity_index),
  max(city_aridity$avg_aridity_index),
  length.out = n_aridity_points
)

# Gaussian versus Gamma family comparison

candidate_registry <- family_info %>%
  filter(eligible) %>%
  select(family_key, family_label, family_name, link) %>%
  crossing(habitat_info)

fit_candidate_row <- function(
    family_key,
    family_label,
    family_name,
    link,
    habitat,
    stub,
    within_term,
    context_term
) {
  message("Family comparison: ", family_label, " - ", habitat)
  
  fit <- fit_model(
    model_data,
    within_term,
    context_term,
    family_key
  )
  
  out <- extract_candidate_summary(
    fit,
    family_key,
    family_label,
    habitat
  )
  
  if (!is.null(fit)) rm(fit)
  gc()
  out
}

candidate_results <- candidate_registry %>%
  pmap_dfr(fit_candidate_row)

family_selection_core <- candidate_results %>%
  group_by(family_key, family_label) %>%
  summarise(
    n_models = n(),
    all_converged = all(converged & pdHess),
    same_n = n_distinct(n_checklists) == 1,
    total_AIC = if (all(is.finite(AIC))) sum(AIC) else NA_real_,
    mean_AIC = if (all(is.finite(AIC))) mean(AIC) else NA_real_,
    total_BIC = if (all(is.finite(BIC))) sum(BIC) else NA_real_,
    .groups = "drop"
  )

family_selection <- family_info %>%
  select(family_key, family_label, family_name, link, eligible) %>%
  left_join(family_selection_core, by = c("family_key", "family_label")) %>%
  mutate(
    reason_not_eligible = case_when(
      eligible ~ NA_character_,
      family_key == "Gamma" ~ "FDis includes one or more zero values",
      TRUE ~ "Not eligible"
    )
  )

valid_families <- family_selection %>%
  filter(
    eligible,
    n_models == nrow(habitat_info),
    all_converged,
    same_n,
    is.finite(total_AIC)
  )

if (nrow(valid_families) == 0) {
  stop("Neither Gaussian nor Gamma produced a valid full set of primary FDis models.")
}

selected_family <- valid_families %>%
  slice_min(total_AIC, n = 1, with_ties = FALSE)

selected_family_key <- selected_family$family_key[[1]]
selected_family_label <- selected_family$family_label[[1]]
selected_family_name <- selected_family$family_name[[1]]
selected_link <- selected_family$link[[1]]

family_selection <- family_selection %>%
  mutate(
    delta_AIC = if_else(
      eligible & is.finite(total_AIC),
      total_AIC - min(valid_families$total_AIC),
      NA_real_
    ),
    selected = family_key == selected_family_key
  )

write_csv(
  candidate_results,
  fs::path(table_dir, "Functional_Dispersion_Family_Comparison_By_Habitat.csv")
)

write_csv(
  family_selection,
  fs::path(table_dir, "Functional_Dispersion_Family_Selection.csv")
)

message(
  "Selected FDis family: ", selected_family_label,
  " (aggregate AIC across Water, Shrub and scrub, and Built)"
)

# DHARMa diagnostics for selected family

diagnostic_stubs <- c(
  habitat_info$stub,
  paste0("Water_min_", water_min_range, "pp")
)

if (!identical(dharma_models, "all")) {
  invalid_names <- setdiff(dharma_models, diagnostic_stubs)
  if (length(invalid_names) > 0) {
    stop(
      "Unknown DHARMa model name(s): ",
      paste(invalid_names, collapse = ", ")
    )
  }
}

run_dharma_fit <- function(fit, habitat, stub, habitat_term) {
  if (!run_dharma) return(NULL)
  
  if (!identical(dharma_models, "all") && !stub %in% dharma_models) {
    return(NULL)
  }
  
  if (length(dharma_sample_n) != 1 || !is.finite(dharma_sample_n) ||
      dharma_sample_n < 100 || dharma_sample_n != as.integer(dharma_sample_n)) {
    stop("dharma_sample_n must be an integer of at least 100.")
  }
  
  if (length(dharma_n) != 1 || !is.finite(dharma_n) ||
      dharma_n < 20 || dharma_n != as.integer(dharma_n)) {
    stop("dharma_n must be an integer of at least 20.")
  }
  
  family_tag <- if (selected_family_key == "Gaussian") "Gaussian" else "Gamma"
  
  test_path <- fs::path(
    table_dir,
    paste0("Functional_Dispersion_DHARMa_", family_tag, "_", stub, "_Tests.csv")
  )
  residual_path <- fs::path(
    data_dir,
    paste0("Functional_Dispersion_DHARMa_", family_tag, "_", stub, "_Residuals.csv")
  )
  figure_path <- fs::path(
    figure_dir,
    paste0("Functional_Dispersion_DHARMa_", family_tag, "_", stub, ".png")
  )
  predictor_figure_path <- fs::path(
    figure_dir,
    paste0("Functional_Dispersion_DHARMa_", family_tag, "_", stub, "_Predictors.png")
  )
  object_path <- fs::path(
    data_dir,
    paste0("Functional_Dispersion_DHARMa_", family_tag, "_", stub, ".rds")
  )
  
  frame <- fit$frame
  if (is.null(frame) || nrow(frame) != nobs(fit)) {
    stop("The fitted model frame is unavailable: ", stub)
  }
  
  if (!all(c("FDis", "aridity_centered", habitat_term, "NAME10") %in%
           names(frame))) {
    stop("Required diagnostic columns are missing: ", stub)
  }
  
  set.seed(dharma_seed)
  sampled_rows <- sort(sample.int(
    nobs(fit),
    size = min(dharma_sample_n, nobs(fit))
  ))
  sampled_frame <- frame[sampled_rows, , drop = FALSE]
  
  predicted_response <- as.numeric(predict(
    fit,
    newdata = sampled_frame,
    type = "response",
    re.form = NULL,
    allow.new.levels = FALSE
  ))
  observed_response <- as.numeric(sampled_frame$FDis)
  
  if (length(predicted_response) != length(observed_response) ||
      any(!is.finite(predicted_response)) ||
      any(!is.finite(observed_response))) {
    stop("Invalid sampled predictions or responses: ", stub)
  }
  
  set.seed(dharma_seed)
  
  if (selected_family_key == "Gaussian") {
    residual_sd <- as.numeric(sigma(fit))
    
    if (length(residual_sd) != 1 || !is.finite(residual_sd) || residual_sd <= 0) {
      stop("The Gaussian model has no valid residual SD: ", stub)
    }
    
    simulated_response <- replicate(
      dharma_n,
      rnorm(
        length(predicted_response),
        mean = predicted_response,
        sd = residual_sd
      )
    )
    
    simulation_label <- "Conditional Gaussian"
    gamma_sigma <- NA_real_
    gamma_shape <- NA_real_
  } else {
    gamma_sigma <- as.numeric(sigma(fit))
    gamma_shape <- 1 / gamma_sigma^2
    
    if (length(gamma_sigma) != 1 || !is.finite(gamma_sigma) ||
        gamma_sigma <= 0 || !is.finite(gamma_shape) || gamma_shape <= 0 ||
        any(predicted_response <= 0)) {
      stop("The Gamma model has invalid conditional simulation parameters: ", stub)
    }
    
    simulated_response <- replicate(
      dharma_n,
      rgamma(
        length(predicted_response),
        shape = gamma_shape,
        scale = predicted_response / gamma_shape
      )
    )
    
    simulation_label <- "Conditional Gamma"
  }
  
  message(
    "DHARMa: ", stub, " - ", length(sampled_rows), " of ", nobs(fit),
    " checklists, ", dharma_n, " simulations"
  )
  
  sim <- DHARMa::createDHARMa(
    simulatedResponse = simulated_response,
    observedResponse = observed_response,
    fittedPredictedResponse = predicted_response,
    integerResponse = FALSE,
    seed = dharma_seed
  )
  
  tests <- list(
    Uniformity = DHARMa::testUniformity(sim, plot = FALSE),
    Dispersion = DHARMa::testDispersion(sim, plot = FALSE),
    Outliers = DHARMa::testOutliers(sim, type = "binomial", plot = FALSE)
  )
  
  dharma_version <- as.character(utils::packageVersion("DHARMa"))
  
  test_results <- imap_dfr(tests, function(result, test_name) {
    tibble(
      habitat = habitat,
      stub = stub,
      selected_family = selected_family_name,
      selected_link = selected_link,
      test = test_name,
      statistic = as.numeric(result$statistic)[1],
      p = as.numeric(result$p.value)[1],
      n_checklists = nobs(fit),
      n_sampled = length(sampled_rows),
      n_simulations = dharma_n,
      simulation = simulation_label,
      gamma_sigma = gamma_sigma,
      gamma_shape = gamma_shape,
      outlier_method = "binomial",
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
  write_rds(sim, object_path, compress = "gz")
  
  grDevices::png(figure_path, width = 1600, height = 800, res = 150)
  tryCatch({
    graphics::par(mfrow = c(1, 2))
    DHARMa::plotQQunif(sim, testOutliers = FALSE)
    DHARMa::plotResiduals(sim, quantreg = FALSE)
  }, finally = grDevices::dev.off())
  
  grDevices::png(
    predictor_figure_path,
    width = 1600,
    height = 800,
    res = 150
  )
  tryCatch({
    graphics::par(mfrow = c(1, 2))
    DHARMa::plotResiduals(
      sim,
      form = diagnostic_data$aridity_centered,
      quantreg = FALSE,
      rank = FALSE,
      xlab = "Centered aridity"
    )
    DHARMa::plotResiduals(
      sim,
      form = diagnostic_data$habitat_within_model_10,
      quantreg = FALSE,
      rank = FALSE,
      xlab = "Within-city habitat (per 10 pp)"
    )
  }, finally = grDevices::dev.off())
  
  rm(frame, sampled_frame, sim, simulated_response, diagnostic_data)
  gc()
  
  test_results
}

# Selected primary models

process_selected_habitat <- function(habitat, stub, within_term, context_term) {
  message("Primary selected family: ", selected_family_label, " - ", habitat)
  
  fit <- fit_model(
    model_data,
    within_term,
    context_term,
    selected_family_key
  )
  
  if (is.null(fit)) {
    stop("Selected-family primary FDis model failed to converge: ", habitat)
  }
  
  model_label <- paste0("Primary ", selected_family_label)
  
  out <- list(
    interaction = extract_interaction(fit, habitat, within_term, model_label),
    qc = extract_qc(fit, habitat, within_term),
    curve = make_curve(fit, habitat, within_term, context_term, aridity_values),
    city_effects = make_city_effects(fit, habitat, within_term, context_term),
    coefficients = extract_coefficients(fit, habitat, model_label),
    summary = extract_model_summary(fit, habitat, model_label),
    dharma = run_dharma_fit(fit, habitat, stub, within_term)
  )
  
  rm(fit)
  gc()
  out
}

primary_outputs <- pmap(habitat_info, process_selected_habitat)

interaction_results <- map_dfr(primary_outputs, "interaction")
model_qc <- map_dfr(primary_outputs, "qc")
resource_curves <- map_dfr(primary_outputs, "curve")
city_effects <- map_dfr(primary_outputs, "city_effects")
model_coefficients <- map_dfr(primary_outputs, "coefficients")
model_summaries <- map_dfr(primary_outputs, "summary")
dharma_results <- bind_rows(map(primary_outputs, "dharma"))

model_registry <- habitat_info %>%
  mutate(
    response = "FDis",
    family = selected_family_name,
    link = selected_link,
    full_model_saved = FALSE,
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
    NAME10 = forcats::fct_drop(NAME10),
    locality_id = forcats::fct_drop(locality_id),
    protocol_name = forcats::fct_drop(protocol_name)
  )

water_info <- habitat_info %>%
  filter(habitat == "Water")

message(
  "Water support sensitivity: ", selected_family_label,
  " - range >= ", water_min_range, " pp"
)

water_support_fit <- fit_model(
  water_data,
  water_info$within_term,
  water_info$context_term,
  selected_family_key
)

if (is.null(water_support_fit)) {
  stop("Selected-family Water support model failed to converge.")
}

water_sensitivity <- bind_rows(
  interaction_results %>%
    filter(habitat == "Water") %>%
    mutate(
      model = paste0("Primary ", selected_family_label, ": all cities"),
      n_cities = n_distinct(city_aridity$NAME10),
      n_localities = sample_summary$n_localities,
      n_checklists = sample_summary$n_checklists,
      .before = family
    ),
  extract_interaction(
    water_support_fit,
    "Water",
    water_info$within_term,
    paste0(
      selected_family_label,
      ": Water range >= ", water_min_range, " pp"
    )
  ) %>%
    mutate(
      n_cities = n_distinct(water_data$NAME10),
      n_localities = n_distinct(water_data$locality_id),
      n_checklists = nrow(water_data),
      .before = family
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

model_coefficients <- bind_rows(
  model_coefficients,
  extract_coefficients(
    water_support_fit,
    "Water support",
    paste0("Water support ", selected_family_label)
  )
)

model_summaries <- bind_rows(
  model_summaries,
  extract_model_summary(
    water_support_fit,
    "Water support",
    paste0("Water support ", selected_family_label)
  )
)

water_dharma <- run_dharma_fit(
  water_support_fit,
  "Water support",
  paste0("Water_min_", water_min_range, "pp"),
  water_info$within_term
)

dharma_results <- bind_rows(dharma_results, water_dharma)

rm(water_support_fit, water_data, model_data)
gc()

# Figures

habitat_levels <- c("Water", "Shrub and scrub", "Built")

plot_curves <- resource_curves %>%
  mutate(habitat = factor(habitat, levels = habitat_levels))

plot_city_effects <- city_effects %>%
  mutate(habitat = factor(habitat, levels = habitat_levels))

if (any(!plot_curves$baseline_positive)) {
  warning(
    "At least one selected-family population baseline prediction was <= 0. ",
    "Relative percentage effects are NA at those aridity values."
  )
}

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
  fs::path(table_dir, "Functional_Dispersion_Habitat_x_Aridity.csv")
)
write_csv(
  model_qc,
  fs::path(table_dir, "Functional_Dispersion_Model_QC.csv")
)
write_csv(
  model_coefficients,
  fs::path(table_dir, "Functional_Dispersion_All_Model_Coefficients.csv")
)
write_csv(
  model_summaries,
  fs::path(table_dir, "Functional_Dispersion_All_Model_Summaries.csv")
)
write_csv(
  resource_curves,
  fs::path(table_dir, "Functional_Dispersion_Effect_Curves.csv")
)
write_csv(
  city_effects,
  fs::path(table_dir, "Functional_Dispersion_City_Conditional_Effects.csv")
)
write_csv(
  water_sensitivity,
  fs::path(table_dir, "Functional_Dispersion_Water_Support_Sensitivity.csv")
)
write_csv(
  water_support_qc,
  fs::path(table_dir, "Functional_Dispersion_Water_Support_Model_QC.csv")
)
write_csv(
  water_support,
  fs::path(table_dir, "Functional_Dispersion_Water_Support_By_City.csv")
)
write_csv(
  sample_summary,
  fs::path(table_dir, "Functional_Dispersion_Sample_Summary.csv")
)

if (run_dharma && nrow(dharma_results) > 0) {
  write_csv(
    dharma_results,
    fs::path(table_dir, "Functional_Dispersion_DHARMa_Tests.csv")
  )
}

model_bundle <- list(
  response = "Functional dispersion",
  raw_response_variable = "FDis",
  modeled_response_variable = "FDis",
  model_family = selected_family_name,
  model_family_label = selected_family_label,
  link = selected_link,
  family_selection_method = paste0(
    "Lowest aggregate AIC across the three primary habitat models; ",
    "candidate family required convergence == 0 and pdHess == TRUE for all three models"
  ),
  gamma_eligible = gamma_eligible,
  family_selection = family_selection,
  family_comparison_by_habitat = candidate_results,
  primary_response_offset = "None",
  city_baseline_structure = "Random intercept for NAME10",
  mean_city_aridity = mean_city_aridity,
  full_model_objects_saved = FALSE,
  primary_model_paths = character(0),
  model_registry = model_registry,
  primary_interactions = interaction_results,
  primary_model_qc = model_qc,
  primary_resource_curves = resource_curves
)

write_rds(
  model_bundle,
  fs::path(data_dir, "Functional_Dispersion_Resource_Importance_Models.rds"),
  compress = "gz"
)

write_rds(
  list(
    selected_family = selected_family_label,
    selected_family_name = selected_family_name,
    selected_link = selected_link,
    gamma_eligible = gamma_eligible,
    family_selection = family_selection,
    family_comparison_by_habitat = candidate_results,
    selected_interactions = interaction_results,
    selected_resource_curves = resource_curves,
    selected_model_qc = model_qc,
    water_support_sensitivity = water_sensitivity,
    water_support_model_qc = water_support_qc,
    dharma_results = dharma_results,
    full_model_objects_saved = FALSE
  ),
  fs::path(data_dir, "Functional_Dispersion_Diagnostic_Bundle.rds"),
  compress = "gz"
)

ggsave(
  fs::path(figure_dir, "Figure_01_Functional_Dispersion_Resource_Importance.png"),
  resource_figure,
  width = 12,
  height = 4.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  fs::path(figure_dir, "Figure_01_Functional_Dispersion_Resource_Importance.pdf"),
  resource_figure,
  width = 12,
  height = 4.5,
  units = "in",
  bg = "white"
)

ggsave(
  fs::path(figure_dir, "Figure_02_Functional_Dispersion_Habitat_x_Aridity.png"),
  interaction_figure,
  width = 7,
  height = 4,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  fs::path(figure_dir, "Figure_02_Functional_Dispersion_Habitat_x_Aridity.pdf"),
  interaction_figure,
  width = 7,
  height = 4,
  units = "in",
  bg = "white"
)

message(
  "Script 04 complete. Selected family: ", selected_family_label,
  ". Full fitted glmmTMB objects were not saved."
)
