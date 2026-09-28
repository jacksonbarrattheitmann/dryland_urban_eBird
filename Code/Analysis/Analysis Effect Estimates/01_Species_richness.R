# ============================================================
# Dryland bird project
# 01_Species_richness habitat collinearity and city effect sizes
#
# Workflow:
# 1. Apply checklist/locality/city eligibility filters.
# 2. Quantify water-shrub-built collinearity within every city.
# 3. Use separate habitat models consistently across all cities.
# 4. Extract one city-specific coefficient for each habitat.
# 5. Relate each city-specific coefficient to city aridity.
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(glmmTMB)
library(ggpubr)

TEMPORARY_DEDUPLICATE <- TRUE

MIN_CHECKLISTS_PER_LOCALITY <- 10
MIN_LOCALITIES_PER_CITY <- 5

DE_MARTONNE_DRYLAND_MAX <- 30

OUTPUT_ROOT <- file.path(
  "Outputs",
  "Analysis_City_Effect_Sizes"
)

dir.create(
  OUTPUT_ROOT,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

data_path <- file.path(
  "Intermediate_data",
  "dfs_for_analysis_combined",
  "combined_data.RDS"
)

if (!file.exists(data_path)) {
  stop("combined_data.RDS was not found at: ", data_path)
}

combined_data <- readRDS(data_path)

required_variables <- c(
  "checklist_id",
  "locality_id",
  "NAME10",
  "state_name",
  "protocol_name",
  "duration_minutes",
  "effort_distance_km",
  "sr",
  "avg_aridity_index",
  "water",
  "shrub_and_scrub",
  "built"
)

missing_variables <- setdiff(
  required_variables,
  names(combined_data)
)

if (length(missing_variables) > 0) {
  stop(
    "Missing required variables: ",
    paste(missing_variables, collapse = ", ")
  )
}


# ------------------------------------------------------------
# Checklist filters
# ------------------------------------------------------------

duplicate_checklists <- combined_data %>%
  count(
    checklist_id,
    name = "n_rows"
  ) %>%
  filter(
    n_rows > 1
  )

if (
  nrow(duplicate_checklists) > 0 &&
  !TEMPORARY_DEDUPLICATE
) {
  stop("Duplicate checklist IDs are present.")
}

if (
  nrow(duplicate_checklists) > 0 &&
  TEMPORARY_DEDUPLICATE
) {
  analysis_base <- combined_data %>%
    distinct(
      checklist_id,
      .keep_all = TRUE
    )
} else {
  analysis_base <- combined_data
}

sr_pool <- analysis_base %>%
  mutate(
    effort_distance_km = if_else(
      protocol_name == "Stationary",
      0,
      effort_distance_km
    )
  ) %>%
  filter(
    protocol_name %in% c(
      "Stationary",
      "Traveling"
    ),
    duration_minutes <= 300,
    effort_distance_km <= 5
  )

# ------------------------------------------------------------
# Restrict analysis to De Martonne dryland conditions
# ------------------------------------------------------------

excluded_humid_cities <- sr_pool %>%
  distinct(
    NAME10,
    state_name,
    avg_aridity_index
  ) %>%
  filter(
    avg_aridity_index >= DE_MARTONNE_DRYLAND_MAX
  ) %>%
  arrange(
    desc(avg_aridity_index)
  )

sr_pool <- sr_pool %>%
  filter(
    avg_aridity_index < DE_MARTONNE_DRYLAND_MAX
  )


model_variables <- c(
  "sr",
  "avg_aridity_index",
  "water",
  "shrub_and_scrub",
  "built",
  "duration_minutes",
  "effort_distance_km",
  "NAME10",
  "locality_id"
)

missing_summary <- data.frame(
  variable = model_variables,
  missing_n = sapply(
    sr_pool[model_variables],
    function(x) sum(is.na(x))
  ),
  stringsAsFactors = FALSE
)

missing_summary$missing_percent <-
  100 *
  missing_summary$missing_n /
  nrow(sr_pool)

sr_pool <- sr_pool %>%
  filter(
    if_all(
      all_of(model_variables),
      ~ !is.na(.)
    )
  )


# ------------------------------------------------------------
# Spatial-grain checks
# ------------------------------------------------------------

aridity_qc <- sr_pool %>%
  group_by(NAME10) %>%
  summarise(
    n_aridity_values =
      n_distinct(avg_aridity_index),
    .groups = "drop"
  )

if (any(aridity_qc$n_aridity_values != 1)) {
  stop(
    "At least one city has more than one aridity value."
  )
}

habitat_qc <- sr_pool %>%
  group_by(
    NAME10,
    locality_id
  ) %>%
  summarise(
    n_water =
      n_distinct(water),
    n_shrub =
      n_distinct(shrub_and_scrub),
    n_built =
      n_distinct(built),
    .groups = "drop"
  )

if (
  any(habitat_qc$n_water > 1) ||
  any(habitat_qc$n_shrub > 1) ||
  any(habitat_qc$n_built > 1)
) {
  stop(
    paste0(
      "At least one locality has multiple values ",
      "for a habitat predictor."
    )
  )
}


# ------------------------------------------------------------
# Scale habitat variables
# ------------------------------------------------------------

maximum_habitat_value <- max(
  c(
    sr_pool$water,
    sr_pool$shrub_and_scrub,
    sr_pool$built
  ),
  na.rm = TRUE
)

if (maximum_habitat_value <= 1.05) {
  habitat_divisor <- 0.10
} else {
  habitat_divisor <- 10
}

sr_pool <- sr_pool %>%
  mutate(
    water_10 =
      water / habitat_divisor,
    shrub_10 =
      shrub_and_scrub / habitat_divisor,
    built_10 =
      built / habitat_divisor,
    duration_10 =
      duration_minutes / 10
  )


# ------------------------------------------------------------
# Locality and city eligibility
# ------------------------------------------------------------

locality_support <- sr_pool %>%
  count(
    NAME10,
    locality_id,
    name = "n_checklists"
  ) %>%
  mutate(
    locality_eligible =
      n_checklists >=
      MIN_CHECKLISTS_PER_LOCALITY
  )

eligible_localities <- locality_support %>%
  filter(
    locality_eligible
  )

city_total_support <- sr_pool %>%
  group_by(NAME10) %>%
  summarise(
    state_name =
      first(state_name),
    avg_aridity_index =
      first(avg_aridity_index),
    n_total_checklists =
      n(),
    n_total_localities =
      n_distinct(locality_id),
    .groups = "drop"
  )

city_eligible_support <- eligible_localities %>%
  group_by(NAME10) %>%
  summarise(
    n_eligible_localities =
      n(),
    n_retained_checklists =
      sum(n_checklists),
    .groups = "drop"
  )

city_qc <- city_total_support %>%
  left_join(
    city_eligible_support,
    by = "NAME10"
  ) %>%
  mutate(
    n_eligible_localities =
      coalesce(
        n_eligible_localities,
        0L
      ),
    n_retained_checklists =
      coalesce(
        n_retained_checklists,
        0L
      ),
    city_eligible =
      n_eligible_localities >=
      MIN_LOCALITIES_PER_CITY
  )

eligible_cities <- city_qc %>%
  filter(
    city_eligible
  )

sr_stage1 <- sr_pool %>%
  semi_join(
    eligible_localities,
    by = c(
      "NAME10",
      "locality_id"
    )
  ) %>%
  semi_join(
    eligible_cities,
    by = "NAME10"
  )

stage1_locality_qc <- sr_stage1 %>%
  count(
    NAME10,
    locality_id,
    name = "n_checklists"
  )

stage1_city_qc <- stage1_locality_qc %>%
  count(
    NAME10,
    name = "n_localities"
  )

if (
  any(
    stage1_locality_qc$n_checklists <
    MIN_CHECKLISTS_PER_LOCALITY
  )
) {
  stop("Locality eligibility QC failed.")
}

if (
  any(
    stage1_city_qc$n_localities <
    MIN_LOCALITIES_PER_CITY
  )
) {
  stop("City eligibility QC failed.")
}


# ------------------------------------------------------------
# Collinearity QC
# Habitat correlations are calculated within each city
# using one row per locality.
# ------------------------------------------------------------

safe_cor <- function(x, y) {
  
  if (
    length(unique(x)) < 2 ||
    length(unique(y)) < 2
  ) {
    return(NA_real_)
  }
  
  cor(
    x,
    y,
    use = "complete.obs",
    method = "pearson"
  )
}

city_habitat_qc_list <- lapply(
  split(
    sr_stage1,
    sr_stage1$NAME10
  ),
  function(city_data) {
    
    locality_data <- city_data %>%
      distinct(
        locality_id,
        water_10,
        shrub_10,
        built_10
      )
    
    r_water_shrub <- safe_cor(
      locality_data$water_10,
      locality_data$shrub_10
    )
    
    r_water_built <- safe_cor(
      locality_data$water_10,
      locality_data$built_10
    )
    
    r_shrub_built <- safe_cor(
      locality_data$shrub_10,
      locality_data$built_10
    )
    
    abs_correlations <- abs(
      c(
        r_water_shrub,
        r_water_built,
        r_shrub_built
      )
    )
    
    max_abs_habitat_cor <- if (
      all(is.na(abs_correlations))
    ) {
      NA_real_
    } else {
      max(
        abs_correlations,
        na.rm = TRUE
      )
    }
    
    data.frame(
      NAME10 =
        first(city_data$NAME10),
      state_name =
        first(city_data$state_name),
      avg_aridity_index =
        first(city_data$avg_aridity_index),
      n_checklists =
        nrow(city_data),
      n_localities =
        nrow(locality_data),
      r_water_shrub =
        r_water_shrub,
      r_water_built =
        r_water_built,
      r_shrub_built =
        r_shrub_built,
      max_abs_habitat_cor =
        max_abs_habitat_cor,
      water_range =
        diff(range(locality_data$water_10)),
      shrub_range =
        diff(range(locality_data$shrub_10)),
      built_range =
        diff(range(locality_data$built_10)),
      water_sd =
        sd(locality_data$water_10),
      shrub_sd =
        sd(locality_data$shrub_10),
      built_sd =
        sd(locality_data$built_10),
      stringsAsFactors = FALSE
    )
  }
)

city_habitat_qc <- bind_rows(
  city_habitat_qc_list
)

correlation_long <- city_habitat_qc %>%
  select(
    NAME10,
    r_water_shrub,
    r_water_built,
    r_shrub_built
  ) %>%
  pivot_longer(
    cols = starts_with("r_"),
    names_to = "comparison",
    values_to = "r"
  ) %>%
  mutate(
    comparison = recode(
      comparison,
      r_water_shrub =
        "Water vs shrub and scrub",
      r_water_built =
        "Water vs built",
      r_shrub_built =
        "Shrub and scrub vs built"
    ),
    abs_r = abs(r)
  )

correlation_pair_summary <- correlation_long %>%
  group_by(comparison) %>%
  summarise(
    n_cities =
      sum(!is.na(r)),
    median_r =
      median(r, na.rm = TRUE),
    median_abs_r =
      median(abs_r, na.rm = TRUE),
    mean_abs_r =
      mean(abs_r, na.rm = TRUE),
    q25_abs_r =
      quantile(
        abs_r,
        0.25,
        na.rm = TRUE
      ),
    q75_abs_r =
      quantile(
        abs_r,
        0.75,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

thresholds <- c(
  0.70,
  0.80,
  0.90,
  0.95
)

threshold_summary <- bind_rows(
  lapply(
    thresholds,
    function(threshold) {
      
      valid <- !is.na(
        city_habitat_qc$max_abs_habitat_cor
      )
      
      n_valid <- sum(valid)
      
      n_above <- sum(
        city_habitat_qc$max_abs_habitat_cor[
          valid
        ] >= threshold
      )
      
      data.frame(
        threshold =
          threshold,
        n_cities =
          n_valid,
        n_cities_at_or_above =
          n_above,
        percent_cities_at_or_above =
          100 * n_above / n_valid,
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  locality_support,
  file.path(
    OUTPUT_ROOT,
    "Table_00_Locality_Eligibility.csv"
  ),
  row.names = FALSE
)

write.csv(
  city_habitat_qc,
  file.path(
    OUTPUT_ROOT,
    "Table_01_City_Habitat_Collinearity.csv"
  ),
  row.names = FALSE
)

write.csv(
  correlation_pair_summary,
  file.path(
    OUTPUT_ROOT,
    "Table_01B_Pairwise_Correlation_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  threshold_summary,
  file.path(
    OUTPUT_ROOT,
    "Table_01C_Collinearity_Threshold_Summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Collinearity figures
# ------------------------------------------------------------

figure_qc1 <- ggplot(
  correlation_long,
  aes(
    x = abs_r
  )
) +
  geom_histogram(
    bins = 30
  ) +
  geom_vline(
    xintercept = c(
      0.7,
      0.8,
      0.9
    ),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  facet_wrap(
    ~ comparison,
    ncol = 1
  ) +
  labs(
    x = "Absolute within-city Pearson correlation",
    y = "Number of cities"
  ) +
  theme_classic(
    base_size = 12
  )

ggsave(
  file.path(
    OUTPUT_ROOT,
    "Figure_QC1_Within_City_Habitat_Correlations.png"
  ),
  figure_qc1,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

figure_qc2 <- ggplot(
  city_habitat_qc,
  aes(
    x = max_abs_habitat_cor
  )
) +
  geom_histogram(
    bins = 30
  ) +
  geom_vline(
    xintercept = c(
      0.7,
      0.8,
      0.9
    ),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  labs(
    x = "Maximum absolute habitat correlation within city",
    y = "Number of cities"
  ) +
  theme_classic(
    base_size = 12
  )

ggsave(
  file.path(
    OUTPUT_ROOT,
    "Figure_QC2_Maximum_Within_City_Habitat_Correlation.png"
  ),
  figure_qc2,
  width = 7,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ============================================================
# Stage 1
# Separate habitat models for every eligible city
# ============================================================

model_is_usable <- function(fit) {
  
  if (inherits(fit, "error")) {
    return(FALSE)
  }
  
  convergence_ok <-
    isTRUE(
      fit$fit$convergence == 0
    )
  
  hessian_ok <- tryCatch(
    isTRUE(fit$sdr$pdHess),
    error = function(e) FALSE
  )
  
  convergence_ok && hessian_ok
}

habitat_info <- data.frame(
  habitat = c(
    "Water",
    "Shrub and scrub",
    "Built"
  ),
  term = c(
    "water_10",
    "shrub_10",
    "built_10"
  ),
  stringsAsFactors = FALSE
)

fit_city_habitat_model <- function(
    city_name,
    habitat_name,
    habitat_term
) {
  
  city_data <- sr_stage1 %>%
    filter(
      NAME10 == city_name
    ) %>%
    mutate(
      locality_id =
        factor(locality_id)
    )
  
  city_meta <- city_habitat_qc %>%
    filter(
      NAME10 == city_name
    )
  
  habitat_vector <- city_data[[habitat_term]]
  
  if (
    length(unique(habitat_vector)) < 2
  ) {
    
    return(
      data.frame(
        NAME10 =
          city_name,
        state_name =
          first(city_data$state_name),
        avg_aridity_index =
          first(city_data$avg_aridity_index),
        habitat =
          habitat_name,
        term =
          habitat_term,
        beta =
          NA_real_,
        SE =
          NA_real_,
        z =
          NA_real_,
        p =
          NA_real_,
        CI_lower =
          NA_real_,
        CI_upper =
          NA_real_,
        n_checklists =
          nrow(city_data),
        n_localities =
          n_distinct(city_data$locality_id),
        habitat_range =
          diff(range(habitat_vector)),
        max_abs_habitat_cor =
          city_meta$max_abs_habitat_cor,
        model_status =
          "No habitat variation",
        fit_attempt =
          NA_character_,
        stringsAsFactors = FALSE
      )
    )
  }
  
  model_formula <- as.formula(
    paste(
      "sr ~",
      habitat_term,
      "+ duration_10 +",
      "effort_distance_km +",
      "(1 | locality_id)"
    )
  )
  
  fit_default <- tryCatch(
    glmmTMB(
      formula =
        model_formula,
      family =
        nbinom2(link = "log"),
      data =
        city_data
    ),
    error = function(e) e
  )
  
  fit <- fit_default
  fit_attempt <- "Default optimizer"
  
  if (!model_is_usable(fit_default)) {
    
    fit_bfgs <- tryCatch(
      glmmTMB(
        formula =
          model_formula,
        family =
          nbinom2(link = "log"),
        data =
          city_data,
        control =
          glmmTMBControl(
            optimizer = optim,
            optArgs = list(
              method = "BFGS"
            )
          )
      ),
      error = function(e) e
    )
    
    if (model_is_usable(fit_bfgs)) {
      fit <- fit_bfgs
      fit_attempt <- "BFGS retry"
    } else if (
      !inherits(fit_bfgs, "error")
    ) {
      fit <- fit_bfgs
      fit_attempt <- "BFGS retry - flagged"
    } else {
      fit_attempt <- "Model failed"
    }
  }
  
  if (inherits(fit, "error")) {
    
    return(
      data.frame(
        NAME10 =
          city_name,
        state_name =
          first(city_data$state_name),
        avg_aridity_index =
          first(city_data$avg_aridity_index),
        habitat =
          habitat_name,
        term =
          habitat_term,
        beta =
          NA_real_,
        SE =
          NA_real_,
        z =
          NA_real_,
        p =
          NA_real_,
        CI_lower =
          NA_real_,
        CI_upper =
          NA_real_,
        n_checklists =
          nrow(city_data),
        n_localities =
          n_distinct(city_data$locality_id),
        habitat_range =
          diff(range(habitat_vector)),
        max_abs_habitat_cor =
          city_meta$max_abs_habitat_cor,
        model_status =
          "Model failed",
        fit_attempt =
          fit_attempt,
        stringsAsFactors = FALSE
      )
    )
  }
  
  convergence_code <-
    fit$fit$convergence
  
  pd_hessian <- tryCatch(
    isTRUE(fit$sdr$pdHess),
    error = function(e) FALSE
  )
  
  if (
    convergence_code != 0 ||
    !pd_hessian
  ) {
    
    return(
      data.frame(
        NAME10 =
          city_name,
        state_name =
          first(city_data$state_name),
        avg_aridity_index =
          first(city_data$avg_aridity_index),
        habitat =
          habitat_name,
        term =
          habitat_term,
        beta =
          NA_real_,
        SE =
          NA_real_,
        z =
          NA_real_,
        p =
          NA_real_,
        CI_lower =
          NA_real_,
        CI_upper =
          NA_real_,
        n_checklists =
          nrow(city_data),
        n_localities =
          n_distinct(city_data$locality_id),
        habitat_range =
          diff(range(habitat_vector)),
        max_abs_habitat_cor =
          city_meta$max_abs_habitat_cor,
        model_status =
          "Convergence problem",
        fit_attempt =
          fit_attempt,
        stringsAsFactors = FALSE
      )
    )
  }
  
  coefficient_table <- as.data.frame(
    summary(fit)$coefficients$cond
  )
  
  coefficient_table$term <-
    rownames(coefficient_table)
  
  row_index <- match(
    habitat_term,
    coefficient_table$term
  )
  
  if (is.na(row_index)) {
    
    return(
      data.frame(
        NAME10 =
          city_name,
        state_name =
          first(city_data$state_name),
        avg_aridity_index =
          first(city_data$avg_aridity_index),
        habitat =
          habitat_name,
        term =
          habitat_term,
        beta =
          NA_real_,
        SE =
          NA_real_,
        z =
          NA_real_,
        p =
          NA_real_,
        CI_lower =
          NA_real_,
        CI_upper =
          NA_real_,
        n_checklists =
          nrow(city_data),
        n_localities =
          n_distinct(city_data$locality_id),
        habitat_range =
          diff(range(habitat_vector)),
        max_abs_habitat_cor =
          city_meta$max_abs_habitat_cor,
        model_status =
          "Coefficient not returned",
        fit_attempt =
          fit_attempt,
        stringsAsFactors = FALSE
      )
    )
  }
  
  beta <-
    coefficient_table$Estimate[row_index]
  
  se <-
    coefficient_table$`Std. Error`[row_index]
  
  z_value <-
    coefficient_table$`z value`[row_index]
  
  p_value <-
    coefficient_table$`Pr(>|z|)`[row_index]
  
  coefficient_ok <-
    is.finite(beta) &&
    is.finite(se) &&
    se > 0
  
  data.frame(
    NAME10 =
      city_name,
    state_name =
      first(city_data$state_name),
    avg_aridity_index =
      first(city_data$avg_aridity_index),
    habitat =
      habitat_name,
    term =
      habitat_term,
    beta =
      beta,
    SE =
      se,
    z =
      z_value,
    p =
      p_value,
    CI_lower =
      beta - 1.96 * se,
    CI_upper =
      beta + 1.96 * se,
    n_checklists =
      nrow(city_data),
    n_localities =
      n_distinct(city_data$locality_id),
    habitat_range =
      diff(range(habitat_vector)),
    max_abs_habitat_cor =
      city_meta$max_abs_habitat_cor,
    model_status =
      ifelse(
        coefficient_ok,
        "Usable",
        "Coefficient problem"
      ),
    fit_attempt =
      fit_attempt,
    stringsAsFactors = FALSE
  )
}


# ------------------------------------------------------------
# Fit every city x habitat combination
# ------------------------------------------------------------

city_names <- sort(
  unique(
    sr_stage1$NAME10
  )
)

model_results <- list()
counter <- 0

for (
  i in seq_along(city_names)
) {
  
  city_name <- city_names[i]
  
  cat(
    "[",
    i,
    "/",
    length(city_names),
    "] ",
    city_name,
    "\n",
    sep = ""
  )
  
  for (
    h in seq_len(
      nrow(habitat_info)
    )
  ) {
    
    counter <- counter + 1
    
    model_results[[counter]] <-
      fit_city_habitat_model(
        city_name =
          city_name,
        habitat_name =
          habitat_info$habitat[h],
        habitat_term =
          habitat_info$term[h]
      )
  }
}

city_effect_sizes_all <- bind_rows(
  model_results
)

city_effect_sizes_usable <- city_effect_sizes_all %>%
  filter(
    model_status == "Usable"
  ) %>%
  mutate(
    habitat = factor(
      habitat,
      levels = c(
        "Water",
        "Shrub and scrub",
        "Built"
      )
    )
  )

stage1_summary <- city_effect_sizes_all %>%
  group_by(habitat) %>%
  summarise(
    n_attempted =
      n(),
    n_usable =
      sum(model_status == "Usable"),
    n_failed =
      sum(model_status != "Usable"),
    median_beta =
      median(
        beta[model_status == "Usable"],
        na.rm = TRUE
      ),
    median_SE =
      median(
        SE[model_status == "Usable"],
        na.rm = TRUE
      ),
    .groups = "drop"
  )

write.csv(
  city_effect_sizes_all,
  file.path(
    OUTPUT_ROOT,
    "Table_02_City_Effect_Sizes_All.csv"
  ),
  row.names = FALSE
)

write.csv(
  city_effect_sizes_usable,
  file.path(
    OUTPUT_ROOT,
    "Table_02B_City_Effect_Sizes_Usable.csv"
  ),
  row.names = FALSE
)

write.csv(
  stage1_summary,
  file.path(
    OUTPUT_ROOT,
    "Table_02C_Stage1_Model_Summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Wide city effect-size table
# ------------------------------------------------------------

city_effect_sizes_wide <- city_effect_sizes_usable %>%
  select(
    NAME10,
    state_name,
    avg_aridity_index,
    habitat,
    beta,
    SE,
    p
  ) %>%
  mutate(
    habitat_short = recode(
      as.character(habitat),
      "Water" = "water",
      "Shrub and scrub" = "shrub",
      "Built" = "built"
    )
  ) %>%
  select(
    -habitat
  ) %>%
  pivot_wider(
    names_from =
      habitat_short,
    values_from = c(
      beta,
      SE,
      p
    ),
    names_glue =
      "{habitat_short}_{.value}"
  )

write.csv(
  city_effect_sizes_wide,
  file.path(
    OUTPUT_ROOT,
    "Table_02D_City_Effect_Sizes_Wide.csv"
  ),
  row.names = FALSE
)


























# ------------------------------------------------------------
# Figure: distributions of Stage-1 city-level beta estimates
# Central 95% shown for visualization only
# NOTE THIS IS VISUALIZING 95% OF THE DATA
# ------------------------------------------------------------

make_beta_histogram <- function(
    model_data,
    habitat_name
) {
  
  plot_data <- model_data %>%
    filter(
      habitat == habitat_name,
      is.finite(beta)
    )
  
  # Central 95% used only for plotting limits
  beta_limits <- quantile(
    plot_data$beta,
    probs = c(
      0.025,
      0.975
    ),
    na.rm = TRUE
  )
  
  # Bin width based on the visible range
  beta_binwidth <- diff(beta_limits) / 20
  
  ggplot(
    plot_data,
    aes(
      x = beta
    )
  ) +
    
    geom_histogram(
      binwidth = beta_binwidth,
      boundary = 0,
      fill = "grey70",
      color = "black",
      linewidth = 0.3
    ) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    
    coord_cartesian(
      xlim = beta_limits
    ) +
    
    labs(
      x = "City-specific beta estimate",
      y = "Number of cities",
      title = habitat_name
    ) +
    
    theme_classic(
      base_size = 12
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )
}


# ------------------------------------------------------------
# Create individual panels
# ------------------------------------------------------------

fig_beta_water <- make_beta_histogram(
  city_effect_sizes_usable,
  "Water"
)

fig_beta_shrub <- make_beta_histogram(
  city_effect_sizes_usable,
  "Shrub and scrub"
)

fig_beta_built <- make_beta_histogram(
  city_effect_sizes_usable,
  "Built"
)


# ------------------------------------------------------------
# Combine panels
# ------------------------------------------------------------

figure_02 <- ggpubr::ggarrange(
  fig_beta_water,
  fig_beta_shrub,
  fig_beta_built,
  ncol = 3,
  nrow = 1,
  labels = c(
    "A",
    "B",
    "C"
  ),
  font.label = list(
    size = 14,
    face = "bold"
  )
)


# ------------------------------------------------------------
# Save figure
# ------------------------------------------------------------

ggsave(
  file.path(
    OUTPUT_ROOT,
    "Figure_02_Beta_Distributions.png"
  ),
  figure_02,
  width = 14,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    OUTPUT_ROOT,
    "Figure_02_Beta_Distributions.pdf"
  ),
  figure_02,
  width = 14,
  height = 5,
  units = "in",
  bg = "white"
)




# ============================================================
# Stage 2
# City-specific habitat coefficient ~ city aridity
# ============================================================

fit_stage2_model <- function(
    habitat_name
) {
  
  model_data <- city_effect_sizes_usable %>%
    filter(
      habitat == habitat_name
    )
  
  fit <- lm(
    beta ~ avg_aridity_index,
    data = model_data
  )
  
  coefficient_table <-
    summary(fit)$coefficients
  
  confint_table <-
    confint(fit)
  
  model_summary <-
    summary(fit)
  
  stats <- data.frame(
    habitat =
      habitat_name,
    n_cities =
      nrow(model_data),
    intercept =
      coefficient_table[
        "(Intercept)",
        "Estimate"
      ],
    aridity_slope =
      coefficient_table[
        "avg_aridity_index",
        "Estimate"
      ],
    aridity_SE =
      coefficient_table[
        "avg_aridity_index",
        "Std. Error"
      ],
    aridity_t =
      coefficient_table[
        "avg_aridity_index",
        "t value"
      ],
    aridity_p =
      coefficient_table[
        "avg_aridity_index",
        "Pr(>|t|)"
      ],
    aridity_CI_lower =
      confint_table[
        "avg_aridity_index",
        1
      ],
    aridity_CI_upper =
      confint_table[
        "avg_aridity_index",
        2
      ],
    R2 =
      model_summary$r.squared,
    adjusted_R2 =
      model_summary$adj.r.squared,
    stringsAsFactors = FALSE
  )
  
  list(
    data = model_data,
    fit = fit,
    stats = stats
  )
}

stage2_water <- fit_stage2_model(
  "Water"
)

stage2_shrub <- fit_stage2_model(
  "Shrub and scrub"
)

stage2_built <- fit_stage2_model(
  "Built"
)

stage2_statistics <- bind_rows(
  stage2_water$stats,
  stage2_shrub$stats,
  stage2_built$stats
)

write.csv(
  stage2_statistics,
  file.path(
    OUTPUT_ROOT,
    "Table_03_Effect_Size_Aridity_Models.csv"
  ),
  row.names = FALSE
)



# ------------------------------------------------------------
# Stage-2 summary table
# Tests whether city-specific habitat effects vary with aridity
# ------------------------------------------------------------

stage2_results_table <- stage2_statistics %>%
  transmute(
    Relationship = case_when(
      habitat == "Water" ~
        "Water effect on richness ~ Aridity",
      
      habitat == "Shrub and scrub" ~
        "Shrub and scrub effect on richness ~ Aridity",
      
      habitat == "Built" ~
        "Built cover effect on richness ~ Aridity"
    ),
    
    n_cities =
      n_cities,
    
    Estimate =
      aridity_slope,
    
    SE =
      aridity_SE,
    
    t =
      aridity_t,
    
    CI_lower =
      aridity_CI_lower,
    
    CI_upper =
      aridity_CI_upper,
    
    p =
      aridity_p,
    
    R2 =
      R2,
    
    adjusted_R2 =
      adjusted_R2
  )

write.csv(
  stage2_results_table,
  file.path(
    OUTPUT_ROOT,
    "Table_04_Stage2_ES_Aridity_Results.csv"
  ),
  row.names = FALSE
)

print(
  stage2_results_table,
  row.names = FALSE
)

# ------------------------------------------------------------
# Stage-2 figures
# Central 95% of effect sizes shown for visualization only
# All cities remain in the fitted models
# ------------------------------------------------------------

make_stage2_panel <- function(
    model_data,
    habitat_name
) {
  
  effect_label <- switch(
    habitat_name,
    "Water" =
      "Effect of Water on Richness",
    "Shrub and scrub" =
      "Effect of Shrub and Scrub on Richness",
    "Built" =
      "Effect of Built Cover on Richness"
  )
  
  # Stage-2 model
  fit <- lm(
    beta ~ avg_aridity_index,
    data = model_data
  )
  
  # Predictions across the observed aridity gradient
  prediction_data <- data.frame(
    avg_aridity_index = seq(
      min(
        model_data$avg_aridity_index,
        na.rm = TRUE
      ),
      max(
        model_data$avg_aridity_index,
        na.rm = TRUE
      ),
      length.out = 200
    )
  )
  
  predictions <- predict(
    fit,
    newdata = prediction_data,
    interval = "confidence",
    level = 0.95
  )
  
  prediction_data <- cbind(
    prediction_data,
    as.data.frame(predictions)
  )
  
  # Central 95% used only for plotting limits
  effect_limits <- quantile(
    model_data$beta,
    probs = c(
      0.025,
      0.975
    ),
    na.rm = TRUE
  )
  
  ggplot() +
    
    geom_ribbon(
      data = prediction_data,
      aes(
        x = avg_aridity_index,
        ymin = lwr,
        ymax = upr
      ),
      alpha = 0.20
    ) +
    
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    
    geom_line(
      data = prediction_data,
      aes(
        x = avg_aridity_index,
        y = fit
      ),
      linewidth = 0.9
    ) +
    
    geom_point(
      data = model_data,
      aes(
        x = avg_aridity_index,
        y = beta
      ),
      alpha = 0.55,
      size = 2
    ) +
    
    coord_cartesian(
      ylim = effect_limits
    ) +
    
    labs(
      x = "Aridity Index (higher = less arid)",
      y = effect_label,
      title = habitat_name
    ) +
    
    theme_classic(
      base_size = 12
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )
}


# ------------------------------------------------------------
# Create individual panels
# ------------------------------------------------------------

fig_water <- make_stage2_panel(
  stage2_water$data,
  "Water"
)

fig_shrub <- make_stage2_panel(
  stage2_shrub$data,
  "Shrub and scrub"
)

fig_built <- make_stage2_panel(
  stage2_built$data,
  "Built"
)


# ------------------------------------------------------------
# Combine panels
# ------------------------------------------------------------

figure_01 <- ggpubr::ggarrange(
  fig_water,
  fig_shrub,
  fig_built,
  ncol = 3,
  nrow = 1,
  labels = c(
    "A",
    "B",
    "C"
  ),
  font.label = list(
    size = 14,
    face = "bold"
  )
)


# ------------------------------------------------------------
# Save figure
# ------------------------------------------------------------

ggsave(
  file.path(
    OUTPUT_ROOT,
    "Figure_01_City_Effect_Size_vs_Aridity.png"
  ),
  figure_01,
  width = 14,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    OUTPUT_ROOT,
    "Figure_01_City_Effect_Size_vs_Aridity.pdf"
  ),
  figure_01,
  width = 14,
  height = 5,
  units = "in",
  bg = "white"
)

# ------------------------------------------------------------
# QC figure: coefficient magnitude vs habitat range
# ------------------------------------------------------------

figure_qc3 <- city_effect_sizes_usable %>%
  mutate(
    abs_beta = abs(beta),
    habitat_range_pp = habitat_range * 10
  ) %>%
  ggplot(
    aes(
      x = habitat_range_pp,
      y = abs_beta
    )
  ) +
  geom_point(
    alpha = 0.4
  ) +
  facet_wrap(
    ~ habitat,
    scales = "free"
  ) +
  labs(
    x = "Within-city habitat range (percentage points)",
    y = "Absolute city-specific coefficient"
  ) +
  theme_classic(
    base_size = 12
  )

ggsave(
  file.path(
    OUTPUT_ROOT,
    "Figure_QC3_Coefficient_Magnitude_vs_Habitat_Range.png"
  ),
  figure_qc3,
  width = 10,
  height = 4,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ------------------------------------------------------------
# Save analysis object
# ------------------------------------------------------------

analysis_bundle <- list(
  settings = list(
    min_checklists_per_locality =
      MIN_CHECKLISTS_PER_LOCALITY,
    min_localities_per_city =
      MIN_LOCALITIES_PER_CITY,
    habitat_divisor =
      habitat_divisor
  ),
  missing_summary =
    missing_summary,
  locality_support =
    locality_support,
  city_qc =
    city_qc,
  city_habitat_qc =
    city_habitat_qc,
  correlation_pair_summary =
    correlation_pair_summary,
  threshold_summary =
    threshold_summary,
  stage1_data =
    sr_stage1,
  city_effect_sizes_all =
    city_effect_sizes_all,
  city_effect_sizes_usable =
    city_effect_sizes_usable,
  city_effect_sizes_wide =
    city_effect_sizes_wide,
  stage1_summary =
    stage1_summary,
  stage2_statistics =
    stage2_statistics,
  stage2_water =
    stage2_water,
  stage2_shrub =
    stage2_shrub,
  stage2_built =
    stage2_built
)

saveRDS(
  analysis_bundle,
  file.path(
    OUTPUT_ROOT,
    "Species_Richness_City_Effect_Sizes.rds"
  )
)


# ------------------------------------------------------------
# Console summary
# ------------------------------------------------------------

cat("\n")
cat("Species richness city-effect analysis complete\n")
cat("----------------------------------------------\n")

cat(
  "Cities after complete cases: ",
  n_distinct(sr_pool$NAME10),
  "\n",
  sep = ""
)

cat(
  "Cities passing 10 x 5 eligibility: ",
  nrow(eligible_cities),
  "\n",
  sep = ""
)

cat(
  "Stage-1 checklists retained: ",
  nrow(sr_stage1),
  "\n",
  sep = ""
)

cat("\nWithin-city habitat collinearity\n")
print(
  correlation_pair_summary,
  row.names = FALSE
)

cat("\nMaximum-correlation thresholds\n")
print(
  threshold_summary,
  row.names = FALSE
)

cat("\nStage-1 separate habitat models\n")
print(
  stage1_summary,
  row.names = FALSE
)

cat("\nStage-2 effect-size ~ aridity models\n")
print(
  stage2_statistics,
  row.names = FALSE
)













# ============================================================
# FINAL ARTIFACTS
# ============================================================

FINAL_OUTPUT_DIR <- file.path(
  OUTPUT_ROOT,
  "Final Artifacts"
)

dir.create(
  FINAL_OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Final Stage-1 effect-size distribution figure
# ------------------------------------------------------------

ggsave(
  file.path(
    FINAL_OUTPUT_DIR,
    "Figure_01_ES_Distributions.png"
  ),
  figure_02,
  width = 14,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    FINAL_OUTPUT_DIR,
    "Figure_01_ES_Distributions.pdf"
  ),
  figure_02,
  width = 14,
  height = 5,
  units = "in",
  bg = "white"
)


# ------------------------------------------------------------
# 2. Final Aridity x effect-size figure
# ------------------------------------------------------------

ggsave(
  file.path(
    FINAL_OUTPUT_DIR,
    "Figure_02_Aridity_x_ES.png"
  ),
  figure_01,
  width = 14,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    FINAL_OUTPUT_DIR,
    "Figure_02_Aridity_x_ES.pdf"
  ),
  figure_01,
  width = 14,
  height = 5,
  units = "in",
  bg = "white"
)


# ------------------------------------------------------------
# 3. Final Stage-2 statistical results table
# ------------------------------------------------------------

final_stats_table <- stage2_statistics %>%
  transmute(
    Relationship = case_when(
      habitat == "Water" ~
        "Water effect on richness ~ Aridity",
      
      habitat == "Shrub and scrub" ~
        "Shrub and scrub effect on richness ~ Aridity",
      
      habitat == "Built" ~
        "Built cover effect on richness ~ Aridity"
    ),
    
    n_cities =
      n_cities,
    
    Estimate =
      aridity_slope,
    
    SE =
      aridity_SE,
    
    t =
      aridity_t,
    
    df =
      n_cities - 2,
    
    CI_lower =
      aridity_CI_lower,
    
    CI_upper =
      aridity_CI_upper,
    
    p =
      aridity_p,
    
    R2 =
      R2,
    
    adjusted_R2 =
      adjusted_R2
  )

write.csv(
  final_stats_table,
  file.path(
    FINAL_OUTPUT_DIR,
    "Table_01_Aridity_x_ES_Statistics.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# Final console output
# ------------------------------------------------------------

cat("\n")
cat("FINAL ARTIFACTS\n")
cat("===============\n")
cat("1. Figure_01_ES_Distributions.png / .pdf\n")
cat("2. Figure_02_Aridity_x_ES.png / .pdf\n")
cat("3. Table_01_Aridity_x_ES_Statistics.csv\n\n")

print(
  final_stats_table,
  row.names = FALSE
)