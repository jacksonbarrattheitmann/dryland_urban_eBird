# ============================================================
# Dryland Bird Project
# 02_Species_Richness.R
#
# Purpose:
# Model checklist-level bird species richness as a function of
# aridity, local habitat, and checklist effort using a GAM/BAM
# framework and a balanced hierarchical subsample.
#
# Sampling design:
#   1. Retain Stationary and Traveling checklists
#   2. Stationary checklist distance = 0 km
#   3. Duration <= 300 minutes
#   4. Distance <= 5 km
#   5. Require >= 5 valid checklists per locality
#   6. Require >= 5 eligible localities per city
#   7. Randomly sample 5 localities per city
#   8. Randomly sample 5 checklists per selected locality
#
# Final model:
#   Negative-binomial BAM with log link
#
# Final response structure:
#   City level:     avg_aridity_index
#   Locality level: water, shrub_and_scrub, built
#   Checklist level: species richness, duration, distance,
#                    protocol
#
# Random-effect smooths:
#   s(NAME10, bs = "re")
#   s(locality_id, bs = "re")
#
# Outputs:
#   Outputs/02_Species_Richness/
#     Figure/
#       Figure_02_Main_Data_Trends.png
#       Figure_02_Main_Data_Trends.pdf
#     Table/
#       Table_02_Main_Model_Results.csv
#     Data/
#       Species_Richness_Manuscript_Input.rds
#
# Figure design:
#   A. Aridity vs mean species richness (city level)
#   B. Water vs mean species richness (locality level)
#   C. Shrub/scrub vs mean species richness (locality level)
#   D. Built vs mean species richness (locality level)
#
# Panels B-D color locality points by the percentile rank of the
# continuous city-level aridity index:
#   warmer colors = more arid
#   cooler colors = less arid
#
# Percentile rank is used only for visualization so that the full
# warm-to-cool color range is visible across the observed data.
# The statistical models continue to use the original continuous
# aridity index.
#
# The descriptive GAM smooths are unadjusted and are not the
# source of inferential tests. Statistical inference comes from
# the final negative-binomial BAM and is exported in Table 02.
#
# Canonical input:
# Intermediate_data/dfs_for_analysis_combined/combined_data.RDS
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(mgcv)
library(ggpubr)


# ------------------------------------------------------------
# 2. Analysis settings
# ------------------------------------------------------------

MIN_CHECKLISTS_PER_LOCALITY <- 5
MIN_LOCALITIES_PER_CITY <- 5
N_LOCALITIES_PER_CITY <- 5
N_CHECKLISTS_PER_LOCALITY <- 5

LOCALITY_SAMPLE_SEED <- 20260911
CHECKLIST_SAMPLE_SEED <- 20260912
K_CHECK_SEED <- 20260912

TEMPORARY_DEDUPLICATE <- TRUE


# ------------------------------------------------------------
# 3. Output directories
# ------------------------------------------------------------

OUTPUT_ROOT <- file.path(
  "Outputs",
  "02_Species_Richness"
)

FIGURE_DIR <- file.path(
  OUTPUT_ROOT,
  "Figure"
)

TABLE_DIR <- file.path(
  OUTPUT_ROOT,
  "Table"
)

DATA_DIR <- file.path(
  OUTPUT_ROOT,
  "Data"
)

for (x in c(FIGURE_DIR, TABLE_DIR, DATA_DIR)) {
  dir.create(
    x,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ------------------------------------------------------------
# 4. Load combined dataset
# ------------------------------------------------------------

data_path <- file.path(
  "Intermediate_data",
  "dfs_for_analysis_combined",
  "combined_data.RDS"
)

if (!file.exists(data_path)) {
  stop(
    "combined_data.RDS was not found at:\n",
    data_path
  )
}

combined_data <- readRDS(data_path)

cat("\n============================================================\n")
cat("SPECIES RICHNESS GAM ANALYSIS\n")
cat("============================================================\n\n")

cat(
  "Rows in combined_data: ",
  nrow(combined_data),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 5. Required-variable check
# ------------------------------------------------------------

required_variables <- c(
  "checklist_id",
  "locality_id",
  "NAME10",
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
    "Required variables are missing:\n",
    paste(missing_variables, collapse = "\n")
  )
}

cat("PASS: All required analysis variables are present.\n")


# ------------------------------------------------------------
# 6. Duplicate checklist handling
# ------------------------------------------------------------

duplicate_checklists <- combined_data %>%
  count(
    checklist_id,
    name = "n_rows"
  ) %>%
  filter(
    n_rows > 1
  )

cat(
  "\nDuplicated checklist IDs in input: ",
  nrow(duplicate_checklists),
  "\n",
  sep = ""
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

cat(
  "Rows after duplicate handling: ",
  nrow(analysis_base),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. Apply protocol and effort filters
# ------------------------------------------------------------

sr_pool <- analysis_base %>%
  mutate(
    effort_distance_km =
      if_else(
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

cat("\n============================================================\n")
cat("EFFORT-FILTERED DATA\n")
cat("============================================================\n\n")

sr_pool %>%
  summarise(
    n_checklists = n(),
    n_localities = n_distinct(locality_id),
    n_cities = n_distinct(NAME10)
  ) %>%
  print()


# ------------------------------------------------------------
# 8. Complete-case restriction
# ------------------------------------------------------------

model_variables <- c(
  "sr",
  "avg_aridity_index",
  "water",
  "shrub_and_scrub",
  "built",
  "duration_minutes",
  "effort_distance_km",
  "protocol_name",
  "NAME10",
  "locality_id"
)

missing_summary <- data.frame(
  variable = model_variables,
  missing_n = sapply(
    sr_pool[model_variables],
    function(x) sum(is.na(x))
  )
)

missing_summary$missing_percent <-
  100 * missing_summary$missing_n / nrow(sr_pool)

cat("\n============================================================\n")
cat("MODEL-VARIABLE MISSINGNESS\n")
cat("============================================================\n\n")

print(missing_summary)

sr_pool <- sr_pool %>%
  filter(
    if_all(
      all_of(model_variables),
      ~ !is.na(.)
    )
  )


# ------------------------------------------------------------
# 9. Identify eligible localities
# ------------------------------------------------------------

locality_sampling <- sr_pool %>%
  count(
    NAME10,
    locality_id,
    name = "n_checklists"
  )

eligible_localities <- locality_sampling %>%
  filter(
    n_checklists >= MIN_CHECKLISTS_PER_LOCALITY
  )


# ------------------------------------------------------------
# 10. Identify eligible cities
# ------------------------------------------------------------

city_sampling <- eligible_localities %>%
  count(
    NAME10,
    name = "n_eligible_localities"
  )

eligible_cities <- city_sampling %>%
  filter(
    n_eligible_localities >= MIN_LOCALITIES_PER_CITY
  )

cat("\n============================================================\n")
cat("ELIGIBILITY\n")
cat("============================================================\n\n")

cat(
  "Eligible localities: ",
  nrow(eligible_localities),
  "\n",
  sep = ""
)

cat(
  "Eligible cities:     ",
  nrow(eligible_cities),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 11. Create eligible sampling pool
# ------------------------------------------------------------

eligible_pool <- sr_pool %>%
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


# ------------------------------------------------------------
# 12. Randomly select five localities per city
# ------------------------------------------------------------

set.seed(LOCALITY_SAMPLE_SEED)

selected_localities <- eligible_localities %>%
  semi_join(
    eligible_cities,
    by = "NAME10"
  ) %>%
  arrange(
    NAME10,
    locality_id
  ) %>%
  group_by(NAME10) %>%
  slice_sample(
    n = N_LOCALITIES_PER_CITY
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 13. Randomly select five checklists per locality
# ------------------------------------------------------------

set.seed(CHECKLIST_SAMPLE_SEED)

sr_data <- eligible_pool %>%
  semi_join(
    selected_localities,
    by = c(
      "NAME10",
      "locality_id"
    )
  ) %>%
  arrange(
    NAME10,
    locality_id,
    checklist_id
  ) %>%
  group_by(
    NAME10,
    locality_id
  ) %>%
  slice_sample(
    n = N_CHECKLISTS_PER_LOCALITY
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 14. Balanced sample QC
# ------------------------------------------------------------

sample_locality_qc <- sr_data %>%
  count(
    NAME10,
    locality_id,
    name = "n_checklists"
  )

sample_city_qc <- sr_data %>%
  group_by(NAME10) %>%
  summarise(
    n_checklists = n(),
    n_localities = n_distinct(locality_id),
    .groups = "drop"
  )

bad_localities <- sum(
  sample_locality_qc$n_checklists !=
    N_CHECKLISTS_PER_LOCALITY
)

bad_cities_localities <- sum(
  sample_city_qc$n_localities !=
    N_LOCALITIES_PER_CITY
)

bad_cities_checklists <- sum(
  sample_city_qc$n_checklists !=
    (
      N_LOCALITIES_PER_CITY *
        N_CHECKLISTS_PER_LOCALITY
    )
)

cat("\n============================================================\n")
cat("BALANCED SAMPLE QC\n")
cat("============================================================\n\n")

cat(
  "Sampled cities:      ",
  n_distinct(sr_data$NAME10),
  "\n",
  sep = ""
)

cat(
  "Sampled localities:  ",
  n_distinct(sr_data$locality_id),
  "\n",
  sep = ""
)

cat(
  "Sampled checklists:  ",
  nrow(sr_data),
  "\n",
  sep = ""
)

cat(
  "Expected checklists: ",
  N_LOCALITIES_PER_CITY *
    N_CHECKLISTS_PER_LOCALITY *
    nrow(eligible_cities),
  "\n",
  sep = ""
)

cat(
  "\nLocalities not containing exactly 5 checklists: ",
  bad_localities,
  "\n",
  sep = ""
)

cat(
  "Cities not containing exactly 5 localities:     ",
  bad_cities_localities,
  "\n",
  sep = ""
)

cat(
  "Cities not containing exactly 25 checklists:    ",
  bad_cities_checklists,
  "\n",
  sep = ""
)

if (
  bad_localities > 0 ||
  bad_cities_localities > 0 ||
  bad_cities_checklists > 0
) {
  stop("Balanced sampling QC failed.")
}

cat("\nPASS: Balanced sampling design verified.\n")


# ------------------------------------------------------------
# 15. Prepare grouping variables
# ------------------------------------------------------------

sr_data <- sr_data %>%
  mutate(
    NAME10 = factor(NAME10),
    locality_id = factor(locality_id),
    protocol_name = factor(
      protocol_name,
      levels = c(
        "Stationary",
        "Traveling"
      )
    )
  )


# ------------------------------------------------------------
# 16. Examine richness distribution
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SPECIES RICHNESS DISTRIBUTION\n")
cat("============================================================\n\n")

print(summary(sr_data$sr))

sr_mean <- mean(sr_data$sr)
sr_variance <- var(sr_data$sr)
sr_var_mean <- sr_variance / sr_mean

cat(
  "\nMean richness: ",
  round(sr_mean, 3),
  "\n",
  sep = ""
)

cat(
  "Variance:      ",
  round(sr_variance, 3),
  "\n",
  sep = ""
)

cat(
  "Variance/mean: ",
  round(sr_var_mean, 3),
  "\n",
  sep = ""
)


# ============================================================
# GAM / BAM MODELS
# ============================================================


# ------------------------------------------------------------
# 17. Global model formula
# ------------------------------------------------------------

sr_formula <- sr ~
  s(
    avg_aridity_index,
    k = 5
  ) +
  s(
    water,
    k = 5
  ) +
  s(
    shrub_and_scrub,
    k = 5
  ) +
  s(
    built,
    k = 5
  ) +
  ti(
    avg_aridity_index,
    water,
    k = c(5, 5)
  ) +
  ti(
    avg_aridity_index,
    shrub_and_scrub,
    k = c(5, 5)
  ) +
  ti(
    avg_aridity_index,
    built,
    k = c(5, 5)
  ) +
  s(
    duration_minutes,
    k = 10
  ) +
  s(
    effort_distance_km,
    k = 5
  ) +
  protocol_name +
  s(
    NAME10,
    bs = "re"
  ) +
  s(
    locality_id,
    bs = "re"
  )


# ------------------------------------------------------------
# 18. Fit Poisson BAM
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FITTING POISSON BAM\n")
cat("============================================================\n\n")

sr_model_poisson <- bam(
  formula = sr_formula,
  family = poisson(link = "log"),
  data = sr_data,
  method = "fREML",
  discrete = TRUE
)


# ------------------------------------------------------------
# 19. Poisson model summary and dispersion
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("POISSON BAM SUMMARY\n")
cat("============================================================\n\n")

print(summary(sr_model_poisson))

poisson_dispersion <- sum(
  residuals(
    sr_model_poisson,
    type = "pearson"
  ) ^ 2
) /
  df.residual(sr_model_poisson)

cat(
  "\nPoisson Pearson dispersion ratio: ",
  round(poisson_dispersion, 3),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 20. Fit negative-binomial BAM
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FITTING NEGATIVE-BINOMIAL BAM\n")
cat("============================================================\n\n")

sr_model_nb <- bam(
  formula = sr_formula,
  family = nb(link = "log"),
  data = sr_data,
  method = "fREML",
  discrete = TRUE
)


# ------------------------------------------------------------
# 21. Negative-binomial summary and dispersion
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("NEGATIVE-BINOMIAL BAM SUMMARY\n")
cat("============================================================\n\n")

print(summary(sr_model_nb))

nb_dispersion <- sum(
  residuals(
    sr_model_nb,
    type = "pearson"
  ) ^ 2
) /
  df.residual(sr_model_nb)

cat(
  "\nNegative-binomial Pearson dispersion ratio: ",
  round(nb_dispersion, 3),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 22. Compare count families
# ------------------------------------------------------------

model_comparison <- data.frame(
  model = c(
    "Poisson BAM",
    "Negative-binomial BAM"
  ),
  AIC = c(
    AIC(sr_model_poisson),
    AIC(sr_model_nb)
  ),
  Pearson_dispersion = c(
    poisson_dispersion,
    nb_dispersion
  )
) %>%
  mutate(
    delta_AIC =
      AIC - min(AIC)
  )

cat("\n============================================================\n")
cat("COUNT-FAMILY COMPARISON\n")
cat("============================================================\n\n")

print(model_comparison)


# ------------------------------------------------------------
# 23. Select final richness model
# ------------------------------------------------------------

sr_model_final <- sr_model_nb
sr_summary <- summary(sr_model_final)

cat("\n============================================================\n")
cat("FINAL SPECIES RICHNESS BAM\n")
cat("============================================================\n\n")

print(sr_summary)


# ============================================================
# FINAL MODEL QC
# ============================================================


# ------------------------------------------------------------
# 24. Convergence and rank
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL MODEL CONVERGENCE AND RANK\n")
cat("============================================================\n\n")

cat(
  "Model converged: ",
  sr_model_final$converged,
  "\n",
  sep = ""
)

cat(
  "Model rank: ",
  sr_model_final$rank,
  " / ",
  length(coef(sr_model_final)),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 25. Stable basis-dimension check
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STABLE K CHECK - 1000 REPLICATES\n")
cat("============================================================\n\n")

set.seed(K_CHECK_SEED)

sr_k_stable <- k.check(
  sr_model_final,
  n.rep = 1000
)

print(sr_k_stable)

k_valid_rows <- !is.na(
  sr_k_stable[, "k-index"]
)

k_followup_flag <- any(
  sr_k_stable[
    k_valid_rows,
    "k-index"
  ] < 1 &
    sr_k_stable[
      k_valid_rows,
      "p-value"
    ] < 0.05 &
    (
      sr_k_stable[
        k_valid_rows,
        "edf"
      ] /
        sr_k_stable[
          k_valid_rows,
          "k'"
        ]
    ) >= 0.80
)

cat(
  "\nBasis-dimension follow-up required: ",
  k_followup_flag,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 26. Standard GAM diagnostics
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GAM CHECK\n")
cat("============================================================\n\n")

set.seed(K_CHECK_SEED)
gam.check(sr_model_final)


# ------------------------------------------------------------
# 27. Concurvity
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FULL CONCURVITY CHECK\n")
cat("============================================================\n\n")

sr_concurvity_full <- concurvity(
  sr_model_final,
  full = TRUE
)

print(sr_concurvity_full)

cat("\n============================================================\n")
cat("PAIRWISE CONCURVITY CHECK\n")
cat("============================================================\n\n")

sr_concurvity_pairwise <- concurvity(
  sr_model_final,
  full = FALSE
)

print(
  round(
    sr_concurvity_pairwise$estimate,
    3
  )
)


# ------------------------------------------------------------
# 28. Variance components
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("VARIANCE COMPONENTS\n")
cat("============================================================\n\n")

sr_vcomp <- gam.vcomp(sr_model_final)
print(sr_vcomp)


# ------------------------------------------------------------
# 29. Residual diagnostics
# ------------------------------------------------------------

sr_pearson_residuals <- residuals(
  sr_model_final,
  type = "pearson"
)

sr_deviance_residuals <- residuals(
  sr_model_final,
  type = "deviance"
)

sr_fitted_values <- fitted(sr_model_final)

cat("\n============================================================\n")
cat("RESIDUAL SUMMARY\n")
cat("============================================================\n\n")

cat("Pearson residuals:\n")
print(summary(sr_pearson_residuals))

cat("\nDeviance residuals:\n")
print(summary(sr_deviance_residuals))

cat(
  "\nAbsolute Pearson residual > 2: ",
  sum(abs(sr_pearson_residuals) > 2),
  " (",
  round(
    100 * mean(abs(sr_pearson_residuals) > 2),
    2
  ),
  "%)\n",
  sep = ""
)

cat(
  "Absolute Pearson residual > 3: ",
  sum(abs(sr_pearson_residuals) > 3),
  " (",
  round(
    100 * mean(abs(sr_pearson_residuals) > 3),
    2
  ),
  "%)\n",
  sep = ""
)

plot(
  sr_fitted_values,
  sr_pearson_residuals,
  pch = 16,
  cex = 0.5,
  xlab = "Fitted species richness",
  ylab = "Pearson residual"
)

abline(
  h = 0,
  lty = 2
)


# ------------------------------------------------------------
# 30. Extract model results
# ------------------------------------------------------------

sr_parametric_results <- as.data.frame(
  sr_summary$p.table
)

sr_parametric_results$term <- rownames(
  sr_parametric_results
)

sr_parametric_results <- sr_parametric_results %>%
  select(
    term,
    everything()
  )

sr_smooth_results <- as.data.frame(
  sr_summary$s.table
)

sr_smooth_results$term <- rownames(
  sr_smooth_results
)

sr_smooth_results <- sr_smooth_results %>%
  select(
    term,
    everything()
  )

cat("\n============================================================\n")
cat("PARAMETRIC TERMS\n")
cat("============================================================\n\n")

print(
  sr_parametric_results,
  row.names = FALSE
)

cat("\n============================================================\n")
cat("SMOOTH TERMS\n")
cat("============================================================\n\n")

print(
  sr_smooth_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 31. Model fit summary
# ------------------------------------------------------------

sr_theta <- sr_model_final$family$getTheta(TRUE)

cat("\n============================================================\n")
cat("MODEL FIT\n")
cat("============================================================\n\n")

cat(
  "Adjusted R-squared: ",
  round(sr_summary$r.sq, 3),
  "\n",
  sep = ""
)

cat(
  "Deviance explained: ",
  round(
    100 * sr_summary$dev.expl,
    2
  ),
  "%\n",
  sep = ""
)

cat(
  "Negative-binomial theta: ",
  round(sr_theta, 3),
  "\n",
  sep = ""
)


# ============================================================
# MAIN PUBLICATION / MEETING FIGURE
# ============================================================
#
# Descriptive empirical trends are shown for the four focal
# ecological predictors.
#
# Panel A:
#   one point per city.
#
# Panels B-D:
#   one point per locality, colored by city-level aridity.
#
# Warm colors represent more-arid cities (lower index values).
# Cool colors represent less-arid cities (higher index values).
#
# Smooths are descriptive and unadjusted.
# Statistical inference is based on the final negative-binomial
# BAM and is reported in Table 02.
# ============================================================


# ------------------------------------------------------------
# 32. Create empirical trend datasets
# ------------------------------------------------------------

sr_city_plot <- sr_data %>%
  group_by(NAME10) %>%
  summarise(
    avg_aridity_index = first(avg_aridity_index),
    mean_sr = mean(
      sr,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

sr_locality_plot <- sr_data %>%
  group_by(
    NAME10,
    locality_id
  ) %>%
  summarise(
    avg_aridity_index = first(avg_aridity_index),
    water = first(water),
    shrub_and_scrub = first(shrub_and_scrub),
    built = first(built),
    mean_sr = mean(
      sr,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    aridity_rank = percent_rank(
      avg_aridity_index
    )
  )


# ------------------------------------------------------------
# 33. Convert habitat cover to percent for plotting
# ------------------------------------------------------------

maximum_habitat_value <- max(
  c(
    sr_locality_plot$water,
    sr_locality_plot$shrub_and_scrub,
    sr_locality_plot$built
  ),
  na.rm = TRUE
)

if (maximum_habitat_value <= 1.05) {
  sr_locality_plot <- sr_locality_plot %>%
    mutate(
      water_plot = water * 100,
      shrub_plot = shrub_and_scrub * 100,
      built_plot = built * 100
    )
} else {
  sr_locality_plot <- sr_locality_plot %>%
    mutate(
      water_plot = water,
      shrub_plot = shrub_and_scrub,
      built_plot = built
    )
}


# ------------------------------------------------------------
# 34. Common figure settings
# ------------------------------------------------------------

main_y_limits <- range(
  c(
    sr_city_plot$mean_sr,
    sr_locality_plot$mean_sr
  ),
  na.rm = TRUE
)

aridity_palette <- c(
  "#7F0000",
  "#B30000",
  "#D7301F",
  "#F46D43",
  "#FDAE61",
  "#FEE08B",
  "#D9D9D9",
  "#E0F3F8",
  "#ABD9E9",
  "#74ADD1",
  "#4575B4",
  "#313695"
)

theme_sr_main <- function() {
  theme_classic(
    base_size = 12
  ) +
    theme(
      axis.text = element_text(
        colour = "black"
      ),
      plot.title = element_text(
        face = "bold",
        size = 12,
        hjust = 0.5
      )
    )
}

aridity_color_scale <- scale_color_gradientn(
  colors = aridity_palette,
  values = c(
    0.00,
    0.07,
    0.15,
    0.25,
    0.36,
    0.46,
    0.50,
    0.56,
    0.66,
    0.76,
    0.88,
    1.00
  ),
  limits = c(
    0,
    1
  ),
  breaks = c(
    0,
    0.5,
    1
  ),
  labels = c(
    "More arid",
    "Intermediate",
    "Less arid"
  ),
  name = "Aridity gradient"
)


# ------------------------------------------------------------
# 35. Panel A: Aridity
# ------------------------------------------------------------

fig_sr_aridity <- ggplot(
  sr_city_plot,
  aes(
    x = avg_aridity_index,
    y = mean_sr
  )
) +
  geom_point(
    alpha = 0.45,
    size = 1.8
  ) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 4),
    se = TRUE,
    color = "black",
    fill = "grey80"
  ) +
  coord_cartesian(
    ylim = main_y_limits
  ) +
  labs(
    x = "Aridity index (higher = less arid)",
    y = "Mean species richness",
    title = "Aridity"
  ) +
  theme_sr_main() +
  theme(
    legend.position = "none"
  )


# ------------------------------------------------------------
# 36. Panel B: Water cover
# ------------------------------------------------------------

fig_sr_water <- ggplot(
  sr_locality_plot,
  aes(
    x = water_plot,
    y = mean_sr
  )
) +
  geom_point(
    aes(
      color = aridity_rank
    ),
    alpha = 0.55,
    size = 1.5
  ) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 4),
    se = TRUE,
    color = "black",
    fill = "grey80"
  ) +
  aridity_color_scale +
  coord_cartesian(
    ylim = main_y_limits
  ) +
  labs(
    x = "Water cover (%)",
    y = "Mean species richness",
    title = "Water"
  ) +
  theme_sr_main()


# ------------------------------------------------------------
# 37. Panel C: Shrub and scrub cover
# ------------------------------------------------------------

fig_sr_shrub <- ggplot(
  sr_locality_plot,
  aes(
    x = shrub_plot,
    y = mean_sr
  )
) +
  geom_point(
    aes(
      color = aridity_rank
    ),
    alpha = 0.55,
    size = 1.5
  ) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 4),
    se = TRUE,
    color = "black",
    fill = "grey80"
  ) +
  aridity_color_scale +
  coord_cartesian(
    ylim = main_y_limits
  ) +
  labs(
    x = "Shrub and scrub cover (%)",
    y = "Mean species richness",
    title = "Shrub and scrub"
  ) +
  theme_sr_main()


# ------------------------------------------------------------
# 38. Panel D: Built cover
# ------------------------------------------------------------

fig_sr_built <- ggplot(
  sr_locality_plot,
  aes(
    x = built_plot,
    y = mean_sr
  )
) +
  geom_point(
    aes(
      color = aridity_rank
    ),
    alpha = 0.55,
    size = 1.5
  ) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x, k = 4),
    se = TRUE,
    color = "black",
    fill = "grey80"
  ) +
  aridity_color_scale +
  coord_cartesian(
    ylim = main_y_limits
  ) +
  labs(
    x = "Built cover (%)",
    y = "Mean species richness",
    title = "Built"
  ) +
  theme_sr_main()


# ------------------------------------------------------------
# 39. Combine empirical panels
# ------------------------------------------------------------

figure_02_main <- ggpubr::ggarrange(
  fig_sr_aridity,
  fig_sr_water,
  fig_sr_shrub,
  fig_sr_built,
  ncol = 2,
  nrow = 2,
  labels = c(
    "A",
    "B",
    "C",
    "D"
  ),
  font.label = list(
    size = 14,
    face = "bold"
  ),
  common.legend = TRUE,
  legend = "bottom"
)

print(figure_02_main)


# ------------------------------------------------------------
# 40. Save main figure
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    FIGURE_DIR,
    "Figure_02_Main_Data_Trends.png"
  ),
  plot = figure_02_main,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    FIGURE_DIR,
    "Figure_02_Main_Data_Trends.pdf"
  ),
  plot = figure_02_main,
  width = 10,
  height = 8,
  units = "in",
  bg = "white"
)


# ============================================================
# MAIN STATISTICAL TABLE
# ============================================================


# ------------------------------------------------------------
# 41. Focal ecological term labels
# ------------------------------------------------------------

focal_term_labels <- c(
  "s(avg_aridity_index)" =
    "Aridity",
  "s(water)" =
    "Water cover",
  "s(shrub_and_scrub)" =
    "Shrub and scrub cover",
  "s(built)" =
    "Built cover",
  "ti(avg_aridity_index,water)" =
    "Aridity x Water",
  "ti(avg_aridity_index,shrub_and_scrub)" =
    "Aridity x Shrub and scrub",
  "ti(avg_aridity_index,built)" =
    "Aridity x Built"
)


# ------------------------------------------------------------
# 42. P-value formatting helper
# ------------------------------------------------------------

format_p_value <- function(p) {
  ifelse(
    p < 0.001,
    "<0.001",
    sprintf("%.3f", p)
  )
}


# ------------------------------------------------------------
# 43. Build main model-results table
# ------------------------------------------------------------

table_02_main <- sr_smooth_results %>%
  filter(
    term %in% names(focal_term_labels)
  ) %>%
  mutate(
    term_order = match(
      term,
      names(focal_term_labels)
    ),
    Predictor = unname(
      focal_term_labels[term]
    ),
    EDF = sprintf(
      "%.2f",
      edf
    ),
    Reference_df = sprintf(
      "%.2f",
      Ref.df
    ),
    F_statistic = sprintf(
      "%.2f",
      F
    ),
    P_value = format_p_value(
      `p-value`
    )
  ) %>%
  arrange(term_order) %>%
  select(
    Predictor,
    EDF,
    Reference_df,
    F_statistic,
    P_value
  )

names(table_02_main) <- c(
  "Predictor",
  "EDF",
  "Ref. df",
  "F statistic",
  "P-value"
)

if (
  nrow(table_02_main) !=
  length(focal_term_labels)
) {
  stop(
    "Table 02 does not contain all seven focal ecological terms."
  )
}

cat("\n============================================================\n")
cat("TABLE 02 - MAIN MODEL RESULTS\n")
cat("============================================================\n\n")

print(
  table_02_main,
  row.names = FALSE
)


# ------------------------------------------------------------
# 44. Save main statistical table
# ------------------------------------------------------------

write.csv(
  table_02_main,
  file = file.path(
    TABLE_DIR,
    "Table_02_Main_Model_Results.csv"
  ),
  row.names = FALSE
)


# ============================================================
# MANUSCRIPT-SYNTHESIS EXPORT
# ============================================================


# ------------------------------------------------------------
# 45. Build standardized manuscript-input bundle
# ------------------------------------------------------------

sr_model_metadata <- data.frame(
  response_id =
    "Species_Richness",
  response_label =
    "Species richness",
  response_variable =
    "sr",
  model_family =
    "Negative binomial",
  link =
    "log",
  n_cities =
    n_distinct(sr_data$NAME10),
  n_localities =
    n_distinct(sr_data$locality_id),
  n_checklists =
    nrow(sr_data),
  adjusted_R2 =
    sr_summary$r.sq,
  deviance_explained =
    sr_summary$dev.expl,
  family_parameter =
    sr_theta,
  family_parameter_name =
    "theta",
  poisson_AIC =
    AIC(sr_model_poisson),
  final_AIC =
    AIC(sr_model_final),
  poisson_dispersion =
    poisson_dispersion,
  final_dispersion =
    nb_dispersion,
  stringsAsFactors = FALSE
)

species_richness_manuscript_input <- list(
  response_id =
    "Species_Richness",
  response_label =
    "Species richness",
  response_variable =
    "sr",
  analysis_data =
    sr_data,
  final_model =
    sr_model_final,
  smooth_results =
    sr_smooth_results,
  parametric_results =
    sr_parametric_results,
  focal_results =
    table_02_main,
  model_metadata =
    sr_model_metadata,
  plot_data = list(
    city = sr_city_plot,
    locality = sr_locality_plot
  )
)

species_richness_manuscript_file <- file.path(
  DATA_DIR,
  "Species_Richness_Manuscript_Input.rds"
)

saveRDS(
  species_richness_manuscript_input,
  species_richness_manuscript_file
)

if (
  !file.exists(
    species_richness_manuscript_file
  )
) {
  stop(
    "Species-richness manuscript-input RDS was not created."
  )
}

cat("\n============================================================\n")
cat("MANUSCRIPT-SYNTHESIS EXPORT\n")
cat("============================================================\n\n")

cat(
  "Saved manuscript input:\n",
  species_richness_manuscript_file,
  "\n",
  sep = ""
)


# ============================================================
# AUTOMATED QC AND OUTPUT CHECK
# ============================================================


# ------------------------------------------------------------
# 46. Automated analysis QC
# ------------------------------------------------------------

qc_balanced_sample <-
  bad_localities == 0 &&
  bad_cities_localities == 0 &&
  bad_cities_checklists == 0

qc_converged <- isTRUE(
  sr_model_final$converged
)

qc_full_rank <-
  sr_model_final$rank ==
  length(coef(sr_model_final))

qc_basis_dimension <-
  !k_followup_flag

sr_qc_summary <- data.frame(
  check = c(
    "Balanced 5 x 5 sampling",
    "Final BAM converged",
    "Final BAM is full rank",
    "Basis dimensions acceptable"
  ),
  pass = c(
    qc_balanced_sample,
    qc_converged,
    qc_full_rank,
    qc_basis_dimension
  )
)

cat("\n============================================================\n")
cat("SPECIES RICHNESS AUTOMATED QC SUMMARY\n")
cat("============================================================\n\n")

print(
  sr_qc_summary,
  row.names = FALSE
)

cat(
  "\nChecks passed: ",
  sum(sr_qc_summary$pass),
  " of ",
  nrow(sr_qc_summary),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 47. Output check
# ------------------------------------------------------------

expected_output_files <- c(
  file.path(
    FIGURE_DIR,
    "Figure_02_Main_Data_Trends.png"
  ),
  file.path(
    FIGURE_DIR,
    "Figure_02_Main_Data_Trends.pdf"
  ),
  file.path(
    TABLE_DIR,
    "Table_02_Main_Model_Results.csv"
  ),
  file.path(
    DATA_DIR,
    "Species_Richness_Manuscript_Input.rds"
  )
)

output_check <- file.exists(
  expected_output_files
)

cat("\n============================================================\n")
cat("OUTPUT FILE CHECK\n")
cat("============================================================\n\n")

print(
  data.frame(
    file = expected_output_files,
    exists = output_check
  ),
  row.names = FALSE
)

if (!all(output_check)) {
  stop(
    "One or more expected Script 02 outputs were not created."
  )
}

cat(
  "\nPASS: All Script 02 outputs were created.\n"
)


# ------------------------------------------------------------
# 48. Final analysis summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SPECIES RICHNESS ANALYSIS COMPLETE\n")
cat("============================================================\n\n")

cat(
  "Eligible cities:               ",
  nrow(eligible_cities),
  "\n",
  sep = ""
)

cat(
  "Sampled cities:                ",
  n_distinct(sr_data$NAME10),
  "\n",
  sep = ""
)

cat(
  "Sampled localities:            ",
  n_distinct(sr_data$locality_id),
  "\n",
  sep = ""
)

cat(
  "Sampled checklists:            ",
  nrow(sr_data),
  "\n",
  sep = ""
)

cat(
  "Poisson dispersion:            ",
  round(poisson_dispersion, 3),
  "\n",
  sep = ""
)

cat(
  "Negative-binomial dispersion:  ",
  round(nb_dispersion, 3),
  "\n",
  sep = ""
)

cat(
  "Poisson AIC:                   ",
  round(AIC(sr_model_poisson), 1),
  "\n",
  sep = ""
)

cat(
  "Negative-binomial AIC:         ",
  round(AIC(sr_model_nb), 1),
  "\n",
  sep = ""
)

cat(
  "Deviance explained:            ",
  round(
    100 * sr_summary$dev.expl,
    2
  ),
  "%\n",
  sep = ""
)

cat(
  "Adjusted R-squared:            ",
  round(sr_summary$r.sq, 3),
  "\n",
  sep = ""
)

cat(
  "Negative-binomial theta:       ",
  round(sr_theta, 3),
  "\n",
  sep = ""
)

cat(
  "Basis-dimension follow-up:     ",
  k_followup_flag,
  "\n",
  sep = ""
)

cat(
  "\nPublication figure:\n",
  file.path(
    FIGURE_DIR,
    "Figure_02_Main_Data_Trends.png"
  ),
  "\n",
  sep = ""
)

cat(
  "\nPublication table:\n",
  file.path(
    TABLE_DIR,
    "Table_02_Main_Model_Results.csv"
  ),
  "\n",
  sep = ""
)

cat(
  "\nManuscript input:\n",
  species_richness_manuscript_file,
  "\n",
  sep = ""
)


# ============================================================
# OVERALL INTERPRETATION
# ============================================================
#
# Model performance:
#
# - The balanced sampling design retained 273 cities,
#   1,365 localities, and 6,825 checklists.
#
# - The negative-binomial BAM was preferred over the Poisson
#   BAM and provided substantially better dispersion behavior.
#
# - The final negative-binomial BAM explained approximately
#   69.5% of deviance, with adjusted R-squared approximately
#   0.65.
#
# Ecological interpretation:
#
# - Aridity alone provides relatively weak evidence for a
#   direct relationship with species richness after accounting
#   for local habitat, sampling effort, and hierarchical
#   structure.
#
# - Water cover shows a supported nonlinear relationship with
#   richness.
#
# - Shrub and scrub cover shows a supported nonlinear
#   relationship with richness.
#
# - Built cover shows a supported nonlinear relationship with
#   richness.
#
# - The aridity x water interaction is supported, indicating
#   that the water-richness relationship varies along the
#   regional aridity gradient.
#
# - The aridity x shrub/scrub interaction is supported,
#   indicating that the shrub-richness relationship varies
#   along the regional aridity gradient.
#
# - The aridity x built interaction is unsupported.
#
# - Checklist duration and distance traveled remain important
#   sampling-effort controls.
#
# - City and locality random effects account for hierarchical
#   variation in checklist-level richness.
#
# Figure interpretation:
#
# - Panel A shows the city-level empirical aridity-richness
#   relationship.
#
# - Panels B-D show locality-level habitat-richness
#   relationships, with point color representing the percentile
#   rank of the city-level aridity index.
#
# - Warm colors indicate more-arid cities and cool colors
#   indicate less-arid cities. Percentile rank is used only for
#   visualization; the BAM uses the original continuous aridity
#   index.
#
# - The descriptive GAM lines are not the source of inferential
#   tests; inference comes from the full negative-binomial BAM.
#
# Overall:
#
# Species richness appears to be structured primarily by local
# habitat and by modification of water and shrub relationships
# across the regional aridity gradient, rather than by a strong
# independent aridity effect alone.
# ============================================================
