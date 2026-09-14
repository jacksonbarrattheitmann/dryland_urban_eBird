# ============================================================
# Dryland Bird Project
# 03_Functional_Richness.R
#
# Purpose:
# Model checklist-level functional richness (FRic) as a
# function of aridity, local habitat, and checklist effort
# using a GAM/BAM framework and balanced hierarchical sampling.
#
# Sampling design:
#
#   1. Retain Stationary and Traveling checklists
#   2. Duration <= 300 minutes
#   3. Distance <= 5 km
#   4. Retain checklists with non-missing FRic
#   5. Require >= 5 FRic checklists per locality
#   6. Require >= 5 eligible localities per city
#   7. Randomly sample 5 localities per city
#   8. Randomly sample 5 checklists per locality
#
# Candidate models:
#
#   A. Additive habitat model
#   B. Aridity x habitat interaction model
#
# The interaction model is retained as the primary/global model
# because the aridity x habitat interactions represent the
# a priori ecological hypotheses and allow direct comparison
# among biodiversity responses.
#
# Publication figures:
#
#   1. FRic response distribution
#   2. Conditional effects of:
#        - Aridity
#        - Water cover
#        - Shrub/scrub cover
#
# Canonical input:
#
# Intermediate_data/dfs_for_analysis_combined/combined_data.RDS
#
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(mgcv)


# ------------------------------------------------------------
# 2. Analysis settings
# ------------------------------------------------------------

MIN_CHECKLISTS_PER_LOCALITY <- 5

MIN_LOCALITIES_PER_CITY <- 5

N_LOCALITIES_PER_CITY <- 5

N_CHECKLISTS_PER_LOCALITY <- 5


# Same sampling seeds used across biodiversity responses.

LOCALITY_SAMPLE_SEED <- 20260911

CHECKLIST_SAMPLE_SEED <- 20260912


# Temporary duplicate handling.

TEMPORARY_DEDUPLICATE <- TRUE


# Figures are returned as R objects by default.
# Change to TRUE only when ready to save final figures.

SAVE_FIGURES <- FALSE


# ------------------------------------------------------------
# 3. Load combined dataset
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


combined_data <- readRDS(
  data_path
)


cat("\n============================================================\n")
cat("FUNCTIONAL RICHNESS ANALYSIS\n")
cat("============================================================\n\n")


cat(
  "Rows in combined_data: ",
  nrow(combined_data),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 4. Required-variable check
# ------------------------------------------------------------

required_variables <- c(
  "checklist_id",
  "locality_id",
  "NAME10",
  "state_name",
  "protocol_name",
  "duration_minutes",
  "effort_distance_km",
  "FRic",
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
    paste(
      missing_variables,
      collapse = "\n"
    )
  )
}


cat(
  "PASS: All required FRic analysis variables are present.\n"
)


# ------------------------------------------------------------
# 5. Duplicate checklist handling
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
  
  stop(
    "Duplicate checklist IDs are present."
  )
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
# 6. Apply protocol and effort filters
# ------------------------------------------------------------

fric_pool <- analysis_base %>%
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


fric_pool %>%
  summarise(
    n_checklists = n(),
    n_localities = n_distinct(locality_id),
    n_cities = n_distinct(NAME10)
  ) %>%
  print()


# ------------------------------------------------------------
# 7. Examine FRic before removing missing values
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FRIC INITIAL QC\n")
cat("============================================================\n\n")


cat(
  "FRic missing: ",
  sum(
    is.na(
      fric_pool$FRic
    )
  ),
  "\n",
  sep = ""
)


cat(
  "FRic missing (%): ",
  round(
    100 *
      mean(
        is.na(
          fric_pool$FRic
        )
      ),
    2
  ),
  "%\n",
  sep = ""
)


cat(
  "FRic zero values: ",
  sum(
    fric_pool$FRic == 0,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "FRic values equal to 1: ",
  sum(
    fric_pool$FRic == 1,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "FRic values < 0: ",
  sum(
    fric_pool$FRic < 0,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "FRic values > 1: ",
  sum(
    fric_pool$FRic > 1,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


print(
  summary(
    fric_pool$FRic
  )
)


# ------------------------------------------------------------
# 8. Remove undefined FRic values
# ------------------------------------------------------------
#
# FRic is undefined for some low-richness checklists.
#
# FRic-specific sampling eligibility is therefore determined
# after removing missing FRic values.
# ------------------------------------------------------------

fric_pool <- fric_pool %>%
  filter(
    !is.na(FRic)
  )


cat(
  "\nRows with defined FRic: ",
  nrow(fric_pool),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Complete-case restriction
# ------------------------------------------------------------

model_variables <- c(
  "FRic",
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
  variable =
    model_variables,
  
  missing_n =
    sapply(
      fric_pool[
        model_variables
      ],
      function(x) {
        sum(
          is.na(x)
        )
      }
    )
)


cat("\n============================================================\n")
cat("MODEL-VARIABLE MISSINGNESS\n")
cat("============================================================\n\n")


print(
  missing_summary
)


fric_pool <- fric_pool %>%
  filter(
    if_all(
      all_of(
        model_variables
      ),
      ~ !is.na(.)
    )
  )


# ------------------------------------------------------------
# 10. Identify eligible localities
# ------------------------------------------------------------

locality_sampling <- fric_pool %>%
  count(
    NAME10,
    locality_id,
    name = "n_checklists"
  )


eligible_localities <- locality_sampling %>%
  filter(
    n_checklists >=
      MIN_CHECKLISTS_PER_LOCALITY
  )


# ------------------------------------------------------------
# 11. Identify eligible cities
# ------------------------------------------------------------

city_sampling <- eligible_localities %>%
  count(
    NAME10,
    name = "n_eligible_localities"
  )


eligible_cities <- city_sampling %>%
  filter(
    n_eligible_localities >=
      MIN_LOCALITIES_PER_CITY
  )


cat("\n============================================================\n")
cat("FRIC ELIGIBILITY\n")
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
# 12. Create eligible sampling pool
# ------------------------------------------------------------

eligible_pool <- fric_pool %>%
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
# 13. Sample five localities per city
# ------------------------------------------------------------

set.seed(
  LOCALITY_SAMPLE_SEED
)


selected_localities <- eligible_localities %>%
  semi_join(
    eligible_cities,
    by = "NAME10"
  ) %>%
  arrange(
    NAME10,
    locality_id
  ) %>%
  group_by(
    NAME10
  ) %>%
  slice_sample(
    n = N_LOCALITIES_PER_CITY
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 14. Sample five FRic checklists per locality
# ------------------------------------------------------------

set.seed(
  CHECKLIST_SAMPLE_SEED
)


fric_data <- eligible_pool %>%
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
# 15. Balanced sample QC
# ------------------------------------------------------------

locality_qc <- fric_data %>%
  count(
    NAME10,
    locality_id,
    name = "n_checklists"
  )


city_qc <- fric_data %>%
  group_by(
    NAME10
  ) %>%
  summarise(
    n_localities =
      n_distinct(
        locality_id
      ),
    
    n_checklists =
      n(),
    
    .groups = "drop"
  )


cat("\n============================================================\n")
cat("FRIC BALANCED SAMPLE QC\n")
cat("============================================================\n\n")


cat(
  "Sampled cities:     ",
  n_distinct(
    fric_data$NAME10
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled localities: ",
  n_distinct(
    fric_data$locality_id
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled checklists: ",
  nrow(
    fric_data
  ),
  "\n",
  sep = ""
)


cat(
  "Localities with != 5 checklists: ",
  sum(
    locality_qc$n_checklists !=
      N_CHECKLISTS_PER_LOCALITY
  ),
  "\n",
  sep = ""
)


cat(
  "Cities with != 5 localities: ",
  sum(
    city_qc$n_localities !=
      N_LOCALITIES_PER_CITY
  ),
  "\n",
  sep = ""
)


cat(
  "Cities with != 25 checklists: ",
  sum(
    city_qc$n_checklists !=
      N_LOCALITIES_PER_CITY *
      N_CHECKLISTS_PER_LOCALITY
  ),
  "\n",
  sep = ""
)


if (
  any(
    locality_qc$n_checklists !=
    N_CHECKLISTS_PER_LOCALITY
  ) ||
  any(
    city_qc$n_localities !=
    N_LOCALITIES_PER_CITY
  ) ||
  any(
    city_qc$n_checklists !=
    N_LOCALITIES_PER_CITY *
    N_CHECKLISTS_PER_LOCALITY
  )
) {
  
  stop(
    "FRic balanced sampling QC failed."
  )
}


cat(
  "\nPASS: FRic balanced sampling verified.\n"
)


# ------------------------------------------------------------
# 16. Prepare factors
# ------------------------------------------------------------

fric_data <- fric_data %>%
  mutate(
    NAME10 =
      factor(
        NAME10
      ),
    
    locality_id =
      factor(
        locality_id
      ),
    
    protocol_name =
      factor(
        protocol_name,
        levels = c(
          "Stationary",
          "Traveling"
        )
      )
  )


# ------------------------------------------------------------
# 17. Examine final FRic response distribution
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL FRIC DISTRIBUTION\n")
cat("============================================================\n\n")


print(
  summary(
    fric_data$FRic
  )
)


cat(
  "\nZeros: ",
  sum(
    fric_data$FRic == 0
  ),
  "\n",
  sep = ""
)


cat(
  "Ones: ",
  sum(
    fric_data$FRic == 1
  ),
  "\n",
  sep = ""
)


cat(
  "Minimum: ",
  min(
    fric_data$FRic
  ),
  "\n",
  sep = ""
)


cat(
  "Maximum: ",
  max(
    fric_data$FRic
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 18. Confirm beta-family requirements
# ------------------------------------------------------------

n_boundary <- sum(
  fric_data$FRic <= 0 |
    fric_data$FRic >= 1
)


cat("\n============================================================\n")
cat("BETA FAMILY QC\n")
cat("============================================================\n\n")


cat(
  "FRic observations outside (0,1): ",
  n_boundary,
  "\n",
  sep = ""
)


if (
  n_boundary > 0
) {
  
  stop(
    paste0(
      "FRic contains boundary values. ",
      "Do not fit the beta BAM until these values are reviewed."
    )
  )
}


cat(
  "PASS: FRic is strictly bounded between 0 and 1.\n"
)


# ============================================================
# CANDIDATE BAM MODELS
# ============================================================


# ------------------------------------------------------------
# 19. Additive model formula
# ------------------------------------------------------------

fric_formula_additive <- FRic ~
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
# 20. Interaction model formula
# ------------------------------------------------------------

fric_formula_interaction <- FRic ~
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
    k = c(
      5,
      5
    )
  ) +
  ti(
    avg_aridity_index,
    shrub_and_scrub,
    k = c(
      5,
      5
    )
  ) +
  ti(
    avg_aridity_index,
    built,
    k = c(
      5,
      5
    )
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
# 21. Fit additive beta BAM
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FITTING ADDITIVE FRIC BAM\n")
cat("============================================================\n\n")


additive_fit_time <- system.time(
  
  fric_model_additive <- bam(
    formula =
      fric_formula_additive,
    
    family =
      betar(
        link = "logit"
      ),
    
    data =
      fric_data,
    
    method =
      "fREML",
    
    discrete =
      TRUE
  )
  
)


cat("\nAdditive BAM fit time:\n")


print(
  additive_fit_time
)


# ------------------------------------------------------------
# 22. Fit interaction beta BAM
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FITTING INTERACTION FRIC BAM\n")
cat("============================================================\n\n")


interaction_fit_time <- system.time(
  
  fric_model_interaction <- bam(
    formula =
      fric_formula_interaction,
    
    family =
      betar(
        link = "logit"
      ),
    
    data =
      fric_data,
    
    method =
      "fREML",
    
    discrete =
      TRUE
  )
  
)


cat("\nInteraction BAM fit time:\n")


print(
  interaction_fit_time
)


# ------------------------------------------------------------
# 23. Candidate-model summaries
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ADDITIVE MODEL SUMMARY\n")
cat("============================================================\n\n")


print(
  summary(
    fric_model_additive
  )
)


cat("\n============================================================\n")
cat("INTERACTION MODEL SUMMARY\n")
cat("============================================================\n\n")


print(
  summary(
    fric_model_interaction
  )
)


# ------------------------------------------------------------
# 24. Compare candidate models
# ------------------------------------------------------------

additive_summary <- summary(
  fric_model_additive
)


interaction_summary <- summary(
  fric_model_interaction
)


fric_model_comparison <- data.frame(
  model = c(
    "Additive",
    "Aridity x habitat interactions"
  ),
  
  AIC = c(
    AIC(
      fric_model_additive
    ),
    
    AIC(
      fric_model_interaction
    )
  ),
  
  deviance_explained = c(
    additive_summary$dev.expl,
    interaction_summary$dev.expl
  ),
  
  adjusted_R2 = c(
    additive_summary$r.sq,
    interaction_summary$r.sq
  )
)


fric_model_comparison <- fric_model_comparison %>%
  mutate(
    delta_AIC =
      AIC -
      min(
        AIC
      )
  )


cat("\n============================================================\n")
cat("FRIC MODEL COMPARISON\n")
cat("============================================================\n\n")


print(
  fric_model_comparison
)


# ------------------------------------------------------------
# 25. Candidate-model convergence
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CANDIDATE MODEL CONVERGENCE\n")
cat("============================================================\n\n")


cat(
  "Additive model converged: ",
  fric_model_additive$converged,
  "\n",
  sep = ""
)


cat(
  "Interaction model converged: ",
  fric_model_interaction$converged,
  "\n",
  sep = ""
)


cat(
  "Additive model rank: ",
  fric_model_additive$rank,
  " / ",
  length(
    coef(
      fric_model_additive
    )
  ),
  "\n",
  sep = ""
)


cat(
  "Interaction model rank: ",
  fric_model_interaction$rank,
  " / ",
  length(
    coef(
      fric_model_interaction
    )
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 26. Preliminary basis checks
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ADDITIVE MODEL K CHECK\n")
cat("============================================================\n\n")


print(
  k.check(
    fric_model_additive
  )
)


cat("\n============================================================\n")
cat("INTERACTION MODEL K CHECK\n")
cat("============================================================\n\n")


print(
  k.check(
    fric_model_interaction
  )
)


# ============================================================
# FINAL FRIC MODEL
# ============================================================


# ------------------------------------------------------------
# 27. Select primary/global FRic model
# ------------------------------------------------------------
#
# The interaction model is retained as the primary FRic model.
#
# The additive and interaction models have essentially
# equivalent overall fit.
#
# The aridity x habitat relationships represent the a priori
# ecological hypotheses for the project. Retaining the same
# global interaction structure across biodiversity responses
# allows direct comparison of whether these interactions are
# supported for species richness, FRic, FDis, and FEve.
#
# Individual interaction terms are interpreted from the global
# fitted model rather than removed solely because they are not
# statistically significant.
# ------------------------------------------------------------

fric_model_final <- fric_model_interaction


fric_summary <- summary(
  fric_model_final
)


cat("\n============================================================\n")
cat("FINAL FRIC MODEL\n")
cat("============================================================\n\n")


print(
  fric_summary
)


# ============================================================
# FINAL FRIC QC AND DIAGNOSTICS
# ============================================================


# ------------------------------------------------------------
# 28. Final model convergence and rank
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FRIC MODEL CONVERGENCE\n")
cat("============================================================\n\n")


cat(
  "Model converged: ",
  fric_model_final$converged,
  "\n",
  sep = ""
)


cat(
  "Model rank: ",
  fric_model_final$rank,
  " / ",
  length(
    coef(
      fric_model_final
    )
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 29. Basis-dimension check
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FRIC BASIS DIMENSION CHECK\n")
cat("============================================================\n\n")


fric_k_check <- k.check(
  fric_model_final
)


print(
  fric_k_check
)


# ------------------------------------------------------------
# 30. GAM variance components
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FRIC VARIANCE COMPONENTS\n")
cat("============================================================\n\n")


fric_vcomp <- gam.vcomp(
  fric_model_final
)


print(
  fric_vcomp
)


# ------------------------------------------------------------
# 31. Full concurvity
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FRIC FULL CONCURVITY\n")
cat("============================================================\n\n")


fric_concurvity_full <- concurvity(
  fric_model_final,
  full = TRUE
)


print(
  fric_concurvity_full
)


# ------------------------------------------------------------
# 32. Pairwise concurvity
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FRIC PAIRWISE CONCURVITY\n")
cat("============================================================\n\n")


fric_concurvity_pairwise <- concurvity(
  fric_model_final,
  full = FALSE
)


fric_pairwise_estimate <-
  fric_concurvity_pairwise$estimate


print(
  round(
    fric_pairwise_estimate,
    3
  )
)


# ------------------------------------------------------------
# 33. Ecological-term concurvity
# ------------------------------------------------------------

fric_ecological_terms <- c(
  "s(avg_aridity_index)",
  "s(water)",
  "s(shrub_and_scrub)",
  "s(built)",
  "ti(avg_aridity_index,water)",
  "ti(avg_aridity_index,shrub_and_scrub)",
  "ti(avg_aridity_index,built)",
  "s(duration_minutes)",
  "s(effort_distance_km)"
)


fric_ecological_concurvity <-
  fric_pairwise_estimate[
    fric_ecological_terms,
    fric_ecological_terms,
    drop = FALSE
  ]


cat("\n============================================================\n")
cat("FRIC ECOLOGICAL-TERM CONCURVITY\n")
cat("============================================================\n\n")


print(
  round(
    fric_ecological_concurvity,
    3
  )
)


# ------------------------------------------------------------
# 34. Residual summaries
# ------------------------------------------------------------

fric_pearson_residuals <- residuals(
  fric_model_final,
  type = "pearson"
)


fric_deviance_residuals <- residuals(
  fric_model_final,
  type = "deviance"
)


fric_fitted_values <- fitted(
  fric_model_final
)


cat("\n============================================================\n")
cat("FRIC RESIDUAL SUMMARY\n")
cat("============================================================\n\n")


cat(
  "Pearson residuals:\n"
)


print(
  summary(
    fric_pearson_residuals
  )
)


cat(
  "\nDeviance residuals:\n"
)


print(
  summary(
    fric_deviance_residuals
  )
)


# ------------------------------------------------------------
# 35. Extreme Pearson residuals
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FRIC EXTREME RESIDUALS\n")
cat("============================================================\n\n")


cat(
  "Absolute Pearson residual > 2: ",
  sum(
    abs(
      fric_pearson_residuals
    ) > 2
  ),
  " (",
  round(
    100 *
      mean(
        abs(
          fric_pearson_residuals
        ) > 2
      ),
    2
  ),
  "%)\n",
  sep = ""
)


cat(
  "Absolute Pearson residual > 3: ",
  sum(
    abs(
      fric_pearson_residuals
    ) > 3
  ),
  " (",
  round(
    100 *
      mean(
        abs(
          fric_pearson_residuals
        ) > 3
      ),
    2
  ),
  "%)\n",
  sep = ""
)


# ------------------------------------------------------------
# 36. Residuals versus fitted
# ------------------------------------------------------------

plot(
  fric_fitted_values,
  fric_pearson_residuals,
  pch = 16,
  cex = 0.5,
  xlab = "Fitted functional richness",
  ylab = "Pearson residuals"
)


abline(
  h = 0,
  lty = 2
)


# ------------------------------------------------------------
# 37. QQ diagnostic
# ------------------------------------------------------------

qq.gam(
  fric_model_final,
  rep = 100,
  level = 0.90
)


# ------------------------------------------------------------
# 38. Standard GAM diagnostics
# ------------------------------------------------------------

gam.check(
  fric_model_final
)


# ------------------------------------------------------------
# 39. Final parametric-term table
# ------------------------------------------------------------

fric_parametric_results <- as.data.frame(
  fric_summary$p.table
)


fric_parametric_results$term <- rownames(
  fric_parametric_results
)


fric_parametric_results <- fric_parametric_results %>%
  select(
    term,
    everything()
  )


cat("\n============================================================\n")
cat("FINAL FRIC PARAMETRIC TERMS\n")
cat("============================================================\n\n")


print(
  fric_parametric_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 40. Final smooth-term table
# ------------------------------------------------------------

fric_smooth_results <- as.data.frame(
  fric_summary$s.table
)


fric_smooth_results$term <- rownames(
  fric_smooth_results
)


fric_smooth_results <- fric_smooth_results %>%
  select(
    term,
    everything()
  )


cat("\n============================================================\n")
cat("FINAL FRIC SMOOTH TERMS\n")
cat("============================================================\n\n")


print(
  fric_smooth_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 41. Final model fit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL FRIC MODEL FIT\n")
cat("============================================================\n\n")


cat(
  "Adjusted R-squared: ",
  round(
    fric_summary$r.sq,
    3
  ),
  "\n",
  sep = ""
)


cat(
  "Deviance explained: ",
  round(
    100 *
      fric_summary$dev.expl,
    2
  ),
  "%\n",
  sep = ""
)


cat(
  "Beta precision parameter: ",
  round(
    fric_model_final$family$getTheta(
      TRUE
    ),
    3
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 42. Diagnostic smooth plots
# ------------------------------------------------------------
#
# Smooth order:
#
#   1  Aridity
#   2  Water
#   3  Shrub/scrub
#   4  Built
#   5  Aridity x water
#   6  Aridity x shrub/scrub
#   7  Aridity x built
#   8  Duration
#   9  Distance
#
# Random-effect smooths are not plotted here.
# ------------------------------------------------------------

plot(
  fric_model_final,
  select = 1,
  shade = TRUE,
  pages = 1
)


plot(
  fric_model_final,
  select = 2,
  shade = TRUE,
  pages = 1
)


plot(
  fric_model_final,
  select = 3,
  shade = TRUE,
  pages = 1
)


plot(
  fric_model_final,
  select = 4,
  shade = TRUE,
  pages = 1
)


plot(
  fric_model_final,
  select = 5,
  pages = 1
)


plot(
  fric_model_final,
  select = 6,
  pages = 1
)


plot(
  fric_model_final,
  select = 7,
  pages = 1
)


plot(
  fric_model_final,
  select = 8,
  shade = TRUE,
  pages = 1
)


plot(
  fric_model_final,
  select = 9,
  shade = TRUE,
  pages = 1
)


# ============================================================
# PUBLICATION-STYLE FIGURES
# ============================================================


# ------------------------------------------------------------
# 43. Reference values for population-level predictions
# ------------------------------------------------------------

reference_aridity <- median(
  fric_data$avg_aridity_index,
  na.rm = TRUE
)


reference_duration <- median(
  fric_data$duration_minutes,
  na.rm = TRUE
)


reference_distance <- median(
  fric_data$effort_distance_km,
  na.rm = TRUE
)


reference_water <- median(
  fric_data$water,
  na.rm = TRUE
)


reference_shrub <- median(
  fric_data$shrub_and_scrub,
  na.rm = TRUE
)


reference_built <- median(
  fric_data$built,
  na.rm = TRUE
)


reference_protocol <- factor(
  "Stationary",
  levels = levels(
    fric_data$protocol_name
  )
)


reference_city <- factor(
  levels(
    fric_data$NAME10
  )[1],
  levels = levels(
    fric_data$NAME10
  )
)


reference_locality <- factor(
  levels(
    fric_data$locality_id
  )[1],
  levels = levels(
    fric_data$locality_id
  )
)


# ------------------------------------------------------------
# 44. Determine habitat scaling
# ------------------------------------------------------------

maximum_habitat_value <- max(
  c(
    fric_data$water,
    fric_data$shrub_and_scrub,
    fric_data$built
  ),
  na.rm = TRUE
)


if (
  maximum_habitat_value <= 1.05
) {
  
  habitat_scale_factor <- 100
  
} else {
  
  habitat_scale_factor <- 1
}


# ------------------------------------------------------------
# 45. Prediction helper
# ------------------------------------------------------------
#
# Predictions and confidence intervals are calculated on the
# beta-model logit scale and transformed back to the original
# FRic scale.
#
# City and locality random effects are excluded so the curves
# represent population-level relationships.
# ------------------------------------------------------------

predict_fric <- function(
    newdata
) {
  
  prediction <- predict(
    fric_model_final,
    newdata = newdata,
    type = "link",
    se.fit = TRUE,
    exclude = c(
      "s(NAME10)",
      "s(locality_id)"
    )
  )
  
  
  link_inverse <-
    fric_model_final$family$linkinv
  
  
  fit_link <-
    as.numeric(
      prediction$fit
    )
  
  
  se_link <-
    as.numeric(
      prediction$se.fit
    )
  
  
  data.frame(
    fit =
      link_inverse(
        fit_link
      ),
    
    lower =
      link_inverse(
        fit_link -
          1.96 *
          se_link
      ),
    
    upper =
      link_inverse(
        fit_link +
          1.96 *
          se_link
      )
  )
}


# ------------------------------------------------------------
# 46. Aridity predictions
# ------------------------------------------------------------

aridity_limits <- quantile(
  fric_data$avg_aridity_index,
  probs = c(
    0.05,
    0.95
  ),
  na.rm = TRUE
)


aridity_sequence <- seq(
  aridity_limits[1],
  aridity_limits[2],
  length.out = 150
)


aridity_prediction_data <- data.frame(
  avg_aridity_index =
    aridity_sequence,
  
  water =
    reference_water,
  
  shrub_and_scrub =
    reference_shrub,
  
  built =
    reference_built,
  
  duration_minutes =
    reference_duration,
  
  effort_distance_km =
    reference_distance,
  
  protocol_name =
    reference_protocol,
  
  NAME10 =
    reference_city,
  
  locality_id =
    reference_locality
)


aridity_prediction_results <- predict_fric(
  aridity_prediction_data
)


aridity_plot_data <- data.frame(
  predictor =
    "Aridity index\n(higher = less arid)",
  
  x_value =
    aridity_prediction_data$avg_aridity_index,
  
  fit =
    aridity_prediction_results$fit,
  
  lower =
    aridity_prediction_results$lower,
  
  upper =
    aridity_prediction_results$upper
)


# ------------------------------------------------------------
# 47. Water-cover predictions
# ------------------------------------------------------------

water_limits <- quantile(
  fric_data$water,
  probs = c(
    0.05,
    0.95
  ),
  na.rm = TRUE
)


water_sequence <- seq(
  water_limits[1],
  water_limits[2],
  length.out = 150
)


water_prediction_data <- data.frame(
  avg_aridity_index =
    reference_aridity,
  
  water =
    water_sequence,
  
  shrub_and_scrub =
    reference_shrub,
  
  built =
    reference_built,
  
  duration_minutes =
    reference_duration,
  
  effort_distance_km =
    reference_distance,
  
  protocol_name =
    reference_protocol,
  
  NAME10 =
    reference_city,
  
  locality_id =
    reference_locality
)


water_prediction_results <- predict_fric(
  water_prediction_data
)


water_plot_data <- data.frame(
  predictor =
    "Water cover (%)",
  
  x_value =
    water_prediction_data$water *
    habitat_scale_factor,
  
  fit =
    water_prediction_results$fit,
  
  lower =
    water_prediction_results$lower,
  
  upper =
    water_prediction_results$upper
)


# ------------------------------------------------------------
# 48. Shrub/scrub predictions
# ------------------------------------------------------------

shrub_limits <- quantile(
  fric_data$shrub_and_scrub,
  probs = c(
    0.05,
    0.95
  ),
  na.rm = TRUE
)


shrub_sequence <- seq(
  shrub_limits[1],
  shrub_limits[2],
  length.out = 150
)


shrub_prediction_data <- data.frame(
  avg_aridity_index =
    reference_aridity,
  
  water =
    reference_water,
  
  shrub_and_scrub =
    shrub_sequence,
  
  built =
    reference_built,
  
  duration_minutes =
    reference_duration,
  
  effort_distance_km =
    reference_distance,
  
  protocol_name =
    reference_protocol,
  
  NAME10 =
    reference_city,
  
  locality_id =
    reference_locality
)


shrub_prediction_results <- predict_fric(
  shrub_prediction_data
)


shrub_plot_data <- data.frame(
  predictor =
    "Shrub and scrub cover (%)",
  
  x_value =
    shrub_prediction_data$shrub_and_scrub *
    habitat_scale_factor,
  
  fit =
    shrub_prediction_results$fit,
  
  lower =
    shrub_prediction_results$lower,
  
  upper =
    shrub_prediction_results$upper
)


# ------------------------------------------------------------
# 49. Combine conditional-effect data
# ------------------------------------------------------------

fric_conditional_effects <- bind_rows(
  aridity_plot_data,
  water_plot_data,
  shrub_plot_data
) %>%
  mutate(
    predictor =
      factor(
        predictor,
        levels = c(
          "Aridity index\n(higher = less arid)",
          "Water cover (%)",
          "Shrub and scrub cover (%)"
        )
      )
  )


# ------------------------------------------------------------
# 50. Main publication figure
# ------------------------------------------------------------

fig_fric_main <- ggplot(
  fric_conditional_effects,
  aes(
    x = x_value,
    y = fit
  )
) +
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    fill = "grey80",
    alpha = 0.7
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  facet_wrap(
    ~ predictor,
    scales = "free_x",
    nrow = 1
  ) +
  labs(
    x = NULL,
    y = "Predicted functional richness"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    strip.background =
      element_blank(),
    
    strip.text =
      element_text(
        face = "bold",
        size = 11
      )
  )


print(
  fig_fric_main
)


# ------------------------------------------------------------
# 51. FRic response-distribution figure
# ------------------------------------------------------------

fig_fric_distribution <- ggplot(
  fric_data,
  aes(
    x = FRic
  )
) +
  geom_histogram(
    bins = 40,
    linewidth = 0.25
  ) +
  labs(
    x = "Functional richness",
    y = "Number of checklists"
  ) +
  theme_classic(
    base_size = 12
  )


print(
  fig_fric_distribution
)


# ------------------------------------------------------------
# 52. Optional figure saving
# ------------------------------------------------------------

if (
  SAVE_FIGURES
) {
  
  figure_directory <- file.path(
    "Outputs",
    "03_Functional_Richness",
    "Figures"
  )
  
  
  dir.create(
    figure_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  ggsave(
    filename =
      file.path(
        figure_directory,
        "FRic_Conditional_Effects.pdf"
      ),
    plot =
      fig_fric_main,
    width = 9,
    height = 3.8
  )
  
  
  ggsave(
    filename =
      file.path(
        figure_directory,
        "FRic_Distribution.pdf"
      ),
    plot =
      fig_fric_distribution,
    width = 5.5,
    height = 4
  )
}


# ============================================================
# FINAL QC AND ANALYSIS SUMMARY
# ============================================================


# ------------------------------------------------------------
# 53. Final automated QC summary
# ------------------------------------------------------------

qc_balanced_sample <-
  all(
    locality_qc$n_checklists ==
      N_CHECKLISTS_PER_LOCALITY
  ) &&
  all(
    city_qc$n_localities ==
      N_LOCALITIES_PER_CITY
  ) &&
  all(
    city_qc$n_checklists ==
      N_LOCALITIES_PER_CITY *
      N_CHECKLISTS_PER_LOCALITY
  )


qc_beta_response <-
  all(
    fric_data$FRic > 0 &
      fric_data$FRic < 1
  )


qc_converged <-
  isTRUE(
    fric_model_final$converged
  )


qc_full_rank <-
  fric_model_final$rank ==
  length(
    coef(
      fric_model_final
    )
  )


fric_qc_summary <- data.frame(
  check = c(
    "Balanced 5 x 5 sampling",
    "FRic values strictly inside (0,1)",
    "Final BAM converged",
    "Final BAM is full rank"
  ),
  
  pass = c(
    qc_balanced_sample,
    qc_beta_response,
    qc_converged,
    qc_full_rank
  )
)


cat("\n============================================================\n")
cat("FRIC AUTOMATED QC SUMMARY\n")
cat("============================================================\n\n")


print(
  fric_qc_summary,
  row.names = FALSE
)


cat(
  "\nChecks passed: ",
  sum(
    fric_qc_summary$pass
  ),
  " of ",
  nrow(
    fric_qc_summary
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 54. Final analysis summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FUNCTIONAL RICHNESS ANALYSIS COMPLETE\n")
cat("============================================================\n\n")


cat(
  "Eligible cities:                ",
  nrow(
    eligible_cities
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled cities:                 ",
  n_distinct(
    fric_data$NAME10
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled localities:             ",
  n_distinct(
    fric_data$locality_id
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled checklists:             ",
  nrow(
    fric_data
  ),
  "\n",
  sep = ""
)


cat(
  "Additive AIC:                   ",
  round(
    AIC(
      fric_model_additive
    ),
    2
  ),
  "\n",
  sep = ""
)


cat(
  "Interaction AIC:                ",
  round(
    AIC(
      fric_model_interaction
    ),
    2
  ),
  "\n",
  sep = ""
)


cat(
  "Final adjusted R-squared:       ",
  round(
    fric_summary$r.sq,
    3
  ),
  "\n",
  sep = ""
)


cat(
  "Final deviance explained:       ",
  round(
    100 *
      fric_summary$dev.expl,
    2
  ),
  "%\n",
  sep = ""
)


cat(
  "Beta precision parameter:       ",
  round(
    fric_model_final$family$getTheta(
      TRUE
    ),
    3
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 55. Overall interpretation
# ------------------------------------------------------------
#
# The primary FRic model uses the same global aridity x habitat
# interaction structure used across the biodiversity-response
# analyses.
#
# The additive and interaction models showed essentially
# equivalent overall model fit. The interaction model is
# retained as the primary/global model because the aridity x
# habitat relationships represent a priori ecological
# hypotheses and maintaining the same global structure
# facilitates comparison among species richness, FRic, FDis,
# and FEve.
#
# Final FRic results indicate:
#
# - Water cover shows the clearest relationship with functional
#   richness. The relationship is nonlinear and strongly
#   supported in the final beta BAM.
#
# - There is little evidence for an independent relationship
#   between aridity and FRic after accounting for local habitat,
#   effort, and the interaction terms.
#
# - Shrub/scrub cover shows comparatively weak or marginal
#   evidence of an association with functional richness.
#
# - Built cover shows little evidence of an independent
#   relationship with FRic.
#
# - None of the three aridity x habitat interaction surfaces
#   show strong statistical support. Thus, there is little
#   evidence that the effects of water, shrub/scrub, or built
#   cover on FRic change strongly along the regional aridity
#   gradient.
#
# - This differs from the species-richness analysis, where
#   aridity x water and aridity x shrub/scrub interactions were
#   supported.
#
# - Checklist duration, distance traveled, and protocol remain
#   important sampling-effort controls.
#
# - City and locality random-effect smooths account for
#   hierarchical variation in the data, with locality-level
#   variation greater than city-level variation.
#
# - The final beta BAM explains approximately 67.9% of
#   deviance, with an adjusted R-squared of approximately 0.56.
#
# - Model convergence, model rank, basis dimensions, response
#   support, concurvity, and residual diagnostics were
#   acceptable, and all automated QC checks passed.
#
# The publication-style conditional-effect figure emphasizes
# aridity, water cover, and shrub/scrub cover.
#
# Population-level predictions exclude city and locality random
# effects. Other habitat and sampling-effort variables are held
# at representative median values, and protocol is held at
# Stationary.1
#
# Built cover is omitted from the primary figure because the
# model provides little evidence for an independent built-cover
# effect. Aridity x habitat surfaces are also omitted from the
# primary visualization because none received strong statistical
# support.
#
# Overall, functional richness appears to be associated more
# strongly with local habitat structure, particularly water
# cover, than with the regional aridity gradient itself.
# Moreover, there is little evidence that local habitat-FRic
# relationships are strongly modified across the aridity
# gradient.
#
# ============================================================