# ============================================================
# Dryland Bird Project
# 04_Functional_Dispersion.R
#
# Purpose:
# Model checklist-level functional dispersion (FDis) as a
# function of aridity, local habitat, and checklist effort
# using a Gaussian BAM and balanced hierarchical sampling.
#
# IMPORTANT:
#
# FDis is restricted to multispecies assemblages (sr >= 2).
#
# Single-species checklists have FDis = 0 by definition because
# no among-species trait dispersion exists when only one
# species is present. Previous Tweedie QC showed that these
# structural zeros were not adequately represented by a
# single Tweedie process.
#
# Therefore this analysis asks:
#
#   Among multispecies assemblages, how does functional
#   dispersion vary with aridity, local habitat, and effort?
#
# Sampling design:
#
#   1. Retain Stationary and Traveling checklists
#   2. Stationary distance set to zero
#   3. Duration <= 300 minutes
#   4. Distance <= 5 km
#   5. Retain multispecies checklists (sr >= 2)
#   6. Retain non-missing FDis
#   7. Require >= 5 eligible checklists per locality
#   8. Require >= 5 eligible localities per city
#   9. Randomly sample 5 localities per city
#  10. Randomly sample 5 checklists per locality
#
# Candidate models:
#
#   A. Additive Gaussian BAM
#   B. Aridity x habitat interaction Gaussian BAM
#
# The interaction model represents the primary/global
# ecological model because the aridity x habitat relationships
# are part of the a priori hypotheses and the same global
# structure is used across biodiversity responses.
#
# Both city and locality random-effect smooths are retained.
#
# Publication figures:
#
#   1. FDis response distribution
#   2. Conditional aridity effect
#   3. Aridity x habitat prediction figure
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

K_CHECK_SEED <- 20260912


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
cat("FUNCTIONAL DISPERSION ANALYSIS\n")
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
  "sr",
  "FDis",
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
  "PASS: All required FDis analysis variables are present.\n"
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

fdis_pool <- analysis_base %>%
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


fdis_pool %>%
  summarise(
    n_checklists = n(),
    n_localities = n_distinct(locality_id),
    n_cities = n_distinct(NAME10)
  ) %>%
  print()


# ------------------------------------------------------------
# 7. Examine original FDis structure
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ORIGINAL FDIS STRUCTURE\n")
cat("============================================================\n\n")


cat(
  "Single-species checklists: ",
  sum(
    fdis_pool$sr == 1,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "Multispecies checklists: ",
  sum(
    fdis_pool$sr >= 2,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "FDis = 0 and sr = 1: ",
  sum(
    fdis_pool$FDis == 0 &
      fdis_pool$sr == 1,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "FDis = 0 and sr >= 2: ",
  sum(
    fdis_pool$FDis == 0 &
      fdis_pool$sr >= 2,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "FDis > 0 and sr = 1: ",
  sum(
    fdis_pool$FDis > 0 &
      fdis_pool$sr == 1,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Restrict analysis to multispecies assemblages
# ------------------------------------------------------------

fdis_pool <- fdis_pool %>%
  filter(
    sr >= 2,
    !is.na(FDis)
  )


cat("\n============================================================\n")
cat("MULTISPECIES FDIS DATA\n")
cat("============================================================\n\n")


cat(
  "Rows after sr >= 2 restriction: ",
  nrow(fdis_pool),
  "\n",
  sep = ""
)


cat(
  "Remaining FDis zeros: ",
  sum(
    fdis_pool$FDis == 0
  ),
  "\n",
  sep = ""
)


cat(
  "Remaining negative FDis values: ",
  sum(
    fdis_pool$FDis < 0
  ),
  "\n",
  sep = ""
)


print(
  summary(
    fdis_pool$FDis
  )
)


if (
  any(
    fdis_pool$FDis <= 0
  )
) {
  
  stop(
    paste0(
      "Non-positive FDis values remain after restricting ",
      "to sr >= 2. Review these observations before modeling."
    )
  )
}


cat(
  "\nPASS: Multispecies FDis values are strictly positive.\n"
)


# ------------------------------------------------------------
# 9. Complete-case restriction
# ------------------------------------------------------------

model_variables <- c(
  "FDis",
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
      fdis_pool[
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


fdis_pool <- fdis_pool %>%
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

locality_sampling <- fdis_pool %>%
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
cat("FDIS ELIGIBILITY\n")
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

eligible_pool <- fdis_pool %>%
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
# 14. Sample five checklists per locality
# ------------------------------------------------------------

set.seed(
  CHECKLIST_SAMPLE_SEED
)


fdis_data <- eligible_pool %>%
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

locality_qc <- fdis_data %>%
  count(
    NAME10,
    locality_id,
    name = "n_checklists"
  )


city_qc <- fdis_data %>%
  group_by(
    NAME10
  ) %>%
  summarise(
    n_localities =
      n_distinct(locality_id),
    
    n_checklists =
      n(),
    
    .groups = "drop"
  )


cat("\n============================================================\n")
cat("FDIS BALANCED SAMPLE QC\n")
cat("============================================================\n\n")


cat(
  "Sampled cities:     ",
  n_distinct(fdis_data$NAME10),
  "\n",
  sep = ""
)


cat(
  "Sampled localities: ",
  n_distinct(fdis_data$locality_id),
  "\n",
  sep = ""
)


cat(
  "Sampled checklists: ",
  nrow(fdis_data),
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
    "FDis balanced sampling QC failed."
  )
}


cat(
  "\nPASS: FDis balanced sampling verified.\n"
)


# ------------------------------------------------------------
# 16. Prepare factors
# ------------------------------------------------------------
#
# IMPORTANT:
# Random-effect smooths require grouping variables to be
# factors. This also fixes the mgcv error encountered in the
# previous draft.
# ------------------------------------------------------------

fdis_data <- fdis_data %>%
  mutate(
    NAME10 =
      factor(NAME10),
    
    locality_id =
      factor(locality_id),
    
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
# 17. Final FDis response distribution
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL MULTISPECIES FDIS DISTRIBUTION\n")
cat("============================================================\n\n")


print(
  summary(
    fdis_data$FDis
  )
)


cat(
  "\nMinimum: ",
  min(
    fdis_data$FDis
  ),
  "\n",
  sep = ""
)


cat(
  "Maximum: ",
  max(
    fdis_data$FDis
  ),
  "\n",
  sep = ""
)


cat(
  "Mean: ",
  mean(
    fdis_data$FDis
  ),
  "\n",
  sep = ""
)


cat(
  "SD: ",
  sd(
    fdis_data$FDis
  ),
  "\n",
  sep = ""
)


# Publication/supplementary response-distribution figure.

fig_fdis_distribution <- ggplot(
  fdis_data,
  aes(
    x = FDis
  )
) +
  geom_histogram(
    bins = 40,
    linewidth = 0.25
  ) +
  labs(
    x = "Functional dispersion",
    y = "Number of checklists"
  ) +
  theme_classic(
    base_size = 12
  )


print(
  fig_fdis_distribution
)


# ============================================================
# CANDIDATE GAUSSIAN BAM MODELS
# ============================================================


# ------------------------------------------------------------
# 18. Additive model formula
# ------------------------------------------------------------

fdis_formula_additive <- FDis ~
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
# 19. Interaction/global model formula
# ------------------------------------------------------------

fdis_formula_interaction <- FDis ~
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
# 20. Fit additive Gaussian BAM
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FITTING ADDITIVE FDIS GAUSSIAN BAM\n")
cat("============================================================\n\n")


additive_fit_time <- system.time(
  
  fdis_model_additive <- bam(
    formula =
      fdis_formula_additive,
    
    family =
      gaussian(
        link = "identity"
      ),
    
    data =
      fdis_data,
    
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
# 21. Fit interaction Gaussian BAM
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FITTING INTERACTION FDIS GAUSSIAN BAM\n")
cat("============================================================\n\n")


interaction_fit_time <- system.time(
  
  fdis_model_interaction <- bam(
    formula =
      fdis_formula_interaction,
    
    family =
      gaussian(
        link = "identity"
      ),
    
    data =
      fdis_data,
    
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
# 22. Candidate-model summaries
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ADDITIVE MODEL SUMMARY\n")
cat("============================================================\n\n")


print(
  summary(
    fdis_model_additive
  )
)


cat("\n============================================================\n")
cat("INTERACTION MODEL SUMMARY\n")
cat("============================================================\n\n")


print(
  summary(
    fdis_model_interaction
  )
)


# ------------------------------------------------------------
# 23. Candidate-model comparison
# ------------------------------------------------------------
#
# AIC is retained as a descriptive comparison.
#
# Because candidate models are fitted using fREML and differ
# in smooth structure, the model-selection decision is not
# based solely on the AIC difference.
#
# The interaction model represents the a priori global
# ecological hypothesis.
# ------------------------------------------------------------

additive_summary <- summary(
  fdis_model_additive
)


interaction_summary <- summary(
  fdis_model_interaction
)


fdis_model_comparison <- data.frame(
  model = c(
    "Additive",
    "Aridity x habitat interactions"
  ),
  
  AIC = c(
    AIC(
      fdis_model_additive
    ),
    
    AIC(
      fdis_model_interaction
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


fdis_model_comparison <- fdis_model_comparison %>%
  mutate(
    delta_AIC =
      AIC -
      min(AIC)
  )


cat("\n============================================================\n")
cat("FDIS MODEL COMPARISON\n")
cat("============================================================\n\n")


print(
  fdis_model_comparison
)


# ------------------------------------------------------------
# 24. Candidate-model convergence
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CANDIDATE MODEL CONVERGENCE\n")
cat("============================================================\n\n")


cat(
  "Additive model converged: ",
  fdis_model_additive$converged,
  "\n",
  sep = ""
)


cat(
  "Interaction model converged: ",
  fdis_model_interaction$converged,
  "\n",
  sep = ""
)


cat(
  "Additive model rank: ",
  fdis_model_additive$rank,
  " / ",
  length(
    coef(
      fdis_model_additive
    )
  ),
  "\n",
  sep = ""
)


cat(
  "Interaction model rank: ",
  fdis_model_interaction$rank,
  " / ",
  length(
    coef(
      fdis_model_interaction
    )
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 25. Preliminary basis checks
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ADDITIVE MODEL K CHECK\n")
cat("============================================================\n\n")


set.seed(
  K_CHECK_SEED
)


print(
  k.check(
    fdis_model_additive
  )
)


cat("\n============================================================\n")
cat("INTERACTION MODEL K CHECK\n")
cat("============================================================\n\n")


set.seed(
  K_CHECK_SEED
)


print(
  k.check(
    fdis_model_interaction
  )
)


# ============================================================
# FINAL FDIS MODEL
# ============================================================


# ------------------------------------------------------------
# 26. Select primary/global model
# ------------------------------------------------------------
#
# The interaction model is retained as the primary/global FDis
# model because the aridity x habitat interactions represent
# the a priori ecological hypotheses and maintaining the same
# global model structure facilitates direct comparison among
# species richness, FRic, FDis, and FEve.
#
# The additive model remains as a simpler comparison.
# ------------------------------------------------------------

fdis_model_final <- fdis_model_interaction


fdis_summary <- summary(
  fdis_model_final
)


cat("\n============================================================\n")
cat("FINAL FDIS MODEL\n")
cat("============================================================\n\n")


print(
  fdis_summary
)


# ============================================================
# FINAL MODEL QC
# ============================================================


# ------------------------------------------------------------
# 27. Convergence and model rank
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FDIS MODEL CONVERGENCE\n")
cat("============================================================\n\n")


cat(
  "Model converged: ",
  fdis_model_final$converged,
  "\n",
  sep = ""
)


cat(
  "Model rank: ",
  fdis_model_final$rank,
  " / ",
  length(
    coef(
      fdis_model_final
    )
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 28. Stable basis-dimension check
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FDIS STABLE K CHECK - 1000 REPLICATES\n")
cat("============================================================\n\n")


set.seed(
  K_CHECK_SEED
)


fdis_k_stable <- k.check(
  fdis_model_final,
  n.rep = 1000
)


print(
  fdis_k_stable
)


# Flag only convincing basis-dimension problems:
#
#   k-index < 1
#   p < 0.05
#   EDF uses >= 80% of available basis dimension

k_valid_rows <-
  !is.na(
    fdis_k_stable[
      ,
      "k-index"
    ]
  )


k_followup_flag <- any(
  fdis_k_stable[
    k_valid_rows,
    "k-index"
  ] < 1 &
    
    fdis_k_stable[
      k_valid_rows,
      "p-value"
    ] < 0.05 &
    
    (
      fdis_k_stable[
        k_valid_rows,
        "edf"
      ] /
        fdis_k_stable[
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
# 29. Variance components
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FDIS VARIANCE COMPONENTS\n")
cat("============================================================\n\n")


fdis_vcomp <- gam.vcomp(
  fdis_model_final
)


print(
  fdis_vcomp
)


# ------------------------------------------------------------
# 30. Full concurvity
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FDIS FULL CONCURVITY\n")
cat("============================================================\n\n")


fdis_concurvity_full <- concurvity(
  fdis_model_final,
  full = TRUE
)


print(
  fdis_concurvity_full
)


# ------------------------------------------------------------
# 31. Pairwise concurvity
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FDIS PAIRWISE CONCURVITY\n")
cat("============================================================\n\n")


fdis_concurvity_pairwise <- concurvity(
  fdis_model_final,
  full = FALSE
)


fdis_pairwise_estimate <-
  fdis_concurvity_pairwise$estimate


print(
  round(
    fdis_pairwise_estimate,
    3
  )
)


# ------------------------------------------------------------
# 32. Ecological-term concurvity
# ------------------------------------------------------------

fdis_ecological_terms <- c(
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


fdis_ecological_concurvity <-
  fdis_pairwise_estimate[
    fdis_ecological_terms,
    fdis_ecological_terms,
    drop = FALSE
  ]


cat("\n============================================================\n")
cat("FDIS ECOLOGICAL-TERM CONCURVITY\n")
cat("============================================================\n\n")


print(
  round(
    fdis_ecological_concurvity,
    3
  )
)


# ------------------------------------------------------------
# 33. Residual diagnostics
# ------------------------------------------------------------

fdis_response_residuals <- residuals(
  fdis_model_final,
  type = "response"
)


fdis_sigma <- sqrt(
  fdis_summary$scale
)


fdis_standardized_residuals <-
  fdis_response_residuals /
  fdis_sigma


fdis_fitted_values <- fitted(
  fdis_model_final
)


cat("\n============================================================\n")
cat("FDIS STANDARDIZED RESIDUAL SUMMARY\n")
cat("============================================================\n\n")


print(
  summary(
    fdis_standardized_residuals
  )
)


cat(
  "\nAbsolute standardized residual > 2: ",
  sum(
    abs(
      fdis_standardized_residuals
    ) > 2
  ),
  " (",
  round(
    100 *
      mean(
        abs(
          fdis_standardized_residuals
        ) > 2
      ),
    2
  ),
  "%)\n",
  sep = ""
)


cat(
  "Absolute standardized residual > 3: ",
  sum(
    abs(
      fdis_standardized_residuals
    ) > 3
  ),
  " (",
  round(
    100 *
      mean(
        abs(
          fdis_standardized_residuals
        ) > 3
      ),
    2
  ),
  "%)\n",
  sep = ""
)


# ------------------------------------------------------------
# 34. Residual-vs-fitted diagnostic figure
# ------------------------------------------------------------

fdis_diagnostic_data <- data.frame(
  fitted =
    fdis_fitted_values,
  
  standardized_residual =
    fdis_standardized_residuals
)


fig_fdis_residuals <- ggplot(
  fdis_diagnostic_data,
  aes(
    x = fitted,
    y = standardized_residual
  )
) +
  geom_point(
    alpha = 0.25,
    size = 0.7
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "Fitted functional dispersion",
    y = "Standardized residual"
  ) +
  theme_classic(
    base_size = 12
  )


print(
  fig_fdis_residuals
)


# ------------------------------------------------------------
# 35. QQ diagnostic figure
# ------------------------------------------------------------

fig_fdis_qq <- ggplot(
  fdis_diagnostic_data,
  aes(
    sample = standardized_residual
  )
) +
  stat_qq(
    alpha = 0.30,
    size = 0.8
  ) +
  stat_qq_line() +
  labs(
    x = "Theoretical quantile",
    y = "Standardized residual quantile"
  ) +
  theme_classic(
    base_size = 12
  )


print(
  fig_fdis_qq
)


# ------------------------------------------------------------
# 36. Standard mgcv diagnostic
# ------------------------------------------------------------

set.seed(
  K_CHECK_SEED
)


gam.check(
  fdis_model_final
)


# ------------------------------------------------------------
# 37. Final parametric terms
# ------------------------------------------------------------

fdis_parametric_results <- as.data.frame(
  fdis_summary$p.table
)


fdis_parametric_results$term <- rownames(
  fdis_parametric_results
)


fdis_parametric_results <- fdis_parametric_results %>%
  select(
    term,
    everything()
  )


cat("\n============================================================\n")
cat("FINAL FDIS PARAMETRIC TERMS\n")
cat("============================================================\n\n")


print(
  fdis_parametric_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 38. Final smooth terms
# ------------------------------------------------------------

fdis_smooth_results <- as.data.frame(
  fdis_summary$s.table
)


fdis_smooth_results$term <- rownames(
  fdis_smooth_results
)


fdis_smooth_results <- fdis_smooth_results %>%
  select(
    term,
    everything()
  )


cat("\n============================================================\n")
cat("FINAL FDIS SMOOTH TERMS\n")
cat("============================================================\n\n")


print(
  fdis_smooth_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 39. Final model fit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL FDIS MODEL FIT\n")
cat("============================================================\n\n")


cat(
  "Adjusted R-squared: ",
  round(
    fdis_summary$r.sq,
    3
  ),
  "\n",
  sep = ""
)


cat(
  "Deviance explained: ",
  round(
    100 *
      fdis_summary$dev.expl,
    2
  ),
  "%\n",
  sep = ""
)


cat(
  "Residual SD: ",
  round(
    fdis_sigma,
    4
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 40. Fitted-value support
# ------------------------------------------------------------
#
# Gaussian models can theoretically predict negative FDis.
# Verify that fitted values remain biologically sensible.
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FDIS FITTED-VALUE QC\n")
cat("============================================================\n\n")


cat(
  "Minimum fitted FDis: ",
  min(
    fdis_fitted_values
  ),
  "\n",
  sep = ""
)


cat(
  "Maximum fitted FDis: ",
  max(
    fdis_fitted_values
  ),
  "\n",
  sep = ""
)


cat(
  "Negative fitted values: ",
  sum(
    fdis_fitted_values < 0
  ),
  "\n",
  sep = ""
)


# ============================================================
# PUBLICATION-STYLE PREDICTIONS
# ============================================================


# ============================================================
# PUBLICATION-STYLE PREDICTIONS
# ============================================================


# ------------------------------------------------------------
# 41. Reference values for population-level predictions
# ------------------------------------------------------------

reference_aridity <- median(
  fdis_data$avg_aridity_index,
  na.rm = TRUE
)


reference_duration <- median(
  fdis_data$duration_minutes,
  na.rm = TRUE
)


reference_distance <- median(
  fdis_data$effort_distance_km,
  na.rm = TRUE
)


reference_water <- median(
  fdis_data$water,
  na.rm = TRUE
)


reference_shrub <- median(
  fdis_data$shrub_and_scrub,
  na.rm = TRUE
)


reference_built <- median(
  fdis_data$built,
  na.rm = TRUE
)


reference_protocol <- factor(
  "Stationary",
  levels = levels(
    fdis_data$protocol_name
  )
)


reference_city <- factor(
  levels(
    fdis_data$NAME10
  )[1],
  levels = levels(
    fdis_data$NAME10
  )
)


reference_locality <- factor(
  levels(
    fdis_data$locality_id
  )[1],
  levels = levels(
    fdis_data$locality_id
  )
)


# ------------------------------------------------------------
# 42. Determine habitat scaling for figures
# ------------------------------------------------------------
#
# Convert habitat variables to percentages on the x-axis if
# they are stored as proportions.
# ------------------------------------------------------------

maximum_habitat_value <- max(
  c(
    fdis_data$water,
    fdis_data$shrub_and_scrub,
    fdis_data$built
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
# 43. Aridity prediction data
# ------------------------------------------------------------

aridity_limits <- quantile(
  fdis_data$avg_aridity_index,
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


aridity_prediction <- predict(
  fdis_model_final,
  newdata =
    aridity_prediction_data,
  type = "response",
  se.fit = TRUE,
  exclude = c(
    "s(NAME10)",
    "s(locality_id)"
  )
)


aridity_prediction_data$fit <-
  as.numeric(
    aridity_prediction$fit
  )


aridity_prediction_data$se <-
  as.numeric(
    aridity_prediction$se.fit
  )


aridity_prediction_data$lower <-
  aridity_prediction_data$fit -
  1.96 *
  aridity_prediction_data$se


aridity_prediction_data$upper <-
  aridity_prediction_data$fit +
  1.96 *
  aridity_prediction_data$se


aridity_plot_data <- data.frame(
  predictor =
    "Aridity index\n(higher = less arid)",
  
  x_value =
    aridity_prediction_data$avg_aridity_index,
  
  fit =
    aridity_prediction_data$fit,
  
  lower =
    aridity_prediction_data$lower,
  
  upper =
    aridity_prediction_data$upper
)


# ------------------------------------------------------------
# 44. Water prediction data
# ------------------------------------------------------------

water_limits <- quantile(
  fdis_data$water,
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


water_prediction <- predict(
  fdis_model_final,
  newdata =
    water_prediction_data,
  type = "response",
  se.fit = TRUE,
  exclude = c(
    "s(NAME10)",
    "s(locality_id)"
  )
)


water_prediction_data$fit <-
  as.numeric(
    water_prediction$fit
  )


water_prediction_data$se <-
  as.numeric(
    water_prediction$se.fit
  )


water_prediction_data$lower <-
  water_prediction_data$fit -
  1.96 *
  water_prediction_data$se


water_prediction_data$upper <-
  water_prediction_data$fit +
  1.96 *
  water_prediction_data$se


water_plot_data <- data.frame(
  predictor =
    "Water cover (%)",
  
  x_value =
    water_prediction_data$water *
    habitat_scale_factor,
  
  fit =
    water_prediction_data$fit,
  
  lower =
    water_prediction_data$lower,
  
  upper =
    water_prediction_data$upper
)


# ------------------------------------------------------------
# 45. Shrub/scrub prediction data
# ------------------------------------------------------------

shrub_limits <- quantile(
  fdis_data$shrub_and_scrub,
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


shrub_prediction <- predict(
  fdis_model_final,
  newdata =
    shrub_prediction_data,
  type = "response",
  se.fit = TRUE,
  exclude = c(
    "s(NAME10)",
    "s(locality_id)"
  )
)


shrub_prediction_data$fit <-
  as.numeric(
    shrub_prediction$fit
  )


shrub_prediction_data$se <-
  as.numeric(
    shrub_prediction$se.fit
  )


shrub_prediction_data$lower <-
  shrub_prediction_data$fit -
  1.96 *
  shrub_prediction_data$se


shrub_prediction_data$upper <-
  shrub_prediction_data$fit +
  1.96 *
  shrub_prediction_data$se


shrub_plot_data <- data.frame(
  predictor =
    "Shrub and scrub cover (%)",
  
  x_value =
    shrub_prediction_data$shrub_and_scrub *
    habitat_scale_factor,
  
  fit =
    shrub_prediction_data$fit,
  
  lower =
    shrub_prediction_data$lower,
  
  upper =
    shrub_prediction_data$upper
)


# ------------------------------------------------------------
# 46. Combine publication prediction data
# ------------------------------------------------------------

fdis_conditional_effects <- bind_rows(
  aridity_plot_data,
  water_plot_data,
  shrub_plot_data
) %>%
  mutate(
    predictor = factor(
      predictor,
      levels = c(
        "Aridity index\n(higher = less arid)",
        "Water cover (%)",
        "Shrub and scrub cover (%)"
      )
    )
  )


# ------------------------------------------------------------
# 47. Main publication figure
# ------------------------------------------------------------
#
# Each panel shows the conditional population-level relationship
# between FDis and one focal predictor.
#
# Other habitat and effort predictors are held at their median
# values, protocol is held at Stationary, and city/locality
# random effects are excluded.
# ------------------------------------------------------------

fig_fdis_main <- ggplot(
  fdis_conditional_effects,
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
    y = "Predicted functional dispersion"
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
  fig_fdis_main
)


# ------------------------------------------------------------
# 48. Automated QC summary
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


qc_positive_response <-
  all(
    fdis_data$FDis > 0
  )


qc_converged <-
  isTRUE(
    fdis_model_final$converged
  )


qc_full_rank <-
  fdis_model_final$rank ==
  length(
    coef(
      fdis_model_final
    )
  )


qc_basis_dimension <-
  !k_followup_flag


qc_fitted_support <-
  all(
    fdis_fitted_values >= 0
  )


fdis_qc_summary <- data.frame(
  check = c(
    "Balanced 5 x 5 sampling",
    "Multispecies FDis is strictly positive",
    "Final BAM converged",
    "Final BAM is full rank",
    "Basis dimensions acceptable",
    "Gaussian fitted values are non-negative"
  ),
  
  pass = c(
    qc_balanced_sample,
    qc_positive_response,
    qc_converged,
    qc_full_rank,
    qc_basis_dimension,
    qc_fitted_support
  )
)


cat("\n============================================================\n")
cat("FDIS AUTOMATED QC SUMMARY\n")
cat("============================================================\n\n")


print(
  fdis_qc_summary,
  row.names = FALSE
)


cat(
  "\nChecks passed: ",
  sum(
    fdis_qc_summary$pass
  ),
  " of ",
  nrow(
    fdis_qc_summary
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 49. Final analysis summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FUNCTIONAL DISPERSION ANALYSIS COMPLETE\n")
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
    fdis_data$NAME10
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled localities:             ",
  n_distinct(
    fdis_data$locality_id
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled checklists:             ",
  nrow(
    fdis_data
  ),
  "\n",
  sep = ""
)


cat(
  "Additive AIC:                   ",
  round(
    AIC(
      fdis_model_additive
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
      fdis_model_interaction
    ),
    2
  ),
  "\n",
  sep = ""
)


cat(
  "Final adjusted R-squared:       ",
  round(
    fdis_summary$r.sq,
    3
  ),
  "\n",
  sep = ""
)


cat(
  "Final deviance explained:       ",
  round(
    100 *
      fdis_summary$dev.expl,
    2
  ),
  "%\n",
  sep = ""
)


cat(
  "Basis-dimension follow-up:      ",
  k_followup_flag,
  "\n",
  sep = ""
)


cat(
  "Negative fitted values:         ",
  sum(
    fdis_fitted_values < 0
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 50. Optional figure saving
# ------------------------------------------------------------

if (
  SAVE_FIGURES
) {
  
  figure_directory <- file.path(
    "Outputs",
    "04_Functional_Dispersion",
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
        "FDis_Distribution.pdf"
      ),
    plot =
      fig_fdis_distribution,
    width = 5.5,
    height = 4
  )
  
  
  ggsave(
    filename =
      file.path(
        figure_directory,
        "FDis_Conditional_Effects.pdf"
      ),
    plot =
      fig_fdis_main,
    width = 9,
    height = 3.8
  )
}


# ------------------------------------------------------------
# 51. Overall interpretation
# ------------------------------------------------------------
#
# The final FDis analysis is restricted to multispecies
# checklists (sr >= 2) because single-species assemblages have
# FDis = 0 by definition and therefore represent structural
# zeros rather than variation in among-species functional
# dispersion.
#
# A Gaussian BAM was selected for multispecies FDis after
# comparison with a Gamma BAM. The Gaussian model provided
# substantially better overall fit, explained more deviance,
# showed acceptable residual behavior, and produced no
# biologically impossible negative fitted values.
#
# The primary FDis model retains the global aridity x habitat
# interaction structure used across biodiversity responses.
# The additive model had slightly lower AIC, but the difference
# was small and AIC comparisons among fREML smooth models are
# treated descriptively rather than as the sole basis for
# model selection. The interaction model is therefore retained
# because the interactions represent the a priori ecological
# hypotheses.
#
# Final model results indicate:
#
# - Functional dispersion increases with increasing aridity
#   index, corresponding to greater FDis under relatively less
#   arid conditions.
#
# - Water cover is strongly associated with functional
#   dispersion, with model-based predictions indicating
#   increasing FDis across the observed water-cover gradient.
#
# - Shrub/scrub cover is also strongly associated with FDis,
#   with predicted functional dispersion generally increasing
#   across the observed shrub/scrub gradient.
#
# - Built cover shows little evidence of an independent
#   relationship with FDis.
#
# - None of the three aridity x habitat interactions show
#   strong statistical support. The aridity x water interaction
#   is comparatively weak/marginal, while the aridity x shrub
#   and aridity x built interactions are unsupported.
#
# - These results therefore provide stronger evidence for
#   additive effects of regional aridity and local habitat than
#   for strong modification of habitat-FDis relationships
#   across the aridity gradient.
#
# - Checklist duration remains an important sampling-effort
#   covariate. Distance traveled shows comparatively weak
#   evidence after accounting for the other model terms.
#
# - City and locality random-effect smooths are retained to
#   account for the hierarchical structure of the eBird data.
#
# - The final Gaussian model explains approximately 44.6% of
#   deviance, with an adjusted R-squared of approximately 0.36.
#
# - Residual diagnostics show mild heteroscedasticity and some
#   departure from normality in the tails, but no major pattern
#   indicating model failure. Basis dimensions are adequate,
#   the model is full rank and converged, and fitted FDis values
#   remain non-negative.
#
# Overall, functional dispersion among multispecies bird
# assemblages appears to respond to both the regional aridity
# gradient and local habitat structure. However, unlike the
# species-richness response, there is comparatively little
# evidence that the effects of local habitat on functional
# dispersion change strongly across the aridity gradient.
#
# ============================================================