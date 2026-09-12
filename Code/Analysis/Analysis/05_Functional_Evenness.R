# ============================================================
# Dryland Bird Project
# 05_Functional_Evenness.R
#
# Purpose:
# Model checklist-level functional evenness (FEve) as a
# function of aridity, local habitat, and checklist effort
# using a beta GAM/BAM framework and balanced hierarchical
# sampling.
#
# Sampling design:
#
#   1. Retain Stationary and Traveling checklists
#   2. Stationary distance set to zero
#   3. Duration <= 300 minutes
#   4. Distance <= 5 km
#   5. Retain checklists with non-missing FEve
#   6. Require >= 5 FEve checklists per locality
#   7. Require >= 5 eligible localities per city
#   8. Randomly sample 5 localities per city
#   9. Randomly sample 5 checklists per locality
#
# Candidate models:
#
#   A. Additive beta BAM
#   B. Aridity x habitat interaction beta BAM
#
# Primary model:
#
# The interaction model is retained as the primary/global model
# because aridity x habitat relationships represent the a priori
# ecological hypotheses and allow direct comparison among
# biodiversity responses.
#
# Figures:
#
#   - FEve response distribution
#
#   Observed/descriptive:
#     - FEve vs aridity
#     - FEve vs water cover
#     - FEve vs shrub/scrub cover
#     - FEve vs built cover
#
#   Adjusted/model-based:
#     - Predicted FEve vs aridity
#     - Predicted FEve vs water cover
#     - Predicted FEve vs shrub/scrub cover
#     - Predicted FEve vs built cover
#
# Observed and predicted relationships are intentionally shown
# in separate figures.
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
cat("FUNCTIONAL EVENNESS ANALYSIS\n")
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
  "FEve",
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
  "PASS: All required FEve analysis variables are present.\n"
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

feve_pool <- analysis_base %>%
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


feve_pool %>%
  summarise(
    n_checklists = n(),
    n_localities = n_distinct(locality_id),
    n_cities = n_distinct(NAME10)
  ) %>%
  print()


# ------------------------------------------------------------
# 7. Initial FEve QC
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FEVE INITIAL QC\n")
cat("============================================================\n\n")


cat(
  "FEve missing: ",
  sum(
    is.na(
      feve_pool$FEve
    )
  ),
  "\n",
  sep = ""
)


cat(
  "FEve missing (%): ",
  round(
    100 *
      mean(
        is.na(
          feve_pool$FEve
        )
      ),
    2
  ),
  "%\n",
  sep = ""
)


cat(
  "FEve zeros: ",
  sum(
    feve_pool$FEve == 0,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "FEve ones: ",
  sum(
    feve_pool$FEve == 1,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "FEve < 0: ",
  sum(
    feve_pool$FEve < 0,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "FEve > 1: ",
  sum(
    feve_pool$FEve > 1,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


print(
  summary(
    feve_pool$FEve
  )
)


# ------------------------------------------------------------
# 8. FEve and species-richness structure
# ------------------------------------------------------------
#
# FEve is undefined for some low-richness assemblages.
#
# Examine the richness structure of defined versus undefined
# FEve observations before response-specific sampling.
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FEVE AND SPECIES-RICHNESS STRUCTURE\n")
cat("============================================================\n\n")


cat(
  "Missing FEve with sr = 1: ",
  sum(
    is.na(feve_pool$FEve) &
      feve_pool$sr == 1,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "Missing FEve with sr >= 2: ",
  sum(
    is.na(feve_pool$FEve) &
      feve_pool$sr >= 2,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "Defined FEve with sr = 1: ",
  sum(
    !is.na(feve_pool$FEve) &
      feve_pool$sr == 1,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "Defined FEve with sr >= 2: ",
  sum(
    !is.na(feve_pool$FEve) &
      feve_pool$sr >= 2,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Remove undefined FEve values
# ------------------------------------------------------------

feve_pool <- feve_pool %>%
  filter(
    !is.na(FEve)
  )


cat(
  "\nRows with defined FEve: ",
  nrow(feve_pool),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Complete-case restriction
# ------------------------------------------------------------

model_variables <- c(
  "FEve",
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
      feve_pool[
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


feve_pool <- feve_pool %>%
  filter(
    if_all(
      all_of(
        model_variables
      ),
      ~ !is.na(.)
    )
  )


# ------------------------------------------------------------
# 11. Identify eligible localities
# ------------------------------------------------------------

locality_sampling <- feve_pool %>%
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
# 12. Identify eligible cities
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
cat("FEVE ELIGIBILITY\n")
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
# 13. Create eligible sampling pool
# ------------------------------------------------------------

eligible_pool <- feve_pool %>%
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
# 14. Sample five localities per city
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
# 15. Sample five FEve checklists per locality
# ------------------------------------------------------------

set.seed(
  CHECKLIST_SAMPLE_SEED
)


feve_data <- eligible_pool %>%
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
# 16. Balanced sample QC
# ------------------------------------------------------------

locality_qc <- feve_data %>%
  count(
    NAME10,
    locality_id,
    name = "n_checklists"
  )


city_qc <- feve_data %>%
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
cat("FEVE BALANCED SAMPLE QC\n")
cat("============================================================\n\n")


cat(
  "Sampled cities:     ",
  n_distinct(
    feve_data$NAME10
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled localities: ",
  n_distinct(
    feve_data$locality_id
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled checklists: ",
  nrow(
    feve_data
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
    "FEve balanced sampling QC failed."
  )
}


cat(
  "\nPASS: FEve balanced sampling verified.\n"
)


# ------------------------------------------------------------
# 17. Prepare factors
# ------------------------------------------------------------

feve_data <- feve_data %>%
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
# 18. Final FEve response distribution
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL FEVE DISTRIBUTION\n")
cat("============================================================\n\n")


print(
  summary(
    feve_data$FEve
  )
)


cat(
  "\nZeros: ",
  sum(
    feve_data$FEve == 0
  ),
  "\n",
  sep = ""
)


cat(
  "Ones: ",
  sum(
    feve_data$FEve == 1
  ),
  "\n",
  sep = ""
)


cat(
  "Values < 0: ",
  sum(
    feve_data$FEve < 0
  ),
  "\n",
  sep = ""
)


cat(
  "Values > 1: ",
  sum(
    feve_data$FEve > 1
  ),
  "\n",
  sep = ""
)


cat(
  "Minimum: ",
  min(
    feve_data$FEve
  ),
  "\n",
  sep = ""
)


cat(
  "Maximum: ",
  max(
    feve_data$FEve
  ),
  "\n",
  sep = ""
)


cat(
  "Mean: ",
  mean(
    feve_data$FEve
  ),
  "\n",
  sep = ""
)


cat(
  "Median: ",
  median(
    feve_data$FEve
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 19. Beta-family gate
# ------------------------------------------------------------

n_boundary <- sum(
  feve_data$FEve <= 0 |
    feve_data$FEve >= 1
)


cat("\n============================================================\n")
cat("FEVE BETA FAMILY QC\n")
cat("============================================================\n\n")


cat(
  "FEve observations outside (0,1): ",
  n_boundary,
  "\n",
  sep = ""
)


if (
  n_boundary > 0
) {
  
  stop(
    paste0(
      "FEve contains exact 0/1 or out-of-range observations. ",
      "Review boundary values before fitting the beta BAM."
    )
  )
}


cat(
  "PASS: FEve is strictly inside (0,1).\n"
)


# ============================================================
# CANDIDATE BETA BAM MODELS
# ============================================================


# ------------------------------------------------------------
# 20. Additive model formula
# ------------------------------------------------------------

feve_formula_additive <- FEve ~
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
# 21. Interaction/global model formula
# ------------------------------------------------------------

feve_formula_interaction <- FEve ~
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
# 22. Fit additive beta BAM
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FITTING ADDITIVE FEVE BETA BAM\n")
cat("============================================================\n\n")


additive_fit_time <- system.time(
  
  feve_model_additive <- bam(
    formula =
      feve_formula_additive,
    
    family =
      betar(
        link = "logit"
      ),
    
    data =
      feve_data,
    
    method =
      "fREML",
    
    discrete =
      TRUE
  )
  
)


cat(
  "\nAdditive BAM fit time:\n"
)


print(
  additive_fit_time
)


# ------------------------------------------------------------
# 23. Fit interaction beta BAM
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FITTING INTERACTION FEVE BETA BAM\n")
cat("============================================================\n\n")


interaction_fit_time <- system.time(
  
  feve_model_interaction <- bam(
    formula =
      feve_formula_interaction,
    
    family =
      betar(
        link = "logit"
      ),
    
    data =
      feve_data,
    
    method =
      "fREML",
    
    discrete =
      TRUE
  )
  
)


cat(
  "\nInteraction BAM fit time:\n"
)


print(
  interaction_fit_time
)


# ------------------------------------------------------------
# 24. Candidate-model summaries
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ADDITIVE MODEL SUMMARY\n")
cat("============================================================\n\n")


print(
  summary(
    feve_model_additive
  )
)


cat("\n============================================================\n")
cat("INTERACTION MODEL SUMMARY\n")
cat("============================================================\n\n")


print(
  summary(
    feve_model_interaction
  )
)


# ------------------------------------------------------------
# 25. Candidate-model comparison
# ------------------------------------------------------------
#
# AIC is treated as descriptive because these models were fit
# with fREML and differ in smooth structure.
# ------------------------------------------------------------

additive_summary <- summary(
  feve_model_additive
)


interaction_summary <- summary(
  feve_model_interaction
)


feve_model_comparison <- data.frame(
  model = c(
    "Additive",
    "Aridity x habitat interactions"
  ),
  
  AIC = c(
    AIC(
      feve_model_additive
    ),
    
    AIC(
      feve_model_interaction
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


feve_model_comparison <- feve_model_comparison %>%
  mutate(
    delta_AIC =
      AIC -
      min(
        AIC
      )
  )


cat("\n============================================================\n")
cat("FEVE MODEL COMPARISON\n")
cat("============================================================\n\n")


print(
  feve_model_comparison
)


# ------------------------------------------------------------
# 26. Candidate-model convergence
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CANDIDATE MODEL CONVERGENCE\n")
cat("============================================================\n\n")


cat(
  "Additive model converged: ",
  feve_model_additive$converged,
  "\n",
  sep = ""
)


cat(
  "Interaction model converged: ",
  feve_model_interaction$converged,
  "\n",
  sep = ""
)


cat(
  "Additive model rank: ",
  feve_model_additive$rank,
  " / ",
  length(
    coef(
      feve_model_additive
    )
  ),
  "\n",
  sep = ""
)


cat(
  "Interaction model rank: ",
  feve_model_interaction$rank,
  " / ",
  length(
    coef(
      feve_model_interaction
    )
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 27. Preliminary basis checks
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ADDITIVE MODEL K CHECK\n")
cat("============================================================\n\n")


set.seed(
  K_CHECK_SEED
)


print(
  k.check(
    feve_model_additive
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
    feve_model_interaction
  )
)


# ============================================================
# FINAL FEVE MODEL
# ============================================================


# ------------------------------------------------------------
# 28. Select primary/global FEve model
# ------------------------------------------------------------
#
# The interaction model is retained as the primary/global FEve
# model.
#
# The additive model has slightly lower AIC, but the difference
# is small and AIC comparisons among models fit with fREML and
# differing smooth structures are treated descriptively.
#
# Aridity x habitat interactions represent the project's
# a priori ecological hypotheses. Maintaining the same global
# structure also facilitates comparison among species richness,
# FRic, FDis, and FEve.
# ------------------------------------------------------------

feve_model_final <- feve_model_interaction


feve_summary <- summary(
  feve_model_final
)


cat("\n============================================================\n")
cat("FINAL FEVE MODEL\n")
cat("============================================================\n\n")


print(
  feve_summary
)


# ============================================================
# FINAL FEVE QC AND DIAGNOSTICS
# ============================================================


# ------------------------------------------------------------
# 29. Final convergence and rank
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FEVE MODEL CONVERGENCE\n")
cat("============================================================\n\n")


cat(
  "Model converged: ",
  feve_model_final$converged,
  "\n",
  sep = ""
)


cat(
  "Model rank: ",
  feve_model_final$rank,
  " / ",
  length(
    coef(
      feve_model_final
    )
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 30. Stable basis-dimension check
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FEVE STABLE K CHECK - 1000 REPLICATES\n")
cat("============================================================\n\n")


set.seed(
  K_CHECK_SEED
)


feve_k_stable <- k.check(
  feve_model_final,
  n.rep = 1000
)


print(
  feve_k_stable
)


k_valid_rows <-
  !is.na(
    feve_k_stable[
      ,
      "k-index"
    ]
  )


k_followup_flag <- any(
  feve_k_stable[
    k_valid_rows,
    "k-index"
  ] < 1 &
    
    feve_k_stable[
      k_valid_rows,
      "p-value"
    ] < 0.05 &
    
    (
      feve_k_stable[
        k_valid_rows,
        "edf"
      ] /
        feve_k_stable[
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
# 31. Variance components
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FEVE VARIANCE COMPONENTS\n")
cat("============================================================\n\n")


feve_vcomp <- gam.vcomp(
  feve_model_final
)


print(
  feve_vcomp
)


# ------------------------------------------------------------
# 32. Full concurvity
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FEVE FULL CONCURVITY\n")
cat("============================================================\n\n")


feve_concurvity_full <- concurvity(
  feve_model_final,
  full = TRUE
)


print(
  feve_concurvity_full
)


# ------------------------------------------------------------
# 33. Pairwise concurvity
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FEVE PAIRWISE CONCURVITY\n")
cat("============================================================\n\n")


feve_concurvity_pairwise <- concurvity(
  feve_model_final,
  full = FALSE
)


feve_pairwise_estimate <-
  feve_concurvity_pairwise$estimate


print(
  round(
    feve_pairwise_estimate,
    3
  )
)


# ------------------------------------------------------------
# 34. Ecological-term concurvity
# ------------------------------------------------------------

feve_ecological_terms <- c(
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


feve_ecological_concurvity <-
  feve_pairwise_estimate[
    feve_ecological_terms,
    feve_ecological_terms,
    drop = FALSE
  ]


cat("\n============================================================\n")
cat("FEVE ECOLOGICAL-TERM CONCURVITY\n")
cat("============================================================\n\n")


print(
  round(
    feve_ecological_concurvity,
    3
  )
)


# ------------------------------------------------------------
# 35. Residual summaries
# ------------------------------------------------------------

feve_pearson_residuals <- residuals(
  feve_model_final,
  type = "pearson"
)


feve_deviance_residuals <- residuals(
  feve_model_final,
  type = "deviance"
)


feve_fitted_values <- fitted(
  feve_model_final
)


cat("\n============================================================\n")
cat("FEVE RESIDUAL SUMMARY\n")
cat("============================================================\n\n")


cat(
  "Pearson residuals:\n"
)


print(
  summary(
    feve_pearson_residuals
  )
)


cat(
  "\nDeviance residuals:\n"
)


print(
  summary(
    feve_deviance_residuals
  )
)


# ------------------------------------------------------------
# 36. Extreme Pearson residuals
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FEVE EXTREME RESIDUALS\n")
cat("============================================================\n\n")


cat(
  "Absolute Pearson residual > 2: ",
  sum(
    abs(
      feve_pearson_residuals
    ) > 2
  ),
  " (",
  round(
    100 *
      mean(
        abs(
          feve_pearson_residuals
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
      feve_pearson_residuals
    ) > 3
  ),
  " (",
  round(
    100 *
      mean(
        abs(
          feve_pearson_residuals
        ) > 3
      ),
    2
  ),
  "%)\n",
  sep = ""
)


# ------------------------------------------------------------
# 37. Residual-vs-fitted diagnostic figure
# ------------------------------------------------------------

feve_diagnostic_data <- data.frame(
  fitted =
    feve_fitted_values,
  
  pearson_residual =
    feve_pearson_residuals
)


fig_feve_residuals <- ggplot(
  feve_diagnostic_data,
  aes(
    x = fitted,
    y = pearson_residual
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
    x = "Fitted functional evenness",
    y = "Pearson residual"
  ) +
  theme_classic(
    base_size = 12
  )


print(
  fig_feve_residuals
)


# ------------------------------------------------------------
# 38. QQ diagnostic
# ------------------------------------------------------------

set.seed(
  K_CHECK_SEED
)


qq.gam(
  feve_model_final,
  rep = 100,
  level = 0.90
)


# ------------------------------------------------------------
# 39. Standard GAM diagnostics
# ------------------------------------------------------------

set.seed(
  K_CHECK_SEED
)


gam.check(
  feve_model_final
)


# ------------------------------------------------------------
# 40. Final parametric-term table
# ------------------------------------------------------------

feve_parametric_results <- as.data.frame(
  feve_summary$p.table
)


feve_parametric_results$term <- rownames(
  feve_parametric_results
)


feve_parametric_results <- feve_parametric_results %>%
  select(
    term,
    everything()
  )


cat("\n============================================================\n")
cat("FINAL FEVE PARAMETRIC TERMS\n")
cat("============================================================\n\n")


print(
  feve_parametric_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 41. Final smooth-term table
# ------------------------------------------------------------

feve_smooth_results <- as.data.frame(
  feve_summary$s.table
)


feve_smooth_results$term <- rownames(
  feve_smooth_results
)


feve_smooth_results <- feve_smooth_results %>%
  select(
    term,
    everything()
  )


cat("\n============================================================\n")
cat("FINAL FEVE SMOOTH TERMS\n")
cat("============================================================\n\n")


print(
  feve_smooth_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 42. Final model fit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL FEVE MODEL FIT\n")
cat("============================================================\n\n")


cat(
  "Adjusted R-squared: ",
  round(
    feve_summary$r.sq,
    3
  ),
  "\n",
  sep = ""
)


cat(
  "Deviance explained: ",
  round(
    100 *
      feve_summary$dev.expl,
    2
  ),
  "%\n",
  sep = ""
)


cat(
  "Beta precision parameter: ",
  round(
    feve_model_final$family$getTheta(
      TRUE
    ),
    3
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 43. Diagnostic smooth plots
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
  feve_model_final,
  select = 1,
  shade = TRUE,
  pages = 1
)


plot(
  feve_model_final,
  select = 2,
  shade = TRUE,
  pages = 1
)


plot(
  feve_model_final,
  select = 3,
  shade = TRUE,
  pages = 1
)


plot(
  feve_model_final,
  select = 4,
  shade = TRUE,
  pages = 1
)


plot(
  feve_model_final,
  select = 5,
  pages = 1
)


plot(
  feve_model_final,
  select = 6,
  pages = 1
)


plot(
  feve_model_final,
  select = 7,
  pages = 1
)


plot(
  feve_model_final,
  select = 8,
  shade = TRUE,
  pages = 1
)


plot(
  feve_model_final,
  select = 9,
  shade = TRUE,
  pages = 1
)


# ============================================================
# PUBLICATION-STYLE FIGURES
# ============================================================
#
# OBSERVED figures:
#
#   Actual checklist-level FEve values plus a descriptive,
#   unadjusted single-predictor beta GAM.
#
# PREDICTED figures:
#
#   Population-level predictions from the final multivariable
#   beta BAM, with other predictors held at representative
#   values and city/locality random effects excluded.
#
# Observed and predicted figures are intentionally kept
# completely separate.
#
# ============================================================


# ------------------------------------------------------------
# 44. FEve response-distribution figure
# ------------------------------------------------------------

fig_feve_distribution <- ggplot(
  feve_data,
  aes(
    x = FEve
  )
) +
  geom_histogram(
    bins = 40,
    linewidth = 0.25
  ) +
  labs(
    x = "Functional evenness",
    y = "Number of checklists"
  ) +
  theme_classic(
    base_size = 12
  )


print(
  fig_feve_distribution
)


# ------------------------------------------------------------
# 45. Determine habitat scaling
# ------------------------------------------------------------

maximum_habitat_value <- max(
  c(
    feve_data$water,
    feve_data$shrub_and_scrub,
    feve_data$built
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


# ============================================================
# OBSERVED RELATIONSHIPS
# ============================================================


# ------------------------------------------------------------
# 46. Helper for observed relationships
# ------------------------------------------------------------
#
# Each observed figure contains the actual checklist-level FEve
# observations and an unadjusted single-predictor beta GAM.
#
# The descriptive smooth is shown only across the central 90%
# of the focal predictor distribution to reduce visual leverage
# from extreme predictor values.
#
# These figures describe the raw empirical relationships.
# They are NOT adjusted effects from the final multivariable
# model.
# ------------------------------------------------------------

make_feve_observed <- function(
    variable,
    x_label,
    scale_factor = 1
) {
  
  observed_data <- data.frame(
    x =
      feve_data[[variable]] *
      scale_factor,
    
    FEve =
      feve_data$FEve
  )
  
  
  raw_formula <- as.formula(
    paste0(
      "FEve ~ s(",
      variable,
      ", k = 5)"
    )
  )
  
  
  raw_model <- gam(
    formula =
      raw_formula,
    
    family =
      betar(
        link = "logit"
      ),
    
    data =
      feve_data,
    
    method =
      "REML"
  )
  
  
  x_limits <- quantile(
    feve_data[[variable]],
    probs = c(
      0.05,
      0.95
    ),
    na.rm = TRUE
  )
  
  
  x_sequence <- seq(
    x_limits[1],
    x_limits[2],
    length.out = 200
  )
  
  
  prediction_data <- data.frame(
    x_sequence
  )
  
  
  names(
    prediction_data
  ) <- variable
  
  
  raw_prediction <- predict(
    raw_model,
    newdata =
      prediction_data,
    type = "link",
    se.fit = TRUE
  )
  
  
  fit_link <- as.numeric(
    raw_prediction$fit
  )
  
  
  se_link <- as.numeric(
    raw_prediction$se.fit
  )
  
  
  link_inverse <-
    raw_model$family$linkinv
  
  
  trend_data <- data.frame(
    x =
      x_sequence *
      scale_factor,
    
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
  
  
  ggplot() +
    geom_point(
      data =
        observed_data,
      aes(
        x = x,
        y = FEve
      ),
      alpha = 0.08,
      size = 0.55
    ) +
    geom_ribbon(
      data =
        trend_data,
      aes(
        x = x,
        ymin = lower,
        ymax = upper
      ),
      fill = "grey80",
      alpha = 0.7
    ) +
    geom_line(
      data =
        trend_data,
      aes(
        x = x,
        y = fit
      ),
      linewidth = 0.9
    ) +
    coord_cartesian(
      xlim =
        x_limits *
        scale_factor
    ) +
    labs(
      x = x_label,
      y = "Functional evenness"
    ) +
    theme_classic(
      base_size = 12
    )
}


# ------------------------------------------------------------
# 47. Observed FEve vs aridity
# ------------------------------------------------------------

fig_feve_aridity_observed <- make_feve_observed(
  variable =
    "avg_aridity_index",
  
  x_label =
    "Aridity index (higher = less arid)"
)


print(
  fig_feve_aridity_observed
)


# ------------------------------------------------------------
# 48. Observed FEve vs water cover
# ------------------------------------------------------------

fig_feve_water_observed <- make_feve_observed(
  variable =
    "water",
  
  x_label =
    "Water cover (%)",
  
  scale_factor =
    habitat_scale_factor
)


print(
  fig_feve_water_observed
)


# ------------------------------------------------------------
# 49. Observed FEve vs shrub/scrub cover
# ------------------------------------------------------------

fig_feve_shrub_observed <- make_feve_observed(
  variable =
    "shrub_and_scrub",
  
  x_label =
    "Shrub and scrub cover (%)",
  
  scale_factor =
    habitat_scale_factor
)


print(
  fig_feve_shrub_observed
)


# ------------------------------------------------------------
# 50. Observed FEve vs built cover
# ------------------------------------------------------------

fig_feve_built_observed <- make_feve_observed(
  variable =
    "built",
  
  x_label =
    "Built cover (%)",
  
  scale_factor =
    habitat_scale_factor
)


print(
  fig_feve_built_observed
)


# ============================================================
# MODEL-ADJUSTED RELATIONSHIPS
# ============================================================


# ------------------------------------------------------------
# 51. Reference values for adjusted predictions
# ------------------------------------------------------------

reference_aridity <- median(
  feve_data$avg_aridity_index,
  na.rm = TRUE
)


reference_duration <- median(
  feve_data$duration_minutes,
  na.rm = TRUE
)


reference_distance <- median(
  feve_data$effort_distance_km,
  na.rm = TRUE
)


reference_water <- median(
  feve_data$water,
  na.rm = TRUE
)


reference_shrub <- median(
  feve_data$shrub_and_scrub,
  na.rm = TRUE
)


reference_built <- median(
  feve_data$built,
  na.rm = TRUE
)


reference_protocol <- factor(
  "Stationary",
  levels =
    levels(
      feve_data$protocol_name
    )
)


reference_city <- factor(
  levels(
    feve_data$NAME10
  )[1],
  levels =
    levels(
      feve_data$NAME10
    )
)


reference_locality <- factor(
  levels(
    feve_data$locality_id
  )[1],
  levels =
    levels(
      feve_data$locality_id
    )
)


# ------------------------------------------------------------
# 52. Helper for adjusted FEve predictions
# ------------------------------------------------------------
#
# Prediction intervals are calculated on the beta-model logit
# scale and transformed back to the original FEve scale.
#
# Other continuous predictors are held at their median values.
# Protocol is held at Stationary.
#
# City and locality random-effect smooths are excluded to
# produce population-level predictions.
# ------------------------------------------------------------

make_feve_prediction <- function(
    variable,
    x_label,
    scale_factor = 1
) {
  
  x_limits <- quantile(
    feve_data[[variable]],
    probs = c(
      0.05,
      0.95
    ),
    na.rm = TRUE
  )
  
  
  x_sequence <- seq(
    x_limits[1],
    x_limits[2],
    length.out = 200
  )
  
  
  n_prediction <- length(
    x_sequence
  )
  
  
  prediction_data <- data.frame(
    avg_aridity_index =
      rep(
        reference_aridity,
        n_prediction
      ),
    
    water =
      rep(
        reference_water,
        n_prediction
      ),
    
    shrub_and_scrub =
      rep(
        reference_shrub,
        n_prediction
      ),
    
    built =
      rep(
        reference_built,
        n_prediction
      ),
    
    duration_minutes =
      rep(
        reference_duration,
        n_prediction
      ),
    
    effort_distance_km =
      rep(
        reference_distance,
        n_prediction
      )
  )
  
  
  prediction_data$protocol_name <- factor(
    rep(
      "Stationary",
      n_prediction
    ),
    levels =
      levels(
        feve_data$protocol_name
      )
  )
  
  
  prediction_data$NAME10 <- factor(
    rep(
      levels(
        feve_data$NAME10
      )[1],
      n_prediction
    ),
    levels =
      levels(
        feve_data$NAME10
      )
  )
  
  
  prediction_data$locality_id <- factor(
    rep(
      levels(
        feve_data$locality_id
      )[1],
      n_prediction
    ),
    levels =
      levels(
        feve_data$locality_id
      )
  )
  
  
  prediction_data[[variable]] <-
    x_sequence
  
  
  model_prediction <- predict(
    feve_model_final,
    newdata =
      prediction_data,
    type = "link",
    se.fit = TRUE,
    exclude = c(
      "s(NAME10)",
      "s(locality_id)"
    )
  )
  
  
  fit_link <- as.numeric(
    model_prediction$fit
  )
  
  
  se_link <- as.numeric(
    model_prediction$se.fit
  )
  
  
  link_inverse <-
    feve_model_final$family$linkinv
  
  
  prediction_plot_data <- data.frame(
    x =
      x_sequence *
      scale_factor,
    
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
  
  
  ggplot(
    prediction_plot_data,
    aes(
      x = x,
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
    labs(
      x = x_label,
      y = "Predicted functional evenness"
    ) +
    theme_classic(
      base_size = 12
    )
}


# ------------------------------------------------------------
# 53. Predicted FEve vs aridity
# ------------------------------------------------------------

fig_feve_aridity_predicted <- make_feve_prediction(
  variable =
    "avg_aridity_index",
  
  x_label =
    "Aridity index (higher = less arid)"
)


print(
  fig_feve_aridity_predicted
)


# ------------------------------------------------------------
# 54. Predicted FEve vs water cover
# ------------------------------------------------------------

fig_feve_water_predicted <- make_feve_prediction(
  variable =
    "water",
  
  x_label =
    "Water cover (%)",
  
  scale_factor =
    habitat_scale_factor
)


print(
  fig_feve_water_predicted
)


# ------------------------------------------------------------
# 55. Predicted FEve vs shrub/scrub cover
# ------------------------------------------------------------

fig_feve_shrub_predicted <- make_feve_prediction(
  variable =
    "shrub_and_scrub",
  
  x_label =
    "Shrub and scrub cover (%)",
  
  scale_factor =
    habitat_scale_factor
)


print(
  fig_feve_shrub_predicted
)


# ------------------------------------------------------------
# 56. Predicted FEve vs built cover
# ------------------------------------------------------------

fig_feve_built_predicted <- make_feve_prediction(
  variable =
    "built",
  
  x_label =
    "Built cover (%)",
  
  scale_factor =
    habitat_scale_factor
)


print(
  fig_feve_built_predicted
)


# ------------------------------------------------------------
# 57. Optional figure saving
# ------------------------------------------------------------

if (
  SAVE_FIGURES
) {
  
  figure_directory <- file.path(
    "Outputs",
    "05_Functional_Evenness",
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
        "FEve_Distribution.pdf"
      ),
    plot =
      fig_feve_distribution,
    width = 5.5,
    height = 4
  )
  
  
  ggsave(
    filename =
      file.path(
        figure_directory,
        "FEve_Aridity_Observed.pdf"
      ),
    plot =
      fig_feve_aridity_observed,
    width = 5.5,
    height = 4
  )
  
  
  ggsave(
    filename =
      file.path(
        figure_directory,
        "FEve_Aridity_Predicted.pdf"
      ),
    plot =
      fig_feve_aridity_predicted,
    width = 5.5,
    height = 4
  )
  
  
  ggsave(
    filename =
      file.path(
        figure_directory,
        "FEve_Water_Observed.pdf"
      ),
    plot =
      fig_feve_water_observed,
    width = 5.5,
    height = 4
  )
  
  
  ggsave(
    filename =
      file.path(
        figure_directory,
        "FEve_Water_Predicted.pdf"
      ),
    plot =
      fig_feve_water_predicted,
    width = 5.5,
    height = 4
  )
  
  
  ggsave(
    filename =
      file.path(
        figure_directory,
        "FEve_Shrub_Observed.pdf"
      ),
    plot =
      fig_feve_shrub_observed,
    width = 5.5,
    height = 4
  )
  
  
  ggsave(
    filename =
      file.path(
        figure_directory,
        "FEve_Shrub_Predicted.pdf"
      ),
    plot =
      fig_feve_shrub_predicted,
    width = 5.5,
    height = 4
  )
  
  
  ggsave(
    filename =
      file.path(
        figure_directory,
        "FEve_Built_Observed.pdf"
      ),
    plot =
      fig_feve_built_observed,
    width = 5.5,
    height = 4
  )
  
  
  ggsave(
    filename =
      file.path(
        figure_directory,
        "FEve_Built_Predicted.pdf"
      ),
    plot =
      fig_feve_built_predicted,
    width = 5.5,
    height = 4
  )
}


# ============================================================
# FINAL QC AND ANALYSIS SUMMARY
# ============================================================


# ------------------------------------------------------------
# 58. Automated QC summary
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
    feve_data$FEve > 0 &
      feve_data$FEve < 1
  )


qc_converged <-
  isTRUE(
    feve_model_final$converged
  )


qc_full_rank <-
  feve_model_final$rank ==
  length(
    coef(
      feve_model_final
    )
  )


qc_basis_dimension <-
  !k_followup_flag


feve_qc_summary <- data.frame(
  check = c(
    "Balanced 5 x 5 sampling",
    "FEve values strictly inside (0,1)",
    "Final BAM converged",
    "Final BAM is full rank",
    "Basis dimensions acceptable"
  ),
  
  pass = c(
    qc_balanced_sample,
    qc_beta_response,
    qc_converged,
    qc_full_rank,
    qc_basis_dimension
  )
)


cat("\n============================================================\n")
cat("FEVE AUTOMATED QC SUMMARY\n")
cat("============================================================\n\n")


print(
  feve_qc_summary,
  row.names = FALSE
)


cat(
  "\nChecks passed: ",
  sum(
    feve_qc_summary$pass
  ),
  " of ",
  nrow(
    feve_qc_summary
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 59. Final analysis summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FUNCTIONAL EVENNESS ANALYSIS COMPLETE\n")
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
    feve_data$NAME10
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled localities:             ",
  n_distinct(
    feve_data$locality_id
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled checklists:             ",
  nrow(
    feve_data
  ),
  "\n",
  sep = ""
)


cat(
  "Additive AIC:                   ",
  round(
    AIC(
      feve_model_additive
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
      feve_model_interaction
    ),
    2
  ),
  "\n",
  sep = ""
)


cat(
  "Final adjusted R-squared:       ",
  round(
    feve_summary$r.sq,
    3
  ),
  "\n",
  sep = ""
)


cat(
  "Final deviance explained:       ",
  round(
    100 *
      feve_summary$dev.expl,
    2
  ),
  "%\n",
  sep = ""
)


cat(
  "Beta precision parameter:       ",
  round(
    feve_model_final$family$getTheta(
      TRUE
    ),
    3
  ),
  "\n",
  sep = ""
)


cat(
  "Basis-dimension follow-up:      ",
  k_followup_flag,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 60. Overall interpretation
# ------------------------------------------------------------
#
# The primary FEve model uses the same global aridity x habitat
# interaction structure used across the biodiversity-response
# analyses.
#
# The additive model has slightly lower AIC than the global
# interaction model. However, the difference is small and AIC
# comparisons among models fitted using fREML and different
# smooth structures are treated descriptively. The interaction
# model is retained because the aridity x habitat relationships
# represent a priori ecological hypotheses and maintaining the
# same global structure facilitates comparison among species
# richness, FRic, FDis, and FEve.
#
# Final FEve results indicate:
#
# - Aridity shows a clear association with functional evenness.
#   The main aridity smooth is approximately linear and is
#   statistically supported in the global model.
#
# - Because the global model also contains aridity x habitat
#   interactions, the full conditional aridity relationship can
#   show some curvature even though the main aridity smooth has
#   EDF approximately equal to 1.
#
# - Shrub/scrub cover shows the clearest local-habitat
#   association with FEve. The relationship is modestly
#   nonlinear and statistically supported.
#
# - Built cover shows weak evidence of an independent
#   relationship with FEve. Its p-value is close to the
#   conventional 0.05 threshold and should therefore be
#   interpreted cautiously.
#
# - Water cover also lies close to the conventional significance
#   threshold and provides only marginal evidence of an
#   independent FEve relationship in the global model.
#
# - None of the three aridity x habitat interaction surfaces
#   show strong statistical support.
#
# - The aridity x shrub/scrub interaction is marginal, while the
#   aridity x built interaction provides weaker suggestive
#   evidence. Neither is treated as a supported interaction.
#
# - Overall, the interaction results provide stronger evidence
#   for regional aridity and local habitat relationships than
#   for strong modification of habitat-FEve relationships along
#   the aridity gradient.
#
# - Checklist duration and distance traveled are both supported
#   sampling-effort controls.
#
# - Traveling versus Stationary protocol shows little
#   independent relationship with FEve after accounting for
#   habitat, effort, and hierarchical structure.
#
# - City and locality random-effect smooths account for
#   hierarchical variation in the eBird data.
#
# - Locality-level variation is substantially greater than
#   city-level variation, indicating important fine-scale
#   heterogeneity among sampling locations.
#
# - The final beta BAM explains approximately 31.8% of deviance,
#   with an adjusted R-squared of approximately 0.21 and an
#   estimated beta precision parameter of approximately 13.
#
# - Model convergence and rank are satisfactory.
#
# - The stable basis-dimension check identifies no smooth terms
#   requiring a larger basis.
#
# - Ecological-term concurvity is generally acceptable.
#   Elevated values involving aridity and its tensor-product
#   interactions are expected to some degree from the global
#   interaction parameterization, while the extreme concurvity
#   values involving city and locality are associated with the
#   hierarchical structure of the predictors and random-effect
#   bases.
#
# - Pearson residuals are generally centered near zero, although
#   a modest negative tail is present. Approximately 4% of
#   observations have absolute Pearson residuals greater than 2
#   and approximately 1% exceed 3.
#
# The publication figures deliberately distinguish descriptive
# raw-data relationships from adjusted model-based effects.
#
# OBSERVED FIGURES:
#
# - These show the actual checklist-level FEve observations and
#   an unadjusted single-predictor beta-GAM trend.
#
# - They describe the empirical relationship present in the
#   sampled data and do not control for other habitat variables,
#   sampling effort, city, or locality.
#
# PREDICTED FIGURES:
#
# - These show population-level predictions from the final
#   multivariable beta BAM.
#
# - Non-focal habitat and sampling-effort predictors are held at
#   representative median values.
#
# - Protocol is held at Stationary.
#
# - City and locality random-effect smooths are excluded.
#
# Keeping these figures separate allows the raw empirical
# patterns to be viewed independently from the adjusted
# relationships estimated by the final inferential model.
#
# Overall, functional evenness appears to respond to both the
# regional aridity gradient and local habitat structure,
# particularly shrub/scrub cover. However, there is relatively
# little evidence that local habitat relationships with
# functional evenness change strongly across the regional
# aridity gradient.
#
# ============================================================