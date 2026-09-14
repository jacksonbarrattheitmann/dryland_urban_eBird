# ============================================================
# Dryland Bird Project
# 07_CWM_Beak_Length.R
#
# Purpose:
# Model checklist-level community-weighted mean beak length
# (CWM.Beak.Length) as a function of aridity, local habitat,
# and checklist effort using a Gamma GAM/BAM framework and
# balanced hierarchical sampling.
#
# CWM.Beak.Length was calculated upstream as the
# abundance-weighted community mean of the AVONET
# Beak.Length_Culmen trait.
#
# Sampling design:
#
#   1. Retain Stationary and Traveling checklists
#   2. Stationary distance set to zero
#   3. Duration <= 300 minutes
#   4. Distance <= 5 km
#   5. Retain checklists with defined CWM.Beak.Length
#   6. Require >= 5 valid checklists per locality
#   7. Require >= 5 eligible localities per city
#   8. Randomly sample 5 localities per city
#   9. Randomly sample 5 checklists per locality
#
# Response family:
#
# A prior family screen compared Gaussian identity and Gamma log
# models using the same global ecological model.
#
# The Gamma log model was preferred because:
#
#   - CWM.Beak.Length was strongly right-skewed
#   - Gamma had substantially lower AIC
#   - Gamma explained more deviance
#   - Gamma maintained positive fitted values
#   - Gamma diagnostics and basis dimensions were acceptable
#
# Therefore, only Gamma models are fit in this final script.
#
# Candidate models:
#
#   A. Additive Gamma BAM
#   B. Aridity x habitat interaction Gamma BAM
#
# The interaction model is retained as the primary/global model
# because aridity x habitat relationships represent the a priori
# ecological hypotheses and facilitate direct comparison among
# biodiversity responses.
#
# Figures:
#
#   - CWM beak-length response distribution
#   - Log-scaled response distribution
#
#   Observed/descriptive:
#     - CWM beak length vs aridity
#     - CWM beak length vs water cover
#     - CWM beak length vs shrub/scrub cover
#     - CWM beak length vs built cover
#
#   Adjusted/model-based:
#     - Predicted CWM beak length vs aridity
#     - Predicted CWM beak length vs water cover
#     - Predicted CWM beak length vs shrub/scrub cover
#     - Predicted CWM beak length vs built cover
#
# Observed and predicted relationships are intentionally shown
# in completely separate figures.
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


LOCALITY_SAMPLE_SEED <- 20260911

CHECKLIST_SAMPLE_SEED <- 20260912

K_CHECK_SEED <- 20260912


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
cat("CWM BEAK LENGTH ANALYSIS\n")
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
  "CWM.Beak.Length",
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
  "PASS: All required CWM beak-length variables are present.\n"
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

beak_pool <- analysis_base %>%
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


beak_pool %>%
  summarise(
    n_checklists =
      n(),
    
    n_localities =
      n_distinct(
        locality_id
      ),
    
    n_cities =
      n_distinct(
        NAME10
      )
  ) %>%
  print()


# ------------------------------------------------------------
# 7. Initial CWM beak-length QC
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CWM BEAK LENGTH INITIAL QC\n")
cat("============================================================\n\n")


cat(
  "CWM.Beak.Length missing: ",
  sum(
    is.na(
      beak_pool$CWM.Beak.Length
    )
  ),
  "\n",
  sep = ""
)


cat(
  "CWM.Beak.Length missing (%): ",
  round(
    100 *
      mean(
        is.na(
          beak_pool$CWM.Beak.Length
        )
      ),
    2
  ),
  "%\n",
  sep = ""
)


cat(
  "CWM.Beak.Length equal to zero: ",
  sum(
    beak_pool$CWM.Beak.Length == 0,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "CWM.Beak.Length below zero: ",
  sum(
    beak_pool$CWM.Beak.Length < 0,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


print(
  summary(
    beak_pool$CWM.Beak.Length
  )
)


# ------------------------------------------------------------
# 8. Detailed response quantiles
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CWM BEAK LENGTH QUANTILES\n")
cat("============================================================\n\n")


beak_quantiles <- quantile(
  beak_pool$CWM.Beak.Length,
  probs = c(
    0,
    0.01,
    0.05,
    0.25,
    0.50,
    0.75,
    0.95,
    0.99,
    1
  ),
  na.rm = TRUE
)


print(
  beak_quantiles
)


# ------------------------------------------------------------
# 9. Remove undefined response values
# ------------------------------------------------------------

beak_pool <- beak_pool %>%
  filter(
    !is.na(CWM.Beak.Length)
  )


cat(
  "\nRows with defined CWM.Beak.Length: ",
  nrow(beak_pool),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Positive-response gate
# ------------------------------------------------------------
#
# Gamma models require a strictly positive response.
# ------------------------------------------------------------

n_nonpositive <- sum(
  beak_pool$CWM.Beak.Length <= 0
)


cat("\n============================================================\n")
cat("POSITIVE-RESPONSE QC\n")
cat("============================================================\n\n")


cat(
  "CWM.Beak.Length observations <= 0: ",
  n_nonpositive,
  "\n",
  sep = ""
)


if (
  n_nonpositive > 0
) {
  
  stop(
    paste0(
      "CWM.Beak.Length contains zero or negative observations. ",
      "Review these values before fitting the Gamma BAM."
    )
  )
}


cat(
  "PASS: CWM.Beak.Length is strictly positive.\n"
)


# ------------------------------------------------------------
# 11. Complete-case restriction
# ------------------------------------------------------------

model_variables <- c(
  "CWM.Beak.Length",
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
      beak_pool[
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


beak_pool <- beak_pool %>%
  filter(
    if_all(
      all_of(
        model_variables
      ),
      ~ !is.na(.)
    )
  )


# ------------------------------------------------------------
# 12. Identify eligible localities
# ------------------------------------------------------------

locality_sampling <- beak_pool %>%
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
# 13. Identify eligible cities
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
cat("CWM BEAK LENGTH ELIGIBILITY\n")
cat("============================================================\n\n")


cat(
  "Eligible localities: ",
  nrow(
    eligible_localities
  ),
  "\n",
  sep = ""
)


cat(
  "Eligible cities:     ",
  nrow(
    eligible_cities
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 14. Create eligible sampling pool
# ------------------------------------------------------------

eligible_pool <- beak_pool %>%
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
# 15. Sample five localities per city
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
# 16. Sample five checklists per locality
# ------------------------------------------------------------

set.seed(
  CHECKLIST_SAMPLE_SEED
)


beak_data <- eligible_pool %>%
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
# 17. Balanced sample QC
# ------------------------------------------------------------

locality_qc <- beak_data %>%
  count(
    NAME10,
    locality_id,
    name = "n_checklists"
  )


city_qc <- beak_data %>%
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
cat("CWM BEAK LENGTH BALANCED SAMPLE QC\n")
cat("============================================================\n\n")


cat(
  "Sampled cities:     ",
  n_distinct(
    beak_data$NAME10
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled localities: ",
  n_distinct(
    beak_data$locality_id
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled checklists: ",
  nrow(
    beak_data
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
    "CWM beak-length balanced sampling QC failed."
  )
}


cat(
  "\nPASS: CWM beak-length balanced sampling verified.\n"
)


# ------------------------------------------------------------
# 18. Prepare factors
# ------------------------------------------------------------

beak_data <- beak_data %>%
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
# 19. Final sampled response distribution
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL CWM BEAK LENGTH DISTRIBUTION\n")
cat("============================================================\n\n")


print(
  summary(
    beak_data$CWM.Beak.Length
  )
)


beak_mean <- mean(
  beak_data$CWM.Beak.Length
)


beak_sd <- sd(
  beak_data$CWM.Beak.Length
)


beak_skewness <- mean(
  (
    (
      beak_data$CWM.Beak.Length -
        beak_mean
    ) /
      beak_sd
  ) ^ 3
)


cat(
  "\nStandard deviation: ",
  beak_sd,
  "\n",
  sep = ""
)


cat(
  "Minimum: ",
  min(
    beak_data$CWM.Beak.Length
  ),
  "\n",
  sep = ""
)


cat(
  "Maximum: ",
  max(
    beak_data$CWM.Beak.Length
  ),
  "\n",
  sep = ""
)


cat(
  "Skewness: ",
  round(
    beak_skewness,
    3
  ),
  "\n",
  sep = ""
)


# ============================================================
# CANDIDATE GAMMA BAM MODELS
# ============================================================


# ------------------------------------------------------------
# 20. Additive model formula
# ------------------------------------------------------------

beak_formula_additive <- CWM.Beak.Length ~
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

beak_formula_interaction <- CWM.Beak.Length ~
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
# 22. Fit additive Gamma BAM
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FITTING ADDITIVE CWM BEAK LENGTH GAMMA BAM\n")
cat("============================================================\n\n")


additive_fit_time <- system.time(
  
  beak_model_additive <- bam(
    formula =
      beak_formula_additive,
    
    family =
      Gamma(
        link = "log"
      ),
    
    data =
      beak_data,
    
    method =
      "fREML",
    
    discrete =
      TRUE
  )
  
)


print(
  additive_fit_time
)


# ------------------------------------------------------------
# 23. Fit interaction Gamma BAM
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FITTING INTERACTION CWM BEAK LENGTH GAMMA BAM\n")
cat("============================================================\n\n")


interaction_fit_time <- system.time(
  
  beak_model_interaction <- bam(
    formula =
      beak_formula_interaction,
    
    family =
      Gamma(
        link = "log"
      ),
    
    data =
      beak_data,
    
    method =
      "fREML",
    
    discrete =
      TRUE
  )
  
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
    beak_model_additive
  )
)


cat("\n============================================================\n")
cat("INTERACTION MODEL SUMMARY\n")
cat("============================================================\n\n")


print(
  summary(
    beak_model_interaction
  )
)


# ------------------------------------------------------------
# 25. Candidate-model comparison
# ------------------------------------------------------------
#
# AIC is treated descriptively because these models are fit
# using fREML and differ in smooth structure.
# ------------------------------------------------------------

additive_summary <- summary(
  beak_model_additive
)


interaction_summary <- summary(
  beak_model_interaction
)


beak_model_comparison <- data.frame(
  model = c(
    "Additive",
    "Aridity x habitat interactions"
  ),
  
  AIC = c(
    AIC(
      beak_model_additive
    ),
    
    AIC(
      beak_model_interaction
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


beak_model_comparison <- beak_model_comparison %>%
  mutate(
    delta_AIC =
      AIC -
      min(
        AIC
      )
  )


cat("\n============================================================\n")
cat("CWM BEAK LENGTH MODEL COMPARISON\n")
cat("============================================================\n\n")


print(
  beak_model_comparison
)


# ------------------------------------------------------------
# 26. Candidate-model convergence and rank
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CANDIDATE MODEL CONVERGENCE\n")
cat("============================================================\n\n")


cat(
  "Additive model converged: ",
  beak_model_additive$converged,
  "\n",
  sep = ""
)


cat(
  "Interaction model converged: ",
  beak_model_interaction$converged,
  "\n",
  sep = ""
)


cat(
  "Additive model rank: ",
  beak_model_additive$rank,
  " / ",
  length(
    coef(
      beak_model_additive
    )
  ),
  "\n",
  sep = ""
)


cat(
  "Interaction model rank: ",
  beak_model_interaction$rank,
  " / ",
  length(
    coef(
      beak_model_interaction
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
    beak_model_additive
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
    beak_model_interaction
  )
)


# ============================================================
# FINAL CWM BEAK LENGTH MODEL
# ============================================================


# ------------------------------------------------------------
# 28. Select primary/global model
# ------------------------------------------------------------
#
# The interaction model is retained as the primary/global
# CWM beak-length model.
#
# The aridity x habitat terms represent the project's a priori
# ecological hypotheses.
#
# Retaining the same global structure facilitates comparison
# with species richness, FRic, FDis, FEve, and CWM body mass.
#
# Individual interactions are interpreted from the global model
# rather than removed solely because they are not statistically
# significant.
# ------------------------------------------------------------

beak_model_final <- beak_model_interaction


beak_summary <- summary(
  beak_model_final
)


cat("\n============================================================\n")
cat("FINAL CWM BEAK LENGTH MODEL\n")
cat("============================================================\n\n")


print(
  beak_summary
)


# ============================================================
# FINAL MODEL QC AND DIAGNOSTICS
# ============================================================


# ------------------------------------------------------------
# 29. Final convergence and rank
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CWM BEAK LENGTH MODEL CONVERGENCE\n")
cat("============================================================\n\n")


cat(
  "Model converged: ",
  beak_model_final$converged,
  "\n",
  sep = ""
)


cat(
  "Model rank: ",
  beak_model_final$rank,
  " / ",
  length(
    coef(
      beak_model_final
    )
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 30. Stable basis-dimension check
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CWM BEAK LENGTH STABLE K CHECK - 1000 REPLICATES\n")
cat("============================================================\n\n")


set.seed(
  K_CHECK_SEED
)


beak_k_stable <- k.check(
  beak_model_final,
  n.rep = 1000
)


print(
  beak_k_stable
)


k_valid_rows <-
  !is.na(
    beak_k_stable[
      ,
      "k-index"
    ]
  )


k_followup_flag <- any(
  beak_k_stable[
    k_valid_rows,
    "k-index"
  ] < 1 &
    
    beak_k_stable[
      k_valid_rows,
      "p-value"
    ] < 0.05 &
    
    (
      beak_k_stable[
        k_valid_rows,
        "edf"
      ] /
        beak_k_stable[
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
cat("CWM BEAK LENGTH VARIANCE COMPONENTS\n")
cat("============================================================\n\n")


beak_vcomp <- gam.vcomp(
  beak_model_final
)


print(
  beak_vcomp
)


# ------------------------------------------------------------
# 32. Full concurvity
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CWM BEAK LENGTH FULL CONCURVITY\n")
cat("============================================================\n\n")


beak_concurvity_full <- concurvity(
  beak_model_final,
  full = TRUE
)


print(
  beak_concurvity_full
)


# ------------------------------------------------------------
# 33. Pairwise concurvity
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CWM BEAK LENGTH PAIRWISE CONCURVITY\n")
cat("============================================================\n\n")


beak_concurvity_pairwise <- concurvity(
  beak_model_final,
  full = FALSE
)


beak_pairwise_estimate <-
  beak_concurvity_pairwise$estimate


print(
  round(
    beak_pairwise_estimate,
    3
  )
)


# ------------------------------------------------------------
# 34. Ecological-term concurvity
# ------------------------------------------------------------

beak_ecological_terms <- c(
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


beak_ecological_concurvity <-
  beak_pairwise_estimate[
    beak_ecological_terms,
    beak_ecological_terms,
    drop = FALSE
  ]


cat("\n============================================================\n")
cat("CWM BEAK LENGTH ECOLOGICAL-TERM CONCURVITY\n")
cat("============================================================\n\n")


print(
  round(
    beak_ecological_concurvity,
    3
  )
)


# ------------------------------------------------------------
# 35. Residual summaries
# ------------------------------------------------------------

beak_pearson_residuals <- residuals(
  beak_model_final,
  type = "pearson"
)


beak_deviance_residuals <- residuals(
  beak_model_final,
  type = "deviance"
)


beak_fitted_values <- fitted(
  beak_model_final
)


beak_standardized_residuals <-
  beak_pearson_residuals /
  sqrt(
    beak_summary$scale
  )


cat("\n============================================================\n")
cat("CWM BEAK LENGTH RESIDUAL SUMMARY\n")
cat("============================================================\n\n")


cat(
  "Pearson residuals:\n"
)


print(
  summary(
    beak_pearson_residuals
  )
)


cat(
  "\nStandardized Pearson residuals:\n"
)


print(
  summary(
    beak_standardized_residuals
  )
)


cat(
  "\nDeviance residuals:\n"
)


print(
  summary(
    beak_deviance_residuals
  )
)


# ------------------------------------------------------------
# 36. Extreme standardized Pearson residuals
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CWM BEAK LENGTH EXTREME RESIDUALS\n")
cat("============================================================\n\n")


cat(
  "Absolute standardized Pearson residual > 2: ",
  sum(
    abs(
      beak_standardized_residuals
    ) > 2
  ),
  " (",
  round(
    100 *
      mean(
        abs(
          beak_standardized_residuals
        ) > 2
      ),
    2
  ),
  "%)\n",
  sep = ""
)


cat(
  "Absolute standardized Pearson residual > 3: ",
  sum(
    abs(
      beak_standardized_residuals
    ) > 3
  ),
  " (",
  round(
    100 *
      mean(
        abs(
          beak_standardized_residuals
        ) > 3
      ),
    2
  ),
  "%)\n",
  sep = ""
)


# ------------------------------------------------------------
# 37. Residual-vs-fitted diagnostic
# ------------------------------------------------------------

beak_diagnostic_data <- data.frame(
  fitted =
    beak_fitted_values,
  
  standardized_residual =
    beak_standardized_residuals
)


fig_beak_residuals <- ggplot(
  beak_diagnostic_data,
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
    x = "Fitted CWM beak length",
    y = "Standardized Pearson residual"
  ) +
  theme_classic(
    base_size = 12
  )


print(
  fig_beak_residuals
)


# ------------------------------------------------------------
# 38. QQ diagnostic
# ------------------------------------------------------------

set.seed(
  K_CHECK_SEED
)


qq.gam(
  beak_model_final,
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
  beak_model_final
)


# ------------------------------------------------------------
# 40. Final parametric-term table
# ------------------------------------------------------------

beak_parametric_results <- as.data.frame(
  beak_summary$p.table
)


beak_parametric_results$term <- rownames(
  beak_parametric_results
)


beak_parametric_results <- beak_parametric_results %>%
  select(
    term,
    everything()
  )


cat("\n============================================================\n")
cat("FINAL CWM BEAK LENGTH PARAMETRIC TERMS\n")
cat("============================================================\n\n")


print(
  beak_parametric_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 41. Final smooth-term table
# ------------------------------------------------------------

beak_smooth_results <- as.data.frame(
  beak_summary$s.table
)


beak_smooth_results$term <- rownames(
  beak_smooth_results
)


beak_smooth_results <- beak_smooth_results %>%
  select(
    term,
    everything()
  )


cat("\n============================================================\n")
cat("FINAL CWM BEAK LENGTH SMOOTH TERMS\n")
cat("============================================================\n\n")


print(
  beak_smooth_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 42. Final model fit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL CWM BEAK LENGTH MODEL FIT\n")
cat("============================================================\n\n")


cat(
  "Adjusted R-squared: ",
  round(
    beak_summary$r.sq,
    3
  ),
  "\n",
  sep = ""
)


cat(
  "Deviance explained: ",
  round(
    100 *
      beak_summary$dev.expl,
    2
  ),
  "%\n",
  sep = ""
)


cat(
  "Gamma scale parameter: ",
  round(
    beak_summary$scale,
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
  beak_model_final,
  select = 1,
  shade = TRUE,
  pages = 1
)


plot(
  beak_model_final,
  select = 2,
  shade = TRUE,
  pages = 1
)


plot(
  beak_model_final,
  select = 3,
  shade = TRUE,
  pages = 1
)


plot(
  beak_model_final,
  select = 4,
  shade = TRUE,
  pages = 1
)


plot(
  beak_model_final,
  select = 5,
  pages = 1
)


plot(
  beak_model_final,
  select = 6,
  pages = 1
)


plot(
  beak_model_final,
  select = 7,
  pages = 1
)


plot(
  beak_model_final,
  select = 8,
  shade = TRUE,
  pages = 1
)


plot(
  beak_model_final,
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
#   Actual checklist-level CWM beak-length observations plus a
#   descriptive, unadjusted single-predictor Gamma GAM.
#
# PREDICTED figures:
#
#   Population-level predictions from the final multivariable
#   Gamma BAM.
#
# Observed and predicted figures are intentionally kept
# completely separate.
#
# ============================================================


# ------------------------------------------------------------
# 44. Response-distribution figure
# ------------------------------------------------------------

fig_beak_distribution <- ggplot(
  beak_data,
  aes(
    x = CWM.Beak.Length
  )
) +
  geom_histogram(
    bins = 40,
    linewidth = 0.25
  ) +
  labs(
    x = "Community-weighted mean beak length",
    y = "Number of checklists"
  ) +
  theme_classic(
    base_size = 12
  )


print(
  fig_beak_distribution
)


# ------------------------------------------------------------
# 45. Log-scaled response-distribution figure
# ------------------------------------------------------------

fig_beak_distribution_logx <- ggplot(
  beak_data,
  aes(
    x = CWM.Beak.Length
  )
) +
  geom_histogram(
    bins = 40,
    linewidth = 0.25
  ) +
  scale_x_log10() +
  labs(
    x = "Community-weighted mean beak length (log-scaled axis)",
    y = "Number of checklists"
  ) +
  theme_classic(
    base_size = 12
  )


print(
  fig_beak_distribution_logx
)


# ------------------------------------------------------------
# 46. Determine habitat scaling
# ------------------------------------------------------------

maximum_habitat_value <- max(
  c(
    beak_data$water,
    beak_data$shrub_and_scrub,
    beak_data$built
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
# 47. Helper for observed relationships
# ------------------------------------------------------------
#
# Each observed figure contains:
#
#   - actual checklist-level CWM.Beak.Length observations
#   - an unadjusted single-predictor Gamma GAM
#   - 95% confidence interval
#
# The trend is displayed across the central 90% of the focal
# predictor distribution.
#
# These figures describe raw empirical relationships.
# They are not adjusted effects from the final model.
# ------------------------------------------------------------

make_beak_observed <- function(
    variable,
    x_label,
    scale_factor = 1
) {
  
  observed_data <- data.frame(
    x =
      beak_data[[variable]] *
      scale_factor,
    
    CWM.Beak.Length =
      beak_data$CWM.Beak.Length
  )
  
  
  raw_formula <- as.formula(
    paste0(
      "CWM.Beak.Length ~ s(",
      variable,
      ", k = 5)"
    )
  )
  
  
  raw_model <- gam(
    formula =
      raw_formula,
    
    family =
      Gamma(
        link = "log"
      ),
    
    data =
      beak_data,
    
    method =
      "REML"
  )
  
  
  x_limits <- quantile(
    beak_data[[variable]],
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
        y = CWM.Beak.Length
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
      y = "Community-weighted mean beak length"
    ) +
    theme_classic(
      base_size = 12
    )
}


# ------------------------------------------------------------
# 48. Observed beak length vs aridity
# ------------------------------------------------------------

fig_beak_aridity_observed <- make_beak_observed(
  variable =
    "avg_aridity_index",
  
  x_label =
    "Aridity index (higher = less arid)"
)


print(
  fig_beak_aridity_observed
)


# ------------------------------------------------------------
# 49. Observed beak length vs water cover
# ------------------------------------------------------------

fig_beak_water_observed <- make_beak_observed(
  variable =
    "water",
  
  x_label =
    "Water cover (%)",
  
  scale_factor =
    habitat_scale_factor
)


print(
  fig_beak_water_observed
)


# ------------------------------------------------------------
# 50. Observed beak length vs shrub/scrub cover
# ------------------------------------------------------------

fig_beak_shrub_observed <- make_beak_observed(
  variable =
    "shrub_and_scrub",
  
  x_label =
    "Shrub and scrub cover (%)",
  
  scale_factor =
    habitat_scale_factor
)


print(
  fig_beak_shrub_observed
)


# ------------------------------------------------------------
# 51. Observed beak length vs built cover
# ------------------------------------------------------------

fig_beak_built_observed <- make_beak_observed(
  variable =
    "built",
  
  x_label =
    "Built cover (%)",
  
  scale_factor =
    habitat_scale_factor
)


print(
  fig_beak_built_observed
)


# ============================================================
# MODEL-ADJUSTED RELATIONSHIPS
# ============================================================


# ------------------------------------------------------------
# 52. Reference values for adjusted predictions
# ------------------------------------------------------------

reference_aridity <- median(
  beak_data$avg_aridity_index,
  na.rm = TRUE
)


reference_duration <- median(
  beak_data$duration_minutes,
  na.rm = TRUE
)


reference_distance <- median(
  beak_data$effort_distance_km,
  na.rm = TRUE
)


reference_water <- median(
  beak_data$water,
  na.rm = TRUE
)


reference_shrub <- median(
  beak_data$shrub_and_scrub,
  na.rm = TRUE
)


reference_built <- median(
  beak_data$built,
  na.rm = TRUE
)


# ------------------------------------------------------------
# 53. Helper for adjusted predictions
# ------------------------------------------------------------
#
# Predictions and confidence intervals are calculated on the
# Gamma log-link scale and transformed back to the original
# CWM beak-length scale.
#
# Non-focal continuous predictors are held at median values.
# Protocol is held at Stationary.
#
# City and locality random-effect smooths are excluded to
# produce population-level relationships.
# ------------------------------------------------------------

make_beak_prediction <- function(
    variable,
    x_label,
    scale_factor = 1
) {
  
  x_limits <- quantile(
    beak_data[[variable]],
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
        beak_data$protocol_name
      )
  )
  
  
  prediction_data$NAME10 <- factor(
    rep(
      levels(
        beak_data$NAME10
      )[1],
      n_prediction
    ),
    levels =
      levels(
        beak_data$NAME10
      )
  )
  
  
  prediction_data$locality_id <- factor(
    rep(
      levels(
        beak_data$locality_id
      )[1],
      n_prediction
    ),
    levels =
      levels(
        beak_data$locality_id
      )
  )
  
  
  prediction_data[[variable]] <-
    x_sequence
  
  
  model_prediction <- predict(
    beak_model_final,
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
    beak_model_final$family$linkinv
  
  
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
      y = "Predicted community-weighted mean beak length"
    ) +
    theme_classic(
      base_size = 12
    )
}


# ------------------------------------------------------------
# 54. Predicted beak length vs aridity
# ------------------------------------------------------------

fig_beak_aridity_predicted <- make_beak_prediction(
  variable =
    "avg_aridity_index",
  
  x_label =
    "Aridity index (higher = less arid)"
)


print(
  fig_beak_aridity_predicted
)


# ------------------------------------------------------------
# 55. Predicted beak length vs water cover
# ------------------------------------------------------------

fig_beak_water_predicted <- make_beak_prediction(
  variable =
    "water",
  
  x_label =
    "Water cover (%)",
  
  scale_factor =
    habitat_scale_factor
)


print(
  fig_beak_water_predicted
)


# ------------------------------------------------------------
# 56. Predicted beak length vs shrub/scrub cover
# ------------------------------------------------------------

fig_beak_shrub_predicted <- make_beak_prediction(
  variable =
    "shrub_and_scrub",
  
  x_label =
    "Shrub and scrub cover (%)",
  
  scale_factor =
    habitat_scale_factor
)


print(
  fig_beak_shrub_predicted
)


# ------------------------------------------------------------
# 57. Predicted beak length vs built cover
# ------------------------------------------------------------

fig_beak_built_predicted <- make_beak_prediction(
  variable =
    "built",
  
  x_label =
    "Built cover (%)",
  
  scale_factor =
    habitat_scale_factor
)


print(
  fig_beak_built_predicted
)


# ------------------------------------------------------------
# 58. Optional figure saving
# ------------------------------------------------------------

if (
  SAVE_FIGURES
) {
  
  figure_directory <- file.path(
    "Outputs",
    "07_CWM_Beak_Length",
    "Figures"
  )
  
  
  dir.create(
    figure_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  ggsave(
    file.path(
      figure_directory,
      "CWM_Beak_Length_Distribution.pdf"
    ),
    fig_beak_distribution,
    width = 5.5,
    height = 4
  )
  
  
  ggsave(
    file.path(
      figure_directory,
      "CWM_Beak_Length_Distribution_LogX.pdf"
    ),
    fig_beak_distribution_logx,
    width = 5.5,
    height = 4
  )
  
  
  ggsave(
    file.path(
      figure_directory,
      "CWM_Beak_Length_Aridity_Observed.pdf"
    ),
    fig_beak_aridity_observed,
    width = 5.5,
    height = 4
  )
  
  
  ggsave(
    file.path(
      figure_directory,
      "CWM_Beak_Length_Aridity_Predicted.pdf"
    ),
    fig_beak_aridity_predicted,
    width = 5.5,
    height = 4
  )
  
  
  ggsave(
    file.path(
      figure_directory,
      "CWM_Beak_Length_Water_Observed.pdf"
    ),
    fig_beak_water_observed,
    width = 5.5,
    height = 4
  )
  
  
  ggsave(
    file.path(
      figure_directory,
      "CWM_Beak_Length_Water_Predicted.pdf"
    ),
    fig_beak_water_predicted,
    width = 5.5,
    height = 4
  )
  
  
  ggsave(
    file.path(
      figure_directory,
      "CWM_Beak_Length_Shrub_Observed.pdf"
    ),
    fig_beak_shrub_observed,
    width = 5.5,
    height = 4
  )
  
  
  ggsave(
    file.path(
      figure_directory,
      "CWM_Beak_Length_Shrub_Predicted.pdf"
    ),
    fig_beak_shrub_predicted,
    width = 5.5,
    height = 4
  )
  
  
  ggsave(
    file.path(
      figure_directory,
      "CWM_Beak_Length_Built_Observed.pdf"
    ),
    fig_beak_built_observed,
    width = 5.5,
    height = 4
  )
  
  
  ggsave(
    file.path(
      figure_directory,
      "CWM_Beak_Length_Built_Predicted.pdf"
    ),
    fig_beak_built_predicted,
    width = 5.5,
    height = 4
  )
}


# ============================================================
# FINAL QC AND ANALYSIS SUMMARY
# ============================================================


# ------------------------------------------------------------
# 59. Automated QC summary
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
    beak_data$CWM.Beak.Length > 0
  )


qc_converged <-
  isTRUE(
    beak_model_final$converged
  )


qc_full_rank <-
  beak_model_final$rank ==
  length(
    coef(
      beak_model_final
    )
  )


qc_basis_dimension <-
  !k_followup_flag


qc_positive_fitted <-
  all(
    beak_fitted_values > 0
  )


beak_qc_summary <- data.frame(
  check = c(
    "Balanced 5 x 5 sampling",
    "CWM beak length strictly positive",
    "Final BAM converged",
    "Final BAM is full rank",
    "Basis dimensions acceptable",
    "Gamma fitted values positive"
  ),
  
  pass = c(
    qc_balanced_sample,
    qc_positive_response,
    qc_converged,
    qc_full_rank,
    qc_basis_dimension,
    qc_positive_fitted
  )
)


cat("\n============================================================\n")
cat("CWM BEAK LENGTH AUTOMATED QC SUMMARY\n")
cat("============================================================\n\n")


print(
  beak_qc_summary,
  row.names = FALSE
)


cat(
  "\nChecks passed: ",
  sum(
    beak_qc_summary$pass
  ),
  " of ",
  nrow(
    beak_qc_summary
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 60. Final analysis summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CWM BEAK LENGTH ANALYSIS COMPLETE\n")
cat("============================================================\n\n")


cat(
  "Eligible cities:                  ",
  nrow(
    eligible_cities
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled cities:                   ",
  n_distinct(
    beak_data$NAME10
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled localities:               ",
  n_distinct(
    beak_data$locality_id
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled checklists:               ",
  nrow(
    beak_data
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled CWM.Beak.Length skewness: ",
  round(
    beak_skewness,
    3
  ),
  "\n",
  sep = ""
)


cat(
  "Additive AIC:                     ",
  round(
    AIC(
      beak_model_additive
    ),
    2
  ),
  "\n",
  sep = ""
)


cat(
  "Interaction AIC:                  ",
  round(
    AIC(
      beak_model_interaction
    ),
    2
  ),
  "\n",
  sep = ""
)


cat(
  "Final adjusted R-squared:         ",
  round(
    beak_summary$r.sq,
    3
  ),
  "\n",
  sep = ""
)


cat(
  "Final deviance explained:         ",
  round(
    100 *
      beak_summary$dev.expl,
    2
  ),
  "%\n",
  sep = ""
)


cat(
  "Gamma scale parameter:            ",
  round(
    beak_summary$scale,
    3
  ),
  "\n",
  sep = ""
)


cat(
  "Basis-dimension follow-up:        ",
  k_followup_flag,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 61. Overall interpretation
# ------------------------------------------------------------
#
# Community-weighted mean beak length is a strictly positive and
# strongly right-skewed response.
#
# A prior working-family screen strongly favored a Gamma model
# with a log link over a Gaussian identity model. Therefore,
# Gamma models are used throughout the final analysis.
#
# The primary model retains the same global aridity x habitat
# interaction structure used across the biodiversity-response
# analyses.
#
# Final CWM beak-length results indicate:
#
# - Water cover shows the clearest ecological relationship with
#   community-weighted mean beak length.
#
# - The water-cover relationship is strongly supported and
#   nonlinear.
#
# - Built cover also shows evidence of an independent nonlinear
#   relationship with CWM beak length.
#
# - There is little evidence for an independent relationship
#   between aridity and CWM beak length after accounting for
#   local habitat, sampling effort, and hierarchical structure.
#
# - Shrub/scrub cover shows little evidence of an independent
#   relationship with CWM beak length.
#
# - The aridity x water interaction is not strongly supported
#   under the final Gamma model.
#
# - The aridity x shrub/scrub interaction is also unsupported.
#
# - The aridity x built interaction provides weak or marginal
#   evidence but should not be treated as a supported
#   interaction.
#
# - Overall, current evidence is stronger for direct local
#   habitat relationships, particularly water and built cover,
#   than for strong modification of habitat relationships along
#   the regional aridity gradient.
#
# - Checklist duration and distance traveled are important
#   sampling-effort controls.
#
# - Traveling checklists are associated with higher CWM beak
#   length than Stationary checklists after accounting for
#   habitat, effort, and hierarchical structure.
#
# - City and locality random-effect smooths account for
#   substantial hierarchical variation in checklist-level CWM
#   beak length.
#
# - The final Gamma model should be interpreted primarily using
#   deviance explained, smooth-term tests, conditional
#   predictions, residual diagnostics, and ecological effect
#   shapes.
#
# - Stable basis-dimension checks are used to confirm that the
#   chosen smooth bases are sufficiently flexible.
#
# - Standardized Pearson residuals are used when evaluating
#   residual thresholds because the Gamma family includes an
#   estimated dispersion parameter.
#
# PUBLICATION FIGURES:
#
# Observed figures show the actual checklist-level CWM
# beak-length observations and an unadjusted single-predictor
# Gamma-GAM trend.
#
# Predicted figures show population-level relationships from the
# final multivariable Gamma BAM.
#
# For adjusted figures:
#
#   - non-focal continuous predictors are held at median values
#   - protocol is held at Stationary
#   - city and locality random-effect smooths are excluded
#
# Keeping observed and adjusted figures separate allows the raw
# empirical pattern to be compared directly with the relationship
# remaining after adjustment for habitat, effort, and
# hierarchical structure.
#
# Overall, community-weighted mean beak length appears to respond
# primarily to local habitat structure, particularly water cover
# and secondarily built cover, rather than directly to the
# regional aridity gradient.
#
# ============================================================