# ============================================================
# Dryland Bird Project
# 06_CWM_Body_Mass.R
#
# Purpose:
# Model checklist-level community-weighted mean body mass
# (CWM.Mass) as a function of aridity, local habitat, and
# checklist effort using a Gamma GAM/BAM framework and balanced
# hierarchical sampling.
#
# CWM.Mass was calculated upstream from AVONET body-mass traits
# using abundance-weighted community means.
#
# Sampling design:
#
#   1. Retain Stationary and Traveling checklists
#   2. Stationary distance set to zero
#   3. Duration <= 300 minutes
#   4. Distance <= 5 km
#   5. Retain checklists with defined CWM.Mass
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
# The Gamma log model was strongly preferred because:
#
#   - CWM.Mass was strongly right-skewed
#   - Gamma had substantially lower AIC
#   - Gamma explained substantially more deviance
#   - Gamma residual behavior was more appropriate
#
# Therefore, only Gamma models are fit in this final script.
#
# Candidate models:
#
#   A. Additive Gamma BAM
#   B. Aridity x habitat interaction Gamma BAM
#
# The interaction model is retained as the primary/global model
# because the aridity x habitat terms represent the project's
# a priori ecological hypotheses and allow direct comparison
# among biodiversity responses.
#
# Publication / meeting outputs:
#
#   Outputs/
#     06_CWM_Body_Mass/
#       Figure/
#         Figure_06_Main_Data_Trends.png
#         Figure_06_Main_Data_Trends.pdf
#
#       Table/
#         Table_06_Main_Model_Results.csv
#
#       Data/
#         CWM_Body_Mass_Manuscript_Input.rds
#
# Main figure:
#
#   A. Aridity vs mean CWM body mass
#   B. Water cover vs mean CWM body mass
#   C. Shrub/scrub cover vs mean CWM body mass
#   D. Built cover vs mean CWM body mass
#
# Aridity is plotted at the city level because aridity is a
# city-level predictor.
#
# Habitat variables are plotted at the locality level because
# water, shrub/scrub, and built cover are locality-level
# predictors.
#
# Panels B-D color locality points by the percentile rank of
# the continuous city-level aridity index:
#   warmer colors = more arid
#   cooler colors = less arid
#
# Percentile rank is used only for visualization so the full
# warm-to-cool color range is represented across the observed
# data. The statistical model continues to use the original
# continuous aridity index.
#
# Figure smooths are descriptive and unadjusted.
# Statistical inference comes from the full Gamma BAM and is
# summarized in Table 06.
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
  "06_CWM_Body_Mass"
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

for (x in c(
  FIGURE_DIR,
  TABLE_DIR,
  DATA_DIR
)) {
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


combined_data <- readRDS(
  data_path
)


cat("\n============================================================\n")
cat("CWM BODY MASS ANALYSIS\n")
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
  "CWM.Mass",
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
  "PASS: All required CWM body-mass variables are present.\n"
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

mass_pool <- analysis_base %>%
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


mass_pool %>%
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
# 7. Initial CWM body-mass QC
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CWM BODY MASS INITIAL QC\n")
cat("============================================================\n\n")


cat(
  "CWM.Mass missing: ",
  sum(
    is.na(
      mass_pool$CWM.Mass
    )
  ),
  "\n",
  sep = ""
)


cat(
  "CWM.Mass missing (%): ",
  round(
    100 *
      mean(
        is.na(
          mass_pool$CWM.Mass
        )
      ),
    2
  ),
  "%\n",
  sep = ""
)


cat(
  "CWM.Mass equal to zero: ",
  sum(
    mass_pool$CWM.Mass == 0,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


cat(
  "CWM.Mass below zero: ",
  sum(
    mass_pool$CWM.Mass < 0,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)


print(
  summary(
    mass_pool$CWM.Mass
  )
)


# ------------------------------------------------------------
# 8. Detailed response quantiles
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CWM BODY MASS QUANTILES\n")
cat("============================================================\n\n")


mass_quantiles <- quantile(
  mass_pool$CWM.Mass,
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
  mass_quantiles
)


# ------------------------------------------------------------
# 9. Remove undefined response values
# ------------------------------------------------------------

mass_pool <- mass_pool %>%
  filter(
    !is.na(CWM.Mass)
  )


cat(
  "\nRows with defined CWM.Mass: ",
  nrow(mass_pool),
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
  mass_pool$CWM.Mass <= 0
)


cat("\n============================================================\n")
cat("POSITIVE-RESPONSE QC\n")
cat("============================================================\n\n")


cat(
  "CWM.Mass observations <= 0: ",
  n_nonpositive,
  "\n",
  sep = ""
)


if (
  n_nonpositive > 0
) {
  
  stop(
    paste0(
      "CWM.Mass contains zero or negative observations. ",
      "Review these values before fitting the Gamma BAM."
    )
  )
}


cat(
  "PASS: CWM.Mass is strictly positive.\n"
)


# ------------------------------------------------------------
# 11. Complete-case restriction
# ------------------------------------------------------------

model_variables <- c(
  "CWM.Mass",
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
      mass_pool[
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


mass_pool <- mass_pool %>%
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

locality_sampling <- mass_pool %>%
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
cat("CWM BODY MASS ELIGIBILITY\n")
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

eligible_pool <- mass_pool %>%
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


mass_data <- eligible_pool %>%
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

locality_qc <- mass_data %>%
  count(
    NAME10,
    locality_id,
    name = "n_checklists"
  )


city_qc <- mass_data %>%
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
cat("CWM BODY MASS BALANCED SAMPLE QC\n")
cat("============================================================\n\n")


cat(
  "Sampled cities:     ",
  n_distinct(
    mass_data$NAME10
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled localities: ",
  n_distinct(
    mass_data$locality_id
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled checklists: ",
  nrow(
    mass_data
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
    "CWM body-mass balanced sampling QC failed."
  )
}


cat(
  "\nPASS: CWM body-mass balanced sampling verified.\n"
)


# ------------------------------------------------------------
# 18. Prepare factors
# ------------------------------------------------------------

mass_data <- mass_data %>%
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
cat("FINAL CWM BODY MASS DISTRIBUTION\n")
cat("============================================================\n\n")


print(
  summary(
    mass_data$CWM.Mass
  )
)


mass_mean <- mean(
  mass_data$CWM.Mass
)


mass_sd <- sd(
  mass_data$CWM.Mass
)


mass_skewness <- mean(
  (
    (
      mass_data$CWM.Mass -
        mass_mean
    ) /
      mass_sd
  ) ^ 3
)


cat(
  "\nStandard deviation: ",
  mass_sd,
  "\n",
  sep = ""
)


cat(
  "Minimum: ",
  min(
    mass_data$CWM.Mass
  ),
  "\n",
  sep = ""
)


cat(
  "Maximum: ",
  max(
    mass_data$CWM.Mass
  ),
  "\n",
  sep = ""
)


cat(
  "Skewness: ",
  round(
    mass_skewness,
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

mass_formula_additive <- CWM.Mass ~
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

mass_formula_interaction <- CWM.Mass ~
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
cat("FITTING ADDITIVE CWM MASS GAMMA BAM\n")
cat("============================================================\n\n")


additive_fit_time <- system.time(
  
  mass_model_additive <- bam(
    formula =
      mass_formula_additive,
    
    family =
      Gamma(
        link = "log"
      ),
    
    data =
      mass_data,
    
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
cat("FITTING INTERACTION CWM MASS GAMMA BAM\n")
cat("============================================================\n\n")


interaction_fit_time <- system.time(
  
  mass_model_interaction <- bam(
    formula =
      mass_formula_interaction,
    
    family =
      Gamma(
        link = "log"
      ),
    
    data =
      mass_data,
    
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
    mass_model_additive
  )
)


cat("\n============================================================\n")
cat("INTERACTION MODEL SUMMARY\n")
cat("============================================================\n\n")


print(
  summary(
    mass_model_interaction
  )
)


# ------------------------------------------------------------
# 25. Candidate-model comparison
# ------------------------------------------------------------
#
# AIC is treated descriptively because these models were fit
# using fREML and differ in smooth structure.
# ------------------------------------------------------------

additive_summary <- summary(
  mass_model_additive
)


interaction_summary <- summary(
  mass_model_interaction
)


mass_model_comparison <- data.frame(
  model = c(
    "Additive",
    "Aridity x habitat interactions"
  ),
  
  AIC = c(
    AIC(
      mass_model_additive
    ),
    
    AIC(
      mass_model_interaction
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


mass_model_comparison <- mass_model_comparison %>%
  mutate(
    delta_AIC =
      AIC -
      min(
        AIC
      )
  )


cat("\n============================================================\n")
cat("CWM BODY MASS MODEL COMPARISON\n")
cat("============================================================\n\n")


print(
  mass_model_comparison
)


# ------------------------------------------------------------
# 26. Candidate-model convergence and rank
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CANDIDATE MODEL CONVERGENCE\n")
cat("============================================================\n\n")


cat(
  "Additive model converged: ",
  mass_model_additive$converged,
  "\n",
  sep = ""
)


cat(
  "Interaction model converged: ",
  mass_model_interaction$converged,
  "\n",
  sep = ""
)


cat(
  "Additive model rank: ",
  mass_model_additive$rank,
  " / ",
  length(
    coef(
      mass_model_additive
    )
  ),
  "\n",
  sep = ""
)


cat(
  "Interaction model rank: ",
  mass_model_interaction$rank,
  " / ",
  length(
    coef(
      mass_model_interaction
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
    mass_model_additive
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
    mass_model_interaction
  )
)


# ============================================================
# FINAL CWM BODY MASS MODEL
# ============================================================


# ------------------------------------------------------------
# 28. Select primary/global model
# ------------------------------------------------------------
#
# The interaction model is retained as the primary/global CWM
# body-mass model.
#
# Aridity x habitat terms represent the project's a priori
# ecological hypotheses. Maintaining the same global structure
# also facilitates comparison with species richness, FRic,
# FDis, and FEve.
#
# Individual interactions are interpreted from the global model
# rather than removed solely because they are not statistically
# significant.
# ------------------------------------------------------------

mass_model_final <- mass_model_interaction


mass_summary <- summary(
  mass_model_final
)


cat("\n============================================================\n")
cat("FINAL CWM BODY MASS MODEL\n")
cat("============================================================\n\n")


print(
  mass_summary
)


# ============================================================
# FINAL MODEL QC AND DIAGNOSTICS
# ============================================================


# ------------------------------------------------------------
# 29. Final convergence and rank
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CWM BODY MASS MODEL CONVERGENCE\n")
cat("============================================================\n\n")


cat(
  "Model converged: ",
  mass_model_final$converged,
  "\n",
  sep = ""
)


cat(
  "Model rank: ",
  mass_model_final$rank,
  " / ",
  length(
    coef(
      mass_model_final
    )
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 30. Stable basis-dimension check
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CWM BODY MASS STABLE K CHECK - 1000 REPLICATES\n")
cat("============================================================\n\n")


set.seed(
  K_CHECK_SEED
)


mass_k_stable <- k.check(
  mass_model_final,
  n.rep = 1000
)


print(
  mass_k_stable
)


k_valid_rows <-
  !is.na(
    mass_k_stable[
      ,
      "k-index"
    ]
  )


k_followup_flag <- any(
  mass_k_stable[
    k_valid_rows,
    "k-index"
  ] < 1 &
    
    mass_k_stable[
      k_valid_rows,
      "p-value"
    ] < 0.05 &
    
    (
      mass_k_stable[
        k_valid_rows,
        "edf"
      ] /
        mass_k_stable[
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
cat("CWM BODY MASS VARIANCE COMPONENTS\n")
cat("============================================================\n\n")


mass_vcomp <- gam.vcomp(
  mass_model_final
)


print(
  mass_vcomp
)


# ------------------------------------------------------------
# 32. Full concurvity
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CWM BODY MASS FULL CONCURVITY\n")
cat("============================================================\n\n")


mass_concurvity_full <- concurvity(
  mass_model_final,
  full = TRUE
)


print(
  mass_concurvity_full
)


# ------------------------------------------------------------
# 33. Pairwise concurvity
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CWM BODY MASS PAIRWISE CONCURVITY\n")
cat("============================================================\n\n")


mass_concurvity_pairwise <- concurvity(
  mass_model_final,
  full = FALSE
)


mass_pairwise_estimate <-
  mass_concurvity_pairwise$estimate


print(
  round(
    mass_pairwise_estimate,
    3
  )
)


# ------------------------------------------------------------
# 34. Ecological-term concurvity
# ------------------------------------------------------------

mass_ecological_terms <- c(
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


mass_ecological_concurvity <-
  mass_pairwise_estimate[
    mass_ecological_terms,
    mass_ecological_terms,
    drop = FALSE
  ]


cat("\n============================================================\n")
cat("CWM BODY MASS ECOLOGICAL-TERM CONCURVITY\n")
cat("============================================================\n\n")


print(
  round(
    mass_ecological_concurvity,
    3
  )
)


# ------------------------------------------------------------
# 35. Residual summaries
# ------------------------------------------------------------

mass_pearson_residuals <- residuals(
  mass_model_final,
  type = "pearson"
)


mass_deviance_residuals <- residuals(
  mass_model_final,
  type = "deviance"
)


mass_fitted_values <- fitted(
  mass_model_final
)


cat("\n============================================================\n")
cat("CWM BODY MASS RESIDUAL SUMMARY\n")
cat("============================================================\n\n")


cat(
  "Pearson residuals:\n"
)


print(
  summary(
    mass_pearson_residuals
  )
)


cat(
  "\nDeviance residuals:\n"
)


print(
  summary(
    mass_deviance_residuals
  )
)


# ------------------------------------------------------------
# 36. Extreme Pearson residuals
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CWM BODY MASS EXTREME RESIDUALS\n")
cat("============================================================\n\n")


cat(
  "Absolute Pearson residual > 2: ",
  sum(
    abs(
      mass_pearson_residuals
    ) > 2
  ),
  " (",
  round(
    100 *
      mean(
        abs(
          mass_pearson_residuals
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
      mass_pearson_residuals
    ) > 3
  ),
  " (",
  round(
    100 *
      mean(
        abs(
          mass_pearson_residuals
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

mass_diagnostic_data <- data.frame(
  fitted =
    mass_fitted_values,
  
  pearson_residual =
    mass_pearson_residuals
)


fig_mass_residuals <- ggplot(
  mass_diagnostic_data,
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
    x = "Fitted CWM body mass",
    y = "Pearson residual"
  ) +
  theme_classic(
    base_size = 12
  )


print(
  fig_mass_residuals
)


# ------------------------------------------------------------
# 38. QQ diagnostic
# ------------------------------------------------------------

set.seed(
  K_CHECK_SEED
)


qq.gam(
  mass_model_final,
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
  mass_model_final
)


# ------------------------------------------------------------
# 40. Final parametric-term table
# ------------------------------------------------------------

mass_parametric_results <- as.data.frame(
  mass_summary$p.table
)


mass_parametric_results$term <- rownames(
  mass_parametric_results
)


mass_parametric_results <- mass_parametric_results %>%
  select(
    term,
    everything()
  )


cat("\n============================================================\n")
cat("FINAL CWM BODY MASS PARAMETRIC TERMS\n")
cat("============================================================\n\n")


print(
  mass_parametric_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 41. Final smooth-term table
# ------------------------------------------------------------

mass_smooth_results <- as.data.frame(
  mass_summary$s.table
)


mass_smooth_results$term <- rownames(
  mass_smooth_results
)


mass_smooth_results <- mass_smooth_results %>%
  select(
    term,
    everything()
  )


cat("\n============================================================\n")
cat("FINAL CWM BODY MASS SMOOTH TERMS\n")
cat("============================================================\n\n")


print(
  mass_smooth_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 42. Final model fit
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL CWM BODY MASS MODEL FIT\n")
cat("============================================================\n\n")


cat(
  "Adjusted R-squared: ",
  round(
    mass_summary$r.sq,
    3
  ),
  "\n",
  sep = ""
)


cat(
  "Deviance explained: ",
  round(
    100 *
      mass_summary$dev.expl,
    2
  ),
  "%\n",
  sep = ""
)


cat(
  "Gamma scale parameter: ",
  round(
    mass_summary$scale,
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
  mass_model_final,
  select = 1,
  shade = TRUE,
  pages = 1
)


plot(
  mass_model_final,
  select = 2,
  shade = TRUE,
  pages = 1
)


plot(
  mass_model_final,
  select = 3,
  shade = TRUE,
  pages = 1
)


plot(
  mass_model_final,
  select = 4,
  shade = TRUE,
  pages = 1
)


plot(
  mass_model_final,
  select = 5,
  pages = 1
)


plot(
  mass_model_final,
  select = 6,
  pages = 1
)


plot(
  mass_model_final,
  select = 7,
  pages = 1
)


plot(
  mass_model_final,
  select = 8,
  shade = TRUE,
  pages = 1
)


plot(
  mass_model_final,
  select = 9,
  shade = TRUE,
  pages = 1
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
# Warm colors represent more-arid cities.
# Cool colors represent less-arid cities.
#
# Smooths are descriptive and unadjusted.
# Statistical inference is based on the final Gamma BAM and is
# reported in Table 06.
#
# ============================================================


# ------------------------------------------------------------
# 44. Create empirical trend datasets
# ------------------------------------------------------------

mass_city_plot <- mass_data %>%
  group_by(
    NAME10
  ) %>%
  summarise(
    avg_aridity_index =
      first(
        avg_aridity_index
      ),
    
    mean_mass =
      mean(
        CWM.Mass,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


mass_locality_plot <- mass_data %>%
  group_by(
    NAME10,
    locality_id
  ) %>%
  summarise(
    avg_aridity_index =
      first(
        avg_aridity_index
      ),
    
    water =
      first(
        water
      ),
    
    shrub_and_scrub =
      first(
        shrub_and_scrub
      ),
    
    built =
      first(
        built
      ),
    
    mean_mass =
      mean(
        CWM.Mass,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  mutate(
    aridity_rank =
      percent_rank(
        avg_aridity_index
      )
  )


# ------------------------------------------------------------
# 45. Convert habitat cover to percent
# ------------------------------------------------------------

maximum_habitat_value <- max(
  c(
    mass_locality_plot$water,
    mass_locality_plot$shrub_and_scrub,
    mass_locality_plot$built
  ),
  na.rm = TRUE
)


if (
  maximum_habitat_value <= 1.05
) {
  
  mass_locality_plot <- mass_locality_plot %>%
    mutate(
      water_plot =
        water * 100,
      
      shrub_plot =
        shrub_and_scrub * 100,
      
      built_plot =
        built * 100
    )
  
} else {
  
  mass_locality_plot <- mass_locality_plot %>%
    mutate(
      water_plot =
        water,
      
      shrub_plot =
        shrub_and_scrub,
      
      built_plot =
        built
    )
}


# ------------------------------------------------------------
# 46. Common figure settings
# ------------------------------------------------------------

main_y_limits <- range(
  c(
    mass_city_plot$mean_mass,
    mass_locality_plot$mean_mass
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


theme_mass_main <- function() {
  
  theme_classic(
    base_size = 12
  ) +
    theme(
      axis.text =
        element_text(
          colour = "black"
        ),
      
      plot.title =
        element_text(
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
# 47. Panel A: Aridity
# ------------------------------------------------------------

fig_mass_aridity <- ggplot(
  mass_city_plot,
  aes(
    x = avg_aridity_index,
    y = mean_mass
  )
) +
  geom_point(
    alpha = 0.45,
    size = 1.8
  ) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(
      x,
      k = 4
    ),
    se = TRUE,
    color = "black",
    fill = "grey80"
  ) +
  coord_cartesian(
    ylim =
      main_y_limits
  ) +
  labs(
    x = "Aridity index (higher = less arid)",
    y = "Mean CWM body mass",
    title = "Aridity"
  ) +
  theme_mass_main() +
  theme(
    legend.position = "none"
  )


# ------------------------------------------------------------
# 48. Panel B: Water cover
# ------------------------------------------------------------

fig_mass_water <- ggplot(
  mass_locality_plot,
  aes(
    x = water_plot,
    y = mean_mass
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
    formula = y ~ s(
      x,
      k = 4
    ),
    se = TRUE,
    color = "black",
    fill = "grey80"
  ) +
  aridity_color_scale +
  coord_cartesian(
    ylim =
      main_y_limits
  ) +
  labs(
    x = "Water cover (%)",
    y = "Mean CWM body mass",
    title = "Water"
  ) +
  theme_mass_main()


# ------------------------------------------------------------
# 49. Panel C: Shrub and scrub cover
# ------------------------------------------------------------

fig_mass_shrub <- ggplot(
  mass_locality_plot,
  aes(
    x = shrub_plot,
    y = mean_mass
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
    formula = y ~ s(
      x,
      k = 4
    ),
    se = TRUE,
    color = "black",
    fill = "grey80"
  ) +
  aridity_color_scale +
  coord_cartesian(
    ylim =
      main_y_limits
  ) +
  labs(
    x = "Shrub and scrub cover (%)",
    y = "Mean CWM body mass",
    title = "Shrub and scrub"
  ) +
  theme_mass_main()


# ------------------------------------------------------------
# 50. Panel D: Built cover
# ------------------------------------------------------------

fig_mass_built <- ggplot(
  mass_locality_plot,
  aes(
    x = built_plot,
    y = mean_mass
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
    formula = y ~ s(
      x,
      k = 4
    ),
    se = TRUE,
    color = "black",
    fill = "grey80"
  ) +
  aridity_color_scale +
  coord_cartesian(
    ylim =
      main_y_limits
  ) +
  labs(
    x = "Built cover (%)",
    y = "Mean CWM body mass",
    title = "Built"
  ) +
  theme_mass_main()


# ------------------------------------------------------------
# 51. Combine four empirical panels
# ------------------------------------------------------------

figure_06_main <- ggpubr::ggarrange(
  fig_mass_aridity,
  fig_mass_water,
  fig_mass_shrub,
  fig_mass_built,
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


print(
  figure_06_main
)


# ------------------------------------------------------------
# 52. Save main figure
# ------------------------------------------------------------

ggsave(
  filename =
    file.path(
      FIGURE_DIR,
      "Figure_06_Main_Data_Trends.png"
    ),
  plot =
    figure_06_main,
  width = 10,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)


ggsave(
  filename =
    file.path(
      FIGURE_DIR,
      "Figure_06_Main_Data_Trends.pdf"
    ),
  plot =
    figure_06_main,
  width = 10,
  height = 8,
  units = "in",
  bg = "white"
)


# ============================================================
# MAIN STATISTICAL TABLE
# ============================================================


# ------------------------------------------------------------
# 53. Focal ecological term labels
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
# 54. P-value formatting helper
# ------------------------------------------------------------

format_p_value <- function(
    p
) {
  
  ifelse(
    p < 0.001,
    "<0.001",
    sprintf(
      "%.3f",
      p
    )
  )
}


# ------------------------------------------------------------
# 55. Build main model-results table
# ------------------------------------------------------------

table_06_main <- mass_smooth_results %>%
  filter(
    term %in%
      names(
        focal_term_labels
      )
  ) %>%
  mutate(
    term_order =
      match(
        term,
        names(
          focal_term_labels
        )
      ),
    
    Predictor =
      unname(
        focal_term_labels[
          term
        ]
      ),
    
    EDF =
      sprintf(
        "%.2f",
        edf
      ),
    
    Reference_df =
      sprintf(
        "%.2f",
        Ref.df
      ),
    
    F_statistic =
      sprintf(
        "%.2f",
        F
      ),
    
    P_value =
      format_p_value(
        `p-value`
      )
  ) %>%
  arrange(
    term_order
  ) %>%
  select(
    Predictor,
    EDF,
    Reference_df,
    F_statistic,
    P_value
  )


names(
  table_06_main
) <- c(
  "Predictor",
  "EDF",
  "Ref. df",
  "F statistic",
  "P-value"
)


if (
  nrow(
    table_06_main
  ) !=
  length(
    focal_term_labels
  )
) {
  
  stop(
    "Table 06 does not contain all seven focal ecological terms."
  )
}


cat("\n============================================================\n")
cat("TABLE 06 - MAIN MODEL RESULTS\n")
cat("============================================================\n\n")


print(
  table_06_main,
  row.names = FALSE
)


# ------------------------------------------------------------
# 56. Save main statistical table
# ------------------------------------------------------------

write.csv(
  table_06_main,
  file =
    file.path(
      TABLE_DIR,
      "Table_06_Main_Model_Results.csv"
    ),
  row.names = FALSE
)


# ============================================================
# MANUSCRIPT-SYNTHESIS EXPORT
# ============================================================


# ------------------------------------------------------------
# 57. Build standardized manuscript-input bundle
# ------------------------------------------------------------

mass_model_metadata <- data.frame(
  response_id =
    "CWM_Body_Mass",
  
  response_label =
    "Community-weighted mean body mass",
  
  response_variable =
    "CWM.Mass",
  
  model_family =
    "Gamma",
  
  link =
    "log",
  
  n_cities =
    n_distinct(
      mass_data$NAME10
    ),
  
  n_localities =
    n_distinct(
      mass_data$locality_id
    ),
  
  n_checklists =
    nrow(
      mass_data
    ),
  
  adjusted_R2 =
    mass_summary$r.sq,
  
  deviance_explained =
    mass_summary$dev.expl,
  
  family_parameter =
    mass_summary$scale,
  
  family_parameter_name =
    "scale",
  
  response_skewness =
    mass_skewness,
  
  additive_AIC =
    AIC(
      mass_model_additive
    ),
  
  final_AIC =
    AIC(
      mass_model_final
    ),
  
  stringsAsFactors = FALSE
)


cwm_body_mass_manuscript_input <- list(
  response_id =
    "CWM_Body_Mass",
  
  response_label =
    "Community-weighted mean body mass",
  
  response_variable =
    "CWM.Mass",
  
  analysis_data =
    mass_data,
  
  final_model =
    mass_model_final,
  
  smooth_results =
    mass_smooth_results,
  
  parametric_results =
    mass_parametric_results,
  
  focal_results =
    table_06_main,
  
  model_metadata =
    mass_model_metadata,
  
  plot_data = list(
    city =
      mass_city_plot,
    
    locality =
      mass_locality_plot
  )
)


cwm_body_mass_manuscript_file <- file.path(
  DATA_DIR,
  "CWM_Body_Mass_Manuscript_Input.rds"
)


saveRDS(
  cwm_body_mass_manuscript_input,
  cwm_body_mass_manuscript_file
)


if (
  !file.exists(
    cwm_body_mass_manuscript_file
  )
) {
  
  stop(
    "CWM body-mass manuscript-input RDS was not created."
  )
}


cat("\n============================================================\n")
cat("MANUSCRIPT-SYNTHESIS EXPORT\n")
cat("============================================================\n\n")


cat(
  "Saved manuscript input:\n",
  cwm_body_mass_manuscript_file,
  "\n",
  sep = ""
)


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


qc_positive_response <-
  all(
    mass_data$CWM.Mass > 0
  )


qc_converged <-
  isTRUE(
    mass_model_final$converged
  )


qc_full_rank <-
  mass_model_final$rank ==
  length(
    coef(
      mass_model_final
    )
  )


qc_basis_dimension <-
  !k_followup_flag


qc_positive_fitted <-
  all(
    mass_fitted_values > 0
  )


mass_qc_summary <- data.frame(
  check = c(
    "Balanced 5 x 5 sampling",
    "CWM body mass strictly positive",
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
cat("CWM BODY MASS AUTOMATED QC SUMMARY\n")
cat("============================================================\n\n")


print(
  mass_qc_summary,
  row.names = FALSE
)


cat(
  "\nChecks passed: ",
  sum(
    mass_qc_summary$pass
  ),
  " of ",
  nrow(
    mass_qc_summary
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 59. Output check
# ------------------------------------------------------------

expected_output_files <- c(
  file.path(
    FIGURE_DIR,
    "Figure_06_Main_Data_Trends.png"
  ),
  
  file.path(
    FIGURE_DIR,
    "Figure_06_Main_Data_Trends.pdf"
  ),
  
  file.path(
    TABLE_DIR,
    "Table_06_Main_Model_Results.csv"
  ),
  
  file.path(
    DATA_DIR,
    "CWM_Body_Mass_Manuscript_Input.rds"
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
    file =
      expected_output_files,
    
    exists =
      output_check
  ),
  row.names = FALSE
)


if (
  !all(
    output_check
  )
) {
  
  stop(
    "One or more expected Script 06 outputs were not created."
  )
}


cat(
  "\nPASS: All Script 06 outputs were created.\n"
)


# ------------------------------------------------------------
# 60. Final analysis summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CWM BODY MASS ANALYSIS COMPLETE\n")
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
    mass_data$NAME10
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled localities:             ",
  n_distinct(
    mass_data$locality_id
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled checklists:             ",
  nrow(
    mass_data
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled CWM.Mass skewness:      ",
  round(
    mass_skewness,
    3
  ),
  "\n",
  sep = ""
)


cat(
  "Additive AIC:                   ",
  round(
    AIC(
      mass_model_additive
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
      mass_model_interaction
    ),
    2
  ),
  "\n",
  sep = ""
)


cat(
  "Final adjusted R-squared:       ",
  round(
    mass_summary$r.sq,
    3
  ),
  "\n",
  sep = ""
)


cat(
  "Final deviance explained:       ",
  round(
    100 *
      mass_summary$dev.expl,
    2
  ),
  "%\n",
  sep = ""
)


cat(
  "Gamma scale parameter:          ",
  round(
    mass_summary$scale,
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


cat(
  "\nPublication figure:\n",
  file.path(
    FIGURE_DIR,
    "Figure_06_Main_Data_Trends.png"
  ),
  "\n",
  sep = ""
)


cat(
  "\nPublication table:\n",
  file.path(
    TABLE_DIR,
    "Table_06_Main_Model_Results.csv"
  ),
  "\n",
  sep = ""
)


cat(
  "\nManuscript input:\n",
  cwm_body_mass_manuscript_file,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 61. Overall interpretation
# ------------------------------------------------------------
#
# Community-weighted mean body mass is a strictly positive and
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
# Based on the initial global Gamma model:
#
# - Water cover shows the clearest ecological relationship with
#   community-weighted mean body mass.
#
# - The water relationship is strongly supported and nonlinear.
#
# - There is little evidence for an independent relationship
#   between aridity and CWM body mass after accounting for local
#   habitat, sampling effort, and hierarchical structure.
#
# - Shrub/scrub cover shows little evidence of an independent
#   relationship with CWM body mass.
#
# - Built cover shows little evidence of an independent
#   relationship with CWM body mass.
#
# - None of the three aridity x habitat interaction surfaces
#   show strong statistical support.
#
# - Thus, current evidence is stronger for a direct local-water
#   relationship than for modification of habitat effects along
#   the regional aridity gradient.
#
# - Checklist duration and distance traveled are important
#   sampling-effort controls.
#
# - City and locality random-effect smooths account for
#   substantial hierarchical variation in checklist-level CWM
#   body mass.
#
# - The final Gamma model should be interpreted primarily using
#   deviance explained, residual diagnostics, smooth-term tests,
#   conditional predictions, and ecological effect shapes.
#
# - Traveling checklists are associated with higher CWM body
#   mass than Stationary checklists after accounting for habitat,
#   effort, and hierarchical structure.
#
# Figure interpretation:
#
# - Panel A shows the city-level empirical aridity-body-mass
#   relationship.
#
# - Panels B-D show locality-level habitat-body-mass
#   relationships, with point color representing percentile
#   rank along the continuous city-level aridity gradient.
#
# - Warm colors indicate more-arid cities and cool colors
#   indicate less-arid cities. Percentile rank is used only for
#   visualization; the BAM uses the original continuous aridity
#   index.
#
# - The descriptive GAM lines are not the source of inferential
#   tests; inference comes from the full Gamma BAM summarized in
#   Table 06.
#
#  Overall, community-weighted mean body mass appears to respond
#  primarily to local water cover rather than directly to the
#  regional aridity gradient or to strong aridity-dependent
#  modification of habitat relationships.
# 
# ============================================================