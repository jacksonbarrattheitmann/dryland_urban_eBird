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
# Publication / meeting outputs:
#
#   Outputs/
#     04_Functional_Dispersion/
#       Figure/
#         Figure_04_Main_Data_Trends.png
#         Figure_04_Main_Data_Trends.pdf
#
#       Table/
#         Table_04_Main_Model_Results.csv
#
#       Data/
#         Functional_Dispersion_Manuscript_Input.rds
#
# Main figure:
#
#   A. Aridity vs mean functional dispersion
#   B. Water cover vs mean functional dispersion
#   C. Shrub/scrub cover vs mean functional dispersion
#   D. Built cover vs mean functional dispersion
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
# Statistical inference comes from the full Gaussian BAM and
# is summarized in Table 04.
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


# Same sampling seeds used across biodiversity responses.

LOCALITY_SAMPLE_SEED <- 20260911

CHECKLIST_SAMPLE_SEED <- 20260912

K_CHECK_SEED <- 20260912


# Temporary duplicate handling.

TEMPORARY_DEDUPLICATE <- TRUE


# ------------------------------------------------------------
# 3. Output directories
# ------------------------------------------------------------

OUTPUT_ROOT <- file.path(
  "Outputs",
  "04_Functional_Dispersion"
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
# Statistical inference is based on the final Gaussian BAM and
# is reported in Table 04.
#
# ============================================================


# ------------------------------------------------------------
# 41. Create empirical trend datasets
# ------------------------------------------------------------

fdis_city_plot <- fdis_data %>%
  group_by(
    NAME10
  ) %>%
  summarise(
    avg_aridity_index =
      first(
        avg_aridity_index
      ),
    
    mean_fdis =
      mean(
        FDis,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


fdis_locality_plot <- fdis_data %>%
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
    
    mean_fdis =
      mean(
        FDis,
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
# 42. Convert habitat cover to percent
# ------------------------------------------------------------

maximum_habitat_value <- max(
  c(
    fdis_locality_plot$water,
    fdis_locality_plot$shrub_and_scrub,
    fdis_locality_plot$built
  ),
  na.rm = TRUE
)


if (
  maximum_habitat_value <= 1.05
) {
  
  fdis_locality_plot <- fdis_locality_plot %>%
    mutate(
      water_plot =
        water * 100,
      
      shrub_plot =
        shrub_and_scrub * 100,
      
      built_plot =
        built * 100
    )
  
} else {
  
  fdis_locality_plot <- fdis_locality_plot %>%
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
# 43. Common figure settings
# ------------------------------------------------------------

main_y_limits <- range(
  c(
    fdis_city_plot$mean_fdis,
    fdis_locality_plot$mean_fdis
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


theme_fdis_main <- function() {
  
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
# 44. Panel A: Aridity
# ------------------------------------------------------------

fig_fdis_aridity <- ggplot(
  fdis_city_plot,
  aes(
    x = avg_aridity_index,
    y = mean_fdis
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
    y = "Mean functional dispersion",
    title = "Aridity"
  ) +
  theme_fdis_main() +
  theme(
    legend.position = "none"
  )


# ------------------------------------------------------------
# 45. Panel B: Water cover
# ------------------------------------------------------------

fig_fdis_water <- ggplot(
  fdis_locality_plot,
  aes(
    x = water_plot,
    y = mean_fdis
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
    y = "Mean functional dispersion",
    title = "Water"
  ) +
  theme_fdis_main()


# ------------------------------------------------------------
# 46. Panel C: Shrub and scrub cover
# ------------------------------------------------------------

fig_fdis_shrub <- ggplot(
  fdis_locality_plot,
  aes(
    x = shrub_plot,
    y = mean_fdis
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
    y = "Mean functional dispersion",
    title = "Shrub and scrub"
  ) +
  theme_fdis_main()


# ------------------------------------------------------------
# 47. Panel D: Built cover
# ------------------------------------------------------------

fig_fdis_built <- ggplot(
  fdis_locality_plot,
  aes(
    x = built_plot,
    y = mean_fdis
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
    y = "Mean functional dispersion",
    title = "Built"
  ) +
  theme_fdis_main()


# ------------------------------------------------------------
# 48. Combine four empirical panels
# ------------------------------------------------------------

figure_04_main <- ggpubr::ggarrange(
  fig_fdis_aridity,
  fig_fdis_water,
  fig_fdis_shrub,
  fig_fdis_built,
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
  figure_04_main
)


# ------------------------------------------------------------
# 49. Save main figure
# ------------------------------------------------------------

ggsave(
  filename =
    file.path(
      FIGURE_DIR,
      "Figure_04_Main_Data_Trends.png"
    ),
  plot =
    figure_04_main,
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
      "Figure_04_Main_Data_Trends.pdf"
    ),
  plot =
    figure_04_main,
  width = 10,
  height = 8,
  units = "in",
  bg = "white"
)


# ============================================================
# MAIN STATISTICAL TABLE
# ============================================================


# ------------------------------------------------------------
# 50. Focal ecological term labels
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
# 51. P-value formatting helper
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
# 52. Build main model-results table
# ------------------------------------------------------------

table_04_main <- fdis_smooth_results %>%
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
  table_04_main
) <- c(
  "Predictor",
  "EDF",
  "Ref. df",
  "F statistic",
  "P-value"
)


if (
  nrow(
    table_04_main
  ) !=
  length(
    focal_term_labels
  )
) {
  
  stop(
    "Table 04 does not contain all seven focal ecological terms."
  )
}


cat("\n============================================================\n")
cat("TABLE 04 - MAIN MODEL RESULTS\n")
cat("============================================================\n\n")


print(
  table_04_main,
  row.names = FALSE
)


# ------------------------------------------------------------
# 53. Save main statistical table
# ------------------------------------------------------------

write.csv(
  table_04_main,
  file =
    file.path(
      TABLE_DIR,
      "Table_04_Main_Model_Results.csv"
    ),
  row.names = FALSE
)


# ============================================================
# MANUSCRIPT-SYNTHESIS EXPORT
# ============================================================


# ------------------------------------------------------------
# 54. Build standardized manuscript-input bundle
# ------------------------------------------------------------

fdis_model_metadata <- data.frame(
  response_id =
    "Functional_Dispersion",
  
  response_label =
    "Functional dispersion",
  
  response_variable =
    "FDis",
  
  model_family =
    "Gaussian",
  
  link =
    "identity",
  
  analysis_subset =
    "Multispecies assemblages (sr >= 2)",
  
  n_cities =
    n_distinct(
      fdis_data$NAME10
    ),
  
  n_localities =
    n_distinct(
      fdis_data$locality_id
    ),
  
  n_checklists =
    nrow(
      fdis_data
    ),
  
  adjusted_R2 =
    fdis_summary$r.sq,
  
  deviance_explained =
    fdis_summary$dev.expl,
  
  family_parameter =
    fdis_sigma,
  
  family_parameter_name =
    "residual_sd",
  
  additive_AIC =
    AIC(
      fdis_model_additive
    ),
  
  final_AIC =
    AIC(
      fdis_model_final
    ),
  
  stringsAsFactors = FALSE
)


functional_dispersion_manuscript_input <- list(
  response_id =
    "Functional_Dispersion",
  
  response_label =
    "Functional dispersion",
  
  response_variable =
    "FDis",
  
  analysis_data =
    fdis_data,
  
  final_model =
    fdis_model_final,
  
  smooth_results =
    fdis_smooth_results,
  
  parametric_results =
    fdis_parametric_results,
  
  focal_results =
    table_04_main,
  
  model_metadata =
    fdis_model_metadata,
  
  plot_data = list(
    city =
      fdis_city_plot,
    
    locality =
      fdis_locality_plot
  )
)


functional_dispersion_manuscript_file <- file.path(
  DATA_DIR,
  "Functional_Dispersion_Manuscript_Input.rds"
)


saveRDS(
  functional_dispersion_manuscript_input,
  functional_dispersion_manuscript_file
)


if (
  !file.exists(
    functional_dispersion_manuscript_file
  )
) {
  
  stop(
    "Functional-dispersion manuscript-input RDS was not created."
  )
}


cat("\n============================================================\n")
cat("MANUSCRIPT-SYNTHESIS EXPORT\n")
cat("============================================================\n\n")


cat(
  "Saved manuscript input:\n",
  functional_dispersion_manuscript_file,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 55. Automated QC summary
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
# 56. Output check
# ------------------------------------------------------------

expected_output_files <- c(
  file.path(
    FIGURE_DIR,
    "Figure_04_Main_Data_Trends.png"
  ),
  
  file.path(
    FIGURE_DIR,
    "Figure_04_Main_Data_Trends.pdf"
  ),
  
  file.path(
    TABLE_DIR,
    "Table_04_Main_Model_Results.csv"
  ),
  
  file.path(
    DATA_DIR,
    "Functional_Dispersion_Manuscript_Input.rds"
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
    "One or more expected Script 04 outputs were not created."
  )
}


cat(
  "\nPASS: All Script 04 outputs were created.\n"
)


# ------------------------------------------------------------
# 57. Final analysis summary
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
  "Residual SD:                    ",
  round(
    fdis_sigma,
    4
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
  "Negative fitted values:         ",
  sum(
    fdis_fitted_values < 0
  ),
  "\n",
  sep = ""
)


cat(
  "\nPublication figure:\n",
  file.path(
    FIGURE_DIR,
    "Figure_04_Main_Data_Trends.png"
  ),
  "\n",
  sep = ""
)


cat(
  "\nPublication table:\n",
  file.path(
    TABLE_DIR,
    "Table_04_Main_Model_Results.csv"
  ),
  "\n",
  sep = ""
)


cat(
  "\nManuscript input:\n",
  functional_dispersion_manuscript_file,
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 58. Overall interpretation
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
# Figure interpretation:
#
# - Panel A shows the city-level empirical aridity-FDis
#   relationship.
#
# - Panels B-D show locality-level habitat-FDis relationships,
#   with point color representing percentile rank along the
#   continuous city-level aridity gradient.
#
# - Warm colors indicate more-arid cities and cool colors
#   indicate less-arid cities. Percentile rank is used only for
#   visualization; the BAM uses the original continuous aridity
#   index.
#
# - The descriptive GAM lines are not the source of inferential
#   tests; inference comes from the full Gaussian BAM summarized
#   in Table 04.
#
# Overall, functional dispersion among multispecies bird
# assemblages appears to respond to both the regional aridity
# gradient and local habitat structure. However, unlike the
# species-richness response, there is comparatively little
# evidence that the effects of local habitat on functional
# dispersion change strongly across the aridity gradient.
#
# ============================================================