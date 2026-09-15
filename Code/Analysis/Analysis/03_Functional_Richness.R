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
#   2. Stationary checklist distance = 0 km
#   3. Duration <= 300 minutes
#   4. Distance <= 5 km
#   5. Retain checklists with defined FRic
#   6. Require >= 5 valid FRic checklists per locality
#   7. Require >= 5 eligible localities per city
#   8. Randomly sample 5 localities per city
#   9. Randomly sample 5 checklists per selected locality
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
# Final response family:
#
#   Beta regression with logit link
#
# Publication / meeting outputs:
#
#   Outputs/
#     03_Functional_Richness/
#       Figure/
#         Figure_03_Main_Data_Trends.png
#         Figure_03_Main_Data_Trends.pdf
#
#       Table/
#         Table_03_Main_Model_Results.csv
#
#       Data/
#         Functional_Richness_Manuscript_Input.rds
#
# Main figure:
#
#   A. Aridity vs mean functional richness
#   B. Water cover vs mean functional richness
#   C. Shrub/scrub cover vs mean functional richness
#   D. Built cover vs mean functional richness
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
# data. The statistical models continue to use the original
# continuous aridity index.
#
# Figure smooths are descriptive and unadjusted.
# Statistical inference comes from the full beta BAM and is
# summarized in Table 03.
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
  "03_Functional_Richness"
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


if (
  !file.exists(
    data_path
  )
) {
  
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
  nrow(
    combined_data
  ),
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
  "FRic",
  "avg_aridity_index",
  "water",
  "shrub_and_scrub",
  "built"
)


missing_variables <- setdiff(
  required_variables,
  names(
    combined_data
  )
)


if (
  length(
    missing_variables
  ) > 0
) {
  
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
  nrow(
    duplicate_checklists
  ),
  "\n",
  sep = ""
)


if (
  nrow(
    duplicate_checklists
  ) > 0 &&
  !TEMPORARY_DEDUPLICATE
) {
  
  stop(
    "Duplicate checklist IDs are present."
  )
}


if (
  nrow(
    duplicate_checklists
  ) > 0 &&
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
  nrow(
    analysis_base
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 7. Apply protocol and effort filters
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
# 8. Examine FRic before removing missing values
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
# 9. Remove undefined FRic values
# ------------------------------------------------------------

fric_pool <- fric_pool %>%
  filter(
    !is.na(
      FRic
    )
  )


cat(
  "\nRows with defined FRic: ",
  nrow(
    fric_pool
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 10. Complete-case restriction
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
          is.na(
            x
          )
        )
      }
    )
)


missing_summary$missing_percent <-
  100 *
  missing_summary$missing_n /
  nrow(
    fric_pool
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
# 11. Identify eligible localities
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
cat("FRIC ELIGIBILITY\n")
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
# 13. Create eligible sampling pool
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
    n =
      N_LOCALITIES_PER_CITY
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 15. Sample five FRic checklists per locality
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
    n =
      N_CHECKLISTS_PER_LOCALITY
  ) %>%
  ungroup()


# ------------------------------------------------------------
# 16. Balanced sample QC
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


bad_localities <- sum(
  locality_qc$n_checklists !=
    N_CHECKLISTS_PER_LOCALITY
)


bad_cities_localities <- sum(
  city_qc$n_localities !=
    N_LOCALITIES_PER_CITY
)


bad_cities_checklists <- sum(
  city_qc$n_checklists !=
    (
      N_LOCALITIES_PER_CITY *
        N_CHECKLISTS_PER_LOCALITY
    )
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
  bad_localities,
  "\n",
  sep = ""
)


cat(
  "Cities with != 5 localities: ",
  bad_cities_localities,
  "\n",
  sep = ""
)


cat(
  "Cities with != 25 checklists: ",
  bad_cities_checklists,
  "\n",
  sep = ""
)


if (
  bad_localities > 0 ||
  bad_cities_localities > 0 ||
  bad_cities_checklists > 0
) {
  
  stop(
    "FRic balanced sampling QC failed."
  )
}


cat(
  "\nPASS: FRic balanced sampling verified.\n"
)


# ------------------------------------------------------------
# 17. Prepare factors
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
# 18. Examine final FRic response distribution
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
# 19. Confirm beta-family requirements
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
# 20. Additive model formula
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
# 21. Interaction model formula
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
# 22. Fit additive beta BAM
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FITTING ADDITIVE FRIC BAM\n")
cat("============================================================\n\n")


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


# ------------------------------------------------------------
# 23. Fit interaction beta BAM
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FITTING INTERACTION FRIC BAM\n")
cat("============================================================\n\n")


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


# ------------------------------------------------------------
# 24. Candidate-model summaries
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
# 25. Compare candidate models
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
# 26. Candidate-model convergence
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
# 29. Stable basis-dimension check
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FRIC STABLE BASIS-DIMENSION CHECK\n")
cat("============================================================\n\n")


set.seed(
  K_CHECK_SEED
)


fric_k_check <- k.check(
  fric_model_final,
  n.rep = 1000
)


print(
  fric_k_check
)


k_valid_rows <- !is.na(
  fric_k_check[
    ,
    "k-index"
  ]
)


k_followup_flag <- any(
  fric_k_check[
    k_valid_rows,
    "k-index"
  ] < 1 &
    
    fric_k_check[
      k_valid_rows,
      "p-value"
    ] < 0.05 &
    
    (
      fric_k_check[
        k_valid_rows,
        "edf"
      ] /
        fric_k_check[
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
# 34. Residual diagnostics
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


cat(
  "\nAbsolute Pearson residual > 2: ",
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


plot(
  fric_fitted_values,
  fric_pearson_residuals,
  pch = 16,
  cex = 0.5,
  xlab = "Fitted functional richness",
  ylab = "Pearson residual"
)


abline(
  h = 0,
  lty = 2
)


qq.gam(
  fric_model_final,
  rep = 100,
  level = 0.90
)


set.seed(
  K_CHECK_SEED
)


gam.check(
  fric_model_final
)


# ------------------------------------------------------------
# 35. Extract model results
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
cat("FINAL FRIC PARAMETRIC TERMS\n")
cat("============================================================\n\n")


print(
  fric_parametric_results,
  row.names = FALSE
)


cat("\n============================================================\n")
cat("FINAL FRIC SMOOTH TERMS\n")
cat("============================================================\n\n")


print(
  fric_smooth_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 36. Final model fit
# ------------------------------------------------------------

fric_precision <- fric_model_final$family$getTheta(
  TRUE
)


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
    fric_precision,
    3
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
#
# Statistical inference is based on the full beta BAM and is
# reported in Table 03.
#
# ============================================================


# ------------------------------------------------------------
# 37. Create empirical trend datasets
# ------------------------------------------------------------

fric_city_plot <- fric_data %>%
  group_by(
    NAME10
  ) %>%
  summarise(
    avg_aridity_index =
      first(
        avg_aridity_index
      ),
    
    mean_fric =
      mean(
        FRic,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


fric_locality_plot <- fric_data %>%
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
    
    mean_fric =
      mean(
        FRic,
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
# 38. Convert habitat cover to percent
# ------------------------------------------------------------

maximum_habitat_value <- max(
  c(
    fric_locality_plot$water,
    fric_locality_plot$shrub_and_scrub,
    fric_locality_plot$built
  ),
  na.rm = TRUE
)


if (
  maximum_habitat_value <= 1.05
) {
  
  fric_locality_plot <- fric_locality_plot %>%
    mutate(
      water_plot =
        water * 100,
      
      shrub_plot =
        shrub_and_scrub * 100,
      
      built_plot =
        built * 100
    )
  
} else {
  
  fric_locality_plot <- fric_locality_plot %>%
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
# 39. Common y-axis range
# ------------------------------------------------------------

main_y_limits <- range(
  c(
    fric_city_plot$mean_fric,
    fric_locality_plot$mean_fric
  ),
  na.rm = TRUE
)


# ------------------------------------------------------------
# 40. Common plotting theme and aridity color scale
# ------------------------------------------------------------

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


theme_fric_main <- function() {
  
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
# 41. Panel A: Aridity
# ------------------------------------------------------------

fig_fric_aridity <- ggplot(
  fric_city_plot,
  aes(
    x = avg_aridity_index,
    y = mean_fric
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
    y = "Mean functional richness",
    title = "Aridity"
  ) +
  theme_fric_main() +
  theme(
    legend.position = "none"
  )


# ------------------------------------------------------------
# 42. Panel B: Water cover
# ------------------------------------------------------------

fig_fric_water <- ggplot(
  fric_locality_plot,
  aes(
    x = water_plot,
    y = mean_fric
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
    y = "Mean functional richness",
    title = "Water"
  ) +
  theme_fric_main()


# ------------------------------------------------------------
# 43. Panel C: Shrub and scrub cover
# ------------------------------------------------------------

fig_fric_shrub <- ggplot(
  fric_locality_plot,
  aes(
    x = shrub_plot,
    y = mean_fric
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
    y = "Mean functional richness",
    title = "Shrub and scrub"
  ) +
  theme_fric_main()


# ------------------------------------------------------------
# 44. Panel D: Built cover
# ------------------------------------------------------------

fig_fric_built <- ggplot(
  fric_locality_plot,
  aes(
    x = built_plot,
    y = mean_fric
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
    y = "Mean functional richness",
    title = "Built"
  ) +
  theme_fric_main()


# ------------------------------------------------------------
# 45. Combine four empirical panels
# ------------------------------------------------------------

figure_03_main <- ggpubr::ggarrange(
  fig_fric_aridity,
  fig_fric_water,
  fig_fric_shrub,
  fig_fric_built,
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
  figure_03_main
)


# ------------------------------------------------------------
# 46. Save main figure
# ------------------------------------------------------------

ggsave(
  filename =
    file.path(
      FIGURE_DIR,
      "Figure_03_Main_Data_Trends.png"
    ),
  plot =
    figure_03_main,
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
      "Figure_03_Main_Data_Trends.pdf"
    ),
  plot =
    figure_03_main,
  width = 10,
  height = 8,
  units = "in",
  bg = "white"
)


# ============================================================
# MAIN STATISTICAL TABLE
# ============================================================
#
# Only focal ecological terms are exported.
#
# Effort covariates, protocol, city, and locality remain in
# the fitted BAM but are not included in the main table.
#
# ============================================================


# ------------------------------------------------------------
# 47. Focal ecological term labels
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
# 48. P-value formatting helper
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
# 49. Build main model-results table
# ------------------------------------------------------------

table_03_main <- fric_smooth_results %>%
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
  table_03_main
) <- c(
  "Predictor",
  "EDF",
  "Ref. df",
  "F statistic",
  "P-value"
)


if (
  nrow(
    table_03_main
  ) !=
  length(
    focal_term_labels
  )
) {
  
  stop(
    "Table 03 does not contain all seven focal ecological terms."
  )
}


cat("\n============================================================\n")
cat("TABLE 03 - MAIN MODEL RESULTS\n")
cat("============================================================\n\n")


print(
  table_03_main,
  row.names = FALSE
)


# ------------------------------------------------------------
# 50. Save main statistical table
# ------------------------------------------------------------

write.csv(
  table_03_main,
  file =
    file.path(
      TABLE_DIR,
      "Table_03_Main_Model_Results.csv"
    ),
  row.names = FALSE
)


# ============================================================
# MANUSCRIPT-SYNTHESIS EXPORT
# ============================================================


# ------------------------------------------------------------
# 51. Build standardized manuscript-input bundle
# ------------------------------------------------------------

fric_model_metadata <- data.frame(
  response_id =
    "Functional_Richness",
  
  response_label =
    "Functional richness",
  
  response_variable =
    "FRic",
  
  model_family =
    "Beta",
  
  link =
    "logit",
  
  n_cities =
    n_distinct(
      fric_data$NAME10
    ),
  
  n_localities =
    n_distinct(
      fric_data$locality_id
    ),
  
  n_checklists =
    nrow(
      fric_data
    ),
  
  adjusted_R2 =
    fric_summary$r.sq,
  
  deviance_explained =
    fric_summary$dev.expl,
  
  family_parameter =
    fric_precision,
  
  family_parameter_name =
    "precision",
  
  additive_AIC =
    AIC(
      fric_model_additive
    ),
  
  final_AIC =
    AIC(
      fric_model_final
    ),
  
  stringsAsFactors = FALSE
)


functional_richness_manuscript_input <- list(
  response_id =
    "Functional_Richness",
  
  response_label =
    "Functional richness",
  
  response_variable =
    "FRic",
  
  analysis_data =
    fric_data,
  
  final_model =
    fric_model_final,
  
  smooth_results =
    fric_smooth_results,
  
  parametric_results =
    fric_parametric_results,
  
  focal_results =
    table_03_main,
  
  model_metadata =
    fric_model_metadata,
  
  plot_data = list(
    city =
      fric_city_plot,
    
    locality =
      fric_locality_plot
  )
)


functional_richness_manuscript_file <- file.path(
  DATA_DIR,
  "Functional_Richness_Manuscript_Input.rds"
)


saveRDS(
  functional_richness_manuscript_input,
  functional_richness_manuscript_file
)


if (
  !file.exists(
    functional_richness_manuscript_file
  )
) {
  
  stop(
    "Functional-richness manuscript-input RDS was not created."
  )
}


cat("\n============================================================\n")
cat("MANUSCRIPT-SYNTHESIS EXPORT\n")
cat("============================================================\n\n")


cat(
  "Saved manuscript input:\n",
  functional_richness_manuscript_file,
  "\n",
  sep = ""
)


# ============================================================
# AUTOMATED QC AND OUTPUT CHECK
# ============================================================


# ------------------------------------------------------------
# 52. Automated analysis QC
# ------------------------------------------------------------

qc_balanced_sample <-
  bad_localities == 0 &&
  bad_cities_localities == 0 &&
  bad_cities_checklists == 0


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


qc_basis_dimension <-
  !k_followup_flag


fric_qc_summary <- data.frame(
  check = c(
    "Balanced 5 x 5 sampling",
    "FRic values strictly inside (0,1)",
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
# 53. Publication-output check
# ------------------------------------------------------------

expected_output_files <- c(
  file.path(
    FIGURE_DIR,
    "Figure_03_Main_Data_Trends.png"
  ),
  
  file.path(
    FIGURE_DIR,
    "Figure_03_Main_Data_Trends.pdf"
  ),
  
  file.path(
    TABLE_DIR,
    "Table_03_Main_Model_Results.csv"
  ),
  
  file.path(
    DATA_DIR,
    "Functional_Richness_Manuscript_Input.rds"
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
    "One or more expected Script 03 outputs were not created."
  )
}


cat(
  "\nPASS: All Script 03 publication outputs were created.\n"
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
    fric_precision,
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
    "Figure_03_Main_Data_Trends.png"
  ),
  "\n",
  sep = ""
)


cat(
  "\nPublication table:\n",
  file.path(
    TABLE_DIR,
    "Table_03_Main_Model_Results.csv"
  ),
  "\n",
  sep = ""
)


cat(
  "\nManuscript input:\n",
  functional_richness_manuscript_file,
  "\n",
  sep = ""
)


# ============================================================
# OVERALL INTERPRETATION
# ============================================================
#
# Model performance:
#
# - The balanced sampling design retains 25 checklists per
#   eligible city: five checklists from each of five sampled
#   localities.
#
# - FRic is modeled with a beta BAM because retained FRic
#   values are strictly between 0 and 1.
#
# - The additive and interaction models have essentially
#   equivalent overall fit. The interaction model is retained
#   as the primary/global model because the aridity x habitat
#   relationships represent the a priori ecological hypotheses
#   and maintaining the same global structure facilitates
#   comparison among biodiversity responses.
#
#
# Ecological interpretation:
#
# - Water cover shows the clearest relationship with functional
#   richness in the final beta BAM.
#
# - There is little evidence for an independent relationship
#   between aridity and functional richness after accounting
#   for local habitat, checklist effort, and hierarchical
#   structure.
#
# - Shrub and scrub cover provides comparatively weak evidence
#   of an independent relationship with functional richness.
#
# - Built cover provides little evidence of an independent
#   relationship with functional richness.
#
# - None of the three aridity x habitat interactions receives
#   strong statistical support.
#
# - Checklist duration, distance traveled, and protocol are
#   retained as sampling-effort controls.
#
# - City and locality random-effect smooths account for
#   hierarchical variation in checklist-level functional
#   richness.
#
#
# Figure interpretation:
#
# - The exported four-panel figure displays descriptive,
#   unadjusted empirical relationships.
#
# - Panel A represents city-level mean functional richness
#   because aridity is defined at the city level.
#
# - Panels B-D represent locality-level mean functional richness
#   because habitat variables are defined at the locality level.
#
# - Point color in panels B-D represents percentile rank along
#   the continuous city-level aridity gradient. Warm colors are
#   more arid and cool colors are less arid. Percentile rank is
#   used only for visualization; the BAM uses the original
#   continuous aridity index.
#
# - The descriptive GAM lines in the figure are not the source
#   of inferential tests.
#
# - Statistical inference for the focal ecological predictors
#   and aridity x habitat interactions comes from the full beta
#   BAM summarized in Table 03.
#
#
# Overall:
#
# Functional richness appears to be associated more strongly
# with local habitat structure, particularly water cover, than
# with the regional aridity gradient itself. There is little
# evidence that local habitat-FRic relationships are strongly
# modified across the regional aridity gradient.
#
# ============================================================
