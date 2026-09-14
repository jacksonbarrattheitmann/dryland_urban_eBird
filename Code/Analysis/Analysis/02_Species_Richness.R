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
#
#   1. Retain Stationary and Traveling checklists
#   2. Duration <= 300 minutes
#   3. Distance <= 5 km
#   4. Require >= 5 checklists per locality
#   5. Require >= 5 eligible localities per city
#   6. Randomly sample 5 localities per city
#   7. Randomly sample 5 checklists per selected locality
#
# Model hierarchy:
#
#   City level:
#     avg_aridity_index
#
#   Locality level:
#     water
#     shrub_and_scrub
#     built
#
#   Checklist level:
#     species richness
#     duration
#     distance
#     protocol
#
# Random-effect smooths:
#
#   s(NAME10, bs = "re")
#   s(locality_id, bs = "re")
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


# Separate fixed seeds for the two sampling stages

LOCALITY_SAMPLE_SEED <- 20260911

CHECKLIST_SAMPLE_SEED <- 20260912


# ------------------------------------------------------------
# 3. Temporary duplicate handling
# ------------------------------------------------------------
#
# Jackson is currently rebuilding combined_data.RDS after
# correcting duplicated GEE locality records.
#
# TRUE:
#   temporarily retain one row per checklist for exploratory
#   analysis.
#
# FALSE:
#   stop if duplicate checklist IDs are detected.
#
# FINAL analysis should use FALSE once the corrected RDS is
# available.
# ------------------------------------------------------------

TEMPORARY_DEDUPLICATE <- TRUE


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
    "Required variables are missing:\n",
    paste(
      missing_variables,
      collapse = "\n"
    )
  )
}


cat(
  "PASS: All required analysis variables are present.\n"
)


# ------------------------------------------------------------
# 6. Duplicate checklist check
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
    paste0(
      "\nDuplicate checklist IDs are present.\n",
      "Use Jackson's corrected combined_data.RDS before ",
      "running the final analysis."
    )
  )
}


if (
  nrow(duplicate_checklists) > 0 &&
  TEMPORARY_DEDUPLICATE
) {
  
  cat(
    "\nWARNING:\n",
    "Temporarily retaining one row per checklist_id.\n",
    "Results are exploratory until the corrected RDS arrives.\n",
    sep = ""
  )
  
  
  analysis_base <- combined_data %>%
    distinct(
      checklist_id,
      .keep_all = TRUE
    )
  
  
} else {
  
  analysis_base <- combined_data
}


# ------------------------------------------------------------
# 7. Apply checklist protocol and effort filters
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
  variable =
    model_variables,
  
  missing_n =
    sapply(
      sr_pool[
        model_variables
      ],
      function(x) {
        sum(
          is.na(x)
        )
      }
    )
)


missing_summary$missing_percent <-
  100 *
  missing_summary$missing_n /
  nrow(sr_pool)


cat("\n============================================================\n")
cat("MODEL-VARIABLE MISSINGNESS\n")
cat("============================================================\n\n")


print(
  missing_summary
)


sr_pool <- sr_pool %>%
  filter(
    if_all(
      all_of(
        model_variables
      ),
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
    n_checklists >=
      MIN_CHECKLISTS_PER_LOCALITY
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
    n_eligible_localities >=
      MIN_LOCALITIES_PER_CITY
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
# 13. Randomly select five checklists per locality
# ------------------------------------------------------------

set.seed(
  CHECKLIST_SAMPLE_SEED
)


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
  group_by(
    NAME10
  ) %>%
  summarise(
    n_checklists = n(),
    
    n_localities =
      n_distinct(
        locality_id
      ),
    
    .groups = "drop"
  )


cat("\n============================================================\n")
cat("BALANCED SAMPLE QC\n")
cat("============================================================\n\n")


cat(
  "Sampled cities:      ",
  n_distinct(
    sr_data$NAME10
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled localities:  ",
  n_distinct(
    sr_data$locality_id
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled checklists:  ",
  nrow(
    sr_data
  ),
  "\n",
  sep = ""
)


cat(
  "Expected checklists: ",
  N_LOCALITIES_PER_CITY *
    N_CHECKLISTS_PER_LOCALITY *
    nrow(
      eligible_cities
    ),
  "\n",
  sep = ""
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
  
  stop(
    "Balanced sampling QC failed."
  )
}


cat(
  "\nPASS: Balanced sampling design verified.\n"
)


# ------------------------------------------------------------
# 15. Prepare grouping variables
# ------------------------------------------------------------

sr_data <- sr_data %>%
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
        protocol_name
      )
  )


# ------------------------------------------------------------
# 16. Examine richness distribution
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SPECIES RICHNESS DISTRIBUTION\n")
cat("============================================================\n\n")


print(
  summary(
    sr_data$sr
  )
)


sr_mean <- mean(
  sr_data$sr
)


sr_variance <- var(
  sr_data$sr
)


sr_var_mean <-
  sr_variance /
  sr_mean


cat(
  "\nMean richness: ",
  round(
    sr_mean,
    3
  ),
  "\n",
  sep = ""
)


cat(
  "Variance:      ",
  round(
    sr_variance,
    3
  ),
  "\n",
  sep = ""
)


cat(
  "Variance/mean: ",
  round(
    sr_var_mean,
    3
  ),
  "\n",
  sep = ""
)


ggplot(
  sr_data,
  aes(
    x = sr
  )
) +
  geom_histogram(
    binwidth = 1
  ) +
  labs(
    x = "Species richness",
    y = "Number of checklists"
  ) +
  theme_classic()


# ============================================================
# GAM / BAM MODELS
# ============================================================


# ------------------------------------------------------------
# 17. Model formula
# ------------------------------------------------------------
#
# ti() terms represent the smooth equivalents of the
# hypothesized aridity x habitat interactions.
#
# Main smooths are included separately so ti() represents the
# additional interaction beyond the marginal relationships.
#
# k values are kept deliberately modest initially.
# gam.check() below will tell us whether any basis dimensions
# need to be increased.
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
# 18. Fit Poisson BAM
# ------------------------------------------------------------
#
# This model is used as the count-family baseline.
#
# We use BAM from the start so the analytical framework is the
# same one we intend to use when testing the full dataset.
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FITTING POISSON BAM\n")
cat("============================================================\n\n")


sr_model_poisson <- bam(
  formula =
    sr_formula,
  
  family =
    poisson(
      link = "log"
    ),
  
  data =
    sr_data,
  
  method =
    "fREML",
  
  discrete =
    TRUE
)


# ------------------------------------------------------------
# 19. Poisson model summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("POISSON BAM SUMMARY\n")
cat("============================================================\n\n")


print(
  summary(
    sr_model_poisson
  )
)


# ------------------------------------------------------------
# 20. Poisson dispersion
# ------------------------------------------------------------

poisson_dispersion <- sum(
  residuals(
    sr_model_poisson,
    type = "pearson"
  )^2
) /
  df.residual(
    sr_model_poisson
  )


cat(
  "\nPoisson Pearson dispersion ratio: ",
  round(
    poisson_dispersion,
    3
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 21. Fit negative-binomial BAM
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FITTING NEGATIVE-BINOMIAL BAM\n")
cat("============================================================\n\n")


sr_model_nb <- bam(
  formula =
    sr_formula,
  
  family =
    nb(
      link = "log"
    ),
  
  data =
    sr_data,
  
  method =
    "fREML",
  
  discrete =
    TRUE
)


# ------------------------------------------------------------
# 22. Negative-binomial model summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("NEGATIVE-BINOMIAL BAM SUMMARY\n")
cat("============================================================\n\n")


print(
  summary(
    sr_model_nb
  )
)


# ------------------------------------------------------------
# 23. Negative-binomial dispersion
# ------------------------------------------------------------

nb_dispersion <- sum(
  residuals(
    sr_model_nb,
    type = "pearson"
  )^2
) /
  df.residual(
    sr_model_nb
  )


cat(
  "\nNegative-binomial Pearson dispersion ratio: ",
  round(
    nb_dispersion,
    3
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 24. Compare count families
# ------------------------------------------------------------

model_comparison <- data.frame(
  model = c(
    "Poisson BAM",
    "Negative-binomial BAM"
  ),
  
  AIC = c(
    AIC(
      sr_model_poisson
    ),
    
    AIC(
      sr_model_nb
    )
  ),
  
  Pearson_dispersion = c(
    poisson_dispersion,
    nb_dispersion
  )
)


cat("\n============================================================\n")
cat("COUNT-FAMILY COMPARISON\n")
cat("============================================================\n\n")


print(
  model_comparison
)


# ------------------------------------------------------------
# 25. Select final richness model
# ------------------------------------------------------------
#
# We expect the negative-binomial family to be more appropriate
# for overdispersed richness data.
#
# The Poisson comparison above provides the empirical check.
# ------------------------------------------------------------

sr_model_final <- sr_model_nb


# ------------------------------------------------------------
# 26. Final model summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL SPECIES RICHNESS BAM\n")
cat("============================================================\n\n")


sr_summary <- summary(
  sr_model_final
)


print(
  sr_summary
)


# ------------------------------------------------------------
# 27. GAM basis-dimension and residual diagnostics
# ------------------------------------------------------------
#
# Pay particular attention to:
#
#   k-index
#   p-values associated with k-index
#
# A low k-index together with a small p-value can indicate that
# a basis dimension should be increased.
#
# Random-effect smooths are not interpreted using the same
# k-index logic as ordinary smooths.
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("GAM CHECK\n")
cat("============================================================\n\n")


gam.check(
  sr_model_final
)


# ------------------------------------------------------------
# 28. Concurvity
# ------------------------------------------------------------
#
# Concurvity is the GAM analogue of multicollinearity.
#
# Values approaching 1 indicate that one model component can be
# approximated strongly by combinations of other components.
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CONCURVITY CHECK\n")
cat("============================================================\n\n")


sr_concurvity <- concurvity(
  sr_model_final,
  full = TRUE
)


print(
  sr_concurvity
)


# ------------------------------------------------------------
# 29. Extract parametric results
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


cat("\n============================================================\n")
cat("PARAMETRIC TERMS\n")
cat("============================================================\n\n")


print(
  sr_parametric_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 30. Extract smooth-term results
# ------------------------------------------------------------

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
cat("SMOOTH TERMS\n")
cat("============================================================\n\n")


print(
  sr_smooth_results,
  row.names = FALSE
)


# ------------------------------------------------------------
# 31. Model fit summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("MODEL FIT\n")
cat("============================================================\n\n")


cat(
  "Adjusted R-squared: ",
  round(
    sr_summary$r.sq,
    3
  ),
  "\n",
  sep = ""
)


cat(
  "Deviance explained: ",
  round(
    100 *
      sr_summary$dev.expl,
    2
  ),
  "%\n",
  sep = ""
)


cat(
  "Negative-binomial theta: ",
  round(
    sr_model_final$family$getTheta(
      TRUE
    ),
    3
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 32. Plot non-random smooth terms
# ------------------------------------------------------------
#
# Smooth terms 1-9 represent:
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
# Random-effect smooths are excluded from these plots.
# ------------------------------------------------------------

plot(
  sr_model_final,
  select = 1,
  shade = TRUE,
  pages = 1
)


plot(
  sr_model_final,
  select = 2,
  shade = TRUE,
  pages = 1
)


plot(
  sr_model_final,
  select = 3,
  shade = TRUE,
  pages = 1
)


plot(
  sr_model_final,
  select = 4,
  shade = TRUE,
  pages = 1
)


plot(
  sr_model_final,
  select = 5,
  pages = 1
)


plot(
  sr_model_final,
  select = 6,
  pages = 1
)


plot(
  sr_model_final,
  select = 7,
  pages = 1
)


plot(
  sr_model_final,
  select = 8,
  shade = TRUE,
  pages = 1
)


plot(
  sr_model_final,
  select = 9,
  shade = TRUE,
  pages = 1
)


# ------------------------------------------------------------
# 33. Descriptive aridity plot
# ------------------------------------------------------------

sr_city_plot <- sr_data %>%
  group_by(
    NAME10
  ) %>%
  summarise(
    avg_aridity_index =
      first(
        avg_aridity_index
      ),
    
    mean_sr =
      mean(
        sr
      ),
    
    .groups = "drop"
  )


ggplot(
  sr_city_plot,
  aes(
    x = avg_aridity_index,
    y = mean_sr
  )
) +
  geom_point(
    alpha = 0.5
  ) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = TRUE
  ) +
  labs(
    x = "Aridity index (higher = less arid)",
    y = "Mean species richness"
  ) +
  theme_classic()


# ------------------------------------------------------------
# 34. Descriptive locality-level habitat plots
# ------------------------------------------------------------

sr_locality_plot <- sr_data %>%
  group_by(
    NAME10,
    locality_id
  ) %>%
  summarise(
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
    
    mean_sr =
      mean(
        sr
      ),
    
    .groups = "drop"
  )


# Water

ggplot(
  sr_locality_plot,
  aes(
    x = water,
    y = mean_sr
  )
) +
  geom_point(
    alpha = 0.3
  ) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = TRUE
  ) +
  scale_x_continuous(
    labels = scales::percent
  ) +
  labs(
    x = "Water cover",
    y = "Mean species richness"
  ) +
  theme_classic()


# Shrub and scrub

ggplot(
  sr_locality_plot,
  aes(
    x = shrub_and_scrub,
    y = mean_sr
  )
) +
  geom_point(
    alpha = 0.3
  ) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = TRUE
  ) +
  scale_x_continuous(
    labels = scales::percent
  ) +
  labs(
    x = "Shrub and scrub cover",
    y = "Mean species richness"
  ) +
  theme_classic()


# Built

ggplot(
  sr_locality_plot,
  aes(
    x = built,
    y = mean_sr
  )
) +
  geom_point(
    alpha = 0.3
  ) +
  geom_smooth(
    method = "gam",
    formula = y ~ s(x),
    se = TRUE
  ) +
  scale_x_continuous(
    labels = scales::percent
  ) +
  labs(
    x = "Built cover",
    y = "Mean species richness"
  ) +
  theme_classic()


# ------------------------------------------------------------
# 35. Final analysis summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SPECIES RICHNESS GAM ANALYSIS COMPLETE\n")
cat("============================================================\n\n")


cat(
  "Eligible cities:               ",
  nrow(
    eligible_cities
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled cities:                ",
  n_distinct(
    sr_data$NAME10
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled localities:            ",
  n_distinct(
    sr_data$locality_id
  ),
  "\n",
  sep = ""
)


cat(
  "Sampled checklists:            ",
  nrow(
    sr_data
  ),
  "\n",
  sep = ""
)


cat(
  "Checklists per city:           ",
  N_LOCALITIES_PER_CITY *
    N_CHECKLISTS_PER_LOCALITY,
  "\n",
  sep = ""
)


cat(
  "Poisson dispersion:            ",
  round(
    poisson_dispersion,
    3
  ),
  "\n",
  sep = ""
)


cat(
  "Negative-binomial dispersion:  ",
  round(
    nb_dispersion,
    3
  ),
  "\n",
  sep = ""
)


cat(
  "Poisson AIC:                   ",
  round(
    AIC(
      sr_model_poisson
    ),
    1
  ),
  "\n",
  sep = ""
)


cat(
  "Negative-binomial AIC:         ",
  round(
    AIC(
      sr_model_nb
    ),
    1
  ),
  "\n",
  sep = ""
)


cat(
  "Deviance explained:            ",
  round(
    100 *
      sr_summary$dev.expl,
    2
  ),
  "%\n",
  sep = ""
)


cat(
  "Adjusted R-squared:            ",
  round(
    sr_summary$r.sq,
    3
  ),
  "\n",
  sep = ""
)








# ============================================================
# OVERALL MODEL INTERPRETATION
# ============================================================
#
# Model performance:
#
# - The balanced sampling design retained 273 cities,
#   1,365 localities, and 6,825 checklists.
#
# - The negative-binomial BAM was clearly preferred over the
#   Poisson BAM. The NB model had substantially lower AIC and
#   reduced Pearson dispersion from 1.53 to 1.16.
#
# - The final negative-binomial BAM converged successfully,
#   was full rank, and explained approximately 69.5% of the
#   deviance (adjusted R-squared = 0.649).
#
# - GAM basis-dimension checks indicated that the selected k
#   values were sufficient for the fitted smooths.
#
# - Residual diagnostics were broadly acceptable, with 4.47%
#   of Pearson residuals exceeding an absolute value of 2 and
#   only 0.62% exceeding an absolute value of 3.
#
#
# Main ecological patterns:
#
# - Aridity alone showed relatively weak evidence of a direct
#   relationship with species richness. The aridity smooth was
#   essentially linear (edf ~ 1) and was not statistically
#   significant at alpha = 0.05. Funny enough, a lot of my 
#   work has actually showed the same thing.
#
# - Local habitat showed stronger relationships with species
#   richness. Water, shrub/scrub, and built cover each showed
#   significant nonlinear associations with richness. Again,
#   some of my old work has show that shrub/scrib presence
#   alone can do this for Vertebrate species. AWESOME!
#
# - The relationship between habitat and species richness also
#   changed across the aridity gradient. AGAIN, a lot of my
#   old work showed that aridity needs to be interacting
#   with another variable and then shows significant.
#   My 2025 paper with fake shrubs showed this aswell.
#
# - The aridity x water interaction was significant, indicating
#   that the relationship between water cover and species
#   richness differs across wetter versus more arid cities.
#   This is an awesome finding.
#
# - The aridity x shrub/scrub interaction was also significant,
#   indicating that the relationship between shrub cover and
#   species richness changes along the aridity gradient.
#   Again, my old work has shown this at a local scale.
#
# - In contrast, there was no evidence for an aridity x built
#   cover interaction, suggesting that the relationship between
#   built cover and richness is comparatively consistent across
#   the aridity gradient. This to me make sense.
#
#
# Checklist effort:
#
# - Checklist duration was strongly associated with observed
#   species richness and showed a clearly nonlinear relationship.
#
# - Distance traveled was also significantly associated with
#   richness.
#
# - These results reinforce the importance of controlling for
#   eBird sampling effort when estimating ecological effects.
#
#
# Hierarchical structure:
#
# - Both city and locality contributed meaningful hierarchical
#   structure to the model.
#
# - The estimated random-effect standard deviation was
#   approximately 0.18 for city and 0.36 for locality.
#
# - The stronger locality-level variation supports retaining
#   locality_id explicitly in the BAM rather than accounting
#   for city alone.
#
#
# Concurvity and model structure:
#
# - Pairwise concurvity among most ecological predictors was
#   low to moderate.
#
# - Higher concurvity involving aridity and the aridity x
#   habitat interactions is expected because the interaction
#   surfaces explicitly contain the aridity gradient.
#
# - No diagnostic evidence indicated that the current model
#   structure required simplification or removal of the
#   locality-level random effect.
#
#
# Overall interpretation:
#
# Species richness appears to be structured primarily by local
# habitat and by the way habitat relationships change across
# the regional aridity gradient, rather than by a strong,
# independent effect of aridity alone.
#
# Water and shrub/scrub relationships with richness vary along
# the aridity gradient, whereas the effect of built cover appears
# more consistent across the gradient.
#
# Together, these results suggest that regional aridity acts in
# part as a modifier of local habitat-richness relationships,
# while local habitat and checklist effort remain important
# determinants of observed bird species richness.
# ============================================================