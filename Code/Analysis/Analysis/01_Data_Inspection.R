# ============================================================
# Dryland Bird Project
# 01_Data_Inspection.R
#
# Purpose:
# Inspect the final combined analysis dataset, apply the basic
# checklist-effort filters, evaluate sampling imbalance among
# localities and cities, and define the pool of observations
# eligible for downstream analyses.
#
# This script DOES NOT randomly subsample the data.
# Random subsampling will be performed in Script 02 using a
# fixed random seed.
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


# ------------------------------------------------------------
# 2. Settings
# ------------------------------------------------------------

# Jackson is currently rebuilding the combined dataset after
# removing duplicated GEE locality records.
#
# TRUE allows us to inspect the current dataset temporarily by
# retaining one row per checklist_id.
#
# When the corrected combined_data.RDS arrives, this should
# have no effect because duplicate checklist IDs should = 0.

TEMPORARY_DEDUPLICATE <- TRUE


# Sampling thresholds suggested for evaluation

MIN_CHECKLISTS_PER_LOCALITY <- 5

MIN_LOCALITIES_PER_CITY <- 5


# ------------------------------------------------------------
# 3. Load final combined dataset
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
cat("DRYLAND BIRD DATA INSPECTION\n")
cat("============================================================\n\n")


# ------------------------------------------------------------
# 4. Dataset dimensions and structure
# ------------------------------------------------------------

cat("Dataset dimensions:\n")
print(
  dim(combined_data)
)


cat("\nColumn names:\n")
print(
  names(combined_data)
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
  "number_observers",
  "sr",
  "FRic",
  "FDis",
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


cat("\n============================================================\n")
cat("REQUIRED VARIABLES\n")
cat("============================================================\n\n")


if (length(missing_variables) == 0) {
  
  cat("PASS: All required analysis variables are present.\n")
  
} else {
  
  cat("Missing variables:\n")
  print(missing_variables)
  
  stop(
    "Required variables are missing from combined_data.RDS."
  )
}


# ------------------------------------------------------------
# 6. Basic sampling structure
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("RAW DATA STRUCTURE\n")
cat("============================================================\n\n")


raw_structure <- combined_data %>%
  summarise(
    n_rows = n(),
    n_unique_checklists = n_distinct(checklist_id),
    n_localities = n_distinct(locality_id),
    n_cities = n_distinct(NAME10),
    n_states = n_distinct(state_name)
  )


print(raw_structure)


# ------------------------------------------------------------
# 7. Duplicate checklist QC
# ------------------------------------------------------------

duplicate_checklists <- combined_data %>%
  count(
    checklist_id,
    name = "n_rows"
  ) %>%
  filter(
    n_rows > 1
  )


cat("\n============================================================\n")
cat("DUPLICATE CHECKLIST QC\n")
cat("============================================================\n\n")


cat(
  "Duplicated checklist IDs: ",
  nrow(duplicate_checklists),
  "\n",
  sep = ""
)


cat(
  "Extra rows associated with duplicated checklist IDs: ",
  sum(
    duplicate_checklists$n_rows - 1
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Temporary duplicate handling
# ------------------------------------------------------------
#
# Jackson has identified and cleaned the duplicate locality
# records upstream and is rebuilding combined_data.RDS.
#
# For temporary exploratory work only, retain one row per
# checklist_id as he suggested.
#
# THE FINAL ANALYSIS MUST BE RERUN WITH HIS CLEANED RDS.
# ------------------------------------------------------------

if (
  nrow(duplicate_checklists) > 0 &&
  TEMPORARY_DEDUPLICATE
) {
  
  cat(
    "\nNOTE: Temporarily retaining one row per checklist_id.\n"
  )
  
  cat(
    "Final analyses must be rerun with the corrected RDS.\n"
  )
  
  
  analysis_base <- combined_data %>%
    distinct(
      checklist_id,
      .keep_all = TRUE
    )
  
  
} else {
  
  analysis_base <- combined_data
}


cat(
  "\nRows after duplicate handling: ",
  nrow(analysis_base),
  "\n",
  sep = ""
)


cat(
  "Unique checklist IDs after duplicate handling: ",
  n_distinct(
    analysis_base$checklist_id
  ),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 9. Examine checklist protocols
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("PROTOCOL STRUCTURE\n")
cat("============================================================\n\n")


protocol_summary <- analysis_base %>%
  count(
    protocol_name,
    sort = TRUE
  )


print(
  protocol_summary,
  n = Inf
)


# ------------------------------------------------------------
# 10. Examine effort variables by protocol
# ------------------------------------------------------------

effort_by_protocol <- analysis_base %>%
  group_by(
    protocol_name
  ) %>%
  summarise(
    n = n(),
    
    missing_duration =
      sum(
        is.na(duration_minutes)
      ),
    
    missing_distance =
      sum(
        is.na(effort_distance_km)
      ),
    
    median_duration =
      median(
        duration_minutes,
        na.rm = TRUE
      ),
    
    median_distance =
      median(
        effort_distance_km,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


print(
  effort_by_protocol,
  n = Inf
)


# ------------------------------------------------------------
# 11. Apply basic analysis filters
# ------------------------------------------------------------
#
# Retain:
#
#   Stationary
#   Traveling
#
# Stationary checklist distance is set to zero.
#
# Effort restrictions:
#
#   duration <= 300 minutes
#   distance <= 5 km
# ------------------------------------------------------------

dryland_data <- analysis_base %>%
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


# ------------------------------------------------------------
# 12. Check filtered analysis dataset
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("EFFORT-FILTERED DATA\n")
cat("============================================================\n\n")


filtered_structure <- dryland_data %>%
  summarise(
    n_checklists = n(),
    n_unique_checklists = n_distinct(checklist_id),
    n_localities = n_distinct(locality_id),
    n_cities = n_distinct(NAME10),
    n_states = n_distinct(state_name),
    
    max_duration =
      max(
        duration_minutes,
        na.rm = TRUE
      ),
    
    max_distance =
      max(
        effort_distance_km,
        na.rm = TRUE
      )
  )


print(filtered_structure)


cat("\nProtocol counts after filtering:\n")


dryland_data %>%
  count(
    protocol_name,
    sort = TRUE
  ) %>%
  print()


# ------------------------------------------------------------
# 13. Quantify filtering
# ------------------------------------------------------------

filter_summary <- data.frame(
  stage = c(
    "Raw combined data",
    "After duplicate handling",
    "After protocol/effort filters"
  ),
  
  n_rows = c(
    nrow(combined_data),
    nrow(analysis_base),
    nrow(dryland_data)
  ),
  
  n_localities = c(
    n_distinct(combined_data$locality_id),
    n_distinct(analysis_base$locality_id),
    n_distinct(dryland_data$locality_id)
  ),
  
  n_cities = c(
    n_distinct(combined_data$NAME10),
    n_distinct(analysis_base$NAME10),
    n_distinct(dryland_data$NAME10)
  )
)


cat("\nFiltering summary:\n")
print(filter_summary)


# ------------------------------------------------------------
# 14. Missingness in analysis variables
# ------------------------------------------------------------

analysis_variables <- c(
  "sr",
  "FRic",
  "FDis",
  "FEve",
  "avg_aridity_index",
  "water",
  "shrub_and_scrub",
  "built",
  "duration_minutes",
  "effort_distance_km",
  "protocol_name",
  "locality_id",
  "NAME10"
)


missing_data <- data.frame(
  variable =
    analysis_variables,
  
  missing_n =
    sapply(
      dryland_data[
        analysis_variables
      ],
      function(x) {
        sum(
          is.na(x)
        )
      }
    )
)


missing_data$missing_percent <-
  100 *
  missing_data$missing_n /
  nrow(dryland_data)


cat("\n============================================================\n")
cat("MISSING DATA\n")
cat("============================================================\n\n")


print(missing_data)


# ------------------------------------------------------------
# 15. Species-richness complete-case check
# ------------------------------------------------------------

sr_model_variables <- c(
  "sr",
  "avg_aridity_index",
  "water",
  "shrub_and_scrub",
  "built",
  "duration_minutes",
  "effort_distance_km",
  "protocol_name",
  "locality_id",
  "NAME10"
)


sr_complete <- complete.cases(
  dryland_data[
    sr_model_variables
  ]
)


cat("\n============================================================\n")
cat("RICHNESS COMPLETE-CASE QC\n")
cat("============================================================\n\n")


cat(
  "Rows before complete-case restriction: ",
  nrow(dryland_data),
  "\n",
  sep = ""
)


cat(
  "Complete richness-model rows: ",
  sum(sr_complete),
  "\n",
  sep = ""
)


cat(
  "Rows missing >=1 richness variable: ",
  sum(!sr_complete),
  "\n",
  sep = ""
)


cat(
  "Percent retained: ",
  round(
    100 *
      mean(sr_complete),
    2
  ),
  "%\n",
  sep = ""
)


# ------------------------------------------------------------
# 16. Check aridity consistency within cities
# ------------------------------------------------------------
#
# Aridity is a city-level predictor.
# Each city should therefore have one aridity value.
# ------------------------------------------------------------

aridity_consistency <- dryland_data %>%
  group_by(
    NAME10
  ) %>%
  summarise(
    n_aridity_values =
      n_distinct(
        avg_aridity_index,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


cat("\n============================================================\n")
cat("CITY-LEVEL ARIDITY QC\n")
cat("============================================================\n\n")


print(
  aridity_consistency %>%
    count(
      n_aridity_values
    )
)


# ------------------------------------------------------------
# 17. Check habitat consistency within localities
# ------------------------------------------------------------
#
# Water, shrub/scrub, and built cover are locality-level
# predictors.
# ------------------------------------------------------------

habitat_consistency <- dryland_data %>%
  group_by(
    locality_id
  ) %>%
  summarise(
    n_water_values =
      n_distinct(
        water,
        na.rm = TRUE
      ),
    
    n_shrub_values =
      n_distinct(
        shrub_and_scrub,
        na.rm = TRUE
      ),
    
    n_built_values =
      n_distinct(
        built,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


habitat_consistency_summary <- habitat_consistency %>%
  summarise(
    multiple_water =
      sum(
        n_water_values > 1
      ),
    
    multiple_shrub =
      sum(
        n_shrub_values > 1
      ),
    
    multiple_built =
      sum(
        n_built_values > 1
      )
  )


cat("\n============================================================\n")
cat("LOCALITY HABITAT QC\n")
cat("============================================================\n\n")


print(
  habitat_consistency_summary
)


# ------------------------------------------------------------
# 18. Confirm locality-city nesting
# ------------------------------------------------------------

locality_city_nesting <- dryland_data %>%
  group_by(
    locality_id
  ) %>%
  summarise(
    n_cities =
      n_distinct(
        NAME10
      ),
    
    .groups = "drop"
  )


cat("\n============================================================\n")
cat("LOCALITY-CITY NESTING\n")
cat("============================================================\n\n")


print(
  locality_city_nesting %>%
    count(
      n_cities
    )
)


# ------------------------------------------------------------
# 19. Land-cover total QC
# ------------------------------------------------------------

landcover_variables <- c(
  "bare",
  "built",
  "crops",
  "flooded_vegetation",
  "grass",
  "shrub_and_scrub",
  "snow_and_ice",
  "trees",
  "water"
)


missing_landcover_columns <- setdiff(
  landcover_variables,
  names(dryland_data)
)


if (
  length(
    missing_landcover_columns
  ) > 0
) {
  
  stop(
    "Missing land-cover variables:\n",
    paste(
      missing_landcover_columns,
      collapse = "\n"
    )
  )
}


landcover_complete <- dryland_data %>%
  filter(
    if_all(
      all_of(
        landcover_variables
      ),
      ~ !is.na(.)
    )
  ) %>%
  mutate(
    landcover_total =
      bare +
      built +
      crops +
      flooded_vegetation +
      grass +
      shrub_and_scrub +
      snow_and_ice +
      trees +
      water
  )


landcover_total_qc <- landcover_complete %>%
  summarise(
    min_total =
      min(landcover_total),
    
    median_total =
      median(landcover_total),
    
    mean_total =
      mean(landcover_total),
    
    max_total =
      max(landcover_total)
  )


cat("\n============================================================\n")
cat("LAND-COVER TOTAL QC\n")
cat("============================================================\n\n")


print(
  landcover_total_qc
)


# ============================================================
# SAMPLING-STRUCTURE ASSESSMENT
# ============================================================


# ------------------------------------------------------------
# 20. Number of checklists per locality
# ------------------------------------------------------------

locality_sampling <- dryland_data %>%
  count(
    NAME10,
    locality_id,
    name = "n_checklists"
  )


cat("\n============================================================\n")
cat("CHECKLISTS PER LOCALITY\n")
cat("============================================================\n\n")


summary(
  locality_sampling$n_checklists
)


quantile(
  locality_sampling$n_checklists,
  probs = c(
    0,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.99,
    1
  )
)


# ------------------------------------------------------------
# 21. Examine locality checklist thresholds
# ------------------------------------------------------------

locality_thresholds <- c(
  1,
  2,
  3,
  5,
  10,
  20
)


locality_threshold_summary <- lapply(
  locality_thresholds,
  function(x) {
    
    tmp <- locality_sampling %>%
      filter(
        n_checklists >= x
      )
    
    
    data.frame(
      min_checklists_per_locality = x,
      n_localities = nrow(tmp),
      n_cities = n_distinct(tmp$NAME10)
    )
  }
) %>%
  bind_rows()


cat("\n============================================================\n")
cat("LOCALITY ELIGIBILITY THRESHOLDS\n")
cat("============================================================\n\n")


print(
  locality_threshold_summary
)


# ------------------------------------------------------------
# 22. Define eligible localities
# ------------------------------------------------------------
#
# Jackson suggested restricting the analysis to localities with
# at least 5 checklists.
# ------------------------------------------------------------

eligible_localities <- locality_sampling %>%
  filter(
    n_checklists >=
      MIN_CHECKLISTS_PER_LOCALITY
  )


cat(
  "\nEligible localities with >= ",
  MIN_CHECKLISTS_PER_LOCALITY,
  " checklists: ",
  nrow(eligible_localities),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 23. Number of eligible localities per city
# ------------------------------------------------------------

city_eligibility <- eligible_localities %>%
  count(
    NAME10,
    name = "n_eligible_localities"
  )


cat("\n============================================================\n")
cat("ELIGIBLE LOCALITIES PER CITY\n")
cat("============================================================\n\n")


summary(
  city_eligibility$n_eligible_localities
)


quantile(
  city_eligibility$n_eligible_localities,
  probs = c(
    0,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    0.99,
    1
  )
)


# ------------------------------------------------------------
# 24. Examine city-locality thresholds
# ------------------------------------------------------------

city_thresholds <- c(
  1,
  3,
  5,
  10,
  20,
  50
)


city_threshold_summary <- lapply(
  city_thresholds,
  function(x) {
    
    data.frame(
      min_eligible_localities_per_city = x,
      
      n_cities =
        sum(
          city_eligibility$n_eligible_localities >= x
        )
    )
  }
) %>%
  bind_rows()


cat("\n============================================================\n")
cat("CITY ELIGIBILITY THRESHOLDS\n")
cat("============================================================\n\n")


print(
  city_threshold_summary
)


# ------------------------------------------------------------
# 25. Define eligible cities
# ------------------------------------------------------------
#
# Jackson suggested retaining cities with at least 5 eligible
# locality IDs.
# ------------------------------------------------------------

eligible_cities <- city_eligibility %>%
  filter(
    n_eligible_localities >=
      MIN_LOCALITIES_PER_CITY
  )


cat(
  "\nEligible cities with >= ",
  MIN_LOCALITIES_PER_CITY,
  " eligible localities: ",
  nrow(eligible_cities),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 26. Create the eligible analysis pool
# ------------------------------------------------------------
#
# This is NOT the random subsample.
#
# This object contains all checklists belonging to:
#
#   - localities with >= 5 checklists
#   - cities with >= 5 eligible localities
#
# Script 02 will sample from this pool.
# ------------------------------------------------------------

eligible_pool <- dryland_data %>%
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


cat("\n============================================================\n")
cat("ELIGIBLE ANALYSIS POOL\n")
cat("============================================================\n\n")


eligible_pool_summary <- eligible_pool %>%
  summarise(
    n_checklists = n(),
    n_localities = n_distinct(locality_id),
    n_cities = n_distinct(NAME10),
    n_states = n_distinct(state_name)
  )


print(
  eligible_pool_summary
)


# ------------------------------------------------------------
# 27. Eligible sampling imbalance among cities
# ------------------------------------------------------------

eligible_city_structure <- eligible_pool %>%
  group_by(
    NAME10
  ) %>%
  summarise(
    n_checklists =
      n(),
    
    n_localities =
      n_distinct(
        locality_id
      ),
    
    .groups = "drop"
  )


cat("\n============================================================\n")
cat("ELIGIBLE CITY SAMPLING STRUCTURE\n")
cat("============================================================\n\n")


eligible_city_structure %>%
  summarise(
    min_checklists =
      min(n_checklists),
    
    median_checklists =
      median(n_checklists),
    
    mean_checklists =
      mean(n_checklists),
    
    max_checklists =
      max(n_checklists),
    
    min_localities =
      min(n_localities),
    
    median_localities =
      median(n_localities),
    
    mean_localities =
      mean(n_localities),
    
    max_localities =
      max(n_localities)
  ) %>%
  print()


# ------------------------------------------------------------
# 28. Compare potential sampling designs
# ------------------------------------------------------------
#
# Jackson's suggested approach:
#
#   Sample 5 checklists from every eligible locality.
#
# Alternative fully balanced design:
#
#   Sample 5 eligible localities from each eligible city,
#   then sample 5 checklists from each locality.
#
# No random sampling is performed here.
# ------------------------------------------------------------

eligible_localities_final <- eligible_localities %>%
  semi_join(
    eligible_cities,
    by = "NAME10"
  )


sampling_design_summary <- data.frame(
  design = c(
    "All eligible data",
    "5 checklists per eligible locality",
    "5 localities per city x 5 checklists"
  ),
  
  expected_checklists = c(
    nrow(eligible_pool),
    
    5 *
      nrow(
        eligible_localities_final
      ),
    
    5 *
      5 *
      nrow(
        eligible_cities
      )
  )
)


cat("\n============================================================\n")
cat("POTENTIAL SAMPLING DESIGNS\n")
cat("============================================================\n\n")


print(
  sampling_design_summary
)


# ------------------------------------------------------------
# 29. Check whether observation date is available
# ------------------------------------------------------------
#
# Jackson also suggested examining seasonality.
#
# We will not define breeding/non-breeding seasons yet.
# First determine whether observation_date is available in
# combined_data.RDS.
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SEASONALITY DATA AVAILABILITY\n")
cat("============================================================\n\n")


if (
  "observation_date" %in%
  names(dryland_data)
) {
  
  cat(
    "observation_date is available.\n"
  )
  
  
  dryland_data <- dryland_data %>%
    mutate(
      observation_date =
        as.Date(
          observation_date
        ),
      
      month =
        as.integer(
          format(
            observation_date,
            "%m"
          )
        )
    )
  
  
  monthly_sampling <- dryland_data %>%
    count(
      month
    )
  
  
  cat("\nChecklist counts by month:\n")
  
  print(
    monthly_sampling
  )
  
  
} else {
  
  cat(
    "observation_date is NOT present in combined_data.RDS.\n"
  )
  
  cat(
    "Breeding/non-breeding analyses cannot yet be evaluated ",
    "from this object.\n",
    sep = ""
  )
}


# ============================================================
# RESPONSE AND PREDICTOR INSPECTION
# ============================================================


# ------------------------------------------------------------
# 30. Response-variable distributions
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("RESPONSE DISTRIBUTIONS\n")
cat("============================================================\n\n")


summary(
  dryland_data$sr
)

summary(
  dryland_data$FRic
)

summary(
  dryland_data$FDis
)

summary(
  dryland_data$FEve
)


dryland_data %>%
  summarise(
    zero_sr =
      sum(
        sr == 0,
        na.rm = TRUE
      ),
    
    zero_FRic =
      sum(
        FRic == 0,
        na.rm = TRUE
      ),
    
    zero_FDis =
      sum(
        FDis == 0,
        na.rm = TRUE
      ),
    
    zero_FEve =
      sum(
        FEve == 0,
        na.rm = TRUE
      )
  ) %>%
  print()


# ------------------------------------------------------------
# 31. Species-richness distribution
# ------------------------------------------------------------

ggplot(
  dryland_data,
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


cat("\nSpecies richness mean:\n")

print(
  mean(
    dryland_data$sr,
    na.rm = TRUE
  )
)


cat("\nSpecies richness variance:\n")

print(
  var(
    dryland_data$sr,
    na.rm = TRUE
  )
)


cat("\nVariance / mean:\n")

print(
  var(
    dryland_data$sr,
    na.rm = TRUE
  ) /
    mean(
      dryland_data$sr,
      na.rm = TRUE
    )
)


# ------------------------------------------------------------
# 32. Functional diversity distributions
# ------------------------------------------------------------

ggplot(
  dryland_data,
  aes(
    x = FRic
  )
) +
  geom_histogram(
    binwidth = 0.02,
    na.rm = TRUE
  ) +
  labs(
    x = "Functional richness",
    y = "Number of checklists"
  ) +
  theme_classic()


ggplot(
  dryland_data,
  aes(
    x = FDis
  )
) +
  geom_histogram(
    binwidth = 0.01,
    na.rm = TRUE
  ) +
  labs(
    x = "Functional dispersion",
    y = "Number of checklists"
  ) +
  theme_classic()


ggplot(
  dryland_data,
  aes(
    x = FEve
  )
) +
  geom_histogram(
    binwidth = 0.02,
    na.rm = TRUE
  ) +
  labs(
    x = "Functional evenness",
    y = "Number of checklists"
  ) +
  theme_classic()


# ------------------------------------------------------------
# 33. FDis zero values
# ------------------------------------------------------------

dryland_data %>%
  group_by(
    FDis_zero =
      FDis == 0
  ) %>%
  summarise(
    n = n(),
    
    min_sr =
      min(
        sr,
        na.rm = TRUE
      ),
    
    median_sr =
      median(
        sr,
        na.rm = TRUE
      ),
    
    mean_sr =
      mean(
        sr,
        na.rm = TRUE
      ),
    
    max_sr =
      max(
        sr,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  ) %>%
  print()


# ------------------------------------------------------------
# 34. Habitat correlations
# ------------------------------------------------------------

habitat_data <- dryland_data %>%
  distinct(
    locality_id,
    water,
    shrub_and_scrub,
    built
  )


habitat_cor <- cor(
  habitat_data %>%
    select(
      water,
      shrub_and_scrub,
      built
    )
)


cat("\n============================================================\n")
cat("HABITAT CORRELATIONS\n")
cat("============================================================\n\n")


print(
  habitat_cor
)


# ------------------------------------------------------------
# 35. Aridity / habitat correlations
# ------------------------------------------------------------

predictor_data <- dryland_data %>%
  distinct(
    locality_id,
    avg_aridity_index,
    water,
    shrub_and_scrub,
    built
  )


predictor_cor <- cor(
  predictor_data %>%
    select(
      avg_aridity_index,
      water,
      shrub_and_scrub,
      built
    )
)


cat("\n============================================================\n")
cat("ARIDITY / HABITAT CORRELATIONS\n")
cat("============================================================\n\n")


print(
  predictor_cor
)


# ------------------------------------------------------------
# 36. Aridity distribution across cities
# ------------------------------------------------------------

aridity_city_data <- dryland_data %>%
  distinct(
    NAME10,
    avg_aridity_index
  )


ggplot(
  aridity_city_data,
  aes(
    x = avg_aridity_index
  )
) +
  geom_histogram(
    binwidth = 2
  ) +
  labs(
    x = "Aridity index",
    y = "Number of cities"
  ) +
  theme_classic()


summary(
  aridity_city_data$avg_aridity_index
)


# IMPORTANT:
#
# Lower De Martonne values = more arid.
# Higher De Martonne values = wetter / less arid.


# ------------------------------------------------------------
# 37. Effort distributions
# ------------------------------------------------------------

ggplot(
  dryland_data,
  aes(
    x = duration_minutes
  )
) +
  geom_histogram(
    binwidth = 10
  ) +
  labs(
    x = "Checklist duration (minutes)",
    y = "Number of checklists"
  ) +
  theme_classic()


ggplot(
  dryland_data,
  aes(
    x = effort_distance_km
  )
) +
  geom_histogram(
    binwidth = 0.25
  ) +
  labs(
    x = "Distance traveled (km)",
    y = "Number of checklists"
  ) +
  theme_classic()


# ------------------------------------------------------------
# 38. Final inspection summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("DATA INSPECTION COMPLETE\n")
cat("============================================================\n\n")


cat(
  "Raw rows:                          ",
  nrow(combined_data),
  "\n",
  sep = ""
)


cat(
  "Raw unique checklist IDs:          ",
  n_distinct(
    combined_data$checklist_id
  ),
  "\n",
  sep = ""
)


cat(
  "Duplicate checklist IDs detected:  ",
  nrow(
    duplicate_checklists
  ),
  "\n",
  sep = ""
)


cat(
  "Rows after effort filters:         ",
  nrow(
    dryland_data
  ),
  "\n",
  sep = ""
)


cat(
  "Eligible analysis-pool rows:       ",
  nrow(
    eligible_pool
  ),
  "\n",
  sep = ""
)


cat(
  "Eligible localities:               ",
  n_distinct(
    eligible_pool$locality_id
  ),
  "\n",
  sep = ""
)


cat(
  "Eligible cities:                   ",
  n_distinct(
    eligible_pool$NAME10
  ),
  "\n",
  sep = ""
)


cat(
  "Richness complete-case retention:  ",
  round(
    100 *
      mean(
        sr_complete
      ),
    2
  ),
  "%\n",
  sep = ""
)


cat("\n")


# ------------------------------------------------------------
# 39. Important note for downstream modeling
# ------------------------------------------------------------

cat("============================================================\n")
cat("NEXT STEP\n")
cat("============================================================\n\n")


cat(
  "Script 02 should draw a reproducible random subsample from\n",
  "eligible_pool using the sampling design selected after\n",
  "reviewing the summaries above.\n",
  sep = ""
)


if (
  nrow(
    duplicate_checklists
  ) > 0
) {
  
  cat(
    "\nWARNING:\n",
    "This run used temporary checklist deduplication.\n",
    "All final models must be rerun using Jackson's corrected\n",
    "combined_data.RDS once it becomes available.\n",
    sep = ""
  )
}