# ============================================================
# Dryland Bird Project
# 01_Data_Inspection.R
#
# Purpose:
# Inspect the final combined analysis dataset to make sure
# variables, sampling structure, response variables, habitat
# predictors, aridity, and checklist-effort variables are
# appropriate for downstream modeling.
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
# 2. Load final combined dataset
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


# ------------------------------------------------------------
# 3. Dataset dimensions
# ------------------------------------------------------------

dim(combined_data)


# ------------------------------------------------------------
# 4. Column names and structure
# ------------------------------------------------------------

names(combined_data)

str(combined_data)


# ------------------------------------------------------------
# 5. Basic sampling structure
# ------------------------------------------------------------

combined_data %>%
  summarise(
    n_rows = n(),
    n_checklists = n_distinct(checklist_id),
    n_localities = n_distinct(locality_id),
    n_cities = n_distinct(NAME10),
    n_states = n_distinct(state_name)
  )


# ------------------------------------------------------------
# 6. Check duplicate checklist IDs
# ------------------------------------------------------------

duplicate_checklists <- combined_data %>%
  count(
    checklist_id
  ) %>%
  filter(
    n > 1
  )


cat(
  "\nDuplicated checklist IDs:",
  nrow(duplicate_checklists),
  "\n"
)


if (nrow(duplicate_checklists) > 0) {
  
  print(
    head(
      duplicate_checklists,
      20
    )
  )
}


# ------------------------------------------------------------
# 7. Examine checklist protocols
# ------------------------------------------------------------

combined_data %>%
  count(
    protocol_name,
    sort = TRUE
  )


# ------------------------------------------------------------
# 8. Examine effort variables by protocol
# ------------------------------------------------------------

combined_data %>%
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


# ------------------------------------------------------------
# 9. Create analysis dataframe
# ------------------------------------------------------------
#
# Stationary checklists have zero distance traveled.
#
# Retain:
#   Stationary
#   Traveling
#
# Effort filters:
#   duration <= 300 minutes
#   distance <= 5 km
# ------------------------------------------------------------

dryland_data <- combined_data %>%
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
# 10. Check analysis dataset
# ------------------------------------------------------------

dryland_data %>%
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


dryland_data %>%
  count(
    protocol_name,
    sort = TRUE
  )


# ------------------------------------------------------------
# 11. Quantify filtering
# ------------------------------------------------------------

filter_summary <- data.frame(
  stage = c(
    "combined_data",
    "analysis_data"
  ),
  
  n_checklists = c(
    nrow(combined_data),
    nrow(dryland_data)
  ),
  
  n_localities = c(
    n_distinct(combined_data$locality_id),
    n_distinct(dryland_data$locality_id)
  ),
  
  n_cities = c(
    n_distinct(combined_data$NAME10),
    n_distinct(dryland_data$NAME10)
  )
)


filter_summary


cat(
  "\nRows removed by protocol/effort filters:",
  nrow(combined_data) - nrow(dryland_data),
  "\n"
)


# ------------------------------------------------------------
# 12. Examine missing data in analysis variables
# ------------------------------------------------------------

model_variables <- c(
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
  "locality_id",
  "NAME10"
)


missing_variables <- setdiff(
  model_variables,
  names(dryland_data)
)


if (length(missing_variables) > 0) {
  
  stop(
    "The combined dataset is missing required variables:\n",
    paste(
      missing_variables,
      collapse = "\n"
    )
  )
}


missing_data <- data.frame(
  variable = model_variables,
  
  missing_n =
    sapply(
      dryland_data[
        model_variables
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


missing_data


# ------------------------------------------------------------
# 13. Richness-model complete-case check
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


cat(
  "\nSpecies-richness rows before complete-case check:",
  nrow(dryland_data),
  "\n"
)


cat(
  "Complete richness-model rows:",
  sum(sr_complete),
  "\n"
)


cat(
  "Rows missing >=1 richness-model variable:",
  sum(!sr_complete),
  "\n"
)


cat(
  "Percent retained:",
  round(
    100 *
      mean(sr_complete),
    2
  ),
  "%\n"
)


# ------------------------------------------------------------
# 14. Examine missing functional diversity values
# ------------------------------------------------------------

dryland_data %>%
  group_by(
    FRic_missing =
      is.na(FRic)
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
  )


# ------------------------------------------------------------
# 15. Compare FRic and FEve missing values
# ------------------------------------------------------------

table(
  FRic_missing =
    is.na(
      dryland_data$FRic
    ),
  
  FEve_missing =
    is.na(
      dryland_data$FEve
    )
)


# ------------------------------------------------------------
# 16. Examine aridity variation within cities
# ------------------------------------------------------------

dryland_data %>%
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
  ) %>%
  count(
    n_aridity_values
  )


# ------------------------------------------------------------
# 17. Examine habitat variation within localities
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


habitat_consistency %>%
  summarise(
    localities_with_multiple_water_values =
      sum(
        n_water_values > 1
      ),
    
    localities_with_multiple_shrub_values =
      sum(
        n_shrub_values > 1
      ),
    
    localities_with_multiple_built_values =
      sum(
        n_built_values > 1
      )
  )


# ------------------------------------------------------------
# 18. Examine sampling structure within cities
# ------------------------------------------------------------

dryland_data %>%
  group_by(
    NAME10
  ) %>%
  summarise(
    n_localities =
      n_distinct(locality_id),
    
    n_checklists =
      n(),
    
    .groups = "drop"
  ) %>%
  summarise(
    min_localities =
      min(n_localities),
    
    median_localities =
      median(n_localities),
    
    mean_localities =
      mean(n_localities),
    
    max_localities =
      max(n_localities),
    
    min_checklists =
      min(n_checklists),
    
    median_checklists =
      median(n_checklists),
    
    mean_checklists =
      mean(n_checklists),
    
    max_checklists =
      max(n_checklists)
  )


# ------------------------------------------------------------
# 19. Examine sampling structure within localities
# ------------------------------------------------------------

dryland_data %>%
  group_by(
    locality_id
  ) %>%
  summarise(
    n_checklists =
      n(),
    
    .groups = "drop"
  ) %>%
  summarise(
    min_checklists =
      min(n_checklists),
    
    median_checklists =
      median(n_checklists),
    
    mean_checklists =
      mean(n_checklists),
    
    max_checklists =
      max(n_checklists)
  )


# ------------------------------------------------------------
# 20. Examine response-variable distributions
# ------------------------------------------------------------

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
  )


# ------------------------------------------------------------
# 21. Examine zero functional dispersion values
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
  )


# ------------------------------------------------------------
# 22. Examine species-richness distribution
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


# Preliminary overdispersion check

mean(
  dryland_data$sr,
  na.rm = TRUE
)


var(
  dryland_data$sr,
  na.rm = TRUE
)


var(
  dryland_data$sr,
  na.rm = TRUE
) /
  mean(
    dryland_data$sr,
    na.rm = TRUE
  )


# ------------------------------------------------------------
# 23. Examine functional richness distribution
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


summary(
  dryland_data$FRic
)


# ------------------------------------------------------------
# 24. Examine functional dispersion distribution
# ------------------------------------------------------------

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


summary(
  dryland_data$FDis
)


# ------------------------------------------------------------
# 25. Examine functional evenness distribution
# ------------------------------------------------------------

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


summary(
  dryland_data$FEve
)


# ------------------------------------------------------------
# 26. Examine boundary values in functional evenness
# ------------------------------------------------------------

dryland_data %>%
  summarise(
    n_FEve_1 =
      sum(
        FEve == 1,
        na.rm = TRUE
      ),
    
    percent_FEve_1 =
      100 *
      sum(
        FEve == 1,
        na.rm = TRUE
      ) /
      sum(
        !is.na(FEve)
      )
  )


# ------------------------------------------------------------
# 27. Examine total land-cover coverage
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
    "Missing land-cover columns:\n",
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


landcover_complete %>%
  summarise(
    n_complete =
      n(),
    
    min_total =
      min(landcover_total),
    
    median_total =
      median(landcover_total),
    
    mean_total =
      mean(landcover_total),
    
    max_total =
      max(landcover_total)
  )


# ------------------------------------------------------------
# 28. Examine correlations among habitat predictors
# ------------------------------------------------------------
#
# Correlations are calculated at the locality level rather than
# the checklist level so heavily sampled localities do not
# receive additional weight.
# ------------------------------------------------------------

habitat_data <- dryland_data %>%
  filter(
    !is.na(water),
    !is.na(shrub_and_scrub),
    !is.na(built)
  ) %>%
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


habitat_cor


# ------------------------------------------------------------
# 29. Visualize habitat correlation matrix
# ------------------------------------------------------------

habitat_cor_df <- as.data.frame(
  as.table(
    habitat_cor
  )
)


ggplot(
  habitat_cor_df,
  aes(
    x = Var1,
    y = Var2,
    fill = Freq
  )
) +
  geom_tile(
    color = "white"
  ) +
  geom_text(
    aes(
      label =
        round(
          Freq,
          2
        )
    ),
    size = 4
  ) +
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  coord_equal() +
  labs(
    x = NULL,
    y = NULL,
    fill = "Correlation"
  ) +
  theme_classic()


# ------------------------------------------------------------
# 30. Examine habitat predictor distributions
# ------------------------------------------------------------

ggplot(
  habitat_data,
  aes(
    x = water
  )
) +
  geom_histogram(
    binwidth = 0.05
  ) +
  labs(
    x = "Proportion water",
    y = "Number of localities"
  ) +
  theme_classic()


ggplot(
  habitat_data,
  aes(
    x = shrub_and_scrub
  )
) +
  geom_histogram(
    binwidth = 0.05
  ) +
  labs(
    x = "Proportion shrub and scrub",
    y = "Number of localities"
  ) +
  theme_classic()


ggplot(
  habitat_data,
  aes(
    x = built
  )
) +
  geom_histogram(
    binwidth = 0.05
  ) +
  labs(
    x = "Proportion built",
    y = "Number of localities"
  ) +
  theme_classic()


# ------------------------------------------------------------
# 31. Examine aridity distribution across cities
# ------------------------------------------------------------

aridity_city_data <- dryland_data %>%
  filter(
    !is.na(
      avg_aridity_index
    )
  ) %>%
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
# 32. Examine relationships between aridity and habitat
# ------------------------------------------------------------

predictor_data <- dryland_data %>%
  filter(
    !is.na(avg_aridity_index),
    !is.na(water),
    !is.na(shrub_and_scrub),
    !is.na(built)
  ) %>%
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


predictor_cor


# ------------------------------------------------------------
# 33. Visualize predictor correlation matrix
# ------------------------------------------------------------

predictor_cor_df <- as.data.frame(
  as.table(
    predictor_cor
  )
)


ggplot(
  predictor_cor_df,
  aes(
    x = Var1,
    y = Var2,
    fill = Freq
  )
) +
  geom_tile(
    color = "white"
  ) +
  geom_text(
    aes(
      label =
        round(
          Freq,
          2
        )
    ),
    size = 4
  ) +
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  scale_x_discrete(
    labels = c(
      avg_aridity_index =
        "Aridity index",
      water =
        "Water",
      shrub_and_scrub =
        "Shrub and scrub",
      built =
        "Built"
    )
  ) +
  scale_y_discrete(
    labels = c(
      avg_aridity_index =
        "Aridity index",
      water =
        "Water",
      shrub_and_scrub =
        "Shrub and scrub",
      built =
        "Built"
    )
  ) +
  coord_equal() +
  labs(
    x = NULL,
    y = NULL,
    fill = "Correlation"
  ) +
  theme_classic()


# ------------------------------------------------------------
# 34. Examine checklist duration
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


summary(
  dryland_data$duration_minutes
)


# ------------------------------------------------------------
# 35. Examine checklist distance traveled
# ------------------------------------------------------------

ggplot(
  dryland_data,
  aes(
    x = effort_distance_km
  )
) +
  geom_histogram(
    binwidth = 0.25,
    na.rm = TRUE
  ) +
  labs(
    x = "Distance traveled (km)",
    y = "Number of checklists"
  ) +
  theme_classic()


summary(
  dryland_data$effort_distance_km
)


# ------------------------------------------------------------
# 36. Examine distance traveled by protocol
# ------------------------------------------------------------

dryland_data %>%
  group_by(
    protocol_name
  ) %>%
  summarise(
    n =
      n(),
    
    median_distance =
      median(
        effort_distance_km,
        na.rm = TRUE
      ),
    
    mean_distance =
      mean(
        effort_distance_km,
        na.rm = TRUE
      ),
    
    max_distance =
      max(
        effort_distance_km,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


dryland_data %>%
  filter(
    protocol_name == "Traveling"
  ) %>%
  summarise(
    min_distance =
      min(
        effort_distance_km,
        na.rm = TRUE
      ),
    
    q25_distance =
      quantile(
        effort_distance_km,
        0.25,
        na.rm = TRUE
      ),
    
    median_distance =
      median(
        effort_distance_km,
        na.rm = TRUE
      ),
    
    q75_distance =
      quantile(
        effort_distance_km,
        0.75,
        na.rm = TRUE
      ),
    
    q95_distance =
      quantile(
        effort_distance_km,
        0.95,
        na.rm = TRUE
      ),
    
    q99_distance =
      quantile(
        effort_distance_km,
        0.99,
        na.rm = TRUE
      ),
    
    max_distance =
      max(
        effort_distance_km,
        na.rm = TRUE
      )
  )


# ------------------------------------------------------------
# 37. Examine long-distance Traveling checklists
# ------------------------------------------------------------
#
# Use combined_data rather than dryland_data because
# dryland_data has already been restricted to <= 5 km.
# ------------------------------------------------------------

combined_data %>%
  filter(
    protocol_name == "Traveling",
    duration_minutes <= 300
  ) %>%
  summarise(
    n_total =
      n(),
    
    n_missing_distance =
      sum(
        is.na(effort_distance_km)
      ),
    
    n_over_5km =
      sum(
        effort_distance_km > 5,
        na.rm = TRUE
      ),
    
    n_over_10km =
      sum(
        effort_distance_km > 10,
        na.rm = TRUE
      ),
    
    n_over_20km =
      sum(
        effort_distance_km > 20,
        na.rm = TRUE
      ),
    
    n_over_50km =
      sum(
        effort_distance_km > 50,
        na.rm = TRUE
      )
  )


# ------------------------------------------------------------
# 38. Confirm locality-city nesting
# ------------------------------------------------------------

dryland_data %>%
  group_by(
    locality_id
  ) %>%
  summarise(
    n_cities =
      n_distinct(
        NAME10
      ),
    
    .groups = "drop"
  ) %>%
  count(
    n_cities
  )


# ------------------------------------------------------------
# 39. Examine number of observers
# ------------------------------------------------------------

table(
  dryland_data$number_observers
)


summary(
  dryland_data$number_observers
)


# ------------------------------------------------------------
# 40. Final inspection summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("DATA INSPECTION COMPLETE\n")
cat("============================================================\n\n")


cat(
  "Combined rows:                     ",
  nrow(combined_data),
  "\n",
  sep = ""
)


cat(
  "Combined unique checklist IDs:     ",
  n_distinct(
    combined_data$checklist_id
  ),
  "\n",
  sep = ""
)


cat(
  "Duplicate checklist IDs:           ",
  nrow(
    duplicate_checklists
  ),
  "\n",
  sep = ""
)


cat(
  "Rows after protocol/effort filters:",
  nrow(dryland_data),
  "\n",
  sep = ""
)


cat(
  "Richness complete-case rows:       ",
  sum(sr_complete),
  "\n",
  sep = ""
)


cat(
  "Richness complete-case retention:  ",
  round(
    100 *
      mean(sr_complete),
    2
  ),
  "%\n",
  sep = ""
)


cat(
  "Cities after effort filtering:     ",
  n_distinct(
    dryland_data$NAME10
  ),
  "\n",
  sep = ""
)


cat(
  "Localities after effort filtering: ",
  n_distinct(
    dryland_data$locality_id
  ),
  "\n",
  sep = ""
)


cat("\n")