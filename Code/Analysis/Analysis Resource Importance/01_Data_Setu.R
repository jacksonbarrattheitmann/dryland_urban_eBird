library(tidyverse)

# Settings

temporary_deduplicate <- TRUE

input_file <- file.path(
  "Intermediate_data",
  "dfs_for_analysis_combined",
  "combined_data.RDS"
)

output_dir <- file.path(
  "Outputs",
  "Analysis_Resource_Importance",
  "01_Data_Setup"
)

data_dir <- file.path(output_dir, "Data")
table_dir <- file.path(output_dir, "Table")

dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

# Load data

if (!file.exists(input_file)) {
  stop("combined_data.RDS was not found.")
}

combined_data <- read_rds(input_file)

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

missing_variables <- setdiff(required_variables, names(combined_data))

if (length(missing_variables) > 0) {
  stop(
    "Missing required variables: ",
    paste(missing_variables, collapse = ", ")
  )
}

# Checklist filtering

duplicate_checklists <- combined_data %>%
  count(checklist_id) %>%
  filter(n > 1)

if (nrow(duplicate_checklists) > 0 && !temporary_deduplicate) {
  stop("Duplicate checklist IDs are present.")
}

analysis_base <- if (temporary_deduplicate) {
  combined_data %>%
    distinct(checklist_id, .keep_all = TRUE)
} else {
  combined_data
}

analysis_base <- analysis_base %>%
  mutate(
    effort_distance_km = if_else(
      protocol_name == "Stationary",
      0,
      effort_distance_km
    )
  ) %>%
  filter(
    protocol_name %in% c("Stationary", "Traveling"),
    duration_minutes <= 300,
    effort_distance_km <= 5
  )

# Spatial checks

aridity_qc <- analysis_base %>%
  filter(!is.na(avg_aridity_index)) %>%
  group_by(NAME10) %>%
  summarise(
    n_aridity_values = n_distinct(avg_aridity_index),
    .groups = "drop"
  )

if (any(aridity_qc$n_aridity_values != 1)) {
  stop("At least one city has more than one non-missing aridity value.")
}

habitat_qc <- analysis_base %>%
  group_by(NAME10, locality_id) %>%
  summarise(
    n_water = n_distinct(water, na.rm = TRUE),
    n_shrub = n_distinct(shrub_and_scrub, na.rm = TRUE),
    n_built = n_distinct(built, na.rm = TRUE),
    .groups = "drop"
  )

if (
  any(habitat_qc$n_water > 1) ||
  any(habitat_qc$n_shrub > 1) ||
  any(habitat_qc$n_built > 1)
) {
  stop("At least one locality has multiple values for a habitat predictor.")
}

locality_city_qc <- analysis_base %>%
  group_by(locality_id) %>%
  summarise(
    n_cities = n_distinct(NAME10),
    .groups = "drop"
  )

if (any(locality_city_qc$n_cities > 1)) {
  stop("At least one locality_id occurs in more than one city.")
}

# Habitat cover

max_habitat_value <- max(
  c(
    analysis_base$water,
    analysis_base$shrub_and_scrub,
    analysis_base$built
  ),
  na.rm = TRUE
)

habitat_scale <- if (max_habitat_value <= 1.05) 100 else 1

analysis_base <- analysis_base %>%
  mutate(
    water_pp = water * habitat_scale,
    shrub_pp = shrub_and_scrub * habitat_scale,
    built_pp = built * habitat_scale
  )

# City habitat means use one row per locality so heavily sampled
# localities do not receive extra weight.

locality_habitat <- analysis_base %>%
  distinct(
    NAME10,
    state_name,
    locality_id,
    water_pp,
    shrub_pp,
    built_pp
  )

city_habitat_means <- locality_habitat %>%
  group_by(NAME10) %>%
  summarise(
    state_name = first(state_name),
    n_localities_habitat = n(),
    city_mean_water_pp = mean(water_pp, na.rm = TRUE),
    city_mean_shrub_pp = mean(shrub_pp, na.rm = TRUE),
    city_mean_built_pp = mean(built_pp, na.rm = TRUE),
    city_sd_water_pp = sd(water_pp, na.rm = TRUE),
    city_sd_shrub_pp = sd(shrub_pp, na.rm = TRUE),
    city_sd_built_pp = sd(built_pp, na.rm = TRUE),
    .groups = "drop"
  )

analysis_base <- analysis_base %>%
  left_join(
    city_habitat_means %>% select(-state_name),
    by = "NAME10"
  ) %>%
  mutate(
    water_within_pp = water_pp - city_mean_water_pp,
    shrub_within_pp = shrub_pp - city_mean_shrub_pp,
    built_within_pp = built_pp - city_mean_built_pp,
    water_within_10 = water_within_pp / 10,
    shrub_within_10 = shrub_within_pp / 10,
    built_within_10 = built_within_pp / 10,
    city_mean_water_10 = city_mean_water_pp / 10,
    city_mean_shrub_10 = city_mean_shrub_pp / 10,
    city_mean_built_10 = city_mean_built_pp / 10,
    duration_10 = duration_minutes / 10
  )

# Check that locality-level habitat deviations are centered within cities.

within_city_centering_qc <- analysis_base %>%
  distinct(
    NAME10,
    locality_id,
    water_within_pp,
    shrub_within_pp,
    built_within_pp
  ) %>%
  group_by(NAME10) %>%
  summarise(
    mean_water_within_pp = mean(water_within_pp, na.rm = TRUE),
    mean_shrub_within_pp = mean(shrub_within_pp, na.rm = TRUE),
    mean_built_within_pp = mean(built_within_pp, na.rm = TRUE),
    .groups = "drop"
  )

centering_tolerance <- 1e-8

if (
  any(abs(within_city_centering_qc$mean_water_within_pp) > centering_tolerance) ||
  any(abs(within_city_centering_qc$mean_shrub_within_pp) > centering_tolerance) ||
  any(abs(within_city_centering_qc$mean_built_within_pp) > centering_tolerance)
) {
  warning("At least one within-city habitat mean differs from zero beyond tolerance.")
}

# Save outputs

write_csv(
  city_habitat_means,
  file.path(table_dir, "Table_01_City_Habitat_Means.csv")
)

write_csv(
  within_city_centering_qc,
  file.path(table_dir, "Table_02_Within_City_Habitat_Centering_QC.csv")
)

write_rds(
  analysis_base,
  file.path(data_dir, "Resource_Importance_Analysis_Base.rds")
)
