# ============================================================
# 04_Filter_eBird_Sampling_Events_to_Dryland_Cities.R
#
# Purpose:
# Filter eBird Sampling Event Data to complete stationary
# checklists from 2015-2025 that fall inside dryland city polygons.
#
# Approach:
# 1. Use auk to filter the large eBird Sampling Event Data file.
# 2. Read the reduced file into R.
# 3. Spatially join checklist points to dryland city polygons.
# 4. Save checklist table, checklist IDs, and city summary.
# ============================================================

library(sf)
library(dplyr)
library(readr)
library(auk)

# ------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------

sed_path <- "Data/Raw/eBird/Unpacked/ebd_US_relApr-2026_sampling.txt.gz"

dryland_cities_path <- "Data/Processed/Shape Files/Cities/dryland_cities_50pct.shp"

ebird_out_dir <- "Data/Processed/eBird"
dir.create(ebird_out_dir, recursive = TRUE, showWarnings = FALSE)

# Intermediate auk-filtered sampling file
sampling_filtered_out <- file.path(
  ebird_out_dir,
  "sampling_events_auk_filtered_2015_2025.txt"
)

# Final outputs
checklists_out <- file.path(
  ebird_out_dir,
  "dryland_city_checklists_2015_2025.csv"
)

checklist_ids_out <- file.path(
  ebird_out_dir,
  "dryland_city_checklist_ids_2015_2025.csv"
)

summary_by_city_out <- file.path(
  ebird_out_dir,
  "dryland_city_checklist_summary_2015_2025.csv"
)

# Remove old outputs if rerunning
file.remove(
  sampling_filtered_out,
  checklists_out,
  checklist_ids_out,
  summary_by_city_out
)

# ------------------------------------------------------------
# 2. Confirm eBird sampling file is valid
# ------------------------------------------------------------

system2("file", sed_path)
system2("gzip", c("-t", sed_path))

# ------------------------------------------------------------
# 3. Read dryland city polygons
# ------------------------------------------------------------

# This is a shapefile, so some field names were shortened.
# We rename them back to readable names here.
dryland_cities <- st_read(dryland_cities_path) %>%
  st_make_valid() %>%
  st_transform(4326) %>%
  dplyr::select(
    city_id,
    NAME10,
    city_ascii = city_sc,
    state_name = stat_nm,
    county_name = cnty_nm,
    dryland_percent = dryl__2,
    geometry
  )

# ------------------------------------------------------------
# 4. Filter Sampling Event Data with auk
# ------------------------------------------------------------

# auk filters the large eBird file without loading the whole thing
# into memory. We filter to:
# - 2015-2025
# - stationary checklists
# - complete checklists
# - duration 5-240 minutes
#
# We do the one-observer filter after reading the reduced file.

sed_filter <- auk_sampling(sed_path) %>%
  auk_date(date = c("2015-01-01", "2025-12-31")) %>%
  auk_protocol(protocol = "Stationary") %>%
  auk_complete() %>%
  auk_duration(duration = c(5, 240))

auk_filter(
  sed_filter,
  file = sampling_filtered_out,
  overwrite = TRUE
)

cat("\nauk-filtered sampling file info:\n")
print(file.info(sampling_filtered_out)[, c("size", "mtime")])

if (!file.exists(sampling_filtered_out) || file.info(sampling_filtered_out)$size == 0) {
  stop("auk did not create a non-empty filtered sampling file.")
}

# ------------------------------------------------------------
# 5. Read filtered Sampling Event Data
# ------------------------------------------------------------

sampling_events <- read_sampling(sampling_filtered_out)

cat("\nSampling events after auk filters:\n")
print(nrow(sampling_events))

cat("\nColumn names in filtered sampling data:\n")
print(names(sampling_events))

# ------------------------------------------------------------
# 6. Apply one-observer and coordinate filters
# ------------------------------------------------------------

sampling_events <- sampling_events %>%
  filter(
    number_observers == 1,
    !is.na(latitude),
    !is.na(longitude)
  )

cat("\nSampling events after one-observer and coordinate filters:\n")
print(nrow(sampling_events))

if (nrow(sampling_events) == 0) {
  stop("No sampling events remained after one-observer and coordinate filters.")
}

# ------------------------------------------------------------
# 7. Spatially join checklists to dryland city polygons
# ------------------------------------------------------------

sampling_points <- sampling_events %>%
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  )

dryland_city_checklists <- sampling_points %>%
  st_join(
    dryland_cities,
    join = st_within,
    left = FALSE
  )

cat("\nChecklists inside dryland city polygons:\n")
print(nrow(dryland_city_checklists))

if (nrow(dryland_city_checklists) == 0) {
  stop("No checklists fell inside the dryland city polygons after spatial join.")
}

# ------------------------------------------------------------
# 8. Save checklist table and checklist IDs
# ------------------------------------------------------------

dryland_city_checklists_table <- dryland_city_checklists %>%
  st_drop_geometry() %>%
  arrange(state_name, city_ascii, observation_date) %>%
  dplyr::select(
    checklist_id = checklist_identifier,
    observation_date,
    protocol_type,
    all_species_reported,
    duration_minutes,
    effort_distance_km,
    number_observers,
    latitude,
    longitude,
    city_id,
    NAME10,
    city_ascii,
    state_name,
    county_name,
    dryland_percent,
    everything()
  )

write_csv(
  dryland_city_checklists_table,
  checklists_out
)

dryland_checklist_ids <- dryland_city_checklists_table %>%
  distinct(checklist_id)

write_csv(
  dryland_checklist_ids,
  checklist_ids_out
)

# ------------------------------------------------------------
# 9. Save summary by city
# ------------------------------------------------------------

dryland_city_summary <- dryland_city_checklists_table %>%
  group_by(
    city_id,
    NAME10,
    city_ascii,
    state_name,
    county_name,
    dryland_percent
  ) %>%
  summarise(
    n_checklists = n_distinct(checklist_id),
    n_unique_locations = n_distinct(latitude, longitude),
    first_date = min(observation_date, na.rm = TRUE),
    last_date = max(observation_date, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_checklists))

write_csv(
  dryland_city_summary,
  summary_by_city_out
)

# ------------------------------------------------------------
# 10. Final checks
# ------------------------------------------------------------

cat("\nDryland city checklists retained:\n")
print(nrow(dryland_city_checklists_table))

cat("\nUnique checklist IDs retained:\n")
print(nrow(dryland_checklist_ids))

cat("\nDryland cities with at least one checklist:\n")
dryland_city_checklists_table %>%
  distinct(city_id) %>%
  nrow() %>%
  print()

cat("\nTop dryland cities by checklist count:\n")
print(head(dryland_city_summary, 25))