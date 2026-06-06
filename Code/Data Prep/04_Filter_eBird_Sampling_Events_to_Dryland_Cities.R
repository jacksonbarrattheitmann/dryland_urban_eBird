# ============================================================
# 04_Filter_eBird_Sampling_Events_to_Dryland_Cities.R
#
# Purpose:
# Filter eBird Sampling Event Data to complete checklists from
# 2015-2025 that fall inside dryland city polygons.
#
# Notes:
# The Sampling Event Data file used here does not include a
# protocol column, so this script does NOT filter to stationary
# checklists. It filters by checklist completeness, observer number,
# duration, valid coordinates, and dryland city polygon location.
#
# Approach:
# 1. Stream the large gzipped Sampling Event Data file with fread().
# 2. Keep only the columns needed for filtering and spatial joining.
# 3. Apply table filters.
# 4. Convert checklist coordinates to spatial points.
# 5. Spatially join checklists to dryland city polygons.
# 6. Save checklist table, checklist IDs, and city-level summary.
# ============================================================

library(sf)
library(dplyr)
library(readr)
library(data.table)

# ------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------

sed_path <- "Data/Raw/eBird/Unpacked/ebd_US_relApr-2026_sampling.txt.gz"

dryland_cities_path <- "Data/Processed/Shape Files/Cities/dryland_cities_50pct.shp"

ebird_out_dir <- "Data/Processed/eBird"
dir.create(ebird_out_dir, recursive = TRUE, showWarnings = FALSE)

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

file.remove(
  checklists_out,
  checklist_ids_out,
  summary_by_city_out
)

# ------------------------------------------------------------
# 2. Read dryland city polygons
# ------------------------------------------------------------

# This is a shapefile, so several field names were shortened.
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

# Broad bounding box around all dryland city polygons.
# This quickly reduces the eBird data before the exact polygon join.
bb <- st_bbox(dryland_cities)

xmin <- as.numeric(bb["xmin"])
xmax <- as.numeric(bb["xmax"])
ymin <- as.numeric(bb["ymin"])
ymax <- as.numeric(bb["ymax"])

cat("\nDryland city bounding box:\n")
print(bb)

# ------------------------------------------------------------
# 3. Stream and read needed Sampling Event Data columns
# ------------------------------------------------------------

# This streams the compressed file instead of manually decompressing it.
# It still scans the large file, so this step may take a while.
cmd <- paste("gzip -cd", shQuote(sed_path))

sampling_events <- fread(
  cmd = cmd,
  sep = "\t",
  quote = "",
  showProgress = TRUE,
  select = c(
    "SAMPLING EVENT IDENTIFIER",
    "OBSERVATION DATE",
    "LATITUDE",
    "LONGITUDE",
    "ALL SPECIES REPORTED",
    "DURATION MINUTES",
    "EFFORT DISTANCE KM",
    "NUMBER OBSERVERS"
  )
)

cat("\nSampling rows read:\n")
print(nrow(sampling_events))

# ------------------------------------------------------------
# 4. Rename columns
# ------------------------------------------------------------

setnames(
  sampling_events,
  old = c(
    "SAMPLING EVENT IDENTIFIER",
    "OBSERVATION DATE",
    "LATITUDE",
    "LONGITUDE",
    "ALL SPECIES REPORTED",
    "DURATION MINUTES",
    "EFFORT DISTANCE KM",
    "NUMBER OBSERVERS"
  ),
  new = c(
    "checklist_id",
    "observation_date",
    "latitude",
    "longitude",
    "all_species_reported",
    "duration_minutes",
    "effort_distance_km",
    "number_observers"
  )
)

# ------------------------------------------------------------
# 5. Apply table filters
# ------------------------------------------------------------

# Filters:
# - date range: 2015-2025
# - complete checklists only
# - one observer
# - duration between 5 and 240 minutes
# - valid coordinates
# - coordinates inside broad dryland-city bounding box
#
# The exact city polygon filter happens in the next section.
sampling_events <- sampling_events[
  observation_date >= "2015-01-01" &
    observation_date <= "2025-12-31" &
    all_species_reported == 1 &
    number_observers == 1 &
    duration_minutes >= 5 &
    duration_minutes <= 240 &
    !is.na(latitude) &
    !is.na(longitude) &
    latitude >= ymin &
    latitude <= ymax &
    longitude >= xmin &
    longitude <= xmax
]

cat("\nSampling events after table filters:\n")
print(nrow(sampling_events))

if (nrow(sampling_events) == 0) {
  stop("No sampling events passed table filters. Check column values or filters.")
}

sampling_events <- as_tibble(sampling_events) %>%
  mutate(
    observation_date = as.Date(observation_date),
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude),
    duration_minutes = as.numeric(duration_minutes),
    effort_distance_km = as.numeric(effort_distance_km),
    number_observers = as.numeric(number_observers)
  )

# ------------------------------------------------------------
# 6. Spatially join checklists to dryland city polygons
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
# 7. Save checklist table and checklist IDs
# ------------------------------------------------------------

dryland_city_checklists_table <- dryland_city_checklists %>%
  st_drop_geometry() %>%
  arrange(state_name, city_ascii, observation_date) %>%
  dplyr::select(
    checklist_id,
    observation_date,
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
    dryland_percent
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
# 8. Save summary by city
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
# 9. Final checks
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