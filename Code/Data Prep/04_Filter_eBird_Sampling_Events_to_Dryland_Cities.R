# ============================================================
# 04_Filter_eBird_Sampling_Events_to_Dryland_Cities_WITH_HOTSPOTS.R
#
# Purpose:
# Filter eBird Sampling Event Data to complete checklists from
# 2015-2025 that fall inside dryland city polygons.
#
# This revised version preserves eBird hotspot/location fields:
#   - locality_id = LOCALITY ID
#   - locality    = LOCALITY
#
# These fields are needed downstream so each checklist can be linked
# to an eBird hotspot/location for modelling and GEE extraction.
#
# Notes:
# The Sampling Event Data file used here does not include a protocol
# column, so this script does NOT filter to stationary checklists.
# It filters by checklist completeness, observer number, duration,
# valid coordinates, and dryland city polygon location.
#
# Approach:
# 1. Stream the large gzipped Sampling Event Data file with fread().
# 2. Keep only the columns needed for filtering and spatial joining.
# 3. Preserve LOCALITY ID and LOCALITY if present.
# 4. Apply table filters.
# 5. Convert checklist coordinates to spatial points.
# 6. Spatially join checklists to dryland city polygons.
# 7. Save checklist table, checklist IDs, hotspot list, and city summary.
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

hotspot_list_out <- file.path(
  ebird_out_dir,
  "dryland_city_hotspot_list_for_GEE_2015_2025.csv"
)

file.remove(
  checklists_out,
  checklist_ids_out,
  summary_by_city_out,
  hotspot_list_out
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
# 3. Inspect Sampling Event Data column names
# ------------------------------------------------------------

# fread(select = ...) errors if a requested column does not exist.
# So we first read the header only, then select only columns that exist.
sed_header <- fread(
  cmd = paste("gzip -cd", shQuote(sed_path)),
  sep = "\t",
  quote = "",
  nrows = 0,
  showProgress = FALSE
)

sed_names <- names(sed_header)

cat("\nSampling Event Data columns related to locality/location:\n")
print(
  grep(
    "LOCALITY|LOCATION|HOTSPOT",
    sed_names,
    ignore.case = TRUE,
    value = TRUE
  )
)

required_sed_cols <- c(
  "SAMPLING EVENT IDENTIFIER",
  "OBSERVATION DATE",
  "LATITUDE",
  "LONGITUDE",
  "ALL SPECIES REPORTED",
  "DURATION MINUTES",
  "EFFORT DISTANCE KM",
  "NUMBER OBSERVERS"
)

optional_sed_cols <- c(
  "LOCALITY ID",
  "LOCALITY"
)

missing_required_sed_cols <- setdiff(required_sed_cols, sed_names)

if (length(missing_required_sed_cols) > 0) {
  stop(
    "The Sampling Event Data file is missing required columns:\n",
    paste(missing_required_sed_cols, collapse = "\n")
  )
}

sed_cols_to_read <- c(
  required_sed_cols,
  intersect(optional_sed_cols, sed_names)
)

cat("\nColumns selected from Sampling Event Data:\n")
print(sed_cols_to_read)

# ------------------------------------------------------------
# 4. Stream and read needed Sampling Event Data columns
# ------------------------------------------------------------

# This streams the compressed file instead of manually decompressing it.
# It still scans the large file, so this step may take a while.
cmd <- paste("gzip -cd", shQuote(sed_path))

sampling_events <- fread(
  cmd = cmd,
  sep = "\t",
  quote = "",
  showProgress = TRUE,
  select = sed_cols_to_read
)

cat("\nSampling rows read:\n")
print(nrow(sampling_events))

# ------------------------------------------------------------
# 5. Rename columns
# ------------------------------------------------------------

rename_old <- c(
  "SAMPLING EVENT IDENTIFIER",
  "OBSERVATION DATE",
  "LATITUDE",
  "LONGITUDE",
  "ALL SPECIES REPORTED",
  "DURATION MINUTES",
  "EFFORT DISTANCE KM",
  "NUMBER OBSERVERS"
)

rename_new <- c(
  "checklist_id",
  "observation_date",
  "latitude",
  "longitude",
  "all_species_reported",
  "duration_minutes",
  "effort_distance_km",
  "number_observers"
)

if ("LOCALITY ID" %in% names(sampling_events)) {
  rename_old <- c(rename_old, "LOCALITY ID")
  rename_new <- c(rename_new, "locality_id")
}

if ("LOCALITY" %in% names(sampling_events)) {
  rename_old <- c(rename_old, "LOCALITY")
  rename_new <- c(rename_new, "locality")
}

setnames(
  sampling_events,
  old = rename_old,
  new = rename_new
)

# If the SED file did not contain locality fields, create them as NA
# so downstream scripts still have stable column names.
if (!"locality_id" %in% names(sampling_events)) {
  sampling_events[, locality_id := NA_character_]
}

if (!"locality" %in% names(sampling_events)) {
  sampling_events[, locality := NA_character_]
}

cat("\nLocality fields after renaming:\n")
print(
  grep(
    "locality|hotspot|location",
    names(sampling_events),
    ignore.case = TRUE,
    value = TRUE
  )
)

# ------------------------------------------------------------
# 6. Apply table filters
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
    number_observers = as.numeric(number_observers),
    locality_id = as.character(locality_id),
    locality = as.character(locality)
  )

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
    checklist_id,
    locality_id,
    locality,
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
# 9. Save hotspot/location list for GEE
# ------------------------------------------------------------

dryland_hotspots_for_GEE <- dryland_city_checklists_table %>%
  filter(!is.na(locality_id), locality_id != "") %>%
  group_by(
    city_id,
    NAME10,
    city_ascii,
    state_name,
    county_name,
    dryland_percent,
    locality_id,
    locality
  ) %>%
  summarise(
    n_checklists = n_distinct(checklist_id),
    first_date = min(observation_date, na.rm = TRUE),
    last_date = max(observation_date, na.rm = TRUE),
    mean_latitude = mean(latitude, na.rm = TRUE),
    mean_longitude = mean(longitude, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(state_name, city_ascii, locality)

write_csv(
  dryland_hotspots_for_GEE,
  hotspot_list_out
)

# ------------------------------------------------------------
# 10. Save summary by city
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
    n_unique_locality_ids = n_distinct(locality_id[!is.na(locality_id) & locality_id != ""]),
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
# 11. Final checks
# ------------------------------------------------------------

cat("\nDryland city checklists retained:\n")
print(nrow(dryland_city_checklists_table))

cat("\nUnique checklist IDs retained:\n")
print(nrow(dryland_checklist_ids))

cat("\nRows with locality_id:\n")
print(sum(!is.na(dryland_city_checklists_table$locality_id) & dryland_city_checklists_table$locality_id != ""))

cat("\nUnique locality IDs retained:\n")
print(
  dryland_city_checklists_table %>%
    filter(!is.na(locality_id), locality_id != "") %>%
    distinct(locality_id) %>%
    nrow()
)

cat("\nDryland cities with at least one checklist:\n")
dryland_city_checklists_table %>%
  distinct(city_id) %>%
  nrow() %>%
  print()

cat("\nTop dryland cities by checklist count:\n")
print(head(dryland_city_summary, 25))

cat("\nOutput files written:\n")
print(checklists_out)
print(checklist_ids_out)
print(summary_by_city_out)
print(hotspot_list_out)
