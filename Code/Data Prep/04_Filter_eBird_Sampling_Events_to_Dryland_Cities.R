# ============================================================
# 04_Filter_eBird_Sampling_Events_to_Dryland_Cities_WITH_HOTSPOTS.R
#
# Purpose:
# Filter the July 2026 eBird Sampling Event Data to complete
# checklists from 2015-2025 that fall inside dryland city polygons.
#
# EFFICIENCY APPROACH:
#
# The full U.S. Sampling Event Data is NOT loaded into R.
#
# Instead:
#   1. gzip streams the compressed national file.
#   2. AWK filters each row immediately.
#   3. Only potentially relevant sampling events are written
#      to a much smaller intermediate TSV.
#   4. R reads that reduced file.
#   5. sf performs the exact dryland-city polygon join.
#
# Filters applied before data enter R:
#   - observation date = 2015-01-01 through 2025-12-31
#   - complete checklists only
#   - one observer
#   - duration = 5-240 minutes
#   - valid coordinates
#   - coordinates inside broad dryland-city bounding box
#
# Exact city membership is determined afterward using sf.
#
# IMPORTANT:
# Spatial geometry comes from the Script 02 shapefile.
# City metadata come from the authoritative Script 02 CSV.
#
# This avoids ESRI Shapefile field-name truncation and ensures
# that dryland_percent is correctly retained.
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(sf)
library(dplyr)
library(readr)
library(data.table)


# ------------------------------------------------------------
# 2. Run settings
# ------------------------------------------------------------

# FALSE:
# Reuse the reduced intermediate Sampling Event Data file if it
# already exists.
#
# TRUE:
# Force a complete re-scan of the compressed U.S. Sampling
# Event Data.
#
# The July 2026 intermediate file has already been successfully
# created, so leave this FALSE.

rebuild_stream_subset <- FALSE


# Number of fread threads used for the reduced file.

fread_threads <- 4


# Floating-point tolerance used for percentage checks.
#
# This prevents values such as 100.00000000000001 from being
# incorrectly treated as biologically impossible.

percent_tolerance <- 1e-8


# ------------------------------------------------------------
# 3. File paths
# ------------------------------------------------------------

# July 2026 Sampling Event Data extracted by Script 03

sed_path <-
  "Data/Raw/eBird/Unpacked/ebd_US_relJul-2026_sampling.txt.gz"


# Spatial dryland-city polygons from Script 02

dryland_cities_path <-
  "Data/Processed/Shape Files/Cities/dryland_cities_50pct.shp"


# Authoritative dryland-city metadata from Script 02

dryland_city_table_path <-
  "Data/Processed/Tables/dryland_cities_50pct.csv"


# Main eBird output folder

ebird_out_dir <- "Data/Processed/eBird"

dir.create(
  ebird_out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Intermediate folder

intermediate_dir <- file.path(
  ebird_out_dir,
  "Intermediate"
)

dir.create(
  intermediate_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Reduced Sampling Event Data file

stream_subset_path <- file.path(
  intermediate_dir,
  "dryland_bbox_sampling_events_2015_2025.tsv"
)


# AWK filtering script

awk_script_path <- file.path(
  intermediate_dir,
  "filter_sampling_events.awk"
)


# Final Script 04 outputs

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


# ------------------------------------------------------------
# 4. Confirm required inputs
# ------------------------------------------------------------

if (!file.exists(sed_path)) {
  
  stop(
    "Sampling Event Data file not found:\n",
    sed_path
  )
}


if (!file.exists(dryland_cities_path)) {
  
  stop(
    "Dryland city shapefile not found:\n",
    dryland_cities_path
  )
}


if (!file.exists(dryland_city_table_path)) {
  
  stop(
    "Dryland city metadata table not found:\n",
    dryland_city_table_path
  )
}


sed_size <- file.info(
  sed_path
)$size


if (
  is.na(sed_size) ||
  sed_size == 0
) {
  
  stop(
    "Sampling Event Data file exists but is empty:\n",
    sed_path
  )
}


cat("\n============================================================\n")
cat("INPUT FILES\n")
cat("============================================================\n\n")


cat("Sampling Event Data:\n")
print(sed_path)


cat("\nSampling Event Data size (GB):\n")

print(
  round(
    sed_size / 1024^3,
    3
  )
)


cat("\nDryland city shapefile:\n")
print(dryland_cities_path)


cat("\nDryland city metadata table:\n")
print(dryland_city_table_path)


# ------------------------------------------------------------
# 5. Read authoritative city metadata
# ------------------------------------------------------------

# All descriptive city attributes come from the CSV created
# directly by Script 02.
#
# This avoids relying on truncated ESRI Shapefile field names.

dryland_city_metadata <- read_csv(
  dryland_city_table_path,
  show_col_types = FALSE
) %>%
  dplyr::select(
    city_id,
    NAME10,
    city_ascii,
    state_name,
    county_name,
    dryland_percent
  )


cat("\n============================================================\n")
cat("DRYLAND CITY METADATA\n")
cat("============================================================\n\n")


cat("Rows in authoritative city table:\n")
print(nrow(dryland_city_metadata))


cat("\nUnique city IDs:\n")
print(n_distinct(dryland_city_metadata$city_id))


cat("\nRaw dryland-percent range:\n")

print(
  range(
    dryland_city_metadata$dryland_percent,
    na.rm = TRUE
  )
)


# ------------------------------------------------------------
# 6. QC and clean dryland_percent
# ------------------------------------------------------------

# Check for missing percentages.

missing_dryland_percent <- sum(
  is.na(
    dryland_city_metadata$dryland_percent
  )
)


cat("\nMissing dryland-percent values:\n")
print(missing_dryland_percent)


if (missing_dryland_percent > 0) {
  
  stop(
    "Missing dryland_percent values detected in authoritative city metadata."
  )
}


# city_id must be unique.

if (
  n_distinct(dryland_city_metadata$city_id) !=
  nrow(dryland_city_metadata)
) {
  
  stop(
    "city_id is not unique in the authoritative city metadata table."
  )
}


# Count values that are genuinely outside the expected
# 50-100% range after allowing tiny floating-point deviations.

invalid_dryland_percent_metadata <- sum(
  
  dryland_city_metadata$dryland_percent <
    (50 - percent_tolerance) |
    
    dryland_city_metadata$dryland_percent >
    (100 + percent_tolerance),
  
  na.rm = TRUE
)


cat("\nTruly invalid dryland-percent values:\n")
print(invalid_dryland_percent_metadata)


if (invalid_dryland_percent_metadata > 0) {
  
  stop(
    "Invalid dryland_percent values detected in authoritative city metadata."
  )
}


# Clean tiny floating-point artifacts around 100%.
#
# Example:
# 100.00000000000001 becomes exactly 100.

dryland_city_metadata <- dryland_city_metadata %>%
  mutate(
    
    dryland_percent = case_when(
      
      dryland_percent > 100 &
        dryland_percent <=
        (100 + percent_tolerance) ~ 100,
      
      TRUE ~ dryland_percent
    )
  )


cat("\nCleaned dryland-percent range:\n")

print(
  range(
    dryland_city_metadata$dryland_percent,
    na.rm = TRUE
  )
)


# ------------------------------------------------------------
# 7. Read city geometry
# ------------------------------------------------------------

# Only city_id and geometry are taken from the shapefile.
#
# All descriptive fields come from the authoritative CSV.

dryland_city_geometry <- st_read(
  dryland_cities_path
) %>%
  st_make_valid() %>%
  st_transform(4326) %>%
  dplyr::select(
    city_id,
    geometry
  )


cat("\n============================================================\n")
cat("DRYLAND CITY GEOMETRY\n")
cat("============================================================\n\n")


cat("Dryland city geometries loaded:\n")
print(nrow(dryland_city_geometry))


cat("\nUnique city IDs in spatial layer:\n")
print(n_distinct(dryland_city_geometry$city_id))


if (
  n_distinct(dryland_city_geometry$city_id) !=
  nrow(dryland_city_geometry)
) {
  
  stop(
    "city_id is not unique in the dryland city spatial layer."
  )
}


# ------------------------------------------------------------
# 8. Confirm spatial and metadata city IDs match
# ------------------------------------------------------------

cities_missing_metadata <- dryland_city_geometry %>%
  st_drop_geometry() %>%
  anti_join(
    dryland_city_metadata %>%
      distinct(city_id),
    by = "city_id"
  )


cities_missing_geometry <- dryland_city_metadata %>%
  distinct(city_id) %>%
  anti_join(
    dryland_city_geometry %>%
      st_drop_geometry() %>%
      distinct(city_id),
    by = "city_id"
  )


cat("\nSpatial cities without metadata:\n")
print(nrow(cities_missing_metadata))


cat("\nMetadata cities without geometry:\n")
print(nrow(cities_missing_geometry))


if (
  nrow(cities_missing_metadata) > 0 ||
  nrow(cities_missing_geometry) > 0
) {
  
  stop(
    "City IDs do not match between the Script 02 shapefile and CSV."
  )
}


# ------------------------------------------------------------
# 9. Join authoritative metadata to city geometry
# ------------------------------------------------------------

dryland_cities <- dryland_city_geometry %>%
  left_join(
    dryland_city_metadata,
    by = "city_id"
  ) %>%
  dplyr::select(
    city_id,
    NAME10,
    city_ascii,
    state_name,
    county_name,
    dryland_percent,
    geometry
  )


cat("\n============================================================\n")
cat("FINAL DRYLAND CITY LAYER\n")
cat("============================================================\n\n")


cat("Dryland city polygons:\n")
print(nrow(dryland_cities))


cat("\nUnique city IDs:\n")
print(n_distinct(dryland_cities$city_id))


cat("\nDryland-percent range:\n")

print(
  range(
    dryland_cities$dryland_percent,
    na.rm = TRUE
  )
)


if (nrow(dryland_cities) != 653) {
  
  stop(
    "Expected 653 dryland cities but found ",
    nrow(dryland_cities),
    "."
  )
}


invalid_dryland_percent_spatial <- sum(
  
  dryland_cities$dryland_percent <
    (50 - percent_tolerance) |
    
    dryland_cities$dryland_percent >
    (100 + percent_tolerance),
  
  na.rm = TRUE
)


cat("\nInvalid dryland-percent values after spatial join:\n")
print(invalid_dryland_percent_spatial)


if (invalid_dryland_percent_spatial > 0) {
  
  stop(
    "Invalid dryland_percent values detected after metadata join."
  )
}


# Geometry object no longer needed separately.

rm(
  dryland_city_geometry
)

gc()


# ------------------------------------------------------------
# 10. Create broad dryland-city bounding box
# ------------------------------------------------------------

# This is used only as an efficient coarse geographic filter.
#
# Exact point-in-polygon membership is tested later.

bb <- st_bbox(
  dryland_cities
)


xmin <- as.numeric(
  bb["xmin"]
)

xmax <- as.numeric(
  bb["xmax"]
)

ymin <- as.numeric(
  bb["ymin"]
)

ymax <- as.numeric(
  bb["ymax"]
)


cat("\nDryland city bounding box:\n")
print(bb)


# ------------------------------------------------------------
# 11. Confirm Sampling Event Data structure
# ------------------------------------------------------------

# Read only the header of the compressed national file.

sed_header <- fread(
  cmd = paste(
    "gzip -cd",
    shQuote(sed_path)
  ),
  sep = "\t",
  quote = "",
  nrows = 0,
  showProgress = FALSE
)


sed_names <- names(
  sed_header
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


missing_required_sed_cols <- setdiff(
  required_sed_cols,
  sed_names
)


if (length(missing_required_sed_cols) > 0) {
  
  stop(
    "Sampling Event Data is missing required columns:\n",
    paste(
      missing_required_sed_cols,
      collapse = "\n"
    )
  )
}


cat("\nRequired Sampling Event Data columns confirmed.\n")


cat("\nLocality-related columns available:\n")

print(
  grep(
    "LOCALITY|LOCATION|HOTSPOT",
    sed_names,
    ignore.case = TRUE,
    value = TRUE
  )
)


# ------------------------------------------------------------
# 12. Stream and filter Sampling Event Data
# ------------------------------------------------------------
#
# This is the expensive national scan.
#
# AWK processes one row at a time and discards irrelevant rows
# before they enter R.
#
# This section is skipped when:
#
# rebuild_stream_subset <- FALSE
#
# AND the intermediate file already exists.
# ------------------------------------------------------------

if (
  rebuild_stream_subset ||
  !file.exists(stream_subset_path)
) {
  
  
  cat("\n============================================================\n")
  cat("STREAMING FILTER OF U.S. SAMPLING EVENT DATA\n")
  cat("============================================================\n\n")
  
  
  # Remove any existing/partial intermediate file if we are
  # explicitly rebuilding it.
  
  if (file.exists(stream_subset_path)) {
    
    cat("Removing previous intermediate subset:\n")
    print(stream_subset_path)
    
    file.remove(
      stream_subset_path
    )
  }
  
  
  # ----------------------------------------------------------
  # 12.1 Write AWK filtering script
  # ----------------------------------------------------------
  
  awk_script <- '
BEGIN {
  FS = "\t"
  OFS = "\t"
}

FNR == 1 {

  for (i = 1; i <= NF; i++) {
    col[$i] = i
  }

  id_col = col["SAMPLING EVENT IDENTIFIER"]
  date_col = col["OBSERVATION DATE"]
  lat_col = col["LATITUDE"]
  lon_col = col["LONGITUDE"]
  complete_col = col["ALL SPECIES REPORTED"]
  duration_col = col["DURATION MINUTES"]
  distance_col = col["EFFORT DISTANCE KM"]
  observer_col = col["NUMBER OBSERVERS"]

  locality_id_col = col["LOCALITY ID"]
  locality_col = col["LOCALITY"]

  print "checklist_id", "observation_date", "latitude", "longitude", "all_species_reported", "duration_minutes", "effort_distance_km", "number_observers", "locality_id", "locality"

  next
}

{

  date_val = $date_col
  lat_val = $lat_col
  lon_val = $lon_col
  duration_val = $duration_col
  observer_val = $observer_col
  complete_val = $complete_col

  if (date_val < start_date || date_val > end_date) next
  if (complete_val != 1) next
  if (observer_val != 1) next
  if (duration_val < 5 || duration_val > 240) next

  if (lat_val == "") next
  if (lon_val == "") next

  if ((lat_val + 0) < ymin || (lat_val + 0) > ymax) next
  if ((lon_val + 0) < xmin || (lon_val + 0) > xmax) next

  locality_id_val = ""
  locality_val = ""
  distance_val = ""

  if (locality_id_col > 0) {
    locality_id_val = $locality_id_col
  }

  if (locality_col > 0) {
    locality_val = $locality_col
  }

  if (distance_col > 0) {
    distance_val = $distance_col
  }

  print $id_col, date_val, lat_val, lon_val, complete_val, duration_val, distance_val, observer_val, locality_id_val, locality_val
}
'


writeLines(
  awk_script,
  awk_script_path
)


# ----------------------------------------------------------
# 12.2 Build streaming command
# ----------------------------------------------------------

# LC_ALL=C reduces unnecessary locale-aware processing.
#
# nice -n 10 lowers CPU scheduling priority so the Mac
# remains more responsive during the scan.

filter_cmd <- paste(
  
  "LC_ALL=C",
  
  "nice -n 10",
  
  "gzip -cd",
  shQuote(sed_path),
  
  "|",
  
  "LC_ALL=C",
  
  "nice -n 10",
  
  "awk",
  
  "-v",
  paste0("xmin=", xmin),
  
  "-v",
  paste0("xmax=", xmax),
  
  "-v",
  paste0("ymin=", ymin),
  
  "-v",
  paste0("ymax=", ymax),
  
  "-v",
  "start_date=2015-01-01",
  
  "-v",
  "end_date=2025-12-31",
  
  "-f",
  shQuote(awk_script_path),
  
  ">",
  shQuote(stream_subset_path)
)


cat("Scanning compressed U.S. Sampling Event Data.\n")
cat("The complete compressed file must be read once.\n")
cat("Irrelevant rows will be discarded before they reach R.\n\n")


stream_start <- Sys.time()


exit_code <- system(
  filter_cmd
)


stream_end <- Sys.time()


cat("\nStreaming filter exit code:\n")
print(exit_code)


cat("\nStreaming filter elapsed time:\n")

print(
  stream_end - stream_start
)


# ----------------------------------------------------------
# 12.3 Confirm streaming filter succeeded
# ----------------------------------------------------------

if (exit_code != 0) {
  
  stop(
    "Streaming Sampling Event Data filter failed."
  )
}


if (!file.exists(stream_subset_path)) {
  
  stop(
    "Filtered Sampling Event Data file was not created."
  )
}


subset_size <- file.info(
  stream_subset_path
)$size


if (
  is.na(subset_size) ||
  subset_size == 0
) {
  
  stop(
    "Filtered Sampling Event Data file was created but is empty."
  )
}


cat("\nStreaming filter completed successfully.\n")


} else {
  
  
  cat("\n============================================================\n")
  cat("REUSING FILTERED SAMPLING EVENT DATA\n")
  cat("============================================================\n\n")
  
  
  cat("Existing intermediate subset:\n")
  print(stream_subset_path)
}


# ------------------------------------------------------------
# 13. Inspect reduced Sampling Event Data file
# ------------------------------------------------------------

if (!file.exists(stream_subset_path)) {
  
  stop(
    "Reduced Sampling Event Data file does not exist."
  )
}


stream_subset_size <- file.info(
  stream_subset_path
)$size


if (
  is.na(stream_subset_size) ||
  stream_subset_size == 0
) {
  
  stop(
    "Reduced Sampling Event Data file is empty."
  )
}


stream_subset_size_mb <-
  stream_subset_size / 1024^2


cat("\nFiltered Sampling Event Data size:\n")

print(
  paste0(
    round(
      stream_subset_size_mb,
      2
    ),
    " MB"
  )
)


# Confirm expected structure.

stream_header <- fread(
  stream_subset_path,
  sep = "\t",
  quote = "",
  nrows = 0,
  showProgress = FALSE
)


expected_stream_cols <- c(
  "checklist_id",
  "observation_date",
  "latitude",
  "longitude",
  "all_species_reported",
  "duration_minutes",
  "effort_distance_km",
  "number_observers",
  "locality_id",
  "locality"
)


missing_stream_cols <- setdiff(
  expected_stream_cols,
  names(stream_header)
)


if (length(missing_stream_cols) > 0) {
  
  stop(
    "Filtered intermediate file is missing expected columns:\n",
    paste(
      missing_stream_cols,
      collapse = "\n"
    )
  )
}


cat("\nFiltered intermediate file structure confirmed.\n")


# ------------------------------------------------------------
# 14. Read reduced Sampling Event Data into R
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("READING REDUCED DATA INTO R\n")
cat("============================================================\n\n")


read_start <- Sys.time()


sampling_events <- fread(
  stream_subset_path,
  sep = "\t",
  quote = "",
  showProgress = TRUE,
  nThread = fread_threads
)


read_end <- Sys.time()


cat("\nReduced file read time:\n")

print(
  read_end - read_start
)


cat("\nFiltered sampling-event rows read into R:\n")
print(nrow(sampling_events))


if (nrow(sampling_events) == 0) {
  
  stop(
    "Filtered Sampling Event Data contains zero rows."
  )
}


# ------------------------------------------------------------
# 15. Convert fields to appropriate types
# ------------------------------------------------------------

sampling_events <- sampling_events %>%
  mutate(
    
    checklist_id = as.character(
      checklist_id
    ),
    
    observation_date = as.Date(
      observation_date
    ),
    
    latitude = as.numeric(
      latitude
    ),
    
    longitude = as.numeric(
      longitude
    ),
    
    all_species_reported = as.numeric(
      all_species_reported
    ),
    
    duration_minutes = as.numeric(
      duration_minutes
    ),
    
    effort_distance_km = as.numeric(
      effort_distance_km
    ),
    
    number_observers = as.numeric(
      number_observers
    ),
    
    locality_id = as.character(
      locality_id
    ),
    
    locality = as.character(
      locality
    )
  )


# ------------------------------------------------------------
# 16. QC reduced Sampling Event Data
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FILTERED SAMPLING EVENT QC\n")
cat("============================================================\n\n")


cat("Rows:\n")
print(nrow(sampling_events))


cat("\nDate range:\n")

print(
  range(
    sampling_events$observation_date,
    na.rm = TRUE
  )
)


cat("\nDuration range:\n")

print(
  range(
    sampling_events$duration_minutes,
    na.rm = TRUE
  )
)


cat("\nObserver values:\n")

print(
  sort(
    unique(
      sampling_events$number_observers
    )
  )
)


cat("\nComplete-checklist values:\n")

print(
  sort(
    unique(
      sampling_events$all_species_reported
    )
  )
)


# ------------------------------------------------------------
# 17. Hard QC checks for filtered Sampling Event Data
# ------------------------------------------------------------

if (
  any(
    is.na(
      sampling_events$observation_date
    )
  )
) {
  
  stop(
    "Missing observation dates detected."
  )
}


if (
  min(
    sampling_events$observation_date,
    na.rm = TRUE
  ) < as.Date("2015-01-01") ||
  
  max(
    sampling_events$observation_date,
    na.rm = TRUE
  ) > as.Date("2025-12-31")
) {
  
  stop(
    "Observation dates outside the intended study period were detected."
  )
}


if (
  any(
    sampling_events$duration_minutes < 5 |
    sampling_events$duration_minutes > 240,
    na.rm = TRUE
  )
) {
  
  stop(
    "Sampling events outside the intended duration range were detected."
  )
}


if (
  any(
    sampling_events$number_observers != 1,
    na.rm = TRUE
  )
) {
  
  stop(
    "Sampling events with observer counts other than 1 were detected."
  )
}


if (
  any(
    sampling_events$all_species_reported != 1,
    na.rm = TRUE
  )
) {
  
  stop(
    "Incomplete checklists were detected."
  )
}


if (
  any(
    is.na(sampling_events$latitude) |
    is.na(sampling_events$longitude)
  )
) {
  
  stop(
    "Missing latitude or longitude values detected."
  )
}


# ------------------------------------------------------------
# 18. Convert sampling events to spatial points
# ------------------------------------------------------------

cat("\nConverting retained sampling events to spatial points...\n")


sampling_points <- st_as_sf(
  sampling_events,
  coords = c(
    "longitude",
    "latitude"
  ),
  crs = 4326,
  remove = FALSE
)


# Non-spatial version is no longer needed.

rm(
  sampling_events
)

gc()


# ------------------------------------------------------------
# 19. Exact spatial join to dryland city polygons
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("EXACT DRYLAND-CITY SPATIAL JOIN\n")
cat("============================================================\n\n")


spatial_start <- Sys.time()


dryland_city_checklists <- st_join(
  sampling_points,
  dryland_cities,
  join = st_within,
  left = FALSE
)


spatial_end <- Sys.time()


cat("Spatial join elapsed time:\n")

print(
  spatial_end - spatial_start
)


cat("\nChecklist rows inside dryland city polygons:\n")
print(nrow(dryland_city_checklists))


if (nrow(dryland_city_checklists) == 0) {
  
  stop(
    "No sampling events fell inside the dryland city polygons."
  )
}


rm(
  sampling_points
)

gc()


# ------------------------------------------------------------
# 20. Check whether any checklist maps to multiple cities
# ------------------------------------------------------------

checklist_city_qc <- dryland_city_checklists %>%
  st_drop_geometry() %>%
  distinct(
    checklist_id,
    city_id
  ) %>%
  count(
    checklist_id,
    name = "n_cities"
  ) %>%
  filter(
    n_cities > 1
  )


cat("\nChecklists assigned to more than one city polygon:\n")
print(nrow(checklist_city_qc))


if (nrow(checklist_city_qc) > 0) {
  
  stop(
    "Some checklists were assigned to more than one city polygon. ",
    "Inspect overlapping city polygons before continuing."
  )
}


# ------------------------------------------------------------
# 21. Create final checklist table
# ------------------------------------------------------------

dryland_city_checklists_table <- dryland_city_checklists %>%
  st_drop_geometry() %>%
  arrange(
    state_name,
    city_ascii,
    observation_date
  ) %>%
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
  ) %>%
  as_tibble()


rm(
  dryland_city_checklists
)

gc()


# ------------------------------------------------------------
# 22. Checklist-level QC
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CHECKLIST-LEVEL QC\n")
cat("============================================================\n\n")


cat("Checklist rows after spatial join:\n")
print(nrow(dryland_city_checklists_table))


cat("\nUnique checklist IDs:\n")

unique_checklists <- n_distinct(
  dryland_city_checklists_table$checklist_id
)

print(
  unique_checklists
)


cat("\nDuplicated checklist IDs:\n")

duplicate_checklists <- sum(
  duplicated(
    dryland_city_checklists_table$checklist_id
  )
)

print(
  duplicate_checklists
)


if (duplicate_checklists > 0) {
  
  stop(
    "Duplicated checklist IDs detected after spatial join."
  )
}


cat("\nDryland-percent range in checklist table:\n")

print(
  range(
    dryland_city_checklists_table$dryland_percent,
    na.rm = TRUE
  )
)


invalid_dryland_percent <- sum(
  
  dryland_city_checklists_table$dryland_percent <
    (50 - percent_tolerance) |
    
    dryland_city_checklists_table$dryland_percent >
    (100 + percent_tolerance),
  
  na.rm = TRUE
)


cat("\nRows with invalid dryland_percent:\n")
print(invalid_dryland_percent)


if (invalid_dryland_percent > 0) {
  
  stop(
    "Invalid dryland_percent values detected in checklist output."
  )
}


# ------------------------------------------------------------
# 23. Check city metadata against authoritative table
# ------------------------------------------------------------

city_metadata_qc <- dryland_city_checklists_table %>%
  distinct(
    city_id,
    NAME10,
    city_ascii,
    state_name,
    county_name,
    dryland_percent
  ) %>%
  left_join(
    
    dryland_city_metadata %>%
      rename(
        NAME10_master = NAME10,
        city_ascii_master = city_ascii,
        state_name_master = state_name,
        county_name_master = county_name,
        dryland_percent_master = dryland_percent
      ),
    
    by = "city_id"
  ) %>%
  filter(
    
    NAME10 != NAME10_master |
      
      city_ascii != city_ascii_master |
      
      state_name != state_name_master |
      
      county_name != county_name_master |
      
      abs(
        dryland_percent -
          dryland_percent_master
      ) > percent_tolerance
  )


cat("\nCities with metadata mismatches:\n")
print(nrow(city_metadata_qc))


if (nrow(city_metadata_qc) > 0) {
  
  print(
    city_metadata_qc,
    n = Inf
  )
  
  stop(
    "City metadata do not match the authoritative Script 02 table."
  )
}


# ------------------------------------------------------------
# 24. Create unique checklist ID table
# ------------------------------------------------------------

dryland_checklist_ids <- dryland_city_checklists_table %>%
  distinct(
    checklist_id
  )


# ------------------------------------------------------------
# 25. Create hotspot/location list for GEE
# ------------------------------------------------------------

dryland_hotspots_for_GEE <- dryland_city_checklists_table %>%
  filter(
    !is.na(locality_id),
    locality_id != ""
  ) %>%
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
    
    n_checklists = n_distinct(
      checklist_id
    ),
    
    first_date = min(
      observation_date,
      na.rm = TRUE
    ),
    
    last_date = max(
      observation_date,
      na.rm = TRUE
    ),
    
    mean_latitude = mean(
      latitude,
      na.rm = TRUE
    ),
    
    mean_longitude = mean(
      longitude,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  arrange(
    state_name,
    city_ascii,
    locality
  )


# ------------------------------------------------------------
# 26. Create summary by city
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
    
    n_checklists = n_distinct(
      checklist_id
    ),
    
    n_unique_locations = n_distinct(
      latitude,
      longitude
    ),
    
    n_unique_locality_ids = n_distinct(
      
      locality_id[
        !is.na(locality_id) &
          locality_id != ""
      ]
    ),
    
    first_date = min(
      observation_date,
      na.rm = TRUE
    ),
    
    last_date = max(
      observation_date,
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  arrange(
    desc(n_checklists)
  )


# ------------------------------------------------------------
# 27. Identify cities without retained eBird checklists
# ------------------------------------------------------------

cities_with_ebird <- dryland_city_summary %>%
  distinct(
    city_id
  )


cities_without_ebird <- dryland_city_metadata %>%
  anti_join(
    cities_with_ebird,
    by = "city_id"
  ) %>%
  arrange(
    state_name,
    city_ascii
  ) %>%
  as_tibble()


# ------------------------------------------------------------
# 28. Cross-output consistency QC
# ------------------------------------------------------------

summary_missing_cities <- dryland_city_checklists_table %>%
  distinct(
    city_id
  ) %>%
  anti_join(
    
    dryland_city_summary %>%
      distinct(
        city_id
      ),
    
    by = "city_id"
  )


hotspot_missing_cities <- dryland_city_checklists_table %>%
  distinct(
    city_id
  ) %>%
  anti_join(
    
    dryland_hotspots_for_GEE %>%
      distinct(
        city_id
      ),
    
    by = "city_id"
  )


# ------------------------------------------------------------
# 29. Final QC before saving
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SCRIPT 04 FINAL CHECKS\n")
cat("============================================================\n\n")


cat("Total dryland cities:\n")
print(nrow(dryland_cities))


cat("\nDryland city checklist rows retained:\n")
print(nrow(dryland_city_checklists_table))


cat("\nUnique checklist IDs retained:\n")
print(unique_checklists)


cat("\nDuplicated checklist IDs:\n")
print(duplicate_checklists)


cat("\nRows with locality_id:\n")

rows_with_locality_id <- sum(
  !is.na(
    dryland_city_checklists_table$locality_id
  ) &
    dryland_city_checklists_table$locality_id != ""
)

print(
  rows_with_locality_id
)


cat("\nMissing locality IDs:\n")

missing_locality_ids <- sum(
  is.na(
    dryland_city_checklists_table$locality_id
  ) |
    dryland_city_checklists_table$locality_id == ""
)

print(
  missing_locality_ids
)


cat("\nUnique locality IDs retained:\n")

unique_locality_ids <- dryland_city_checklists_table %>%
  filter(
    !is.na(locality_id),
    locality_id != ""
  ) %>%
  distinct(
    locality_id
  ) %>%
  nrow()

print(
  unique_locality_ids
)


cat("\nDryland-percent range:\n")

print(
  range(
    dryland_city_checklists_table$dryland_percent,
    na.rm = TRUE
  )
)


cat("\nRows with invalid dryland_percent:\n")
print(invalid_dryland_percent)


cat("\nCities with metadata mismatches:\n")
print(nrow(city_metadata_qc))


cat("\nDryland cities with at least one retained checklist:\n")
print(nrow(cities_with_ebird))


cat("\nDryland cities with no retained checklists:\n")
print(nrow(cities_without_ebird))


cat("\nCities with no retained checklists:\n")

print(
  cities_without_ebird,
  n = Inf
)


cat("\nCities in checklist table:\n")

checklist_city_count <- n_distinct(
  dryland_city_checklists_table$city_id
)

print(
  checklist_city_count
)


cat("\nCities in city-summary table:\n")

summary_city_count <- n_distinct(
  dryland_city_summary$city_id
)

print(
  summary_city_count
)


cat("\nCities in hotspot table:\n")

hotspot_city_count <- n_distinct(
  dryland_hotspots_for_GEE$city_id
)

print(
  hotspot_city_count
)


cat("\nChecklist cities absent from city summary:\n")
print(nrow(summary_missing_cities))


cat("\nChecklist cities absent from hotspot output:\n")
print(nrow(hotspot_missing_cities))


cat("\nTop 25 dryland cities by checklist count:\n")

print(
  as_tibble(
    head(
      dryland_city_summary,
      25
    )
  ),
  n = 25
)


# ------------------------------------------------------------
# 30. Hard final QC checks
# ------------------------------------------------------------

if (
  nrow(cities_with_ebird) +
  nrow(cities_without_ebird) !=
  nrow(dryland_city_metadata)
) {
  
  stop(
    "City accounting does not equal the master dryland-city dataset."
  )
}


if (
  checklist_city_count !=
  nrow(cities_with_ebird)
) {
  
  stop(
    "Checklist city count does not match cities_with_ebird."
  )
}


if (
  summary_city_count !=
  nrow(cities_with_ebird)
) {
  
  stop(
    "City-summary count does not match cities_with_ebird."
  )
}


if (nrow(summary_missing_cities) > 0) {
  
  stop(
    "Some checklist cities are absent from the city summary."
  )
}


if (nrow(hotspot_missing_cities) > 0) {
  
  stop(
    "Some checklist cities are absent from the hotspot output."
  )
}


if (
  nrow(dryland_checklist_ids) !=
  nrow(dryland_city_checklists_table)
) {
  
  stop(
    "Checklist row count and unique checklist count do not match."
  )
}


# ------------------------------------------------------------
# 31. Save final Script 04 outputs
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SAVING SCRIPT 04 OUTPUTS\n")
cat("============================================================\n\n")


write_csv(
  dryland_city_checklists_table,
  checklists_out
)


write_csv(
  dryland_checklist_ids,
  checklist_ids_out
)


write_csv(
  dryland_city_summary,
  summary_by_city_out
)


write_csv(
  dryland_hotspots_for_GEE,
  hotspot_list_out
)


# ------------------------------------------------------------
# 32. Confirm output files
# ------------------------------------------------------------

final_outputs <- c(
  checklists_out,
  checklist_ids_out,
  summary_by_city_out,
  hotspot_list_out
)


output_check <- data.frame(
  
  file = final_outputs,
  
  exists = file.exists(
    final_outputs
  ),
  
  size_mb = round(
    
    file.info(
      final_outputs
    )$size / 1024^2,
    
    2
  )
)


cat("\nFinal output files:\n")
print(output_check)


if (!all(output_check$exists)) {
  
  stop(
    "One or more final Script 04 outputs were not written successfully."
  )
}


# ------------------------------------------------------------
# 33. Final summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL SCRIPT 04 SUMMARY\n")
cat("============================================================\n\n")


cat("Master dryland cities:\n")
print(nrow(dryland_city_metadata))


cat("\nCities with retained eBird checklists:\n")
print(nrow(cities_with_ebird))


cat("\nCities without retained eBird checklists:\n")
print(nrow(cities_without_ebird))


cat("\nRetained checklist rows:\n")
print(nrow(dryland_city_checklists_table))


cat("\nUnique checklist IDs:\n")
print(unique_checklists)


cat("\nDuplicate checklist IDs:\n")
print(duplicate_checklists)


cat("\nUnique locality IDs:\n")
print(unique_locality_ids)


cat("\nInvalid dryland-percent rows:\n")
print(invalid_dryland_percent)


cat("\nCity metadata mismatches:\n")
print(nrow(city_metadata_qc))


cat("\nChecklist cities absent from city summary:\n")
print(nrow(summary_missing_cities))


cat("\nChecklist cities absent from hotspot output:\n")
print(nrow(hotspot_missing_cities))


# ------------------------------------------------------------
# 34. Complete
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SCRIPT 04 COMPLETE\n")
cat("============================================================\n\n")


cat("Reusable filtered Sampling Event Data:\n")
cat(stream_subset_path, "\n")


cat("\nCurrent setting:\n")
cat("rebuild_stream_subset <- FALSE\n")


cat("\nThe national Sampling Event Data will only be rescanned if:\n")
cat("- rebuild_stream_subset is changed to TRUE, or\n")
cat("- the intermediate filtered file no longer exists.\n")