# ============================================================
# 05_Filter_eBird_Observations_to_Dryland_Cities.R
#
# Purpose:
# Use the dryland-city checklist IDs from File 04 to filter the
# full observation-level eBird Basic Dataset.
#
# Outputs:
# 1. Observation records for dryland-city checklists
# 2. Observation records joined to city/checklist metadata
# 3. Species richness per checklist
# 4. Species richness per city per year
#
# Important:
# This script streams the full EBD observation file directly from
# the .tar archive. It does NOT extract the 129 GB .gz file first.
# ============================================================

library(data.table)
library(dplyr)
library(readr)

# ------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------

# Original eBird TAR file
tar_path <- "Data/Raw/eBird/ebd_US_smp_relApr-2026.tar"

# Full observation file inside the TAR archive
ebd_member <- "ebd_US_relApr-2026.txt.gz"

# Outputs from File 04
checklists_path <- "Data/Processed/eBird/dryland_city_checklists_2015_2025.csv"

checklist_ids_path <- "Data/Processed/eBird/dryland_city_checklist_ids_2015_2025.csv"

# Output folder
ebird_out_dir <- "Data/Processed/eBird"
dir.create(ebird_out_dir, recursive = TRUE, showWarnings = FALSE)

# Observation-level outputs
observations_tsv_out <- file.path(
  ebird_out_dir,
  "dryland_city_observations_2015_2025.tsv"
)

observations_with_city_out <- file.path(
  ebird_out_dir,
  "dryland_city_observations_with_city_2015_2025.tsv"
)

# Richness outputs
checklist_richness_out <- file.path(
  ebird_out_dir,
  "dryland_city_checklist_species_richness_2015_2025.csv"
)

city_year_richness_out <- file.path(
  ebird_out_dir,
  "dryland_city_year_species_richness_2015_2025.csv"
)

# AWK script and log
awk_script_path <- file.path(
  ebird_out_dir,
  "filter_ebd_observations_by_checklist.awk"
)

awk_log_path <- file.path(
  ebird_out_dir,
  "filter_ebd_observations_by_checklist_log.txt"
)

# Remove old outputs if rerunning
unlink(
  c(
    observations_tsv_out,
    observations_with_city_out,
    checklist_richness_out,
    city_year_richness_out,
    awk_script_path,
    awk_log_path
  ),
  force = TRUE
)

# ------------------------------------------------------------
# 2. Read checklist IDs and checklist metadata from File 04
# ------------------------------------------------------------

dryland_checklist_ids <- fread(checklist_ids_path)
dryland_checklists <- fread(checklists_path)

cat("\nDryland checklist IDs loaded:\n")
print(nrow(dryland_checklist_ids))

cat("\nDryland checklist metadata rows loaded:\n")
print(nrow(dryland_checklists))

# ------------------------------------------------------------
# 3. Confirm EBD file exists inside TAR
# ------------------------------------------------------------

cat("\nFiles inside eBird TAR:\n")

tar_list_cmd <- paste(
  "/usr/bin/tar -tf",
  shQuote(tar_path)
)

tar_files <- system(
  tar_list_cmd,
  intern = TRUE
)

print(tar_files)

if (!ebd_member %in% tar_files) {
  stop("The EBD member file was not found inside the TAR archive.")
}

# ------------------------------------------------------------
# 4. Write AWK script to filter EBD observations by checklist ID
# ------------------------------------------------------------

# This AWK script:
# - reads checklist IDs from File 04
# - streams the full observation-level EBD
# - keeps only observation rows whose checklist ID is in the dryland-city list
# - writes a smaller observation-level TSV file
#
# Notes:
# EBD uses PROTOCOL NAME, not PROTOCOL TYPE.

awk_script <- '
BEGIN {
  FS = "\t"
  OFS = "\t"
}

NR == FNR {
  if (FNR == 1 && $1 == "checklist_id") {
    next
  }
  ids[$1] = 1
  next
}

FNR == 1 {
  for (i = 1; i <= NF; i++) {
    col[$i] = i
  }

  id_col       = col["SAMPLING EVENT IDENTIFIER"]
  date_col     = col["OBSERVATION DATE"]
  common_col   = col["COMMON NAME"]
  sci_col      = col["SCIENTIFIC NAME"]
  tax_col      = col["TAXONOMIC ORDER"]
  category_col = col["CATEGORY"]
  count_col    = col["OBSERVATION COUNT"]
  lat_col      = col["LATITUDE"]
  lon_col      = col["LONGITUDE"]
  protocol_col = col["PROTOCOL NAME"]
  approved_col = col["APPROVED"]
  reviewed_col = col["REVIEWED"]

  missing = 0

  if (id_col == 0) {
    print "Missing SAMPLING EVENT IDENTIFIER" > "/dev/stderr"
    missing = 1
  }
  if (date_col == 0) {
    print "Missing OBSERVATION DATE" > "/dev/stderr"
    missing = 1
  }
  if (common_col == 0) {
    print "Missing COMMON NAME" > "/dev/stderr"
    missing = 1
  }
  if (sci_col == 0) {
    print "Missing SCIENTIFIC NAME" > "/dev/stderr"
    missing = 1
  }
  if (count_col == 0) {
    print "Missing OBSERVATION COUNT" > "/dev/stderr"
    missing = 1
  }

  if (missing == 1) {
    print "Columns found in EBD file:" > "/dev/stderr"
    for (i = 1; i <= NF; i++) {
      print i, $i > "/dev/stderr"
    }
    exit 1
  }

  print "checklist_id",
        "observation_date",
        "common_name",
        "scientific_name",
        "taxonomic_order",
        "category",
        "observation_count",
        "latitude",
        "longitude",
        "protocol_name",
        "approved",
        "reviewed"

  next
}

{
  checklist_id = $id_col

  if (checklist_id in ids) {
    print checklist_id,
          $date_col,
          $common_col,
          $sci_col,
          (tax_col > 0 ? $tax_col : ""),
          (category_col > 0 ? $category_col : ""),
          $count_col,
          (lat_col > 0 ? $lat_col : ""),
          (lon_col > 0 ? $lon_col : ""),
          (protocol_col > 0 ? $protocol_col : ""),
          (approved_col > 0 ? $approved_col : ""),
          (reviewed_col > 0 ? $reviewed_col : "")
  }
}
'

writeLines(awk_script, awk_script_path)

# ------------------------------------------------------------
# 5. Stream full EBD from TAR and filter by checklist ID
# ------------------------------------------------------------

# Correct macOS TAR syntax:
# /usr/bin/tar -O -x -f archive.tar member.gz
#
# This streams the .gz file from inside the TAR, decompresses it,
# and filters rows with AWK.

filter_cmd <- paste(
  "/usr/bin/tar -O -x -f",
  shQuote(tar_path),
  shQuote(ebd_member),
  "| gzip -cd",
  "| awk -f",
  shQuote(awk_script_path),
  shQuote(checklist_ids_path),
  "-",
  ">",
  shQuote(observations_tsv_out),
  "2>",
  shQuote(awk_log_path)
)

cat("\nFiltering full EBD observation file. This may take several hours...\n")
cat(filter_cmd, "\n")

exit_code <- system(filter_cmd)

cat("\nAWK exit code:\n")
print(exit_code)

cat("\nAWK log:\n")
if (file.exists(awk_log_path)) {
  awk_log <- readLines(awk_log_path, warn = FALSE)
  if (length(awk_log) > 0) {
    print(awk_log)
  } else {
    cat("No AWK log messages.\n")
  }
}

if (!file.exists(observations_tsv_out)) {
  stop("The filtered observation file was not created.")
}

obs_file_info <- file.info(observations_tsv_out)

cat("\nFiltered observation file info:\n")
print(obs_file_info[, c("size", "mtime")])

if (is.na(obs_file_info$size) || obs_file_info$size == 0) {
  stop("The filtered observation file exists but is empty. Check the AWK log above.")
}

# ------------------------------------------------------------
# 6. Read filtered observations
# ------------------------------------------------------------

observations <- fread(
  observations_tsv_out,
  sep = "\t",
  quote = "",
  showProgress = TRUE
)

cat("\nFiltered observation rows loaded:\n")
print(nrow(observations))

cat("\nUnique checklists with observations:\n")
print(uniqueN(observations$checklist_id))

# ------------------------------------------------------------
# 7. Join observations to city/checklist metadata
# ------------------------------------------------------------

checklist_metadata <- dryland_checklists %>%
  as_tibble() %>%
  select(
    checklist_id,
    city_id,
    NAME10,
    city_ascii,
    state_name,
    county_name,
    dryland_percent,
    checklist_latitude = latitude,
    checklist_longitude = longitude,
    duration_minutes,
    effort_distance_km,
    number_observers
  ) %>%
  as.data.table()

observations <- as.data.table(observations)

observations_with_city <- merge(
  observations,
  checklist_metadata,
  by = "checklist_id",
  all.x = TRUE,
  all.y = FALSE
)

observations_with_city[
  ,
  observation_year := as.integer(substr(observation_date, 1, 4))
]

cat("\nObservation rows after city join:\n")
print(nrow(observations_with_city))

cat("\nUnique cities with observations:\n")
print(uniqueN(observations_with_city$city_id))

fwrite(
  observations_with_city,
  observations_with_city_out,
  sep = "\t"
)

# ------------------------------------------------------------
# 8. Species richness per checklist
# ------------------------------------------------------------

# checklist_species_richness:
# number of unique scientific names where EBD category == "species"
#
# checklist_taxa_richness_all:
# number of unique scientific names across all EBD categories,
# including hybrids, slash taxa, spuhs, etc.

checklist_richness <- observations_with_city[
  ,
  .(
    checklist_species_richness = uniqueN(scientific_name[category == "species"]),
    checklist_taxa_richness_all = uniqueN(scientific_name),
    n_observation_rows = .N
  ),
  by = .(
    checklist_id,
    observation_date,
    observation_year,
    city_id,
    NAME10,
    city_ascii,
    state_name,
    county_name,
    dryland_percent,
    checklist_latitude,
    checklist_longitude,
    duration_minutes,
    effort_distance_km,
    number_observers,
    protocol_name
  )
]

fwrite(
  checklist_richness,
  checklist_richness_out
)

# ------------------------------------------------------------
# 9. Species richness per city per year
# ------------------------------------------------------------

# City-year species richness:
# total number of unique species observed across all retained
# checklists in a city during a given year.

city_year_species <- observations_with_city[
  category == "species",
  .(
    city_year_species_richness = as.integer(uniqueN(scientific_name)),
    n_species_observation_rows = as.integer(.N)
  ),
  by = .(
    city_id,
    NAME10,
    city_ascii,
    state_name,
    county_name,
    dryland_percent,
    observation_year
  )
]

city_year_taxa_all <- observations_with_city[
  ,
  .(
    city_year_taxa_richness_all = as.integer(uniqueN(scientific_name)),
    n_total_observation_rows = as.integer(.N)
  ),
  by = .(
    city_id,
    NAME10,
    city_ascii,
    state_name,
    county_name,
    dryland_percent,
    observation_year
  )
]

city_year_effort <- checklist_richness[
  ,
  .(
    n_checklists = as.integer(uniqueN(checklist_id)),
    n_unique_locations = as.integer(uniqueN(paste(checklist_latitude, checklist_longitude))),
    mean_checklist_species_richness = as.numeric(mean(checklist_species_richness, na.rm = TRUE)),
    median_checklist_species_richness = as.numeric(median(checklist_species_richness, na.rm = TRUE)),
    mean_duration_minutes = as.numeric(mean(duration_minutes, na.rm = TRUE)),
    mean_effort_distance_km = as.numeric(mean(effort_distance_km, na.rm = TRUE))
  ),
  by = .(
    city_id,
    NAME10,
    city_ascii,
    state_name,
    county_name,
    dryland_percent,
    observation_year
  )
]

city_year_summary <- Reduce(
  function(x, y) merge(
    x,
    y,
    by = c(
      "city_id",
      "NAME10",
      "city_ascii",
      "state_name",
      "county_name",
      "dryland_percent",
      "observation_year"
    ),
    all = TRUE
  ),
  list(
    city_year_effort,
    city_year_species,
    city_year_taxa_all
  )
)

setorder(
  city_year_summary,
  state_name,
  city_ascii,
  observation_year
)

fwrite(
  city_year_summary,
  city_year_richness_out
)

# ------------------------------------------------------------
# 10. Final checks
# ------------------------------------------------------------

cat("\nFinal outputs created:\n")
print(observations_tsv_out)
print(observations_with_city_out)
print(checklist_richness_out)
print(city_year_richness_out)

cat("\nObservation rows retained:\n")
print(nrow(observations_with_city))

cat("\nUnique checklists retained with observations:\n")
print(uniqueN(observations_with_city$checklist_id))

cat("\nUnique cities retained with observations:\n")
print(uniqueN(observations_with_city$city_id))

cat("\nCity-year richness rows:\n")
print(nrow(city_year_summary))

cat("\nTop city-years by species richness:\n")
print(
  head(
    city_year_summary[
      order(-city_year_species_richness),
      .(
        NAME10,
        city_ascii,
        state_name,
        observation_year,
        city_year_species_richness,
        n_checklists,
        n_unique_locations,
        mean_checklist_species_richness
      )
    ],
    25
  )
)