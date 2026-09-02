# ============================================================
# 05_Filter_eBird_Observations_to_Dryland_Cities_WITH_HOTSPOTS.R
#
# Purpose:
# Use the QC-approved dryland-city checklist metadata from
# Script 04 to filter the July 2026 eBird Basic Dataset (EBD)
# and calculate checklist-level and city-year species richness.
#
# EFFICIENCY REVISION:
#
# The national EBD is streamed only once. Matching observations
# are written to ONE combined filtered file rather than hundreds
# of city files during the national scan.
#
# This removes the major macOS bottleneck from repeatedly doing:
#
#   print row >> city_file
#   close(city_file)
#
# for millions of EBD rows.
#
# Main workflow:
#
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(data.table)
library(stringr)


# ------------------------------------------------------------
# 2. macOS / RStudio warning control
# ------------------------------------------------------------

Sys.unsetenv("MallocStackLogging")
Sys.unsetenv("MallocStackLoggingNoCompact")

Sys.setenv(
  MallocStackLogging = "0",
  MallocStackLoggingNoCompact = "0"
)


# ------------------------------------------------------------
# 3. Run settings
# ------------------------------------------------------------

# eBird release used by this Script 05 run.

ebird_release <- "relJul-2026"


# Normally FALSE.
#
# TRUE forces the expensive U.S.-wide EBD filter to run again.

force_rebuild_filtered_ebd <- FALSE


# Normally FALSE.
#
# TRUE rebuilds the compact checklist/city-year summary files
# from the already-filtered EBD without rescanning the national
# EBD.

force_rebuild_ebd_summaries <- FALSE


# Normally FALSE.
#
# TRUE rewrites the per-city checklist-richness and city-year
# summary CSVs from the combined outputs.

rebuild_city_summaries <- FALSE


# ------------------------------------------------------------
# ONE-TIME LEGACY RECOVERY SETTING FOR THE CURRENT RUN
# ------------------------------------------------------------
#
# The previous Script 05 run shown in the console completed the
# national tar -> gzip -> awk pipeline with exit code 0 but then
# stalled before its R-side completion marker was written.
#
# TRUE allows this revised script to combine those completed 647
# raw city TSVs into the new single filtered EBD file instead of
# repeating the ~3.5-hour national scan.
#
# IMPORTANT:
# Use TRUE only because the previous console explicitly showed:
#
#   Raw EBD split exit code:
#   [1] 0
#
# Once this revised script has successfully created the new
# filtered-file completion marker, this recovery option is
# ignored on later runs.

recover_confirmed_legacy_split <- FALSE


# Raw per-city observation files are no longer required for the
# analysis. Leave FALSE for the efficient workflow.
#
# The old files can remain on disk for inspection/recovery, but
# the revised pipeline does not depend on them after the new
# filtered file has been created.

export_raw_observations_by_city <- FALSE


# ------------------------------------------------------------
# 4. File paths
# ------------------------------------------------------------

# July 2026 eBird TAR and full EBD member.

tar_path <-
  "Data/Raw/eBird/ebd_US_smp_relJul-2026.tar"

ebd_member <-
  "ebd_US_relJul-2026.txt.gz"


# QC-approved Script 04 checklist metadata.

checklists_path <-
  "Data/Processed/eBird/dryland_city_checklists_2015_2025.csv"


# Main eBird output directory.

ebird_out_dir <-
  "Data/Processed/eBird"


dir.create(
  ebird_out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Existing / compatibility folders.

city_raw_out_dir <- file.path(
  ebird_out_dir,
  "Observations_By_City_Raw"
)

checklist_richness_by_city_dir <- file.path(
  ebird_out_dir,
  "Checklist_Richness_By_City"
)

city_year_unique_by_city_dir <- file.path(
  ebird_out_dir,
  "City_Year_Unique_Species_By_City"
)


dir.create(
  city_raw_out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  checklist_richness_by_city_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  city_year_unique_by_city_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Release marker.

release_marker_path <- file.path(
  ebird_out_dir,
  "script05_ebird_release.txt"
)


# NEW primary filtered EBD product + completion marker.

filtered_ebd_out <- file.path(
  ebird_out_dir,
  "dryland_city_filtered_ebd_observations_2015_2025.tsv"
)

filtered_ebd_complete_marker <- file.path(
  ebird_out_dir,
  "script05_filtered_ebd_complete.txt"
)


# Compact EBD-derived summary intermediates.

checklist_species_counts_ebd_out <- file.path(
  ebird_out_dir,
  "dryland_city_checklist_species_counts_from_ebd_2015_2025.tsv"
)

city_year_species_counts_ebd_out <- file.path(
  ebird_out_dir,
  "dryland_city_year_species_counts_from_ebd_2015_2025.tsv"
)

protocol_lookup_ebd_out <- file.path(
  ebird_out_dir,
  "dryland_city_checklist_protocol_lookup_from_ebd_2015_2025.tsv"
)

ebd_summaries_complete_marker <- file.path(
  ebird_out_dir,
  "script05_ebd_summaries_complete.txt"
)


# Script 05 intermediate / output files.

city_manifest_out <- file.path(
  ebird_out_dir,
  "dryland_city_observation_file_manifest_2015_2025.csv"
)

checklist_city_map_out <- file.path(
  ebird_out_dir,
  "dryland_city_checklist_to_city_file_map_2015_2025.tsv"
)

checklist_metadata_out <- file.path(
  ebird_out_dir,
  "dryland_city_checklist_metadata_for_joins_2015_2025.csv"
)

checklist_richness_out <- file.path(
  ebird_out_dir,
  "dryland_city_checklist_species_richness_2015_2025.csv"
)

city_year_unique_species_out <- file.path(
  ebird_out_dir,
  "dryland_city_year_unique_species_from_observations_2015_2025.csv"
)

city_year_richness_out <- file.path(
  ebird_out_dir,
  "dryland_city_year_species_richness_2015_2025.csv"
)

processing_status_out <- file.path(
  ebird_out_dir,
  "dryland_city_processing_status_2015_2025.csv"
)


# Script 04 already produces the first of these.
# Script 05 overwrites it only after successful richness
# processing.

hotspot_list_out <- file.path(
  ebird_out_dir,
  "dryland_city_hotspot_list_for_GEE_2015_2025.csv"
)

hotspot_latlong_out <- file.path(
  ebird_out_dir,
  "dryland_city_hotspot_latlong_for_GEE_2015_2025.csv"
)


# AWK / shell helper files.

filter_awk_path <- file.path(
  ebird_out_dir,
  "filter_ebd_to_dryland_checklists.awk"
)

filter_log_path <- file.path(
  ebird_out_dir,
  "filter_ebd_to_dryland_checklists_log.txt"
)

summary_awk_path <- file.path(
  ebird_out_dir,
  "summarize_filtered_dryland_ebd.awk"
)

summary_log_path <- file.path(
  ebird_out_dir,
  "summarize_filtered_dryland_ebd_log.txt"
)

legacy_recovery_map_path <- file.path(
  ebird_out_dir,
  "legacy_city_file_recovery_map.tsv"
)


# Legacy completion marker from the previous architecture.
# It is not used as the new recovery checkpoint.

legacy_raw_split_complete_marker <- file.path(
  ebird_out_dir,
  "script05_raw_split_complete.txt"
)


# ------------------------------------------------------------
# 5. Utility functions
# ------------------------------------------------------------

atomic_fwrite <- function(x, path, sep = ",") {
  
  tmp <- paste0(path, ".tmp")
  
  if (file.exists(tmp)) {
    unlink(tmp, force = TRUE)
  }
  
  fwrite(
    x,
    tmp,
    sep = sep
  )
  
  if (!file.exists(tmp) || file.info(tmp)$size == 0) {
    stop("Temporary output was not written:\n", tmp)
  }
  
  if (file.exists(path)) {
    unlink(path, force = TRUE)
  }
  
  if (!file.rename(tmp, path)) {
    stop("Could not finalize output:\n", path)
  }
  
  invisible(path)
}


safe_read_marker_release <- function(marker_path) {
  
  if (!file.exists(marker_path)) {
    return(NA_character_)
  }
  
  x <- readLines(
    marker_path,
    warn = FALSE
  )
  
  release_line <- x[
    grepl("^release=", x)
  ]
  
  if (length(release_line) == 0) {
    return(NA_character_)
  }
  
  sub(
    "^release=",
    "",
    release_line[1]
  )
}


# ------------------------------------------------------------
# 6. Release-safety check
# ------------------------------------------------------------

existing_release <- NA_character_

if (file.exists(release_marker_path)) {
  
  release_lines <- readLines(
    release_marker_path,
    warn = FALSE
  )
  
  if (length(release_lines) > 0) {
    existing_release <- trimws(release_lines[1])
  }
}


new_release_detected <-
  is.na(existing_release) ||
  existing_release != ebird_release


if (new_release_detected) {
  
  cat("\n============================================================\n")
  cat("NEW SCRIPT 05 eBIRD RELEASE DETECTED\n")
  cat("============================================================\n\n")
  
  cat("Previous release marker:\n")
  print(existing_release)
  
  cat("\nCurrent release:\n")
  print(ebird_release)
  
  cat(
    "\nRemoving old Script 05 derived outputs so eBird releases ",
    "cannot be mixed.\n"
  )
  
  
  # Remove old per-city raw observations because they are tied
  # to the previous eBird release.
  
  if (dir.exists(city_raw_out_dir)) {
    unlink(
      city_raw_out_dir,
      recursive = TRUE,
      force = TRUE
    )
  }
  
  
  # Remove per-city summary directories.
  
  if (dir.exists(checklist_richness_by_city_dir)) {
    unlink(
      checklist_richness_by_city_dir,
      recursive = TRUE,
      force = TRUE
    )
  }
  
  if (dir.exists(city_year_unique_by_city_dir)) {
    unlink(
      city_year_unique_by_city_dir,
      recursive = TRUE,
      force = TRUE
    )
  }
  
  
  script05_files_to_clear <- c(
    filtered_ebd_out,
    filtered_ebd_complete_marker,
    checklist_species_counts_ebd_out,
    city_year_species_counts_ebd_out,
    protocol_lookup_ebd_out,
    ebd_summaries_complete_marker,
    city_manifest_out,
    checklist_city_map_out,
    checklist_metadata_out,
    checklist_richness_out,
    city_year_unique_species_out,
    city_year_richness_out,
    processing_status_out,
    hotspot_latlong_out,
    filter_awk_path,
    filter_log_path,
    summary_awk_path,
    summary_log_path,
    legacy_recovery_map_path,
    legacy_raw_split_complete_marker
  )
  
  unlink(
    script05_files_to_clear[
      file.exists(script05_files_to_clear)
    ],
    force = TRUE
  )
  
  writeLines(
    ebird_release,
    release_marker_path
  )
  
  cat("\nCurrent eBird release marker written.\n")
  
} else {
  
  cat("\n============================================================\n")
  cat("SCRIPT 05 RELEASE CHECK\n")
  cat("============================================================\n\n")
  
  cat("Existing outputs belong to:\n")
  print(existing_release)
  
  cat("\nRelease matches current run.\n")
}


# Re-create directories after release cleanup.

dir.create(
  city_raw_out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  checklist_richness_by_city_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  city_year_unique_by_city_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Remove abandoned temporary files from interrupted writes.

for (d in c(
  ebird_out_dir,
  checklist_richness_by_city_dir,
  city_year_unique_by_city_dir
)) {
  
  unlink(
    list.files(
      d,
      pattern = "\\.tmp$",
      full.names = TRUE
    ),
    force = TRUE
  )
}


# ------------------------------------------------------------
# 7. Confirm source files
# ------------------------------------------------------------

if (!file.exists(tar_path)) {
  stop(
    "July 2026 eBird TAR file not found:\n",
    tar_path
  )
}

if (!file.exists(checklists_path)) {
  stop(
    "Script 04 checklist metadata file not found:\n",
    checklists_path
  )
}


cat("\n============================================================\n")
cat("SOURCE FILES\n")
cat("============================================================\n\n")

cat("eBird TAR:\n")
print(tar_path)

cat("\nTAR size (GB):\n")
print(
  round(
    file.info(tar_path)$size / 1024^3,
    3
  )
)

cat("\nScript 04 checklist file:\n")
print(checklists_path)


# ------------------------------------------------------------
# 8. Read and QC Script 04 checklist metadata
# ------------------------------------------------------------

dryland_checklists <- fread(
  checklists_path
)


cat("\n============================================================\n")
cat("SCRIPT 04 INPUT QC\n")
cat("============================================================\n\n")

cat("Checklist rows loaded:\n")
print(nrow(dryland_checklists))


required_checklist_cols <- c(
  "checklist_id",
  "locality_id",
  "locality",
  "observation_date",
  "city_id",
  "NAME10",
  "city_ascii",
  "state_name",
  "county_name",
  "dryland_percent",
  "latitude",
  "longitude",
  "duration_minutes",
  "effort_distance_km",
  "number_observers"
)


missing_checklist_cols <- setdiff(
  required_checklist_cols,
  names(dryland_checklists)
)


if (length(missing_checklist_cols) > 0) {
  stop(
    "Script 04 checklist file is missing required columns:\n",
    paste(
      missing_checklist_cols,
      collapse = "\n"
    )
  )
}


# Standardize key field types.

dryland_checklists[
  ,
  checklist_id := as.character(checklist_id)
]


dryland_checklists[
  ,
  locality_id := as.character(locality_id)
]


dryland_checklists[
  ,
  observation_date := as.character(observation_date)
]


input_checklist_rows <- nrow(dryland_checklists)
input_unique_checklists <- uniqueN(dryland_checklists$checklist_id)
input_city_count <- uniqueN(dryland_checklists$city_id)

input_locality_count <- uniqueN(
  dryland_checklists$locality_id[
    !is.na(dryland_checklists$locality_id) &
      dryland_checklists$locality_id != ""
  ]
)


cat("\nUnique checklist IDs:\n")
print(input_unique_checklists)

cat("\nDryland cities:\n")
print(input_city_count)

cat("\nUnique locality IDs:\n")
print(input_locality_count)

cat("\nObservation-date range:\n")
print(
  range(
    as.Date(dryland_checklists$observation_date),
    na.rm = TRUE
  )
)

cat("\nDryland-percent range:\n")
print(
  range(
    dryland_checklists$dryland_percent,
    na.rm = TRUE
  )
)

cat("\nMissing locality IDs:\n")
print(
  sum(
    is.na(dryland_checklists$locality_id) |
      dryland_checklists$locality_id == ""
  )
)


# Hard input checks.

if (input_checklist_rows != input_unique_checklists) {
  stop("Script 04 input contains duplicate checklist IDs.")
}

if (any(is.na(dryland_checklists$city_id))) {
  stop("Missing city_id values detected in Script 04 input.")
}

if (
  any(
    is.na(dryland_checklists$locality_id) |
    dryland_checklists$locality_id == ""
  )
) {
  stop("Missing locality_id values detected in Script 04 input.")
}

if (
  min(
    as.Date(dryland_checklists$observation_date),
    na.rm = TRUE
  ) < as.Date("2015-01-01") ||
  max(
    as.Date(dryland_checklists$observation_date),
    na.rm = TRUE
  ) > as.Date("2025-12-31")
) {
  stop("Checklist dates fall outside 2015-2025.")
}

if (
  any(
    dryland_checklists$dryland_percent < 50 - 1e-8 |
    dryland_checklists$dryland_percent > 100 + 1e-8,
    na.rm = TRUE
  )
) {
  stop("Invalid dryland_percent values detected in Script 04 input.")
}


# ------------------------------------------------------------
# 9. Create safe city file names
# ------------------------------------------------------------

dryland_checklists[
  ,
  city_file_stub := paste0(
    "city_",
    stringr::str_pad(
      city_id,
      width = 4,
      pad = "0"
    ),
    "_",
    stringr::str_to_lower(
      stringr::str_replace_all(
        state_name,
        "[^A-Za-z0-9]+",
        "_"
      )
    ),
    "_",
    stringr::str_to_lower(
      stringr::str_replace_all(
        city_ascii,
        "[^A-Za-z0-9]+",
        "_"
      )
    )
  )
]


dryland_checklists[
  ,
  city_file_stub := stringr::str_replace_all(
    city_file_stub,
    "_+$",
    ""
  )
]


# ------------------------------------------------------------
# 10. Create city manifest
# ------------------------------------------------------------

city_manifest <- unique(
  dryland_checklists[
    ,
    .(
      city_id,
      NAME10,
      city_ascii,
      state_name,
      county_name,
      dryland_percent,
      city_file_stub
    )
  ]
)

setorder(
  city_manifest,
  state_name,
  city_ascii
)


city_manifest[
  ,
  raw_observation_file := file.path(
    city_raw_out_dir,
    paste0(
      city_file_stub,
      ".tsv"
    )
  )
]


city_manifest[
  ,
  checklist_richness_file := file.path(
    checklist_richness_by_city_dir,
    paste0(
      city_file_stub,
      "_checklist_richness.csv"
    )
  )
]


city_manifest[
  ,
  city_year_unique_file := file.path(
    city_year_unique_by_city_dir,
    paste0(
      city_file_stub,
      "_city_year_unique_species.csv"
    )
  )
]


if (nrow(city_manifest) != input_city_count) {
  stop(
    "City manifest row count does not equal unique Script 04 city count."
  )
}

if (anyDuplicated(city_manifest$city_file_stub) > 0) {
  stop("Duplicate city_file_stub values detected.")
}


atomic_fwrite(
  city_manifest,
  city_manifest_out
)

cat("\nCity manifest rows:\n")
print(nrow(city_manifest))


# ------------------------------------------------------------
# 11. Create checklist-to-city map
# ------------------------------------------------------------

checklist_city_map <- unique(
  dryland_checklists[
    ,
    .(
      checklist_id,
      city_file_stub
    )
  ]
)


if (nrow(checklist_city_map) != input_unique_checklists) {
  stop(
    "Checklist-to-city map does not contain exactly one row per checklist."
  )
}

if (anyDuplicated(checklist_city_map$checklist_id) > 0) {
  stop("Duplicate checklist IDs detected in checklist-to-city map.")
}


atomic_fwrite(
  checklist_city_map,
  checklist_city_map_out,
  sep = "\t"
)

cat("\nChecklist IDs in EBD filter map:\n")
print(nrow(checklist_city_map))


# ------------------------------------------------------------
# 12. Create checklist metadata crosswalk
# ------------------------------------------------------------

checklist_metadata_small <- unique(
  dryland_checklists[
    ,
    .(
      checklist_id,
      locality_id,
      locality,
      observation_date,
      city_id,
      NAME10,
      city_ascii,
      state_name,
      county_name,
      dryland_percent,
      latitude_checklist = latitude,
      longitude_checklist = longitude,
      duration_minutes,
      effort_distance_km,
      number_observers,
      city_file_stub
    )
  ]
)


checklist_metadata_small[
  ,
  checklist_id := as.character(checklist_id)
]


if (nrow(checklist_metadata_small) != input_unique_checklists) {
  stop(
    "Checklist metadata crosswalk does not contain exactly one row per checklist."
  )
}

if (anyDuplicated(checklist_metadata_small$checklist_id) > 0) {
  stop("Duplicate checklist IDs detected in checklist metadata crosswalk.")
}


atomic_fwrite(
  checklist_metadata_small,
  checklist_metadata_out
)

cat("\nChecklist metadata crosswalk written to:\n")
cat(checklist_metadata_out, "\n")


# ------------------------------------------------------------
# 13. Confirm July EBD exists inside TAR
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("JULY 2026 EBD ARCHIVE CHECK\n")
cat("============================================================\n\n")


tar_list_cmd <- paste(
  "/usr/bin/tar -tf",
  shQuote(tar_path)
)


tar_files <- system(
  tar_list_cmd,
  intern = TRUE
)


cat("Requested EBD member:\n")
print(ebd_member)

cat("\nEBD member present in TAR:\n")
print(ebd_member %in% tar_files)


if (!ebd_member %in% tar_files) {
  stop(
    "The July 2026 EBD member was not found inside the TAR:\n",
    ebd_member
  )
}


# ------------------------------------------------------------
# 14. Determine whether the new combined filtered EBD is ready
# ------------------------------------------------------------

filtered_marker_release <- safe_read_marker_release(
  filtered_ebd_complete_marker
)


filtered_ebd_ready <-
  file.exists(filtered_ebd_complete_marker) &&
  identical(filtered_marker_release, ebird_release) &&
  file.exists(filtered_ebd_out) &&
  file.info(filtered_ebd_out)$size > 0


# Track whether upstream products changed during this run so
# dependent summary products cannot be reused accidentally.

filtered_rebuilt_this_run <- FALSE
summaries_rebuilt_this_run <- FALSE


if (force_rebuild_filtered_ebd) {
  filtered_ebd_ready <- FALSE
}


cat("\n============================================================\n")
cat("FILTERED EBD STATUS\n")
cat("============================================================\n\n")

cat("Filtered EBD file exists:\n")
print(file.exists(filtered_ebd_out))

cat("\nFiltered EBD completion marker exists:\n")
print(file.exists(filtered_ebd_complete_marker))

cat("\nFiltered EBD marker release:\n")
print(filtered_marker_release)

cat("\nCompleted filtered EBD can be reused:\n")
print(filtered_ebd_ready)


# ------------------------------------------------------------
# 15. One-time recovery from the completed legacy 647-file split
# ------------------------------------------------------------
#
# This section exists specifically so the successful ~3.5-hour
# split from the previous Script 05 run does not need to be
# repeated.
# ------------------------------------------------------------

if (
  !filtered_ebd_ready &&
  recover_confirmed_legacy_split &&
  !force_rebuild_filtered_ebd
) {
  
  cat("\n============================================================\n")
  cat("CHECKING CONFIRMED LEGACY SPLIT FOR RECOVERY\n")
  cat("============================================================\n\n")
  
  
  expected_legacy_files <- city_manifest$raw_observation_file
  
  legacy_all_exist <- all(
    file.exists(expected_legacy_files)
  )
  
  cat("Expected legacy raw city files:\n")
  print(length(expected_legacy_files))
  
  cat("\nAll legacy raw city files exist:\n")
  print(legacy_all_exist)
  
  
  legacy_base_cols <- c(
    "checklist_id",
    "observation_date",
    "locality_id",
    "locality",
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
  )
  
  
  legacy_headers_ok <- FALSE
  
  
  if (legacy_all_exist) {
    
    legacy_header_bad <- character(0)
    
    for (f in expected_legacy_files) {
      
      h <- names(
        fread(
          f,
          nrows = 0,
          showProgress = FALSE
        )
      )
      
      if (!identical(h, legacy_base_cols)) {
        legacy_header_bad <- c(
          legacy_header_bad,
          f
        )
      }
    }
    
    legacy_headers_ok <- length(legacy_header_bad) == 0
    
    cat("\nLegacy raw-file headers valid:\n")
    print(legacy_headers_ok)
    
    if (!legacy_headers_ok) {
      cat("\nFiles with unexpected legacy headers:\n")
      print(head(legacy_header_bad, 25))
    }
  }
  
  
  if (legacy_all_exist && legacy_headers_ok) {
    
    cat(
      "\nRecovering the completed legacy split into ONE combined ",
      "filtered EBD file.\n"
    )
    cat(
      "This scans only the already-filtered city TSVs; it does NOT ",
      "rescan the national eBird archive.\n\n"
    )
    
    
    legacy_recovery_map <- city_manifest[
      ,
      .(
        city_file_stub,
        raw_observation_file
      )
    ]
    
    atomic_fwrite(
      legacy_recovery_map,
      legacy_recovery_map_path,
      sep = "\t"
    )
    
    
    filtered_tmp <- paste0(
      filtered_ebd_out,
      ".tmp"
    )
    
    unlink(
      c(
        filtered_tmp,
        filtered_ebd_complete_marker
      ),
      force = TRUE
    )
    
    
    legacy_recovery_cmd <- paste(
      "set -o pipefail;",
      "rm -f", shQuote(filtered_tmp), ";",
      "printf '%s\\n'",
      shQuote(
        paste(
          c(
            "city_file_stub",
            legacy_base_cols
          ),
          collapse = "\t"
        )
      ),
      ">", shQuote(filtered_tmp), ";",
      "while IFS=$(printf '\\t') read -r stub filepath; do",
      "if [ \"$stub\" = \"city_file_stub\" ]; then continue; fi;",
      "tail -n +2 \"$filepath\" |",
      "awk -v stub=\"$stub\" 'BEGIN { OFS=\"\\t\" } { print stub, $0 }'",
      ">>", shQuote(filtered_tmp), ";",
      "done <", shQuote(legacy_recovery_map_path), ";",
      "mv", shQuote(filtered_tmp), shQuote(filtered_ebd_out), ";",
      "printf '%s\\n'",
      shQuote(paste0("release=", ebird_release)),
      shQuote(paste0("completed=", Sys.time())),
      shQuote("source=confirmed_legacy_city_split"),
      ">", shQuote(filtered_ebd_complete_marker)
    )
    
    
    legacy_recovery_start <- Sys.time()
    
    legacy_recovery_exit <- system(
      paste(
        "bash -lc",
        shQuote(legacy_recovery_cmd)
      )
    )
    
    legacy_recovery_end <- Sys.time()
    
    
    cat("Legacy recovery exit code:\n")
    print(legacy_recovery_exit)
    
    cat("\nLegacy recovery elapsed time:\n")
    print(
      legacy_recovery_end - legacy_recovery_start
    )
    
    
    if (legacy_recovery_exit != 0) {
      stop(
        "Legacy city-file recovery failed. The national EBD has NOT been rescanned."
      )
    }
    
    
    filtered_marker_release <- safe_read_marker_release(
      filtered_ebd_complete_marker
    )
    
    filtered_ebd_ready <-
      file.exists(filtered_ebd_complete_marker) &&
      identical(filtered_marker_release, ebird_release) &&
      file.exists(filtered_ebd_out) &&
      file.info(filtered_ebd_out)$size > 0
    
    filtered_rebuilt_this_run <- TRUE
    
    
    if (!filtered_ebd_ready) {
      stop(
        "Legacy recovery command completed but the new filtered EBD checkpoint is invalid."
      )
    }
    
    cat("\nLegacy recovery completed successfully.\n")
    cat("New combined filtered EBD:\n")
    cat(filtered_ebd_out, "\n")
    
  } else {
    
    cat(
      "\nA safe legacy recovery source was not available. ",
      "The revised national filter will be used instead.\n"
    )
  }
}


# ------------------------------------------------------------
# 16. Build the new combined filtered EBD if needed
# ------------------------------------------------------------

if (!filtered_ebd_ready) {
  
  cat("\n============================================================\n")
  cat("FILTERING JULY 2026 NATIONAL EBD\n")
  cat("============================================================\n\n")
  
  cat(
    "The national EBD will be streamed once. Matching rows will ",
    "be written to ONE combined filtered file.\n"
  )
  
  cat(
    "This replaces the previous 647-file AWK write/close pattern.\n\n"
  )
  
  
  filtered_tmp <- paste0(
    filtered_ebd_out,
    ".tmp"
  )
  
  unlink(
    c(
      filtered_tmp,
      filtered_ebd_complete_marker,
      filter_log_path
    ),
    force = TRUE
  )
  
  
  filter_awk <- '
BEGIN {
  FS = "\\t"
  OFS = "\\t"
}

NR == FNR {

  if (FNR == 1 && $1 == "checklist_id") {
    next
  }

  checklist_to_city[$1] = $2
  next
}

FNR == 1 {

  for (i = 1; i <= NF; i++) {
    col[$i] = i
  }

  id_col = col["SAMPLING EVENT IDENTIFIER"]
  date_col = col["OBSERVATION DATE"]
  locality_id_col = col["LOCALITY ID"]
  locality_col = col["LOCALITY"]
  common_col = col["COMMON NAME"]
  sci_col = col["SCIENTIFIC NAME"]
  tax_col = col["TAXONOMIC ORDER"]
  category_col = col["CATEGORY"]
  count_col = col["OBSERVATION COUNT"]
  lat_col = col["LATITUDE"]
  lon_col = col["LONGITUDE"]
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

  if (locality_id_col == 0) {
    print "Missing LOCALITY ID" > "/dev/stderr"
    missing = 1
  }

  if (locality_col == 0) {
    print "Missing LOCALITY" > "/dev/stderr"
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

  if (category_col == 0) {
    print "Missing CATEGORY" > "/dev/stderr"
    missing = 1
  }

  if (count_col == 0) {
    print "Missing OBSERVATION COUNT" > "/dev/stderr"
    missing = 1
  }

  if (lat_col == 0) {
    print "Missing LATITUDE" > "/dev/stderr"
    missing = 1
  }

  if (lon_col == 0) {
    print "Missing LONGITUDE" > "/dev/stderr"
    missing = 1
  }

  if (protocol_col == 0) {
    print "Missing PROTOCOL NAME" > "/dev/stderr"
    missing = 1
  }

  if (missing == 1) {

    print "Columns found in EBD file:" > "/dev/stderr"

    for (i = 1; i <= NF; i++) {
      print i, $i > "/dev/stderr"
    }

    exit 1
  }

  header = "city_file_stub"
  header = header OFS "checklist_id"
  header = header OFS "observation_date"
  header = header OFS "locality_id"
  header = header OFS "locality"
  header = header OFS "common_name"
  header = header OFS "scientific_name"
  header = header OFS "taxonomic_order"
  header = header OFS "category"
  header = header OFS "observation_count"
  header = header OFS "latitude"
  header = header OFS "longitude"
  header = header OFS "protocol_name"
  header = header OFS "approved"
  header = header OFS "reviewed"

  print header > outtmp

  next
}

{

  checklist_id = $id_col

  if (checklist_id in checklist_to_city) {

    city_stub = checklist_to_city[checklist_id]

    locality_id_val = ""
    locality_val = ""
    tax_val = ""
    category_val = ""
    lat_val = ""
    lon_val = ""
    protocol_val = ""
    approved_val = ""
    reviewed_val = ""

    if (locality_id_col > 0) locality_id_val = $locality_id_col
    if (locality_col > 0) locality_val = $locality_col
    if (tax_col > 0) tax_val = $tax_col
    if (category_col > 0) category_val = $category_col
    if (lat_col > 0) lat_val = $lat_col
    if (lon_col > 0) lon_val = $lon_col
    if (protocol_col > 0) protocol_val = $protocol_col
    if (approved_col > 0) approved_val = $approved_col
    if (reviewed_col > 0) reviewed_val = $reviewed_col

    row = city_stub
    row = row OFS checklist_id
    row = row OFS $date_col
    row = row OFS locality_id_val
    row = row OFS locality_val
    row = row OFS $common_col
    row = row OFS $sci_col
    row = row OFS tax_val
    row = row OFS category_val
    row = row OFS $count_col
    row = row OFS lat_val
    row = row OFS lon_val
    row = row OFS protocol_val
    row = row OFS approved_val
    row = row OFS reviewed_val

    print row > outtmp
  }
}
'


writeLines(
  filter_awk,
  filter_awk_path,
  useBytes = TRUE
)


filter_pipeline_cmd <- paste(
  "set -o pipefail;",
  "rm -f", shQuote(filtered_tmp), shQuote(filter_log_path), ";",
  "{",
  "LC_ALL=C env -u MallocStackLogging -u MallocStackLoggingNoCompact",
  "nice -n 10 /usr/bin/tar -O -x -f",
  shQuote(tar_path),
  shQuote(ebd_member),
  "|",
  "LC_ALL=C nice -n 10 gzip -cd",
  "|",
  "LC_ALL=C env -u MallocStackLogging -u MallocStackLoggingNoCompact",
  "nice -n 10 awk",
  "-v", shQuote(paste0("outtmp=", filtered_tmp)),
  "-f", shQuote(filter_awk_path),
  shQuote(checklist_city_map_out),
  "-",
  ";",
  "}",
  "2>", shQuote(filter_log_path),
  "&&",
  "mv", shQuote(filtered_tmp), shQuote(filtered_ebd_out),
  "&&",
  "printf '%s\\n'",
  shQuote(paste0("release=", ebird_release)),
  shQuote(paste0("completed=", Sys.time())),
  shQuote("source=national_ebd_stream"),
  ">", shQuote(filtered_ebd_complete_marker)
)


cat("Starting national EBD stream.\n")
cat("Filter log:\n")
cat(filter_log_path, "\n\n")

filter_start <- Sys.time()

filter_exit <- system(
  paste(
    "bash -lc",
    shQuote(filter_pipeline_cmd)
  )
)

filter_end <- Sys.time()


cat("\nFiltered EBD exit code:\n")
print(filter_exit)

cat("\nFiltered EBD elapsed time:\n")
print(filter_end - filter_start)

cat("\nFiltered EBD log:\n")

if (file.exists(filter_log_path)) {
  
  filter_log <- readLines(
    filter_log_path,
    warn = FALSE
  )
  
  if (length(filter_log) > 0) {
    print(
      tail(
        filter_log,
        50
      )
    )
  } else {
    cat("No log messages.\n")
  }
}


if (filter_exit != 0) {
  stop(
    "National July 2026 EBD filter failed. The completion marker was not written."
  )
}


filtered_marker_release <- safe_read_marker_release(
  filtered_ebd_complete_marker
)

filtered_ebd_ready <-
  file.exists(filtered_ebd_complete_marker) &&
  identical(filtered_marker_release, ebird_release) &&
  file.exists(filtered_ebd_out) &&
  file.info(filtered_ebd_out)$size > 0

filtered_rebuilt_this_run <- TRUE


if (!filtered_ebd_ready) {
  stop(
    "National filter returned success but the filtered EBD checkpoint is invalid."
  )
}

cat("\nFiltered EBD checkpoint created successfully.\n")

} else {
  
  cat("\n============================================================\n")
  cat("REUSING COMPLETED FILTERED EBD\n")
  cat("============================================================\n\n")
  
  cat("Filtered EBD:\n")
  cat(filtered_ebd_out, "\n")
}


# ------------------------------------------------------------
# 17. Filtered EBD QC
# ------------------------------------------------------------

filtered_header <- names(
  fread(
    filtered_ebd_out,
    nrows = 0,
    sep = "\t",
    showProgress = FALSE
  )
)


required_filtered_cols <- c(
  "city_file_stub",
  "checklist_id",
  "observation_date",
  "locality_id",
  "locality",
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
)


missing_filtered_cols <- setdiff(
  required_filtered_cols,
  filtered_header
)


if (length(missing_filtered_cols) > 0) {
  stop(
    "Combined filtered EBD is missing required columns:\n",
    paste(
      missing_filtered_cols,
      collapse = "\n"
    )
  )
}


cat("\n============================================================\n")
cat("FILTERED EBD QC\n")
cat("============================================================\n\n")

cat("Filtered EBD size (GB):\n")
print(
  round(
    file.info(filtered_ebd_out)$size / 1024^3,
    3
  )
)

cat("\nFiltered EBD header columns:\n")
print(filtered_header)


# ------------------------------------------------------------
# 18. Build compact EBD summaries using ONE pass + disk sort
# ------------------------------------------------------------
#
# The filtered observation file can still be large, so we do
# NOT load it into R. Instead:
#
#   1. AWK scans it once and writes three compact pair files.
#   2. macOS sort groups those pairs on disk.
#   3. AWK converts the sorted pairs into small count tables.
#
# This preserves unique-species logic without requiring the
# full observation table in R memory.
# ------------------------------------------------------------

if (filtered_rebuilt_this_run) {
  
  unlink(
    c(
      checklist_species_counts_ebd_out,
      city_year_species_counts_ebd_out,
      protocol_lookup_ebd_out,
      ebd_summaries_complete_marker
    ),
    force = TRUE
  )
}


summary_marker_release <- safe_read_marker_release(
  ebd_summaries_complete_marker
)


summaries_ready <-
  file.exists(ebd_summaries_complete_marker) &&
  identical(summary_marker_release, ebird_release) &&
  file.exists(checklist_species_counts_ebd_out) &&
  file.exists(city_year_species_counts_ebd_out) &&
  file.exists(protocol_lookup_ebd_out) &&
  file.info(checklist_species_counts_ebd_out)$size > 0 &&
  file.info(city_year_species_counts_ebd_out)$size > 0 &&
  file.info(protocol_lookup_ebd_out)$size > 0


if (force_rebuild_ebd_summaries) {
  summaries_ready <- FALSE
}


if (!summaries_ready) {
  
  cat("\n============================================================\n")
  cat("BUILDING COMPACT EBD SUMMARIES\n")
  cat("============================================================\n\n")
  
  
  checklist_pairs_tmp <- file.path(
    ebird_out_dir,
    "script05_checklist_species_pairs.tmp"
  )
  
  city_year_pairs_tmp <- file.path(
    ebird_out_dir,
    "script05_city_year_species_pairs.tmp"
  )
  
  protocol_pairs_tmp <- file.path(
    ebird_out_dir,
    "script05_protocol_pairs.tmp"
  )
  
  checklist_counts_tmp <- paste0(
    checklist_species_counts_ebd_out,
    ".tmp"
  )
  
  city_year_counts_tmp <- paste0(
    city_year_species_counts_ebd_out,
    ".tmp"
  )
  
  protocol_lookup_tmp <- paste0(
    protocol_lookup_ebd_out,
    ".tmp"
  )
  
  
  unlink(
    c(
      checklist_pairs_tmp,
      city_year_pairs_tmp,
      protocol_pairs_tmp,
      checklist_counts_tmp,
      city_year_counts_tmp,
      protocol_lookup_tmp,
      ebd_summaries_complete_marker,
      summary_log_path
    ),
    force = TRUE
  )
  
  
  summary_awk <- '
BEGIN {
  FS = "\\t"
  OFS = "\\t"
}

FNR == 1 {

  for (i = 1; i <= NF; i++) {
    col[$i] = i
  }

  city_col = col["city_file_stub"]
  id_col = col["checklist_id"]
  date_col = col["observation_date"]
  sci_col = col["scientific_name"]
  category_col = col["category"]
  protocol_col = col["protocol_name"]

  # Keep the full condition on one line for macOS / BSD awk.
  # BSD awk can reject a newline immediately after `if (`.
  if (city_col == 0 || id_col == 0 || date_col == 0 || sci_col == 0 || category_col == 0 || protocol_col == 0) {
    print "Filtered EBD is missing required summary columns." > "/dev/stderr"
    exit 1
  }

  next
}

{

  checklist_id = $id_col
  protocol_name = $protocol_col

  if (protocol_name != "") {
    print checklist_id, protocol_name >> protocol_pairs
  }

  if ($category_col == "species" && $sci_col != "") {

    scientific_name = $sci_col
    city_stub = $city_col
    year = substr($date_col, 1, 4)

    print checklist_id, scientific_name >> checklist_pairs
    print city_stub, year, scientific_name >> city_year_pairs
  }
}
'


writeLines(
  summary_awk,
  summary_awk_path,
  useBytes = TRUE
)


checklist_count_awk <- paste0(
  "BEGIN { FS=OFS=\"\\t\"; ",
  "print \"checklist_id\",\"species_richness\",\"total_observation_records\"; ",
  "prev_id=\"\"; prev_sp=\"\" } ",
  "{ ",
  "if ($1 != prev_id) { ",
  "if (prev_id != \"\") print prev_id, richness, total; ",
  "prev_id=$1; prev_sp=\"\"; richness=0; total=0; ",
  "} ",
  "total++; ",
  "if ($2 != prev_sp) { richness++; prev_sp=$2; } ",
  "} ",
  "END { if (prev_id != \"\") print prev_id, richness, total }"
)


city_year_count_awk <- paste0(
  "BEGIN { FS=OFS=\"\\t\"; ",
  "print \"city_file_stub\",\"year\",\"true_city_year_species_observed\",\"total_species_observation_records\"; ",
  "prev_city=\"\"; prev_year=\"\"; prev_sp=\"\" } ",
  "{ ",
  "if ($1 != prev_city || $2 != prev_year) { ",
  "if (prev_city != \"\") print prev_city, prev_year, richness, total; ",
  "prev_city=$1; prev_year=$2; prev_sp=\"\"; richness=0; total=0; ",
  "} ",
  "total++; ",
  "if ($3 != prev_sp) { richness++; prev_sp=$3; } ",
  "} ",
  "END { if (prev_city != \"\") print prev_city, prev_year, richness, total }"
)


# Portable tab delimiter for both GNU and macOS / BSD sort.
# This expands to one literal tab character in the shell.
sort_tab_shell <- "\"$(printf '\\t')\""


summary_pipeline_cmd <- paste(
  "set -o pipefail;",
  "rm -f",
  shQuote(checklist_pairs_tmp),
  shQuote(city_year_pairs_tmp),
  shQuote(protocol_pairs_tmp),
  shQuote(checklist_counts_tmp),
  shQuote(city_year_counts_tmp),
  shQuote(protocol_lookup_tmp),
  shQuote(summary_log_path),
  ";",
  "LC_ALL=C nice -n 10 awk",
  "-v", shQuote(paste0("checklist_pairs=", checklist_pairs_tmp)),
  "-v", shQuote(paste0("city_year_pairs=", city_year_pairs_tmp)),
  "-v", shQuote(paste0("protocol_pairs=", protocol_pairs_tmp)),
  "-f", shQuote(summary_awk_path),
  shQuote(filtered_ebd_out),
  "2>", shQuote(summary_log_path),
  "&&",
  "LC_ALL=C nice -n 10 sort -t", sort_tab_shell,
  "-k1,1 -k2,2",
  shQuote(checklist_pairs_tmp),
  "| awk", shQuote(checklist_count_awk),
  ">", shQuote(checklist_counts_tmp),
  "&&",
  "LC_ALL=C nice -n 10 sort -t", sort_tab_shell,
  "-k1,1 -k2,2 -k3,3",
  shQuote(city_year_pairs_tmp),
  "| awk", shQuote(city_year_count_awk),
  ">", shQuote(city_year_counts_tmp),
  "&&",
  "{ printf 'checklist_id\\tprotocol_name\\n';",
  "LC_ALL=C nice -n 10 sort -t", sort_tab_shell,
  "-k1,1 -k2,2 -u",
  shQuote(protocol_pairs_tmp),
  "; } >", shQuote(protocol_lookup_tmp),
  "&&",
  "mv", shQuote(checklist_counts_tmp), shQuote(checklist_species_counts_ebd_out),
  "&&",
  "mv", shQuote(city_year_counts_tmp), shQuote(city_year_species_counts_ebd_out),
  "&&",
  "mv", shQuote(protocol_lookup_tmp), shQuote(protocol_lookup_ebd_out),
  "&&",
  "rm -f",
  shQuote(checklist_pairs_tmp),
  shQuote(city_year_pairs_tmp),
  shQuote(protocol_pairs_tmp),
  "&&",
  "printf '%s\\n'",
  shQuote(paste0("release=", ebird_release)),
  shQuote(paste0("completed=", Sys.time())),
  ">", shQuote(ebd_summaries_complete_marker)
)


summary_start <- Sys.time()

summary_exit <- system(
  paste(
    "bash -lc",
    shQuote(summary_pipeline_cmd)
  )
)

summary_end <- Sys.time()


cat("Compact-summary exit code:\n")
print(summary_exit)

cat("\nCompact-summary elapsed time:\n")
print(summary_end - summary_start)


if (file.exists(summary_log_path)) {
  
  summary_log <- readLines(
    summary_log_path,
    warn = FALSE
  )
  
  cat("\nCompact-summary log:\n")
  
  if (length(summary_log) > 0) {
    print(
      tail(
        summary_log,
        50
      )
    )
  } else {
    cat("No log messages.\n")
  }
}


if (summary_exit != 0) {
  stop(
    "Compact EBD summary stage failed. The filtered national result remains reusable."
  )
}


summary_marker_release <- safe_read_marker_release(
  ebd_summaries_complete_marker
)

summaries_ready <-
  file.exists(ebd_summaries_complete_marker) &&
  identical(summary_marker_release, ebird_release) &&
  file.exists(checklist_species_counts_ebd_out) &&
  file.exists(city_year_species_counts_ebd_out) &&
  file.exists(protocol_lookup_ebd_out)


if (!summaries_ready) {
  stop(
    "Compact summary command returned success but summary checkpoint is invalid."
  )
}

summaries_rebuilt_this_run <- TRUE

cat("\nCompact EBD summary checkpoint created successfully.\n")

} else {
  
  cat("\n============================================================\n")
  cat("REUSING COMPLETED COMPACT EBD SUMMARIES\n")
  cat("============================================================\n\n")
}


# ------------------------------------------------------------
# 19. Read and QC compact EBD summaries
# ------------------------------------------------------------

checklist_counts <- fread(
  checklist_species_counts_ebd_out,
  sep = "\t",
  showProgress = FALSE
)


city_year_species_counts <- fread(
  city_year_species_counts_ebd_out,
  sep = "\t",
  showProgress = FALSE
)


protocol_lookup_pairs <- fread(
  protocol_lookup_ebd_out,
  sep = "\t",
  showProgress = FALSE
)


checklist_counts[
  ,
  checklist_id := as.character(checklist_id)
]

protocol_lookup_pairs[
  ,
  checklist_id := as.character(checklist_id)
]

city_year_species_counts[
  ,
  year := as.integer(year)
]


cat("\n============================================================\n")
cat("COMPACT EBD SUMMARY QC\n")
cat("============================================================\n\n")

cat("Checklists with >=1 species-level EBD record:\n")
print(nrow(checklist_counts))

cat("\nCity-year count rows from EBD:\n")
print(nrow(city_year_species_counts))

cat("\nUnique checklist/protocol pairs:\n")
print(nrow(protocol_lookup_pairs))


if (anyDuplicated(checklist_counts$checklist_id) > 0) {
  stop("Duplicate checklist IDs detected in checklist species counts.")
}


unexpected_count_ids <- setdiff(
  checklist_counts$checklist_id,
  checklist_metadata_small$checklist_id
)

unexpected_protocol_ids <- setdiff(
  protocol_lookup_pairs$checklist_id,
  checklist_metadata_small$checklist_id
)

if (length(unexpected_count_ids) > 0) {
  stop(
    "Species-count summaries contain checklist IDs not present in Script 04."
  )
}

if (length(unexpected_protocol_ids) > 0) {
  stop(
    "Protocol summaries contain checklist IDs not present in Script 04."
  )
}


# Because protocol_lookup_ebd_out was sorted with -u, duplicate
# checklist IDs now mean the same checklist had >1 distinct
# non-missing protocol value in the EBD.

protocol_multi <- protocol_lookup_pairs[
  ,
  .N,
  by = checklist_id
][
  N > 1
]


if (nrow(protocol_multi) > 0) {
  stop(
    "More than one distinct protocol_name was detected within ",
    "one or more checklists. Example checklist IDs:\n",
    paste(
      head(protocol_multi$checklist_id, 25),
      collapse = "\n"
    )
  )
}


protocol_lookup <- protocol_lookup_pairs[
  ,
  .(
    checklist_id,
    protocol_name
  )
]


# ------------------------------------------------------------
# 20. Build combined checklist richness table
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("BUILDING CHECKLIST RICHNESS\n")
cat("============================================================\n\n")


checklist_richness <- merge(
  copy(checklist_metadata_small),
  checklist_counts,
  by = "checklist_id",
  all.x = TRUE,
  sort = FALSE
)


checklist_richness <- merge(
  checklist_richness,
  protocol_lookup,
  by = "checklist_id",
  all.x = TRUE,
  sort = FALSE
)


checklist_richness[
  is.na(species_richness),
  species_richness := 0
]

checklist_richness[
  is.na(total_observation_records),
  total_observation_records := 0
]


checklist_richness[
  ,
  year := as.integer(
    substr(
      observation_date,
      1,
      4
    )
  )
]


if (nrow(checklist_richness) != input_checklist_rows) {
  stop(
    "Checklist richness row count does not equal Script 04 input."
  )
}

if (uniqueN(checklist_richness$checklist_id) != input_unique_checklists) {
  stop(
    "Checklist richness does not contain every unique Script 04 checklist."
  )
}

if (anyDuplicated(checklist_richness$checklist_id) > 0) {
  stop("Duplicate checklist IDs were created during richness processing.")
}

if (any(is.na(checklist_richness$species_richness))) {
  stop("Missing species_richness values detected.")
}

if (any(checklist_richness$species_richness < 0)) {
  stop("Negative species richness detected.")
}


setcolorder(
  checklist_richness,
  c(
    "checklist_id",
    "locality_id",
    "locality",
    "observation_date",
    "year",
    "city_id",
    "NAME10",
    "city_ascii",
    "state_name",
    "county_name",
    "dryland_percent",
    "duration_minutes",
    "effort_distance_km",
    "number_observers",
    "protocol_name",
    "latitude_checklist",
    "longitude_checklist",
    "species_richness",
    "total_observation_records",
    "city_file_stub"
  )
)


missing_checklist_ids <- setdiff(
  dryland_checklists$checklist_id,
  checklist_richness$checklist_id
)

extra_checklist_ids <- setdiff(
  checklist_richness$checklist_id,
  dryland_checklists$checklist_id
)


if (
  length(missing_checklist_ids) > 0 ||
  length(extra_checklist_ids) > 0
) {
  stop(
    "Checklist ID sets do not match between Script 04 and Script 05."
  )
}


atomic_fwrite(
  checklist_richness,
  checklist_richness_out
)

cat("Checklist richness rows:\n")
print(nrow(checklist_richness))


# ------------------------------------------------------------
# 21. Build combined city-year unique species table
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("BUILDING CITY-YEAR UNIQUE SPECIES\n")
cat("============================================================\n\n")


city_year_base <- unique(
  checklist_metadata_small[
    ,
    .(
      city_file_stub,
      year = as.integer(
        substr(
          observation_date,
          1,
          4
        )
      )
    )
  ]
)


city_year_base_keys <- paste(
  city_year_base$city_file_stub,
  city_year_base$year,
  sep = "__"
)

city_year_count_keys <- paste(
  city_year_species_counts$city_file_stub,
  city_year_species_counts$year,
  sep = "__"
)

unexpected_city_year_keys <- setdiff(
  city_year_count_keys,
  city_year_base_keys
)

if (length(unexpected_city_year_keys) > 0) {
  stop(
    "EBD city-year summaries contain combinations not present in Script 04."
  )
}


city_year_unique_species <- merge(
  city_year_base,
  city_year_species_counts,
  by = c(
    "city_file_stub",
    "year"
  ),
  all.x = TRUE,
  sort = FALSE
)


city_year_unique_species[
  is.na(true_city_year_species_observed),
  true_city_year_species_observed := 0
]

city_year_unique_species[
  is.na(total_species_observation_records),
  total_species_observation_records := 0
]


city_year_unique_species <- merge(
  city_year_unique_species,
  city_manifest[
    ,
    .(
      city_id,
      NAME10,
      city_ascii,
      state_name,
      county_name,
      dryland_percent,
      city_file_stub
    )
  ],
  by = "city_file_stub",
  all.x = TRUE,
  sort = FALSE
)


setcolorder(
  city_year_unique_species,
  c(
    "city_id",
    "NAME10",
    "city_ascii",
    "state_name",
    "county_name",
    "year",
    "dryland_percent",
    "true_city_year_species_observed",
    "total_species_observation_records",
    "city_file_stub"
  )
)


setorder(
  city_year_unique_species,
  state_name,
  city_ascii,
  year
)


if (any(is.na(city_year_unique_species$true_city_year_species_observed))) {
  stop("Missing city-year unique species values detected.")
}


atomic_fwrite(
  city_year_unique_species,
  city_year_unique_species_out
)

cat("Combined city-year unique species rows:\n")
print(nrow(city_year_unique_species))


# ------------------------------------------------------------
# 22. Write per-city summary files from compact combined tables
# ------------------------------------------------------------
#
# These outputs are retained for compatibility, but the expensive
# city-by-city RAW EBD processing loop has been removed.
# ------------------------------------------------------------

if (rebuild_city_summaries) {
  
  unlink(
    list.files(
      checklist_richness_by_city_dir,
      pattern = "\\.csv$",
      full.names = TRUE
    ),
    force = TRUE
  )
  
  unlink(
    list.files(
      city_year_unique_by_city_dir,
      pattern = "\\.csv$",
      full.names = TRUE
    ),
    force = TRUE
  )
}


cat("\n============================================================\n")
cat("WRITING / CHECKING PER-CITY SUMMARY FILES\n")
cat("============================================================\n\n")


status_rows <- vector(
  "list",
  nrow(city_manifest)
)


for (i in seq_len(nrow(city_manifest))) {
  
  this_stub <- city_manifest$city_file_stub[i]
  checklist_city_out <- city_manifest$checklist_richness_file[i]
  city_year_city_out <- city_manifest$city_year_unique_file[i]
  
  
  checklist_done <-
    file.exists(checklist_city_out) &&
    file.info(checklist_city_out)$size > 0
  
  city_year_done <-
    file.exists(city_year_city_out) &&
    file.info(city_year_city_out)$size > 0
  
  
  if (
    !checklist_done ||
    rebuild_city_summaries ||
    summaries_rebuilt_this_run
  ) {
    
    checklist_i <- checklist_richness[
      city_file_stub == this_stub
    ]
    
    if (nrow(checklist_i) == 0) {
      stop(
        "No checklist richness rows found for city:\n",
        this_stub
      )
    }
    
    atomic_fwrite(
      checklist_i,
      checklist_city_out
    )
    
    checklist_done <- TRUE
  }
  
  
  if (
    !city_year_done ||
    rebuild_city_summaries ||
    summaries_rebuilt_this_run
  ) {
    
    city_year_i <- city_year_unique_species[
      city_file_stub == this_stub
    ]
    
    if (nrow(city_year_i) == 0) {
      stop(
        "No city-year rows found for city:\n",
        this_stub
      )
    }
    
    atomic_fwrite(
      city_year_i,
      city_year_city_out
    )
    
    city_year_done <- TRUE
  }
  
  
  status_rows[[i]] <- data.table(
    city_file_stub = this_stub,
    checklist_richness_done = checklist_done,
    city_year_unique_done = city_year_done,
    status = "complete_from_combined_filtered_ebd"
  )
  
  
  if (i %% 50 == 0 || i == nrow(city_manifest)) {
    cat(
      "Per-city summaries checked/written:",
      i,
      "of",
      nrow(city_manifest),
      "\n"
    )
  }
}


processing_status <- rbindlist(
  status_rows,
  fill = TRUE
)


atomic_fwrite(
  processing_status,
  processing_status_out
)


# ------------------------------------------------------------
# 23. Optional raw per-city observation export
# ------------------------------------------------------------
#
# The revised analysis does not need these files. This block is
# intentionally disabled by default because recreating them adds
# substantial disk I/O with no analytical benefit.
#
# If a later task specifically requires city-level raw EBD files,
# they can be exported from filtered_ebd_out without rescanning
# the national archive.
# ------------------------------------------------------------

if (export_raw_observations_by_city) {
  
  warning(
    "export_raw_observations_by_city = TRUE was requested. ",
    "For efficiency, the revised core Script 05 does not perform ",
    "this optional export automatically. Use filtered_ebd_out as ",
    "the analysis-ready raw observation source, or create a ",
    "separate export-only helper script if city TSVs are required."
  )
}


# ------------------------------------------------------------
# 24. Build hotspot/location outputs for GEE
# ------------------------------------------------------------

hotspot_list <- checklist_richness[
  
  !is.na(locality_id) &
    locality_id != "",
  
  .(
    n_checklists = as.integer(
      uniqueN(checklist_id)
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
      latitude_checklist,
      na.rm = TRUE
    ),
    mean_longitude = mean(
      longitude_checklist,
      na.rm = TRUE
    )
  ),
  
  by = .(
    city_id,
    NAME10,
    city_ascii,
    state_name,
    county_name,
    dryland_percent,
    locality_id,
    locality
  )
]


setorder(
  hotspot_list,
  state_name,
  city_ascii,
  locality
)


atomic_fwrite(
  hotspot_list,
  hotspot_list_out
)


hotspot_latlong <- hotspot_list[
  ,
  .(
    latitude = mean(
      mean_latitude,
      na.rm = TRUE
    ),
    longitude = mean(
      mean_longitude,
      na.rm = TRUE
    ),
    n_checklists = sum(
      n_checklists,
      na.rm = TRUE
    )
  ),
  by = .(
    locality_id,
    locality
  )
]


setorder(
  hotspot_latlong,
  locality_id
)


atomic_fwrite(
  hotspot_latlong,
  hotspot_latlong_out
)


cat("\nHotspot/locality rows:\n")
print(nrow(hotspot_list))

cat("\nUnique hotspot/locality IDs:\n")
print(uniqueN(hotspot_list$locality_id))

cat("\nSimple hotspot lat/long rows:\n")
print(nrow(hotspot_latlong))


# ------------------------------------------------------------
# 25. Build final city-year richness summary
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("BUILDING FINAL CITY-YEAR RICHNESS SUMMARY\n")
cat("============================================================\n\n")


checklist_richness[
  ,
  species_richness := as.numeric(species_richness)
]

checklist_richness[
  ,
  total_observation_records := as.numeric(total_observation_records)
]

checklist_richness[
  ,
  duration_minutes := as.numeric(duration_minutes)
]

checklist_richness[
  ,
  effort_distance_km := as.numeric(effort_distance_km)
]

checklist_richness[
  ,
  number_observers := as.numeric(number_observers)
]


city_year_effort_summary <- checklist_richness[
  ,
  .(
    n_checklists = as.integer(
      uniqueN(checklist_id)
    ),
    n_localities = as.integer(
      uniqueN(
        locality_id[
          !is.na(locality_id) &
            locality_id != ""
        ]
      )
    ),
    mean_checklist_richness = as.numeric(
      mean(
        species_richness,
        na.rm = TRUE
      )
    ),
    median_checklist_richness = as.numeric(
      median(
        species_richness,
        na.rm = TRUE
      )
    ),
    summed_checklist_richness = as.numeric(
      sum(
        species_richness,
        na.rm = TRUE
      )
    ),
    total_observation_records = as.numeric(
      sum(
        total_observation_records,
        na.rm = TRUE
      )
    ),
    mean_duration_minutes = as.numeric(
      mean(
        duration_minutes,
        na.rm = TRUE
      )
    ),
    mean_effort_distance_km = as.numeric(
      mean(
        effort_distance_km,
        na.rm = TRUE
      )
    ),
    mean_number_observers = as.numeric(
      mean(
        number_observers,
        na.rm = TRUE
      )
    )
  ),
  
  by = .(
    city_id,
    NAME10,
    city_ascii,
    state_name,
    county_name,
    year,
    dryland_percent
  )
]


city_year_richness <- merge(
  city_year_effort_summary,
  city_year_unique_species[
    ,
    .(
      city_id,
      NAME10,
      city_ascii,
      state_name,
      county_name,
      year,
      dryland_percent,
      true_city_year_species_observed,
      total_species_observation_records
    )
  ],
  by = c(
    "city_id",
    "NAME10",
    "city_ascii",
    "state_name",
    "county_name",
    "year",
    "dryland_percent"
  ),
  all.x = TRUE,
  sort = FALSE
)


setorder(
  city_year_richness,
  state_name,
  city_ascii,
  year
)


if (any(is.na(city_year_richness$true_city_year_species_observed))) {
  stop(
    "Missing true city-year species richness after final join."
  )
}


atomic_fwrite(
  city_year_richness,
  city_year_richness_out
)


cat("Final city-year richness rows:\n")
print(nrow(city_year_richness))


# ------------------------------------------------------------
# 26. Final Script 05 QC
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SCRIPT 05 FINAL QC\n")
cat("============================================================\n\n")


final_checklist_rows <- nrow(checklist_richness)
final_unique_checklists <- uniqueN(checklist_richness$checklist_id)
final_city_count <- uniqueN(checklist_richness$city_id)

final_locality_count <- uniqueN(
  checklist_richness$locality_id[
    !is.na(checklist_richness$locality_id) &
      checklist_richness$locality_id != ""
  ]
)

zero_richness_checklists <- sum(
  checklist_richness$species_richness == 0,
  na.rm = TRUE
)

missing_protocol_checklists <- sum(
  is.na(checklist_richness$protocol_name) |
    checklist_richness$protocol_name == ""
)

expected_city_year_rows <- uniqueN(
  dryland_checklists[
    ,
    paste(
      city_id,
      substr(
        observation_date,
        1,
        4
      ),
      sep = "_"
    )
  ]
)

expected_checklist_files <- city_manifest$checklist_richness_file
expected_city_year_files <- city_manifest$city_year_unique_file

missing_checklist_files <- expected_checklist_files[
  !file.exists(expected_checklist_files)
]

missing_city_year_files <- expected_city_year_files[
  !file.exists(expected_city_year_files)
]


cat("Script 04 input checklist rows:\n")
print(input_checklist_rows)

cat("\nScript 05 checklist richness rows:\n")
print(final_checklist_rows)

cat("\nUnique checklist IDs:\n")
print(final_unique_checklists)

cat("\nCities:\n")
print(final_city_count)

cat("\nUnique locality IDs:\n")
print(final_locality_count)

cat("\nZero-richness complete checklists:\n")
print(zero_richness_checklists)

cat("\nChecklists missing protocol_name:\n")
print(missing_protocol_checklists)

cat("\nExpected city-year combinations:\n")
print(expected_city_year_rows)

cat("\nFinal city-year rows:\n")
print(nrow(city_year_richness))

cat("\nPer-city checklist richness files:\n")
print(length(expected_checklist_files))

cat("\nMissing per-city checklist files:\n")
print(length(missing_checklist_files))

cat("\nPer-city city-year files:\n")
print(length(expected_city_year_files))

cat("\nMissing per-city city-year files:\n")
print(length(missing_city_year_files))

cat("\nUnique locality IDs in GEE hotspot output:\n")
print(uniqueN(hotspot_list$locality_id))

cat("\nScript 04 checklist IDs missing from Script 05:\n")
print(length(missing_checklist_ids))

cat("\nUnexpected Script 05 checklist IDs:\n")
print(length(extra_checklist_ids))


# ------------------------------------------------------------
# 27. Hard final QC checks
# ------------------------------------------------------------

stopifnot(
  file.exists(filtered_ebd_complete_marker),
  file.exists(filtered_ebd_out),
  file.exists(ebd_summaries_complete_marker),
  file.exists(checklist_species_counts_ebd_out),
  file.exists(city_year_species_counts_ebd_out),
  file.exists(protocol_lookup_ebd_out),
  length(expected_checklist_files) == input_city_count,
  all(file.exists(expected_checklist_files)),
  length(expected_city_year_files) == input_city_count,
  all(file.exists(expected_city_year_files)),
  final_checklist_rows == input_checklist_rows,
  final_unique_checklists == input_unique_checklists,
  final_city_count == input_city_count,
  final_locality_count == input_locality_count,
  length(missing_checklist_ids) == 0,
  length(extra_checklist_ids) == 0,
  all(checklist_richness$species_richness >= 0),
  nrow(city_year_richness) == expected_city_year_rows,
  uniqueN(hotspot_list$locality_id) == input_locality_count
)


# ------------------------------------------------------------
# 28. Confirm final outputs exist
# ------------------------------------------------------------

final_outputs <- c(
  filtered_ebd_out,
  filtered_ebd_complete_marker,
  checklist_species_counts_ebd_out,
  city_year_species_counts_ebd_out,
  protocol_lookup_ebd_out,
  ebd_summaries_complete_marker,
  city_manifest_out,
  checklist_city_map_out,
  checklist_metadata_out,
  checklist_richness_out,
  city_year_unique_species_out,
  city_year_richness_out,
  processing_status_out,
  hotspot_list_out,
  hotspot_latlong_out
)


output_check <- data.table(
  file = final_outputs,
  exists = file.exists(final_outputs),
  size_mb = round(
    file.info(final_outputs)$size / 1024^2,
    2
  )
)


cat("\n============================================================\n")
cat("FINAL OUTPUT FILES\n")
cat("============================================================\n\n")

print(output_check)


if (!all(output_check$exists)) {
  stop("One or more final Script 05 outputs were not written.")
}


# ------------------------------------------------------------
# 29. Final result
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SCRIPT 05 QC PASSED\n")
cat("============================================================\n\n")

cat("eBird release:                 ", ebird_release, "\n")
cat("Input Script 04 checklists:    ", input_checklist_rows, "\n")
cat("Final checklist rows:          ", final_checklist_rows, "\n")
cat("Unique checklist IDs:          ", final_unique_checklists, "\n")
cat("Cities:                        ", final_city_count, "\n")
cat("Unique locality IDs:           ", final_locality_count, "\n")
cat("Zero-richness checklists:      ", zero_richness_checklists, "\n")
cat("Missing protocol values:       ", missing_protocol_checklists, "\n")
cat("Final city-year rows:          ", nrow(city_year_richness), "\n")
cat("Missing checklist IDs:         ", length(missing_checklist_ids), "\n")
cat("Unexpected checklist IDs:      ", length(extra_checklist_ids), "\n")


cat("\nMain outputs:\n")

cat(
  "Combined filtered EBD:\n",
  filtered_ebd_out,
  "\n\n"
)

cat(
  "Checklist richness:\n",
  checklist_richness_out,
  "\n\n"
)

cat(
  "City-year richness:\n",
  city_year_richness_out,
  "\n\n"
)

cat(
  "Checklist metadata crosswalk:\n",
  checklist_metadata_out,
  "\n\n"
)

cat(
  "GEE hotspot lat/long:\n",
  hotspot_latlong_out,
  "\n"
)
