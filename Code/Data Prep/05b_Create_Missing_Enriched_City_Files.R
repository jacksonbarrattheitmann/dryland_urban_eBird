# ============================================================
# 05b_Create_Missing_Enriched_City_Files.R
#
# Purpose:
# 1. Create any city-level eBird observation files that are
#    present in Observations_By_City_Raw but missing from
#    Observations_By_City_With_Checklist_Metadata.
#
# 2. Preserve the .x / .y structure expected by Jackson's
#    downstream processing script.
#
# 3. Standardize Stationary checklist effort distance:
#       missing effort_distance_km -> 0
#
# 4. Run final QC across all 647 enriched city files.
#
#
# IMPORTANT:
#
# - Does NOT rerun the national eBird extraction.
# - Does NOT modify Script 05.
# - Does NOT recreate enriched files that already exist.
# - Uses current July 2026 raw city files.
# - Uses current Script 05 checklist metadata.
# - Joins by checklist_id.
#
#
# Jackson-compatible duplicate fields:
#
#   locality.x
#   locality.y
#   locality_id.x
#   locality_id.y
#   observation_date.x
#   observation_date.y
#
#
# Structural QC:
#
# - checklist_id must match metadata.
# - locality_id.x and locality_id.y must agree.
# - observation_date.x and observation_date.y must agree.
#
# locality.x and locality.y may differ in descriptive text.
# These differences are reported but are not treated as errors.
#
#
# Effort-distance convention:
#
# - Missing distance for Stationary protocols is converted to 0.
# - Missing distance for other protocols is left unchanged.
#
#
# Expected final result:
#
#   647 raw city files
#   647 enriched city files
#   0 missing enriched city files
#   0 Stationary checklists with missing distance
#
# ============================================================


# ------------------------------------------------------------
# 1. Package
# ------------------------------------------------------------

library(data.table)


# ------------------------------------------------------------
# 2. File paths
# ------------------------------------------------------------

raw_dir <- file.path(
  "Data",
  "Processed",
  "eBird",
  "Observations_By_City_Raw"
)


enriched_dir <- file.path(
  "Data",
  "Processed",
  "eBird",
  "Observations_By_City_With_Checklist_Metadata"
)


metadata_path <- file.path(
  "Data",
  "Processed",
  "eBird",
  "dryland_city_checklist_metadata_for_joins_2015_2025.csv"
)


manifest_path <- file.path(
  "Data",
  "Processed",
  "eBird",
  "dryland_city_observation_file_manifest_2015_2025.csv"
)


enrichment_status_out <- file.path(
  "Data",
  "Processed",
  "eBird",
  "missing_city_enrichment_status_2015_2025.csv"
)


stationary_status_out <- file.path(
  "Data",
  "Processed",
  "eBird",
  "stationary_effort_distance_correction_status_2015_2025.csv"
)


effort_qc_out <- file.path(
  "Data",
  "Processed",
  "eBird",
  "enriched_city_effort_distance_qc_2015_2025.csv"
)


dir.create(
  enriched_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Helper function for safe replacement
# ------------------------------------------------------------

safe_replace_tsv <- function(x, final_path) {
  
  temp_path <- paste0(
    final_path,
    ".tmp"
  )
  
  
  backup_path <- paste0(
    final_path,
    ".bak_tmp"
  )
  
  
  if (file.exists(temp_path)) {
    unlink(temp_path, force = TRUE)
  }
  
  
  if (file.exists(backup_path)) {
    unlink(backup_path, force = TRUE)
  }
  
  
  fwrite(
    x,
    temp_path,
    sep = "\t"
  )
  
  
  if (
    !file.exists(temp_path) ||
    file.info(temp_path)$size == 0
  ) {
    
    stop(
      "Temporary TSV was not written correctly:\n",
      temp_path
    )
  }
  
  
  # If replacing an existing file, temporarily move the
  # original file out of the way.
  
  if (file.exists(final_path)) {
    
    if (
      !file.rename(
        final_path,
        backup_path
      )
    ) {
      
      unlink(
        temp_path,
        force = TRUE
      )
      
      stop(
        "Could not create temporary backup of:\n",
        final_path
      )
    }
  }
  
  
  # Move completed temporary file into final location.
  
  if (
    !file.rename(
      temp_path,
      final_path
    )
  ) {
    
    # Restore backup if replacement failed.
    
    if (file.exists(backup_path)) {
      
      file.rename(
        backup_path,
        final_path
      )
    }
    
    
    stop(
      "Could not finalize TSV:\n",
      final_path
    )
  }
  
  
  # Remove backup after successful replacement.
  
  if (file.exists(backup_path)) {
    
    unlink(
      backup_path,
      force = TRUE
    )
  }
}


# ------------------------------------------------------------
# 4. Confirm required inputs exist
# ------------------------------------------------------------

if (!dir.exists(raw_dir)) {
  
  stop(
    "Raw city folder does not exist:\n",
    raw_dir
  )
}


if (!file.exists(metadata_path)) {
  
  stop(
    "Checklist metadata crosswalk not found:\n",
    metadata_path
  )
}


if (!file.exists(manifest_path)) {
  
  stop(
    "Current Script 05 city manifest not found:\n",
    manifest_path
  )
}


# ------------------------------------------------------------
# 5. Read current city manifest
# ------------------------------------------------------------

manifest <- fread(
  manifest_path,
  showProgress = FALSE
)


required_manifest_cols <- c(
  "city_id",
  "city_file_stub"
)


missing_manifest_cols <- setdiff(
  required_manifest_cols,
  names(manifest)
)


if (length(missing_manifest_cols) > 0) {
  
  stop(
    "Manifest is missing required columns:\n",
    paste(
      missing_manifest_cols,
      collapse = "\n"
    )
  )
}


if (
  anyDuplicated(
    manifest$city_file_stub
  ) > 0
) {
  
  stop(
    "Duplicate city_file_stub values detected in manifest."
  )
}


expected_files <- paste0(
  manifest$city_file_stub,
  ".tsv"
)


cat("\n============================================================\n")
cat("CURRENT CITY FILE INVENTORY\n")
cat("============================================================\n\n")


cat("Expected cities from Script 05 manifest:\n")
print(length(expected_files))


# ------------------------------------------------------------
# 6. Inventory raw and enriched city files
# ------------------------------------------------------------

raw_files <- basename(
  list.files(
    raw_dir,
    pattern = "\\.tsv$",
    full.names = TRUE
  )
)


enriched_files <- basename(
  list.files(
    enriched_dir,
    pattern = "\\.tsv$",
    full.names = TRUE
  )
)


missing_raw <- setdiff(
  expected_files,
  raw_files
)


missing_enriched <- setdiff(
  expected_files,
  enriched_files
)


extra_enriched <- setdiff(
  enriched_files,
  expected_files
)


cat("\nRaw city files present:\n")
print(length(raw_files))


cat("\nExisting enriched city files:\n")
print(length(enriched_files))


cat("\nMissing raw city files:\n")
print(length(missing_raw))


cat("\nMissing enriched city files:\n")
print(length(missing_enriched))


cat("\nUnexpected enriched city files:\n")
print(length(extra_enriched))


if (length(missing_raw) > 0) {
  
  cat("\nMissing raw city files:\n")
  print(missing_raw)
  
  stop(
    "One or more expected raw city files are missing."
  )
}


if (length(extra_enriched) > 0) {
  
  cat("\nUnexpected enriched city files:\n")
  print(extra_enriched)
  
  stop(
    "Unexpected city files exist in enriched folder."
  )
}


if (length(missing_enriched) > 0) {
  
  cat("\nCities requiring enrichment:\n")
  print(missing_enriched)
  
} else {
  
  cat(
    "\nAll expected enriched city files already exist.\n"
  )
}


# ------------------------------------------------------------
# 7. Read current Script 05 checklist metadata
# ------------------------------------------------------------

metadata <- fread(
  metadata_path,
  showProgress = TRUE
)


required_metadata_cols <- c(
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
  "latitude_checklist",
  "longitude_checklist",
  "duration_minutes",
  "effort_distance_km",
  "number_observers",
  "city_file_stub"
)


missing_metadata_cols <- setdiff(
  required_metadata_cols,
  names(metadata)
)


if (length(missing_metadata_cols) > 0) {
  
  stop(
    "Checklist metadata crosswalk is missing required columns:\n",
    paste(
      missing_metadata_cols,
      collapse = "\n"
    )
  )
}


metadata[
  ,
  checklist_id := as.character(checklist_id)
]


cat("\n============================================================\n")
cat("CHECKLIST METADATA QC\n")
cat("============================================================\n\n")


cat("Metadata rows:\n")
print(nrow(metadata))


cat("\nUnique checklist IDs:\n")
print(uniqueN(metadata$checklist_id))


cat("\nCities in metadata:\n")
print(uniqueN(metadata$city_file_stub))


cat("\nDryland-percent range:\n")

print(
  range(
    metadata$dryland_percent,
    na.rm = TRUE
  )
)


if (
  anyDuplicated(
    metadata$checklist_id
  ) > 0
) {
  
  stop(
    "Checklist metadata contains duplicate checklist IDs."
  )
}


if (
  uniqueN(
    metadata$city_file_stub
  ) != 647
) {
  
  stop(
    "Expected metadata for exactly 647 cities."
  )
}


if (
  any(
    metadata$dryland_percent < 50 - 1e-8 |
    metadata$dryland_percent > 100 + 1e-8,
    na.rm = TRUE
  )
) {
  
  stop(
    "Invalid dryland_percent values detected."
  )
}


# ------------------------------------------------------------
# 8. Required raw observation columns
# ------------------------------------------------------------

required_raw_cols <- c(
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


# ------------------------------------------------------------
# 9. Required final columns for Jackson
# ------------------------------------------------------------

required_jackson_cols <- c(
  "checklist_id",
  "observation_date.x",
  "locality_id.x",
  "locality.x",
  "common_name",
  "scientific_name",
  "category",
  "observation_count",
  "protocol_name",
  "observation_date.y",
  "locality_id.y",
  "locality.y",
  "city_id",
  "NAME10",
  "city_ascii",
  "state_name",
  "county_name",
  "dryland_percent",
  "duration_minutes",
  "effort_distance_km",
  "number_observers"
)


# ------------------------------------------------------------
# 10. Create missing enriched city files
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("CREATING MISSING ENRICHED CITY FILES\n")
cat("============================================================\n\n")


if (length(missing_enriched) == 0) {
  
  cat(
    "No enriched city files are missing.\n"
  )
  
  
  enrichment_status <- data.table()
  
  
} else {
  
  status_list <- vector(
    "list",
    length(missing_enriched)
  )
  
  
  for (i in seq_along(missing_enriched)) {
    
    file_name <- missing_enriched[i]
    
    
    city_stub <- tools::file_path_sans_ext(
      file_name
    )
    
    
    raw_path <- file.path(
      raw_dir,
      file_name
    )
    
    
    enriched_path <- file.path(
      enriched_dir,
      file_name
    )
    
    
    cat(
      "\n------------------------------------------------------------\n"
    )
    
    
    cat(
      "Processing ",
      i,
      " of ",
      length(missing_enriched),
      ": ",
      file_name,
      "\n",
      sep = ""
    )
    
    
    # --------------------------------------------------------
    # 10.1 Read raw city file
    # --------------------------------------------------------
    
    raw <- fread(
      raw_path,
      showProgress = TRUE
    )
    
    
    missing_raw_cols_i <- setdiff(
      required_raw_cols,
      names(raw)
    )
    
    
    if (length(missing_raw_cols_i) > 0) {
      
      stop(
        "Raw city file is missing required columns:\n",
        file_name,
        "\n\n",
        paste(
          missing_raw_cols_i,
          collapse = "\n"
        )
      )
    }
    
    
    raw[
      ,
      checklist_id := as.character(checklist_id)
    ]
    
    
    raw_rows_before <- nrow(raw)
    
    
    raw_unique_checklists <- uniqueN(
      raw$checklist_id
    )
    
    
    # --------------------------------------------------------
    # 10.2 Get metadata for city
    # --------------------------------------------------------
    
    metadata_city <- metadata[
      city_file_stub == city_stub
    ]
    
    
    if (nrow(metadata_city) == 0) {
      
      stop(
        "No checklist metadata found for:\n",
        city_stub
      )
    }
    
    
    if (
      anyDuplicated(
        metadata_city$checklist_id
      ) > 0
    ) {
      
      stop(
        "Duplicate checklist IDs in metadata for:\n",
        city_stub
      )
    }
    
    
    raw_ids_missing_metadata <- setdiff(
      unique(raw$checklist_id),
      metadata_city$checklist_id
    )
    
    
    if (
      length(
        raw_ids_missing_metadata
      ) > 0
    ) {
      
      stop(
        "Raw observation checklist IDs missing metadata for:\n",
        city_stub,
        "\nNumber missing: ",
        length(raw_ids_missing_metadata)
      )
    }
    
    
    metadata_ids_without_raw <- setdiff(
      metadata_city$checklist_id,
      unique(raw$checklist_id)
    )
    
    
    # --------------------------------------------------------
    # 10.3 Join raw observations to checklist metadata
    # --------------------------------------------------------
    
    enriched <- merge(
      raw,
      metadata_city,
      by = "checklist_id",
      all.x = TRUE,
      sort = FALSE,
      suffixes = c(
        ".x",
        ".y"
      )
    )
    
    
    # --------------------------------------------------------
    # 10.4 Row-count QC
    # --------------------------------------------------------
    
    if (
      nrow(enriched) !=
      raw_rows_before
    ) {
      
      stop(
        "Observation row count changed during join for:\n",
        city_stub,
        "\nRaw rows: ",
        raw_rows_before,
        "\nEnriched rows: ",
        nrow(enriched)
      )
    }
    
    
    # --------------------------------------------------------
    # 10.5 Required-column QC
    # --------------------------------------------------------
    
    missing_jackson_cols <- setdiff(
      required_jackson_cols,
      names(enriched)
    )
    
    
    if (
      length(
        missing_jackson_cols
      ) > 0
    ) {
      
      stop(
        "Enriched file is missing Jackson-required columns:\n",
        city_stub,
        "\n\n",
        paste(
          missing_jackson_cols,
          collapse = "\n"
        )
      )
    }
    
    
    # --------------------------------------------------------
    # 10.6 Standardize Stationary distance
    # --------------------------------------------------------
    #
    # Stationary checklists have zero distance traveled.
    #
    # Only missing values for Stationary protocols are changed.
    # Other protocol types are left exactly as recorded.
    # --------------------------------------------------------
    
    stationary_missing <- (
      !is.na(enriched$protocol_name) &
        grepl(
          "^Stationary",
          enriched$protocol_name
        ) &
        is.na(enriched$effort_distance_km)
    )
    
    
    stationary_rows_corrected <- sum(
      stationary_missing
    )
    
    
    stationary_checklists_corrected <- uniqueN(
      enriched$checklist_id[
        stationary_missing
      ]
    )
    
    
    enriched[
      stationary_missing,
      effort_distance_km := 0
    ]
    
    
    # --------------------------------------------------------
    # 10.7 Metadata completeness
    # --------------------------------------------------------
    
    missing_duration <- sum(
      is.na(
        enriched$duration_minutes
      )
    )
    
    
    missing_observers <- sum(
      is.na(
        enriched$number_observers
      )
    )
    
    
    missing_city_name <- sum(
      is.na(enriched$NAME10) |
        enriched$NAME10 == ""
    )
    
    
    missing_dryland <- sum(
      is.na(
        enriched$dryland_percent
      )
    )
    
    
    missing_distance_after <- sum(
      is.na(
        enriched$effort_distance_km
      )
    )
    
    
    if (
      missing_duration > 0 ||
      missing_observers > 0 ||
      missing_city_name > 0 ||
      missing_dryland > 0
    ) {
      
      stop(
        "Required checklist metadata failed to join for:\n",
        city_stub,
        "\nMissing duration: ",
        missing_duration,
        "\nMissing observers: ",
        missing_observers,
        "\nMissing city name: ",
        missing_city_name,
        "\nMissing dryland percent: ",
        missing_dryland
      )
    }
    
    
    # --------------------------------------------------------
    # 10.8 Validate duplicate structural fields
    # --------------------------------------------------------
    
    locality_id_mismatch <- sum(
      
      !is.na(enriched$locality_id.x) &
        !is.na(enriched$locality_id.y) &
        enriched$locality_id.x !=
        enriched$locality_id.y,
      
      na.rm = TRUE
    )
    
    
    locality_name_difference <- sum(
      
      !is.na(enriched$locality.x) &
        !is.na(enriched$locality.y) &
        enriched$locality.x !=
        enriched$locality.y,
      
      na.rm = TRUE
    )
    
    
    date_mismatch <- sum(
      
      !is.na(enriched$observation_date.x) &
        !is.na(enriched$observation_date.y) &
        enriched$observation_date.x !=
        enriched$observation_date.y,
      
      na.rm = TRUE
    )
    
    
    if (
      locality_id_mismatch > 0 ||
      date_mismatch > 0
    ) {
      
      stop(
        "Raw and metadata structural fields disagree for:\n",
        city_stub,
        "\nlocality_id mismatches: ",
        locality_id_mismatch,
        "\ndate mismatches: ",
        date_mismatch
      )
    }
    
    
    # --------------------------------------------------------
    # 10.9 Write file
    # --------------------------------------------------------
    
    if (file.exists(enriched_path)) {
      
      stop(
        "Refusing to overwrite existing enriched file:\n",
        enriched_path
      )
    }
    
    
    safe_replace_tsv(
      enriched,
      enriched_path
    )
    
    
    # --------------------------------------------------------
    # 10.10 Save city status
    # --------------------------------------------------------
    
    status_list[[i]] <- data.table(
      
      city_file_stub =
        city_stub,
      
      raw_rows =
        raw_rows_before,
      
      enriched_rows =
        nrow(enriched),
      
      raw_unique_checklists =
        raw_unique_checklists,
      
      metadata_checklists =
        nrow(metadata_city),
      
      metadata_checklists_without_raw_rows =
        length(metadata_ids_without_raw),
      
      missing_duration =
        missing_duration,
      
      missing_number_observers =
        missing_observers,
      
      missing_effort_distance_after =
        missing_distance_after,
      
      stationary_rows_corrected =
        stationary_rows_corrected,
      
      stationary_checklists_corrected =
        stationary_checklists_corrected,
      
      locality_id_mismatches =
        locality_id_mismatch,
      
      locality_name_differences =
        locality_name_difference,
      
      date_mismatches =
        date_mismatch,
      
      enriched_file =
        enriched_path
    )
    
    
    cat(
      "  Raw rows: ",
      raw_rows_before,
      "\n",
      sep = ""
    )
    
    
    cat(
      "  Unique raw checklists: ",
      raw_unique_checklists,
      "\n",
      sep = ""
    )
    
    
    cat(
      "  Stationary checklists corrected: ",
      stationary_checklists_corrected,
      "\n",
      sep = ""
    )
    
    
    cat(
      "  locality_id mismatches: ",
      locality_id_mismatch,
      "\n",
      sep = ""
    )
    
    
    cat(
      "  Locality-name differences: ",
      locality_name_difference,
      "\n",
      sep = ""
    )
    
    
    cat(
      "  Date mismatches: ",
      date_mismatch,
      "\n",
      sep = ""
    )
    
    
    cat(
      "  Written successfully.\n"
    )
    
    
    rm(
      raw,
      metadata_city,
      enriched,
      metadata_ids_without_raw
    )
    
    
    gc()
  }
  
  
  enrichment_status <- rbindlist(
    status_list,
    fill = TRUE
  )
  
  
  fwrite(
    enrichment_status,
    enrichment_status_out
  )
}


# ------------------------------------------------------------
# 11. Confirm all 647 enriched files are now present
# ------------------------------------------------------------

final_enriched_files <- basename(
  list.files(
    enriched_dir,
    pattern = "\\.tsv$",
    full.names = TRUE
  )
)


final_missing <- setdiff(
  expected_files,
  final_enriched_files
)


final_extra <- setdiff(
  final_enriched_files,
  expected_files
)


if (
  length(final_missing) > 0 ||
  length(final_extra) > 0
) {
  
  cat("\nStill missing:\n")
  print(final_missing)
  
  
  cat("\nUnexpected:\n")
  print(final_extra)
  
  
  stop(
    "Final enriched city inventory failed."
  )
}


# ------------------------------------------------------------
# 12. Standardize Stationary distance across ALL 647 files
# ------------------------------------------------------------
#
# This makes the script rerunnable.
#
# It also checks/fixes the older 631 enriched files so the
# entire Jackson handoff uses the same distance convention.
#
# Files are only rewritten if they actually contain a
# Stationary checklist with missing effort_distance_km.
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("STANDARDIZING STATIONARY EFFORT DISTANCE\n")
cat("============================================================\n\n")


all_enriched_paths <- file.path(
  enriched_dir,
  expected_files
)


stationary_status_list <- vector(
  "list",
  length(all_enriched_paths)
)


for (i in seq_along(all_enriched_paths)) {
  
  f <- all_enriched_paths[i]
  
  
  cat(
    "Checking ",
    i,
    " of ",
    length(all_enriched_paths),
    ": ",
    basename(f),
    "\n",
    sep = ""
  )
  
  
  # Read only columns needed to determine whether
  # a rewrite is necessary.
  
  check <- fread(
    f,
    select = c(
      "checklist_id",
      "protocol_name",
      "effort_distance_km"
    ),
    showProgress = FALSE
  )
  
  
  stationary_missing <- (
    !is.na(check$protocol_name) &
      grepl(
        "^Stationary",
        check$protocol_name
      ) &
      is.na(check$effort_distance_km)
  )
  
  
  rows_to_fix <- sum(
    stationary_missing
  )
  
  
  checklists_to_fix <- uniqueN(
    check$checklist_id[
      stationary_missing
    ]
  )
  
  
  rewritten <- FALSE
  
  
  if (rows_to_fix > 0) {
    
    cat(
      "  Stationary checklists requiring correction: ",
      checklists_to_fix,
      "\n",
      sep = ""
    )
    
    
    x <- fread(
      f,
      showProgress = TRUE
    )
    
    
    idx <- (
      !is.na(x$protocol_name) &
        grepl(
          "^Stationary",
          x$protocol_name
        ) &
        is.na(x$effort_distance_km)
    )
    
    
    x[
      idx,
      effort_distance_km := 0
    ]
    
    
    # Verify correction before writing.
    
    remaining_stationary_missing <- sum(
      
      !is.na(x$protocol_name) &
        grepl(
          "^Stationary",
          x$protocol_name
        ) &
        is.na(x$effort_distance_km)
    )
    
    
    if (
      remaining_stationary_missing > 0
    ) {
      
      stop(
        "Stationary correction failed for:\n",
        f
      )
    }
    
    
    safe_replace_tsv(
      x,
      f
    )
    
    
    rewritten <- TRUE
    
    
    rm(x)
    gc()
  }
  
  
  stationary_status_list[[i]] <- data.table(
    
    file =
      basename(f),
    
    stationary_rows_corrected =
      rows_to_fix,
    
    stationary_checklists_corrected =
      checklists_to_fix,
    
    file_rewritten =
      rewritten
  )
  
  
  rm(check)
}


stationary_status <- rbindlist(
  stationary_status_list,
  fill = TRUE
)


fwrite(
  stationary_status,
  stationary_status_out
)


cat("\nStationary correction summary:\n")


cat(
  "Files rewritten: ",
  sum(stationary_status$file_rewritten),
  "\n",
  sep = ""
)


cat(
  "Stationary checklists corrected: ",
  sum(stationary_status$stationary_checklists_corrected),
  "\n",
  sep = ""
)


cat(
  "Stationary observation rows corrected: ",
  sum(stationary_status$stationary_rows_corrected),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 13. Header QC across all 647 enriched files
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ALL-CITY JACKSON HEADER QC\n")
cat("============================================================\n\n")


bad_headers_list <- list()


for (i in seq_along(all_enriched_paths)) {
  
  f <- all_enriched_paths[i]
  
  
  h <- names(
    fread(
      f,
      nrows = 0,
      showProgress = FALSE
    )
  )
  
  
  missing_cols_i <- setdiff(
    required_jackson_cols,
    h
  )
  
  
  if (
    length(
      missing_cols_i
    ) > 0
  ) {
    
    bad_headers_list[[
      length(bad_headers_list) + 1L
    ]] <- data.table(
      
      file =
        basename(f),
      
      problem =
        paste(
          missing_cols_i,
          collapse = "; "
        )
    )
  }
  
  
  if (
    i %% 100 == 0 ||
    i == length(all_enriched_paths)
  ) {
    
    cat(
      "Headers checked: ",
      i,
      " of ",
      length(all_enriched_paths),
      "\n",
      sep = ""
    )
  }
}


if (
  length(
    bad_headers_list
  ) > 0
) {
  
  bad_headers <- rbindlist(
    bad_headers_list,
    fill = TRUE
  )
  
  
} else {
  
  bad_headers <- data.table(
    file = character(),
    problem = character()
  )
}


cat(
  "\nFiles failing Jackson header QC: ",
  nrow(bad_headers),
  "\n",
  sep = ""
)


if (
  nrow(
    bad_headers
  ) > 0
) {
  
  print(
    bad_headers,
    nrows = Inf
  )
}


# ------------------------------------------------------------
# 14. Final effort-distance QC across ALL 647 cities
# ------------------------------------------------------------
#
# Work at checklist level rather than observation-row level.
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("ALL-CITY EFFORT-DISTANCE QC\n")
cat("============================================================\n\n")


effort_qc_list <- vector(
  "list",
  length(all_enriched_paths)
)


for (i in seq_along(all_enriched_paths)) {
  
  f <- all_enriched_paths[i]
  
  
  x <- fread(
    f,
    select = c(
      "checklist_id",
      "protocol_name",
      "effort_distance_km"
    ),
    showProgress = FALSE
  )
  
  
  x <- unique(
    x,
    by = "checklist_id"
  )
  
  
  effort_qc_list[[i]] <- x[
    ,
    .(
      n_checklists = .N,
      n_missing_distance =
        sum(
          is.na(
            effort_distance_km
          )
        )
    ),
    by = protocol_name
  ]
  
  
  rm(x)
  
  
  if (
    i %% 100 == 0 ||
    i == length(all_enriched_paths)
  ) {
    
    cat(
      "Effort QC checked: ",
      i,
      " of ",
      length(all_enriched_paths),
      "\n",
      sep = ""
    )
  }
}


effort_qc <- rbindlist(
  effort_qc_list,
  fill = TRUE
)


effort_qc_summary <- effort_qc[
  ,
  .(
    n_checklists =
      sum(n_checklists),
    
    n_missing_distance =
      sum(n_missing_distance)
  ),
  by = protocol_name
]


effort_qc_summary[
  ,
  pct_missing_distance :=
    round(
      100 *
        n_missing_distance /
        n_checklists,
      2
    )
]


setorder(
  effort_qc_summary,
  -n_missing_distance
)


print(
  effort_qc_summary,
  nrows = Inf
)


fwrite(
  effort_qc_summary,
  effort_qc_out
)


# ------------------------------------------------------------
# 15. Hard QC for Stationary protocols
# ------------------------------------------------------------

stationary_final_qc <- effort_qc_summary[
  
  !is.na(protocol_name) &
    grepl(
      "^Stationary",
      protocol_name
    )
]


stationary_missing_final <- sum(
  stationary_final_qc$n_missing_distance
)


cat(
  "\nStationary checklists with missing distance after correction: ",
  stationary_missing_final,
  "\n",
  sep = ""
)


if (
  stationary_missing_final > 0
) {
  
  stop(
    "One or more Stationary checklists still has missing distance."
  )
}


# Traveling missing values are reported but NOT altered.

traveling_qc <- effort_qc_summary[
  protocol_name == "Traveling"
]


if (
  nrow(
    traveling_qc
  ) > 0
) {
  
  cat(
    "\nTraveling checklists with missing distance: ",
    traveling_qc$n_missing_distance,
    " of ",
    traveling_qc$n_checklists,
    " (",
    traveling_qc$pct_missing_distance,
    "%)\n",
    sep = ""
  )
}


# ------------------------------------------------------------
# 16. Final inventory
# ------------------------------------------------------------

final_enriched_files <- basename(
  list.files(
    enriched_dir,
    pattern = "\\.tsv$",
    full.names = TRUE
  )
)


final_missing <- setdiff(
  expected_files,
  final_enriched_files
)


final_extra <- setdiff(
  final_enriched_files,
  expected_files
)


# ------------------------------------------------------------
# 17. Final hard QC
# ------------------------------------------------------------

stopifnot(
  
  length(expected_files) == 647,
  
  length(raw_files) == 647,
  
  length(final_enriched_files) == 647,
  
  length(final_missing) == 0,
  
  length(final_extra) == 0,
  
  nrow(bad_headers) == 0,
  
  stationary_missing_final == 0
)


# ------------------------------------------------------------
# 18. Final result
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SCRIPT 05b QC PASSED\n")
cat("============================================================\n\n")


cat(
  "Expected cities:                    ",
  length(expected_files),
  "\n",
  sep = ""
)


cat(
  "Raw city files:                     ",
  length(raw_files),
  "\n",
  sep = ""
)


cat(
  "Enriched city files:                ",
  length(final_enriched_files),
  "\n",
  sep = ""
)


cat(
  "Missing enriched files:             ",
  length(final_missing),
  "\n",
  sep = ""
)


cat(
  "Unexpected enriched files:          ",
  length(final_extra),
  "\n",
  sep = ""
)


cat(
  "Metadata-header failures:           ",
  nrow(bad_headers),
  "\n",
  sep = ""
)


cat(
  "Stationary missing distance:        ",
  stationary_missing_final,
  "\n",
  sep = ""
)


cat(
  "Files corrected for Stationary NA:  ",
  sum(stationary_status$file_rewritten),
  "\n",
  sep = ""
)


cat(
  "\nJackson input folder:\n",
  enriched_dir,
  "\n",
  sep = ""
)