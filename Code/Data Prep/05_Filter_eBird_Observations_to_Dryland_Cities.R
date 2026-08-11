# ============================================================
# 05_Filter_eBird_Observations_to_Dryland_Cities_WITH_HOTSPOTS.R
#
# Purpose:
# Use the dryland-city checklist metadata from Script 04 to filter
# the full eBird Basic Dataset and split retained observation records
# into one raw observation file per dryland city.
#
# This version is designed for long runs and safe recovery:
#   1. Raw city observation files are written by AWK.
#   2. Checklist richness is written as one file per city.
#   3. City-year unique species richness is written as one file per city.
#   4. If the script stops halfway through, rerun it and it skips completed cities.
#   5. locality_id and locality are retained for each checklist.
#   6. Combined outputs are rebuilt at the end from the per-city files.
#   7. Hotspot/locality lists are written for GEE extraction.
# ============================================================

library(data.table)
library(dplyr)
library(readr)
library(stringr)

# ------------------------------------------------------------
# macOS / RStudio warning control
# ------------------------------------------------------------

# These warnings can appear repeatedly on Apple Silicon Macs when RStudio
# launches external commands. They are usually harmless, but can flood the
# console and make long runs difficult to monitor.
Sys.unsetenv("MallocStackLogging")
Sys.unsetenv("MallocStackLoggingNoCompact")
Sys.setenv(MallocStackLogging = "0")
Sys.setenv(MallocStackLoggingNoCompact = "0")

# ------------------------------------------------------------
# 0. Run settings
# ------------------------------------------------------------

# Keep FALSE unless you intentionally want to delete raw city files
# and rerun the full expensive AWK split.
reset_raw_city_outputs <- FALSE

# TRUE deletes only summary outputs, not raw city files.
# Use TRUE if you changed the summary logic and want to rebuild city summaries.
reset_summary_outputs <- FALSE

# If raw city files already exist, skip the AWK split and go straight
# to city-level summary processing.
skip_awk_if_raw_city_files_exist <- TRUE

# Creating metadata-expanded observation files duplicates checklist metadata
# onto every observation row and is usually not needed.
write_metadata_observation_files <- TRUE

# Optional helper: run this in a separate R session while Script 05 is running
# to check whether raw files are still being written.
check_raw_split_progress <- function() {
  raw_files <- list.files(
    "Data/Processed/eBird/Observations_By_City_Raw",
    pattern = "\\.tsv$",
    full.names = TRUE
  )
  
  info <- file.info(raw_files)
  
  cat("Number of raw city files:", length(raw_files), "\n")
  cat("Total size GB:", round(sum(info$size, na.rm = TRUE) / 1024^3, 2), "\n")
  cat("Most recent modified:", max(info$mtime, na.rm = TRUE), "\n")
  
  invisible(
    data.frame(
      file = basename(raw_files),
      size_mb = round(info$size / 1024^2, 2),
      modified = info$mtime
    )[order(-info$mtime), ]
  )
}

# ------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------

tar_path <- "Data/Raw/eBird/ebd_US_smp_relApr-2026.tar"
ebd_member <- "ebd_US_relApr-2026.txt.gz"
checklists_path <- "Data/Processed/eBird/dryland_city_checklists_2015_2025.csv"

ebird_out_dir <- "Data/Processed/eBird"
dir.create(ebird_out_dir, recursive = TRUE, showWarnings = FALSE)

city_raw_out_dir <- file.path(ebird_out_dir, "Observations_By_City_Raw")
city_metadata_out_dir <- file.path(ebird_out_dir, "Observations_By_City_With_Checklist_Metadata")
checklist_richness_by_city_dir <- file.path(ebird_out_dir, "Checklist_Richness_By_City")
city_year_unique_by_city_dir <- file.path(ebird_out_dir, "City_Year_Unique_Species_By_City")

dir.create(city_raw_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(city_metadata_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checklist_richness_by_city_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(city_year_unique_by_city_dir, recursive = TRUE, showWarnings = FALSE)

city_manifest_out <- file.path(ebird_out_dir, "dryland_city_observation_file_manifest_2015_2025.csv")
checklist_city_map_out <- file.path(ebird_out_dir, "dryland_city_checklist_to_city_file_map_2015_2025.tsv")
checklist_metadata_out <- file.path(ebird_out_dir, "dryland_city_checklist_metadata_for_joins_2015_2025.csv")

checklist_richness_out <- file.path(ebird_out_dir, "dryland_city_checklist_species_richness_2015_2025.csv")
city_year_unique_species_out <- file.path(ebird_out_dir, "dryland_city_year_unique_species_from_observations_2015_2025.csv")
city_year_richness_out <- file.path(ebird_out_dir, "dryland_city_year_species_richness_2015_2025.csv")
processing_status_out <- file.path(ebird_out_dir, "dryland_city_processing_status_2015_2025.csv")
hotspot_list_out <- file.path(ebird_out_dir, "dryland_city_hotspot_list_for_GEE_2015_2025.csv")
hotspot_latlong_out <- file.path(ebird_out_dir, "dryland_city_hotspot_latlong_for_GEE_2015_2025.csv")

awk_script_path <- file.path(ebird_out_dir, "split_ebd_observations_by_city.awk")
awk_log_path <- file.path(ebird_out_dir, "split_ebd_observations_by_city_log.txt")
awk_quiet_log_path <- file.path(ebird_out_dir, "split_ebd_observations_by_city_console_messages.txt")

# ------------------------------------------------------------
# 1.1 Optional resets
# ------------------------------------------------------------

if (reset_raw_city_outputs) {
  cat("\nResetting raw city outputs and all dependent outputs...\n")
  unlink(
    c(
      list.files(city_raw_out_dir, full.names = TRUE, pattern = "\\.tsv$"),
      list.files(city_metadata_out_dir, full.names = TRUE, pattern = "\\.tsv$"),
      list.files(checklist_richness_by_city_dir, full.names = TRUE, pattern = "\\.csv$"),
      list.files(city_year_unique_by_city_dir, full.names = TRUE, pattern = "\\.csv$"),
      city_manifest_out,
      checklist_city_map_out,
      checklist_metadata_out,
      checklist_richness_out,
      city_year_unique_species_out,
      city_year_richness_out,
      processing_status_out,
      hotspot_list_out,
      hotspot_latlong_out,
      awk_script_path,
      awk_log_path,
      awk_quiet_log_path
    ),
    force = TRUE
  )
}

if (reset_summary_outputs && !reset_raw_city_outputs) {
  cat("\nResetting summary outputs only. Raw city files are kept.\n")
  unlink(
    c(
      list.files(checklist_richness_by_city_dir, full.names = TRUE, pattern = "\\.csv$"),
      list.files(city_year_unique_by_city_dir, full.names = TRUE, pattern = "\\.csv$"),
      checklist_richness_out,
      city_year_unique_species_out,
      city_year_richness_out,
      processing_status_out,
      hotspot_list_out,
      hotspot_latlong_out
    ),
    force = TRUE
  )
}

# ------------------------------------------------------------
# 2. Read checklist metadata from Script 04
# ------------------------------------------------------------

if (!file.exists(checklists_path)) {
  stop("Checklist metadata file not found: ", checklists_path)
}

dryland_checklists <- fread(checklists_path)

cat("\nDryland checklist metadata rows loaded:\n")
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

missing_checklist_cols <- setdiff(required_checklist_cols, names(dryland_checklists))

if (length(missing_checklist_cols) > 0) {
  stop(
    "The checklist metadata file is missing required columns: ",
    paste(missing_checklist_cols, collapse = ", ")
  )
}

cat("\nRows with locality_id in checklist metadata:\n")
print(sum(!is.na(dryland_checklists$locality_id) & dryland_checklists$locality_id != ""))

cat("\nUnique locality IDs in checklist metadata:\n")
print(uniqueN(dryland_checklists$locality_id[!is.na(dryland_checklists$locality_id) & dryland_checklists$locality_id != ""]))

# ------------------------------------------------------------
# 3. Create safe city file names and checklist-to-city map
# ------------------------------------------------------------

clean_file_part <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("[^a-z0-9]+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

dryland_checklists <- dryland_checklists %>%
  as_tibble() %>%
  mutate(
    city_file_stub = paste0(
      "city_",
      stringr::str_pad(city_id, width = 4, pad = "0"),
      "_",
      clean_file_part(state_name),
      "_",
      clean_file_part(city_ascii)
    )
  ) %>%
  as.data.table()

city_manifest <- dryland_checklists %>%
  as_tibble() %>%
  distinct(city_id, NAME10, city_ascii, state_name, county_name, dryland_percent, city_file_stub) %>%
  arrange(state_name, city_ascii) %>%
  mutate(
    raw_observation_file = file.path(city_raw_out_dir, paste0(city_file_stub, ".tsv")),
    metadata_observation_file = file.path(city_metadata_out_dir, paste0(city_file_stub, ".tsv")),
    checklist_richness_file = file.path(checklist_richness_by_city_dir, paste0(city_file_stub, "_checklist_richness.csv")),
    city_year_unique_file = file.path(city_year_unique_by_city_dir, paste0(city_file_stub, "_city_year_unique_species.csv"))
  )

write_csv(city_manifest, city_manifest_out)

checklist_city_map <- dryland_checklists %>%
  as_tibble() %>%
  distinct(checklist_id, city_file_stub)

write_tsv(checklist_city_map, checklist_city_map_out)

checklist_metadata_small <- dryland_checklists %>%
  as_tibble() %>%
  select(
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
  ) %>%
  distinct() %>%
  as.data.table()

fwrite(checklist_metadata_small, checklist_metadata_out)

cat("\nDryland cities in checklist metadata:\n")
print(nrow(city_manifest))

cat("\nUnique checklist IDs to filter:\n")
print(n_distinct(checklist_city_map$checklist_id))

cat("\nChecklist metadata crosswalk written to:\n")
cat(checklist_metadata_out, "\n")

# ------------------------------------------------------------
# 4. Confirm EBD file exists inside TAR
# ------------------------------------------------------------

if (!file.exists(tar_path)) {
  stop("eBird TAR file not found: ", tar_path)
}

cat("\nFiles inside eBird TAR:\n")

tar_list_cmd <- paste("/usr/bin/tar -tf", shQuote(tar_path))
tar_files <- system(tar_list_cmd, intern = TRUE)

print(tar_files)

if (!ebd_member %in% tar_files) {
  stop("The EBD member file was not found inside the TAR archive: ", ebd_member)
}

# ------------------------------------------------------------
# 5. Raw city observation split using AWK
# ------------------------------------------------------------

existing_raw_city_files <- list.files(city_raw_out_dir, full.names = TRUE, pattern = "\\.tsv$")

if (
  skip_awk_if_raw_city_files_exist &&
  length(existing_raw_city_files) > 0 &&
  !reset_raw_city_outputs
) {
  
  cat("\nSkipping AWK split because raw city files already exist.\n")
  cat("Existing raw city files:\n")
  print(length(existing_raw_city_files))
  
} else {
  
  cat("\nWriting AWK script for EBD filtering and city splitting...\n")
  
  awk_script <- '
BEGIN {
  FS = "\t"
  OFS = "\t"
}

NR == FNR {
  if (FNR == 1 && $1 == "checklist_id") {
    next
  }

  checklist_to_city_file[$1] = $2
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

  header = "checklist_id"
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

  next
}

{
  checklist_id = $id_col

  if (checklist_id in checklist_to_city_file) {
    outfile = outdir "/" checklist_to_city_file[checklist_id] ".tsv"

    if (!(outfile in wrote_header)) {
      print header > outfile
      wrote_header[outfile] = 1
      close(outfile)
    }

    locality_id_val = ""
    locality_val = ""
    tax_val = ""
    category_val = ""
    lat_val = ""
    lon_val = ""
    protocol_val = ""
    approved_val = ""
    reviewed_val = ""

    if (locality_id_col > 0) {
      locality_id_val = $locality_id_col
    }

    if (locality_col > 0) {
      locality_val = $locality_col
    }

    if (tax_col > 0) {
      tax_val = $tax_col
    }

    if (category_col > 0) {
      category_val = $category_col
    }

    if (lat_col > 0) {
      lat_val = $lat_col
    }

    if (lon_col > 0) {
      lon_val = $lon_col
    }

    if (protocol_col > 0) {
      protocol_val = $protocol_col
    }

    if (approved_col > 0) {
      approved_val = $approved_col
    }

    if (reviewed_col > 0) {
      reviewed_val = $reviewed_col
    }

    row = checklist_id
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

    print row >> outfile
    close(outfile)
  }
}
'

writeLines(awk_script, awk_script_path)

filter_cmd <- paste(
  "env -u MallocStackLogging -u MallocStackLoggingNoCompact",
  "/usr/bin/tar -O -x -f",
  shQuote(tar_path),
  shQuote(ebd_member),
  "| gzip -cd",
  paste0("| env -u MallocStackLogging -u MallocStackLoggingNoCompact awk -v outdir=", shQuote(city_raw_out_dir)),
  "-f",
  shQuote(awk_script_path),
  shQuote(checklist_city_map_out),
  "-",
  "2>",
  shQuote(awk_log_path)
)

# Wrap the external command in bash and redirect console chatter.
# The AWK-specific stderr log is still written to awk_log_path above.
shell_cmd <- paste(
  "bash -lc",
  shQuote(paste0(filter_cmd, " > ", shQuote(awk_quiet_log_path), " 2>&1"))
)

cat("\nFiltering and splitting full EBD observation file by city.\n")
cat("This may take a long time because it streams the full EBD file.\n")
cat("\nExternal command console messages will be written to:\n")
cat(awk_quiet_log_path, "\n")
cat("\nAWK diagnostic log will be written to:\n")
cat(awk_log_path, "\n")

exit_code <- system(shell_cmd)

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
} else {
  cat("No AWK log file found.\n")
}

cat("\nExternal command console log preview:\n")
if (file.exists(awk_quiet_log_path)) {
  quiet_log <- readLines(awk_quiet_log_path, warn = FALSE)
  if (length(quiet_log) > 0) {
    print(tail(quiet_log, 25))
  } else {
    cat("No external console messages.\n")
  }
} else {
  cat("No external console log file found.\n")
}

if (!identical(exit_code, 0L)) {
  stop("AWK filtering did not complete successfully. Check the logs above.")
}
}

# ------------------------------------------------------------
# 6. Hard checkpoint after AWK split
# ------------------------------------------------------------

city_raw_files <- list.files(city_raw_out_dir, full.names = TRUE, pattern = "\\.tsv$")

cat("\nCHECKPOINT: Raw city files available:\n")
print(length(city_raw_files))

if (length(city_raw_files) == 0) {
  stop("Zero raw city files were found. Do not continue. Check checklist IDs, EBD member, and AWK log.")
}

raw_test <- fread(city_raw_files[1], nrows = 5)
cat("\nRaw city file columns from first file:\n")
print(names(raw_test))

if (!all(c("locality_id", "locality") %in% names(raw_test))) {
  stop("Raw city files do not contain locality_id/locality. Do not continue.")
}

city_file_sizes <- data.table(
  raw_observation_file = city_raw_files,
  raw_file_size_mb = file.info(city_raw_files)$size / 1024^2
)[order(-raw_file_size_mb)]

cat("\nLargest raw city files:\n")
print(head(city_file_sizes, 20))

city_manifest_dt <- as.data.table(city_manifest)

city_manifest_dt[, raw_file_exists := file.exists(raw_observation_file)]

city_manifest_dt[
  ,
  raw_file_size_mb := fifelse(
    raw_file_exists,
    file.info(raw_observation_file)$size / 1024^2,
    as.numeric(NA)
  )
]

fwrite(city_manifest_dt, city_manifest_out)

cat("\nUpdated city manifest written to:\n")
cat(city_manifest_out, "\n")

# ------------------------------------------------------------
# 7. Optional metadata-expanded observation files
# ------------------------------------------------------------

if (write_metadata_observation_files) {
  
  cat("\nCreating metadata-expanded city observation files...\n")
  cat("This can be very slow and can create very large files.\n")
  
  for (i in seq_len(nrow(city_manifest_dt))) {
    
    this_stub <- city_manifest_dt$city_file_stub[i]
    
    raw_file <- file.path(city_raw_out_dir, paste0(this_stub, ".tsv"))
    metadata_file <- file.path(city_metadata_out_dir, paste0(this_stub, ".tsv"))
    
    if (!file.exists(raw_file)) {
      next
    }
    
    obs <- fread(raw_file)
    
    meta <- checklist_metadata_small[city_file_stub == this_stub][, city_file_stub := NULL]
    
    obs_with_meta <- merge(obs, meta, by = "checklist_id", all.x = TRUE, sort = FALSE)
    
    fwrite(obs_with_meta, metadata_file, sep = "\t")
    
    if (i %% 25 == 0) {
      cat("Metadata-expanded files processed:", i, "of", nrow(city_manifest_dt), "\n")
    }
    
    rm(obs, meta, obs_with_meta)
    gc()
  }
}

if (!write_metadata_observation_files) {
  cat("\nSkipping metadata-expanded observation files because write_metadata_observation_files = FALSE.\n")
}

# ------------------------------------------------------------
# 8. Process each city and write one summary file per city
# ------------------------------------------------------------

cat("\nProcessing city summaries incrementally.\n")
cat("Each completed city writes its own files and can be skipped on rerun.\n")

checklist_metadata_small <- as.data.table(checklist_metadata_small)
setkey(checklist_metadata_small, checklist_id)

status_rows <- vector("list", nrow(city_manifest_dt))

for (i in seq_len(nrow(city_manifest_dt))) {
  
  this_stub <- city_manifest_dt$city_file_stub[i]
  
  raw_file <- city_manifest_dt$raw_observation_file[i]
  checklist_city_out <- city_manifest_dt$checklist_richness_file[i]
  city_year_city_out <- city_manifest_dt$city_year_unique_file[i]
  
  raw_exists <- file.exists(raw_file)
  
  checklist_done <- file.exists(checklist_city_out) &&
    file.info(checklist_city_out)$size > 0
  
  city_year_done <- file.exists(city_year_city_out) &&
    file.info(city_year_city_out)$size > 0
  
  if (!raw_exists) {
    
    cat("\nSkipping city", i, "of", nrow(city_manifest_dt), "because no raw file exists:", this_stub, "\n")
    
    status_rows[[i]] <- data.table(
      city_file_stub = this_stub,
      raw_file_exists = FALSE,
      checklist_richness_done = FALSE,
      city_year_unique_done = FALSE,
      status = "missing_raw_file"
    )
    
    next
  }
  
  if (checklist_done && city_year_done) {
    
    cat("\nSkipping completed city", i, "of", nrow(city_manifest_dt), ":", this_stub, "\n")
    
    status_rows[[i]] <- data.table(
      city_file_stub = this_stub,
      raw_file_exists = TRUE,
      checklist_richness_done = TRUE,
      city_year_unique_done = TRUE,
      status = "already_complete"
    )
    
    next
  }
  
  cat("\nProcessing city", i, "of", nrow(city_manifest_dt), ":", this_stub, "\n")
  
  obs <- fread(raw_file)
  
  required_obs_cols <- c("checklist_id", "observation_date", "locality_id", "locality", "scientific_name", "category")
  missing_obs_cols <- setdiff(required_obs_cols, names(obs))
  
  if (length(missing_obs_cols) > 0) {
    stop(
      "Raw city observation file is missing columns: ",
      paste(missing_obs_cols, collapse = ", "),
      "\nFile: ",
      raw_file
    )
  }
  
  obs[, year := as.integer(substr(observation_date, 1, 4))]
  
  obs_species <- obs[
    category == "species" &
      !is.na(scientific_name) &
      scientific_name != ""
  ]
  
  if (nrow(obs_species) == 0) {
    
    cat("No species-level records in this file. Writing empty status and skipping.\n")
    
    status_rows[[i]] <- data.table(
      city_file_stub = this_stub,
      raw_file_exists = TRUE,
      checklist_richness_done = FALSE,
      city_year_unique_done = FALSE,
      status = "no_species_records"
    )
    
    rm(obs, obs_species)
    gc()
    next
  }
  
  checklist_counts <- obs_species[
    ,
    .(
      species_richness = as.numeric(uniqueN(scientific_name)),
      total_observation_records = as.numeric(.N)
    ),
    by = .(
      checklist_id
    )
  ]
  
  checklist_richness_i <- merge(
    checklist_counts,
    checklist_metadata_small,
    by = "checklist_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  checklist_richness_i[, year := as.integer(substr(observation_date, 1, 4))]
  
  checklist_richness_i <- checklist_richness_i[
    ,
    .(
      checklist_id,
      locality_id,
      locality,
      observation_date,
      year,
      city_id,
      NAME10,
      city_ascii,
      state_name,
      county_name,
      dryland_percent,
      duration_minutes,
      effort_distance_km,
      number_observers,
      latitude_checklist,
      longitude_checklist,
      species_richness,
      total_observation_records,
      city_file_stub
    )
  ]
  
  fwrite(checklist_richness_i, checklist_city_out)
  
  city_year_unique_i <- obs_species[
    ,
    .(
      true_city_year_species_observed = as.numeric(uniqueN(scientific_name)),
      total_species_observation_records = as.numeric(.N)
    ),
    by = .(
      year
    )
  ]
  
  city_info <- city_manifest_dt[
    city_file_stub == this_stub,
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
  
  if (nrow(city_info) != 1) {
    stop("Could not uniquely match city manifest row for: ", this_stub)
  }
  
  city_year_unique_i[, city_id := city_info$city_id]
  city_year_unique_i[, NAME10 := city_info$NAME10]
  city_year_unique_i[, city_ascii := city_info$city_ascii]
  city_year_unique_i[, state_name := city_info$state_name]
  city_year_unique_i[, county_name := city_info$county_name]
  city_year_unique_i[, dryland_percent := city_info$dryland_percent]
  city_year_unique_i[, city_file_stub := city_info$city_file_stub]
  
  setcolorder(
    city_year_unique_i,
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
  
  fwrite(city_year_unique_i, city_year_city_out)
  
  cat("Checklist rows written for this city:", nrow(checklist_richness_i), "\n")
  cat("City-year rows written for this city:", nrow(city_year_unique_i), "\n")
  
  status_rows[[i]] <- data.table(
    city_file_stub = this_stub,
    raw_file_exists = TRUE,
    checklist_richness_done = TRUE,
    city_year_unique_done = TRUE,
    status = "processed_now"
  )
  
  # Save status after every city so progress is visible and recoverable.
  status_now <- rbindlist(status_rows, fill = TRUE)
  fwrite(status_now, processing_status_out)
  
  rm(obs, obs_species, checklist_counts, checklist_richness_i, city_year_unique_i, city_info)
  gc()
}

processing_status <- rbindlist(status_rows, fill = TRUE)
fwrite(processing_status, processing_status_out)

cat("\nCity processing status written to:\n")
cat(processing_status_out, "\n")

cat("\nStatus counts:\n")
print(processing_status[, .N, by = status])

# ------------------------------------------------------------
# 9. Combine per-city checklist richness outputs
# ------------------------------------------------------------

checklist_city_files <- list.files(checklist_richness_by_city_dir, full.names = TRUE, pattern = "\\.csv$")

cat("\nCombining per-city checklist richness files:\n")
print(length(checklist_city_files))

if (length(checklist_city_files) == 0) {
  stop("No per-city checklist richness files found.")
}

checklist_richness_list <- lapply(checklist_city_files, fread)
checklist_richness <- rbindlist(checklist_richness_list, fill = TRUE)

fwrite(checklist_richness, checklist_richness_out)

cat("\nCombined checklist richness rows written:\n")
print(nrow(checklist_richness))

cat("\nCombined checklist richness written to:\n")
cat(checklist_richness_out, "\n")

# ------------------------------------------------------------
# 9.1 Build hotspot/location list for GEE
# ------------------------------------------------------------

hotspot_list <- checklist_richness[
  !is.na(locality_id) & locality_id != "",
  .(
    n_checklists = as.integer(uniqueN(checklist_id)),
    first_date = min(observation_date, na.rm = TRUE),
    last_date = max(observation_date, na.rm = TRUE),
    mean_latitude = mean(latitude_checklist, na.rm = TRUE),
    mean_longitude = mean(longitude_checklist, na.rm = TRUE)
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

setorder(hotspot_list, state_name, city_ascii, locality)
fwrite(hotspot_list, hotspot_list_out)

hotspot_latlong <- hotspot_list[
  ,
  .(
    latitude = mean(mean_latitude, na.rm = TRUE),
    longitude = mean(mean_longitude, na.rm = TRUE),
    n_checklists = sum(n_checklists, na.rm = TRUE)
  ),
  by = .(locality_id, locality)
]

setorder(hotspot_latlong, locality_id)
fwrite(hotspot_latlong, hotspot_latlong_out)

cat("\nHotspot list written to:\n")
cat(hotspot_list_out, "\n")

cat("\nSimple hotspot lat/long file written to:\n")
cat(hotspot_latlong_out, "\n")

cat("\nRows with locality_id in combined checklist richness:\n")
print(sum(!is.na(checklist_richness$locality_id) & checklist_richness$locality_id != ""))

cat("\nUnique locality IDs in combined checklist richness:\n")
print(uniqueN(checklist_richness$locality_id[!is.na(checklist_richness$locality_id) & checklist_richness$locality_id != ""]))

rm(checklist_richness_list)
gc()

# ------------------------------------------------------------
# 10. Combine per-city city-year unique species outputs
# ------------------------------------------------------------

city_year_city_files <- list.files(city_year_unique_by_city_dir, full.names = TRUE, pattern = "\\.csv$")

cat("\nCombining per-city city-year unique species files:\n")
print(length(city_year_city_files))

if (length(city_year_city_files) == 0) {
  stop("No per-city city-year unique species files found.")
}

city_year_unique_list <- lapply(city_year_city_files, fread)
city_year_unique_species <- rbindlist(city_year_unique_list, fill = TRUE)

fwrite(city_year_unique_species, city_year_unique_species_out)

cat("\nCombined city-year unique species rows written:\n")
print(nrow(city_year_unique_species))

cat("\nCombined city-year unique species written to:\n")
cat(city_year_unique_species_out, "\n")

rm(city_year_unique_list)
gc()

# ------------------------------------------------------------
# 11. Build final city-year richness summary
# ------------------------------------------------------------

cat("\nBuilding final city-year richness summary...\n")

checklist_richness[, species_richness := as.numeric(species_richness)]
checklist_richness[, total_observation_records := as.numeric(total_observation_records)]
checklist_richness[, duration_minutes := as.numeric(duration_minutes)]
checklist_richness[, effort_distance_km := as.numeric(effort_distance_km)]
checklist_richness[, number_observers := as.numeric(number_observers)]

city_year_effort_summary <- checklist_richness[
  ,
  .(
    n_checklists = as.integer(uniqueN(checklist_id)),
    n_localities = as.integer(uniqueN(locality_id[!is.na(locality_id) & locality_id != ""])),
    mean_checklist_richness = as.numeric(mean(species_richness, na.rm = TRUE)),
    median_checklist_richness = as.numeric(median(species_richness, na.rm = TRUE)),
    summed_checklist_richness = as.numeric(sum(species_richness, na.rm = TRUE)),
    total_observation_records = as.numeric(sum(total_observation_records, na.rm = TRUE)),
    mean_duration_minutes = as.numeric(mean(duration_minutes, na.rm = TRUE)),
    mean_effort_distance_km = as.numeric(mean(effort_distance_km, na.rm = TRUE)),
    mean_number_observers = as.numeric(mean(number_observers, na.rm = TRUE))
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

setorder(city_year_richness, state_name, city_ascii, year)

fwrite(city_year_richness, city_year_richness_out)

cat("\nFinal city-year richness rows written:\n")
print(nrow(city_year_richness))

cat("\nFinal city-year richness written to:\n")
cat(city_year_richness_out, "\n")

# ------------------------------------------------------------
# 12. Final summary
# ------------------------------------------------------------

cat("\nScript 05 complete.\n")
cat("\nMain outputs:\n")
cat("City raw observation folder:\n")
cat(city_raw_out_dir, "\n")
cat("Per-city checklist richness folder:\n")
cat(checklist_richness_by_city_dir, "\n")
cat("Per-city city-year unique species folder:\n")
cat(city_year_unique_by_city_dir, "\n")
cat("City file manifest:\n")
cat(city_manifest_out, "\n")
cat("Checklist metadata crosswalk:\n")
cat(checklist_metadata_out, "\n")
cat("Processing status:\n")
cat(processing_status_out, "\n")
cat("Combined checklist richness summary:\n")
cat(checklist_richness_out, "\n")
cat("Combined city-year unique species intermediate:\n")
cat(city_year_unique_species_out, "\n")
cat("Final city-year richness summary:\n")
cat(city_year_richness_out, "\n")
cat("Hotspot list for GEE:\n")
cat(hotspot_list_out, "\n")
cat("Simple hotspot lat/long file for GEE:\n")
cat(hotspot_latlong_out, "\n")

if (write_metadata_observation_files) {
  cat("City metadata observation folder:\n")
  cat(city_metadata_out_dir, "\n")
} else {
  cat("\nMetadata-expanded observation files were skipped by design.\n")
  cat("Use raw city files plus checklist metadata crosswalk for downstream joins.\n")
}
