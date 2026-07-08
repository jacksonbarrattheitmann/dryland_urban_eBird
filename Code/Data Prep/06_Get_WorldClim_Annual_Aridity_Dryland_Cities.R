# ============================================================
# 06_Get_WorldClim_Annual_Aridity_Dryland_Cities_WITH_HOTSPOTS.R
#
# Purpose:
# Download/extract WorldClim historical monthly weather data,
# calculate annual De Martonne aridity for dryland cities, and join
# aridity to the new incremental eBird outputs from Script 05.
#
# This version matches the revised Script 05 structure:
#   1. It joins aridity to the combined city-year richness file.
#   2. It joins aridity to the combined checklist-level richness file.
#   3. It joins aridity to each per-city checklist richness file.
#   4. It joins aridity to each per-city city-year richness file.
#   5. It copies the hotspot/locality files created by Script 05
#      into the For Analysis folder.
#   6. It does NOT require the old giant observation-level file:
#      dryland_city_observations_with_city_2015_2025.tsv
#
# Aridity metric:
# De Martonne Aridity Index
#
#   AI = total annual precipitation / (mean annual temperature + 10)
#
# where:
#   total annual precipitation = sum of monthly precipitation in mm
#   mean annual temperature    = mean of monthly average temperature in Celsius
#   monthly average temp       = (monthly tmin + monthly tmax) / 2
#
# Notes:
# WorldClim historical monthly weather is currently available through 2024,
# so this script calculates aridity for 2015-2024. eBird rows for 2025
# will remain unmatched/NA unless another climate product is added later.
# ============================================================

library(sf)
library(dplyr)
library(readr)
library(data.table)
library(terra)

# ------------------------------------------------------------
# 0. Run settings
# ------------------------------------------------------------

# TRUE = redownload/recalculate WorldClim-derived city-year aridity.
# FALSE = reuse existing aridity table if it already exists.
recalculate_aridity <- FALSE

# TRUE = rebuild aridity-enhanced per-city files.
# FALSE = skip per-city files that already exist.
#
# Use TRUE for this first run after adding locality_id/locality,
# because older per-city aridity files may not include the new fields.
overwrite_per_city_aridity_outputs <- TRUE

# ------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------

dryland_cities_path <- "Data/Processed/Shape Files/Cities/dryland_cities_50pct.shp"

# eBird outputs from revised incremental Script 05
city_year_richness_path <- "Data/Processed/eBird/dryland_city_year_species_richness_2015_2025.csv"

checklist_richness_path <- "Data/Processed/eBird/dryland_city_checklist_species_richness_2015_2025.csv"

city_manifest_path <- "Data/Processed/eBird/dryland_city_observation_file_manifest_2015_2025.csv"

hotspot_list_path <- "Data/Processed/eBird/dryland_city_hotspot_list_for_GEE_2015_2025.csv"

hotspot_latlong_path <- "Data/Processed/eBird/dryland_city_hotspot_latlong_for_GEE_2015_2025.csv"

checklist_richness_by_city_dir <- "Data/Processed/eBird/Checklist_Richness_By_City"

city_year_unique_by_city_dir <- "Data/Processed/eBird/City_Year_Unique_Species_By_City"

# Raw WorldClim folder
raw_worldclim_dir <- "Data/Raw/WorldClim"

# Processed climate outputs
processed_climate_dir <- "Data/Processed/Climate"

# Final analysis-ready outputs
for_analysis_dir <- "Data/Processed/For Analysis"

# Per-city aridity-enhanced output folders
checklist_richness_aridity_by_city_dir <- file.path(
  for_analysis_dir,
  "Checklist_Richness_Aridity_By_City"
)

city_year_aridity_by_city_dir <- file.path(
  for_analysis_dir,
  "City_Year_Richness_Aridity_By_City"
)

dir.create(raw_worldclim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_climate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(for_analysis_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checklist_richness_aridity_by_city_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(city_year_aridity_by_city_dir, recursive = TRUE, showWarnings = FALSE)

# Climate-only intermediate outputs
city_monthly_climate_out <- file.path(
  processed_climate_dir,
  "dryland_city_monthly_worldclim_climate_2015_2024.csv"
)

city_year_aridity_out <- file.path(
  processed_climate_dir,
  "dryland_city_year_worldclim_aridity_2015_2024.csv"
)

# Copies / final analysis-ready outputs
city_year_without_aridity_for_analysis_out <- file.path(
  for_analysis_dir,
  "dryland_city_year_species_richness_2015_2025.csv"
)

checklist_without_aridity_for_analysis_out <- file.path(
  for_analysis_dir,
  "dryland_city_checklist_species_richness_2015_2025.csv"
)

hotspot_list_for_analysis_out <- file.path(
  for_analysis_dir,
  "dryland_city_hotspot_list_for_GEE_2015_2025.csv"
)

hotspot_latlong_for_analysis_out <- file.path(
  for_analysis_dir,
  "dryland_city_hotspot_latlong_for_GEE_2015_2025.csv"
)

city_year_richness_aridity_out <- file.path(
  for_analysis_dir,
  "dryland_city_year_species_richness_aridity_2015_2025.csv"
)

checklist_richness_aridity_out <- file.path(
  for_analysis_dir,
  "dryland_city_checklist_species_richness_aridity_2015_2025.csv"
)

city_year_aridity_for_analysis_out <- file.path(
  for_analysis_dir,
  "dryland_city_year_worldclim_aridity_2015_2024.csv"
)

per_city_aridity_status_out <- file.path(
  for_analysis_dir,
  "dryland_city_per_city_aridity_join_status_2015_2025.csv"
)

# ------------------------------------------------------------
# 2. Confirm required eBird outputs exist
# ------------------------------------------------------------

required_files <- c(
  city_year_richness_path,
  checklist_richness_path,
  city_manifest_path,
  hotspot_list_path,
  hotspot_latlong_path
)

missing_required_files <- required_files[!file.exists(required_files)]

if (length(missing_required_files) > 0) {
  stop(
    "Missing required eBird output files from Script 05:\n",
    paste(missing_required_files, collapse = "\n")
  )
}

if (!dir.exists(checklist_richness_by_city_dir)) {
  stop("Missing per-city checklist richness folder: ", checklist_richness_by_city_dir)
}

if (!dir.exists(city_year_unique_by_city_dir)) {
  stop("Missing per-city city-year folder: ", city_year_unique_by_city_dir)
}

# Confirm the revised Script 05 outputs include hotspot/locality fields.
checklist_header <- fread(checklist_richness_path, nrows = 0)

required_checklist_hotspot_cols <- c(
  "checklist_id",
  "locality_id",
  "locality",
  "latitude_checklist",
  "longitude_checklist"
)

missing_checklist_hotspot_cols <- setdiff(
  required_checklist_hotspot_cols,
  names(checklist_header)
)

if (length(missing_checklist_hotspot_cols) > 0) {
  stop(
    "The checklist richness file is missing hotspot/locality columns. ",
    "Rerun revised Script 05 first. Missing columns: ",
    paste(missing_checklist_hotspot_cols, collapse = ", ")
  )
}

hotspot_latlong_header <- fread(hotspot_latlong_path, nrows = 0)

required_hotspot_latlong_cols <- c(
  "locality_id",
  "locality",
  "latitude",
  "longitude",
  "n_checklists"
)

missing_hotspot_latlong_cols <- setdiff(
  required_hotspot_latlong_cols,
  names(hotspot_latlong_header)
)

if (length(missing_hotspot_latlong_cols) > 0) {
  stop(
    "The simple hotspot lat/long file is missing required columns: ",
    paste(missing_hotspot_latlong_cols, collapse = ", ")
  )
}

cat("\nHotspot/locality outputs from Script 05 confirmed.\n")

# ------------------------------------------------------------
# 3. Build or reuse city-year aridity
# ------------------------------------------------------------

if (file.exists(city_year_aridity_out) && !recalculate_aridity) {
  
  cat("\nReusing existing city-year aridity table:\n")
  cat(city_year_aridity_out, "\n")
  
  city_year_aridity <- fread(city_year_aridity_out) %>%
    as_tibble()
  
} else {
  
  # ------------------------------------------------------------
  # 3.1 Read dryland city polygons
  # ------------------------------------------------------------
  
  cat("\nReading dryland city polygons:\n")
  cat(dryland_cities_path, "\n")
  
  dryland_cities <- st_read(dryland_cities_path) %>%
    st_make_valid() %>%
    st_transform(4326) %>%
    dplyr::select(
      city_id,
      NAME10,
      city_ascii = city_sc,
      state_name = stat_nm,
      county_name = cnty_nm,
      geometry
    )
  
  # Use points on surface for climate extraction.
  city_points <- dryland_cities %>%
    st_transform(5070) %>%
    st_point_on_surface() %>%
    st_transform(4326)
  
  city_points_vect <- terra::vect(city_points)
  
  # ------------------------------------------------------------
  # 3.2 WorldClim settings
  # ------------------------------------------------------------
  
  resolution <- "10m"
  
  periods <- c(
    "2010-2019",
    "2020-2024"
  )
  
  variables <- c(
    "tmin",
    "tmax",
    "prec"
  )
  
  years <- 2015:2024
  months <- 1:12
  
  make_worldclim_url <- function(var, period, resolution = "10m") {
    paste0(
      "https://geodata.ucdavis.edu/climate/worldclim/2_1/hist/cts4.09/",
      "wc2.1_cruts4.09_",
      resolution,
      "_",
      var,
      "_",
      period,
      ".zip"
    )
  }
  
  # ------------------------------------------------------------
  # 3.3 Download raw WorldClim files
  # ------------------------------------------------------------
  
  for (period in periods) {
    for (var in variables) {
      
      url <- make_worldclim_url(
        var = var,
        period = period,
        resolution = resolution
      )
      
      zip_path <- file.path(
        raw_worldclim_dir,
        basename(url)
      )
      
      unzip_dir <- file.path(
        raw_worldclim_dir,
        tools::file_path_sans_ext(basename(url))
      )
      
      if (file.exists(zip_path)) {
        zip_size <- file.info(zip_path)$size
        
        if (is.na(zip_size) || zip_size < 100000) {
          cat("\nRemoving likely failed WorldClim ZIP:\n")
          cat(zip_path, "\n")
          file.remove(zip_path)
        }
      }
      
      if (!file.exists(zip_path)) {
        cat("\nDownloading WorldClim raw data:\n")
        cat(url, "\n")
        cat("Saving to:\n")
        cat(zip_path, "\n")
        
        download.file(
          url = url,
          destfile = zip_path,
          mode = "wb"
        )
      } else {
        cat("\nWorldClim ZIP already exists:\n")
        cat(zip_path, "\n")
      }
      
      if (!dir.exists(unzip_dir)) {
        cat("Extracting WorldClim rasters to:\n")
        cat(unzip_dir, "\n")
        
        dir.create(unzip_dir, recursive = TRUE, showWarnings = FALSE)
        unzip(zip_path, exdir = unzip_dir)
      } else {
        cat("WorldClim rasters already extracted:\n")
        cat(unzip_dir, "\n")
      }
    }
  }
  
  # ------------------------------------------------------------
  # 3.4 Helper function to find monthly raster files
  # ------------------------------------------------------------
  
  find_monthly_file <- function(var, year, month, resolution = "10m") {
    
    period <- ifelse(year <= 2019, "2010-2019", "2020-2024")
    
    unzip_dir <- file.path(
      raw_worldclim_dir,
      paste0("wc2.1_cruts4.09_", resolution, "_", var, "_", period)
    )
    
    files <- list.files(
      unzip_dir,
      pattern = "\\.tif$",
      recursive = TRUE,
      full.names = TRUE
    )
    
    if (length(files) == 0) {
      stop(
        "No .tif files found in ",
        unzip_dir,
        ". Check whether the WorldClim ZIP was downloaded and extracted correctly."
      )
    }
    
    year_month_patterns <- c(
      paste0(year, "-", sprintf("%02d", month)),
      paste0(year, "_", sprintf("%02d", month)),
      paste0(year, ".", sprintf("%02d", month)),
      paste0(year, sprintf("%02d", month))
    )
    
    match_files <- files[
      Reduce(
        `|`,
        lapply(
          year_month_patterns,
          function(p) grepl(p, basename(files))
        )
      )
    ]
    
    if (length(match_files) == 0) {
      cat("\nAvailable files in folder:\n")
      print(head(basename(files), 30))
      
      stop(
        "Could not find raster for ",
        var,
        " ",
        year,
        "-",
        sprintf("%02d", month),
        " in ",
        unzip_dir
      )
    }
    
    match_files[1]
  }
  
  # ------------------------------------------------------------
  # 3.5 Extract monthly climate values at city points
  # ------------------------------------------------------------
  
  city_monthly_climate <- list()
  counter <- 1
  
  for (year in years) {
    
    cat("\nProcessing WorldClim year:", year, "\n")
    
    for (month in months) {
      
      cat("  Month:", sprintf("%02d", month), "\n")
      
      tmin_file <- find_monthly_file("tmin", year, month, resolution)
      tmax_file <- find_monthly_file("tmax", year, month, resolution)
      prec_file <- find_monthly_file("prec", year, month, resolution)
      
      tmin_r <- terra::rast(tmin_file)
      tmax_r <- terra::rast(tmax_file)
      prec_r <- terra::rast(prec_file)
      
      tmin_vals <- terra::extract(tmin_r, city_points_vect)
      tmax_vals <- terra::extract(tmax_r, city_points_vect)
      prec_vals <- terra::extract(prec_r, city_points_vect)
      
      month_df <- city_points %>%
        st_drop_geometry() %>%
        mutate(
          year = year,
          month = month,
          tmin_c = as.numeric(tmin_vals[[2]]),
          tmax_c = as.numeric(tmax_vals[[2]]),
          monthly_precip_mm = as.numeric(prec_vals[[2]]),
          monthly_mean_temp_c = (tmin_c + tmax_c) / 2
        )
      
      city_monthly_climate[[counter]] <- month_df
      counter <- counter + 1
    }
  }
  
  city_monthly_climate <- bind_rows(city_monthly_climate)
  
  write_csv(
    city_monthly_climate,
    city_monthly_climate_out
  )
  
  # ------------------------------------------------------------
  # 3.6 Calculate annual aridity by city-year
  # ------------------------------------------------------------
  
  city_year_aridity <- city_monthly_climate %>%
    group_by(
      city_id,
      NAME10,
      city_ascii,
      state_name,
      county_name,
      year
    ) %>%
    summarise(
      total_annual_precip_mm = sum(monthly_precip_mm, na.rm = TRUE),
      mean_annual_temp_c = mean(monthly_mean_temp_c, na.rm = TRUE),
      n_months_with_climate = sum(
        !is.na(monthly_precip_mm) &
          !is.na(monthly_mean_temp_c)
      ),
      .groups = "drop"
    ) %>%
    mutate(
      aridity_index = total_annual_precip_mm / (mean_annual_temp_c + 10),
      aridity_class_de_martonne = case_when(
        aridity_index < 5 ~ "extremely arid",
        aridity_index < 10 ~ "arid",
        aridity_index < 20 ~ "semiarid",
        aridity_index < 30 ~ "dry subhumid",
        aridity_index < 60 ~ "humid",
        TRUE ~ "very humid"
      )
    )
  
  write_csv(
    city_year_aridity,
    city_year_aridity_out
  )
}

# Save a copy of the climate-only aridity table in For Analysis.
write_csv(
  city_year_aridity,
  city_year_aridity_for_analysis_out
)

cat("\nCity-year aridity available:\n")
cat(city_year_aridity_out, "\n")

cat("\nCity-year aridity copy written to For Analysis:\n")
cat(city_year_aridity_for_analysis_out, "\n")

# ------------------------------------------------------------
# 4. Prepare aridity lookup table
# ------------------------------------------------------------

aridity_lookup <- city_year_aridity %>%
  dplyr::select(
    city_id,
    year,
    total_annual_precip_mm,
    mean_annual_temp_c,
    aridity_index,
    aridity_class_de_martonne,
    n_months_with_climate
  )

# ------------------------------------------------------------
# 5. Copy non-aridity combined outputs to For Analysis
# ------------------------------------------------------------

file.copy(
  from = city_year_richness_path,
  to = city_year_without_aridity_for_analysis_out,
  overwrite = TRUE
)

file.copy(
  from = checklist_richness_path,
  to = checklist_without_aridity_for_analysis_out,
  overwrite = TRUE
)

file.copy(
  from = hotspot_list_path,
  to = hotspot_list_for_analysis_out,
  overwrite = TRUE
)

file.copy(
  from = hotspot_latlong_path,
  to = hotspot_latlong_for_analysis_out,
  overwrite = TRUE
)

cat("\nCopied non-aridity combined eBird files to For Analysis:\n")
cat(city_year_without_aridity_for_analysis_out, "\n")
cat(checklist_without_aridity_for_analysis_out, "\n")

cat("\nCopied hotspot/locality files to For Analysis:\n")
cat(hotspot_list_for_analysis_out, "\n")
cat(hotspot_latlong_for_analysis_out, "\n")

# ------------------------------------------------------------
# 6. Join aridity to combined city-year eBird richness
# ------------------------------------------------------------

city_year_richness <- fread(city_year_richness_path) %>%
  as_tibble()

if (!"year" %in% names(city_year_richness)) {
  stop(
    "The combined city-year richness file does not have a year column. ",
    "Available columns are: ",
    paste(names(city_year_richness), collapse = ", ")
  )
}

city_year_richness_aridity <- city_year_richness %>%
  left_join(
    aridity_lookup,
    by = c(
      "city_id" = "city_id",
      "year" = "year"
    )
  )

write_csv(
  city_year_richness_aridity,
  city_year_richness_aridity_out
)

cat("\nCity-year richness + aridity output written:\n")
cat(city_year_richness_aridity_out, "\n")

# ------------------------------------------------------------
# 7. Join aridity to combined checklist-level richness
# ------------------------------------------------------------

checklist_richness <- fread(checklist_richness_path) %>%
  as_tibble()

if (!"year" %in% names(checklist_richness)) {
  stop(
    "The combined checklist richness file does not have a year column. ",
    "Available columns are: ",
    paste(names(checklist_richness), collapse = ", ")
  )
}

checklist_richness_aridity <- checklist_richness %>%
  left_join(
    aridity_lookup,
    by = c(
      "city_id" = "city_id",
      "year" = "year"
    )
  )

write_csv(
  checklist_richness_aridity,
  checklist_richness_aridity_out
)

cat("\nChecklist richness + aridity output written:\n")
cat(checklist_richness_aridity_out, "\n")

cat("\nLocality columns retained in checklist richness + aridity:\n")
print(
  grep(
    "locality|hotspot|location",
    names(checklist_richness_aridity),
    ignore.case = TRUE,
    value = TRUE
  )
)

cat("\nRows with locality_id in checklist richness + aridity:\n")
print(sum(!is.na(checklist_richness_aridity$locality_id) & checklist_richness_aridity$locality_id != ""))

# ------------------------------------------------------------
# 8. Join aridity to each per-city checklist richness file
# ------------------------------------------------------------

cat("\nJoining aridity to per-city checklist richness files...\n")

checklist_city_files <- list.files(
  checklist_richness_by_city_dir,
  full.names = TRUE,
  pattern = "\\.csv$"
)

city_year_city_files <- list.files(
  city_year_unique_by_city_dir,
  full.names = TRUE,
  pattern = "\\.csv$"
)

if (length(checklist_city_files) == 0) {
  stop("No per-city checklist richness files found in: ", checklist_richness_by_city_dir)
}

if (length(city_year_city_files) == 0) {
  stop("No per-city city-year files found in: ", city_year_unique_by_city_dir)
}

status_rows <- list()
status_counter <- 1

for (i in seq_along(checklist_city_files)) {
  
  in_file <- checklist_city_files[i]
  out_file <- file.path(
    checklist_richness_aridity_by_city_dir,
    basename(in_file)
  )
  
  if (
    file.exists(out_file) &&
    file.info(out_file)$size > 0 &&
    !overwrite_per_city_aridity_outputs
  ) {
    
    cat(
      "\nSkipping existing per-city checklist aridity file",
      i,
      "of",
      length(checklist_city_files),
      ":",
      basename(out_file),
      "\n"
    )
    
    status_rows[[status_counter]] <- data.table(
      input_file = basename(in_file),
      output_file = basename(out_file),
      output_type = "checklist_richness_aridity_by_city",
      status = "already_complete"
    )
    
    status_counter <- status_counter + 1
    next
  }
  
  cat(
    "\nProcessing per-city checklist file",
    i,
    "of",
    length(checklist_city_files),
    ":",
    basename(in_file),
    "\n"
  )
  
  x <- fread(in_file) %>%
    as_tibble()
  
  if (!"year" %in% names(x)) {
    stop("Missing year column in per-city checklist file: ", in_file)
  }
  
  x_aridity <- x %>%
    left_join(
      aridity_lookup,
      by = c(
        "city_id" = "city_id",
        "year" = "year"
      )
    )
  
  write_csv(x_aridity, out_file)
  
  status_rows[[status_counter]] <- data.table(
    input_file = basename(in_file),
    output_file = basename(out_file),
    output_type = "checklist_richness_aridity_by_city",
    status = "processed_now"
  )
  
  status_counter <- status_counter + 1
}

# ------------------------------------------------------------
# 9. Join aridity to each per-city city-year file
# ------------------------------------------------------------

cat("\nJoining aridity to per-city city-year files...\n")

for (i in seq_along(city_year_city_files)) {
  
  in_file <- city_year_city_files[i]
  out_file <- file.path(
    city_year_aridity_by_city_dir,
    basename(in_file)
  )
  
  if (
    file.exists(out_file) &&
    file.info(out_file)$size > 0 &&
    !overwrite_per_city_aridity_outputs
  ) {
    
    cat(
      "\nSkipping existing per-city city-year aridity file",
      i,
      "of",
      length(city_year_city_files),
      ":",
      basename(out_file),
      "\n"
    )
    
    status_rows[[status_counter]] <- data.table(
      input_file = basename(in_file),
      output_file = basename(out_file),
      output_type = "city_year_aridity_by_city",
      status = "already_complete"
    )
    
    status_counter <- status_counter + 1
    next
  }
  
  cat(
    "\nProcessing per-city city-year file",
    i,
    "of",
    length(city_year_city_files),
    ":",
    basename(in_file),
    "\n"
  )
  
  x <- fread(in_file) %>%
    as_tibble()
  
  if (!"year" %in% names(x)) {
    stop("Missing year column in per-city city-year file: ", in_file)
  }
  
  x_aridity <- x %>%
    left_join(
      aridity_lookup,
      by = c(
        "city_id" = "city_id",
        "year" = "year"
      )
    )
  
  write_csv(x_aridity, out_file)
  
  status_rows[[status_counter]] <- data.table(
    input_file = basename(in_file),
    output_file = basename(out_file),
    output_type = "city_year_aridity_by_city",
    status = "processed_now"
  )
  
  status_counter <- status_counter + 1
}

per_city_aridity_status <- rbindlist(status_rows, fill = TRUE)

fwrite(
  per_city_aridity_status,
  per_city_aridity_status_out
)

cat("\nPer-city aridity join status written:\n")
cat(per_city_aridity_status_out, "\n")

cat("\nPer-city aridity join status counts:\n")
print(per_city_aridity_status[, .N, by = .(output_type, status)])

# ------------------------------------------------------------
# 10. Final checks
# ------------------------------------------------------------

cat("\nAridity years available:\n")
print(sort(unique(city_year_aridity$year)))

cat("\nRows in city-year richness + aridity:\n")
print(nrow(city_year_richness_aridity))

cat("\nRows with aridity matched in city-year richness:\n")
print(sum(!is.na(city_year_richness_aridity$aridity_index)))

cat("\nRows without aridity matched in city-year richness, expected for 2025 if present:\n")
print(sum(is.na(city_year_richness_aridity$aridity_index)))

cat("\nRows in checklist richness + aridity:\n")
print(nrow(checklist_richness_aridity))

cat("\nRows with aridity matched in checklist richness:\n")
print(sum(!is.na(checklist_richness_aridity$aridity_index)))

cat("\nRows without aridity matched in checklist richness, expected for 2025 if present:\n")
print(sum(is.na(checklist_richness_aridity$aridity_index)))

cat("\nPer-city checklist aridity files written:\n")
print(length(list.files(
  checklist_richness_aridity_by_city_dir,
  pattern = "\\.csv$"
)))

cat("\nPer-city city-year aridity files written:\n")
print(length(list.files(
  city_year_aridity_by_city_dir,
  pattern = "\\.csv$"
)))

cat("\nChecking final files exist:\n")

final_files <- c(
  city_year_aridity_for_analysis_out,
  city_year_without_aridity_for_analysis_out,
  checklist_without_aridity_for_analysis_out,
  hotspot_list_for_analysis_out,
  hotspot_latlong_for_analysis_out,
  city_year_richness_aridity_out,
  checklist_richness_aridity_out,
  per_city_aridity_status_out
)

print(
  data.frame(
    file = final_files,
    exists = file.exists(final_files),
    size_mb = round(file.info(final_files)$size / 1024^2, 2)
  )
)

cat("\nPreview of hotspot lat/long file for GEE:\n")

hotspot_latlong_preview <- fread(hotspot_latlong_for_analysis_out)

print(head(hotspot_latlong_preview, 20))

cat("\nRows in hotspot lat/long file for GEE:\n")
print(nrow(hotspot_latlong_preview))

cat("\nUnique locality IDs in hotspot lat/long file for GEE:\n")
print(uniqueN(hotspot_latlong_preview$locality_id))

cat("\nPreview of city-year richness + aridity:\n")

preview_cols <- c(
  "NAME10",
  "city_ascii",
  "state_name",
  "year",
  "true_city_year_species_observed",
  "n_checklists",
  "total_annual_precip_mm",
  "mean_annual_temp_c",
  "aridity_index",
  "aridity_class_de_martonne"
)

available_preview_cols <- intersect(preview_cols, names(city_year_richness_aridity))

print(
  city_year_richness_aridity %>%
    arrange(desc(true_city_year_species_observed)) %>%
    select(all_of(available_preview_cols)) %>%
    head(25)
)

cat("\nScript 06 complete.\n")
cat("\nMain aridity-enhanced outputs:\n")
cat(city_year_richness_aridity_out, "\n")
cat(checklist_richness_aridity_out, "\n")
cat(checklist_richness_aridity_by_city_dir, "\n")
cat(city_year_aridity_by_city_dir, "\n")

cat("\nHotspot/locality outputs copied to For Analysis:\n")
cat(hotspot_list_for_analysis_out, "\n")
cat(hotspot_latlong_for_analysis_out, "\n")
