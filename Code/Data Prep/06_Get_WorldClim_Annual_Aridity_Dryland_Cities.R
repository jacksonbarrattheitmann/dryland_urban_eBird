# ============================================================
# 06_Get_WorldClim_Annual_Aridity_Dryland_Cities.R
#
# Purpose:
# Download WorldClim historical monthly weather data, calculate
# annual aridity for dryland cities, and join aridity to the
# eBird city-year, checklist-level, and observation-level outputs.
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
# WorldClim historical monthly weather is currently available
# through 2024, so this script calculates aridity for 2015-2024.
# eBird rows for 2025 will remain unmatched unless another climate
# product is added later.
# ============================================================

library(sf)
library(dplyr)
library(readr)
library(data.table)
library(terra)

# ------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------

dryland_cities_path <- "Data/Processed/Shape Files/Cities/dryland_cities_50pct.shp"

# eBird outputs from File 05
city_year_richness_path <- "Data/Processed/eBird/dryland_city_year_species_richness_2015_2025.csv"

checklist_richness_path <- "Data/Processed/eBird/dryland_city_checklist_species_richness_2015_2025.csv"

observations_with_city_path <- "Data/Processed/eBird/dryland_city_observations_with_city_2015_2025.tsv"

# Raw WorldClim folder:
# downloaded ZIPs and extracted rasters will both be stored here.
raw_worldclim_dir <- "Data/Raw/WorldClim"

# Intermediate processed climate outputs
processed_climate_dir <- "Data/Processed/Climate"

# Intermediate processed eBird outputs
processed_ebird_dir <- "Data/Processed/eBird"

# Final analysis-ready outputs
for_analysis_dir <- "Data/Processed/For Analysis"

dir.create(raw_worldclim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_climate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_ebird_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(for_analysis_dir, recursive = TRUE, showWarnings = FALSE)

# Climate-only intermediate outputs
city_monthly_climate_out <- file.path(
  processed_climate_dir,
  "dryland_city_monthly_worldclim_climate_2015_2024.csv"
)

city_year_aridity_out <- file.path(
  processed_climate_dir,
  "dryland_city_year_worldclim_aridity_2015_2024.csv"
)

# Final analysis-ready outputs
city_year_richness_aridity_out <- file.path(
  for_analysis_dir,
  "dryland_city_year_species_richness_aridity_2015_2025.csv"
)

checklist_richness_aridity_out <- file.path(
  for_analysis_dir,
  "dryland_city_checklist_species_richness_aridity_2015_2025.csv"
)

observations_with_city_aridity_out <- file.path(
  for_analysis_dir,
  "dryland_city_observations_with_city_aridity_2015_2025.tsv"
)

# Optional copy of the aridity-only table in For Analysis
city_year_aridity_for_analysis_out <- file.path(
  for_analysis_dir,
  "dryland_city_year_worldclim_aridity_2015_2024.csv"
)

# ------------------------------------------------------------
# 2. Read dryland city polygons
# ------------------------------------------------------------

# The city output is currently a shapefile, so several fields were
# shortened. We rename them back to readable names here.
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
# This gives one climate value per city-year.
# We do this in an equal-area CRS first to avoid lon/lat warnings,
# then transform back to WGS84 for raster extraction.
city_points <- dryland_cities %>%
  st_transform(5070) %>%
  st_point_on_surface() %>%
  st_transform(4326)

city_points_vect <- terra::vect(city_points)

# ------------------------------------------------------------
# 3. WorldClim settings
# ------------------------------------------------------------

# WorldClim historical monthly weather files are grouped by period.
# For 2015-2024 we need:
# - 2010-2019
# - 2020-2024
#
# Variables:
# - tmin = monthly minimum temperature in Celsius
# - tmax = monthly maximum temperature in Celsius
# - prec = monthly total precipitation in mm

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
# 4. Download raw WorldClim files
# ------------------------------------------------------------

# Raw downloads are saved to:
# Data/Raw/WorldClim/
#
# Each ZIP is also extracted into a matching folder under:
# Data/Raw/WorldClim/

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
    
    # Remove failed partial downloads, such as 0-byte or tiny 404 files.
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
# 5. Helper function to find monthly raster files
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
  
  # WorldClim historical monthly files usually include year and month.
  # Use flexible matching to protect against minor filename differences.
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
# 6. Extract monthly climate values at city points
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
        
        # Monthly average temperature in Celsius
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
# 7. Calculate annual aridity by city-year
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
    # Total annual precipitation in mm
    total_annual_precip_mm = sum(monthly_precip_mm, na.rm = TRUE),
    
    # Mean annual temperature in Celsius
    mean_annual_temp_c = mean(monthly_mean_temp_c, na.rm = TRUE),
    
    n_months_with_climate = sum(
      !is.na(monthly_precip_mm) &
        !is.na(monthly_mean_temp_c)
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    # De Martonne Aridity Index:
    # AI = total annual precipitation / (mean annual temp + 10)
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

# Also save a copy of the climate-only aridity table in For Analysis.
write_csv(
  city_year_aridity,
  city_year_aridity_for_analysis_out
)

cat("\nCity-year aridity output written:\n")
print(city_year_aridity_out)

cat("\nCity-year aridity copy written to For Analysis:\n")
print(city_year_aridity_for_analysis_out)

# ------------------------------------------------------------
# 8. Join aridity to city-year eBird richness
# ------------------------------------------------------------

city_year_richness <- fread(city_year_richness_path) %>%
  as_tibble()

city_year_richness_aridity <- city_year_richness %>%
  left_join(
    city_year_aridity %>%
      select(
        city_id,
        year,
        total_annual_precip_mm,
        mean_annual_temp_c,
        aridity_index,
        aridity_class_de_martonne,
        n_months_with_climate
      ),
    by = c(
      "city_id" = "city_id",
      "observation_year" = "year"
    )
  )

write_csv(
  city_year_richness_aridity,
  city_year_richness_aridity_out
)

cat("\nCity-year richness + aridity output written:\n")
print(city_year_richness_aridity_out)

# ------------------------------------------------------------
# 9. Join aridity to checklist-level richness
# ------------------------------------------------------------

checklist_richness <- fread(checklist_richness_path) %>%
  as_tibble()

checklist_richness_aridity <- checklist_richness %>%
  left_join(
    city_year_aridity %>%
      select(
        city_id,
        year,
        total_annual_precip_mm,
        mean_annual_temp_c,
        aridity_index,
        aridity_class_de_martonne,
        n_months_with_climate
      ),
    by = c(
      "city_id" = "city_id",
      "observation_year" = "year"
    )
  )

write_csv(
  checklist_richness_aridity,
  checklist_richness_aridity_out
)

cat("\nChecklist richness + aridity output written:\n")
print(checklist_richness_aridity_out)

# ------------------------------------------------------------
# 10. Join aridity to full observation-level eBird file
# ------------------------------------------------------------

# This joins annual city-level aridity to the full observation-level
# eBird dataset produced in File 05.
#
# This file has ~27 million rows, so we avoid unnecessary copies.
# fread() already returns a data.table, so do NOT call as.data.table()
# on the full observations object.
#
# Aridity is joined by:
#   city_id + observation_year
#
# WorldClim aridity is available for 2015-2024. Observation rows from
# 2025 will remain NA unless another climate product is added later.

cat("\nPreparing memory for observation-level aridity join...\n")

rm(
  city_monthly_climate,
  city_year_richness,
  city_year_richness_aridity,
  checklist_richness,
  checklist_richness_aridity,
  dryland_cities,
  city_points,
  city_points_vect
)

gc()

cat("\nReading full observation-level eBird file for aridity join:\n")
print(observations_with_city_path)

observations_with_city <- fread(
  observations_with_city_path,
  sep = "\t",
  quote = "",
  showProgress = TRUE
)

cat("\nObservation rows loaded:\n")
print(nrow(observations_with_city))

# fread() already returns a data.table.
# Do not run as.data.table(observations_with_city), because that can
# create an unnecessary full copy and hit the memory limit.

if (!"observation_year" %in% names(observations_with_city)) {
  observations_with_city[
    ,
    observation_year := as.integer(substr(observation_date, 1, 4))
  ]
}

aridity_lookup <- city_year_aridity %>%
  dplyr::select(
    city_id,
    year,
    total_annual_precip_mm,
    mean_annual_temp_c,
    aridity_index,
    aridity_class_de_martonne,
    n_months_with_climate
  ) %>%
  as.data.table()

setnames(
  aridity_lookup,
  old = "year",
  new = "observation_year"
)

# Update join by reference: adds aridity columns without creating a
# second 27-million-row merged copy.
observations_with_city[
  aridity_lookup,
  on = .(city_id, observation_year),
  `:=`(
    total_annual_precip_mm = i.total_annual_precip_mm,
    mean_annual_temp_c = i.mean_annual_temp_c,
    aridity_index = i.aridity_index,
    aridity_class_de_martonne = i.aridity_class_de_martonne,
    n_months_with_climate = i.n_months_with_climate
  )
]

cat("\nObservation rows after aridity join:\n")
print(nrow(observations_with_city))

cat("\nRows with aridity matched:\n")
print(sum(!is.na(observations_with_city$aridity_index)))

cat("\nRows without aridity matched, expected mostly/entirely for 2025:\n")
print(sum(is.na(observations_with_city$aridity_index)))

cat("\nWriting full observation-level eBird + aridity file:\n")
print(observations_with_city_aridity_out)

fwrite(
  observations_with_city,
  observations_with_city_aridity_out,
  sep = "\t"
)

cat("\nFull observation-level eBird + aridity output written:\n")
print(observations_with_city_aridity_out)

# ------------------------------------------------------------
# 11. Final checks because this is a lot all together lmao
# ------------------------------------------------------------

# This section is written so it can be run even after restarting R.
# If the objects are not currently in memory, it reloads the final
# files from disk.

# ------------------------------------------------------------
# Reload objects if needed
# ------------------------------------------------------------

if (!exists("city_year_aridity")) {
  city_year_aridity <- fread(city_year_aridity_for_analysis_out) %>%
    as_tibble()
}

if (!exists("city_year_richness_aridity")) {
  city_year_richness_aridity <- fread(city_year_richness_aridity_out) %>%
    as_tibble()
}

if (!exists("checklist_richness_aridity")) {
  checklist_richness_aridity <- fread(checklist_richness_aridity_out) %>%
    as_tibble()
}

# For the huge observation-level file, only read it if it is not already
# in memory. This may take a little while.
if (!exists("observations_with_city")) {
  observations_with_city <- fread(
    observations_with_city_aridity_out,
    sep = "\t",
    quote = "",
    showProgress = TRUE
  )
}

# ------------------------------------------------------------
# Aridity-year checks
# ------------------------------------------------------------

cat("\nAridity years available:\n")
print(sort(unique(city_year_aridity$year)))

# ------------------------------------------------------------
# City-year richness + aridity checks
# ------------------------------------------------------------

cat("\nRows in city-year richness + aridity:\n")
print(nrow(city_year_richness_aridity))

cat("\nRows with aridity matched in city-year richness:\n")
print(sum(!is.na(city_year_richness_aridity$aridity_index)))

cat("\nRows without aridity matched in city-year richness, expected for 2025 if present:\n")
print(sum(is.na(city_year_richness_aridity$aridity_index)))

# ------------------------------------------------------------
# Checklist richness + aridity checks
# ------------------------------------------------------------

cat("\nRows in checklist richness + aridity:\n")
print(nrow(checklist_richness_aridity))

cat("\nRows with aridity matched in checklist richness:\n")
print(sum(!is.na(checklist_richness_aridity$aridity_index)))

cat("\nRows without aridity matched in checklist richness, expected for 2025 if present:\n")
print(sum(is.na(checklist_richness_aridity$aridity_index)))

# ------------------------------------------------------------
# Observation-level eBird + aridity checks
# ------------------------------------------------------------

cat("\nRows in observation-level eBird + aridity:\n")
print(nrow(observations_with_city))

cat("\nRows with aridity matched in observation-level eBird:\n")
print(sum(!is.na(observations_with_city$aridity_index)))

cat("\nRows without aridity matched in observation-level eBird, expected for 2025 if present:\n")
print(sum(is.na(observations_with_city$aridity_index)))

# ------------------------------------------------------------
# Final output paths
# ------------------------------------------------------------

cat("\nFinal analysis-ready outputs:\n")
print(city_year_aridity_for_analysis_out)
print(city_year_richness_aridity_out)
print(checklist_richness_aridity_out)
print(observations_with_city_aridity_out)

# ------------------------------------------------------------
# Preview city-year output
# ------------------------------------------------------------

cat("\nPreview of city-year richness + aridity:\n")
print(
  city_year_richness_aridity %>%
    arrange(desc(city_year_species_richness)) %>%
    select(
      NAME10,
      city_ascii,
      state_name,
      observation_year,
      city_year_species_richness,
      n_checklists,
      total_annual_precip_mm,
      mean_annual_temp_c,
      aridity_index,
      aridity_class_de_martonne
    ) %>%
    head(25)
)

# ------------------------------------------------------------
# Preview observation-level output
# ------------------------------------------------------------

cat("\nPreview of observation-level eBird + aridity:\n")

obs_preview_cols <- c(
  "checklist_id",
  "observation_date",
  "observation_year",
  "city_id",
  "NAME10",
  "common_name",
  "scientific_name",
  "total_annual_precip_mm",
  "mean_annual_temp_c",
  "aridity_index",
  "aridity_class_de_martonne"
)

missing_preview_cols <- setdiff(obs_preview_cols, names(observations_with_city))

if (length(missing_preview_cols) > 0) {
  cat("\nWarning: These preview columns are missing from observations_with_city:\n")
  print(missing_preview_cols)
  cat("\nAvailable columns are:\n")
  print(names(observations_with_city))
} else {
  print(
    observations_with_city[
      ,
      ..obs_preview_cols
    ][1:10]
  )
}

# ------------------------------------------------------------
# File existence checks
# ------------------------------------------------------------

cat("\nChecking that final files exist:\n")

final_files <- c(
  city_year_aridity_for_analysis_out,
  city_year_richness_aridity_out,
  checklist_richness_aridity_out,
  observations_with_city_aridity_out
)

print(
  data.frame(
    file = final_files,
    exists = file.exists(final_files),
    size_mb = round(file.info(final_files)$size / 1024^2, 2)
  )
)