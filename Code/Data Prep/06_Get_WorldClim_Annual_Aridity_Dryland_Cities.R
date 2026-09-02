# ============================================================
# 06_Get_WorldClim_Annual_Aridity_Dryland_Cities_WITH_HOTSPOTS.R
#
# Purpose:
# Calculate or reuse annual WorldClim-derived De Martonne aridity
# for the 653 dryland cities and join those climate values to the
# QC-approved Script 05 eBird outputs.
#
# IMPORTANT:
#
# - Script 05 uses the July 2026 eBird release.
# - Script 05 has independently passed QC.
# - City geometry comes from the Script 02 shapefile.
# - Authoritative city metadata come from the Script 02 CSV.
# - WorldClim historical monthly climate is available here for
#   2015-2024.
# - Therefore 2025 eBird rows intentionally retain NA climate.
#
# De Martonne Aridity Index:
#
#   AI = annual precipitation / (mean annual temperature + 10)
#
# Higher AI = wetter / less arid conditions.
#
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(sf)
library(data.table)
library(terra)


# ------------------------------------------------------------
# 2. Run settings
# ------------------------------------------------------------

# FALSE:
# Reuse the existing 2015-2024 city-year WorldClim table if it
# exists and passes QC.
#
# TRUE:
# Recalculate climate from the raw WorldClim rasters.

recalculate_aridity <- FALSE


# TRUE for this first run after rebuilding Script 05.
#
# This clears/rebuilds the per-city climate-enhanced outputs so
# no older April-derived products can remain.

overwrite_per_city_aridity_outputs <- TRUE


# Expected data structure.

expected_master_city_count <- 653L

expected_ebird_city_count <- 647L

climate_years <- 2015:2024

ebird_years <- 2015:2025

months <- 1:12


# ------------------------------------------------------------
# 3. File paths
# ------------------------------------------------------------

# Authoritative Script 02 city metadata.

dryland_city_metadata_path <-
  "Data/Processed/Tables/dryland_cities_50pct.csv"


# Script 02 city geometry.

dryland_cities_path <-
  "Data/Processed/Shape Files/Cities/dryland_cities_50pct.shp"


# Script 05 completion / release information.

script05_release_marker_path <- file.path(
  "Data/Processed/eBird",
  "script05_ebird_release.txt"
)

script05_filtered_marker_path <- file.path(
  "Data/Processed/eBird",
  "script05_filtered_ebd_complete.txt"
)

script05_summary_marker_path <- file.path(
  "Data/Processed/eBird",
  "script05_ebd_summaries_complete.txt"
)


# Script 05 combined products.

city_year_richness_path <- file.path(
  "Data/Processed/eBird",
  "dryland_city_year_species_richness_2015_2025.csv"
)

checklist_richness_path <- file.path(
  "Data/Processed/eBird",
  "dryland_city_checklist_species_richness_2015_2025.csv"
)

city_manifest_path <- file.path(
  "Data/Processed/eBird",
  "dryland_city_observation_file_manifest_2015_2025.csv"
)

hotspot_list_path <- file.path(
  "Data/Processed/eBird",
  "dryland_city_hotspot_list_for_GEE_2015_2025.csv"
)

hotspot_latlong_path <- file.path(
  "Data/Processed/eBird",
  "dryland_city_hotspot_latlong_for_GEE_2015_2025.csv"
)


# Script 05 per-city folders.

checklist_richness_by_city_dir <- file.path(
  "Data/Processed/eBird",
  "Checklist_Richness_By_City"
)

city_year_unique_by_city_dir <- file.path(
  "Data/Processed/eBird",
  "City_Year_Unique_Species_By_City"
)


# WorldClim folders.

raw_worldclim_dir <-
  "Data/Raw/WorldClim"

processed_climate_dir <-
  "Data/Processed/Climate"


# Final analysis-ready folder.

for_analysis_dir <-
  "Data/Processed/For Analysis"


# Per-city aridity-enhanced output folders.

checklist_richness_aridity_by_city_dir <- file.path(
  for_analysis_dir,
  "Checklist_Richness_Aridity_By_City"
)

city_year_aridity_by_city_dir <- file.path(
  for_analysis_dir,
  "City_Year_Richness_Aridity_By_City"
)


dir.create(
  raw_worldclim_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  processed_climate_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  for_analysis_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  checklist_richness_aridity_by_city_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  city_year_aridity_by_city_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 4. Output paths
# ------------------------------------------------------------

city_monthly_climate_out <- file.path(
  processed_climate_dir,
  "dryland_city_monthly_worldclim_climate_2015_2024.csv"
)

city_year_aridity_out <- file.path(
  processed_climate_dir,
  "dryland_city_year_worldclim_aridity_2015_2024.csv"
)


# Copies of non-aridity Script 05 products.

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


# Climate-enhanced products.

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
# 5. Small atomic-write helper
# ------------------------------------------------------------

atomic_fwrite <- function(x, path) {
  
  tmp <- paste0(
    path,
    ".tmp"
  )
  
  
  if (file.exists(tmp)) {
    
    unlink(
      tmp,
      force = TRUE
    )
  }
  
  
  fwrite(
    x,
    tmp
  )
  
  
  if (
    !file.exists(tmp) ||
    file.info(tmp)$size == 0
  ) {
    
    stop(
      "Temporary output was not written:\n",
      tmp
    )
  }
  
  
  if (file.exists(path)) {
    
    unlink(
      path,
      force = TRUE
    )
  }
  
  
  if (
    !file.rename(
      tmp,
      path
    )
  ) {
    
    stop(
      "Could not finalize output:\n",
      path
    )
  }
  
  
  invisible(path)
}


# ------------------------------------------------------------
# 6. Confirm required upstream products
# ------------------------------------------------------------

required_files <- c(
  dryland_city_metadata_path,
  city_year_richness_path,
  checklist_richness_path,
  city_manifest_path,
  hotspot_list_path,
  hotspot_latlong_path,
  script05_release_marker_path,
  script05_filtered_marker_path,
  script05_summary_marker_path
)


missing_required_files <- required_files[
  !file.exists(required_files)
]


if (
  length(missing_required_files) > 0
) {
  
  stop(
    "Missing required upstream files:\n",
    paste(
      missing_required_files,
      collapse = "\n"
    )
  )
}


if (
  !dir.exists(
    checklist_richness_by_city_dir
  )
) {
  
  stop(
    "Missing Script 05 per-city checklist folder:\n",
    checklist_richness_by_city_dir
  )
}


if (
  !dir.exists(
    city_year_unique_by_city_dir
  )
) {
  
  stop(
    "Missing Script 05 per-city city-year folder:\n",
    city_year_unique_by_city_dir
  )
}


# ------------------------------------------------------------
# 7. Confirm Script 05 belongs to July 2026
# ------------------------------------------------------------

script05_release <- readLines(
  script05_release_marker_path,
  warn = FALSE
)


filtered_marker <- readLines(
  script05_filtered_marker_path,
  warn = FALSE
)


summary_marker <- readLines(
  script05_summary_marker_path,
  warn = FALSE
)


script05_release_ok <-
  length(script05_release) > 0 &&
  trimws(script05_release[1]) == "relJul-2026"


filtered_marker_ok <- any(
  filtered_marker == "release=relJul-2026"
)


summary_marker_ok <- any(
  summary_marker == "release=relJul-2026"
)


cat("\n============================================================\n")
cat("SCRIPT 05 RELEASE QC\n")
cat("============================================================\n\n")


cat("Script 05 release:\n")
print(script05_release)


cat("\nFiltered EBD marker = July 2026:\n")
print(filtered_marker_ok)


cat("\nCompact-summary marker = July 2026:\n")
print(summary_marker_ok)


if (!script05_release_ok) {
  
  stop(
    "Script 05 release marker is not relJul-2026."
  )
}


if (!filtered_marker_ok) {
  
  stop(
    "Script 05 filtered-EBD marker is not relJul-2026."
  )
}


if (!summary_marker_ok) {
  
  stop(
    "Script 05 summary marker is not relJul-2026."
  )
}


# ------------------------------------------------------------
# 8. Read authoritative 653-city metadata
# ------------------------------------------------------------

city_master <- fread(
  dryland_city_metadata_path,
  showProgress = FALSE
)


required_city_cols <- c(
  "city_id",
  "NAME10",
  "city_ascii",
  "state_name",
  "county_name",
  "dryland_percent"
)


missing_city_cols <- setdiff(
  required_city_cols,
  names(city_master)
)


if (
  length(missing_city_cols) > 0
) {
  
  stop(
    "Authoritative city table is missing required columns:\n",
    paste(
      missing_city_cols,
      collapse = "\n"
    )
  )
}


city_master[
  ,
  city_id := as.integer(city_id)
]


# Clean only tiny floating-point artifacts around 100%.

city_master[
  dryland_percent > 100 &
    dryland_percent <= 100 + 1e-8,
  dryland_percent := 100
]


cat("\n============================================================\n")
cat("MASTER CITY QC\n")
cat("============================================================\n\n")


cat("Master city rows:\n")
print(nrow(city_master))


cat("\nUnique city IDs:\n")
print(uniqueN(city_master$city_id))


cat("\nDryland-percent range:\n")

print(
  range(
    city_master$dryland_percent,
    na.rm = TRUE
  )
)


if (
  nrow(city_master) !=
  expected_master_city_count
) {
  
  stop(
    "Expected 653 master cities."
  )
}


if (
  uniqueN(city_master$city_id) !=
  expected_master_city_count
) {
  
  stop(
    "city_id is not unique in master city table."
  )
}


if (
  any(
    city_master$dryland_percent < 50 - 1e-8 |
    city_master$dryland_percent > 100 + 1e-8,
    na.rm = TRUE
  )
) {
  
  stop(
    "Invalid dryland_percent values detected."
  )
}


# ------------------------------------------------------------
# 9. Read and QC Script 05 combined products
# ------------------------------------------------------------

checklist_richness <- fread(
  checklist_richness_path,
  showProgress = TRUE
)


city_year_richness <- fread(
  city_year_richness_path,
  showProgress = FALSE
)


city_manifest <- fread(
  city_manifest_path,
  showProgress = FALSE
)


checklist_richness[
  ,
  checklist_id := as.character(checklist_id)
]


checklist_richness[
  ,
  city_id := as.integer(city_id)
]


checklist_richness[
  ,
  year := as.integer(year)
]


city_year_richness[
  ,
  city_id := as.integer(city_id)
]


city_year_richness[
  ,
  year := as.integer(year)
]


city_manifest[
  ,
  city_id := as.integer(city_id)
]


required_checklist_cols <- c(
  "checklist_id",
  "locality_id",
  "locality",
  "year",
  "city_id",
  "species_richness",
  "latitude_checklist",
  "longitude_checklist"
)


missing_checklist_cols <- setdiff(
  required_checklist_cols,
  names(checklist_richness)
)


if (
  length(missing_checklist_cols) > 0
) {
  
  stop(
    "Script 05 checklist richness file is missing columns:\n",
    paste(
      missing_checklist_cols,
      collapse = "\n"
    )
  )
}


if (
  !"year" %in%
  names(city_year_richness)
) {
  
  stop(
    "Script 05 city-year richness file is missing year."
  )
}


input_checklist_rows <- nrow(
  checklist_richness
)


input_unique_checklists <- uniqueN(
  checklist_richness$checklist_id
)


input_ebird_city_count <- uniqueN(
  checklist_richness$city_id
)


input_locality_count <- uniqueN(
  checklist_richness$locality_id[
    !is.na(checklist_richness$locality_id) &
      checklist_richness$locality_id != ""
  ]
)


input_city_year_rows <- nrow(
  city_year_richness
)


cat("\n============================================================\n")
cat("SCRIPT 05 INPUT QC\n")
cat("============================================================\n\n")


cat("Checklist rows:\n")
print(input_checklist_rows)


cat("\nUnique checklist IDs:\n")
print(input_unique_checklists)


cat("\nCities:\n")
print(input_ebird_city_count)


cat("\nUnique locality IDs:\n")
print(input_locality_count)


cat("\nCity-year rows:\n")
print(input_city_year_rows)


cat("\nCity manifest rows:\n")
print(nrow(city_manifest))


if (
  input_checklist_rows !=
  input_unique_checklists
) {
  
  stop(
    "Duplicate checklist IDs detected in Script 05 input."
  )
}


if (
  input_ebird_city_count !=
  expected_ebird_city_count
) {
  
  stop(
    "Expected 647 eBird cities from Script 05."
  )
}


if (
  nrow(city_manifest) !=
  expected_ebird_city_count
) {
  
  stop(
    "Expected 647 rows in Script 05 city manifest."
  )
}


unknown_ebird_city_ids <- setdiff(
  unique(
    checklist_richness$city_id
  ),
  city_master$city_id
)


if (
  length(
    unknown_ebird_city_ids
  ) > 0
) {
  
  stop(
    "Script 05 contains city IDs not found in master city table."
  )
}


if (
  any(
    !checklist_richness$year %in%
    ebird_years
  )
) {
  
  stop(
    "Unexpected years detected in checklist richness."
  )
}


if (
  any(
    !city_year_richness$year %in%
    ebird_years
  )
) {
  
  stop(
    "Unexpected years detected in city-year richness."
  )
}


# ------------------------------------------------------------
# 10. Build or reuse WorldClim city-year aridity
# ------------------------------------------------------------

if (
  file.exists(city_year_aridity_out) &&
  !recalculate_aridity
) {
  
  cat("\n============================================================\n")
  cat("REUSING EXISTING WORLDCLIM ARIDITY\n")
  cat("============================================================\n\n")
  
  
  cat("Existing climate table:\n")
  cat(city_year_aridity_out, "\n")
  
  
  city_year_aridity_raw <- fread(
    city_year_aridity_out,
    showProgress = FALSE
  )
  
  
  required_existing_climate_cols <- c(
    "city_id",
    "year",
    "total_annual_precip_mm",
    "mean_annual_temp_c",
    "n_months_with_climate",
    "aridity_index",
    "aridity_class_de_martonne"
  )
  
  
  missing_existing_climate_cols <- setdiff(
    required_existing_climate_cols,
    names(city_year_aridity_raw)
  )
  
  
  if (
    length(
      missing_existing_climate_cols
    ) > 0
  ) {
    
    stop(
      "Existing aridity table is missing required columns:\n",
      paste(
        missing_existing_climate_cols,
        collapse = "\n"
      ),
      "\nSet recalculate_aridity <- TRUE after resolving this."
    )
  }
  
  
  city_year_aridity_raw[
    ,
    city_id := as.integer(city_id)
  ]
  
  
  city_year_aridity_raw[
    ,
    year := as.integer(year)
  ]
  
  
  # Retain the climate values only.
  #
  # City names/state/county are reattached below from the
  # authoritative Script 02 city table.
  
  climate_numeric <- unique(
    city_year_aridity_raw[
      ,
      .(
        city_id,
        year,
        total_annual_precip_mm,
        mean_annual_temp_c,
        n_months_with_climate,
        aridity_index,
        aridity_class_de_martonne
      )
    ]
  )
  
  
  rm(
    city_year_aridity_raw
  )
  
  
} else {
  
  cat("\n============================================================\n")
  cat("RECALCULATING WORLDCLIM ARIDITY\n")
  cat("============================================================\n\n")
  
  
  # ----------------------------------------------------------
  # 10.1 Read geometry only
  # ----------------------------------------------------------
  
  if (
    !file.exists(
      dryland_cities_path
    )
  ) {
    
    stop(
      "Dryland city shapefile not found:\n",
      dryland_cities_path
    )
  }
  
  
  dryland_city_geometry <- st_read(
    dryland_cities_path,
    quiet = TRUE
  )
  
  
  if (
    !"city_id" %in%
    names(dryland_city_geometry)
  ) {
    
    stop(
      "Dryland city shapefile is missing city_id."
    )
  }
  
  
  dryland_city_geometry <- dryland_city_geometry[
    "city_id"
  ]
  
  
  dryland_city_geometry$city_id <- as.integer(
    dryland_city_geometry$city_id
  )
  
  
  dryland_city_geometry <- st_make_valid(
    dryland_city_geometry
  )
  
  
  dryland_city_geometry <- st_transform(
    dryland_city_geometry,
    4326
  )
  
  
  geometry_ids_missing_from_master <- setdiff(
    dryland_city_geometry$city_id,
    city_master$city_id
  )
  
  
  master_ids_missing_from_geometry <- setdiff(
    city_master$city_id,
    dryland_city_geometry$city_id
  )
  
  
  if (
    length(
      geometry_ids_missing_from_master
    ) > 0 ||
    
    length(
      master_ids_missing_from_geometry
    ) > 0
  ) {
    
    stop(
      "City IDs differ between Script 02 geometry and metadata."
    )
  }
  
  
  if (
    nrow(
      dryland_city_geometry
    ) !=
    expected_master_city_count
  ) {
    
    stop(
      "Expected 653 dryland city geometries."
    )
  }
  
  
  # ----------------------------------------------------------
  # 10.2 City climate extraction points
  # ----------------------------------------------------------
  
  city_points <- st_transform(
    dryland_city_geometry,
    5070
  )
  
  
  city_points <- st_point_on_surface(
    city_points
  )
  
  
  city_points <- st_transform(
    city_points,
    4326
  )
  
  
  city_points_vect <- terra::vect(
    city_points
  )
  
  
  # ----------------------------------------------------------
  # 10.3 WorldClim settings
  # ----------------------------------------------------------
  
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
  
  
  make_worldclim_url <- function(
    var,
    period,
    resolution = "10m"
  ) {
    
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
  
  
  # ----------------------------------------------------------
  # 10.4 Download / extract raw WorldClim files
  # ----------------------------------------------------------
  
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
        tools::file_path_sans_ext(
          basename(url)
        )
      )
      
      
      if (file.exists(zip_path)) {
        
        zip_size <- file.info(
          zip_path
        )$size
        
        
        if (
          is.na(zip_size) ||
          zip_size < 100000
        ) {
          
          cat(
            "\nRemoving likely failed WorldClim ZIP:\n"
          )
          
          cat(
            zip_path,
            "\n"
          )
          
          file.remove(
            zip_path
          )
        }
      }
      
      
      if (!file.exists(zip_path)) {
        
        cat(
          "\nDownloading WorldClim data:\n"
        )
        
        cat(
          url,
          "\n"
        )
        
        
        download.file(
          url = url,
          destfile = zip_path,
          mode = "wb"
        )
        
      } else {
        
        cat(
          "\nWorldClim ZIP already exists:\n"
        )
        
        cat(
          zip_path,
          "\n"
        )
      }
      
      
      if (!dir.exists(unzip_dir)) {
        
        cat(
          "Extracting WorldClim rasters to:\n"
        )
        
        cat(
          unzip_dir,
          "\n"
        )
        
        
        dir.create(
          unzip_dir,
          recursive = TRUE,
          showWarnings = FALSE
        )
        
        
        unzip(
          zip_path,
          exdir = unzip_dir
        )
        
      } else {
        
        cat(
          "WorldClim rasters already extracted:\n"
        )
        
        cat(
          unzip_dir,
          "\n"
        )
      }
    }
  }
  
  
  # ----------------------------------------------------------
  # 10.5 Locate monthly raster files
  # ----------------------------------------------------------
  
  find_monthly_file <- function(
    var,
    year,
    month,
    resolution = "10m"
  ) {
    
    period <- ifelse(
      year <= 2019,
      "2010-2019",
      "2020-2024"
    )
    
    
    unzip_dir <- file.path(
      raw_worldclim_dir,
      paste0(
        "wc2.1_cruts4.09_",
        resolution,
        "_",
        var,
        "_",
        period
      )
    )
    
    
    files <- list.files(
      unzip_dir,
      pattern = "\\.tif$",
      recursive = TRUE,
      full.names = TRUE
    )
    
    
    if (
      length(files) == 0
    ) {
      
      stop(
        "No WorldClim .tif files found in:\n",
        unzip_dir
      )
    }
    
    
    year_month_patterns <- c(
      paste0(
        year,
        "-",
        sprintf("%02d", month)
      ),
      
      paste0(
        year,
        "_",
        sprintf("%02d", month)
      ),
      
      paste0(
        year,
        ".",
        sprintf("%02d", month)
      ),
      
      paste0(
        year,
        sprintf("%02d", month)
      )
    )
    
    
    matched <- rep(
      FALSE,
      length(files)
    )
    
    
    for (p in year_month_patterns) {
      
      matched <- matched |
        grepl(
          p,
          basename(files)
        )
    }
    
    
    match_files <- files[
      matched
    ]
    
    
    if (
      length(match_files) == 0
    ) {
      
      stop(
        "Could not locate WorldClim raster for ",
        var,
        " ",
        year,
        "-",
        sprintf("%02d", month)
      )
    }
    
    
    match_files[1]
  }
  
  
  # ----------------------------------------------------------
  # 10.6 Extract monthly climate
  # ----------------------------------------------------------
  
  city_monthly_list <- vector(
    "list",
    length(climate_years) *
      length(months)
  )
  
  
  counter <- 1L
  
  
  for (year in climate_years) {
    
    cat(
      "\nProcessing WorldClim year:",
      year,
      "\n"
    )
    
    
    for (month in months) {
      
      cat(
        "  Month:",
        sprintf("%02d", month),
        "\n"
      )
      
      
      tmin_file <- find_monthly_file(
        "tmin",
        year,
        month,
        resolution
      )
      
      
      tmax_file <- find_monthly_file(
        "tmax",
        year,
        month,
        resolution
      )
      
      
      prec_file <- find_monthly_file(
        "prec",
        year,
        month,
        resolution
      )
      
      
      tmin_r <- terra::rast(
        tmin_file
      )
      
      
      tmax_r <- terra::rast(
        tmax_file
      )
      
      
      prec_r <- terra::rast(
        prec_file
      )
      
      
      tmin_vals <- terra::extract(
        tmin_r,
        city_points_vect
      )
      
      
      tmax_vals <- terra::extract(
        tmax_r,
        city_points_vect
      )
      
      
      prec_vals <- terra::extract(
        prec_r,
        city_points_vect
      )
      
      
      month_dt <- data.table(
        
        city_id =
          city_points$city_id,
        
        year =
          as.integer(year),
        
        month =
          as.integer(month),
        
        tmin_c =
          as.numeric(
            tmin_vals[[2]]
          ),
        
        tmax_c =
          as.numeric(
            tmax_vals[[2]]
          ),
        
        monthly_precip_mm =
          as.numeric(
            prec_vals[[2]]
          )
      )
      
      
      month_dt[
        ,
        monthly_mean_temp_c :=
          (
            tmin_c +
              tmax_c
          ) / 2
      ]
      
      
      city_monthly_list[[counter]] <-
        month_dt
      
      
      counter <- counter + 1L
      
      
      rm(
        tmin_r,
        tmax_r,
        prec_r,
        tmin_vals,
        tmax_vals,
        prec_vals,
        month_dt
      )
    }
  }
  
  
  city_monthly_numeric <- rbindlist(
    city_monthly_list
  )
  
  
  rm(
    city_monthly_list
  )
  
  
  gc()
  
  
  # ----------------------------------------------------------
  # 10.7 Monthly climate QC
  # ----------------------------------------------------------
  
  expected_monthly_rows <-
    expected_master_city_count *
    length(climate_years) *
    length(months)
  
  
  monthly_unique_keys <- uniqueN(
    paste(
      city_monthly_numeric$city_id,
      city_monthly_numeric$year,
      city_monthly_numeric$month,
      sep = "__"
    )
  )
  
  
  missing_monthly_climate <- sum(
    is.na(
      city_monthly_numeric$tmin_c
    ) |
      is.na(
        city_monthly_numeric$tmax_c
      ) |
      is.na(
        city_monthly_numeric$monthly_precip_mm
      )
  )
  
  
  cat("\n============================================================\n")
  cat("MONTHLY WORLDCLIM QC\n")
  cat("============================================================\n\n")
  
  
  cat("Expected monthly rows:\n")
  print(expected_monthly_rows)
  
  
  cat("\nActual monthly rows:\n")
  print(nrow(city_monthly_numeric))
  
  
  cat("\nUnique city-year-month combinations:\n")
  print(monthly_unique_keys)
  
  
  cat("\nRows with missing climate:\n")
  print(missing_monthly_climate)
  
  
  if (
    nrow(city_monthly_numeric) !=
    expected_monthly_rows
  ) {
    
    stop(
      "Monthly climate row count is incorrect."
    )
  }
  
  
  if (
    monthly_unique_keys !=
    expected_monthly_rows
  ) {
    
    stop(
      "Duplicated or missing city-year-month climate combinations detected."
    )
  }
  
  
  if (
    missing_monthly_climate > 0
  ) {
    
    stop(
      "Missing WorldClim climate values detected."
    )
  }
  
  
  # Attach authoritative metadata before saving monthly table.
  
  city_monthly_climate <- merge(
    
    city_monthly_numeric,
    
    city_master[
      ,
      .(
        city_id,
        NAME10,
        city_ascii,
        state_name,
        county_name
      )
    ],
    
    by = "city_id",
    all.x = TRUE,
    sort = FALSE
  )
  
  
  setcolorder(
    city_monthly_climate,
    c(
      "city_id",
      "NAME10",
      "city_ascii",
      "state_name",
      "county_name",
      "year",
      "month",
      "tmin_c",
      "tmax_c",
      "monthly_mean_temp_c",
      "monthly_precip_mm"
    )
  )
  
  
  setorder(
    city_monthly_climate,
    city_id,
    year,
    month
  )
  
  
  atomic_fwrite(
    city_monthly_climate,
    city_monthly_climate_out
  )
  
  
  # ----------------------------------------------------------
  # 10.8 Annual climate / De Martonne aridity
  # ----------------------------------------------------------
  
  climate_numeric <- city_monthly_numeric[
    ,
    .(
      
      total_annual_precip_mm =
        sum(
          monthly_precip_mm
        ),
      
      mean_annual_temp_c =
        mean(
          monthly_mean_temp_c
        ),
      
      n_months_with_climate =
        .N
    ),
    
    by = .(
      city_id,
      year
    )
  ]
  
  
  climate_numeric[
    ,
    aridity_index :=
      total_annual_precip_mm /
      (
        mean_annual_temp_c +
          10
      )
  ]
  
  
  climate_numeric[
    ,
    aridity_class_de_martonne :=
      fcase(
        
        aridity_index < 5,
        "extremely arid",
        
        aridity_index < 10,
        "arid",
        
        aridity_index < 20,
        "semiarid",
        
        aridity_index < 30,
        "dry subhumid",
        
        aridity_index < 60,
        "humid",
        
        default =
          "very humid"
      )
  ]
  
  
  rm(
    city_monthly_numeric
  )
  
  
  gc()
}


# ------------------------------------------------------------
# 11. Standardize and QC city-year climate table
# ------------------------------------------------------------

climate_numeric[
  ,
  city_id := as.integer(city_id)
]


climate_numeric[
  ,
  year := as.integer(year)
]


expected_climate_rows <-
  expected_master_city_count *
  length(climate_years)


climate_unique_keys <- uniqueN(
  paste(
    climate_numeric$city_id,
    climate_numeric$year,
    sep = "__"
  )
)


unknown_climate_city_ids <- setdiff(
  unique(
    climate_numeric$city_id
  ),
  city_master$city_id
)


master_city_ids_missing_climate <- setdiff(
  city_master$city_id,
  unique(
    climate_numeric$city_id
  )
)


unexpected_climate_years <- setdiff(
  unique(
    climate_numeric$year
  ),
  climate_years
)


missing_expected_climate_years <- setdiff(
  climate_years,
  unique(
    climate_numeric$year
  )
)


invalid_climate_rows <- climate_numeric[
  
  is.na(total_annual_precip_mm) |
    is.na(mean_annual_temp_c) |
    is.na(aridity_index) |
    
    !is.finite(total_annual_precip_mm) |
    !is.finite(mean_annual_temp_c) |
    !is.finite(aridity_index) |
    
    n_months_with_climate != 12
]


formula_check <- climate_numeric[
  ,
  total_annual_precip_mm /
    (
      mean_annual_temp_c +
        10
    )
]


formula_difference <- max(
  abs(
    climate_numeric$aridity_index -
      formula_check
  ),
  na.rm = TRUE
)


cat("\n============================================================\n")
cat("CITY-YEAR CLIMATE QC\n")
cat("============================================================\n\n")


cat("Expected climate rows:\n")
print(expected_climate_rows)


cat("\nActual climate rows:\n")
print(nrow(climate_numeric))


cat("\nUnique city-year keys:\n")
print(climate_unique_keys)


cat("\nClimate cities:\n")
print(uniqueN(climate_numeric$city_id))


cat("\nClimate years:\n")
print(sort(unique(climate_numeric$year)))


cat("\nInvalid climate rows:\n")
print(nrow(invalid_climate_rows))


cat("\nMaximum De Martonne formula difference:\n")
print(formula_difference)


if (
  nrow(climate_numeric) !=
  expected_climate_rows
) {
  
  stop(
    "Expected exactly 6,530 city-year climate rows."
  )
}


if (
  climate_unique_keys !=
  expected_climate_rows
) {
  
  stop(
    "Climate table contains duplicate or missing city-year combinations."
  )
}


if (
  uniqueN(climate_numeric$city_id) !=
  expected_master_city_count
) {
  
  stop(
    "Climate table does not contain all 653 cities."
  )
}


if (
  length(
    unknown_climate_city_ids
  ) > 0
) {
  
  stop(
    "Climate table contains unknown city IDs."
  )
}


if (
  length(
    master_city_ids_missing_climate
  ) > 0
) {
  
  stop(
    "Some master cities are missing from climate table."
  )
}


if (
  length(
    unexpected_climate_years
  ) > 0 ||
  
  length(
    missing_expected_climate_years
  ) > 0
) {
  
  stop(
    "Climate years do not exactly equal 2015-2024."
  )
}


if (
  nrow(
    invalid_climate_rows
  ) > 0
) {
  
  stop(
    "Invalid or incomplete annual climate rows detected."
  )
}


if (
  formula_difference > 1e-8
) {
  
  stop(
    "Stored De Martonne values do not match the expected formula."
  )
}


# ------------------------------------------------------------
# 12. Reattach authoritative city metadata
# ------------------------------------------------------------

city_year_aridity <- merge(
  
  climate_numeric,
  
  city_master[
    ,
    .(
      city_id,
      NAME10,
      city_ascii,
      state_name,
      county_name
    )
  ],
  
  by = "city_id",
  all.x = TRUE,
  sort = FALSE
)


setcolorder(
  city_year_aridity,
  c(
    "city_id",
    "NAME10",
    "city_ascii",
    "state_name",
    "county_name",
    "year",
    "total_annual_precip_mm",
    "mean_annual_temp_c",
    "n_months_with_climate",
    "aridity_index",
    "aridity_class_de_martonne"
  )
)


setorder(
  city_year_aridity,
  state_name,
  city_ascii,
  year
)


# Refresh the processed climate file even when an older valid
# climate table was reused. This guarantees authoritative city
# metadata in the final version.

atomic_fwrite(
  city_year_aridity,
  city_year_aridity_out
)


atomic_fwrite(
  city_year_aridity,
  city_year_aridity_for_analysis_out
)


cat("\nCity-year WorldClim aridity validated.\n")
cat("Rows:\n")
print(nrow(city_year_aridity))


# ------------------------------------------------------------
# 13. Build unique aridity lookup
# ------------------------------------------------------------

aridity_lookup <- city_year_aridity[
  ,
  .(
    city_id,
    year,
    total_annual_precip_mm,
    mean_annual_temp_c,
    aridity_index,
    aridity_class_de_martonne,
    n_months_with_climate
  )
]


if (
  anyDuplicated(
    paste(
      aridity_lookup$city_id,
      aridity_lookup$year,
      sep = "__"
    )
  ) > 0
) {
  
  stop(
    "Aridity lookup contains duplicate city-year keys."
  )
}


# ------------------------------------------------------------
# 14. Copy non-aridity Script 05 products to For Analysis
# ------------------------------------------------------------

copy_results <- c(
  
  file.copy(
    city_year_richness_path,
    city_year_without_aridity_for_analysis_out,
    overwrite = TRUE
  ),
  
  file.copy(
    checklist_richness_path,
    checklist_without_aridity_for_analysis_out,
    overwrite = TRUE
  ),
  
  file.copy(
    hotspot_list_path,
    hotspot_list_for_analysis_out,
    overwrite = TRUE
  ),
  
  file.copy(
    hotspot_latlong_path,
    hotspot_latlong_for_analysis_out,
    overwrite = TRUE
  )
)


if (
  !all(
    copy_results
  )
) {
  
  stop(
    "One or more Script 05 files could not be copied to For Analysis."
  )
}


cat("\nScript 05 outputs copied to For Analysis.\n")


# ------------------------------------------------------------
# 15. Join aridity to combined checklist richness
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("JOINING ARIDITY TO CHECKLIST RICHNESS\n")
cat("============================================================\n\n")


# Protect against accidental duplicate climate columns.

climate_cols_to_add <- c(
  "total_annual_precip_mm",
  "mean_annual_temp_c",
  "aridity_index",
  "aridity_class_de_martonne",
  "n_months_with_climate"
)


existing_climate_cols <- intersect(
  climate_cols_to_add,
  names(checklist_richness)
)


if (
  length(existing_climate_cols) > 0
) {
  
  stop(
    "Checklist richness already contains climate columns:\n",
    paste(
      existing_climate_cols,
      collapse = "\n"
    )
  )
}


checklist_richness[
  ,
  script06_row_order := .I
]


checklist_richness_aridity <- merge(
  
  checklist_richness,
  
  aridity_lookup,
  
  by = c(
    "city_id",
    "year"
  ),
  
  all.x = TRUE,
  sort = FALSE
)


setorder(
  checklist_richness_aridity,
  script06_row_order
)


checklist_richness_aridity[
  ,
  script06_row_order := NULL
]


if (
  nrow(checklist_richness_aridity) !=
  input_checklist_rows
) {
  
  stop(
    "Checklist row count changed during aridity join."
  )
}


if (
  uniqueN(
    checklist_richness_aridity$checklist_id
  ) !=
  input_unique_checklists
) {
  
  stop(
    "Checklist IDs changed during aridity join."
  )
}


historical_checklist_missing <- checklist_richness_aridity[
  year %in% climate_years &
    is.na(aridity_index)
]


checklist_2025_rows <- checklist_richness_aridity[
  year == 2025,
  .N
]


checklist_2025_with_climate <- checklist_richness_aridity[
  year == 2025 &
    !is.na(aridity_index),
  .N
]


cat("Checklist rows after join:\n")
print(nrow(checklist_richness_aridity))


cat("\n2015-2024 checklist rows missing climate:\n")
print(nrow(historical_checklist_missing))


cat("\n2025 checklist rows:\n")
print(checklist_2025_rows)


cat("\n2025 checklist rows unexpectedly matched to climate:\n")
print(checklist_2025_with_climate)


if (
  nrow(
    historical_checklist_missing
  ) > 0
) {
  
  stop(
    "One or more 2015-2024 checklist rows failed to match climate."
  )
}


if (
  checklist_2025_with_climate > 0
) {
  
  stop(
    "2025 checklist rows unexpectedly received WorldClim climate."
  )
}


atomic_fwrite(
  checklist_richness_aridity,
  checklist_richness_aridity_out
)


# ------------------------------------------------------------
# 16. Join aridity to combined city-year richness
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("JOINING ARIDITY TO CITY-YEAR RICHNESS\n")
cat("============================================================\n\n")


existing_city_year_climate_cols <- intersect(
  climate_cols_to_add,
  names(city_year_richness)
)


if (
  length(
    existing_city_year_climate_cols
  ) > 0
) {
  
  stop(
    "City-year richness already contains climate columns."
  )
}


city_year_richness[
  ,
  script06_row_order := .I
]


city_year_richness_aridity <- merge(
  
  city_year_richness,
  
  aridity_lookup,
  
  by = c(
    "city_id",
    "year"
  ),
  
  all.x = TRUE,
  sort = FALSE
)


setorder(
  city_year_richness_aridity,
  script06_row_order
)


city_year_richness_aridity[
  ,
  script06_row_order := NULL
]


if (
  nrow(city_year_richness_aridity) !=
  input_city_year_rows
) {
  
  stop(
    "City-year row count changed during aridity join."
  )
}


historical_city_year_missing <- city_year_richness_aridity[
  year %in% climate_years &
    is.na(aridity_index)
]


city_year_2025_rows <- city_year_richness_aridity[
  year == 2025,
  .N
]


city_year_2025_with_climate <- city_year_richness_aridity[
  year == 2025 &
    !is.na(aridity_index),
  .N
]


cat("City-year rows after join:\n")
print(nrow(city_year_richness_aridity))


cat("\n2015-2024 city-year rows missing climate:\n")
print(nrow(historical_city_year_missing))


cat("\n2025 city-year rows:\n")
print(city_year_2025_rows)


cat("\n2025 city-year rows unexpectedly matched to climate:\n")
print(city_year_2025_with_climate)


if (
  nrow(
    historical_city_year_missing
  ) > 0
) {
  
  stop(
    "One or more 2015-2024 city-year rows failed to match climate."
  )
}


if (
  city_year_2025_with_climate > 0
) {
  
  stop(
    "2025 city-year rows unexpectedly received WorldClim climate."
  )
}


atomic_fwrite(
  city_year_richness_aridity,
  city_year_richness_aridity_out
)


# ------------------------------------------------------------
# 17. Prepare per-city output folders
# ------------------------------------------------------------

if (
  overwrite_per_city_aridity_outputs
) {
  
  cat("\n============================================================\n")
  cat("RESETTING PER-CITY ARIDITY OUTPUTS\n")
  cat("============================================================\n\n")
  
  
  unlink(
    list.files(
      checklist_richness_aridity_by_city_dir,
      pattern = "\\.csv$",
      full.names = TRUE
    ),
    force = TRUE
  )
  
  
  unlink(
    list.files(
      city_year_aridity_by_city_dir,
      pattern = "\\.csv$",
      full.names = TRUE
    ),
    force = TRUE
  )
}


# ------------------------------------------------------------
# 18. Validate Script 05 per-city source files
# ------------------------------------------------------------

expected_checklist_city_files <-
  city_manifest$checklist_richness_file


expected_city_year_files <-
  city_manifest$city_year_unique_file


missing_source_checklist_files <- expected_checklist_city_files[
  !file.exists(
    expected_checklist_city_files
  )
]


missing_source_city_year_files <- expected_city_year_files[
  !file.exists(
    expected_city_year_files
  )
]


cat("\n============================================================\n")
cat("PER-CITY SOURCE FILE QC\n")
cat("============================================================\n\n")


cat("Expected per-city checklist files:\n")
print(length(expected_checklist_city_files))


cat("\nMissing per-city checklist files:\n")
print(length(missing_source_checklist_files))


cat("\nExpected per-city city-year files:\n")
print(length(expected_city_year_files))


cat("\nMissing per-city city-year files:\n")
print(length(missing_source_city_year_files))


if (
  length(
    missing_source_checklist_files
  ) > 0
) {
  
  stop(
    "Some Script 05 per-city checklist files are missing."
  )
}


if (
  length(
    missing_source_city_year_files
  ) > 0
) {
  
  stop(
    "Some Script 05 per-city city-year files are missing."
  )
}


# ------------------------------------------------------------
# 19. Join aridity to per-city checklist files
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("PER-CITY CHECKLIST + ARIDITY\n")
cat("============================================================\n\n")


status_rows <- vector(
  "list",
  expected_ebird_city_count * 2L
)


status_counter <- 1L


for (
  i in seq_len(
    nrow(city_manifest)
  )
) {
  
  in_file <-
    city_manifest$checklist_richness_file[i]
  
  
  out_file <- file.path(
    checklist_richness_aridity_by_city_dir,
    basename(in_file)
  )
  
  
  if (
    file.exists(out_file) &&
    file.info(out_file)$size > 0 &&
    !overwrite_per_city_aridity_outputs
  ) {
    
    x_existing <- fread(
      out_file,
      showProgress = FALSE
    )
    
    
    status_rows[[status_counter]] <- data.table(
      city_id = city_manifest$city_id[i],
      input_file = basename(in_file),
      output_file = basename(out_file),
      output_type =
        "checklist_richness_aridity_by_city",
      input_rows = NA_integer_,
      output_rows = nrow(x_existing),
      historical_missing_climate = NA_integer_,
      unexpected_2025_match = NA_integer_,
      status = "already_complete"
    )
    
    
    status_counter <- status_counter + 1L
    
    
    rm(
      x_existing
    )
    
    
    next
  }
  
  
  x <- fread(
    in_file,
    showProgress = FALSE
  )
  
  
  x[
    ,
    city_id := as.integer(city_id)
  ]
  
  
  x[
    ,
    year := as.integer(year)
  ]
  
  
  input_rows_i <- nrow(
    x
  )
  
  
  x[
    ,
    script06_row_order := .I
  ]
  
  
  x_aridity <- merge(
    
    x,
    
    aridity_lookup,
    
    by = c(
      "city_id",
      "year"
    ),
    
    all.x = TRUE,
    sort = FALSE
  )
  
  
  setorder(
    x_aridity,
    script06_row_order
  )
  
  
  x_aridity[
    ,
    script06_row_order := NULL
  ]
  
  
  historical_missing_i <- x_aridity[
    year %in% climate_years &
      is.na(aridity_index),
    .N
  ]
  
  
  unexpected_2025_i <- x_aridity[
    year == 2025 &
      !is.na(aridity_index),
    .N
  ]
  
  
  if (
    nrow(x_aridity) !=
    input_rows_i
  ) {
    
    stop(
      "Row count changed for per-city checklist file:\n",
      in_file
    )
  }
  
  
  if (
    historical_missing_i > 0
  ) {
    
    stop(
      "Historical climate missing in per-city checklist file:\n",
      in_file
    )
  }
  
  
  if (
    unexpected_2025_i > 0
  ) {
    
    stop(
      "2025 unexpectedly matched climate in:\n",
      in_file
    )
  }
  
  
  atomic_fwrite(
    x_aridity,
    out_file
  )
  
  
  status_rows[[status_counter]] <- data.table(
    city_id = city_manifest$city_id[i],
    input_file = basename(in_file),
    output_file = basename(out_file),
    output_type =
      "checklist_richness_aridity_by_city",
    input_rows = input_rows_i,
    output_rows = nrow(x_aridity),
    historical_missing_climate =
      historical_missing_i,
    unexpected_2025_match =
      unexpected_2025_i,
    status = "processed_now"
  )
  
  
  status_counter <- status_counter + 1L
  
  
  if (
    i %% 50 == 0 ||
    i == nrow(city_manifest)
  ) {
    
    cat(
      "Checklist files processed:",
      i,
      "of",
      nrow(city_manifest),
      "\n"
    )
  }
  
  
  rm(
    x,
    x_aridity
  )
}


# ------------------------------------------------------------
# 20. Join aridity to per-city city-year files
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("PER-CITY CITY-YEAR + ARIDITY\n")
cat("============================================================\n\n")


for (
  i in seq_len(
    nrow(city_manifest)
  )
) {
  
  in_file <-
    city_manifest$city_year_unique_file[i]
  
  
  out_file <- file.path(
    city_year_aridity_by_city_dir,
    basename(in_file)
  )
  
  
  if (
    file.exists(out_file) &&
    file.info(out_file)$size > 0 &&
    !overwrite_per_city_aridity_outputs
  ) {
    
    x_existing <- fread(
      out_file,
      showProgress = FALSE
    )
    
    
    status_rows[[status_counter]] <- data.table(
      city_id = city_manifest$city_id[i],
      input_file = basename(in_file),
      output_file = basename(out_file),
      output_type =
        "city_year_aridity_by_city",
      input_rows = NA_integer_,
      output_rows = nrow(x_existing),
      historical_missing_climate = NA_integer_,
      unexpected_2025_match = NA_integer_,
      status = "already_complete"
    )
    
    
    status_counter <- status_counter + 1L
    
    
    rm(
      x_existing
    )
    
    
    next
  }
  
  
  x <- fread(
    in_file,
    showProgress = FALSE
  )
  
  
  x[
    ,
    city_id := as.integer(city_id)
  ]
  
  
  x[
    ,
    year := as.integer(year)
  ]
  
  
  input_rows_i <- nrow(
    x
  )
  
  
  x[
    ,
    script06_row_order := .I
  ]
  
  
  x_aridity <- merge(
    
    x,
    
    aridity_lookup,
    
    by = c(
      "city_id",
      "year"
    ),
    
    all.x = TRUE,
    sort = FALSE
  )
  
  
  setorder(
    x_aridity,
    script06_row_order
  )
  
  
  x_aridity[
    ,
    script06_row_order := NULL
  ]
  
  
  historical_missing_i <- x_aridity[
    year %in% climate_years &
      is.na(aridity_index),
    .N
  ]
  
  
  unexpected_2025_i <- x_aridity[
    year == 2025 &
      !is.na(aridity_index),
    .N
  ]
  
  
  if (
    nrow(x_aridity) !=
    input_rows_i
  ) {
    
    stop(
      "Row count changed for per-city city-year file:\n",
      in_file
    )
  }
  
  
  if (
    historical_missing_i > 0
  ) {
    
    stop(
      "Historical climate missing in per-city city-year file:\n",
      in_file
    )
  }
  
  
  if (
    unexpected_2025_i > 0
  ) {
    
    stop(
      "2025 unexpectedly matched climate in:\n",
      in_file
    )
  }
  
  
  atomic_fwrite(
    x_aridity,
    out_file
  )
  
  
  status_rows[[status_counter]] <- data.table(
    city_id = city_manifest$city_id[i],
    input_file = basename(in_file),
    output_file = basename(out_file),
    output_type =
      "city_year_aridity_by_city",
    input_rows = input_rows_i,
    output_rows = nrow(x_aridity),
    historical_missing_climate =
      historical_missing_i,
    unexpected_2025_match =
      unexpected_2025_i,
    status = "processed_now"
  )
  
  
  status_counter <- status_counter + 1L
  
  
  if (
    i %% 50 == 0 ||
    i == nrow(city_manifest)
  ) {
    
    cat(
      "City-year files processed:",
      i,
      "of",
      nrow(city_manifest),
      "\n"
    )
  }
  
  
  rm(
    x,
    x_aridity
  )
}


# ------------------------------------------------------------
# 21. Save per-city processing status
# ------------------------------------------------------------

per_city_aridity_status <- rbindlist(
  status_rows,
  fill = TRUE
)


atomic_fwrite(
  per_city_aridity_status,
  per_city_aridity_status_out
)


cat("\nPer-city processing status:\n")

print(
  per_city_aridity_status[
    ,
    .N,
    by = .(
      output_type,
      status
    )
  ]
)


# ------------------------------------------------------------
# 22. Per-city output QC
# ------------------------------------------------------------

actual_checklist_aridity_files <- list.files(
  checklist_richness_aridity_by_city_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)


actual_city_year_aridity_files <- list.files(
  city_year_aridity_by_city_dir,
  pattern = "\\.csv$",
  full.names = TRUE
)


expected_checklist_aridity_basenames <- basename(
  expected_checklist_city_files
)


expected_city_year_aridity_basenames <- basename(
  expected_city_year_files
)


actual_checklist_aridity_basenames <- basename(
  actual_checklist_aridity_files
)


actual_city_year_aridity_basenames <- basename(
  actual_city_year_aridity_files
)


missing_checklist_aridity_files <- setdiff(
  expected_checklist_aridity_basenames,
  actual_checklist_aridity_basenames
)


extra_checklist_aridity_files <- setdiff(
  actual_checklist_aridity_basenames,
  expected_checklist_aridity_basenames
)


missing_city_year_aridity_files <- setdiff(
  expected_city_year_aridity_basenames,
  actual_city_year_aridity_basenames
)


extra_city_year_aridity_files <- setdiff(
  actual_city_year_aridity_basenames,
  expected_city_year_aridity_basenames
)


cat("\n============================================================\n")
cat("PER-CITY OUTPUT QC\n")
cat("============================================================\n\n")


cat("Checklist + aridity files:\n")
print(length(actual_checklist_aridity_files))


cat("\nMissing checklist + aridity files:\n")
print(length(missing_checklist_aridity_files))


cat("\nUnexpected checklist + aridity files:\n")
print(length(extra_checklist_aridity_files))


cat("\nCity-year + aridity files:\n")
print(length(actual_city_year_aridity_files))


cat("\nMissing city-year + aridity files:\n")
print(length(missing_city_year_aridity_files))


cat("\nUnexpected city-year + aridity files:\n")
print(length(extra_city_year_aridity_files))


if (
  length(
    actual_checklist_aridity_files
  ) != expected_ebird_city_count
) {
  
  stop(
    "Expected 647 checklist+aridity files."
  )
}


if (
  length(
    actual_city_year_aridity_files
  ) != expected_ebird_city_count
) {
  
  stop(
    "Expected 647 city-year+aridity files."
  )
}


if (
  length(
    missing_checklist_aridity_files
  ) > 0 ||
  
  length(
    extra_checklist_aridity_files
  ) > 0 ||
  
  length(
    missing_city_year_aridity_files
  ) > 0 ||
  
  length(
    extra_city_year_aridity_files
  ) > 0
) {
  
  stop(
    "Per-city climate-enhanced file sets do not match Script 05."
  )
}


if (
  nrow(
    per_city_aridity_status
  ) !=
  expected_ebird_city_count * 2L
) {
  
  stop(
    "Expected 1,294 per-city processing-status rows."
  )
}


if (
  any(
    per_city_aridity_status$input_rows !=
    per_city_aridity_status$output_rows,
    na.rm = TRUE
  )
) {
  
  stop(
    "Per-city row counts changed during climate joins."
  )
}


if (
  any(
    per_city_aridity_status$historical_missing_climate > 0,
    na.rm = TRUE
  )
) {
  
  stop(
    "Historical climate is missing from one or more per-city outputs."
  )
}


if (
  any(
    per_city_aridity_status$unexpected_2025_match > 0,
    na.rm = TRUE
  )
) {
  
  stop(
    "2025 climate was unexpectedly matched in a per-city output."
  )
}


# ------------------------------------------------------------
# 23. Hotspot / locality QC
# ------------------------------------------------------------

hotspot_list <- fread(
  hotspot_list_path,
  showProgress = FALSE
)


hotspot_latlong <- fread(
  hotspot_latlong_path,
  showProgress = FALSE
)


hotspot_list_localities <- uniqueN(
  hotspot_list$locality_id
)


hotspot_latlong_localities <- uniqueN(
  hotspot_latlong$locality_id
)


cat("\n============================================================\n")
cat("HOTSPOT / LOCALITY QC\n")
cat("============================================================\n\n")


cat("Checklist locality IDs:\n")
print(input_locality_count)


cat("\nHotspot-list locality IDs:\n")
print(hotspot_list_localities)


cat("\nHotspot-latlong locality IDs:\n")
print(hotspot_latlong_localities)


cat("\nHotspot-list rows:\n")
print(nrow(hotspot_list))


cat("\nHotspot-latlong rows:\n")
print(nrow(hotspot_latlong))


cat("\nSum of hotspot checklist counts:\n")

hotspot_checklist_total <- sum(
  hotspot_list$n_checklists,
  na.rm = TRUE
)

print(
  hotspot_checklist_total
)


if (
  hotspot_list_localities !=
  input_locality_count
) {
  
  stop(
    "Hotspot list does not preserve all Script 05 locality IDs."
  )
}


if (
  hotspot_latlong_localities !=
  input_locality_count
) {
  
  stop(
    "Hotspot lat/long file does not preserve all Script 05 locality IDs."
  )
}


if (
  hotspot_checklist_total !=
  input_checklist_rows
) {
  
  stop(
    "Hotspot checklist total does not equal Script 05 checklist count."
  )
}


# ------------------------------------------------------------
# 24. Confirm final analysis files
# ------------------------------------------------------------

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


final_file_check <- data.table(
  
  file = final_files,
  
  exists = file.exists(
    final_files
  ),
  
  size_mb = round(
    file.info(
      final_files
    )$size / 1024^2,
    2
  )
)


cat("\n============================================================\n")
cat("FINAL ANALYSIS FILES\n")
cat("============================================================\n\n")


print(
  final_file_check
)


if (
  !all(
    final_file_check$exists
  )
) {
  
  stop(
    "One or more final Script 06 outputs are missing."
  )
}


# ------------------------------------------------------------
# 25. Final hard QC
# ------------------------------------------------------------

final_checklist_rows <- nrow(
  checklist_richness_aridity
)


final_checklist_ids <- uniqueN(
  checklist_richness_aridity$checklist_id
)


final_ebird_cities <- uniqueN(
  checklist_richness_aridity$city_id
)


final_localities <- uniqueN(
  checklist_richness_aridity$locality_id[
    !is.na(checklist_richness_aridity$locality_id) &
      checklist_richness_aridity$locality_id != ""
  ]
)


final_city_year_rows <- nrow(
  city_year_richness_aridity
)


checklist_historical_matched <- checklist_richness_aridity[
  year %in% climate_years &
    !is.na(aridity_index),
  .N
]


checklist_historical_total <- checklist_richness_aridity[
  year %in% climate_years,
  .N
]


city_year_historical_matched <- city_year_richness_aridity[
  year %in% climate_years &
    !is.na(aridity_index),
  .N
]


city_year_historical_total <- city_year_richness_aridity[
  year %in% climate_years,
  .N
]


cat("\n============================================================\n")
cat("SCRIPT 06 FINAL QC\n")
cat("============================================================\n\n")


cat("Climate city-year rows:\n")
print(nrow(city_year_aridity))


cat("\nChecklist rows before join:\n")
print(input_checklist_rows)


cat("\nChecklist rows after join:\n")
print(final_checklist_rows)


cat("\nUnique checklist IDs:\n")
print(final_checklist_ids)


cat("\nCities:\n")
print(final_ebird_cities)


cat("\nLocalities:\n")
print(final_localities)


cat("\nCity-year rows before join:\n")
print(input_city_year_rows)


cat("\nCity-year rows after join:\n")
print(final_city_year_rows)


cat("\n2015-2024 checklist rows:\n")
print(checklist_historical_total)


cat("\n2015-2024 checklist rows with climate:\n")
print(checklist_historical_matched)


cat("\n2025 checklist rows:\n")
print(checklist_2025_rows)


cat("\n2015-2024 city-year rows:\n")
print(city_year_historical_total)


cat("\n2015-2024 city-year rows with climate:\n")
print(city_year_historical_matched)


cat("\n2025 city-year rows:\n")
print(city_year_2025_rows)


cat("\nPer-city checklist+aridity files:\n")
print(length(actual_checklist_aridity_files))


cat("\nPer-city city-year+aridity files:\n")
print(length(actual_city_year_aridity_files))


stopifnot(
  
  # Climate table
  
  nrow(city_year_aridity) ==
    expected_master_city_count *
    length(climate_years),
  
  uniqueN(city_year_aridity$city_id) ==
    expected_master_city_count,
  
  all(
    city_year_aridity$n_months_with_climate ==
      12
  ),
  
  
  # Checklist table
  
  final_checklist_rows ==
    input_checklist_rows,
  
  final_checklist_ids ==
    input_unique_checklists,
  
  final_ebird_cities ==
    input_ebird_city_count,
  
  final_localities ==
    input_locality_count,
  
  checklist_historical_matched ==
    checklist_historical_total,
  
  checklist_2025_with_climate ==
    0,
  
  
  # City-year table
  
  final_city_year_rows ==
    input_city_year_rows,
  
  city_year_historical_matched ==
    city_year_historical_total,
  
  city_year_2025_with_climate ==
    0,
  
  
  # Hotspots
  
  hotspot_list_localities ==
    input_locality_count,
  
  hotspot_latlong_localities ==
    input_locality_count,
  
  hotspot_checklist_total ==
    input_checklist_rows,
  
  
  # Per-city outputs
  
  length(actual_checklist_aridity_files) ==
    expected_ebird_city_count,
  
  length(actual_city_year_aridity_files) ==
    expected_ebird_city_count,
  
  nrow(per_city_aridity_status) ==
    expected_ebird_city_count * 2L
)


# ------------------------------------------------------------
# 26. Final result
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("SCRIPT 06 QC PASSED\n")
cat("============================================================\n\n")


cat("WorldClim years:                2015-2024\n")

cat(
  "Climate city-year rows:       ",
  nrow(city_year_aridity),
  "\n"
)

cat(
  "eBird checklist rows:          ",
  final_checklist_rows,
  "\n"
)

cat(
  "Unique checklist IDs:          ",
  final_checklist_ids,
  "\n"
)

cat(
  "eBird cities:                  ",
  final_ebird_cities,
  "\n"
)

cat(
  "Unique locality IDs:           ",
  final_localities,
  "\n"
)

cat(
  "eBird city-year rows:          ",
  final_city_year_rows,
  "\n"
)

cat(
  "2015-2024 checklist climate:   ",
  checklist_historical_matched,
  " / ",
  checklist_historical_total,
  "\n",
  sep = ""
)

cat(
  "2025 checklist rows with AI:   ",
  checklist_2025_with_climate,
  "\n"
)

cat(
  "2015-2024 city-year climate:   ",
  city_year_historical_matched,
  " / ",
  city_year_historical_total,
  "\n",
  sep = ""
)

cat(
  "2025 city-year rows with AI:   ",
  city_year_2025_with_climate,
  "\n"
)

cat(
  "Per-city checklist files:      ",
  length(actual_checklist_aridity_files),
  "\n"
)

cat(
  "Per-city city-year files:      ",
  length(actual_city_year_aridity_files),
  "\n"
)


cat("\nMain analysis-ready outputs:\n")

cat(
  checklist_richness_aridity_out,
  "\n"
)

cat(
  city_year_richness_aridity_out,
  "\n"
)

cat(
  city_year_aridity_for_analysis_out,
  "\n"
)

cat(
  hotspot_latlong_for_analysis_out,
  "\n"
)