### Prepare data for analysis script ###

# This is a first attempt at trying to build some GAM/GLM
# models for our urban dryland birds dataset. 

library(dplyr)
library(ggplot2)
library(readr)
library(tidyverse)
library(FD)
library(purrr)

# ---------------------------------------------------------
# Read in data that is the same for every city

land <- readRDS(
  "Intermediate_data/GEE_land_cover/full_gee_combined.RDS"
)

trait <- readRDS(
  "Intermediate_data/functional_traits_data_AVONET/AVONET_cleaned.RDS"
) %>%
  rename(scientific_name = new_scientific_name,
         Hand.Wing.Index = 'Hand-Wing.Index') %>%
  as.data.frame()


# ---------------------------------------------------------
# Function
# ---------------------------------------------------------

process_city <- function(tsv_file,
                         land,
                         trait,
                         output_dir = "Intermediate_data/dfs_for_analysis_by_city") {
  
  # -------------------------------------------------------
  # Get city name from file name
  # -------------------------------------------------------
  
  city_name <- tools::file_path_sans_ext(
    basename(tsv_file)
  )
  
  message("Processing: ", city_name)
  
  
  # -------------------------------------------------------
  # Read eBird data
  # -------------------------------------------------------
  
  bird_dat <- read_tsv(tsv_file, show_col_types = FALSE) %>%
    select(-locality.y, -locality_id.y, -observation_date.y) %>%
    rename(
      locality = locality.x,
      locality_id = locality_id.x,
      observation_date = observation_date.x
    )
  
  
  # -------------------------------------------------------
  # Clean eBird data
  # Remove checklists with X and spuh/slash observations
  # -------------------------------------------------------
  
  bird_clean <- bird_dat %>%
    group_by(checklist_id) %>%
    filter(
      !any(observation_count == "X"),
      !category == "spuh",
      !category == "slash",
      !category == "domestic", 
      !category == "hybrid"
    ) %>%
    ungroup()
  
  bird_clean$observation_count <- as.numeric(
    bird_clean$observation_count
  )
  
  
  # -------------------------------------------------------
  # Total abundance and richness per checklist
  # -------------------------------------------------------
  
  tot <- bird_clean %>%
    group_by(checklist_id) %>%
    summarize(
      sr = length(unique(common_name)),
      n = sum(observation_count),
      locality_id = first(locality_id),
      dryland_percent = first(dryland_percent)/100000000,
      NAME10 = first(NAME10),
      duration_minutes = first(duration_minutes),
      number_observers = first(number_observers),
      effort_distance_km = first(effort_distance_km),
      state_name = first(state_name),
      protocol_name = first(protocol_name),
      .groups = "drop"
    )
  
  
  # -------------------------------------------------------
  # Checklist x species matrix
  # -------------------------------------------------------
  
  bird_mat <- bird_clean %>%
    group_by(checklist_id, common_name) %>%
    summarise(
      abundance = sum(observation_count, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(common_name) %>%
    pivot_wider(
      names_from = common_name,
      values_from = abundance,
      values_fill = 0
    ) %>%
    as.data.frame()
  
  rownames(bird_mat) <- bird_mat$checklist_id
  bird_mat$checklist_id <- NULL
  
  
  # -------------------------------------------------------
  # Filter AVONET to species found in this city
  # -------------------------------------------------------
  
  trait_filt <- trait %>%
    filter(
      scientific_name %in% bird_clean$scientific_name
    )
  
  
  # -------------------------------------------------------
  # Check for taxonomy mismatches
  # -------------------------------------------------------
  
  check <- bird_clean %>%
    filter(
      !scientific_name %in% trait_filt$scientific_name
    )
  
  if (nrow(check) > 0) {
    warning(
      city_name,
      ": taxonomy check failed. ",
      nrow(check),
      " observations were not found in AVONET."
    )
  }
  
  
  # -------------------------------------------------------
  # Join common names to AVONET
  # -------------------------------------------------------
  
  bird_names <- bird_clean %>%
    select(common_name, scientific_name) %>%
    distinct()
  
  
  # -------------------------------------------------------
  # Prepare trait matrix
  # -------------------------------------------------------
  
  traits_fr <- trait_filt %>%
    select(
      scientific_name,
      Beak.Length_Culmen,
      Beak.Width,
      Tarsus.Length,
      Hand.Wing.Index,
      Tail.Length,
      Trophic.Level,
      Trophic.Niche,
      Mass
    ) %>%
    left_join(
      bird_names,
      by = "scientific_name"
    ) %>%
    arrange(common_name)
  
  
  rownames(traits_fr) <- traits_fr$common_name
  
  traits_fr$scientific_name <- NULL
  traits_fr$common_name <- NULL
  
  
  # -------------------------------------------------------
  # Calculate functional diversity
  # -------------------------------------------------------
  
  d <- gowdis(traits_fr)
  
  fd <- FD::dbFD(
    d,
    bird_mat,
    calc.FRic = TRUE,
    w.abun = TRUE,
    corr = "cailliez"
  )
  
  results <- data.frame(
    FRic = fd$FRic,
    FDis = fd$FDis,
    FEve = fd$FEve
  ) %>%
    tibble::rownames_to_column("checklist_id")
  
  
  # -------------------------------------------------------
  # Community weighted means
  # -------------------------------------------------------
  
  traits_for_join <- trait_filt %>%
    select(
      scientific_name,
      Beak.Length_Culmen,
      Beak.Width,
      Tarsus.Length,
      Hand.Wing.Index,
      Tail.Length,
      Trophic.Level,
      Trophic.Niche,
      Mass
    )
  
  
  bird_trait <- bird_clean %>%
    left_join(
      traits_for_join,
      by = "scientific_name"
    )
  
  bird_trait$observation_count <- as.numeric(
    bird_trait$observation_count
  )
  
  
  comm_weights <- bird_trait %>%
    group_by(checklist_id, common_name) %>%
    reframe(
      n = sum(observation_count),
      total.Beak.Length =
        observation_count * Beak.Length_Culmen,
      total.Beak.Width =
        observation_count * Beak.Width,
      total.Tarsus.Length =
        observation_count * Tarsus.Length,
      total.HWI =
        observation_count * Hand.Wing.Index,
      total.Tail.Length =
        observation_count * Tail.Length,
      total.Mass =
        observation_count * Mass
    )
  
  
  comm_weights_final <- comm_weights %>%
    group_by(checklist_id) %>%
    summarize(
      CWM.Beak.Length =
        sum(total.Beak.Length) / sum(n),
      CWM.Beak.Width =
        sum(total.Beak.Width) / sum(n),
      CWM.Tarsus.Length =
        sum(total.Tarsus.Length) / sum(n),
      CWM.HWI =
        sum(total.HWI) / sum(n),
      CWM.Tail.Length =
        sum(total.Tail.Length) / sum(n),
      CWM.Mass =
        sum(total.Mass) / sum(n),
      .groups = "drop"
    )
  
  
  # -------------------------------------------------------
  # Filter land-cover data to localities in this city
  # -------------------------------------------------------
  
  land_filt <- land %>%
    rename(locality_id = LOCALITY_ID) %>%
    filter(
      locality_id %in% tot$locality_id
    )
  
  
  # -------------------------------------------------------
  # Combine everything
  # -------------------------------------------------------
  
  big_df <- tot %>%
    left_join(
      comm_weights_final,
      by = "checklist_id"
    ) %>%
    left_join(
      results,
      by = "checklist_id"
    ) %>%
    left_join(
      land_filt,
      by = "locality_id"
    )
  
  
  # -------------------------------------------------------
  # Save
  # -------------------------------------------------------
  
  output_file <- file.path(
    output_dir,
    paste0(city_name, ".RDS")
  )
  
  saveRDS(
    big_df,
    output_file
  )
  
  message(
    "Finished: ", city_name,
    " | ",
    nrow(big_df),
    " rows | saved to ",
    output_file
  )
  
  # Return the dataframe
  return(big_df)
}



tsv_files <- list.files(
  "Data/Processed/eBird/Observations_By_City_With_Checklist_Metadata/",
  pattern = "\\.tsv$",
  full.names = TRUE
)

for (file in tsv_files[509:647]) {
  
  process_city(
    tsv_file = file,
    land = land,
    trait = trait
  )
  
}


### Now how about rbinding all of these together
### for Brittany's analysis script

rds_folder <- "Intermediate_data/dfs_for_analysis_by_city/"

rds_files <- list.files(
  path = rds_folder,
  pattern = "\\.RDS$",
  full.names = TRUE
)

combined_data <- map_dfr(
  rds_files,
  ~ {
    message("Reading: ", basename(.x))
    readRDS(.x)
  }
)

## read in the aridity data
arid <- read.csv("Data/Processed/Climate/dryland_city_year_worldclim_aridity_2015_2024.csv")

arid_sum <- arid %>%
  group_by(NAME10) %>%
  summarize(avg_aridity_index = mean(aridity_index),
            avg_precp = mean(total_annual_precip_mm),
            avg_annual_temp = mean(mean_annual_temp_c))

## now join this to the big_data_frame

combined_data <- combined_data %>%
  left_join(arid_sum, by = "NAME10")

saveRDS(combined_data,"Intermediate_data/dfs_for_analysis_combined/combined_data.RDS")

## let's do another quick export with a filtered subset, because our land cover 
## is still not done for every LOCALITY_ID

combined_data_filt <- combined_data %>%
  filter(locality_id %in% land$LOCALITY_ID)

saveRDS(combined_data_filt, "Intermediate_data/dfs_for_analysis_combined/combined_data_filt.RDS")





## Code to check a single city for the broken name issue with AVONET ##



bird_dat <- read_tsv("Data/Processed/eBird/Observations_By_City_With_Checklist_Metadata/city_2311_california_san_diego.tsv") %>%
  select(-locality.y, -locality_id.y, -observation_date.y) %>%
  rename(
    locality = locality.x,
    locality_id = locality_id.x,
    observation_date = observation_date.x)

bird_clean <- bird_dat %>%
  group_by(checklist_id) %>%
  filter(
    !any(observation_count == "X"),
    !category == "spuh",
    !category == "slash",
    !category == "domestic", 
    !category == "hybrid"
  ) %>%
  ungroup()

bird_clean$observation_count <- as.numeric(
  bird_clean$observation_count
)

trait_filt <- trait %>%
  filter(
    scientific_name %in% bird_clean$scientific_name
  )

check <- bird_clean %>%
  filter(
    !scientific_name %in% trait_filt$scientific_name
  )
