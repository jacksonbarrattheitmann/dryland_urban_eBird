## This script was used to assign batch numbers to all of the possible locality_ids 
## that we had to get the land cover data for from GEE

library(tidyverse)
library(dplyr)
library(readr)
library(sf)

dat <- read.csv("Data/Processed/For Analysis/dryland_city_hotspot_latlong_for_GEE_2015_2025.csv")

dat_sf <- st_as_sf(dat, coords = c("longitude", "latitude"), crs = 4326)

# now let's just get the locality_id and the geo

dat_sf <- dat_sf %>%
  select(locality_id, geometry)

dat_sf <- dat_sf %>%
  mutate(batch = ceiling(row_number() / 2500))

## now let's export this as a shapefile for GEE

write_sf(dat_sf, "Intermediate_data/locality_id_shapefile/eBird_localities_drylands.shp")

#### need to probably subset these to make them manageable

batch_size <- 2500
n <- nrow(dat_sf)

# Starting indices for each batch
starts <- seq(1, n, by = batch_size)

for (i in seq_along(starts)) {
  
  start_row <- starts[i]
  end_row <- min(start_row + batch_size - 1, n)
  
  batch <- dat_sf[start_row:end_row, ]
  
  st_write(
    batch,
    paste0(
      "Intermediate_data/locality_id_shapefile/ebird_batch_",
      i,
      ".shp"
    ),
    delete_dsn = TRUE,
    quiet = TRUE
  )
  
  message("Exported batch ", i, ": rows ", start_row, "–", end_row)
}

## Need to make these even smaller for the batches that failed the first time
## read in the final dbf that we made

ebird_localities <- read_sf("Intermediate_data/locality_id_shapefile/eBird_localities_drylands.shp")

missed_batches <- c("11", "14", "18", "28", "31", "33", "36", "37", "41", "43", "44", "45", "51", "56", "62", "66", "69", "70", "72")

## now filter for the missed ones

ebird_locals_missed <- ebird_localities %>%
  filter(batch %in% missed_batches)

## we need to make the batches smaller
ebird_locals_missed <- ebird_locals_missed %>%
  mutate(new_batch = rep(1:46, each = 1000, length.out = n()))

## save this dataframe
st_write(ebird_locals_missed, "Intermediate_data/locality_id_shapefile/ebird_localities_missed.shp")

## read it in to deal with another failure

ebird_missed <- read_sf("Intermediate_data/locality_id_shapefile/ebird_localities_missed.shp")

#### We had another failure ####
missed_batches_2nd_attempt <- c("25", "29", "33", "36")

## now filter again, and chop it up again into smaller pieces
## since GEE is abyssmal at dealing with large amounts of data

ebird_locals_missed2 <- ebird_missed %>%
  filter(new_batch %in% missed_batches_2nd_attempt)

ebird_locals_missed2 <- ebird_locals_missed2 %>%
  mutate(new_batch2 = rep(1:8, each = 500, length.out = n()))

## save this dataframe, and try again to make it work on GEE
st_write(ebird_locals_missed2, "Intermediate_data/locality_id_shapefile/ebird_localities_missed2.shp", append = FALSE)
