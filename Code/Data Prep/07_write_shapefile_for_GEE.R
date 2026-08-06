### let's try reading in a city

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


