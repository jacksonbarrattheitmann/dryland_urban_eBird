## This script combines all of the .csv files
## that contain the land cover estimates for
## the eBird hotspots in 1km buffers

library(dplyr)
library(readr)

# creating a list of file names and reading them in
# total we have 427,963 LOCALITY_IDs, so that's how many rows we have
file_list <- list.files(path = "Intermediate_data/GEE_land_cover/batch_1st_attempt/", pattern = "\\.csv", full.names = TRUE)
r <- lapply(file_list, read.csv)
s <- bind_rows(r)
s <- as.data.frame(s)

# getting rid of some not so useful columns
s <- s %>%
  select(-.geo, -system.index, -batch) %>%
  rename(LOCALITY_ID = lclty_d)


# saving an RDS of the new combined file
saveRDS(s, file = "Intermediate_data/GEE_land_cover/gee_land_cover_combined.RDS")

gee_land_cover <- readRDS("Intermediate_data/GEE_land_cover/gee_land_cover_combined1.RDS")

### creating a list of the missed localities
file_list <- list.files(path = "Intermediate_data/GEE_land_cover/missed_localities/", pattern = "\\.csv", full.names = TRUE)
r <- lapply(file_list, read.csv)
s <- bind_rows(r)
s <- as.data.frame(s)

# getting rid of some not so useful columns
s <- s %>%
  select(-.geo, -system.index, -batch, -new_batch) %>%
  rename(LOCALITY_ID = lclty_d)

## combine these two, and then figure out what
## we're still missing from GEE

full <- gee_land_cover %>%
  rbind(s)


## Read in the second batch of missed localities in GEE
file_list <- list.files(path = "Intermediate_data/GEE_land_cover/missed_localities2/", pattern = "\\.csv", full.names = TRUE)
r <- lapply(file_list, read.csv)
s <- bind_rows(r)
s <- as.data.frame(s)

# getting rid of some not so useful columns
s <- s %>%
  select(-.geo, -system.index, -batch, -new_batch, -new_batch2) %>%
  rename(LOCALITY_ID = lclty_d)

### Now append this to the 'full' df

full <- full %>%
  rbind(s)


## now save this as the final GEE data frame

saveRDS(full, "Intermediate_data/GEE_land_cover/full_gee_combined.RDS")



## CODE WE HAD TO USE TO GRAB THE MISSING LOCALITY_IDs##

## Now intersect with our newest lat/long data file, and see what we're missing
##
# dat <- read.csv("Data/Processed/For Analysis/dryland_city_hotspot_latlong_for_GEE_2015_2025.csv")
# 
# still_missing <- dat %>%
#   filter(!locality_id %in% full$LOCALITY_ID)
# 
# saveRDS(still_missing, "Intermediate_data/GEE_land_cover/missed_localities2/missed_localities2.RDS")
