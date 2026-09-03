## This script combines all of the .csv files
## that contain the land cover estimates for
## the eBird hotspots in 1km buffers

library(dplyr)
library(readr)

test <- read.csv("Intermediate_data/GEE_land_cover/batch_1st_attempt/DynamicWorld_2021_batch_1.csv")

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

gee_land_cover <- readRDS("Intermediate_data/GEE_land_cover/gee_land_cover_combined.RDS")
