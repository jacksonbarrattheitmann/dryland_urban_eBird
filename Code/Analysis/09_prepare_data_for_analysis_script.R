### Prepare data for analysis script ###

# This is a first attempt at trying to build some GAM/GLM
# models for our urban dryland birds dataset. 

library(dplyr)
library(ggplot2)
library(readr)

# first let's grab some data
bird_dat <- read_tsv("Data/Processed/eBird/")

head(bird_dat)

## land cover data
land <- readRDS("Intermediate_data/GEE_land_cover/gee_land_cover_combined.RDS")

