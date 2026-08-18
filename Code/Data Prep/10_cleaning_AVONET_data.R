## Cleaning the AVONET Data

library(readr)
library(tidyverse)

## read in the AVONET data
trait <- readxl::read_xlsx("Intermediate_data/functional_traits_data_AVONET/AVONET Supplementary dataset 1.xlsx",
                           sheet = "AVONET2_eBird") %>%
  rename(old_scientific_name = Species2)

## table of taxonomy with old and new scientific names for harmonization
tax_bird <- readxl::read_xlsx("Intermediate_data/functional_traits_data_AVONET/AVONET_eBird_crosswalk_ALL_SPECIES.xlsx") %>%
  select(common_name, old_scientific_name, new_scientific_name)

## let's just join them, so we have all the possible species in eBird

trait_filtered <- trait %>%
  filter(old_scientific_name %in% tax_bird$old_scientific_name) %>%
  left_join(tax_bird, by = "old_scientific_name", relationship = "many-to-many")


## quick check for what we lost, only extinct species, so that's fine
lost_spp <- tax_bird %>%
  filter(!old_scientific_name %in% trait_filtered$old_scientific_name)


saveRDS(trait_filtered, "Intermediate_data/functional_traits_data_AVONET/AVONET_cleaned.RDS")




