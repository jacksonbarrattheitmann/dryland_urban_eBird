# ============================================================
# 02_Find_Dryland_Cities.R
#
# Purpose:
# Identify which city polygons fall primarily within U.S. drylands.
#
# Inputs:
# 1. City polygons from the previous urban biodiversity paper (Callaghan et all 2020)
# 2. Cleaned U.S. dryland shapefile from 01_Set_Up_Dryland_Shapefil.R
#
# Rule:
# A city is classified as a "dryland city" if at least 50% of its
# polygon area overlaps the cleaned dryland shapefile.
# ============================================================

library(sf)
library(dplyr)
library(ggplot2)
library(readr)

# ------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------

cities_path <- "Data/Raw/Shape Files/Cities/all_cities.geojson"

drylands_path <- "Data/Processed/Shape Files/Dryland/drylands_US_lower48_UNCCD_clean.shp"

city_shape_out_dir <- "Data/Processed/Shape Files/Cities"
table_out_dir <- "Data/Processed/Tables"
map_out_dir <- "Outputs/Maps"

dir.create(city_shape_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(map_out_dir, recursive = TRUE, showWarnings = FALSE)

# Output files
all_cities_shp_out <- file.path(
  city_shape_out_dir,
  "cities_with_dryland_overlap.shp"
)

dryland_cities_shp_out <- file.path(
  city_shape_out_dir,
  "dryland_cities_50pct.shp"
)

all_cities_csv_out <- file.path(
  table_out_dir,
  "cities_with_dryland_overlap.csv"
)

dryland_cities_csv_out <- file.path(
  table_out_dir,
  "dryland_cities_50pct.csv"
)

dryland_cities_map_out <- file.path(
  map_out_dir,
  "dryland_cities_50pct_map.png"
)

# ------------------------------------------------------------
# 2. Read spatial data
# ------------------------------------------------------------

cities <- st_read(cities_path)
drylands <- st_read(drylands_path)

# Make geometries valid before spatial operations
cities <- st_make_valid(cities)
drylands <- st_make_valid(drylands)

# ------------------------------------------------------------
# 3. Project to equal-area CRS
# ------------------------------------------------------------

# EPSG:5070 = NAD83 / Conus Albers.
# This is used so area calculations are in square meters and
# are appropriate across the contiguous United States.
cities_aea <- st_transform(cities, 5070)
drylands_aea <- st_transform(drylands, 5070)

# ------------------------------------------------------------
# 4. Prepare city and dryland geometries
# ------------------------------------------------------------

# Add a simple unique ID to each city. This avoids relying on
# city names, which can be duplicated across states.
cities_aea <- cities_aea %>%
  mutate(
    city_id = row_number(),
    city_area_m2 = as.numeric(st_area(geometry))
  )

# Combine all dryland polygons into one dryland boundary.
# We only need the dryland geometry here, not the dryland classes.
drylands_union <- drylands_aea %>%
  summarise(geometry = st_union(geometry)) %>%
  st_make_valid()

# ------------------------------------------------------------
# 5. Calculate dryland overlap for each city
# ------------------------------------------------------------

# Intersect city polygons with drylands. The resulting geometry
# represents only the portion of each city that falls inside drylands.
city_dryland_intersection <- st_intersection(
  cities_aea,
  drylands_union
) %>%
  mutate(
    dryland_overlap_m2 = as.numeric(st_area(geometry))
  )

# Sum dryland overlap area for each city
city_dryland_summary <- city_dryland_intersection %>%
  st_drop_geometry() %>%
  group_by(city_id) %>%
  summarise(
    dryland_overlap_m2 = sum(dryland_overlap_m2, na.rm = TRUE),
    .groups = "drop"
  )

# Join dryland overlap back to the full city layer
cities_with_dryland <- cities_aea %>%
  left_join(city_dryland_summary, by = "city_id") %>%
  mutate(
    dryland_overlap_m2 = if_else(
      is.na(dryland_overlap_m2),
      0,
      dryland_overlap_m2
    ),
    dryland_percent = 100 * dryland_overlap_m2 / city_area_m2,
    
    # 50% rule:
    # A city is considered a dryland city if at least half
    # of its polygon area overlaps the dryland shapefile.
    dryland_city_50pct = dryland_percent >= 50
  )

# ------------------------------------------------------------
# 6. Filter to dryland cities
# ------------------------------------------------------------

dryland_cities_50pct <- cities_with_dryland %>%
  filter(dryland_city_50pct)

cat("\nTotal cities:\n")
print(nrow(cities_with_dryland))

cat("\nCities with >=50% dryland overlap:\n")
print(nrow(dryland_cities_50pct))

cat("\nTop dryland-overlap cities:\n")

top_dryland_cities <- dryland_cities_50pct %>%
  st_drop_geometry() %>%
  arrange(desc(dryland_percent)) %>%
  select(
    city_id,
    NAME10,
    city_ascii,
    state_name,
    county_name,
    lat,
    lng,
    dryland_percent
  ) %>%
  head(50)

print(top_dryland_cities)

# ------------------------------------------------------------
# 7. Save shapefiles
# ------------------------------------------------------------

st_write(
  cities_with_dryland,
  all_cities_shp_out,
  delete_dsn = TRUE
)

st_write(
  dryland_cities_50pct,
  dryland_cities_shp_out,
  delete_dsn = TRUE
)

# ------------------------------------------------------------
# 8. Save tables
# ------------------------------------------------------------

all_cities_table <- cities_with_dryland %>%
  st_drop_geometry() %>%
  arrange(state_name, city_ascii) %>%
  select(
    city_id,
    GEOID10,
    UACE10,
    NAME10,
    city_ascii,
    state_name,
    county_name,
    lat,
    lng,
    population,
    population_proper,
    density,
    ALAND10,
    AWATER10,
    city_area_m2,
    dryland_overlap_m2,
    dryland_percent,
    dryland_city_50pct
  )

dryland_cities_table <- dryland_cities_50pct %>%
  st_drop_geometry() %>%
  arrange(state_name, city_ascii) %>%
  select(
    city_id,
    GEOID10,
    UACE10,
    NAME10,
    city_ascii,
    state_name,
    county_name,
    lat,
    lng,
    population,
    population_proper,
    density,
    ALAND10,
    AWATER10,
    city_area_m2,
    dryland_overlap_m2,
    dryland_percent
  )

write_csv(all_cities_table, all_cities_csv_out)
write_csv(dryland_cities_table, dryland_cities_csv_out)

# ------------------------------------------------------------
# 9. Map dryland cities
# ------------------------------------------------------------

dryland_city_map <- ggplot() +
  geom_sf(
    data = drylands_union,
    fill = "orange",
    color = NA,
    alpha = 0.35
  ) +
  geom_sf(
    data = cities_aea,
    fill = NA,
    color = "grey75",
    linewidth = 0.05
  ) +
  geom_sf(
    data = dryland_cities_50pct,
    fill = "red",
    color = "black",
    linewidth = 0.1,
    alpha = 0.75
  ) +
  theme_minimal() +
  labs(
    x = NULL,
    y = NULL
  )

dryland_city_map

ggsave(
  filename = dryland_cities_map_out,
  plot = dryland_city_map,
  width = 10,
  height = 6,
  dpi = 300
)