# ============================================================
# 01_Set_Up_Dryland_Shapefil.R
# ============================================================

library(sf)
library(dplyr)
library(tigris)
library(ggplot2)

options(tigris_use_cache = TRUE)

# ------------------------------------------------------------
# 1. File paths
# ------------------------------------------------------------

# All raw shapefile pieces must be in this folder:
# .shp, .dbf, .shx, .prj
drylands_path <- "Data/Raw/Shape Files/drylands_UNCCD_CBD_july2014.shp"

# Processed shapefile output folder
shape_out_dir <- "Data/Processed/Shape Files"
dir.create(shape_out_dir, recursive = TRUE, showWarnings = FALSE)

# Map output folder
map_out_dir <- "Outputs/Maps"
dir.create(map_out_dir, recursive = TRUE, showWarnings = FALSE)

# Raw clipped output
out_path_raw <- file.path(shape_out_dir, "drylands_US_lower48_raw.shp")

# Cleaned UNCCD drylands output
out_path_clean <- file.path(shape_out_dir, "drylands_US_lower48_UNCCD_clean.shp")

# Saved map output
map_out_path <- file.path(map_out_dir, "drylands_US_lower48_UNCCD_clean_map.png")

# ------------------------------------------------------------
# 2. Read global drylands shapefile
# ------------------------------------------------------------

drylands_global <- st_read(drylands_path)

# Inspect fields
names(drylands_global)

# Inspect original dryland classes
drylands_global %>%
  st_drop_geometry() %>%
  count(HIX_NAME, HIX_DESC, UNCCD, UNCCDDESC, CBDDRYLAND, CBDDRYDEF, sort = TRUE)

# ------------------------------------------------------------
# 3. Download U.S. state boundaries
# ------------------------------------------------------------

us_states <- states(cb = TRUE, year = 2023, class = "sf")

# Keep contiguous U.S. only
lower48_states <- us_states %>%
  filter(!STUSPS %in% c(
    "AK", "HI", "PR", "GU", "VI", "MP", "AS"
  ))

# ------------------------------------------------------------
# 4. Reproject states to match drylands CRS
# ------------------------------------------------------------

lower48_states <- st_transform(lower48_states, st_crs(drylands_global))

# ------------------------------------------------------------
# 5. Dissolve lower-48 states into one boundary
# ------------------------------------------------------------

lower48_boundary <- lower48_states %>%
  summarise(geometry = st_union(geometry)) %>%
  st_make_valid()

drylands_global <- st_make_valid(drylands_global)

# ------------------------------------------------------------
# 6. Clip drylands to lower 48 U.S.
# ------------------------------------------------------------

drylands_US_lower48_raw <- st_intersection(
  drylands_global,
  lower48_boundary
)

# ------------------------------------------------------------
# 7. Create clean dryland class field
# ------------------------------------------------------------

drylands_US_lower48_clean <- drylands_US_lower48_raw %>%
  mutate(
    dryland_class = case_when(
      HIX_DESC == "Hyperarid" ~ "hyperarid: P/PET < 0.05",
      !is.na(UNCCDDESC) ~ UNCCDDESC,
      TRUE ~ NA_character_
    )
  )

# Inspect clipped classes before filtering
drylands_US_lower48_clean %>%
  st_drop_geometry() %>%
  count(HIX_NAME, HIX_DESC, UNCCD, UNCCDDESC, CBDDRYLAND, CBDDRYDEF, dryland_class, sort = TRUE)

# ------------------------------------------------------------
# 8. Keep only aridity-defined dryland classes
# ------------------------------------------------------------

# This keeps:
# - hyperarid
# - arid
# - semiarid
# - dry subhumid
#
# This removes:
# - CBD-only presumed included polygons with P/PET >= 0.65
#   where no aridity class is assigned.
drylands_US_lower48_UNCCD <- drylands_US_lower48_clean %>%
  filter(!is.na(dryland_class))

# Confirm final classes
drylands_US_lower48_UNCCD %>%
  st_drop_geometry() %>%
  count(dryland_class, sort = TRUE)

# ------------------------------------------------------------
# 9. Save clipped shapefiles
# ------------------------------------------------------------

# Save raw lower-48 clip
st_write(
  drylands_US_lower48_raw,
  out_path_raw,
  delete_dsn = TRUE
)

# Save cleaned UNCCD-only lower-48 drylands
st_write(
  drylands_US_lower48_UNCCD,
  out_path_clean,
  delete_dsn = TRUE
)

# ------------------------------------------------------------
# 10. Quick base R plot
# ------------------------------------------------------------

plot(st_geometry(lower48_boundary), border = "grey70")
plot(st_geometry(drylands_US_lower48_UNCCD), add = TRUE, col = "orange")

# Save orange base R plot
orange_map_out_path <- file.path(map_out_dir, "drylands_US_lower48_orange_map.png")

png(
  filename = orange_map_out_path,
  width = 10,
  height = 6,
  units = "in",
  res = 300
)

plot(st_geometry(lower48_boundary), border = "grey70")
plot(st_geometry(drylands_US_lower48_UNCCD), add = TRUE, col = "orange")

dev.off()

# ------------------------------------------------------------
# 11. Class-colored ggplot map
# ------------------------------------------------------------

drylands_map <- ggplot() +
  geom_sf(data = lower48_boundary, fill = "white", color = "grey70") +
  geom_sf(
    data = drylands_US_lower48_UNCCD,
    aes(fill = dryland_class),
    color = "black",
    linewidth = 0.1
  ) +
  theme_minimal() +
  labs(
    title = "UNCCD dryland classes in the contiguous United States",
    fill = "Dryland class"
  )

drylands_map

# ------------------------------------------------------------
# 12. Save map
# ------------------------------------------------------------

ggsave(
  filename = map_out_path,
  plot = drylands_map,
  width = 10,
  height = 6,
  dpi = 300
)