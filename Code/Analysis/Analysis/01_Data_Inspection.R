# ============================================================
# Dryland Bird Project
# 01_Data_Inspection.R
# So all I am doing here is inspecting all the data to make sure everything makes sense,
# that the variables line up, nothing is missing, and such.
# Basically inspecting everything to make sure the structure of our models will make sense.
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(dplyr)


# ------------------------------------------------------------
# 2. Load datasets
# ------------------------------------------------------------

combined_data <- readRDS(
  "Intermediate_data/dfs_for_analysis_combined/combined_data.RDS"
)

combined_data_filt <- readRDS(
  "Intermediate_data/dfs_for_analysis_combined/combined_data_filt.RDS"
)


# ------------------------------------------------------------
# 3. Dataset dimensions
# ------------------------------------------------------------

dim(combined_data)

dim(combined_data_filt)


# ------------------------------------------------------------
# 4. Column names
# ------------------------------------------------------------

names(combined_data)

names(combined_data_filt)


# ------------------------------------------------------------
# 5. Compare columns between datasets
# ------------------------------------------------------------

# Variables only found in combined_data

setdiff(
  names(combined_data),
  names(combined_data_filt)
)


# Variables only found in combined_data_filt

setdiff(
  names(combined_data_filt),
  names(combined_data)
)


# Variables found in both datasets

intersect(
  names(combined_data),
  names(combined_data_filt)
)



# ------------------------------------------------------------
# 6. Compare raw and filtered datasets
# ------------------------------------------------------------

data_comparison <- data.frame(
  dataset = c("combined_data", "combined_data_filt"),
  
  n_checklists = c(
    nrow(combined_data),
    nrow(combined_data_filt)
  ),
  
  n_localities = c(
    n_distinct(combined_data$locality_id),
    n_distinct(combined_data_filt$locality_id)
  ),
  
  n_cities = c(
    n_distinct(combined_data$NAME10),
    n_distinct(combined_data_filt$NAME10)
  ),
  
  n_states = c(
    n_distinct(combined_data$state_name),
    n_distinct(combined_data_filt$state_name)
  )
)

data_comparison



# ------------------------------------------------------------
# 7. Identify removed checklists
# ------------------------------------------------------------

removed_checklists <- combined_data %>%
  anti_join(
    combined_data_filt,
    by = "checklist_id"
  )

dim(removed_checklists)

combined_data %>%
  count(checklist_id) %>%
  filter(n > 1)

combined_data_filt %>%
  count(checklist_id) %>%
  filter(n > 1)



# ------------------------------------------------------------
# 8. Compare retained and removed checklists
# ------------------------------------------------------------

retained_summary <- combined_data_filt %>%
  summarise(
    mean_sr = mean(sr, na.rm = TRUE),
    median_sr = median(sr, na.rm = TRUE),
    mean_duration = mean(duration_minutes, na.rm = TRUE),
    median_duration = median(duration_minutes, na.rm = TRUE),
    mean_distance = mean(effort_distance_km, na.rm = TRUE),
    median_distance = median(effort_distance_km, na.rm = TRUE),
    mean_observers = mean(number_observers, na.rm = TRUE)
  )

removed_summary <- removed_checklists %>%
  summarise(
    mean_sr = mean(sr, na.rm = TRUE),
    median_sr = median(sr, na.rm = TRUE),
    mean_duration = mean(duration_minutes, na.rm = TRUE),
    median_duration = median(duration_minutes, na.rm = TRUE),
    mean_distance = mean(effort_distance_km, na.rm = TRUE),
    median_distance = median(effort_distance_km, na.rm = TRUE),
    mean_observers = mean(number_observers, na.rm = TRUE)
  )

retained_summary
removed_summary


# ------------------------------------------------------------
# 9. Compare checklist protocols
# ------------------------------------------------------------

combined_data_filt %>%
  count(protocol_name, sort = TRUE)

removed_checklists %>%
  count(protocol_name, sort = TRUE)



# ------------------------------------------------------------
# 10. Examine effort variables by protocol
# ------------------------------------------------------------

combined_data_filt %>%
  group_by(protocol_name) %>%
  summarise(
    n = n(),
    missing_duration = sum(is.na(duration_minutes)),
    missing_distance = sum(is.na(effort_distance_km)),
    median_duration = median(duration_minutes, na.rm = TRUE),
    median_distance = median(effort_distance_km, na.rm = TRUE),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 11. Create analysis dataframe
# ------------------------------------------------------------

dryland_data <- combined_data_filt %>%
  mutate(
    effort_distance_km = if_else(
      protocol_name == "Stationary",
      0,
      effort_distance_km
    )
  ) %>%
  filter(
    protocol_name %in% c("Stationary", "Traveling"),
    duration_minutes <= 300,
    effort_distance_km <= 5
  )


# Check analysis dataset

dryland_data %>%
  summarise(
    n_checklists = n(),
    n_localities = n_distinct(locality_id),
    n_cities = n_distinct(NAME10),
    n_states = n_distinct(state_name),
    max_duration = max(duration_minutes),
    max_distance = max(effort_distance_km)
  )

dryland_data %>%
  count(protocol_name, sort = TRUE)

### FILTER DISTANCE TRAVEL TO 5KM 

# ------------------------------------------------------------
# 12. Examine missing data in analysis variables
# ------------------------------------------------------------

model_variables <- c(
  "sr",
  "FRic",
  "FDis",
  "FEve",
  "avg_aridity_index",
  "water",
  "shrub_and_scrub",
  "built",
  "duration_minutes",
  "effort_distance_km",
  "locality_id",
  "NAME10"
)

missing_data <- data.frame(
  variable = model_variables,
  missing_n = sapply(
    dryland_data[model_variables],
    function(x) sum(is.na(x))
  )
)

missing_data$missing_percent <-
  missing_data$missing_n / nrow(dryland_data) * 100

missing_data


# ------------------------------------------------------------
# 13. Examine missing functional diversity values
# ------------------------------------------------------------

dryland_data %>%
  group_by(
    FRic_missing = is.na(FRic)
  ) %>%
  summarise(
    n = n(),
    min_sr = min(sr),
    median_sr = median(sr),
    mean_sr = mean(sr),
    max_sr = max(sr),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 14. Compare FRic and FEve missing values
# ------------------------------------------------------------

table(
  FRic_missing = is.na(dryland_data$FRic),
  FEve_missing = is.na(dryland_data$FEve)
)



# ------------------------------------------------------------
# 15. Examine aridity variation within cities
# ------------------------------------------------------------

dryland_data %>%
  group_by(NAME10) %>%
  summarise(
    n_aridity_values = n_distinct(avg_aridity_index),
    .groups = "drop"
  ) %>%
  count(n_aridity_values)


# ------------------------------------------------------------
# 16. Examine habitat variation within localities
# ------------------------------------------------------------

dryland_data %>%
  group_by(locality_id) %>%
  summarise(
    n_water_values = n_distinct(water),
    n_shrub_values = n_distinct(shrub_and_scrub),
    n_built_values = n_distinct(built),
    .groups = "drop"
  ) %>%
  summarise(
    localities_with_multiple_water_values = sum(n_water_values > 1),
    localities_with_multiple_shrub_values = sum(n_shrub_values > 1),
    localities_with_multiple_built_values = sum(n_built_values > 1)
  )

# ------------------------------------------------------------
# 17. Examine sampling structure within cities and localities
# ------------------------------------------------------------

dryland_data %>%
  group_by(NAME10) %>%
  summarise(
    n_localities = n_distinct(locality_id),
    n_checklists = n(),
    .groups = "drop"
  ) %>%
  summarise(
    min_localities = min(n_localities),
    median_localities = median(n_localities),
    mean_localities = mean(n_localities),
    max_localities = max(n_localities),
    min_checklists = min(n_checklists),
    median_checklists = median(n_checklists),
    mean_checklists = mean(n_checklists),
    max_checklists = max(n_checklists)
  )


dryland_data %>%
  group_by(locality_id) %>%
  summarise(
    n_checklists = n(),
    .groups = "drop"
  ) %>%
  summarise(
    min_checklists = min(n_checklists),
    median_checklists = median(n_checklists),
    mean_checklists = mean(n_checklists),
    max_checklists = max(n_checklists)
  )



# ------------------------------------------------------------
# 18. Examine response variable distributions
# ------------------------------------------------------------

summary(dryland_data$sr)
summary(dryland_data$FRic)
summary(dryland_data$FDis)
summary(dryland_data$FEve)

dryland_data %>%
  summarise(
    zero_sr = sum(sr == 0, na.rm = TRUE),
    zero_FRic = sum(FRic == 0, na.rm = TRUE),
    zero_FDis = sum(FDis == 0, na.rm = TRUE),
    zero_FEve = sum(FEve == 0, na.rm = TRUE)
  )


# ------------------------------------------------------------
# 19. Examine zero functional dispersion values
# ------------------------------------------------------------

dryland_data %>%
  group_by(
    FDis_zero = FDis == 0
  ) %>%
  summarise(
    n = n(),
    min_sr = min(sr),
    median_sr = median(sr),
    mean_sr = mean(sr),
    max_sr = max(sr),
    .groups = "drop"
  )



# ------------------------------------------------------------
# Check dispresion of response variables now to make sure the proposed models
# will work.
# ------------------------------------------------------------

# ------------------------------------------------------------
# 20. Examine species richness distribution
# ------------------------------------------------------------

ggplot(dryland_data, aes(x = sr)) +
  geom_histogram(binwidth = 1) +
  labs(
    x = "Species richness",
    y = "Number of checklists"
  ) +
  theme_classic() # Righ Skewed like i figured it would be for SR

# Negative binomial is already the more plausible candidate for species richness, but we should still verify that using the fitted model diagnostics


# Preliminary check for overdispersion

mean(dryland_data$sr)

var(dryland_data$sr)

var(dryland_data$sr) / mean(dryland_data$sr)



# ------------------------------------------------------------
# 21. Examine functional richness distribution
# ------------------------------------------------------------

ggplot(dryland_data, aes(x = FRic)) +
  geom_histogram(binwidth = 0.02, na.rm = TRUE) +
  labs(
    x = "Functional richness",
    y = "Number of checklists"
  ) +
  theme_classic()

summary(dryland_data$FRic)

# We cannot pre-select Gaussian for this just yet. Bounded and Skewed shape means a BETA GAM might be appropriate. Will test in next script (02)



# ------------------------------------------------------------
# 22. Examine functional dispersion distribution
# ------------------------------------------------------------

ggplot(dryland_data, aes(x = FDis)) +
  geom_histogram(binwidth = 0.01) +
  labs(
    x = "Functional dispersion",
    y = "Number of checklists"
  ) +
  theme_classic()

summary(dryland_data$FDis)

# This one most likely is just a Guassian model.



# ------------------------------------------------------------
# 23. Examine functional evenness distribution
# ------------------------------------------------------------

ggplot(dryland_data, aes(x = FEve)) +
  geom_histogram(binwidth = 0.02, na.rm = TRUE) +
  labs(
    x = "Functional evenness",
    y = "Number of checklists"
  ) +
  theme_classic()

summary(dryland_data$FEve)

# Modst likely Gaussian model.


# ------------------------------------------------------------
# 24. Examine boundary values in functional evenness
# ------------------------------------------------------------

dryland_data %>%
  summarise(
    n_FEve_1 = sum(FEve == 1, na.rm = TRUE),
    percent_FEve_1 = sum(FEve == 1, na.rm = TRUE) /
      sum(!is.na(FEve)) * 100
  )
# No FEve are actually equal to 1


# ------------------------------------------------------------
# 25. Examine total land-cover coverage
# ------------------------------------------------------------

dryland_data %>%
  mutate(
    landcover_total =
      bare +
      built +
      crops +
      flooded_vegetation +
      grass +
      shrub_and_scrub +
      snow_and_ice +
      trees +
      water
  ) %>%
  summarise(
    min_total = min(landcover_total),
    median_total = median(landcover_total),
    mean_total = mean(landcover_total),
    max_total = max(landcover_total)
  )
# Just making sure land cover is all at 100% or 1.00


# ------------------------------------------------------------
# 26. Examine correlations among habitat predictors
# This is checking the Locality Level Correlations
# ------------------------------------------------------------

habitat_data <- dryland_data %>%
  distinct(
    locality_id,
    water,
    shrub_and_scrub,
    built
  )

habitat_cor <- cor(
  habitat_data %>%
    select(
      water,
      shrub_and_scrub,
      built
    )
)

habitat_cor


# Visualize correlation matrix

habitat_cor_df <- as.data.frame(as.table(habitat_cor))

ggplot(habitat_cor_df, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Freq, 2)), size = 4) +
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  coord_equal() +
  labs(
    x = NULL,
    y = NULL,
    fill = "Correlation"
  ) +
  theme_classic()


# ------------------------------------------------------------
# 27. Examine habitat predictor distributions
# ------------------------------------------------------------

ggplot(habitat_data, aes(x = water)) +
  geom_histogram(binwidth = 0.05) +
  labs(
    x = "Proportion water",
    y = "Number of localities"
  ) +
  theme_classic()

# Very skewed.


ggplot(habitat_data, aes(x = shrub_and_scrub)) +
  geom_histogram(binwidth = 0.05) +
  labs(
    x = "Proportion shrub and scrub",
    y = "Number of localities"
  ) +
  theme_classic()

# Also Skewed


ggplot(habitat_data, aes(x = built)) +
  geom_histogram(binwidth = 0.05) +
  labs(
    x = "Proportion built",
    y = "Number of localities"
  ) +
  theme_classic()

# Lots of variation so GOOD!



# ------------------------------------------------------------
# 28. Examine aridity distribution across cities
# ------------------------------------------------------------

dryland_data %>%
  distinct(
    NAME10,
    avg_aridity_index
  ) %>%
  ggplot(aes(x = avg_aridity_index)) +
  geom_histogram(binwidth = 2) +
  labs(
    x = "Aridity index",
    y = "Number of cities"
  ) +
  theme_classic()


dryland_data %>%
  distinct(
    NAME10,
    avg_aridity_index
  ) %>%
  pull(avg_aridity_index) %>%
  summary()

# REMEMBER!!!!! THE LOWER THE VALUE THE MORE ARID!!!!! SO IF THE MODELS LATER SHOW A POSITIVE EFFECT THAT ACTUALLY MEANS LESS ARID


# ------------------------------------------------------------
# 29. Examine relationships between aridity and habitat
# ------------------------------------------------------------

predictor_data <- dryland_data %>%
  distinct(
    locality_id,
    avg_aridity_index,
    water,
    shrub_and_scrub,
    built
  )


# Correlation matrix

predictor_cor <- cor(
  predictor_data %>%
    select(
      avg_aridity_index,
      water,
      shrub_and_scrub,
      built
    )
)

predictor_cor


# Visualize correlation matrix

predictor_cor_df <- as.data.frame(as.table(predictor_cor))

ggplot(predictor_cor_df, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = round(Freq, 2)),
    size = 4
  ) +
  scale_fill_gradient2(
    low = "steelblue",
    mid = "white",
    high = "firebrick",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  scale_x_discrete(
    labels = c(
      avg_aridity_index = "Aridity index",
      water = "Water",
      shrub_and_scrub = "Shrub and scrub",
      built = "Built"
    )
  ) +
  scale_y_discrete(
    labels = c(
      avg_aridity_index = "Aridity index",
      water = "Water",
      shrub_and_scrub = "Shrub and scrub",
      built = "Built"
    )
  ) +
  coord_equal() +
  labs(
    x = NULL,
    y = NULL,
    fill = "Correlation"
  ) +
  theme_classic()

# looks good. I want to next examine the effort covariates since we have looked at all of the habitat and locality ones.



# ------------------------------------------------------------
# 30. Examine checklist duration
# ------------------------------------------------------------

ggplot(dryland_data, aes(x = duration_minutes)) +
  geom_histogram(binwidth = 10) +
  labs(
    x = "Checklist duration (minutes)",
    y = "Number of checklists"
  ) +
  theme_classic()

summary(dryland_data$duration_minutes)

# So from the looks of this, for a GAM it might be s(duration_minutes) INSTEAD of just duration_minutes


# ------------------------------------------------------------
# 31. Examine checklist distance traveled
# ------------------------------------------------------------

ggplot(dryland_data, aes(x = effort_distance_km)) +
  geom_histogram(binwidth = 0.25, na.rm = TRUE) +
  labs(
    x = "Distance traveled (km)",
    y = "Number of checklists"
  ) +
  theme_classic()

summary(dryland_data$effort_distance_km)

# So because there are just so many damn stationary ones we MIGHT need to include protocol into the model to account for this.


# ------------------------------------------------------------
# 32. Examine distance traveled by protocol
# ------------------------------------------------------------

dryland_data %>%
  group_by(protocol_name) %>%
  summarise(
    n = n(),
    median_distance = median(effort_distance_km, na.rm = TRUE),
    mean_distance = mean(effort_distance_km, na.rm = TRUE),
    max_distance = max(effort_distance_km, na.rm = TRUE),
    .groups = "drop"
  )

dryland_data %>%
  filter(protocol_name == "Traveling") %>%
  summarise(
    min_distance = min(effort_distance_km, na.rm = TRUE),
    q25_distance = quantile(effort_distance_km, 0.25, na.rm = TRUE),
    median_distance = median(effort_distance_km, na.rm = TRUE),
    q75_distance = quantile(effort_distance_km, 0.75, na.rm = TRUE),
    q95_distance = quantile(effort_distance_km, 0.95, na.rm = TRUE),
    q99_distance = quantile(effort_distance_km, 0.99, na.rm = TRUE),
    max_distance = max(effort_distance_km, na.rm = TRUE)
  )

# Someone did an 80km checklist? Yeah thats gotta get filtered...


# ------------------------------------------------------------
# 33. Examine long-distance Traveling checklists
# ------------------------------------------------------------

dryland_data %>%
  filter(protocol_name == "Traveling") %>%
  summarise(
    n_total = n(),
    n_over_5km = sum(effort_distance_km > 5, na.rm = TRUE),
    n_over_10km = sum(effort_distance_km > 10, na.rm = TRUE),
    n_over_20km = sum(effort_distance_km > 20, na.rm = TRUE),
    n_over_50km = sum(effort_distance_km > 50, na.rm = TRUE)
  )

# Ok there are a few other checklists that are way too long. From the lit we usually do not use them when they are that long for distance so lets see when we reduce it to 5


# ------------------------------------------------------------
# 34. Confirm locality-city nesting
# ------------------------------------------------------------

dryland_data %>%
  group_by(locality_id) %>%
  summarise(
    n_cities = n_distinct(NAME10),
    .groups = "drop"
  ) %>%
  count(n_cities)



# ------------------------------------------------------------
# 35. Examine number of observers
# ------------------------------------------------------------

table(dryland_data$number_observers)

summary(dryland_data$number_observers)