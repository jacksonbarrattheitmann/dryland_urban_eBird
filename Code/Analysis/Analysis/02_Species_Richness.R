# ============================================================
# Dryland Bird Project
# 02_Species_Richness.R
#
# Purpose:
# Model checklist-level bird species richness as a function of
# aridity, local habitat, and checklist effort.
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(mgcv)
library(glmmTMB) # Added bc of the BAM issues we had at part 4 (see my notes)
library(DHARMa)
library(performance)

# ------------------------------------------------------------
# 2. Load and prepare analysis data
# ------------------------------------------------------------

sr_data <- readRDS(
  "Intermediate_data/dfs_for_analysis_combined/combined_data_filt.RDS"
)

sr_data <- sr_data %>%
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
# Again we filter effort_distance_km to <= 5


# ------------------------------------------------------------
# 3. Prepare grouping variables
# ------------------------------------------------------------

sr_data <- sr_data %>%
  mutate(
    NAME10 = factor(NAME10),
    locality_id = factor(locality_id),
    protocol_name = factor(protocol_name)
  )

# ------------------------------------------------------------
# 4. Fit initial Poisson richness model
# ------------------------------------------------------------

#sr_model_poisson <- bam(
#  sr ~
#    avg_aridity_index * water +
#    avg_aridity_index * shrub_and_scrub +
#    avg_aridity_index * built +
#    s(duration_minutes, k = 5) +
#    s(effort_distance_km, k = 5) +
#    protocol_name +
#    s(NAME10, bs = "re") +
#    s(locality_id, bs = "re"),
#  family = poisson(link = "log"),
#  data = sr_data,
#  method = "fREML",
#  discrete = TRUE
#)


# The locality random-effect smooth is the primary computational
# bottleneck because locality_id contains >50,000 levels. Representing
# that many random-effect coefficients within BAM exceeded the available
# 16 GB vector-memory limit.
#
# Therefore, two alternatives are evaluated:
# 1) retain the GAM/BAM framework with city as the random effect, or
# 2) use glmmTMB, which can represent the high-cardinality locality
#    random intercept more efficiently.


# ------------------------------------------------------------
# 5. Fit alternative Poisson richness models
# ------------------------------------------------------------
# Here im going to test the 2 differnt options to see what works.

# ------------------------------------------------------------
# Option A: GAM with city random effect
# ------------------------------------------------------------

sr_model_A_poisson <- bam(
  sr ~
    avg_aridity_index * water +
    avg_aridity_index * shrub_and_scrub +
    avg_aridity_index * built +
    s(duration_minutes, k = 5) +
    s(effort_distance_km, k = 5) +
    protocol_name +
    s(NAME10, bs = "re"),
  family = poisson(link = "log"),
  data = sr_data,
  method = "fREML",
  discrete = TRUE
)


# ------------------------------------------------------------
# Option B: Mixed model with city and locality random effects
# ------------------------------------------------------------

sr_model_B_poisson <- glmmTMB(
  sr ~
    avg_aridity_index * water +
    avg_aridity_index * shrub_and_scrub +
    avg_aridity_index * built +
    splines::ns(duration_minutes, df = 4) +
    effort_distance_km +
    protocol_name +
    (1 | NAME10) +
    (1 | locality_id),
  family = poisson(link = "log"),
  data = sr_data
)

# ------------------------------------------------------------
# 6. Examine initial Poisson models
# ------------------------------------------------------------

summary(sr_model_A_poisson)

summary(sr_model_B_poisson)


#### Let's just use the performance package to check the results

check_model(sr_model_B_poisson)

# ------------------------------------------------------------
# 7. Examine overdispersion in Poisson models
# ------------------------------------------------------------


# ------------------------------------------------------------
# Option A
# ------------------------------------------------------------

# Pearson dispersion ratio

sum(
  residuals(sr_model_A_poisson, type = "pearson")^2
) / df.residual(sr_model_A_poisson)


# Simulation-based dispersion test

library(DHARMa)

sr_residuals_A <- simulateResiduals(
  sr_model_A_poisson,
  n = 100
)

testDispersion(sr_residuals_A)


# ------------------------------------------------------------
# Option B
# ------------------------------------------------------------

# Pearson dispersion ratio

sum(
  residuals(sr_model_B_poisson, type = "pearson")^2
) / df.residual(sr_model_B_poisson)


# DHARMa simulation was also attempted for Option B, but the
# large dataset and >50,000 locality random effects exceeded
# the available 16 GB vector memory.
#
# sr_residuals_B <- simulateResiduals(
#   sr_model_B_poisson,
#   n = 100
# )
#
# Error:
# vector memory limit of 16.0 Gb reached


# ------------------------------------------------------------
# 8. Fit negative-binomial richness models
# ------------------------------------------------------------


# ------------------------------------------------------------
# Option A: GAM with city random effect
# ------------------------------------------------------------

sr_model_A_nb <- bam(
  sr ~
    avg_aridity_index * water +
    avg_aridity_index * shrub_and_scrub +
    avg_aridity_index * built +
    s(duration_minutes, k = 5) +
    s(effort_distance_km, k = 5) +
    protocol_name +
    s(NAME10, bs = "re"),
  family = nb(link = "log"),
  data = sr_data,
  method = "fREML",
  discrete = TRUE
)


# ------------------------------------------------------------
# Option B: Mixed model with city and locality random effects
# ------------------------------------------------------------

sr_model_B_nb <- glmmTMB(
  sr ~
    avg_aridity_index * water +
    avg_aridity_index * shrub_and_scrub +
    avg_aridity_index * built +
    splines::ns(duration_minutes, df = 4) +
    effort_distance_km +
    protocol_name +
    (1 | NAME10) +
    (1 | locality_id),
  family = nbinom2(link = "log"),
  data = sr_data
)


# ------------------------------------------------------------
# 9. Examine negative-binomial models
# ------------------------------------------------------------


# ------------------------------------------------------------
# Model summaries
# ------------------------------------------------------------

summary(sr_model_A_nb)

summary(sr_model_B_nb)


# ------------------------------------------------------------
# Pearson dispersion ratios
# ------------------------------------------------------------

sum(
  residuals(sr_model_A_nb, type = "pearson")^2
) / df.residual(sr_model_A_nb)


sum(
  residuals(sr_model_B_nb, type = "pearson")^2
) / df.residual(sr_model_B_nb)



# ------------------------------------------------------------
# 10. Simulation-based check of NB dispersion
# ------------------------------------------------------------

sr_residuals_A_nb <- simulateResiduals(
  sr_model_A_nb,
  n = 100
)

testDispersion(sr_residuals_A_nb)


# Poisson clearly failed; negative binomial provides a substantially better variance structure and is appropriate for the richness analysis.
# DHARMa indicates slight underdispersion (0.96). Although the
# test is significant, the magnitude of departure from 1 is small
# and significance is strongly influenced by the very large sample size.

# ------------------------------------------------------------
# 11. Check GAM basis dimensions
# ------------------------------------------------------------

gam.check(sr_model_A_nb)

# The GAM converged successfully and showed no major issues with model rank
# or optimization. The basis-dimension check indicated that k = 5 was
# sufficient for effort_distance_km, but duration_minutes reached the
# maximum available EDF and had a low k-index (0.94, p = 0.005).
# Therefore, the duration smooth will be refit with a larger basis dimension
# while keeping the distance smooth unchanged.


# ------------------------------------------------------------
# 12. Refit Option A with increased duration basis
# ------------------------------------------------------------

sr_model_A_nb2 <- bam(
  sr ~
    avg_aridity_index * water +
    avg_aridity_index * shrub_and_scrub +
    avg_aridity_index * built +
    s(duration_minutes, k = 10) +
    s(effort_distance_km, k = 5) +
    protocol_name +
    s(NAME10, bs = "re"),
  family = nb(link = "log"),
  data = sr_data,
  method = "fREML",
  discrete = TRUE
)


# Check basis dimensions again

gam.check(sr_model_A_nb2)

# Increasing the duration basis dimension from k = 5 to k = 10
# resolved the basis-dimension issue (k-index = 0.97, p = 0.26).
# The distance smooth also showed no evidence of an insufficient
# basis dimension (k-index = 0.97, p = 0.12). Therefore, no
# further increase in smooth complexity is needed.



# ------------------------------------------------------------
# 13. Center predictors and fit final species richness model
# ------------------------------------------------------------

# Center focal predictors so that main effects are evaluated at
# average conditions rather than at zero.

sr_data <- sr_data %>%
  mutate(
    aridity_c = avg_aridity_index - mean(avg_aridity_index),
    water_c = water - mean(water),
    shrub_c = shrub_and_scrub - mean(shrub_and_scrub),
    built_c = built - mean(built)
  )


# Fit final negative-binomial mixed model

sr_model_final <- glmmTMB(
  sr ~
    aridity_c * water_c +
    aridity_c * shrub_c +
    aridity_c * built_c +
    splines::ns(duration_minutes, df = 4) +
    effort_distance_km +
    protocol_name +
    (1 | NAME10) +
    (1 | locality_id),
  family = nbinom2(link = "log"),
  data = sr_data)

# Examine final model

summary(sr_model_final)





# ------------------------------------------------------------
# 14. Prepare data for species richness figures
# ------------------------------------------------------------

# City-level richness for the city-level aridity predictor

sr_city_plot <- sr_data %>%
  group_by(NAME10) %>%
  summarise(
    avg_aridity_index = first(avg_aridity_index),
    mean_sr = mean(sr),
    .groups = "drop"
  )


# Locality-level richness for locality-level habitat predictors

sr_locality_plot <- sr_data %>%
  group_by(locality_id) %>%
  summarise(
    water = first(water),
    shrub_and_scrub = first(shrub_and_scrub),
    built = first(built),
    mean_sr = mean(sr),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 15. Visualize species richness relationships
# ------------------------------------------------------------


# Aridity

ggplot(
  sr_city_plot,
  aes(
    x = avg_aridity_index,
    y = mean_sr
  )
) +
  geom_point(
    alpha = 0.4,
    size = 1.5
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 1
  ) +
  labs(
    x = "Aridity index (higher = less arid)",
    y = "Mean species richness"
  ) +
  theme_classic()


# Water

ggplot(
  sr_locality_plot,
  aes(
    x = water,
    y = mean_sr
  )
) +
  geom_point(
    alpha = 0.1,
    size = 0.8
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 1
  ) +
  scale_x_continuous(
    labels = scales::percent
  ) +
  labs(
    x = "Water cover",
    y = "Mean species richness"
  ) +
  theme_classic()


# Shrub and scrub

ggplot(
  sr_locality_plot,
  aes(
    x = shrub_and_scrub,
    y = mean_sr
  )
) +
  geom_point(
    alpha = 0.1,
    size = 0.8
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 1
  ) +
  scale_x_continuous(
    labels = scales::percent
  ) +
  labs(
    x = "Shrub and scrub cover",
    y = "Mean species richness"
  ) +
  theme_classic()


# Built

ggplot(
  sr_locality_plot,
  aes(
    x = built,
    y = mean_sr
  )
) +
  geom_point(
    alpha = 0.1,
    size = 0.8
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 1
  ) +
  scale_x_continuous(
    labels = scales::percent
  ) +
  labs(
    x = "Built cover",
    y = "Mean species richness"
  ) +
  theme_classic()


# ------------------------------------------------------------
# 16. Create species richness model results table
# ------------------------------------------------------------

sr_results <- as.data.frame(
  summary(sr_model_final)$coefficients$cond
)

sr_results$term <- rownames(sr_results)

sr_results <- sr_results %>%
  select(
    term,
    Estimate,
    `Std. Error`,
    `z value`,
    `Pr(>|z|)`
  )

sr_results

ggplot(sr_data, aes(x = avg_aridity_index, y = water)) + 
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()
