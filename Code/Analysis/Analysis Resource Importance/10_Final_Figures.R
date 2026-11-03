library(dplyr)
library(ggplot2)
library(readr)

# Settings

figure_dpi <- 600
response_order <- c(
  "Species richness",
  "Functional richness",
  "Functional dispersion",
  "Functional evenness",
  "CWM body mass",
  "CWM beak length"
)

main_responses <- response_order[1:4]
functional_responses <- response_order[2:4]
trait_responses <- response_order[5:6]
habitat_order <- c("Water", "Shrub and scrub", "Built")

analysis_root <- fs::path("Outputs", "Analysis_Resource_Importance")
synthesis_data_dir <- fs::path(
  analysis_root, "08_Resource_Importance_Synthesis", "Data"
)
prediction_root <- fs::path(
  analysis_root, "09_Resource_Importance_Predictions"
)
output_root <- fs::path(analysis_root, "10_Final_Figures")
main_figure_dir <- fs::path(output_root, "Main_Figures")
supp_figure_dir <- fs::path(output_root, "Supplementary_Figures")

fs::dir_create(main_figure_dir)
fs::dir_create(supp_figure_dir)

curve_file <- fs::path(
  synthesis_data_dir, "Final_Primary_Response_Effect_Curves.csv"
)
city_effect_file <- fs::path(
  synthesis_data_dir, "Figure_01_City_Specific_Conditional_Effect_Points.csv"
)
water_prediction_file <- fs::path(
  prediction_root, "Data", "Figure_S3_City_Specific_Water_Predictions.csv"
)
water_summary_file <- fs::path(
  prediction_root, "Table", "Table_05_Water_Prediction_Summary.csv"
)

input_files <- c(
  curve_file, city_effect_file, water_prediction_file, water_summary_file
)
missing_files <- input_files[!fs::file_exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing saved figure input(s). Run Scripts 08 and 09 first:\n",
    paste(missing_files, collapse = "\n")
  )
}

# Load saved figure data

primary_curves <- read_csv(curve_file, show_col_types = FALSE)
city_effect_points <- read_csv(city_effect_file, show_col_types = FALSE)
water_prediction_data <- read_csv(
  water_prediction_file, show_col_types = FALSE
)
water_prediction_summary <- read_csv(
  water_summary_file, show_col_types = FALSE
)

check_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop("Missing columns in ", label, ": ", paste(missing, collapse = ", "))
  }
}

check_columns(
  primary_curves,
  c(
    "response", "habitat", "avg_aridity_index", "relative_change_percent",
    "CI_lower_percent", "CI_upper_percent"
  ),
  "effect curves"
)
check_columns(
  city_effect_points,
  c("response", "habitat", "NAME10", "avg_aridity_index", "relative_change_percent"),
  "city effects"
)
check_columns(
  water_prediction_data,
  c(
    "Response", "Habitat", "NAME10", "avg_aridity_index",
    "habitat_deviation_pp", "predicted_relative_change_percent"
  ),
  "Water predictions"
)
check_columns(
  water_prediction_summary,
  c("Response", "Habitat", "retained_cities"),
  "Water prediction summary"
)

# Input checks

expected_pairs <- expand.grid(
  response = response_order,
  habitat = habitat_order,
  stringsAsFactors = FALSE
)
observed_pairs <- primary_curves %>%
  distinct(response, habitat)

pair_key <- function(x) paste(x$response, x$habitat, sep = "::")

if (!setequal(pair_key(expected_pairs), pair_key(observed_pairs)) ||
    nrow(observed_pairs) != nrow(expected_pairs)) {
  stop("Effect curves must contain all 18 response-habitat combinations.")
}

if (any(!is.finite(as.matrix(primary_curves[, c(
  "avg_aridity_index", "relative_change_percent",
  "CI_lower_percent", "CI_upper_percent"
)]))) ||
any(primary_curves$CI_lower_percent > primary_curves$CI_upper_percent)) {
  stop("Effect curves contain invalid values or confidence intervals.")
}

curve_counts <- primary_curves %>%
  group_by(response, habitat) %>%
  summarise(n_aridity = n_distinct(avg_aridity_index), .groups = "drop")

if (any(curve_counts$n_aridity < 2)) {
  stop("At least two aridity values are required per effect curve.")
}

if (any(!city_effect_points$response %in% main_responses) ||
    any(!city_effect_points$habitat %in% habitat_order) ||
    any(!is.finite(city_effect_points$avg_aridity_index)) ||
    any(!is.finite(city_effect_points$relative_change_percent))) {
  stop("City-effect points contain unexpected responses or invalid values.")
}

water_prediction_data <- water_prediction_data %>%
  filter(Habitat == "Water")
water_prediction_summary <- water_prediction_summary %>%
  filter(Habitat == "Water")

if (nrow(water_prediction_data) == 0 ||
    nrow(water_prediction_summary) != length(main_responses) ||
    !setequal(unique(water_prediction_data$Response), main_responses) ||
    !setequal(water_prediction_summary$Response, main_responses) ||
    anyDuplicated(water_prediction_summary$Response) > 0 ||
    anyNA(water_prediction_data[, c(
      "Response", "NAME10", "habitat_deviation_pp",
      "predicted_relative_change_percent", "avg_aridity_index"
    )]) ||
    any(!is.finite(water_prediction_data$habitat_deviation_pp)) ||
    any(!is.finite(water_prediction_data$predicted_relative_change_percent)) ||
    any(!is.finite(water_prediction_data$avg_aridity_index))) {
  stop("Water predictions must contain four complete biodiversity response panels.")
}

observed_city_counts <- water_prediction_data %>%
  distinct(Response, NAME10) %>%
  count(Response, name = "n_cities_in_csv") %>%
  left_join(
    water_prediction_summary %>% select(Response, retained_cities),
    by = "Response", relationship = "one-to-one"
  )

if (anyNA(observed_city_counts$retained_cities) ||
    any(observed_city_counts$n_cities_in_csv !=
        observed_city_counts$retained_cities)) {
  stop("Figure S3 city counts do not match Script 09's prediction summary.")
}

baseline_check <- water_prediction_data %>%
  group_by(Response, NAME10) %>%
  summarise(
    has_zero_baseline = any(
      abs(habitat_deviation_pp) < 1e-8 &
        abs(predicted_relative_change_percent) < 1e-8
    ),
    .groups = "drop"
  )

if (any(!baseline_check$has_zero_baseline)) {
  stop("Figure S3 is missing a normalized zero baseline for at least one city.")
}

# Figures 1, 2, S1, and S2

make_effect_figure <- function(responses, show_city_points = FALSE) {
  plot_curves <- primary_curves %>%
    filter(response %in% responses) %>%
    mutate(
      response = factor(response, levels = responses),
      habitat = factor(habitat, levels = habitat_order)
    )
  
  plot_figure <- ggplot(
    plot_curves,
    aes(avg_aridity_index, relative_change_percent)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
    geom_ribbon(
      aes(ymin = CI_lower_percent, ymax = CI_upper_percent),
      alpha = 0.15
    )
  
  if (show_city_points) {
    plot_points <- city_effect_points %>%
      filter(response %in% responses) %>%
      mutate(
        response = factor(response, levels = responses),
        habitat = factor(habitat, levels = habitat_order)
      )
    
    plot_figure <- plot_figure +
      geom_point(
        data = plot_points,
        aes(avg_aridity_index, relative_change_percent),
        inherit.aes = FALSE,
        size = 0.85,
        alpha = 0.35
      )
  }
  
  plot_figure <- plot_figure +
    geom_line(linewidth = 0.9, lineend = "round")
  
  if (length(responses) == 1) {
    plot_figure <- plot_figure +
      facet_grid(cols = vars(habitat))
  } else {
    plot_figure <- plot_figure +
      facet_grid(
        rows = vars(response),
        cols = vars(habitat),
        scales = "free_y",
        switch = "y"
      )
  }
  
  plot_figure +
    scale_x_reverse() +
    labs(
      x = "De Martonne aridity index",
      y = "Expected change in response (%)",
      caption = "Less arid  <-                                   ->  More arid"
    ) +
    theme_classic(base_size = 11) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_text(face = "bold", size = 11.5),
      strip.text.y.left = element_text(face = "bold", angle = 0, size = 10.5),
      strip.placement = "outside",
      panel.spacing = grid::unit(0.8, "lines"),
      axis.text = element_text(size = 9.5),
      axis.title = element_text(size = 11),
      plot.caption = element_text(hjust = 0.5, size = 9.5)
    )
}

figure_01 <- make_effect_figure(response_order[1], show_city_points = TRUE)
figure_02 <- make_effect_figure(functional_responses, show_city_points = TRUE)
figure_s1 <- make_effect_figure(trait_responses[1])
figure_s2 <- make_effect_figure(trait_responses[2])

# Figure S3

water_prediction_data <- water_prediction_data %>%
  mutate(Response = factor(Response, levels = main_responses))

figure_s3 <- ggplot(
  water_prediction_data,
  aes(
    x = habitat_deviation_pp,
    y = predicted_relative_change_percent,
    group = NAME10,
    colour = avg_aridity_index
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.4) +
  geom_line(alpha = 0.35, linewidth = 0.5, lineend = "round") +
  scale_colour_viridis_c(
    option = "plasma",
    direction = -1,
    name = "De Martonne\naridity index\n(lower = more arid)"
  ) +
  scale_x_continuous(
    limits = c(-0.5, 0.5),
    breaks = seq(-0.5, 0.5, length.out = 5),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  facet_wrap(vars(Response), ncol = 2, scales = "free_y") +
  labs(
    x = "Water cover relative to city mean (percentage points)",
    y = "Predicted change relative to city-mean water cover (%)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9.5),
    panel.spacing = grid::unit(0.9, "lines")
  )

# Save manuscript figures

figure_outputs <- data.frame(
  figure_id = c("Figure 1", "Figure 2", "Figure S1", "Figure S2", "Figure S3"),
  directory = c(
    main_figure_dir, main_figure_dir,
    supp_figure_dir, supp_figure_dir, supp_figure_dir
  ),
  filename = c(
    "Figure_01_Species_Richness_Resource_Importance",
    "Figure_02_Functional_Diversity_Resource_Importance",
    "Figure_S1_CWM_Body_Mass_Resource_Importance",
    "Figure_S2_CWM_Beak_Length_Resource_Importance",
    "Figure_S3_City_Specific_Water_Predictions"
  ),
  width_in = c(12.5, 12.5, 12.5, 12.5, 12.5),
  height_in = c(4.5, 10, 4.5, 4.5, 8.5),
  stringsAsFactors = FALSE
)

plots <- list(figure_01, figure_02, figure_s1, figure_s2, figure_s3)

for (i in seq_along(plots)) {
  path_png <- fs::path(
    figure_outputs$directory[i], paste0(figure_outputs$filename[i], ".png")
  )
  path_pdf <- fs::path(
    figure_outputs$directory[i], paste0(figure_outputs$filename[i], ".pdf")
  )
  
  ggsave(
    path_png, plots[[i]],
    width = figure_outputs$width_in[i],
    height = figure_outputs$height_in[i],
    units = "in", dpi = figure_dpi, bg = "white"
  )
  
  ggsave(
    path_pdf, plots[[i]],
    width = figure_outputs$width_in[i],
    height = figure_outputs$height_in[i],
    units = "in", bg = "white"
  )
  
  message("Saved ", figure_outputs$figure_id[i], ": ", path_png)
}

# Reproducibility manifest

input_manifest <- data.frame(
  input_file = input_files,
  size_bytes = as.numeric(file.info(input_files)$size),
  modified = as.character(file.info(input_files)$mtime),
  stringsAsFactors = FALSE
)

write_csv(
  input_manifest,
  fs::path(output_root, "Final_Figure_Input_Manifest.csv")
)
write_csv(
  figure_outputs,
  fs::path(output_root, "Final_Figure_Output_Manifest.csv")
)

message("Script 10 complete: five figures regenerated from saved CSV outputs.")
