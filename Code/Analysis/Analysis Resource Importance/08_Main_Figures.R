library(tidyverse)

# Settings

response_order <- c(
  "Species richness",
  "Functional richness",
  "Functional dispersion",
  "Functional evenness",
  "CWM body mass",
  "CWM beak length"
)

main_responses <- response_order[1:4]
trait_responses <- response_order[5:6]
habitat_order <- c("Water", "Shrub and scrub", "Built")

output_root <- fs::path(
  "Outputs", "Analysis_Resource_Importance", "08_Resource_Importance_Synthesis"
)

main_figure_dir <- fs::path(output_root, "Main_Figures")
supp_figure_dir <- fs::path(output_root, "Supplementary_Figures")
table_dir <- fs::path(output_root, "Tables")
data_dir <- fs::path(output_root, "Data")

walk(
  c(main_figure_dir, supp_figure_dir, table_dir, data_dir),
  fs::dir_create
)

# Final response-specific outputs

response_info <- tribble(
  ~response_id,            ~response,                ~folder,                    ~prefix,
  "Species_Richness",      "Species richness",       "02_Species_Richness",      "Species_Richness",
  "Functional_Richness",   "Functional richness",    "03_Functional_Richness",   "Functional_Richness",
  "Functional_Dispersion", "Functional dispersion",  "04_Functional_Dispersion", "Functional_Dispersion",
  "Functional_Evenness",   "Functional evenness",    "05_Functional_Evenness",   "Functional_Evenness",
  "CWM_Body_Mass",         "CWM body mass",          "06_CWM_Body_Mass",         "CWM_Body_Mass",
  "CWM_Beak_Length",       "CWM beak length",        "07_CWM_Beak_Length",       "CWM_Beak_Length"
) %>%
  mutate(
    root = map_chr(
      folder,
      ~ fs::path("Outputs", "Analysis_Resource_Importance", .x)
    ),
    data_file = map2_chr(
      root, prefix,
      ~ fs::path(.x, "Data", str_c(.y, "_Resource_Importance_Data.rds"))
    ),
    interaction_file = map2_chr(
      root, prefix,
      ~ fs::path(.x, "Table", str_c(.y, "_Habitat_x_Aridity.csv"))
    ),
    curve_file = map2_chr(
      root, prefix,
      ~ fs::path(.x, "Table", str_c(.y, "_Effect_Curves.csv"))
    ),
    sample_file = map2_chr(
      root, prefix,
      ~ fs::path(.x, "Table", str_c(.y, "_Sample_Summary.csv"))
    ),
    city_effect_file = if_else(
      response %in% main_responses,
      map2_chr(
        root, prefix,
        ~ fs::path(.x, "Table", str_c(.y, "_City_Conditional_Effects.csv"))
      ),
      NA_character_
    )
  )

required_files <- c(
  response_info$data_file,
  response_info$interaction_file,
  response_info$curve_file,
  response_info$sample_file,
  response_info$city_effect_file[!is.na(response_info$city_effect_file)]
)

missing_files <- required_files[!fs::file_exists(required_files)]

if (length(missing_files) > 0) {
  abort(str_c("Missing synthesis input(s):\n", str_c(missing_files, collapse = "\n")))
}

# Helpers

pick_column <- function(data, choices) {
  found <- intersect(choices, names(data))
  
  if (length(found) != 1) {
    abort(str_c("Could not uniquely identify one of: ", str_c(choices, collapse = ", ")))
  }
  
  found
}

normalize_curve <- function(data, response_id, response) {
  effect_col <- pick_column(
    data,
    c(
      "relative_change_percent",
      "relative_change_percent_at_typical_baseline"
    )
  )
  
  lower_col <- pick_column(
    data,
    c(
      "relative_CI_lower_percent",
      "relative_CI_lower_percent_at_typical_baseline"
    )
  )
  
  upper_col <- pick_column(
    data,
    c(
      "relative_CI_upper_percent",
      "relative_CI_upper_percent_at_typical_baseline"
    )
  )
  
  data %>%
    transmute(
      response_id,
      response,
      habitat = as.character(habitat),
      avg_aridity_index,
      relative_change_percent = .data[[effect_col]],
      CI_lower_percent = .data[[lower_col]],
      CI_upper_percent = .data[[upper_col]]
    )
}

interpolate_curve <- function(x, y, at) {
  approx(x, y, xout = at, rule = 2)$y
}

# Sample and aridity reference

sample_parts <- response_info %>%
  select(response_id, response, data_file, sample_file) %>%
  pmap(function(response_id, response, data_file, sample_file) {
    sample_summary <- read_csv(sample_file, show_col_types = FALSE)
    
    city_aridity <- read_rds(data_file) %>%
      distinct(NAME10, avg_aridity_index)
    
    if (n_distinct(city_aridity$NAME10) != nrow(city_aridity)) {
      abort(str_c("Aridity is not unique by city for ", response, "."))
    }
    
    aridity_quantiles <- quantile(
      city_aridity$avg_aridity_index,
      c(0.10, 0.50, 0.90),
      names = FALSE
    )
    
    list(
      summary = tibble(
        response_id,
        response,
        n_checklists = sample_summary$n_checklists[[1]],
        n_localities = sample_summary$n_localities[[1]],
        n_cities = sample_summary$n_cities[[1]],
        mean_city_aridity = mean(city_aridity$avg_aridity_index),
        min_city_aridity = min(city_aridity$avg_aridity_index),
        aridity_q10 = aridity_quantiles[[1]],
        aridity_q50 = aridity_quantiles[[2]],
        aridity_q90 = aridity_quantiles[[3]],
        max_city_aridity = max(city_aridity$avg_aridity_index)
      ),
      aridity = city_aridity %>%
        transmute(response_id, response, NAME10, avg_aridity_index)
    )
  })

sample_reference <- map_dfr(sample_parts, "summary")
aridity_support <- map_dfr(sample_parts, "aridity")

rm(sample_parts)

# Final interaction results

interaction_results <- response_info %>%
  select(response_id, response, interaction_file) %>%
  pmap_dfr(function(response_id, response, interaction_file) {
    read_csv(interaction_file, show_col_types = FALSE) %>%
      transmute(
        response_id,
        response,
        habitat = as.character(habitat),
        Estimate = Estimate / 10,
        SE = SE / 10,
        z,
        `P-value` = p,
        CI_lower = CI_lower / 10,
        CI_upper = CI_upper / 10
      )
  }) %>%
  mutate(
    response = factor(response, levels = response_order),
    habitat = factor(habitat, levels = habitat_order)
  ) %>%
  arrange(response, habitat)

# Final response-scale curves

primary_curves <- response_info %>%
  select(response_id, response, curve_file) %>%
  pmap_dfr(function(response_id, response, curve_file) {
    read_csv(curve_file, show_col_types = FALSE) %>%
      normalize_curve(response_id, response)
  }) %>%
  mutate(
    response = factor(response, levels = response_order),
    habitat = factor(habitat, levels = habitat_order)
  )

if (nrow(distinct(primary_curves, response, habitat)) != 18) {
  abort("Final response curves are incomplete.")
}

# City-specific conditional effects for the four biodiversity responses

city_effect_points <- response_info %>%
  filter(response %in% main_responses) %>%
  select(response_id, response, city_effect_file) %>%
  pmap_dfr(function(response_id, response, city_effect_file) {
    read_csv(city_effect_file, show_col_types = FALSE) %>%
      transmute(
        response_id,
        response,
        habitat = as.character(habitat),
        NAME10 = as.character(NAME10),
        avg_aridity_index,
        relative_change_percent
      )
  }) %>%
  mutate(
    response = factor(response, levels = main_responses),
    habitat = factor(habitat, levels = habitat_order)
  )

# Effects at mean and representative aridity values

mean_effects <- interaction_results %>%
  select(response_id, response, habitat) %>%
  left_join(
    sample_reference %>% select(response_id, mean_city_aridity),
    by = "response_id"
  ) %>%
  pmap_dfr(function(response_id, response, habitat, mean_city_aridity) {
    curve <- primary_curves %>%
      filter(
        .data$response_id == .env$response_id,
        as.character(.data$habitat) == as.character(.env$habitat)
      ) %>%
      arrange(avg_aridity_index)
    
    tibble(
      response_id,
      response = as.character(response),
      habitat = as.character(habitat),
      effect_at_mean_aridity_percent = interpolate_curve(
        curve$avg_aridity_index,
        curve$relative_change_percent,
        mean_city_aridity
      ),
      effect_at_mean_aridity_CI_lower = interpolate_curve(
        curve$avg_aridity_index,
        curve$CI_lower_percent,
        mean_city_aridity
      ),
      effect_at_mean_aridity_CI_upper = interpolate_curve(
        curve$avg_aridity_index,
        curve$CI_upper_percent,
        mean_city_aridity
      )
    )
  })

aridity_quantiles <- sample_reference %>%
  select(response_id, aridity_q10, aridity_q50, aridity_q90) %>%
  pivot_longer(
    starts_with("aridity_q"),
    names_to = "quantile",
    values_to = "aridity_value"
  ) %>%
  mutate(
    aridity_position = recode(
      quantile,
      aridity_q10 = "More arid (10th percentile)",
      aridity_q50 = "Median aridity (50th percentile)",
      aridity_q90 = "Less arid (90th percentile)"
    )
  )

quantile_effects <- interaction_results %>%
  select(response_id, response, habitat) %>%
  left_join(aridity_quantiles, by = "response_id") %>%
  pmap_dfr(function(response_id, response, habitat, quantile,
                    aridity_value, aridity_position) {
    curve <- primary_curves %>%
      filter(
        .data$response_id == .env$response_id,
        as.character(.data$habitat) == as.character(.env$habitat)
      ) %>%
      arrange(avg_aridity_index)
    
    tibble(
      response_id,
      response = as.character(response),
      habitat = as.character(habitat),
      aridity_position,
      aridity_value,
      expected_relative_change_percent = interpolate_curve(
        curve$avg_aridity_index,
        curve$relative_change_percent,
        aridity_value
      ),
      CI_lower_percent = interpolate_curve(
        curve$avg_aridity_index,
        curve$CI_lower_percent,
        aridity_value
      ),
      CI_upper_percent = interpolate_curve(
        curve$avg_aridity_index,
        curve$CI_upper_percent,
        aridity_value
      )
    )
  })

# Main figures

main_curves <- primary_curves %>%
  filter(as.character(response) %in% main_responses) %>%
  mutate(
    response = factor(as.character(response), levels = main_responses),
    habitat = factor(as.character(habitat), levels = habitat_order)
  )

figure_01 <- ggplot(
  main_curves,
  aes(avg_aridity_index, relative_change_percent)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_ribbon(
    aes(ymin = CI_lower_percent, ymax = CI_upper_percent),
    alpha = 0.15
  ) +
  geom_point(
    data = city_effect_points,
    aes(avg_aridity_index, relative_change_percent),
    inherit.aes = FALSE,
    size = 0.85,
    alpha = 0.35
  ) +
  geom_line(linewidth = 0.9, lineend = "round") +
  facet_grid(
    rows = vars(response),
    cols = vars(habitat),
    scales = "free_y",
    switch = "y"
  ) +
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

figure_02_data <- interaction_results %>%
  filter(as.character(response) %in% main_responses) %>%
  mutate(
    response = factor(as.character(response), levels = main_responses),
    habitat = factor(as.character(habitat), levels = rev(habitat_order))
  )

figure_02 <- ggplot(
  figure_02_data,
  aes(Estimate, habitat)
) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45) +
  geom_errorbar(
    aes(xmin = CI_lower, xmax = CI_upper),
    orientation = "y",
    width = 0,
    linewidth = 0.75
  ) +
  geom_point(size = 2.8) +
  facet_wrap(vars(response), ncol = 2, scales = "free_x") +
  labs(
    x = "Habitat x aridity interaction estimate per +1 percentage point (95% CI)",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 10),
    axis.title.x = element_text(size = 11),
    panel.spacing = grid::unit(1.1, "lines")
  )

# Supplementary CWM figures

make_trait_figure <- function(response_name, y_label) {
  curve <- primary_curves %>%
    filter(as.character(response) == response_name) %>%
    mutate(habitat = factor(as.character(habitat), levels = habitat_order))
  
  support <- aridity_support %>%
    filter(response == response_name) %>%
    distinct(NAME10, avg_aridity_index)
  
  ggplot(curve, aes(avg_aridity_index, relative_change_percent)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.45) +
    geom_ribbon(
      aes(ymin = CI_lower_percent, ymax = CI_upper_percent),
      alpha = 0.15
    ) +
    geom_line(linewidth = 1.05, lineend = "round") +
    geom_rug(
      data = support,
      aes(x = avg_aridity_index),
      inherit.aes = FALSE,
      sides = "b",
      alpha = 0.12,
      linewidth = 0.2
    ) +
    facet_wrap(vars(habitat), nrow = 1) +
    scale_x_reverse() +
    labs(
      x = "De Martonne aridity index",
      y = y_label,
      caption = "Less arid  <-                                   ->  More arid"
    ) +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 10.5),
      axis.title = element_text(size = 12),
      plot.caption = element_text(hjust = 0.5, size = 9.5)
    )
}

figure_s1 <- make_trait_figure(
  "CWM body mass",
  "Expected change in CWM body mass (%)"
)

figure_s2 <- make_trait_figure(
  "CWM beak length",
  "Expected change in CWM beak length (%)"
)

# Save figures

ggsave(
  fs::path(main_figure_dir, "Figure_01_Combined_Biodiversity_Resource_Importance.png"),
  figure_01,
  width = 12.5,
  height = 12,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  fs::path(main_figure_dir, "Figure_01_Combined_Biodiversity_Resource_Importance.pdf"),
  figure_01,
  width = 12.5,
  height = 12,
  units = "in",
  bg = "white"
)

ggsave(
  fs::path(main_figure_dir, "Figure_02_Habitat_x_Aridity_Effect_Sizes.png"),
  figure_02,
  width = 9,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  fs::path(main_figure_dir, "Figure_02_Habitat_x_Aridity_Effect_Sizes.pdf"),
  figure_02,
  width = 9,
  height = 7,
  units = "in",
  bg = "white"
)

ggsave(
  fs::path(supp_figure_dir, "Figure_S1_Expected_Change_in_CWM_Body_Mass.png"),
  figure_s1,
  width = 14,
  height = 5.4,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  fs::path(supp_figure_dir, "Figure_S1_Expected_Change_in_CWM_Body_Mass.pdf"),
  figure_s1,
  width = 14,
  height = 5.4,
  units = "in",
  bg = "white"
)

ggsave(
  fs::path(supp_figure_dir, "Figure_S2_Expected_Change_in_CWM_Beak_Length.png"),
  figure_s2,
  width = 14,
  height = 5.4,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  fs::path(supp_figure_dir, "Figure_S2_Expected_Change_in_CWM_Beak_Length.pdf"),
  figure_s2,
  width = 14,
  height = 5.4,
  units = "in",
  bg = "white"
)

# Tables

publication_stats <- interaction_results %>%
  transmute(
    Response = as.character(response),
    Habitat = as.character(habitat),
    Estimate,
    SE,
    z,
    `P-value`
  )

main_stats <- publication_stats %>%
  filter(Response %in% main_responses) %>%
  mutate(
    Response = factor(Response, levels = main_responses),
    Habitat = factor(Habitat, levels = habitat_order)
  ) %>%
  arrange(Response, Habitat) %>%
  mutate(
    Response = as.character(Response),
    Habitat = as.character(Habitat)
  )

trait_stats <- publication_stats %>%
  filter(Response %in% trait_responses) %>%
  mutate(
    Response = factor(Response, levels = trait_responses),
    Habitat = factor(Habitat, levels = habitat_order)
  ) %>%
  arrange(Response, Habitat) %>%
  mutate(
    Response = as.character(Response),
    Habitat = as.character(Habitat)
  )

full_stats <- interaction_results %>%
  left_join(
    sample_reference %>%
      select(response_id, n_checklists, n_cities, n_localities),
    by = "response_id"
  ) %>%
  left_join(
    mean_effects,
    by = c("response_id", "response", "habitat")
  ) %>%
  transmute(
    Response = as.character(response),
    Habitat = as.character(habitat),
    Checklists = n_checklists,
    Cities = n_cities,
    Localities = n_localities,
    Effect_at_mean_aridity_percent = effect_at_mean_aridity_percent,
    Effect_at_mean_aridity_CI_lower_percent = effect_at_mean_aridity_CI_lower,
    Effect_at_mean_aridity_CI_upper_percent = effect_at_mean_aridity_CI_upper,
    Habitat_x_Aridity_Estimate_per_1pp = Estimate,
    Habitat_x_Aridity_SE_per_1pp = SE,
    Habitat_x_Aridity_CI_lower_per_1pp = CI_lower,
    Habitat_x_Aridity_CI_upper_per_1pp = CI_upper,
    z,
    `P-value`
  )

main_quantile_effects <- quantile_effects %>%
  filter(response %in% main_responses)

trait_quantile_effects <- quantile_effects %>%
  filter(response %in% trait_responses)

write_csv(
  main_stats,
  fs::path(table_dir, "Table_S1_Biodiversity_Habitat_x_Aridity_Statistics.csv")
)

write_csv(
  trait_stats,
  fs::path(table_dir, "Table_S2_CWM_Trait_Habitat_x_Aridity_Statistics.csv")
)

write_csv(
  main_quantile_effects,
  fs::path(table_dir, "Table_S3_Biodiversity_Effects_at_Aridity_Quantiles.csv")
)

write_csv(
  trait_quantile_effects,
  fs::path(table_dir, "Table_S4_CWM_Trait_Effects_at_Aridity_Quantiles.csv")
)

write_csv(
  full_stats,
  fs::path(table_dir, "Table_01_Full_Habitat_x_Aridity_Statistics.csv")
)

write_csv(
  sample_reference,
  fs::path(table_dir, "Table_02_Response_Sample_and_Aridity_Reference.csv")
)

# Figure data

write_csv(
  primary_curves,
  fs::path(data_dir, "Final_Primary_Response_Effect_Curves.csv")
)

write_csv(
  city_effect_points,
  fs::path(data_dir, "Figure_01_City_Specific_Conditional_Effect_Points.csv")
)
