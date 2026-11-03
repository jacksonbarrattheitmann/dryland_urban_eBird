library(dplyr)
library(ggplot2)
library(purrr)
library(readr)
library(tibble)
library(tidyr)

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
functional_responses <- response_order[2:4]
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

# Response-specific outputs

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
    interaction_file = map2_chr(
      root, prefix,
      ~ fs::path(.x, "Table", paste0(.y, "_Habitat_x_Aridity.csv"))
    ),
    curve_file = map2_chr(
      root, prefix,
      ~ fs::path(.x, "Table", paste0(.y, "_Effect_Curves.csv"))
    ),
    sample_file = map2_chr(
      root, prefix,
      ~ fs::path(.x, "Table", paste0(.y, "_Sample_Summary.csv"))
    ),
    city_effect_file = map2_chr(
      root, prefix,
      ~ fs::path(.x, "Table", paste0(.y, "_City_Conditional_Effects.csv"))
    )
  )

required_files <- c(
  response_info$interaction_file,
  response_info$curve_file,
  response_info$sample_file,
  response_info$city_effect_file
)

missing_files <- required_files[!fs::file_exists(required_files)]

if (length(missing_files) > 0) {
  stop("Missing synthesis input(s):\n", paste(missing_files, collapse = "\n"))
}

# Helpers

pick_column <- function(data, choices) {
  found <- intersect(choices, names(data))
  
  if (length(found) != 1) {
    stop("Could not uniquely identify one of: ", paste(choices, collapse = ", "))
  }
  
  found
}

normalize_curve <- function(data, response_id, response) {
  if ("habitat_increment_percent_cover" %in% names(data) &&
      any(abs(data$habitat_increment_percent_cover - 1) > 1e-8, na.rm = TRUE)) {
    stop("Unexpected habitat increment in effect curves for ", response, ".")
  }
  
  effect_col <- pick_column(
    data,
    c("relative_change_percent", "relative_change_percent_at_typical_baseline")
  )
  
  lower_col <- pick_column(
    data,
    c("relative_CI_lower_percent", "relative_CI_lower_percent_at_typical_baseline")
  )
  
  upper_col <- pick_column(
    data,
    c("relative_CI_upper_percent", "relative_CI_upper_percent_at_typical_baseline")
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
  select(response_id, response, sample_file, city_effect_file) %>%
  pmap(function(response_id, response, sample_file, city_effect_file) {
    sample_summary <- read_csv(sample_file, show_col_types = FALSE)
    city_effects <- read_csv(city_effect_file, show_col_types = FALSE)
    
    required_sample_cols <- c("n_checklists", "n_localities", "n_cities")
    required_city_cols <- c(
      "NAME10", "avg_aridity_index", "habitat", "relative_change_percent"
    )
    
    if (nrow(sample_summary) != 1 ||
        !all(required_sample_cols %in% names(sample_summary))) {
      stop("Invalid sample summary for ", response, ".")
    }
    
    if (!all(required_city_cols %in% names(city_effects))) {
      stop("Missing city-conditional effect fields for ", response, ".")
    }
    
    if (anyNA(city_effects$NAME10) ||
        any(!is.finite(city_effects$avg_aridity_index)) ||
        !setequal(as.character(city_effects$habitat), habitat_order)) {
      stop("Invalid city-conditional effect data for ", response, ".")
    }
    
    if (nrow(city_effects) !=
        nrow(distinct(city_effects, habitat, NAME10))) {
      stop("Duplicated city-habitat effects for ", response, ".")
    }
    
    city_aridity <- city_effects %>%
      distinct(NAME10, avg_aridity_index)
    
    if (nrow(city_aridity) != n_distinct(city_effects$NAME10) ||
        nrow(city_aridity) != sample_summary$n_cities[[1]]) {
      stop("City/aridity reference mismatch for ", response, ".")
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
        transmute(response_id, response, NAME10, avg_aridity_index),
      city_effects = city_effects %>%
        transmute(
          response_id,
          response,
          habitat = as.character(habitat),
          NAME10 = as.character(NAME10),
          avg_aridity_index,
          relative_change_percent
        )
    )
  })

sample_reference <- map_dfr(sample_parts, "summary")
aridity_support <- map_dfr(sample_parts, "aridity")

city_effect_points <- map_dfr(sample_parts, "city_effects") %>%
  filter(response %in% main_responses) %>%
  mutate(
    response = factor(response, levels = main_responses),
    habitat = factor(habitat, levels = habitat_order)
  )

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

if (nrow(interaction_results) != 18 ||
    nrow(distinct(interaction_results, response, habitat)) != 18) {
  stop("Expected exactly 18 habitat-by-aridity interaction results.")
}

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

if (nrow(distinct(primary_curves, response, habitat)) != 18 ||
    any(!is.finite(primary_curves$relative_change_percent)) ||
    any(!is.finite(primary_curves$CI_lower_percent)) ||
    any(!is.finite(primary_curves$CI_upper_percent))) {
  stop("Final response curves are incomplete or invalid.")
}

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

# Figures

make_effect_figure <- function(responses, show_city_points = FALSE) {
  plot_curves <- primary_curves %>%
    filter(as.character(response) %in% responses) %>%
    mutate(
      response = factor(as.character(response), levels = responses),
      habitat = factor(as.character(habitat), levels = habitat_order)
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
      filter(as.character(response) %in% responses) %>%
      mutate(response = factor(as.character(response), levels = responses))
    
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

# Save figures

figure_outputs <- tibble(
  figure = list(figure_01, figure_02, figure_s1, figure_s2),
  directory = c(
    main_figure_dir, main_figure_dir, supp_figure_dir, supp_figure_dir
  ),
  filename = c(
    "Figure_01_Species_Richness_Resource_Importance",
    "Figure_02_Functional_Diversity_Resource_Importance",
    "Figure_S1_CWM_Body_Mass_Resource_Importance",
    "Figure_S2_CWM_Beak_Length_Resource_Importance"
  ),
  width = c(12.5, 12.5, 12.5, 12.5),
  height = c(4.5, 10, 4.5, 4.5)
)

pwalk(
  figure_outputs,
  function(figure, directory, filename, width, height) {
    ggsave(
      fs::path(directory, paste0(filename, ".png")),
      figure,
      width = width,
      height = height,
      units = "in",
      dpi = 600,
      bg = "white"
    )
    
    ggsave(
      fs::path(directory, paste0(filename, ".pdf")),
      figure,
      width = width,
      height = height,
      units = "in",
      bg = "white"
    )
  }
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

write_csv(
  aridity_support,
  fs::path(data_dir, "City_Aridity_Reference.csv")
)

write_csv(
  response_info %>%
    select(response, folder, interaction_file, curve_file, sample_file, city_effect_file),
  fs::path(data_dir, "Synthesis_Input_Manifest.csv")
)
