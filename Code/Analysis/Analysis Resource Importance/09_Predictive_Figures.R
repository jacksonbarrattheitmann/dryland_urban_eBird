library(tidyverse)
library(glmmTMB)

# Settings

alpha_threshold <- 0.05
min_habitat_range_pp <- 1
display_half_range_pp <- 0.5
n_habitat_points <- 61

synthesis_table <- fs::path(
  "Outputs", "Analysis_Resource_Importance",
  "08_Resource_Importance_Synthesis", "Tables",
  "Table_01_Full_Habitat_x_Aridity_Statistics.csv"
)

output_root <- fs::path(
  "Outputs", "Analysis_Resource_Importance",
  "09_Resource_Importance_Predictions"
)

data_dir <- fs::path(output_root, "Data")
table_dir <- fs::path(output_root, "Table")
figure_dir <- fs::path(output_root, "Figure")

walk(c(data_dir, table_dir, figure_dir), fs::dir_create)

if (!fs::file_exists(synthesis_table)) {
  stop("Run 08_Resource_Importance_Synthesis.R first.")
}

response_info <- tribble(
  ~Response, ~folder, ~data_name,
  "Species richness", "02_Species_Richness", "Species_Richness_Resource_Importance_Data.rds",
  "Functional richness", "03_Functional_Richness", "Functional_Richness_Resource_Importance_Data.rds",
  "Functional dispersion", "04_Functional_Dispersion", "Functional_Dispersion_Resource_Importance_Data.rds",
  "Functional evenness", "05_Functional_Evenness", "Functional_Evenness_Resource_Importance_Data.rds",
  "CWM body mass", "06_CWM_Body_Mass", "CWM_Body_Mass_Resource_Importance_Data.rds",
  "CWM beak length", "07_CWM_Beak_Length", "CWM_Beak_Length_Resource_Importance_Data.rds"
) %>%
  mutate(
    data_file = fs::path(
      "Outputs", "Analysis_Resource_Importance",
      folder, "Data", data_name
    )
  )

habitat_info <- tribble(
  ~Habitat, ~stub, ~within_term,
  "Water", "Water", "water_within_model_10",
  "Shrub and scrub", "Shrub_and_scrub", "shrub_within_model_10",
  "Built", "Built", "built_within_model_10"
)

significant_interactions <- read_csv(
  synthesis_table,
  show_col_types = FALSE
) %>%
  filter(`P-value` < alpha_threshold) %>%
  select(
    Response,
    Habitat,
    interaction_estimate_1pp = Habitat_x_Aridity_Estimate_per_1pp,
    interaction_p = `P-value`
  ) %>%
  left_join(response_info, by = "Response") %>%
  left_join(habitat_info, by = "Habitat") %>%
  mutate(
    model_file = fs::path(
      "Outputs", "Analysis_Resource_Importance",
      folder, "Models", paste0(stub, "_model.rds")
    )
  )

if (nrow(significant_interactions) == 0) {
  stop("No Habitat x Aridity interactions have P < 0.05.")
}

required_files <- c(
  significant_interactions$data_file,
  significant_interactions$model_file
)

if (any(!fs::file_exists(required_files))) {
  stop("At least one prediction input file is missing.")
}

interaction_term <- function(fit, within_term) {
  terms <- names(fixef(fit)$cond)
  out <- terms[
    str_detect(terms, fixed(within_term)) &
      str_detect(terms, fixed("aridity_centered"))
  ]
  
  if (length(out) != 1) {
    stop("Could not identify the Habitat x Aridity term for ", within_term, ".")
  }
  
  out
}

inverse_link <- function(link, x) {
  if (link == "log") {
    return(exp(x))
  }
  
  if (link == "logit") {
    return(plogis(x))
  }
  
  stop("Unsupported model link: ", link)
}

make_prediction <- function(
    Response,
    Habitat,
    interaction_estimate_1pp,
    interaction_p,
    folder,
    data_name,
    data_file,
    stub,
    within_term,
    model_file
) {
  analysis_data <- read_rds(data_file)
  fit <- read_rds(model_file)
  fit_frame <- tryCatch(model.frame(fit), error = function(e) fit$frame)
  
  if (is.null(fit_frame)) {
    stop("Could not recover the model frame for ", Response, " / ", Habitat, ".")
  }
  
  model_link <- family(fit)$link
  
  if (!model_link %in% c("log", "logit")) {
    stop("Predictions currently support only log and logit links.")
  }
  
  beta <- fixef(fit)$cond
  int_term <- interaction_term(fit, within_term)
  
  city_re <- ranef(fit)$cond$NAME10 %>%
    as.data.frame() %>%
    rownames_to_column("NAME10") %>%
    as_tibble()
  
  if (!all(c("(Intercept)", within_term) %in% names(city_re))) {
    stop("Could not recover city random intercept/slope for ", Response, " / ", Habitat, ".")
  }
  
  city_re <- city_re %>%
    transmute(
      NAME10,
      city_random_intercept = .data[["(Intercept)"]],
      city_random_slope = .data[[within_term]]
    )
  
  support <- analysis_data %>%
    distinct(
      NAME10,
      locality_id,
      habitat_deviation_pp = .data[[within_term]] * 10
    ) %>%
    group_by(NAME10) %>%
    summarise(
      n_localities = n(),
      habitat_min_pp = min(habitat_deviation_pp),
      habitat_max_pp = max(habitat_deviation_pp),
      habitat_range_pp = habitat_max_pp - habitat_min_pp,
      .groups = "drop"
    ) %>%
    mutate(
      left_support_pp = pmax(0, -habitat_min_pp),
      right_support_pp = pmax(0, habitat_max_pp),
      retained =
        habitat_range_pp >= min_habitat_range_pp &
        left_support_pp > 0 &
        right_support_pp > 0,
      plot_min_pp = pmax(habitat_min_pp, -display_half_range_pp),
      plot_max_pp = pmin(habitat_max_pp, display_half_range_pp),
      retained = retained & plot_max_pp > plot_min_pp
    )
  
  reference_protocol <- analysis_data %>%
    count(protocol_name, sort = TRUE) %>%
    slice(1) %>%
    pull(protocol_name) %>%
    as.character()
  
  reference_duration <- median(analysis_data$duration_10, na.rm = TRUE)
  reference_distance <- median(analysis_data$effort_distance_km, na.rm = TRUE)
  
  baseline_data <- analysis_data %>%
    distinct(NAME10, .keep_all = TRUE) %>%
    mutate(
      "{within_term}" := 0,
      duration_10 = reference_duration,
      effort_distance_km = reference_distance,
      protocol_name = factor(
        reference_protocol,
        levels = levels(fit_frame$protocol_name)
      ),
      NAME10 = factor(
        as.character(NAME10),
        levels = levels(fit_frame$NAME10)
      ),
      locality_id = factor(
        levels(fit_frame$locality_id)[1],
        levels = levels(fit_frame$locality_id)
      )
    )
  
  fixed_baseline_link <- as.numeric(
    predict(
      fit,
      newdata = baseline_data,
      type = "link",
      re.form = NA,
      allow.new.levels = FALSE
    )
  )
  
  city_reference <- baseline_data %>%
    transmute(
      NAME10 = as.character(NAME10),
      avg_aridity_index,
      aridity_centered,
      fixed_baseline_link = fixed_baseline_link
    ) %>%
    left_join(city_re, by = "NAME10") %>%
    left_join(support, by = "NAME10") %>%
    mutate(
      fixed_slope_link_per_10pp =
        unname(beta[[within_term]]) +
        unname(beta[[int_term]]) * aridity_centered,
      city_slope_link_per_10pp =
        fixed_slope_link_per_10pp + city_random_slope,
      baseline_link = fixed_baseline_link + city_random_intercept,
      baseline_response = inverse_link(model_link, baseline_link)
    )
  
  city_audit <- city_reference %>%
    transmute(
      Response,
      Habitat,
      NAME10,
      avg_aridity_index,
      n_localities,
      habitat_min_pp,
      habitat_max_pp,
      habitat_range_pp,
      left_support_pp,
      right_support_pp,
      retained,
      city_slope_link_per_10pp,
      interaction_estimate_1pp,
      interaction_p
    )
  
  prediction_data <- city_reference %>%
    filter(retained) %>%
    mutate(
      habitat_deviation_pp = map2(
        plot_min_pp,
        plot_max_pp,
        ~ sort(unique(c(seq(.x, .y, length.out = n_habitat_points), 0)))
      )
    ) %>%
    select(
      NAME10,
      avg_aridity_index,
      city_slope_link_per_10pp,
      baseline_link,
      baseline_response,
      habitat_deviation_pp
    ) %>%
    unnest(habitat_deviation_pp) %>%
    mutate(
      prediction_link =
        baseline_link +
        city_slope_link_per_10pp * habitat_deviation_pp / 10,
      predicted_response = inverse_link(model_link, prediction_link),
      predicted_relative_change_percent =
        100 * (predicted_response / baseline_response - 1),
      Response,
      Habitat,
      interaction_estimate_1pp,
      interaction_p,
      .before = 1
    )
  
  baseline_qc <- prediction_data %>%
    group_by(Response, Habitat, NAME10) %>%
    slice_min(abs(habitat_deviation_pp), n = 1, with_ties = FALSE) %>%
    ungroup()
  
  if (max(abs(baseline_qc$predicted_relative_change_percent), na.rm = TRUE) > 1e-8) {
    stop("Prediction normalization failed for ", Response, " / ", Habitat, ".")
  }
  
  prediction_summary <- tibble(
    Response,
    Habitat,
    interaction_estimate_1pp,
    interaction_p,
    fitted_cities = n_distinct(city_reference$NAME10),
    retained_cities = sum(city_reference$retained),
    minimum_habitat_range_pp = min_habitat_range_pp,
    display_lower_pp = -display_half_range_pp,
    display_upper_pp = display_half_range_pp,
    prediction_points = nrow(prediction_data)
  )
  
  city_aridity <- analysis_data %>%
    distinct(NAME10, avg_aridity_index)
  
  aridity_reference <- tibble(
    aridity_group = factor(
      c("More arid", "Median", "Less arid"),
      levels = c("More arid", "Median", "Less arid")
    ),
    avg_aridity_index = quantile(
      city_aridity$avg_aridity_index,
      c(0.10, 0.50, 0.90),
      names = FALSE
    )
  )
  
  mean_city_aridity <- mean(
    analysis_data$avg_aridity_index - analysis_data$aridity_centered
  )
  
  representative_grid <- aridity_reference %>%
    crossing(
      habitat_deviation_pp = seq(
        -display_half_range_pp,
        display_half_range_pp,
        length.out = n_habitat_points
      )
    )
  
  representative_data <- analysis_data %>%
    slice(rep(1, nrow(representative_grid))) %>%
    mutate(
      avg_aridity_index = representative_grid$avg_aridity_index,
      aridity_centered = avg_aridity_index - mean_city_aridity,
      "{within_term}" := representative_grid$habitat_deviation_pp / 10,
      across(starts_with("city_context_"), ~ 0),
      duration_10 = reference_duration,
      effort_distance_km = reference_distance,
      protocol_name = factor(
        reference_protocol,
        levels = levels(fit_frame$protocol_name)
      ),
      NAME10 = factor(
        levels(fit_frame$NAME10)[1],
        levels = levels(fit_frame$NAME10)
      ),
      locality_id = factor(
        levels(fit_frame$locality_id)[1],
        levels = levels(fit_frame$locality_id)
      )
    )
  
  representative_predictions <- representative_grid %>%
    mutate(
      predicted_response = as.numeric(
        predict(
          fit,
          newdata = representative_data,
          type = "response",
          re.form = NA,
          allow.new.levels = FALSE
        )
      )
    ) %>%
    group_by(aridity_group, avg_aridity_index) %>%
    mutate(
      baseline_response = predicted_response[
        which.min(abs(habitat_deviation_pp))
      ],
      predicted_relative_change_percent =
        100 * (predicted_response / baseline_response - 1)
    ) %>%
    ungroup() %>%
    mutate(
      Response,
      Habitat,
      interaction_estimate_1pp,
      interaction_p,
      .before = 1
    )
  
  out <- list(
    predictions = prediction_data,
    representative_predictions = representative_predictions,
    audit = city_audit,
    summary = prediction_summary
  )
  
  rm(
    analysis_data,
    fit,
    fit_frame,
    city_re,
    baseline_data,
    city_reference,
    representative_data
  )
  gc()
  
  out
}

prediction_parts <- significant_interactions %>%
  pmap(make_prediction)

prediction_data <- map_dfr(prediction_parts, "predictions")
representative_prediction_data <- map_dfr(
  prediction_parts,
  "representative_predictions"
)
city_support_audit <- map_dfr(prediction_parts, "audit")
prediction_summary <- map_dfr(prediction_parts, "summary")

panel_order <- significant_interactions %>%
  mutate(
    habitat_order = match(Habitat, c("Water", "Shrub and scrub", "Built")),
    response_order = match(
      Response,
      c(
        "Species richness",
        "Functional richness",
        "Functional dispersion",
        "Functional evenness",
        "CWM body mass",
        "CWM beak length"
      )
    )
  ) %>%
  arrange(habitat_order, response_order) %>%
  transmute(panel = str_c(Response, "\n", Habitat)) %>%
  pull(panel)

prediction_data <- prediction_data %>%
  mutate(
    panel = factor(
      str_c(Response, "\n", Habitat),
      levels = panel_order
    )
  )

representative_prediction_data <- representative_prediction_data %>%
  mutate(
    panel = factor(
      str_c(Response, "\n", Habitat),
      levels = panel_order
    )
  )

# Current city-specific figure

figure_03 <- ggplot(
  prediction_data,
  aes(
    x = habitat_deviation_pp,
    y = predicted_relative_change_percent,
    group = NAME10,
    colour = avg_aridity_index
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_line(
    alpha = 0.32,
    linewidth = 0.45,
    lineend = "round"
  ) +
  scale_colour_viridis_c(
    option = "plasma",
    direction = -1,
    name = "De Martonne\naridity index\n(lower = more arid)"
  ) +
  scale_x_continuous(
    limits = c(-display_half_range_pp, display_half_range_pp),
    breaks = c(-0.5, -0.25, 0, 0.25, 0.5),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  facet_wrap(
    vars(panel),
    ncol = 3,
    scales = "free_y"
  ) +
  labs(
    x = "Local habitat relative to city mean (percentage points)",
    y = "Predicted change in response relative to city-mean habitat (%)"
  ) +
  theme_classic(base_size = 12) +
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

# Candidate figure using representative aridity levels

figure_03_candidate <- ggplot(
  representative_prediction_data,
  aes(
    x = habitat_deviation_pp,
    y = predicted_relative_change_percent,
    group = aridity_group,
    colour = aridity_group
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    linewidth = 0.4
  ) +
  geom_line(
    linewidth = 1,
    lineend = "round"
  ) +
  scale_colour_viridis_d(
    option = "plasma",
    direction = -1,
    name = "Aridity"
  ) +
  scale_x_continuous(
    limits = c(-display_half_range_pp, display_half_range_pp),
    breaks = c(-0.5, -0.25, 0, 0.25, 0.5),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  facet_wrap(
    vars(panel),
    ncol = 3,
    scales = "free_y"
  ) +
  labs(
    x = "Local habitat cover relative to city mean (%)",
    y = "Predicted change in response relative to city-mean habitat (%)"
  ) +
  theme_classic(base_size = 12) +
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

# Save outputs

write_csv(
  prediction_data %>% select(-panel),
  fs::path(data_dir, "Significant_Interaction_City_Predictions.csv")
)

write_csv(
  representative_prediction_data %>% select(-panel),
  fs::path(data_dir, "Representative_Aridity_Predictions.csv")
)

write_csv(
  city_support_audit,
  fs::path(table_dir, "Table_01_City_Habitat_Support_Audit.csv")
)

write_csv(
  prediction_summary,
  fs::path(table_dir, "Table_02_Prediction_Summary.csv")
)

write_csv(
  significant_interactions %>%
    select(Response, Habitat, interaction_estimate_1pp, interaction_p),
  fs::path(table_dir, "Table_03_Significant_Interactions_Predicted.csv")
)

ggsave(
  fs::path(
    figure_dir,
    "Figure_03_Significant_Habitat_x_Aridity_Predictions.png"
  ),
  figure_03,
  width = 14,
  height = 8.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  fs::path(
    figure_dir,
    "Figure_03_Significant_Habitat_x_Aridity_Predictions.pdf"
  ),
  figure_03,
  width = 14,
  height = 8.5,
  units = "in",
  bg = "white"
)

ggsave(
  fs::path(
    figure_dir,
    "Figure_03_Candidate_Representative_Aridity_Predictions.png"
  ),
  figure_03_candidate,
  width = 14,
  height = 8.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  fs::path(
    figure_dir,
    "Figure_03_Candidate_Representative_Aridity_Predictions.pdf"
  ),
  figure_03_candidate,
  width = 14,
  height = 8.5,
  units = "in",
  bg = "white"
)
