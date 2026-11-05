library(tidyverse)
library(glmmTMB)

# Settings

alpha_threshold <- 0.05
min_habitat_range_pp <- 1
display_half_range_pp <- 0.5
n_habitat_points <- 61

main_responses <- c(
  "Species richness",
  "Functional richness",
  "Functional dispersion",
  "Functional evenness"
)

analysis_root <- fs::path(
  "Outputs",
  "Analysis_Resource_Importance"
)

output_root <- fs::path(
  analysis_root,
  "09_Resource_Importance_Predictions"
)

data_dir <- fs::path(output_root, "Data")
table_dir <- fs::path(output_root, "Table")
figure_dir <- fs::path(output_root, "Figure")

walk(
  c(data_dir, table_dir, figure_dir),
  fs::dir_create
)

# Response and habitat inputs

response_info <- tribble(
  ~Response, ~folder, ~prefix, ~data_name, ~prediction_source,
  "Species richness",
  "02_Species_Richness",
  "Species_Richness",
  "Species_Richness_Resource_Importance_Data.rds",
  "saved_model",
  "Functional richness",
  "03_Functional_Richness",
  "Functional_Richness",
  "Functional_Richness_Resource_Importance_Data.rds",
  "saved_model",
  "Functional dispersion",
  "04_Functional_Dispersion",
  "Functional_Dispersion",
  "Functional_Dispersion_Resource_Importance_Data.rds",
  "saved_city_effects",
  "Functional evenness",
  "05_Functional_Evenness",
  "Functional_Evenness",
  "Functional_Evenness_Resource_Importance_Data.rds",
  "saved_model"
) %>%
  mutate(
    data_file = fs::path(
      analysis_root,
      folder,
      "Data",
      data_name
    ),
    interaction_file = fs::path(
      analysis_root,
      folder,
      "Table",
      paste0(
        prefix,
        "_Habitat_x_Aridity.csv"
      )
    ),
    city_effect_file = if_else(
      Response == "Functional dispersion",
      fs::path(
        analysis_root,
        folder,
        "Table",
        "Functional_Dispersion_City_Conditional_Effects.csv"
      ),
      NA_character_
    )
  )

habitat_info <- tribble(
  ~Habitat, ~stub, ~within_term,
  "Water",
  "Water",
  "water_within_model_10",
  "Shrub and scrub",
  "Shrub_and_scrub",
  "shrub_within_model_10",
  "Built",
  "Built",
  "built_within_model_10"
)

required_files <- c(
  response_info$data_file,
  response_info$interaction_file,
  response_info$city_effect_file[
    !is.na(response_info$city_effect_file)
  ]
)

missing_files <- required_files[
  !fs::file_exists(required_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing prediction input(s):\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

# Final supported interactions

interaction_results <- response_info %>%
  select(
    Response,
    folder,
    prediction_source,
    data_file,
    interaction_file,
    city_effect_file
  ) %>%
  pmap_dfr(
    function(
    Response,
    folder,
    prediction_source,
    data_file,
    interaction_file,
    city_effect_file
    ) {
      read_csv(
        interaction_file,
        show_col_types = FALSE
      ) %>%
        transmute(
          Response = Response,
          Habitat = habitat,
          interaction_estimate_1pp =
            if (
              "Estimate_1pp" %in%
              names(.)
            ) {
              Estimate_1pp
            } else {
              Estimate / 10
            },
          interaction_p = p,
          folder = folder,
          prediction_source =
            prediction_source,
          data_file = data_file,
          city_effect_file =
            city_effect_file
        )
    }
  ) %>%
  filter(
    interaction_p <
      alpha_threshold
  ) %>%
  left_join(
    habitat_info,
    by = "Habitat"
  ) %>%
  mutate(
    model_file = if_else(
      prediction_source ==
        "saved_model",
      fs::path(
        analysis_root,
        folder,
        "Models",
        paste0(
          stub,
          "_model.rds"
        )
      ),
      NA_character_
    )
  )

if (nrow(interaction_results) == 0) {
  stop(
    "No supported Habitat x Aridity interactions have P < 0.05."
  )
}

model_files <- interaction_results$model_file[
  !is.na(
    interaction_results$model_file
  )
]

missing_models <- model_files[
  !fs::file_exists(model_files)
]

if (length(missing_models) > 0) {
  stop(
    "Missing fitted model(s) for supported interactions:\n",
    paste(
      missing_models,
      collapse = "\n"
    )
  )
}

# Helpers

interaction_term <- function(
    fit,
    within_term
) {
  terms <- names(
    fixef(fit)$cond
  )
  
  out <- terms[
    str_detect(
      terms,
      fixed(within_term)
    ) &
      str_detect(
        terms,
        fixed(
          "aridity_centered"
        )
      )
  ]
  
  if (length(out) != 1) {
    stop(
      "Could not identify the Habitat x Aridity term for ",
      within_term,
      "."
    )
  }
  
  out
}

inverse_link <- function(
    link,
    x
) {
  if (link == "identity") {
    return(x)
  }
  
  if (link == "log") {
    return(exp(x))
  }
  
  if (link == "logit") {
    return(plogis(x))
  }
  
  stop(
    "Unsupported model link: ",
    link
  )
}

make_support <- function(
    analysis_data,
    within_term
) {
  analysis_data %>%
    distinct(
      NAME10,
      locality_id,
      habitat_deviation_pp =
        .data[[within_term]] * 10
    ) %>%
    group_by(NAME10) %>%
    summarise(
      n_localities = n(),
      habitat_min_pp =
        min(
          habitat_deviation_pp
        ),
      habitat_max_pp =
        max(
          habitat_deviation_pp
        ),
      habitat_range_pp =
        habitat_max_pp -
        habitat_min_pp,
      .groups = "drop"
    ) %>%
    mutate(
      left_support_pp =
        pmax(
          0,
          -habitat_min_pp
        ),
      right_support_pp =
        pmax(
          0,
          habitat_max_pp
        ),
      retained =
        habitat_range_pp >=
        min_habitat_range_pp &
        left_support_pp > 0 &
        right_support_pp > 0,
      plot_min_pp =
        pmax(
          habitat_min_pp,
          -display_half_range_pp
        ),
      plot_max_pp =
        pmin(
          habitat_max_pp,
          display_half_range_pp
        ),
      retained =
        retained &
        plot_max_pp >
        plot_min_pp
    )
}

make_prediction_grid <- function(
    city_reference
) {
  city_reference %>%
    filter(retained) %>%
    mutate(
      habitat_deviation_pp =
        map2(
          plot_min_pp,
          plot_max_pp,
          ~ sort(
            unique(
              c(
                seq(
                  .x,
                  .y,
                  length.out =
                    n_habitat_points
                ),
                0
              )
            )
          )
        )
    ) %>%
    unnest(
      habitat_deviation_pp
    )
}

check_baseline <- function(
    prediction_data,
    Response,
    Habitat
) {
  baseline_qc <- prediction_data %>%
    group_by(
      Response,
      Habitat,
      NAME10
    ) %>%
    slice_min(
      abs(
        habitat_deviation_pp
      ),
      n = 1,
      with_ties = FALSE
    ) %>%
    ungroup()
  
  if (max(
    abs(
      baseline_qc$
      predicted_relative_change_percent
    ),
    na.rm = TRUE
  ) > 1e-8) {
    stop(
      "Prediction normalization failed for ",
      Response,
      " / ",
      Habitat,
      "."
    )
  }
}

# Predictions from saved fitted models

make_model_prediction <- function(
    Response,
    Habitat,
    interaction_estimate_1pp,
    interaction_p,
    data_file,
    model_file,
    within_term
) {
  message(
    "Supported interaction prediction: ",
    Response,
    " / ",
    Habitat
  )
  
  analysis_data <- read_rds(
    data_file
  )
  
  fit <- read_rds(
    model_file
  )
  
  fit_frame <- tryCatch(
    model.frame(fit),
    error = function(e) {
      fit$frame
    }
  )
  
  if (is.null(fit_frame)) {
    stop(
      "Could not recover the model frame for ",
      Response,
      " / ",
      Habitat,
      "."
    )
  }
  
  model_family <-
    family(fit)$family
  
  model_link <-
    family(fit)$link
  
  if (!model_link %in%
      c(
        "identity",
        "log",
        "logit"
      )) {
    stop(
      "Unsupported model link for ",
      Response,
      " / ",
      Habitat,
      ": ",
      model_link
    )
  }
  
  beta <- fixef(fit)$cond
  
  int_term <- interaction_term(
    fit,
    within_term
  )
  
  city_re <- ranef(fit)$cond$NAME10 %>%
    as.data.frame() %>%
    rownames_to_column(
      "NAME10"
    ) %>%
    as_tibble()
  
  if (!all(
    c(
      "(Intercept)",
      within_term
    ) %in%
    names(city_re)
  )) {
    stop(
      "Could not recover city random intercept/slope for ",
      Response,
      " / ",
      Habitat,
      "."
    )
  }
  
  city_re <- city_re %>%
    transmute(
      NAME10,
      city_random_intercept =
        .data[["(Intercept)"]],
      city_random_slope =
        .data[[within_term]]
    )
  
  support <- make_support(
    analysis_data,
    within_term
  )
  
  reference_protocol <-
    analysis_data %>%
    count(
      protocol_name,
      sort = TRUE
    ) %>%
    slice(1) %>%
    pull(protocol_name) %>%
    as.character()
  
  reference_duration <-
    median(
      analysis_data$duration_10,
      na.rm = TRUE
    )
  
  reference_distance <-
    median(
      analysis_data$
        effort_distance_km,
      na.rm = TRUE
    )
  
  baseline_data <-
    analysis_data %>%
    distinct(
      NAME10,
      .keep_all = TRUE
    ) %>%
    mutate(
      "{within_term}" := 0,
      duration_10 =
        reference_duration,
      effort_distance_km =
        reference_distance,
      protocol_name =
        factor(
          reference_protocol,
          levels =
            levels(
              fit_frame$
                protocol_name
            )
        ),
      NAME10 =
        factor(
          as.character(
            NAME10
          ),
          levels =
            levels(
              fit_frame$NAME10
            )
        ),
      locality_id =
        factor(
          levels(
            fit_frame$
              locality_id
          )[1],
          levels =
            levels(
              fit_frame$
                locality_id
            )
        )
    )
  
  fixed_baseline_link <-
    as.numeric(
      predict(
        fit,
        newdata =
          baseline_data,
        type = "link",
        re.form = NA,
        allow.new.levels =
          FALSE
      )
    )
  
  city_reference <-
    baseline_data %>%
    transmute(
      NAME10 =
        as.character(
          NAME10
        ),
      avg_aridity_index,
      aridity_centered,
      fixed_baseline_link =
        fixed_baseline_link
    ) %>%
    left_join(
      city_re,
      by = "NAME10"
    ) %>%
    left_join(
      support %>%
        mutate(
          NAME10 =
            as.character(
              NAME10
            )
        ),
      by = "NAME10"
    ) %>%
    mutate(
      fixed_slope =
        unname(
          beta[[within_term]]
        ) +
        unname(
          beta[[int_term]]
        ) *
        aridity_centered,
      city_slope =
        fixed_slope +
        city_random_slope,
      baseline_link =
        fixed_baseline_link +
        city_random_intercept,
      baseline_response =
        inverse_link(
          model_link,
          baseline_link
        )
    )
  
  if (anyNA(
    city_reference$retained
  ) ||
  any(
    !is.finite(
      city_reference$city_slope
    )
  ) ||
  any(
    !is.finite(
      city_reference$
      baseline_response
    )
  ) ||
  any(
    city_reference$
    baseline_response <= 0
  )) {
    stop(
      "Invalid city reference for ",
      Response,
      " / ",
      Habitat,
      "."
    )
  }
  
  prediction_data <-
    make_prediction_grid(
      city_reference
    ) %>%
    mutate(
      prediction_link =
        baseline_link +
        city_slope *
        habitat_deviation_pp /
        10,
      predicted_response =
        inverse_link(
          model_link,
          prediction_link
        ),
      predicted_relative_change_percent =
        100 * (
          predicted_response /
            baseline_response - 1
        ),
      Response = Response,
      Habitat = Habitat,
      family = model_family,
      link = model_link,
      interaction_estimate_1pp =
        interaction_estimate_1pp,
      interaction_p =
        interaction_p,
      .before = 1
    ) %>%
    select(
      Response,
      Habitat,
      NAME10,
      avg_aridity_index,
      family,
      link,
      habitat_deviation_pp,
      predicted_response,
      baseline_response,
      predicted_relative_change_percent,
      interaction_estimate_1pp,
      interaction_p
    )
  
  check_baseline(
    prediction_data,
    Response,
    Habitat
  )
  
  audit <- city_reference %>%
    transmute(
      Response,
      Habitat,
      NAME10,
      avg_aridity_index,
      family =
        model_family,
      link =
        model_link,
      n_localities,
      habitat_min_pp,
      habitat_max_pp,
      habitat_range_pp,
      left_support_pp,
      right_support_pp,
      retained,
      city_slope,
      baseline_response,
      interaction_estimate_1pp =
        interaction_estimate_1pp,
      interaction_p =
        interaction_p
    )
  
  summary <- tibble(
    Response,
    Habitat,
    family =
      model_family,
    link =
      model_link,
    interaction_estimate_1pp =
      interaction_estimate_1pp,
    interaction_p =
      interaction_p,
    fitted_cities =
      n_distinct(
        city_reference$NAME10
      ),
    retained_cities =
      sum(
        city_reference$retained
      ),
    minimum_habitat_range_pp =
      min_habitat_range_pp,
    display_lower_pp =
      -display_half_range_pp,
    display_upper_pp =
      display_half_range_pp,
    prediction_points =
      nrow(
        prediction_data
      )
  )
  
  rm(
    analysis_data,
    fit,
    fit_frame,
    city_re,
    baseline_data,
    city_reference
  )
  
  gc()
  
  list(
    predictions =
      prediction_data,
    audit =
      audit,
    summary =
      summary
  )
}

# FDis prediction from saved city-specific effects

make_fdis_prediction <- function(
    Response,
    Habitat,
    interaction_estimate_1pp,
    interaction_p,
    data_file,
    city_effect_file,
    within_term
) {
  message(
    "Supported FDis interaction prediction: ",
    Habitat
  )
  
  analysis_data <- read_rds(
    data_file
  )
  
  city_effects <- read_csv(
    city_effect_file,
    show_col_types = FALSE
  ) %>%
    filter(
      habitat == Habitat
    )
  
  required_cols <- c(
    "NAME10",
    "avg_aridity_index",
    "family",
    "link",
    "city_slope",
    "baseline_response"
  )
  
  if (!all(
    required_cols %in%
    names(city_effects)
  )) {
    stop(
      "Final FDis city-effect output is missing required prediction fields."
    )
  }
  
  support <- make_support(
    analysis_data,
    within_term
  )
  
  city_reference <-
    city_effects %>%
    transmute(
      NAME10 =
        as.character(
          NAME10
        ),
      avg_aridity_index,
      family,
      link,
      city_slope,
      baseline_response
    ) %>%
    left_join(
      support %>%
        mutate(
          NAME10 =
            as.character(
              NAME10
            )
        ),
      by = "NAME10"
    )
  
  model_family <-
    unique(
      city_reference$family
    )
  
  model_link <-
    unique(
      city_reference$link
    )
  
  if (length(model_family) != 1 ||
      length(model_link) != 1 ||
      model_link != "identity") {
    stop(
      "Final FDis prediction input must contain one Gaussian identity-link model."
    )
  }
  
  prediction_data <-
    make_prediction_grid(
      city_reference
    ) %>%
    mutate(
      predicted_response =
        baseline_response +
        city_slope *
        habitat_deviation_pp /
        10,
      predicted_relative_change_percent =
        100 * (
          predicted_response /
            baseline_response - 1
        ),
      Response = Response,
      Habitat = Habitat,
      interaction_estimate_1pp =
        interaction_estimate_1pp,
      interaction_p =
        interaction_p,
      .before = 1
    ) %>%
    select(
      Response,
      Habitat,
      NAME10,
      avg_aridity_index,
      family,
      link,
      habitat_deviation_pp,
      predicted_response,
      baseline_response,
      predicted_relative_change_percent,
      interaction_estimate_1pp,
      interaction_p
    )
  
  check_baseline(
    prediction_data,
    Response,
    Habitat
  )
  
  audit <- city_reference %>%
    transmute(
      Response,
      Habitat,
      NAME10,
      avg_aridity_index,
      family,
      link,
      n_localities,
      habitat_min_pp,
      habitat_max_pp,
      habitat_range_pp,
      left_support_pp,
      right_support_pp,
      retained,
      city_slope,
      baseline_response,
      interaction_estimate_1pp =
        interaction_estimate_1pp,
      interaction_p =
        interaction_p
    )
  
  summary <- tibble(
    Response,
    Habitat,
    family =
      model_family,
    link =
      model_link,
    interaction_estimate_1pp =
      interaction_estimate_1pp,
    interaction_p =
      interaction_p,
    fitted_cities =
      n_distinct(
        city_reference$NAME10
      ),
    retained_cities =
      sum(
        city_reference$retained
      ),
    minimum_habitat_range_pp =
      min_habitat_range_pp,
    display_lower_pp =
      -display_half_range_pp,
    display_upper_pp =
      display_half_range_pp,
    prediction_points =
      nrow(
        prediction_data
      )
  )
  
  rm(
    analysis_data,
    city_effects,
    city_reference
  )
  
  gc()
  
  list(
    predictions =
      prediction_data,
    audit =
      audit,
    summary =
      summary
  )
}

# Generate supported-interaction predictions

prediction_parts <-
  interaction_results %>%
  pmap(
    function(
    Response,
    Habitat,
    interaction_estimate_1pp,
    interaction_p,
    folder,
    prediction_source,
    data_file,
    city_effect_file,
    stub,
    within_term,
    model_file
    ) {
      if (
        prediction_source ==
        "saved_model"
      ) {
        return(
          make_model_prediction(
            Response =
              Response,
            Habitat =
              Habitat,
            interaction_estimate_1pp =
              interaction_estimate_1pp,
            interaction_p =
              interaction_p,
            data_file =
              data_file,
            model_file =
              model_file,
            within_term =
              within_term
          )
        )
      }
      
      make_fdis_prediction(
        Response =
          Response,
        Habitat =
          Habitat,
        interaction_estimate_1pp =
          interaction_estimate_1pp,
        interaction_p =
          interaction_p,
        data_file =
          data_file,
        city_effect_file =
          city_effect_file,
        within_term =
          within_term
      )
    }
  )

prediction_data <- map_dfr(
  prediction_parts,
  "predictions"
)

city_support_audit <- map_dfr(
  prediction_parts,
  "audit"
)

prediction_summary <- map_dfr(
  prediction_parts,
  "summary"
)

# Panel order follows the original figure:
# Water first, then Shrub and scrub, then Built;
# responses follow manuscript order.

panel_order <-
  interaction_results %>%
  mutate(
    habitat_order =
      match(
        Habitat,
        c(
          "Water",
          "Shrub and scrub",
          "Built"
        )
      ),
    response_order =
      match(
        Response,
        main_responses
      )
  ) %>%
  arrange(
    habitat_order,
    response_order
  ) %>%
  transmute(
    panel =
      str_c(
        Response,
        "\n",
        Habitat
      )
  ) %>%
  pull(panel)

prediction_data <-
  prediction_data %>%
  mutate(
    panel = factor(
      str_c(
        Response,
        "\n",
        Habitat
      ),
      levels =
        panel_order
    )
  )

# Preview

figure_s1_preview <- ggplot(
  prediction_data,
  aes(
    x =
      habitat_deviation_pp,
    y =
      predicted_relative_change_percent,
    group =
      NAME10,
    colour =
      avg_aridity_index
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
    alpha = 0.35,
    linewidth = 0.5,
    lineend = "round"
  ) +
  scale_colour_viridis_c(
    option = "plasma",
    direction = -1,
    name =
      "De Martonne\naridity index\n(lower = more arid)"
  ) +
  scale_x_continuous(
    limits = c(
      -display_half_range_pp,
      display_half_range_pp
    ),
    breaks = seq(
      -display_half_range_pp,
      display_half_range_pp,
      length.out = 5
    ),
    expand = expansion(
      mult = c(
        0.02,
        0.02
      )
    )
  ) +
  facet_wrap(
    vars(panel),
    ncol = 3,
    scales = "free_y"
  ) +
  labs(
    x =
      "Local habitat relative to city mean (percentage points)",
    y =
      "Predicted change in response relative to city-mean habitat (%)"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    strip.background =
      element_blank(),
    strip.text =
      element_text(
        face = "bold",
        size = 10.5
      ),
    legend.position =
      "right",
    legend.title =
      element_text(
        size = 9.5
      ),
    legend.text =
      element_text(
        size = 9
      ),
    axis.title =
      element_text(
        size = 10.5
      ),
    axis.text =
      element_text(
        size = 9
      ),
    panel.spacing =
      grid::unit(
        0.9,
        "lines"
      )
  )

# Save outputs

write_csv(
  prediction_data %>%
    select(-panel),
  fs::path(
    data_dir,
    "Significant_Interaction_City_Predictions.csv"
  )
)

write_csv(
  city_support_audit,
  fs::path(
    table_dir,
    "Table_01_City_Habitat_Support_Audit.csv"
  )
)

write_csv(
  prediction_summary,
  fs::path(
    table_dir,
    "Table_02_Prediction_Summary.csv"
  )
)

write_csv(
  interaction_results %>%
    select(
      Response,
      Habitat,
      interaction_estimate_1pp,
      interaction_p
    ),
  fs::path(
    table_dir,
    "Table_03_Significant_Interactions_Predicted.csv"
  )
)

ggsave(
  fs::path(
    figure_dir,
    "Figure_S1_Supported_Habitat_x_Aridity_Predictions_PREVIEW.png"
  ),
  figure_s1_preview,
  width = 12.5,
  height = 8.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

message(
  "Script 09 complete: supported Habitat x Aridity predictions regenerated from final response-specific results."
)
