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
trait_responses <- response_order[5:6]

habitat_order <- c(
  "Water",
  "Shrub and scrub",
  "Built"
)

analysis_root <- fs::path(
  "Outputs",
  "Analysis_Resource_Importance"
)

synthesis_root <- fs::path(
  analysis_root,
  "08_Resource_Importance_Synthesis"
)

prediction_root <- fs::path(
  analysis_root,
  "09_Resource_Importance_Predictions"
)

output_root <- fs::path(
  analysis_root,
  "10_Final_Figures"
)

main_figure_dir <- fs::path(
  output_root,
  "Main_Figures"
)

supp_figure_dir <- fs::path(
  output_root,
  "Supplementary_Figures"
)

walk(
  c(
    main_figure_dir,
    supp_figure_dir
  ),
  fs::dir_create
)

# Final inputs

curve_file <- fs::path(
  synthesis_root,
  "Data",
  "Final_Primary_Response_Effect_Curves.csv"
)

main_city_effect_file <- fs::path(
  synthesis_root,
  "Data",
  "Figure_01_City_Specific_Conditional_Effect_Points.csv"
)

main_interaction_file <- fs::path(
  synthesis_root,
  "Tables",
  "Table_S1_Biodiversity_Habitat_x_Aridity_Statistics.csv"
)

trait_interaction_file <- fs::path(
  synthesis_root,
  "Tables",
  "Table_S2_CWM_Trait_Habitat_x_Aridity_Statistics.csv"
)

prediction_file <- fs::path(
  prediction_root,
  "Data",
  "Significant_Interaction_City_Predictions.csv"
)

prediction_summary_file <- fs::path(
  prediction_root,
  "Table",
  "Table_02_Prediction_Summary.csv"
)

body_mass_city_file <- fs::path(
  analysis_root,
  "06_CWM_Body_Mass",
  "Table",
  "CWM_Body_Mass_City_Conditional_Effects.csv"
)

beak_length_city_file <- fs::path(
  analysis_root,
  "07_CWM_Beak_Length",
  "Table",
  "CWM_Beak_Length_City_Conditional_Effects.csv"
)

input_files <- c(
  curve_file,
  main_city_effect_file,
  main_interaction_file,
  trait_interaction_file,
  prediction_file,
  prediction_summary_file,
  body_mass_city_file,
  beak_length_city_file
)

missing_files <- input_files[
  !fs::file_exists(
    input_files
  )
]

if (length(missing_files) > 0) {
  stop(
    "Missing final-figure input(s). Run Scripts 08 and 09 first:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

# Load inputs

primary_curves <- read_csv(
  curve_file,
  show_col_types = FALSE
)

main_city_effects <- read_csv(
  main_city_effect_file,
  show_col_types = FALSE
)

main_interactions <- read_csv(
  main_interaction_file,
  show_col_types = FALSE
)

trait_interactions <- read_csv(
  trait_interaction_file,
  show_col_types = FALSE
)

prediction_data <- read_csv(
  prediction_file,
  show_col_types = FALSE
)

prediction_summary <- read_csv(
  prediction_summary_file,
  show_col_types = FALSE
)

body_mass_city <- read_csv(
  body_mass_city_file,
  show_col_types = FALSE
) %>%
  transmute(
    response =
      "CWM body mass",
    habitat,
    NAME10,
    avg_aridity_index,
    relative_change_percent
  )

beak_length_city <- read_csv(
  beak_length_city_file,
  show_col_types = FALSE
) %>%
  transmute(
    response =
      "CWM beak length",
    habitat,
    NAME10,
    avg_aridity_index,
    relative_change_percent
  )

trait_city_effects <- bind_rows(
  body_mass_city,
  beak_length_city
)

# Input checks

expected_curve_pairs <- expand.grid(
  response = response_order,
  habitat = habitat_order,
  stringsAsFactors = FALSE
)

observed_curve_pairs <-
  primary_curves %>%
  distinct(
    response,
    habitat
  )

pair_key <- function(x) {
  paste(
    x$response,
    x$habitat,
    sep = "::"
  )
}

if (!setequal(
  pair_key(
    expected_curve_pairs
  ),
  pair_key(
    observed_curve_pairs
  )
) ||
nrow(
  observed_curve_pairs
) != 18) {
  stop(
    "Effect curves must contain all 18 response-habitat combinations."
  )
}

if (nrow(main_interactions) != 12 ||
    nrow(
      distinct(
        main_interactions,
        Response,
        Habitat
      )
    ) != 12) {
  stop(
    "Main interaction table must contain 12 biodiversity response-habitat rows."
  )
}

if (nrow(trait_interactions) != 6 ||
    nrow(
      distinct(
        trait_interactions,
        Response,
        Habitat
      )
    ) != 6) {
  stop(
    "Trait interaction table must contain six response-habitat rows."
  )
}

if (nrow(prediction_data) == 0 ||
    nrow(prediction_summary) == 0) {
  stop(
    "Supported-interaction prediction outputs are empty."
  )
}

if (any(
  prediction_summary$
  interaction_p >= 0.05
)) {
  stop(
    "Figure S1 contains an interaction that is no longer supported at P < 0.05."
  )
}

# Figure helper for response-scale curves

make_effect_figure <- function(
    responses,
    city_points,
    n_rows
) {
  plot_curves <-
    primary_curves %>%
    filter(
      response %in%
        responses
    ) %>%
    mutate(
      response =
        factor(
          response,
          levels =
            responses
        ),
      habitat =
        factor(
          habitat,
          levels =
            habitat_order
        )
    )
  
  plot_points <-
    city_points %>%
    filter(
      response %in%
        responses
    ) %>%
    mutate(
      response =
        factor(
          response,
          levels =
            responses
        ),
      habitat =
        factor(
          habitat,
          levels =
            habitat_order
        )
    )
  
  ggplot(
    plot_curves,
    aes(
      x =
        avg_aridity_index,
      y =
        relative_change_percent
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype =
        "dashed",
      linewidth = 0.4
    ) +
    geom_ribbon(
      aes(
        ymin =
          CI_lower_percent,
        ymax =
          CI_upper_percent
      ),
      alpha = 0.15
    ) +
    geom_point(
      data =
        plot_points,
      aes(
        x =
          avg_aridity_index,
        y =
          relative_change_percent
      ),
      inherit.aes =
        FALSE,
      size = 0.8,
      alpha = 0.35
    ) +
    geom_line(
      linewidth = 0.9,
      lineend = "round"
    ) +
    facet_grid(
      rows =
        vars(response),
      cols =
        vars(habitat),
      scales =
        "free_y",
      switch = "y"
    ) +
    scale_x_reverse() +
    labs(
      x =
        "De Martonne aridity index",
      y =
        "Expected change in response (%)",
      caption =
        "Less arid  <-                                   ->  More arid"
    ) +
    theme_classic(
      base_size = 11
    ) +
    theme(
      strip.background =
        element_blank(),
      strip.text.x =
        element_text(
          face = "bold",
          size = 11.5
        ),
      strip.text.y.left =
        element_text(
          face = "bold",
          angle = 0,
          size = 10.5
        ),
      strip.placement =
        "outside",
      panel.spacing =
        grid::unit(
          0.8,
          "lines"
        ),
      axis.text =
        element_text(
          size = 9.5
        ),
      axis.title =
        element_text(
          size = 11
        ),
      plot.caption =
        element_text(
          hjust = 0.5,
          size = 9.5
        )
    )
}

# Figure 2

figure_02 <- make_effect_figure(
  responses =
    main_responses,
  city_points =
    main_city_effects,
  n_rows = 4
)

# Figure 3

figure_03_data <-
  main_interactions %>%
  mutate(
    Response =
      factor(
        Response,
        levels =
          main_responses
      ),
    Habitat =
      factor(
        Habitat,
        levels =
          habitat_order
      ),
    CI_lower =
      Estimate -
      1.96 * SE,
    CI_upper =
      Estimate +
      1.96 * SE
  )

figure_03 <- ggplot(
  figure_03_data,
  aes(
    x = Estimate,
    y = Habitat
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_errorbar(
    aes(
      xmin = CI_lower,
      xmax = CI_upper
    ),
    orientation = "y",
    width = 0
  ) +
  geom_point(
    size = 2.2
  ) +
  facet_wrap(
    vars(Response),
    ncol = 2,
    scales = "free_x"
  ) +
  labs(
    x =
      "Habitat x aridity interaction estimate per +1 percentage point (95% CI)",
    y = NULL
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
        size = 11
      ),
    panel.spacing =
      grid::unit(
        1.1,
        "lines"
      ),
    axis.text =
      element_text(
        size = 9.5
      ),
    axis.title.x =
      element_text(
        size = 11
      )
  )

# Figure S1

prediction_data <-
  prediction_data %>%
  mutate(
    panel =
      str_c(
        Response,
        "\n",
        Habitat
      )
  )

panel_order <-
  prediction_summary %>%
  mutate(
    habitat_order =
      match(
        Habitat,
        habitat_order
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
    panel =
      factor(
        panel,
        levels =
          panel_order
      )
  )

figure_s1 <- ggplot(
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
    limits =
      c(
        -0.5,
        0.5
      ),
    breaks =
      seq(
        -0.5,
        0.5,
        length.out = 5
      ),
    expand =
      expansion(
        mult =
          c(
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

# Figure S2

figure_s2 <- make_effect_figure(
  responses =
    trait_responses,
  city_points =
    trait_city_effects,
  n_rows = 2
)

# Figure S3

figure_s3_data <-
  trait_interactions %>%
  mutate(
    Response =
      factor(
        Response,
        levels =
          trait_responses
      ),
    Habitat =
      factor(
        Habitat,
        levels =
          habitat_order
      ),
    CI_lower =
      Estimate -
      1.96 * SE,
    CI_upper =
      Estimate +
      1.96 * SE
  )

figure_s3 <- ggplot(
  figure_s3_data,
  aes(
    x = Estimate,
    y = Habitat
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_errorbar(
    aes(
      xmin = CI_lower,
      xmax = CI_upper
    ),
    orientation = "y",
    width = 0
  ) +
  geom_point(
    size = 2.2
  ) +
  facet_wrap(
    vars(Response),
    nrow = 1,
    scales = "free_x"
  ) +
  labs(
    x =
      "Habitat x aridity interaction estimate per +1 percentage point (95% CI)",
    y = NULL
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
        size = 11
      ),
    panel.spacing =
      grid::unit(
        1.1,
        "lines"
      ),
    axis.text =
      element_text(
        size = 9.5
      ),
    axis.title.x =
      element_text(
        size = 11
      )
  )

# Remove stale outputs from prior layouts

stale_files <- c(
  fs::path(
    main_figure_dir,
    "Figure_01_Biodiversity_Resource_Importance.png"
  ),
  fs::path(
    main_figure_dir,
    "Figure_01_Biodiversity_Resource_Importance.pdf"
  ),
  fs::path(
    main_figure_dir,
    "Figure_01_Species_Richness_Resource_Importance.png"
  ),
  fs::path(
    main_figure_dir,
    "Figure_01_Species_Richness_Resource_Importance.pdf"
  ),
  fs::path(
    main_figure_dir,
    "Figure_02_Functional_Diversity_Resource_Importance.png"
  ),
  fs::path(
    main_figure_dir,
    "Figure_02_Functional_Diversity_Resource_Importance.pdf"
  ),
  fs::path(
    supp_figure_dir,
    "Figure_S3_City_Specific_Water_Predictions.png"
  ),
  fs::path(
    supp_figure_dir,
    "Figure_S3_City_Specific_Water_Predictions.pdf"
  )
)

stale_files <- stale_files[
  fs::file_exists(
    stale_files
  )
]

if (length(stale_files) > 0) {
  fs::file_delete(
    stale_files
  )
}

# Save original final-report figure suite

figure_outputs <- tibble(
  figure_id = c(
    "Figure 2",
    "Figure 3",
    "Figure S1",
    "Figure S2",
    "Figure S3"
  ),
  figure = list(
    figure_02,
    figure_03,
    figure_s1,
    figure_s2,
    figure_s3
  ),
  directory = c(
    main_figure_dir,
    main_figure_dir,
    supp_figure_dir,
    supp_figure_dir,
    supp_figure_dir
  ),
  filename = c(
    "Figure_02_Biodiversity_Resource_Importance",
    "Figure_03_Biodiversity_Habitat_x_Aridity_Interactions",
    "Figure_S1_Supported_Habitat_x_Aridity_Predictions",
    "Figure_S2_CWM_Trait_Resource_Importance",
    "Figure_S3_CWM_Trait_Habitat_x_Aridity_Interactions"
  ),
  width_in = c(
    12.5,
    10,
    12.5,
    12.5,
    10
  ),
  height_in = c(
    12.5,
    8,
    8.5,
    7,
    5.5
  )
)

pwalk(
  figure_outputs,
  function(
    figure_id,
    figure,
    directory,
    filename,
    width_in,
    height_in
  ) {
    path_png <- fs::path(
      directory,
      paste0(
        filename,
        ".png"
      )
    )
    
    path_pdf <- fs::path(
      directory,
      paste0(
        filename,
        ".pdf"
      )
    )
    
    ggsave(
      path_png,
      figure,
      width = width_in,
      height = height_in,
      units = "in",
      dpi = figure_dpi,
      bg = "white"
    )
    
    ggsave(
      path_pdf,
      figure,
      width = width_in,
      height = height_in,
      units = "in",
      bg = "white"
    )
    
    message(
      "Saved ",
      figure_id,
      ": ",
      path_png
    )
  }
)

# Reproducibility manifests

input_manifest <- tibble(
  input_file =
    input_files,
  size_bytes =
    as.numeric(
      file.info(
        input_files
      )$size
    ),
  modified =
    as.character(
      file.info(
        input_files
      )$mtime
    )
)

write_csv(
  input_manifest,
  fs::path(
    output_root,
    "Final_Figure_Input_Manifest.csv"
  )
)

write_csv(
  figure_outputs %>%
    select(
      -figure
    ),
  fs::path(
    output_root,
    "Final_Figure_Output_Manifest.csv"
  )
)

message(
  "Script 10 complete: original final-report figure suite regenerated with current model results."
)
