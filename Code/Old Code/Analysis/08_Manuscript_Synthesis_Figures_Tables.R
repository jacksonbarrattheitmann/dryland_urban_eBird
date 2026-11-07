# ============================================================
# Dryland Bird Project
# 08_Manuscript_Synthesis_Figures_Tables.R
#
# Purpose:
# Assemble manuscript-ready synthesis figures and a combined
# statistical table from the standardized outputs of Scripts
# 02-05.
#
# This script DOES NOT refit any models.
#
# Included responses:
#
#   1. Species richness
#   2. Functional richness
#   3. Functional dispersion
#   4. Functional evenness
#
# CWM body mass and CWM beak length are intentionally excluded
# from the synthesis figures and tables for now.
#
# Figures are organized BY PREDICTOR:
#
#   1. Aridity
#   2. Water cover
#   3. Shrub and scrub cover
#   4. Built cover
#
# For habitat figures, locality-level points are colored along
# the same warm-to-cool aridity gradient used in Scripts 02-05.
#
# Warm colors = more arid
# Cool colors = less arid
#
# Aridity percentile rank is used ONLY for visualization.
# All inferential models used the original continuous aridity
# index.
#
# A single combined statistical table is also created using the
# focal ecological terms from each of the four final models.
#
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggpubr)


# ------------------------------------------------------------
# 2. Input manuscript bundles
# ------------------------------------------------------------

input_files <- c(
  Species_Richness =
    file.path(
      "Outputs",
      "02_Species_Richness",
      "Data",
      "Species_Richness_Manuscript_Input.rds"
    ),
  
  Functional_Richness =
    file.path(
      "Outputs",
      "03_Functional_Richness",
      "Data",
      "Functional_Richness_Manuscript_Input.rds"
    ),
  
  Functional_Dispersion =
    file.path(
      "Outputs",
      "04_Functional_Dispersion",
      "Data",
      "Functional_Dispersion_Manuscript_Input.rds"
    ),
  
  Functional_Evenness =
    file.path(
      "Outputs",
      "05_Functional_Evenness",
      "Data",
      "Functional_Evenness_Manuscript_Input.rds"
    )
)


missing_input_files <- input_files[
  !file.exists(
    input_files
  )
]


if (
  length(
    missing_input_files
  ) > 0
) {
  
  stop(
    paste0(
      "The following manuscript-input RDS files are missing:\n\n",
      paste(
        missing_input_files,
        collapse = "\n"
      ),
      "\n\nRun the corresponding analysis scripts before running Script 08."
    )
  )
}


cat("\n============================================================\n")
cat("MANUSCRIPT SYNTHESIS\n")
cat("============================================================\n\n")


cat(
  "PASS: All four manuscript-input files were found.\n"
)


# ------------------------------------------------------------
# 3. Output directories
# ------------------------------------------------------------

OUTPUT_ROOT <- file.path(
  "Outputs",
  "08_Manuscript_Synthesis"
)


FIGURE_DIR <- file.path(
  OUTPUT_ROOT,
  "Figure"
)


TABLE_DIR <- file.path(
  OUTPUT_ROOT,
  "Table"
)


for (
  x in c(
    FIGURE_DIR,
    TABLE_DIR
  )
) {
  
  dir.create(
    x,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ------------------------------------------------------------
# 4. Read manuscript bundles
# ------------------------------------------------------------

bundles <- lapply(
  input_files,
  readRDS
)


# ------------------------------------------------------------
# 5. Validate manuscript bundles
# ------------------------------------------------------------

required_bundle_elements <- c(
  "response_id",
  "response_label",
  "response_variable",
  "analysis_data",
  "final_model",
  "smooth_results",
  "parametric_results",
  "focal_results",
  "model_metadata",
  "plot_data"
)


bundle_qc <- lapply(
  names(
    bundles
  ),
  function(bundle_name) {
    
    bundle <- bundles[[bundle_name]]
    
    
    missing_elements <- setdiff(
      required_bundle_elements,
      names(
        bundle
      )
    )
    
    
    data.frame(
      response =
        bundle_name,
      
      complete =
        length(
          missing_elements
        ) == 0,
      
      missing_elements =
        ifelse(
          length(
            missing_elements
          ) == 0,
          "",
          paste(
            missing_elements,
            collapse = "; "
          )
        ),
      
      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows()


print(
  bundle_qc,
  row.names = FALSE
)


if (
  any(
    !bundle_qc$complete
  )
) {
  
  stop(
    "One or more manuscript-input bundles are incomplete."
  )
}


cat(
  "\nPASS: All manuscript-input bundles have the expected structure.\n"
)


# ------------------------------------------------------------
# 6. Response metadata for plotting
# ------------------------------------------------------------

response_info <- data.frame(
  response_id = c(
    "Species_Richness",
    "Functional_Richness",
    "Functional_Dispersion",
    "Functional_Evenness"
  ),
  
  response_label = c(
    "Species richness",
    "Functional richness",
    "Functional dispersion",
    "Functional evenness"
  ),
  
  city_y = c(
    "mean_sr",
    "mean_fric",
    "mean_fdis",
    "mean_feve"
  ),
  
  locality_y = c(
    "mean_sr",
    "mean_fric",
    "mean_fdis",
    "mean_feve"
  ),
  
  stringsAsFactors = FALSE
)


response_ids <- response_info$response_id


# ------------------------------------------------------------
# 7. Response metadata QC
# ------------------------------------------------------------

if (
  anyDuplicated(
    response_info$response_id
  ) > 0
) {
  
  stop(
    "Duplicate response_id values are present in response_info."
  )
}


if (
  any(
    is.na(
      response_info$city_y
    )
  ) ||
  any(
    is.na(
      response_info$locality_y
    )
  )
) {
  
  stop(
    "Missing plotting-variable names are present in response_info."
  )
}


cat(
  "\nPASS: Response metadata are unique and complete.\n"
)


# ------------------------------------------------------------
# 8. Standard aridity palette
# ------------------------------------------------------------

aridity_palette <- c(
  "#7F0000",
  "#B30000",
  "#D7301F",
  "#F46D43",
  "#FDAE61",
  "#FEE08B",
  "#D9D9D9",
  "#E0F3F8",
  "#ABD9E9",
  "#74ADD1",
  "#4575B4",
  "#313695"
)


aridity_color_scale <- scale_color_gradientn(
  colors =
    aridity_palette,
  
  values = c(
    0.00,
    0.07,
    0.15,
    0.25,
    0.36,
    0.46,
    0.50,
    0.56,
    0.66,
    0.76,
    0.88,
    1.00
  ),
  
  limits = c(
    0,
    1
  ),
  
  breaks = c(
    0,
    0.5,
    1
  ),
  
  labels = c(
    "More arid",
    "Intermediate",
    "Less arid"
  ),
  
  name =
    "Aridity gradient"
)


# ------------------------------------------------------------
# 9. Common manuscript theme
# ------------------------------------------------------------

theme_synthesis <- function() {
  
  theme_classic(
    base_size = 12
  ) +
    theme(
      axis.text =
        element_text(
          colour = "black"
        ),
      
      plot.title =
        element_text(
          face = "bold",
          size = 12,
          hjust = 0.5
        ),
      
      legend.title =
        element_text(
          size = 11
        ),
      
      legend.text =
        element_text(
          size = 10
        )
    )
}


# ------------------------------------------------------------
# 10. Helper: get response plotting information
# ------------------------------------------------------------

get_response_info <- function(
    response_id
) {
  
  out <- response_info %>%
    filter(
      .data$response_id ==
        .env$response_id
    )
  
  
  if (
    nrow(
      out
    ) != 1
  ) {
    
    stop(
      paste0(
        "Could not uniquely identify response metadata for: ",
        response_id
      )
    )
  }
  
  
  out
}


# ------------------------------------------------------------
# 11. Helper: obtain a bundle by response ID
# ------------------------------------------------------------

get_bundle <- function(
    response_id
) {
  
  matches <- which(
    vapply(
      bundles,
      function(x) {
        identical(
          x$response_id,
          response_id
        )
      },
      logical(
        1
      )
    )
  )
  
  
  if (
    length(
      matches
    ) != 1
  ) {
    
    stop(
      paste0(
        "Could not uniquely identify bundle for response: ",
        response_id
      )
    )
  }
  
  
  bundles[[matches]]
}


# ------------------------------------------------------------
# 12. Helper: aridity panel
# ------------------------------------------------------------

make_aridity_panel <- function(
    response_id
) {
  
  info <- get_response_info(
    response_id
  )
  
  
  bundle <- get_bundle(
    response_id
  )
  
  
  plot_data <- bundle$plot_data$city
  
  
  y_variable <- info$city_y[
    1
  ]
  
  
  if (
    !y_variable %in%
    names(
      plot_data
    )
  ) {
    
    stop(
      paste0(
        "City-level plotting variable missing for ",
        response_id,
        ": ",
        y_variable
      )
    )
  }
  
  
  ggplot(
    plot_data,
    aes(
      x = avg_aridity_index,
      y = .data[[y_variable]]
    )
  ) +
    geom_point(
      alpha = 0.45,
      size = 1.7
    ) +
    geom_smooth(
      method = "gam",
      formula = y ~ s(
        x,
        k = 4
      ),
      se = TRUE,
      color = "black",
      fill = "grey80"
    ) +
    labs(
      x =
        "Aridity index (higher = less arid)",
      
      y =
        info$response_label[
          1
        ],
      
      title =
        info$response_label[
          1
        ]
    ) +
    theme_synthesis() +
    theme(
      legend.position =
        "none"
    )
}


# ------------------------------------------------------------
# 13. Helper: habitat panel
# ------------------------------------------------------------

make_habitat_panel <- function(
    response_id,
    habitat
) {
  
  info <- get_response_info(
    response_id
  )
  
  
  bundle <- get_bundle(
    response_id
  )
  
  
  plot_data <- bundle$plot_data$locality
  
  
  y_variable <- info$locality_y[
    1
  ]
  
  
  habitat_column <- switch(
    habitat,
    
    water =
      "water_plot",
    
    shrub =
      "shrub_plot",
    
    built =
      "built_plot",
    
    stop(
      "Unknown habitat: ",
      habitat
    )
  )
  
  
  x_label <- switch(
    habitat,
    
    water =
      "Water cover (%)",
    
    shrub =
      "Shrub and scrub cover (%)",
    
    built =
      "Built cover (%)"
  )
  
  
  required_columns <- c(
    y_variable,
    habitat_column,
    "aridity_rank"
  )
  
  
  missing_columns <- setdiff(
    required_columns,
    names(
      plot_data
    )
  )
  
  
  if (
    length(
      missing_columns
    ) > 0
  ) {
    
    stop(
      paste0(
        "Missing locality plotting columns for ",
        response_id,
        ": ",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
  }
  
  
  ggplot(
    plot_data,
    aes(
      x = .data[[habitat_column]],
      y = .data[[y_variable]]
    )
  ) +
    geom_point(
      aes(
        color =
          aridity_rank
      ),
      alpha = 0.55,
      size = 1.45
    ) +
    geom_smooth(
      method = "gam",
      formula = y ~ s(
        x,
        k = 4
      ),
      se = TRUE,
      color = "black",
      fill = "grey80"
    ) +
    aridity_color_scale +
    labs(
      x =
        x_label,
      
      y =
        info$response_label[
          1
        ],
      
      title =
        info$response_label[
          1
        ]
    ) +
    theme_synthesis()
}


# ------------------------------------------------------------
# 14. Helper: create four-response predictor figure
# ------------------------------------------------------------

make_synthesis_figure <- function(
    predictor
) {
  
  plot_list <- lapply(
    response_ids,
    function(response_id) {
      
      if (
        predictor ==
        "aridity"
      ) {
        
        make_aridity_panel(
          response_id
        )
        
      } else {
        
        make_habitat_panel(
          response_id =
            response_id,
          
          habitat =
            predictor
        )
      }
    }
  )
  
  
  ggpubr::ggarrange(
    plotlist =
      plot_list,
    
    ncol =
      2,
    
    nrow =
      2,
    
    labels =
      LETTERS[
        seq_along(
          plot_list
        )
      ],
    
    font.label = list(
      size = 14,
      face = "bold"
    ),
    
    common.legend =
      predictor !=
      "aridity",
    
    legend =
      "bottom"
  )
}


# ============================================================
# SYNTHESIS FIGURES
# ============================================================


# ------------------------------------------------------------
# 15. Aridity
# ------------------------------------------------------------

fig_aridity <- make_synthesis_figure(
  predictor =
    "aridity"
)


print(
  fig_aridity
)


# ------------------------------------------------------------
# 16. Water
# ------------------------------------------------------------

fig_water <- make_synthesis_figure(
  predictor =
    "water"
)


print(
  fig_water
)


# ------------------------------------------------------------
# 17. Shrub and scrub
# ------------------------------------------------------------

fig_shrub <- make_synthesis_figure(
  predictor =
    "shrub"
)


print(
  fig_shrub
)


# ------------------------------------------------------------
# 18. Built
# ------------------------------------------------------------

fig_built <- make_synthesis_figure(
  predictor =
    "built"
)


print(
  fig_built
)


# ============================================================
# SAVE SYNTHESIS FIGURES
# ============================================================


# ------------------------------------------------------------
# 19. Helper for figure saving
# ------------------------------------------------------------

save_synthesis_figure <- function(
    figure,
    filename_stub
) {
  
  ggsave(
    filename =
      file.path(
        FIGURE_DIR,
        paste0(
          filename_stub,
          ".png"
        )
      ),
    
    plot =
      figure,
    
    width =
      10,
    
    height =
      8,
    
    units =
      "in",
    
    dpi =
      600,
    
    bg =
      "white"
  )
  
  
  ggsave(
    filename =
      file.path(
        FIGURE_DIR,
        paste0(
          filename_stub,
          ".pdf"
        )
      ),
    
    plot =
      figure,
    
    width =
      10,
    
    height =
      8,
    
    units =
      "in",
    
    bg =
      "white"
  )
}


# ------------------------------------------------------------
# 20. Save four predictor-centered figures
# ------------------------------------------------------------

save_synthesis_figure(
  figure =
    fig_aridity,
  
  filename_stub =
    "Figure_Aridity_Core_Responses"
)


save_synthesis_figure(
  figure =
    fig_water,
  
  filename_stub =
    "Figure_Water_Core_Responses"
)


save_synthesis_figure(
  figure =
    fig_shrub,
  
  filename_stub =
    "Figure_Shrub_Core_Responses"
)


save_synthesis_figure(
  figure =
    fig_built,
  
  filename_stub =
    "Figure_Built_Core_Responses"
)


# ============================================================
# MASTER STATISTICAL TABLE
# ============================================================


# ------------------------------------------------------------
# 21. Standardize and combine focal model statistics
# ------------------------------------------------------------

master_statistics_list <- lapply(
  response_ids,
  function(response_id) {
    
    bundle <- get_bundle(
      response_id
    )
    
    
    info <- get_response_info(
      response_id
    )
    
    
    stats <- bundle$focal_results
    
    
    required_stats_columns <- c(
      "Predictor",
      "EDF",
      "Ref. df",
      "F statistic",
      "P-value"
    )
    
    
    missing_stats_columns <- setdiff(
      required_stats_columns,
      names(
        stats
      )
    )
    
    
    if (
      length(
        missing_stats_columns
      ) > 0
    ) {
      
      stop(
        paste0(
          "Missing focal-results columns for ",
          response_id,
          ": ",
          paste(
            missing_stats_columns,
            collapse = ", "
          )
        )
      )
    }
    
    
    family_name <- bundle$model_metadata$model_family[
      1
    ]
    
    
    stats %>%
      mutate(
        Response =
          info$response_label[
            1
          ],
        
        Family =
          family_name,
        
        .before =
          1
      )
  }
)


table_all_model_statistics <- bind_rows(
  master_statistics_list
)


# ------------------------------------------------------------
# 22. Order master statistics table
# ------------------------------------------------------------

predictor_order <- c(
  "Aridity",
  "Water cover",
  "Shrub and scrub cover",
  "Built cover",
  "Aridity x Water",
  "Aridity x Shrub and scrub",
  "Aridity x Built"
)


response_order <- response_info$response_label


table_all_model_statistics <- table_all_model_statistics %>%
  mutate(
    Response =
      factor(
        Response,
        levels =
          response_order
      ),
    
    Predictor =
      factor(
        Predictor,
        levels =
          predictor_order
      )
  ) %>%
  arrange(
    Response,
    Predictor
  ) %>%
  mutate(
    Response =
      as.character(
        Response
      ),
    
    Predictor =
      as.character(
        Predictor
      )
  )


# ------------------------------------------------------------
# 23. Print master statistics table
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("MASTER MODEL STATISTICS TABLE\n")
cat("============================================================\n\n")


print(
  table_all_model_statistics,
  row.names = FALSE
)


# ------------------------------------------------------------
# 24. Save master statistics table
# ------------------------------------------------------------

write.csv(
  table_all_model_statistics,
  file =
    file.path(
      TABLE_DIR,
      "Table_Core_Model_Statistics.csv"
    ),
  row.names = FALSE
)


write.table(
  table_all_model_statistics,
  file =
    file.path(
      TABLE_DIR,
      "Table_Core_Model_Statistics.txt"
    ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


# ============================================================
# MODEL-SUMMARY TABLE
# ============================================================


# ------------------------------------------------------------
# 25. Combine model metadata
# ------------------------------------------------------------

model_summary_list <- lapply(
  response_ids,
  function(response_id) {
    
    bundle <- get_bundle(
      response_id
    )
    
    
    metadata <- bundle$model_metadata
    
    
    data.frame(
      Response =
        bundle$response_label,
      
      Family =
        metadata$model_family[
          1
        ],
      
      Link =
        metadata$link[
          1
        ],
      
      Cities =
        metadata$n_cities[
          1
        ],
      
      Localities =
        metadata$n_localities[
          1
        ],
      
      Checklists =
        metadata$n_checklists[
          1
        ],
      
      Adjusted_R2 =
        round(
          metadata$adjusted_R2[
            1
          ],
          3
        ),
      
      Deviance_explained_percent =
        round(
          100 *
            metadata$deviance_explained[
              1
            ],
          2
        ),
      
      Family_parameter =
        round(
          metadata$family_parameter[
            1
          ],
          3
        ),
      
      Family_parameter_name =
        metadata$family_parameter_name[
          1
        ],
      
      stringsAsFactors = FALSE
    )
  }
)


table_model_summary <- bind_rows(
  model_summary_list
)


print(
  table_model_summary,
  row.names = FALSE
)


write.csv(
  table_model_summary,
  file =
    file.path(
      TABLE_DIR,
      "Table_Core_Model_Summary.csv"
    ),
  row.names = FALSE
)


# ============================================================
# OUTPUT QC
# ============================================================


# ------------------------------------------------------------
# 26. Expected outputs
# ------------------------------------------------------------

expected_outputs <- c(
  file.path(
    FIGURE_DIR,
    "Figure_Aridity_Core_Responses.png"
  ),
  
  file.path(
    FIGURE_DIR,
    "Figure_Water_Core_Responses.png"
  ),
  
  file.path(
    FIGURE_DIR,
    "Figure_Shrub_Core_Responses.png"
  ),
  
  file.path(
    FIGURE_DIR,
    "Figure_Built_Core_Responses.png"
  ),
  
  file.path(
    TABLE_DIR,
    "Table_Core_Model_Statistics.csv"
  ),
  
  file.path(
    TABLE_DIR,
    "Table_Core_Model_Summary.csv"
  )
)


output_qc <- data.frame(
  file =
    expected_outputs,
  
  exists =
    file.exists(
      expected_outputs
    ),
  
  stringsAsFactors = FALSE
)


cat("\n============================================================\n")
cat("SYNTHESIS OUTPUT QC\n")
cat("============================================================\n\n")


print(
  output_qc,
  row.names = FALSE
)


if (
  !all(
    output_qc$exists
  )
) {
  
  stop(
    "One or more expected Script 08 outputs were not created."
  )
}


cat(
  "\nPASS: All expected synthesis outputs were created.\n"
)


# ============================================================
# FINAL SUMMARY
# ============================================================


cat("\n============================================================\n")
cat("MANUSCRIPT SYNTHESIS COMPLETE\n")
cat("============================================================\n\n")


cat(
  "Synthesis figure directory:\n",
  FIGURE_DIR,
  "\n\n",
  sep = ""
)


cat(
  "Master statistics table:\n",
  file.path(
    TABLE_DIR,
    "Table_Core_Model_Statistics.csv"
  ),
  "\n\n",
  sep = ""
)


cat(
  "Model summary table:\n",
  file.path(
    TABLE_DIR,
    "Table_Core_Model_Summary.csv"
  ),
  "\n",
  sep = ""
)


# ============================================================
# INTERPRETATION OF SCRIPT PURPOSE
# ============================================================
#
# Script 08 is a presentation/synthesis script only.
#
# It does not change:
#
#   - response definitions
#   - sampling
#   - model families
#   - model formulas
#   - fitted models
#   - statistical inference
#
# The four synthesis figures compare species richness, functional
# richness, functional dispersion, and functional evenness across
# the same focal environmental predictor.
#
# CWM body mass and CWM beak length are intentionally excluded
# from this version and can be added later if needed.
#
# ============================================================
