# **Aridity alters biodiversity responses to resource availability in dryland cities** 

## Summary

This repository contains the data, R scripts, and analytical outputs associated with a study examining how local habitat cover and aridity influence bird biodiversity in dryland cities across the western United States. Using eBird observations from 2015–2025, AVONET functional trait data, Google Earth Engine land-cover estimates, and WorldClim climate data, the study evaluates relationships between water, shrub and scrub, and built habitat cover and multiple measures of avian community diversity. Response variables include species richness, functional richness (FRic), functional dispersion (FDis), functional evenness (FEve), and community-weighted mean body mass and beak length. Mixed-effects models were used to assess how local habitat relationships vary across the De Martonne aridity gradient while accounting for sampling effort and variation among cities and localities. The repository provides the datasets, statistical analyses, and figures supporting the study.

## Code

The `Code` folder contains R scripts used for data preparation, statistical modelling, and figure generation. Scripts are organized into two subfolders: `Data Prep` and `Analysis`.

### Data Prep

Located in `Code/Data Prep/`.

-   **01_Set_Up_Dryland_Shapefile.R:** Clips the global drylands shapefile to the contiguous United States and retains classified hyperarid, arid, semiarid, and dry subhumid regions. Exports processed shapefiles and maps of the study area.

-   **02_Find_Dryland_Cities.R:** Identifies dryland cities by calculating the proportion of each urban polygon overlapping the dryland study area. Retains cities with at least 50% dryland coverage and exports the corresponding spatial layers and city metadata.

-   **03_Unpack_eBird_TAR.R:** Extracts and validates the compressed eBird Sampling Event Data from the July 2026 United States eBird archive. The full eBird Basic Dataset remains compressed within the original archive for subsequent processing.

-   **04_Filter_eBird_Sampling_Events_to_Dryland_Cities.R:** Filters eBird sampling events from 2015–2025 to complete, single-observer checklists lasting 5–240 minutes. Retains checklists located within the selected dryland city polygons and generates checklist metadata, city-level summaries, and locality information.

-   **05_Filter_eBird_Observations_to_Dryland_Cities.R:** Filters the eBird Basic Dataset using the retained checklist identifiers and calculates species richness at the checklist and city-year levels. Produces processed observation summaries and locality information for subsequent analyses and land-cover extraction.

-   **06_Get_WorldClim_Annual_Aridity_Dryland_Cities.R:** Extracts historical WorldClim temperature and precipitation data for 2015–2024 and calculates annual De Martonne Aridity Index values for each dryland city. Integrates climate information with the processed eBird datasets for subsequent analyses.

-   **07_write_shapefile_for_GEE.R:** Converts eBird locality coordinates into spatial point data and exports shapefiles in batches for land-cover extraction using Google Earth Engine (GEE). Also prepares additional batches for localities requiring repeated extraction.

-   **08_combining_GEE_landCover_csvs.R:** Combines GEE-derived land-cover estimates within 1-km buffers around eBird localities, including records recovered through additional extraction attempts. Removes duplicate locality records and exports the consolidated land-cover dataset.

-   **10_cleaning_AVONET_data.R:** Processes AVONET bird trait data and harmonizes species names with eBird taxonomy using a taxonomic crosswalk. Exports a cleaned trait dataset for functional diversity and community-weighted mean analyses.

### Analysis

Located in `Code/Analysis/Analysis Resource Importance/`.

The scripts conduct statistical analyses examining how local habitat cover (water, shrub and scrub, and built environments) interacts with city-level aridity to influence bird species richness, functional diversity, and community-weighted mean traits. Analyses use mixed-effects models that account for variation among cities and sampling localities.

-   **01_Data_Setup.R:** Prepares the combined dataset for statistical analyses by filtering eBird checklists, standardizing sampling effort variables, calculating city-level habitat means, and centering local habitat cover relative to city means. Exports the standardized analysis dataset and quality-control summaries.

-   **02_Species_Richness.R:** Fits negative binomial mixed-effects models examining relationships between habitat cover, aridity, and bird species richness. Evaluates habitat-by-aridity interactions, performs a water-cover sensitivity analysis, and generates model diagnostics, statistical summaries, and conditional effect figures.

-   **03_Functional_Richness.R:** Fits beta mixed-effects models examining relationships between habitat cover, aridity, and functional richness (FRic). Generates habitat-by-aridity interaction estimates, relative percentage-change predictions, sensitivity analyses, and model diagnostic outputs.

-   **04_Functional_Dispersion.R:** Evaluates relationships between habitat cover, aridity, and functional dispersion (FDis) using mixed-effects models. Compares Gaussian and Gamma model families where appropriate, selects the supported model structure, and generates interaction statistics, conditional effects, sensitivity analyses, and diagnostic summaries.

-   **05_Functional_Evenness.R:** Fits beta mixed-effects models examining relationships between habitat cover, aridity, and functional evenness (FEve). Produces habitat-by-aridity interaction estimates, relative percentage-change predictions, sensitivity analyses, and model diagnostic outputs.

-   **06_CWM_Body_Mass.R:** Fits Gamma mixed-effects models to examine relationships between habitat cover, aridity, and community-weighted mean bird body mass. Produces model statistics, conditional effects, sensitivity analyses, and figures used in the supplementary analyses.

-   **07_CWM_Beak_Length.R:** Fits Gamma mixed-effects models to examine relationships between habitat cover, aridity, and community-weighted mean bird beak length. Produces model statistics, conditional effects, sensitivity analyses, and figures used in the supplementary analyses.

-   **08_Main_Figures.R:** Combines statistical outputs from the six response-specific analyses to produce summary tables, habitat-by-aridity interaction statistics, response-scale conditional effects, and initial main and supplementary figures. Exports the consolidated datasets required for final figure generation.

-   **09_Predictive_Figures.R:** Generates city-specific predictions for statistically supported habitat-by-aridity interactions (P \< 0.05) across species richness and functional diversity responses. Restricts predictions to the observed range of local habitat variation and produces prediction summaries, support assessments, and a supplementary figure preview.

-   **10_Final_Figures.R:** Combines the finalized statistical summaries and model predictions to generate Figures 2 and 3 and Supplementary Figures S1–S3. Exports publication-ready figures in PNG and PDF formats, together with input and output manifests for reproducibility.

## Data

The `Data` folder contains raw and processed datasets used to identify dryland cities, characterize bird communities, and calculate environmental predictors. Data are organized into `Raw` and `Processed` subfolders.

### Raw

Located in `Data/Raw/`.

-   **eBird:** Directory designated for the original eBird Basic Dataset (EBD) and Sampling Event Data (SED) from the July 2026 United States release. These datasets provide bird observations and checklist-level sampling information for 2015–2025. The raw archives are not included in the current folder and must be obtained separately from eBird.

-   **Shape Files/Cities:** Contains the original urban boundary polygons used to identify dryland cities. The `all_cities` dataset is provided in ESRI shapefile format (`.shp`, `.shx`, `.dbf`, `.prj`) and GeoJSON format (`.geojson`). These polygons define the spatial boundaries of urban areas used in the analysis.

-   **Shape Files/Dryland:** Contains the global UNEP-WCMC/UNCCD/CBD drylands spatial dataset (`drylands_UNCCD_CBD_july2014`), stored in ESRI shapefile format. The dataset provides dryland boundaries and classifications used to identify dryland regions and calculate the proportion of each city overlapping dryland areas.

-   **WorldClim:** Contains historical monthly climate rasters from WorldClim 2.1, downscaled using CRU-TS 4.09, at 10-arc-minute resolution. Data include minimum temperature (`tmin`), maximum temperature (`tmax`), and precipitation (`prec`), organized into the 2010–2019 and 2020–2024 periods. Files are organized by climate variable and period, with downloaded ZIP archives and extracted raster directories. Monthly climate data for 2015–2024 were used to calculate city-level annual De Martonne Aridity Index values.

### Processed

Located in `Data/Processed/`.

Contains datasets generated during the processing of spatial, eBird, and climatic data. These files are used to identify study cities, calculate biodiversity metrics, and prepare datasets for subsequent analyses.

-   **Climate:** Contains two CSV datasets generated from historical WorldClim data. `dryland_city_monthly_worldclim_climate_2015_2024.csv` provides monthly minimum and maximum temperatures, mean temperature, and precipitation for each city. `dryland_city_year_worldclim_aridity_2015_2024.csv` summarizes annual precipitation, mean annual temperature, and the calculated De Martonne Aridity Index. Records are organized by city, year, and, for monthly data, month.

-   **eBird:** Contains filtered eBird checklist and observation datasets for 2015–2025, including checklist identifiers, locality information, observation metadata, species richness estimates, and city-year summaries. The `Checklist_Richness_By_City` folder contains individual city-level CSV files with checklist information and species richness, while `City_Year_Unique_Species_By_City` contains annual unique-species summaries for each city. The `Intermediate` folder stores temporary processing files. Additional files include processing manifests, quality-control summaries, and scripts and logs associated with eBird filtering.

-   **For Analysis:** Contains processed eBird datasets combined with annual climatic information. The `Checklist_Richness_Aridity_By_City` and `City_Year_Richness_Aridity_By_City` folders contain city-specific CSV files integrating biodiversity summaries with precipitation, temperature, and De Martonne aridity measurements. Additional CSV files provide combined city-year richness summaries, climate data, eBird locality coordinates for Google Earth Engine, and quality-control information. Climate values are available for 2015–2024; records from 2025 retain missing climate values.

-   **Shape Files:** Contains processed spatial datasets in ESRI shapefile format. The `Cities` folder includes urban boundaries with calculated dryland overlap and the selected cities meeting the 50% dryland coverage threshold. The `Dryland` folder contains spatial layers representing the dryland regions of the contiguous United States, including the cleaned and classified dryland polygons used to identify study cities.

-   **Tables:** Contains CSV summaries of urban and dryland spatial intersections. `cities_with_dryland_overlap.csv` provides city-level dryland overlap information, while `dryland_cities_50pct.csv` contains the selected cities with at least 50% dryland coverage and their associated geographic metadata.

## Intermediate Data

Located in `Intermediate_data/`.

Contains intermediate datasets used to combine eBird observations, environmental predictors, and bird functional traits into the datasets required for statistical analyses. Data are organized into five subfolders.

-   **dfs_for_analysis_by_city:** Contains individual city-level datasets stored as RDS files. Files follow the naming convention `city_[ID]_[state]_[city].RDS`, allowing datasets to be identified by city and geographic location. These files represent intermediate city-specific analysis datasets prior to compilation into the combined analysis dataset.

-   **dfs_for_analysis_combined:** Contains `combined_data.RDS`, the consolidated dataset used for the statistical analyses, and an accompanying `metadata.xlsx` workbook. The combined dataset includes checklist and locality identifiers, city information, species richness, functional diversity metrics, community-weighted mean traits, sampling effort variables, aridity, and habitat cover. This dataset serves as the primary input for `01_Data_Setup.R`.

-   **functional_traits_data_AVONET:** Contains original AVONET bird functional trait datasets, taxonomic crosswalks, and processed trait information. Files include Excel workbooks and RDS objects used to harmonize AVONET species names with eBird taxonomy. `AVONET_cleaned.RDS` contains the cleaned trait dataset generated by `10_cleaning_AVONET_data.R` for subsequent functional diversity analyses.

-   **GEE_land_cover:** Contains land-cover estimates extracted from Google Earth Engine (GEE) for eBird localities using 1-km buffers. The `batch_1st_attempt`, `missed_localities`, and `missed_localities2` folders contain CSV exports from successive extraction attempts. These datasets are combined and deduplicated using `08_combining_GEE_landCover_csvs.R`, producing `full_gee_combined.RDS`, which contains the consolidated locality-level land-cover estimates.

-   **locality_id_shapefile:** Contains spatial point datasets representing eBird sampling localities, prepared for land-cover extraction in GEE. Files are stored in ESRI shapefile format (`.shp`, `.shx`, `.dbf`, `.prj`) and include locality identifiers, geographic coordinates, and batch assignments. Localities were initially divided into batches of up to 2,500 points, with additional shapefiles generated to recover missing land-cover estimates.

## Outputs

Located in `Outputs/Analysis_Resource_Importance/`.

Contains the processed analysis datasets, fitted models, statistical summaries, diagnostic results, and figures generated by the resource importance analyses. Outputs are organized into ten numbered folders corresponding to the R scripts in `Code/Analysis/Analysis Resource Importance/`.

### Data Setup

-   **01_Data_Setup:** Contains the combined analysis-ready dataset (`Data`) and associated data preparation and quality-control summaries (`Table`). These outputs provide the standardized inputs used by the response-specific statistical analyses.

### Response-Specific Analyses

Folders `02_Species_Richness` through `07_CWM_Beak_Length` contain the outputs for six response variables: species richness, functional richness (FRic), functional dispersion (FDis), functional evenness (FEve), community-weighted mean body mass, and community-weighted mean beak length.

These folders share a common structure:

-   **Data:** Contains response-specific analysis datasets and associated model or diagnostic data stored in RDS format.

-   **Models:** Contains saved fitted model objects where retained. The functional dispersion analysis differs by storing model summaries and diagnostic bundles rather than full fitted model objects.

-   **Table:** Contains CSV files summarizing model coefficients, habitat-by-aridity interactions, conditional effects, model diagnostics, sensitivity analyses, and sample characteristics.

-   **Figure:** Contains response-specific statistical figures and diagnostic visualizations, primarily exported in PNG and PDF formats.

### Synthesis and Final Figures

-   **08_Resource_Importance_Synthesis:** Combines results across the six response variables and produces consolidated statistical tables, response-scale conditional effects, and preliminary main and supplementary figures. Outputs are organized into `Data`, `Tables`, `Main_Figures`, and `Supplementary_Figures`.

-   **09_Resource_Importance_Predictions:** Contains city-specific predictions for statistically supported habitat-by-aridity interactions, including prediction datasets (`Data`), summaries of habitat support and predicted effects (`Table`), and a supplementary figure preview (`Figure`).

-   **10_Final_Figures:** Contains the final consolidated figure suite, organized into `Main_Figures` and `Supplementary_Figures`. The script produces Figures 2 and 3 and Supplementary Figures S1–S3 in PNG and PDF formats. Input and output manifest CSV files document the datasets and figure outputs used in the final visualization workflow.
