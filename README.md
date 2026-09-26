# drisdiagnostics

<!-- badges: start -->
<!-- badges: end -->

The goal of drisdiagnostics is to provide tools for plant nutrient
diagnostics using the Diagnosis and Recommendation Integrated System (DRIS),
Modified DRIS (MDRIS) and the PASS approach. It supports splitting crop
populations into high- and low-yielding subgroups, calculating diagnostic
norms and indices, and visualizing nutrient balance through diverging bar
charts and tri-axial DRIS wheel plots. It also provide MDRIS and PASS indices.

## Installation


``` r
install.packages("drisdiagnostics")
```

## Example: Standard DRIS workflow

`my_crop_data` below is a placeholder for your own data frame of plant
tissue nutrient concentrations plus a yield column. Substitute it with your
own dataset before running.

``` r
library(drisdiagnostics)

# Split populations using calculated mean + a custom SD
processed_data <- split_populations(
  data = my_crop_data,
  yield_col = "Yield",
  method = "mean_sd",
  value = 2/3
)

# Calculate DRIS norms from the split populations
dris_norms <- calculate_dris_norms(
  data = processed_data,
  nutrient_cols = c("N", "P", "K", "S",  "Ca", "Mg")
)

# Generate the optimal HYP ratio reference table
hyp_ratios_results <- generate_hyp_ratios_table(
  data = processed_data,
  dris_norms = dris_norms,
  sample_col = "Sample"
)

# Run diagnostics on the sample set
diagnostics_results <- calculate_dris_diagnostics(
  data = processed_data,
  dris_norms = dris_norms,
  nutrient_cols = c("N", "P", "K", "S",  "Ca", "Mg"),
  sample_col = "Sample",
  hyp_only = TRUE
)

# Run diagnostics separately for each level of a grouping factor (e.g. leaf age)
grouped_results <- calculate_dris_by_group(
  data = processed_data,
  dris_norms = dris_norms,
  nutrient_cols = c("N", "P", "K", "S",  "Ca", "Mg"),
  filter_col = "Age",   # The column to filter and group by
  sample_col = "Sample",
  hyp_only = FALSE      # TRUE to analyze only HYP; FALSE to analyze everyone
)

# Generate the nutrient sufficiency ranges table
ranges_table <- sufficiency_ranges(
  data = processed_data,
  nutrient_cols = c("N", "P", "K", "Zn", "B")
)
```

## Visualizing results

``` r
# Diverging bar chart of nutrient indices for a single sample
indices <- c(
  N_index  = -15.4,
  P_index  = -8.2,
  K_index  = 12.6,
  Ca_index = 5.3
)

plot_dris_diverging_bar(
  indices = indices,
  sample_name = "Sample 1"
)

# Tri-axial DRIS wheel using three custom ratios
plot_dris_norms_wheel(
  dris_norms = dris_norms,
  ratios_to_plot = c("P/N", "N/K", "P/K")
)
```

## Example: MDRIS workflow

``` r
mdris_norms <- calculate_mdris_norms(
  data = my_crop_data,
  nutrient_cols = c("N", "P", "K", "Zn", "B")
)

mdris_results <- calculate_mdris_diagnostics(
  data = my_crop_data,
  mdris_norms = mdris_norms,
  nutrient_cols = c("N", "P", "K", "Zn", "B"),
  sample_col = "Sample"
)
```

## Example: PASS workflow

Unlike the workflows above, this example is fully self-contained and can be
copied and run as-is.

``` r
# Define the pre-established INI norms (critical levels & SDs)
ini_norms <- data.frame(
  nutrient = c("N", "P", "K", "S", "Zn", "Ca", "Mg", "B", "Mn", "Cu", "Fe"),
  CL       = c(25.0, 2.5, 17.5, 1.6, 19.0, 3.0, 1.6, 6.0, 19.0, 3.0, 21.0),
  SD       = c(4.159, 0.464, 4.394, 0.426, 7.490, 1.254, 0.900, 3.426, 24.296, 2.640, 74.753),
  group    = c("common", "common", "common", "common", "common",
               "rare", "rare", "rare", "rare", "rare", "rare"),
  stringsAsFactors = FALSE
)

# Define the DNI norms (pairwise ratio means & SDs)
dni_norms <- data.frame(
  ratio = c("N/P", "N/K", "N/S", "N/Zn", "P/K", "P/S", "P/Zn", "K/S", "K/Zn", "S/Zn"),
  mean  = c(9.030, 1.460, 11.900, 0.118, 0.169, 1.420, 0.0113, 8.770, 0.0714, 0.0095),
  SD    = c(2.140, 0.426, 2.610, 0.0446, 0.054, 0.3449, 0.0056, 2.990, 0.0476, 0.0036),
  stringsAsFactors = FALSE
)

# Input tissue concentrations for your samples (Macros in g/kg, Micros in mg/kg)
my_corn_samples <- data.frame(
  Sample = c("Field_1", "Field_2"),
  N  = c(22.6, 26.5),
  P  = c(3.1, 2.8),
  K  = c(16.2, 19.5),
  S  = c(2.0, 1.8),
  Zn = c(21.0, 24.0),
  Ca = c(4.4, 4.8),
  Mg = c(2.1, 1.9),
  B  = c(7.0, 8.5),
  Cu = c(3.0, 4.2),
  Mn = c(22.0, 31.0),
  Fe = c(150.0, 185.0),
  stringsAsFactors = FALSE
)

# Run the PASS analysis
pass_report <- calculate_pass_analysis(
  data = my_corn_samples,
  ini_norms = ini_norms,
  dni_norms = dni_norms,
  sample_col = "Sample"
)

pass_report
```
