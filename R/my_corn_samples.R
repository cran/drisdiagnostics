#'Leaf Samples for PASS Analysis
#'
#' A dataset containing nutrient concentrations from leaf samples
#' used to demonstratenPASS calculation.
#'
#' @format A data frame with 3 rows and 12 columns:
#' \describe{
#'   \item{Sample}{Unique identifier for each leaf sample}
#'   \item{N}{Nitrogen concentration (%)}
#'   \item{P}{Phosphorus concentration (%)}
#'   \item{K}{Potassium concentration (%)}
#'   \item{Ca}{Calcium concentration (%)}
#'   \item{Mg}{Magnesium concentration (%)}
#'   \item{S}{Sulfur concentration (%)}
#'   \item{Zn}{Zinc concentration (mg/kg)}
#'   \item{B}{Boron concentration (mg/kg)}
#'   \item{Cu}{Copper concentration (mg/kg)}
#'   \item{Mn}{Manganese concentration (mg/kg)}
#'   \item{Fe}{Iron concentration (mg/kg)}
#' }
#'
#' @details
#' These leaf samples are used as example data for running PASS analyses.
#'
#' @examples
#' data(my_corn_samples)
#' data(ini_norms)
#' data(dni_norms)
#'
#' # Run the PASS analysis
#' pass_report <- calculate_pass_analysis(
#'   data = my_corn_samples,
#'   ini_norms = ini_norms,
#'   dni_norms = dni_norms,
#'   sample_col = "Sample"
#' )
#'
#' # View the master diagnostic compiled table
#' print(pass_report$compiled_table)
"my_corn_samples"
