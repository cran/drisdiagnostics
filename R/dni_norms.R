#' DRIS Norms for PASS Analysis
#'
#' A dataset containing DRIS nutrient ratio norms used in the PASS
#' method.
#'
#' @format A data frame with 11 rows and 3 columns:
#' \describe{
#'   \item{ratio}{Nutrient ratio (e.g., N/P, K/Ca)}
#'   \item{mean}{Mean value of the ratio}
#'   \item{SD}{Standard deviation of the ratio}
#' }
#'
#' @details
#' These DRIS norms are taken from Baldock and Schulte (1996) and
#' are used together with initial norms to generate PASS diagnostics.
#'
#' @examples
#' data(dni_norms)
#' head(dni_norms)
"dni_norms"
