#' Initial Norms for PASS Analysis
#'
#' A dataset containing initial nutrient norms used in the PASS
#' (Plant Analysis with Standardized Scores) method.
#'
#' @format A data frame with 12 rows and 4 columns:
#' \describe{
#'   \item{nutrient}{Nutrient name}
#'   \item{CL}{Critical level value}
#'   \item{SD}{Standard deviation}
#'   \item{group}{Grouping category}
#' }
#'
#' @details
#' These initial norms are taken from Baldock and Schulte (1996) and
#' are used in combination with DRIS norms to compute PASS scores.
#'
#' @examples
#' data(ini_norms)
#' head(ini_norms)
"ini_norms"
