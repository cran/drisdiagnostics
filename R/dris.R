#' Split Plant Populations into High-Yielding (HYP) and Low-Yielding (LYP) Subpopulations
#'
#' Categorizes crop sample tissue records into high-yielding (HYP) and low-yielding (LYP)
#' subpopulations based on a designated economic or yield threshold. This partitioning is a fundamental
#' preprocessing step in DRIS allowing the calculation of baseline
#' nutritional norms and standard deviation diagnostics exclusively from the healthy,
#' high-yielding reference cohort.
#'
#' @param data A data frame containing plant tissue nutrient concentrations and crop yield values or economic parameter.
#' @param yield_col Character string representing the column name of the yield variable in \code{data}.
#' @param method Character string specifying the mathematical cutoff method. Must be one of
#'   \code{"mean_sd"} (mean plus standard deviation threshold), \code{"percentile"} (quantile-based
#'   cutoff), or \code{"absolute"} (fixed yield threshold).
#' @param value Numeric threshold value matching the selected \code{method}. For \code{"mean_sd"},
#'   it represents a custom standard deviation multiplier added to the mean (defaults to 1.0 if \code{NULL}).
#'   For \code{"percentile"}, it is the percentile ranking cutoff (e.g., 75 for the 75th percentile).
#'   For \code{"absolute"}, it is the raw absolute yield value.
#'
#' @return The original data frame \code{data} with an appended factor column named
#'   \code{"population_category"} containing levels \code{"HYP"} (High-Yielding Population)
#'   and \code{"LYP"} (Low-Yielding Population).
#' @examples
#' \donttest{
#' set.seed(123)
#' my_crop_data <- data.frame(
#'   Sample = paste0("S", 1:20),
#'   Yield = c(
#'     2.0, 2.1, 2.2, 2.3, 2.4,
#'     4.0, 4.1, 4.2, 4.3, 4.4,
#'     2.0, 2.1, 2.2, 2.3, 2.4,
#'     4.0, 4.1, 4.2, 4.3, 4.4
#'   ),
#'   N  = runif(20, 2, 4),
#'   P  = runif(20, 0.2, 0.5),
#'   K  = runif(20, 1.5, 3),
#'   S  = runif(20, 0.1, 0.4),
#'   Ca = runif(20, 0.5, 2),
#'   Mg = runif(20, 0.2, 0.8)
#' )
#'
#' # Split populations using a custom percentile value
#' processed_data_per <- split_populations(
#'   data = my_crop_data,
#'   yield_col = "Yield",
#'   method = "percentile",
#'   value = 75
#' )
#'
#' # Split populations using an absolute yield threshold of 3.2
#' processed_data_abs <- split_populations(
#'   data = my_crop_data,
#'   yield_col = "Yield",
#'   method = "absolute",
#'   value = 3.2
#' )
#'
#' # Split populations using a custom standard deviation multiplier
#' processed_data <- split_populations(
#'   data = my_crop_data,
#'   yield_col = "Yield",
#'   method = "mean_sd",
#'   value = 2/3
#' )
#' }
#' @export
split_populations <- function(data, yield_col, method = c("mean_sd", "percentile", "absolute"), value = NULL) {
  method <- match.arg(method)

  if (!yield_col %in% colnames(data)) {
    stop(paste("Yield column", yield_col, "not found in the dataset."))
  }

  yields <- data[[yield_col]]

  cutoff <- switch(method,
                   "mean_sd" = {
                     mean_val <- mean(yields, na.rm = TRUE)
                     sd_val <- stats::sd(yields, na.rm = TRUE)
                     multiplier <- if (!is.null(value)) value else 1
                     mean_val + multiplier * sd_val
                   },
                   "percentile" = {
                     if (is.null(value)) {
                       stop("For 'percentile' method, you must provide a 'value' (e.g., 80 for 80th percentile).")
                     }
                     prob <- if (value > 1) value / 100 else value
                     stats::quantile(yields, probs = prob, na.rm = TRUE)
                   },
                   "absolute" = {
                     if (is.null(value)) {
                       stop("For 'absolute' method, you must provide a numeric 'value' as a cutoff threshold.")
                     }
                     value
                   }
  )

  data$population_category <- ifelse(yields >= cutoff, "HYP", "LYP")
  data$population_category <- factor(data$population_category, levels = c("HYP", "LYP"))

  return(data)
}

#' Calculate DRIS Ratio Norms and Select Optimal Expressions
#'
#' Computes descriptive statistics for both high-yielding (HYP) and low-yielding (LYP) crop
#' subpopulations across all possible bi-variate nutrient expressions, including ratio forms (A/B, B/A)
#' and product forms (A*B). To determine which expression is the most agronomically sensitive,
#' the function compares the variances of the two subpopulations using a variance ratio test (Sa/Sb)
#' and evaluates statistical significance via a standard F-test. The nutrient expression that
#' maximizes this variance ratio is flagged as optimal, establishing the baseline diagnostic norm.
#' @param data A data frame containing nutrient columns and a "population_category" column.
#' @param nutrient_cols A character vector containing the names of the nutrient columns.
#' @return A data frame containing calculated population statistics (means, variances, standard deviations,
#'   and coefficients of variation) for both HYP and LYP groups across all combinations, alongside
#'   variance ratios (Sa_Sb), F-test p-values, and a selection indicator (\code{"Optimal"})
#'   flagging the expression with the largest variance ratio.
#' @examples
#' \donttest{
#' set.seed(123)
#'
#' my_crop_data <- data.frame(
#'   Sample = paste0("S", 1:20),
#'   Yield = c(
#'     2.0, 2.1, 2.2, 2.3, 2.4,
#'     4.0, 4.1, 4.2, 4.3, 4.4,
#'     2.0, 2.1, 2.2, 2.3, 2.4,
#'     4.0, 4.1, 4.2, 4.3, 4.4
#'   ),
#'   N = runif(20, 2, 4),
#'   P = runif(20, 0.2, 0.5),
#'   K = runif(20, 1.5, 3),
#'   S = runif(20, 0.1, 0.4),
#'   Ca = runif(20, 0.5, 2),
#'   Mg = runif(20, 0.2, 0.8)
#' )
#'
#' processed_data <- split_populations(
#'   data = my_crop_data,
#'   yield_col = "Yield",
#'   method = "percentile",
#'   value = 50
#' )
#'
#' dris_norms <- calculate_dris_norms(
#'   data = processed_data,
#'   nutrient_cols = c("N", "P", "K", "S", "Ca", "Mg")
#' )
#'
#' print(dris_norms)
#' }
#' @export
calculate_dris_norms <- function(data, nutrient_cols) {
  if (!"population_category" %in% colnames(data)) {
    stop("The dataset must contain a 'population_category' column. Please run split_populations() first.")
  }

  hyp_data <- data[data$population_category == "HYP", ]
  lyp_data <- data[data$population_category == "LYP", ]

  n_hyp <- nrow(hyp_data)
  n_lyp <- nrow(lyp_data)

  if (n_hyp < 2 || n_lyp < 2) {
    stop("Both HYP and LYP must have at least 2 samples to calculate variances.")
  }

  results <- list()
  idx <- 1
  pairs <- utils::combn(nutrient_cols, 2, simplify = FALSE)

  for (pair in pairs) {
    nut_a <- pair[1]
    nut_b <- pair[2]

    expressions <- list(
      list(name = paste0(nut_a, "/", nut_b), calc = function(a, b) a / b),
      list(name = paste0(nut_b, "/", nut_a), calc = function(a, b) b / a),
      list(name = paste0(nut_a, "*", nut_b), calc = function(a, b) a * b)
    )

    pair_results <- list()
    for (expr in expressions) {
      hyp_vals <- expr$calc(hyp_data[[nut_a]], hyp_data[[nut_b]])
      lyp_vals <- expr$calc(lyp_data[[nut_a]], lyp_data[[nut_b]])

      # HYP Statistics (Sb)
      mean_b <- mean(hyp_vals, na.rm = TRUE)
      var_b  <- stats::var(hyp_vals, na.rm = TRUE)
      sd_b   <- stats::sd(hyp_vals, na.rm = TRUE)
      cv_b   <- (sd_b / mean_b) * 100

      # LYP Statistics (Sa)
      mean_a <- mean(lyp_vals, na.rm = TRUE)
      var_a  <- stats::var(lyp_vals, na.rm = TRUE)
      sd_a   <- stats::sd(lyp_vals, na.rm = TRUE)
      cv_a   <- (sd_a / mean_a) * 100

      # Variance Ratio (Sa / Sb)
      sa_sb_ratio <- var_a / var_b
      p_val <- stats::pf(sa_sb_ratio, df1 = n_lyp - 1, df2 = n_hyp - 1, lower.tail = FALSE)

      pair_results[[length(pair_results) + 1]] <- data.frame(
        Expression = expr$name,
        HYP_mean = round(mean_b, 2),
        HYP_variance = round(var_b, 2),
        HYP_SD = round(sd_b, 2),
        HYP_CV = round(cv_b, 2),
        LYP_mean = round(mean_a, 2),
        LYP_variance = round(var_a, 2),
        LYP_SD = round(sd_a, 2),
        LYP_CV = round(cv_a, 2),
        Sa_Sb = round(sa_sb_ratio, 2),
        p_value = round(p_val, 2),
        stringsAsFactors = FALSE
      )
    }

    pair_df <- do.call(rbind, pair_results)
    best_idx <- which.max(pair_df$Sa_Sb)
    pair_df$Selected <- ""
    pair_df$Selected[best_idx] <- "Optimal"

    results[[idx]] <- pair_df
    idx <- idx + 1
  }

  final_table <- do.call(rbind, results)
  rownames(final_table) <- NULL
  return(final_table)
}

#' Calculate Optimal Nutrient Ratios for the High-Yielding Population (HYP)
#'
#' Extracts the optimal diagnostic expressions identified in the norms calibration step and calculates
#' their exact ratio values for each sample in the high-yielding reference subpopulation (HYP).
#' This function provides a structured reference table summarizing baseline healthy ratios, appending
#' critical summary statistics (arithmetic means, variances, standard deviations, and coefficients of
#' variation) as summary rows at the bottom of the table.
#' @param data A data frame containing the crop tissue and population category columns.
#' @param dris_norms A data frame returned by `calculate_dris_norms()`.
#' @param sample_col Character string representing the column name of the sample identifier.
#' @return A character data frame where rows represent individual high-yielding samples and columns
#' represent the calculated optimal ratio expressions, with appended summary rows labeled
#' \code{"Mean"}, \code{"Variance"}, \code{"SD"}, and \code{"CV"} at the bottom.
#' @examples
#' \donttest{
#' set.seed(123)
#'
#' my_crop_data <- data.frame(
#'   Sample = paste0("S", 1:20),
#'   Yield = c(
#'     2.0, 2.1, 2.2, 2.3, 2.4,
#'     4.0, 4.1, 4.2, 4.3, 4.4,
#'     2.0, 2.1, 2.2, 2.3, 2.4,
#'     4.0, 4.1, 4.2, 4.3, 4.4
#'   ),
#'   N = runif(20, 2, 4),
#'   P = runif(20, 0.2, 0.5),
#'   K = runif(20, 1.5, 3),
#'   S = runif(20, 0.1, 0.4),
#'   Ca = runif(20, 0.5, 2),
#'   Mg = runif(20, 0.2, 0.8)
#' )
#'
#' processed_data <- split_populations(
#'   data = my_crop_data,
#'   yield_col = "Yield",
#'   method = "percentile",
#'   value = 50
#' )
#'
#' dris_norms <- calculate_dris_norms(
#'   data = processed_data,
#'   nutrient_cols = c("N", "P", "K", "S", "Ca", "Mg")
#' )
#'
#' hyp_ratios_results <- generate_hyp_ratios_table(
#'   data = processed_data,
#'   dris_norms = dris_norms,
#'   sample_col = "Sample"
#' )
#' print(hyp_ratios_results)
#' }
#' @export
generate_hyp_ratios_table <- function(data, dris_norms, sample_col = "Sample") {
  # Filter for HYP population
  hyp_data <- data[data$population_category == "HYP", ]
  if (nrow(hyp_data) == 0) {
    stop("No samples found in the High-Yielding Population (HYP).")
  }

  # Extract only the optimal expressions
  optimal_norms <- dris_norms[dris_norms$Selected == "Optimal", ]
  if (nrow(optimal_norms) == 0) {
    stop("No optimal expressions found in the dris_norms table. Ensure 'Selected' column has 'Optimal'.")
  }

  # Initialize the results table with the sample column
  ratios_table <- data.frame(Sample = hyp_data[[sample_col]], stringsAsFactors = FALSE)
  colnames(ratios_table) <- sample_col

  # Compute each optimal ratio for the HYP samples
  for (i in 1:nrow(optimal_norms)) {
    expr_name <- optimal_norms$Expression[i]

    # Parse and compute the expression
    if (grepl("/", expr_name)) {
      parts <- unlist(strsplit(expr_name, "/"))
      nut_num <- parts[1]
      nut_den <- parts[2]
      ratios_table[[expr_name]] <- round(hyp_data[[nut_num]] / hyp_data[[nut_den]], 2)
    } else if (grepl("\\*", expr_name)) {
      parts <- unlist(strsplit(expr_name, "\\*"))
      nut1 <- parts[1]
      nut2 <- parts[2]
      ratios_table[[expr_name]] <- round(hyp_data[[nut1]] * hyp_data[[nut2]], 2)
    }
  }

  # Calculate summary statistics for the bottom of the table
  num_cols <- colnames(ratios_table)[-1]

  means <- c("Mean", sapply(num_cols, function(col) sprintf("%.2f", mean(ratios_table[[col]], na.rm = TRUE))))
  vars  <- c("Variance", sapply(num_cols, function(col) sprintf("%.2f", stats::var(ratios_table[[col]], na.rm = TRUE))))
  sds   <- c("SD", sapply(num_cols, function(col) sprintf("%.2f", stats::sd(ratios_table[[col]], na.rm = TRUE))))
  cvs   <- c("CV", sapply(num_cols, function(col) {
    m <- mean(ratios_table[[col]], na.rm = TRUE)
    s <- stats::sd(ratios_table[[col]], na.rm = TRUE)
    sprintf("%.2f", (s / m) * 100)
  }))

  # Combine data and summaries into a single printable data frame
  summary_df <- as.data.frame(rbind(means, vars, sds, cvs), stringsAsFactors = FALSE)
  colnames(summary_df) <- colnames(ratios_table)

  ratios_table_char <- as.data.frame(lapply(ratios_table, function(x) {
    if (is.numeric(x)) sprintf("%.2f", x) else as.character(x)
  }), stringsAsFactors = FALSE)
  colnames(ratios_table_char) <- colnames(ratios_table)

  final_table <- rbind(ratios_table_char, summary_df)

  return(final_table)
}

#' Calculate DRIS Functions, Individual Indices, and NBI
#'
#' Performs standard Diagnosis and Recommendation Integrated System (DRIS) calculations for plant tissue
#' samples. It computes individual function values for each optimal ratio using the symmetrical
#' Beaufils formulas, derives individual relative nutrient indices by averaging these function
#' values with appropriate sign conventions and computes the overall Nutrient Balance Index (NBI)
#' as the sum of absolute index values to measure general nutritional balance.
#'
#' @param data A data frame containing the original plant sample dataset to be diagnosed.
#' @param dris_norms A data frame of calibrated norms returned by \code{calculate_dris_norms()}.
#' @param nutrient_cols A character vector specifying the names of the selected nutrients to evaluate.
#' @param sample_col A character string specifying the column name containing the sample identifiers.
#'   Defaults to \code{"Sample"}.
#' @param hyp_only A logical value. If \code{TRUE} (default), restricts the diagnostic calculations
#'   strictly to samples labeled \code{"HYP"}. If \code{FALSE}, runs the diagnosis over the entire
#'   dataset (both HYP and LYP).
#'
#' @return A list containing four diagnostic components:
#'   \item{nutrient_tables}{A list of individual data frames (one for each nutrient) detailing sample
#'     ratios, calculated function values, and final indices with a summary Mean row.}
#'   \item{formulas}{A list of character strings representing the final algebraic formulas used to compile
#'     the index for each nutrient.}
#'   \item{compiled_table}{A master summary data frame containing the final relative index values and
#'     overall Nutrient Balance Index (NBI) for each sample, with an overall Mean row appended at the bottom.}
#'   \item{raw_indices}{A numeric matrix containing unrounded, raw index values for subsequent statistical
#'     modeling or graphing.}
#' @examples
#' \donttest{
#' set.seed(123)
#'
#' my_crop_data <- data.frame(
#'   Sample = paste0("S", 1:20),
#'   Yield = c(
#'     2.0, 2.1, 2.2, 2.3, 2.4,
#'     4.0, 4.1, 4.2, 4.3, 4.4,
#'     2.0, 2.1, 2.2, 2.3, 2.4,
#'     4.0, 4.1, 4.2, 4.3, 4.4
#'   ),
#'   N = runif(20, 2, 4),
#'   P = runif(20, 0.2, 0.5),
#'   K = runif(20, 1.5, 3),
#'   S = runif(20, 0.1, 0.4),
#'   Ca = runif(20, 0.5, 2),
#'   Mg = runif(20, 0.2, 0.8)
#' )
#'
#' processed_data <- split_populations(
#'   data = my_crop_data,
#'   yield_col = "Yield",
#'   method = "percentile",
#'   value = 50
#' )
#'
#' dris_norms <- calculate_dris_norms(
#'   data = processed_data,
#'   nutrient_cols = c("N", "P", "K", "S", "Ca", "Mg")
#' )
#'
#' diagnostics_results <- calculate_dris_diagnostics(
#'   data = processed_data,
#'   dris_norms = dris_norms,
#'   nutrient_cols = c("N", "P", "K", "S", "Ca", "Mg"),
#'   sample_col = "Sample"
#' )
#' print(diagnostics_results$compiled_table)
#' }
#' @export
calculate_dris_diagnostics <- function(data, dris_norms, nutrient_cols, sample_col = "Sample", hyp_only = TRUE) {

  # Check if population categories have been assigned
  if (!"population_category" %in% colnames(data)) {
    stop("The dataset must contain a 'population_category' column. Please run split_populations() first.")
  }

  # Condition: Filter for HYP if hyp_only is TRUE
  if (hyp_only) {
    data <- data[data$population_category == "HYP", ]
    message("Diagnostic Condition Applied: Analyzing only the High-Yielding Population (HYP).")
  } else {
    message("Diagnostic Condition Applied: Analyzing all populations (HYP & LYP).")
  }

  if (nrow(data) == 0) {
    stop("No samples found matching the diagnostic selection.")
  }

  # Extract only optimal norms
  optimal_norms <- dris_norms[dris_norms$Selected == "Optimal", ]
  if (nrow(optimal_norms) == 0) {
    stop("No optimal expressions found in the dris_norms table. Ensure 'Selected' column has 'Optimal'.")
  }

  n_samples <- nrow(data)
  samples <- data[[sample_col]]

  # Step 6a: Calculate the f(X/Y) function values for all samples
  functions_df <- data.frame(Sample = samples, stringsAsFactors = FALSE)
  colnames(functions_df) <- sample_col

  for (i in 1:nrow(optimal_norms)) {
    expr_name <- optimal_norms$Expression[i]
    norm_mean <- optimal_norms$HYP_mean[i]
    norm_cv   <- optimal_norms$HYP_CV[i]

    # Calculate sample ratio
    if (grepl("/", expr_name)) {
      parts <- unlist(strsplit(expr_name, "/"))
      nut_num <- parts[1]
      nut_den <- parts[2]
      sample_ratio <- data[[nut_num]] / data[[nut_den]]
    } else if (grepl("\\*", expr_name)) {
      parts <- unlist(strsplit(expr_name, "\\*"))
      nut1 <- parts[1]
      nut2 <- parts[2]
      sample_ratio <- data[[nut1]] * data[[nut2]]
    }

    # Calculate f(X/Y) based on Beaufils formulas
    f_vals <- ifelse(
      sample_ratio >= norm_mean,
      ((sample_ratio / norm_mean) - 1) * (1000 / norm_cv),
      (1 - (norm_mean / sample_ratio)) * (1000 / norm_cv)
    )

    functions_df[[expr_name]] <- f_vals
  }

  # Step 6b: Create individual nutrient tables and build formulas dynamically
  nutrient_tables <- list()
  formulas_list <- list()
  indices_matrix <- matrix(0, nrow = n_samples, ncol = length(nutrient_cols))
  colnames(indices_matrix) <- paste0(nutrient_cols, "_index")

  for (j in seq_along(nutrient_cols)) {
    nut <- nutrient_cols[j]

    # Find all expressions containing this nutrient
    expr_indices <- which(sapply(optimal_norms$Expression, function(expr) {
      if (grepl("/", expr)) {
        parts <- unlist(strsplit(expr, "/"))
        return(nut %in% parts)
      } else if (grepl("\\*", expr)) {
        parts <- unlist(strsplit(expr, "\\*"))
        return(nut %in% parts)
      }
      return(FALSE)
    }))

    matching_exprs <- optimal_norms$Expression[expr_indices]
    z <- length(matching_exprs)

    # Determine the sign (+ or -) for each function value in the index equation
    signs <- sapply(matching_exprs, function(expr) {
      if (grepl("/", expr)) {
        parts <- unlist(strsplit(expr, "/"))
        if (parts[1] == nut) return("+")
        if (parts[2] == nut) return("-")
      }
      return("+") # Default fallback
    })

    # Build the algebra formula string for documentation (Corrected in v3)
    formula_terms <- paste0(signs, "f(", matching_exprs, ")")
    formula_str_terms <- paste(formula_terms, collapse = " ")
    formula_str_terms <- sub("^\\+\\s*", "", formula_str_terms)
    formula_str_terms <- gsub(" -f\\(", " - f\\(", formula_str_terms)
    formula_str_terms <- gsub(" \\+f\\(", " + f\\(", formula_str_terms)
    formula_str <- paste0(nut, " Index = { ", formula_str_terms, " } / ", z)
    formulas_list[[nut]] <- formula_str

    # Calculate the actual index values for each sample
    index_vals <- rep(0, n_samples)
    for (k in seq_along(matching_exprs)) {
      expr <- matching_exprs[k]
      sign_mult <- if (signs[k] == "+") 1 else -1
      index_vals <- index_vals + (sign_mult * functions_df[[expr]])
    }
    index_vals <- index_vals / z
    indices_matrix[, j] <- index_vals # Keep numeric for compiling

    # Assemble the individual nutrient table
    nut_table <- data.frame(Sample = samples, stringsAsFactors = FALSE)
    colnames(nut_table) <- sample_col

    # Add matching f(X/Y) columns
    for (expr in matching_exprs) {
      nut_table[[paste0("f(", expr, ")")]] <- functions_df[[expr]]
    }

    # Add final index column (rounded to 2 decimals nearest)
    index_col_name <- paste0(nut, "_index")
    nut_table[[index_col_name]] <- indices_matrix[, j]

    # Append Mean row at bottom
    col_means <- sapply(colnames(nut_table)[-1], function(col) mean(as.numeric(nut_table[[col]]), na.rm = TRUE))
    mean_row <- c("Mean", sprintf("%.2f", col_means))

    # Protect column names from mangling and force 2 decimal places on numeric columns
    nut_table_char <- as.data.frame(lapply(nut_table, function(col) {
      if (is.numeric(col)) sprintf("%.2f", col) else as.character(col)
    }), stringsAsFactors = FALSE)
    colnames(nut_table_char) <- colnames(nut_table)

    nut_table_final <- rbind(nut_table_char, mean_row)

    nutrient_tables[[nut]] <- nut_table_final
  }

  # Step 7: Compile all indices into a master summary table
  compiled_df <- data.frame(Sample = samples, stringsAsFactors = FALSE)
  colnames(compiled_df) <- sample_col

  for (nut in nutrient_cols) {
    compiled_df[[paste0(nut, "_index")]] <- indices_matrix[, paste0(nut, "_index")]
  }

  # Step 8: Calculate Nutrient Balance Index (NBI)
  nbi_vals <- apply(indices_matrix, 1, function(row) sum(abs(row)))
  compiled_df$NBI <- nbi_vals

  # Add Mean row at bottom of compiled table
  compiled_means <- sapply(colnames(compiled_df)[-1], function(col) mean(compiled_df[[col]], na.rm = TRUE))
  mean_row_compiled <- c("Mean", sprintf("%.2f", compiled_means))

  # Protect names and format to 2 decimals
  compiled_table_char <- as.data.frame(lapply(compiled_df, function(col) {
    if (is.numeric(col)) sprintf("%.2f", col) else as.character(col)
  }), stringsAsFactors = FALSE)
  colnames(compiled_table_char) <- colnames(compiled_df)

  compiled_table_final <- rbind(compiled_table_char, mean_row_compiled)

  return(list(
    nutrient_tables = nutrient_tables,
    formulas = formulas_list,
    compiled_table = compiled_table_final,
    raw_indices = indices_matrix
  ))
}

#' Perform Grouped DRIS Diagnostics and Multi-Level Factor Filtering
#'
#' Subsets a plant tissue dataset by a specified grouping factor (such as leaf age, soil type,
#' crop variety, or geographic region) and runs the standard Beaufils DRIS diagnostic workflow
#' independently for each factor level. It compiles individual diagnostic tables for each sub-group
#' and generates a master summary table containing mean nutrient indices, the overall Nutrient
#' Balance Index (NBI) computed directly from the group means, and a sorted nutrient requirement
#' order to facilitate comparative agronomic evaluations across varying environmental or management conditions.
#'
#' @param data A data frame containing crop nutrient concentration columns, the grouping factor column,
#'   and the \code{"population_category"} column (if \code{hyp_only = TRUE}) as generated by
#'   \code{split_populations()}.
#' @param dris_norms A data frame of calibrated reference norms returned by \code{calculate_dris_norms()}.
#' @param nutrient_cols A character vector specifying the exact names of the nutrient columns in \code{data}
#'   to be evaluated.
#' @param filter_col A character string representing the name of the column in \code{data} to filter
#'   and group by (e.g., \code{"Age"}, \code{"pH"}, or \code{"Variety"}).
#' @param sample_col A character string specifying the column name containing the sample identifiers
#'   in \code{data}. Defaults to \code{"Sample"}.
#' @param hyp_only A logical value. If \code{TRUE} (default), restricts the sub-group diagnostic
#'   evaluations strictly to high-yielding reference samples (\code{"HYP"}). If \code{FALSE}, runs the
#'   diagnostics over all samples in each group.
#'
#' @return A list containing two diagnostic components:
#'   \item{grouped_tables}{A named list of individual compiled diagnostic tables (data frames) for
#'     each unique level of the grouping factor.}
#'   \item{mean_summary_table}{A master summary data frame containing the mean index values for each
#'     nutrient, the overall group NBI (sum of absolute values of the mean indices), and the sorted
#'     \code{Requirement_Order} for each factor level.}
#' @examples
#' \donttest{
#' set.seed(123)
#' my_crop_data <- data.frame(
#'   Sample = paste0("S", 1:20),
#'   Age = rep(c("Young", "Medium", "Old"), length.out = 20),
#'   Yield = runif(20, 2, 5),
#'   N = runif(20, 2, 4),
#'   P = runif(20, 0.2, 0.5),
#'   K = runif(20, 1.5, 3),
#'   S = runif(20, 0.1, 0.4),
#'   Ca = runif(20, 0.5, 2),
#'   Mg = runif(20, 0.2, 0.8)
#' )
#' #  Split populations
#' processed_data <- split_populations(
#'   data = my_crop_data,
#'   yield_col = "Yield",
#'   method = "percentile",
#'   value = 50
#' )
#'
#' #  Calculate DRIS norms
#' dris_norms <- calculate_dris_norms(
#'   data = processed_data,
#'   nutrient_cols = c("N", "P", "K", "S", "Ca", "Mg")
#' )
#'
#' #  Grouped DRIS diagnostics
#' grouped_results <- calculate_dris_by_group(
#'   data = processed_data,
#'   dris_norms = dris_norms,
#'   nutrient_cols = c("N", "P", "K", "S", "Ca", "Mg"),
#'   filter_col = "Age",
#'   sample_col = "Sample",
#'   hyp_only = TRUE
#' )
#' }
#' @export
calculate_dris_by_group <- function(data, dris_norms, nutrient_cols, filter_col, sample_col = "Sample", hyp_only = TRUE) {
  if (!filter_col %in% colnames(data)) {
    stop(paste("Filter column", filter_col, "not found in the dataset."))
  }

  # Extract unique, non-missing levels of the grouping factor
  levels_list <- unique(data[[filter_col]])
  levels_list <- levels_list[!is.na(levels_list)]

  grouped_tables <- list()
  summary_rows <- list()

  for (level in levels_list) {
    # Filter dataset for the specific level (e.g., "6-7 old")
    sub_data <- data[data[[filter_col]] == level, ]

    # Run diagnostics on this subset
    diag_res <- calculate_dris_diagnostics(
      data = sub_data,
      dris_norms = dris_norms,
      nutrient_cols = nutrient_cols,
      sample_col = sample_col,
      hyp_only = hyp_only
    )

    # Store individual compiled tables
    grouped_tables[[as.character(level)]] <- diag_res$compiled_table

    # FIX: Extract strictly numeric raw indices to calculate means cleanly (prevents NA values)
    raw_indices <- diag_res$raw_indices
    group_means <- colMeans(raw_indices, na.rm = TRUE)
    names(group_means) <- nutrient_cols

    # REQ FIX: Calculate NBI directly from the row means themselves (Sum of absolute values of mean indices)
    group_nbi <- sum(abs(group_means))

    # Sort indices in ascending order (most negative to most positive)
    sorted_idx <- order(group_means)
    sorted_nuts <- nutrient_cols[sorted_idx]
    sorted_vals <- group_means[sorted_idx]

    # Build requirement order string with '=' or '>' signs
    req_order_str <- sorted_nuts[1]
    if (length(sorted_nuts) > 1) {
      for (k in 2:length(sorted_nuts)) {
        sep <- if (abs(sorted_vals[k] - sorted_vals[k-1]) < 1e-5) " = " else " > "
        req_order_str <- paste0(req_order_str, sep, sorted_nuts[k])
      }
    }

    # Assemble row data
    row_data <- data.frame(Factor_Level = as.character(level), stringsAsFactors = FALSE)
    for (nut in nutrient_cols) {
      row_data[[paste0("Mean_", nut, "_index")]] <- group_means[nut]
    }
    row_data$NBI <- group_nbi
    row_data$Requirement_Order <- req_order_str

    summary_rows[[as.character(level)]] <- row_data
  }

  # Bind all rows into the final summary table
  mean_summary_table <- do.call(rbind, summary_rows)
  rownames(mean_summary_table) <- NULL

  # Rename first column to match user's factor column
  colnames(mean_summary_table)[1] <- filter_col

  # Format all numeric columns to exactly 2 decimal places nearest
  num_cols <- sapply(mean_summary_table, is.numeric)
  mean_summary_table[num_cols] <- lapply(mean_summary_table[num_cols], function(col) {
    round(col, 2)
  })

  return(list(
    grouped_tables = grouped_tables,
    mean_summary_table = mean_summary_table
  ))
}



#' Create Sufficiency Ranges Table for Nutrients based on HYP
#'
#' This function calculates the sufficiency ranges for each nutrient in the High-Yielding
#' Population (HYP) based on standard deviation thresholds (4/3 SD and 8/3 SD).
#'
#' Classification rules:
#' - Deficient: < mean - (8/3) * SD
#' - Low: mean - (8/3) * SD to mean - (4/3) * SD
#' - Optimum: mean - (4/3) * SD to mean + (4/3) * SD
#' - High: mean + (4/3) * SD to mean + (8/3) * SD
#' - Excessive: > mean + (8/3) * SD
#'
#' @param data A data frame containing the crop tissue and population category columns (after running split_populations()).
#' @param nutrient_cols A character vector containing the names of the nutrient columns.
#' @return A data frame containing calculated sufficiency ranges for each nutrient with columns:
#'   Nutrient, Mean, SD, Deficient, low, optimum, high and excessive values.
#' @examples
#' \donttest{
#' set.seed(123)
#' my_crop_data <- data.frame(
#'   Sample = paste0("S", 1:20),
#'   Yield = c(
#'     2.0, 2.1, 2.2, 2.3, 2.4,
#'     4.0, 4.1, 4.2, 4.3, 4.4,
#'     2.0, 2.1, 2.2, 2.3, 2.4,
#'     4.0, 4.1, 4.2, 4.3, 4.4
#'   ),
#'   N  = runif(20, 2, 4),
#'   P  = runif(20, 0.2, 0.5),
#'   K  = runif(20, 1.5, 3),
#'   Zn = runif(20, 0.001, 0.01),
#'   B  = runif(20, 0.0005, 0.002)
#' )
#'
#' processed_data <- split_populations(
#'   data = my_crop_data,
#'   yield_col = "Yield",
#'   method = "percentile",
#'   value = 50
#' )
#'
#' # Generate the nutrient sufficiency ranges table
#' ranges_table <- sufficiency_ranges(
#'   data = processed_data,
#'   nutrient_cols = c("N", "P", "K", "Zn", "B")
#' )
#' print(ranges_table)
#' }
#' @export
sufficiency_ranges <- function(data, nutrient_cols) {
  # Check if population categories have been assigned
  if (!"population_category" %in% colnames(data)) {
    stop("The dataset must contain a 'population_category' column. Please run split_populations() first.")
  }

  # Filter for High-Yielding Population (HYP)
  hyp_data <- data[data$population_category == "HYP", ]
  if (nrow(hyp_data) == 0) {
    stop("No samples found in the High-Yielding Population (HYP).")
  }

  results <- data.frame(
    Nutrient = character(),
    Mean = numeric(),
    SD = numeric(),
    Deficient = character(),
    low = character(),
    optimum = character(),
    high = character(),
    excessive = character(),
    stringsAsFactors = FALSE
  )

  for (nut in nutrient_cols) {
    if (!nut %in% colnames(hyp_data)) {
      stop(paste("Nutrient column", nut, "not found in the dataset."))
    }

    vals <- as.numeric(hyp_data[[nut]])
    mean_val <- mean(vals, na.rm = TRUE)
    sd_val <- stats::sd(vals, na.rm = TRUE)

    # Calculate limits
    t1 <- mean_val - (8/3) * sd_val
    t2 <- mean_val - (4/3) * sd_val
    t3 <- mean_val + (4/3) * sd_val
    t4 <- mean_val + (8/3) * sd_val

    # Format to 2 decimal places
    mean_fmt <- round(mean_val, 2)
    sd_fmt <- round(sd_val, 2)
    t1_fmt <- round(t1, 2)
    t2_fmt <- round(t2, 2)
    t3_fmt <- round(t3, 2)
    t4_fmt <- round(t4, 2)

    deficient_str <- paste0("< ", sprintf("%.2f", t1_fmt))
    low_str <- paste0(sprintf("%.2f", t1_fmt), " - ", sprintf("%.2f", t2_fmt))
    opt_str <- paste0(sprintf("%.2f", t2_fmt), " - ", sprintf("%.2f", t3_fmt))
    high_str <- paste0(sprintf("%.2f", t3_fmt), " - ", sprintf("%.2f", t4_fmt))
    excess_str <- paste0("> ", sprintf("%.2f", t4_fmt))

    results <- rbind(results, data.frame(
      Nutrient = nut,
      Mean = mean_fmt,
      SD = sd_fmt,
      Deficient = deficient_str,
      low = low_str,
      optimum = opt_str,
      high = high_str,
      excessive = excess_str,
      stringsAsFactors = FALSE
    ))
  }

  return(results)
}


# ==============================================================================
# VISUALIZATION FUNCTIONS (ggplot2)
# ==============================================================================

#' Plot DRIS Indices as a Diverging Bar Chart
#'
#' Visualizes DRIS index values using a horizontal diverging bar chart.
#' Negative index values indicate relative nutrient deficiency, whereas
#' positive index values indicate relative nutrient excess.
#'
#' @param indices A named numeric vector containing DRIS index values.
#'   The names of the vector should correspond to nutrient names. Names ending
#'   in \code{"_index"} are automatically simplified for plotting.
#'
#' @param sample_name A character string specifying the sample name to display
#'   in the plot title. Defaults to \code{"Sample"}.
#'
#' @return A \code{ggplot} object representing a horizontal diverging bar chart
#'   of DRIS indices, with nutrients classified as deficient (negative values)
#'   or excessive (positive values).
#'
#' @import ggplot2
#'
#' @examples
#' \donttest{
#' indices <- c(
#'   N_index = -15.4,
#'   P_index = -8.2,
#'   K_index = 12.6,
#'   Ca_index = 5.3
#' )
#'
#' plot_dris_diverging_bar(
#'   indices = indices,
#'   sample_name = "Sample 1"
#' )
#' }
#' @export
plot_dris_diverging_bar <- function(indices, sample_name = "Sample") {

  # Convert named vector to data frame
  df <- data.frame(
    Nutrient = names(indices),
    Index = as.numeric(indices),
    stringsAsFactors = FALSE
  )

  # Strip "_index" suffix if present
  df$Nutrient <- gsub("_index", "", df$Nutrient)

  # Sort nutrients from most deficient to most excessive
  df$Nutrient <- factor(
    df$Nutrient,
    levels = df$Nutrient[order(df$Index)]
  )

  # Classify nutrient status
  df$Status <- ifelse(
    df$Index < 0,
    "Deficient",
    "Excess"
  )

  # Create labels
  df$Index_Label <- sprintf("%.2f", df$Index)

  df$Label_Hjust <- ifelse(
    df$Index >= 0,
    -0.25,
    1.25
  )

  # Create symbols for ggplot aesthetics
  index_var <- rlang::sym("Index")
  nutrient_var <- rlang::sym("Nutrient")
  status_var <- rlang::sym("Status")
  label_var <- rlang::sym("Index_Label")
  hjust_var <- rlang::sym("Label_Hjust")

  # Build diverging bar chart
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = !!index_var,
      y = !!nutrient_var,
      fill = !!status_var
    )
  ) +

    # Draw bars
    ggplot2::geom_col(
      width = 0.6,
      color = "black",
      linewidth = 0.3,
      alpha = 0.85
    ) +

    # Central reference line
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey30",
      linewidth = 0.8
    ) +

    # Index value labels
    ggplot2::geom_text(
      ggplot2::aes(
        label = !!label_var,
        hjust = !!hjust_var
      ),
      size = 4,
      fontface = "bold"
    ) +

    # Custom colors
    ggplot2::scale_fill_manual(
      values = c(
        "Deficient" = "#e41a1c",
        "Excess" = "#377eb8"
      )
    ) +

    # Axis expansion
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.18, 0.18))
    ) +

    ggplot2::theme_minimal(base_size = 14) +

    ggplot2::labs(
      title = paste("Diverging DRIS Indices for", sample_name),
      subtitle =
        "Negative values indicate deficiency | Positive values indicate excess",
      x = "DRIS Index Value",
      y = "Nutrient Category"
    ) +

    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none",

      plot.title = ggplot2::element_text(
        face = "bold",
        size = 16,
        margin = ggplot2::margin(b = 5)
      ),

      plot.subtitle = ggplot2::element_text(
        color = "grey40",
        size = 11,
        margin = ggplot2::margin(b = 15)
      ),

      axis.title.x = ggplot2::element_text(
        face = "bold",
        margin = ggplot2::margin(t = 10)
      ),

      axis.title.y = ggplot2::element_text(
        face = "bold",
        margin = ggplot2::margin(r = 10)
      ),

      axis.text.y = ggplot2::element_text(
        face = "bold",
        size = 12
      )
    )

  return(p)
}

# Declare global variables used in ggplot2 aesthetics
utils::globalVariables(
  c(
    "Index",
    "Nutrient",
    "Status",
    "x",
    "y",
    "xend",
    "yend",
    "x_label",
    "y_label",
    "Label",
    "Text",
    "Color",
    "Size"
  )
)

#' Plot Tri-Axial DRIS Wheel Chart from Calibration Norms
#'
#' Projects the calculated Z-score deviations of the low-yielding population (LYP) relative to
#' the high-yielding population (HYP) norms on a classic concentric tri-axial balance wheel diagram.
#' The function extracts descriptive statistics for exactly three selected ratio expressions, projects
#' their Cartesian coordinate points at 120-degree radial intervals. The visual boundaries (representing the balanced zone at 2/3 SD radius and the
#' imbalanced zone at 4/3 SD radius) allow agronomists to instantly diagnose nutrient balance and predict
#' qualitative deficiency or excess dynamics.
#'
#' @param dris_norms A data frame of calibrated reference norms and statistics returned by
#'   \code{calculate_dris_norms()}.
#' @param ratios_to_plot A character vector of exactly three ratio expressions present in
#'   \code{dris_norms$Expression} (e.g., \code{c("N/K", "Ca/K", "Ca/N")}).
#'
#' @return A \code{ggplot} object representing a fully annotated, concentric tri-axial
#'   DRIS balance wheel with directional nutrient deficiency and excess indicators.
#' @import ggplot2
#' @examples
#' \donttest{
#' set.seed(123)
#' my_crop_data <- data.frame(
#'   Sample = paste0("S", 1:20),
#'   Yield = c(
#'     2.0, 2.1, 2.2, 2.3, 2.4,
#'     4.0, 4.1, 4.2, 4.3, 4.4,
#'     2.0, 2.1, 2.2, 2.3, 2.4,
#'     4.0, 4.1, 4.2, 4.3, 4.4
#'   ),
#'   N  = runif(20, 2, 4),
#'   P  = runif(20, 0.2, 0.5),
#'   K  = runif(20, 1.5, 3),
#'   S  = runif(20, 0.1, 0.4),
#'   Ca = runif(20, 0.5, 2),
#'   Mg = runif(20, 0.2, 0.8)
#' )
#'
#' processed_data <- split_populations(
#'   data = my_crop_data,
#'   yield_col = "Yield",
#'   method = "percentile",
#'   value = 50
#' )
#'
#' dris_norms <- calculate_dris_norms(
#'   data = processed_data,
#'   nutrient_cols = c("N", "P", "K", "S", "Ca", "Mg")
#' )
#'
#' plot_dris_norms_wheel(
#'   dris_norms = dris_norms,
#'   ratios_to_plot = c("P/N", "N/K", "P/K")
#' )
#' }
#' @export
plot_dris_norms_wheel <- function(dris_norms, ratios_to_plot) {
  # Internal helper to generate circle boundary coordinates
  make_circle <- function(r, x_center = 0, y_center = 0, n_points = 100) {
    theta <- seq(0, 2 * pi, length.out = n_points)
    data.frame(
      x = x_center + r * cos(theta),
      y = y_center + r * sin(theta)
    )
  }


  # 1. Validation Checks
  if (length(ratios_to_plot) != 3) {
    stop("You must specify exactly three ratio expressions to plot.")
  }

  missing_ratios <- setdiff(ratios_to_plot, dris_norms$Expression)
  if (length(missing_ratios) > 0) {
    stop(paste("The following selected ratios were not found in dris_norms$Expression:",
               paste(missing_ratios, collapse = ", ")))
  }

  selected_norms <- dris_norms[dris_norms$Expression %in% ratios_to_plot, ]
  selected_norms <- selected_norms[match(ratios_to_plot, selected_norms$Expression), ]

  # 2. Extract statistics and compute HYP_SD & Z-scores
  selected_norms$HYP_SD <- (selected_norms$HYP_mean * selected_norms$HYP_CV) / 100

  # Calculate Z-score for the low-yielding population (LYP)
  selected_norms$Z <- (selected_norms$LYP_mean - selected_norms$HYP_mean) / selected_norms$HYP_SD

  # 3. Project three axes at 120-degree intervals (90, 30, 150 degrees in radians)
  angles <- c(90, 30, 150) * (pi / 180)

  # Convert polar coordinates (Z, angle) to Cartesian coordinates (X, Y) for plotting
  vertices <- data.frame(
    Ratio = selected_norms$Expression,
    Z = selected_norms$Z,
    x = selected_norms$Z * cos(angles),
    y = selected_norms$Z * sin(angles),
    stringsAsFactors = FALSE
  )

  # 4. Define concentric circle boundaries (Inner = 2/3 SD, Outer = 4/3 SD)
  r_inner <- 2/3
  r_outer <- 4/3
  r_max   <- 3.0   # Fixed Maximum Axis limit to maintain aesthetic scale

  inner_circle <- make_circle(r_inner)
  outer_circle <- make_circle(r_outer)

  # Set up the intersecting axes lines
  axes_df <- data.frame(
    x    = 0,
    y    = 0,
    xend = r_max * cos(angles),
    yend = r_max * sin(angles),
    Label = ratios_to_plot
  )

  # Opposing lines to complete the full diameters across the chart
  opp_axes_df <- data.frame(
    x    = 0,
    y    = 0,
    xend = r_max * cos(angles + pi),
    yend = r_max * sin(angles + pi)
  )

  # 5. Dynamic Tick-Mark Generation Along Axes
  tick_steps <- c(-2, -1, 1, 2)
  ticks_list <- list()

  for (i in 1:3) {
    mean_val <- selected_norms$HYP_mean[i]
    sd_val   <- selected_norms$HYP_SD[i]
    angle    <- angles[i]

    # Calculate values at each step
    step_values <- mean_val + (tick_steps * sd_val)

    # Project positions in Cartesian space
    tick_df <- data.frame(
      Z_score = tick_steps,
      Value   = step_values,
      x       = tick_steps * cos(angle),
      y       = tick_steps * sin(angle),
      stringsAsFactors = FALSE
    )

    # Offset labels to prevent overlapping the axis lines
    if (i == 1) { # Vertical Axis
      tick_df$x_label <- tick_df$x + 0.18
      tick_df$y_label <- tick_df$y
    } else if (i == 2) { # Top-Right Axis
      tick_df$x_label <- tick_df$x + 0.15 * cos(120 * pi / 180)
      tick_df$y_label <- tick_df$y + 0.15 * sin(120 * pi / 180)
    } else { # Top-Left Axis
      tick_df$x_label <- tick_df$x + 0.15 * cos(60 * pi / 180)
      tick_df$y_label <- tick_df$y + 0.15 * sin(60 * pi / 180)
    }

    # Format labels to 2 decimal places
    tick_df$Label <- sprintf("%.2f", tick_df$Value)
    ticks_list[[i]] = tick_df
  }

  ticks_df <- do.call(rbind, ticks_list)

  # 6. Build Qualitative Annotation Map (Arrows) - Corrected vector indexing
  el1 <- unlist(strsplit(ratios_to_plot[1], "/"))
  el2 <- unlist(strsplit(ratios_to_plot[2], "/"))
  el3 <- unlist(strsplit(ratios_to_plot[3], "/"))

  labels_df <- data.frame(
    x = c(
      # Center Balanced Sector (Optimum Zone - Inside Inner Circle)
      -0.30, 0.30,
      0.35, -0.15,
      -0.35, 0.15,

      # Transition Zone (Ring - Slight Imbalance)
      0.50, 1.10, -1.10,

      # Extreme Outer Zone (Imbalance - Deficient/Excess Arrows)
      -1.40, 1.40,
      2.40, -2.40,
      -2.40, 2.40
    ),
    y = c(
      # Center Balanced
      0.40, 0.40,
      -0.20, -0.40,
      -0.20, -0.40,

      # Transition Zone
      1.40, -0.60, -0.60,

      # Extreme Outer
      2.30, -2.30,
      1.30, -1.30,
      1.30, -1.30
    ),
    Text = c(
      # Center Balanced (Element)
      paste0("-> ", el1[1]), paste0("-> ", el1[2]),
      paste0("-> ", el2[1]), paste0("-> ", el2[2]),
      paste0("-> ", el3[1]), paste0("-> ", el3[2]),

      # Transition (Diagonal Arrows)
      paste0(el1[1], " \\"), paste0(el1[2], " \\"), paste0(el2[2], " \\"),

      # Extreme Outer (Vertical & Diagonal Imbalance Arrows)
      paste0("^ ", el1[1]), paste0("v ", el1[1]),
      paste0("^ ", el2[1]), paste0("v ", el2[1]),
      paste0("^ ", el3[1]), paste0("v ", el3[1])
    ),
    Color = c(
      # Center = Green (Balanced)
      rep("#2ca02c", 6),
      # Transition = Orange (Slight Imbalance)
      rep("#ff7f0e", 3),
      # Outer = Red (Severe Imbalance)
      rep("#d62728", 6)
    ),
    Size = c(
      rep(4.0, 6),
      rep(4.5, 3),
      rep(5.0, 6)
    )
  )

  # 7. Build the ggplot2 visual wheel
  p <- ggplot() +

    # Concentric circles representing Zone of Balance and Imbalance
    geom_polygon(
      data = outer_circle,
      aes(x, y),
      fill = "#fff5f0",
      color = "#fc9272",
      linetype = "solid",
      linewidth = 0.8
    ) +

    geom_polygon(
      data = inner_circle,
      aes(x, y),
      fill = "#f7fcf5",
      color = "#a1d99b",
      linetype = "solid",
      linewidth = 0.8
    ) +

    # Central origin point
    geom_point(
      aes(x = 0, y = 0),
      color = "grey30",
      size = 2
    ) +

    # Grid axes
    geom_segment(
      data = axes_df,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      arrow = arrow(length = unit(0.25, "cm")),
      color = "black",
      linewidth = 0.7
    ) +

    geom_segment(
      data = opp_axes_df,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      color = "black",
      linewidth = 0.5,
      linetype = "solid"
    ) +

    # Axis labels
    annotate(
      "text",
      x = 0,
      y = r_max + 0.20,
      label = ratios_to_plot[1],
      size = 5.5,
      fontface = "bold"
    ) +

    annotate(
      "text",
      x = r_max * cos(30 * pi / 180) + 0.25,
      y = r_max * sin(30 * pi / 180) + 0.15,
      label = ratios_to_plot[2],
      size = 5.5,
      fontface = "bold"
    ) +

    annotate(
      "text",
      x = r_max * cos(150 * pi / 180) - 0.25,
      y = r_max * sin(150 * pi / 180) + 0.15,
      label = ratios_to_plot[3],
      size = 5.5,
      fontface = "bold"
    ) +

    # Axis ticks
    geom_text(
      data = ticks_df,
      aes(
        x = x_label,
        y = y_label,
        label = Label
      ),
      size = 3.5,
      fontface = "bold",
      color = "grey20"
    ) +

    # Segment labels
    geom_text(
      data = labels_df,
      aes(
        x = x,
        y = y,
        label = Text,
        color = Color,
        size = Size
      ),
      fontface = "bold"
    ) +

    # Coordinate cleanup
    scale_color_identity() +
    scale_size_identity() +

    coord_fixed(
      xlim = c(-3.4, 3.4),
      ylim = c(-3.4, 3.4)
    ) +

    theme_void() +

    labs(
      title = "Tri-Axial DRIS Physiological Diagnosis Chart",
      subtitle = paste(
        "Plotting LYP deviation Z-scores relative to HYP norms\n",
        "Inner Circle =", round(r_inner, 2),
        "SD | Outer Circle =", round(r_outer, 2), "SD"
      )
    ) +

    theme(
      plot.title = element_text(
        face = "bold",
        size = 16,
        hjust = 0.5,
        margin = margin(b = 5)
      ),
      plot.subtitle = element_text(
        color = "grey40",
        size = 10,
        hjust = 0.5,
        margin = margin(b = 20)
      )
    )

  return(p)
}

#' Calculate Beverly's Modified DRIS (MDRIS) Norms
#'
#' This function computes both log-concentration norms and log-ratio norms for the entire
#' population as proposed by Beverly (1987). Unlike traditional DRIS, MDRIS norms are computed
#' for the entire population database without splitting into high/low-yielding groups, as the
#' log-transformation stabilizes variances across different yield levels.
#'
#' @param data A data frame containing raw nutrient concentrations (dry matter percentage).
#' @param nutrient_cols A character vector containing the names of the nutrient columns.
#' @return A list containing:
#'   1. `nutrient_norms`: Data frame of individual nutrient statistics (Raw Mean, Raw SD, Log Mean, Log SD).
#'   2. `ratio_norms`: Data frame of log-transformed ratio statistics (Expression, Log_Mean, Log_SD).
#' @examples
#' \donttest{
#' set.seed(123)
#' my_crop_data <- data.frame(
#'   Sample = paste0("S", 1:20),
#'   N  = runif(20, 2, 4),
#'   P  = runif(20, 0.2, 0.5),
#'   K  = runif(20, 1.5, 3),
#'   Zn = runif(20, 0.001, 0.01),
#'   B  = runif(20, 0.0005, 0.002)
#' )
#'
#' mdris_norms <- calculate_mdris_norms(
#'   data = my_crop_data,
#'   nutrient_cols = c("N", "P", "K", "Zn", "B")
#' )
#' print(mdris_norms)
#' }
#' @export
calculate_mdris_norms <- function(data, nutrient_cols) {
  # 1. Individual Nutrient Norms
  nut_norms_list <- list()
  for (nut in nutrient_cols) {
    raw_vals <- as.numeric(data[[nut]])
    log_vals <- log(raw_vals)

    nut_norms_list[[nut]] <- data.frame(
      Nutrient = nut,
      Raw_Mean = round(mean(raw_vals, na.rm = TRUE), 4),
      Raw_SD   = round(stats::sd(raw_vals, na.rm = TRUE), 4),
      Log_Mean = round(mean(log_vals, na.rm = TRUE), 4),
      Log_SD   = round(stats::sd(log_vals, na.rm = TRUE), 4),
      stringsAsFactors = FALSE
    )
  }
  nut_norms_df <- do.call(rbind, nut_norms_list)
  rownames(nut_norms_df) <- NULL

  # 2. Log-Ratio Norms (for all possible pairs)
  ratio_norms_list <- list()
  idx <- 1
  pairs <- utils::combn(nutrient_cols, 2, simplify = FALSE)

  for (pair in pairs) {
    nut_a <- pair[1]
    nut_b <- pair[2]

    # Calculate log(A/B)
    log_ratio_vals <- log(data[[nut_a]] / data[[nut_b]])

    ratio_norms_list[[idx]] <- data.frame(
      Expression = paste0(nut_a, "/", nut_b),
      Log_Mean   = round(mean(log_ratio_vals, na.rm = TRUE), 4),
      Log_SD     = round(stats::sd(log_ratio_vals, na.rm = TRUE), 4),
      stringsAsFactors = FALSE
    )
    idx <- idx + 1
  }
  ratio_norms_df <- do.call(rbind, ratio_norms_list)
  rownames(ratio_norms_df) <- NULL

  return(list(
    nutrient_norms = nut_norms_df,
    ratio_norms = ratio_norms_df
  ))
}

#' Calculate Beverly's Modified DRIS (MDRIS) Diagnostics & Concentration Indices
#'
#' This function takes a tissue sample dataset and pre-calculated MDRIS norms, and computes:
#' 1. Ratio-Based MDRIS indices (I_R) using the unified symmetrical logarithmic equation.
#' 2. Log-Concentration Indices (J_log_r) evaluating nutrients on an absolute log scale.
#' 3. Direct Concentration Indices (K_r) evaluating nutrients directly on raw scales.
#'
#' @param data A data frame of plant tissue analysis to diagnose.
#' @param mdris_norms The list returned by `calculate_mdris_norms()`.
#' @param nutrient_cols A character vector containing the names of the selected nutrients.
#' @param sample_col Character string representing the column name of the sample identifier.
#' @return A list containing:
#'   1. `weighted_deviations`: A data frame of calculated WDR/S values for each sample.
#'   2. `compiled_table`: A master data frame containing MDRIS indices, J_log_r values, K_r values and corresponding diagnostic rankings.
#' @examples
#' \donttest{
#' set.seed(123)
#' my_crop_data <- data.frame(
#'   Sample = paste0("S", 1:20),
#'   N  = runif(20, 2, 4),
#'   P  = runif(20, 0.2, 0.5),
#'   K  = runif(20, 1.5, 3),
#'   Zn = runif(20, 0.001, 0.01),
#'   B  = runif(20, 0.0005, 0.002)
#' )
#'
#' mdris_norms <- calculate_mdris_norms(
#'   data = my_crop_data,
#'   nutrient_cols = c("N", "P", "K", "Zn", "B")
#' )
#'
#' mdris_results <- calculate_mdris_diagnostics(
#'   data = my_crop_data,
#'   mdris_norms = mdris_norms,
#'   nutrient_cols = c("N", "P", "K", "Zn", "B"),
#'   sample_col = "Sample"
#' )
#' print(mdris_results$compiled_table)
#' }
#' @export
calculate_mdris_diagnostics <- function(data, mdris_norms, nutrient_cols, sample_col = "Sample") {
  n_samples <- nrow(data)
  samples <- data[[sample_col]]

  nutrient_norms <- mdris_norms$nutrient_norms
  ratio_norms <- mdris_norms$ratio_norms

  # Step 1: Calculate Weighted Deviations (WDR/S)
  wdr_df <- data.frame(Sample = samples, stringsAsFactors = FALSE)
  colnames(wdr_df) <- sample_col

  for (i in 1:nrow(ratio_norms)) {
    expr <- ratio_norms$Expression[i]
    norm_mean <- ratio_norms$Log_Mean[i]
    norm_sd   <- ratio_norms$Log_SD[i]

    parts <- unlist(strsplit(expr, "/"))
    nut_num <- parts[1]
    nut_den <- parts[2]

    # Calculate sample log-ratio
    sample_log_ratio <- log(data[[nut_num]] / data[[nut_den]])

    # Weighted Symmetrical Deviation Formula (Equation 1 of Beverly 1987)
    wdr_vals <- 10 * (sample_log_ratio - norm_mean) / norm_sd
    wdr_df[[expr]] <- wdr_vals
  }

  # Step 2: Calculate Ratio-Based Modified DRIS Indices (I_R)
  indices_matrix <- matrix(0, nrow = n_samples, ncol = length(nutrient_cols))
  colnames(indices_matrix) <- paste0("I_", nutrient_cols)

  for (j in seq_along(nutrient_cols)) {
    nut <- nutrient_cols[j]

    # Find all ratio norms containing this nutrient
    matching_ratios_idx <- which(sapply(ratio_norms$Expression, function(expr) {
      parts <- unlist(strsplit(expr, "/"))
      return(nut %in% parts)
    }))

    matching_ratios <- ratio_norms$Expression[matching_ratios_idx]
    z <- length(matching_ratios)

    if (z == 0) {
      stop(paste("No ratio norms found involving nutrient:", nut))
    }

    # Accumulate weighted deviations with correct signs
    index_vals <- rep(0, n_samples)
    for (expr in matching_ratios) {
      parts <- unlist(strsplit(expr, "/"))
      sign_mult <- if (parts[1] == nut) 1 else -1
      index_vals <- index_vals + (sign_mult * wdr_df[[expr]])
    }

    # Average the deviations (Equation 2 of Beverly 1987)
    indices_matrix[, j] <- index_vals / z
  }

  # Step 3: Calculate Log-Concentration Indices (J_log_r)
  j_matrix <- matrix(0, nrow = n_samples, ncol = length(nutrient_cols))
  colnames(j_matrix) <- paste0("J_log_", nutrient_cols)

  for (j in seq_along(nutrient_cols)) {
    nut <- nutrient_cols[j]
    norm_row <- nutrient_norms[nutrient_norms$Nutrient == nut, ]

    log_mean <- norm_row$Log_Mean
    log_sd   <- norm_row$Log_SD

    sample_log_conc <- log(as.numeric(data[[nut]]))

    # J_log_r Formula (Equation 3 of Beverly 1987)
    j_matrix[, j] <- (sample_log_conc - log_mean) / log_sd
  }

  # Step 4: Calculate Direct Concentration Indices (K_r)
  k_matrix <- matrix(0, nrow = n_samples, ncol = length(nutrient_cols))
  colnames(k_matrix) <- paste0("K_", nutrient_cols)

  for (j in seq_along(nutrient_cols)) {
    nut <- nutrient_cols[j]
    norm_row <- nutrient_norms[nutrient_norms$Nutrient == nut, ]

    raw_mean <- norm_row$Raw_Mean
    raw_sd   <- norm_row$Raw_SD

    sample_raw_conc <- as.numeric(data[[nut]])

    # K_r Formula (Equation 4 of Beverly 1987)
    k_matrix[, j] <- (sample_raw_conc - raw_mean) / raw_sd
  }

  # Step 5: Assemble Compiled Master Table
  compiled_df <- data.frame(Sample = samples, stringsAsFactors = FALSE)
  colnames(compiled_df) <- sample_col

  # Append I_R columns
  for (nut in nutrient_cols) {
    compiled_df[[paste0("I_", nut)]] <- round(indices_matrix[, paste0("I_", nut)], 2)
  }

  # Append J_log_r columns
  for (nut in nutrient_cols) {
    compiled_df[[paste0("J_log_", nut)]] <- round(j_matrix[, paste0("J_log_", nut)], 2)
  }

  # Append K_r columns
  for (nut in nutrient_cols) {
    compiled_df[[paste0("K_", nut)]] <- round(k_matrix[, paste0("K_", nut)], 2)
  }

  # Add Requirement Order based on Ratio Indices I_R
  req_orders <- sapply(1:n_samples, function(i) {
    row_vals <- indices_matrix[i, ]
    names(row_vals) <- nutrient_cols

    # Sort in ascending order
    sorted_idx <- order(row_vals)
    sorted_nuts <- nutrient_cols[sorted_idx]
    sorted_vals <- row_vals[sorted_idx]

    order_str <- sorted_nuts[1]
    if (length(sorted_nuts) > 1) {
      for (k in 2:length(sorted_nuts)) {
        sep <- if (abs(sorted_vals[k] - sorted_vals[k-1]) < 1e-5) " = " else " > "
        order_str <- paste0(order_str, sep, sorted_nuts[k])
      }
    }
    return(order_str)
  })
  compiled_df$MDRIS_Ratio_Order <- req_orders

  # Format numeric fields in compiled_df
  compiled_df_char <- as.data.frame(lapply(compiled_df, function(col) {
    if (is.numeric(col)) sprintf("%.2f", col) else as.character(col)
  }), stringsAsFactors = FALSE)
  colnames(compiled_df_char) <- colnames(compiled_df)

  # Create Row 1 (Method)
  cols <- colnames(compiled_df)
  I_cols <- grep("^I_", cols)
  J_cols <- grep("^J_log_", cols)
  K_cols <- grep("^K_", cols)

  header1 <- rep("", length(cols))
  header1[1] <- "Method"
  if (length(I_cols) > 0) header1[I_cols] <- "MDRIS Ratio"
  if (length(J_cols) > 0) header1[J_cols] <- "Log Concentration Index"
  if (length(K_cols) > 0) header1[K_cols] <- "Direct Concentration Index"

  # Create Row 2 (Ratios/Indicators)
  header2 <- cols

  # Convert to data frames to prepend
  header1_df <- as.data.frame(t(header1), stringsAsFactors = FALSE)
  colnames(header1_df) <- cols

  header2_df <- as.data.frame(t(header2), stringsAsFactors = FALSE)
  colnames(header2_df) <- cols

  # Combine to form the final table (Row 1: Method, Row 2: Ratios, Row 3: Samples)
  compiled_table_final <- rbind(header1_df, header2_df, compiled_df_char)

  return(list(
    weighted_deviations = wdr_df,
    compiled_table = compiled_table_final,
    raw_indices_I = indices_matrix,
    raw_indices_J = j_matrix,
    raw_indices_K = k_matrix
  ))
}



# ==============================================================================
# NEW MODULE: PLANT ANALYSIS WITH STANDARDIZED SCORES (PASS) DIAGNOSTICS
# ==============================================================================

#' Perform Plant Analysis with Standardized Scores (PASS) with Custom Norms
#'
#' This function implements the hybrid PASS diagnostic engine of Baldock and Schulte (1996)
#' using data frame inputs for INI and DNI norms. It is highly robust and avoids R index-stripping issues.
#'
#' @param data A data frame containing tissue concentration columns for the crop.
#' @param ini_norms A data frame containing columns: `nutrient`, `CL`, `SD`, and `group`.
#' @param dni_norms A data frame containing columns: `ratio`, `mean`, and `SD`.
#' @param sample_col Character string representing the column name of the sample identifier.
#' @return A list containing:
#'   1. `compiled_table`: A master data frame containing calculated indices, PASS YI, and recommendations.
#'   2. `raw_ini`: Numeric matrix of calculated raw INI values.
#'   3. `raw_dni`: Numeric matrix of calculated raw DNI values.
#' @examples
#' \donttest{
#' set.seed(123)
#'
#' # Tissue samples to be diagnosed
#' my_corn_samples <- data.frame(
#'   Sample = paste0("Corn_", 1:10),
#'   N  = runif(10, 2.5, 3.5),
#'   P  = runif(10, 0.20, 0.40),
#'   K  = runif(10, 1.80, 2.60),
#'   S  = runif(10, 0.15, 0.30),
#'   Ca = runif(10, 0.40, 0.90),
#'   Mg = runif(10, 0.20, 0.50)
#' )
#'
#' # Independent Nutrient Index (INI) norms: critical level (CL), SD, and
#' # whether the nutrient participates in the "common" DNI calculations
#' ini_norms <- data.frame(
#'   nutrient = c("N", "P", "K", "S", "Ca", "Mg"),
#'   CL       = c(3.00, 0.30, 2.20, 0.22, 0.65, 0.35),
#'   SD       = c(0.25, 0.04, 0.20, 0.03, 0.10, 0.05),
#'   group    = c("common", "common", "common",
#'                "rare", "rare", "rare"),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Dependent Nutrient Index (DNI) norms: mean and SD for each nutrient ratio
#' dni_norms <- data.frame(
#'   ratio = c("N/P", "P/N", "N/K", "K/N", "P/K", "K/P"),
#'   mean  = c(10.00, 0.10, 1.36, 0.74, 0.136, 7.33),
#'   SD    = c(1.20, 0.02, 0.15, 0.08, 0.020, 0.90),
#'   stringsAsFactors = FALSE
#' )
#'
#' # Run the PASS analysis
#' pass_report <- calculate_pass_analysis(
#'   data = my_corn_samples,
#'   ini_norms = ini_norms,
#'   dni_norms = dni_norms,
#'   sample_col = "Sample"
#' )
#' print(pass_report$compiled_table)
#' }
#' @export
calculate_pass_analysis <- function(data, ini_norms, dni_norms, sample_col = "Sample") {

  # ============================================================================
  # STEP 1: VALIDATION AND INITIALIZATION
  # ============================================================================
  if (!sample_col %in% colnames(data)) {
    stop(paste("Sample column '", sample_col, "' not found in the dataset."))
  }

  required_ini_cols <- c("nutrient", "CL", "SD", "group")
  if (!all(required_ini_cols %in% colnames(ini_norms))) {
    stop(paste("ini_norms must contain columns:", paste(required_ini_cols, collapse = ", ")))
  }

  required_dni_cols <- c("ratio", "mean", "SD")
  if (!all(required_dni_cols %in% colnames(dni_norms))) {
    stop(paste("dni_norms must contain columns:", paste(required_dni_cols, collapse = ", ")))
  }

  n_samples <- nrow(data)
  samples   <- data[[sample_col]]
  nutrients <- unique(ini_norms$nutrient)

  # Ensure all required nutrient columns are present in the sample data
  missing_data_cols <- setdiff(nutrients, colnames(data))
  if (length(missing_data_cols) > 0) {
    stop(paste("The following nutrient columns are missing from 'data':",
               paste(missing_data_cols, collapse = ", ")))
  }

  # ============================================================================
  # STEP 2: CALCULATE INDEPENDENT NUTRIENT INDICES (INI)
  # ============================================================================
  ini_matrix <- matrix(0, nrow = n_samples, ncol = length(nutrients))
  colnames(ini_matrix) <- paste0(nutrients, "_INI")

  for (j in seq_along(nutrients)) {
    nut <- nutrients[j]
    norm_row <- ini_norms[ini_norms$nutrient == nut, ]

    cl_val <- norm_row$CL
    sd_val <- norm_row$SD
    sample_concs <- as.numeric(data[[nut]])

    ini_matrix[, j] <- 10 * (sample_concs - cl_val) / sd_val - 10
  }

  # ============================================================================
  # STEP 3: CALCULATE DEPENDENT NUTRIENT INDICES (DNI)
  # ============================================================================
  common_nuts <- ini_norms$nutrient[ini_norms$group == "common"]
  dni_matrix  <- matrix(0, nrow = n_samples, ncol = length(common_nuts))
  colnames(dni_matrix) <- paste0(common_nuts, "_DNI")

  # Index constants to bypass citation filter bracket-stripping
  first_idx  <- 1
  second_idx <- 2

  # Calculate Standardized Scores (SS) for each ratio
  ss_df <- data.frame(Sample = samples, stringsAsFactors = FALSE)
  for (i in 1:nrow(dni_norms)) {
    ratio_expr <- dni_norms$ratio[i]
    norm_mean  <- dni_norms$mean[i]
    norm_sd    <- dni_norms$SD[i]

    parts   <- unlist(strsplit(ratio_expr, "/"))
    nut_num <- parts[first_idx]
    nut_den <- parts[second_idx]

    # Calculate sample ratio
    sample_ratio <- data[[nut_num]] / data[[nut_den]]

    # Standardized ratio score
    ss_df[[ratio_expr]] <- 10 * (sample_ratio - norm_mean) / norm_sd
  }

  # Accumulate ratio scores to compile DNI for each common response nutrient
  for (j in seq_along(common_nuts)) {
    nut <- common_nuts[j]

    # Find all ratio rows in dni_norms involving this nutrient
    matching_idx <- which(sapply(dni_norms$ratio, function(ratio) {
      parts <- unlist(strsplit(ratio, "/"))
      return(nut %in% parts)
    }))

    matching_ratios <- dni_norms$ratio[matching_idx]
    z <- length(matching_ratios)

    if (z > 0) {
      dni_vals <- rep(0, n_samples)
      for (expr in matching_ratios) {
        parts <- unlist(strsplit(expr, "/"))
        sign_mult <- if (parts[first_idx] == nut) 1 else -1
        dni_vals <- dni_vals + (sign_mult * ss_df[[expr]])
      }
      dni_matrix[, j] <- dni_vals / z
    }
  }

  # ============================================================================
  # STEP 4: IDENTIFY DEFICIENCIES & ACTIONABLE RECOMMENDATIONS
  # ============================================================================
  rec_matrix <- matrix("", nrow = n_samples, ncol = length(nutrients))
  colnames(rec_matrix) <- paste0(nutrients, "_Status")

  for (i in 1:n_samples) {
    for (j in seq_along(nutrients)) {
      nut <- nutrients[j]
      norm_row <- ini_norms[ini_norms$nutrient == nut, ]

      ini_val <- ini_matrix[i, j]
      grp_type <- norm_row$group

      dni_col <- paste0(nut, "_DNI")
      dni_val <- if (dni_col %in% colnames(dni_matrix)) dni_matrix[i, dni_col] else NA

      if (grp_type == "common") {
        if (ini_val < -10) {
          rec_matrix[i, j] <- "Probable Deficiency"
        } else if (ini_val < 0 && !is.na(dni_val) && (ini_val + dni_val) < -10) {
          rec_matrix[i, j] <- "Slightly Possible Deficiency"
        } else {
          rec_matrix[i, j] <- "Sufficient"
        }
      } else {
        if (ini_val < -10) {
          rec_matrix[i, j] <- "Slightly Possible Deficiency"
        } else {
          rec_matrix[i, j] <- "Sufficient"
        }
      }
    }
  }

  # ============================================================================
  # STEP 5: CALCULATE PASS YIELD INDEX (PASS YI)
  # ============================================================================
  pass_yi <- rep(0, n_samples)
  for (i in 1:n_samples) {
    y_vals <- rep(0, length(nutrients))

    for (j in seq_along(nutrients)) {
      nut <- nutrients[j]
      norm_row <- ini_norms[ini_norms$nutrient == nut, ]

      ini_val <- ini_matrix[i, j]
      grp_type <- norm_row$group

      y_val <- if (ini_val < -10) (ini_val + 10)^2 else 0
      y_vals[j] <- if (grp_type == "common") 2 * y_val else y_val
    }

    pass_yi[i] <- sum(y_vals)
  }

  # ============================================================================
  # STEP 6: COMPILE MASTER TABLE
  # ============================================================================
  compiled_df <- data.frame(Sample = samples, stringsAsFactors = FALSE)
  colnames(compiled_df) <- sample_col

  # Format and append INI columns
  for (nut in nutrients) {
    compiled_df[[paste0("INI_", nut)]] <- round(ini_matrix[, paste0(nut, "_INI")], 2)
  }

  # Format and append DNI columns
  for (nut in common_nuts) {
    compiled_df[[paste0("DNI_", nut)]] <- round(dni_matrix[, paste0(nut, "_DNI")], 2)
  }

  # Append overall Yield Index
  compiled_df$PASS_YI <- round(pass_yi, 2)

  # Append qualitative statuses
  for (nut in nutrients) {
    compiled_df[[paste0("Status_", nut)]] <- rec_matrix[, paste0(nut, "_Status")]
  }

  return(list(
    compiled_table = compiled_df,
    raw_ini = ini_matrix,
    raw_dni = dni_matrix
  ))
}

