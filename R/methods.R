#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### latentFactoR S3Methods ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Conventions (shared with {psychNetwoRk}): a title, a rule of dashes,
# left-justified labels (`%-24s`), values with three decimals, and
# `print` and `summary` both return the input invisibly. `summary` adds
# details (matrices, fit indices) to the `print` overview

# Updated 20.11.2026

# lf_simulate ----

# Print `lf_simulate`
# Handles output from `simulate_factors`, `simulate_hierarchical_factors`,
# and the functions that modify them (`add_cross_loadings`,
# `add_local_dependence`, `add_method_factors`, `add_wording_effects`,
# `add_population_error`, and `data_to_zipfs`), which append their own class
# Updated 20.11.2026
#' @export
print.lf_simulate <- function(x, ...)
{

  # Print compact overview
  lf_simulate_describe(x, detailed = FALSE)

  # Return input invisibly
  invisible(x)

}

# Summary `lf_simulate`
# Updated 20.11.2026
#' @export
summary.lf_simulate <- function(object, ...)
{

  # Print detailed overview
  lf_simulate_describe(object, detailed = TRUE)

  # Return input invisibly
  invisible(object)

}

# lf_estimate ----

# Print `lf_estimate` (`estimate_dimensions`)
# Updated 20.11.2026
#' @export
print.lf_estimate <- function(x, ...)
{

  # Header
  lf_header("Estimated Dimensions")

  # Dimensions by method
  lf_print_dimensions(x$dimensions)

  # Return input invisibly
  invisible(x)

}

# Summary `lf_estimate`
# Updated 20.11.2026
#' @export
summary.lf_estimate <- function(object, ...)
{

  # Overview
  print.lf_estimate(object)

  # Agreement across methods
  dimensions <- object$dimensions[!is.na(object$dimensions)]
  if(length(dimensions) > 0){
    cat("\n")
    cat(sprintf("%-24s%s\n", "Range:", paste(range(dimensions), collapse = " to ")))
    cat(sprintf("%-24s%s\n", "Most frequent:", paste(lf_mode(dimensions), collapse = ", ")))
  }

  # Non-estimated methods
  missing_methods <- names(object$dimensions)[is.na(object$dimensions)]
  if(length(missing_methods) > 0){
    cat(sprintf("%-24s%s\n", "Not estimated:", paste(missing_methods, collapse = ", ")))
  }

  # NEST convergence
  if(!is.null(object$NEST$converged)){
    cat(sprintf("%-24s%s\n", "NEST converged:", object$NEST$converged))
  }

  # Return input invisibly
  invisible(object)

}

# Estimation functions ----

# Print `lf_ekc` (`EKC`)
# Updated 20.11.2026
#' @export
print.lf_ekc <- function(x, ...)
{
  lf_header("Empirical Kaiser Criterion")
  cat(sprintf("%-24s%d\n", "Dimensions:", as.integer(x$dimensions)))
  invisible(x)
}

# Summary `lf_ekc`
# Updated 20.11.2026
#' @export
summary.lf_ekc <- function(object, ...)
{
  print.lf_ekc(object)
  shown <- seq_len(min(length(object$eigenvalues), max(object$dimensions + 2, 5)))
  cat("\nEigenvalues and reference values:\n")
  print(
    round(
      data.frame(
        eigenvalue = object$eigenvalues[shown],
        reference = object$reference[shown],
        row.names = shown
      ), 3
    )
  )
  invisible(object)
}

# Print `lf_nest` (`NEST`)
# Updated 20.11.2026
#' @export
print.lf_nest <- function(x, ...)
{
  lf_header("Next Eigenvalue Sufficiency Test")
  cat(sprintf("%-24s%d\n", "Dimensions:", as.integer(x$dimensions)))
  cat(sprintf("%-24s%s\n", "Converged:", x$converged))
  invisible(x)
}

# Summary `lf_nest`
# Updated 20.11.2026
#' @export
summary.lf_nest <- function(object, ...)
{
  print.lf_nest(object)
  if(length(object$loadings) > 0 && ncol(object$loadings) > 0){
    cat("\nLoadings:\n")
    print(round(object$loadings, 3))
  }
  invisible(object)
}

# Print `lf_ff` (`factor_forest`)
# Updated 20.11.2026
#' @export
print.lf_ff <- function(x, ...)
{
  lf_header("Factor Forest")
  cat(sprintf("%-24s%d\n", "Dimensions:", as.integer(x$dimensions)))
  invisible(x)
}

# Summary `lf_ff`
# Updated 20.11.2026
#' @export
summary.lf_ff <- function(object, ...)
{
  print.lf_ff(object)
  cat("\nProbabilities:\n")
  print(round(unlist(object$probabilities), 3))
  invisible(object)
}

# Print `lf_esem` (`ESEM`)
# Updated 20.11.2026
#' @export
print.lf_esem <- function(x, ...)
{

  # Header
  lf_header("Exploratory Structural Equation Model")

  # Check for convergence
  if(!is(x$model, "lavaan")){
    cat("Model could not converge.\n")
    return(invisible(x))
  }

  # Main fit measures
  fit <- x$fit
  names(fit) <- gsub(".scaled", "", names(fit), fixed = TRUE)
  main <- intersect(c("chisq", "df", "cfi", "tli", "rmsea", "srmr"), names(fit))
  labels <- c(
    chisq = "Chi-square:", df = "Degrees of freedom:", cfi = "CFI:",
    tli = "TLI:", rmsea = "RMSEA:", srmr = "SRMR:"
  )
  for(measure in main){
    value <- fit[[measure]]
    cat(
      sprintf(
        "%-24s%s\n", labels[[measure]],
        swiftelse(measure == "df", format(round(value)), format(round(value, 3), nsmall = 3))
      )
    )
  }

  # Return input invisibly
  invisible(x)

}

# Summary `lf_esem`
# Updated 20.11.2026
#' @export
summary.lf_esem <- function(object, ...)
{

  # Overview
  print.lf_esem(object)

  # All fit measures
  if(is(object$model, "lavaan")){
    cat("\nFit measures:\n")
    print(unclass(object$fit))
    cat("\nUse `lavaan::summary(<object>$model)` for parameter estimates\n")
  }

  # Return input invisibly
  invisible(object)

}

# Print `lf_zipfs_parameters` (`obtain_zipfs_parameters`)
# Updated 20.11.2026
#' @export
print.lf_zipfs_parameters <- function(x, ...)
{
  lf_header("Zipf's Distribution Parameters")
  cat(sprintf("%-24s%.3f\n", "Alpha:", x[["alpha"]]))
  cat(sprintf("%-24s%.3f\n", "Beta:", x[["beta"]]))
  cat(sprintf("%-24s%.5f\n", "Sum of squared error:", x[["zipfs_sse"]]))
  invisible(x)
}

# Summary `lf_zipfs_parameters`
# Updated 20.11.2026
#' @export
summary.lf_zipfs_parameters <- function(object, ...)
{
  print.lf_zipfs_parameters(object)
  invisible(object)
}

# Shared helpers ----

#' @noRd
# Print a title with a rule of dashes
# Updated 20.11.2026
lf_header <- function(title)
{
  cat(title, "\n", sep = "")
  cat(strrep("-", 40), "\n", sep = "")
}

#' @noRd
# Print a named vector of dimensions (one per line)
# Updated 20.11.2026
lf_print_dimensions <- function(dimensions)
{
  for(method in names(dimensions)){
    cat(
      sprintf(
        "%-24s%s\n", paste0(gsub("_", " ", method), ":"),
        swiftelse(is.na(dimensions[[method]]), "NA", as.character(dimensions[[method]]))
      )
    )
  }
}

#' @noRd
# Most frequent value(s)
# Updated 20.11.2026
lf_mode <- function(values)
{
  counts <- table(values)
  return(as.numeric(names(counts)[counts == max(counts)]))
}

#' @noRd
# Formats the range and mean of a numeric vector for print/summary
# Updated 12.07.2026
lf_simulate_format_range <- function(values)
{

  # Single value
  if(length(unique(values)) == 1){
    return(sprintf("%.3f (all)", values[1]))
  }

  # Range and mean
  sprintf(
    "%.3f to %.3f (mean = %.3f)",
    min(values), max(values), mean(values)
  )

}

#' @noRd
# Splits a loading matrix into dominant (within-block) and cross (outside-block)
# loadings, where `block_sizes` gives the number of rows belonging to each column
# (e.g. variables per factor, or lower-order factors per higher-order factor)
# Updated 12.07.2026
lf_simulate_summarize_loadings <- function(loading_matrix, block_sizes)
{

  sequence <- factor_variable_sequence(block_sizes)
  dominant <- numeric(0)
  cross <- numeric(0)

  for(i in seq_along(block_sizes)){

    block <- sequence$start[i]:sequence$end[i]
    column <- loading_matrix[, i]

    dominant <- c(dominant, column[block])
    cross <- c(cross, column[-block])

  }

  list(dominant = dominant, cross = cross)

}

#' @noRd
# Shared print/summary body for `lf_simulate` objects
# Describes the original simulated data (the root of the `original_results`
# chain) and then each modification in the order it was added
# Updated 20.11.2026
lf_simulate_describe <- function(x, detailed = FALSE)
{

  # Obtain chain of objects (root first)
  chain <- list()
  object <- x
  while(!is.null(object)){
    chain <- c(list(object), chain)
    object <- object$original_results
  }

  # Describe the original simulated data
  lf_simulate_base_describe(chain[[1]], detailed = detailed)

  # Describe each modification
  for(object in chain[-1]){

    # Class added by the modification
    added <- class(object)[length(class(object))]

    # Describe
    switch(
      added,
      "lf_cl" = lf_describe_cl(object, detailed),
      "lf_ld" = lf_describe_ld(object, detailed),
      "lf_me" = lf_describe_me(object, detailed),
      "lf_we" = lf_describe_we(object, detailed),
      "lf_pe" = lf_describe_pe(object, detailed),
      "lf_zipfs" = lf_describe_zipfs(object, detailed)
    )

  }

}

#' @noRd
# Description of the simulated factor data (flat or hierarchical)
# Updated 20.11.2026
lf_simulate_base_describe <- function(x, detailed = FALSE)
{

  parameters <- x$parameters
  hierarchical <- !is.null(parameters$lower$factors)
  dimensions <- dim(x$data)

  # Header
  cat(
    if(isTRUE(hierarchical)){
      "Simulated Hierarchical Factor Data\n"
    }else{
      "Simulated Factor Data\n"
    }
  )
  cat(strrep("-", 40), "\n", sep = "")

  cat(sprintf("%-24s%d\n", "Sample size:", dimensions[1]))
  cat(sprintf("%-24s%d\n", "Variables:", dimensions[2]))

  if(isTRUE(hierarchical)){

    lower_parameters <- parameters$lower
    higher_parameters <- parameters$higher

    cat(sprintf("%-24s%d\n", "Lower-order factors:", lower_parameters$factors))
    cat(sprintf("%-24s%d\n", "Higher-order factors:", higher_parameters$factors))
    cat(
      sprintf("%-24s", "Factors per higher:"),
      paste(higher_parameters$variables, collapse = ", "), "\n", sep = ""
    )

    lower <- lf_simulate_summarize_loadings(lower_parameters$loadings, lower_parameters$variables)
    higher <- lf_simulate_summarize_loadings(higher_parameters$loadings, higher_parameters$variables)

    cat("\nVariables onto lower-order factors:\n")
    cat("  Loadings:       ", lf_simulate_format_range(lower$dominant), "\n", sep = "")
    if(any(lower$cross != 0)){
      cat("  Cross-loadings: ", lf_simulate_format_range(lower$cross), "\n", sep = "")
    }

    cat("\nLower-order factors onto higher-order factors:\n")
    cat("  Loadings:       ", lf_simulate_format_range(higher$dominant), "\n", sep = "")
    if(any(higher$cross != 0)){
      cat("  Cross-loadings: ", lf_simulate_format_range(higher$cross), "\n", sep = "")
    }

    lower_correlations_offdiag <- lower_parameters$correlations[lower.tri(lower_parameters$correlations)]
    if(length(lower_correlations_offdiag) != 0){
      cat(
        "\nLower-order factor correlations: ",
        lf_simulate_format_range(lower_correlations_offdiag), "\n", sep = ""
      )
    }

    higher_correlations_offdiag <- higher_parameters$correlations[lower.tri(higher_parameters$correlations)]
    if(length(higher_correlations_offdiag) != 0){
      cat(
        "Higher-order factor correlations: ",
        lf_simulate_format_range(higher_correlations_offdiag), "\n", sep = ""
      )
    }

    disturbance_correlations <- disturbance_correlation_matrix(higher_parameters$disturbances)
    disturbance_offdiag <- disturbance_correlations[lower.tri(disturbance_correlations)]
    if(length(disturbance_offdiag) != 0){
      cat(
        "\nLower-order disturbance correlations: ",
        lf_simulate_format_range(disturbance_offdiag), "\n", sep = ""
      )
    }

    if(isTRUE(detailed)){

      cat("\nImplied lower-order factor correlation matrix:\n")
      print(round(lower_parameters$correlations, 3))

      cat("\nLower-order disturbance (residual) correlation matrix:\n")
      print(round(disturbance_correlations, 3))

    }

  }else{

    cat(sprintf("%-24s%d\n", "Factors:", parameters$factors))

    loadings <- lf_simulate_summarize_loadings(parameters$loadings, parameters$variables)

    cat("\nLoadings:       ", lf_simulate_format_range(loadings$dominant), "\n", sep = "")
    if(any(loadings$cross != 0)){
      cat("Cross-loadings: ", lf_simulate_format_range(loadings$cross), "\n", sep = "")
    }

    correlations_offdiag <- parameters$factor_correlations[lower.tri(parameters$factor_correlations)]
    if(length(correlations_offdiag) != 0){
      cat("Factor correlations: ", lf_simulate_format_range(correlations_offdiag), "\n", sep = "")
    }

    if(isTRUE(detailed) && parameters$factors > 1){
      cat("\nFactor correlation matrix:\n")
      print(round(parameters$factor_correlations, 3))
    }

  }

  # Shared: categorical variables and skew
  categories <- parameters$categories
  continuous <- sum(is.infinite(categories))
  categorical <- length(categories) - continuous

  cat("\n")
  if(categorical == 0){
    cat("Categories:  all variables continuous\n")
  }else if(continuous == 0){
    cat(sprintf(
      "Categories:  all variables categorical, %s\n",
      lf_simulate_format_range(categories)
    ))
  }else{
    cat(sprintf(
      "Categories:  %d continuous, %d categorical (categories: %s)\n",
      continuous, categorical, lf_simulate_format_range(categories[!is.infinite(categories)])
    ))
  }

  if(any(parameters$skew != 0)){
    cat(sprintf("Skew:        %s\n", lf_simulate_format_range(parameters$skew)))
  }else{
    cat("Skew:        none\n")
  }

  if(isTRUE(detailed) && categorical != 0){
    cat("\nCategory counts:\n")
    print(table(categories[!is.infinite(categories)], dnn = NULL))
  }
}

#' @noRd
# Disturbance correlations from the disturbance covariance matrix
# (disturbances with zero variance have zero correlations)
# Updated 20.11.2026
disturbance_correlation_matrix <- function(disturbances)
{
  sds <- sqrt(pmax(diag(disturbances), 0))
  scaling <- swiftelse(sds > 0, 1 / sds, 0)
  correlations <- disturbances * tcrossprod(scaling)
  diag(correlations) <- 1
  return(correlations)
}

#' @noRd
# Section header for a modification
# Updated 20.11.2026
lf_section <- function(title)
{
  cat("\n", title, ":\n", sep = "")
}

#' @noRd
# Cross-loadings (`add_cross_loadings`)
# Updated 20.11.2026
lf_describe_cl <- function(x, detailed = FALSE)
{

  # Obtain loadings
  parameters <- x$parameters
  loadings <- lf_simulate_summarize_loadings(parameters$loadings, parameters$variables)
  cross <- loadings$cross[loadings$cross != 0]

  # Describe
  lf_section("Cross-loadings added")
  cat(sprintf("  %-24s%d\n", "Cross-loadings:", length(cross)))
  if(length(cross) > 0){
    cat(sprintf("  %-24s%s\n", "Magnitude (absolute):", lf_simulate_format_range(abs(cross))))
  }

  # Detailed loadings
  if(isTRUE(detailed)){
    cat("\n  Loadings:\n")
    print(round(parameters$loadings, 3))
  }

}

#' @noRd
# Local dependence (`add_local_dependence`)
# Updated 20.11.2026
lf_describe_ld <- function(x, detailed = FALSE)
{

  # Obtain residuals
  residuals <- x$correlated_residuals

  # Describe
  lf_section("Local dependence added")
  cat(sprintf("  %-24s%d\n", "Correlated residuals:", nrow(residuals)))
  if(nrow(residuals) > 0){
    cat(
      sprintf(
        "  %-24s%s\n", "Residuals:",
        lf_simulate_format_range(residuals$residual_correlation)
      )
    )
  }

  # Detailed pairs
  if(isTRUE(detailed) && nrow(residuals) > 0){
    cat("\n  Locally dependent pairs:\n")
    print(
      data.frame(
        variable_1 = residuals$item_rows, variable_2 = residuals$item_columns,
        residual = round(residuals$residual_correlation, 3)
      ), row.names = FALSE
    )
  }

}

#' @noRd
# Method factors (`add_method_factors`)
# Updated 20.11.2026
lf_describe_me <- function(x, detailed = FALSE)
{

  # Obtain method loadings
  parameters <- x$parameters
  methods <- seq.int(parameters$factors + 1, ncol(parameters$loadings))
  method_loadings <- parameters$loadings[, methods, drop = FALSE]
  active <- colSums(method_loadings != 0) > 0
  values <- method_loadings[method_loadings != 0]

  # Describe
  lf_section("Method factors added")
  cat(sprintf("  %-24s%d\n", "Method factors:", sum(active)))
  if(length(values) > 0){
    cat(sprintf("  %-24s%s\n", "Method loadings:", lf_simulate_format_range(values)))
  }
  negative <- sum(dominant_signs(parameters$loadings, parameters$variables) < 0)
  cat(sprintf("  %-24s%d\n", "Negative loadings:", negative))

  # Detailed correlations
  if(isTRUE(detailed) && sum(active) > 1){
    cat("\n  Method factor correlations:\n")
    print(round(parameters$factor_correlations[methods[active], methods[active]], 3))
  }

}

#' @noRd
# Wording effects (`add_wording_effects`)
# Updated 20.11.2026
lf_describe_we <- function(x, detailed = FALSE)
{

  # Describe
  lf_section("Wording effects added")
  effects <- table(x$replaced_sample_effects)
  cat(sprintf("  %-24s%s\n", "Method(s):", paste(names(effects), collapse = ", ")))
  cat(
    sprintf(
      "  %-24s%d (%.1f%%)\n", "Biased cases:", as.integer(x$biased_sample_size),
      100 * x$biased_sample_size / nrow(x$data)
    )
  )
  adjusted <- x$adjusted_results$parameters
  negative <- sum(dominant_signs(adjusted$loadings, adjusted$variables) < 0)
  cat(sprintf("  %-24s%d\n", "Negative loadings:", negative))

  # Detailed counts
  if(isTRUE(detailed) && length(effects) > 1){
    cat("\n  Biased cases by method:\n")
    print(effects)
  }

}

#' @noRd
# Population error (`add_population_error`)
# Updated 20.11.2026
lf_describe_pe <- function(x, detailed = FALSE)
{

  # Obtain population error
  error <- x$population_error

  # Describe
  lf_section("Population error added")
  cat(sprintf("  %-24s%s\n", "Method:", toupper(error$method)))
  cat(sprintf("  %-24s%.3f\n", paste0("Target ", error$metric, ":"), error$fit_value))
  cat(sprintf("  %-24s%.3f\n", paste0("Achieved ", error$metric, ":"), error$achieved))
  cat(sprintf("  %-24s%.5f\n", "Delta:", error$delta))
  cat(sprintf("  %-24s%.3f\n", "Maximum error:", error$max_error))

  # Detailed fit
  if(isTRUE(detailed)){
    error_values <- error$error[upper.tri(error$error)]
    cat(sprintf("  %-24s%s\n", "Error (absolute):", lf_simulate_format_range(abs(error_values))))
    cat(
      sprintf(
        "  %-24s%s\n", "Fit indices:",
        paste(names(error$fit), sprintf("%.3f", error$fit), collapse = ", ")
      )
    )
    cat(sprintf("  %-24s%.2e\n", "Orthogonality (max):", error$R_difference))
    cat(sprintf("  %-24s%d\n", "Optimization attempts:", as.integer(error$iterations)))
    cat("\n  Loadings of the true model fit to the population with error:\n")
    print(round(error$loadings, 3))
  }

}

#' @noRd
# Zipf's distribution (`data_to_zipfs`)
# Updated 20.11.2026
lf_describe_zipfs <- function(x, detailed = FALSE)
{

  # Describe
  lf_section("Data transformed to Zipf's distribution")
  RMSE <- unlist(x$RMSE)
  for(name in names(RMSE)){
    cat(sprintf("  %-24s%.3f\n", paste0("RMSE (", name, "):"), RMSE[[name]]))
  }

  # Detailed correlations
  if(isTRUE(detailed)){
    cat("\n  Spearman correlations:\n")
    print(round(x$spearman_correlation, 3))
  }

}
