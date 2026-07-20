#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### latentFactoR S3Methods ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Updated 30.09.2022

# print() Methods ----

# Print `lf_estimate`
# Updated 30.09.2022
#' @export
print.lf_estimate <- function(x, ...)
{

  # Print dimensions
  print(x$dimensions)

}

# Print `lf_simulate`
# Handles output from both `simulate_factors` (flat) and
# `simulate_hierarchical_factors` (hierarchical), which share this class
# Updated 12.07.2026
#' @export
print.lf_simulate <- function(x, ...)
{

  # Print compact overview
  lf_simulate_describe(x, detailed = FALSE)

  # Return input invisibly
  invisible(x)

}

# summary() Methods ----

# Summary `lf_estimate`
# Updated 30.09.2022
#' @export
summary.lf_estimate <- function(object, ...)
{

  # Print dimensions
  print(object$dimensions)

}

# Summary `lf_simulate`
# Updated 12.07.2026
#' @export
summary.lf_simulate <- function(object, ...)
{

  # Print detailed overview
  lf_simulate_describe(object, detailed = TRUE)

  # Return input invisibly
  invisible(object)

}

# Shared helpers ----

#' @noRd
# Formats the range and mean of a numeric vector for `lf_simulate` print/summary
# Updated 12.07.2026
lf_simulate_format_range <- function(values)
{
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
# Shared print/summary body for `lf_simulate` objects (flat or hierarchical)
# Updated 12.07.2026
lf_simulate_describe <- function(x, detailed = FALSE)
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
  cat(strrep("-", 35), "\n", sep = "")

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
    if(length(lower$cross) != 0){
      cat("  Cross-loadings: ", lf_simulate_format_range(lower$cross), "\n", sep = "")
    }

    cat("\nLower-order factors onto higher-order factors:\n")
    cat("  Loadings:       ", lf_simulate_format_range(higher$dominant), "\n", sep = "")
    if(length(higher$cross) != 0){
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

    disturbance_offdiag <- higher_parameters$disturbances[lower.tri(higher_parameters$disturbances)]
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
      print(round(higher_parameters$disturbances, 3))

    }

  }else{

    cat(sprintf("%-24s%d\n", "Factors:", parameters$factors))

    loadings <- lf_simulate_summarize_loadings(parameters$loadings, parameters$variables)

    cat("\nLoadings:       ", lf_simulate_format_range(loadings$dominant), "\n", sep = "")
    if(length(loadings$cross) != 0){
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

# predict() Methods ----
# {xgboost}: predictLearner
# Updated 30.09.2022
#' @export
predictLearner.classif.xgboost.earlystop = function(.learner, .model, .newdata, ...) {
  td = .model$task.desc
  m = .model$learner.model
  cls = td$class.levels
  nc = length(cls)
  obj = .learner$par.vals$objective

  if (is.null(obj))
    .learner$par.vals$objective = ifelse(nc == 2L, "binary:logistic", "multi:softprob")

  p = predict(m, newdata = data.matrix(BBmisc::convertDataFrameCols(.newdata, ints.as.num = TRUE)), ...)

  if (nc == 2L) { #binaryclass
    if (.learner$par.vals$objective == "multi:softprob") {
      y = matrix(p, nrow = length(p) / nc, ncol = nc, byrow = TRUE)
      colnames(y) = cls
    } else {
      y = matrix(0, ncol = 2, nrow = nrow(.newdata))
      colnames(y) = cls
      y[, 1L] = 1 - p
      y[, 2L] = p
    }
    if (.learner$predict.type == "prob") {
      return(y)
    } else {
      p = colnames(y)[max.col(y)]
      names(p) = NULL
      p = factor(p, levels = colnames(y))
      return(p)
    }
  } else { #multiclass
    if (.learner$par.vals$objective  == "multi:softmax") {
      return(factor(p, levels = cls)) #special handling for multi:softmax which directly predicts class levels
    } else {
      p = matrix(p, nrow = length(p) / nc, ncol = nc, byrow = TRUE)
      colnames(p) = cls
      if (.learner$predict.type == "prob") {
        return(p)
      } else {
        ind = max.col(p)
        cns = colnames(p)
        return(factor(cns[ind], levels = cns))
      }
    }
  }
}