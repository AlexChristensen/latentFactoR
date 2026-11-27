#' Obtain Fit Indices
#'
#' Computes fit indices of a model-implied correlation matrix relative to a
#' correlation matrix, either for a sample (finite \code{sample_size}) or for
#' the population (\code{sample_size = Inf})
#'
#' @param R_model Numeric matrix or data frame.
#' Model-implied correlation matrix.
#' Must be positive definite
#'
#' @param R Numeric matrix or data frame.
#' Sample or population correlation matrix with the same variables (in the same
#' order) as \code{R_model}.
#' Must be positive definite
#'
#' @param sample_size Numeric (length = 1).
#' Number of cases used to compute \code{R}.
#' Use \code{Inf} for population fit indices
#'
#' @param df Numeric (length = 1).
#' Degrees of freedom of the model. For a model of a correlation matrix with
#' \eqn{p} variables, this is \eqn{p (p - 1) / 2} minus the number of free
#' parameters (excluding uniquenesses, which are determined by the unit
#' diagonal). For example, a model with simple structure has \eqn{p} loadings
#' and \eqn{m (m - 1) / 2} factor correlations for \eqn{m} factors
#'
#' @details
#' The maximum likelihood discrepancy of the model is
#' \eqn{F = \log|R_{model}| + tr(R R_{model}^{-1}) - \log|R| - p} and that of
#' the independence model is \eqn{F_0 = -\log|R|} with
#' \eqn{df_0 = p (p - 1) / 2}.
#'
#' For a finite sample size \eqn{N}, the chi-square statistics are
#' \eqn{\chi^2 = N F} and \eqn{\chi^2_0 = N F_0} (as in
#' \code{\link[lavaan]{lavaan}} with maximum likelihood estimation):
#'
#' \itemize{
#'
#' \item \eqn{CFI = 1 - \max(\chi^2 - df, 0) /
#' \max(\chi^2_0 - df_0, \chi^2 - df, 0)}
#'
#' \item \eqn{TLI = 1 - ((\chi^2 - df) / df) / ((\chi^2_0 - df_0) / df_0)}
#'
#' \item \eqn{RMSEA = \sqrt{\max(\chi^2 - df, 0) / (df N)}}
#'
#' }
#'
#' CFI and TLI are one when their denominators are zero. For the population
#' (the limit as \eqn{N} goes to infinity):
#'
#' \itemize{
#'
#' \item \eqn{CFI = 1 - F / F_0} (values below zero are set to zero)
#'
#' \item \eqn{TLI = 1 - (F / df) / (F_0 / df_0)}
#'
#' \item \eqn{RMSEA = \sqrt{F / df}}
#'
#' }
#'
#' SRMR does not depend on the sample size. It is the root mean square of the
#' \eqn{p (p + 1) / 2} unique elements (including the diagonal) of
#' \eqn{R - R_{model}}
#'
#' No model is estimated: the indices describe \code{R_model} as it is given.
#' CFI, TLI, and RMSEA have their usual interpretation when \code{R_model} is
#' the maximum likelihood solution of the model for \code{R}
#'
#' @return Returns a named numeric vector containing \code{CFI}, \code{TLI},
#' \code{RMSEA}, and \code{SRMR}
#'
#' @examples
#' # Generate factor data
#' two_factor <- simulate_factors(
#'   factors = 2, # factors = 2
#'   variables = 6, # variables per factor = 6
#'   loadings = 0.55, # loadings between = 0.45 to 0.65
#'   cross_loadings = 0.05, # cross-loadings N(0, 0.05)
#'   correlations = 0.30, # correlation between factors = 0.30
#'   sample_size = 1000 # number of cases = 1000
#' )
#'
#' # Add population error (CFI = 0.95 using Cudeck and Browne)
#' two_factor_cb <- add_population_error(lf_object = two_factor, method = "cb")
#'
#' # Population fit of the true model
#' # (df = 66 correlations - 12 loadings - 1 factor correlation)
#' obtain_fit_indices(
#'   R_model = two_factor_cb$original_correlation,
#'   R = two_factor_cb$population_correlation,
#'   sample_size = Inf, df = 53
#' )
#'
#' # Fit of the true model to the sample data
#' obtain_fit_indices(
#'   R_model = two_factor_cb$original_correlation,
#'   R = cor(two_factor_cb$data),
#'   sample_size = 1000, df = 53
#' )
#'
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#'
# Obtain fit indices
# Updated 27.11.2026
obtain_fit_indices <- function(R_model, R, sample_size, df)
{

  # Ensure matrices
  R_model <- as.matrix(R_model)
  R <- as.matrix(R)

  # Check for appropriate correlation matrices
  # (see below for function; returns Cholesky decomposition)
  model_chol <- correlation_chol(R_model, "R_model")
  R_chol <- correlation_chol(R, "R")

  # Check for matching variables
  if(ncol(R_model) != ncol(R)){
    stop(
      "'R_model' and 'R' must have the same number of variables",
      call. = FALSE
    )
  }
  if(
    !is.null(colnames(R_model)) && !is.null(colnames(R)) &&
    !identical(colnames(R_model), colnames(R))
  ){
    stop(
      "'R_model' and 'R' must have the same variable names in the same order",
      call. = FALSE
    )
  }

  # Check for appropriate types, lengths, and ranges
  object_error(sample_size, c("numeric", "integer"))
  object_error(df, c("numeric", "integer"))
  length_error(sample_size, 1); length_error(df, 1)
  if(is.na(sample_size) || sample_size <= 1){
    stop(
      "Input into 'sample_size' argument must be greater than one (or `Inf`)",
      call. = FALSE
    )
  }
  if(!is.finite(df) || df <= 0){
    stop("Input into 'df' argument must be a positive number", call. = FALSE)
  }

  # Obtain population fit (see `error_optimizations` for function)
  fit <- fit_indices(
    R, reference_model(R_model, model_chol), df, population_chol = R_chol
  )

  # Check for degenerate independence model
  if(!fit$valid){
    stop(
      paste(
        "Fit indices could not be computed.",
        "CFI and TLI are undefined when 'R' is an identity matrix"
      ), call. = FALSE
    )
  }

  # Check for population or sample
  if(is.infinite(sample_size)){

    CFI <- max(fit$CFI, 0)
    TLI <- fit$TLI
    RMSEA <- fit$RMSEA

  }else{

    # Noncentrality (chi-square minus degrees of freedom)
    noncentrality <- sample_size * fit$discrepancy - df
    baseline_noncentrality <- sample_size * fit$baseline - fit$baseline_df

    # Comparative fit index
    denominator <- max(baseline_noncentrality, noncentrality, 0)
    CFI <- swiftelse(
      denominator > 0, 1 - max(noncentrality, 0) / denominator, 1
    )

    # Tucker-Lewis index
    TLI <- swiftelse(
      baseline_noncentrality != 0,
      1 - (noncentrality / df) / (baseline_noncentrality / fit$baseline_df), 1
    )

    # Root mean square error of approximation
    RMSEA <- sqrt(max(noncentrality, 0) / (df * sample_size))

  }

  # Return fit indices
  return(c(CFI = CFI, TLI = TLI, RMSEA = RMSEA, SRMR = fit$SRMR))

}

#' @noRd
# Check for a positive definite correlation matrix
# (returns its Cholesky decomposition)
# Updated 27.11.2026
correlation_chol <- function(R, name)
{

  # Check for numeric, finite, square matrix
  if(!is.numeric(R) || any(!is.finite(R)) || nrow(R) != ncol(R) || ncol(R) < 2){
    stop(
      paste0(
        "Input into '", name, "' argument must be a square numeric matrix ",
        "(at least 2 variables) without missing or infinite values"
      ), call. = FALSE
    )
  }

  # Check for symmetry and unit diagonal (within rounding error)
  if(any(abs(R - t(R)) > 1e-08) || any(abs(diag(R) - 1) > 1e-08)){
    stop(
      paste0(
        "Input into '", name, "' argument must be a correlation matrix ",
        "(symmetric with ones on the diagonal)"
      ), call. = FALSE
    )
  }

  # Check for positive definiteness
  R_chol <- try(chol(R), silent = TRUE)
  if(is(R_chol, "try-error")){
    stop(
      paste0("Input into '", name, "' argument is not positive definite"),
      call. = FALSE
    )
  }

  # Return Cholesky decomposition
  return(R_chol)

}
