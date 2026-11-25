#' Adds Population Error to \code{\link[latentFactoR]{simulate_factors}} Data
#'
#' Adds population (model) error to simulated factor model data.
#' The population correlation matrix is perturbed so that the true factor model
#' (known loading pattern, factor correlations, and residuals) no longer fits
#' perfectly but instead fits at a specified value of a fit index (CFI, TLI,
#' RMSEA, or SRMR). Data are then re-generated from the perturbed population
#' correlation matrix.
#' See examples to get started
#'
#' @details
#' The error is generated with one of four methods (see \code{method}) and
#' optimized so that the \code{metric} of the true factor model, fit to the
#' population correlation matrix with error using maximum likelihood (i.e., the
#' known model), equals \code{fit_value}. The fit indices are computed at the
#' population level using the maximum likelihood discrepancy function between
#' the population correlation matrix with error and the correlation matrix
#' implied by the true model. The error added to any cell of the original
#' population correlation matrix is limited to the range -1 to 1 and the
#' population correlation matrix with error is always positive definite with a
#' unit diagonal.
#'
#' The true model includes the non-zero loadings (cross-loadings are set to
#' zero unless \code{leave_cross_loadings = TRUE}), the non-zero factor
#' correlations, the uniquenesses, and any correlated residuals from
#' \code{\link[latentFactoR]{add_local_dependence}}
#'
#' The fit indices are population values (infinite sample size) computed from
#' the maximum likelihood discrepancy \eqn{F} between the population
#' correlation matrix with error, \eqn{\Sigma}, and the correlation matrix
#' implied by the true model fit to it (\eqn{df} are the degrees of freedom of
#' the true model): \eqn{RMSEA = \sqrt{F / df}} and
#' \eqn{CFI = 1 - F / F_0}, where \eqn{F_0 = -\log|\Sigma|} is the
#' discrepancy of the independence model
#' (\eqn{TLI = 1 - (F / df) / (F_0 / df_0)} with
#' \eqn{df_0 = p (p - 1) / 2} for \eqn{p} variables)
#'
#' Each optimization attempt uses new random draws of the error. Attempts
#' continue (up to \code{max_iterations}) until \code{fit_value} is matched
#' (to within \eqn{10^{-6}}) and the attempt closest to \code{fit_value} is
#' returned. If \code{fit_value} cannot be reached within \code{tolerance}
#' (e.g., because the target requires more error than can be added while
#' remaining positive definite), a warning is given and the closest
#' solution is used
#'
#' @param lf_object Data object from \code{\link[latentFactoR]{simulate_factors}}
#' (or \code{\link[latentFactoR]{add_cross_loadings}},
#' \code{\link[latentFactoR]{add_local_dependence}}, or
#' \code{\link[latentFactoR]{add_method_factors}} applied to it)
#'
#' @param method Character (length = 1).
#' Method used to generate the population error.
#' There is no default; \code{method} must be specified.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"cb"} --- Cudeck and Browne (1992). Error is a random direction
#' that is orthogonal to the maximum likelihood score of every free parameter of
#' the true model, so the true model remains the best-fitting model of the
#' population correlation matrix with error (the true parameters are recovered
#' exactly)
#'
#' \item \code{"tkl"} --- Tucker, Koopman, and Linn (1969). Error is
#' the contribution of many random minor factors whose loadings decay
#' geometrically. The minor factors take a proportion of each variable's
#' uniqueness (the major factors are unchanged) and this proportion is
#' optimized
#'
#' \item \code{"wb"} --- Wu and Browne (2015). The population correlation
#' matrix with error is the inverse of a Wishart draw for the precision matrix
#' (i.e., an inverted Wishart draw for the covariance matrix), where the
#' Wishart is centered on the original population correlation matrix and the
#' result is standardized to a correlation matrix. The dispersion of the
#' distribution is optimized. Draws are not rescaled; draws with an
#' error outside of -1 to 1 are treated as infeasible
#'
#' \item \code{"hb"} --- Hu and Bentler (1998). Underparameterized
#' misspecification: a random subset of the loadings (and, optionally, factor
#' covariances) that are zero in the true model are non-zero in the population
#' (omitted cross-loadings and factor covariances). The omitted parameters have
#' random signs and magnitudes that vary uniformly between one half and one
#' times a common size, and the common size is optimized
#'
#' }
#'
#' @param metric Character (length = 1).
#' Fit index used as the objective criterion. The true model is fit to the
#' population correlation matrix with error (maximum likelihood)
#' and the error is optimized so that \code{metric} equals \code{fit_value}.
#' Defaults to the fit index that is consistent with \code{method}:
#' \code{"CFI"} for \code{"cb"} and \code{"tkl"}, \code{"RMSEA"} for
#' \code{"wb"} (the dispersion is approximately RMSEA\eqn{^2}), and
#' \code{"SRMR"} for \code{"hb"} (error is defined through the size of
#' omitted cross-loadings).
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"CFI"} --- Comparative fit index
#'
#' \item \code{"TLI"} --- Tucker-Lewis index
#'
#' \item \code{"RMSEA"} --- Root mean square error of approximation
#'
#' \item \code{"SRMR"} --- Standardized root mean square residual
#' (computed over the lower triangle of the correlation matrices including
#' the diagonal, following \pkg{lavaan})
#'
#' }
#'
#' @param fit_value Numeric (length = 1).
#' Target value of \code{metric} for the true model fit to the population
#' correlation matrix with error. Must be between 0 and 1 (exclusive).
#' Defaults to the conventional cutoff for the \code{metric}:
#' \code{0.95} for \code{"CFI"} and \code{"TLI"}, \code{0.06} for
#' \code{"RMSEA"}, and \code{0.05} for \code{"SRMR"}
#'
#' @param tolerance Numeric (length = 1).
#' Maximum absolute difference between the target \code{fit_value} and the
#' \code{metric} of the population correlation matrix with error that is
#' accepted without a warning. Ensures that the appropriate amount of
#' population error was added.
#' Defaults to \code{0.01}
#'
#' @param max_iterations Numeric (length = 1).
#' Maximum number of optimization attempts (random draws of the error)
#' to reach the target \code{fit_value} within the specified
#' \code{tolerance}.
#' Defaults to \code{10}
#'
#' @param control List.
#' Method-specific options. Options that do not apply to \code{method} are
#' ignored. Available options:
#'
#' \itemize{
#'
#' \item \code{minor_factors} --- Number of minor factors for
#' \code{method = "tkl"}. Defaults to \code{max(50, 2 * variables)}
#'
#' \item \code{decay} --- Geometric decay (\eqn{\epsilon}) of the minor factor
#' loadings for \code{method = "tkl"}. Must be between 0 and 1.
#' The loading of minor factor \eqn{k} is multiplied by
#' \eqn{(1 - \epsilon)^{k - 1}}, so the variance it contributes is proportional
#' to \eqn{(1 - \epsilon)^{2(k - 1)}}. Smaller values (e.g., \code{0.02}) spread the
#' error over many minor factors (diffuse, unstructured noise) while larger
#' values (e.g., \code{0.30}, the default) concentrate the error in a few
#' minor factors (more structured error). Very large values (\code{> 0.60}) can
#' make the target \code{fit_value} harder to reach within the error bounds
#' (\eqn{[-1, 1]}).
#' Defaults to \code{0.30}
#'
#' \item \code{omitted} --- Number of loadings that are zero in the true
#' model that are omitted cross-loadings for \code{method = "hb"}. If not
#' specified, the number is limited to the number of zero loadings. Defaults to
#' \code{2}
#'
#' \item \code{omitted_covariances} --- Number of factor covariances that are
#' zero in the true model (e.g., between substantive and method factors) that
#' are omitted factor covariances for \code{method = "hb"}. Defaults to
#' \code{0}
#'
#' }
#'
#' @param leave_cross_loadings Boolean.
#' Should cross-loadings be kept in the true model?
#' Defaults to \code{FALSE}.
#' Convergence problems can arise if cross-loadings are kept,
#' so setting them to zero is the default. Only set to \code{TRUE}
#' with careful consideration of the structure. Make sure to perform
#' additional checks that the data are adequate
#'
#' @param seed Numeric (length = 1).
#' Seed for reproducible results. The same \code{seed} gives the same results
#' (and does not change the random number stream of the R session).
#' Defaults to \code{NULL} (random seed)
#'
#' @return Returns a list containing:
#'
#' \item{data}{Simulated data from the population correlation matrix with error
#' (same dimensions, column names, categories, and skew as the input data)}
#'
#' \item{population_correlation}{Population correlation matrix with error added}
#'
#' \item{original_correlation}{Population correlation matrix of the true model
#' \emph{before} population error was added}
#'
#' \item{population_error}{
#' A list containing the parameters used to generate population error:
#'
#' \itemize{
#'
#' \item \code{method} --- Method used to generate the population error
#'
#' \item \code{metric} --- Fit index used as the objective criterion
#'
#' \item \code{fit_value} --- Target value of \code{metric}
#'
#' \item \code{achieved} --- Value of \code{metric} achieved by the true model
#' fit to the population correlation matrix with error
#'
#' \item \code{fit} --- Named vector of all fit indices (\code{CFI},
#' \code{TLI}, \code{RMSEA}, and \code{SRMR}) achieved
#'
#' \item \code{delta} --- Minimum of the (maximum likelihood) discrepancy
#' function corresponding to the achieved fit
#'
#' \item \code{error} --- Matrix of the error that was added to the
#' original population correlation matrix
#'
#' \item \code{max_error} --- Largest absolute error that was added
#'
#' \item \code{R_difference} --- Largest absolute difference between the
#' original population correlation matrix and the true model (maximum
#' likelihood) fit to the population correlation matrix with error.
#' Values near zero indicate that the error is orthogonal to the true model
#' (expected for \code{method = "cb"} only)
#'
#' \item \code{loadings} --- Estimated loadings of the true model fit to the
#' population correlation matrix with error (signs aligned to the original
#' loadings)
#'
#' \item \code{factor_correlations} --- Estimated factor correlations of the
#' true model fit to the population correlation matrix with error
#'
#' \item \code{iterations} --- Number of optimization attempts used
#'
#' \item \code{omitted} --- Only for \code{method = "hb"}. A list of the
#' omitted parameters that were added to the population (their values in the
#' population; true model values are zero): \code{loadings} (matrix of omitted
#' cross-loadings), \code{factor_covariances} (matrix of omitted factor
#' covariances), and \code{size} (common size; the magnitude of each omitted
#' parameter is between one half and one times \code{size})
#'
#' }
#'
#' }
#'
#' \item{original_results}{Original \code{lf_object} input into function}
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
#' \donttest{
#' # Use Wu and Browne to target TLI = 0.90
#' two_factor_wb <- add_population_error(
#'   lf_object = two_factor, method = "wb",
#'   metric = "TLI", fit_value = 0.90
#' )
#'
#' # Target RMSEA = 0.04
#' two_factor_rmsea <- add_population_error(
#'   lf_object = two_factor, method = "cb", metric = "RMSEA",
#'   fit_value = 0.04
#' )
#'
#' # Use Hu and Bentler to target SRMR = 0.05
#' two_factor_srmr <- add_population_error(
#'   lf_object = two_factor, method = "hb", fit_value = 0.05
#' )}
#'
#' @author
#' {\code{\link{latentFactoR}}} authors \cr
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Hudson Golino <hfg9s@virginia.edu>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>,
#' Marcos Jimenez,
#' Francisco J. Abad,
#' Eduardo Garcia-Garzon,
#' Vithor R. Franco
#'
#' @references
#' \strong{Cudeck and Browne (\code{"cb"})} \cr
#' Cudeck, R., & Browne, M.W. (1992).
#' Constructing a covariance matrix that yields a specified minimizer and a specified minimum discrepancy function value.
#' \emph{Psychometrika}, \emph{57}, 357--369.
#'
#' \strong{Tucker, Koopman, and Linn (\code{"tkl"})} \cr
#' Tucker, L.R., Koopman, R.F., & Linn, R.L. (1969).
#' Evaluation of factor analytic research procedures by means of simulated correlation matrices.
#' \emph{Psychometrika}, \emph{34}, 421--459.
#'
#' \strong{Wu and Browne (\code{"wb"})} \cr
#' Wu, H., & Browne, M.W. (2015).
#' Quantifying adventitious error in a covariance structure as a random effect.
#' \emph{Psychometrika}, \emph{80}, 571--600.
#'
#' \strong{Hu and Bentler (\code{"hb"}) and fit indices} \cr
#' Hu, L., & Bentler, P.M. (1998).
#' Fit indices in covariance structure modeling: Sensitivity to underparameterized model misspecification.
#' \emph{Psychological Methods}, \emph{3}, 424--453.
#'
#' @export
#'
# Add population error to simulated data ----
# Updated 25.11.2026
add_population_error <- function(
    lf_object,
    method = c("cb", "tkl", "wb", "hb"),
    metric = c("cfi", "tli", "rmsea", "srmr"),
    fit_value = NULL,
    tolerance = 0.01, max_iterations = 10,
    control = list(),
    leave_cross_loadings = FALSE,
    seed = NULL
)
{

  # Set seed (reproducible results)
  seed_token <- lf_seed_start(seed)
  on.exit(lf_seed_end(seed_token), add = TRUE)


  # Obtain options specified by the user
  specified <- names(control)

  # Check inputs
  inputs <- add_population_error_errors(
    lf_object, method, metric, fit_value,
    tolerance, max_iterations, control, leave_cross_loadings,
    method_missing = missing(method), metric_missing = missing(metric)
  )
  method <- inputs$method
  metric <- inputs$metric
  target <- inputs$fit_value
  control <- inputs$control

  # Obtain true model
  model <- true_model(lf_object, leave_cross_loadings)
  parts <- factor_model_parts(model$loadings, model$phi, model$psi)
  R <- model$R
  nodes <- ncol(R)

  # Check for identification (e.g., cross-loadings on every variable
  # kept with `leave_cross_loadings = TRUE` allow rotation)
  if(!model_identified(parts)){
    stop(
      paste(
        "The true model is not identified (its parameters are not unique),",
        "so its degrees of freedom and fit are not defined.",
        "\n\nUse `leave_cross_loadings = FALSE` or a factor structure with",
        "fewer free loadings"
      ), call. = FALSE
    )
  }

  # Check for degrees of freedom
  if(parts$df < 1){
    stop(
      paste(
        "The true model has no degrees of freedom (the model is saturated or",
        "not identified) so population error cannot be added.",
        "\n\nUse a model with fewer parameters (e.g., more variables per factor)"
      ), call. = FALSE
    )
  }

  # Method-specific preparation
  if(method == "cb"){

    # Space of errors that leave the true model as the minimizer
    preparation <- CB_prepare(R, parts)

    # Check for room for error
    if(preparation$dimension < 1){
      stop(
        paste(
          "There is no room for error orthogonal to the true model.",
          "\n\nUse a model with more degrees of freedom"
        ), call. = FALSE
      )
    }

  }else if(method == "hb"){

    # Number of omitted parameters (must be available in the true model)
    control <- HB_counts(model, control, specified)

  }

  # Initialize best solution
  best <- NULL
  best_difference <- Inf

  # Loop through attempts
  for(iteration in seq_len(max_iterations)){

    # Optimize error
    solution <- switch(
      method,
      "cb" = CB_optimize(
        metric, target, preparation, parts, tolerance = tolerance
      ),
      "tkl" = TKL_optimize(
        metric, target, R, parts, diag(model$psi),
        control$minor_factors, control$decay
      ),
      "wb" = WB_optimize(metric, target, R, parts),
      "hb" = HB_optimize(
        metric, target, model, R, parts,
        control$omitted, control$omitted_covariances
      )
    )

    # Check for invalid fit
    if(!isTRUE(solution$valid)){
      next
    }

    # Obtain difference from target
    difference <- abs(solution[[metric]] - target)

    # Check for improvement
    if(difference < best_difference){

      # Store results
      best <- list(fit = solution, iterations = iteration)
      best_difference <- difference

    }

    # Check for an exact solution (target is matched)
    # (attempts continue past `tolerance` until the target is matched so that
    # the achieved fit is as close to `fit_value` as the error allows)
    if(best_difference <= 1e-06){
      break
    }

  }

  # Check for no solution
  if(is.null(best)){
    stop(
      paste(
        "Population error could not be generated with a positive",
        "definite matrix. \n\n",
        "Try increasing the number of iterations using the",
        "`max_iterations` argument."
      ), call. = FALSE
    )
  }

  # Obtain population correlation matrix with error
  R_population <- best$fit$R_population
  dimnames(R_population) <- dimnames(R)

  # Check for convergence
  if(best_difference > tolerance){

    warning(
      paste0(
        "Convergence counter has exceeded its limit. ",
        "There were issues converging to the target ", metric, " = ", target,
        " within `tolerance` = ", tolerance, ". \n\n",
        "Using the solution with the closest ", metric, " = ",
        round(best$fit[[metric]], 4), "."
      ), call. = FALSE
    )

  }

  # Obtain the true model fit to the population correlation matrix with error
  # (from the optimization; the discrepancy is the minimum of the true model)
  known <- best$fit$model

  # Check for improper solutions
  if(isTRUE(known$heywood)){
    warning(
      paste(
        "The true model fit to the population correlation matrix with error",
        "has negative uniquenesses (Heywood case). The amount of error may",
        "be too large for this factor structure"
      ), call. = FALSE
    )
  }

  # Sign indeterminacy is per factor (a factor's loadings can be estimated as
  # a whole with reflected sign) so align each factor to the original loadings
  # (factor correlations are reflected accordingly)
  signs <- sign(colSums(known$lambda * model$loadings))
  signs[signs == 0] <- 1
  estimated_loadings <- known$lambda %*% diag(signs, nrow = length(signs))
  estimated_phi <- diag(signs, nrow = length(signs)) %*% known$phi %*%
    diag(signs, nrow = length(signs))

  # Add names
  factor_names <- paste0("F", seq_len(ncol(model$loadings)))
  dimnames(estimated_loadings) <- list(rownames(R), factor_names)
  dimnames(estimated_phi) <- list(factor_names, factor_names)

  # Obtain fit indices and discrepancy
  fit_vector <- unlist(best$fit[c("CFI", "TLI", "RMSEA", "SRMR")])
  names(fit_vector) <- c("CFI", "TLI", "RMSEA", "SRMR")

  # Obtain error
  error <- R_population - R

  # Generate data
  data <- simulate_error_data(
    R_population, lf_object$parameters, nrow(lf_object$data)
  )

  # Populate results
  results <- list(
    data = data,
    population_correlation = R_population,
    original_correlation = R,
    population_error = list(
      method = method,
      metric = metric,
      fit_value = target,
      achieved = unname(fit_vector[metric]),
      fit = fit_vector,
      delta = best$fit$discrepancy,
      error = error,
      max_error = max(abs(error)),
      R_difference = max(abs(known$sigma - R)),
      loadings = estimated_loadings,
      factor_correlations = estimated_phi,
      iterations = best$iterations
    ),
    original_results = lf_object
  )

  # Add omitted parameters (Hu and Bentler)
  if(method == "hb"){

    omitted <- best$fit$omitted
    dimnames(omitted$loadings) <- list(rownames(R), factor_names)
    dimnames(omitted$factor_covariances) <- list(factor_names, factor_names)
    results$population_error$omitted <- omitted

  }

  # Add class
  class(results) <- c(class(lf_object), "lf_pe")

  # Return results
  return(results)

}

# Input checking ----
#' @noRd
# Updated 20.11.2026
add_population_error_errors <- function(
    lf_object, method, metric, fit_value,
    tolerance, max_iterations, control, leave_cross_loadings,
    method_missing, metric_missing
)
{

  # Check for appropriate class
  # (population error must be added last)
  check_lf_object(
    lf_object, "add_population_error",
    forbid = c("lf_we", "lf_pe", "lf_zipfs"), flat = TRUE
  )

  # Check for population correlation matrix
  if(is.null(lf_object$population_correlation)){
    stop(
      "`lf_object` does not contain a population correlation matrix",
      call. = FALSE
    )
  }

  # Check for method (no default)
  if(method_missing || is.null(method)){
    stop(
      paste0(
        "`method` must be specified. ",
        "\n\nUse `method = \"cb\"`, `method = \"tkl\"`, ",
        "`method = \"wb\"`, or `method = \"hb\"`. ",
        "See `?add_population_error` for details"
      ), call. = FALSE
    )
  }
  type_error(method, "character")
  length_error(method, 1)
  method <- tolower(method)

  # Check for available methods
  available_methods <- c("cb", "tkl", "wb", "hb")
  if(!method %in% available_methods){
    stop(
      paste0(
        "`method` = \"", method, "\" is not available.",
        "\n\nUse ", paste0("`method = \"", available_methods, "\"`", collapse = ", ")
      ), call. = FALSE
    )
  }

  # Set default metric (consistent with the method) and check for metric
  if(metric_missing){
    metric <- switch(
      method,
      "cb" = "CFI", "tkl" = "CFI", "wb" = "RMSEA", "hb" = "SRMR"
    )
  }else{

    # Check for appropriate metric
    type_error(metric, "character")
    length_error(metric, 1)
    metric <- toupper(metric)

    # Check for available metrics
    available_metrics <- c("CFI", "TLI", "RMSEA", "SRMR")
    if(!metric %in% available_metrics){
      stop(
        paste0(
          "`metric` = \"", metric, "\" is not available.",
          "\n\nUse ", paste0("`metric = \"", tolower(available_metrics), "\"`", collapse = ", ")
        ), call. = FALSE
      )
    }

  }

  # Set default fit_value (conventional cutoffs)
  if(is.null(fit_value)){
    fit_value <- switch(
      metric,
      "CFI" = 0.95,
      "TLI" = 0.95,
      "RMSEA" = 0.06,
      "SRMR" = 0.05
    )
  }

  # Check for appropriate fit_value
  type_error(fit_value, "numeric")
  length_error(fit_value, 1)

  # Check for range
  if(!is.finite(fit_value) || fit_value <= 0 || fit_value >= 1){
    stop(
      paste0(
        "`fit_value` = ", fit_value, " is not a valid ", metric, ".",
        "\n\nNumeric `fit_value` must be between 0 and 1 (exclusive)"
      ), call. = FALSE
    )
  }

  # Check tolerance and maximum iterations
  type_error(tolerance, "numeric")
  length_error(tolerance, 1)
  range_error(tolerance, c(0, 1))
  type_error(max_iterations, "numeric")
  length_error(max_iterations, 1)
  range_error(max_iterations, c(1, Inf))

  # Check leave cross-loadings
  type_error(leave_cross_loadings, "logical")
  length_error(leave_cross_loadings, 1)

  # Check control
  if(!is.list(control)){
    stop("`control` must be a list. See `?add_population_error` for options")
  }

  # Set control defaults
  control_defaults <- list(
    minor_factors = max(50, 2 * ncol(lf_object$data)),
    decay = 0.30, omitted = 2, omitted_covariances = 0
  )

  # Check for unknown options
  unknown <- setdiff(names(control), names(control_defaults))
  if(length(unknown) > 0 || (length(control) > 0 && is.null(names(control)))){
    stop(
      paste0(
        "`control` option(s) not available: ",
        paste0("\"", unknown, "\"", collapse = ", "),
        "\n\nUse ", paste0("\"", names(control_defaults), "\"", collapse = ", ")
      )
    )
  }

  # Update defaults
  control_defaults[names(control)] <- control
  control <- control_defaults

  # Check control options
  for(option in names(control)){
    type_error(control[[option]], "numeric")
    length_error(control[[option]], 1)
  }
  range_error(control$minor_factors, c(1, Inf))
  range_error(control$decay, c(0, 1))
  range_error(control$omitted, c(0, Inf))
  range_error(control$omitted_covariances, c(0, Inf))
  for(option in c("omitted", "omitted_covariances")){
    if(control[[option]] %% 1 != 0){
      stop(
        paste0(
          "`control$", option, "` = ", control[[option]],
          " must be a whole number"
        ), call. = FALSE
      )
    }
  }

  # Return inputs
  return(
    list(
      method = method, metric = metric,
      fit_value = fit_value, control = control
    )
  )

}

# Debugging ----
# two_factor <- simulate_factors(
#   factors = 2, variables = 6,
#   loadings = 0.55, cross_loadings = 0.05,
#   correlations = 0.30, sample_size = 1000
# )
# lf_object = two_factor; method = "cb"; metric = "CFI"; fit_value = NULL
# tolerance = 0.01; max_iterations = 10; control = list()
# leave_cross_loadings = FALSE
