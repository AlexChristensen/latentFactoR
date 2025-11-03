#' Adds Population Error to \code{\link[latentFactoR]{simulate_factors}} Data
#'
#' Adds population error to simulated data from \code{\link[latentFactoR]{simulate_factors}}.
#' See examples to get started
#'
#' @param lf_object Data object from \code{\link[latentFactoR]{simulate_factors}}
#'
#' @param cfa_method Character (length = 1).
#' Method to generate population error.
#' Defaults to \code{"ml"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"minres"} --- Minimum residual
#'
#' \item \code{"ml"} --- Maximum likelihood
#'
#' }
#'
#' @param fit Character (length = 1).
#' Fit index to control population error.
#' Defaults to \code{"rmsr"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"cfi"} --- Comparative fit index
#'
#' \item \code{"rmsea"} --- Root mean square error of approximation
#'
#' \item \code{"rmsr"} --- Root mean square residuals
#'
#' \item \code{"raw"} --- Direct application of error
#'
#' }
#'
#' @param misfit Character or numeric (length = 1).
#' Magnitude of error to add.
#' Defaults to \code{"close"}.
#' Available options:
#'
#' \itemize{
#'
#' \item \code{"close"} --- Slight deviations from original population correlation matrix
#'
#' \item \code{"acceptable"} --- Moderate deviations from original population correlation matrix
#'
#' }
#'
#' While numbers can be used, they are \strong{not} recommended. They can be
#' used to specify misfit but the level of misfit will vary depending
#' on the factor structure
#'
#' @param error_method Character (length = 1).
#' Method to control population error.
#' Defaults to \code{"cudeck"}.
#' Description of methods:
#'
#' \itemize{
#'
#' \item \code{"cudeck"} --- Description coming soon... see Cudeck & Browne, 1992
#' for more details
#'
#' \item \code{"yuan"} --- Description coming soon...
#'
#' }
#'
#' @param tolerance Numeric (length = 1).
#' Tolerance of SRMR difference between population error
#' correlation matrix and the original population correlation
#' matrix. Ensures that appropriate population error
#' was added. Similarly, verifies that the MAE of the
#' loadings are not greater than the specified amount,
#' ensuring proper convergence.
#' Defaults to \code{0.01}
#'
#' @param convergence_iterations Numeric (length = 1).
#' Number of iterations to reach parameter convergence
#' within the specified `tolerance`.
#' Defaults to \code{10}
#'
#' @param leave_cross_loadings Boolean.
#' Should cross-loadings be kept?
#' Defaults to \code{FALSE}.
#' Convergence problems can arise if cross-loadings are kept,
#' so setting them to zero is the default. Only set to \code{TRUE}
#' with careful consideration of the structure. Make sure to perform
#' additional checks that the data are adequate
#'
#' @return Returns a list containing:
#'
#' \item{data}{Simulated data from the specified factor model}
#'
#' \item{population_correlation}{Population correlation matrix with local dependence added}
#'
#' \item{population_error}{
#' A list containing the parameters used to generate population error:
#'
#' \itemize{
#'
#' \item \code{error_correlation} --- Correlation matrix with population
#' error added (same as \code{population_correlation})
#'
#' \item \code{fit} --- Fit measure used to control population error
#'
#' \item \code{delta} --- Minimum of the objective function corresponding
#' to the misfit value
#'
#' \item \code{misfit} --- Specified misfit value
#'
#' \item \code{loadings} --- Estiamted CFA loadings after error has been added
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
#' # Add small population error using Cudeck method
#' two_factor_Cudeck <- add_population_error(
#'   lf_object = two_factor,
#'   cfa_method = "minres",
#'   fit = "rmsr", misfit = "close",
#'   error_method = "cudeck"
#' )
#'
#' # Add small population error using Yuan method
#' two_factor_Yuan <- add_population_error(
#'   lf_object = two_factor,
#'   cfa_method = "minres",
#'   fit = "rmsr", misfit = "close",
#'   error_method = "yuan"
#' )
#'
#' @author
#' {\code{bifactor}} authors \cr
#' Marcos Jimenez,
#' Francisco J. Abad,
#' Eduardo Garcia-Garzon,
#' Vithor R. Franco,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
#'
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
#' Cudeck, R., & Browne, M.W. (1992).
#' Constructing a covariance matrix that yields a specified minimizer and a specified minimum discrepancy function value.
#' \emph{Psychometrika}, \emph{57}, 357–369.
#'
#' @export
#'
# Add population error to simulated data
# Updated 03.11.2025
add_population_error <- function(
    lf_object,
    cfa_method = c("minres", "ml"),
    fit = c("cfi", "rmsea", "rmsr", "raw"),
    misfit = c("close", "acceptable"),
    error_method = c("cudeck", "yuan"),
    tolerance = 0.01,
    convergence_iterations = 10,
    leave_cross_loadings = FALSE
)
{

  # Check for appropriate class
  if(!is(lf_object, "lf_simulate")){

    # Produce error
    stop(
      paste(
        "`lf_object` input is not class \"lf_simulate\" from the `simulate_factors` function.",
        "\n\nInput class(es) of current `lf_object`:",
        paste0("\"", class(lf_object), "\"", collapse = ", "),
        "\n\nUse `simulate_factors` to generate your data to input into this function"
      )
    )

  }

  # Obtain parameters from simulated data
  parameters <- lf_object$parameters

  # Check for missing CFA method
  if(missing(cfa_method)){
    cfa_method <- "ml"
  }else{cfa_method <- tolower(match.arg(cfa_method))}

  # Check for missing fit
  if(missing(fit)){
    fit <- "rmsr"
  }else{fit <- tolower(match.arg(fit))}

  # Check for missing misfit
  if(missing(misfit)){
    misfit <- "close"
  }

  # Check for missing error method
  if(missing(error_method)){
    error_method <- "cudeck"
  }else{error_method <- tolower(match.arg(error_method))}

  # Check for appropriate misfit
  length_error(misfit, 1);

  # Obtain loadings
  loadings <- parameters$loadings

  # Determine error structure
  if(is(lf_object, "lf_ld")){

    # Obtain correlated residuals
    errors <- lf_object$population_correlation -
      lf_object$original_results$population_correlation

  }else{

    # Initialize correlated residuals
    errors <- matrix(
      0, nrow = ncol(lf_object$data),
      ncol = ncol(lf_object$data)
    )

  }

  # Check for whether cross-loadings should remain
  if(!isTRUE(leave_cross_loadings)){

    # Set sequence of variables for each factor
    end_variables <- cumsum(parameters$variables)
    start_variables <- (end_variables + 1) - parameters$variables

    # Loop through loadings
    for(i in 1:ncol(loadings)){

      # Set cross-loadings to zero
      loadings[
        start_variables[i]:end_variables[i],
        -i
      ] <- 0

    }

  }else if(is(lf_object, "lf_cl")){

    # Set factor correlations
    factor_correlations <- parameters$factor_correlations

    # Check communalities
    communalities <- diag(
      loadings %*%
        factor_correlations %*%
        t(loadings)
    )

    # Initialize break count
    break_count <- 0

    # Loop through until all communalities < 0.80
    while(any(communalities >= 0.80)){

      # Increase break count
      break_count <- break_count + 1

      # Message about adjustment
      if(break_count == 1){

        message(
          paste(
            "Communalities for the following variable(s) were >= 0.80:",
            paste0(
              which(communalities >= 0.80),
              collapse = ", "
            ),
            "\nThe dominant loadings on these variable(s) were decreased",
            "incrementally by 0.01 until their communalities were < 0.80"
          )
        )

      }

      # Identify loadings with communalities greater than 0.90
      target_loadings <- matrix(
        loadings[which(communalities >= 0.80),],
        ncol = ncol(loadings),
        byrow = FALSE
      )

      # Decrease maximum loadings by 0.01
      replace_loadings <- matrix(
        apply(target_loadings, 1, function(x){

          # Obtain signs
          signs <- sign(x)

          # Compute absolute max
          x <- abs(x)

          # Decrease by 0.01
          x[which.max(x)] <- x[which.max(x)] - 0.01

          # Add back signs
          x <- x * signs

          # Return loadings
          return(x)

        }),
        ncol = ncol(loadings),
        byrow = TRUE
      )

      # Replace loadings
      loadings[which(communalities >= 0.80),] <- replace_loadings

      # Check communalities
      communalities <- diag(
        loadings %*%
          factor_correlations %*%
          t(loadings)
      )

    }

  }

  # Re-compute population correlation matrix
  lf_object$population_correlation <- loadings %*%
    parameters$factor_correlations %*%
    t(loadings)

  # Obtain uniquenesses and put them on the diagonal of errors
  diag(errors) <- 1 - diag(lf_object$population_correlation)

  # Reset diagonal of population correlation matrix
  diag(lf_object$population_correlation) <- 1

  # Initialize positive definite and convergence
  positive_definite <- FALSE
  convergence <- FALSE

  # Initialize counts so it doesn't get stuck
  pd_stuck_count <- 0
  convergence_stuck_count <- 0

  # Initialize maximum absolute residual vector
  residual <- rep(NA, length = convergence_iterations)

  # Initialize previous minimum
  previous_minimum <- Inf

  # Ensure proper convergence
  while(!convergence){

    # Try to get positive definite matrix
    while(!positive_definite){

      # Obtain population error
      if(error_method == "cudeck"){

        # Using Cudeck method
        # From {bifactor} version 0.1.0
        # See `utils-latentFactoR`
        population_error <- try(
          cudeck(
            R = lf_object$population_correlation,
            lambda = loadings,
            Phi = parameters$factor_correlations,
            Psi = errors, fit = fit,
            misfit = misfit, method = cfa_method
          ),
          silent = TRUE
        )

      }else if(error_method == "yuan"){

        # Using Yuan method
        # From {bifactor} version 0.1.0
        # See `utils-latentFactoR`
        population_error <- try(
          yuan(
            R = lf_object$population_correlation,
            lambda = loadings,
            Phi = parameters$factor_correlations,
            Psi = errors, fit = fit,
            misfit = misfit, method = cfa_method
          ),
          silent = TRUE
        )

      }

      # Check for error
      if(is(population_error, "try-error")){

        # Increase positive definite stuck count
        pd_stuck_count <- pd_stuck_count + 1

        # Check if a break is necessary
        if(pd_stuck_count >= convergence_iterations){

          # Stop and tell user to increase convergence iterations
          stop(
            paste(
              "Convergence counter has exceeded its limit.",
              "There were issues converging the model with proper",
              "population error. \n\n",
              "Population error could not converge with a positive",
              "definite matrix. \n\n",
              "Try increasing the number of iterations for convergence",
              "using the `convergence_iterations` argument."
            )
          )

        }

      }else if(
        any(
          eigen(
            x = population_error$R_error,
            symmetric = TRUE,
            only.values = TRUE
          )$values < .Machine$double.eps
        )
      ){

        # Increase positive definite stuck count
        pd_stuck_count <- pd_stuck_count + 1

        # Check if a break is necessary
        if(pd_stuck_count >= convergence_iterations){

          # Stop and tell user to increase convergence iterations
          stop(
            paste(
              "Convergence counter has exceeded its limit.",
              "There were issues converging the model with proper",
              "population error. \n\n",
              "Population error could not converge with a positive",
              "definite matrix. \n\n",
              "Try increasing the number of iterations for convergence",
              "using the `convergence_iterations` argument."
            )
          )

        }

      }else{

        # Set positive definite to be TRUE
        positive_definite <- TRUE

      }

    }

    # Obtain population error correlation matrix
    error_correlation <- population_error$R_error

    # Add row and column names to population error correlation matrix
    colnames(error_correlation) <- paste0("V", 1:ncol(error_correlation))
    row.names(error_correlation) <- colnames(error_correlation)

    # Specify the CFA model
    # (set 1 to estimate the parameter and 0 to fix it to zero)
    target <- ifelse(loadings != 0, 1, 0) # For loadings
    target_phi <- ifelse(parameters$factor_correlations != 0, 1, 0)  # For factor correlations
    target_psi <- ifelse(errors != 0, 1, 0)  # For uniquenesses and errors

    # Fit the CFA
    cfa <- CFA(
      S = error_correlation,
      target = target,
      targetphi = target_phi,
      targetpsi = target_psi,
      method = cfa_method
    )

    # Check SRMR
    SRMR <- sqrt(mean(cfa$residuals[lower.tri(cfa$residuals)]^2))

    # Obtain SRMR difference
    SRMR_difference <- abs(SRMR - population_error$misfit)

    # Compute the largest absolute residual
    max_abs_res <- max(abs(cfa$residuals))

    # Cutoff for the maximum absolute residual
    max_res <- switch(
      as.character(misfit),
      "close" = 0.10,
      "acceptable" = 0.15,
      as.numeric(misfit) + 0.05
    )

    # Ensure same order of loadings
    error_loadings <- cfa$lambda

    # Sometimes the loadings can be in opposite directions
    # Check that...
    # Get signs
    error_signs <- sign(error_loadings)
    loading_signs <- sign(loadings)
    if(any(error_signs != loading_signs)){

      # Get non-zero signs
      non_zero <- error_signs != 0

      # Flip signs
      error_loadings[non_zero] <- -error_loadings[non_zero]

    }

    # Obtain difference between error and population loadings
    MAE <- mean(abs(error_loadings - loadings))

    # Check for convergence
    if(
      SRMR_difference <= tolerance &
      MAE <= tolerance &
      max_abs_res < max_res
    ){convergence <- TRUE}

    # Increase convergence stuck count
    convergence_stuck_count <- convergence_stuck_count + 1

    # Add residual
    residual[convergence_stuck_count] <- max_abs_res

    # Check for minimum
    if(min(residual, na.rm = TRUE) < previous_minimum){

      # Update minimum
      previous_minimum <- round(min(residual, na.rm = TRUE), 3)

      # Store results
      pe_stored <- population_error
      loadings_stored <- error_loadings
      SRMR_stored <- SRMR_difference
      MAE_loadings <- MAE

    }

    # Check if a break is necessary
    if(convergence_stuck_count >= convergence_iterations){

      warning(
        paste(
          "Convergence counter has exceeded its limit.",
          "There were issues converging the model with proper",
          "population error. \n\n",
          "Using the solution with the minimum maximum absolute residual =",
          previous_minimum
        )
      )

      # Restore stored values
      population_error <- pe_stored
      error_loadings <- loadings_stored
      SRMR_difference <- SRMR_stored
      MAE <- MAE_loadings

      # Break out of loop
      break

    }

  }

  # Re-estimate data
  ## Cholesky decomposition
  cholesky <- chol(population_error$R_error)

  ## Obtain sample size and total variables
  sample_size <- nrow(lf_object$data)
  total_variables <- ncol(lf_object$data)

  ## Generate data
  data <- mvtnorm::rmvnorm(sample_size, sigma = diag(total_variables))

  ## Make data based on factor structure
  data <- data %*% cholesky

  ## Set variable categories, limit, and skew
  variable_categories <- parameters$categories
  categorical_limit <- parameters$categorical_limit
  skew <- parameters$skew

  ## Ensure appropriate type and length for categories
  type_error(variable_categories, "numeric")
  length_error(variable_categories, c(1, total_variables))

  ## Identify categories to variables
  if(length(variable_categories) == 1){
    variable_categories <- rep(variable_categories, total_variables)
  }

  ## Check for categories greater than categorical limit and not infinite
  if(any(variable_categories > categorical_limit & !is.infinite(variable_categories))){

    ## Make variables with categories greater than 7 (or categorical_limit) continuous
    variable_categories[
      variable_categories > categorical_limit & !is.infinite(variable_categories)
    ] <- Inf

  }

  ## Find categories
  if(any(variable_categories <= categorical_limit)){

    ## Target columns to categorize
    columns <- which(variable_categories <= categorical_limit)

    ## Set skew
    if(length(skew) != length(columns)){
      skew <- sample(skew, length(columns), replace = TRUE)
    }

    ## Loop through columns
    for(i in columns){

      data[,i] <- categorize(
        data = data[,i],
        categories = variable_categories[i],
        skew_value = skew[i]
      )

    }

  }

  ## Add column names to data
  colnames(data) <- paste0(
    "V", formatC(
      x = 1:total_variables,
      digits = floor(log10(total_variables)),
      flag = "0", format = "d"
    )
  )

  ## Populate results
  results <- list(
    data = data,
    population_correlation = population_error$R_error,
    original_correlation = lf_object$population_correlation
  )

  # Add parameters from population error
  error_parameters <- list(
    error_correlation = population_error$R_error,
    fit = population_error$fit,
    delta = population_error$delta,
    misfit = population_error$misfit,
    loadings = error_loadings,
    MAX_residual = max_abs_res,
    SRMR_difference = SRMR_difference,
    MAE_loadings = MAE
  )

  # Add population error parameters to results
  results$population_error <- error_parameters

  # Add original results
  results$original_results <- lf_object

  # Add class
  class(results) <- c(class(lf_object), "lf_pe")

  # Message to cite {bifactor}
  message(
    paste0(
      "Please cite the {bifactor} package:\n\n",
      "Jimenez, M., Abad, F. J., Garcia-Garzon, E., Garrido, L. E., Franco, V. R. (2022). ",
      styletext("bifactor: Exploratory factor and bi-factor modeling with multiple general factors. ", defaults = "italics"),
      "R package version 0.1.0. ",
      "Retrieved from https://github.com/Marcosjnez/bifactor"
    )
  )

  # Return results
  return(results)

}

# Debugging ----
# two_factor <- simulate_factors(
#   factors = 2, variables = 6,
#   loadings = 0.55, cross_loadings = 0.05,
#   correlations = 0.30, sample_size = 1000
# )
# lf_object = two_factor; cfa_method = "minres";
# fit = "cfi"; misfit = "acceptable"; error_method = "cudeck";
# tolerance = 0.01; convergence_iterations = 10;
# leave_cross_loadings = FALSE

