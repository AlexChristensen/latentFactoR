#' Adds Methods Factors to \code{\link[latentFactoR]{simulate_factors}} Data
#'
#' Adds methods factors to simulated data from \code{\link[latentFactoR]{simulate_factors}}.
#' See examples to get started
#'
#' @param lf_object Data object from \code{\link[latentFactoR]{simulate_factors}}.
#' Data \strong{must} be categorical. If data are not categorical, then
#' there function with throw an error
#'
#' @param proportion_negative Numeric (length = 1 or \code{factors}).
#' Proportion of variables that should have negative (or flipped) dominant loadings across all
#' or each factor. Accepts number of variables as well.
#' The first variables on each factor, up to the corresponding proportion, will be
#' flipped. Set to \code{0} to not have any loadings flipped.
#' Defaults to \code{0.50}
#'
#' @param proportion_negative_range Numeric (length = 2).
#' Range of proportion of variables that are randomly selected from
#' a uniform distribution. Accepts number of number of variables as well.
#' Defaults to \code{NULL}
#'
#' @param methods_factors Numeric
#'
#' @param methods_loadings Numeric
#'
#' @param methods_loadings_range Numeric
#'
#' @param methods_correlations Numeric
#'
#' @param methods_correlations_range Numeric
#'
#' @param seed Numeric (length = 1).
#' Seed for reproducible results. The same \code{seed} gives the same results
#' (and does not change the random number stream of the R session).
#' Defaults to \code{NULL} (random seed)
#'
#' @return Returns a list containing:
#'
#' \item{data}{Biased data simulated data from the specified factor model}
#'
#' \item{unbiased_data}{The corresponding unbiased data prior to replacing values
#' to generate the (biased) \code{data}}
#'
#' \item{population_correlation}{Population correlation matrix with the
#' methods factors}
#'
#' \item{parameters}{Bias-adjusted parameters of the \code{lf_object} input into function}
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
#'   sample_size = 1000, # number of cases = 1000
#'   variable_categories = 5 # 5-point Likert scale
#' )
#'
#' # Add methods factors
#' two_factor_methods_effect <- add_method_factors(
#'   lf_object = two_factor,
#'   proportion_negative = 0.50,
#'   methods_loadings = 0.20,
#'   methods_loadings_range = 0.10
#' )
#'
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
#'
#' @references
#' Garcia-Pardina, A., Abad, F. J., Christensen, A. P., Golino, H., & Garrido, L. E. (2024).
#' Dimensionality assessment in the presence of wording effects: A network psychometric and factorial approach.
#' \emph{Behavior Research Methods}.
#'
#' @export
#'
# Add methods factors to simulated data ----
# Updated 12.07.2026
add_method_factors <- function(
    lf_object,
    proportion_negative = 0.50,
    proportion_negative_range = NULL,
    methods_factors,
    methods_loadings,
    methods_loadings_range = 0.00,
    methods_correlations,
    methods_correlations_range = NULL,
    seed = NULL
)
{

  # Set seed (reproducible results)
  seed_token <- lf_seed_start(seed)
  on.exit(lf_seed_end(seed_token), add = TRUE)


  # Check for appropriate class
  # (method factors are added to the substantive factors, so they must be
  # added before local dependence, wording effects, and population error)
  check_lf_object(
    lf_object, "add_method_factors",
    forbid = c("lf_me", "lf_ld", "lf_we", "lf_pe", "lf_zipfs"), flat = TRUE
  )

  # Obtain parameters from simulated data
  parameters <- lf_object$parameters

  # Check for missing methods factors
  if(missing(methods_factors)){methods_factors <- parameters$factors}
  if(missing(methods_correlations)){
    methods_correlations <- parameters$factor_correlations
    methods_correlations[] <- 1
  }

  # Check inputs
  inputs <- add_method_factors_errors(
    parameters,
    proportion_negative, proportion_negative_range,
    methods_correlations, methods_correlations_range
  )

  # Collect inputs
  proportion_negative <- inputs$proportion_negative
  methods_correlations <- inputs$methods_correlations

  # Obtain substantive loadings
  loadings <- parameters$loadings

  # Obtain variables
  variables <- parameters$variables

  # Set sequence of variables for each factor
  variable_sequence <- factor_variable_sequence(parameters$variables)
  start_variables <- variable_sequence$start
  end_variables <- variable_sequence$end

  # Flip dominant loadings
  loadings <- flip_dominant_loadings(loadings, variables, proportion_negative)

  # Ensure proper signs for skew
  signs <- dominant_signs(loadings, variables)

  # Obtain skews
  skews <- parameters$skew

  # Handle skew signs and re-assign
  parameters$skew <- handle_skew_signs(
    skews = skews, signs = signs
  )

  # Determine number of wording factors
  if(length(methods_factors) == 1){

    # Determine whether proportion was used
    if(methods_factors < 1){

      # Check for error in range
      range_error(methods_factors, c(0, 1)) # object range error

      # Determine number of methods factors
      methods_factors <- round(methods_factors * parameters$factors)

    }

    # Determine whether methods factors is not one
    if(methods_factors != 1){
      methods_factors <- shuffle(seq_len(parameters$factors), size = methods_factors)
    }

  }

  # Check for error in length and range
  length_error(methods_factors, seq_len(parameters$factors)) # object length error
  range_error(methods_factors, c(0, parameters$factors)) # object range error

  # Initialize checks
  check_eigenvalues <- TRUE
  check_communalities <- TRUE

  # Initialize iteration counter (guards against an infinite loop when the
  # supplied structure leaves no randomness to escape an inadmissible
  # solution, e.g. `methods_loadings`/`methods_correlations` are fixed matrices)
  iterations <- 0
  max_iterations <- 1000

  # Run through loop
  while(isTRUE(check_eigenvalues) | isTRUE(check_communalities)){

    # Increment iteration counter and check against maximum
    iterations <- iterations + 1
    if(iterations > max_iterations){
      stop(
        paste0(
          "Could not find an admissible (positive semi-definite) solution after ",
          max_iterations, " attempts. The supplied `methods_loadings` and/or ",
          "`methods_correlations` may be jointly inadmissible. Try relaxing fixed ",
          "matrices into ranges or adjusting values"
        )
      )
    }

    # Populate methods loadings matrix
    if(!is(methods_loadings, "matrix")){

      # Get dimensions of original loadings matrix
      loading_dimensions <- dim(loadings)

      # Create methods loadings matrix
      methods_loadings_matrix <- matrix(
        0, nrow = loading_dimensions[1], ncol = loading_dimensions[2]
      )

      # Identify dominant loadings
      if(length(methods_loadings) == 1){methods_loadings <- rep(methods_loadings, parameters$factors)}
      if(length(methods_loadings_range) == 1){methods_loadings_range <- rep(methods_loadings_range, parameters$factors)}

      # Populate factor loadings
      for(i in methods_factors){

        # Generate loadings from uniform distribution
        methods_loadings_matrix[
          start_variables[i]:end_variables[i], i # dominant loadings
        ] <- runif_xoshiro(
          variables[i], # target variables
          min = methods_loadings[i] - methods_loadings_range[i] * 0.50,
          max = methods_loadings[i] + methods_loadings_range[i] * 0.50
        )

      }

    }else{# Input is already loadings matrix
      methods_loadings_matrix <- methods_loadings
    }

    # Methods factor correlations
    if(length(methods_correlations) == 1){

      # Generate correlation matrix
      methods_correlation_matrix <- matrix(
        data = methods_correlations,
        nrow = parameters$factors, ncol = parameters$factors
      )

    }else{# Input is already correlation matrix
      methods_correlation_matrix <- methods_correlations
    }

    # Ensure diagonal of methods correlation matrix is 1
    diag(methods_correlation_matrix) <- 1

    # Combine loadings matrices
    combined_loadings <- cbind(loadings, methods_loadings_matrix)

    # Combine correlation matrices
    combined_correlations <- matrix(
      0, nrow = parameters$factors * 2,
      ncol = parameters$factors * 2
    )

    # Populate combined correlation matrix
    ## Substantive factors
    combined_correlations[
      seq_len(parameters$factors), seq_len(parameters$factors)
    ] <- parameters$factor_correlations
    ## Methods factors
    combined_correlations[
      seq_len(parameters$factors) + parameters$factors,
      seq_len(parameters$factors) + parameters$factors
    ] <- methods_correlation_matrix

    # Create population correlation matrix
    population_correlation <- combined_loadings %*%
      combined_correlations %*%
      t(combined_loadings)

    # Check communalities
    check_communalities <- any(diag(population_correlation) > 0.90)

    # Ensure diagonal of correlation matrix is 1
    diag(population_correlation) <- 1

    # Check eigenvalues
    check_eigenvalues <- any(matrix_eigenvalues(population_correlation) <= 0)

  }

  # Generate data (categorize and skew variables)
  generated <- generate_lf_data(
    population_correlation, nrow(lf_object$data), parameters$categories,
    parameters$categorical_limit, parameters$skew
  )
  data <- generated$data

  # Add names to population correlation matrix
  dimnames(population_correlation) <- list(colnames(data), colnames(data))

  # Update parameters
  parameters$loadings <- combined_loadings
  parameters$cross_loadings <- c(
    parameters$cross_loadings, rep(0, length(parameters$cross_loadings))
  )
  parameters$factor_correlations <- combined_correlations
  parameters$skew <- generated$skew

  # Populate results
  results <- list(
    data = data,
    unbiased_data = lf_object$data,
    population_correlation = population_correlation,
    parameters = parameters,
    original_results = lf_object
  )

  # Add class
  class(results) <- c(class(lf_object), "lf_me")

  # Return results
  return(results)

}

# Input checking ----
#' @noRd
# Updated 12.07.2026
add_method_factors_errors <- function(
    parameters,
    proportion_negative, proportion_negative_range = NULL,
    methods_correlations, methods_correlations_range = NULL
)
{

  # Check for methods correlations range
  if(!is.null(methods_correlations_range)){
    type_error(methods_correlations_range, "numeric") # object type error
    length_error(methods_correlations_range, 2) # object length error
    range_error(methods_correlations_range, c(-1, 1)) # object range error
    methods_correlations <- draw_correlation_matrix(
      parameters$factors, methods_correlations_range
    )
  }

  # Resolve proportions (or numbers of variables)
  proportion_negative <- resolve_proportion(
    proportion_negative, proportion_negative_range,
    parameters$variables, "proportion_negative"
  )

  # Return checked input
  return(
    list(
      proportion_negative = proportion_negative,
      methods_correlations = methods_correlations
    )
  )

}
