#' Simulates Latent Factor Data
#'
#' @description
#' Simulates data from a latent factor model based on many
#' manipulable parameters. Parameters do not have default values and
#' must each be set. See examples to get started
#'
#' The population correlation matrix (variable-level) implied by the factor model is:
#'
#' \deqn{\boldsymbol{\Sigma} = \boldsymbol{\Lambda} \boldsymbol{\Phi} \boldsymbol{\Lambda}' + \boldsymbol{\Theta}}{Sigma = Lambda * Phi * Lambda' + Theta}
#'
#' where \eqn{\boldsymbol{\Lambda}}{Lambda} is the loading matrix (variables x factors; \code{loadings}
#' and \code{cross_loadings}), \eqn{\boldsymbol{\Phi}}{Phi} is the factor correlation matrix
#' (\code{correlations}), and \eqn{\boldsymbol{\Theta}}{Theta} is a diagonal matrix of uniquenesses
#' that ensures unit variance for each variable
#'
#' @param factors Numeric (length = 1).
#' Number of factors
#'
#' @param variables Numeric (length = 1 or \code{factors}).
#' Number of variables per factor.
#' Can be a single value or as many values as there are factors.
#' Minimum three variables per factor
#'
#' @param variables_range Numeric (length = 2).
#' Range of variables to randomly select from a random uniform distribution.
#' Minimum three variables per factor
#'
#' @param loadings Numeric or matrix (length = 1, \code{factors}, total number of variables (\code{factors} x \code{variables}), or total number of variables x \code{factors}).
#' Loadings drawn from a random uniform distribution using +/- 0.10 of value input.
#' Can be a single value or as many values as there are factors (corresponding to the factors).
#' Can also be a loading matrix. Columns must match factors and rows must match total variables (\code{factors} x \code{variables})
#' General effect sizes range from small (0.40), moderate (0.55), to large (0.70)
#'
#' @param loadings_range Numeric (length = 2).
#' Range of loadings to randomly select from a random uniform distribution.
#' General effect sizes range from small (0.40), moderate (0.55), to large (0.70)
#'
#' @param cross_loadings Numeric or matrix (length = 1, \code{factors}, or total number of variables x \code{factors}).
#' Cross-loadings drawn from a random normal distribution with a mean of 0 and standard deviation of value input.
#' Can be a single value or as many values as there are factors (corresponding to the factors).
#' Can also be a loading matrix. Columns must match factors and rows must match total variables (\code{factors} x \code{variables})
#'
#' @param cross_loadings_range Numeric (length = 2).
#' Range of cross-loadings to randomly select from a random uniform distribution
#'
#' @param correlations Numeric (length = 1 or \code{factors} x \code{factors}).
#' Can be a single value that will be used for all correlations between factors.
#' Can also be a square matrix (\code{factors} x \code{factors}).
#' General effect sizes range from orthogonal (0.00), small (0.30), moderate (0.50), to large (0.70)
#'
#' @param correlations_range Numeric (length = 2).
#' Range of correlations to randomly select from a random uniform distribution.
#' General effect sizes range from orthogonal (0.00), small (0.30), moderate (0.50), to large (0.70)
#'
#' @param sample_size Numeric (length = 1).
#' Number of cases to generate from a random multivariate normal distribution using
#' \code{\link[mvtnorm]{rmvnorm}}
#'
#' @param variable_categories Numeric (length = 1 or total variables (\code{factors} x \code{variables})).
#' Number of categories for each variable. \code{Inf} is used for continuous variables; otherwise,
#' values reflect number of categories
#'
#' @param categorical_limit Numeric (length = 1).
#' Values greater than input value are considered continuous.
#' Defaults to \code{7} meaning that 8 or more categories are considered continuous
#' (i.e., variables are \emph{not} categorized from continuous to categorical)
#'
#' @param skew Numeric (length = 1 or categorical variables).
#' Skew to be included in categorical variables. It is randomly sampled from provided values.
#' Can be a single value or as many values as there are (total) variables.
#' Current skew implementation is between -2 and 2 in increments of 0.05.
#' Skews that are not in this sequence will be converted to their nearest
#' value in the sequence. Not recommended to use with \code{variables_range}.
#' Future versions will incorporate finer skews
#'
#' @param skew_range Numeric (length = 2).
#' Randomly selects skews within in the range.
#' Somewhat redundant with \code{skew} but more flexible.
#' Compatible with \code{variables_range}
#'
#' @return Returns a list containing:
#'
#' \item{data}{Simulated data from the specified factor model}
#'
#' \item{population_correlation}{Population correlation matrix}
#'
#' \item{continuous_data}{Simulated data before categories and skew were applied
#' (i.e., the underlying continuous data)}
#'
#' \item{parameters}{
#' A list containing the parameters used to generate the data:
#'
#' \itemize{
#'
#' \item \code{factors} --- Number of factors
#'
#' \item \code{variables} --- Variables on each factor
#'
#' \item \code{loadings} --- Loading matrix
#'
#' \item \code{cross_loadings} --- Cross-loadings used for each factor
#'
#' \item \code{factor_correlations} --- Correlations between factors
#'
#' \item \code{categories} --- Categories for each variable
#'
#' \item \code{categorical_limit} --- Category limit used to determine continuous variables
#'
#' \item \code{skew} --- Skew for each variable
#'
#' }
#'
#' }
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
#' # Randomly vary loadings
#' two_factor_loadings <- simulate_factors(
#'   factors = 2, # factors = 2
#'   variables = 6, # variables per factor = 6
#'   loadings_range = c(0.30, 0.80), # loadings between = 0.30 to 0.80
#'   cross_loadings = 0.05, # cross-loadings N(0, 0.05)
#'   correlations = 0.30, # correlation between factors = 0.30
#'   sample_size = 1000 # number of cases = 1000
#' )
#'
#' # Generate dichotomous data
#' two_factor_dichotomous <- simulate_factors(
#'   factors = 2, # factors = 2
#'   variables = 6, # variables per factor = 6
#'   loadings = 0.55, # loadings between = 0.45 to 0.65
#'   cross_loadings = 0.05, # cross-loadings N(0, 0.05)
#'   correlations = 0.30, # correlation between factors = 0.30
#'   sample_size = 1000, # number of cases = 1000
#'   variable_categories = 2 # dichotomous data
#' )
#'
#' # Generate dichotomous data with skew
#' two_factor_dichotomous_skew <- simulate_factors(
#'   factors = 2, # factors = 2
#'   variables = 6, # variables per factor = 6
#'   loadings = 0.55, # loadings between = 0.45 to 0.65
#'   cross_loadings = 0.05, # cross-loadings N(0, 0.05)
#'   correlations = 0.30, # correlation between factors = 0.30
#'   sample_size = 1000, # number of cases = 1000
#'   variable_categories = 2, # dichotomous data
#'   skew = 1 # all variables with have a skew of 1
#' )
#'
#' # Generate dichotomous data with variable skew
#' two_factor_dichotomous_skew <- simulate_factors(
#'   factors = 2, # factors = 2
#'   variables = 6, # variables per factor = 6
#'   loadings = 0.55, # loadings between = 0.45 to 0.65
#'   cross_loadings = 0.05, # cross-loadings N(0, 0.05)
#'   correlations = 0.30, # correlation between factors = 0.30
#'   sample_size = 1000, # number of cases = 1000
#'   variable_categories = 2, # dichotomous data
#'   skew_range = c(-2, 2) # skew = -2 to 2 (increments of 0.05)
#' )
#'
#' @author
#' Maria Dolores Nieto Canaveras <mnietoca@nebrija.es>,
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Hudson Golino <hfg9s@virginia.edu>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
#'
#' @references
#' Garrido, L. E., Abad, F. J., & Ponsoda, V. (2011). \cr
#' Performance of Velicer’s minimum average partial factor retention method with categorical variables. \cr
#' \emph{Educational and Psychological Measurement}, \emph{71}(3), 551-570.
#'
#' Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., ... & Martinez-Molina, A. (2020).
#' Investigating the performance of exploratory graph analysis and traditional techniques to identify the number of latent factors: A simulation and tutorial.
#' \emph{Psychological Methods}, \emph{25}(3), 292-320.
#'
#' @export
#'
# Main factor simulation function
# Updated 12.07.2026
simulate_factors <- function(
    factors, variables, variables_range = NULL,
    loadings, loadings_range = NULL,
    cross_loadings, cross_loadings_range = NULL,
    correlations, correlations_range = NULL,
    sample_size, variable_categories = Inf,
    categorical_limit = 7,
    skew = 0, skew_range = NULL
)
{

  # Check inputs
  inputs <- simulate_factors_errors(
    factors, variables, variables_range,
    loadings, loadings_range,
    cross_loadings, cross_loadings_range,
    correlations, correlations_range,
    sample_size, variable_categories,
    categorical_limit, skew, skew_range
  )

  # Collect inputs
  variables <- inputs$variables; total_variables <- inputs$total_variables
  loadings <- inputs$loadings; cross_loadings <- inputs$cross_loadings
  correlations <- inputs$correlations; skew <- inputs$skew

  # Ensure appropriate lengths
  length_error(correlations, c(1, factors * factors))

  # Ensure appropriate ranges
  range_error(loadings, c(-1, 1))
  range_error(cross_loadings, c(-1, 1));
  range_error(correlations, c(-1, 1))

  # Ensure appropriate types
  if(!is(loadings, "matrix")){type_error(loadings, "numeric")}
  if(!is(cross_loadings, "matrix")){type_error(cross_loadings, "numeric")}
  if(!is(correlations, "matrix")){type_error(correlations, "numeric")}

  # Set factor sequence
  factor_sequence <- seq_len(factors)

  # Identify variables structure on factors
  if(length(variables) == 1){variables <- rep(variables, factors)}

  # Create starting and ending of variable sequences
  variable_sequence <- factor_variable_sequence(variables)
  start_variables <- variable_sequence$start
  end_variables <- variable_sequence$end

  # Verify generalizable lengths
  if(length(loadings) == 1){loadings <- rep(loadings, factors)}
  if(length(cross_loadings) == 1){cross_loadings <- rep(cross_loadings, factors)}
  if(length(variable_categories) == 1){
    variable_categories <- rep(variable_categories, total_variables)
  }

  # Initialize checks
  check_eigenvalues <- TRUE
  check_communalities <- TRUE

  # Initialize iteration counter (guards against an infinite loop when the supplied
  # structure leaves no randomness to escape an inadmissible solution, e.g. `loadings`
  # and `correlations` are both fixed matrices rather than drawn from ranges)
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
          max_iterations, " attempts. The supplied `loadings` and `correlations` may be ",
          "jointly inadmissible (e.g. communalities exceeding .90 or a non-positive-definite ",
          "correlation matrix). Try relaxing fixed matrices into ranges or adjusting values"
        )
      )
    }

    # Generate loadings matrix
    if(!is(loadings, "matrix")){

      # Initialize loading matrix
      loading_matrix <- matrix(data = 0, nrow = total_variables, ncol = factors)

      # Populate factor loadings
      for(i in factor_sequence){

        # Populate dominant loadings
        if(is.null(loadings_range)){

          # Generate loadings from uniform distribution
          loading_matrix[start_variables[i]:end_variables[i], i] <-
          runif_xoshiro(variables[i], min = loadings[i] - 0.10, max = loadings[i] + 0.10)

        }else{

          # Accept loadings from range (generated earlier)
          loading_matrix[start_variables[i]:end_variables[i], i] <-
          loadings[start_variables[i]:end_variables[i]]

        }

        # Add cross-loadings (if more than one factor)
        if(factors > 1){

          # Loop through factors
          for(j in factor_sequence){

            # Ignore dominant factor
            if(i != j){

              # Check for range of cross-loadings
              if(!is.null(cross_loadings_range)){
                loading_matrix[start_variables[j]:end_variables[j], i] <-
                runif_xoshiro(
                  variables[j], min = min(cross_loadings_range),
                  max = max(cross_loadings_range)
                )
              }else{
                loading_matrix[start_variables[j]:end_variables[j], i] <-
                rnorm_ziggurat(variables[j]) * cross_loadings[j]
              }

            }

          }

        }

      }

    }else{# Input is already loadings matrix
      loading_matrix <- loadings
    }

    # Factor correlations
    if(length(correlations) == 1){

      # Generate correlation matrix
      correlation_matrix <- matrix(data = correlations, nrow = factors, ncol = factors)

    }else{# Input is already correlation matrix
      correlation_matrix <- correlations
    }

    # Ensure diagonal of correlation matrix is 1
    diag(correlation_matrix) <- 1

    # Assume always positive correlations (for now)
    correlation_matrix <- abs(correlation_matrix)

    # Create population correlation matrix
    population_correlation <- loading_matrix %*% tcrossprod(correlation_matrix, loading_matrix)

    # Check communalities
    check_communalities <- any(diag(population_correlation) > 0.90)

    # Ensure diagonal of correlation matrix is 1
    diag(population_correlation) <- 1

    # Check eigenvalues
    check_eigenvalues <- any(matrix_eigenvalues(population_correlation) <= 0)

  }

  # Cholesky decomposition
  cholesky <- chol(population_correlation)

  # Generate data
  data <- mvtnorm::rmvnorm(sample_size, sigma = diag(total_variables))

  # Make data based on factor structure
  data <- data %*% cholesky

  # Set data dimensions
  dimensions <- dim(data)

  # Store continuous data
  continuous_data <- data

  # Check for categories greater than categorical limit and not infinite
  variable_categories <- mark_continuous_categories(variable_categories, categorical_limit)

  # Set skew/categories
  categorize_columns <- which(variable_categories <= categorical_limit)
  continuous_columns <- setdiff(seq_len(dimensions[2]), categorize_columns)
  n_continuous <- length(continuous_columns)
  n_categorize <- length(categorize_columns)

  # Store skew (bug fix for later use of `skew`)
  skew_stored <- skew

  # Initialize final skew
  final_skew <- numeric(dimensions[2])

  ## Check for categories
  if(length(categorize_columns) != 0){

    # Set skew
    if(length(skew_stored) == 1){
      skew <- rep(skew_stored, n_categorize)
    }else if(length(skew_stored) != dimensions[2]){
      skew <- shuffle_replace(skew_stored, n_categorize)
    }

    # Check for categories greater than 6
    six_categories <- variable_categories[categorize_columns] > 6
    if(any(six_categories)){

      # Check for skew not equal to zero
      if(any(skew[six_categories] != 0)){

        # Set them equal to zero (overwrite all skews)
        skew[six_categories] <- 0

        # Send warning
        warning("Variables with categories > 6 are not available to add skew. Skew for these variables were set to zero.")

      }

    }

    # Loop through columns
    for(i in seq_along(categorize_columns)){

      data[,categorize_columns[i]] <- categorize(
        data = data[,categorize_columns[i]],
        categories = variable_categories[categorize_columns[i]],
        skew_value = skew[i]
      )

    }

    # Add to final skew
    final_skew[categorize_columns] <- skew
  }

  ## Check for continuous
  if(n_continuous != 0){

    # Set skew
    if(length(skew_stored) == 1){
      skew <- rep(skew_stored, n_continuous)
    }else if(length(skew_stored) != dimensions[2]){
      skew <- shuffle_replace(skew_stored, n_continuous)
    }

    # Loop through columns
    for(i in seq_len(n_continuous)){

      data[,continuous_columns[i]] <- skew_continuous( # function in `utils-latentFactoR`
        skew = skew[i], data = data[,continuous_columns[i]]
      )

    }

    # Add to final skew
    final_skew[continuous_columns] <- skew

  }

  # Add column names to data
  colnames(data) <- paste0(
    "V", format_integer(seq_len(total_variables), digits(total_variables) - 1)
  )

  # Populate parameters
  parameters <- list(
    factors = factors,
    variables = variables,
    loadings = loading_matrix,
    cross_loadings = cross_loadings,
    factor_correlations = correlation_matrix,
    categories = variable_categories,
    categorical_limit = categorical_limit,
    skew = final_skew
  )

  # Populate results
  results <- list(
    data = data,
    population_correlation = population_correlation,
    continuous_data = continuous_data,
    parameters = parameters
  )

  # Add class
  class(results) <- "lf_simulate"

  # Return results
  return(results)

}

# Bug checking ----
# factors = 3; variables = 6;
# variables_range = NULL; loadings = 0.55;
# loadings_range = NULL; cross_loadings = 0.10;
# cross_loadings_range = NULL; correlations = 0.30;
# correlations_range = NULL; sample_size = 1000;
# variable_categories = Inf; categorical_limit = 7;
# skew = 0; skew_range = NULL;
#
# source("D:/R Packages/latentFactoR/R/simulation_helpers.R")
# source("D:/R Packages/latentFactoR/R/utils-latentFactoR.R")

# Input checking ----
#' @noRd
# Updated 12.07.2026
simulate_factors_errors <- function(
    factors, variables, variables_range = NULL,
    loadings, loadings_range = NULL,
    cross_loadings, cross_loadings_range = NULL,
    correlations, correlations_range = NULL,
    sample_size, variable_categories = Inf,
    categorical_limit = 7,
    skew = 0, skew_range = NULL
)
{

  # Check factors first
  type_error(factors, "numeric")
  length_error(factors, 1)
  range_error(factors, c(1, Inf))

  # Check sample size second
  type_error(sample_size, "numeric")
  length_error(sample_size, 1)
  range_error(sample_size, c(1, Inf))

  # Check for variables range
  if(!is.null(variables_range)){
    type_error(variables_range, "numeric") # object type error
    length_error(variables_range, 2) # object length error
    range_error(variables_range, c(3, Inf)) # object range error
    variables <- round(runif_xoshiro(
      factors,
      min = min(variables_range),
      max = max(variables_range)
    ))
  }

  # Determine total number of variables
  if(length(variables) == 1){
    total_variables <- factors * variables
  }else if(length(variables) == factors){
    total_variables <- sum(variables)
  }

  # Check for loadings range
  if(!is.null(loadings_range)){
    type_error(loadings_range, "numeric") # object type error
    length_error(loadings_range, 2)  # object length error
    range_error(loadings_range, c(-1, 1)) # object range error
    loadings <- runif_xoshiro(
      total_variables,
      min = min(loadings_range),
      max = max(loadings_range)
    )
  }

  # Check for cross-loadings range
  if(!is.null(cross_loadings_range)){
    type_error(cross_loadings_range, "numeric") # object type error
    length_error(cross_loadings_range, 2) # object length error
    range_error(cross_loadings_range, c(-1, 1)) # object range error
    cross_loadings <- runif_xoshiro(
      factors,
      min = min(cross_loadings_range),
      max = max(cross_loadings_range)
    )
  }

  # Check for correlations range
  if(!is.null(correlations_range)){
    type_error(correlations_range, "numeric") # object type error
    length_error(correlations_range, 2) # object length error
    range_error(correlations_range, c(-1, 1)) # object range error

    # Initialize correlation matrix
    correlation_matrix <- matrix(data = 0, nrow = factors, ncol = factors)

    # Population correlation matrix
    correlation_matrix[
      lower.tri(correlation_matrix)
    ] <- runif_xoshiro(
      sum(lower.tri(correlation_matrix)),
      min = min(correlations_range),
      max = max(correlations_range)
    )

    # Make correlation matrix symmetric
    correlations <- correlation_matrix + t(correlation_matrix)

  }

  # Check for skew range
  if(!is.null(skew_range)){
    type_error(skew_range, "numeric") # object type error
    length_error(skew_range, 2) # object length error
    range_error(skew_range, c(-2, 2)) # object range error
    possible_skews <- seq(-2, 2, 0.05) # possible skews
    skew_range <- round(skew_range, 2) # get to hundredths digit
    min_range <- abs(min(skew_range) - possible_skews) # difference for minimum
    min_skew <- possible_skews[which.min(min_range)] # get minimum skew
    max_range <- abs(max(skew_range) - possible_skews) # difference for maximum
    max_skew <- possible_skews[which.min(max_range)] # get maximum skew
    skew <- seq(min_skew, max_skew, 0.05) # obtain skews
  }

  # Ensure appropriate type and length for categories
  type_error(variable_categories, "numeric")
  length_error(variable_categories, c(1, total_variables))

  # Ensure appropriate types
  type_error(variables, "numeric"); type_error(variable_categories, "numeric")
  type_error(categorical_limit, "numeric"); type_error(skew, "numeric")

  # Ensure appropriate lengths
  length_error(variables, c(1, factors)); length_error(categorical_limit, 1)

  # Ensure appropriate ranges
  range_error(variables, c(2, Inf)); range_error(skew, c(-2, 2))
  range_error(variable_categories, c(2, Inf))

  # Return checked input
  return(
    list(
      variables = variables, total_variables = total_variables, loadings = loadings,
      cross_loadings = cross_loadings, correlations = correlations, skew = skew
    )
  )

}

# Shared helpers ----
#' @noRd
# Sets starting and ending variable indices for each factor ----
# Used across `simulate_factors`, `add_cross_loadings`, `add_local_dependence`,
# `add_method_factors`, `add_population_error`, and `add_wording_effects`
# Updated 12.07.2026
factor_variable_sequence <- function(variables)
{

  # Determine ending variable index for each factor
  end_variables <- cumsum(variables)

  # Determine starting variable index for each factor
  start_variables <- (end_variables + 1) - variables

  # Return start and end sequences
  return(
    list(
      start = start_variables,
      end = end_variables
    )
  )

}

#' @noRd
# Converts categories greater than the categorical limit to continuous ----
# Used across `simulate_factors`, `add_local_dependence`,
# `add_method_factors`, and `add_population_error`
# Updated 12.07.2026
mark_continuous_categories <- function(variable_categories, categorical_limit)
{

  # Check for categories greater than categorical limit and not infinite
  categories_check <- (variable_categories > categorical_limit) & (!is.infinite(variable_categories))

  # Make variables with categories greater than the limit continuous
  if(any(categories_check)){
    variable_categories[categories_check] <- Inf
  }

  # Return adjusted categories
  return(variable_categories)

}
