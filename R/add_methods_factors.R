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
#' @return Returns a list containing:
#'
#' \item{data}{Biased data simulated data from the specified factor model}
#'
#' \item{unbiased_data}{The corresponding unbiased data prior to replacing values
#' to generate the (biased) \code{data}}
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
#' @export
#'
# Add methods factors to simulated data ----
# Updated 17.04.2024
add_method_factors <- function(
    lf_object,
    proportion_negative = 0.50,
    proportion_negative_range = NULL,
    methods_factors,
    methods_loadings,
    methods_loadings_range = 0.00,
    methods_correlations,
    methods_correlations_range = NULL
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

  # Check for missing methods factors
  if(missing(methods_factors)){methods_factors <- parameters$factors}
  if(missing(methods_correlations)){
    methods_correlations <- parameters$factor_correlations
    methods_correlations[] <- 1
   }

  # Check for methods correlations range
  if(!is.null(methods_correlations_range)){
    type_error(methods_correlations_range, "numeric") # object type error
    length_error(methods_correlations_range, 2) # object length error
    range_error(methods_correlations_range, c(-1, 1)) # object range error

    # Initialize correlation matrix
    methods_correlation_matrix <- matrix(
      data = 0, nrow = parameters$factors, ncol = parameters$factors
    )

    # Population correlation matrix
    methods_correlation_matrix[
      lower.tri(methods_correlation_matrix)
    ] <- runif(
      sum(lower.tri(methods_correlation_matrix)),
      min = min(methods_correlations_range),
      max = max(methods_correlations_range)
    )

    # Make correlation matrix symmetric
    methods_correlations <- methods_correlation_matrix + t(methods_correlation_matrix)

  }

  # Ensure appropriate types
  type_error(proportion_negative, "numeric")

  # Ensure appropriate lengths
  length_error(proportion_negative, c(1, parameters$factors))

  # Set proportions
  if(length(proportion_negative) == 1){
    proportion_negative <- rep(proportion_negative, parameters$factors)
  }

  # Convert negative wording proportions to proportions
  if(any(proportion_negative > 1)){

    # Target values
    target_negative <- which(proportion_negative > 1)

    # Ensure proportions
    proportion_negative[target_negative] <-
      proportion_negative[target_negative] / parameters$variables[target_negative]

  }

  # Ensure appropriate ranges
  range_error(proportion_negative, c(0, 1));

  # Obtain substantive loadings
  loadings <- parameters$loadings

  # Obtain variables
  variables <- parameters$variables

  # Set sequence of variables for each factor
  end_variables <- cumsum(parameters$variables)
  start_variables <- (end_variables + 1) - parameters$variables

  # Flip dominant loadings
  for(i in 1:ncol(loadings)){

    # Obtain number of flipped variables
    negative_variables <- round(proportion_negative[i] * variables[i])

    # Check for zero negative variables
    if(negative_variables != 0){

      # Target dominant loadings
      target_loadings <- start_variables[i]:end_variables[i]

      # Make loadings absolute
      loadings[target_loadings, i] <- abs(loadings[target_loadings, i])

      # Set dominant loadings to inverse
      loadings[
        target_loadings[1:negative_variables],
        i
      ] <- -loadings[
        target_loadings[1:negative_variables],
        i
      ]

    }

  }

  # Initialize signs
  signs <- numeric(nrow(loadings))

  # Ensure proper signs for skew
  for(i in 1:ncol(loadings)){

    # Target dominant loadings
    target_loadings <- start_variables[i]:end_variables[i]

    # Determine sign
    signs[target_loadings] <- sign(loadings[target_loadings, i])

  }

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
      methods_factors <- sample(parameters$factors, size = methods_factors)
    }

  }

  # Check for error in length and range
  length_error(methods_factors, seq_len(parameters$factors)) # object length error
  range_error(methods_factors, c(0, parameters$factors)) # object range error

  # Initialize checks
  check_eigenvalues <- TRUE
  check_communalities <- TRUE

  # Run through loop
  while(isTRUE(check_eigenvalues) | isTRUE(check_communalities)){

    # Populate methods loadings matrix
    if(!is(methods_loadings, "matrix")){

      # Get dimensions of original loadings matrix
      loading_dimensions <- dim(loadings)

      # Create methods loadings matrix
      methods_loadings_matrix <- matrix(
        0, nrow = loading_dimensions[1], ncol = loading_dimensions[2]
      )

      # Create starting and ending of variable sequences
      end_variables <- cumsum(variables)
      start_variables <- (end_variables + 1) - variables

      # Identify dominant loadings
      if(length(methods_loadings) == 1){methods_loadings <- rep(methods_loadings, parameters$factors)}
      if(length(methods_loadings_range) == 1){methods_loadings_range <- rep(methods_loadings_range, parameters$factors)}

      # Populate factor loadings
      for(i in methods_factors){

        # Generate loadings from uniform distribution
        methods_loadings_matrix[
          start_variables[i]:end_variables[i], i # dominant loadings
        ] <- runif(
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
    check_eigenvalues <- any(eigen(population_correlation)$values <= 0)

  }

  # Cholesky decomposition
  cholesky <- chol(population_correlation)

  # Total variables
  total_variables <- nrow(combined_loadings)

  # Generate data
  data <- mvtnorm::rmvnorm(nrow(lf_object$data), sigma = diag(total_variables))

  # Make data based on factor structure
  data <- data %*% cholesky

  # Store continuous data
  continuous_data <- data

  # Variable categories
  variable_categories <- parameters$categories

  # Ensure appropriate type and length for categories
  type_error(variable_categories, "numeric")
  length_error(variable_categories, c(1, total_variables))

  # Identify categories to variables
  if(length(variable_categories) == 1){
    variable_categories <- rep(variable_categories, total_variables)
  }

  # Check for categories greater than categorical limit and not infinite
  if(any(variable_categories > parameters$categorical_limit & !is.infinite(variable_categories))){

    # Make variables with categories greater than 7 (or categorical_limit) continuous
    variable_categories[
      variable_categories > parameters$categorical_limit & !is.infinite(variable_categories)
    ] <- Inf

  }

  # Set skew/categories
  ## Target columns to categorize and/or add skew
  categorize_columns <- which(variable_categories <= parameters$categorical_limit)
  continuous_columns <- setdiff(1:ncol(data), categorize_columns)

  # Store skew (bug fix for later use of `skew`)
  skew_stored <- parameters$skew

  # Initialize final skew
  final_skew <- numeric(ncol(lf_object$data))

  ## Check for categories
  if(length(categorize_columns) != 0){

    # Set skew
    if(length(skew_stored) == 1){
      skew <- rep(skew_stored, length(categorize_columns))
    }else if(length(skew_stored) != ncol(data)){
      skew <- sample(skew_stored, length(categorize_columns), replace = TRUE)
    }else{
      skew <- skew_stored
    }

    # Check for categories greater than 6
    if(any(variable_categories[categorize_columns] > 6)){

      # Obtain indices
      categories_greater <- which(variable_categories[categorize_columns] > 6)

      # Check for skew not equal to zero
      if(any(skew[categories_greater] != 0)){

        # Set them equal to zero (overwrite all skews)
        skew[categories_greater] <- 0

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
  if(length(continuous_columns) != 0){

    # Set skew
    if(length(skew_stored) == 1){
      skew <- rep(skew_stored, length(continuous_columns))
    }else if(length(skew_stored) != ncol(lf_object$data)){
      skew <- sample(skew_stored, length(continuous_columns), replace = TRUE)
    }else{
      skew <- skew_stored
    }

    # Loop through columns
    for(i in seq_along(continuous_columns)){

      data[,continuous_columns[i]] <- skew_continuous( # function in `utils-latentFactoR`
        skewness = skew[i],
        data = data[,continuous_columns[i]]
      )

    }

    # Add to final skew
    final_skew[continuous_columns] <- skew

  }

  # Add column names to data
  colnames(data) <- paste0(
    "V", formatC(
      x = 1:total_variables,
      digits = floor(log10(total_variables)),
      flag = "0", format = "d"
    )
  )

  # Update parameters
  parameters$loadings <- combined_loadings
  parameters$cross_loadings <- c(
    parameters$cross_loadings, rep(0, length(parameters$cross_loadings))
  )
  parameters$factor_correlations <- combined_correlations

  # Populate results
  results <- list(
    data = data,
    unbiased_data = lf_object$data,
    parameters = parameters,
    original_results = lf_object
  )

  # Add class
  class(results) <- c(class(lf_object), "lf_me")

  # Return results
  return(results)

}


