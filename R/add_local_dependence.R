#' Adds Local Dependence to \code{\link[latentFactoR]{simulate_factors}} Data
#'
#' Adds local dependence to simulated data from \code{\link[latentFactoR]{simulate_factors}}.
#' See examples to get started
#'
#' @param lf_object Data object from \code{\link[latentFactoR]{simulate_factors}}
#'
#' @param method Character (length = 1).
#' Method to generate local dependence between variables.
#' Only \code{"correlate_residuals"} at the moment.
#' Future developments will include minor factor
#' and threshold-shift methods. Description of methods:
#'
#' \itemize{
#'
#' \item \code{"correlate_residuals"} --- Adds residuals directly to the population
#' correlation matrix prior to data generation (uses population correlation matrix
#' from \code{\link[latentFactoR]{simulate_factors}})
#'
#' \item \code{"minor_factors"} --- Coming soon...
#'
#' \item \code{"threshold_shifts"} --- Coming soon...
#'
#' }
#'
#' @param proportion_LD Numeric (length = 1 or \code{factors}).
#' Proportion of variables that should be locally dependent across all
#' or each factor. Accepts number of locally dependent values as well
#'
#' @param proportion_LD_range Numeric (length = 2).
#' Range of proportion of variables that are randomly selected from
#' a random uniform distribution. Accepts number of locally dependent values as well.
#' Defaults to \code{NULL}
#'
#' @param add_residuals Numeric (length = 1, \code{factors}, or total number of locally dependent variables).
#' Amount of residual to add to the population correlation matrix between two variables.
#' Only used when \code{method = "correlated_residuals"}. Magnitudes are drawn from
#' a random uniform distribution using +/- 0.05 of value input.
#' Can also be specified directly (same length as total number of locally dependent variables).
#' General effect sizes range from small (0.20), moderate (0.30), to large (0.40)
#'
#' @param add_residuals_range Numeric (length = 2).
#' Range of the residuals to add to the correlation matrix are randomly selected from
#' a random uniform distribution.
#' Defaults to \code{NULL}
#'
#' @param allow_multiple Boolean.
#' Whether a variable should be allowed to be locally dependent with
#' more than one other variable.
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} for more complex locally dependence patterns
#'
#' @return Returns a list containing:
#'
#' \item{data}{Simulated data from the specified factor model}
#'
#' \item{population_correlation}{Population correlation matrix with local dependence added}
#'
#' \item{original_correlation}{Original population correlation matrix \emph{before}
#' local dependence was added}
#'
#' \item{correlated_residuals}{A data frame with the first two columns specifying
#' the variables that are locally dependent and the third column specifying the
#' magnitude of the added residual for each locally dependent pair}
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
#' # Add local dependence
#' two_factor_LD <- add_local_dependence(
#'   lf_object = two_factor,
#'   proportion_LD = 0.25,
#'   add_residuals = 0.20,
#'   allow_multiple = FALSE
#' )
#'
#' # Randomly vary proportions
#' two_factor_LD <- add_local_dependence(
#'   lf_object = two_factor,
#'   proportion_LD_range = c(0.10, 0.50),
#'   add_residuals = 0.20,
#'   allow_multiple = FALSE
#' )
#'
#' # Randomly vary residuals
#' two_factor_LD <- add_local_dependence(
#'   lf_object = two_factor,
#'   proportion_LD = 0.25,
#'   add_residuals_range = c(0.20, 0.40),
#'   allow_multiple = FALSE
#' )
#'
#' # Randomly vary proportions, residuals, and allow multiple
#' two_factor_LD <- add_local_dependence(
#'   lf_object = two_factor,
#'   proportion_LD_range = c(0.10, 0.50),
#'   add_residuals_range = c(0.20, 0.40),
#'   allow_multiple = TRUE
#' )
#'
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Hudson Golino <hfg9s@virginia.edu>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
#'
#' @references
#' Christensen, A. P., Garrido, L. E., & Golino, H. (2023).
#' Unique variable analysis: A network psychometrics method to detect local dependence.
#' \emph{Multivariate Behavioral Research}, 1–18.
#'
#' @export
#'
# Add local dependence to simulated data
# Updated 12.07.2026
add_local_dependence <- function(
    lf_object,
    method = c(
      "correlate_residuals",
      "minor_factors",
      "threshold_shifts"
    ),
    proportion_LD, proportion_LD_range = NULL,
    add_residuals = NULL, add_residuals_range = NULL,
    allow_multiple = FALSE
)
{

  # Check for missing method
  if(missing(method)){
    method <- "correlate_residuals"
  }else{method <- tolower(match.arg(method))}

  # Check inputs
  inputs <- add_local_dependence_errors(
    lf_object, proportion_LD, proportion_LD_range
  )

  # Collect inputs
  parameters <- inputs$parameters
  proportion_LD <- inputs$proportion_LD

  # Add local dependence
  if(method == "correlate_residuals"){

    # Obtain results
    results <- correlate_residuals(
      lf_object = lf_object,
      proportion_LD = proportion_LD,
      allow_multiple = allow_multiple,
      add_residuals = add_residuals,
      add_residuals_range = add_residuals_range
    )

  }else{

    stop(
      paste0(
        "'", method, "' is not available yet and will be implemented in the future.",
        "\nFor now, use `method = \"correlate_residuals\""
      )
    )

  }

  # Add class
  class(results) <- c(class(lf_object), "lf_ld")

  # Return results
  return(results)

}

# Input checking ----
#' @noRd
# Updated 12.07.2026
add_local_dependence_errors <- function(
    lf_object, proportion_LD, proportion_LD_range = NULL
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

  # Check that population error has not yet been added
  if(is(lf_object, "lf_pe")){

    # Produce error
    stop(
      paste(
        "`lf_object` input is class \"lf_pe\" from the `add_population_error` function.",
        "Population error must be added after any local dependence to ensure proper simulation.",
        "\n\nUse `simulate_factors` to generate your data to input into this function first,",
        "then use this output in the `add_population_error` function."
      )
    )

  }

  # Obtain parameters from simulated data
  parameters <- lf_object$parameters

  # Check for proportion local dependence range
  if(!is.null(proportion_LD_range)){
    type_error(proportion_LD_range, "numeric") # object type error
    length_error(proportion_LD_range, 2) # object length error

    # Check for number of variables in range
    if(any(proportion_LD_range > 1)){

      # Target values (positions within the range vector,
      # not factor indices, so use a single representative
      # variable count rather than indexing `parameters$variables`)
      target_LD <- which(proportion_LD_range > 1)

      # Ensure proportions
      proportion_LD_range[target_LD] <-
        proportion_LD_range[target_LD] / mean(parameters$variables)

    }

    # Check for error in range
    range_error(proportion_LD_range, c(0, 1)) # object range error
    proportion_LD <- runif_xoshiro(
      parameters$factors,
      min = min(proportion_LD_range),
      max = max(proportion_LD_range)
    )
  }

  # Ensure appropriate types
  type_error(proportion_LD, "numeric")

  # Ensure appropriate lengths
  length_error(proportion_LD, c(1, parameters$factors))

  # Set proportions
  if(length(proportion_LD) == 1){
    proportion_LD <- rep(proportion_LD, parameters$factors)
  }

  # Convert local dependence proportions to proportions
  if(any(proportion_LD >= 1)){

    # Target values
    target_LD <- which(proportion_LD >= 1)

    # Ensure proportions
    proportion_LD[target_LD] <-
      proportion_LD[target_LD] / parameters$variables[target_LD]

  }

  # Ensure appropriate ranges
  range_error(proportion_LD, c(0, 1))

  # Return checked input
  return(
    list(
      parameters = parameters,
      proportion_LD = proportion_LD
    )
  )

}

#' @noRd
# Adds correlated residuals to generated data
# Updated 12.07.2026
correlate_residuals <- function(
    lf_object,
    proportion_LD, allow_multiple = FALSE,
    add_residuals, add_residuals_range
)
{

  # Obtain parameters
  parameters <- lf_object$parameters

  # Check inputs
  inputs <- correlate_residuals_errors(
    parameters, proportion_LD, allow_multiple,
    add_residuals, add_residuals_range
  )

  # Collect inputs
  variables_LD <- inputs$variables_LD
  add_residuals <- inputs$add_residuals

  # Set parameters
  factors <- parameters$factors
  variables <- parameters$variables
  loadings <- parameters$loadings
  total_variables <- sum(variables)
  sample_size <- nrow(lf_object$data)
  variable_categories <- parameters$categories
  categorical_limit <- parameters$categorical_limit
  skew <- parameters$skew
  original_correlation <- lf_object$population_correlation
  population_correlation <- original_correlation

  # Set start and end points for variables
  variable_sequence <- factor_variable_sequence(variables)
  start_variables <- variable_sequence$start
  end_variables <- variable_sequence$end

  # Set factor sequence
  factor_sequence <- seq_len(factors)

  # Initialize checks
  check_eigenvalues <- TRUE

  # Initialize iteration counter (guards against an infinite loop when the
  # supplied structure leaves no randomness to escape an inadmissible solution)
  iterations <- 0
  max_iterations <- 1000

  # Run through loop
  while(isTRUE(check_eigenvalues)){

    # Increment iteration counter and check against maximum
    iterations <- iterations + 1
    if(iterations > max_iterations){
      stop(
        paste0(
          "Could not find an admissible (positive semi-definite) solution after ",
          max_iterations, " attempts. The supplied local dependence structure may ",
          "be inadmissible. Try adjusting `proportion_LD`/`add_residuals` values"
        )
      )
    }

    # Initialize correlated residual matrix
    correlated_residuals <- matrix(
      0, nrow = 0, ncol = 2
    )

    # Loop through factors and add local dependence
    for(f in factor_sequence){

      # Determine available variables
      available_variables <- start_variables[f]:end_variables[f]

      # Check for cross-loading class
      if(is(lf_object, "lf_cl")){

        # Obtain loadings for available variables (ensure matrix)
        available_loadings <- matrix(
          loadings[available_variables, -f],
          ncol = factors - 1
        )

        # Collect sum of cross-loadings
        sum_cross_loadings <- rowSums(available_loadings)

        # Remove variables with non-zero cross-loadings
        available_variables <- available_variables[
          sum_cross_loadings == 0
        ]

        # Cross-loadings may have depleted the pool of available
        # variables for this factor; each local dependence pair
        # needs two distinct variables, so cap the number of pairs
        # to what the (possibly reduced) pool can actually supply
        variables_LD[f] <- min(variables_LD[f], floor(length(available_variables) / 2))

      }

      # Item rows
      item_rows <- if(isTRUE(allow_multiple)){
        shuffle_replace(available_variables, variables_LD[f])
      }else{
        shuffle(available_variables, variables_LD[f])
      }

      # Set remaining variables
      if(isTRUE(allow_multiple)){

        # Do not remove variables
        remaining_variables <- available_variables

      }else{

        # Remove already included variables
        remaining_variables <- setdiff(
          available_variables, item_rows
        )

      }

      # Item columns
      item_columns <- if(isTRUE(allow_multiple)){
        shuffle_replace(remaining_variables, variables_LD[f])
      }else{
        shuffle(remaining_variables, variables_LD[f])
      }

      # Bind to correlated residual matrix
      correlated_residuals <- rbind(
        correlated_residuals,
        cbind(item_rows, item_columns)
      )

      # Obtain same variable rows
      same_variable_rows <- which(
        correlated_residuals[,"item_rows"] ==
          correlated_residuals[,"item_columns"]
      )

      # Replace until there are no duplicate rows
      while(any(same_variable_rows)){

        # Replace second column with new variable
        correlated_residuals[same_variable_rows, 2] <- if(isTRUE(allow_multiple)){
          shuffle_replace(remaining_variables, length(same_variable_rows))
        }else{
          shuffle(remaining_variables, length(same_variable_rows))
        }

        # Re-check for duplicate rows
        same_variable_rows <- which(
          correlated_residuals[,"item_rows"] ==
            correlated_residuals[,"item_columns"]
        )

      }

      # Obtain duplicate rows
      duplicate_rows <- match_row(correlated_residuals)

      # Replace until there are no duplicate rows
      while(any(duplicate_rows)){

        # Replace second column with new variable
        correlated_residuals[duplicate_rows, 2] <- if(isTRUE(allow_multiple)){
          shuffle_replace(remaining_variables, length(duplicate_rows))
        }else{
          shuffle(remaining_variables, length(duplicate_rows))
        }

        # Re-check for duplicate rows
        duplicate_rows <- match_row(correlated_residuals)

      }

    }

    # Add column for residual
    correlated_residuals <- cbind(
      correlated_residuals, rep(0, nrow(correlated_residuals))
    )

    # Add column name for residual
    colnames(correlated_residuals)[3] <- "residual_correlation"

    # Loop through correlated residuals
    if(nrow(correlated_residuals) != 0){

      # Check if add residuals length equals 1 or number of factors
      if(length(add_residuals) != sum(variables_LD)){

        # If only one value
        if(length(add_residuals) == 1){
          add_residuals <- rep(add_residuals, sum(variables_LD))
        }else if(length(add_residuals) == factors){# Length of factors

          # Loop through number of local dependence variables
          add_residuals <- unlist(lapply(factor_sequence, function(i){
            rep(add_residuals[i], variables_LD[i])
          }))

        }

      }

      # Loop through correlated residuals
      for(i in seq_len(nrow(correlated_residuals))){

        # Insert or compute random residual
        if(!is.null(add_residuals_range)){
          random_residual <- add_residuals[i]
        }else{
          random_residual <- runif_xoshiro(
            1,
            min = add_residuals[i] - 0.05,
            max = add_residuals[i] + 0.05
          )
        }

        # Obtain sign
        original_sign <- sign(original_correlation[
          correlated_residuals[i,1],
          correlated_residuals[i,2]
        ])

        # Add residuals to correlation matrix
        population_correlation[
          correlated_residuals[i,1],
          correlated_residuals[i,2]
        ] <- (abs(original_correlation[
          correlated_residuals[i,1],
          correlated_residuals[i,2]
        ]) + random_residual) * original_sign

        # Ensure symmetric
        population_correlation[
          correlated_residuals[i,2],
          correlated_residuals[i,1]
        ] <- population_correlation[
          correlated_residuals[i,1],
          correlated_residuals[i,2]
        ]

        # Add random residual to output
        correlated_residuals[i,"residual_correlation"] <-
          random_residual * original_sign

      }

    }

    # Check eigenvalues
    check_eigenvalues <- any(matrix_eigenvalues(population_correlation) <= 0)

    # Return population correlation to original state (if necessary)
    if(isTRUE(check_eigenvalues)){
      population_correlation <- original_correlation
    }

  }

  # Cholesky decomposition
  cholesky <- chol(population_correlation)

  # Generate data
  data <- mvtnorm::rmvnorm(sample_size, sigma = diag(total_variables))

  # Make data based on factor structure
  data <- data %*% cholesky

  # Set data dimensions
  dimensions <- dim(data)

  # Ensure appropriate type and length for categories
  type_error(variable_categories, "numeric")
  length_error(variable_categories, c(1, total_variables))

  # Identify categories to variables
  if(length(variable_categories) == 1){
    variable_categories <- rep(variable_categories, total_variables)
  }

  # Check for categories greater than categorical limit and not infinite
  variable_categories <- mark_continuous_categories(variable_categories, categorical_limit)

  # Set skew/categories
  ## Target columns to categorize and/or add skew
  categorize_columns <- which(variable_categories <= categorical_limit)
  continuous_columns <- setdiff(seq_len(dimensions[2]), categorize_columns)

  # Ensure skew is in the appropriate direction for correlated residuals
  if(nrow(correlated_residuals) != 0){

    # 1. Obtain loading signs
    signs <- numeric(nrow(loadings))

    # Ensure proper signs for skew
    for(i in factor_sequence){

      # Target dominant loadings
      target_loadings <- start_variables[i]:end_variables[i]

      # Determine sign
      signs[target_loadings] <- sign(loadings[target_loadings, i])

    }

    # 2. Obtain correlated residual chains

    # Obtain residual variables
    residual_variables <- correlated_residuals[,c(
      "item_rows", "item_columns"
    )]

    # Ensure matrix
    if(!is.matrix(residual_variables)){
      residual_variables <- t(as.matrix(residual_variables))
    }

    # Initialize residual chain list
    residual_chain <- vector("list", nrow(residual_variables))

    # Create residual chain
    while(nrow(residual_variables) != 0){

      # Determine if either variable exists in residual chain
      combine_residuals <- unlist(
        lapply(residual_chain, function(x){

          if(any(residual_variables[1,] %in% x)){
            return(TRUE)
          }else{
            return(FALSE)
          }

        })
      )

      # Check for whether to combine
      if(!any(combine_residuals)){

        # Find NULL lists
        null_list <- unlist(lapply(residual_chain, is.null))

        # First NULL list
        residual_chain[[which(null_list)[1]]] <-
          unique(as.vector(residual_variables[1,]))

      }else if(sum(combine_residuals) == 1){

        # Find residual chain to add to
        add_to_chain <- which(combine_residuals)

        # Add variables to chain
        residual_chain[[add_to_chain]] <- unique( # keep unique
          c(
            residual_chain[[add_to_chain]], # existing chain
            as.vector(residual_variables[1,]) # new to add
          )
        )

      }else{

        # Find residual chains to merge
        merge_chains <- which(combine_residuals)

        # Merge into first merge
        residual_chain[[merge_chains[1]]] <- unique(unlist(residual_chain[merge_chains]))

        # Add to first merge chain
        residual_chain[[merge_chains[1]]] <- unique( # keep unique
          c(
            residual_chain[[merge_chains[1]]], # existing chain
            as.vector(residual_variables[1,]) # new to add
          )
        )

        # Set second merge chain to NULL
        residual_chain[[merge_chains[2]]] <- NULL

      }

      # Remove residual variables from consideration
      residual_variables <- matrix(residual_variables[-1,], ncol = 2)

    }

    # Find NULL lists
    null_list <- unlist(lapply(residual_chain, is.null))

    # Keep non-NULL chains
    residual_chain <- residual_chain[!null_list]

    # Ensure skew in the same direction
    for(i in seq_along(residual_chain)){

      # Handle skew
      skew[residual_chain[[i]]] <- handle_skew_signs(
        skews = skew[residual_chain[[i]]],
        signs = signs[residual_chain[[i]]]
      )

    }

  }

  ## Check for categories
  if(length(categorize_columns) != 0){

    # Set skew
    if(length(skew) == 1){
      skew <- rep(skew, length(categorize_columns))
    }else if(length(skew) != dimensions[2]){
      skew <- shuffle_replace(skew, length(categorize_columns))
    }

    # Loop through columns
    for(i in categorize_columns){

      data[,i] <- categorize(
        data = data[,i],
        categories = variable_categories[i],
        skew_value = skew[i]
      )

    }

  }

  ## Check for continuous
  if(length(continuous_columns) != 0){

    # Set skew
    if(length(skew) == 1){
      skew <- rep(skew, length(continuous_columns))
    }else if(length(skew) != dimensions[2]){
      skew <- shuffle_replace(skew, length(continuous_columns))
    }

    # Loop through columns
    for(i in continuous_columns){

      data[,i] <- skew_continuous( # function in `utils-latentFactoR`
        skew = skew[i], data = data[,i]
      )

    }

  }

  # Add column names to data
  colnames(data) <- paste0(
    "V", format_integer(seq_len(total_variables), digits(total_variables) - 1)
  )

  # Update skews
  parameters$skew <- skew

  # Populate results
  results <- list(
    data = data,
    population_correlation = population_correlation,
    parameters = parameters,
    correlated_residuals = as.data.frame(correlated_residuals),
    original_results = lf_object
  )

  # Return results
  return(results)

}

# Input checking ----
#' @noRd
# Updated 12.07.2026
correlate_residuals_errors <- function(
    parameters, proportion_LD, allow_multiple,
    add_residuals, add_residuals_range
)
{

  # Obtain number of local dependencies
  variables_LD <- round(proportion_LD * parameters$variables)

  # If variables cannot have multiple local dependencies,
  # then number of local dependencies needs to be cut in half (per factor)
  if(!isTRUE(allow_multiple)){
    variables_LD <- ifelse(
      variables_LD == 1, variables_LD,
      floor(variables_LD / 2)
    )
  }

  # Check for add residual range
  if(!is.null(add_residuals_range)){
    type_error(add_residuals_range, "numeric") # object type error
    length_error(add_residuals_range, 2) # object length error
    range_error(add_residuals_range, c(0, 1)) # object range error
    add_residuals <- runif_xoshiro(
      sum(variables_LD),
      min = min(add_residuals_range),
      max = max(add_residuals_range)
    )
  }

  # Ensure appropriate types
  type_error(add_residuals, "numeric")

  # Ensure appropriate lengths
  length_error(add_residuals, c(1, parameters$factors, sum(variables_LD)))

  # Ensure appropriate ranges
  range_error(add_residuals, c(0, 1))

  # Return checked input
  return(
    list(
      variables_LD = variables_LD,
      add_residuals = add_residuals
    )
  )

}
