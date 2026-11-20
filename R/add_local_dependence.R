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
#' @param seed Numeric (length = 1).
#' Seed for reproducible results. The same \code{seed} gives the same results
#' (and does not change the random number stream of the R session).
#' Defaults to \code{NULL} (random seed)
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
    allow_multiple = FALSE,
    seed = NULL
)
{

  # Set seed (reproducible results)
  seed_token <- lf_seed_start(seed)
  on.exit(lf_seed_end(seed_token), add = TRUE)


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
  # (population error must be added after any local dependence)
  check_lf_object(
    lf_object, "add_local_dependence",
    forbid = c("lf_we", "lf_pe", "lf_zipfs"), flat = TRUE
  )

  # Obtain parameters from simulated data
  parameters <- lf_object$parameters

  # Resolve proportions (or numbers of variables)
  proportion_LD <- resolve_proportion(
    proportion_LD, proportion_LD_range, parameters$variables,
    "proportion_LD", one_is_count = TRUE
  )

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

      # Replace pairs that are the same variable or duplicate another
      # pair (in either order) until all pairs are unique
      bad_rows <- invalid_residual_pairs(correlated_residuals)
      replacement_count <- 0
      while(length(bad_rows) != 0){

        # Guard against impossible configurations
        replacement_count <- replacement_count + 1
        if(replacement_count > 1000){
          stop(
            paste(
              "Could not find unique pairs of locally dependent variables.",
              "Try decreasing `proportion_LD`"
            ), call. = FALSE
          )
        }

        # Replace second column with new variable
        correlated_residuals[bad_rows, 2] <- if(isTRUE(allow_multiple)){
          shuffle_replace(remaining_variables, length(bad_rows))
        }else{
          shuffle(remaining_variables, length(bad_rows))
        }

        # Re-check pairs
        bad_rows <- invalid_residual_pairs(correlated_residuals)

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

      # Residual for each pair: one value, one value per factor, or one
      # value per pair (kept separate from `add_residuals` so that retries
      # and capped numbers of pairs are expanded correctly)
      residual_values <- if(length(add_residuals) == nrow(correlated_residuals)){
        add_residuals
      }else if(length(add_residuals) == 1){
        rep(add_residuals, nrow(correlated_residuals))
      }else if(length(add_residuals) == factors){
        rep(add_residuals, times = variables_LD)
      }else{
        add_residuals[seq_len(nrow(correlated_residuals))]
      }

      # Loop through correlated residuals
      for(i in seq_len(nrow(correlated_residuals))){

        # Insert or compute random residual
        if(!is.null(add_residuals_range)){
          random_residual <- residual_values[i]
        }else{
          random_residual <- runif_xoshiro(
            1,
            min = residual_values[i] - 0.05,
            max = residual_values[i] + 0.05
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

  # Per-variable skew
  if(length(skew) == 1){
    skew <- rep(skew, total_variables)
  }

  # Ensure skew is in the appropriate direction for correlated residuals
  # (variables in a chain of correlated residuals are skewed in the same
  # direction as their dominant loadings)
  if(nrow(correlated_residuals) != 0 && length(skew) == total_variables){

    # Obtain loading signs
    signs <- dominant_signs(loadings, variables)

    # Obtain correlated residual chains (connected variables)
    residual_chain <- residual_chains(
      correlated_residuals[, c("item_rows", "item_columns"), drop = FALSE]
    )

    # Ensure skew in the same direction
    for(i in seq_along(residual_chain)){

      # Handle skew
      skew[residual_chain[[i]]] <- handle_skew_signs(
        skews = skew[residual_chain[[i]]],
        signs = signs[residual_chain[[i]]]
      )

    }

  }

  # Generate data (categorize and skew variables)
  generated <- generate_lf_data(
    population_correlation, sample_size, variable_categories,
    categorical_limit, skew
  )
  data <- generated$data

  # Update skews
  parameters$skew <- generated$skew

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

#' @noRd
# Rows of correlated residual pairs that are the same variable or
# that duplicate another pair (in either order)
# Updated 20.11.2026
invalid_residual_pairs <- function(pairs)
{

  # Same variable
  same <- pairs[, 1] == pairs[, 2]

  # Duplicated pairs (unordered)
  keys <- paste(
    pmin(pairs[, 1], pairs[, 2]), pmax(pairs[, 1], pairs[, 2])
  )
  duplicate <- duplicated(keys)

  # Return rows
  return(which(same | duplicate))

}

#' @noRd
# Chains of variables connected by correlated residuals
# (connected components of the residual pairs)
# Updated 20.11.2026
residual_chains <- function(pairs)
{

  # Initialize chains
  chains <- list()

  # Loop through pairs
  for(k in seq_len(nrow(pairs))){

    # Variables in pair
    members <- unique(as.vector(pairs[k, ]))

    # Chains that share a variable with the pair
    connected <- which(
      vapply(chains, function(x){any(members %in% x)}, logical(1))
    )

    # Merge all connected chains with the pair
    if(length(connected) == 0){
      chains[[length(chains) + 1]] <- members
    }else{
      merged <- unique(c(unlist(chains[connected]), members))
      chains <- chains[-connected]
      chains[[length(chains) + 1]] <- merged
    }

  }

  # Return chains
  return(chains)

}
