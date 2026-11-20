#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Shared Simulation Helpers ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#' @noRd
# Variable names ("V1", ..., or zero-padded "V01", ...)
# Used across the simulation and estimation functions
# Updated 20.11.2026
variable_names <- function(total_variables)
{
  return(
    paste0(
      "V", format_integer(seq_len(total_variables), digits(total_variables) - 1)
    )
  )
}

#' @noRd
# Generate data from a population correlation matrix ----
# Draws multivariate normal data, categorizes variables with categories
# at or below the categorical limit (with skew), and skews continuous variables.
# Skew can be a single value, one value per variable, or a set of
# values to sample from (e.g., from `skew_range`).
# Used across `simulate_factors`, `simulate_hierarchical_factors`,
# `add_local_dependence`, `add_method_factors`, and `add_population_error`
# Returns the data, the skew applied to each variable, and the
# (continuous-marked) categories
# Updated 20.11.2026
generate_lf_data <- function(
    R, sample_size, variable_categories,
    categorical_limit, skew
)
{

  # Obtain total variables
  total_variables <- ncol(R)

  # Identify categories to variables
  if(length(variable_categories) == 1){
    variable_categories <- rep(variable_categories, total_variables)
  }

  # Cholesky decomposition
  cholesky <- chol(R)

  # Generate data (standard normal draws from the C random number generator)
  data <- matrix(
    rnorm_ziggurat(sample_size * total_variables),
    nrow = sample_size, ncol = total_variables
  )

  # Make data based on factor structure
  data <- data %*% cholesky

  # Check for categories greater than categorical limit and not infinite
  variable_categories <- mark_continuous_categories(variable_categories, categorical_limit)

  # Set skew/categories
  categorize_columns <- which(variable_categories <= categorical_limit)
  continuous_columns <- setdiff(seq_len(total_variables), categorize_columns)

  # Initialize final skew
  final_skew <- numeric(total_variables)

  # Select skew for a set of columns
  select_skew <- function(columns){

    if(length(skew) == 1){
      rep(skew, length(columns))
    }else if(length(skew) == total_variables){
      # Full-length (per-variable) skew: select this subset's own values
      skew[columns]
    }else{
      shuffle_replace(skew, length(columns))
    }

  }

  ## Check for categories
  if(length(categorize_columns) != 0){

    # Set skew
    categorical_skew <- select_skew(categorize_columns)

    # Check for categories greater than 6
    six_categories <- variable_categories[categorize_columns] > 6
    if(any(six_categories) && any(categorical_skew[six_categories] != 0)){

      # Set them equal to zero (overwrite all skews)
      categorical_skew[six_categories] <- 0

      # Send warning
      warning("Variables with categories > 6 are not available to add skew. Skew for these variables were set to zero.")

    }

    # Loop through columns
    for(i in seq_along(categorize_columns)){

      data[,categorize_columns[i]] <- categorize(
        data = data[,categorize_columns[i]],
        categories = variable_categories[categorize_columns[i]],
        skew_value = categorical_skew[i]
      )

    }

    # Add to final skew
    final_skew[categorize_columns] <- categorical_skew

  }

  ## Check for continuous
  if(length(continuous_columns) != 0){

    # Set skew
    continuous_skew <- select_skew(continuous_columns)

    # Loop through columns
    for(i in seq_along(continuous_columns)){

      data[,continuous_columns[i]] <- skew_continuous( # function in `utils-latentFactoR`
        skew = continuous_skew[i], data = data[,continuous_columns[i]]
      )

    }

    # Add to final skew
    final_skew[continuous_columns] <- continuous_skew

  }

  # Add column names to data
  colnames(data) <- variable_names(total_variables)

  # Return results
  return(
    list(
      data = data, skew = final_skew,
      categories = variable_categories
    )
  )

}

#' @noRd
# Check {latentFactoR} simulated objects ----
# Ensures `lf_object` is from `simulate_factors` (class "lf_simulate"),
# optionally forbids classes (e.g., population error must be added last),
# and optionally requires a flat (non-hierarchical) factor structure
# Used across the `add_*` functions and `data_to_zipfs`
# Updated 20.11.2026
check_lf_object <- function(
    lf_object, function_name, forbid = NULL, flat = FALSE
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
      ), call. = FALSE
    )

  }

  # Check for forbidden classes
  forbidden <- forbid[vapply(forbid, function(x){is(lf_object, x)}, logical(1))]
  if(length(forbidden) > 0){

    # Function that produces each class
    producer <- c(
      lf_cl = "add_cross_loadings", lf_ld = "add_local_dependence",
      lf_me = "add_method_factors", lf_we = "add_wording_effects",
      lf_pe = "add_population_error", lf_zipfs = "data_to_zipfs"
    )

    # Produce error
    stop(
      paste0(
        "`lf_object` input is class \"", forbidden[1], "\" from the `",
        producer[forbidden[1]], "` function, which cannot be used as input into `",
        function_name, "`.",
        "\n\nApply `", function_name, "` before `", producer[forbidden[1]],
        "` (or use the `original_results` element of the object)"
      ), call. = FALSE
    )

  }

  # Check for flat factor structure (e.g., not hierarchical)
  if(
    isTRUE(flat) && (
      is.null(lf_object$parameters$loadings) ||
      !is.matrix(lf_object$parameters$loadings)
    )
  ){

    # Produce error
    stop(
      paste0(
        "`lf_object` does not have a flat factor structure ",
        "(e.g., output from `simulate_hierarchical_factors`), which is not ",
        "supported by `", function_name, "`.",
        "\n\nUse output from `simulate_factors`"
      ), call. = FALSE
    )

  }

}

#' @noRd
# Symmetric matrix with off-diagonal values drawn uniformly from a range ----
# (diagonal is zero; set by the caller)
# Updated 20.11.2026
draw_correlation_matrix <- function(dimension, range)
{

  # Initialize matrix
  correlation_matrix <- matrix(data = 0, nrow = dimension, ncol = dimension)

  # Draw lower triangle
  correlation_matrix[lower.tri(correlation_matrix)] <- runif_xoshiro(
    sum(lower.tri(correlation_matrix)),
    min = min(range), max = max(range)
  )

  # Return symmetric matrix
  return(correlation_matrix + t(correlation_matrix))

}

#' @noRd
# Snap a skew range onto the available skews (-10 to 10 in 0.05) ----
# Updated 20.11.2026
snap_skew_range <- function(skew_range)
{

  # Possible skews
  possible_skews <- seq(-10, 10, 0.05)

  # Get to hundredths digit
  skew_range <- round(skew_range, 2)

  # Minimum and maximum skews
  min_skew <- possible_skews[which.min(abs(min(skew_range) - possible_skews))]
  max_skew <- possible_skews[which.min(abs(max(skew_range) - possible_skews))]

  # Return skews
  return(seq(min_skew, max_skew, 0.05))

}

#' @noRd
# Signs of the dominant loadings of each variable ----
# Updated 20.11.2026
dominant_signs <- function(loadings, variables)
{

  # Set sequence of variables for each factor
  variable_sequence <- factor_variable_sequence(variables)

  # Initialize signs
  signs <- numeric(nrow(loadings))

  # Determine sign of the dominant loadings
  for(i in seq_along(variables)){
    target <- variable_sequence$start[i]:variable_sequence$end[i]
    signs[target] <- sign(loadings[target, i])
  }

  # Return signs
  return(signs)

}

#' @noRd
# Flip the first dominant loadings on each factor ----
# (the first `round(proportion * variables)` dominant loadings
# on each factor are made negative)
# Updated 20.11.2026
flip_dominant_loadings <- function(loadings, variables, proportion)
{

  # Set sequence of variables for each factor
  variable_sequence <- factor_variable_sequence(variables)

  # Loop through factors
  for(i in seq_along(variables)){

    # Obtain number of flipped variables
    negative_variables <- round(proportion[i] * variables[i])

    # Check for zero negative variables
    if(negative_variables != 0){

      # Target dominant loadings
      target <- variable_sequence$start[i]:variable_sequence$end[i]

      # Make loadings absolute
      loadings[target, i] <- abs(loadings[target, i])

      # Set dominant loadings to inverse
      flipped <- target[seq_len(negative_variables)]
      loadings[flipped, i] <- -loadings[flipped, i]

    }

  }

  # Return loadings
  return(loadings)

}

#' @noRd
# Resolve proportions of variables per factor ----
# Accepts a proportion (0 to 1) or a number of variables (> 1) per factor,
# either as values (length 1 or `factors`) or a range (length 2) from which
# a proportion is drawn for each factor. Values greater than one are
# treated as numbers of variables; a value of exactly one is a number of
# variables when `one_is_count = TRUE` (`add_cross_loadings`,
# `add_local_dependence`) and a proportion (all variables) otherwise
# Updated 20.11.2026
resolve_proportion <- function(
    proportion, proportion_range, variables, name, one_is_count = FALSE
)
{

  # Identify numbers of variables
  is_count <- function(x){
    if(isTRUE(one_is_count)){x >= 1}else{x > 1}
  }

  # Number of factors
  factors <- length(variables)

  # Check for range
  if(!is.null(proportion_range)){

    # Check range
    type_error(proportion_range, "numeric")
    length_error(proportion_range, 2)

    # Convert numbers of variables to proportions
    # (positions within the range vector, not factor indices, so use a
    # single representative variable count)
    counts <- is_count(proportion_range)
    proportion_range[counts] <- proportion_range[counts] / mean(variables)

    # Check for error in range
    range_error(proportion_range, c(0, 1))

    # Draw proportions
    proportion <- runif_xoshiro(
      factors, min = min(proportion_range), max = max(proportion_range)
    )

  }

  # Check values
  if(!is.numeric(proportion)){
    stop(
      paste0(
        "Input into '", name, "' argument is not 'numeric'. Input is ",
        paste0("'", class(proportion), "'", collapse = ", ")
      ), call. = FALSE
    )
  }
  if(!length(proportion) %in% c(1, factors)){
    stop(
      paste0(
        "Length of '", name, "' (", length(proportion), ") does not match ",
        "expected length(s). Length must be: '1' or '", factors, "'"
      ), call. = FALSE
    )
  }

  # Set proportions
  if(length(proportion) == 1){
    proportion <- rep(proportion, factors)
  }

  # Convert numbers of variables to proportions
  counts <- is_count(proportion)
  proportion[counts] <- proportion[counts] / variables[counts]

  # Check range
  if(any(proportion < 0 | proportion > 1)){
    stop(
      paste0(
        "'", name, "' must be between 0 and 1 (proportion) or between 1 and ",
        "the number of variables on each factor (number of variables)"
      ), call. = FALSE
    )
  }

  # Return proportions
  return(proportion)

}
