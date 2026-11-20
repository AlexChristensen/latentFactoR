#%%%%%%%%%%%%%%%%%%%%%%%
# GENERAL FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Determine number of categories in data ----
# Updated 02.02.2023
data_categories <- function(data)
{

  # Ensure data is matrix
  data <- as.matrix(data)

  # Loop over columns
  categories <- apply(
    data, 2, function(x){
      length(na.omit(unique(x)))
    }
  )

  # Return categories
  return(categories)

}

#' @noRd
# Convert version to number ----
# Updated 02.02.2023
version_conversion <- function(version)
{

  # Convert to character
  version <- as.character(version)

  # Remove periods
  version <- gsub("\\.", "", version)

  # Convert to numeric
  version <- as.numeric(version)

  # Return version
  return(version)

}

#' @noRd
# All-purpose symmetric checker ----
# Updated 03.02.2023
is_symmetric <- function(data){

  # Check for whether rows equal columns
  if(nrow(data) == ncol(data)){

    # Convert to matrix
    data_matrix <- as.matrix(data)

    # Remove names
    data_matrix <- unname(data_matrix)

    # Check that lower triangle equal upper triangle
    lower_triangle <- data_matrix[lower.tri(data_matrix)]
    transpose_matrix <- t(data_matrix) # ensures similar orientation
    upper_triangle <- transpose_matrix[lower.tri(transpose_matrix)]

    # Check that all are equal
    all_equal <- all(lower_triangle == upper_triangle, na.rm = TRUE)

  }else{

    # Not a matrix
    return(FALSE)

  }

  # Return whether all are equal
  return(all_equal)

}

#' @noRd
# Ensure data has dimension names ----
# Updated 03.02.2023
ensure_dimension_names <- function(data)
{

  # Check for column names
  if(is.null(colnames(data))){

    # Standardize names
    colnames(data) <- paste0(
      "V", formatC(
        x = 1:ncol(data),
        digits = (digits(ncol(data)) - 1),
        format = "d", flag = "0"
      )
    )

  }

  # Check for matrix
  if(nrow(data) == ncol(data)){

    # Check for row names
    if(is.null(data) | any(row.names(data) != colnames(data))){

      # Assign column names to row names
      row.names(data) <- colnames(data)

    }

  }

  # Return named data
  return(data)

}

#' @noRd
# No names print ----
# Updated 03.02.2023
no_name_print <- function(object){

  # Convert object to data frame
  df <- as.data.frame(object)

  # Remove column names
  colnames(df) <- NULL

  # Print with no quotes or row names
  print(df, quote = FALSE, row.names = FALSE)

}

#' @noRd
# Wrapper for `eigen` ----
# Extracts eigenvalues only
# Updated 11.03.2026
matrix_eigenvalues <- function(X)
{
  return(eigen(X, symmetric = TRUE, only.values = TRUE)$values)
}

#' @noRd
# Determine number of digits in a number ----
# Updated 24.07.2023
digits <- function(number)
{
  # Obtain the lowest value of log base 10 and add 1
  return(floor(log10(abs(number))) + 1)

}

#' @noRd
# Format number with certain decimals ----
# Mainly for naming and printing
# Updated 14.06.2024
format_decimal <- function(numbers, places)
{

  return(
    formatC(
      x = numbers,
      digits = places,
      format = "f", flag = "0"
    )
  )

}

#' @noRd
# Format number with certain integer ----
# Mainly for naming and printing
# Updated 14.06.2024
format_integer <- function(numbers, places)
{

  return(
    formatC(
      x = numbers,
      digits = places,
      format = "d", flag = "0"
    )
  )

}

#' @noRd
# Faster `ifelse` ----
# For single value replacements:
# 1.5x faster with 1 value
# 2.5x faster with 10 values
# >= 18x faster with >= 100 values
# Updated 26.02.2026
swiftelse <- function(condition, true, false)
{
  # Get condition length
  condition_length <- length(condition)

  # Check for single value
  if(condition_length == 1){

    # If TRUE
    if(condition){
      return(true)
    }else{ # Otherwise, FALSE
      return(false)
    }

  }

  # Initialize result
  result <- vector(typeof(false), condition_length)

  # Set TRUE condition
  if(length(true) == 1){
    result[condition] <- true
  }else{
    result[condition] <- true[condition]
  }

  # Set FALSE condition
  if(length(false) == 1){
    result[!condition] <- false
  }else{

    # Get opposite condition (slightly faster than repeated calls)
    opposite_condition <- !condition
    result[opposite_condition] <- false[opposite_condition]
  }


  # Return result
  return(result)

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Random number generation ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Reproducible random number generation ----
# Exported functions with a `seed` argument open a seed scope: while it is
# active, every call to the C random number generators (without its own
# seed) uses a deterministic seed derived from the scope's seed and a call
# counter. All random numbers in the package come from the C generators; R's
# random number generator is also set with `set.seed(seed)` so that
# functions from other packages used within (e.g., `fspe::fspe` in
# `estimate_dimensions`) are reproducible as well. Nested calls (e.g., `add_*`
# functions that call `simulate_factors`) share the outer scope. When the
# scope ends, R's random number generator state is restored, so a seeded call
# does not change the user's random stream
# Updated 20.11.2026
lf_rng <- new.env(parent = emptyenv())

#' @noRd
# Start a seed scope (returns a token for `lf_seed_end`; NULL when no scope
# was started)
# Updated 20.11.2026
lf_seed_start <- function(seed)
{

  # No seed or a scope is already active
  if(is.null(seed) || isTRUE(lf_rng$active)){
    return(NULL)
  }

  # Check seed
  if(
    !is.numeric(seed) || length(seed) != 1 || !is.finite(seed) ||
    seed != round(seed) || abs(seed) > .Machine$integer.max
  ){
    stop("`seed` must be a single whole number (or NULL for a random seed)", call. = FALSE)
  }

  # Store user's random number generator state
  global <- globalenv()
  token <- list(
    exists = exists(".Random.seed", envir = global, inherits = FALSE)
  )
  if(token$exists){
    token$state <- get(".Random.seed", envir = global, inherits = FALSE)
  }

  # Start scope
  lf_rng$active <- TRUE
  lf_rng$base <- abs(seed) + 1
  lf_rng$counter <- 0

  # Set R's random number generator
  set.seed(seed)

  # Return token
  return(token)

}

#' @noRd
# End a seed scope (restores the user's random number generator state)
# Updated 20.11.2026
lf_seed_end <- function(token)
{

  # Check for scope
  if(is.null(token)){
    return(invisible(NULL))
  }

  # End scope
  lf_rng$active <- FALSE

  # Restore user's random number generator state
  global <- globalenv()
  if(isTRUE(token$exists)){
    assign(".Random.seed", token$state, envir = global)
  }else if(exists(".Random.seed", envir = global, inherits = FALSE)){
    rm(".Random.seed", envir = global)
  }

  # Return invisibly
  return(invisible(NULL))

}

#' @noRd
# Seed for a call to the C random number generators
# (0 = random seed from the clock, process ID, and a counter)
# Updated 20.11.2026
lf_call_seed <- function(seed)
{

  # Supplied seed
  if(!is.null(seed)){
    return(as.numeric(seed))
  }

  # No active scope
  if(!isTRUE(lf_rng$active)){
    return(0)
  }

  # Deterministic seed from the scope's seed and the call counter
  # (non-zero and exactly representable as a double)
  lf_rng$counter <- lf_rng$counter + 1
  return(lf_rng$base * 2097152 + lf_rng$counter)

}

#' @noRd
# Generate uniform data ----
# Allows adjustment of range
# Updated 26.10.2023
runif_xoshiro <- function(n, min = 0, max = 1, seed = NULL)
{

  # Get values
  values <- .Call(
    "r_xoshiro_uniform",
    as.integer(n),
    lf_call_seed(seed),
    PACKAGE = "latentFactoR"
  )

  # Check for changes to minimum and maximum
  if(min != 0 || max != 1){ # transform
    values <- min + (max - min) * values
  }


  # Return call from C
  return(values)

}

#' @noRd
# Shuffle (without replacement) ----
# Uses xoshiro256++ random number generation: https://prng.di.unimi.it/
# Updated 30.07.2023
shuffle <- function(x, size = length(x), seed = NULL)
{

  # Return call from C
  return(
    x[.Call(
      "r_xoshiro_shuffle",
      as.integer(seq_along(x)),
      lf_call_seed(seed),
      PACKAGE = "latentFactoR"
    )][seq_len(size)]
  )

}

#' @noRd
# Weighted Shuffle (without replacement) ----
# Uses xoshiro256++ random number generation: https://prng.di.unimi.it/
# Updated 14.04.2026
weighted_shuffle <- function(x, size = length(x), prob, seed = NULL)
{

  # Return call from C
  return(
    x[.Call(
      "r_xoshiro_weighted_shuffle",
      as.integer(seq_along(x)), prob,
      lf_call_seed(seed),
      PACKAGE = "latentFactoR"
    )][seq_len(size)]
  )

}

#' @noRd
# Shuffle (with replacement) ----
# Uses xoshiro256++ random number generation: https://prng.di.unimi.it/
# Updated 10.03.2026
shuffle_replace <- function(x, size = length(x), seed = NULL)
{

  # No draws (the C routine treats a size of zero as `length(x)`)
  if(size == 0){
    return(x[0])
  }

  # Return call from C
  return(
    x[.Call(
      "r_xoshiro_shuffle_replace",
      x, size, lf_call_seed(seed),
      PACKAGE = "latentFactoR"
    )]
  )

}

#' @noRd
# Random normal generation with Ziggurat ----
# https://people.sc.fsu.edu/~jburkardt/cpp_src/ziggurat/ziggurat.html
# Updated 27.07.2023
rnorm_ziggurat <- function(n, seed = NULL)
{

  # Return call from C
  return(
    .Call(
      "r_ziggurat",
      as.integer(n),
      lf_call_seed(seed),
      PACKAGE = "latentFactoR"
    )
  )

}