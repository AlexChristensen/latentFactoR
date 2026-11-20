#%%%%%%%%%%%%%%%%
# categorize ----
#%%%%%%%%%%%%%%%%

#' @noRd
# Ensures skew is in an appropriate direction based on loadings
# Updated 12.12.2022
handle_skew_signs <- function(skews, signs){

  # Check if all skew signs are in the same direction
  if(all(sign(skews) == -1)){ # All in negative direction

    # Make all skew for negative variables positive
    skews[signs == -1] <- abs(skews[signs == -1])

  }else if(all(sign(skews) == 1)){ # All in negative direction

    # Make all skew for positive variables negative
    skews[signs == 1] <- -abs(skews[signs == 1])

  }else{

    # If mixed, base signs on mode of positive variables
    # with preference for negative skew (i.e., ties go to negative skew)

    # Obtain skew signs for positive variables
    positive_skew_signs <- sign(skews[signs == 1])

    # Determine whether mode with ties going to negative skew
    if(
      sum(positive_skew_signs == -1) >=
      sum(positive_skew_signs == 1)
    ){

      # Make all skew for positive variables negative
      skews[signs == 1] <- -abs(skews[signs == 1])

      # Make all skew for negative variables positive
      skews[signs == -1] <- abs(skews[signs == -1])

    }else{ # Do positive skew for positive variables

      # Make all skew for positive variables positive
      skews[signs == 1] <- abs(skews[signs == 1])

      # Make all skew for negative variables negative
      skews[signs == -1] <- -abs(skews[signs == -1])

    }

  }

  # Return skews
  return(skews)

}

#' @noRd
# Adds skew to a single variable based on threshold (skew) values
# Updated 30.11.2022
skew_single_variable <- function(data, skew_values){

  # Categorize biased data with updated thresholds
  for(i in (length(skew_values) + 1):1){

    # First category
    if(i == 1){
      data[data < skew_values[i]] <- i
    }else if(i == length(skew_values) + 1){ # Last category
      data[data >= skew_values[i-1]] <- i
    }else{ # Middle category
      data[data >= skew_values[i-1] & data < skew_values[i]] <- i
    }

  }

  # Return skewed data
  return(data)

}

# Based on
# Garrido, L. E., Abad, F. J., & Ponsoda, V. (2011).
# Performance of Velicer’s minimum average partial factor retention
# method with categorical variables.
# Educational and Psychological Measurement, 71(3), 551-570.
# https://doi.org/10.1177/0013164410389489
#
#' @noRd
# Generates skewed data for continuous data ----
# Updated 13.03.2026
skew_continuous <- function(skew, data = NULL, sample_size = 1e06, tolerance = 1e-05)
{

  # Generate data
  if(is.null(data)){
    data <- rnorm_ziggurat(sample_size)
  }

  # Check for zero skew (add minimal skew)
  if(skew == 0){
    return(data)
  }

  # Obtain sign and absolute skew
  skew_sign <- sign(skew)
  absolute <- abs(skew)

  # Compute arcsinh once
  arcsinh_data <- asinh(data)

  # Initialize kurtosis
  kurtosis <- 1

  # Optimize for result
  objective <- function(kurtosis){

    # Optimize for skew
    optimize(
      f = function(x){
        abs(absolute - psych::skew(sinh(kurtosis * (arcsinh_data + x))))
      }, interval = c(0, absolute + 2), tol = tolerance
    )

  }

  # Repeat until convergence
  repeat{

    # Get result
    result <- silent_call(objective(kurtosis))

    # Check for minimum too low
    if(is.na(result$objective)){
      return(data)
    }

    # Check for tolerance
    if(result$objective < tolerance){
      break
    }

    # Update kurtosis
    kurtosis <- kurtosis + result$objective * 0.25 # learning rate

  }

  # Return result
  return(skew_sign * scale(sinh(kurtosis * (arcsinh_data + result$minimum))))

}

# Based on
# Garrido, L. E., Abad, F. J., & Ponsoda, V. (2011).
# Performance of Velicer’s minimum average partial factor retention
# method with categorical variables.
# Educational and Psychological Measurement, 71(3), 551-570.
# https://doi.org/10.1177/0013164410389489
#
#' @noRd
# Skewness of the categorical distribution implied by `proportion`
# (the cumulative share taken up by the lowest-value category),
# `reduction_factor`, `categories`, and `sample_size` -- computed
# analytically from the category counts instead of materializing a
# `sample_size`-length vector and calling `psych::skew()` on it.
# Mathematically identical to that (the type = 3 moment-skewness formula
# only depends on the counts per value, not on realized data order), just
# computed directly from the block-boundary indices, which are the same
# ones the original vector-based version rounded to when assigning blocks
# Updated 12.08.2026
categorical_skew <- function(
    categories, proportion, divided_proportion,
    reduction_factor, sample_size
)
{

  propinf <- 1 / sample_size
  propsup <- proportion
  E <- divided_proportion / reduction_factor
  counts <- numeric(categories)
  previous_index <- round(sample_size * propinf)

  # Convert cumulative proportions to per-category counts
  for(j in 1:categories){

    current_index <- round(sample_size * propsup)

    # The last category keeps its final index; earlier categories lose
    # theirs to the next category's starting index (the two ranges share
    # a boundary index, which the next category's assignment overwrites)
    counts[j] <- if(j < categories){
      current_index - previous_index
    }else{
      current_index - previous_index + 1
    }

    E <- E * reduction_factor
    propinf <- propsup
    propsup <- propinf + E
    previous_index <- current_index

  }

  # Moments, equivalent to `psych::skew(x, type = 3)`
  values <- seq_len(categories)
  n <- sum(counts)
  mean_value <- sum(counts * values) / n
  deviations <- values - mean_value
  variance <- sum(counts * deviations^2) / (n - 1)
  third_moment <- sum(counts * deviations^3)

  return(third_moment / (n * variance^1.5))

}

#' @noRd
# Sum of the geometric allocation weights (ratio = `reduction_factor`)
# used to split the categories after the first one
# Updated 12.08.2026
skew_allocation <- function(categories, reduction_factor)
{

  allocation <- 1
  allocation_1 <- 1

  if(categories > 2){
    for(j in 1:(categories - 2)){
      allocation_1 <- allocation_1 * reduction_factor
      allocation <- allocation + allocation_1
    }
  }

  return(allocation)

}

#' @noRd
# Finds `proportion` (the share of the lowest-value category) whose
# implied skewness matches `skewness`, for a fixed `reduction_factor`
# and `sample_size`. `skew(proportion)` is *not* monotonic over the full
# (0, 1) range -- it dips to an interior minimum before rising -- so this
# locates that minimum first (`optimize()`) and searches only the
# monotonic-increasing branch to its right. Bisecting from a fixed start
# (e.g., 0.5) can fail to converge whenever that minimum sits between the
# start and the true root, even when `skewness` is technically reachable.
# Returns NULL when `skewness` is unreachable for this `reduction_factor`/
# `sample_size` combination
# Updated 12.08.2026
skew_bisection <- function(
    skewness, categories, reduction_factor, sample_size,
    tolerance = 0.00001, maximum_iterations = 300
)
{

  allocation <- skew_allocation(categories, reduction_factor)

  skew_at <- function(proportion){
    divided_proportion <- (1 - proportion) / allocation
    categorical_skew(categories, proportion, divided_proportion, reduction_factor, sample_size)
  }

  floor <- optimize(skew_at, interval = c(0.00005, 0.99995))
  ceiling_value <- skew_at(0.99995)

  if(skewness < floor$objective - tolerance || skewness > ceiling_value + tolerance){
    return(NULL)
  }

  limitsinf <- floor$minimum
  limitsup <- 0.99995
  proportion <- (limitsinf + limitsup) / 2
  skew_actual <- skew_at(proportion)

  iteration <- 0
  while(abs(skew_actual - skewness) > tolerance){

    iteration <- iteration + 1
    if(iteration > maximum_iterations || (limitsup - limitsinf) < .Machine$double.eps * 4){
      return(NULL)
    }

    if(skew_actual < skewness){
      limitsinf <- proportion
      proportion <- (proportion + limitsup) / 2
    }else{
      limitsup <- proportion
      proportion <- (proportion + limitsinf) / 2
    }

    skew_actual <- skew_at(proportion)

  }

  return(list(
    proportion = proportion, allocation = allocation,
    reduction_factor = reduction_factor, sample_size = sample_size
  ))

}

#' @noRd
# Searches over `reduction_factor` and `sample_size` for a combination
# that can reach `skewness` for `categories`. The default
# `reduction_factor` at the smallest `sample_size` reaches the large
# majority of targets directly; `sample_size` only needs to grow when the
# target sits in a region where `skew(proportion)` is so steep that
# adjacent representable proportions jump past it; `reduction_factor` only
# needs to move off its default when `skewness` is outside the achievable
# range entirely (e.g., skew near zero for categories >= 4, which the
# default's shape cannot reach at any sample_size). When a search over
# `reduction_factor` is needed, the value closest to the default is kept,
# to minimize distortion from the "standard" shape
# Updated 12.08.2026
skew_search <- function(
    skewness, categories, tolerance = 0.00001,
    default_reduction_factor = 0.75,
    reduction_factor_grid = seq(0.02, 0.98, by = 0.02),
    sample_sizes = c(1e6, 1e8, 1e10, 1e12)
)
{

  # Fast path: default reduction_factor, escalating resolution
  for(sample_size in sample_sizes){
    fit <- skew_bisection(skewness, categories, default_reduction_factor, sample_size, tolerance)
    if(!is.null(fit)){
      return(fit)
    }
  }

  # Fallback: search reduction_factor too, preferring the smallest
  # sample_size and the reduction_factor closest to the default
  best <- NULL
  for(sample_size in sample_sizes){
    for(reduction_factor in reduction_factor_grid){
      fit <- skew_bisection(skewness, categories, reduction_factor, sample_size, tolerance)
      if(!is.null(fit)){
        if(is.null(best) || abs(fit$reduction_factor - default_reduction_factor) <
           abs(best$reduction_factor - default_reduction_factor)){
          best <- fit
        }
      }
    }
    if(!is.null(best)){
      break
    }
  }

  return(best)

}

#' @noRd
# Generates the `categories - 1` z-score breakpoints that split a
# `categories`-category variable to match `skewness`, on the fly rather
# than looked up from a precomputed table. Negative targets are generated
# by solving for `abs(skewness)` and mirroring (negate + reverse category
# order): this family (mass concentrated toward the lowest-value category,
# geometrically decaying across the rest) cannot reach deep negative skew
# for categories >= 3 directly, at any `reduction_factor` -- mirroring is
# the only robust way to get there. Returns NULL (with a warning) if
# `skewness` is not achievable for `categories`
# Updated 12.08.2026
skew_generator <- function(
    skewness, categories,
    reduction_factor = 0.75,
    sample_size = 1000000,
    tolerance = 0.00001
)
{

  fit <- skew_search(
    skewness = abs(skewness), categories = categories, tolerance = tolerance,
    default_reduction_factor = reduction_factor
  )

  if(is.null(fit)){
    warning("`skewness` is not achievable for this `categories`", call. = FALSE)
    return(NULL)
  }

  # Cumulative probabilities for categories 1..categories
  E <- (1 - fit$proportion) / fit$allocation / fit$reduction_factor
  cumulative_probability <- fit$proportion
  cumulative <- numeric(categories)

  for(j in 1:categories){
    cumulative[j] <- cumulative_probability
    E <- E * fit$reduction_factor
    cumulative_probability <- cumulative_probability + E
  }

  # Only `categories - 1` breakpoints are needed to split `categories` groups
  breakpoints <- qnorm(cumulative[1:(categories - 1)])
  breakpoints[is.infinite(breakpoints)] <- 0

  if(skewness < 0){
    breakpoints <- -rev(breakpoints)
  }

  return(list(
    breakpoints = breakpoints,
    reduction_factor = fit$reduction_factor,
    sample_size = fit$sample_size
  ))

}

#%%%%%%%%%%%%%%%%%%%
# data_to_zipfs ----
#%%%%%%%%%%%%%%%%%%%

#' @noRd
# Finds nearest non-zero decimal
# Updated 23.09.2022
nearest_decimal <- function(vec)
{

  # Obtain minimum
  minimum <- min(vec[vec!=0])

  # Zap digit
  zap_digit <- 0

  # Count
  digit <- -1

  # Find zap digit where it is not zero
  while(zap_digit == 0){

    # Increase digit
    digit <- digit + 1

    # Set zap
    zap_digit <- round(minimum, digits = digit)

  }

  # Return digit
  return(digit)

}

#%%%%%%%%%%%%%%%%%%%
# factor_forest ----
#%%%%%%%%%%%%%%%%%%%

# Comparison data approach (Ruscio & Roche, 2012)
#
# Optimized implementation of the original `EFA.Comp.Data()` (adapted from the
# Non-hierarchical-PRISM project): same algorithm and defaults, but much faster
# (random numbers are from the package's C random number generators so that
# results are reproducible with `seed`):
#   - Correlations of resampled data use BLAS (`crossprod`) rather than
#     `cor(use = "pairwise.complete.obs")`. The comparison population is built
#     by resampling observed non-missing values, so it never has missing values
#   - Only eigenvalues are computed where eigenvectors are not needed
#   - The shared + unique component matrix is computed once per iteration
#     rather than once per variable
#   - Imposing the observed distributions tracks the row permutation as an
#     index vector instead of re-sorting the whole population per variable
#   - Values that do not change across iterations are computed once
# With `cor()` the generated populations are identical to the original; the
# BLAS correlations differ from `cor()` at the ~1e-17 level

#' @noRd
# Correlation matrix of data that has no missing values
# Updated 20.11.2026
cd_cor <- function(x, cor_type)
{

  # Spearman needs the ranking done by `cor()`
  if(cor_type == "spearman"){
    return(cor(x, method = "spearman"))
  }

  # Pearson via BLAS
  x <- x - rep(colMeans(x), each = nrow(x))
  cross <- crossprod(x)
  scale <- sqrt(diag(cross))

  # Return correlations
  return(cross / tcrossprod(scale))

}

#' @noRd
# Principal axis factor analysis on a correlation matrix (loadings only)
# Updated 20.11.2026
cd_factor_analysis <- function(cor_matrix, n_factors, max_iter = 50, criterion = 0.001)
{

  # Initialize
  k <- ncol(cor_matrix)
  old_h2 <- rep(99, k)
  change <- 1
  iter <- 0

  # Iterate communalities
  while((change >= criterion) && (iter < max_iter)){

    # Next iteration
    iter <- iter + 1

    # Eigendecomposition
    eig <- eigen(cor_matrix, symmetric = TRUE)

    # Loadings and communalities
    loadings <- eig$vectors[, seq_len(n_factors), drop = FALSE] *
      rep(sqrt(eig$values[seq_len(n_factors)]), each = k)
    h2 <- rowSums(loadings * loadings)

    # Update
    change <- max(abs(old_h2 - h2))
    old_h2 <- h2
    diag(cor_matrix) <- h2

  }

  # Return loadings
  return(loadings)

}

#' @noRd
# Generate comparison population
# Updated 20.11.2026
cd_generate <- function(
    n_factors, N, cor_type, target_cor, distribution_values,
    max_trials = 5, initial_multiplier = 1
)
{

  # Initialize
  k <- ncol(target_cor)
  best_rmsr <- 1
  trials_without_improvement <- 0
  lower <- lower.tri(target_cor)
  n_lower <- sum(lower)

  # Distribution for each variable (sorted bootstrap of observed values)
  distributions <- matrix(0, nrow = N, ncol = k)
  for(i in seq_len(k)){
    distributions[, i] <- sort(
      shuffle_replace(distribution_values[[i]], size = N)
    )
  }

  # Shared and unique components
  shared_comp <- matrix(rnorm_ziggurat(N * n_factors), nrow = N, ncol = n_factors)
  unique_comp <- matrix(rnorm_ziggurat(N * k), nrow = N, ncol = k)

  # Loadings -> data with the observed distributions imposed
  # (rows are returned in the same order as the original implementation,
  # which keeps the order of floating point summation identical)
  build <- function(cor_matrix, standardize = FALSE){

    # Shared loadings (bounded, positive first loading)
    shared_load <- cd_factor_analysis(cor_matrix, n_factors)
    shared_load[shared_load > 1] <- 1
    shared_load[shared_load < -1] <- -1
    if(shared_load[1, 1] < 0){
      shared_load <- shared_load * -1
    }

    # Unique loadings
    h2 <- rowSums(shared_load * shared_load)
    unique_load <- sqrt(ifelse(h2 < 1, 1 - h2, 0))

    # Normally distributed data
    normal <- tcrossprod(shared_comp, shared_load) +
      unique_comp * rep(unique_load, each = N)
    if(standardize){
      normal <- scale(normal)
    }

    # Impose observed distributions: for each variable in turn, order the
    # rows by that variable and assign its sorted distribution (only the
    # row order is tracked)
    order_rows <- seq_len(N)
    imposed <- normal
    for(i in seq_len(k)){
      order_rows <- order_rows[sort.list(normal[order_rows, i])]
      imposed[order_rows, i] <- distributions[, i]
    }

    # Return population (rows in the original implementation's order)
    return(imposed[order_rows, , drop = FALSE])

  }

  # Search for the intermediate correlation matrix
  intermediate_cor <- target_cor
  while(trials_without_improvement < max_trials){

    # Population from current intermediate correlations
    population <- build(intermediate_cor)

    # Residual correlations
    reproduced_cor <- cd_cor(population, cor_type)
    residual_cor <- target_cor - reproduced_cor
    rmsr <- sqrt(sum(residual_cor[lower]^2) / n_lower)

    # Update
    if(rmsr < best_rmsr){
      best_rmsr <- rmsr
      best_cor <- intermediate_cor
      best_res <- residual_cor
      intermediate_cor <- intermediate_cor + initial_multiplier * residual_cor
      trials_without_improvement <- 0
    }else{
      trials_without_improvement <- trials_without_improvement + 1
      intermediate_cor <- best_cor +
        initial_multiplier * 0.5^trials_without_improvement * best_res
    }

  }

  # Final population (standardized as in the original)
  return(build(best_cor, standardize = TRUE))

}

#' @noRd
# Comparison data approach (Ruscio & Roche, 2012)
# Data = N x k data matrix; F.Max = largest number of factors considered;
# N.Pop = population size; N.Samples = samples per population;
# Alpha = significance level for adding a factor; Spearman = rank correlations;
# use = handling of missing values for the observed correlations
# Updated 20.11.2026
EFA.Comp.Data <- function(
    Data, F.Max = 8, N.Pop = 10000, N.Samples = 500, Alpha = 0.30,
    Spearman = FALSE, use = "pairwise.complete.obs"
)
{

  # Dimensions
  Data <- as.matrix(Data)
  N <- nrow(Data)
  k <- ncol(Data)
  cor_type <- swiftelse(isTRUE(Spearman), "spearman", "pearson")

  # A k-factor solution on k variables is not identified, so cap the search
  # at k - 1 factors (also avoids requesting more eigenvectors than exist)
  F.Max <- max(1, min(F.Max, k - 1))

  # Quantities that do not depend on the number of factors
  target_cor <- cor(Data, method = cor_type, use = use)
  eigs_data <- eigen(target_cor, symmetric = TRUE, only.values = TRUE)$values
  distribution_values <- lapply(seq_len(k), function(i){na.omit(Data[, i])})

  # Rows of the population (for resampling)
  population_rows <- seq_len(N.Pop)

  # Search over number of factors
  rmsr_eigs <- matrix(0, nrow = N.Samples, ncol = F.Max)
  significant <- TRUE
  factors <- 1

  while((factors <= F.Max) && significant){

    # Comparison population
    population <- cd_generate(
      n_factors = factors, N = N.Pop, cor_type = cor_type,
      target_cor = target_cor, distribution_values = distribution_values
    )

    # Fit of eigenvalues from samples of the population
    for(j in seq_len(N.Samples)){
      sample_cor <- cd_cor(
        population[shuffle_replace(population_rows, size = N), , drop = FALSE],
        cor_type
      )
      eigs_sample <- eigen(sample_cor, symmetric = TRUE, only.values = TRUE)$values
      rmsr_eigs[j, factors] <- sqrt(sum((eigs_sample - eigs_data)^2) / k)
    }

    # Is the added factor a significant improvement?
    if(factors > 1){
      significant <- wilcox.test(
        rmsr_eigs[, factors], rmsr_eigs[, factors - 1], "less"
      )$p.value < Alpha
    }
    if(significant){
      factors <- factors + 1
    }

  }

  # Return number of factors
  return(factors - 1)

}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# obtain_zipfs_parameters ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Obtains Zipf SSE
# Updated 27.09.2022
zipf_sse <- function(values, zipfs){
  sum((zipfs - values)^2, na.rm = TRUE)
}


#' @noRd
# Obtains Zipf's values
# Updated 27.09.2022
zipf_values <- function(alpha, beta, rank_order)
{1 / (rank_order + beta)^alpha}

#' @noRd
# Estimate parameters
# Updated 28.09.2022
estimate_parameters <- function(
    alpha_sequence,
    beta_sequence,
    zipfs,
    rank_order
)
{

  # All possible combinations
  sequences <- expand.grid(
    alpha = alpha_sequence,
    beta = beta_sequence
  )

  # Initialize RMSE vector
  sse <- numeric(nrow(sequences))

  # Possible values and differences
  for(i in 1:nrow(sequences)){

    # Progress message
    cat(
      colortext(
        text =  paste(
          "\r Estimating alpha...",
          formatC(
            sequences[i,1], digits = 2,
            format = "f", flag = "0"
          ), "  ",
          "Estimating beta...",
          formatC(
            sequences[i,2], digits = 2,
            format = "f", flag = "0"
          ), "  "
        ),
        defaults = "message"
      )
    )

    # Values based on parameters
    values <- zipf_values(
      alpha = sequences[i,1], beta = sequences[i,2],
      rank_order = rank_order
    )

    # Difference
    sse[i] <- zipf_sse(
      values = values,
      zipfs = zipfs
    )


  }

  # Update message
  cat(
    colortext(
      text =  paste(
        "\r Estimating alpha...",
        formatC(
          sequences[which.min(sse),1], digits = 2,
          format = "f", flag = "0"
        ), "  ",
        "Estimating beta...",
        formatC(
          sequences[which.min(sse),2], digits = 2,
          format = "f", flag = "0"
        ), "  "
      ),
      defaults = "message"
    )
  )

  # Return parameters
  return(sequences[which.min(sse),])

}

#' @noRd
# General function to silently obtain output ----
# Updated 01.07.2023
silent_call <- function(...)
{

  # Make call
  sink <- capture.output(
    result <- suppressWarnings(
      suppressMessages(
        ...
      )
    )
  )

  # Return result
  return(result)

}

#%%%%%%%%%%%%%%%%%%%%%
# ERROR FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Error for object type
# Updated 30.09.2022
object_error <- function(input, expected_type){

  # Check for possible object types
  possible_types <- sapply(
    X = expected_type,
    FUN = is,
    object = input
  )

  # Check for object types
  if(all(!possible_types)){
    stop(
      paste(
        "Input into '", deparse(substitute(input)),
        "' argument is not ", paste("'", expected_type, "'", sep = "", collapse = ", "),
        ". Input is ", paste("'", class(input), "'", sep = "", collapse = ", "),
        sep = ""
      )
    )
  }

}

#' @noRd
# Error for input type
# Updated 08.08.2022
type_error <- function(input, expected_type){

  # Check for type
  if(!is(input, expected_type)){
    stop(
      paste(
        "Input into '", deparse(substitute(input)),
        "' argument is not '", expected_type,
        "'. Input is ", paste("'", class(input), "'", sep = "", collapse = ", "),
        sep = ""
      )
    )
  }

}

#' @noRd
# Error for input length
# Updated 08.08.2022
length_error <- function(input, expected_lengths){

  # Check for length of input in expected length
  if(!length(input) %in% expected_lengths){
    stop(
      paste(
        "Length of '", deparse(substitute(input)),
        "' (", length(input),") does not match expected length(s). Length must be: ",
        paste("'", expected_lengths, "'", collapse = " or ", sep = ""),
        sep = ""
      )
    )
  }

}

#' @noRd
# Error for input range
# Updated 05.09.2022
range_error <- function(input, expected_ranges){

  # Obtain expected maximum and minimum values
  expected_maximum <- max(expected_ranges)
  expected_minimum <- min(expected_ranges)

  # Obtain maximum and minimum values
  # (compared unrounded, with a small tolerance for floating point error;
  # rounded only for the message)
  actual_maximum <- max(input)
  actual_minimum <- min(input)
  tolerance <- 1e-12

  # Check for maximum of input in expected range
  if(actual_maximum > expected_maximum + tolerance){
    actual_maximum <- signif(actual_maximum, 6)
    stop(
      paste(
        "Maximum of '", deparse(substitute(input)),
        "' (", actual_maximum,") does not match expected range(s). Range must be between: ",
        paste0("'", expected_ranges, "'", collapse = " and "),
        sep = ""
      )
    )
  }

  # Check for minimum of input in expected range
  if(actual_minimum < expected_minimum - tolerance){
    actual_minimum <- signif(actual_minimum, 6)
    stop(
      paste(
        "Minimum of '", deparse(substitute(input)),
        "' (", actual_minimum,") does not match expected range(s). Range must be between: ",
        paste0("'", expected_ranges, "'", collapse = " and "),
        sep = ""
      )
    )
  }

}

#%%%%%%%%%%%%%%%%%%%%%%
# SYSTEM FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%

#' System check for OS and RSTUDIO
#'
#' @description Checks for whether text options are available
#'
#' @param ... Additional arguments
#'
#' @return \code{TRUE} if text options are available and \code{FALSE} if not
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
# System Check
# Updated 08.09.2020
system.check <- function (...)
{
  OS <- unname(tolower(Sys.info()["sysname"]))

  RSTUDIO <- ifelse(Sys.getenv("RSTUDIO") == "1", TRUE, FALSE)

  TEXT <- TRUE

  if(!RSTUDIO){if(OS != "linux"){TEXT <- FALSE}}

  res <- list()

  res$OS <- OS
  res$RSTUDIO <- RSTUDIO
  res$TEXT <- TEXT

  return(res)
}

#' Colorfies Text
#'
#' Makes text a wide range of colors (8-bit color codes)
#'
#' @param text Character.
#' Text to color
#'
#' @return Colorfied text
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
#'
# Color text
# Updated 08.09.2020
colortext <- function(text, number = NULL, defaults = NULL)
{
  # Check system
  sys.check <- system.check()

  if(sys.check$TEXT)
  {
    # Defaults for number (white text)
    if(is.null(number) || number < 0 || number > 231)
    {number <- 15}

    # Check for default color
    if(!is.null(defaults))
    {
      # Adjust highlight color based on background color
      if(defaults == "highlight")
      {
        if(sys.check$RSTUDIO)
        {

          if(rstudioapi::getThemeInfo()$dark)
          {number <- 226
          }else{number <- 208}

        }else{number <- 208}
      }else{

        number <- switch(defaults,
                         message = 204,
                         red = 9,
                         orange = 208,
                         yellow = 11,
                         "light green" = 10,
                         green = 34,
                         cyan = 14,
                         blue = 12,
                         magenta = 13,
                         pink = 211,
        )

      }

    }

    return(paste("\033[38;5;", number, "m", text, "\033[0m", sep = ""))

  }else{return(text)}
}

#' Stylizes Text
#'
#' Makes text bold, italics, underlined, and strikethrough
#'
#' @param text Character.
#' Text to stylized
#'
#' @return Sytlized text
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
# Style text
# Updated 08.09.2020
styletext <- function(text, defaults = c("bold", "italics", "highlight",
                                         "underline", "strikethrough"))
{
  # Check system
  sys.check <- system.check()

  if(sys.check$TEXT)
  {
    if(missing(defaults))
    {number <- 0
    }else{

      # Get number code
      number <- switch(defaults,
                       bold = 1,
                       italics = 3,
                       underline = 4,
                       highlight = 7,
                       strikethrough = 9
      )

    }

    return(paste("\033[", number, ";m", text, "\033[0m", sep = ""))
  }else{return(text)}
}

#' Text Symbols
#'
#' Makes text symbols (star, checkmark, square root)
#'
#' @param symbol Character.
#'
#' @return Outputs symbol
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
# Symbols
# Updated 24.04.2020
textsymbol <- function(symbol = c("alpha", "beta", "chi", "delta",
                                  "eta", "gamma", "lambda", "omega",
                                  "phi", "pi", "rho", "sigma", "tau",
                                  "theta", "square root", "infinity",
                                  "check mark", "x", "bullet")
)
{
  # Get number code
  sym <- switch(symbol,
                alpha = "\u03B1",
                beta = "\u03B2",
                chi = "\u03C7",
                delta = "\u03B4",
                eta = "\u03B7",
                gamma = "\u03B3",
                lambda = "\u03BB,",
                omega = "\u03C9",
                phi = "\u03C6",
                pi = "\u03C0",
                rho = "\u03C1",
                sigma = "\u03C3",
                tau = "\u03C4",
                theta = "\u03B8",
                "square root" = "\u221A",
                infinity = "\u221E",
                "check mark" = "\u2713",
                x = "\u2717",
                bullet = "\u2022"
  )

  return(sym)
}

#%%%%%%%%%%%%%%%%%%%%%%%
# UTILITY FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Ensures column names to data
# Updated 02.12.2022
ensure_column_names <- function(data)
{

  # Check for whether column names exist
  if(is.null(colnames(data))){

    # Add column names
    colnames(data) <- paste0(
      "V", formatC(
        1:ncol(data), digits = floor(log10(ncol(data))),
        format = "d", flag = "0"
      )
    )

  }

  # Return data
  return(data)

}

#' @noRd
# Checks for duplicated rows
# Updated 05.09.2022
match_row <- function(data)
{
  # Make data frame
  df <- as.data.frame(data)

  # Obtain duplicate indices
  dupe_ind <- duplicated(df)

  # Return rows
  return(which(dupe_ind))



}

#' @noRd
# Function update default arguments
# Updated 12.12.2022
update_defaults <- function(FUN, FUN_args)
{

  # Exploratory Graph Analysis
  if(FUN == "EGA"){

    # Defaults
    defaults <- list(
      corr = "cor_auto", uni.method = "louvain",
      model = "glasso", consensus.method = "most_common",
      plot.EGA = FALSE
    )

  }

  # Factor Forest
  if(FUN == "FF"){

    # Defaults
    defaults <- list(
      maximum_factors = 8,
      PA_correlation = "cor"
    )

  }

  # Out-of-sample Factor Analysis
  if(FUN == "FSPE"){

    # Defaults
    defaults <- list(
      maxK = 8, rep = 1, method = "PE", pbar = FALSE
    )

  }

  # Next Eigenvalue Sufficiency Test
  if(FUN == "NEST"){

    # Defaults
    defaults <- list(
      iterations = 1000,
      maximum_iterations = 500,
      alpha = 0.05,
      convergence = 0.00001
    )

  }

  # Parallel Analysis
  if(FUN == "PA"){

    # Defaults
    defaults <- list(
      fm = "minres", fa = "both", cor = "cor",
      n.iter = 20, sim = FALSE, plot = FALSE
    )

  }

  # Check for any function arguments
  if(any(names(FUN_args) %in% names(defaults))){

    # Obtain indices
    args_index <- which(names(FUN_args) %in% names(defaults))

    # Replace arguments
    defaults[names(FUN_args)[args_index]] <- FUN_args[args_index]

  }

  # Return arguments
  return(defaults)

}

#' @noRd
# Function to obtain arguments
# Updated 30.09.2022
obtain_arguments <- function(FUN, FUN_args)
{

  # Obtain formal arguments
  FUN_formals <- formals(FUN)

  # Check for input arguments
  if(length(FUN_args) != 0){

    ## Check for matching arguments
    if(any(names(FUN_args) %in% names(FUN_formals))){

      replace_args <- FUN_args[na.omit(match(names(FUN_formals), names(FUN_args)))]

      FUN_formals[names(replace_args)] <- replace_args
    }

  }

  # Remove ellipses
  if("..." %in% names(FUN_formals)){
    FUN_formals[which(names(FUN_formals) == "...")] <- NULL
  }

  # Return agrumnets
  return(FUN_formals)

}











