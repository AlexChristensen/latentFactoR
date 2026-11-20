#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Population Error Optimization Functions ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Adapted from {psychNetwoRk} for factor models. The "true model" is the
# population factor model (Lambda Phi Lambda' + Psi) with a known pattern of
# free (non-zero) parameters. Population error is added to the population
# correlation matrix such that the fit of the true model (constrained maximum
# likelihood) to the correlation matrix with error equals a target value of a
# fit index. Fit indices are computed from the maximum likelihood discrepancy
# at the size of the human population (so they are at their large-sample limits)

#' @noRd
# Population size (approximately the human population)
# Updated 20.11.2026
population_size <- 8.3e09

# True model ----

#' @noRd
# Obtain the true factor model (loadings, factor correlations,
# uniquenesses/residuals, and population correlation matrix)
# from a {latentFactoR} object. Correlated residuals (e.g., from
# `add_local_dependence`) are the part of the population correlation matrix
# that is not explained by the factors. Method factors (from
# `add_method_factors`) are part of the true model
# Updated 20.11.2026
true_model <- function(lf_object, leave_cross_loadings = FALSE)
{

  # Obtain parameters
  parameters <- lf_object$parameters
  loadings <- parameters$loadings
  phi <- parameters$factor_correlations
  substantive <- parameters$factors

  # Correlated residuals: population correlations not explained by the
  # factors (computed with the generating loadings, before any
  # cross-loadings are removed from the true model)
  psi <- lf_object$population_correlation -
    loadings %*% tcrossprod(phi, loadings)
  diag(psi) <- 0
  psi[abs(psi) < 1e-10] <- 0
  psi <- unname(psi)

  # Check for whether cross-loadings should remain
  if(!isTRUE(leave_cross_loadings)){

    # Set sequence of variables for each factor
    variable_sequence <- factor_variable_sequence(parameters$variables)
    start_variables <- variable_sequence$start
    end_variables <- variable_sequence$end

    # Loop through substantive factors (method factors only load
    # on their own factor's variables)
    for(i in seq_len(substantive)){

      # Set cross-loadings to zero
      loadings[start_variables[i]:end_variables[i], setdiff(seq_len(substantive), i)] <- 0

    }

  }

  # Method factors
  if(ncol(loadings) > substantive){

    # Method factor columns
    methods <- seq.int(substantive + 1, ncol(loadings))

    # Remove method factors without loadings
    empty <- methods[colSums(loadings[, methods, drop = FALSE] != 0) == 0]
    if(length(empty) > 0){
      loadings <- loadings[, -empty, drop = FALSE]
      phi <- phi[-empty, -empty, drop = FALSE]
      methods <- seq.int(substantive + 1, length.out = ncol(loadings) - substantive)
    }

    # Perfectly correlated method factors are a single method factor
    # (Lambda_m 1 1' Lambda_m' = (Lambda_m 1)(Lambda_m 1)')
    if(
      length(methods) > 1 &&
      all(abs(phi[methods, methods] - 1) < 1e-12) &&
      all(phi[methods, -methods] == 0)
    ){
      combined <- rowSums(loadings[, methods, drop = FALSE])
      loadings <- cbind(loadings[, -methods, drop = FALSE], combined)
      phi <- rbind(
        cbind(phi[-methods, -methods, drop = FALSE], 0),
        c(rep(0, substantive), 1)
      )
    }

  }

  # Obtain common part of the population correlation matrix
  common <- loadings %*% tcrossprod(phi, loadings)

  # Uniquenesses (diagonal of the residuals)
  diag(psi) <- 1 - diag(common)

  # Check for inadmissible uniquenesses
  if(any(diag(psi) <= 0)){
    stop(
      paste(
        "The true model has communalities greater than or equal to one",
        "(non-positive uniquenesses).",
        "\n\nPopulation error cannot be added to an inadmissible model"
      ), call. = FALSE
    )
  }

  # Population correlation matrix (diagonal is exactly one)
  R <- common + psi
  diag(R) <- 1

  # Add names
  names <- variable_names(nrow(R))
  dimnames(R) <- list(names, names)
  loadings <- unname(loadings); phi <- unname(phi)

  # Return true model
  return(list(loadings = loadings, phi = phi, psi = psi, R = R))

}

#' @noRd
# Parameter structure of a factor model
# Free parameters are the non-zero loadings, non-zero factor correlations
# (lower triangle), and the uniquenesses and non-zero residual correlations
# (lower triangle including diagonal). Factor variances are fixed to one
# Updated 20.11.2026
factor_model_parts <- function(loadings, phi, psi)
{

  # Obtain dimensions
  nodes <- nrow(loadings)
  factors <- ncol(loadings)

  # Identify free parameters
  lambda_free <- which(loadings != 0)
  phi_free <- which(lower.tri(phi) & phi != 0)
  psi_free <- which(lower.tri(psi, diag = TRUE) & (psi != 0 | diag(nodes) == 1))

  # Count parameters
  parameters <- length(lambda_free) + length(phi_free) + length(psi_free)

  # Return parts
  return(
    list(
      nodes = nodes, factors = factors,
      lambda_free = lambda_free, phi_free = phi_free, psi_free = psi_free,
      parameters = parameters,
      df = nodes * (nodes + 1) / 2 - parameters,
      start = c(loadings[lambda_free], phi[phi_free], psi[psi_free])
    )
  )

}

#' @noRd
# Model matrices from the parameter vector
# Updated 20.11.2026
model_matrices <- function(theta, parts)
{

  # Obtain counts
  nodes <- parts$nodes
  factors <- parts$factors
  n_lambda <- length(parts$lambda_free)
  n_phi <- length(parts$phi_free)

  # Loadings
  lambda <- matrix(0, nrow = nodes, ncol = factors)
  lambda[parts$lambda_free] <- theta[seq_len(n_lambda)]

  # Factor correlations (unit diagonal)
  phi <- diag(factors)
  phi[parts$phi_free] <- theta[n_lambda + seq_len(n_phi)]
  phi <- phi + t(phi) - diag(factors)

  # Residuals
  psi <- matrix(0, nrow = nodes, ncol = nodes)
  psi[parts$psi_free] <- theta[n_lambda + n_phi + seq_len(length(parts$psi_free))]
  psi <- psi + t(psi) - diag(diag(psi), nodes)

  # Return matrices
  return(
    list(
      lambda = lambda, phi = phi, psi = psi,
      sigma = lambda %*% tcrossprod(phi, lambda) + psi
    )
  )

}

#' @noRd
# Every free parameter has a derivative of the form u v' + v u'
# (loadings: u = e_i, v = (Lambda Phi)[, k]; factor correlations: u = lambda_k,
# v = lambda_l; residuals: u = e_i, v = e_j and u = v = e_i / sqrt(2) for the
# diagonal). This returns U and V (nodes by parameters)
# Updated 20.11.2026
model_uv <- function(matrices, parts)
{

  # Obtain dimensions and matrices
  nodes <- parts$nodes
  factors <- parts$factors
  lambda <- matrices$lambda
  lambda_phi <- lambda %*% matrices$phi

  # Initialize
  U <- V <- matrix(0, nrow = nodes, ncol = parts$parameters)

  # Loadings
  n_lambda <- length(parts$lambda_free)
  if(n_lambda > 0){

    # Obtain row and column
    columns <- seq_len(n_lambda)
    i <- (parts$lambda_free - 1) %% nodes + 1
    k <- (parts$lambda_free - 1) %/% nodes + 1

    # Set U and V
    U[cbind(i, columns)] <- 1
    V[, columns] <- lambda_phi[, k, drop = FALSE]

  }

  # Factor correlations
  n_phi <- length(parts$phi_free)
  if(n_phi > 0){

    # Obtain row and column
    columns <- n_lambda + seq_len(n_phi)
    k <- (parts$phi_free - 1) %% factors + 1
    l <- (parts$phi_free - 1) %/% factors + 1

    # Set U and V
    U[, columns] <- lambda[, k, drop = FALSE]
    V[, columns] <- lambda[, l, drop = FALSE]

  }

  # Residuals
  n_psi <- length(parts$psi_free)
  columns <- n_lambda + n_phi + seq_len(n_psi)
  i <- (parts$psi_free - 1) %% nodes + 1
  j <- (parts$psi_free - 1) %/% nodes + 1
  value <- swiftelse(i == j, 1 / sqrt(2), 1)
  U[cbind(i, columns)] <- value
  V[cbind(j, columns)] <- value

  # Return U and V
  return(list(U = U, V = V))

}

#' @noRd
# Jacobian of vec(Sigma) with respect to the parameters
# (nodes^2 rows by number of parameters columns; used for validation)
# Updated 20.11.2026
model_jacobian <- function(matrices, parts)
{

  # Obtain U and V
  uv <- model_uv(matrices, parts)

  # Derivative of each parameter: u v' + v u'
  jacobian <- vapply(
    seq_len(parts$parameters),
    function(column){
      derivative <- tcrossprod(uv$U[, column], uv$V[, column])
      as.vector(derivative + t(derivative))
    },
    numeric(parts$nodes^2)
  )

  # Return Jacobian
  return(matrix(jacobian, nrow = parts$nodes^2))

}

#' @noRd
# Expected (Fisher) information of the parameters
#   tr(Sigma^-1 D_k Sigma^-1 D_l) with D = u v' + v u'
# Updated 20.11.2026
model_information <- function(uv, inverse)
{

  # Obtain products
  inverse_U <- inverse %*% uv$U
  inverse_V <- inverse %*% uv$V
  A <- crossprod(uv$U, inverse_U)
  B <- crossprod(uv$V, inverse_V)
  C <- crossprod(uv$U, inverse_V)

  # Return information
  return(2 * (A * B + C * t(C)))

}

#' @noRd
# Check that the true model is identified
# (expected information at the true parameters must be full rank)
# Updated 20.11.2026
model_identified <- function(parts)
{

  # Obtain information at the true parameters
  matrices <- model_matrices(parts$start, parts)
  information <- model_information(
    model_uv(matrices, parts), solve(matrices$sigma)
  )

  # Scale to a correlation-like matrix (rank is scale-free)
  scaling <- 1 / sqrt(diag(information))
  information <- information * tcrossprod(scaling)

  # Return whether full rank
  return(qr(information, tol = 1e-08)$rank == parts$parameters)

}

#' @noRd
# Reference (model-implied) matrix with inverse and log-determinant
# Updated 20.11.2026
reference_model <- function(R, R_chol = chol(R))
{

  # Return reference
  return(
    list(
      R = R, inverse = chol2inv(R_chol),
      log_det = 2 * sum(log(diag(R_chol)))
    )
  )

}

# Known model ----

#' @noRd
# Known model (constrained maximum likelihood) fit to a correlation matrix
# Fisher scoring (analytic score and expected information) with step halving,
# started at the true parameters (or a warm start); falls back to a
# quasi-Newton optimizer (followed by Fisher scoring) if Fisher scoring
# does not converge. Convergence is declared when the largest score is below
# 1e-09 (the discrepancy is then within ~1e-18 of its minimum) or when the
# discrepancy no longer decreases beyond rounding error
# Updated 20.11.2026
known_cfa <- function(S, parts, start = parts$start, max.iter = 100, S_chol = chol(S))
{

  # Obtain nodes and log-determinant of S
  nodes <- parts$nodes
  log_det_S <- 2 * sum(log(diag(S_chol)))

  # Maximum likelihood discrepancy (and optionally score and information)
  # (invalid when not positive definite)
  state <- function(theta, derivatives = TRUE){

    # Obtain model-implied matrix
    matrices <- model_matrices(theta, parts)
    sigma_chol <- try(chol(matrices$sigma), silent = TRUE)
    if(is(sigma_chol, "try-error")){
      return(list(valid = FALSE))
    }

    # Obtain discrepancy
    inverse <- chol2inv(sigma_chol)
    log_det <- 2 * sum(log(diag(sigma_chol)))
    discrepancy <- sum(S * inverse) + log_det - log_det_S - nodes
    result <- list(
      valid = is.finite(discrepancy), discrepancy = discrepancy,
      matrices = matrices, inverse = inverse, log_det = log_det
    )
    if(!derivatives || !result$valid){
      return(result)
    }

    # Score: 2 * u_k' M v_k with M = Sigma^-1 (Sigma - S) Sigma^-1
    uv <- model_uv(matrices, parts)
    M <- inverse %*% (matrices$sigma - S) %*% inverse
    result$score <- 2 * colSums(uv$U * (M %*% uv$V))

    # Expected information
    result$information <- model_information(uv, inverse)

    # Return state
    return(result)

  }

  # Fisher scoring
  scoring <- function(theta){

    # Initial state
    current <- state(theta)
    if(!current$valid){
      return(list(theta = theta, state = current))
    }

    # Loop through iterations
    for(iteration in seq_len(max.iter)){

      # Check for convergence
      if(max(abs(current$score)) < 1e-09){
        break
      }

      # Obtain step (add a small ridge if the information is singular)
      step <- try(solve(current$information, current$score), silent = TRUE)
      if(is(step, "try-error")){
        ridge <- 1e-08 * mean(diag(current$information))
        step <- try(
          solve(current$information + diag(ridge, length(theta)), current$score),
          silent = TRUE
        )
        if(is(step, "try-error")){
          break
        }
      }

      # Step halving (candidates are evaluated without derivatives)
      improved <- FALSE
      for(half in 0:30){

        # Candidate
        candidate <- theta - step / 2^half
        proposal <- state(candidate, derivatives = FALSE)

        # Check for a decrease in the discrepancy
        if(proposal$valid && proposal$discrepancy < current$discrepancy){
          improved <- TRUE
          break
        }

      }

      # Check for no improvement (minimum is reached within rounding error)
      if(!improved){
        break
      }

      # Update
      decrease <- current$discrepancy - proposal$discrepancy
      theta <- candidate
      current <- state(theta)

      # Check for convergence (decrease is at the level of rounding error)
      if(
        decrease < 8 * .Machine$double.eps * max(1, abs(current$discrepancy)) &&
        max(abs(current$score)) < 1e-08
      ){
        break
      }

    }

    # Return solution
    return(list(theta = theta, state = current))

  }

  # Try Fisher scoring first
  solution <- scoring(start)
  converged <- isTRUE(solution$state$valid) && max(abs(solution$state$score)) < 1e-06

  # Otherwise, quasi-Newton optimization followed by Fisher scoring
  if(!converged){

    # Objective and gradient
    optimization <- try(
      suppressWarnings(
        nlminb(
          start = start,
          objective = function(theta){
            proposal <- state(theta, derivatives = FALSE)
            if(proposal$valid){proposal$discrepancy}else{1e10}
          },
          gradient = function(theta){
            proposal <- state(theta)
            if(proposal$valid){proposal$score}else{numeric(length(theta))}
          },
          control = list(
            iter.max = 1000, eval.max = 2000, rel.tol = 1e-14, x.tol = 1e-12
          )
        )
      ), silent = TRUE
    )

    # Polish with Fisher scoring
    if(!is(optimization, "try-error")){
      polished <- scoring(optimization$par)
      if(
        isTRUE(polished$state$valid) && (
          !isTRUE(solution$state$valid) ||
          polished$state$discrepancy <= solution$state$discrepancy
        )
      ){
        solution <- polished
      }
    }

  }

  # Obtain final state
  final <- solution$state
  valid <- isTRUE(final$valid)
  max_gradient <- if(valid){max(abs(final$score))}else{Inf}

  # Return results
  return(
    list(
      lambda = final$matrices$lambda, phi = final$matrices$phi,
      psi = final$matrices$psi, sigma = final$matrices$sigma,
      theta = solution$theta,
      discrepancy = if(valid){final$discrepancy}else{Inf},
      reference = if(valid){
        list(R = final$matrices$sigma, inverse = final$inverse, log_det = final$log_det)
      },
      gradient = max_gradient,
      heywood = valid && any(diag(final$matrices$psi) < 0),
      converged = valid && max_gradient < 1e-06
    )
  )

}

#' @noRd
# Fit indices of a model-implied matrix to a population matrix
# Updated 20.11.2026
fit_indices <- function(
    R_population, reference, df, n = population_size,
    population_chol = try(chol(R_population), silent = TRUE)
)
{

  # Obtain number of nodes
  nodes <- ncol(R_population)

  # Check positive definiteness
  if(is(population_chol, "try-error")){
    return(list(valid = FALSE))
  }

  # Log-determinant of population matrix
  log_det_population <- 2 * sum(log(diag(population_chol)))

  # Return fit
  return(
    fit_from_discrepancy(
      discrepancy = reference$log_det + sum(R_population * reference$inverse) -
        log_det_population - nodes,
      log_det_population = log_det_population,
      residual = R_population - reference$R,
      nodes = nodes, df = df, n = n,
      R_population = R_population
    )
  )

}

#' @noRd
# Fit indices from the maximum likelihood discrepancy
#   F = log|Sigma_model| + tr(Sigma_pop Sigma_model^-1) - log|Sigma_pop| - p
# and the log-determinant of the population matrix (independence baseline)
# (SRMR is from the residuals or, if supplied, `SRMR` directly)
# Updated 20.11.2026
fit_from_discrepancy <- function(
    discrepancy, log_det_population, residual = NULL, nodes, df,
    n = population_size, R_population = NULL, SRMR = NULL
)
{

  # Chi-square
  chi_square <- n * discrepancy

  # Baseline chi-square (independence model): -log|Sigma_pop|
  baseline_chi_square <- -n * log_det_population
  baseline_df <- nodes * (nodes - 1) / 2

  # Guard against non-finite results
  if(!is.finite(chi_square) || !is.finite(baseline_chi_square)){
    return(list(valid = FALSE))
  }

  # Chi-square above degrees of freedom
  model_excess <- chi_square - df
  baseline_excess <- baseline_chi_square - baseline_df

  # Check for degenerate baseline
  denominator <- max(baseline_excess, model_excess, 0)
  baseline_ratio <- baseline_chi_square / baseline_df
  if(denominator == 0 || baseline_ratio <= 1){
    return(list(valid = FALSE))
  }

  # Standardized root mean square residual (lower triangle including the
  # diagonal, p * (p + 1) / 2 unique elements; same as {lavaan} for correlations)
  if(is.null(SRMR)){
    SRMR <- sqrt(mean(residual[lower.tri(residual, diag = TRUE)]^2))
  }

  # Return fit
  return(
    list(
      valid = TRUE, R_population = R_population,
      discrepancy = discrepancy, chi_square = chi_square,
      baseline_chi_square = baseline_chi_square,
      baseline_df = baseline_df, df = df,
      CFI = 1 - max(model_excess, 0) / denominator,
      TLI = (baseline_ratio - chi_square / df) / (baseline_ratio - 1),
      RMSEA = sqrt(max(model_excess, 0) / (df * n)),
      SRMR = SRMR
    )
  )

}

#' @noRd
# Fit indices of the known model (constrained maximum likelihood)
# Updated 20.11.2026
refit_indices <- function(R_population, parts, start = parts$start, n = population_size)
{

  # Check positive definiteness
  population_chol <- try(chol(R_population), silent = TRUE)
  if(is(population_chol, "try-error")){
    return(list(valid = FALSE))
  }

  # Obtain known model
  known <- try(
    known_cfa(R_population, parts, start = start, S_chol = population_chol),
    silent = TRUE
  )

  # Check for errors and non-convergence
  # (a non-converged solution overstates misfit)
  if(is(known, "try-error") || !isTRUE(known$converged)){
    return(list(valid = FALSE))
  }

  # Obtain fit of the known model (decompositions are reused)
  fit <- fit_indices(
    R_population, known$reference, parts$df, n,
    population_chol = population_chol
  )

  # Add known model
  if(fit$valid){
    fit$model <- known
  }

  # Return fit
  return(fit)

}

#' @noRd
# Refit function with warm starts
# Each refit starts from the solution at the nearest previously evaluated
# parameter (falling back to the true parameters). Returns a function of the
# scalar parameter and the population matrix function
# Updated 20.11.2026
warm_refit <- function(population, parts, n = population_size)
{

  # Store of evaluated parameters and solutions
  evaluated <- numeric(0)
  solutions <- list()

  # Return function
  return(
    function(parameter, cold = FALSE){

      # Population correlation matrix with error
      R_population <- population(parameter)
      if(is.null(R_population)){
        return(list(valid = FALSE))
      }

      # Starting values from the nearest evaluated parameter
      start <- parts$start
      if(!isTRUE(cold) && length(evaluated) > 0){
        start <- solutions[[which.min(abs(evaluated - parameter))]]
      }

      # Fit known model
      fit <- refit_indices(R_population, parts, start = start, n = n)

      # Store solution
      if(fit$valid){
        evaluated <<- c(evaluated, parameter)
        solutions[[length(solutions) + 1]] <<- fit$model$theta
      }

      # Return fit
      return(fit)

    }
  )

}

# Scalar solver ----

#' @noRd
# Solve for a scalar parameter that reaches the target metric
# The parameter increases with the amount of error and `lower` corresponds to
# (essentially) no error, where the metric is at its best value (CFI/TLI near
# one, RMSEA/SRMR near zero). The smallest error that reaches the target is
# found by searching upward from `lower` (geometric steps) to the first
# change of sign, then refining with Brent's method. Evaluations are memoized
# Updated 20.11.2026
solve_scalar <- function(
    FUN, lower, upper, metric, target, steps = 10, tolerance = 1e-10,
    spacing = c("geometric", "linear")
)
{

  # Set spacing of the search
  spacing <- match.arg(spacing)

  # Sign of the difference at no error
  # (CFI/TLI are above and RMSEA/SRMR are below their targets)
  zero_sign <- swiftelse(metric %in% c("CFI", "TLI"), 1, -1)

  # Memoized evaluations
  parameters <- numeric(0)
  fits <- list()
  evaluate <- function(parameter){

    # Check for previous evaluation
    index <- which(parameters == parameter)
    if(length(index) > 0){
      return(fits[[index[1]]])
    }

    # Evaluate
    fit <- FUN(parameter)
    parameters <<- c(parameters, parameter)
    fits[[length(fits) + 1]] <<- fit

    # Return fit
    return(fit)

  }

  # Difference between metric and target (missing when invalid)
  difference <- function(parameter){
    fit <- evaluate(parameter)
    return(swiftelse(isTRUE(fit$valid), fit[[metric]] - target, NA))
  }

  # Search upward from the smallest error: geometric steps
  # (lower + (upper - lower) / 2^k) or equally spaced steps
  grid <- switch(
    spacing,
    "geometric" = lower + (upper - lower) / 2^(rev(seq_len(steps)) - 1),
    "linear" = seq(lower, upper, length.out = steps + 1)[-1]
  )
  previous <- lower; previous_value <- zero_sign * abs(target)
  bracket <- NULL; first_invalid <- NULL
  for(point in grid){

    # Obtain difference
    value <- difference(point)

    # Skip invalid points (bracket only between valid points)
    if(is.na(value)){
      if(is.null(first_invalid)){first_invalid <- point}
      next
    }
    first_invalid <- NULL

    # Check for exact solution
    if(value == 0){
      return(c(evaluate(point), list(parameter = point)))
    }

    # Check for change of sign
    if(sign(value) != zero_sign){
      bracket <- c(previous, point)
      bracket_values <- c(previous_value, value)
      break
    }

    # Update previous valid point
    previous <- point; previous_value <- value

  }

  # The search ended in invalid points (e.g., errors outside of -1 to 1):
  # refine toward the boundary of valid errors looking for the target
  if(is.null(bracket) && !is.null(first_invalid)){

    # Interval between the last valid and the first invalid point
    a <- previous; b <- first_invalid
    for(i in seq_len(30)){

      # Middle point
      middle <- (a + b) / 2
      value <- difference(middle)

      # Update interval
      if(is.na(value)){
        b <- middle
      }else if(value == 0){
        return(c(evaluate(middle), list(parameter = middle)))
      }else if(sign(value) != zero_sign){
        bracket <- c(a, middle)
        bracket_values <- c(previous_value, value)
        break
      }else{
        a <- middle; previous_value <- value
      }

    }

  }

  # Target is not reachable so use the closest evaluated parameter
  if(is.null(bracket)){

    # Valid evaluations
    valid <- vapply(fits, function(x){isTRUE(x$valid)}, logical(1))
    if(!any(valid)){
      return(list(valid = FALSE))
    }

    # Closest
    distances <- vapply(
      fits[valid], function(x){abs(x[[metric]] - target)}, numeric(1)
    )
    parameter <- parameters[valid][which.min(distances)]
    return(c(evaluate(parameter), list(parameter = parameter)))

  }

  # Refine with Brent's method (invalid points inside the bracket are
  # treated as beyond the target so the bracket shrinks toward smaller error)
  root <- try(
    uniroot(
      function(parameter){
        value <- difference(parameter)
        swiftelse(is.na(value), bracket_values[2], value)
      },
      interval = bracket, f.lower = bracket_values[1], f.upper = bracket_values[2],
      tol = tolerance * max(1, abs(upper - lower)), maxiter = 200
    ), silent = TRUE
  )

  # Obtain parameter (closest valid bracket end if root finding failed)
  if(is(root, "try-error")){
    parameter <- bracket[2]
  }else{

    # Choose the evaluated point closest to the target near the root
    # (the root itself is always evaluated)
    parameter <- root$root
    if(is.na(difference(parameter))){
      parameter <- bracket[2]
    }

  }

  # Return fit
  return(c(evaluate(parameter), list(parameter = parameter)))

}

#' @noRd
# Solve with warm-started refits and confirm with a refit from the true
# parameters (guards against the warm start following a different minimum)
# Updated 20.11.2026
solve_refit <- function(
    population, parts, lower, upper, metric, target,
    n = population_size, spacing = "geometric"
)
{

  # Solve with warm starts
  refit <- warm_refit(population, parts, n)
  solution <- solve_scalar(refit, lower, upper, metric, target, spacing = spacing)

  # Check for no solution
  if(!isTRUE(solution$valid)){
    return(solution)
  }

  # Refit from the true parameters at the solution
  cold <- refit(solution$parameter, cold = TRUE)

  # Same minimum: return solution
  if(
    isTRUE(cold$valid) &&
    abs(cold$discrepancy - solution$discrepancy) <= 1e-10 * max(1, solution$discrepancy)
  ){
    return(solution)
  }

  # Otherwise, solve again starting every refit from the true parameters
  return(
    solve_scalar(
      function(parameter){refit(parameter, cold = TRUE)},
      lower, upper, metric, target, spacing = spacing
    )
  )

}

# Cudeck and Browne (1992) ----

#' @noRd
# Cudeck and Browne (1992) preparation
# Error must be orthogonal to the maximum likelihood score of every free
# parameter of the true model so that the true parameters remain the
# minimizer of the discrepancy (i.e., the true model is the best-fitting
# model of the population correlation matrix with error)
# Updated 20.11.2026
CB_prepare <- function(R, parts)
{

  # Obtain nodes
  nodes <- parts$nodes

  # Cells of the error (upper triangle; diagonal remains one)
  cells <- which(upper.tri(R))
  ij <- arrayInd(cells, c(nodes, nodes))

  # Reference
  R_chol <- chol(R)
  reference <- reference_model(R, R_chol)

  # Score of each parameter in the error cells:
  #   R^-1 (u v' + v u') R^-1 = a b' + b a', a = R^-1 u, b = R^-1 v
  uv <- model_uv(model_matrices(parts$start, parts), parts)
  A <- reference$inverse %*% uv$U
  B <- reference$inverse %*% uv$V
  scores <- A[ij[, 1], , drop = FALSE] * B[ij[, 2], , drop = FALSE] +
    B[ij[, 1], , drop = FALSE] * A[ij[, 2], , drop = FALSE]

  # Orthonormal basis of the scores
  decomposition <- qr(scores, tol = 1e-08)
  basis <- qr.Q(decomposition)[, seq_len(decomposition$rank), drop = FALSE]

  # Return preparation
  return(
    list(
      cells = cells, reference = reference, basis = basis,
      inverse_chol = backsolve(R_chol, diag(nodes)),
      dimension = length(cells) - decomposition$rank
    )
  )

}

#' @noRd
# Project vector onto the space orthogonal to the scores
# Updated 20.11.2026
CB_project <- function(vector, preparation)
{
  return(vector - preparation$basis %*% crossprod(preparation$basis, vector))
}

#' @noRd
# Cudeck and Browne (1992) population correlation matrix
# Updated 20.11.2026
CB_population <- function(error_vector, preparation)
{

  # Set up error matrix (upper triangle)
  nodes <- ncol(preparation$reference$R)
  error_matrix <- matrix(0, nodes, nodes)
  error_matrix[preparation$cells] <- error_vector

  # Population correlation matrix with error
  return(preparation$reference$R + (error_matrix + t(error_matrix)))

}

#' @noRd
# Cudeck and Browne (1992) error function
# (the true model is the minimizer so the reference is the original matrix)
# Updated 20.11.2026
CB_error <- function(error_vector, preparation, df, n = population_size)
{
  return(
    fit_indices(
      CB_population(error_vector, preparation), preparation$reference, df, n
    )
  )
}

#' @noRd
# Cudeck and Browne (1992) closed-form fit along a direction
# For R_population = R + s D with eigenvalues l of U^-T D U^-1 (R = U'U):
#   log|R_population| = log|R| + sum(log(1 + s l))
#   F(s) = s sum(l) - sum(log(1 + s l))
#   SRMR(s) = s * SRMR(D)
# Updated 20.11.2026
CB_direction_fit <- function(scale, eigenvalues, unit_SRMR, preparation, df, n = population_size)
{

  # Positive definiteness
  shifted <- 1 + scale * eigenvalues
  if(any(shifted <= 0)){
    return(list(valid = FALSE))
  }

  # Return fit
  return(
    fit_from_discrepancy(
      discrepancy = scale * sum(eigenvalues) - sum(log(shifted)),
      log_det_population = preparation$reference$log_det + sum(log(shifted)),
      nodes = ncol(preparation$reference$R), df = df, n = n,
      SRMR = scale * unit_SRMR
    )
  )

}

#' @noRd
# Cudeck and Browne (1992) cost function (error is projected)
# Updated 20.11.2026
CB_error_cost <- function(vector, metric, target, preparation, df, n = population_size)
{

  # Obtain fit
  fit <- CB_error(as.vector(CB_project(vector, preparation)), preparation, df, n)

  # Check for invalid results
  if(!fit$valid){
    return(1000)
  }

  # Smooth cost
  return((fit[[metric]] - target)^2)

}

#' @noRd
# Cudeck and Browne (1992) gradient function
# Updated 20.11.2026
CB_error_gradient <- function(vector, metric, target, preparation, df, n = population_size)
{

  # Obtain fit
  reference <- preparation$reference
  nodes <- ncol(reference$R)
  fit <- CB_error(as.vector(CB_project(vector, preparation)), preparation, df, n)

  # Cost is flat when invalid
  if(!fit$valid){
    return(numeric(length(vector)))
  }

  # Inverse of the population matrix
  R_population_inverse <- chol2inv(chol(fit$R_population))

  # Parameter e enters Sigma_pop at (i,j) and (j,i), so
  #   d chi-square / de = 2 * n * (R^-1 - Sigma_pop^-1)_ij
  #   d baseline / de = -2 * n * (Sigma_pop^-1)_ij
  chi_gradient <- 2 * n * (reference$inverse - R_population_inverse)
  baseline_gradient <- -2 * n * R_population_inverse

  # Chi-square above degrees of freedom
  model_excess <- fit$chi_square - df
  baseline_excess <- fit$baseline_chi_square - fit$baseline_df

  # Gradient of the metric
  metric_gradient <- switch(
    metric,
    "CFI" = {

      # Flat when chi-square is below degrees of freedom (CFI = 1)
      # or when it exceeds the baseline (CFI = 0)
      if(model_excess <= 0 || model_excess >= baseline_excess){
        matrix(0, nrow = nodes, ncol = nodes)
      }else{
        -(chi_gradient * baseline_excess - model_excess * baseline_gradient) /
          baseline_excess^2
      }

    },
    "TLI" = {

      # B = baseline chi-square / df, M = model chi-square / df
      B <- fit$baseline_chi_square / fit$baseline_df
      M <- fit$chi_square / df
      (-(chi_gradient / df) * (B - 1) + (baseline_gradient / fit$baseline_df) * (M - 1)) /
        (B - 1)^2

    },
    "RMSEA" = {

      # Flat when chi-square is below degrees of freedom (RMSEA = 0)
      if(model_excess <= 0){
        matrix(0, nrow = nodes, ncol = nodes)
      }else{
        chi_gradient / (2 * df * n * fit$RMSEA)
      }

    },
    "SRMR" = {

      # Flat at zero residual
      if(fit$SRMR == 0){
        matrix(0, nrow = nodes, ncol = nodes)
      }else{
        (fit$R_population - reference$R) /
          (fit$SRMR * nodes * (nodes + 1) / 2)
      }

    }
  )

  # d cost / de = 2 * (metric - target) * d metric / de
  # (the projection is symmetric so the gradient is projected)
  return(
    as.vector(
      CB_project(
        2 * (fit[[metric]] - target) * metric_gradient[preparation$cells],
        preparation
      )
    )
  )

}

#' @noRd
# Cudeck and Browne (1992) optimization
# Updated 20.11.2026
CB_optimize <- function(
    metric, target, preparation, parts,
    n = population_size, tolerance = 0.01
)
{

  # Obtain reference and degrees of freedom
  reference <- preparation$reference
  nodes <- parts$nodes
  df <- parts$df
  cells <- length(preparation$cells)

  # Check for room for error
  if(preparation$dimension < 1){
    return(list(valid = FALSE))
  }

  # Random direction over error cells
  # (random signs with magnitudes in [0.5, 1] before projection)
  direction <- as.vector(
    CB_project(
      shuffle_replace(c(-1, 1), cells) * runif_xoshiro(cells, min = 0.5, max = 1),
      preparation
    )
  )

  # Check for degenerate direction
  if(max(abs(direction)) < 1e-10){
    return(list(valid = FALSE))
  }

  # Error matrix of the direction
  error_matrix <- CB_population(direction, preparation) - reference$R

  # Eigenvalues of U^-T D U^-1 (closed-form fit along the direction)
  eigenvalues <- eigen(
    crossprod(preparation$inverse_chol, error_matrix %*% preparation$inverse_chol),
    symmetric = TRUE, only.values = TRUE
  )$values
  unit_SRMR <- sqrt(mean(error_matrix[lower.tri(error_matrix, diag = TRUE)]^2))

  # Largest scale that keeps the errors within [-1, 1]
  maximum_scale <- 1 / max(abs(direction))

  # Largest scale that keeps the population matrix positive definite
  # (R + scale * D is positive definite when 1 + scale * lambda_min > 0)
  smallest <- min(eigenvalues)
  if(smallest < 0){
    maximum_scale <- min(maximum_scale, 0.999 * -1 / smallest)
  }

  # Scale the direction to reach the target metric (closed form)
  scaled <- solve_scalar(
    FUN = function(scale){
      CB_direction_fit(scale, eigenvalues, unit_SRMR, preparation, df, n)
    },
    lower = 0, upper = maximum_scale,
    metric = metric, target = target, tolerance = 1e-14
  )

  # Exact fit at the solution
  if(isTRUE(scaled$valid)){
    scaled <- c(
      CB_error(direction * scaled$parameter, preparation, df, n),
      list(parameter = scaled$parameter)
    )
  }

  # Otherwise, optimize error directly
  if(!isTRUE(scaled$valid) || abs(scaled[[metric]] - target) > tolerance){

    # Start from the closest valid scaled direction (zero is flat)
    start_scale <- swiftelse(
      isTRUE(scaled$valid), max(scaled$parameter, 0.01 * maximum_scale),
      0.01 * maximum_scale
    )

    optimization <- try(
      optim(
        par = direction * start_scale,
        fn = CB_error_cost, gr = CB_error_gradient,
        metric = metric, target = target, preparation = preparation,
        df = df, n = n, method = "L-BFGS-B"
      ), silent = TRUE
    )

    # Check for errors
    if(!is(optimization, "try-error")){

      # Obtain fit
      projected <- as.vector(CB_project(optimization$par, preparation))
      fit <- CB_error(projected, preparation, df, n)

      # Keep the solution closest to the target
      if(
        fit$valid && (
          !isTRUE(scaled$valid) ||
          abs(fit[[metric]] - target) < abs(scaled[[metric]] - target)
        )
      ){
        scaled <- c(fit, list(parameter = NA_real_))
      }

    }

  }

  # Check for no solution
  if(!isTRUE(scaled$valid)){
    return(list(valid = FALSE))
  }

  # True model at the population matrix (stationary at the true parameters)
  known <- try(known_cfa(scaled$R_population, parts), silent = TRUE)
  if(is(known, "try-error") || !isTRUE(known$converged)){
    return(list(valid = FALSE))
  }

  # Confirm that the true parameters are the minimizer: a refit from
  # perturbed (deterministic) starting values must not find a lower
  # discrepancy (otherwise the true model is not the best-fitting model
  # and the direction is rejected)
  perturbation <- 1 + 0.10 * sin(seq_len(parts$parameters))
  confirmed <- try(
    known_cfa(scaled$R_population, parts, start = parts$start * perturbation),
    silent = TRUE
  )
  if(
    !is(confirmed, "try-error") && isTRUE(confirmed$converged) &&
    confirmed$discrepancy < scaled$discrepancy - 1e-10 * max(1, scaled$discrepancy)
  ){
    return(list(valid = FALSE))
  }

  # Return exact fit with the true model
  scaled$model <- known
  return(scaled)

}

# Tucker, Koopman, and Linn (1969) ----

#' @noRd
# Tucker, Koopman, and Linn (1969) population correlation matrix
# (minor factors contribute `weight` of the variance; diagonal remains one)
# Updated 20.11.2026
TKL_population <- function(weight, minor, R)
{
  return((1 - weight) * R + weight * minor)
}

#' @noRd
# Tucker, Koopman, and Linn (1969) optimization
# Updated 20.11.2026
TKL_optimize <- function(
    metric, target, R, parts,
    minor_factors, decay, n = population_size
)
{

  # Obtain nodes
  nodes <- parts$nodes
  minor_factors <- round(minor_factors)

  # Random loadings on minor factors
  loadings <- matrix(
    rnorm_ziggurat(nodes * minor_factors),
    nrow = nodes, ncol = minor_factors
  )

  # Geometric decay of loadings across minor factors
  loadings <- loadings * rep((1 - decay)^(seq_len(minor_factors) - 1), each = nodes)

  # Scale each variable's loadings to unit length
  loadings <- loadings / sqrt(rowSums(loadings^2))

  # Correlations from minor factors
  minor <- tcrossprod(loadings)

  # Largest weight that keeps the errors within [-1, 1]
  # (error is linear in the weight: weight * (minor - R))
  maximum_weight <- min(
    0.999, 1 / max(abs((minor - R)[upper.tri(R)]))
  )

  # Return solution
  return(
    solve_refit(
      function(weight){TKL_population(weight, minor, R)},
      parts, lower = 0, upper = maximum_weight,
      metric = metric, target = target, n = n
    )
  )

}

# Wu and Browne (2015) ----

#' @noRd
# Wu and Browne (2015) population correlation matrix
# Updated 20.11.2026
WB_population <- function(log_dispersion, draws, R)
{

  # Obtain precision parameter (m = 1 / v)
  nodes <- ncol(R)
  precision <- exp(-log_dispersion)

  # Bartlett decomposition (fixed random draws, continuous precision)
  bartlett <- draws$normals
  diag(bartlett) <- sqrt(qchisq(draws$uniforms, df = precision - seq_len(nodes) + 1))

  # Wishart draw on the precision matrix, centered on the original precision
  # matrix (E[W] = R^-1), which is an inverted Wishart draw for the covariance
  #   Sigma = W^-1, W ~ Wishart(m, R^-1 / m)
  factor_inverse <- try(
    forwardsolve(draws$cholesky %*% bartlett, diag(nodes)), silent = TRUE
  )
  if(is(factor_inverse, "try-error") || any(!is.finite(factor_inverse))){
    return(NULL)
  }
  R_draw <- cov2cor(crossprod(factor_inverse))

  # Check that the errors are within [-1, 1]
  # (the draw is kept as is so that it remains a pure inverted Wishart draw;
  # draws outside the bounds are treated as invalid)
  if(max(abs(R_draw - R)) > 1){
    return(NULL)
  }

  # Return population correlation matrix
  return(R_draw)

}

#' @noRd
# Wu and Browne (2015) optimization
# Updated 20.11.2026
WB_optimize <- function(metric, target, R, parts, n = population_size)
{

  # Obtain nodes
  nodes <- parts$nodes

  # Random draws (held fixed while the dispersion is optimized)
  normals <- matrix(0, nrow = nodes, ncol = nodes)
  normals[lower.tri(normals)] <- rnorm_ziggurat(nodes * (nodes - 1) / 2)
  draws <- list(
    normals = normals, uniforms = runif_xoshiro(nodes),
    cholesky = t(chol(chol2inv(chol(R))))
  )

  # Return solution
  # (dispersion v = 1 / m, with m > nodes - 1; v is approximately RMSEA^2)
  return(
    solve_refit(
      function(log_dispersion){WB_population(log_dispersion, draws, R)},
      parts, lower = log(1e-06), upper = log(1 / nodes),
      metric = metric, target = target, n = n, spacing = "linear"
    )
  )

}

# Hu and Bentler (1998) ----

#' @noRd
# Hu and Bentler (1998) population correlation matrix
# Omitted cross-loadings of magnitude `size` are added to the population
# loadings; uniquenesses are reduced so that the diagonal remains one
# Updated 20.11.2026
HB_population <- function(size, omissions)
{

  # Population loadings with omitted cross-loadings
  loadings <- omissions$loadings + size * omissions$omitted

  # Common part
  common <- loadings %*% tcrossprod(omissions$phi, loadings)

  # Uniquenesses must remain positive
  uniquenesses <- 1 - diag(common)
  if(any(uniquenesses <= 1e-06)){
    return(NULL)
  }

  # Population correlation matrix (correlated residuals remain)
  residuals <- omissions$residuals
  diag(residuals) <- uniquenesses
  R_population <- common + residuals

  # Check positive definiteness
  if(is(try(chol(R_population), silent = TRUE), "try-error")){
    return(NULL)
  }

  # Return population correlation matrix
  return(R_population)

}

#' @noRd
# Hu and Bentler (1998) optimization
# (underparameterized misspecification: omitted cross-loadings)
# Updated 20.11.2026
HB_optimize <- function(
    metric, target, model, R, parts, omitted, n = population_size
)
{

  # Zero loadings in the true model
  zero_cells <- which(model$loadings == 0)

  # Check for omittable loadings
  if(length(zero_cells) == 0){
    return(list(valid = FALSE))
  }

  # Randomly select omitted cross-loadings (at least one)
  omitted_cells <- shuffle(
    zero_cells, size = max(1, round(omitted * length(zero_cells)))
  )

  # Random signs for omitted cross-loadings
  omitted_loadings <- matrix(0, nrow = nrow(model$loadings), ncol = ncol(model$loadings))
  omitted_loadings[omitted_cells] <- swiftelse(
    runif_xoshiro(length(omitted_cells)) < 0.50, -1, 1
  )

  # Correlated residuals of the true model (diagonal is replaced)
  residuals <- model$psi
  diag(residuals) <- 0

  # Population components
  omissions <- list(
    loadings = model$loadings, phi = model$phi,
    omitted = omitted_loadings, residuals = residuals
  )

  # Population correlation matrix (errors within [-1, 1])
  population <- function(size){

    # Obtain population correlation matrix
    R_population <- HB_population(size, omissions)

    # Return if allowed
    if(is.null(R_population) || max(abs(R_population - R)) > 1){
      return(NULL)
    }else{
      return(R_population)
    }

  }

  # Find largest size (magnitude of omitted cross-loadings)
  maximum_size <- 1
  if(is.null(population(maximum_size))){

    # Bisection (smallest size is always allowed)
    lower <- 0; upper <- maximum_size
    for(i in seq_len(50)){

      # Update bounds
      middle <- (lower + upper) / 2
      if(!is.null(population(middle))){
        lower <- middle
      }else{
        upper <- middle
      }

    }

    # Stay away from the boundary
    maximum_size <- lower * 0.999

  }

  # Return solution
  return(
    solve_refit(
      population, parts, lower = 0, upper = maximum_size,
      metric = metric, target = target, n = n
    )
  )

}

# Data ----

#' @noRd
# Generate data from a population correlation matrix
# (same categorization and skew handling as `simulate_factors`)
# Updated 20.11.2026
simulate_error_data <- function(R, parameters, sample_size)
{
  return(
    generate_lf_data(
      R, sample_size, parameters$categories,
      parameters$categorical_limit, parameters$skew
    )$data
  )
}
