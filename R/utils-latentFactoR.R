#%%%%%%%%%%%%%%%%%%%%%%%%%%
# add_population_error ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Grid search for Cudeck function
# Updated 17.09.2022
grid_search <- function(delta, G) {
  
  n <- 1000
  x <- seq(-1e3, 1e3, length.out = n)
  y <- vector(length = n)
  for(i in 1:n) y[i] <- root_ml(x[i], delta, G)
  
  index <- which.min(y)
  x <- seq(x[index-1], x[index+1], length.out = n)
  for(i in 1:n) y[i] <- root_ml(x[i], delta, G)
  
  index <- which.min(y)
  x <- seq(x[index-1], x[index+1], length.out = n)
  for(i in 1:n) y[i] <- root_ml(x[i], delta, G)
  
  return(x[which.min(y)])
  
}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Maximum likelhood for grid search in Cudeck function
# Updated 17.09.2022
root_ml <- function(x, delta, G) {
  
  I <- diag(nrow(G))
  f <- x*sum(diag(G)) - log(det(I + x*G)) - delta
  
  return(f^2)
  
}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Maximum likelhood for grid search in Cudeck function
# Updated 17.09.2022
groot_ml <- function(x, delta, G) {
  
  I <- diag(nrow(G))
  f <- x*sum(diag(G)) - log(det(I + x*G)) - delta
  g <- sum(diag(G)) - sum(diag(solve(I + x*G) %*% G))
  g <- 2*g*f
  
  return(g)
  
}

#' @noRd
#' @importFrom stats optim
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Optimization function
# Updated 17.09.2022
opt <- function(x, delta, G) {
  
  # x <- stats::runif(1, 0, 1)
  # det(diag(nrow(G)) + x*G)
  root <- optim(x, fn = root_ml, gr = groot_ml, method = "L-BFGS-B",
                       lower = -Inf, upper = Inf, G = G, delta = delta)
  k <- root$par
  
  return(k)
  
}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Error catch for optimization function
# Updated 17.09.2022
opt_error <- function(x, delta, G) {
  
  x <- tryCatch({opt(x, delta, G)}, error = return(x))
  return(x)
  
}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Derivative function for Rhat methods
# Updated 17.09.2022
dxt <- function(X) {
  
  # derivative wrt transpose (just a permutation matrix)
  
  p <- nrow(X)
  q <- ncol(X)
  pq <- p*q
  
  res <- array(0, dim = c(pq, pq))
  null <- matrix(0, p, q)
  
  for(i in 1:pq) {
    temp <- null
    temp[i] <- 1
    res[, i] <- c(t(temp))
  }
  
  return(res)
  
}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# LRhat function for Cudeck method
# Updated 17.09.2022
gLRhat <- function(Lambda, Phi) {
  
  # derivative of Lambda wrt Rhat
  
  p <- nrow(Lambda)
  g1 <- (Lambda %*% Phi) %x% diag(p)
  g21 <- diag(p) %x% (Lambda %*% Phi)
  g2 <- g21 %*% dxt(Lambda)
  g <- g1 + g2
  
  return(g)
  
}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# PRhat function for Cudeck method
# Updated 17.09.2022
gPRhat <- function(Lambda, Phi) {
  
  g1 <- Lambda %x% Lambda
  g2 <- g1 %*% dxt(Phi)
  g <- g1 + g2
  g <- g[, which(lower.tri(Phi))]
  
  return(g)
  
}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Rhat function for Cudeck method
# Updated 17.09.2022
guRhat <- function(p) {
  
  gu <- matrix(0, p*p, p)
  
  for(i in 1:p) {
    
    index <- (i-1)*p + i
    gu[index, i] <- 1
    
  }
  
  return(gu)
  
}

#' @noRd
#' @importFrom stats lm
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Adds population error using Cudeck method to generated data
# Updated 17.09.2022
cudeck <- function(R, lambda, Phi, uniquenesses,
                   fit = "rmsr", misfit = "close",
                   method = "ols", confirmatory = TRUE) {
  
  # Method of Cudeck and Browne (1992):
  
  p <- nrow(R)
  q <- ncol(lambda)
  
  tdiag <- TRUE
  
  if(confirmatory) {
    
    # Select the columns corresponding to estimated loadings (only works when
    # not estimating correlations)
    # if(!correlation) dS_dL <- dS_dL[, which(lambda > 0)]
    npars <- sum(lambda > 0) + p + sum(abs(Phi[lower.tri(Phi)]) > 0)
    dS_du <- guRhat(p)
    dS_dL <- gLRhat(lambda, Phi)[, which(lambda != 0)]
    dS_dP <- gPRhat(lambda, Phi)[, which(Phi[lower.tri(Phi)] != 0)]
    gS <- cbind(dS_dL, dS_dP, dS_du) # matrix of derivatives wrt the correlation model
    
  } else {
    
    # dS_dP <- gPRhat(lambda, Phi)
    # gS <- cbind(dS_dL, dS_dP, dS_du)
    npars <- p*q + p - 0.5*q*(q-1) # Number of model parameters
    dS_du <- guRhat(p)
    dS_dL <- gLRhat(lambda, Phi)
    # dS_dP <- gPRhat(lambda, Phi)
    # gS <- cbind(dS_dL, dS_dP, dS_du) # matrix of derivatives wrt the correlation model
    gS <- cbind(dS_dL, dS_du) # matrix of derivatives wrt the correlation model
    
  }
  
  df <- p*(p+1)/2 - npars # Degrees of freedom
  
  if(method == "minres" || method == "ols") {
    
    B <- -2*gS[lower.tri(R, diag = tdiag), ]
    # B <- -2*lambda %*%
    
  } else if(method == "ml") {
    
    # K <- transition(p)
    # MP_inv <- solve(t(K) %*% K) %*% t(K)
    # D <- MP_inv %*% t(MP_inv)
    indexes <- vector(length = p)
    indexes[1] <- 1
    for(i in 2:p) {
      increment <- i
      indexes[i] <- indexes[i-1]+increment
    }
    D <- matrix(0, p*(p+1)/2, p*(p+1)/2)
    diag(D) <- 2
    diag(D)[indexes] <- 1
    R_inv <- solve(R)
    vecs <- apply(gS, 2, FUN = function(x) -t(R_inv %*% matrix(x, p, p) %*% R_inv))
    B <- t(vecs[which(upper.tri(R, diag = tdiag)), ]) %*% D
    B <- t(B)
    
  }
  
  # BtB <- t(B) %*% B
  m <- p+1
  U <- replicate(p, runif(m, 0, 1))
  A1 <- t(U) %*% U
  sq <- diag(1/sqrt(diag(A1)))
  A2 <- sq %*% A1 %*% sq
  diag_u <- diag(sqrt(uniquenesses))
  y <- diag_u %*% A2 %*% diag_u
  y <- y[lower.tri(y, diag = tdiag)]
  # y <- A2[lower.tri(A2, diag = tdiag)]
  # y <- stats::runif(p*(p+1)/2, 0, 1)
  # e <- qr.Q(qr(cbind(B, y)))[, ncol(B)+1]
  # v <- MASS::ginv(BtB) %*% t(B) %*% y
  # e <- y - B %*% v # equation 7
  # B.qr <- qr(B)
  # e <- qr.resid(B.qr, y)
  e <- unname(lm(y ~ B, model = FALSE, qr = TRUE)$residuals)
  
  if(fit == "rmsr") {
    if(misfit == "close") {
      r2 <- mean(1-uniquenesses)
      misfit <- 0.05*r2
    } else if(misfit == "acceptable") {
      r2 <- mean(1-uniquenesses)
      misfit <- 0.10*r2
    }
    delta <- misfit^2*0.5*p*(p-1)
    # delta <- (1-misfit2)*(0.5*(sum(R_error^2) - p))
  } else if(fit == "cfi") {
    null_f <- 0.5*(sum(R^2) - p)
    delta <- (1-misfit)*null_f
  } else if(fit == "rmsea") {
    delta <- misfit^2 * df
  } else if(fit == "raw") {
    delta <- misfit
  }
  
  if(method == "minres" || method == "ols") {
    
    E <- matrix(0, p, p)
    E[lower.tri(E, diag = tdiag)] <- e
    E <- t(E) + E
    diag(E) <- 0
    k <- sqrt(2*delta/sum(E*E))
    E <- k*E
    
  } else if(method == "ml") {
    
    E <- matrix(0, p, p)
    E[upper.tri(R, diag = tdiag)] <- e
    E <- t(E) + E
    diag(E) <- 0
    constant <- 1e-04 / sqrt(mean(E*E))
    E <- constant*E # Fix this to avoid NAs
    R_inv <- solve(R)
    G <- R_inv %*% E
    x <- suppressWarnings(grid_search(delta, G))
    k <- opt(x, delta, G)
    # limits <- c(-1e05, 1e05)
    # k <- GSS(delta, G, limits)
    # k <- grad_descend(delta, G)
    
    E <- k*E
    
  }
  
  R_error <- R + E
  
  # check for positiveness:
  minimum_eigval <- min(eigen(R_error, symmetric = TRUE, only.values = TRUE)$values)
  if(minimum_eigval <= 0) warning("The matrix was not positive-definite. The amount of error may be too big.")
  
  return(list(R_error = R_error, fit = fit, delta = delta, misfit = misfit))
  
}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Minimum residual function for CFA in Yuan method
# Updated 17.09.2022
f_minres <- function(x, S, ldetS, q, indexes_lambda, indexes_phi) {
  
  p <- nrow(S)
  lambda_parameters <- length(indexes_lambda)
  Lambda <- matrix(0, p, q)
  Lambda[indexes_lambda] <- x[1:lambda_parameters]
  Phi <- matrix(0, q, q)
  Phi[indexes_phi] <- x[-(1:lambda_parameters)]
  Phi <- t(Phi) + Phi
  diag(Phi) <- 1
  Rhat <- Lambda %*% Phi %*% t(Lambda)
  diag(Rhat) <- 1
  res <- S - Rhat
  f <- 0.5*sum(res*res)
  
  return(f)
  
}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Minimum residual function for CFA in Yuan method
# Updated 17.09.2022
g_minres <- function(x, S, ldetS, q, indexes_lambda, indexes_phi) {
  
  p <- nrow(S)
  lambda_parameters <- length(indexes_lambda)
  Lambda <- matrix(0, p, q)
  Lambda[indexes_lambda] <- x[1:lambda_parameters]
  Phi <- matrix(0, q, q)
  Phi[indexes_phi] <- x[-(1:lambda_parameters)]
  Phi <- t(Phi) + Phi
  diag(Phi) <- 1
  Rhat <- Lambda %*% Phi %*% t(Lambda)
  diag(Rhat) <- 1
  res <- S - Rhat
  
  g1 <- (res %*% Lambda %*% Phi)[indexes_lambda]
  g2 <- (t(Lambda) %*% res %*% Lambda)[indexes_phi]
  g <- -2*c(g1, g2)
  
  return(g)
  
}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Maximum likelihood function for CFA in Yuan method
# Updated 17.09.2022
f_ml <- function(x, S, ldetS, q, indexes_lambda, indexes_phi) {
  
  p <- nrow(S)
  lambda_p <- length(indexes_lambda)
  Lambda <- matrix(0, p, q)
  Lambda[indexes_lambda] <- x[1:lambda_p]
  phi_p <- length(indexes_phi)
  Phi <- matrix(0, q, q)
  Phi[indexes_phi] <- x[(lambda_p+1):(lambda_p + phi_p)]
  Phi <- t(Phi) + Phi
  diag(Phi) <- 1
  u <- x[-(1:(lambda_p + phi_p))]
  Rhat <- Lambda %*% Phi %*% t(Lambda) + diag(u)
  f <- log(det(Rhat)) - ldetS + sum(S*solve(Rhat)) - p
  
  return(f)
  
}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Maximum likelihood function for CFA in Yuan method
# Updated 17.09.2022
g_ml <- function(x, S, ldetS, q, indexes_lambda, indexes_phi) {
  
  p <- nrow(S)
  lambda_p <- length(indexes_lambda)
  Lambda <- matrix(0, p, q)
  Lambda[indexes_lambda] <- x[1:lambda_p]
  phi_p <- length(indexes_phi)
  Phi <- matrix(0, q, q)
  Phi[indexes_phi] <- x[(lambda_p+1):(lambda_p + phi_p)]
  Phi <- t(Phi) + Phi
  diag(Phi) <- 1
  u <- x[-(1:(lambda_p + phi_p))]
  
  Rhat <- Lambda %*% Phi %*% t(Lambda) + diag(u)
  Rhat_inv <- solve(Rhat)
  Ri_res_Ri <- 2*Rhat_inv %*% (Rhat - S) %*% Rhat_inv
  
  # Joreskog (page 10; 1965) Testing a simple structure in factor analysis
  g <- c(c(Ri_res_Ri %*% Lambda %*% Phi)[indexes_lambda],
         c(t(Lambda) %*% Ri_res_Ri %*% Lambda)[indexes_phi],
         diag(Ri_res_Ri)*0.5) # omitting the sample size constant
  
  return(g)
  
}

#' @noRd
#' @importFrom stats runif nlminb
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# CFA function for Yuan method
# Updated 17.09.2022
CFA <- function(S, target, targetphi, method = "minres") {
  
  p <- nrow(target)
  q <- ncol(target)
  indexes_lambda <- which(target != 0)
  indexes_phi <- which(targetphi != 0 & lower.tri(targetphi))
  lambda_p <- length(indexes_lambda)
  phi_p <- length(indexes_phi)
  init <- 1/diag(solve(S))
  
  if(method == "minres") {
    
    x <- c(runif(lambda_p), rep(0, phi_p))
    ldetS <- NULL
    f <- f_minres
    g <- g_minres
    lower <- c(rep(-Inf, lambda_p), rep(-1, phi_p))
    upper <- c(rep(Inf, lambda_p), rep(1, phi_p))
    
  } else if(method == "ml") {
    
    x <- c(runif(lambda_p), rep(0, phi_p), init)
    ldetS <- log(det(S))
    f <- f_ml
    g <- g_ml
    lower <- c(rep(-Inf, lambda_p), rep(-1, phi_p), rep(0, p))
    upper <- c(rep(Inf, lambda_p), rep(1, phi_p + p))
    
  }
  
  cfa <- nlminb(x, objective = f, gradient = g,
                       lower = lower, upper = upper,
                       S = S, ldetS = ldetS, q = q, indexes_lambda = indexes_lambda,
                       indexes_phi = indexes_phi,
                       control = list(iter.max = 1e4, eval.max = 1e4))
  
  lambda_hat <- matrix(0, p, q)
  lambda_hat[indexes_lambda] <- cfa$par[1:lambda_p]
  phi_hat <- matrix(0, q, q)
  
  if(method == "minres") {
    
    phi_hat[indexes_phi] <- cfa$par[-(1:lambda_p)]
    phi_hat <- phi_hat + t(phi_hat)
    diag(phi_hat) <- 1
    
  } else if(method == "ml") {
    
    phi_hat[indexes_phi] <- cfa$par[(lambda_p+1):(lambda_p + phi_p)]
    phi_hat <- t(phi_hat) + phi_hat
    diag(phi_hat) <- 1
    u <- cfa$par[-(1:(lambda_p + phi_p))]
    
  }
  
  S_hat <- lambda_hat %*% phi_hat %*% t(lambda_hat)
  uniquenesses_hat <- 1 - diag(S_hat)
  diag(S_hat) <- 1
  residuals <- S - S_hat
  
  df <- p*(p+1)/2 - (length(indexes_lambda) + length(indexes_phi) + p)
  
  results <- list(f = cfa$objective, convergence = cfa$convergence,
                  iterations = cfa$iterations, df = df,
                  lambda = lambda_hat, phi = phi_hat,
                  uniquenesses = uniquenesses_hat,
                  model = S_hat, residuals = residuals)
  
  return(results)
  
}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Adds population error using Yuan method to generated data
# Updated 17.09.2022
yuan <- function(
    R, lambda, Phi, uniquenesses,
    fit = "rmsr", misfit = "close",
    method = "minres", confirmatory = TRUE
)
{
  
  p <- nrow(R)
  q <- ncol(lambda)
  
  tdiag <- TRUE
  
  if(confirmatory) {
    
    npars <- sum(lambda > 0) + p + sum(abs(Phi[lower.tri(Phi)]) > 0)
    
  } else {
    
    npars <- p*q + p - 0.5*q*(q-1) # Number of model parameters
    
  }
  
  df <- p*(p+1)/2 - npars # Degrees of freedom
  
  if(fit == "rmsr") {
    if(misfit == "close") {
      r2 <- mean(1-uniquenesses)
      misfit <- 0.05*r2
    } else if(misfit == "acceptable") {
      r2 <- mean(1-uniquenesses)
      misfit <- 0.10*r2
    }
    delta <- misfit^2*0.5*p*(p-1)
    # delta <- (1-misfit2)*(0.5*(sum(R_error^2) - p))
  } else if(fit == "cfi") {
    null_f <- 0.5*(sum(R^2) - p)
    delta <- (1-misfit)*null_f
  } else if(fit == "rmsea") {
    delta <- misfit^2 * df
  } else if(fit == "raw") {
    delta <- misfit
  }
  
  p <- nrow(R)
  L <- lambda + 1e-06
  R1 <- R
  R <- L %*% Phi %*% t(L); diag(R) <- 1
  target <- ifelse(lambda != 0, 1, 0)
  targetphi <- ifelse(Phi != 0, 1, 0)
  fit <- CFA(R, target, targetphi, method = method)
  Phat <- fit$model
  
  # from delta to tau:
  E <- R - Phat
  
  if(method == "minres" || method == "ols") {
    
    tau <- sqrt(2*delta/sum(E*E))
    R_error <- Phat + tau*E
    
  } else if(method == "ml") {
    
    R_inv <- solve(R1)
    constant <- 1e-04 / sqrt(mean(E*E))
    E <- constant*E # Fix this to avoid NAs
    G <- R_inv %*% E
    
    tau <- suppressWarnings(grid_search(delta, G))
    tau <- opt_error(tau, delta, G)
    # limits <- c(-1e3, 1e3)
    # tau <- GSS(delta, G, limits)
    # tau <- grad_descend(delta, G)
    
    R_error <- Phat + tau*E
    
  }
  
  # check for positiveness:
  minimum_eigval <- min(eigen(R_error, symmetric = TRUE, only.values = TRUE)$values)
  if(minimum_eigval <= 0) warning("The matrix was not positive-definite. The amount of error may be too big.")
  
  return(list(R_error = R_error, fit = fit, delta = delta, misfit = misfit))
  
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%
# add_local_dependence ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Adds correlated residuals to generated data
# Updated 23.09.2022
correlate_residuals <- function(
    lf_object,
    proportion_LD, allow_multiple = FALSE,
    add_residuals, add_residuals_range
)
{
  
  # Obtain parameters
  parameters <- lf_object$parameters
  
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
  
  # Obtain number of local dependencies
  variables_LD <- round(proportion_LD * variables)
  
  # If variables cannot have multiple local dependencies,
  # then number of local dependencies needs to be cut in half (per factor)
  if(!isTRUE(allow_multiple)){
    variables_LD <- floor(variables_LD / 2) 
  }
  
  # Check for add residual range
  if(!is.null(add_residuals_range)){
    type_error(add_residuals_range, "numeric") # object type error
    length_error(add_residuals_range, 2) # object length error
    range_error(add_residuals_range, c(0, 1)) # object range error
    add_residuals <- runif(
      sum(variables_LD),
      min = min(add_residuals_range),
      max = max(add_residuals_range)
    )
  }
  
  # Ensure appropriate types
  type_error(add_residuals, "numeric");
  
  # Ensure appropriate lengths
  length_error(add_residuals, c(1, parameters$factors, sum(variables_LD)));
  
  # Ensure appropriate ranges
  range_error(add_residuals, c(0, 1));
  
  # Set start and end points for variables
  end_variables <- cumsum(variables)
  start_variables <- end_variables + 1 - variables
  
  # Initialize checks
  check_eigenvalues <- TRUE
  
  # Run through loop
  while(isTRUE(check_eigenvalues)){
    
    # Initialize correlated residual matrix
    correlated_residuals <- matrix(
      0, nrow = 0, ncol = 2
    )
    
    # Loop through factors and add local dependence
    for(f in 1:factors){
      
      # Item rows
      item_rows <- sample(
        start_variables[f]:end_variables[f],
        variables_LD[f],
        replace = allow_multiple
      )
      
      # Set remaining variables
      if(isTRUE(allow_multiple)){
        
        # Do not remove variables
        remaining_variables <- start_variables[f]:end_variables[f]
        
      }else{
        
        # Remove already included variables
        remaining_variables <- setdiff(
          start_variables[f]:end_variables[f], item_rows
        )
        
      }
      
      # Item columns
      item_columns <- sample(
        remaining_variables,
        variables_LD[f],
        replace = allow_multiple
      )
      
      # Bind to correlated residual matrix
      correlated_residuals <- rbind(
        correlated_residuals,
        cbind(item_rows, item_columns)
      )
      
      # Obtain duplicate rows
      duplicate_rows <- match_row(correlated_residuals)
      
      # Replace until there are no duplicate rows
      while(any(duplicate_rows)){
        
        # Replace second column with new variable
        correlated_residuals[duplicate_rows, 2] <- sample(
          remaining_variables,
          length(duplicate_rows),
          replace = allow_multiple
        )
        
        # Re-check for duplicate rows
        duplicate_rows <- match_row(correlated_residuals)
        
      }
      
    }
    
    # Obtain amount residual to add
    if(length(add_residuals) == length(variables_LD)){
      
      # Loop through correlated_residuals
      for(i in 1:nrow(correlated_residuals)){
        
        # Add residuals to correlation matrix
        original_correlation[
          correlated_residuals[i,1],
          correlated_residuals[i,2]
        ] <- add_residuals[i]
        
        # Ensure symmetric
        original_correlation[
          correlated_residuals[i,2],
          correlated_residuals[i,1]
        ] <- original_correlation[
          correlated_residuals[i,1],
          correlated_residuals[i,2]
        ]
        
      }
      
      
    }else{
      
      # Check if add residuals length equals 1 or number of factors
      if(length(add_residuals) != sum(variables_LD)){
        
        # If only one value
        if(length(add_residuals) == 1){
          add_residuals <- rep(add_residuals, sum(variables_LD))
        }else if(length(add_residuals) == parameters$factors){# Length of factors
          
          # Loop through number of local dependence variables
          add_residuals <- unlist(lapply(1:parameters$factors, function(i){
            rep(add_residuals[i], variables_LD[i])
          }))
          
        }
        
      }
      
      # Loop through correlated residuals
      if(nrow(correlated_residuals) != 0){
        
        # Loop through correlated residuals
        for(i in 1:nrow(correlated_residuals)){
          
          # Compute random residual
          random_residual <- runif(
            1,
            min = add_residuals[i] - 0.05,
            max = add_residuals[i] + 0.05
          )
          
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
          
        }
        
      }
      
    }
    
    # Check eigenvalues
    check_eigenvalues <- any(eigen(population_correlation)$values <= 0)
    
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
  
  # Ensure appropriate type and length for categories
  type_error(variable_categories, "numeric")
  length_error(variable_categories, c(1, total_variables))
  
  # Identify categories to variables
  if(length(variable_categories) == 1){
    variable_categories <- rep(variable_categories, total_variables)
  }
  
  # Check for categories greater than categorical limit and not infinite
  if(any(variable_categories > categorical_limit & !is.infinite(variable_categories))){
    
    # Make variables with categories greater than 7 (or categorical_limit) continuous
    variable_categories[
      variable_categories > categorical_limit & !is.infinite(variable_categories)
    ] <- Inf
    
  }
  
  # Find categories
  if(any(variable_categories <= categorical_limit)){
    
    # Target columns to categorize
    columns <- which(variable_categories <= categorical_limit)
    
    # Set skew
    if(length(skew) != length(columns)){
      skew <- sample(skew, length(columns), replace = TRUE)
    }
    
    # Loop through columns
    for(i in columns){
      
      data[,i] <- categorize(
        data = data[,i],
        categories = variable_categories[i],
        skew_value = skew[i]
      )
      
    }
    
  }
  
  # Add column names to data
  colnames(data) <- paste0(
    "V", formatC(
      x = 1:total_variables,
      digits = floor(log10(total_variables)),
      flag = "0", format = "d"
    )
  )
  
  # Update correlated residuals
  correlated_residuals_df <- data.frame(
    V1 = correlated_residuals[,1],
    V2 = correlated_residuals[,2],
    added_residual = add_residuals
  )
  
  # Populate results
  results <- list(
    correlated_residuals = correlated_residuals_df,
    data = data,
    population_correlation = population_correlation,
    original_correlation = original_correlation,
    original_results = lf_object
  )
  
  # Return results
  return(results)
  
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

#%%%%%%%%%%%%%%%%%%%%%
# ERROR FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%

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
  actual_maximum <- round(max(input), 3)
  actual_minimum <- round(min(input), 3)
  
  # Check for maximum of input in expected range
  if(actual_maximum > expected_maximum){
    stop(
      paste(
        "Maximum of '", deparse(substitute(input)),
        "' (", actual_maximum,") does not match expected range(s). Range must be between: ",
        paste0("'", expected_ranges, "'", collapse = " and "),
        sep = ""
      )
    )
  }
  
  # Check for maximum of input in expected range
  if(actual_minimum < expected_minimum){
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

#%%%%%%%%%%%%%%%%%%%%%%%%%%
# GENERATION FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%

# Based on 
# Garrido, L. E., Abad, F. J., & Ponsoda, V. (2011).
# Performance of Velicer’s minimum average partial factor retention
# method with categorical variables.
# Educational and Psychological Measurement, 71(3), 551-570.
# https://doi.org/10.1177/0013164410389489
#
#' @noRd
# Generates skewed data for continuous data
# Updated 09.08.2022
skew_continuous <- function(
    skewness,
    data = NULL,
    sample_size = 1000000,
    tolerance = 0.00001
)
{
  
  # Original skewness
  original_skewness <- skewness
  
  # Generate data
  if(is.null(data)){
    data <- rnorm(sample_size)
  }
  
  # Kurtosis
  kurtosis <- 1
  
  # Skew data
  skew_data <- sinh(
    kurtosis * (asinh(data) + skewness) 
  )
  
  # Observed skew in data
  observed_skew <- psych::skew(skew_data)
  
  # Minimize difference
  while(abs(observed_skew - original_skewness) > tolerance){
    
    # Obtain difference
    difference <- observed_skew - original_skewness
    
    # Decrease skewness
    # Negative values increase
    # Positive values decrease
    skewness <- skewness - difference
    
    # Skew data
    skew_data <- sinh(
      kurtosis * (asinh(data) + skewness) 
    )
    
    # Observed skew in data
    observed_skew <- psych::skew(skew_data)
    
  }
  
  # Return skewed data
  return(skew_data)
  
}

# Based on 
# Garrido, L. E., Abad, F. J., & Ponsoda, V. (2011).
# Performance of Velicer’s minimum average partial factor retention
# method with categorical variables.
# Educational and Psychological Measurement, 71(3), 551-570.
# https://doi.org/10.1177/0013164410389489
#
#' @noRd
# Generates skew
# Updated 09.08.2022
skew_generator <- function(
    skewness, categories,
    reduction_factor = 0.75,
    sample_size = 1000000,
    initial_proportion = 0.50,
    tolerance = 0.00001
)
{
  
  # Initialize skew matrix
  skew_matrix <- matrix(
    0, nrow = categories, ncol = categories
  )
  
  # Initialize cases
  cases <- numeric(sample_size)
  
  # Loop through categories
  for(i in 2:categories){
    
    # Initialize (largest category) proportion
    proportion <- initial_proportion
    
    # Current proportion for rest of categories
    remaining_proportion <- 1 - proportion
    
    # Initialize categories allocations
    allocation <- 1
    allocation_1 <- 1
    
    # Loop through with reduction factor
    if(i > 2){
      for(j in 1:(i-2)){
        allocation_1 <- allocation_1 * reduction_factor
        allocation <- allocation + allocation_1
      }
    }
    
    # Divide remaining proportion by allocations
    divided_proportion <- remaining_proportion / allocation
    
    # Undefined objects
    propinf <- 1 / sample_size
    propsup <- initial_proportion
    E <- divided_proportion / reduction_factor
    
    # Loop through
    for(j in 1:i){
      cases[round(sample_size * propinf):round(sample_size * propsup)] <- j
      E <- E * reduction_factor
      propinf <- propsup
      propsup <- propinf + E
    }
    
    # Compute skewness
    skew_actual <- psych::skew(cases)
    
    # Limits
    limitsup <- 1
    limitsinf <- 0
    
    # Ensure skew within tolerance
    while(abs(skew_actual - skewness) > tolerance){
      
      # Skew greater than
      if(skew_actual < skewness){
        limitinf <- proportion
        proportion <- (proportion + limitsup) / 2
      }else{
        limitsup <- proportion
        proportion <- (proportion + limitsinf) / 2
      }
      
      # Update
      # Current proportion for rest of categories
      remaining_proportion <- 1 - proportion
      
      # Divide remaining proportion by allocations
      divided_proportion <- remaining_proportion / allocation
      
      # Undefined objects
      propinf <- 1 / sample_size
      propsup <- proportion
      E <- divided_proportion / reduction_factor
      
      # Loop through
      for(j in 1:i){
        cases[round(sample_size * propinf):round(sample_size * propsup)] <- j
        E <- E * reduction_factor
        propinf <- propsup
        propsup <- propinf + E
      }
      
      # Compute skewness
      skew_actual <- psych::skew(cases)
      
    }
    
    # Set E
    E <- divided_proportion / reduction_factor
    cumulative_probability <- proportion
    
    # Update matrix
    for(j in 1:i){
      
      skew_matrix[i,j] <- cumulative_probability
      E <- E * reduction_factor
      cumulative_probability <- cumulative_probability + E
      
    }
    
  }
  
  # Normal inverse
  norm_inv_matrix <- qnorm(skew_matrix)
  
  # Set infinite values to zero
  norm_inv_matrix[is.infinite(norm_inv_matrix)] <- 0
  
  # Category probability
  category_probability <- matrix(
    0, nrow = categories, ncol = categories
  )
  
  # Fill first column
  category_probability[,1] <- skew_matrix[,1]
  
  # Loop through
  for(i in 2:categories){
    category_probability[,i] <- skew_matrix[,i] - skew_matrix[,i-1]
  }
  
  # Make -1 = 0
  category_probability[category_probability == -1] <- 0
  
  # Return skew
  result <- list(
    skew_matrix = norm_inv_matrix,
    probability = category_probability
  )
  return(result)
  
}

#%%%%%%%%%%%%%%%%%%%%%%
# SYSTEM FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%

#' Error report
#' 
#' @description Gives necessary information for user reporting error
#' 
#' @param result Character.
#' The error from the result
#' 
#' @param SUB_FUN Character.
#' Sub-routine the error occurred in
#' 
#' @param FUN Character.
#' Main function the error occurred in
#' 
#' @return Error and message to send to GitHub
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @noRd
#' 
#' @importFrom utils packageVersion
#' 
# Error Report
# Updated 08.08.2022
error.report <- function(result, SUB_FUN, FUN)
{
  # Let user know that an error has occurred
  message(paste("\nAn error has occurred in the '", SUB_FUN, "' function of '", FUN, "':\n", sep =""))
  
  # Give them the error to send to you
  cat(paste(result))
  
  # Tell them where to send it
  message("\nPlease open a new issue on GitHub (bug report): https://github.com/hfgolino/EGAnet/issues/new/choose")
  
  # Give them information to fill out the issue
  OS <- as.character(Sys.info()["sysname"])
  OSversion <- paste(as.character(Sys.info()[c("release", "version")]), collapse = " ")
  Rversion <- paste(R.version$major, R.version$minor, sep = ".")
  latentFactoRversion <- paste(unlist(packageVersion("latentFactoR")), collapse = ".")
  
  # Let them know to provide this information
  message(paste("\nBe sure to provide the following information:\n"))
  
  # To reproduce
  message(styletext("To Reproduce:", defaults = "bold"))
  message(paste(" ", textsymbol("bullet"), " Function error occurred in: ", SUB_FUN, " function of ", FUN, sep = ""))
  
  # R, SemNetCleaner, and SemNetDictionaries
  message(styletext("\nR and latentFactoR versions:", defaults = "bold"))
  message(paste(" ", textsymbol("bullet"), " R version: ", Rversion, sep = ""))
  message(paste(" ", textsymbol("bullet"), " latentFactoR version: ", latentFactoRversion, sep = ""))
  
  # Desktop
  message(styletext("\nOperating System:", defaults = "bold"))
  message(paste(" ", textsymbol("bullet"), " OS: ", OS, sep = ""))
  message(paste(" ", textsymbol("bullet"), " Version: ", OSversion, sep = ""))
}

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