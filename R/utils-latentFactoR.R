#%%%%%%%%%%%%%%%%%%%%%%%%%%
# add_population_error ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Specify CFA model
# Updated 01.10.2022
model_CFA <- function(variables, loadings)
{
  
  # Initialize model
  model <- ""
  
  # Loop through factors
  for(i in 1:ncol(loadings)){
    
    # Append model
    if(i != ncol(loadings)){
      model <- paste0(
        model,
        "F", i, " =~ ",
        paste0(
          "V", which(loadings[,i] != 0),
          collapse = " + "
        ), " \n "
      )
    }else{
      model <- paste0(
        model,
        "F", i, " =~ ",
        paste0(
          "V", which(loadings[,i] != 0),
          collapse = " + "
        )
      )
    }
    
  }
  
  # Return model
  return(model)
  
}

#' @noRd
# Standardized Root Mean Resdiual
# Updated 28.09.2022
srmr <- function(population, error)
{
  
  # Obtain lower triangles
  population_lower <- population[lower.tri(population)]
  error_lower <- error[lower.tri(error)]
  
  # Compute SRMR
  SRMR <- sqrt(
    mean(
      (population_lower - error_lower)^2
    )
  )
  
  # Return SRMR
  return(SRMR)
  
}

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
  
  # Ensure matrix
  if(!is.matrix(g)){
    g <- matrix(g, ncol = 1)
  }
  
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
  
  # Ensure matrix
  if(!is.matrix(g)){
    g <- matrix(g, ncol = 1)
  }
  
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
# Rhat function for Cudeck method
# Updated 07.10.2022
# For Psi
gURhat <- function(p) {
  
  pcov <- p*(p+1)*0.5
  
  Psi <- diag(p)
  gPsi <- diag(p) %x% diag(p)
  gPsi <- gPsi + dxt(Psi) %*% gPsi
  gPsi <- gPsi[, lower.tri(Psi, diag = TRUE)]
  gPsi[gPsi != 0] <- 1
  
  return(gPsi)
  
}

#' @noRd
#' @importFrom stats lm
# Adds population error using Cudeck method to generated data
# Updated 13.10.2022 -- Marcos
cudeck <- function(R, lambda, Phi, Psi,
                   fit = "rmsr", misfit = "close",
                   method = "minres") {
  
  # Method of Cudeck and Browne (1992):
  
  p <- nrow(lambda)
  q <- ncol(lambda)
  uniquenesses <- diag(Psi)
  
  # Count the number of parameters
  nlambda <- sum(lambda != 0)
  nphi <- sum(Phi[lower.tri(Phi)] != 0)
  npsi <- sum(Psi[lower.tri(Psi, diag = TRUE)] != 0)
  npars <- nlambda + nphi + npsi
  df <- p*(p+1)/2 - npars # Degrees of freedom
  
  if(nlambda + nphi > p*q - 0.5*q*(q-1)) {
    warning("The model is not identified. There exists infinite solutions for the model parameters.")
  }
  
  if(nlambda + nphi + npsi > p*(p+1)/2) {
    warning("The true model has negative degrees of freedom.")
  }
  
  # Create the matrix of derivatives wrt the correlation model:
  
  dS_dL <- gLRhat(lambda, Phi)[, which(lambda != 0)]
  dS_dP <- gPRhat(lambda, Phi)[, which(Phi[lower.tri(Phi)] != 0)]
  dS_dU <- gURhat(p)[, which(Psi[lower.tri(Psi, diag = TRUE)] != 0)]
  gS <- cbind(dS_dL, dS_dP, dS_dU)
  
  if(method == "minres" || method == "ols") {
    
    # Select the nonduplicated elements of the correlation matrix wrt each parameter
    B <- -gS[lower.tri(R, diag = TRUE), ]
    
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
    # vecs <- apply(gS, 2, FUN = function(x) -t(R_inv %*% matrix(x, p, p) %*% R_inv))
    # B <- t(vecs[which(upper.tri(R, diag = TRUE)), ]) %*% D
    # B <- t(B)
    B <- -apply(gS, 2, FUN = function(x) t((R_inv %*% matrix(x, p, p) %*% R_inv)[which(upper.tri(R, diag = TRUE))]) %*% D)
    # The error must be orthogonal to the derivative of each parameter derivative wrt the correlation model
    
  }
  
  # Generate random error:
  
  m <- p+1
  U <- replicate(p, stats::runif(m, 0, 1))
  A1 <- t(U) %*% U
  sq <- diag(1/sqrt(diag(A1)))
  A2 <- sq %*% A1 %*% sq
  diag_u <- diag(sqrt(uniquenesses))
  y <- diag_u %*% A2 %*% diag_u
  y <- y[lower.tri(y, diag = TRUE)]
  # y <- A2[lower.tri(A2, diag = TRUE)]
  # e <- y - B %*% v # equation 7 from Cudeck and Browne (1992)
  Q <- qr.Q(qr(B))
  e <- y - Q %*% t(Q) %*% y
  
  # Get the error matrix:
  
  # Adjust the error to satisfy the desired amount of misfit:
  
  if(method == "minres" || method == "ols") {
    
    E <- matrix(0, p, p)
    E[lower.tri(E, diag = TRUE)] <- e
    E <- t(E) + E
    diag(E) <- 0
    
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
    
    k <- sqrt(2*delta/sum(E*E))
    E <- k*E
    
  } else if(method == "ml") {
    
    E <- matrix(0, p, p)
    E[upper.tri(R, diag = TRUE)] <- e
    E <- t(E) + E
    diag(E) <- 0
    
    if(fit == "rmsr") {
      delta <- "A given RMSR is compatible with multiple maximum likelihood discrepancy values and is not provided"
    } else if(fit == "cfi") {
      null_f <- -log(det(R))
      delta <- (1-misfit)*null_f
    } else if(fit == "rmsea") {
      delta <- misfit^2 * df
    } else if(fit == "raw") {
      delta <- misfit
    }
    
    if(fit == "rmsr") {
      
      k <- sqrt((0.5*p*(p-1))*2*misfit^2/sum(E*E))
      E <- k*E
      
    } else {
      
      constant <- 1e-04 / sqrt(mean(E*E))
      E <- constant*E
      R_inv <- solve(R)
      G <- R_inv %*% E
      x <- suppressWarnings(grid_search(delta, G))
      # x <- sqrt(2*delta/sum(G*G)) # Initial value suggested by Cudeck
      k <- opt(x, delta, G)
      # limits <- c(-1e05, 1e05)
      # k <- GSS(delta, G, limits)
      # k <- grad_descend(delta, G)
      E <- k*E
      
    }
  }
  
  R_error <- R + E
  
  # check for positiveness:
  minimum_eigval <- min(eigen(R_error, symmetric = TRUE, only.values = TRUE)$values)
  if(minimum_eigval <= 0) warning("The matrix was not positive-definite. The amount of misfit may be too big.")
  
  return(list(R_error = R_error, fit = fit, delta = delta, misfit = misfit))
  
}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Minimum residual function for CFA in Yuan method
# Updated 13.10.2022
f_minres <- function(
    x, S, ldetS, q,
    indexes_lambda, lambda_p,
    indexes_phi, phi_p,
    indexes_psi
)
{
  
  p <- nrow(S)
  lambda_p <- length(indexes_lambda)
  Lambda <- matrix(0, p, q)
  Lambda[indexes_lambda] <- x[1:lambda_p]
  phi_p <- length(indexes_phi)
  Phi <- matrix(0, q, q)
  Phi[indexes_phi] <- x[(lambda_p+1):(lambda_p + phi_p)]
  Phi <- t(Phi) + Phi
  diag(Phi) <- 1
  # Psi added 07.10.2022 -- Marcos
  Psi <- matrix(0, p, p)
  Psi[indexes_psi] <- x[-(1:(lambda_p + phi_p))]
  Psi[upper.tri(Psi)] <- t(Psi)[upper.tri(Psi)]
  Rhat <- Lambda %*% Phi %*% t(Lambda) + Psi
  res <- S - Rhat
  f <- 0.5*sum(res*res)
  
  return(f)
  
}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Minimum residual function for CFA in Yuan method
# Updated 07.10.2022
g_minres <- function(
    x, S, ldetS, q,
    indexes_lambda, lambda_p,
    indexes_phi, phi_p,
    indexes_psi
)
{
  
  p <- nrow(S)
  lambda_p <- length(indexes_lambda)
  Lambda <- matrix(0, p, q)
  Lambda[indexes_lambda] <- x[1:lambda_p]
  phi_p <- length(indexes_phi)
  Phi <- matrix(0, q, q)
  Phi[indexes_phi] <- x[(lambda_p+1):(lambda_p + phi_p)]
  Phi <- t(Phi) + Phi
  diag(Phi) <- 1
  # Psi added 07.10.2022 -- Marcos
  Psi <- matrix(0, p, p)
  Psi[indexes_psi] <- x[-(1:(lambda_p + phi_p))]
  Psi[upper.tri(Psi)] <- t(Psi)[upper.tri(Psi)]
  Rhat <- Lambda %*% Phi %*% t(Lambda) + Psi
  res <- S - Rhat
  
  # Change 07.10.2022 -- Marcos
  g1 <- (res %*% Lambda %*% Phi)[indexes_lambda]
  g2 <- (t(Lambda) %*% res %*% Lambda)[indexes_phi]
  # g <- -2*c(g1, g2, 0.5*diag(res))
  res2 <- res
  res2[lower.tri(res2)] <- 2*res[lower.tri(res)]
  g <- -2*c(g1, g2, 0.5*res2[indexes_psi])
  
  return(g)
  
}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Maximum likelihood function for CFA in Yuan method
# Updated 07.10.2022
f_ml <- function(
    x, S, ldetS, q,
    indexes_lambda, lambda_p,
    indexes_phi, phi_p,
    indexes_psi
)
{
  
  p <- nrow(S)
  lambda_p <- length(indexes_lambda)
  Lambda <- matrix(0, p, q)
  Lambda[indexes_lambda] <- x[1:lambda_p]
  phi_p <- length(indexes_phi)
  Phi <- matrix(0, q, q)
  Phi[indexes_phi] <- x[(lambda_p+1):(lambda_p + phi_p)]
  Phi <- t(Phi) + Phi
  diag(Phi) <- 1
  # Psi added 07.10.2022 -- Marcos
  Psi <- matrix(0, p, p)
  Psi[indexes_psi] <- x[-(1:(lambda_p + phi_p))]
  Psi[upper.tri(Psi)] <- t(Psi)[upper.tri(Psi)]
  Rhat <- Lambda %*% Phi %*% t(Lambda) + Psi
  f <- log(det(Rhat)) - ldetS + sum(S*solve(Rhat)) - p
  
  return(f)
  
}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Maximum likelihood function for CFA in Yuan method
# Updated 07.10.2022
g_ml <- function(
    x, S, ldetS, q,
    indexes_lambda, lambda_p,
    indexes_phi, phi_p,
    indexes_psi
)
{
  
  p <- nrow(S)
  lambda_p <- length(indexes_lambda)
  Lambda <- matrix(0, p, q)
  Lambda[indexes_lambda] <- x[1:lambda_p]
  phi_p <- length(indexes_phi)
  Phi <- matrix(0, q, q)
  Phi[indexes_phi] <- x[(lambda_p+1):(lambda_p + phi_p)]
  Phi <- t(Phi) + Phi
  diag(Phi) <- 1
  Psi <- matrix(0, p, p)
  Psi[indexes_psi] <- x[-(1:(lambda_p + phi_p))]
  Psi[upper.tri(Psi)] <- t(Psi)[upper.tri(Psi)]
  
  Rhat <- Lambda %*% Phi %*% t(Lambda) + Psi
  Rhat_inv <- solve(Rhat)
  Ri_res_Ri <- 2*Rhat_inv %*% (Rhat - S) %*% Rhat_inv
  Ri_res_Ri2 <- Ri_res_Ri
  Ri_res_Ri2[lower.tri(Ri_res_Ri2)] <- 2*Ri_res_Ri[lower.tri(Ri_res_Ri)]
  
  # Joreskog (page 10; 1965) Testing a simple structure in factor analysis
  # g <- c(c(Ri_res_Ri %*% Lambda %*% Phi)[indexes_lambda],
  #        c(t(Lambda) %*% Ri_res_Ri %*% Lambda)[indexes_phi],
  #        diag(Ri_res_Ri)*0.5)
  g <- c(c(Ri_res_Ri %*% Lambda %*% Phi)[indexes_lambda],
         c(t(Lambda) %*% Ri_res_Ri %*% Lambda)[indexes_phi],
         Ri_res_Ri2[indexes_psi]*0.5)
  
  return(g)
  
}

#' @noRd
#' @importFrom stats runif nlminb
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# CFA function
# Updated 07.10.2022 -- Marcos
CFA <- function(S, target, targetphi, targetpsi = diag(nrow(target)), method = "minres") {
  
  p <- nrow(target)
  q <- ncol(target)
  indexes_lambda <- which(target != 0) # Which lambdas are estimated
  indexes_phi <- which(targetphi != 0 & lower.tri(targetphi)) # Which phis are estimated
  indexes_psi <- which(targetpsi != 0 & lower.tri(targetpsi, diag = TRUE)) # Which psies are estimated
  lambda_p <- length(indexes_lambda) # Number of lambda parameters
  phi_p <- length(indexes_phi) # Number of phi parameters
  psi_p <- length(indexes_psi) # Number of psi parameters
  
  init_diag_psi <- 1/diag(solve(S)) # Initial diagonal psi parameter values
  init_psi <- rep(0, times = psi_p)
  diag_indexes <- (p+1)*0:(p-1)+1 # Indexes for the diagonal of Psi
  offdiag_indexes <- which(targetpsi != 0 & lower.tri(targetpsi)) # Indexes for the off-diagonal of Psi
  cor_res_indexes <- which(indexes_psi %in% offdiag_indexes) # Indexes for correlated residuals
  # Allocate init_diag_psi in the positions of the vector corresponding to the diagonal of Psi:
  init_psi[-cor_res_indexes] <- init_diag_psi
  
  lower_psi <- rep(0.005, psi_p) # Lower bounds for the uniquenessess
  lower_psi[cor_res_indexes] <- -0.995 # Lower bounds for correlated residuals
  upper_psi <- rep(0.995, psi_p) # Upper bounds for correlated residuals
  lower <- c(rep(-Inf, lambda_p), rep(-1, phi_p), lower_psi)
  upper <- c(rep(Inf, lambda_p), rep(1, phi_p), upper_psi)
  
  x <- c(runif(lambda_p), rep(0, phi_p), init_psi)
  
  if(method == "minres") {
    
    ldetS <- NULL
    f <- f_minres
    g <- g_minres
    
  } else if(method == "ml") {
    
    ldetS <- log(det(S))
    f <- f_ml
    g <- g_ml
    
  }
  
  cfa <- nlminb(
    start = x, objective = f, gradient = g,
    lower = lower, upper = upper,
    S = S, ldetS = ldetS, q = q,
    indexes_lambda = indexes_lambda, lambda_p = lambda_p,
    indexes_phi = indexes_phi, phi_p = phi_p,
    indexes_psi = indexes_psi,
    control = list(iter.max = 1e4, eval.max = 1e4)
  )
  
  # Arrange lambda parameter estimates:
  lambda_hat <- matrix(0, p, q)
  lambda_hat[indexes_lambda] <- cfa$par[1:lambda_p]
  
  # Arrange phi parameter estimates:
  phi_hat <- matrix(0, q, q)
  phi_hat[indexes_phi] <- cfa$par[(lambda_p+1):(lambda_p + phi_p)]
  phi_hat <- t(phi_hat) + phi_hat
  diag(phi_hat) <- 1
  
  # Arrange psi parameter estimates:
  psi_hat <- matrix(0, p, p)
  psi_hat[indexes_psi] <- cfa$par[-(1:(lambda_p + phi_p))]
  psi_hat[upper.tri(psi_hat)] <- t(psi_hat)[upper.tri(psi_hat)]
  
  # Model matrix:
  S_hat <- lambda_hat %*% phi_hat %*% t(lambda_hat) + psi_hat
  uniquenesses_hat <- diag(psi_hat)
  diag(S_hat) <- 1 # Fix rounding errors from the optimization
  residuals <- S - S_hat
  
  # Degrees of freedom:
  df <- p*(p+1)/2 - (lambda_p + phi_p + psi_p)
  
  results <- list(f = cfa$objective, convergence = cfa$convergence,
                  iterations = cfa$iterations, df = df,
                  lambda = lambda_hat, phi = phi_hat,
                  psi = psi_hat, uniquenesses = uniquenesses_hat,
                  model = S_hat, residuals = residuals)
  
  return(results)
  
}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Adds population error using Yuan method to generated data
# Updated 07.10.2022 -- Marcos
yuan <- function(R, lambda, Phi, Psi,
                 fit = "rmsr", misfit = "close",
                 method = "minres") {
  
  p <- nrow(R)
  q <- ncol(lambda)
  uniquenesses <- diag(Psi)
  
  # Count the number of parameters
  nlambda <- sum(lambda != 0)
  nphi <- sum(Phi[lower.tri(Phi)] != 0)
  npsi <- sum(Psi[lower.tri(Psi, diag = TRUE)] != 0)
  npars <- nlambda + nphi + npsi
  df <- p*(p+1)/2 - npars # Degrees of freedom
  
  if(nlambda + nphi > p*q - 0.5*q*(q-1)) {
    warning("The population model is not identified. There exists infinite solutions for the model parameters.")
  }
  
  if(nlambda + nphi + npsi > p*(p+1)/2) {
    warning("The population model has negative degrees of freedom.")
  }
  
  # Add an small error to the population parameters
  # lambda_error <- lambda - 1e-04
  # Rerror <- lambda_error %*% Phi %*% t(lambda_error) + Psi; diag(Rerror) <- 1
  Rerror <- R
  Rerror[lower.tri(R)] <- Rerror[lower.tri(R)] + runif(0.5*p*(p-1), -1e-06, 1e-06)
  Rerror[upper.tri(R)] <- t(Rerror)[upper.tri(R)]
  
  # Create the FA model
  target <- ifelse(lambda != 0, 1, 0)
  targetphi <- ifelse(Phi != 0, 1, 0)
  targetpsi <- ifelse(Psi != 0, 1, 0)
  cfa <- CFA(Rerror, target, targetphi, targetpsi, method = method)
  Phat <- cfa$model
  
  # Get the error matrix:
  E <- Rerror - Phat
  # Hopefully, the error is orthogonal to the derivative of each parameter derivative wrt the discrepancy function
  
  # Adjust the error to satisfy the desired amount of misfit:
  
  if(method == "minres" || method == "ols") {
    
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
    
    k <- sqrt(2*delta/sum(E*E))
    E <- k*E
    
  } else if(method == "ml") {
    
    if(fit == "rmsr") {
      delta <- "A given RMSR is compatible with multiple maximum likelihood discrepancy values and is not provided"
    } else if(fit == "cfi") {
      null_f <- -log(det(R))
      delta <- (1-misfit)*null_f
    } else if(fit == "rmsea") {
      delta <- misfit^2 * df
    } else if(fit == "raw") {
      delta <- misfit
    }
    
    if(fit == "rmsr") {
      
      k <- sqrt((0.5*p*(p-1))*2*misfit^2/sum(E*E))
      E <- k*E
      
    } else {
      
      constant <- 1e-04 / sqrt(mean(E*E))
      E <- constant*E # Fix this to avoid NAs
      R_inv <- solve(R)
      G <- R_inv %*% E
      x <- suppressWarnings(grid_search(delta, G))
      # x <- sqrt(2*delta/sum(G*G)) # Initial value suggested by Cudeck
      k <- opt(x, delta, G)
      # limits <- c(-1e05, 1e05)
      # k <- GSS(delta, G, limits)
      # k <- grad_descend(delta, G)
      E <- k*E
      
    }
  }
  
  R_error <- Phat + E
  
  # check for positiveness:
  minimum_eigval <- min(eigen(R_error, symmetric = TRUE, only.values = TRUE)$values)
  if(minimum_eigval <= 0) warning("The matrix was not positive-definite. The amount of misfit may be too big.")
  
  return(list(R_error = R_error, fit = fit, delta = delta, misfit = misfit))
  
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%
# add_local_dependence ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Adds correlated residuals to generated data
# Updated 01.11.2022
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
        population_correlation[
          correlated_residuals[i,1],
          correlated_residuals[i,2]
        ] <- original_correlation[
          correlated_residuals[i,1],
          correlated_residuals[i,2]
        ] + add_residuals[i]
        
        # Ensure symmetric
        population_correlation[
          correlated_residuals[i,2],
          correlated_residuals[i,1]
        ] <- population_correlation[
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
    data = data,
    population_correlation = population_correlation,
    parameters = parameters,
    correlated_residuals = correlated_residuals_df,
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

#%%%%%%%%%%%%%%%%%%%
# factor_forest ----
#%%%%%%%%%%%%%%%%%%%

#' @noRd
#' @importFrom stats wilcox.test
#' @importFrom graphics abline
# EFA comparison
# Updated 30.09.2022
EFA.Comp.Data <- function(Data, F.Max, N.Pop = 10000, N.Samples = 500, Alpha = .30, Graph = F, Spearman = F, use)
{
  # Data = N (sample size) x k (number of variables) data matrix
  # F.Max = largest number of factors to consider
  # N.Pop = size of finite populations of comparison data (default = 10,000 cases)
  # N.Samples = number of samples drawn from each population (default = 500)
  # Alpha = alpha level when testing statistical significance of improvement with add'l factor (default = .30) 
  # Graph = whether to plot the fit of eigenvalues to those for comparison data (default = F)
  # Spearman = whether to use Spearman rank-order correlations rather than Pearson correlations (default = F)
  
  N <- dim(Data)[1]
  k <- dim(Data)[2]
  if (Spearman) Cor.Type <- "spearman" else Cor.Type <- "pearson"
  cor.Data <- cor(Data, method = Cor.Type, use = use)
  Eigs.Data <- eigen(cor.Data)$values
  RMSR.Eigs <- matrix(0, nrow = N.Samples, ncol = F.Max)
  Sig <- T
  F.CD <- 1
  while ((F.CD <= F.Max) & (Sig))
  {
    Pop <- GenData(Data, N.Factors = F.CD, N = N.Pop, Cor.Type = Cor.Type, use = use)
    for (j in 1:N.Samples)
    {
      Samp <- Pop[sample(1:N.Pop, size = N, replace = T),]
      cor.Samp <- cor(Samp, method = Cor.Type, use = use)
      Eigs.Samp <- eigen(cor.Samp)$values
      RMSR.Eigs[j,F.CD] <- sqrt(sum((Eigs.Samp - Eigs.Data) * (Eigs.Samp - Eigs.Data)) / k)
    }
    if (F.CD > 1) Sig <- (wilcox.test(RMSR.Eigs[,F.CD], RMSR.Eigs[,(F.CD - 1)], "less")$p.value < Alpha)
    if (Sig) F.CD <- F.CD + 1
  }
  if (Graph)
  {
    if (Sig) x.max <- F.CD - 1
    else x.max <- F.CD
    ys <- apply(RMSR.Eigs[,1:x.max], 2, mean)
    plot(x = 1:x.max, y = ys, ylim = c(0, max(ys)), xlab = "Factor", ylab = "RMSR Eigenvalue", type = "b", 
         main = "Fit to Comparison Data")
    abline(v = F.CD - 1, lty = 3)
  }
  return(F.CD - 1)
}

#' @noRd
# Data generation
# Updated 30.09.2022
GenData <- function(Supplied.Data, N.Factors, N, Max.Trials = 5, Initial.Multiplier = 1, Cor.Type, use){
  
  k <- dim(Supplied.Data)[2]
  Data <- matrix(0, nrow = N, ncol = k)            # Matrix to store the simulated data
  Distributions <- matrix(0, nrow = N, ncol = k)   # Matrix to store each variable's score distribution
  Iteration <- 0                                   # Iteration counter
  Best.RMSR <- 1                                   # Lowest RMSR correlation
  Trials.Without.Improvement <- 0                  # Trial counter
  
  # Generate distribution for each variable (step 2) -------------------------------------------------------------
  
  for (i in 1:k)
    Distributions[,i] <- sort(sample(na.omit(Supplied.Data[,i]), size = N, replace = T))
  
  # Calculate and store a copy of the target correlation matrix (step 3) -----------------------------------------
  
  Target.Corr <- cor(Supplied.Data, method = Cor.Type, use = use)
  Intermediate.Corr <- Target.Corr
  
  # Generate random normal data for shared and unique components, initialize factor loadings (steps 5, 6) --------
  
  Shared.Comp <- matrix(rnorm(N * N.Factors, 0, 1), nrow = N, ncol = N.Factors)
  Unique.Comp <- matrix(rnorm(N * k, 0, 1), nrow = N, ncol = k)
  Shared.Load <- matrix(0, nrow = k, ncol = N.Factors)
  Unique.Load <- matrix(0, nrow = k, ncol = 1)
  
  # Begin loop that ends when specified number of iterations pass without improvement in RMSR correlation --------
  
  while (Trials.Without.Improvement < Max.Trials)
  {
    Iteration <- Iteration + 1
    
    # Calculate factor loadings and apply to reproduce desired correlations (steps 7, 8) ---------------------------
    
    Fact.Anal <- Factor.Analysis(Intermediate.Corr, Corr.Matrix = TRUE, N.Factors = N.Factors, Cor.Type = Cor.Type, use = use)
    if (N.Factors == 1) Shared.Load[,1] <- Fact.Anal$loadings
    else 
      for (i in 1:N.Factors)
        Shared.Load[,i] <- Fact.Anal$loadings[,i]
      Shared.Load[Shared.Load > 1] <- 1
      Shared.Load[Shared.Load < -1] <- -1
      if (Shared.Load[1,1] < 0) Shared.Load <- Shared.Load * -1
      for (i in 1:k)
        if (sum(Shared.Load[i,] * Shared.Load[i,]) < 1) Unique.Load[i,1] <- 
        (1 - sum(Shared.Load[i,] * Shared.Load[i,]))
      else Unique.Load[i,1] <- 0
      Unique.Load <- sqrt(Unique.Load)
      for (i in 1:k)
        Data[,i] <- (Shared.Comp %*% t(Shared.Load))[,i] + Unique.Comp[,i] * Unique.Load[i,1]
      
      # Replace normal with nonnormal distributions (step 9) ---------------------------------------------------------
      
      for (i in 1:k)
      {
        Data <- Data[sort.list(Data[,i]),]
        Data[,i] <- Distributions[,i]
      }
      
      # Calculate RMSR correlation, compare to lowest value, take appropriate action (steps 10, 11, 12) --------------
      
      Reproduced.Corr <- cor(Data, method = Cor.Type, use = use)
      Residual.Corr <- Target.Corr - Reproduced.Corr
      RMSR <- sqrt(sum(Residual.Corr[lower.tri(Residual.Corr)] * Residual.Corr[lower.tri(Residual.Corr)]) / 
                     (.5 * (k * k - k)))
      if (RMSR < Best.RMSR)
      {
        Best.RMSR <- RMSR
        Best.Corr <- Intermediate.Corr
        Best.Res <- Residual.Corr
        Intermediate.Corr <- Intermediate.Corr + Initial.Multiplier * Residual.Corr
        Trials.Without.Improvement <- 0
      }
      else 
      {
        Trials.Without.Improvement <- Trials.Without.Improvement + 1
        Current.Multiplier <- Initial.Multiplier * .5 ^ Trials.Without.Improvement
        Intermediate.Corr <- Best.Corr + Current.Multiplier * Best.Res
      }
  }
  
  Fact.Anal <- Factor.Analysis(Best.Corr, Corr.Matrix = TRUE, N.Factors = N.Factors, Cor.Type = Cor.Type, use = use)
  if (N.Factors == 1) Shared.Load[,1] <- Fact.Anal$loadings
  else
    for (i in 1:N.Factors)
      Shared.Load[,i] <- Fact.Anal$loadings[,i]
  Shared.Load[Shared.Load > 1] <- 1
  Shared.Load[Shared.Load < -1] <- -1
  if (Shared.Load[1,1] < 0) Shared.Load <- Shared.Load * -1
  for (i in 1:k)
    if (sum(Shared.Load[i,] * Shared.Load[i,]) < 1) Unique.Load[i,1] <-
    (1 - sum(Shared.Load[i,] * Shared.Load[i,]))
  else Unique.Load[i,1] <- 0
  Unique.Load <- sqrt(Unique.Load)
  for (i in 1:k)
    Data[,i] <- (Shared.Comp %*% t(Shared.Load))[,i] + Unique.Comp[,i] * Unique.Load[i,1]
  Data <- apply(Data, 2, scale) # standardizes each variable in the matrix
  for (i in 1:k)
  {
    Data <- Data[sort.list(Data[,i]),]
    Data[,i] <- Distributions[,i]
  }
  
  return(Data)
}

#' @noRd
# Factor analysis
# Updated 30.09.2022
Factor.Analysis <- function(Data, Corr.Matrix = FALSE, Max.Iter = 50, N.Factors = 0, Cor.Type, use)
{
  Data <- as.matrix(Data)
  k <- dim(Data)[2]
  if (N.Factors == 0)
  {
    N.Factors <- k
    Determine <- T
  }
  else Determine <- F
  if (!Corr.Matrix) Cor.Matrix <- cor(Data, method = Cor.Type, use = use)
  else Cor.Matrix <- Data
  Criterion <- .001
  Old.H2 <- rep(99, k)
  H2 <- rep(0, k)
  Change <- 1
  Iter <- 0
  Factor.Loadings <- matrix(nrow = k, ncol = N.Factors)
  while ((Change >= Criterion) & (Iter < Max.Iter))
  {
    Iter <- Iter + 1
    Eig <- eigen(Cor.Matrix)
    L <- sqrt(Eig$values[1:N.Factors])
    for (i in 1:N.Factors)
      Factor.Loadings[,i] <- Eig$vectors[,i] * L[i]
    for (i in 1:k)
      H2[i] <- sum(Factor.Loadings[i,] * Factor.Loadings[i,])
    Change <- max(abs(Old.H2 - H2))
    Old.H2 <- H2
    diag(Cor.Matrix) <- H2
  }
  if (Determine) N.Factors <- sum(Eig$values > 1)
  return(list(loadings = Factor.Loadings[,1:N.Factors], factors = N.Factors))
}

#' @noRd
# {xgboost}: createDMatrixFromTask
# Updated 30.09.2022
createDMatrixFromTask = function(task, weights = NULL) {

  data = mlr::getTaskData(task, target.extra = TRUE)
  data$data = BBmisc::convertDataFrameCols(data$data, ints.as.num = TRUE)
  if (mlr::getTaskType(task) == "classif")  {
    cl = mlr::getTaskClassLevels(task)
    data$target =  match(as.character(data$target), cl) - 1
  }

  if (!is.null(weights))
    xgboost::xgb.DMatrix(data = data.matrix(data$data), label = data$target, weight = weights)
  else if (!is.null(task$weights))
    xgboost::xgb.DMatrix(data = data.matrix(data$data), label = data$target, weight = task$weights)
  else
    xgboost::xgb.DMatrix(data = data.matrix(data$data), label = data$target)
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
# Updated 22.11.2022
skew_continuous <- function(
    skewness,
    data = NULL,
    sample_size = 1000000,
    tolerance = 0.00001
)
{
  
  # Check for zero skew (skip adding skew)
  if(skewness == 0){
    return(data)
  }
  
  # Obtain absolute skewness
  if(sign(skewness) == -1){
    skewness <- abs(skewness)
    flip <- TRUE
  }else{
    flip <- FALSE
  }
  
  # Generate data
  if(is.null(data)){
    data <- rnorm(sample_size)
  }
  
  # Kurtosis
  kurtosis <- 1
  
  # Initialize increments
  increments <- 0.01
  
  # Seek along a range of skews
  skew_values <- seq(
    -2, 2, increments
  )
  
  # Compute skews
  skews <- unlist(lapply(skew_values, function(x){
    # Skew data
    skew_data <- sinh(
      kurtosis * (asinh(data) + x) 
    )
    
    # Observed skew in data
    psych::skew(skew_data)
  }))
  
  # Compute minimum index
  minimum <- which.min(abs(skewness - skews))
  
  # Check for whether skewness is found
  while(abs(skewness - skews[minimum]) > tolerance){
    
    # Check for minimum value
    if(minimum == 1){
      kurtosis <- kurtosis - 0.1
    }else if(minimum == length(skews)){
      kurtosis <- kurtosis + 0.1
    }else{
      
      # Decrease increments
      increments <- 0.01 * 0.1
      
      # Seek along a range of skews
      skew_values <- seq(
        skew_values[minimum - 1],
        skew_values[minimum + 1],
        length.out = 100
      )

    }
    
    # Compute skews
    skews <- unlist(lapply(skew_values, function(x){
      # Skew data
      skew_data <- sinh(
        kurtosis * (asinh(data) + x) 
      )
      
      # Observed skew in data
      psych::skew(skew_data)
    }))
    
    # Compute minimum index
    minimum <- which.min(abs(skewness - skews))
    
  }
  
  # Compute final skew data
  skew_data <- sinh(
    kurtosis * (asinh(data) + skew_values[minimum]) 
  )
  
  # Re-scale
  skew_data <- scale(skew_data)
  
  # Flip skew?
  if(isTRUE(flip)){
    skew_data <- -skew_data
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

#' @noRd
#' @importFrom stats na.omit
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











