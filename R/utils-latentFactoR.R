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
# Updated 24.11.2022
srmr <- function(population, error)
{

  # Obtain lower triangles
  lower_triangle <- lower.tri(population)

  # Return SRMR
  return(sqrt(mean((population[lower_triangle] - error[lower_triangle])^2)))

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
# Updated 26.05.2024
dxt <- function(X) {

  # derivative wrt transpose (just a permutation matrix)

  dimensions <- dim(X)

  # p <- nrow(X)
  # q <- ncol(X)
  # pq <- p*q

  pq <- dimensions[1] * dimensions[2]

  # res <- array(0, dim = c(pq, pq))
  res <- matrix(0, pq, pq)
  # null <- matrix(0, p, q)
  null <- matrix(0, dimensions[1], dimensions[2])

  # for(i in 1:pq) {
  for(i in seq_len(pq)){
    temp <- null
    temp[i] <- 1
    # res[,i] <- c(t(temp))
    res[,i] <- t(temp)
  }

  return(res)

}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# LRhat function for Cudeck method
# Updated 26.05.2024
gLRhat <- function(Lambda, Phi) {

  # derivative of Lambda wrt Rhat

  I <- diag(dim(Lambda)[1])

  # pre-compute Lambda %*% Phi
  LP <- Lambda %*% Phi

  # p <- nrow(Lambda)
  # g1 <- (Lambda %*% Phi) %x% diag(p)
  g1 <- LP %x% I
  # g21 <- diag(p) %x% (Lambda %*% Phi)
  g21 <- I %x% LP
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
# Updated 26.05.2024
gPRhat <- function(Lambda, Phi) {

  g1 <- Lambda %x% Lambda
  g2 <- g1 %*% dxt(Phi)
  g <- g1 + g2
  # g <- g[,which(lower.tri(Phi))]
  return(g[, lower.tri(Phi), drop = FALSE])

  # Ensure matrix
  # if(!is.matrix(g)){
  #   g <- matrix(g, ncol = 1)
  # }
  #
  # return(g)

}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Rhat function for Cudeck method
# Updated 26.05.2024
guRhat <- function(p) {

  gu <- matrix(0, p*p, p)

  # for(i in 1:p) {
  for(i in seq_len(p)){
    # index <- (i-1) * p + i
    gu[(i-1) * p + i, i] <- 1
  }

  return(gu)

}

#' @noRd
# Rhat function for Cudeck method
# Updated 07.10.2022
# For Psi
gURhat <- function(p) {

  # pcov <- p*(p+1)*0.5

  Psi <- diag(p)
  # gPsi <- diag(p) %x% diag(p)
  gPsi <- Psi %x% Psi
  gPsi <- gPsi + dxt(Psi) %*% gPsi
  gPsi <- gPsi[, lower.tri(Psi, diag = TRUE)]
  gPsi[gPsi != 0] <- 1

  return(gPsi)

}

#' @noRd
# Adds population error using Cudeck method to generated data
# Updated 26.05.2024
cudeck <- function(R, lambda, Phi, Psi,
                   fit = "rmsr", misfit = "close",
                   method = "minres", efa = FALSE) {

  # Method of Cudeck and Browne (1992):

  p <- nrow(lambda)
  q <- ncol(lambda)
  uniquenesses <- diag(Psi)

  # Pre-compute indices
  lambda_idx <- lambda != 0
  phi_idx <- Phi[lower.tri(Phi)] != 0
  psi_idx <- Psi[lower.tri(Psi, diag = TRUE)] != 0

  # Count the number of parameters
  # nlambda <- sum(lambda != 0)
  nlambda <- sum(lambda_idx)
  # nphi <- sum(Phi[lower.tri(Phi)] != 0)
  nphi <- sum(phi_idx)
  # npsi <- sum(Psi[lower.tri(Psi, diag = TRUE)] != 0)
  npsi <- sum(psi_idx)
  npars <- nlambda + nphi + npsi
  # df <- p*(p+1)/2 - npars # Degrees of freedom

  # pre-compute df
  pre_df <- p * (p + 1) / 2
  df <- pre_df - npars

  if(nlambda + nphi > p*q - 0.5*q*(q-1)) {
    warning("The model is not identified. There exists infinite solutions for the model parameters.")
  }

  if(nlambda + nphi + npsi > pre_df) {
    warning("The true model has negative degrees of freedom.")
  }

  # Create the matrix of derivatives wrt the correlation model:

  if(efa) {

    dS_dL <- gLRhat(lambda, Phi)
    dS_dU <- gURhat(p)[, which(Psi[lower.tri(Psi, diag = TRUE)] != 0)]
    gS <- cbind(dS_dL, dS_dU)

  } else {

    # dS_dL <- gLRhat(lambda, Phi)[, which(lambda != 0)]
    dS_dL <- gLRhat(lambda, Phi)[, lambda_idx]
    # dS_dP <- gPRhat(lambda, Phi)[, which(Phi[lower.tri(Phi)] != 0)]
    dS_dP <- gPRhat(lambda, Phi)[, phi_idx]
    # dS_dU <- gURhat(p)[, which(Psi[lower.tri(Psi, diag = TRUE)] != 0)]
    dS_dU <- gURhat(p)[, psi_idx]
    gS <- cbind(dS_dL, dS_dP, dS_dU)

  }

  if(method == "minres" || method == "ols") {

    # Select the nonduplicated elements of the correlation matrix wrt each parameter
    B <- -gS[lower.tri(R, diag = TRUE), ]

  } else if(method == "ml") {

    # K <- transition(p)
    # MP_inv <- solve(t(K) %*% K) %*% t(K)
    # D <- MP_inv %*% t(MP_inv)
    # indexes <- vector(length = p)
    # indexes[1] <- 1
    # for(i in 2:p) {
    #   increment <- i
    #   indexes[i] <- indexes[i-1]+increment
    # }
    indexes <- cumsum(seq_len(p))
    # D <- matrix(0, p*(p+1)/2, p*(p+1)/2)
    # diag(D) <- 2
    D <- diag(2, pre_df, pre_df)
    diag(D)[indexes] <- 1
    R_inv <- solve(R)
    # vecs <- apply(gS, 2, FUN = function(x) -t(R_inv %*% matrix(x, p, p) %*% R_inv))
    # B <- t(vecs[which(upper.tri(R, diag = TRUE)), ]) %*% D
    # B <- t(B)
    # B <- -apply(gS, 2, FUN = function(x) t((R_inv %*% matrix(x, p, p) %*% R_inv)[which(upper.tri(R, diag = TRUE))]) %*% D)

    # Pre-compute upper triangle
    upper_triangle <- upper.tri(R, diag = TRUE)

    B <- -apply(
      gS, 2,
      FUN = function(x){
          crossprod((R_inv %*% matrix(x, p, p) %*% R_inv)[upper_triangle], D)
      }
    )
    # The error must be orthogonal to the derivative of each parameter wrt correlation model

  }

  # Generate random error:

  m <- p + 1
  U <- replicate(p, runif_xoshiro(m, min = 1, max = 1))
  A1 <- crossprod(U) # t(U) %*% U
  sq <- diag(1 / sqrt(diag(A1)))
  A2 <- sq %*% A1 %*% sq
  diag_u <- diag(sqrt(uniquenesses))
  y <- diag_u %*% A2 %*% diag_u
  y <- y[lower.tri(y, diag = TRUE)]
  # y <- A2[lower.tri(A2, diag = TRUE)]
  # v <- MASS::ginv(t(B) %*% B) %*% t(B) %*% y
  # e <- y - B %*% v # equation 7 from Cudeck and Browne (1992)
  # Q <- qr.Q(qr(B))
  e <- y - tcrossprod(qr.Q(qr(B))) %*% y # y - Q %*% t(Q) %*% y

  # Adjust the error to satisfy the desired amount of misfit:

  if(method == "minres" || method == "ols") {

    E <- matrix(0, p, p)
    E[lower.tri(E, diag = TRUE)] <- e
    E <- t(E) + E
    # diag(E) <- diag(E)/2
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

      # if(misfit == "close") {
      #   r2 <- mean(1-uniquenesses)
      #   misfit <- 0.05*r2
      # } else if(misfit == "acceptable") {
      #   r2 <- mean(1-uniquenesses)
      #   misfit <- 0.10*r2
      # }

      # Switch out misfit with a value
      misfit_value <- switch(
        as.character(misfit),
        "close" = 0.05,
        "acceptable" = 0.10,
        as.numeric(misfit)
      )

      # Compute misfit
      misfit <- misfit_value * (mean(1 - uniquenesses))

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

  R_error <- R + E

  if(method == "ml" & fit == "rmsr") {
    delta <- log(det(R)) - log(det(R_error)) + sum(R_error*solve(R)) - nrow(R)
  }

  # Check for positiveness:
  minimum_eigval <- min(eigen(R_error, symmetric = TRUE, only.values = TRUE)$values)
  positive <- TRUE
  if(minimum_eigval <= 0) {
    warning("The matrix was not positive-definite. The amount of misfit may be too big.")
    positive <- FALSE
  }

  return(list(R_error = R_error, fit = fit, delta = delta, misfit = misfit, positive = positive))

}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Minimum residual function for CFA in Yuan method
# Updated 26.05.2024
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
  Rhat <- Lambda %*% tcrossprod(Phi, Lambda) + Psi # Lambda %*% Phi %*% t(Lambda) + Psi
  res <- S - Rhat
  f <- 0.5*sum(res*res)

  return(f)

}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Minimum residual function for CFA in Yuan method
# Updated 26.05.2024
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
  Rhat <- Lambda %*% tcrossprod(Phi, Lambda) + Psi # Lambda %*% Phi %*% t(Lambda) + Psi
  res <- S - Rhat

  # Change 07.10.2022 -- Marcos
  g1 <- (res %*% Lambda %*% Phi)[indexes_lambda]
  g2 <- (crossprod(Lambda, res) %*% Lambda)[indexes_phi] # (t(Lambda) %*% res %*% Lambda)[indexes_phi]
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
# Updated 26.05.2024
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
  Rhat <- Lambda %*% tcrossprod(Phi, Lambda) + Psi # Lambda %*% Phi %*% t(Lambda) + Psi
  f <- log(det(Rhat)) - ldetS + sum(S*solve(Rhat)) - p

  return(f)

}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# Maximum likelihood function for CFA in Yuan method
# Updated 26.05.2024
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

  Rhat <- Lambda %*% tcrossprod(Phi, Lambda) + Psi # Lambda %*% Phi %*% t(Lambda) + Psi
  Rhat_inv <- solve(Rhat)
  Ri_res_Ri <- 2*Rhat_inv %*% (Rhat - S) %*% Rhat_inv
  Ri_res_Ri2 <- Ri_res_Ri
  Ri_res_Ri2[lower.tri(Ri_res_Ri2)] <- 2*Ri_res_Ri[lower.tri(Ri_res_Ri)]

  # Joreskog (page 10; 1965) Testing a simple structure in factor analysis
  # g <- c(c(Ri_res_Ri %*% Lambda %*% Phi)[indexes_lambda],
  #        c(t(Lambda) %*% Ri_res_Ri %*% Lambda)[indexes_phi],
  #        diag(Ri_res_Ri)*0.5)
  g <- c(c(Ri_res_Ri %*% Lambda %*% Phi)[indexes_lambda],
         c(crossprod(Lambda, Ri_res_Ri) %*% Lambda)[indexes_phi], # c(t(Lambda) %*% Ri_res_Ri %*% Lambda)[indexes_phi],
         Ri_res_Ri2[indexes_psi]*0.5)

  return(g)

}

#' @noRd
# From {bifactor} version 0.1.0
# Accessed on 17.09.2022
# CFA function
# Updated 07.10.2022 -- Marcos
CFA <- function(S, target, targetphi, targetpsi = diag(nrow(target)), method = "minres", W = NULL) {

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

  x <- c(runif_xoshiro(lambda_p), rep(0, phi_p), init_psi)

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
  S_hat <- lambda_hat %*% tcrossprod(phi_hat, lambda_hat) + psi_hat # lambda_hat %*% phi_hat %*% t(lambda_hat) + psi_hat
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
# Accessed on 26.05.2024
# Adds population error using Yuan method to generated data
# Updated 26.05.2024
yuan <- function(R, lambda, Phi, Psi,
                 fit = "rmsr", misfit = "close",
                 method = "minres", W = NULL) {

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
  Rerror[lower.tri(R)] <- Rerror[lower.tri(R)] + runif_xoshiro(0.5*p*(p-1), min = -1e-06, max = 1e-06)
  Rerror[upper.tri(R)] <- t(Rerror)[upper.tri(R)]

  # Create the FA model
  target <- ifelse(lambda != 0, 1, 0)
  targetphi <- ifelse(Phi != 0, 1, 0)
  targetpsi <- ifelse(Psi != 0, 1, 0)
  cfa <- CFA(Rerror, target, targetphi, targetpsi, method = method, W = NULL)
  Phat <- cfa$model

  # Get the error matrix:
  E <- Rerror - Phat
  # Hopefully, the error is orthogonal to the derivative of each parameter wrt correlation model

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

  if(method == "ml" & fit == "rmsr") {
    delta <- log(det(R)) - log(det(R_error)) + sum(R_error*solve(R)) - nrow(R)
  }

  # Check for positiveness:
  minimum_eigval <- min(eigen(R_error, symmetric = TRUE, only.values = TRUE)$values)
  positive <- TRUE
  if(minimum_eigval <= 0) {
    warning("The matrix was not positive-definite. The amount of misfit may be too big.")
    positive <- FALSE
  }

  return(list(R_error = R_error, fit = fit, delta = delta, misfit = misfit, positive = positive))

}

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

#' @noRd
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

  Shared.Comp <- matrix(rnorm_ziggurat(N * N.Factors), nrow = N, ncol = N.Factors)
  Unique.Comp <- matrix(rnorm_ziggurat(N * k), nrow = N, ncol = k)
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











