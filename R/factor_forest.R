#' Estimate Number of Dimensions using Factor Forest
#'
#' Estimates the number of dimensions in data using the
#' pre-trained XGBoost model from Goretzko and Buhner
#' (2020, 2022). See examples to get started
#'
#' @details
#' The model's 184 features are computed as in the authors' evaluation code:
#' the sample size and number of variables, summary features of the
#' correlation matrix and its eigenvalues, parallel analysis, the empirical
#' Kaiser criterion, the comparison data approach, and the first 80
#' eigenvalues of the correlation matrix and of the reduced correlation matrix.
#'
#' Parallel analysis follows \code{\link[psych]{fa.parallel}} with
#' \code{fa = "fa"}: the one-factor minimum residual eigenvalues of the data
#' are compared to the 95th percentile of those from 20 sets of random normal
#' data with the same dimensions. Because the reference only depends on the
#' sample size and number of variables, no correlations of resampled data are
#' computed. The comparison data approach (Ruscio & Roche, 2012) uses an
#' optimized implementation of the original algorithm.
#'
#' The pre-trained model is included with the package (an XGBoost booster
#' saved in JSON format) and is loaded once per R session
#'
#' @param data Matrix or data frame.
#' A dataset with all numeric values (rows = cases, columns = variables).
#' Correlation matrices are not supported because the comparison data
#' approach requires the data
#'
#' @param sample_size Numeric (length = 1).
#' Deprecated and ignored (the sample size is the number of rows of \code{data})
#'
#' @param maximum_factors Numeric (length = 1).
#' Maximum number of factors to search over in the comparison data approach.
#' Defaults to \code{8}
#'
#' @param PA_correlation Character (length = 1).
#' Type of correlation used in parallel analysis.
#' Defaults to \code{"cor"}:
#'
#' \itemize{
#'
#' \item \code{"cor"} --- Pearson's correlation
#'
#' \item \code{"poly"} --- Polychoric correlation
#'
#' \item \code{"tet"} --- Tetrachoric correlation
#'
#' }
#'
#' @param seed Numeric (length = 1).
#' Seed for reproducible results. The same \code{seed} gives the same results
#' (and does not change the random number stream of the R session).
#' Defaults to \code{NULL} (random seed)
#'
#' @return Returns a list (class \code{"lf_ff"}) containing:
#'
#' \item{dimensions}{Number of dimensions identified}
#'
#' \item{probabilities}{Probability that the number of dimensions
#' is most likely}
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
#' \donttest{
#' # Perform Factor Forest
#' factor_forest(two_factor$data)}
#'
#' @author
#' # Authors of Factor Forest \cr
#' David Goretzko and Markus Buhner
#'
#' # Authors of \{latentFactoR\} \cr
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Hudson Golino <hfg9s@virginia.edu>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
#'
#' @references
#' Goretzko, D., & Buhner, M. (2022).
#' Factor retention using machine learning with ordinal data.
#' \emph{Applied Psychological Measurement}, 01466216221089345.
#'
#' Goretzko, D., & Buhner, M. (2020).
#' One model to rule them all? Using machine learning
#' algorithms to determine the number of factors
#' in exploratory factor analysis.
#' \emph{Psychological Methods}, \emph{25}(6), 776-786.
#'
#' Ruscio, J., & Roche, B. (2012).
#' Determining the number of factors to retain in an exploratory factor
#' analysis using comparison data of known factorial structure.
#' \emph{Psychological Assessment}, \emph{24}(2), 282-292.
#'
#' @export
#'
# Factor Forest
# Updated 20.11.2026
factor_forest <- function(
    data, sample_size,
    maximum_factors = 8,
    PA_correlation = c("cor", "poly", "tet"),
    seed = NULL
)
{

  # Set seed (reproducible results)
  seed_token <- lf_seed_start(seed)
  on.exit(lf_seed_end(seed_token), add = TRUE)


  # Check for correlation type
  PA_correlation <- match.arg(PA_correlation)

  # Check for appropriate data
  object_error(data, c("matrix", "data.frame", "array"));
  sink <- apply(data, 2, type_error, expected_type = "numeric");

  # Ensure data is matrix
  data <- as.matrix(data)

  # Correlation matrices are not supported
  if(is_symmetric(data)){
    stop(
      paste(
        "Input for 'data' is a correlation (symmetric) matrix. `factor_forest`",
        "requires data because the comparison data approach (a feature of the",
        "model) is computed from the data"
      ), call. = FALSE
    )
  }

  # Check for appropriate maximum factors
  type_error(maximum_factors, "numeric");
  length_error(maximum_factors, 1);
  range_error(maximum_factors, c(1, Inf));

  # Check for variable names
  if(is.null(colnames(data))){
    colnames(data) <- variable_names(ncol(data))
  }

  # Compute features
  features <- factor_forest_features(
    data, maximum_factors = maximum_factors, PA_correlation = PA_correlation
  )

  # Ensure features are in the order the model was trained on
  # (the model matches features by position)
  if(!identical(names(features), factor_forest_feature_names)){
    stop("Factor Forest features do not match the model's training features", call. = FALSE)
  }

  # Predict class probabilities (uses the best early-stopping iteration)
  probabilities <- as.vector(
    suppressWarnings(predict(factor_forest_model(), as.matrix(features)))
  )

  # Organize output
  names(probabilities) <- seq_along(probabilities)

  # Set up return list
  results <- list(
    dimensions = which.max(probabilities)[[1]],
    probabilities = round(probabilities, 4)
  )

  # Add class
  class(results) <- "lf_ff"

  # Return results
  return(results)

}

#' @noRd
# Cache so the model is only read once per R session
# Updated 20.11.2026
factor_forest_cache <- new.env(parent = emptyenv())

#' @noRd
# Load (and cache) the pre-trained model
# (`tunedxgb.json` is the booster of Goretzko & Buhner, 2020, re-saved in the
# JSON format that current versions of {xgboost} can read; see inst/COPYRIGHTS)
# Updated 20.11.2026
factor_forest_model <- function()
{

  # Read once
  if(is.null(factor_forest_cache$model)){

    # Locate model
    path <- system.file("extdata", "tunedxgb.json", package = "latentFactoR")
    if(!nzchar(path)){
      stop("Could not find the Factor Forest model (`tunedxgb.json`)", call. = FALSE)
    }

    # Load model
    factor_forest_cache$model <- xgboost::xgb.load(path)

  }

  # Return model
  return(factor_forest_cache$model)

}

#' @noRd
# Feature names in the order the model was trained on
# (from the task of the original {mlr} model, `tunedxgb.rds`)
# Updated 20.11.2026
factor_forest_feature_names <- c(
  "N", "p", "eiggreater1", "releig1", "releig2", "releig3", "eiggreater07",
  "sdeigval", "var50", "var75", "onenorm", "frobnorm", "maxnorm",
  "avgcor", "specnorm", "smlcor", "avgcom", "det", "KMO", "Gini",
  "Kolm", "pa_solution", "ekc", "cd",
  paste0("eigval", 1:80), paste0("fa_eigval", 1:80)
)

#' @noRd
# Pad (with -1000, as in training) or truncate to the model's expected length
# Updated 20.11.2026
factor_forest_pad <- function(x, length = 80)
{

  # Initialize padded
  padded <- rep(-1000, length)

  # Fill
  n <- min(base::length(x), length)
  padded[seq_len(n)] <- x[seq_len(n)]

  # Return padded
  return(padded)

}

#' @noRd
# Parallel analysis (common factor model), as `psych::fa.parallel(fa = "fa")`
# The one-factor minimum residual eigenvalues of the correlation matrix are
# compared to the `quant` quantile of those from random normal data
# (the reference only depends on the sample size and number of variables)
# Returns the number of factors and the empirical eigenvalues
# Updated 20.11.2026
factor_forest_parallel <- function(correlation, sample_size, n.iter = 20, quant = 0.95)
{

  # Number of variables
  variables <- ncol(correlation)

  # Minimum residual eigenvalues (one factor, unrotated)
  fa_values <- function(R, SMC = TRUE){
    suppressWarnings(suppressMessages(
      psych::fa(
        R, nfactors = 1, rotate = "none", fm = "minres",
        SMC = SMC, warnings = FALSE
      )$values
    ))
  }

  # Empirical eigenvalues
  empirical <- fa_values(correlation)

  # Eigenvalues of random normal data
  random <- vapply(seq_len(n.iter), function(i){
    fa_values(
      cor(matrix(rnorm_ziggurat(sample_size * variables), nrow = sample_size)),
      SMC = FALSE
    )
  }, numeric(variables))

  # Quantiles of random eigenvalues
  reference <- apply(random, 1, quantile, probs = quant)

  # Number of consecutive eigenvalues greater than the reference
  above <- empirical > reference
  factors <- swiftelse(all(above), variables, which(!above)[1] - 1)

  # Return results
  return(list(factors = factors, eigenvalues = empirical))

}

#' @noRd
# Compute the model's features (one-row data frame)
# Feature definitions follow the authors' evaluation code (column order is
# the order the model was trained on)
# Updated 20.11.2026
factor_forest_features <- function(data, maximum_factors = 8, PA_correlation = "cor")
{

  # Obtain dimensions
  N <- nrow(data)
  p <- ncol(data)

  # Correlation matrix
  dat_cor <- EGAnet::auto.correlate(data, verbose = FALSE)

  # Eigenvalues
  eigval <- eigen(dat_cor, symmetric = TRUE, only.values = TRUE)$values
  vareig <- cumsum(eigval) / p

  # Empirical Kaiser criterion (as in training: eigenvalues before the first
  # one not exceeding the reference)
  reference <- numeric(p)
  for(i in seq_len(p)){
    reference[i] <- max(
      ((1 + sqrt(p / N))^2) * (p - sum(reference)) / (p - i + 1), 1
    )
  }
  ekc <- swiftelse(
    any(eigval <= reference), which(eigval <= reference)[1] - 1, p
  )

  # Eigenvalue based features
  eiggreater1 <- sum(eigval > 1)
  releig1 <- eigval[1] / p
  releig2 <- sum(eigval[1:2]) / p
  releig3 <- sum(eigval[1:3]) / p
  eiggreater07 <- sum(eigval > 0.7)
  sdeigval <- sd(eigval)
  var50 <- min(which(vareig > 0.50))
  var75 <- min(which(vareig > 0.75))

  # Matrix norm features
  onenorm <- norm(dat_cor, "O")
  frobnorm <- norm(dat_cor, "F")
  maxnorm <- norm(dat_cor - diag(p), "M")
  avgcor <- sum(abs(dat_cor - diag(p))) / (p * (p - 1))
  specnorm <- sqrt(eigen(crossprod(dat_cor), symmetric = TRUE, only.values = TRUE)$values[1])

  smlcor <- sum(dat_cor <= 0.1)
  avgcom <- mean(suppressWarnings(psych::smc(dat_cor)))
  det <- det(dat_cor)

  # "Inequality" features
  # (the model was trained with these computed on `lower.tri(dat_cor)`, a
  # logical matrix, rather than on the correlations; kept so that the
  # features match training)
  KMO <- suppressWarnings(psych::KMO(dat_cor)$MSA)
  Gini <- ineq::ineq(lower.tri(dat_cor), type = "Gini")
  Kolm <- ineq::ineq(lower.tri(dat_cor), type = "Kolm")

  # Parallel analysis
  pa_correlation <- switch(
    PA_correlation,
    "cor" = cor(data, use = "pairwise.complete.obs"),
    "poly" = suppressWarnings(psych::polychoric(data)$rho),
    "tet" = suppressWarnings(psych::tetrachoric(data)$rho)
  )
  pa <- factor_forest_parallel(pa_correlation, N)
  pa_solution <- pa$factors

  # Pad (or truncate) eigenvalues to the model's expected length
  eigval <- factor_forest_pad(eigval)
  fa_eigval <- factor_forest_pad(pa$eigenvalues)
  names(eigval) <- paste0("eigval", seq_along(eigval))
  names(fa_eigval) <- paste0("fa_eigval", seq_along(fa_eigval))

  # Comparison data approach
  cd <- EFA.Comp.Data(Data = data, F.Max = maximum_factors, use = "pairwise.complete.obs")

  # Return features
  return(
    cbind(
      data.frame(
        N = N, p = p, eiggreater1, releig1, releig2, releig3, eiggreater07,
        sdeigval, var50, var75, onenorm, frobnorm, maxnorm, avgcor, specnorm,
        smlcor, avgcom, det, KMO, Gini, Kolm, pa_solution, ekc, cd
      ),
      t(eigval), t(fa_eigval)
    )
  )

}
