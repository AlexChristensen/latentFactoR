#' Obtain Zipf's Distribution Parameters from Data
#'
#' Zipf's distribution is commonly found for text data. Closely related to the
#' Pareto and power-law distributions, the Zipf's distribution produces
#' highly skewed data. This function obtains the best fitting parameters
#' to Zipf's distribution 
#' 
#' 
#' @details The best parameters are optimized by minimizing the aboslute
#' difference between the original frequencies and the frequencies obtained
#' by the \emph{beta} and \emph{alpha} parameters in the following
#' formula (Piantadosi, 2014):
#' 
#' \emph{f(r) proportional to 1 / (r + beta)^alpha}
#' 
#' where \emph{f(r)} is the \emph{r}th most frequency,
#' \emph{r} is the rank-order of the data, \emph{beta}
#' is a shift in the rank (following Mandelbrot, 1953, 1962),
#' and \emph{alpha} is the power of the rank with greater
#' values suggesting greater differences between the largest
#' frequency to the next, and so forth.
#' 
#' @param data Numeric vector, matrix, or data frame.
#' Numeric data to determine Zipf's distribution parameters
#'
#' @return Returns a vector containing the estimated \code{beta} and
#' \code{alpha} parameters. Also contains \code{zipfs_rmse} which corresponds
#' to the root mean square error between frequencies based
#' on the parameter values estimated and the original data frequencies
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
#' # Transform data to Mandelbrot's Zipf's
#' two_factor_zipfs <- data_to_zipfs(
#'   lf_object = two_factor,
#'   beta = 2.7,
#'   alpha = 1
#' )
#' 
#' # Obtain Zipf's distribution parameters
#' obtain_zipfs_parameters(two_factor_zipfs$data)
#' 
#' @references
#' Mandelbrot, B. (1953).
#' An informational theory of the statistical structure of language.
#' \emph{Communication Theory}, \emph{84}, 486–502.
#' 
#' Mandelbrot, B. (1962).
#' On the theory of word frequencies and on related Markovian models of discourse.
#' \emph{Structure of Language and its Mathematical Aspects}, 190–219.
#' 
#' Piantadosi, S. T. (2014).
#' Zipf’s word frequency law in natural language: A critical review and future directions.
#' \emph{Psychonomic Bulletin & Review}, \emph{21}(5), 1112-1130.
#' 
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Hudson Golino <hfg9s@virginia.edu>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
#' 
#' @export
#'
# Obtain Zipf's distribution parameters
# Updated 26.09.2022
obtain_zipfs_parameters <- function(data)
{
  
  # Check if data is vector
  if(is.vector(data)){
    
    # Obtain frequencies
    frequencies <- data
    
    # Obtain Zipf's values
    zipfs <- frequencies / length(data)
    
  }else{
    
    # Ensure data is matrix
    data <- as.matrix(data)
    
    # Obtain frequencies
    frequencies <- as.vector(data)
    
    # Obtain Zipf's values
    zipfs <- frequencies / nrow(data)
    
  }
  
  # Set nearest decimal
  digit <- nearest_decimal(zipfs)
  zipfs <- round(
    zipfs, digits = digit
  )
  
  # Non-zero zipfs
  non_zero_zipfs <- zipfs != 0
  
  # Assume rank order based on frequencies
  rank_order <- rank(-frequencies)
  
  # For large datasets, remove frequencies
  rm(frequencies); gc()
  
  # Equation to solve
  # zipfs = 1 / (rank_order + beta)^alpha
  
  # Set beta (assume alpha = 1)
  beta <- estimate_beta(
    beta_sequence = seq(0, 100, 1),
    zipfs = zipfs,
    non_zero_zipfs = non_zero_zipfs,
    rank_order = rank_order,
    digit = digit
  )
  
  # Fine-tune beta
  est_beta <- seq(
    ifelse(beta == 0, 0, beta - 1),
    ifelse(beta == 0, 2, beta + 1),
    0.1
  )
  
  # Set beta (assume alpha = 1)
  beta <- estimate_beta(
    beta_sequence = est_beta,
    zipfs = zipfs,
    non_zero_zipfs = non_zero_zipfs,
    rank_order = rank_order,
    digit = digit
  )
  
  # Further fine-tune beta
  est_beta <- seq(
    ifelse(beta == 0, 0, beta - 0.1),
    ifelse(beta == 0, 2, beta + 0.1),
    0.01
  )
  
  # Set beta (assume alpha = 1)
  beta <- estimate_beta(
    beta_sequence = est_beta,
    zipfs = zipfs,
    non_zero_zipfs = non_zero_zipfs,
    rank_order = rank_order,
    digit = digit
  )
  
  # New line
  cat("\n")
  
  # Set alpha
  alpha <- estimate_alpha(
    beta = beta,
    alpha_sequence = seq(0, 20, 1),
    zipfs = zipfs,
    non_zero_zipfs = non_zero_zipfs,
    rank_order = rank_order,
    digit = digit
  )
  
  # Fine-tune alpha
  est_alpha <- seq(
    ifelse(alpha == 0, 0, alpha - 1),
    ifelse(alpha == 0, 2, alpha + 1),
    0.1
  )
  
  # Set alpha
  alpha <- estimate_alpha(
    beta = beta,
    alpha_sequence = est_alpha,
    zipfs = zipfs,
    non_zero_zipfs = non_zero_zipfs,
    rank_order = rank_order,
    digit = digit
  )
  
  # Further fine-tune alpha
  est_alpha <- seq(
    ifelse(alpha == 0, 0, alpha - 0.1),
    ifelse(alpha == 0, 2, alpha + 0.1),
    0.01
  )
  
  # Set alpha
  alpha <- estimate_alpha(
    beta = beta,
    alpha_sequence = est_alpha,
    zipfs = zipfs,
    non_zero_zipfs = non_zero_zipfs,
    rank_order = rank_order,
    digit = digit
  )
  
  # New line
  cat("\n")
  
  # Final values
  final_values <- round(1 / (rank_order + beta)^alpha, digit)
  
  # Compute root mean square error
  rmse <- sqrt(
    mean(
      (
        zipfs[non_zero_zipfs] - final_values[non_zero_zipfs]
      )^2, na.rm = TRUE
    )
  )
  
  # Return parameters
  return(
    zapsmall(c(beta = beta, alpha = alpha, zipfs_rmse = rmse))
  )
  
}
