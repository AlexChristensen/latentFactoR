#' Transforms Data to Zipf's Distribution
#'
#' Zipf's distribution is commonly found for text data. Closely related to the
#' Pareto and power-law distributions, the Zipf's distribution produces
#' highly skewed data. This transformation is intended to mirror the data
#' generating process of Zipf's law seen in semantic network and topic
#' modeling data.
#' 
#' @details The formula used to transform data is (Piantadosi, 2014):
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
#' The function will transform continuous data output from \code{\link[latentFactoR]{simulate_factors}}. 
#' See examples to get started
#' 
#' @param lf_object Data object from \code{\link[latentFactoR]{simulate_factors}}
#' 
#' @param beta Numeric (length = 1).
#' Sets the shift in rank.
#' Defaults to \code{2.7}
#' 
#' @param alpha Numeric (length = 1).
#' Sets the power of the rank.
#' Defaults to \code{1}
#' 
#' @param dichotomous Boolean (length = 1).
#' Whether data should be dichotomized rather
#' than frequencies (e.g., semantic network analysis).
#' Defaults to \code{FALSE}
#' 
#' @return Returns a list containing:
#' 
#' \item{data}{Simulated data that has been transform to follow Zipf's distribution}
#' 
#' \item{RMSE}{A vector of root mean square errors for transformed data and data
#' assumed to follow theoretical Zipf's distribution and Spearman's correlation
#' matrix of the transformed data compared to the original population correlation
#' matrix}
#' 
#' \item{spearman_correlation}{Spearman's correlation matrix of the transformed data}
#' 
#' \item{original_correlation}{Original population correlation matrix \emph{before}
#' the data were transformed}
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
#' # Transform data to Mandelbrot's Zipf's
#' two_factor_zipfs <- data_to_zipfs(
#'   lf_object = two_factor,
#'   beta = 2.7,
#'   alpha = 1
#' )
#' 
#' # Transform data to Mandelbrot's Zipf's (dichotomous)
#' two_factor_zipfs_binary <- data_to_zipfs(
#'   lf_object = two_factor,
#'   beta = 2.7,
#'   alpha = 1,
#'   dichotomous = TRUE
#' )
#' 
#' # Transform data to Pareto distribution
#' two_factor_pareto <- data_to_zipfs(
#'   lf_object = two_factor,
#'   beta = 2.7,
#'   alpha = 0.4 # gets very close to Pareto
#' )
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
#' Zipf, G. (1936).
#' \emph{The psychobiology of language}.
#' London, UK: Routledge.
#' 
#' Zipf, G. (1949).
#' \emph{Human behavior and the principle of least effort}. 
#' New York, NY: Addison-Wesley.
#' 
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Hudson Golino <hfg9s@virginia.edu>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
#' 
#' @importFrom stats cor
#' 
#' @export
#'
# Transform data to Zipf's distribution
# Updated 23.09.2022
data_to_zipfs <- function(
    lf_object,
    beta = 2.7,
    alpha = 1,
    dichotomous = FALSE
)
{
  
  # Check for appropriate class
  if(!is(lf_object, "lf-simulate")){
    
    # Produce error
    stop(
      paste(
        "`lf_object` input is not class \"lf-simulate\" from the `simulate_factors` function.",
        "\n\nInput class(es) of current `lf_object`:", 
        paste0("\"", class(lf_object), "\"", collapse = ", "),
        "\n\nUse `simulate_factors` to generate your data to input into this function"
      )
    )
    
  }
  
  # Obtain parameters from simulated data
  parameters <- lf_object$parameters
  
  # Ensure appropriate input
  type_error(beta, "numeric"); length_error(beta, 1);
  type_error(alpha, "numeric"); length_error(alpha, 1);
  range_error(alpha, c(0.001, Inf))
  
  # Obtain original data
  original_data <- lf_object$data
  
  # Obtain correlations
  correlations <- lf_object$population_correlation
  diag(correlations) <- 0 # ensures finding maximum later
  
  # Reverse values and obtain ranked data
  rank_order <- rank(-original_data)
  
  # Transformation based on Mandelbrot
  zipfs <- 1 / (rank_order + beta)^alpha
  
  # Estimate frequencies
  frequencies <- round(zipfs * nrow(original_data))
  
  # Initialize and populate data
  data <- original_data; data[] <- frequencies;
  
  # Check for zero sum totals
  sum_totals <- colSums(data)
  
  # Check for zeros
  if(any(sum_totals == 0)){
    
    # Which?
    targets <- which(sum_totals == 0)
    
    # Loop through
    for(i in targets){
      
      # Find max correlation value
      max_target <- which.max(abs(correlations[i,]))
      
      # Find non-zero values
      non_zero <- which(data[,max_target] != 0)
      
      # Sample from non_zero values
      index <- non_zero[sample(1:length(non_zero), 1)]
      
      # Assign 1 to target index
      data[index, i] <- 1
      
    }
    
  }
  
  # Determine minimum
  minimum <- min(data)
  
  # Check for setting variables to zero
  if(minimum > 0){
    
    # Substract minimum from data
    data <- data - minimum
    
  }
  
  # Check for whether data should be binarized
  if(isTRUE(dichotomous)){
    data[data != 0] <- 1
  }
  
  # Check parameters
  frequencies <- as.vector(data)
  
  # Obtain Zipf's values
  zipfs <- frequencies / nrow(data)
  
  # Set nearest decimal
  digit <- nearest_decimal(zipfs)
  zipfs <- round(zipfs, digits = digit)
  
  # Non-zero zipfs
  non_zero_zipfs <- zipfs != 0
  
  # Assume rank order based on frequencies
  rank_order <- rank(-frequencies)
  
  # Simulated values
  simulated_values <- round(
    1 / (rank_order + beta)^alpha,
    digit
  )
  
  # Compute RMSE for Zipfs parameters
  zipfs_rmse <- sqrt(
    mean(
      (
        zipfs[non_zero_zipfs] - simulated_values[non_zero_zipfs]
      )^2, na.rm = TRUE
    )
  )
  
  # Obtain Spearman correlation matrix
  spearman_correlation <- cor(data, method = "spearman")
  
  # Compute RMSE for correlation matrix
  correlation_rmse <- sqrt(
    mean(
      (spearman_correlation - lf_object$population_correlation)^2, 
      na.rm = TRUE
    )
  )
  
  # Populate results
  results <- list(
    data = data,
    RMSE = c(
      zipfs = round(zipfs_rmse, 4),
      correlation = round(correlation_rmse, 4)
    ),
    spearman_correlation = spearman_correlation,
    original_correlation = lf_object$population_correlation,
    original_results = lf_object
  )
  
  # Add class
  class(results) <- c(class(lf_object), "lf-zipfs")
  
  # Return results
  return(results)
  
}