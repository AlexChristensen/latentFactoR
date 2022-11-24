#' Estimate Number of Dimensions using Empirical Kaiser Criterion
#'
#' Estimates the number of dimensions in data using 
#' Empirical Kaiser Criterion (Braeken & Van Assen, 2017).
#' See examples to get started
#' 
#' @param data Matrix or data frame.
#' Either a dataset with all numeric values
#' (rows = cases, columns = variables) or
#' a symmetric correlation matrix
#' 
#' @param sample_size Numeric (length = 1).
#' If input into \code{data} is a correlation matrix,
#' then specifying the sample size is required
#' 
#' @return Returns a list containing:
#' 
#' \item{dimensions}{Number of dimensions identified}
#' 
#' \item{eigenvalues}{Eigenvalues}
#' 
#' \item{reference}{Reference values compared against eigenvalues}
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
#' # Perform Empirical Kaiser Criterion
#' EKC(two_factor$data)
#' 
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Hudson Golino <hfg9s@virginia.edu>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
#' 
#' @references
#' Braeken, J., & Van Assen, M. A. (2017).
#' An empirical Kaiser criterion.
#' \emph{Psychological Methods}, \emph{22}(3), 450–466.
#'
#' @export
#'
# Empirical Kaiser Criterion
# Updated 24.11.2022
EKC <- function(
    data, sample_size
)
{
  
  # Check for appropriate data
  object_error(data, c("matrix", "data.frame", "array"));
  sink <- apply(data, 2, type_error, expected_type = "numeric");
  
  # Ensure data is matrix
  data <- as.matrix(data)
  
  # Check for variable names
  if(is.null(colnames(data))){
    colnames(data) <- paste0("V", 1:ncol(data))
  }
  
  # Obtain correlation matrix (if not already)
  if(!isSymmetric(data)){
    
    # Compute correlations
    correlation <- qgraph::cor_auto(
      data, forcePD = TRUE, verbose = FALSE
    )
    
    # Set sample size
    sample_size <- nrow(data)
    
  }else{
    
    # Set data as correlations
    correlation <- data
    
    # Check for sample size
    if(missing(sample_size)){
      stop("Input for 'sample_size' is required when input for 'data' is a correlation (symmetric) matrix.")
    }
    
  }
  
  # Check for appropriate sample size
  type_error(sample_size, "numeric"); length_error(sample_size, 1);
  range_error(sample_size, c(2, Inf))
  
  # Obtain number of variables
  variables <- ncol(correlation)
  
  # Obtain eigenvalues
  eigenvalues <- eigen(correlation)$values
  
  # Set null model
  l_up <- (1 + sqrt(variables / sample_size)^2)
  
  # Set eigenvalue criterion
  V <- c(0, eigenvalues[-variables])
  
  # Set order of variables
  W <- variables:1
  
  # Initialize reference criterion
  reference <- sapply(
    (variables - V) / W * l_up, function(x){
      max(x, 1)
    }
  )
  
  # Identify last eigenvalue greater than reference
  dimensions <- which(eigenvalues < reference)[1] - 1

  # Set up results list
  results <- list(
    dimensions = dimensions,
    eigenvalues = eigenvalues,
    reference = reference
  )
  
  # Return results
  return(results)
  
}
