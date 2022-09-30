#' Estimate Number of Dimensions using Next Eigenvalue Sufficiency Test
#'
#' Estimates the number of dimensions in data using NEST (Achim, 2017).
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
#' @param iterations Numeric (length = 1).
#' Number of iterations to estimate rank.
#' Defaults to \code{1000}
#' 
#' @param maximum_iterations Numeric (length = 1).
#' Maximum umber of iterations to obtain convergence
#' of eigenvalues.
#' Defaults to \code{500}
#' 
#' @param alpha Numeric (length = 1).
#' Significance level for determine sufficient eigenvalues.
#' Defaults to \code{0.05}
#' 
#' @param convergence Numeric (length = 1).
#' Value necessary to be less than or equal to 
#' when establishing convergence of eigenvalues
#' 
#' @return Returns a list containing:
#' 
#' \item{dimensions}{Number of dimensions identified}
#' 
#' \item{loadings}{Loading matrix}
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
#' \dontrun{
#' # Perform NEST
#' NEST(two_factor$data)}
#' 
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Hudson Golino <hfg9s@virginia.edu>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
#' 
#' @references
#' Achim, A. (2017).
#' Testing the number of required dimensions in
#' exploratory factor analysis.
#' \emph{The Quantitative Methods for Psychology}, \emph{13}(1), 64–74.
#' 
#' Brandenburg, N., & Papenberg, M. (2022).
#' Reassessment of innovative methods to determine the number 
#' of factors: A simulation-Based comparison of Exploratory 
#' Graph Analysis and Next Eigenvalue Sufficiency Test.
#' \emph{Psychological Methods}. 
#'
#' @export
#'
# Next Eigenvalue Sufficiency Test
# Updated 30.09.2022
NEST <- function(
    data, sample_size,
    iterations = 1000,
    maximum_iterations = 500,
    alpha = 0.05,
    convergence = 0.00001
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
  
  # Check for appropriate input for the rest of the arguments
  ## Type errors
  type_error(iterations, "numeric"); type_error(maximum_iterations, "numeric");
  type_error(alpha, "numeric"); type_error(convergence, "numeric");
  
  ## Length errors
  length_error(iterations, 1); length_error(maximum_iterations, 1);
  length_error(alpha, 1); length_error(convergence, 1);
  
  ## Range errors
  range_error(iterations, c(1, Inf)); range_error(maximum_iterations, c(1, Inf));
  range_error(alpha, c(0, 1)); range_error(convergence, c(0, 1));
  
  # Obtain eigenvalues
  eigenvalues <- eigen(correlation, symmetric = TRUE, only.values = TRUE)$values
  
  # Factor Limit
  factor_limit <- floor(0.80 * ncol(correlation))
  
  # Loop through number of factors
  for(factors in 0:factor_limit){
    
    # Increase factors by 1
    factors_1 <- factors + 1
    
    # Rank
    rank <- rep(1, factors_1)
    
    # Set up model
    if(factors == 0){
      model <- diag(rep(1, ncol(correlation)))
    }else{
      
      # Copy correlation matrix
      R <- correlation
      
      # Obtain diagonal
      diagonal <- diag(R)
      
      # Loop through iterations
      for(i in 1:maximum_iterations){
        
        # Obtain eigenvalues and eigenvectors
        eigens <- eigen(R, symmetric = TRUE)
        
        # Check for unidimensional structure
        if(factors == 1){
          
          # Check LD
          LD <- eigens$vectors[,1] * sqrt(eigens$values[1])
          
          # Obtain communalities
          communalities <- LD^2
          
        }else{
          
          # Check eigenvalues and eigenvectors across factors
          current_factor <- factors
          
          # Loop through factors
          while(eigens$values[current_factor] <= 0){
            current_factor <- current_factor - 1
          }
          
          # Check LD
          LD <- eigens$vectors[,1:current_factor] %*%
            diag(sqrt(eigens$values[1:current_factor]))
          
          # Obtain communalities
          communalities <- rowSums(LD^2)
          
        }
        
        # Check for communalities greater than 1
        if(max(communalities) > 1){
          R <- R + diag(diagonal - diag(R))
          warning(
            paste(
              "Communalities greater than 1 with",
              factors, "factors. Stopping estimation..."
            )
          )
          break
        }
        
        # Compute absolute difference between diagonal and communalities
        difference <- max(abs(diagonal - communalities))
        
        # Break if difference is less than convergence
        if(difference < convergence){
          break
        }
        
        # Re-set correlation matrix with communalities
        R <- R + diag(communalities - diagonal)
        
        # Re-set diagonal
        diagonal <- communalities
        
      }
      
      # Check for whether maximum iterations were reached
      if(i >= maximum_iterations){
        warning(
          paste(
            "Convergence not found with",
            factors, "factors..."
          )
        )
      }
      
      # Set the model
      model <- t(
        cbind(
          LD, diag(sqrt(1 - communalities))
        )
      )
      
    }
    
    # Sum of factors and variables
    factors_variables <- factors + ncol(data)
    
    # Compute rank
    for(j in 1:iterations){
      
      # Generate data
      random_data <- matrix(
        rnorm(sample_size * factors_variables),
        nrow = sample_size, ncol = factors_variables
      )
      
      # Multiply data by model
      new_data <- random_data %*% model
      
      # Compute new correlations
      new_correlation <- cor(new_data)
      
      # Compute eigenvalues
      eigens <- eigen(
        new_correlation,
        symmetric = TRUE,
        only.values = TRUE
      )
      
      # Compute rank
      rank <- rank + (
        eigens$values[1:factors_1] >=
          eigenvalues[1:factors_1]
      )
      
    }
    
    # Check for rank significance
    if(rank[factors_1] > alpha * (iterations + 1)){
      
      # Obtain loadings
      loadings <- as.matrix(
        data.frame(
          t(model[1:factors,])
        )
      )
      
      # Remove loading names
      loadings <- unname(loadings)
      
      # Break out of loop
      break
    }
    
  }
  
  
  # Set up results list
  results <- list(
    dimensions = factors,
    loadings = loadings
  )
  
  # Return results list
  return(results)
  
}
