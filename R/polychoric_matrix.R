#' Computes polychoric correlations
#'
#' @param data Matrix or data frame.
#' A dataset with all ordinal values
#' (rows = cases, columns = variables)
#' 
#' @return Returns a polychoric correlation matrix
#' 
#' @examples
#' # Generate polytomous data
#' two_factor_polytomous <- simulate_factors(
#'   factors = 2, # factors = 2
#'   variables = 6, # variables per factor = 6
#'   loadings = 0.55, # loadings between = 0.45 to 0.65
#'   cross_loadings = 0.05, # cross-loadings N(0, 0.05)
#'   correlations = 0.30, # correlation between factors = 0.30
#'   sample_size = 1000, # number of cases = 1000
#'   variable_categories = 5 # polytomous data
#' )
#' 
#' # Compute polychoric correlation matrix
#' correlations <- polychoric_matrix(two_factor_polytomous$data)
#' 
#' 
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com> with assistance from GPT-4
#'
#' @export
#'
# Compute polychoric correlation matrix
# Updated 03.04.2023
polychoric_matrix <- function(data)
{
  
  # Ensure data is an integer matrix
  data <- apply(as.matrix(data), 2, as.integer)
  
  # Call from C
  correlations <- .Call(
    "r_polychoric_correlation_matrix",
    data,
    PACKAGE = "latentFactoR"
  )

  # Return
  return(correlations)
  
  
}
