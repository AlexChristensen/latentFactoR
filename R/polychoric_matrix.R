#' Computes Polychoric Correlations
#' 
#' A fast implementation of polychoric correlations in C.
#' Uses the Beasley-Springer-Moro algorithm (Boro & Springer, 1977; Moro, 1995)
#' to estimate the inverse univariate normal CDF, the Drezner-Wesolosky 
#' approximation (Drezner & Wesolosky, 1990) to estimate the bivariate normal
#' CDF, and Brent's method (Brent, 2013) for optimization of rho
#'
#' @param data Matrix or data frame.
#' A dataset with all ordinal values
#' (rows = cases, columns = variables)
#' 
#' @param na_data Character (length = 1).
#' How should missing data be handled?
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"pairwise"}}
#' {Computes correlation for all available cases between
#' two variables}
#' 
#' \item{\code{"listwise"}}
#' {Computes correlation for all complete cases in the dataset}
#' 
#' }
#' 
#' @param empty_method Character (length = 1).
#' Method for empty cell correction.
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"none"}}
#' {Adds no value (\code{empty_value = "none"})
#' to the empirical joint frequency table between two variables}
#' 
#' \item{\code{"zero"}}
#' {Adds \code{empty_value} to the cells with zero
#' in the joint frequency table between two variables}
#' 
#' \item{\code{"all"}}
#' {Adds \code{empty_value} to all
#' in the joint frequency table between two variables}
#' 
#' }
#' 
#' @param empty_value Character (length = 1).
#' Value to add to the joint frequency table cells.
#' Accepts numeric values between 0 and 1 or
#' specific methods:
#' 
#' \itemize{
#' 
#' \item{\code{"none"}}
#' {Adds no value (\code{0}) to the empirical joint
#' frequency table between two variables}
#' 
#' \item{\code{"point_five"}}
#' {Adds \code{0.5} to the cells defined by \code{empty_method}}
#' 
#' \item{\code{"one_over"}}
#' {Adds \code{1 / n} where \code{n} equals the number of cells
#' based on \code{empty_method}. For \code{empty_method = "zero"},
#' }
#' 
#' }
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
#' # Obtain data (ensure matrix)
#' simulated_data <- as.matrix(two_factor_polytomous$data)
#' 
#' # Compute polychoric correlation matrix
#' correlations <- polychoric_matrix(simulated_data)
#' 
#' # Randomly assign missing data
#' simulated_data[sample(1:length(simulated_data), 1000)] <- NA
#' 
#' # Compute polychoric correlation matrix with pairwise method
#' na_correlations <- polychoric_matrix(
#'   simulated_data, na_data = "pairwise"
#'  )
#' 
#' @references 
#' Beasley, J. D., & Springer, S. G. (1977).
#' Algorithm AS 111: The percentage points of the normal distribution.
#' \emph{Journal of the Royal Statistical Society. Series C (Applied Statistics)}, \emph{26}(1), 118-121.
#' 
#' Brent, R. P. (2013). 
#' Algorithms for minimization without derivatives.
#' Mineola, NY: Dover Publications, Inc.
#' 
#' Drezner, Z., & Wesolowsky, G. O. (1990).
#' On the computation of the bivariate normal integral.
#' \emph{Journal of Statistical Computation and Simulation}, \emph{35}(1-2), 101-107.
#' 
#' Moro, B. (1995).
#' The full monte.
#' \emph{Risk 8 (February)}, 57-58.
#' 
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com> with assistance from GPT-4
#'
#' @export
#'
# Compute polychoric correlation matrix
# Updated 31.05.2023
polychoric_matrix <- function(
    data, na_data = c("pairwise", "listwise"),
    empty_method = c("none", "zero", "all"),
    empty_value = c("none", "point_five", "one_over")
)
{
  
  # Check for 'missing' argument
  if(missing(na_data)){
    na_data <- "pairwise"
  }else{
    na_data <- tolower(match.arg(na_data))
  }
  
  # Ensure data is a matrix
  data <- as.matrix(data)
  
  # Check for missing data
  if(na_data == "pairwise"){
    data[is.na(data)] <- 99 # "pairwise" is performed in C
  }else if(na_data == "listwise"){
    data <- na.omit(data) # no performance difference with C
  }
  
  # Ensure data is an integer matrix
  data <- apply(data, 2, as.integer)
  
  # Check for 'empty_method' argument
  if(missing(empty_method)){
    empty_method <- "none"
  }else{
    empty_method <- tolower(match.arg(empty_method))
  }
  
  # Check for 'empty_value' argument
  if(missing(empty_value)){
    empty_value <- "none"
  }else{
    
    # Check for numeric 
    if(!is.character(empty_value)){
      
      # Ensure proper range
      range_error(empty_value, c(0, 1))
      
    }else{ # Character input
      empty_value <- tolower(match.arg(empty_value))
    }
    
  }
  
  # Set up 'empty_method' and 'empty_value' for C
  if(empty_method == "none"){
    empty_method <- 0 # Set no value
    empty_value <- 0 # Set no value
  }else{
    
    # Set 'empty_method'
    empty_method <- ifelse(empty_method == "zero", 1, 2)
    
    # Set 'empty_value'
    if(is.character(empty_value)){
      empty_value <- ifelse(empty_value == "point_five", 0.50, 2)
    }
    
  } 
  
  # Call from C
  correlations <- .Call(
    "r_polychoric_correlation_matrix",
    data,
    as.integer(empty_method),
    as.double(empty_value),
    PACKAGE = "latentFactoR"
  )
  
  # Check for variable names
  if(!is.null(colnames(data))){
    
    # Add names to rows and columns
    colnames(correlations) <- 
      row.names(correlations) <- 
      colnames(data)
  }
  
  # Check for estimation issue
  ## Number of variables
  variables <- ncol(correlations)
  ## Get indices of diagonal
  diagonal_index <- seq(1, variables * variables, variables + 1)
  ## Check for all zero off-diagonal elements
  if(all(correlations[-diagonal_index] == 0)){
    
    # Check for this issue with polychoric empty cells when using:
    # empty_method = "zero" or "all"
    # empty_value >= 0.10
    # No idea why < 0.10 works (usually < 0.05) but greater than doesn't
    warning("Correlations did not converge. All correlations are zero.")
    
  }

  # Return
  return(correlations)
  
  
}
