#' Adds Local Dependence to Factor Model Data
#'
#' Adds local dependence to simulated data from \code{\link[latentFactoR]{simulate_factors}}. 
#' See examples to get started
#' 
#' @param lf_object Data object from \code{\link[latentFactoR]{simulate_factors}}
#' 
#' @param method Character (length = 1).
#' Method to generate local dependence between variables.
#' Only \code{"correlate_residuals"} at the moment.
#' Future developments will include minor factor
#' and threshold-shift methods. Description of methods:
#' 
#' \itemize{
#' 
#' \item{\code{"correlate_residuals"}}
#' {Adds residuals directly to the population correlation matrix
#' prior to data generation (uses population correlation matrix
#' from \code{\link[latentFactoR]{simulate_factors}})}
#' 
#' \item{\code{"minor_factors"}}
#' {Coming soon...}
#' 
#' \item{\code{"threshold_shifts"}}
#' {Coming soon...} 
#' 
#' }
#' 
#' @param proportion_LD Numeric (length = 1 or \code{factors}).
#' Proportion of variables that should be locally dependent across all
#' or each factor. Accepts number of locally dependent values as well
#' 
#' @param proportion_LD_range Numeric (length = 2).
#' Range of proportion of variables that are randomly selected from
#' a random uniform distribution. Accepts number of locally dependent values as well
#' 
#' @param add_residuals Numeric (length = 1, \code{factors}, or total number of locally dependent variables).
#' Amount of residual to add to the population correlation matrix between two variables.
#' Only used when \code{method = "correlated_residuals"}. Magnitudes are drawn from
#' a random uniform distribution using +/- 0.05 of value input.
#' Can also be specified directly (same length as total number of locally dependent variables).
#' General effect sizes range from small (0.20), moderate (0.30), to large (0.40)
#' 
#' @param add_residuals_range Numeric (length = 2).
#' Range of the residuals to add to the correlation matrix are randomly selected from
#' a random uniform distribution
#' 
#' @param allow_multiple Boolean.
#' Whether a variable should be allowed to be locally dependent with
#' more than one other variable.
#' Defaults to \code{FALSE}.
#' Set to \code{TRUE} for more complex locally dependence patterns
#' 
#' @return Returns a list containing:
#' 
#' \item{correlated_residuals}{A data frame with the first two columns specifying
#' the variables that are locally dependent and the third column specifying the
#' magnitude of the added residual for each locally dependent pair}
#' 
#' \item{data}{Simulated data from the specified factor model}
#' 
#' \item{population_correlation}{Population correlation matrix with local dependence added}
#' 
#' \item{original_correlation}{Original population correlation matrix \emph{before}
#' local dependence was added}
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
#' # Add local dependence
#' two_factor_LD <- add_local_dependence(
#'   lf_object = two_factor,
#'   proportion_LD = 0.25,
#'   add_residuals = 0.20,
#'   allow_multiple = FALSE 
#' )
#' 
#' # Randomly vary proportions
#' two_factor_LD <- add_local_dependence(
#'   lf_object = two_factor,
#'   proportion_LD_range = c(0.10, 0.50),
#'   add_residuals = 0.20,
#'   allow_multiple = FALSE 
#' )
#' 
#' # Randomly vary residuals
#' two_factor_LD <- add_local_dependence(
#'   lf_object = two_factor,
#'   proportion_LD = 0.25,
#'   add_residuals_range = c(0.20, 0.40),
#'   allow_multiple = FALSE 
#' )
#' 
#' # Randomly vary proportions, residuals, and allow multiple 
#' two_factor_LD <- add_local_dependence(
#'   lf_object = two_factor,
#'   proportion_LD_range = c(0.10, 0.50),
#'   add_residuals_range = c(0.20, 0.40),
#'   allow_multiple = TRUE
#' )
#' 
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Hudson Golino <hfg9s@virginia.edu>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
#' 
#' @references
#' Christensen, A. P., Garrido, L. E., & Golino, H. (2022).
#' Unique variable analysis: A network psychometrics method to detect local dependence.
#' \emph{PsyArXiv}
#'
#' @export
#'
# Add local dependence to simulated data
# Updated 17.09.2022
add_local_dependence <- function(
    lf_object,
    method = c(
      "correlate_residuals",
      "minor_factors",
      "threshold_shifts"
    ),
    proportion_LD, proportion_LD_range = NULL,
    add_residuals = NULL, add_residuals_range = NULL,
    allow_multiple = FALSE
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
  
  # Check for missing method
  if(missing(method)){
    method <- "correlate_residuals"
  }else{method <- tolower(match.arg(method))}
  
  # Check for proportion local dependence range
  if(!is.null(proportion_LD_range)){
    type_error(proportion_LD_range, "numeric") # object type error
    length_error(proportion_LD_range, 2) # object length error
    
    # Check for number of variables in range
    if(any(proportion_LD_range >= 1)){
      
      # Target values
      target_LD <- which(proportion_LD_range >= 1)
      
      # Ensure proportions
      proportion_LD_range[target_LD] <-
        proportion_LD_range[target_LD] / parameters$variables[target_LD]
      
    }
    
    # Check for error in range
    range_error(proportion_LD_range, c(0, 1)) # object range error
    proportion_LD <- runif(
      parameters$factors,
      min = min(proportion_LD_range),
      max = max(proportion_LD_range)
    )
  }
  
  # Ensure appropriate types
  type_error(proportion_LD, "numeric");
  
  # Ensure appropriate lengths
  length_error(proportion_LD, c(1, parameters$factors));
  
  # Set proportions
  if(length(proportion_LD) == 1){
    proportion_LD <- rep(proportion_LD, parameters$factors)
  }
  
  # Convert local dependence proportions to proportions
  if(any(proportion_LD >= 1)){
    
    # Target values
    target_LD <- which(proportion_LD >= 1)
    
    # Ensure proportions
    proportion_LD[target_LD] <-
      proportion_LD[target_LD] / parameters$variables[target_LD]
    
  }
  
  # Ensure appropriate ranges
  range_error(proportion_LD, c(0, 1));
  
  # Add local dependence
  if(method == "correlate_residuals"){
    
    # Obtain results
    results <- correlate_residuals(
      lf_object = lf_object,
      proportion_LD = proportion_LD,
      allow_multiple = allow_multiple,
      add_residuals = add_residuals,
      add_residuals_range = add_residuals_range
    )
    
  }else{
    
    stop(
      paste0(
        "'", method, "' is not available yet and will be implemented in the future.",
        "\nFor now, use `method = \"correlate_residuals\""
      )
    )
    
  }
  
  # Add class
  class(results) <- c(class(lf_object), "lf-ld")
  
  # Return results
  return(results)
  
}