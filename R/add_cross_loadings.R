#' Adds (Substantial) Cross-loadings to \code{\link[latentFactoR]{simulate_factors}} Data
#'
#' Intended to add substantial cross-loadings to simulated data from \code{\link[latentFactoR]{simulate_factors}}. 
#' See examples to get started
#' 
#' @param lf_object Data object from \code{\link[latentFactoR]{simulate_factors}}
#' 
#' @param proportion_cross_loadings Numeric (length = 1 or \code{factors}).
#' Proportion of variables that should be cross-loaded randomly onto
#' one other factor. Accepts number of variables to 
#' cross-load onto one other factor as well
#' 
#' @param proportion_cross_loadings_range Numeric (length = 2).
#' Range of proportion of variables that should be cross-loaded randomly onto
#' one other factor. Accepts number of variables to 
#' cross-load onto one other factor as well
#' 
#' @param magnitude_cross_loadings Numeric (length = 1, \code{factors}, or total number of variables to cross-load across all factors).
#' The magnitude or size of the cross-loadings.
#' Must range between \code{-1} and \code{1}.
#' 
#' @param magnitude_cross_loadings_range Numeric (length = 2).
#' The range of the magnitude or size of the cross-loadings.
#' Defaults to \code{NULL}
#' 
#' @param leave_cross_loadings Boolean.
#' Should cross-loadings be kept?
#' Defaults to \code{FALSE}.
#' Convergence problems can arise if cross-loadings are kept,
#' so setting them to zero is the default. Only set to \code{TRUE}
#' with careful consideration of the structure. Make sure to perform
#' additional checks that the data are adequate
#' 
#' @return Returns a list containing the same parameters as the original
#' \code{lf_object} but with updated \code{data}, \code{population_correlation},
#' and \code{parameters} (specifically, \code{loadings} matrix). Also returns
#' original \code{lf_object} in \code{original_results}
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
#' # Add substantial cross-loadings
#' two_factor_CL <- add_cross_loadings(
#'   lf_object = two_factor,
#'   proportion_cross_loadings = 0.25,
#'   magnitude_cross_loadings = 0.35
#' )
#' 
#' # Randomly vary proportions
#' two_factor_CL <- add_cross_loadings(
#'   lf_object = two_factor,
#'   proportion_cross_loadings_range = c(0, 0.25),
#'   magnitude_cross_loadings = 0.35
#' )
#' 
#' # Randomly vary magnitudes
#' two_factor_CL <- add_cross_loadings(
#'   lf_object = two_factor,
#'   proportion_cross_loadings = 0.25,
#'   magnitude_cross_loadings_range = c(0.35, 0.45)
#' )
#' 
#' # Set number of cross-loadings per factor (rather than proportion)
#' two_factor_CL <- add_cross_loadings(
#'   lf_object = two_factor,
#'   proportion_cross_loadings = 2,
#'   magnitude_cross_loadings = 0.35
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
# Add substantial cross-loadings to simulated data
# Updated 24.10.2022
add_cross_loadings <- function(
    lf_object,
    proportion_cross_loadings,
    proportion_cross_loadings_range = NULL,
    magnitude_cross_loadings,
    magnitude_cross_loadings_range = NULL,
    leave_cross_loadings = FALSE
)
{
  
  # Check for appropriate class
  if(!is(lf_object, "lf_simulate")){
    
    # Produce error
    stop(
      paste(
        "`lf_object` input is not class \"lf_simulate\" from the `simulate_factors` function.",
        "\n\nInput class(es) of current `lf_object`:", 
        paste0("\"", class(lf_object), "\"", collapse = ", "),
        "\n\nUse `simulate_factors` to generate your data to input into this function"
      )
    )
    
  }
  
  # Obtain parameters from simulated data
  parameters <- lf_object$parameters
  
  # Check for proportion cross-loadings range
  if(!is.null(proportion_cross_loadings_range)){
    type_error(proportion_cross_loadings_range, "numeric") # object type error
    length_error(proportion_cross_loadings_range, 2) # object length error
    
    # Check for number of variables in range
    if(any(proportion_cross_loadings_range >= 1)){
      
      # Target values
      target_LD <- which(proportion_cross_loadings_range >= 1)
      
      # Ensure proportions
      proportion_cross_loadings_range[target_LD] <-
        proportion_cross_loadings_range[target_LD] / parameters$variables[target_LD]
      
    }
    
    # Check for error in range
    range_error(proportion_cross_loadings_range, c(0, 1)) # object range error
    proportion_cross_loadings <- runif(
      parameters$factors,
      min = min(proportion_cross_loadings_range),
      max = max(proportion_cross_loadings_range)
    )
  }
  
  # Ensure appropriate types
  type_error(proportion_cross_loadings, "numeric");
  
  # Ensure appropriate lengths
  length_error(proportion_cross_loadings, c(1, parameters$factors));
  
  # Set proportions
  if(length(proportion_cross_loadings) == 1){
    proportion_cross_loadings <- rep(proportion_cross_loadings, parameters$factors)
  }
  
  # Convert local dependence proportions to proportions
  if(any(proportion_cross_loadings >= 1)){
    
    # Target values
    target_LD <- which(proportion_cross_loadings >= 1)
    
    # Ensure proportions
    proportion_cross_loadings[target_LD] <-
      proportion_cross_loadings[target_LD] / parameters$variables[target_LD]
    
  }
  
  # Ensure appropriate ranges
  range_error(proportion_cross_loadings, c(0, 1));
  
  # Obtain integer values of proportions
  integer_cross_loadings <- round(proportion_cross_loadings * parameters$variables)
  
  # Check for magnitude cross-loadings range
  if(!is.null(magnitude_cross_loadings_range)){
    type_error(magnitude_cross_loadings_range, "numeric") # object type error
    length_error(magnitude_cross_loadings_range, 2) # object length error
    range_error(magnitude_cross_loadings_range, c(0, 1)) # object range error
    magnitude_cross_loadings <- runif(
      sum(integer_cross_loadings),
      min = min(magnitude_cross_loadings_range),
      max = max(magnitude_cross_loadings_range)
    )
  }
  
  # Ensure appropriate types
  type_error(magnitude_cross_loadings, "numeric");
  
  # Ensure appropriate lengths
  length_error(magnitude_cross_loadings, c(1, parameters$factors, sum(integer_cross_loadings)));
  
  # Ensure appropriate ranges
  range_error(magnitude_cross_loadings, c(-1, 1));
  
  # Set magnitudes
  if(length(magnitude_cross_loadings) == 1){
    magnitude_cross_loadings <- rep(magnitude_cross_loadings, sum(integer_cross_loadings))
  }
  
  # Set sequence of variables for each factor
  end_variables <- cumsum(parameters$variables)
  start_variables <- (end_variables + 1) - parameters$variables
  
  # Obtain loadings
  loadings <- parameters$loadings
  
  # Check whether original cross-loadings should remain
  if(!isTRUE(leave_cross_loadings)){
    
    for(i in 1:ncol(loadings)){
      
      # Set cross-loadings to zero
      loadings[
        start_variables[i]:end_variables[i],
        -i
      ] <- 0
      
    }
    
  }
  
  # Set sequence of cross-loadings for each factor
  end_cross_loadings <- cumsum(integer_cross_loadings)
  start_cross_loadings <- (end_cross_loadings + 1) - integer_cross_loadings
  
  # Make cross-loading list (easier to manage in loop)
  cross_loading_list <- lapply(seq_along(start_cross_loadings), function(i){
    
    # Return each variables cross-loading value
    magnitude_cross_loadings[
      start_cross_loadings[i]:end_cross_loadings[i]
    ]
    
  })
  
  # Obtain factor correlations
  factor_correlations <- parameters$factor_correlations
  
  # Set cross-loadings on the factors
  for(i in 1:ncol(loadings)){
    
    # Obtain cross-loading structure
    cross_loading_structure <- loadings[start_variables[i]:end_variables[i],-i]
    
    # Ensure cross-loading structure is a matrix
    if(!is.matrix(cross_loading_structure)){
      cross_loading_structure <- matrix(
        cross_loading_structure,
        ncol = 1
      )
    }
    
    # Sample row
    random_row <- sample(
      1:nrow(cross_loading_structure),
      length(start_cross_loadings[i]:end_cross_loadings[i]),
      replace = FALSE
    )
    
    # Sample column
    random_column <- sample(
      1:ncol(cross_loading_structure),
      length(start_cross_loadings[i]:end_cross_loadings[i]),
      replace = TRUE
    )
    
    # Loop through rows and columns
    for(j in 1:length(start_cross_loadings[i]:end_cross_loadings[i])){
      
      # Populate cross-loading structure
      cross_loading_structure[
        random_row[j], random_column[j]
      ] <- cross_loading_list[[i]][j]
      
    }
    
    # Populate cross-loading structure
    loadings[start_variables[i]:end_variables[i],-i] <- cross_loading_structure

    # Check communalities
    communalities <- diag(
      loadings %*%
        factor_correlations %*%
        t(loadings)
    )
    
    # Initialize break count
    break_count <- 0
    
    # Loop through until all communalities < 0.90
    while(any(communalities >= 0.90)){
      
      # Increase break count
      break_count <- break_count + 1
      
      # Message about adjustment
      if(break_count == 1){
        
        message(
          paste(
            "Communalities for the following variable(s) were >= 0.90:",
            paste0(
              which(communalities >= 0.90),
              collapse = ", "
            ),
            "\nThe dominant loadings on these variable(s) were decreased",
            "incrementally by 0.01 until their communalities were < 0.90"
          )
        )
        
      }
      
      # Identify loadings with communalities greater than 0.90
      target_loadings <- matrix(
        loadings[which(communalities >= 0.90),],
        ncol = ncol(loadings),
        byrow = FALSE
      )
      
      # Decrease maximum loadings by 0.01
      replace_loadings <- matrix(
        apply(target_loadings, 1, function(x){
          
          # Obtain signs
          signs <- sign(x)
          
          # Compute absolute max
          x <- abs(x)
          
          # Decrease by 0.01
          x[which.max(x)] <- x[which.max(x)] - 0.01
          
          # Add back signs
          x <- x * signs
          
          # Return loadings
          return(x)
          
        }),
        ncol = ncol(loadings),
        byrow = TRUE
      )
      
      # Replace loadings
      loadings[which(communalities >= 0.90),] <- replace_loadings
      
      # Check communalities
      communalities <- diag(
        loadings %*%
          factor_correlations %*%
          t(loadings)
      )
      
    }
    
  }

  # Re-simulate data
  results <- simulate_factors(
    factors = parameters$factors,
    variables = parameters$variables,
    loadings = loadings, # use loading matrix created in this function
    cross_loadings = 0, # doesn't use
    correlations = parameters$factor_correlations,
    sample_size = nrow(lf_object$data),
    variable_categories = parameters$categories,
    categorical_limit = parameters$categorical_limit,
    skew = parameters$skew
  )
  
  # Add original results to current results
  results$original_results <- lf_object
  
  # Add class
  class(results) <- c(class(lf_object), "lf_cl")
  
  # Return results
  return(results)
  
}
