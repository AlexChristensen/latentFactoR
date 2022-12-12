#%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Simulation Helpers ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%#

# These functions are used internally to facilitate
# the use of functions that are often used in 
# simulations we run. They are not intended for users
# although they may be useful to some.

# Creates an effects table for all methods
effect_table <- function(
    formula, data, method,
    minimum_effect = c("small", "moderate", "large"),
    return_all = FALSE
)
{
  
  # Check that method exists
  if(!method %in% colnames(data)){
    stop(
      paste0(
        "'", method, "' was not found in the `data`. ",
        "Make sure the `method` argument is a column name in data"
      )
    )
  }
  
  # Identify unique methods
  unique_methods <- unique(
    data[,method]
  )
  
  # Check for only one method
  if(length(unique_methods) == 1){
    
    message(
      paste0(
        "Only one method, '", unique_methods, "', was found ",
        "in the data. Results will only include this method."
      )
    )
    
  }
  
  # Separate data into lists
  separate_data <- lapply(unique_methods, function(x){
    data[data[,method] == x,]
  })
  
  # Perform ANOVAs on separate data
  results <- lapply(
    separate_data, function(x){
      
      # Obtain effects
      obtain_effect_sizes(
        formula = formula,
        data = x,
        minimum_effect = minimum_effect,
        return_all = TRUE
      )
      
    }
  )
  
  # Stack results
  stacked_results <- t(do.call(
    cbind.data.frame, results
  ))
  
  # Rename rows
  row.names(stacked_results) <- unique_methods
  
  # Check for returning all
  if(isTRUE(return_all)){
    return(t(stacked_results))
  }else{
    
    # If minimum effect is missing, set to large
    if(missing(minimum_effect)){
      minimum_effect <- "large"
    }else{minimum_effect <- tolower(match.arg(minimum_effect))}
    
    # Check for minimum effect
    if(is.character(minimum_effect)){
      
      effect_size <- switch(
        minimum_effect,
        "small" = 0.01,
        "moderate" = 0.06,
        "large" = 0.14
      )
      
    }else if(is.numeric(minimum_effect)){
      effect_size <- minimum_effect
    }
    
    # Ensure that at least one method has
    # larger than minimum effect
    conditions_to_return <- apply(
      stacked_results,
      2, function(x){
        any(x >= effect_size)
      }
    )

    # Check for no results
    if(sum(conditions_to_return) == 0){
      
      # Return message
      message(
        paste0(
          "No effect sizes were greater than ",
          minimum_effect, ".", "\n",
          "Try setting `return_all = TRUE` to see all effect sizes."
        )
      )
      
      # Return NULL
      return(NULL)
      
    }else{
      
      # Return rounded larger etas
      return(t(stacked_results[,conditions_to_return]))
      
    }
    
  }
  
}

#' All-in-one ANOVA and partial eta squared
#' @noRd
# Updated 12.12.2022
obtain_effect_sizes <- function(
  formula, data,
  minimum_effect = c("small", "moderate", "large"),
  return_all = FALSE
)
{
  
  # If minimum effect is missing, set to large
  if(missing(minimum_effect)){
    minimum_effect <- "large"
  }else{minimum_effect <- tolower(match.arg(minimum_effect))}
  
  # Check for minimum effect
  if(is.character(minimum_effect)){
    
    effect_size <- switch(
      minimum_effect,
      "small" = 0.01,
      "moderate" = 0.06,
      "large" = 0.14
    )
    
  }else if(is.numeric(minimum_effect)){
    effect_size <- minimum_effect
  }
  
  # Obtain ANOVA formula
  formula <- as.formula(formula)
  
  # Perform ANOVA
  aov_object <- aov(
    formula = formula,
    data = data
  )
  
  # Obtain partial etas
  etas <- partial_eta(aov_object)
  
  # Check for returning all
  if(isTRUE(return_all)){
    return(round(etas, 3))
  }else{
    
    # Determine etas greater than or equal to minimum effect
    larger_etas <- etas[which(round(etas, 2) >= effect_size)]
    
    # Check for no results
    if(length(larger_etas) == 0){
      
      # Return message
      message(
        paste0(
          "No effect sizes were greater than ",
          minimum_effect, ".", "\n",
          "Try setting `return_all = TRUE` to see all effect sizes."
        )
      )
      
      # Return NULL
      return(NULL)
      
    }else{
      
      # Return rounded larger etas
      return(round(larger_etas, 3))
      
    }
    
  }

}

#' Partial eta squared
#' @noRd
# Updated 12.12.2022
partial_eta <- function(aov_object)
{
  
  # Obtain summary
  aov_summary <- summary(aov_object)
  
  # Obtain variable names
  variable_names <- trimws(row.names(aov_summary[[1]]))
  
  # Obtain sum of squares
  sum_of_squares <- aov_summary[[1]]$`Sum Sq`
  
  # Obtain residuals
  residuals <- sum_of_squares[length(sum_of_squares)]
  
  # Obtain partial eta squares
  etas <- sum_of_squares / residuals
  
  # Add names
  names(etas) <- variable_names
  
  # Remove residuals
  etas <- etas[-length(etas)]
  
  # Return partial eta squares
  return(etas)
  
}



