#%%%%%%%%%%%%%%%%%%%%%%%%%
# add_wording_effects ----
#%%%%%%%%%%%%%%%%%%%%%%%%%

#' Adds acquiescence effects to simulated data from \code{\link[latentFactoR]{simulate_factors}}
#' 
#' @param wording_data Matrix or data frame.
#' \code{\link{latentFactoR}} data that has been manipulated
#' to have wording effects with positive and negative loadings
#' 
#' @param variables Numeric (length = \code{factors}).
#' Number of variables per factor
#' 
#' @param loadings Matrix or data frame.
#' Loadings from the manipulated \code{wording_data}
#' 
#' @param categories Numeric (length = \code{variables} x \code{factors}).
#' Number of categories for each variable in the \code{wording_data}
#' 
#' @param proportion_biased_variables Numeric (length = 1 or \code{factors}).
#' Proportion of variables that should be biased with wording effects.
#' Proportion of biased variables will only count for variables below the mid-point of the \code{variable_categories}.
#' Defaults to \code{1} or all possible variables
#' 
#' @param proportion_biased_person Numeric (length = 1 or \code{proportion_biased_cases} x \code{sample_size}).
#' Person-specific parameter of how many much bias the \code{proportion_biased_cases} will
#' have over the possible biased variables. This parameter interacts with
#' \code{proportion_biased_variables}. Parameter specifies the proportion of variables 
#' that should have bias per person.
#' If one value is provided, then all biased cases will have the same proportion of variables biased.
#' Individual values are possible by providing values for each biased case
#' (\code{round(nrow(lf_object$data) * proportion_biased_cases)}). Setting individual
#' values for each biased case is not recommended
#' (use \code{proportion_biased_person_range} instead).
#' Defaults to \code{1} or all possible biased variables for all biased cases
#' 
#' @param replacement_index Numeric.
#' Indices for the cases to be replaced using the acquiescence method
#' 
#' @return Returns matrix with acquiescence wording effects added
#' 
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
#'
#' @noRd
#' 
# Applies acquiescence wording effects
# Updated 05.12.2022
add_wording_acquiescence <- function(
    wording_data, variables, loadings, categories,
    proportion_biased_variables,
    proportion_biased_person,
    replacement_index
)
{
  
  # Set sequence of variables for each factor
  end_variables <- cumsum(variables)
  start_variables <- (end_variables + 1) - variables
  
  # Initialize signs
  candidate_variables <- rep(FALSE, nrow(loadings))
  
  # Loop through each factor
  for(i in 1:ncol(loadings)){
    
    # Target variables
    target_variables <- start_variables[i]:end_variables[i]
    
    # Number of variables for re-coding
    number_recode <- round(proportion_biased_variables[i] * variables[i])
    
    # Check for whether recoding is necessary
    if(number_recode != 0){
      
      # Determine variables for re-coding
      recode_variables <- sample(
        target_variables,
        number_recode
      )
      
      # Set candidate variables for re-coding
      candidate_variables[recode_variables] <- TRUE
      
    }
    
  }
  
  # Obtain replacement data
  replacement_data <- wording_data$data
  
  # Loop through biased participants
  for(i in replacement_index){
    
    # Skip if bias is zero
    if(proportion_biased_person[i] != 0){
      
      # Target participant
      participant <- replacement_data[i,]
      
      # Participant candidate variables
      participant_candidate_variables <- candidate_variables &
        participant <= ceiling(categories / 2)
      
      # Determine person bias
      person_bias <- round(sum(participant_candidate_variables) * proportion_biased_person[i])
      
      # Sample candidate variables
      target_variables <- sample(
        which(participant_candidate_variables),
        person_bias,
        replace = FALSE
      )
      
      # Target variable categories
      target_categories <- categories[target_variables]
      
      # Target participant's variables (to nearest agreement point)
      replacement_data[i, target_variables] <- ceiling((target_categories / 2) + 1)
      
    }
    
  }
  
  # Return replacement data
  return(replacement_data[replacement_index,])
  
}

#' Adds difficulty effects to simulated data from \code{\link[latentFactoR]{simulate_factors}}
#' 
#' 
#' @param wording_data Matrix or data frame.
#' \code{\link{latentFactoR}} data that has been manipulated
#' to have wording effects with positive and negative loadings
#' 
#' @param variables Numeric (length = \code{factors}).
#' Number of variables per factor
#' 
#' @param loadings Matrix or data frame.
#' Loadings from the manipulated \code{wording_data}
#' 
#' @param categories Numeric (length = \code{variables} x \code{factors}).
#' Number of categories for each variable in the \code{wording_data}
#' 
#' @param proportion_biased_variables Numeric (length = 1 or \code{factors}).
#' Proportion of variables that should be biased with wording effects.
#' Proportion of biased variables will only count for the negative variables.
#' Defaults to \code{1} or all possible variables
#' 
#' @param proportion_biased_person Numeric (length = 1 or \code{proportion_biased_cases} x \code{sample_size}).
#' Person-specific parameter of how many much bias the \code{proportion_biased_cases} will
#' have over the possible biased variables. This parameter interacts with
#' \code{proportion_biased_variables}. Parameter specifies the proportion of variables 
#' that should have bias per person.
#' If one value is provided, then all biased cases will have the same proportion of variables biased.
#' Individual values are possible by providing values for each biased case
#' (\code{round(nrow(lf_object$data) * proportion_biased_cases)}). Setting individual
#' values for each biased case is not recommended
#' (use \code{proportion_biased_person_range} instead).
#' Defaults to \code{1} or all possible biased variables for all biased cases
#' 
#' @param replacement_index Numeric.
#' Indices for the cases to be replaced using the difficulty method
#' 
#' @return Returns matrix with difficulty wording effects added
#' 
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
#'
#' @noRd
# Applies difficulty wording effects
# Updated 05.12.2022
add_wording_difficulty <- function(
    wording_data, variables, loadings, categories,
    proportion_biased_variables,
    proportion_biased_person,
    replacement_index
)
{
  
  # Set sequence of variables for each factor
  end_variables <- cumsum(variables)
  start_variables <- (end_variables + 1) - variables
  
  # Initialize signs
  signs <- numeric(nrow(loadings))
  
  # Make all dominant loadings positive
  for(i in 1:ncol(loadings)){
    
    # Target dominant loadings
    target_loadings <- start_variables[i]:end_variables[i]
    
    # Determine sign
    signs[target_loadings] <- sign(loadings[target_loadings, i])
    
    # Check for proportion of biased variables
    if(proportion_biased_variables[i] != 1){
      
      # Set signed loadings
      signed_loadings <- signs[target_loadings]
      
      # Modify the signs (candidate variables)
      modify_signs <- round(sum(signed_loadings == -1) * proportion_biased_variables[i])
      
      # Replace some negatives with 1s
      signed_loadings[
        (modify_signs + 1):length(signed_loadings)
      ] <- 1
      
      # Return to signs vector
      signs[target_loadings] <- signed_loadings
      
    }
    
  }
  
  # Candidate variables for re-coding
  candidate_variables <- as.logical(ifelse(signs == -1, 1, 0))
  
  # Obtain replacement data
  replacement_data <- wording_data$data
  
  # Obtain person bias on candidate variables
  person_bias <- round(
    proportion_biased_person * sum(candidate_variables)
  )
  
  # Loop through biased participants
  for(i in replacement_index){
    
    # Skip if bias is zero
    if(person_bias[i] != 0){
      
      # Target participant
      participant <- replacement_data[i,]
      
      # Sample candidate variables
      target_variables <- sample(
        which(candidate_variables),
        person_bias[i],
        replace = FALSE
      )
      
      # Target variable categories
      target_categories <- categories[target_variables]
      
      # Target participant's variables
      replacement_data[i, target_variables] <- (target_categories + 1) -
        participant[target_variables]
      
    }
    
  }
  
  # Return replacement data
  return(replacement_data[replacement_index,])
  
}

#' Adds random careless effects to simulated data from \code{\link[latentFactoR]{simulate_factors}}
#' 
#' 
#' @param wording_data Matrix or data frame.
#' \code{\link{latentFactoR}} data that has been manipulated
#' to have wording effects with positive and negative loadings
#' 
#' @param variables Numeric (length = \code{factors}).
#' Number of variables per factor
#' 
#' @param loadings Matrix or data frame.
#' Loadings from the manipulated \code{wording_data}
#' 
#' @param categories Numeric (length = \code{variables} x \code{factors}).
#' Number of categories for each variable in the \code{wording_data}
#' 
#' @param proportion_biased_variables Numeric (length = 1 or \code{factors}).
#' Proportion of variables that should be biased with wording effects.
#' Defaults to \code{1} or all possible variables
#' 
#' @param proportion_biased_person Numeric (length = 1 or \code{proportion_biased_cases} x \code{sample_size}).
#' Person-specific parameter of how many much bias the \code{proportion_biased_cases} will
#' have over the possible biased variables. This parameter interacts with
#' \code{proportion_biased_variables}. Parameter specifies the proportion of variables 
#' that should have bias per person.
#' If one value is provided, then all biased cases will have the same proportion of variables biased.
#' Individual values are possible by providing values for each biased case
#' (\code{round(nrow(lf_object$data) * proportion_biased_cases)}). Setting individual
#' values for each biased case is not recommended
#' (use \code{proportion_biased_person_range} instead).
#' Defaults to \code{1} or all possible biased variables for all biased cases
#' 
#' @param replacement_index Numeric.
#' Indices for the cases to be replaced using the random careless method
#' 
#' @return Returns matrix with random careless wording effects added
#' 
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
#'
#' @noRd
# Applies random careless wording effects
# Updated 05.12.2022
add_wording_random_careless <- function(
    wording_data, variables, loadings, categories,
    proportion_biased_variables,
    proportion_biased_person,
    replacement_index
)
{
  
  # Set sequence of variables for each factor
  end_variables <- cumsum(variables)
  start_variables <- (end_variables + 1) - variables
  
  # Initialize signs
  candidate_variables <- rep(FALSE, nrow(loadings))
  
  # Loop through each factor
  for(i in 1:ncol(loadings)){
    
    # Target variables
    target_variables <- start_variables[i]:end_variables[i]
    
    # Number of variables for re-coding
    number_recode <- round(proportion_biased_variables[i] * variables[i])
    
    # Check for whether recoding is necessary
    if(number_recode != 0){
      
      # Determine variables for re-coding
      recode_variables <- sample(
        target_variables,
        number_recode
      )
      
      # Set candidate variables for re-coding
      candidate_variables[recode_variables] <- TRUE
      
    }
    
  }
  
  # Obtain replacement data
  replacement_data <- wording_data$data
  
  # Loop through biased participants
  for(i in replacement_index){
    
    # Skip if bias is zero
    if(proportion_biased_person[i] != 0){
      
      # Target participant
      participant <- replacement_data[i,]
      
      # Determine person bias
      person_bias <- round(sum(candidate_variables) * proportion_biased_person[i])
      
      # Sample candidate variables
      target_variables <- sample(
        which(candidate_variables),
        person_bias,
        replace = FALSE
      )
      
      # Loop through variables
      for(recode in target_variables){
        
        # Insert random values
        replacement_data[i, recode] <- sample(
          1:categories[recode], 1
        )
        
      }
      
    }
    
  }
  
  # Return replacement data
  return(replacement_data[replacement_index,])
  
}

#' Adds straight line effects to simulated data from \code{\link[latentFactoR]{simulate_factors}}
#' 
#' 
#' @param wording_data Matrix or data frame.
#' \code{\link{latentFactoR}} data that has been manipulated
#' to have wording effects with positive and negative loadings
#' 
#' @param variables Numeric (length = \code{factors}).
#' Number of variables per factor
#' 
#' @param loadings Matrix or data frame.
#' Loadings from the manipulated \code{wording_data}
#' 
#' @param categories Numeric (length = \code{variables} x \code{factors}).
#' Number of categories for each variable in the \code{wording_data}
#' 
#' @param proportion_biased_variables Numeric (length = 1 or \code{factors}).
#' Proportion of variables that should be biased with wording effects.
#' Defaults to \code{1} or all possible variables
#' 
#' @param proportion_biased_person Numeric (length = 1 or \code{proportion_biased_cases} x \code{sample_size}).
#' Person-specific parameter of how many much bias the \code{proportion_biased_cases} will
#' have over the possible biased variables. This parameter interacts with
#' \code{proportion_biased_variables}. Parameter specifies the proportion of variables 
#' that should have bias per person.
#' If one value is provided, then all biased cases will have the same proportion of variables biased.
#' Individual values are possible by providing values for each biased case
#' (\code{round(nrow(lf_object$data) * proportion_biased_cases)}). Setting individual
#' values for each biased case is not recommended
#' (use \code{proportion_biased_person_range} instead).
#' Defaults to \code{1} or all possible biased variables for all biased cases
#' 
#' @param replacement_index Numeric.
#' Indices for the cases to be replaced using the straight line method
#' 
#' @return Returns matrix with straight line wording effects added
#' 
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
#'
#' @noRd
# Applies straight line wording effects
# Updated 05.12.2022
add_wording_straight_line <- function(
    wording_data, variables, loadings, categories,
    proportion_biased_variables,
    proportion_biased_person,
    replacement_index
)
{
  
  # Set sequence of variables for each factor
  end_variables <- cumsum(variables)
  start_variables <- (end_variables + 1) - variables
  
  # Initialize signs
  candidate_variables <- rep(FALSE, nrow(loadings))
  
  # Loop through each factor
  for(i in 1:ncol(loadings)){
    
    # Target variables
    target_variables <- start_variables[i]:end_variables[i]
    
    # Number of variables for re-coding
    number_recode <- round(proportion_biased_variables[i] * variables[i])
    
    # Check for whether recoding is necessary
    if(number_recode != 0){
      
      # Determine variables for re-coding
      recode_variables <- sample(
        target_variables,
        number_recode
      )
      
      # Set candidate variables for re-coding
      candidate_variables[recode_variables] <- TRUE
      
    }
    
  }
  
  # Obtain replacement data
  replacement_data <- wording_data$data
  
  # Loop through biased participants
  for(i in replacement_index){
    
    # Skip if bias is zero
    if(proportion_biased_person[i] != 0){
      
      # Target participant
      participant <- replacement_data[i,]
      
      # Determine person bias
      person_bias <- round(sum(candidate_variables) * proportion_biased_person[i])
      
      # Sample candidate variables
      target_variables <- sample(
        which(candidate_variables),
        person_bias,
        replace = FALSE
      )
      
      # Obtain max category
      maximum_category <- max(categories)
      
      # Draw random number from maximum category
      random_category <- sample(
        1:maximum_category, 1
      )
      
      # Proportion of category
      proportion_category <- random_category / maximum_category
      
      # Determine replacement category value
      replacement_category <- round(categories * proportion_category)
      
      # Can't have zeros
      replacement_category <- ifelse(
        replacement_category == 0, 1, replacement_category
      )
      
      # Replace target variables
      replacement_data[i, target_variables] <- replacement_category[target_variables]
      
    }
    
  }
  
  # Return replacement data
  return(replacement_data[replacement_index,])
  
}