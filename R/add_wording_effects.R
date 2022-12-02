#' Adds Wording Effects to \code{\link[latentFactoR]{simulate_factors}} Data
#'
#' Adds wording effects to simulated data from \code{\link[latentFactoR]{simulate_factors}}. 
#' See examples to get started
#' 
#' @param lf_object Data object from \code{\link[latentFactoR]{simulate_factors}}.
#' Data \strong{must} be categorical. If data are not categorical, then
#' there function with throw an error
#' 
#' @param method Character (length = 1).
#' Method to generate wording effect to add to the data.
#' Description of methods:
#' 
#' \itemize{
#' 
#' \item{\code{"acquiescence"}}
#' {Generates new data with flipped dominant loadings 
#' (based on \code{proportion_negative}) and sets
#' thresholds below the mid-point to \code{-Inf} for
#' all variables creating a restricted range of responding
#' (e.g., only 3s, 4s, and 5s on a 5-point Likert scale)}
#' 
#' \item{\code{"difficulty"}}
#' {Generates new data with flipped dominant loadings 
#' (based on \code{proportion_negative}) and uses this data
#' as the data without wording effects. Then, the signs of the
#' dominant loadings are obtained and the dominant loadings are 
#' made to be absolute. Finally, the skews are multiplied by
#' the signs of the original dominant loadings when generating
#' the data with the wording effects}
#' 
#' \item{\code{"random_careless"}}
#' {Number of cases up to \code{proportion_biased_cases} are sampled
#' and replaced by values from a random uniform distribution ranging
#' between the lowest and highest response category for each variable.
#' These values then replace the values in the original data} 
#' 
#' \item{\code{"straight_line"}}
#' {Coming soon...} 
#' 
#' }
#' 
#' @param proportion_negative Numeric (length = 1 or \code{factors}).
#' Proportion of variables that should have negative (or flipped) dominant loadings across all
#' or each factor. Accepts number of variables as well.
#' The first variables on each factor, up to the corresponding proportion, will be
#' flipped. Set to \code{0} to not have any loadings flipped.
#' Defaults to \code{0.50}
#' 
#' @param proportion_negative_range Numeric (length = 2).
#' Range of proportion of variables that are randomly selected from
#' a uniform distribution. Accepts number of number of variables as well.
#' Defaults to \code{NULL}
#' 
#' @param proportion_biased_cases Numeric (length = 1).
#' Proportion of cases that should be biased with wording effects.
#' Also accepts number of cases to be biased. The first \emph{n} number of cases,
#' up to the corresponding proportion, will be biased.
#' Defaults to \code{0.10} or 10 percent of cases.
#' 
#' @param proportion_biased_variables Numeric (length = 1 or \code{factors}).
#' Proportion of variables that should be biased with wording effects.
#' For \code{method = "difficulty"}, proportion of biased variables will only
#' count for the negative variables.
#' For \code{method = "acquiescence"}, proportion of biased variables will only
#' count for variables below the mid-point of the \code{variable_categories}.
#' Defaults to \code{1} or all possible variables
#' 
#' @param proportion_biased_variables_range Numeric (length = 2).
#' Range of proportion of variables that should be biased with wording effects.
#' Values are drawn randomly from a uniform distribution.
#' Defaults to \code{NULL}
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
#' @param proportion_biased_person_range Numeric (length = 2).
#' Range to randomly draw bias from a uniform distribution. Allows for random
#' person-specific bias to be obtained.
#' Defaults to \code{NULL}
#' 
#' @return Returns a list containing:
#' 
#' \item{data}{Biased data simulated data from the specified factor model}
#' 
#' \item{unbiased_data}{The corresponding unbiased data prior to replacing values
#' to generate the (biased) \code{data}}
#' 
#' \item{biased_sample_size}{The number of cases that have biased data}
#' 
#' \item{adjusted_results}{Bias-adjusted \code{lf_object} input into function} 
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
#'   sample_size = 1000, # number of cases = 1000
#'   variable_categories = 5 # 5-point Likert scale
#' )
#' 
#' # Add wording effects using acquiescence method
#' two_factor_acquiescence <- add_wording_effects(
#'   lf_object = two_factor,
#'   proportion_negative = 0.50,
#'   proportion_biased_cases = 0.10,
#'   method = "acquiescence"
#' )
#' 
#' # Add wording effects using difficulty method
#' two_factor_difficulty <- add_wording_effects(
#'   lf_object = two_factor,
#'   proportion_negative = 0.50,
#'   proportion_biased_cases = 0.10,
#'   method = "difficulty"
#' )
#' 
#' # Add wording effects using random careless method
#' two_factor_random_careless <- add_wording_effects(
#'   lf_object = two_factor,
#'   proportion_negative = 0.50,
#'   proportion_biased_cases = 0.10,
#'   method = "random_careless"
#' )
#' 
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
#' 
#' @references
#' Garcia-Pardina, A., Abad, F. J., Christensen, A. P., Golino, H., & Garrido, L. E. (2022).
#' Dimensionality assessment in the presence of wording effects: A network psychometric and factorial approach.
#' \emph{PsyArXiv}.
#' 
#' Garrido, L. E., Golino, H., Christensen, A. P., Martinez-Molina, A., Arias, V. B., Guerra-Pena, K., ... & Abad, F. J. (2022).
#' A systematic evaluation of wording effects modeling under the exploratory structural equation modeling framework.
#' \emph{PsyArXiv}.
#' 
#' @importFrom utils data
#'
#' @export
#'
# Add wording effects to simulated data
# Updated 01.12.2022
add_wording_effects <- function(
    lf_object,
    method = c(
      "acquiescence", "difficulty",
      "random_careless", "straight_line"
    ),
    proportion_negative = 0.50,
    proportion_negative_range = NULL,
    proportion_biased_cases = 0.10,
    proportion_biased_variables = 1,
    proportion_biased_variables_range = NULL,
    proportion_biased_person = 1,
    proportion_biased_person_range = NULL
)
{
  
  # Match `method` argument (no default)
  if(missing(method)){
    stop("The `method` argument must be set.")
  }else{method <- match.arg(method)}
  
  # Ensure `method` is lowercase
  method <- tolower(method)
  
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
  
  # Ensure data is categorical
  if(any(lf_object$parameters$categories > 6)){
    
    # Produce error
    stop(
      paste(
        "Data input into `lf_object` must all be categorical (6 categories or less).",
        "These variables were found to be continuous:",
        paste(which(lf_object$parameters$categories > 6), collapse = ", ")
      )
    )
    
  }
  
  # Obtain number of cases
  sample_size <- nrow(lf_object$data)
  
  # Ensure appropriate methods
  type_error(proportion_biased_cases, "numeric");
  
  # Ensure appropriate lengths
  length_error(proportion_biased_cases, 1);
  
  # Convert biased cases to proportions
  if(proportion_biased_cases > 1){
    proportion_biased_cases <- proportion_biased_cases / sample_size
  }
  
  # Ensure appropriate ranges
  range_error(proportion_biased_cases, c(0, 1));
  
  # Obtain sample size
  biased_sample_size <- round(
    proportion_biased_cases * sample_size
  )
  
  # Obtain parameters from simulated data
  parameters <- lf_object$parameters
  
  # Check for percentage negative range
  if(!is.null(proportion_negative_range)){
    type_error(proportion_negative_range, "numeric") # object type error
    length_error(proportion_negative_range, 2) # object length error
    
    # Check for number of variables in range
    if(any(proportion_negative_range > 1)){
      
      # Target values
      target_negative <- which(proportion_negative_range > 1)
      
      # Ensure proportions
      proportion_negative_range[target_negative] <-
        proportion_negative_range[target_negative] / parameters$variables[target_negative]
      
    }
    
    # Check for error in range
    range_error(proportion_negative_range, c(0, 1)) # object range error
    proportion_negative <- runif(
      parameters$factors,
      min = min(proportion_negative_range),
      max = max(proportion_negative_range)
    )
  }
  
  # Ensure appropriate types
  type_error(proportion_negative, "numeric");
  
  # Ensure appropriate lengths
  length_error(proportion_negative, c(1, parameters$factors));
  
  # Set proportions
  if(length(proportion_negative) == 1){
    proportion_negative <- rep(proportion_negative, parameters$factors)
  }
  
  # Convert negative wording proportions to proportions
  if(any(proportion_negative > 1)){
    
    # Target values
    target_negative <- which(proportion_negative > 1)
    
    # Ensure proportions
    proportion_negative[target_negative] <-
      proportion_negative[target_negative] / parameters$variables[target_negative]
    
  }
  
  # Ensure appropriate ranges
  range_error(proportion_negative, c(0, 1));
  
  # Check for proportion variable range
  if(!is.null(proportion_biased_variables_range)){
    type_error(proportion_biased_variables_range, "numeric") # object type error
    length_error(proportion_biased_variables_range, 2) # object length error
    
    # Check for number of variables in range
    if(any(proportion_biased_variables_range > 1)){
      
      # Target values
      target_variables <- which(proportion_biased_variables_range > 1)
      
      # Ensure proportions
      proportion_biased_variables_range[target_variables] <-
        proportion_biased_variables_range[target_variables] / parameters$variables[target_variables]
      
    }
    
    # Check for error in range
    range_error(proportion_biased_variables_range, c(0, 1)) # object range error
    proportion_biased_variables <- runif(
      parameters$factors,
      min = min(proportion_biased_variables_range),
      max = max(proportion_biased_variables_range)
    )
    
  }
  
  # Ensure appropriate types
  type_error(proportion_biased_variables, "numeric");
  
  # Ensure appropriate length
  length_error(proportion_biased_variables, c(1, parameters$factors))
  
  # Set proportions
  if(length(proportion_biased_variables) == 1){
    proportion_biased_variables <- rep(proportion_biased_variables, parameters$factors)
  }
  
  # Convert negative wording proportions to proportions
  if(any(proportion_biased_variables > 1)){
    
    # Target values
    target_variables <- which(proportion_biased_variables > 1)
    
    # Ensure proportions
    proportion_biased_variables[target_variables] <-
      proportion_biased_variables[target_variables] / parameters$variables[target_variables]
    
  }
  
  # Ensure appropriate ranges
  range_error(proportion_biased_variables, c(0, 1));
  
  # Check for proportion variable range
  if(!is.null(proportion_biased_person_range)){
    type_error(proportion_biased_person_range, "numeric") # object type error
    length_error(proportion_biased_person_range, 2) # object length error
    range_error(proportion_biased_person_range, c(0, 1)) # object range error
    proportion_biased_person <- runif(
      biased_sample_size,
      min = min(proportion_biased_person_range),
      max = max(proportion_biased_person_range)
    )
  }
  
  # Ensure appropriate types
  type_error(proportion_biased_person, "numeric");
  
  # Ensure appropriate length
  length_error(proportion_biased_person, c(1, biased_sample_size))
  
  # Set proportions
  if(length(proportion_biased_person) == 1){
    proportion_biased_person <- rep(proportion_biased_person, biased_sample_size)
  }
  
  # Ensure appropriate ranges
  range_error(proportion_biased_person, c(0, 1));
  
  # Obtain loadings
  loadings <- parameters$loadings
  
  # Obtain variables
  variables <- parameters$variables
  
  # Set sequence of variables for each factor
  end_variables <- cumsum(parameters$variables)
  start_variables <- (end_variables + 1) - parameters$variables
  
  # Flip dominant loadings
  for(i in 1:ncol(loadings)){
    
    # Obtain number of flipped variables
    negative_variables <- round(proportion_negative[i] * variables[i])
    
    # Check for zero negative variables
    if(negative_variables != 0){
      
      # Target dominant loadings
      target_loadings <- start_variables[i]:end_variables[i]
      
      # Make loadings absolute
      loadings[target_loadings, i] <- abs(loadings[target_loadings, i])
      
      # Set dominant loadings to inverse
      loadings[
        target_loadings[1:negative_variables],
        i
      ] <- -loadings[
        target_loadings[1:negative_variables],
        i
      ]
      
    }
    
  }
  
  # Ensure proper signs for skew
  for(i in 1:ncol(loadings)){
    
    # Target dominant loadings
    target_loadings <- start_variables[i]:end_variables[i]
    
    # Determine sign
    target_signs <- sign(loadings[target_loadings, i])
    
    # Target skews
    skews <- parameters$skew[target_loadings]
    
    # Check if all skew signs are in the same direction
    if(all(sign(skews) == -1)){ # All in negative direction
      
      # Make all skew for negative variables positive
      skews[target_signs == -1] <- abs(skews[target_signs == -1])
      
    }else if(all(sign(skews) == 1)){ # All in negative direction
      
      # Make all skew for positive variables negative
      skews[target_signs == 1] <- -abs(skews[target_signs == 1])
      
    }else{ 
      
      # If mixed, base signs on mode of positive variables
      # with preference for negative skew (i.e., ties go to negative skew)
      
      # Obtain skew signs for positive variables
      positive_skew_signs <- sign(skews[target_signs == 1])
      
      # Determine whether mode with ties going to negative skew
      if(
        sum(positive_skew_signs == -1) >=
        sum(positive_skew_signs == 1)
      ){
        
        # Make all skew for positive variables negative
        skews[target_signs == 1] <- -abs(skews[target_signs == 1])
        
        # Make all skew for negative variables positive
        skews[target_signs == -1] <- abs(skews[target_signs == -1])
        
      }else{ # Do positive skew for positive variables
        
        # Make all skew for positive variables positive
        skews[target_signs == 1] <- abs(skews[target_signs == 1])
        
        # Make all skew for negative variables negative
        skews[target_signs == -1] <- -abs(skews[target_signs == -1])
        
      }
      
    }
    
    # Re-assign skew to target variables
    parameters$skew[target_loadings] <- skews
    
  }
  
  # Re-generate data
  wording_data <- simulate_factors(
    factors = parameters$factors,
    variables = parameters$variables,
    loadings = loadings,
    cross_loadings = loadings,
    correlations = parameters$factor_correlations,
    sample_size = nrow(lf_object$data),
    variable_categories = parameters$categories,
    categorical_limit = parameters$categorical_limit,
    skew = parameters$skew
  )
  
  # Update parameters
  parameters <- wording_data$parameters
  
  # Obtain loadings
  loadings <- parameters$loadings
  
  # Obtain variables
  variables <- parameters$variables
  
  # Set sequence of variables for each factor
  end_variables <- cumsum(parameters$variables)
  start_variables <- (end_variables + 1) - parameters$variables
  
  # Check for difficulty
  if(method == "difficulty"){
    
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
    replacement_data <- wording_data$data[1:biased_sample_size,]
    
    # Obtain person bias on candidate variables
    person_bias <- round(
      proportion_biased_person * sum(candidate_variables)
    )
    
    # Loop through biased participants
    for(i in 1:nrow(replacement_data)){
      
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
        target_categories <- parameters$categories[target_variables]
        
        # Target participant's variables
        replacement_data[i, target_variables] <- (target_categories + 1) -
          participant[target_variables]
        
      }
      
    }
    
  }else if(method == "acquiescence"){
    
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
    replacement_data <- wording_data$data[1:biased_sample_size,]
    
    # Loop through biased participants
    for(i in 1:nrow(replacement_data)){
      
      # Skip if bias is zero
      if(proportion_biased_person[i] != 0){
        
        # Target participant
        participant <- replacement_data[i,]
        
        # Participant candidate variables
        participant_candidate_variables <- candidate_variables &
          participant <= (parameters$categories / 2)
        
        # Determine person bias
        person_bias <- round(sum(participant_candidate_variables) * proportion_biased_person[i])
        
        # Sample candidate variables
        target_variables <- sample(
          which(participant_candidate_variables),
          person_bias,
          replace = FALSE
        )
        
        # Target variable categories
        target_categories <- parameters$categories[target_variables]
        
        # Target participant's variables (to mid-point or nearest agreement point)
        replacement_data[i, target_variables] <- round((target_categories + 1) / 2)

      }
      
    }
    
  }else if(method == "random_careless"){
    
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
    replacement_data <- wording_data$data[1:biased_sample_size,]
    
    # Loop through biased participants
    for(i in 1:nrow(replacement_data)){
      
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
            1:parameters$categories[recode], 1
          )
          
        }
        
      }
      
    }
    
  }else if(method == "straight_line"){
    
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
    replacement_data <- wording_data$data[1:biased_sample_size,]
    
    # Loop through biased participants
    for(i in 1:nrow(replacement_data)){
      
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
        maximum_category <- max(parameters$categories)
        
        # Draw random number from maximum category
        random_category <- sample(
          1:maximum_category, 1
        )
        
        # Proportion of category
        proportion_category <- random_category / maximum_category
        
        # Determine replacement category value
        replacement_category <- round(parameters$categories * proportion_category)
        
        # Can't have zeros
        replacement_category <- ifelse(
          replacement_category == 0, 1, replacement_category
        )
        
        # Replace target variables
        replacement_data[i, target_variables] <- replacement_category[target_variables]
        
      }
      
    }
    
  }
  
  # Replace original data with replacement data
  new_data <- wording_data$data
  new_data[1:biased_sample_size,] <- replacement_data
  
  # Populate results
  results <- list(
    data = new_data,
    unbiased_data = wording_data$data,
    biased_sample_size = biased_sample_size,
    adjusted_results = wording_data,
    original_results = lf_object
  )
  
  # Add class
  class(results) <- c(class(lf_object), "lf_we")
  
  # Return results
  return(results)
  
}

# NOTES FOR FUTURE IMPLEMENTATION

# Person parameter
# different proportion for bias
# 100% = all replacement
# 50% = 50% of the items

# Item parameter
# different proportion for bias
# all equal = then all the same for person
# not all the same: sampling weight

# sample(
#   1:ncol(original_data), # items
#   ncol(original_data) * 0.50, # person parameter
#   prob = c(  # item parameter
#     rep(0, 3), rep(0.30, 3), rep(0.70, 3)
#   )
# )