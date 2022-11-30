
# Add wording effect
add_wording_effects <- function(
    lf_object,
    proportion_negative = 0.00,
    proportion_negative_range = NULL,
    proportion_biased_cases = 0.20,
    type = c(
      "acquiescence", "difficulty",
      "random_careless", "straight_line"
    )
)
{
  
  # Not yet implemented
  if(type == "straight_line"){
    
    # Produce error
    stop(
      "The \"straight_line\" argument for `type` is not yet implemented."
    )
    
  }
  
  # Match `type` argument (no default)
  if(missing(type)){
    stop("The `type` argument must be set.")
  }else{type <- match.arg(type)}
  
  # Ensure `type` is lowercase
  type <- tolower(type)
  
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
  
  # Ensure appropriate types
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
  
  # Random careless responding
  if(type == "random_careless"){
    
    # Initialize biased data
    biased_data <- lf_object$data
    
    # Loop through variables
    for(i in 1:ncol(biased_data)){
      
      # Insert random values
      biased_data[1:biased_sample_size,i] <- as.numeric(
        cut( # cuts evenly distributed categories
          runif(biased_sample_size), # generates random numbers
          breaks = parameters$categories[i]
        )
      )

    }
    
    # Populate results
    results <- list(
      biased_data = new_data,
      unbiased_data = lf_object$data,
      biased_sample_size = biased_sample_size,
      original_results = lf_object
    )
    
    # Add class
    class(results) <- c(class(lf_object), "lf_we")
    
    # Return data
    return(results)
  
  }
  
  
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
  
  # Check for difficulty
  if(type == "difficulty"){
    
    # Update parameters
    parameters <- wording_data$parameters
    
    # Obtain loadings
    loadings <- parameters$loadings
    
    # Obtain variables
    variables <- parameters$variables
    
    # Set sequence of variables for each factor
    end_variables <- cumsum(parameters$variables)
    start_variables <- (end_variables + 1) - parameters$variables
    
    # Initialize signs
    signs <- numeric(nrow(loadings))
    
    # Make all dominant loadings positive
    for(i in 1:ncol(loadings)){
      
      # Target dominant loadings
      target_loadings <- start_variables[i]:end_variables[i]
      
      # Determine sign
      signs[target_loadings] <- sign(loadings[target_loadings, i])
      
      # Make loadings absolute
      loadings[target_loadings, i] <- abs(loadings[target_loadings, i])
      
    }

    # Just using `parameters$skew * signs`
    # in the argument "skew" below does not work
    # However, inputting the output into `parameters$skew` does work
    # I don't know, R is weird sometimes
    if(all(parameters$skew != 0)){
      parameters$skew <- parameters$skew * signs
    }
    
    # Re-generate data
    difficulty_data <- simulate_factors(
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
    
  }
  
  # Re-set parameters
  parameters <- wording_data$parameters
  
  # Obtain skews
  skews <- parameters$skew
  
  # Obtain categories
  categories <- parameters$categories
  
  # Load skew tables
  skew_tables <- get(data(
    "skew_tables",
    package = "latentFactoR",
    envir = environment()
  ))
  
  # Initialize threshold list
  threshold_list <- vector("list")
  
  # Loop through each variable
  for(i in 1:length(skews)){
    
    # Obtain category name
    category_name <- switch(
      as.character(categories[i]),
      "2" = "two",
      "3" = "three",
      "4" = "four",
      "5" = "five",
      "6" = "six"
    )
    
    # Obtain skew thresholds
    skew_thresholds <- skew_tables[[category_name]]
    
    # Format skew name
    skew_name <- formatC(
      x = skews[i], digits = 2,
      format = "f", flag = "0"
    )
    
    # Obtain thresholds
    threshold_list[[i]] <- unname(skew_thresholds[,skew_name])
    
  }
  
  # Acquiescence
  if(type == "acquiescence"){
    
    # Loop through threshold list
    update_thresholds <- lapply(threshold_list, function(x){
      
      # Determine mid-points
      mid_point <- switch(
        as.character(length(x)), # number of thresholds
        "1" = 1, # 2 categories
        "2" = 1, # 3 categories
        "3" = 2, # 4 categories
        "4" = 2, # 5 categories
        "5" = 3 # 6 categories
      )
      
      # Change up to mid-point
      x[1:mid_point] <- -Inf
      
      # Return updated threshold
      return(x)
      
    })
    
    # Obtain biased data
    biased_data <- lf_object$continuous_data[1:biased_sample_size,]
    
    # Apply skew
    skew_biased_data <- lapply(seq_along(update_thresholds), function(i){
      
      # Obtain skew values
      skew_values <- update_thresholds[[i]]
      
      # Apply skew to data
      skew_biased_data <- skew_single_variable(
        data = biased_data[,i],
        skew_values = skew_values
      )
      
      # Return data
      return(skew_biased_data)
      
    })
    
    # Combine skewed data
    replacement_data <- t(as.matrix(
      unname(do.call(rbind.data.frame, skew_biased_data))
    ))
    
  }else if(type == "difficulty"){
    
    # Obtain replacement data
    replacement_data <- difficulty_data$data[1:biased_sample_size,]
    
  }
  
  
  # Replace original data with replacement data
  new_data <- wording_data$data
  new_data[1:biased_sample_size,] <- replacement_data
  
  # Populate results
  results <- list(
    biased_data = new_data,
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
