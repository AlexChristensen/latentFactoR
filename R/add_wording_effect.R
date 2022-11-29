# Load {latentFactoR}
library(latentFactoR)

# Simulate regular data
simulated <- simulate_factors(
  factors = 3,
  variables = 3,
  loadings = 0.50,
  cross_loadings = 0,
  correlations = 0.30,
  sample_size = 500,
  variable_categories = c(2, 2, 3, 3, 4, 4, 5, 5, 5),
  skew_range = c(-1, 1)
)

# Ensure data are categorical
# ERROR if data are continuous
# or ALLOW number of categories to be set


# Positive and negative worded items
lf_object <- simulated

# Number of percent negatively worded
# Default zero
proportion_negative <- 0.333 # 33%

# First X number of variables

# Obtain parameters from simulated data
parameters <- lf_object$parameters

# Set range null
proportion_negative_range <- c(0, 1)

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

# Proportion of cases that will be biased
proportion_biased_cases <- 0.20 # 20%

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

# Change N cases

# bias to all; bias to some
#
# acquiescence
#
# factor model is same as original model
#
# thresholds are changed (manipulate for each variable)
#
#
#
# parameters$skew # good model

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










































