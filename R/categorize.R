#' Categorize Continuous Data
#'
#' Categorizes continuous data based on Garrido, Abad and Ponsoda (2011; see references).
#' Categorical data with 2 to 6 categories can include skew between -2 to 2 in
#' increments of 0.05
#'
#' @param data Numeric (length = n).
#' A vector of continuous data with \emph{n} values.
#' For matrices, use \code{apply}
#' 
#' @param categories Numeric (length = 1).
#' Number of categories to create.
#' Between 2 and 6 categories can be used with skew
#' 
#' @param skew_value Numeric (length = 1).
#' Value of skew.
#' Ranges between -2 to 2 in increments of 0.05.
#' Skews not in this sequence will be converted to
#' the nearest value in this sequence.
#' Defaults to \code{0} or no skew
#' 
#' @return Returns a numeric vector of the categorize data
#'
#' @examples
#' # Dichotomous data (no skew)
#' dichotomous <- categorize(
#'   data = rnorm(1000),
#'   categories = 2
#' )
#' 
#' # Dichotomous data (with positive skew)
#' dichotomous_skew <- categorize(
#'   data = rnorm(1000),
#'   categories = 2,
#'   skew_value = 1.25
#' )
#' 
#' # 5-point Likert scale (no skew)
#' five_likert <- categorize(
#'   data = rnorm(1000),
#'   categories = 5 
#' )
#' 
#' # 5-point Likert scale (negative skew)
#' five_likert <- categorize(
#'   data = rnorm(1000),
#'   categories = 5,
#'   skew_value = -0.45
#' )
#' 
#' @author
#' Maria Dolores Nieto Canaveras <mnietoca@nebrija.es>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>,
#' Hudson Golino <hfg9s@virginia.edu>,
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>
#' 
#' @references
#' Garrido, L. E., Abad, F. J., & Ponsoda, V. (2011). \cr
#' Performance of Velicer’s minimum average partial factor retention method with categorical variables. \cr
#' \emph{Educational and Psychological Measurement}, \emph{71}(3), 551-570.
#' 
#' Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., ... & Martinez-Molina, A. (2020).
#' Investigating the performance of exploratory graph analysis and traditional techniques to identify the number of latent factors: A simulation and tutorial.
#' \emph{Psychological Methods}, \emph{25}(3), 292-320.
#'
#' @export
#'
# Categorization function
# Updated 30.11.2022
categorize <- function(data, categories, skew_value = 0)
{
  
  # Possible skew values
  possible_skews <- seq(-2, 2, 0.05)
  
  # Check if skew is in possible values
  if(!skew_value %in% possible_skews){
    
    # Round to hundredths digit
    skew <- round(skew_value, 2)
    
    # Get differences
    skew_difference <- abs(skew - possible_skews)
    
    # Get skew
    skew_value <- possible_skews[which.min(skew_difference)]
    
  }
  
  # Load skew tables
  skew_tables <- get(data(
    "skew_tables",
    package = "latentFactoR",
    envir = environment()
  ))
  
  # Switch categories to character
  categories <- switch(
    as.character(categories),
    "2" = "two",
    "3" = "three",
    "4" = "four",
    "5" = "five",
    "6" = "six"
  )
  
  # Obtain skew table
  skew_table <- skew_tables[[categories]]
  
  # Obtain skew values
  skew_values <- skew_table[,formatC(
    skew_value, digits = 2,
    format = "f", flag = "0"
  )]
  
  # Add skew to data (see `utils-latentFactoR` for function)
  skewed_data <- skew_single_variable(
    data = data, skew_values = skew_values
  )
  
  # Return categorized data
  return(skewed_data)
  
}  
  
  
  