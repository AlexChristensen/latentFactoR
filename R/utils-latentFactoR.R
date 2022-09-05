#%%%%%%%%%%%%%%%%%%%%%%%%%%
# GENERATION FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%

# Based on 
# Garrido, L. E., Abad, F. J., & Ponsoda, V. (2011).
# Performance of Velicer’s minimum average partial factor retention
# method with categorical variables.
# Educational and Psychological Measurement, 71(3), 551-570.
# https://doi.org/10.1177/0013164410389489
#
#' @noRd
# Generates skewed data for continuous data
# Updated 09.08.2022
skew_continuous <- function(
    skewness,
    data = NULL,
    sample_size = 1000000,
    tolerance = 0.00001
)
{
  
  # Original skewness
  original_skewness <- skewness
  
  # Generate data
  if(is.null(data)){
    data <- rnorm(sample_size)
  }
  
  # Kurtosis
  kurtosis <- 1
  
  # Skew data
  skew_data <- sinh(
    kurtosis * (asinh(data) + skewness) 
  )
  
  # Observed skew in data
  observed_skew <- psych::skew(skew_data)
  
  # Minimize difference
  while(abs(observed_skew - original_skewness) > tolerance){
    
    # Obtain difference
    difference <- observed_skew - original_skewness
    
    # Decrease skewness
    # Negative values increase
    # Positive values decrease
    skewness <- skewness - difference
    
    # Skew data
    skew_data <- sinh(
      kurtosis * (asinh(data) + skewness) 
    )
    
    # Observed skew in data
    observed_skew <- psych::skew(skew_data)
    
  }
  
  # Return skewed data
  return(skew_data)
  
}

# Based on 
# Garrido, L. E., Abad, F. J., & Ponsoda, V. (2011).
# Performance of Velicer’s minimum average partial factor retention
# method with categorical variables.
# Educational and Psychological Measurement, 71(3), 551-570.
# https://doi.org/10.1177/0013164410389489
#
#' @noRd
# Generates skew
# Updated 09.08.2022
skew_generator <- function(
    skewness, categories,
    reduction_factor = 0.75,
    sample_size = 1000000,
    initial_proportion = 0.50,
    tolerance = 0.00001
)
{
  
  # Initialize skew matrix
  skew_matrix <- matrix(
    0, nrow = categories, ncol = categories
  )
  
  # Initialize cases
  cases <- numeric(sample_size)
  
  # Loop through categories
  for(i in 2:categories){
    
    # Initialize (largest category) proportion
    proportion <- initial_proportion
    
    # Current proportion for rest of categories
    remaining_proportion <- 1 - proportion
    
    # Initialize categories allocations
    allocation <- 1
    allocation_1 <- 1
    
    # Loop through with reduction factor
    if(i > 2){
      for(j in 1:(i-2)){
        allocation_1 <- allocation_1 * reduction_factor
        allocation <- allocation + allocation_1
      }
    }
    
    # Divide remaining proportion by allocations
    divided_proportion <- remaining_proportion / allocation
    
    # Undefined objects
    propinf <- 1 / sample_size
    propsup <- initial_proportion
    E <- divided_proportion / reduction_factor
    
    # Loop through
    for(j in 1:i){
      cases[round(sample_size * propinf):round(sample_size * propsup)] <- j
      E <- E * reduction_factor
      propinf <- propsup
      propsup <- propinf + E
    }
    
    # Compute skewness
    skew_actual <- psych::skew(cases)
    
    # Limits
    limitsup <- 1
    limitsinf <- 0
    
    # Ensure skew within tolerance
    while(abs(skew_actual - skewness) > tolerance){
      
      # Skew greater than
      if(skew_actual < skewness){
        limitinf <- proportion
        proportion <- (proportion + limitsup) / 2
      }else{
        limitsup <- proportion
        proportion <- (proportion + limitsinf) / 2
      }
      
      # Update
      # Current proportion for rest of categories
      remaining_proportion <- 1 - proportion
      
      # Divide remaining proportion by allocations
      divided_proportion <- remaining_proportion / allocation
      
      # Undefined objects
      propinf <- 1 / sample_size
      propsup <- proportion
      E <- divided_proportion / reduction_factor
      
      # Loop through
      for(j in 1:i){
        cases[round(sample_size * propinf):round(sample_size * propsup)] <- j
        E <- E * reduction_factor
        propinf <- propsup
        propsup <- propinf + E
      }
      
      # Compute skewness
      skew_actual <- psych::skew(cases)
      
    }
    
    # Set E
    E <- divided_proportion / reduction_factor
    cumulative_probability <- proportion
    
    # Update matrix
    for(j in 1:i){
      
      skew_matrix[i,j] <- cumulative_probability
      E <- E * reduction_factor
      cumulative_probability <- cumulative_probability + E
      
    }
    
  }
  
  # Normal inverse
  norm_inv_matrix <- qnorm(skew_matrix)
  
  # Set infinite values to zero
  norm_inv_matrix[is.infinite(norm_inv_matrix)] <- 0
  
  # Category probability
  category_probability <- matrix(
    0, nrow = categories, ncol = categories
  )
  
  # Fill first column
  category_probability[,1] <- skew_matrix[,1]
  
  # Loop through
  for(i in 2:categories){
    category_probability[,i] <- skew_matrix[,i] - skew_matrix[,i-1]
  }
  
  # Make -1 = 0
  category_probability[category_probability == -1] <- 0
  
  # Return skew
  result <- list(
    skew_matrix = norm_inv_matrix,
    probability = category_probability
  )
  return(result)
  
}

#%%%%%%%%%%%%%%%%%%%%%
# ERROR FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%

#' @noRd
# Error for input type
# Updated 08.08.2022
type_error <- function(input, expected_type){
  
  # Check for type
  if(!is(input, expected_type)){
    stop(
      paste(
        "Input into '", deparse(substitute(input)),
        "' argument is not '", expected_type,
        "'. Input is ", paste("'", class(input), "'", sep = "", collapse = ", "),
        sep = ""
      )
    )
  }

}

#' @noRd
# Error for input length
# Updated 08.08.2022
length_error <- function(input, expected_lengths){
  
  # Check for length of input in expected length
  if(!length(input) %in% expected_lengths){
    stop(
      paste(
        "Length of '", deparse(substitute(input)),
        "' (", length(input),") does not match expected length(s). Length must be: ",
        paste("'", expected_lengths, "'", collapse = " or ", sep = ""),
        sep = ""
      )
    )
  }
  
}

#' @noRd
# Error for input range
# Updated 05.09.2022
range_error <- function(input, expected_ranges){
  
  # Obtain expected maximum and minimum values
  expected_maximum <- max(expected_ranges)
  expected_minimum <- min(expected_ranges)
  
  # Obtain maximum and minimum values
  actual_maximum <- round(max(input), 3)
  actual_minimum <- round(min(input), 3)
  
  # Check for maximum of input in expected range
  if(actual_maximum > expected_maximum){
    stop(
      paste(
        "Maximum of '", deparse(substitute(input)),
        "' (", actual_maximum,") does not match expected range(s). Range must be between: ",
        paste0("'", expected_ranges, "'", collapse = " and "),
        sep = ""
      )
    )
  }
  
  # Check for maximum of input in expected range
  if(actual_minimum < expected_minimum){
    stop(
      paste(
        "Minimum of '", deparse(substitute(input)),
        "' (", actual_minimum,") does not match expected range(s). Range must be between: ",
        paste0("'", expected_ranges, "'", collapse = " and "),
        sep = ""
      )
    )
  }

}


#%%%%%%%%%%%%%%%%%%%%%%
# SYSTEM FUNCTIONS ----
#%%%%%%%%%%%%%%%%%%%%%%

#' Error report
#' 
#' @description Gives necessary information for user reporting error
#' 
#' @param result Character.
#' The error from the result
#' 
#' @param SUB_FUN Character.
#' Sub-routine the error occurred in
#' 
#' @param FUN Character.
#' Main function the error occurred in
#' 
#' @return Error and message to send to GitHub
#' 
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#' 
#' @noRd
#' 
#' @importFrom utils packageVersion
#' 
# Error Report
# Updated 08.08.2022
error.report <- function(result, SUB_FUN, FUN)
{
  # Let user know that an error has occurred
  message(paste("\nAn error has occurred in the '", SUB_FUN, "' function of '", FUN, "':\n", sep =""))
  
  # Give them the error to send to you
  cat(paste(result))
  
  # Tell them where to send it
  message("\nPlease open a new issue on GitHub (bug report): https://github.com/hfgolino/EGAnet/issues/new/choose")
  
  # Give them information to fill out the issue
  OS <- as.character(Sys.info()["sysname"])
  OSversion <- paste(as.character(Sys.info()[c("release", "version")]), collapse = " ")
  Rversion <- paste(R.version$major, R.version$minor, sep = ".")
  latentFactoRversion <- paste(unlist(packageVersion("latentFactoR")), collapse = ".")
  
  # Let them know to provide this information
  message(paste("\nBe sure to provide the following information:\n"))
  
  # To reproduce
  message(styletext("To Reproduce:", defaults = "bold"))
  message(paste(" ", textsymbol("bullet"), " Function error occurred in: ", SUB_FUN, " function of ", FUN, sep = ""))
  
  # R, SemNetCleaner, and SemNetDictionaries
  message(styletext("\nR and latentFactoR versions:", defaults = "bold"))
  message(paste(" ", textsymbol("bullet"), " R version: ", Rversion, sep = ""))
  message(paste(" ", textsymbol("bullet"), " latentFactoR version: ", latentFactoRversion, sep = ""))
  
  # Desktop
  message(styletext("\nOperating System:", defaults = "bold"))
  message(paste(" ", textsymbol("bullet"), " OS: ", OS, sep = ""))
  message(paste(" ", textsymbol("bullet"), " Version: ", OSversion, sep = ""))
}

#' System check for OS and RSTUDIO
#'
#' @description Checks for whether text options are available
#'
#' @param ... Additional arguments
#'
#' @return \code{TRUE} if text options are available and \code{FALSE} if not
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
# System Check
# Updated 08.09.2020
system.check <- function (...)
{
  OS <- unname(tolower(Sys.info()["sysname"]))
  
  RSTUDIO <- ifelse(Sys.getenv("RSTUDIO") == "1", TRUE, FALSE)
  
  TEXT <- TRUE
  
  if(!RSTUDIO){if(OS != "linux"){TEXT <- FALSE}}
  
  res <- list()
  
  res$OS <- OS
  res$RSTUDIO <- RSTUDIO
  res$TEXT <- TEXT
  
  return(res)
}

#' Colorfies Text
#'
#' Makes text a wide range of colors (8-bit color codes)
#'
#' @param text Character.
#' Text to color
#'
#' @return Colorfied text
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
#'
# Color text
# Updated 08.09.2020
colortext <- function(text, number = NULL, defaults = NULL)
{
  # Check system
  sys.check <- system.check()
  
  if(sys.check$TEXT)
  {
    # Defaults for number (white text)
    if(is.null(number) || number < 0 || number > 231)
    {number <- 15}
    
    # Check for default color
    if(!is.null(defaults))
    {
      # Adjust highlight color based on background color
      if(defaults == "highlight")
      {
        if(sys.check$RSTUDIO)
        {
          
          if(rstudioapi::getThemeInfo()$dark)
          {number <- 226
          }else{number <- 208}
          
        }else{number <- 208}
      }else{
        
        number <- switch(defaults,
                         message = 204,
                         red = 9,
                         orange = 208,
                         yellow = 11,
                         "light green" = 10,
                         green = 34,
                         cyan = 14,
                         blue = 12,
                         magenta = 13,
                         pink = 211,
        )
        
      }
      
    }
    
    return(paste("\033[38;5;", number, "m", text, "\033[0m", sep = ""))
    
  }else{return(text)}
}

#' Stylizes Text
#'
#' Makes text bold, italics, underlined, and strikethrough
#'
#' @param text Character.
#' Text to stylized
#'
#' @return Sytlized text
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
# Style text
# Updated 08.09.2020
styletext <- function(text, defaults = c("bold", "italics", "highlight",
                                         "underline", "strikethrough"))
{
  # Check system
  sys.check <- system.check()
  
  if(sys.check$TEXT)
  {
    if(missing(defaults))
    {number <- 0
    }else{
      
      # Get number code
      number <- switch(defaults,
                       bold = 1,
                       italics = 3,
                       underline = 4,
                       highlight = 7,
                       strikethrough = 9
      )
      
    }
    
    return(paste("\033[", number, ";m", text, "\033[0m", sep = ""))
  }else{return(text)}
}

#' Text Symbols
#'
#' Makes text symbols (star, checkmark, square root)
#'
#' @param symbol Character.
#'
#' @return Outputs symbol
#'
#' @author Alexander Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
# Symbols
# Updated 24.04.2020
textsymbol <- function(symbol = c("alpha", "beta", "chi", "delta",
                                  "eta", "gamma", "lambda", "omega",
                                  "phi", "pi", "rho", "sigma", "tau",
                                  "theta", "square root", "infinity",
                                  "check mark", "x", "bullet")
)
{
  # Get number code
  sym <- switch(symbol,
                alpha = "\u03B1",
                beta = "\u03B2",
                chi = "\u03C7",
                delta = "\u03B4",
                eta = "\u03B7",
                gamma = "\u03B3",
                lambda = "\u03BB,",
                omega = "\u03C9",
                phi = "\u03C6",
                pi = "\u03C0",
                rho = "\u03C1",
                sigma = "\u03C3",
                tau = "\u03C4",
                theta = "\u03B8",
                "square root" = "\u221A",
                infinity = "\u221E",
                "check mark" = "\u2713",
                x = "\u2717",
                bullet = "\u2022"
  )
  
  return(sym)
}

