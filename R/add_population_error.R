#' Adds Population Error to Factor Model Data
#'
#' Adds population error to simulated data from \code{\link[latentFactoR]{simulate_factors}}. 
#' See examples to get started
#' 
#' @param lf_object Data object from \code{\link[latentFactoR]{simulate_factors}}
#' 
#' @param cfa_method Character (length = 1).
#' Method to generate population error.
#' Defaults to \code{"minres"}.
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"minres"}}
#' {Minimum residual}
#' 
#' \item{\code{"ml"}}
#' {Maximum likelihood}
#' 
#' }
#' 
#' @param fit Character (length = 1).
#' Fit index to control population error.
#' Defaults to \code{"rmsr"}.
#' Available options:
#' 
#' \itemize{
#' 
#' \item{\code{"cfi"}}
#' {Comparative fit index}
#' 
#' \item{\code{"rmsea"}}
#' {Root mean square error of approximation}
#' 
#' \item{\code{"rmsr"}}
#' {Root mean square residuals}
#' 
#' \item{\code{"raw"}}
#' {Direct application of error}
#' 
#' }
#' 
#' @param misfit Numeric (length = 1).
#' Magnitude of error to add.
#' Defaults to \code{0}.
#' General effect sizes range from small (0.10), moderate (0.20), to large (0.30)
#' 
#' @param error_method Character (length = 1).
#' Method to control population error.
#' Description of methods:
#' 
#' \itemize{
#' 
#' \item{\code{"cudeck"}}
#' {Description coming soon... see Cudeck & Browne, 1992
#' for more details}
#' 
#' \item{\code{"yuan"}}
#' {Description coming soon...}
#' 
#' }
#' 
#' @return Returns a list containing:
#' 
#' \item{data}{Simulated data from the specified factor model}
#' 
#' \item{population_correlation}{Population correlation matrix with local dependence added}
#' 
#' \item{parameters}{
#' A list containing the parameters used to generate the data:
#' 
#' \itemize{
#' 
#' \item{\code{factors}}
#' {Number of factors}
#' 
#' \item{\code{variables}}
#' {Variables on each factor}
#' 
#' \item{\code{loadings}}
#' {Loading matrix}
#' 
#' \item{\code{factor_correlations}}
#' {Correlations between factors}
#' 
#' \item{\code{categories}}
#' {Categories for each variable}
#' 
#' \item{\code{skew}}
#' {Skew for each variable}
#' 
#' }
#' 
#' }
#' 
#' \item{population_error}{
#' A list containing the parameters used to generate population error:
#' 
#' \itemize{
#' 
#' \item{\code{error_correlation}}
#' {Correlation matrix with population error added
#' (same as \code{population_correlation})}
#' 
#' \item{\code{fit}}
#' {Fit measure used to control population error}
#' 
#' \item{\code{delta}}
#' {Minimum of the objective function corresponding to the misfit value}
#' 
#' \item{\code{misfit}}
#' {Specified misfit value}
#' 
#' }
#' 
#' }
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
#' # Add small population error using Yuan method
#' two_factor_Yuan <- add_population_error(
#'   lf_object = two_factor,
#'   cfa_method = "minres",
#'   fit = "rmsr", misfit = 0.10,
#'   error_method = "yuan"
#' )
#' 
#' # Add small population error using Cudeck method
#' two_factor_Yuan <- add_population_error(
#'   lf_object = two_factor,
#'   cfa_method = "minres",
#'   fit = "rmsr", misfit = 0.10,
#'   error_method = "cudeck"
#' )
#' 
#' @author
#' {\code{\link{bifactor}}} authors \cr
#' Marcos Jimenez,
#' Francisco J. Abad,
#' Eduardo Garcia-Garzon,
#' Vithor R. Franco,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
#' 
#' {\code{\link{latentFactoR}}} authors \cr
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Hudson Golino <hfg9s@virginia.edu>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
#' 
#' @references
#' Christensen, A. P., Garrido, L. E., & Golino, H. (2022). \cr
#' Unique variable analysis: A network psychometrics method to detect local dependence. \cr
#' \emph{PsyArXiv}
#' 
#' Cudeck, R., & Browne, M.W. (1992). \cr
#' Constructing a covariance matrix that yields a specified minimizer and a specified minimum discrepancy function value. \cr
#' \emph{Psychometrika}, \emph{57}, 357–369.
#' 
#' Jimenez, M., Abad, F. J., Garcia-Garzon, E., Golino, H., Christensen, A. P., & Garrido, L. E. (2022). \cr
#' Dimensionality assessment in generalized bi-factor structures: A network psychometrics approach. \cr
#' \emph{PsyArXiv}
#'
#' @export
#'
# Add population to simulated data
# Updated 17.09.2022
add_population_error <- function(
    lf_object,
    cfa_method = c("minres", "ml"),
    fit = c("cfi", "rmsea", "rmsr", "raw"),
    misfit = 0,
    error_method = c("cudeck", "yuan")
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
  
  # Check for missing CFA method
  if(missing(cfa_method)){
    cfa_method <- "minres"
  }else{cfa_method <- tolower(match.arg(cfa_method))}
  
  # Check for missing fit
  if(missing(fit)){
    fit <- "rmsr"
  }else{fit <- tolower(match.arg(fit))}
  
  # Check for missing error method
  if(missing(error_method)){
    error_method <- "yuan"
  }else{error_method <- tolower(match.arg(error_method))}
  
  # Check for appropriate misfit
  type_error(misfit, "numeric"); length_error(misfit, 1);
  range_error(misfit, c(0, 1));
  
  # Obtain population error
  if(error_method == "cudeck"){
    
    # Using Cudeck method
    # From {bifactor} version 0.1.0
    # See `utils-latentFactoR`
    population_error <- cudeck(
      R = lf_object$population_correlation,
      lambda = parameters$loadings,
      Phi = parameters$factor_correlations,
      uniquenesses = 1 - rowSums(parameters$loadings^2),
      misfit = misfit, method = cfa_method,
      confirmatory = FALSE
    )
    
  }else if(error_method == "yuan"){
    
    # Using Yuan method
    # From {bifactor} version 0.1.0
    # See `utils-latentFactoR`
    population_error <- yuan(
      R = lf_object$population_correlation,
      lambda = parameters$loadings,
      Phi = parameters$factor_correlations,
      uniquenesses = 1 - rowSums(parameters$loadings^2),
      fit = fit, misfit = misfit, method = cfa_method,
      confirmatory = FALSE
    )
    
  }
  
  # Re-estimate data
  ## Cholesky decomposition
  cholesky <- chol(population_error$R_error)
  
  ## Obtain sample size and total variables
  sample_size <- nrow(lf_object$data)
  total_variables <- ncol(lf_object$data)
  
  ## Generate data
  data <- mvtnorm::rmvnorm(sample_size, sigma = diag(total_variables))
  
  ## Make data based on factor structure
  data <- data %*% cholesky
  
  ## Set variable categories, limit, and skew
  variable_categories <- parameters$categories
  categorical_limit <- parameters$categorical_limit
  skew <- parameters$skew
  
  ## Ensure appropriate type and length for categories
  type_error(variable_categories, "numeric")
  length_error(variable_categories, c(1, total_variables))
  
  ## Identify categories to variables
  if(length(variable_categories) == 1){
    variable_categories <- rep(variable_categories, total_variables)
  }
  
  ## Check for categories greater than categorical limit and not infinite
  if(any(variable_categories > categorical_limit & !is.infinite(variable_categories))){
    
    ## Make variables with categories greater than 7 (or categorical_limit) continuous
    variable_categories[
      variable_categories > categorical_limit & !is.infinite(variable_categories)
    ] <- Inf
    
  }
  
  ## Find categories
  if(any(variable_categories <= categorical_limit)){
    
    ## Target columns to categorize
    columns <- which(variable_categories <= categorical_limit)
    
    ## Set skew
    if(length(skew) != length(columns)){
      skew <- sample(skew, length(columns), replace = TRUE)
    }
    
    ## Loop through columns
    for(i in columns){
      
      data[,i] <- categorize(
        data = data[,i],
        categories = variable_categories[i],
        skew_value = skew[i]
      )
      
    }
    
  }
  
  ## Add column names to data
  colnames(data) <- paste0(
    "V", formatC(
      x = 1:total_variables,
      digits = floor(log10(total_variables)),
      flag = "0", format = "d"
    )
  )
  
  ## Populate results
  results <- list(
    data = data,
    population_correlation = population_error$R_error,
    original_correlation = lf_object$population_correlation
  )
  
  # Add parameters from population error
  error_parameters <- list(
    error_correlation = population_error$R_error,
    fit = population_error$fit,
    delta = population_error$delta,
    misfit = population_error$misfit
  )
  
  # Add population error parameters to results
  results$population_error <- error_parameters
  
  # Add original results
  results$original_results <- lf_object
  
  # Add class
  class(results) <- c(class(lf_object), "lf-pe")
  
  # Return results
  return(results)
  
}