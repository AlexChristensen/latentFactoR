#' Estimates Dimensions using Several State-of-the-art Methods
#'
#' Estimates dimensions using Exploratory Graph Analysis
#' (\code{\link[EGAnet]{EGA}}), Exploratory Factor Analysis
#' with out-of-sample prediction (\code{\link[fspe]{fspe}}),
#' Next Eigenvalue Sufficiency Test (\code{\link[latentFactoR]{NEST}}),
#' parallel analysis (\code{\link[psych]{fa.parallel}})
#' 
#' @param data Matrix or data frame.
#' Either a dataset with all numeric values
#' (rows = cases, columns = variables) or
#' a symmetric correlation matrix
#' 
#' @param sample_size Numeric (length = 1).
#' If input into \code{data} is a correlation matrix,
#' then specifying the sample size is required
#' 
#' @param EGA_args List.
#' List of arguments to be passed along to
#' \code{\link[EGAnet]{EGA}}.
#' Defaults are listed
#' 
#' @param FSPE_args List.
#' List of arguments to be passed along to
#' \code{\link[fspe]{fspe}}.
#' Defaults are listed
#' 
#' @param NEST_args List.
#' List of arguments to be passed along to
#' \code{\link[latentFactoR]{NEST}}.
#' Defaults are listed
#' 
#' @param PA_args List.
#' List of arguments to be passed along to
#' \code{\link[psych]{fa.parallel}}.
#' Defaults are listed
#' 
#' @return Returns a list containing:
#' 
#' \item{dimensions}{Dimensions estimated from each method}
#' 
#' A list of each methods output (see their respective functions for their outputs)
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
#' \dontrun{
#' # Estimate dimensions
#' estimate_dimensions(two_factor$data)
#' }
#' 
#' @author
#' Maria Dolores Nieto Canaveras <mnietoca@nebrija.es>,
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Hudson Golino <hfg9s@virginia.edu>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
#' 
#' @importFrom utils capture.output
#'
#' @export
#'
# Main factor simulation function
# Updated 05.09.2022
estimate_dimensions <- function(
    data, sample_size,
    EGA_args = list(
      corr = "cor_auto", uni.method = "louvain",
      model = "glasso", algorithm = "walktrap",
      consensus.method = "most_common",
      plot.EGA = FALSE
    ),
    FSPE_args = list(
      maxK = 8, rep = 1, method = "PE", pbar = FALSE
    ),
    NEST_args = list(
      iterations = 1000,
      maximum_iterations = 500,
      alpha = 0.05,
      convergence = 0.00001
    ),
    PA_args = list(
      fm = "minres", fa = "both", n.iter = 20,
      sim = TRUE, plot = FALSE
    )
)
{
  
  # Check for appropriate data
  object_error(data, c("matrix", "data.frame", "array"));
  sink <- apply(data, 2, type_error, expected_type = "numeric");
  
  # Ensure data is matrix
  data <- as.matrix(data)
  
  # Check for variable names
  if(is.null(colnames(data))){
    colnames(data) <- paste0("V", 1:ncol(data))
  }
  
  # Obtain correlation matrix (if not already)
  if(!isSymmetric(data)){
    
    # Compute correlations
    correlation <- qgraph::cor_auto(
      data, forcePD = TRUE, verbose = FALSE
    )
    
    # Set sample size
    sample_size <- nrow(data)
    
  }else{
    
    # Set data as correlations
    correlation <- data
    
    # Check for sample size
    if(missing(sample_size)){
      stop("Input for 'sample_size' is required when input for 'data' is a correlation (symmetric) matrix.")
    }
    
  }
  
  # Check for appropriate sample size
  type_error(sample_size, "numeric"); length_error(sample_size, 1);
  range_error(sample_size, c(2, Inf))
  
  # Estimate dimensions
  ## EGA
  message("Estimating EGA...", appendLF = FALSE)
  
  ## Obtain EGA arguments
  EGA_args <- obtain_arguments(
    FUN = EGAnet::EGA,
    FUN_args = EGA_args
  )
  
  ## Set data and sample size
  EGA_args$data <- correlation
  EGA_args$n <- sample_size
  
  ## Estimate EGA
  ega_results <- suppressWarnings(
    do.call(
      what = EGAnet::EGA,
      args = EGA_args
    )
  )
  
  ## EGA finished
  message("done.")
  
  ## Check for whether FSPE can be estimated
  if(isSymmetric(data)){
    
    ## FSPE
    warning("FSPE cannot be estimated. FSPE requires data rather than a correlation (symmetric) matrix")
    
    ## Set results to NULL
    fspe_results <- list(
      nfactor <- NA
    )
    
  }else{
    
    ## FSPE
    message("Estimating FSPE...", appendLF = FALSE)
    
    ## Obtain EGA arguments
    FSPE_args <- obtain_arguments(
      FUN = fspe::fspe,
      FUN_args = FSPE_args
    )
    
    ## Set data and sample size
    FSPE_args$data <- data
    
    ## Estimate FSPE
    fspe_results <- suppressWarnings(
      do.call(
        what = fspe::fspe,
        args = FSPE_args
      )
    )
    
    ## FSPE finished
    message("done.")
    
  }
  
  ## NEST
  message("Estimating NEST...", appendLF = FALSE)
  
  ## Obtain NEST arguments
  NEST_args <- obtain_arguments(
    FUN = NEST,
    FUN_args = NEST_args
  )
  
  ## Set data and sample size
  NEST_args$data <- correlation
  NEST_args$sample_size <- sample_size
  
  ## Estimate NEST
  nest_results <- suppressWarnings(
    do.call(
      what = NEST,
      args = NEST_args
    )
  )
  
  ## NEST finished
  message("done.")
  
  ## Parallel Analysis
  message("Estimating parallel analysis...", appendLF = FALSE)
  
  ## Obtain Parallel Analysis arguments
  PA_args <- obtain_arguments(
    FUN = psych::fa.parallel,
    FUN_args = PA_args
  )
  
  ## Set data and sample size
  PA_args$x <- correlation
  PA_args$n.obs <- sample_size
  
  ## Estimate Parallel Analysis
  sink <- capture.output(
    pa_results <- suppressWarnings(
      do.call(
        what = psych::fa.parallel,
        args = PA_args
      )
    )
  )
  
  ## Check for dimensions
  ## PAF
  if("nfact" %in% names(pa_results)){
    pa_nfact <- pa_results$nfact
  }
  ## PCA
  if("ncomp" %in% names(pa_results)){
    pa_ncomp <- pa_results$ncomp
  }
  
  ## Parallel Analysis finished
  message("done.")
  
  # Set up results vector
  dimension_results <- c(
    EGA = ega_results$n.dim,
    FSPE = fspe_results$nfactor,
    NEST = nest_results$dimensions,
    PA_PAF = ifelse(
      exists("pa_nfact"), pa_nfact, NA
    ),
    PA_PCA = ifelse(
      exists("pa_ncomp"), pa_ncomp, NA
    )
  )
  
  # Set up results list
  results <- list(
    dimensions = dimension_results,
    EGA = ega_results,
    FSPE = fspe_results,
    NEST = nest_results,
    PA = pa_results
  )
  
  # Set class
  class(results) <- "lf_estimate"
  
  # Return
  return(results)
  
  
}
