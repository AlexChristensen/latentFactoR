#' Estimate Number of Dimensions using Factor Forest
#'
#' Estimates the number of dimensions in data using the
#' pre-trained Random Forest model from Goretzko and Buhner
#' (2020, 2022). See examples to get started
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
#' @param maximum_factors Numeric (length = 1).
#' Maximum number of factors to search over.
#' Defaults to \code{8}
#' 
#' @param PA_correlation Character (length = 1).
#' Type of correlation used in \code{\link[psych]{fa.parallel}}.
#' Must be set:
#' 
#' \itemize{
#' 
#' \item{\code{"cor"}}
#' {Pearson's correlation}
#' 
#' \item{\code{"poly"}}
#' {Polychoric correlation}
#' 
#' \item{\code{"tet"}}
#' {Tetrachoric correlation}
#' 
#' }
#' 
#' @return Returns a list containing:
#' 
#' \item{dimensions}{Number of dimensions identified}
#' 
#' \item{probabilities}{Probability that the number of dimensions
#' is most likely}
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
#' # Perform Factor Forest
#' factor_forest(two_factor$data)}
#' 
#' @author
#' # Authors of Factor Forest \cr
#' David Goretzko and Markus Buhner
#' 
#' # Authors of {latentFactoR} \cr
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Hudson Golino <hfg9s@virginia.edu>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
#' 
#' @references
#' Goretzko, D., & Buhner, M. (2022).
#' Factor retention using machine learning with ordinal data.
#' \emph{Applied Psychological Measurement}, 01466216221089345.
#' 
#' Goretzko, D., & Buhner, M. (2020).
#' One model to rule them all? Using machine learning
#' algorithms to determine the number of factors
#' in exploratory factor analysis.
#' \emph{Psychological Methods}, \emph{25}(6), 776-786.
#'
#' @import xgboost
#' @import mlr
#' @importFrom stats sd predict
#'
#' @export
#'
# Factor Forest
# Updated 19.12.2022
factor_forest <- function(
    data, sample_size,
    maximum_factors = 8,
    PA_correlation = c("cor", "poly", "tet")
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
  
  # Check for appropriate maximum factors
  type_error(maximum_factors, "numeric");
  length_error(maximum_factors, 1);
  range_error(maximum_factors, c(1, Inf));
  
  # REORGANIZE EVENTUALLY...
  
  newdata <- data
  
  # Message user
  message("Generating features...")
  
  # calculate all necessary features:
  
  N <- sample_size
  p <- ncol(newdata)
  # dat_cor <- cor(newdata, use = use)
  dat_cor <- correlation
  eigval <- eigen(dat_cor)$values
  vareig <- cumsum(eigval)/p
  
  # eigenvalue features
  
  eiggreater1 <- sum(eigval > 1)  
  releig1 <- eigval[1]/p
  releig2 <- sum(eigval[1:2])/p
  releig3 <- sum(eigval[1:3])/p
  eiggreater07 <- sum(eigval > 0.7)
  sdeigval <- sd(eigval)
  var50 <- min(which(vareig > 0.50))
  var75 <- min(which(vareig > 0.75))
  
  # matrix norm features
  
  onenorm <- norm(dat_cor,"O")
  frobnorm <- norm(dat_cor,"F")
  maxnorm <- norm(dat_cor-diag(p),"M")
  avgcor <- sum(abs(dat_cor-diag(p)))/(p*(p-1))
  specnorm <- sqrt(eigen(t(dat_cor)%*%dat_cor)$values[1])
  
  smlcor <- sum(dat_cor <= 0.1)
  avgcom <- mean(psych::smc(dat_cor))
  det <- det(dat_cor)
  
  KMO <- psych::KMO(dat_cor)$MSA
  Gini <- ineq::ineq(lower.tri(dat_cor), type = "Gini")
  Kolm <- ineq::ineq(lower.tri(dat_cor), type = "Kolm")
  
  # parallel analysis
  
  sink <- capture.output(
    pa <- psych::fa.parallel(
      x = newdata,
      fa = "fa", cor = PA_correlation,
      plot = FALSE
    )
  )
  pa_solution <- pa$nfact
  fa_eigval <- pa$fa.values
  
  
  # empirical kaiser criterion

  # Obtain number of variables
  variables <- ncol(correlation)
  
  # Obtain eigenvalues
  eigenvalues <- eigen(correlation)$values
  
  # Initialize reference criterion
  reference <- numeric(length = variables)
  
  # Loop over variables
  for(i in 1:variables){
    reference[i] <- max(
      ((1 + sqrt(variables / sample_size))^2) *
        (variables - sum(reference)) / # note difference with `EKC`
        (variables - i + 1), 1
    )
  }
  # Model was trained with `reference` but proper uses `eigenvalues`
  
  
  # Identify last eigenvalue greater than reference
  ekc <- max(which(eigenvalues >= reference))
  
  # setting missing eigenvalues to -1000
  
  eigval[(length(eigval)+1):80] <- -1000
  fa_eigval[(length(fa_eigval)+1):80] <- -1000
  names(eigval) <- paste("eigval", 1:80, sep = "")
  names(fa_eigval) <- paste("fa_eigval", 1:80, sep = "")
  
  
  cd <- EFA.Comp.Data(Data = newdata, F.Max = maximum_factors, use = "pairwise.complete.obs")
  
  # combination of features
  
  features <- cbind(data.frame(N,p,eiggreater1,releig1,releig2,releig3,eiggreater07,sdeigval,var50,var75,onenorm,frobnorm,
                               maxnorm, avgcor, specnorm, smlcor, avgcom,det, KMO, Gini, Kolm, pa_solution, ekc, cd), t(eigval), t(fa_eigval))
  
  # Check if model is stored in the environment
  if(!exists("factor_forest_model", envir = globalenv())){
    
    # Download Factor Forest model from Google Drive
    drive_link <- "1bK-lMOh2lO7sGIVxHy1jI4kr3LjOURDD"
    
    # Check if Factor Forest model exists
    if(
      !"factor_forest_model.rdata" %in%
      tolower(list.files(tempdir()))
    ){
      
      # Let user know Factor Forest model is downloading
      message("Downloading Factor Forest model...", appendLF = FALSE)
      
      # Download Factor Forest model
      model_file <- suppressMessages(
        googledrive::drive_download(
          googledrive::as_id(drive_link),
          path = paste(tempdir(), "factor_forest_model.Rdata", sep = "\\"),
          overwrite = TRUE
        )
      )
      
    }else{
      
      # Create dummy space file list
      model_file <- list()
      model_file$local_path <- paste(
        tempdir(), "\\",
        "factor_forest_model.RData",
        sep = ""
      )
      
    }
    
    # Let user know Factor Forest model is loading
    message("Loading Factor Forest model...", appendLF = FALSE)
    
    # Initialize object
    factor_forest_model <- NULL
    
    # Load Factor Forest model
    load(model_file$local_path)
    
    # Let user know loading is finished
    message("done")
    
  }

  # Obtain output
  sink <- capture.output(
    suppressWarnings(
      out <- predict(factor_forest_model, newdata = features)$data
    )
  )
  
  # Organize output
  ## Dimensions
  dimensions <- as.numeric(
    as.character(
      out$response
    )
  )
  
  ## Obtain probabilities
  probabilities <- out[,grep("prob.", colnames(out))]
  probabilities <- round(
    as.vector(as.matrix(probabilities)), 4
  )
  names(probabilities) <- 1:length(probabilities)
  
  # Set up return list
  results <- list(
    dimensions = dimensions,
    probabilities = probabilities
  )
  
  # Let user know that they can store model in workspace
  if(!exists("factor_forest_model", envir = globalenv())){
    
    # Adjust file path
    adjusted_path <- gsub("\\\\", "/", model_file$local_path)
    
    # Let user know how to get the model
    message(
      paste0(
        "\nFor most efficient use, you can load the model ",
        "into your workspace. Run the following code to ",
        "load the model:\n\n",
        "load(\"", adjusted_path, "\")\n" 
      )
    )
  }
  
  # Return results
  return(results)
  
}
