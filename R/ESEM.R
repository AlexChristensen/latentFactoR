#' Estimates Exploratory Structural Equation Model
#'
#' A general function to estimate an Exploratory Structural
#' Equation Model (ESEM) using the \code{\link{lavaan}} package.
#' With \code{\link{latentFactoR}} objects,
#' the function requires fewer inputs
#'
#' @param data Numeric matrix, data frame, or \code{\link{latentFactoR}} object
#'
#' @param factors Numeric (length = 1).
#' Number of ESEM factors to estimate
#'
#' @param variables Numeric (length = 1 or \code{factors}).
#' Number of variables per factor. A vector the length of the
#' number of factors can be specified to allow varying
#' number of variables on each factor (necessary for some
#' \code{wording_factor} arguments)
#'
#' @param estimator Character.
#' Estimator to be used in \code{\link[lavaan]{cfa}}.
#' Default options are \code{"MLR"} for continuous data
#' and \code{"WLSMV"} for categorical data
#'
#' @param fit_measures Character.
#' Fit measures to be computed using \code{\link[lavaan]{fitMeasures}}.
#' Defaults to: \code{"chisq"}, \code{"df"}, \code{"pvalue"}, \code{"cfi"},
#' \code{"tli"}, \code{"rmsea"}, \code{"rmsea.ci.lower"},
#' \code{"rmsea.ci.upper"}, \code{"rmsea.pvalue"}, and \code{"srmr"}.
#' Other measures can be added but these measures will always be produced.
#'
#' If scaled values are available (not \code{NA}), then scaled fit measures
#' will be used.
#'
#' @param variable_polarity Numeric/character (length = 1 or total variables).
#' Whether all (length = 1) or each variable (length = total variables) are
#' positive (\code{1}, \code{"p"}, \code{"pos"}, \code{"positive"}) or
#' negative (\code{-1}, \code{"n"}, \code{"neg"}, \code{"negative"})
#' polarity on the factor
#'
#' @param wording_factor Character (length = 1).
#' Whether wording factor(s) should be estimated.
#' Defaults to \code{"none"}.
#' Options include:
#'
#' \itemize{
#'
#' \item{"CTCM1"}
#' {Description coming soon...}
#'
#' \item{"CTCM1_each"}
#' {Description coming soon...}
#'
#' \item{"RI"}
#' {Description coming soon...}
#'
#' \item{"RI_each"}
#' {Description coming soon...}
#'
#' }
#'
#' @param CTCM1_polarity Character.
#' Polarity of the CTCM1 wording factor(s).
#' Defaults to \code{"negative"} for negative
#' polarity variables
#'
#' @param ... Additional arguments to be passed on to
#' \code{\link[lavaan]{cfa}}
#'
#' @return Returns a list containing:
#'
#' \item{model}{Estimated ESEM model}
#'
#' \item{fit}{Fit measures of estimated ESEM model}
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
#' \dontrun{
#' # Estimate ESEM model with no wording effects
#' esem_no_wording_effects <- ESEM(
#'   data = two_factor,
#'   estimator = "WLSMV"
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
#' # Estimate ESEM model with wording effects
#' esem_wording_effects <- ESEM(
#'   data = two_factor_acquiescence,
#'   estimator = "WLSMV",
#'   wording_factor = "RI"
#' )}
#'
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
#'
#' @export
#'
# Exploratory Structural Equation Modeling
# Updated 06.12.2022
ESEM <- function(
    data, factors, variables,
    estimator = c("MLR", "WLSMV"),
    fit_measures = NULL,
    variable_polarity = NULL,
    wording_factor = c(
      "none", "CTCM1", "CTCM1_each", "RI", "RI_each"
    ),
    CTCM1_polarity = c("negative", "positive"),
    ...
)
{

  # Check for {latentFactoR} class
  if(is(data, "lf_simulate")){

    # Switch data to `lf_object`
    lf_object <- data

    # Obtain values
    data <- lf_object$data

    # Obtain number of variables and factors
    if(is(lf_object, "lf_we")){
      parameters <- lf_object$adjusted_results$parameters
    }else{
      parameters <- lf_object$parameters
    }

    # Obtain number of variables
    variables <- parameters$variables

    # Obtain number of factors
    factors <- length(variables)

  }

  # Check for missing wording factor
  if(missing(wording_factor)){
    wording_factor <- "none"
  }else{
    wording_factor <- tolower(match.arg(wording_factor))
  }

  # Ensure data is a matrix
  data <- as.matrix(data)

  # Ensure data has column names (see `utils-latentFactoR` for function)
  data <- ensure_column_names(data)

  # Check for appropriate types
  apply(data, 2, type_error, "numeric"); type_error(factors, "numeric");
  type_error(variables, "numeric"); type_error(wording_factor, "character");

  # Check for appropriate lengths
  length_error(factors, 1); length_error(variables, c(1, factors));

  # Check for appropriate ranges
  range_error(factors, c(1, Inf)); range_error(variables, c(3, Inf))

  # Total number of variables
  total_variables <- ncol(data)

  # Set up lavaan model
  for(i in 1:factors){

    # Factor set up
    factor_syntax <- paste(
      'efa("block1")*',
      "Fac", i, " =~", sep = ""
    )

    # Variable set up
    variable_syntax <- paste0(
      "V", formatC(
        1:total_variables, format = "d",
        flag = "0", digits = 1
      ), sep = "", collapse = " + "
    )

    # Combine set up
    if(i == 1){
      model <- paste(factor_syntax, variable_syntax, "\n")
    }else{
      model <- paste(model, factor_syntax, variable_syntax, "\n")
    }

  }

  # Check for CTCM1
  if(wording_factor == "ctcm1" | wording_factor == "ctcm1_each"){

    # Check for CTCM1 polarity
    if(missing(CTCM1_polarity)){
      CTCM1_polarity <- "negative"
    }else{
      CTCM1_polarity <- match.arg(CTCM1_polarity)
    }

    # Ensure proper polarity lengths
    length_error(variable_polarity, c(0, 1, sum(variables)));

    # Check for variable polarity
    if(is.null(variable_polarity)){

      # Check for parameters ({latentFactoR} class)
      if(exists("parameters")){

        # Obtain loadings
        loadings <- parameters$loadings

        # Set sequence of variables for each factor
        end_variables <- cumsum(variables)
        start_variables <- (end_variables + 1) - variables

        # Initialize signs
        variable_polarity <- numeric(nrow(loadings))

        # Make all dominant loadings positive
        for(i in 1:ncol(loadings)){

          # Target dominant loadings
          target_loadings <- start_variables[i]:end_variables[i]

          # Determine sign
          variable_polarity[target_loadings] <- sign(loadings[target_loadings, i])

        }

      }else{

        stop(
          paste(
            "CTCM1 methods require that each variable's polarity",
            '(e.g., 1/-1, "p"/"n", "pos"/"neg", "positive"/"negative")',
            "be specified."
          )
        )

      }

    }

    # Set variable polarity as characters
    variable_polarity <- tolower(as.character(variable_polarity))

    # Homogenize variable polarity
    ## Positive
    positive_possible_polarity <- c("1", "p", "pos", "positive")

    # Loop through possible polarities
    for(possible in positive_possible_polarity){
      variable_polarity <- ifelse(
        variable_polarity == possible,
        "1", variable_polarity
      )
    }

    ## Negative
    negative_possible_polarity <- c("-1", "n", "neg", "negative")

    # Loop through possible polarities
    for(possible in negative_possible_polarity){
      variable_polarity <- ifelse(
        variable_polarity == possible,
        "-1", variable_polarity
      )
    }

    # Convert variable polarity to numeric
    variable_polarity <- as.numeric(variable_polarity)

    # Set target sign
    target_sign <- ifelse(CTCM1_polarity == "negative", -1, 1)

    # Check for whether each factor gets a CTCM1 factor
    if(wording_factor == "ctcm1"){

      # Configure CTCM1 factor
      ctcm1_model <- paste(
        'CTCM1 =~',
        paste(
          "1*",
          colnames(data)[which(variable_polarity == target_sign)],
          sep = "", collapse = " + "
        ), "\n"
      )

      # Add CTCM1 factor
      model <- paste(
        model,
        ctcm1_model,
        collapse = ""
      )

      # Make factors orthogonal
      orthogonal_correlations <- paste0(
        "CTCM1 ~~ 0*Fac", 1:factors, collapse = " \n "
      )

      # Add orthogonal factors
      model <- paste(
        model,
        orthogonal_correlations,
        collapse = ""
      )

    }else if(wording_factor == "ctcm1_each"){

      # Set sequence of variables for each factor
      end_variables <- cumsum(variables)
      start_variables <- (end_variables + 1) - variables

      # Set up lavaan model
      for(i in 1:factors){

        # Factor set up
        factor_syntax <- paste(
          "CTCM1", i, " =~", sep = ""
        )

        # Target variables
        target_variables <- start_variables[i]:end_variables[i]

        # Check that target sign exists in variables
        if(any(variable_polarity[target_variables] == target_sign)){

          # Variable set up
          variable_syntax <- paste0(
            "1*V", formatC(
              target_variables[which(variable_polarity[target_variables] == target_sign)],
              format = "d",
              flag = "0", digits = 1
            ), sep = "", collapse = " + "
          )

        }

        # Combine set up
        if(i == 1){
          ctcm1_model <- paste(factor_syntax, variable_syntax, "\n")
        }else{
          ctcm1_model <- paste(ctcm1_model, factor_syntax, variable_syntax, "\n")
        }

      }

      # Add CTCM1 factors
      model <- paste(
        model,
        ctcm1_model,
        collapse = ""
      )

      # Make factors orthogonal
      for(i in 1:factors){

        # Make factor orthogonal
        orthogonal_correlations <- paste0(
          paste0("CTCM1", i, " ~~ 0*Fac", collapse = ""),
          1:factors, collapse = " \n "
        )

        # Add orthogonal factors
        model <- paste(
          model,
          orthogonal_correlations,
          sep = " \n "
        )

      }

    }

  }else if(wording_factor == "ri" | wording_factor == "ri_each"){

    # Check for whether each factor gets a random-intercept factor
    if(wording_factor == "ri"){

      # Configure random-intercept factor
      ri_model <- paste(
        'RI =~',
        paste(
          "1*", colnames(data), sep = "", collapse = " + "
        ), "\n"
      )

      # Add random-intercept factor
      model <- paste(
        model,
        ri_model,
        collapse = ""
      )

      # Make factors orthogonal
      orthogonal_correlations <- paste0(
        "RI ~~ 0*Fac", 1:factors, collapse = " \n "
      )

      # Add orthogonal factors
      model <- paste(
        model,
        orthogonal_correlations,
        collapse = ""
      )

    }else if(wording_factor == "ri_each"){

      # Set sequence of variables for each factor
      end_variables <- cumsum(variables)
      start_variables <- (end_variables + 1) - variables

      # Set up lavaan model
      for(i in 1:factors){

        # Factor set up
        factor_syntax <- paste(
          "RI", i, " =~", sep = ""
        )

        # Variable set up
        variable_syntax <- paste0(
          "1*V", formatC(
            start_variables[i]:end_variables[i], format = "d",
            flag = "0", digits = 1
          ), sep = "", collapse = " + "
        )

        # Combine set up
        if(i == 1){
          ri_model <- paste(factor_syntax, variable_syntax, "\n")
        }else{
          ri_model <- paste(ri_model, factor_syntax, variable_syntax, "\n")
        }

      }

      # Add random-intercept factor
      model <- paste(
        model,
        ri_model,
        collapse = ""
      )

      # Make factors orthogonal
      for(i in 1:factors){

        # Make factor orthogonal
        orthogonal_correlations <- paste0(
          paste0("RI", i, " ~~ 0*Fac", collapse = ""),
          1:factors, collapse = " \n "
        )

        # Add orthogonal factors
        model <- paste(
          model,
          orthogonal_correlations,
          sep = " \n "
        )

      }

    }

  }

  # Set "ordered" argument
  ordered <- estimator == "WLSMV"
  # In the future, set ordered variables
  # to be only categorical variables
  # (make function work with mixed
  # [continuous and categorical] data)

  # Set "std.lv" argument
  std.lv <- wording_factor == "none"

  # Perform ESEM
  esem_mod <- try(
    lavaan::cfa(
      model = model, data = data,
      estimator = estimator, std.lv = std.lv,
      ordered = ordered
    ),
    silent = TRUE
  )

  # Check for error
  if(is(esem_mod, "try-error")){

    # Return results
    return(
      list(
        model = "Model could not converge.",
        fit = NA
      )
    )

  }

  # Summary of fit
  esem_fit <- lavaan::fitMeasures(esem_mod)

  # Set fit measures
  if(is.null(fit_measures)){

    fit_measures <- c(
      "chisq", "df", "pvalue",
      "cfi", "tli",
      "rmsea", "rmsea.ci.lower",
      "rmsea.ci.upper", "rmsea.pvalue"
    )

  }else{

    fit_measures <- unique(c(
      fit_measures,
      "chisq", "df", "pvalue",
      "cfi", "tli",
      "rmsea", "rmsea.ci.lower",
      "rmsea.ci.upper", "rmsea.pvalue"
    ))

  }

  # Determine whether scaled measures are available
  scaled <- !is.na(esem_fit["chisq.scaled"])

  # Obtain scaled fit measures
  if(isTRUE(scaled)){

    # Add scaled to measures
    fit_measures <- paste0(fit_measures, ".scaled")

    # Check if scaled exist
    no_scaled_exists <- is.na(match(fit_measures, names(esem_fit)))

    # Remove scale from measures that don't have it
    if(any(no_scaled_exists)){

      fit_measures[which(no_scaled_exists)] <-
        gsub(
          ".scaled", "", fit_measures[which(no_scaled_exists)]
        )

    }

  }

  # Add SRMR
  fit_measures <- unique(c(fit_measures, "srmr"))

  # Obtain fit measures
  final_fit <- round(esem_fit[fit_measures], 3)

  # Return results
  return(
    list(
      model = esem_mod,
      fit = final_fit
    )
  )

}
