#' Simulates Hierarchical Latent Factor Data
#'
#' Simulates data from a hierarchical (second-order) latent factor model based
#' on many manipulable parameters. Variables load onto lower-order factors, and
#' lower-order factors load onto higher-order factors. Parameters do not have
#' default values and must each be set, with the exception of \code{higher_disturbances}
#' and \code{lower_correlations} (and their \code{_range} variants), which are mutually
#' exclusive alternatives for specifying the same underlying structure (exactly one must
#' be supplied). See examples to get started
#'
#' The population correlation matrix (variable-level) implied by the hierarchical factor model is:
#'
#' \deqn{\boldsymbol{\Sigma} = \boldsymbol{\Lambda}_L \boldsymbol{\Phi}_L \boldsymbol{\Lambda}_L' + \boldsymbol{\Theta}}{Sigma = Lambda_L * Phi_L * Lambda_L' + Theta}
#'
#' where \eqn{\boldsymbol{\Lambda}_L}{Lambda_L} is the lower-order loading matrix (variables x lower-order
#' factors; \code{lower_loadings} and \code{lower_cross_loadings}), \eqn{\boldsymbol{\Phi}_L}{Phi_L} is the
#' implied lower-order factor correlation matrix, and \eqn{\boldsymbol{\Theta}}{Theta} is a diagonal matrix
#' of uniquenesses that ensures unit variance for each variable
#'
#' \eqn{\boldsymbol{\Phi}_L}{Phi_L} is, in turn, implied by the higher-order factor structure:
#'
#' \deqn{\boldsymbol{\Phi}_L = \boldsymbol{\Lambda}_H \boldsymbol{\Phi}_H \boldsymbol{\Lambda}_H' + \boldsymbol{\Psi}}{Phi_L = Lambda_H * Phi_H * Lambda_H' + Psi}
#'
#' where \eqn{\boldsymbol{\Lambda}_H}{Lambda_H} is the higher-order loading matrix (lower-order factors x
#' higher-order factors; \code{higher_loadings} and \code{higher_cross_loadings}), \eqn{\boldsymbol{\Phi}_H}{Phi_H}
#' is the higher-order factor correlation matrix (\code{higher_correlations}), and \eqn{\boldsymbol{\Psi}}{Psi}
#' is the disturbance (residual) covariance matrix of the lower-order factors (\code{higher_disturbances}) —
#' i.e., the variance in each lower-order factor left unexplained by the higher-order factor(s). When
#' \code{lower_correlations} is supplied instead of \code{higher_disturbances}, \eqn{\boldsymbol{\Phi}_L}{Phi_L}
#' is fixed to the supplied target matrix and \eqn{\boldsymbol{\Psi}}{Psi} is solved for instead
#'
#' @param lower_factors Numeric (length = 1).
#' Number of lower-order (first-order) factors
#'
#' @param variables Numeric (length = 1 or \code{lower_factors}).
#' Number of variables per lower-order factor.
#' Can be a single value or as many values as there are lower-order factors.
#' Minimum two variables per factor
#'
#' @param variables_range Numeric (length = 2).
#' Range of variables to randomly select from a random uniform distribution.
#' Minimum three variables per factor
#'
#' @param lower_loadings Numeric or matrix (length = 1, \code{lower_factors}, or total variables x \code{lower_factors}).
#' Loadings of variables onto their lower-order factor, drawn from a random uniform distribution using +/- 0.10 of value input.
#' Can be a single value or as many values as there are lower-order factors (corresponding to the factors);
#' a vector as long as the total number of variables is \emph{not} applied per variable unless generated
#' automatically via \code{lower_loadings_range}.
#' Can also be a loading matrix. Columns must match \code{lower_factors} and rows must match total variables (\code{lower_factors} x \code{variables})
#' General effect sizes range from small (0.40), moderate (0.55), to large (0.70)
#'
#' @param lower_loadings_range Numeric (length = 2).
#' Range of lower-order loadings to randomly select from a random uniform distribution.
#' General effect sizes range from small (0.40), moderate (0.55), to large (0.70)
#'
#' @param lower_cross_loadings Numeric (length = 1 or \code{lower_factors}).
#' Cross-loadings of variables onto non-dominant lower-order factors, drawn from a random normal distribution with a mean of 0 and standard deviation of value input.
#' Can be a single value or as many values as there are lower-order factors (corresponding to the factors).
#' Unlike \code{lower_loadings}, a full loading matrix is \emph{not} supported: if a matrix is supplied,
#' only its first \code{lower_factors} elements (taken column-wise) are used as the per-factor standard deviations
#'
#' @param lower_cross_loadings_range Numeric (length = 2).
#' Range of lower-order cross-loadings to randomly select from a random uniform distribution
#'
#' @param lower_correlations Numeric or matrix (length = 1 or \code{lower_factors} x \code{lower_factors}).
#' Target correlation matrix between lower-order factors, supplied directly as an alternative to
#' \code{higher_disturbances}. Can be a single value (applied to all off-diagonal correlations) or a full
#' \code{lower_factors} x \code{lower_factors} correlation matrix. When supplied, \code{higher_disturbances}
#' is \emph{implied} rather than set directly: it is solved for as the residual needed to reproduce this
#' target correlation matrix given the higher-order loadings and correlations (\code{higher_loadings} and
#' \code{higher_correlations}), and may end up with off-diagonal values (i.e., correlated disturbances
#' between lower-order factors are permitted). The implied disturbances must still form a valid
#' (positive semi-definite) covariance matrix; if the target correlation matrix is incompatible with the
#' supplied higher-order structure (e.g., it implies negative unique variances), values are redrawn when
#' \code{higher_loadings_range} and/or \code{higher_correlations_range} are used, or an error is raised
#' when the higher-order structure is entirely fixed.
#' Mutually exclusive with \code{higher_disturbances}/\code{higher_disturbances_range}: exactly one of
#' \code{higher_disturbances}/\code{higher_disturbances_range} or
#' \code{lower_correlations}/\code{lower_correlations_range} must be supplied
#'
#' @param lower_correlations_range Numeric (length = 2).
#' Range of lower-order correlations to randomly select (per pair of lower-order factors) from a random
#' uniform distribution. Somewhat redundant with \code{lower_correlations} but more flexible.
#' Mutually exclusive with \code{higher_disturbances}/\code{higher_disturbances_range}
#'
#' @param higher_factors Numeric (length = 1).
#' Number of higher-order (second-order) factors.
#' \code{lower_factors} must be evenly divisible by \code{higher_factors}
#' unless \code{higher_loadings} is supplied as a matrix or \code{higher_variables} is supplied
#'
#' @param higher_variables Numeric (length = 1 or \code{higher_factors}).
#' Number of lower-order factors nested within each higher-order factor.
#' Can be a single value (lower-order factors are split evenly across higher-order factors,
#' equivalent to leaving \code{higher_variables} as \code{NULL}) or as many values as there are
#' higher-order factors (allowing a different number of lower-order factors to load onto each
#' higher-order factor). Values must sum to \code{lower_factors}.
#' Ignored when \code{higher_loadings} is supplied as a matrix
#'
#' @param higher_loadings Numeric or matrix (length = 1, \code{higher_factors}, or \code{lower_factors} x \code{higher_factors}).
#' Loadings of lower-order factors onto their higher-order factor, drawn from a random uniform distribution using +/- 0.10 of value input.
#' Can be a single value or as many values as there are higher-order factors (corresponding to the factors);
#' a vector as long as \code{lower_factors} is \emph{not} applied per lower-order factor unless generated
#' automatically via \code{higher_loadings_range}.
#' Can also be a loading matrix. Columns must match \code{higher_factors} and rows must match \code{lower_factors}
#' General effect sizes range from small (0.40), moderate (0.55), to large (0.70)
#'
#' @param higher_loadings_range Numeric (length = 2).
#' Range of higher-order loadings to randomly select from a random uniform distribution.
#' General effect sizes range from small (0.40), moderate (0.55), to large (0.70)
#'
#' @param higher_cross_loadings Numeric (length = 1 or \code{higher_factors}).
#' Cross-loadings of lower-order factors onto non-dominant higher-order factors, drawn from a random normal distribution with a mean of 0 and standard deviation of value input.
#' Can be a single value or as many values as there are higher-order factors (corresponding to the factors).
#' Unlike \code{higher_loadings}, a full loading matrix is \emph{not} supported: if a matrix is supplied,
#' only its first \code{higher_factors} elements (taken column-wise) are used as the per-factor standard deviations
#'
#' @param higher_cross_loadings_range Numeric (length = 2).
#' Range of higher-order cross-loadings to randomly select from a random uniform distribution
#'
#' @param higher_correlations Numeric (length = 1 or \code{higher_factors} x \code{higher_factors}).
#' Can be a single value that will be used for all correlations between higher-order factors.
#' Can also be a square matrix (\code{higher_factors} x \code{higher_factors}).
#' Negative values are currently converted to their absolute value (only positive
#' correlations between higher-order factors are supported).
#' General effect sizes range from orthogonal (0.00), small (0.30), moderate (0.50), to large (0.70)
#'
#' @param higher_correlations_range Numeric (length = 2).
#' Range of correlations between higher-order factors to randomly select from a random uniform distribution.
#' General effect sizes range from orthogonal (0.00), small (0.30), moderate (0.50), to large (0.70)
#'
#' @param higher_disturbances Numeric or matrix (length = 1, \code{lower_factors}, or \code{lower_factors} x \code{lower_factors}).
#' Disturbance (residual) variance of each lower-order factor left unexplained by the higher-order factor(s).
#' A single value or vector of length \code{lower_factors} is placed on the diagonal of an otherwise
#' zero (uncorrelated disturbances) matrix. Can also be a full \code{lower_factors} x \code{lower_factors}
#' covariance matrix to allow for correlated disturbances between lower-order factors.
#' Diagonal values must be between 0 and 1 and the matrix must be positive semi-definite.
#' Mutually exclusive with \code{lower_correlations}/\code{lower_correlations_range}: exactly one of
#' \code{higher_disturbances}/\code{higher_disturbances_range} or
#' \code{lower_correlations}/\code{lower_correlations_range} must be supplied
#'
#' @param higher_disturbances_range Numeric (length = 2).
#' Range of disturbance variances to randomly and independently select (per lower-order factor)
#' from a random uniform distribution. Resulting disturbances are uncorrelated.
#' Mutually exclusive with \code{lower_correlations}/\code{lower_correlations_range}
#'
#' @param sample_size Numeric (length = 1).
#' Number of cases to generate from a random multivariate normal distribution using
#' \code{\link[mvtnorm]{rmvnorm}}
#'
#' @param variable_categories Numeric (length = 1 or total variables (\code{lower_factors} x \code{variables})).
#' Number of categories for each variable. \code{Inf} is used for continuous variables; otherwise,
#' values reflect number of categories
#'
#' @param categorical_limit Numeric (length = 1).
#' Values greater than input value are considered continuous.
#' Defaults to \code{7} meaning that 8 or more categories are considered continuous
#' (i.e., variables are \emph{not} categorized from continuous to categorical)
#'
#' @param skew Numeric (length = 1 or categorical variables).
#' Skew to be included in categorical variables. It is randomly sampled from provided values.
#' Can be a single value or as many values as there are (total) variables.
#' Current skew implementation is between -2 and 2 in increments of 0.05.
#' Skews that are not in this sequence will be converted to their nearest
#' value in the sequence. Not recommended to use with \code{variables_range}.
#' Future versions will incorporate finer skews
#'
#' @param skew_range Numeric (length = 2).
#' Randomly selects skews within in the range.
#' Somewhat redundant with \code{skew} but more flexible.
#' Compatible with \code{variables_range}
#'
#' @return Returns a list containing:
#'
#' \item{data}{Simulated data from the specified hierarchical factor model}
#'
#' \item{population_correlation}{Population correlation matrix (variable-level)}
#'
#' \item{lower_correlation}{Population correlation matrix between lower-order factors,
#' implied by the higher-order factor structure and disturbances (or, when \code{lower_correlations}
#' is supplied, equal to that target matrix)}
#'
#' \item{continuous_data}{Simulated data before categories and skew were applied
#' (i.e., the underlying continuous data)}
#'
#' \item{parameters}{
#' A list containing the parameters used to generate the data:
#'
#' \itemize{
#'
#' \item \code{lower} --- A list of lower-order parameters:
#' \code{factors} (number of lower-order factors),
#' \code{variables} (variables on each lower-order factor),
#' \code{loadings} (loading matrix of variables onto lower-order factors),
#' \code{cross_loadings} (cross-loadings used for each lower-order factor), and
#' \code{correlations} (population correlation matrix between lower-order factors;
#' same matrix as the top-level \code{lower_correlation})
#'
#' \item \code{higher} --- A list of higher-order parameters:
#' \code{factors} (number of higher-order factors),
#' \code{variables} (lower-order factors nested within each higher-order factor),
#' \code{loadings} (loading matrix of lower-order factors onto higher-order factors),
#' \code{cross_loadings} (cross-loadings used for each higher-order factor),
#' \code{correlations} (correlations between higher-order factors), and
#' \code{disturbances} (disturbance (residual) covariance matrix of lower-order factors,
#' as supplied directly, or implied when \code{lower_correlations} is used instead)
#'
#' \item \code{categories} --- Categories for each variable
#'
#' \item \code{categorical_limit} --- Category limit used to determine continuous variables
#'
#' \item \code{skew} --- Skew for each variable
#'
#' }
#'
#' }
#'
#' @examples
#' # Generate hierarchical factor data
#' # 4 lower-order factors (2 per higher-order factor), 2 higher-order factors
#' hierarchical <- simulate_hierarchical_factors(
#'   lower_factors = 4, # lower-order factors = 4
#'   variables = 6, # variables per lower-order factor = 6
#'   lower_loadings = 0.55, # lower-order loadings = 0.45 to 0.65
#'   lower_cross_loadings = 0.05, # lower-order cross-loadings N(0, 0.05)
#'   higher_factors = 2, # higher-order factors = 2
#'   higher_loadings = 0.60, # higher-order loadings = 0.50 to 0.70
#'   higher_cross_loadings = 0.05, # higher-order cross-loadings N(0, 0.05)
#'   higher_correlations = 0.30, # correlation between higher-order factors = 0.30
#'   higher_disturbances = 0.30, # disturbance variance for each lower-order factor = 0.30
#'   sample_size = 1000 # number of cases = 1000
#' )
#'
#' # Different number of lower-order factors per higher-order factor
#' # 6 lower-order factors (2 on higher-order factor A, 4 on higher-order factor B)
#' hierarchical_uneven <- simulate_hierarchical_factors(
#'   lower_factors = 6, variables = 6,
#'   lower_loadings = 0.55, lower_cross_loadings = 0.05,
#'   higher_factors = 2, higher_variables = c(2, 4), # A = 2 factors, B = 4 factors
#'   higher_loadings = 0.60, higher_cross_loadings = 0.05,
#'   higher_correlations = 0.30, higher_disturbances = 0.30,
#'   sample_size = 1000
#' )
#'
#' # Randomly vary higher-order loadings
#' hierarchical_loadings <- simulate_hierarchical_factors(
#'   lower_factors = 4, variables = 6,
#'   lower_loadings = 0.55, lower_cross_loadings = 0.05,
#'   higher_factors = 2,
#'   higher_loadings_range = c(0.40, 0.80), # higher-order loadings = 0.40 to 0.80
#'   higher_cross_loadings = 0.05,
#'   higher_correlations = 0.30,
#'   higher_disturbances_range = c(0.20, 0.40), # disturbances = 0.20 to 0.40
#'   sample_size = 1000
#' )
#'
#' # Generate dichotomous data
#' hierarchical_dichotomous <- simulate_hierarchical_factors(
#'   lower_factors = 4, variables = 6,
#'   lower_loadings = 0.55, lower_cross_loadings = 0.05,
#'   higher_factors = 2, higher_loadings = 0.60,
#'   higher_cross_loadings = 0.05, higher_correlations = 0.30,
#'   higher_disturbances = 0.30, sample_size = 1000,
#'   variable_categories = 2 # dichotomous data
#' )
#'
#' # Supply a target lower-order correlation matrix directly instead of `higher_disturbances`;
#' # the disturbances are implied (solved for) and may end up correlated
#' hierarchical_lower_correlations <- simulate_hierarchical_factors(
#'   lower_factors = 4, variables = 6,
#'   lower_loadings = 0.55, lower_cross_loadings = 0.05,
#'   higher_factors = 2, higher_loadings = 0.60,
#'   higher_cross_loadings = 0.05, higher_correlations = 0.30,
#'   lower_correlations = 0.40, # target correlation of 0.40 between all lower-order factors
#'   sample_size = 1000
#' )
#'
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>,
#' Maria Dolores Nieto Canaveras <mnietoca@nebrija.es>,
#' Hudson Golino <hfg9s@virginia.edu>,
#' Luis Eduardo Garrido <luisgarrido@pucmm.edu>
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
#' @importFrom stats qnorm rnorm runif
#' @importFrom methods is
#'
#' @export
#'
# Main hierarchical factor simulation function
# Updated 12.07.2026
simulate_hierarchical_factors <- function(
    lower_factors, variables, variables_range = NULL,
    lower_loadings, lower_loadings_range = NULL,
    lower_cross_loadings, lower_cross_loadings_range = NULL,
    lower_correlations = NULL, lower_correlations_range = NULL,
    higher_factors, higher_variables = NULL,
    higher_loadings, higher_loadings_range = NULL,
    higher_cross_loadings, higher_cross_loadings_range = NULL,
    higher_correlations, higher_correlations_range = NULL,
    higher_disturbances = NULL, higher_disturbances_range = NULL,
    sample_size, variable_categories = Inf,
    categorical_limit = 7,
    skew = 0, skew_range = NULL
)
{

  # Check inputs
  inputs <- simulate_hierarchical_factors_errors(
    lower_factors, variables, variables_range,
    lower_loadings, lower_loadings_range,
    lower_cross_loadings, lower_cross_loadings_range,
    lower_correlations, lower_correlations_range,
    higher_factors, higher_variables, higher_loadings, higher_loadings_range,
    higher_cross_loadings, higher_cross_loadings_range,
    higher_correlations, higher_correlations_range,
    higher_disturbances, higher_disturbances_range,
    sample_size, variable_categories,
    categorical_limit, skew, skew_range
  )

  # Collect inputs
  variables <- inputs$variables; total_variables <- inputs$total_variables
  lower_loadings <- inputs$lower_loadings; lower_cross_loadings <- inputs$lower_cross_loadings
  lower_correlations <- inputs$lower_correlations
  higher_variables <- inputs$higher_variables
  higher_loadings <- inputs$higher_loadings; higher_cross_loadings <- inputs$higher_cross_loadings
  higher_correlations <- inputs$higher_correlations; higher_disturbances <- inputs$higher_disturbances
  skew <- inputs$skew

  # Ensure appropriate lengths
  if(!is.null(lower_correlations)){
    length_error(lower_correlations, c(1, lower_factors * lower_factors))
  }
  length_error(higher_correlations, c(1, higher_factors * higher_factors))
  if(!is.null(higher_disturbances)){
    length_error(higher_disturbances, c(1, lower_factors, lower_factors * lower_factors))
  }

  # Ensure appropriate ranges
  range_error(lower_loadings, c(-1, 1)); range_error(lower_cross_loadings, c(-1, 1))
  if(!is.null(lower_correlations)){range_error(lower_correlations, c(-1, 1))}
  range_error(higher_loadings, c(-1, 1)); range_error(higher_cross_loadings, c(-1, 1))
  range_error(higher_correlations, c(-1, 1))
  if(!is.null(higher_disturbances)){
    range_error(diag(as.matrix(higher_disturbances), names = FALSE), c(0, 1))
  }

  # Ensure appropriate types
  if(!is(lower_loadings, "matrix")){type_error(lower_loadings, "numeric")}
  if(!is(lower_cross_loadings, "matrix")){type_error(lower_cross_loadings, "numeric")}
  if(!is.null(lower_correlations) && !is(lower_correlations, "matrix")){
    type_error(lower_correlations, "numeric")
  }
  if(!is(higher_loadings, "matrix")){type_error(higher_loadings, "numeric")}
  if(!is(higher_cross_loadings, "matrix")){type_error(higher_cross_loadings, "numeric")}
  if(!is(higher_correlations, "matrix")){type_error(higher_correlations, "numeric")}
  if(!is.null(higher_disturbances) && !is(higher_disturbances, "matrix")){
    type_error(higher_disturbances, "numeric")
  }

  # Set lower-order factor sequence
  lower_factor_sequence <- seq_len(lower_factors)

  # Set higher-order factor sequence
  higher_factor_sequence <- seq_len(higher_factors)

  # Identify variables structure on lower-order factors
  if(length(variables) == 1){variables <- rep(variables, lower_factors)}

  # Create starting and ending of variable sequences
  variable_sequence <- factor_variable_sequence(variables)
  start_variables <- variable_sequence$start
  end_variables <- variable_sequence$end

  # Verify generalizable lengths (lower-order level)
  if(length(lower_loadings) == 1){lower_loadings <- rep(lower_loadings, lower_factors)}
  if(length(lower_cross_loadings) == 1){lower_cross_loadings <- rep(lower_cross_loadings, lower_factors)}
  if(length(variable_categories) == 1){variable_categories <- rep(variable_categories, total_variables)}

  # Identify lower-order factor structure on higher-order factors
  if(!is(higher_loadings, "matrix")){

    if(is.null(higher_variables)){

      # Ensure lower-order factors divide evenly across higher-order factors
      if(lower_factors %% higher_factors != 0){
        stop(
          paste0(
            "`lower_factors` (", lower_factors, ") must be evenly divisible by `higher_factors` (",
            higher_factors, ") when `higher_loadings` is not supplied as a matrix and ",
            "`higher_variables` is not supplied"
          )
        )
      }

      # Set number of lower-order factors nested within each higher-order factor
      higher_variables <- rep(lower_factors / higher_factors, higher_factors)

    }else if(length(higher_variables) == 1){

      # Repeat single value across higher-order factors
      higher_variables <- rep(higher_variables, higher_factors)

    }

    # Ensure supplied lower-order factor counts account for all lower-order factors
    if(sum(higher_variables) != lower_factors){
      stop(
        paste0(
          "Sum of `higher_variables` (", sum(higher_variables), ") must equal `lower_factors` (",
          lower_factors, ")"
        )
      )
    }

    # Create starting and ending of lower-order factor sequences on higher-order factors
    higher_variable_sequence <- factor_variable_sequence(higher_variables)
    start_higher <- higher_variable_sequence$start
    end_higher <- higher_variable_sequence$end

  }

  # Verify generalizable lengths (higher-order level)
  if(length(higher_loadings) == 1){higher_loadings <- rep(higher_loadings, higher_factors)}
  if(length(higher_cross_loadings) == 1){higher_cross_loadings <- rep(higher_cross_loadings, higher_factors)}

  # `higher_disturbances` and `lower_correlations` are mutually exclusive (enforced in
  # `simulate_hierarchical_factors_errors`); exactly one of the two branches below applies.
  # Recorded now because `higher_disturbances` is reassigned each loop iteration below
  imply_disturbances <- is.null(higher_disturbances)

  if(!isTRUE(imply_disturbances)){

    # Expand higher-order disturbances to a full covariance matrix
    if(!is(higher_disturbances, "matrix")){

      if(length(higher_disturbances) == 1){
        higher_disturbances <- diag(rep(higher_disturbances, lower_factors), nrow = lower_factors)
      }else{
        higher_disturbances <- diag(higher_disturbances, nrow = lower_factors)
      }

    }

    # Ensure higher-order disturbances form a valid (positive semi-definite) covariance matrix
    if(any(matrix_eigenvalues(higher_disturbances) < 0)){
      stop("`higher_disturbances` must be positive semi-definite (a valid covariance matrix)")
    }

  }else{

    # Expand target lower-order correlations to a full matrix. `higher_disturbances` is
    # implied (solved for) inside the loop below as the residual needed to reproduce this
    # target correlation matrix given the (possibly randomly drawn) higher-order structure,
    # and may end up with off-diagonal values (correlated disturbances)
    if(!is(lower_correlations, "matrix")){
      # A single value is recycled to fill the whole matrix;
      # a flattened vector of length `lower_factors^2` is reshaped column-wise
      lower_correlation_target <- matrix(data = lower_correlations, nrow = lower_factors, ncol = lower_factors)
    }else{
      lower_correlation_target <- lower_correlations
    }

    # Ensure diagonal of target lower-order correlation matrix is 1
    diag(lower_correlation_target) <- 1

  }

  # Initialize checks
  check_eigenvalues <- TRUE
  check_communalities <- TRUE

  # Initialize iteration counter (guards against an infinite loop when the supplied
  # structure leaves no randomness to escape an inadmissible solution, e.g. all matrices
  # of loadings/correlations/disturbances are fixed rather than drawn from ranges)
  iterations <- 0
  max_iterations <- 1000

  # Run through loop
  while(isTRUE(check_eigenvalues) | isTRUE(check_communalities)){

    # Increment iteration counter and check against maximum
    iterations <- iterations + 1
    if(iterations > max_iterations){
      stop(
        paste0(
          "Could not find an admissible (positive semi-definite) solution after ",
          max_iterations, " attempts. The supplied `higher_loadings`, `higher_correlations`, ",
          "and `higher_disturbances`/`lower_correlations` may be jointly inadmissible ",
          "(e.g. implying negative unique variances). Try relaxing fixed matrices into ",
          "ranges or adjusting values"
        )
      )
    }

    # Generate higher-order loadings matrix (lower-order factors onto higher-order factors)
    if(!is(higher_loadings, "matrix")){

      # Initialize higher-order loading matrix
      higher_loading_matrix <- matrix(data = 0, nrow = lower_factors, ncol = higher_factors)

      # Populate higher-order factor loadings
      for(i in higher_factor_sequence){

        # Populate dominant loadings
        if(is.null(higher_loadings_range)){

          # Generate loadings from uniform distribution
          higher_loading_matrix[start_higher[i]:end_higher[i], i] <-
          runif_xoshiro(higher_variables[i], min = higher_loadings[i] - 0.10, max = higher_loadings[i] + 0.10)

        }else{

          # Accept loadings from range (generated earlier)
          higher_loading_matrix[start_higher[i]:end_higher[i], i] <-
          higher_loadings[start_higher[i]:end_higher[i]]

        }

        # Add cross-loadings (if more than one higher-order factor)
        if(higher_factors > 1){

          # Loop through higher-order factors
          for(j in higher_factor_sequence){

            # Ignore dominant factor
            if(i != j){

              # Check for range of cross-loadings
              if(!is.null(higher_cross_loadings_range)){
                higher_loading_matrix[start_higher[j]:end_higher[j], i] <-
                runif_xoshiro(
                  higher_variables[j], min = min(higher_cross_loadings_range),
                  max = max(higher_cross_loadings_range)
                )
              }else{
                higher_loading_matrix[start_higher[j]:end_higher[j], i] <-
                rnorm_ziggurat(higher_variables[j]) * higher_cross_loadings[j]
              }

            }

          }

        }

      }

    }else{# Input is already loadings matrix
      higher_loading_matrix <- higher_loadings
    }

    # Higher-order factor correlations
    if(length(higher_correlations) == 1){

      # Generate correlation matrix
      higher_correlation_matrix <- matrix(data = higher_correlations, nrow = higher_factors, ncol = higher_factors)

    }else{# Input is already correlation matrix
      higher_correlation_matrix <- higher_correlations
    }

    # Ensure diagonal of higher-order correlation matrix is 1
    diag(higher_correlation_matrix) <- 1

    # Assume always positive correlations (for now)
    higher_correlation_matrix <- abs(higher_correlation_matrix)

    # Create lower-order factor covariance implied by higher-order structure
    implied_lower_covariance <- higher_loading_matrix %*% tcrossprod(higher_correlation_matrix, higher_loading_matrix)

    # Check higher-order communalities (variance of lower-order factors explained by higher-order factors)
    check_higher_communalities <- any(diag(implied_lower_covariance) > 0.90)

    # When a target lower-order correlation matrix was supplied, imply the higher-order
    # disturbances needed to reproduce it given this iteration's (possibly randomly drawn)
    # higher-order loadings and correlations; may contain off-diagonal (correlated) values
    if(isTRUE(imply_disturbances)){
      higher_disturbances <- lower_correlation_target - implied_lower_covariance
    }

    # Ensure higher-order disturbances form a valid (positive semi-definite) covariance matrix
    check_disturbances <- any(matrix_eigenvalues(higher_disturbances) < 0)

    # Add disturbances (residual/unique covariance not explained by higher-order factors)
    lower_correlation <- implied_lower_covariance + higher_disturbances

    # Ensure diagonal of lower-order correlation matrix is 1
    diag(lower_correlation) <- 1

    # Generate loadings matrix (variables onto lower-order factors)
    if(!is(lower_loadings, "matrix")){

      # Initialize loading matrix
      loading_matrix <- matrix(data = 0, nrow = total_variables, ncol = lower_factors)

      # Populate factor loadings
      for(i in lower_factor_sequence){

        # Populate dominant loadings
        if(is.null(lower_loadings_range)){

          # Generate loadings from uniform distribution
          loading_matrix[start_variables[i]:end_variables[i], i] <-
          runif_xoshiro(variables[i], min = lower_loadings[i] - 0.10, max = lower_loadings[i] + 0.10)

        }else{

          # Accept loadings from range (generated earlier)
          loading_matrix[start_variables[i]:end_variables[i], i] <-
          lower_loadings[start_variables[i]:end_variables[i]]

        }

        # Add cross-loadings (if more than one lower-order factor)
        if(lower_factors > 1){

          # Loop through lower-order factors
          for(j in lower_factor_sequence){

            # Ignore dominant factor
            if(i != j){

              # Check for range of cross-loadings
              if(!is.null(lower_cross_loadings_range)){
                loading_matrix[start_variables[j]:end_variables[j], i] <-
                runif_xoshiro(
                  variables[j], min = min(lower_cross_loadings_range),
                  max = max(lower_cross_loadings_range)
                )
              }else{
                loading_matrix[start_variables[j]:end_variables[j], i] <-
                rnorm_ziggurat(variables[j]) * lower_cross_loadings[j]
              }

            }

          }

        }

      }

    }else{# Input is already loadings matrix
      loading_matrix <- lower_loadings
    }

    # Create population correlation matrix (variable-level)
    population_correlation <- loading_matrix %*% tcrossprod(lower_correlation, loading_matrix)

    # Check communalities (either level failing is enough to trigger a retry)
    check_communalities <- isTRUE(check_higher_communalities) | any(diag(population_correlation) > 0.90)

    # Ensure diagonal of population correlation matrix is 1
    diag(population_correlation) <- 1

    # Check eigenvalues (disturbance, factor-level, and variable-level matrices)
    check_eigenvalues <- isTRUE(check_disturbances) |
      any(matrix_eigenvalues(lower_correlation) <= 0) |
      any(matrix_eigenvalues(population_correlation) <= 0)

  }

  # Cholesky decomposition
  cholesky <- chol(population_correlation)

  # Generate data
  data <- mvtnorm::rmvnorm(sample_size, sigma = diag(total_variables))

  # Make data based on factor structure
  data <- data %*% cholesky

  # Set data dimensions
  dimensions <- dim(data)

  # Store continuous data
  continuous_data <- data

  # Check for categories greater than categorical limit and not infinite
  variable_categories <- mark_continuous_categories(variable_categories, categorical_limit)

  # Set skew/categories
  categorize_columns <- which(variable_categories <= categorical_limit)
  continuous_columns <- setdiff(seq_len(dimensions[2]), categorize_columns)
  n_continuous <- length(continuous_columns)
  n_categorize <- length(categorize_columns)

  # Store skew (bug fix for later use of `skew`)
  skew_stored <- skew

  # Initialize final skew
  final_skew <- numeric(dimensions[2])

  ## Check for categories
  if(length(categorize_columns) != 0){

    # Set skew
    if(length(skew_stored) == 1){
      skew <- rep(skew_stored, n_categorize)
    }else if(length(skew_stored) != dimensions[2]){
      skew <- shuffle_replace(skew_stored, n_categorize)
    }

    # Check for categories greater than 6
    six_categories <- variable_categories[categorize_columns] > 6
    if(any(six_categories)){

      # Check for skew not equal to zero
      if(any(skew[six_categories] != 0)){

        # Set them equal to zero (overwrite all skews)
        skew[six_categories] <- 0

        # Send warning
        warning("Variables with categories > 6 are not available to add skew. Skew for these variables were set to zero.")

      }

    }

    # Loop through columns
    for(i in seq_along(categorize_columns)){

      data[,categorize_columns[i]] <- categorize(
        data = data[,categorize_columns[i]],
        categories = variable_categories[categorize_columns[i]],
        skew_value = skew[i]
      )

    }

    # Add to final skew
    final_skew[categorize_columns] <- skew
  }

  ## Check for continuous
  if(n_continuous != 0){

    # Set skew
    if(length(skew_stored) == 1){
      skew <- rep(skew_stored, n_continuous)
    }else if(length(skew_stored) != dimensions[2]){
      skew <- shuffle_replace(skew_stored, n_continuous)
    }

    # Loop through columns
    for(i in seq_len(n_continuous)){

      data[,continuous_columns[i]] <- skew_continuous( # function in `utils-latentFactoR`
        skew = skew[i], data = data[,continuous_columns[i]]
      )

    }

    # Add to final skew
    final_skew[continuous_columns] <- skew

  }

  # Add column names to data
  colnames(data) <- paste0(
    "V", format_integer(seq_len(total_variables), digits(total_variables) - 1)
  )

  # Populate parameters
  parameters <- list(
    lower = list(
      factors = lower_factors,
      variables = variables,
      loadings = loading_matrix,
      cross_loadings = lower_cross_loadings,
      correlations = lower_correlation
    ),
    higher = list(
      factors = higher_factors,
      variables = higher_variables,
      loadings = higher_loading_matrix,
      cross_loadings = higher_cross_loadings,
      correlations = higher_correlation_matrix,
      disturbances = higher_disturbances
    ),
    categories = variable_categories,
    categorical_limit = categorical_limit,
    skew = final_skew
  )

  # Populate results
  results <- list(
    data = data,
    population_correlation = population_correlation,
    lower_correlation = lower_correlation,
    continuous_data = continuous_data,
    parameters = parameters
  )

  # Add class
  class(results) <- "lf_simulate"

  # Return results
  return(results)

}

# Bug checking
# lower_factors = 4; variables = 6; variables_range = NULL
# lower_loadings = 0.55; lower_loadings_range = NULL
# lower_cross_loadings = 0.05; lower_cross_loadings_range = NULL
# higher_factors = 1; higher_variables = NULL; higher_loadings = 0.55; higher_loadings_range = NULL
# higher_cross_loadings = 0.05; higher_cross_loadings_range = NULL
# higher_correlations = 0.30; higher_correlations_range = NULL
# higher_disturbances = 0.10; higher_disturbances_range = c(0, 0.20)
# sample_size = 1000; variable_categories = Inf
# categorical_limit = 7; skew = 0; skew_range = NULL

# Input checking ----
#' @noRd
# Updated 12.07.2026
simulate_hierarchical_factors_errors <- function(
    lower_factors, variables, variables_range = NULL,
    lower_loadings, lower_loadings_range = NULL,
    lower_cross_loadings, lower_cross_loadings_range = NULL,
    lower_correlations = NULL, lower_correlations_range = NULL,
    higher_factors, higher_variables = NULL,
    higher_loadings, higher_loadings_range = NULL,
    higher_cross_loadings, higher_cross_loadings_range = NULL,
    higher_correlations, higher_correlations_range = NULL,
    higher_disturbances = NULL, higher_disturbances_range = NULL,
    sample_size, variable_categories = Inf,
    categorical_limit = 7,
    skew = 0, skew_range = NULL
)
{

  # Check lower-order factors first
  type_error(lower_factors, "numeric")
  length_error(lower_factors, 1)
  range_error(lower_factors, c(1, Inf))

  # Check higher-order factors second
  type_error(higher_factors, "numeric")
  length_error(higher_factors, 1)
  range_error(higher_factors, c(1, Inf))

  # Check for number of lower-order factors nested within each higher-order factor
  if(!is.null(higher_variables)){
    type_error(higher_variables, "numeric") # object type error
    length_error(higher_variables, c(1, higher_factors)) # object length error
    range_error(higher_variables, c(1, Inf)) # object range error
  }

  # Check sample size third
  type_error(sample_size, "numeric")
  length_error(sample_size, 1)
  range_error(sample_size, c(1, Inf))

  # Disturbances and a target lower-order correlation matrix are two alternate ways
  # to arrive at the same lower-order correlation structure; exactly one is required
  disturbances_supplied <- !is.null(higher_disturbances) || !is.null(higher_disturbances_range)
  correlations_supplied <- !is.null(lower_correlations) || !is.null(lower_correlations_range)

  if(disturbances_supplied && correlations_supplied){
    stop(
      "Only one of `higher_disturbances`/`higher_disturbances_range` or ",
      "`lower_correlations`/`lower_correlations_range` can be supplied"
    )
  }else if(!disturbances_supplied && !correlations_supplied){
    stop(
      "One of `higher_disturbances`/`higher_disturbances_range` or ",
      "`lower_correlations`/`lower_correlations_range` must be supplied"
    )
  }

  # Check for variables range
  if(!is.null(variables_range)){
    type_error(variables_range, "numeric") # object type error
    length_error(variables_range, 2) # object length error
    range_error(variables_range, c(3, Inf)) # object range error
    variables <- round(runif(
      lower_factors,
      min = min(variables_range),
      max = max(variables_range)
    ))
  }

  # Determine total number of variables
  if(length(variables) == 1){
    total_variables <- lower_factors * variables
  }else if(length(variables) == lower_factors){
    total_variables <- sum(variables)
  }

  # Check for lower-order loadings range
  if(!is.null(lower_loadings_range)){
    type_error(lower_loadings_range, "numeric") # object type error
    length_error(lower_loadings_range, 2)  # object length error
    range_error(lower_loadings_range, c(-1, 1)) # object range error
    lower_loadings <- runif(
      total_variables,
      min = min(lower_loadings_range),
      max = max(lower_loadings_range)
    )
  }

  # Check for lower-order cross-loadings range
  if(!is.null(lower_cross_loadings_range)){
    type_error(lower_cross_loadings_range, "numeric") # object type error
    length_error(lower_cross_loadings_range, 2) # object length error
    range_error(lower_cross_loadings_range, c(-1, 1)) # object range error
    lower_cross_loadings <- runif(
      lower_factors,
      min = min(lower_cross_loadings_range),
      max = max(lower_cross_loadings_range)
    )
  }

  # Check for lower-order correlations range (target lower-order correlation matrix,
  # used to imply `higher_disturbances` rather than supplying it directly)
  if(!is.null(lower_correlations_range)){
    type_error(lower_correlations_range, "numeric") # object type error
    length_error(lower_correlations_range, 2) # object length error
    range_error(lower_correlations_range, c(-1, 1)) # object range error

    # Initialize target lower-order correlation matrix
    lower_correlation_matrix <- matrix(data = 0, nrow = lower_factors, ncol = lower_factors)

    # Populate lower-order correlation matrix
    lower_correlation_matrix[
      lower.tri(lower_correlation_matrix)
    ] <- runif(
      sum(lower.tri(lower_correlation_matrix)),
      min = min(lower_correlations_range),
      max = max(lower_correlations_range)
    )

    # Make correlation matrix symmetric
    lower_correlations <- lower_correlation_matrix + t(lower_correlation_matrix)

  }

  # Check for higher-order loadings range
  if(!is.null(higher_loadings_range)){
    type_error(higher_loadings_range, "numeric") # object type error
    length_error(higher_loadings_range, 2)  # object length error
    range_error(higher_loadings_range, c(-1, 1)) # object range error
    higher_loadings <- runif(
      lower_factors,
      min = min(higher_loadings_range),
      max = max(higher_loadings_range)
    )
  }

  # Check for higher-order cross-loadings range
  if(!is.null(higher_cross_loadings_range)){
    type_error(higher_cross_loadings_range, "numeric") # object type error
    length_error(higher_cross_loadings_range, 2) # object length error
    range_error(higher_cross_loadings_range, c(-1, 1)) # object range error
    higher_cross_loadings <- runif(
      higher_factors,
      min = min(higher_cross_loadings_range),
      max = max(higher_cross_loadings_range)
    )
  }

  # Check for higher-order correlations range
  if(!is.null(higher_correlations_range)){
    type_error(higher_correlations_range, "numeric") # object type error
    length_error(higher_correlations_range, 2) # object length error
    range_error(higher_correlations_range, c(-1, 1)) # object range error

    # Initialize correlation matrix
    higher_correlation_matrix <- matrix(data = 0, nrow = higher_factors, ncol = higher_factors)

    # Population correlation matrix
    higher_correlation_matrix[
      lower.tri(higher_correlation_matrix)
    ] <- runif(
      sum(lower.tri(higher_correlation_matrix)),
      min = min(higher_correlations_range),
      max = max(higher_correlations_range)
    )

    # Make correlation matrix symmetric
    higher_correlations <- higher_correlation_matrix + t(higher_correlation_matrix)

  }

  # Check for higher-order disturbances range
  if(!is.null(higher_disturbances_range)){
    type_error(higher_disturbances_range, "numeric") # object type error
    length_error(higher_disturbances_range, 2) # object length error
    range_error(higher_disturbances_range, c(0, 1)) # object range error

    # Independently draw a disturbance variance for each lower-order factor
    higher_disturbances <- diag(
      runif(
        lower_factors,
        min = min(higher_disturbances_range),
        max = max(higher_disturbances_range)
      ),
      nrow = lower_factors
    )
  }

  # Check for skew range
  if(!is.null(skew_range)){
    type_error(skew_range, "numeric") # object type error
    length_error(skew_range, 2) # object length error
    range_error(skew_range, c(-2, 2)) # object range error
    possible_skews <- seq(-2, 2, 0.05) # possible skews
    skew_range <- round(skew_range, 2) # get to hundredths digit
    min_range <- abs(min(skew_range) - possible_skews) # difference for minimum
    min_skew <- possible_skews[which.min(min_range)] # get minimum skew
    max_range <- abs(max(skew_range) - possible_skews) # difference for maximum
    max_skew <- possible_skews[which.min(max_range)] # get maximum skew
    skew <- seq(min_skew, max_skew, 0.05) # obtain skews
  }

  # Ensure appropriate type and length for categories
  type_error(variable_categories, "numeric")
  length_error(variable_categories, c(1, total_variables))

  # Ensure appropriate types
  type_error(variables, "numeric"); type_error(variable_categories, "numeric")
  type_error(categorical_limit, "numeric"); type_error(skew, "numeric")

  # Ensure appropriate lengths
  length_error(variables, c(1, lower_factors)); length_error(categorical_limit, 1)

  # Ensure appropriate ranges
  range_error(variables, c(2, Inf)); range_error(skew, c(-2, 2))
  range_error(variable_categories, c(2, Inf))

  # Return checked input
  return(
    list(
      variables = variables, total_variables = total_variables,
      lower_loadings = lower_loadings, lower_cross_loadings = lower_cross_loadings,
      lower_correlations = lower_correlations,
      higher_variables = higher_variables,
      higher_loadings = higher_loadings, higher_cross_loadings = higher_cross_loadings,
      higher_correlations = higher_correlations, higher_disturbances = higher_disturbances,
      skew = skew
    )
  )

}
