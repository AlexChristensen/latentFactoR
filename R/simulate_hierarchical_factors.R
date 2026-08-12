#' Simulates Hierarchical Latent Factor Data
#'
#' @description
#' Simulates data from a hierarchical (second-order) latent factor model based
#' on many manipulable parameters. Variables load onto lower-order factors, and
#' lower-order factors load onto higher-order factors. Parameters do not have
#' default values and must each be set, with the exception of \code{off_disturbances}
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
#' is the disturbance (residual) correlation matrix of the lower-order factors, with its \emph{off-diagonal}
#' (correlated disturbance) elements set directly via \code{off_disturbances} — its diagonal is not
#' controllable: whatever variance \eqn{\boldsymbol{\Phi}_H}{Phi_H} and \eqn{\boldsymbol{\Lambda}_H}{Lambda_H}
#' leave unexplained in each lower-order factor is always rescaled back to unit variance, so only the
#' correlation between disturbances (not their magnitude) can be set
#'
#' @details
#' How \eqn{\boldsymbol{\Phi}_H}{Phi_H} and \eqn{\boldsymbol{\Psi}}{Psi} are obtained depends on which
#' parameters are supplied. Exactly one of \code{off_disturbances}/\code{off_disturbances_range} or
#' \code{lower_correlations}/\code{lower_correlations_range} must be supplied;
#' \code{higher_correlations}/\code{higher_correlations_range} is required alongside the former, but
#' optional alongside the latter:
#'
#' \tabular{ll}{
#'   \strong{Supplied} \tab \strong{Implied} \cr
#'   \code{higher_correlations} + \code{off_disturbances} \tab \code{lower_correlations} (i.e., \eqn{\boldsymbol{\Phi}_L}{Phi_L}) \cr
#'   \code{higher_correlations} + \code{lower_correlations} \tab \code{off_disturbances} (i.e., \eqn{\boldsymbol{\Psi}}{Psi}) \cr
#'   \code{lower_correlations} (alone) \tab \code{higher_correlations} and \code{off_disturbances} (i.e., \eqn{\boldsymbol{\Phi}_H}{Phi_H} and \eqn{\boldsymbol{\Psi}}{Psi}) \cr
#' }
#'
#' When \eqn{\boldsymbol{\Psi}}{Psi} is implied, it is solved for as the residual needed to reproduce the
#' target \code{lower_correlations} matrix given \eqn{\boldsymbol{\Phi}_H}{Phi_H} and
#' \eqn{\boldsymbol{\Lambda}_H}{Lambda_H}, and may end up with a non-unit diagonal and/or off-diagonal
#' values (i.e., correlated disturbances between lower-order factors are permitted); unlike a directly
#' supplied \code{off_disturbances}, this implied \eqn{\boldsymbol{\Psi}}{Psi} keeps whatever diagonal
#' falls out of the subtraction, since it is not passed back through the unit-variance rescaling that a
#' supplied \code{off_disturbances} is. When \eqn{\boldsymbol{\Phi}_H}{Phi_H} is also
#' implied (i.e., \code{higher_correlations} is omitted), it is obtained first, from the target
#' \code{lower_correlations} matrix and \eqn{\boldsymbol{\Lambda}_H}{Lambda_H}, via a generalized-inverse
#' (least-squares) regression:
#'
#' \deqn{\boldsymbol{\Phi}_H = (\boldsymbol{\Lambda}_H'\boldsymbol{\Lambda}_H)^{-1}\boldsymbol{\Lambda}_H' \boldsymbol{\Phi}_L \boldsymbol{\Lambda}_H(\boldsymbol{\Lambda}_H'\boldsymbol{\Lambda}_H)^{-1}}{Phi_H = solve(crossprod(Lambda_H)) \%*\% t(Lambda_H) \%*\% Phi_L \%*\% Lambda_H \%*\% solve(crossprod(Lambda_H))}
#'
#' This raw least-squares result is not guaranteed to be a valid correlation matrix (unit diagonal,
#' \code{[-1, 1]} off-diagonals), so it is standardized with \code{\link[stats]{cov2cor}}, after which its
#' diagonal is reset to 1 and off-diagonal elements are taken in absolute value (consistent with how a
#' directly supplied \code{higher_correlations} is handled); \eqn{\boldsymbol{\Psi}}{Psi} is then solved
#' for as described above
#'
#' In every case, \eqn{\boldsymbol{\Psi}}{Psi} (whether fixed or implied) must form a valid (positive
#' semi-definite) covariance matrix; if a target correlation matrix or fixed higher-order structure is
#' incompatible (e.g., it implies negative unique variances), values are redrawn when
#' \code{higher_loadings_range}, \code{higher_correlations_range}, and/or \code{higher_cross_loadings_range}
#' are used, or an error is raised when the higher-order structure is entirely fixed
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
#' \code{off_disturbances}. Can be a single value (applied to all off-diagonal correlations) or a full
#' \code{lower_factors} x \code{lower_factors} correlation matrix. When supplied, \code{off_disturbances}
#' (and, if \code{higher_correlations} is also omitted, \code{higher_correlations} too) is \emph{implied}
#' rather than set directly (see Details).
#' Mutually exclusive with \code{off_disturbances}/\code{off_disturbances_range}: exactly one of
#' \code{off_disturbances}/\code{off_disturbances_range} or
#' \code{lower_correlations}/\code{lower_correlations_range} must be supplied
#'
#' @param lower_correlations_range Numeric (length = 2).
#' Range of lower-order correlations to randomly select (per pair of lower-order factors) from a random
#' uniform distribution. Somewhat redundant with \code{lower_correlations} but more flexible.
#' Mutually exclusive with \code{off_disturbances}/\code{off_disturbances_range}
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
#' General effect sizes range from orthogonal (0.00), small (0.30), moderate (0.50), to large (0.70).
#' Required when \code{off_disturbances}/\code{off_disturbances_range} is used. Optional when
#' \code{lower_correlations}/\code{lower_correlations_range} is used instead, in which case omitting it
#' (and \code{higher_correlations_range}) implies it instead of supplying it directly (see Details)
#'
#' @param higher_correlations_range Numeric (length = 2).
#' Range of correlations between higher-order factors to randomly select from a random uniform distribution.
#' General effect sizes range from orthogonal (0.00), small (0.30), moderate (0.50), to large (0.70).
#' Required when \code{off_disturbances}/\code{off_disturbances_range} is used, unless
#' \code{higher_correlations} is supplied instead. Optional when
#' \code{lower_correlations}/\code{lower_correlations_range} is used (see \code{higher_correlations})
#'
#' @param off_disturbances Numeric or matrix (length = 1 or \code{lower_factors} x \code{lower_factors}).
#' \emph{Off-diagonal} (correlated) disturbances between lower-order factors' residuals, left unexplained
#' by the higher-order factor(s); the diagonal is not settable here (see Details for why). A single value
#' is applied uniformly (no jitter) to every off-diagonal element; a full
#' \code{lower_factors} x \code{lower_factors} matrix is used exactly as supplied instead (its diagonal is
#' ignored and reset to 1). Values must be between -1 and 1 and the resulting matrix (with unit diagonal)
#' must be positive semi-definite.
#' Mutually exclusive with \code{lower_correlations}/\code{lower_correlations_range}: exactly one of
#' \code{off_disturbances}/\code{off_disturbances_range} or
#' \code{lower_correlations}/\code{lower_correlations_range} must be supplied
#'
#' @param off_disturbances_range Numeric (length = 2).
#' Range of off-diagonal (correlated) disturbances to randomly and independently select (per pair of
#' lower-order factors) from a random uniform distribution.
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
#' Current skew implementation is between -10 and 10 in increments of 0.05.
#' Skews that are not in this sequence will be converted to their nearest
#' value in the sequence. Not recommended to use with \code{variables_range}
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
#' \code{disturbances} (disturbance (residual) correlation matrix of lower-order factors;
#' when \code{off_disturbances} is supplied directly, only its off-diagonal elements are
#' meaningful, since its diagonal is always reset to 1; when implied instead, via
#' \code{lower_correlations}, the diagonal reflects the residual variance actually solved for)
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
#'   off_disturbances = 0.10, # off-diagonal (correlated) disturbances = 0.10
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
#'   higher_correlations = 0.30, off_disturbances = 0.10,
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
#'   off_disturbances_range = c(0.00, 0.20), # off-diagonal disturbances = 0.00 to 0.20
#'   sample_size = 1000
#' )
#'
#' # Generate dichotomous data
#' hierarchical_dichotomous <- simulate_hierarchical_factors(
#'   lower_factors = 4, variables = 6,
#'   lower_loadings = 0.55, lower_cross_loadings = 0.05,
#'   higher_factors = 2, higher_loadings = 0.60,
#'   higher_cross_loadings = 0.05, higher_correlations = 0.30,
#'   off_disturbances = 0.10, sample_size = 1000,
#'   variable_categories = 2 # dichotomous data
#' )
#'
#' # Supply a target lower-order correlation matrix directly instead of `off_disturbances`;
#' # the disturbances are implied (solved for) and may end up correlated
#' hierarchical_lower_correlations <- simulate_hierarchical_factors(
#'   lower_factors = 4, variables = 6,
#'   lower_loadings = 0.55, lower_cross_loadings = 0.05,
#'   lower_correlations = 0.40, # target correlation of 0.40 between all lower-order factors
#'   higher_factors = 2, higher_loadings = 0.60,
#'   higher_cross_loadings = 0.05, higher_correlations = 0.30,
#'   sample_size = 1000
#' )
#'
#' # Omit `higher_correlations` when `lower_correlations` is supplied: it is implied
#' # (rather than the disturbances) from the target lower-order correlations and `higher_loadings`
#' hierarchical_implied_correlations <- simulate_hierarchical_factors(
#'   lower_factors = 4, variables = 6,
#'   lower_loadings = 0.55, lower_cross_loadings = 0.05,
#'   lower_correlations = 0.40, # target correlation of 0.40 between all lower-order factors
#'   higher_factors = 2, higher_loadings = 0.60,
#'   higher_cross_loadings = 0.05,
#'   sample_size = 1000
#' )
#'
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @export
#'
# Main hierarchical factor simulation function
# Updated 13.07.2026
simulate_hierarchical_factors <- function(
    lower_factors, variables, variables_range = NULL,
    lower_loadings, lower_loadings_range = NULL,
    lower_cross_loadings, lower_cross_loadings_range = NULL,
    lower_correlations = NULL, lower_correlations_range = NULL,
    higher_factors, higher_variables = NULL,
    higher_loadings, higher_loadings_range = NULL,
    higher_cross_loadings, higher_cross_loadings_range = NULL,
    higher_correlations = NULL, higher_correlations_range = NULL,
    off_disturbances = NULL, off_disturbances_range = NULL,
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
    off_disturbances, off_disturbances_range,
    sample_size, variable_categories,
    categorical_limit, skew, skew_range
  )

  # Collect inputs
  variables <- inputs$variables; total_variables <- inputs$total_variables
  lower_loadings <- inputs$lower_loadings; lower_cross_loadings <- inputs$lower_cross_loadings
  lower_correlations <- inputs$lower_correlations
  higher_variables <- inputs$higher_variables
  higher_loadings <- inputs$higher_loadings; higher_cross_loadings <- inputs$higher_cross_loadings
  higher_correlations <- inputs$higher_correlations; off_disturbances <- inputs$off_disturbances
  skew <- inputs$skew

  # Ensure appropriate lengths
  if(!is.null(lower_correlations)){
    length_error(lower_correlations, c(1, lower_factors * lower_factors))
  }
  if(!is.null(higher_correlations)){
    length_error(higher_correlations, c(1, higher_factors * higher_factors))
  }
  if(!is.null(off_disturbances)){
    length_error(off_disturbances, c(1, lower_factors * lower_factors))
  }

  # Ensure appropriate ranges
  range_error(lower_loadings, c(-1, 1)); range_error(lower_cross_loadings, c(-1, 1))
  if(!is.null(lower_correlations)){range_error(lower_correlations, c(-1, 1))}
  range_error(higher_loadings, c(-1, 1)); range_error(higher_cross_loadings, c(-1, 1))
  if(!is.null(higher_correlations)){range_error(higher_correlations, c(-1, 1))}
  if(!is.null(off_disturbances)){range_error(off_disturbances, c(-1, 1))}

  # Ensure appropriate types
  if(!is(lower_loadings, "matrix")){type_error(lower_loadings, "numeric")}
  if(!is(lower_cross_loadings, "matrix")){type_error(lower_cross_loadings, "numeric")}
  if(!is.null(lower_correlations) && !is(lower_correlations, "matrix")){
    type_error(lower_correlations, "numeric")
  }
  if(!is(higher_loadings, "matrix")){type_error(higher_loadings, "numeric")}
  if(!is(higher_cross_loadings, "matrix")){type_error(higher_cross_loadings, "numeric")}
  if(!is.null(higher_correlations) && !is(higher_correlations, "matrix")){
    type_error(higher_correlations, "numeric")
  }
  if(!is.null(off_disturbances) && !is(off_disturbances, "matrix")){
    type_error(off_disturbances, "numeric")
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

  # `off_disturbances` and `lower_correlations` are mutually exclusive (enforced in
  # `simulate_hierarchical_factors_errors`); exactly one of the two branches below applies.
  imply_disturbances <- is.null(off_disturbances)

  # When `lower_correlations` is supplied, `higher_correlations` is optional: if the user
  # also omitted it, it is implied (each loop iteration) from the target lower-order
  # correlations and that iteration's higher-order loadings, rather than being fixed
  imply_higher_correlations <- imply_disturbances && is.null(higher_correlations)

  if(!imply_disturbances){

    # A single value is applied to every off-diagonal element (uniformly, no jitter); a full
    # matrix is used as-is instead. Either way, only the off-diagonal elements matter
    # (correlated disturbances), so the diagonal is reset to 1. This is fixed (not redrawn
    # each loop iteration below), so its validity is checked once, up front
    if(!is(off_disturbances, "matrix")){
      off_disturbances <- matrix(data = off_disturbances, nrow = lower_factors, ncol = lower_factors)
    }

    diag(off_disturbances) <- 1

    if(any(matrix_eigenvalues(off_disturbances) < 0)){
      stop("`off_disturbances` must be positive semi-definite (a valid correlation matrix)")
    }

  }

  if(imply_disturbances){

    # Expand target lower-order correlations to a full matrix. `off_disturbances` is
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
  while(check_eigenvalues | check_communalities){

    # Increment iteration counter and check against maximum
    iterations <- iterations + 1
    if(iterations > max_iterations){
      stop(
        paste0(
          "Could not find an admissible (positive semi-definite) solution after ",
          max_iterations, " attempts. The supplied `higher_loadings`, `higher_correlations`, ",
          "and `off_disturbances`/`lower_correlations` may be jointly inadmissible ",
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

          # Draw fresh loadings from the range on every retry
          # attempt (rather than reusing a single fixed draw),
          # so retries can actually escape an inadmissible solution
          higher_loading_matrix[start_higher[i]:end_higher[i], i] <-
            runif_xoshiro(
              higher_variables[i], min = min(higher_loadings_range), max = max(higher_loadings_range)
            )

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
    if(imply_higher_correlations){

      # Imply higher-order correlations from the target lower-order correlation matrix and
      # this iteration's higher-order loadings via a generalized-inverse (least-squares)
      # regression: Phi_H = pinv(Lambda_H) * Phi_L * pinv(Lambda_H)'. The raw result is not
      # guaranteed to be a valid correlation matrix (unit diagonal, [-1, 1] off-diagonals),
      # so it is standardized into one with `cov2cor`
      higher_inverse <- solve(crossprod(higher_loading_matrix))
      higher_correlation_matrix <- cov2cor(
        tcrossprod(higher_inverse, higher_loading_matrix) %*%
          lower_correlation_target %*% higher_loading_matrix %*% higher_inverse
      )

    }else if(!is.null(higher_correlations_range)){

      # Draw fresh correlations from the range on every retry
      # attempt (rather than reusing a single fixed draw),
      # so retries can actually escape an inadmissible solution
      higher_correlation_matrix <- matrix(data = 0, nrow = higher_factors, ncol = higher_factors)
      higher_correlation_matrix[lower.tri(higher_correlation_matrix)] <- runif_xoshiro(
        sum(lower.tri(higher_correlation_matrix)),
        min = min(higher_correlations_range), max = max(higher_correlations_range)
      )
      higher_correlation_matrix <- higher_correlation_matrix + t(higher_correlation_matrix)

    }else if(length(higher_correlations) == 1){

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

    # Off-diagonal (correlated) disturbances between lower-order factors' residuals. When a
    # target lower-order correlation matrix was supplied instead, the disturbances are fully
    # implied as the residual needed to reproduce it given this iteration's (possibly randomly
    # drawn) higher-order loadings and correlations, and may contain a non-unit diagonal
    # and/or off-diagonal (correlated) values. Otherwise, `off_disturbances` is already a fixed,
    # validated matrix (expanded up front) and is used as-is
    if(isTRUE(imply_disturbances)){

      disturbance_matrix <- lower_correlation_target - implied_lower_covariance

    }else{

      disturbance_matrix <- off_disturbances

    }

    # Ensure disturbances form a valid (positive semi-definite) correlation matrix
    check_disturbances <- any(matrix_eigenvalues(disturbance_matrix) < 0)

    # Add disturbances (residual/correlated covariance not explained by higher-order factors)
    lower_correlation <- implied_lower_covariance + disturbance_matrix

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

          # Draw fresh loadings from the range on every retry
          # attempt (rather than reusing a single fixed draw),
          # so retries can actually escape an inadmissible solution
          loading_matrix[start_variables[i]:end_variables[i], i] <-
            runif_xoshiro(
              variables[i], min = min(lower_loadings_range), max = max(lower_loadings_range)
            )

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
    check_communalities <- check_higher_communalities | any(diag(population_correlation) > 0.90)

    # Ensure diagonal of population correlation matrix is 1
    diag(population_correlation) <- 1

    # Check eigenvalues (disturbance, factor-level, and variable-level matrices)
    check_eigenvalues <- check_disturbances |
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
    }else if(length(skew_stored) == dimensions[2]){
      # Full-length (per-variable) skew: select this subset's own values
      skew <- skew_stored[categorize_columns]
    }else{
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
    }else if(length(skew_stored) == dimensions[2]){
      # Full-length (per-variable) skew: select this subset's own values
      skew <- skew_stored[continuous_columns]
    }else{
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
      disturbances = disturbance_matrix
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
# lower_correlations = 0.30; lower_correlations_range = NULL
# higher_factors = 2; higher_variables = 2;
# higher_loadings = 0.55; higher_loadings_range = NULL
# higher_cross_loadings = 0.05; higher_cross_loadings_range = NULL
# higher_correlations = NULL; higher_correlations_range = NULL
# off_disturbances = NULL; off_disturbances_range = NULL
# sample_size = 1000; variable_categories = Inf
# categorical_limit = 7; skew = 0; skew_range = NULL

# Input checking ----
#' @noRd
# Updated 13.07.2026
simulate_hierarchical_factors_errors <- function(
    lower_factors, variables, variables_range = NULL,
    lower_loadings, lower_loadings_range = NULL,
    lower_cross_loadings, lower_cross_loadings_range = NULL,
    lower_correlations = NULL, lower_correlations_range = NULL,
    higher_factors, higher_variables = NULL,
    higher_loadings, higher_loadings_range = NULL,
    higher_cross_loadings, higher_cross_loadings_range = NULL,
    higher_correlations = NULL, higher_correlations_range = NULL,
    off_disturbances = NULL, off_disturbances_range = NULL,
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
  disturbances_supplied <- !is.null(off_disturbances) || !is.null(off_disturbances_range)
  correlations_supplied <- !is.null(lower_correlations) || !is.null(lower_correlations_range)

  if(disturbances_supplied && correlations_supplied){
    stop(
      "Only one of `off_disturbances`/`off_disturbances_range` or ",
      "`lower_correlations`/`lower_correlations_range` can be supplied"
    )
  }else if(!disturbances_supplied && !correlations_supplied){
    stop(
      "One of `off_disturbances`/`off_disturbances_range` or ",
      "`lower_correlations`/`lower_correlations_range` must be supplied"
    )
  }

  # `higher_correlations` is only optional (implied) when `lower_correlations` is supplied
  # instead of `off_disturbances`; in the `off_disturbances` branch it must be supplied
  if(
    isTRUE(disturbances_supplied) &&
    is.null(higher_correlations) && is.null(higher_correlations_range)
  ){
    stop(
      "`higher_correlations`/`higher_correlations_range` must be supplied when using ",
      "`off_disturbances`/`off_disturbances_range` (they are only optional when ",
      "`lower_correlations`/`lower_correlations_range` is supplied instead)"
    )
  }

  # Check for variables range
  if(!is.null(variables_range)){
    type_error(variables_range, "numeric") # object type error
    length_error(variables_range, 2) # object length error
    range_error(variables_range, c(3, Inf)) # object range error
    variables <- round(runif_xoshiro(
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
    lower_loadings <- runif_xoshiro(
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
    lower_cross_loadings <- runif_xoshiro(
      lower_factors,
      min = min(lower_cross_loadings_range),
      max = max(lower_cross_loadings_range)
    )
  }

  # Check for lower-order correlations range (target lower-order correlation matrix,
  # used to imply `off_disturbances` rather than supplying it directly)
  if(!is.null(lower_correlations_range)){
    type_error(lower_correlations_range, "numeric") # object type error
    length_error(lower_correlations_range, 2) # object length error
    range_error(lower_correlations_range, c(-1, 1)) # object range error

    # Initialize target lower-order correlation matrix
    lower_correlation_matrix <- matrix(data = 0, nrow = lower_factors, ncol = lower_factors)

    # Populate lower-order correlation matrix
    lower_correlation_matrix[
      lower.tri(lower_correlation_matrix)
    ] <- runif_xoshiro(
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
    higher_loadings <- runif_xoshiro(
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
    higher_cross_loadings <- runif_xoshiro(
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
    ] <- runif_xoshiro(
      sum(lower.tri(higher_correlation_matrix)),
      min = min(higher_correlations_range),
      max = max(higher_correlations_range)
    )

    # Make correlation matrix symmetric
    higher_correlations <- higher_correlation_matrix + t(higher_correlation_matrix)

  }

  # Check for off-diagonal disturbances range
  if(!is.null(off_disturbances_range)){
    type_error(off_disturbances_range, "numeric") # object type error
    length_error(off_disturbances_range, 2) # object length error
    range_error(off_disturbances_range, c(-1, 1)) # object range error

    # Initialize off-diagonal disturbance matrix
    off_disturbance_matrix <- matrix(data = 0, nrow = lower_factors, ncol = lower_factors)

    # Independently draw an off-diagonal (correlated) disturbance for each pair of lower-order factors
    off_disturbance_matrix[
      lower.tri(off_disturbance_matrix)
    ] <- runif_xoshiro(
      sum(lower.tri(off_disturbance_matrix)),
      min = min(off_disturbances_range),
      max = max(off_disturbances_range)
    )

    # Make matrix symmetric
    off_disturbances <- off_disturbance_matrix + t(off_disturbance_matrix)
  }

  # Check for skew range
  if(!is.null(skew_range)){
    type_error(skew_range, "numeric") # object type error
    length_error(skew_range, 2) # object length error
    range_error(skew_range, c(-10, 10)) # object range error
    possible_skews <- seq(-10, 10, 0.05) # possible skews
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
  range_error(variables, c(2, Inf)); range_error(skew, c(-10, 10))
  range_error(variable_categories, c(2, Inf))

  # Return checked input
  return(
    list(
      variables = variables, total_variables = total_variables,
      lower_loadings = lower_loadings, lower_cross_loadings = lower_cross_loadings,
      lower_correlations = lower_correlations,
      higher_variables = higher_variables,
      higher_loadings = higher_loadings, higher_cross_loadings = higher_cross_loadings,
      higher_correlations = higher_correlations, off_disturbances = off_disturbances,
      skew = skew
    )
  )

}