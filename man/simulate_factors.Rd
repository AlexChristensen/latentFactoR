% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_factors.R
\name{simulate_factors}
\alias{simulate_factors}
\title{Simulates Latent Factor Data}
\usage{
simulate_factors(
  factors,
  variables,
  variables_range = NULL,
  loadings,
  loadings_range = NULL,
  cross_loadings,
  cross_loadings_range = NULL,
  correlations,
  correlations_range = NULL,
  sample_size,
  variable_categories = Inf,
  categorical_limit = 7,
  skew = 0,
  skew_range = NULL
)
}
\arguments{
\item{factors}{Numeric (length = 1).
Number of factors}

\item{variables}{Numeric (length = 1 or \code{factors}).
Number of variables per factor.
Can be a single value or as many values as there are factors.
Minimum three variables per factor}

\item{variables_range}{Numeric (length = 2).
Range of variables to randomly select from a random uniform distribution.
Minimum three variables per factor}

\item{loadings}{Numeric or matrix (length = 1, \code{factors}, total number of variables (\code{factors} x \code{variables}), or \code{factors} x total number of variables.
Loadings drawn from a random uniform distribution using +/- 0.10 of value input.
Can be a single value or as many values as there are factors (corresponding to the factors).
Can also be a loading matrix. Columns must match factors and rows must match total variables (\code{factors} x \code{variables})
General effect sizes range from small (0.40), moderate (0.55), to large (0.70)}

\item{loadings_range}{Numeric (length = 2).
Range of loadings to randomly select from a random uniform distribution.
General effect sizes range from small (0.40), moderate (0.55), to large (0.70)}

\item{cross_loadings}{Numeric or matrix(length = 1, \code{factors}, or \code{factors} x total number of variables.
Cross-loadings drawn from a random normal distribution with a mean of 0 and standard deviation of value input.
Can be a single value or as many values as there are factors (corresponding to the factors).
Can also be a loading matrix. Columns must match factors and rows must match total variables (\code{factors} x \code{variables})}

\item{cross_loadings_range}{Numeric (length = 2).
Range of cross-loadings to randomly select from a random uniform distribution}

\item{correlations}{Numeric (length = 1 or \code{factors} x \code{factors}).
Can be a single value that will be used for all correlations between factors.
Can also be a square matrix (\code{factors} x \code{factors}).
General effect sizes range from orthogonal (0.00), small (0.30), moderate (0.50), to large (0.70)}

\item{correlations_range}{Numeric (length = 2).
Range of correlations to randomly select from a random uniform distribution.
General effect sizes range from orthogonal (0.00), small (0.30), moderate (0.50), to large (0.70)}

\item{sample_size}{Numeric (length = 1).
Number of cases to generate from a random multivariate normal distribution using
\code{\link[mvtnorm]{rmvnorm}}}

\item{variable_categories}{Numeric (length = 1 or total variables (\code{factors} x \code{variables})).
Number of categories for each variable. \code{Inf} is used for continuous variables; otherwise,
values reflect number of categories}

\item{categorical_limit}{Numeric (length = 1).
Values greater than input value are considered continuous.
Defaults to \code{7} meaning that 8 or more categories are considered continuous
(i.e., variables are \emph{not} categorized from continuous to categorical)}

\item{skew}{Numeric (length = 1 or categorical variables).
Skew to be included in categorical variables. It is randomly sampled from provided values.
Can be a single value or as many values as there are (total) variables.
Current skew implementation is between -2 and 2 in increments of 0.05.
Skews that are not in this sequence will be converted to their nearest
value in the sequence. Not recommended to use with \code{variables_range}.
Future versions will incorporate finer skews}

\item{skew_range}{Numeric (length = 2).
Randomly selects skews within in the range.
Somewhat redundant with \code{skew} but more flexible.
Compatible with \code{variables_range}}
}
\value{
Returns a list containing:

\item{data}{Simulated data from the specified factor model}

\item{population_correlation}{Population correlation matrix}

\item{parameters}{
A list containing the parameters used to generate the data:

\itemize{

\item \code{factors} --- Number of factors

\item \code{variables} --- Variables on each factor

\item \code{loadings} --- Loading matrix

\item \code{factor_correlations} --- Correlations between factors

\item \code{categories} --- Categories for each variable

\item \code{skew} --- Skew for each variable

}

}
}
\description{
Simulates data from a latent factor model based on many
manipulable parameters. Parameters do not have default values and
must each be set. See examples to get started
}
\examples{
# Generate factor data
two_factor <- simulate_factors(
  factors = 2, # factors = 2
  variables = 6, # variables per factor = 6
  loadings = 0.55, # loadings between = 0.45 to 0.65
  cross_loadings = 0.05, # cross-loadings N(0, 0.05)
  correlations = 0.30, # correlation between factors = 0.30
  sample_size = 1000 # number of cases = 1000
)

# Randomly vary loadings
two_factor_loadings <- simulate_factors(
  factors = 2, # factors = 2
  variables = 6, # variables per factor = 6
  loadings_range = c(0.30, 0.80), # loadings between = 0.30 to 0.80
  cross_loadings = 0.05, # cross-loadings N(0, 0.05)
  correlations = 0.30, # correlation between factors = 0.30
  sample_size = 1000 # number of cases = 1000
)

# Generate dichotomous data
two_factor_dichotomous <- simulate_factors(
  factors = 2, # factors = 2
  variables = 6, # variables per factor = 6
  loadings = 0.55, # loadings between = 0.45 to 0.65
  cross_loadings = 0.05, # cross-loadings N(0, 0.05)
  correlations = 0.30, # correlation between factors = 0.30
  sample_size = 1000, # number of cases = 1000
  variable_categories = 2 # dichotomous data
)

# Generate dichotomous data with skew
two_factor_dichotomous_skew <- simulate_factors(
  factors = 2, # factors = 2
  variables = 6, # variables per factor = 6
  loadings = 0.55, # loadings between = 0.45 to 0.65
  cross_loadings = 0.05, # cross-loadings N(0, 0.05)
  correlations = 0.30, # correlation between factors = 0.30
  sample_size = 1000, # number of cases = 1000
  variable_categories = 2, # dichotomous data
  skew = 1 # all variables with have a skew of 1
)

# Generate dichotomous data with variable skew
two_factor_dichotomous_skew <- simulate_factors(
  factors = 2, # factors = 2
  variables = 6, # variables per factor = 6
  loadings = 0.55, # loadings between = 0.45 to 0.65
  cross_loadings = 0.05, # cross-loadings N(0, 0.05)
  correlations = 0.30, # correlation between factors = 0.30
  sample_size = 1000, # number of cases = 1000
  variable_categories = 2, # dichotomous data
  skew_range = c(-2, 2) # skew = -2 to 2 (increments of 0.05)
)

}
\references{
Garrido, L. E., Abad, F. J., & Ponsoda, V. (2011). \cr
Performance of Velicer’s minimum average partial factor retention method with categorical variables. \cr
\emph{Educational and Psychological Measurement}, \emph{71}(3), 551-570.

Golino, H., Shi, D., Christensen, A. P., Garrido, L. E., Nieto, M. D., Sadana, R., ... & Martinez-Molina, A. (2020).
Investigating the performance of exploratory graph analysis and traditional techniques to identify the number of latent factors: A simulation and tutorial.
\emph{Psychological Methods}, \emph{25}(3), 292-320.
}
\author{
Maria Dolores Nieto Canaveras <mnietoca@nebrija.es>,
Alexander P. Christensen <alexpaulchristensen@gmail.com>,
Hudson Golino <hfg9s@virginia.edu>,
Luis Eduardo Garrido <luisgarrido@pucmm.edu>
}
