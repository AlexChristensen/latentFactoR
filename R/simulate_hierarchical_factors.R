simulate_hierarchical_factors <- function(
    lower_factors, variables, variables_range = NULL,
    lower_loadings, lower_loadings_range = NULL,
    lower_cross_loadings, lower_cross_loadings_range = NULL,
    higher_factors, higher_loadings, higher_loadings_range = NULL,
    higher_cross_loadings, higher_cross_loadings_range = NULL,
    higher_correlations, higher_correlations_range = NULL,
    higher_disturbances, higher_disturbances_range = NULL,
    sample_size, variable_categories = Inf,
    categorical_limit = 7,
    skew = 0, skew_range = NULL
)
{

  # Check for admissibility
  if(
    any((diag(higher_disturbances) < 0) | (diag(higher_disturbances) < 1)) |
    any(eigen(higher_disturbances, symmetric = TRUE, only.values = TRUE) < 0)
  ){

    # Send error about inadmissible

  }

  # Compute lower order correlations
  lower_correlations <- higher_loadings %*% higher_correlations %*% t(higher_loadings) + higher_disturbances

  # Obtain correlation matrix
  population_correlation <- lower_loadings %*% lower_correlations %*% lower_loadings

}