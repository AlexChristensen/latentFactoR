#' Checks Implied Lower-Order Correlations Across a Grid of Higher-Order Values
#'
#' @description
#' Computes the lower-order factor correlations implied by
#' \code{\link[latentFactoR]{simulate_hierarchical_factors}}'s higher-order factor structure
#' across a grid of \code{higher_correlations} values (and, optionally, \code{higher_loadings}
#' and/or off-diagonal \code{off_disturbances} values), without generating any sample data.
#' Implied correlations are split into pairs of lower-order factors that share their dominant
#' higher-order factor and pairs that do not, so their diverging response to
#' \code{higher_correlations} can be compared directly. Intended as a lightweight way to explore
#' this relationship before committing to a full (and comparatively costly) call to
#' \code{\link[latentFactoR]{simulate_hierarchical_factors}}. See examples to get started
#'
#' @details
#' For each value of \code{higher_correlations_grid} (crossed with each value of
#' \code{higher_loadings_grid} and/or \code{off_disturbances_grid}, when supplied), the implied
#' lower-order correlation matrix is obtained the same way it is inside
#' \code{\link[latentFactoR]{simulate_hierarchical_factors}}:
#'
#' \deqn{\boldsymbol{\Phi}_L = \boldsymbol{\Lambda}_H \boldsymbol{\Phi}_H \boldsymbol{\Lambda}_H' + \boldsymbol{\Psi}}{Phi_L = Lambda_H * Phi_H * Lambda_H' + Psi}
#'
#' \code{higher_loadings} (and \code{higher_variables}, if supplied) is used to fix which
#' higher-order factor each lower-order factor dominantly loads on (\code{which.max} of its row
#' in the higher-order loading matrix). That pattern is held fixed throughout. When
#' \code{higher_loadings_grid} is left \code{NULL}, the dominant loading magnitudes from
#' \code{higher_loadings} are also held fixed (at their point value, with no random jitter or
#' cross-loadings) so that \code{higher_correlations} is the only thing varying across the grid.
#' When \code{higher_loadings_grid} is supplied instead, each of its values is applied uniformly,
#' in turn, to every dominant loading (overriding \code{higher_loadings}' magnitudes, though its
#' pattern is still used), crossing loading magnitude with \code{higher_correlations}
#'
#' \eqn{\boldsymbol{\Psi}}{Psi} (the disturbance matrix) is only included when
#' \code{off_disturbances_grid} is supplied. When it is left \code{NULL} (the default),
#' \eqn{\boldsymbol{\Psi}}{Psi} is omitted entirely: as an uncorrelated/diagonal variance, as it is
#' whenever a single value is supplied to \code{\link[latentFactoR]{simulate_hierarchical_factors}}'s
#' \code{off_disturbances}, it only ever adds to the diagonal of
#' \eqn{\boldsymbol{\Phi}_L}{Phi_L}, which is then rescaled back to 1 regardless of its value,
#' leaving the off-diagonal --- i.e., the actual lower-order correlations --- entirely unaffected.
#' Each value of \code{off_disturbances_grid}, by contrast, is used as an \emph{off-diagonal}
#' (correlated) disturbance value: it is applied uniformly to every off-diagonal element of
#' \eqn{\boldsymbol{\Psi}}{Psi} (with a unit diagonal), so it \emph{does} shift the implied
#' lower-order correlations directly, and can push them outside of what \code{higher_loadings} and
#' \code{higher_correlations} alone would imply (including, for large enough values, outside of
#' \code{[-1, 1]} or into inadmissible territory)
#'
#' Every pair of lower-order factors is classified as sharing their dominant higher-order factor
#' or not, and the implied correlation for each grid point is averaged separately within each of
#' those two groups. Each grid point is also checked for admissibility: whether the implied
#' lower-order correlation matrix is positive semi-definite (the same check performed inside
#' \code{\link[latentFactoR]{simulate_hierarchical_factors}}), and, when
#' \code{off_disturbances_grid} is supplied, whether \eqn{\boldsymbol{\Psi}}{Psi} itself is a
#' valid (positive semi-definite) covariance matrix
#'
#' @param lower_factors Numeric (length = 1).
#' Number of lower-order (first-order) factors.
#' Minimum two (implied lower-order correlations are not defined for a single lower-order factor)
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
#' higher-order factors. Values must sum to \code{lower_factors}.
#' Ignored when \code{higher_loadings} is supplied as a matrix
#'
#' @param higher_loadings Numeric or matrix (length = 1, \code{higher_factors}, or \code{lower_factors} x \code{higher_factors}).
#' Loadings of lower-order factors onto their higher-order factor, used to fix which higher-order
#' factor each lower-order factor dominantly loads on.
#' Can be a single value or as many values as there are higher-order factors (corresponding to
#' the factors). Can also be a loading matrix: columns must match \code{higher_factors} and rows
#' must match \code{lower_factors}. Unlike \code{\link[latentFactoR]{simulate_hierarchical_factors}},
#' values are used exactly as supplied (no random jitter and no cross-loadings). When
#' \code{higher_loadings_grid} is supplied, only the \emph{pattern} of dominant loadings implied
#' by \code{higher_loadings} is retained --- the magnitudes are overridden by the grid.
#' Defaults to \code{NULL}, which is only valid when \code{higher_loadings_grid} is supplied
#' instead (its magnitude would have been overridden by the grid anyway); one of
#' \code{higher_loadings} or \code{higher_loadings_grid} must be supplied
#'
#' @param higher_loadings_grid Numeric.
#' Values of \code{higher_loadings} to check, crossed with \code{higher_correlations_grid} (and
#' \code{off_disturbances_grid}, when supplied).
#' Each value is applied uniformly to every lower-order factor's dominant loading (the pattern of
#' which higher-order factor is dominant for each lower-order factor still comes from
#' \code{higher_loadings}/\code{higher_variables}).
#' Defaults to \code{NULL}, which holds loadings fixed at \code{higher_loadings} and does not
#' facet the plot by loading value
#'
#' @param higher_correlations_grid Numeric.
#' Values of \code{higher_correlations} to check.
#' Each value is applied uniformly to all off-diagonal correlations between higher-order factors
#' (as when a single value is supplied to \code{\link[latentFactoR]{simulate_hierarchical_factors}}'s
#' \code{higher_correlations})
#'
#' @param off_disturbances_grid Numeric.
#' Values of \code{off_disturbances} to check, crossed with \code{higher_correlations_grid}
#' (and \code{higher_loadings_grid}, when supplied). Each value is applied uniformly to every
#' \emph{off-diagonal} element of the lower-order disturbance matrix (with a unit diagonal), i.e.,
#' it represents a correlation between lower-order factors' disturbances, and is the only way this
#' function's grid can shift the implied lower-order correlations directly (see Details).
#' Defaults to \code{NULL}, which omits disturbances entirely (equivalent to uncorrelated
#' disturbances, which do not affect the implied lower-order correlations) and does not facet the
#' plot by disturbance value
#'
#' @param plot Boolean (length = 1).
#' Whether a \{ggplot2\} line plot of the implied lower-order correlations across the grid should
#' be produced and printed. When \code{higher_loadings_grid} and/or \code{off_disturbances_grid}
#' are supplied, the plot is faceted by loading and/or disturbance value.
#' Defaults to \code{TRUE}
#'
#' @param title Character (length = 1).
#' Title of the plot.
#' Defaults to \code{""}
#'
#' @param subtitle Character (length = 1).
#' Subtitle of the plot.
#' Defaults to \code{""}
#'
#' @param within_color Character (length = 1).
#' Color of the line for pairs of lower-order factors that share their dominant higher-order
#' factor.
#' Defaults to \code{"#A3C743"}
#'
#' @param between_color Character (length = 1).
#' Color of the line for pairs of lower-order factors that do \emph{not} share their dominant
#' higher-order factor.
#' Defaults to \code{"#F1AD41"}
#'
#' @return Returns a list containing:
#'
#' \item{grid}{Data frame with one row per \code{higher_correlations_grid} value (crossed with
#' \code{higher_loadings_grid} and/or \code{off_disturbances_grid}, when supplied), containing
#' the grid value(s) (\code{higher_correlation} and, when supplied, \code{higher_loading} and/or
#' \code{off_disturbance}), the mean implied correlation among lower-order factor pairs that
#' share their dominant higher-order factor (\code{within_correlation}), the mean implied
#' correlation among pairs that do not (\code{between_correlation}), the smallest eigenvalue of the
#' implied lower-order correlation matrix (\code{minimum_eigenvalue}), and whether that matrix ---
#' and, when \code{off_disturbances_grid} is supplied, the disturbance matrix itself --- is
#' admissible (positive semi-definite; \code{admissible})}
#'
#' \item{plot}{\{ggplot2\} line plot of \code{within_correlation} and \code{between_correlation}
#' across the grid, faceted by \code{higher_loading} and/or \code{off_disturbance} when
#' supplied (\code{NULL} if \code{plot = FALSE})}
#'
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' # Check implied lower-order correlations across a grid of higher-order correlations
#' implied <- implied_hierarchical_values(
#'   lower_factors = 4, higher_factors = 2,
#'   higher_loadings = 0.60, # higher-order loadings fixed at 0.60
#'   higher_correlations_grid = seq(0, 0.80, 0.10)
#' )
#'
#' # Also grid over higher-order loadings, crossed with higher-order correlations
#' implied_loadings <- implied_hierarchical_values(
#'   lower_factors = 4, higher_factors = 2,
#'   higher_loadings = 0.60, # only its dominant-loading pattern is used
#'   higher_loadings_grid = seq(0.40, 0.80, 0.20),
#'   higher_correlations_grid = seq(0, 0.80, 0.20)
#' )
#'
#' # `higher_loadings` can be omitted entirely when `higher_loadings_grid` is supplied
#' implied_loadings_only <- implied_hierarchical_values(
#'   lower_factors = 4, higher_factors = 2,
#'   higher_loadings_grid = seq(0.40, 0.80, 0.10),
#'   higher_correlations_grid = seq(0, 0.80, 0.20)
#' )
#'
#' # Also grid over off-diagonal (correlated) higher-order disturbances
#' implied_disturbances <- implied_hierarchical_values(
#'   lower_factors = 4, higher_factors = 2,
#'   higher_loadings = 0.60,
#'   higher_correlations_grid = seq(0, 0.80, 0.20),
#'   off_disturbances_grid = seq(-0.20, 0.20, 0.20)
#' )
#'
#' # Uneven number of lower-order factors nested within each higher-order factor
#' implied_uneven <- implied_hierarchical_values(
#'   lower_factors = 6, higher_factors = 2, higher_variables = c(2, 4),
#'   higher_loadings = 0.60,
#'   higher_correlations_grid = seq(0, 0.80, 0.10)
#' )
#'
#' @author
#' Alexander P. Christensen <alexpaulchristensen@gmail.com>
#'
#' @noRd
#'
# Grid check for implied hierarchical values
# Updated 13.07.2026
implied_hierarchical_values <- function(
    lower_factors, higher_factors, higher_variables = NULL,
    higher_loadings = NULL, higher_loadings_grid = NULL,
    higher_correlations_grid,
    off_disturbances_grid = NULL,
    plot = TRUE, title = "", subtitle = "",
    within_color = "#A3C743", between_color = "#F1AD41"
)
{

  # Check inputs
  inputs <- implied_hierarchical_values_errors(
    lower_factors, higher_factors, higher_variables, higher_loadings, higher_loadings_grid,
    higher_correlations_grid, off_disturbances_grid, plot, title, subtitle,
    within_color, between_color
  )

  # Collect inputs
  higher_variables <- inputs$higher_variables

  # When `higher_loadings` is not supplied, `higher_loadings_grid` must be (enforced in
  # `implied_hierarchical_values_errors`); its magnitude never matters in that case (only
  # the pattern of which higher-order factor is dominant for each lower-order factor does,
  # and that comes entirely from `higher_variables`), so a placeholder value stands in for it
  if(is.null(higher_loadings)){higher_loadings <- 1}

  # Build higher-order loading matrix (lower-order factors onto higher-order factors),
  # fixed at its point value (no jitter, no cross-loadings). Used to fix the pattern of
  # dominant loadings, and, when `higher_loadings_grid` is not supplied, their magnitude too
  if(is(higher_loadings, "matrix")){

    higher_loading_matrix <- higher_loadings

  }else{

    # Repeat single value across higher-order factors
    if(length(higher_loadings) == 1){higher_loadings <- rep(higher_loadings, higher_factors)}

    # Create starting and ending of lower-order factor sequences on higher-order factors
    higher_variable_sequence <- factor_variable_sequence(higher_variables)
    start_higher <- higher_variable_sequence$start
    end_higher <- higher_variable_sequence$end

    # Initialize higher-order loading matrix
    higher_loading_matrix <- matrix(data = 0, nrow = lower_factors, ncol = higher_factors)

    # Populate dominant loadings (exact point values; no jitter, no cross-loadings)
    for(i in seq_len(higher_factors)){
      higher_loading_matrix[start_higher[i]:end_higher[i], i] <- higher_loadings[i]
    }

  }

  # Identify each lower-order factor's dominant higher-order factor (largest loading)
  dominant_factor <- apply(higher_loading_matrix, 1, which.max)

  # Set up every unique pair of lower-order factors
  lower_factor_pairs <- combn(lower_factors, 2)

  # Identify whether each pair shares its dominant higher-order factor
  within_pair <- dominant_factor[lower_factor_pairs[1,]] == dominant_factor[lower_factor_pairs[2,]]

  # Warn when a group has no pairs to average over
  if(all(within_pair)){
    warning(
      "Every pair of lower-order factors shares the same dominant higher-order factor. ",
      "`between_correlation` cannot be computed and will be `NaN`"
    )
  }else if(all(!within_pair)){
    warning(
      "No pair of lower-order factors shares a dominant higher-order factor (each higher-order ",
      "factor has only one lower-order factor nested within it). `within_correlation` cannot be ",
      "computed and will be `NaN`"
    )
  }

  # Whether loading magnitude and/or off-diagonal disturbances are also being gridded
  grid_loadings <- !is.null(higher_loadings_grid)
  grid_disturbances <- !is.null(off_disturbances_grid)

  # Build grid (only cross in dimensions that were actually supplied)
  grid_dimensions <- list()
  if(isTRUE(grid_loadings)){grid_dimensions$higher_loading <- higher_loadings_grid}
  if(isTRUE(grid_disturbances)){grid_dimensions$off_disturbance <- off_disturbances_grid}
  grid_dimensions$higher_correlation <- higher_correlations_grid
  grid <- do.call(expand.grid, grid_dimensions)

  # Initialize columns for implied values
  grid$within_correlation <- NA_real_
  grid$between_correlation <- NA_real_
  grid$minimum_eigenvalue <- NA_real_
  grid$admissible <- NA

  # Loop through grid
  for(i in seq_len(nrow(grid))){

    # Loading matrix for this grid point: fixed pattern, magnitude either fixed or gridded
    if(isTRUE(grid_loadings)){
      loading_matrix <- matrix(data = 0, nrow = lower_factors, ncol = higher_factors)
      loading_matrix[cbind(seq_len(lower_factors), dominant_factor)] <- grid$higher_loading[i]
    }else{
      loading_matrix <- higher_loading_matrix
    }

    # Higher-order correlation matrix (single value applied to all off-diagonal correlations)
    higher_correlation_matrix <- matrix(
      data = grid$higher_correlation[i], nrow = higher_factors, ncol = higher_factors
    )
    diag(higher_correlation_matrix) <- 1

    # Assume always positive correlations (consistent with `simulate_hierarchical_factors`)
    higher_correlation_matrix <- abs(higher_correlation_matrix)

    # Lower-order disturbance matrix: uniform off-diagonal (correlated) value with a unit
    # diagonal, or an all-zero (i.e., no-op) matrix when disturbances are not being gridded
    if(isTRUE(grid_disturbances)){
      disturbance_matrix <- matrix(
        data = grid$off_disturbance[i], nrow = lower_factors, ncol = lower_factors
      )
      diag(disturbance_matrix) <- 1
      disturbance_admissible <- all(matrix_eigenvalues(disturbance_matrix) >= 0)
    }else{
      disturbance_matrix <- matrix(data = 0, nrow = lower_factors, ncol = lower_factors)
      disturbance_admissible <- TRUE
    }

    # Create lower-order factor correlations implied by the higher-order structure
    lower_correlation <- loading_matrix %*%
      tcrossprod(higher_correlation_matrix, loading_matrix) + disturbance_matrix

    # Ensure diagonal of lower-order correlation matrix is 1
    diag(lower_correlation) <- 1

    # Check eigenvalues for admissibility (positive semi-definite)
    eigenvalues <- matrix_eigenvalues(lower_correlation)

    # Pull implied correlation for every pair of lower-order factors
    pair_correlations <- lower_correlation[
      cbind(lower_factor_pairs[1,], lower_factor_pairs[2,])
    ]

    # Store implied values
    grid$within_correlation[i] <- mean(pair_correlations[within_pair])
    grid$between_correlation[i] <- mean(pair_correlations[!within_pair])
    grid$minimum_eigenvalue[i] <- min(eigenvalues)
    grid$admissible[i] <- disturbance_admissible && all(eigenvalues > 0)

  }

  # Initialize plot
  implied_plot <- NULL

  # Produce plot
  if(isTRUE(plot)){

    implied_plot <- implied_hierarchical_values_plot(
      grid, facet_loadings = grid_loadings, facet_disturbances = grid_disturbances,
      title = title, subtitle = subtitle,
      within_color = within_color, between_color = between_color
    )

  }

  # Return results
  return(
    list(
      grid = grid,
      plot = implied_plot
    )
  )

}

# Input checking ----
#' @noRd
# Updated 13.07.2026
implied_hierarchical_values_errors <- function(
    lower_factors, higher_factors, higher_variables, higher_loadings, higher_loadings_grid,
    higher_correlations_grid, off_disturbances_grid, plot, title, subtitle,
    within_color, between_color
)
{

  # Check lower-order factors first
  type_error(lower_factors, "numeric")
  length_error(lower_factors, 1)
  range_error(lower_factors, c(2, Inf))

  # Check higher-order factors second
  type_error(higher_factors, "numeric")
  length_error(higher_factors, 1)
  range_error(higher_factors, c(1, Inf))

  # Check for number of lower-order factors nested within each higher-order factor
  if(!is.null(higher_variables)){

    type_error(higher_variables, "numeric") # object type error
    length_error(higher_variables, c(1, higher_factors)) # object length error
    range_error(higher_variables, c(1, Inf)) # object range error

    # Repeat single value across higher-order factors
    if(length(higher_variables) == 1){higher_variables <- rep(higher_variables, higher_factors)}

  }else if(!is(higher_loadings, "matrix")){

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

  }

  # Ensure supplied lower-order factor counts account for all lower-order factors
  if(!is(higher_loadings, "matrix") && sum(higher_variables) != lower_factors){
    stop(
      paste0(
        "Sum of `higher_variables` (", sum(higher_variables), ") must equal `lower_factors` (",
        lower_factors, ")"
      )
    )
  }

  # `higher_loadings` is only optional when `higher_loadings_grid` is supplied instead
  # (its magnitude is irrelevant then; only its pattern would have mattered, and that comes
  # from `higher_variables`)
  if(is.null(higher_loadings) && is.null(higher_loadings_grid)){
    stop("Either `higher_loadings` or `higher_loadings_grid` must be supplied")
  }

  # Ensure appropriate type, length, and range for higher-order loadings
  if(!is.null(higher_loadings)){
    if(!is(higher_loadings, "matrix")){
      type_error(higher_loadings, "numeric")
      length_error(higher_loadings, c(1, higher_factors))
    }
    range_error(higher_loadings, c(-1, 1))
  }

  # Ensure appropriate type and range for higher-order loadings grid
  if(!is.null(higher_loadings_grid)){
    type_error(higher_loadings_grid, "numeric")
    range_error(higher_loadings_grid, c(-1, 1))
  }

  # Ensure appropriate type and range for higher-order correlations grid
  type_error(higher_correlations_grid, "numeric")
  range_error(higher_correlations_grid, c(-1, 1))

  # Ensure appropriate type and range for higher-order disturbances grid
  if(!is.null(off_disturbances_grid)){
    type_error(off_disturbances_grid, "numeric")
    range_error(off_disturbances_grid, c(-1, 1))
  }

  # Ensure appropriate types and lengths for plotting arguments
  type_error(plot, "logical"); length_error(plot, 1)
  type_error(title, "character"); length_error(title, 1)
  type_error(subtitle, "character"); length_error(subtitle, 1)
  type_error(within_color, "character"); length_error(within_color, 1)
  type_error(between_color, "character"); length_error(between_color, 1)

  # Return checked input
  return(
    list(higher_variables = higher_variables)
  )

}

# Plotting ----
#' @noRd
# Updated 13.07.2026
implied_hierarchical_values_plot <- function(
    grid, facet_loadings, facet_disturbances, title, subtitle, within_color, between_color
)
{

  # Initialize variables for {ggplot2}
  higher_correlation <- higher_loading <- off_disturbance <- implied_correlation <- group <- NULL

  # Set group labels
  within_label <- "Same Higher-Order Factor"
  between_label <- "Different Higher-Order Factors"

  # Create long format data frame
  long_df <- rbind(
    data.frame(
      higher_correlation = grid$higher_correlation,
      group = within_label,
      implied_correlation = grid$within_correlation
    ),
    data.frame(
      higher_correlation = grid$higher_correlation,
      group = between_label,
      implied_correlation = grid$between_correlation
    )
  )

  # Carry loading and/or disturbance grid values along for faceting
  if(isTRUE(facet_loadings)){long_df$higher_loading <- rep(grid$higher_loading, 2)}
  if(isTRUE(facet_disturbances)){long_df$off_disturbance <- rep(grid$off_disturbance, 2)}

  # Set factor levels (controls legend and drawing order)
  long_df$group <- factor(long_df$group, levels = c(within_label, between_label))

  # Plot
  plot_to_return <- ggplot2::ggplot(
    data = long_df,
    ggplot2::aes(x = higher_correlation, y = implied_correlation, color = group)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::geom_line(linewidth = 1, na.rm = TRUE) +
    ggplot2::geom_point(size = 2, na.rm = TRUE) +
    ggplot2::geom_vline(
      data = grid[!grid$admissible,],
      ggplot2::aes(xintercept = higher_correlation),
      linetype = "dashed", color = "grey50"
    ) +
    ggplot2::scale_color_manual(
      name = "Lower-Order Factor Pairs",
      values = setNames(c(within_color, between_color), c(within_label, between_label)),
    ) +
    ggplot2::labs(
      title = title, subtitle = subtitle,
      x = "Higher-Order Correlation", y = "Implied Lower-Order Correlation"
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = 8),
      labels = function(breaks){format_decimal(breaks, 2)}
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::pretty_breaks(n = 8),
      limits = c(swiftelse(isTRUE(facet_disturbances), min(long_df$off_disturbance), 0), 1),
      labels = function(breaks){format_decimal(breaks, 2)}
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 14, hjust = 0.5),
      panel.background = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10),
      axis.line = ggplot2::element_line(color = "black"),
      legend.position = "bottom",
      legend.title.position = "bottom",
      legend.title = ggplot2::element_text(size = 12, hjust = 0.5),
      legend.text = ggplot2::element_text(size = 10),
      strip.text = ggplot2::element_text(size = 10)
    )

  # Loading facet label
  loading_labeller <- function(values){
    paste0("Higher-Order Loading = ", format_decimal(as.numeric(values), 2))
  }

  # Disturbance facet label
  disturbance_labeller <- function(values){
    paste0("Disturbance = ", format_decimal(as.numeric(values), 2))
  }

  # Facet by whichever of loading magnitude and/or off-diagonal disturbances were also gridded
  if(isTRUE(facet_loadings) && isTRUE(facet_disturbances)){

    plot_to_return <- plot_to_return +
      ggplot2::facet_grid(
        off_disturbance ~ higher_loading,
        labeller = ggplot2::labeller(
          higher_loading = loading_labeller, off_disturbance = disturbance_labeller
        )
      )

  }else if(isTRUE(facet_loadings)){

    plot_to_return <- plot_to_return +
      ggplot2::facet_wrap(~ higher_loading, labeller = ggplot2::labeller(higher_loading = loading_labeller))

  }else if(isTRUE(facet_disturbances)){

    plot_to_return <- plot_to_return +
      ggplot2::facet_wrap(
        ~ off_disturbance, labeller = ggplot2::labeller(off_disturbance = disturbance_labeller)
      )

  }

  # Print plot
  plot(plot_to_return)

  # Return plot
  return(plot_to_return)

}
