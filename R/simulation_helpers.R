#%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### Simulation Helpers ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Testing code
#
# # NEO Openness to Experience (n = 800)
# neo <- NetworkToolbox::neoOpen[1:800,]
# 
# # Create groups
# neo <- cbind(
#   rep(1:4, 200),
#   neo
# )
# 
# # Change column name
# colnames(neo)[1] <- "group"
# 
# # Obtain effects
# effect_table_object <- effect_table(
#   formula = Val8 ~ Val7 * Act4 * group,
#   data = neo, method = "group"
# )
# 
# # Plot effects
# effect_plot(effect_table_object)
#

# These functions are used internally to facilitate
# the use of functions that are often used in 
# simulations we run. They are not intended for users
# although they may be useful to some.

#' Creates a plot of the effects table
#' @noRd
# Updated 20.12.2022
effect_plot <- function(
    effect_table_object,
    title = "", subtitle = "",
    fill_color = "blue",
    produce = TRUE
)
{
  
  # Create long format data frame
  long_df <- data.frame(
    Method = rep(
      colnames(effect_table_object), # Methods
      each = nrow(effect_table_object) # Effects
    ),
    Term = rep(
      row.names(effect_table_object), # Effects
      times = ncol(effect_table_object) # Methods
    ),
    Eta = as.vector(
      as.matrix(effect_table_object) # Values
    )
  )
  
  # Set modes
  long_df$Method <- as.factor(long_df$Method)
  long_df$Term <- factor(
    row.names(effect_table_object),
    levels = rev(row.names(effect_table_object))
  )
  long_df$Eta <- as.numeric(
    as.character(long_df$Eta)
  )
  
  # Plot
  plot_to_return <- ggplot2::ggplot(
    data = long_df,
    ggplot2::aes(
      x = Method, y = Term,
      fill = Eta, label = Eta
    )
  ) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(
      low = "white", high = fill_color,
      limits = c(0, 1), 
      name = expression(
        "\u03B7"[p]^2
      )
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::scale_x_discrete(
      position = "bottom", expand = c(0,0)
    ) +
    ggplot2::scale_y_discrete(expand = c(0,0)) + 
    ggplot2::geom_text(size = 6) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 16, hjust = 0.5),
      panel.background = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = 16),
      axis.text = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.ticks = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = 14, hjust = 0.5),
      legend.text = ggplot2::element_text(size = 12),
      legend.key.height = ggplot2::unit(1, "cm")
    )
  
  # Check if plot should be produced
  if(isTRUE(produce)){
    plot(plot_to_return)
  }
  
  # Message to user
  message(
    paste(
      "Use standard {ggplot2} functionality to manipulate the plot\n\n",
      "`aes()` settings:\n",
      "x = Method\n", "y = Term\n", "fill = Eta\n", "label = Eta\n"
    )
  )
  
  # Return plot
  return(plot_to_return)
  
}

#' Creates an effects table for all methods
#' @noRd
# Updated 20.12.2022
effect_table <- function(
    formula, data, method,
    minimum_effect = c("small", "moderate", "large"),
    interactions_limit = 3,
    progress = TRUE,
    return_all = FALSE
)
{
  
  # Check that method exists
  if(!method %in% colnames(data)){
    stop(
      paste0(
        "'", method, "' was not found in the `data`. ",
        "Make sure the `method` argument is a column name in data"
      )
    )
  }
  
  # Identify unique methods
  unique_methods <- unique(
    data[,method]
  )
  
  # Determine character length (for progress)
  if(isTRUE(progress)){
    method_characters <- nchar(as.character(unique_methods))
  }
  
  # Check for only one method
  if(length(unique_methods) == 1){
    
    message(
      paste0(
        "Only one method, '", unique_methods, "', was found ",
        "in the data. Results will only include this method."
      )
    )
    
  }
  
  # Separate data into lists
  separate_data <- lapply(unique_methods, function(x){
    data[data[,method] == x,]
  })
  
  # Perform ANOVAs on separate data
  results <- lapply(
    seq_along(separate_data), function(i){
      
      # Check for progress
      if(isTRUE(progress)){
        
        # Progress message
        cat(
          
          colortext(
            paste0(
              "\r Obtaining effects for '",
              unique_methods[i], "' (",
              i,
              " of ", length(separate_data), ")...",
              paste0(
                rep(" ", abs(diff(range(method_characters)))),
                collapse = ""
              )
            ),
            defaults = "message"
          )
          
        )
        
      }
      
      # Obtain effects
      obtain_effect_sizes(
        formula = formula,
        data = separate_data[[i]],
        minimum_effect = minimum_effect,
        return_all = TRUE
      )
      
    }
  )
  
  # Stack results
  stacked_results <- t(do.call(
    cbind.data.frame, results
  ))
  
  # Rename rows
  row.names(stacked_results) <- unique_methods
  
  # Check for returning all
  if(isTRUE(return_all)){
    return(t(stacked_results))
  }else{
    
    # If minimum effect is missing, set to large
    if(missing(minimum_effect)){
      minimum_effect <- "large"
    }else{minimum_effect <- tolower(match.arg(minimum_effect))}
    
    # Check for minimum effect
    if(is.character(minimum_effect)){
      
      effect_size <- switch(
        minimum_effect,
        "small" = 0.01,
        "moderate" = 0.06,
        "large" = 0.14
      )
      
    }else if(is.numeric(minimum_effect)){
      effect_size <- minimum_effect
    }
    
    # Ensure that at least one method has
    # larger than minimum effect
    conditions_to_return <- apply(
      stacked_results,
      2, function(x){
        any(x >= effect_size)
      }
    )

    # Check for no results
    if(sum(conditions_to_return) == 0){
      
      # Return message
      message(
        paste0(
          "No effect sizes were greater than ",
          minimum_effect, ".", "\n",
          "Try setting `return_all = TRUE` to see all effect sizes."
        )
      )
      
      # Return NULL
      return(NULL)
      
    }else{
      
      # Obtain terms
      terms <- names(conditions_to_return[conditions_to_return])
      
      # Obtain length of terms
      term_lengths <- unlist(
        lapply(strsplit(terms, split = ":"), length)
      )
      
      # Check for interaction limit
      if(any(term_lengths > interactions_limit)){
        conditions_to_return[
          terms[
            term_lengths > interactions_limit
          ]
        ] <- FALSE
      }
      
      # Return rounded larger etas
      return(t(stacked_results[,conditions_to_return]))
      
    }
    
  }
  
}

#' All-in-one ANOVA and partial eta squared
#' @importFrom stats aov as.formula
#' @noRd
# Updated 12.12.2022
obtain_effect_sizes <- function(
  formula, data,
  minimum_effect = c("small", "moderate", "large"),
  return_all = FALSE
)
{
  
  # If minimum effect is missing, set to large
  if(missing(minimum_effect)){
    minimum_effect <- "large"
  }else{minimum_effect <- tolower(match.arg(minimum_effect))}
  
  # Check for minimum effect
  if(is.character(minimum_effect)){
    
    effect_size <- switch(
      minimum_effect,
      "small" = 0.01,
      "moderate" = 0.06,
      "large" = 0.14
    )
    
  }else if(is.numeric(minimum_effect)){
    effect_size <- minimum_effect
  }
  
  # Obtain ANOVA formula
  formula <- as.formula(formula)
  
  # Perform ANOVA
  aov_object <- aov(
    formula = formula,
    data = data
  )
  
  # Obtain partial etas
  etas <- partial_eta(aov_object)
  
  # Check for returning all
  if(isTRUE(return_all)){
    return(round(etas, 3))
  }else{
    
    # Determine etas greater than or equal to minimum effect
    larger_etas <- etas[which(round(etas, 2) >= effect_size)]
    
    # Check for no results
    if(length(larger_etas) == 0){
      
      # Return message
      message(
        paste0(
          "No effect sizes were greater than ",
          minimum_effect, ".", "\n",
          "Try setting `return_all = TRUE` to see all effect sizes."
        )
      )
      
      # Return NULL
      return(NULL)
      
    }else{
      
      # Return rounded larger etas
      return(round(larger_etas, 3))
      
    }
    
  }

}

#' Partial eta squared
#' @noRd
# Updated 13.12.2022
partial_eta <- function(aov_object)
{
  
  # Obtain summary
  aov_summary <- summary(aov_object)
  
  # Obtain variable names
  variable_names <- trimws(row.names(aov_summary[[1]]))
  
  # Obtain sum of squares
  sum_of_squares <- aov_summary[[1]]$`Sum Sq`
  
  # Obtain residuals
  residuals <- sum_of_squares[length(sum_of_squares)]
  
  # Obtain partial eta squares
  etas <- sum_of_squares / (sum_of_squares + residuals)
  
  # Add names
  names(etas) <- variable_names
  
  # Remove residuals
  etas <- etas[-length(etas)]
  
  # Return partial eta squares
  return(etas)
  
}

#' Check levels available in a given dataset
#' @noRd
# Updated 29.12.2022
check_levels <- function(data, limit = 10)
{
  
  # Ensure data is a data frame
  data <- as.data.frame(data)
  
  # Determine variables within a reasonable limit
  variable_lengths <- apply(data, 2, function(x){
    length(na.omit(unique(x)))
  })
  
  # Check levels under or equal to limit
  variable_levels <- apply(
    data[,variable_lengths <= limit],
    2, function(x){
      unique(x)
    }
  )
  
  # Return levels
  return(variable_levels)
  
}
