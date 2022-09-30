#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#### latentFactoR S3Methods ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Updated 30.09.2022

# print() Methods ----

# Print `lf_estimate`
# Updated 30.09.2022
#' @export
print.lf_estimate <- function(x, ...)
{
  
  # Print dimensions
  x$dimensions
  
}

# summary() Methods ----

# Summary `lf_estimate`
# Updated 30.09.2022
#' @export
summary.lf_estimate <- function(object, ...)
{
  
  # Print dimensions
  object$dimensions
  
}