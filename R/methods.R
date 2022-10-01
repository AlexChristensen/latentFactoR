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
  print(x$dimensions)
  
}

# summary() Methods ----

# Summary `lf_estimate`
# Updated 30.09.2022
#' @export
summary.lf_estimate <- function(object, ...)
{
  
  # Print dimensions
  print(object$dimensions)
  
}

# predict() Methods ----
# {xgboost}: predictLearner
# Updated 30.09.2022
#' @export
predictLearner.classif.xgboost.earlystop = function(.learner, .model, .newdata, ...) {
  td = .model$task.desc
  m = .model$learner.model
  cls = td$class.levels
  nc = length(cls)
  obj = .learner$par.vals$objective
  
  if (is.null(obj))
    .learner$par.vals$objective = ifelse(nc == 2L, "binary:logistic", "multi:softprob")
  
  p = predict(m, newdata = data.matrix(BBmisc::convertDataFrameCols(.newdata, ints.as.num = TRUE)), ...)
  
  if (nc == 2L) { #binaryclass
    if (.learner$par.vals$objective == "multi:softprob") {
      y = matrix(p, nrow = length(p) / nc, ncol = nc, byrow = TRUE)
      colnames(y) = cls
    } else {
      y = matrix(0, ncol = 2, nrow = nrow(.newdata))
      colnames(y) = cls
      y[, 1L] = 1 - p
      y[, 2L] = p
    }
    if (.learner$predict.type == "prob") {
      return(y)
    } else {
      p = colnames(y)[max.col(y)]
      names(p) = NULL
      p = factor(p, levels = colnames(y))
      return(p)
    }
  } else { #multiclass
    if (.learner$par.vals$objective  == "multi:softmax") {
      return(factor(p, levels = cls)) #special handling for multi:softmax which directly predicts class levels
    } else {
      p = matrix(p, nrow = length(p) / nc, ncol = nc, byrow = TRUE)
      colnames(p) = cls
      if (.learner$predict.type == "prob") {
        return(p)
      } else {
        ind = max.col(p)
        cns = colnames(p)
        return(factor(cns[ind], levels = cns))
      }
    }
  }
}