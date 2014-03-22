# L1NN_1.R

############################################################
library(linprog)

#' An optimal L1 Large Margin reweighting function.
#'
#' This function computes the optimal weighting of variables to minimize the L1
#' distance for k-nearest neighbors algorithms.
#'
#' @param x a matrix of numeric variables
#' @param y a character vector of class labels 
#' @param eql.wt.vec determines which columns should share weights
#' @param mu determines what share of minimization vector is on the weights and
#' what share is on the errors.
#' @param margin is size of margin
#' @param lambda induces additional sparsity/variable selection
#' @param k number of neighbors to compare
#' @param min.nb logical, determines whether target neighbor distance
#' is minimized
#' @param wt.vec current weighting vector 
#' @keywords metric learning
#' @export
#'
#' 
L1NN <- function(x = x, y = y, wt.vec = rep(1, ncol(x)), min.len = 2,
                 mu = .5, margin = 1, lambda = 0, k = 1,
#                  eql.wt.vec = NULL,
                 min.nb = T, output = F) {
  
#   Catching eql.wt.vec error.
#   if(length(eql.wt.vec)!=ncol(x)){
#     cat("eql.wt.vec under development.\n")
#     eql.wt.vec <- c(1:ncol(x))
#   }
  
  clean.data <- .CleanData(x, y, min.len, k)
  x <- clean.data$x
  y <- clean.data$y

  if(length(table(y)) < 2){
    stop("y has only one labeled class.")
  }
  
  # Calculate optimal weights
  wts.errors <- .CalcWts(x, y, k, wt.vec, min.nb, lambda, margin, mu)# , eql.wt.vec)
  wts <- wts.errors[[1]]
  errors <- wts.errors[[2]]
  # Format output
  stats <- c(summary(errors), sd(errors))
  names(stats) <- c("Min", "1Q", "Mean", "Median", "3Q", "Max", "Sd")
  list(weights = wts, stats = stats)
}

