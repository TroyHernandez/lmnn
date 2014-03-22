# L1NNEquilibrium.R

#' Runs the L1NN algorithm repeatedly to find optimal parameters following
#' possible change of target neighbors
#'
#' @param x a matrix of numeric variables
#' @param y a character vector of class labels
#' @param mu determines what share of minimization vector is on the weights and
#' what share is on the errors.
#' @param margin is size of margin
#' @param lambda induces additional sparsity/variable selection
#' @param k number of neighbors to compare
#' @param min.nb logical, determines whether target neighbor distance
#' is minimized
#' @param equilibrium should the L1NNEquilibrium function be used to compute
#' final reweighting
#' @param output logical indicating if intermediate output should printed 
#' @param num.its integer indicating maximum number of times to rerun L1NN 
#' @param take.mean logical indicating if mean of reweightings should be taken
#' or FALSE to use only the last computed L1NN 
#' @param cv should the L1NN.cv function be used to compute final reweighting
#' @keywords metric learning
#' @export
#'
L1NNEquilibrium <- function(x = x, y = y, mu = .5, margin = 1,
                            lambda = 0, k = 1, min.nb = T, output = F,
                            num.its = 10, take.mean = TRUE, cv = TRUE) {
  if(cv == FALSE){
    knn.wts <- L1NN(x, y, wt.vec = rep(1, ncol(x)), min.len = 2, mu, margin,
                    lambda, k, min.nb, output)$weights

  }else{
    l1nn.fit <- L1NN.cv(x, y, equilibrium = FALSE, n.fold = 10)
    l1nn.pars <- l1nn.fit$pars
    knn.wts <- L1NN(x, y, wt.vec = rep(1, ncol(x)), min.len = 2,
                    mu = l1nn.pars["mu"],
                    margin = l1nn.pars["margin"],
                    lambda = l1nn.pars["lambda"],
                    k = l1nn.pars["k"],
                    min.nb = l1nn.pars["min.nb"])$weights
  }
  
  knn.wts2 <- rep(1, ncol(x))
  
  wt.mat <- matrix(0, nrow = num.its, ncol = ncol(x))
  wt.mat[1,] <- knn.wts
  l1.it <- 1
  
  while(sum(knn.wts!=knn.wts2) > 0 & l1.it < num.its){
    
    knn.wts2 <- knn.wts
    if(cv == FALSE){
      knn.wts <- L1NN(x, y, wt.vec = rep(1, ncol(x)), min.len = 2, mu, margin,
                      lambda, k, min.nb, output)$weights
    }else{
      knn.wts <- L1NN(x, y, wt.vec = rep(1, ncol(x)), min.len = 2,
                      mu = l1nn.pars["mu"],
                      margin = l1nn.pars["margin"],
                      lambda = l1nn.pars["lambda"],
                      k = l1nn.pars["k"],
                      min.nb = l1nn.pars["min.nb"])$weights
    }
    l1.it <- l1.it+1
    if(sum(knn.wts == 1) == ncol(x)){
      knn.wts <- knn.wts2
    }
    wt.mat[l1.it,] <- knn.wts
  }
  
  if(take.mean == TRUE && l1.it > 1){
    # take the mean of the last 2/3rds of calculated wts
    knn.wts <- apply(wt.mat[ceiling(l1.it * 1 / 3):l1.it, ], 2, mean)
  }

  knn.wts
}