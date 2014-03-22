# L1NNCV.R

#' Runs cross validation on the L1NN algorithm to find optimal parameters
#'
#' This function uses random cross validation (as opposed to grid) to find the 
#' optimal hyper-parameters for the L1NN function.
#'
#' @param x a matrix of numeric variables
#' @param y a character vector of class labels
#' @param num.cv integer indicated the number of random cv parameter sets to try
#' @param seed for reproducibility
#' @param nfold integer for the number of folds to used to compute cv errors
#' @param mu determines what share of minimization vector is on the weights and
#' what share is on the errors.
#' @param margin is size of margin
#' @param lambda induces additional sparsity/variable selection
#' @param k number of neighbors to compare
#' @param min.nb logical, determines whether target neighbor distance
#' is minimized
#' @param equilibrium should the L1NNEquilibrium function be used to compute
#' final reweighting
#' @keywords metric learning
#' @export
#'
L1NN.cv <- function(x, y, num.cv = 10, seed = 1, n.fold = 10,
                    mu = runif(num.cv, 0, 1),
                    margin = 10 ^ runif(num.cv, -3, 3),
                    lambda = 10 ^ runif(num.cv, -3, 0),
                    k = 1,
                    min.nb = sample(c(0, 1, 1), num.cv, replace = T),
                    equilibrium = TRUE) {
  
  n.fold <- min(n.fold, min(table(y)))
  parmat <- cbind(mu, margin, lambda, k, min.nb)
  colnames(parmat) <- c("mu", "margin", "lambda", "k", "min.nb")
  
  #Partition data here.
  library(caret)
  set.seed(seed)
  train.ind <- createFolds(y, k = n.fold, list = T)
  #   stop("here")
  errmat <- matrix(0, n.fold, num.cv)
  
  for(i in 1:n.fold){
    for(j in 1:num.cv){
      #       cat(dim(errmat),i,j,"\n")
      training <- unlist(train.ind[-i])
      testing <- train.ind[i][[1]]
      if(equilibrium == TRUE){
        temp.wts <- L1NNEquilibrium(x = x[training, ], y = y[training],
                                    mu = parmat[j, "mu"],
                                    margin = parmat[j, "margin"],
                                    lambda = parmat[j, "lambda"],
                                    k = parmat[j, "k"],
                                    min.nb = parmat[j, "min.nb"],
                                    output = F, cv = FALSE)
      }else{
        temp.wts <- L1NN(x = x[training, ], y = y[training], min.len = 1,
                         mu = parmat[j, "mu"],
                         margin = parmat[j, "margin"],
                         lambda = parmat[j, "lambda"],
                         k = parmat[j, "k"],
                         min.nb = parmat[j, "min.nb"],
                         output = F)$weights
      }
      fit <- knn(train = t(t(x[training, ]) * temp.wts),
                 test = t(t(x[testing, ]) * temp.wts),
                 cl = y[training], k = as.numeric(parmat[j, "k"]))
      errmat[i, j] <- sum(fit != y[testing]) / length(training)
      #       , levels = newlevels))
    }
  }
  
  #Determine best parameters
  errmat2 <- apply(errmat, 1, mean)
  
  best <- which(errmat2 == min(errmat2), arr.ind = T)
  if(length(best) == 2){
    best <- sample(best,1)
  }
  if(length(best) > 2){
    cands <- parmat[best,]
    #Check to remove uniform columns
    if(sum(apply(cands, 2, sd) == 0) > 0){
      cands <- cands[, -which(apply(cands, 2, sd) == 0)]
      mahal <- which.min(as.matrix(dist(rbind(apply(cands, 2, mean), cands),
                                        diag = T, upper = T))[1, -1])
    }else{
      mahal <- which.min(as.matrix(dist(rbind(apply(cands, 2, mean), cands),
                                        diag = T, upper = T))[1, -1])
    }
    best <- best[mahal]
  }
  
  list(num.cv = num.cv, pars = parmat[best, ])
}
