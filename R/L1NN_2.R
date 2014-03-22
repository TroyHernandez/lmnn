# L1NN_2.R

.CleanData <- function(x, y, min.len, k){
  
  rownames(x) <- 1:nrow(x)
  
  if(k > nrow(x)){
    stop("More neighbors than samples.\n")
  }

  #Removing sparse ys
  badlab <- names(which(table(y) < min.len))
  
  if(length(badlab) > 0){
    cat("Removing infrequent classes.")
    for(i in badlab){
      badind <- which(y == i)
      y <- y[-badind]
      x <- x[-badind, ]
    }
  }

  #Removing empty ys
  if(sum(y == "") > 0){
    cat("Removing \"\" class.")
    
    badind <- which(y == "")
    y <- y[-badind]
    x <- x[-badind, ]
  }
  
  if(sum(is.na(y)) > 0){
    cat("Removing NA's.")
    badind <- which(is.na(y))
    y <- y[-badind]
    x <- x[-badind, ]
  }
  list(x = x, y = y)
}

.CalcWts <- function(x, y, k, wt.vec, min.nb, lambda, margin, mu, output = F){ #, eql.wt.vec){
  
  constraint.vars <- .ConstraintVars(x, y, k, wt.vec)
  order.mat <- constraint.vars$order.mat
  target.mat <- constraint.vars$target.mat
  num.equations <- constraint.vars$num.equations
  num.target.total <- constraint.vars$num.target.total
  num.imposter <- constraint.vars$num.imposter
  
  if(num.equations > 0){
    ###############################
    constraints <- .CreateConstraints(x, target.mat, order.mat, min.nb, margin,
                                     num.equations, num.target.total,
                                     num.imposter)
    A <- constraints$A
    b0 <- constraints$b0
    cvec <- constraints$cvec
    ###################################################################
    #Redefining constraints for equal weights on across moments
#     uniq.wts <- unique(eql.wt.vec)
    
#     if(length(uniq.wts) != ncol(x)){
#       uniq.wt.vars <- UnequalWtVars()
#     }
    
    # Adding dummy variable for each equation.
    Adiag <- matrix(0, nrow = num.equations, ncol = num.equations)
    diag(Adiag) <- 1
    A <- cbind(A, Adiag)
    
    # Adding dummy opt variables
    cvec <- scale(cvec, center = FALSE)
    cvec <- c(mu * cvec, rep((1 - mu), num.equations))
    # Sparsity induction
    cvec[1:ncol(x)] <- cvec[1:ncol(x)] + lambda
    
    lin.prog <- solveLP(cvec = cvec, bvec = b0, Amat = A,
                       const.dir = rep(">=", length(b0)), verbose = 1,
                       maxiter = 2000, lpSolve = TRUE)
    ###################################################################
    #Equal weights
    errors <- lin.prog$solution[(ncol(x) + 1):length(lin.prog$solution)]
    weights <- lin.prog$solution[1:ncol(x)]
    ###################################################################
    #Producing output
    if(is.null(weights) || sum(weights) == 0 || sum(is.na(weights)) > 0){
      if(output == T){
        cat("All weights are zero or NA! \n")
      }
      wts <- rep(1, ncol(x))
      errors <- NaN
    }else{
      wts <-  weights / sum(weights) * ncol(x)
      # Under development.
#       wts <- UnequalWtAdj
    }
  }else{
    if(output == TRUE){
      cat("No imposters!\n")
    }
    errors <- 0
    wts <- wt.vec
  }
  list(wts = wts, errors = errors)
}

