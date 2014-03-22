# L1NN_3.R

.ConstraintVars <- function(x, y, k, wt.vec){
  #Find k neighbors for each sample
  dist.mat <- as.matrix(dist(t(t(x) * wt.vec), method = "manhattan",
                             upper = TRUE, diag = TRUE))
  #matrix each row a sample of 'k' closest indices
  order.mat <- as.matrix(t(apply(dist.mat, 1, order))[, 2:(k + 1)],
                         ncol = k)
  #matrix telling if each neighbor is a target
  target.mat <- matrix(y[order.mat] == rep(y, k)
                       , ncol = k)
  
  #vector telling how many imposters are around
  num.imposter <- k - apply(target.mat, 1, sum)
  num.equations <- sum(num.imposter) #*apply(target.mat,1,sum))
  num.target <- apply(target.mat,1,sum)
  num.target.total <- sum(num.target)

  list(order.mat = order.mat, target.mat = target.mat,
       num.imposter = num.imposter, num.equations = num.equations,
       num.target.total = num.target.total)
}

###############################################################

.CreateConstraints <- function(x, target.mat, order.mat, min.nb, margin,
                              num.equations, num.target.total, num.imposter){
  #b0 length depends on the number of mismatches
  b0 <- rep(margin, num.equations)
  
  ###################################################################
  #Minimization matrix, to be flattened into minimization vector
  cvec <- matrix(0, nrow = num.target.total, ncol = ncol(x))
  
  #This is the matrix that measures imposter penetrations
  A <- matrix(0,nrow=num.equations,ncol=ncol(x))
  
  temp.ind <- 1
  temp.ind2 <- 1
  
  for(i in 1:nrow(target.mat)){
    #imposter and target locations
    imposter.loc <- order.mat[i,which(target.mat[i,]!=T)]
    target.loc <- order.mat[i,which(target.mat[i,]==T)]
    
    #If imposters exist, we use them to create A matrix.
    if(num.imposter[i] > 0){
      #Iterate through each pair of targets and i
      for(j in 1:length(imposter.loc)){
        #Creating constraints
        A[temp.ind, ] <- abs(x[i, ] - x[imposter.loc[j], ])
        temp.ind <- temp.ind + 1
      }
    }
    
    #If good neighbors exist, we use them to create cvec.
    if(length(target.loc) > 0){
      for(j in 1:length(target.loc)){
        cvec[temp.ind2, 1:ncol(x)] <- abs(x[i, ] - x[target.loc[j], ])
        temp.ind2 <- temp.ind2 + 1
      }
    }
  }
  
  if(min.nb == T){
    cvec <- apply(cvec, 2, sum)
  }else{
    cvec <- rep(1, ncol(x))
  }
  
  list(A = A, b0 = b0, cvec = cvec)
}


###########################################
#Underconstruction!!!
# UnequalWtVars <- function(A, cvec, eql.wt.vec, uniq.wts, ){
#   
#   #First iteration of grouping weights together
#   temp.ind <- which(eql.wt.vec == uniq.wts[1])
#   if(length(temp.ind) > 1){
#     A.temp[, temp.ind] <- apply(A[, temp.ind], 1, mean)
#     cvec.temp[temp.ind] <- mean(cvec[temp.ind])
#   }else{
#     A.temp <- t(t(A[,temp.ind]))
#     cvec.temp <- cvec[temp.ind]
#   }
#   
#   if(length(uniq.wts)>2){
#     for(i in 2:length(uniq.wts)){
#       temp.ind <- which(eql.wt.vec==uniq.wts[i])
#       if(length(temp.ind)>1){
#         A.temp <- cbind(A.temp,apply(A[,temp.ind],1,sum))
#         cvec.temp <- c(cvec.temp,sum(cvec[temp.ind]))
#       }else{
#         A.temp <- cbind(A.temp,A[,temp.ind])
#         cvec.temp <- c(cvec.temp,cvec[temp.ind])
#       }
#     }
#   }
#   
#   list(A = A.temp, cvec = cvec.temp)
# }

#############################################
# Under construction
# UnequalWtAdj <- function(){
#   if(length(uniq.wts) != ncol(x)){
#     wtstemp <- rep(0, ncol(x))
#     numind <- length(which(eql.wt.vec == uniq.wts[1]))
#     wtstemp[which(eql.wt.vec == uniq.wts[1])] <-
#       rep(wts[1] / numind * ncol(x), numind)
#     
#     if(length(uniq.wts)>2){
#       for(i in 2:length(uniq.wts)){
#         numind <- length(which(eql.wt.vec==uniq.wts[i]))
#         wtstemp[which(eql.wt.vec==uniq.wts[i])] <- rep(wts[i]/numind*ncol(x),numind)
#       }
#     }
#     wts <- wtstemp
#   }
#   wts
# }