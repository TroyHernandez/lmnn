\name{L1NN.cv}
\alias{L1NN.cv}
\title{Runs cross validation on the L1NN algorithm to find optimal parameters}
\usage{
L1NN.cv(x, y, num.cv = 10, seed = 1, n.fold = 10, mu = runif(num.cv, 0,
  1), margin = 10^runif(num.cv, -3, 3), lambda = 10^runif(num.cv, -3, 0),
  k = 1, min.nb = sample(c(0, 1, 1), num.cv, replace = T),
  equilibrium = TRUE)
}
\arguments{
  \item{x}{a matrix of numeric variables}

  \item{y}{a character vector of class labels}

  \item{num.cv}{integer indicated the number of random cv
  parameter sets to try}

  \item{seed}{for reproducibility}

  \item{nfold}{integer for the number of folds to used to
  compute cv errors}

  \item{mu}{determines what share of minimization vector is
  on the weights and what share is on the errors.}

  \item{margin}{is size of margin}

  \item{lambda}{induces additional sparsity/variable
  selection}

  \item{k}{number of neighbors to compare}

  \item{min.nb}{logical, determines whether target neighbor
  distance is minimized}

  \item{equilibrium}{should the L1NNEquilibrium function be
  used to compute final reweighting}
}
\description{
This function uses random cross validation (as opposed to
grid) to find the optimal hyper-parameters for the L1NN
function.
}
\keyword{learning}
\keyword{metric}

