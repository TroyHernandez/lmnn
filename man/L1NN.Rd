\name{L1NN}
\alias{L1NN}
\title{An optimal L1 Large Margin reweighting function.}
\usage{
L1NN(x = x, y = y, wt.vec = rep(1, ncol(x)), min.len = 2, mu = 0.5,
  margin = 1, lambda = 0, k = 1, min.nb = T, output = F)
}
\arguments{
  \item{x}{a matrix of numeric variables}

  \item{y}{a character vector of class labels}

  \item{eql.wt.vec}{determines which columns should share
  weights}

  \item{mu}{determines what share of minimization vector is
  on the weights and what share is on the errors.}

  \item{margin}{is size of margin}

  \item{lambda}{induces additional sparsity/variable
  selection}

  \item{k}{number of neighbors to compare}

  \item{min.nb}{logical, determines whether target neighbor
  distance is minimized}

  \item{wt.vec}{current weighting vector}
}
\description{
This function computes the optimal weighting of variables
to minimize the L1 distance for k-nearest neighbors
algorithms.
}
\keyword{learning}
\keyword{metric}

