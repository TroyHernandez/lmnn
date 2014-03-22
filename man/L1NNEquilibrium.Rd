\name{L1NNEquilibrium}
\alias{L1NNEquilibrium}
\title{Runs the L1NN algorithm repeatedly to find optimal parameters following
possible change of target neighbors}
\usage{
L1NNEquilibrium(x = x, y = y, mu = 0.5, margin = 1, lambda = 0,
  k = 1, min.nb = T, output = F, num.its = 10, take.mean = TRUE,
  cv = TRUE)
}
\arguments{
  \item{x}{a matrix of numeric variables}

  \item{y}{a character vector of class labels}

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

  \item{output}{logical indicating if intermediate output
  should printed}

  \item{num.its}{integer indicating maximum number of times
  to rerun L1NN}

  \item{take.mean}{logical indicating if mean of
  reweightings should be taken or FALSE to use only the
  last computed L1NN}

  \item{cv}{should the L1NN.cv function be used to compute
  final reweighting}
}
\description{
Runs the L1NN algorithm repeatedly to find optimal
parameters following possible change of target neighbors
}
\keyword{learning}
\keyword{metric}

