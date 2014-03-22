lmnn
====

This is an R package that implements the large-margin nearest neighbor algorithm
as found in Weinberger, 2006. Currently, the implementation only computes a
linear transformation matrix, L, that is diagonal.  This results in a linear
programming problem that minimizes/maximizes the L1 distance between
target/imposter neighbors. An additional function is provided that repeatedly implements the LMNN algorithm due to the fact that target neighbors can change depending on the calculated
transformation.
