# Build the sparse indicator matrix Z from a leaf node matrix

Returns the sparse \\n \times L\\ indicator matrix \\Z\\ such that the
proximity kernel is \\K = Z Z^\top / B\\. This factored representation
can be passed to
[`kernel_balance`](https://jaredhuling.github.io/forestBalance/reference/kernel_balance.md)
to enable the CG solver, which avoids forming the full \\n \times n\\
kernel.

## Usage

``` r
leaf_node_kernel_Z(leaf_matrix)
```

## Arguments

- leaf_matrix:

  An integer matrix of dimension \\n \times B\\, where
  `leaf_matrix[i, b]` is the leaf node ID for observation `i` in tree
  `b`. Typically produced by
  [`get_leaf_node_matrix`](https://jaredhuling.github.io/forestBalance/reference/get_leaf_node_matrix.md).

## Value

A sparse `dgCMatrix` of dimension \\n \times L\\, where \\L = \sum_b
L_b\\ is the total number of leaves across all trees. Each row has
exactly \\B\\ nonzero entries (one per tree).

## Examples

``` r
# \donttest{
library(grf)
n <- 100
p <- 5
X <- matrix(rnorm(n * p), n, p)
Y <- cbind(X[, 1] + rnorm(n), X[, 2] + rnorm(n))
forest <- multi_regression_forest(X, Y, num.trees = 50)
leaf_mat <- get_leaf_node_matrix(forest)
Z <- leaf_node_kernel_Z(leaf_mat)
# }
```
