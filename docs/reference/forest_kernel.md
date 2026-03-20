# Compute random forest proximity kernel from a GRF forest

Extracts the leaf node assignments for each observation across all trees
in a trained GRF forest, then computes the \\n \times n\\ proximity
kernel matrix where entry \\(i,j)\\ is the proportion of trees in which
observations \\i\\ and \\j\\ share a leaf node.

## Usage

``` r
forest_kernel(forest, newdata = NULL)
```

## Arguments

- forest:

  A trained forest object from the grf package.

- newdata:

  A numeric matrix of observations. If `NULL` (default), uses the
  original training data.

## Value

A symmetric numeric matrix of dimension \\n \times n\\.

## Details

This is a convenience function that calls
[`get_leaf_node_matrix`](https://jaredhuling.github.io/forestBalance/reference/get_leaf_node_matrix.md)
followed by
[`leaf_node_kernel`](https://jaredhuling.github.io/forestBalance/reference/leaf_node_kernel.md).
If you need both the leaf matrix and the kernel, it is more efficient to
call them separately.

## Examples

``` r
# \donttest{
library(grf)
n <- 100
p <- 5
X <- matrix(rnorm(n * p), n, p)
Y <- cbind(X[, 1] + rnorm(n), X[, 2] + rnorm(n))
forest <- multi_regression_forest(X, Y, num.trees = 50)

K <- forest_kernel(forest)
# }
```
