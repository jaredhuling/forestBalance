# Compute random forest proximity kernel from a leaf node matrix

Given a matrix of leaf node assignments (observations x trees), computes
the \\n \times n\\ kernel matrix where entry \\(i,j)\\ is the proportion
of trees in which observations \\i\\ and \\j\\ fall into the same leaf
node.

## Usage

``` r
leaf_node_kernel(leaf_matrix, sparse = TRUE)
```

## Arguments

- leaf_matrix:

  An integer matrix of dimension \\n \times B\\, where
  `leaf_matrix[i, b]` is the leaf node ID for observation `i` in tree
  `b`. Typically produced by
  [`get_leaf_node_matrix`](http://jaredhuling.org/forestBalance/reference/get_leaf_node_matrix.md).

- sparse:

  Logical; if `TRUE` (default), returns a sparse `dgCMatrix` from the
  Matrix package. If `FALSE`, returns a dense base R matrix.

## Value

A symmetric \\n \times n\\ matrix (sparse or dense depending on
`sparse`) where entry \\(i,j)\\ equals \$\$K(i,j) = \frac{1}{B}
\sum\_{b=1}^{B} \mathbf{1}(L_b(i) = L_b(j)),\$\$ i.e., the fraction of
trees where \\i\\ and \\j\\ share a leaf.

## Details

For each tree \\b\\, the leaf assignment defines a sparse \\n \times
L_b\\ indicator matrix \\Z_b\\. This function stacks all \\B\\ indicator
matrices column-wise into a single sparse matrix \\Z = \[Z_1 \| Z_2 \|
\cdots \| Z_B\]\\ of dimension \\n \times \sum_b L_b\\. The kernel is
then obtained via a single sparse cross-product: \\K = Z Z^\top / B\\.
Leaf ID remapping is done in C++ for speed.

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
K <- leaf_node_kernel(leaf_mat)               # sparse
K_dense <- leaf_node_kernel(leaf_mat, sparse = FALSE)  # dense
# }
```
