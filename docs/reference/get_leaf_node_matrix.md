# Extract leaf node membership matrix from a GRF forest

For each observation and each tree in a trained GRF forest, determines
which leaf node the observation falls into. The tree traversal is
implemented in C++ for speed, directly reading the forest's internal
tree structure and avoiding the overhead of
[`get_tree`](https://rdrr.io/pkg/grf/man/get_tree.html) and
[`get_leaf_node`](https://rdrr.io/pkg/grf/man/get_leaf_node.html).

## Usage

``` r
get_leaf_node_matrix(forest, newdata = NULL)
```

## Arguments

- forest:

  A trained forest object from the grf package (e.g., from
  [`multi_regression_forest`](https://rdrr.io/pkg/grf/man/multi_regression_forest.html),
  [`regression_forest`](https://rdrr.io/pkg/grf/man/regression_forest.html),
  etc.).

- newdata:

  A numeric matrix of observations to predict leaf membership for. Must
  have the same number of columns as the training data. If `NULL`
  (default), uses the original training data stored in the forest.

## Value

An integer matrix of dimension `nrow(newdata)` by `num.trees`, where
entry `[i, b]` is the internal node index (1-based) of the leaf that
observation `i` falls into in tree `b`.

## Details

This function reads the internal tree vectors stored in a GRF forest
object (`_child_nodes`, `_split_vars`, `_split_values`, `_root_nodes`,
`_send_missing_left`) and traverses each tree in compiled C++ code. This
is dramatically faster than the per-observation R-level loop in
[`grf::get_leaf_node`](https://rdrr.io/pkg/grf/man/get_leaf_node.html).

## Examples

``` r
# \donttest{
library(grf)
n <- 200
p <- 5
X <- matrix(rnorm(n * p), n, p)
Y <- cbind(X[, 1] + rnorm(n), X[, 2] + rnorm(n))
forest <- multi_regression_forest(X, Y, num.trees = 100)

# Leaf membership for training data
leaf_mat <- get_leaf_node_matrix(forest)

# Leaf membership for new data
X.test <- matrix(rnorm(50 * p), 50, p)
leaf_mat_test <- get_leaf_node_matrix(forest, X.test)
# }
```
