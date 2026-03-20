# forestBalance: Forest Kernel Energy Balancing for Causal Inference

Estimates average treatment effects (ATE) using kernel energy balancing
with random forest similarity kernels. A multivariate random forest
jointly models covariates, treatment, and outcome to build a proximity
kernel that captures confounding structure. Balancing weights are
obtained via a closed-form kernel energy distance solution.

## Main function

[`forest_balance`](https://jaredhuling.github.io/forestBalance/reference/forest_balance.md)
is the primary interface. It fits the forest, constructs the kernel,
computes balancing weights, and estimates the ATE. By default it uses
K-fold cross-fitting and an adaptive leaf size to minimize overfitting
bias.

## Key features

- Adaptive `min.node.size` that scales with \\n\\ and \\p\\

- K-fold cross-fitting to reduce kernel overfitting bias

- Rcpp-accelerated leaf node extraction

- Sparse kernel construction via single `tcrossprod`

- Conjugate gradient solver for large \\n\\ (avoids forming the kernel
  matrix entirely)

## Lower-level interface

For more control, the pipeline can be run step by step:
[`get_leaf_node_matrix`](https://jaredhuling.github.io/forestBalance/reference/get_leaf_node_matrix.md),
[`leaf_node_kernel`](https://jaredhuling.github.io/forestBalance/reference/leaf_node_kernel.md),
[`kernel_balance`](https://jaredhuling.github.io/forestBalance/reference/kernel_balance.md).

## References

De, S. and Huling, J.D. (2025). Data adaptive covariate balancing for
causal effect estimation for high dimensional data. *arXiv preprint
arXiv:2512.18069*.

## See also

Useful links:

- <https://github.com/jaredhuling/forestBalance>

- Report bugs at <https://github.com/jaredhuling/forestBalance/issues>

## Author

**Maintainer**: Jared Huling <jaredhuling@gmail.com>

Authors:

- Simion De
