# forestBalance 0.1.0

Initial release.

## Features

* `forest_balance()`: High-level ATE estimation with forest kernel energy
  balancing. Supports K-fold cross-fitting (default) and adaptive
  `min.node.size` selection.

* Adaptive `min.node.size` heuristic: `max(20, min(floor(n/200) + p, floor(n/50)))`.
  Scales leaf size with sample size and dimension to balance bias and variance.

* Two linear solvers with automatic selection:
    - **Direct**: sparse Cholesky on treated/control sub-blocks. Used for
      moderate sample sizes.
    - **Conjugate gradient (CG)**: matrix-free solver using the factored
      kernel representation. Avoids forming the full kernel matrix for large n.

* Rcpp-accelerated leaf node extraction from GRF forests (C++ tree traversal).

* Sparse kernel construction via single `tcrossprod` on the leaf indicator
  matrix.

* `compute_balance()`: Covariate balance diagnostics (SMD, ESS, energy
  distance) with support for nonlinear transformations.

* `simulate_data()`: Synthetic data generation with nonlinear confounding for
  demonstrations and testing.
