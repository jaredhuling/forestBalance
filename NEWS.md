# forestBalance 0.1.1

## New features

* `forest_balance()` and `kernel_balance()` gain `estimand = c("ATE", "ATT", "ATC")`.

* Optional ridge penalty `lambda` (default 0) in all solver paths.

* `kernel.response = c("joint", "treatment", "outcome")` selects the responses used to
  grow the kernel forest.

* `crossfit.balance = c("fold", "full")` chooses whether cross-fitted weights are solved on
  the evaluation fold or on the full sample with the held-out forest.

* New `maxiter` argument for the conjugate gradient solver. Fits now return leaf counts,
  solver iteration counts, and diagnostics for every linear solve (residuals, arm totals,
  and projected KKT residuals).

## Changes in behavior

* The default ridge penalty is now `lambda = 0.1` in `forest_balance()` and
  `kernel_balance()`. Without a ridge, the kernel block of an arm is singular when units
  share leaves in nearly every tree (for example, with a binary outcome and many
  covariates), and the balancing weights are then not determined.

* A numerically singular kernel block now raises an informative error instead of
  returning arbitrary weights or failing with a generic message, and an ill-conditioned
  block (condition number above about 10^6) gives a warning that the weights can be
  extreme. Because kernel entries lie in [0, 1], the condition number of an arm block of
  size `n_a` is at most `1 + n_a / lambda`, so with the default ridge the warning requires
  an arm of more than about 100,000 units.

* The target linear solve now uses the right-hand side `n_a^2 b_a` directly instead of
  rescaling its solution by `n_a^2`, so solver error is no longer amplified by the squared
  arm size. A small positive `lambda` with the iterative solvers at large `n` no longer
  fails the final accuracy check. Weights are mathematically unchanged.

* The block Jacobi solver (`solver = "bj"`) now stops on the same absolute residual
  criterion as the conjugate gradient solver. Its relative criterion left errors that
  grew with the target right-hand side and failed the final accuracy check under a
  positive `lambda` at large `n`.

* Failed linear solves or failed weight accuracy checks now raise an error carrying the
  diagnostics; version 0.1.0 returned weights without these checks. Cross-fitted
  estimates no longer drop failed folds when averaging.

## Compatibility

* `src/cg_solve.cpp` uses `Map<SparseMatrix>` for compatibility with Eigen 5 (thanks to
  Dirk Eddelbuettel, #1).

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
