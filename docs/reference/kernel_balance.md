# Kernel energy balancing weights via closed-form solution

Computes balancing weights that minimize a kernelized energy distance
between the weighted treated and control distributions and the overall
sample. The weights are obtained via a closed-form solution to a linear
system derived from the kernel energy distance objective.

## Usage

``` r
kernel_balance(
  trt,
  kern = NULL,
  Z = NULL,
  num.trees = NULL,
  solver = c("auto", "direct", "cg"),
  tol = 1e-08,
  maxiter = 2000L
)
```

## Arguments

- trt:

  A binary (0/1) integer or numeric vector indicating treatment
  assignment (`1` = treated, `0` = control).

- kern:

  A symmetric \\n \times n\\ kernel matrix (dense or sparse), or `NULL`
  if `Z` is provided.

- Z:

  Optional sparse indicator matrix from
  [`leaf_node_kernel_Z`](https://jaredhuling.github.io/forestBalance/reference/leaf_node_kernel_Z.md)
  such that \\K = Z Z^\top / B\\. When supplied, the solver can avoid
  forming the full kernel matrix. If both `kern` and `Z` are given, `Z`
  takes priority when the CG solver is selected.

- num.trees:

  Number of trees \\B\\. Required when `Z` is provided.

- solver:

  Which linear solver to use. `"auto"` (default) selects `"direct"` for
  \\n \le 5000\\ and `"cg"` for \\n \> 5000\\. `"direct"` uses sparse
  Cholesky on the treated and control sub-blocks of the kernel. `"cg"`
  uses conjugate gradient iterations with the factored \\Z\\
  representation, avoiding formation of any kernel matrix.

- tol:

  Convergence tolerance for the CG solver. Default is `1e-8`. Ignored
  when `solver = "direct"`.

- maxiter:

  Maximum CG iterations. Default is 1000.

## Value

A list with the following elements:

- weights:

  A numeric vector of length \\n\\ containing the balancing weights.
  Treated weights sum to \\n_1\\ and control weights sum to \\n_0\\.

- solver:

  The solver that was used (`"direct"` or `"cg"`).

## Details

The modified kernel \\K_q\\ used in the optimization is block-diagonal:
the treated–control cross-blocks are zero because \\K_q(i,j) = 0\\
whenever \\A_i \neq A_j\\. Both solvers exploit this structure by
working on the treated and control blocks independently.

The **direct** solver extracts the sub-blocks \\K\_{tt}\\ and
\\K\_{cc}\\ and solves via sparse Cholesky. This gives exact solutions
but requires forming (at least sub-blocks of) the kernel matrix.

The **CG** solver uses the factored representation \\K = Z Z^\top / B\\
to perform matrix–vector products without forming any kernel matrix, via
\\K v = Z (Z^\top v) / B\\. This is much faster and more
memory-efficient at large \\n\\ (e.g., \\n \> 5000\\). The CG iterates
converge to the exact solution; the default tolerance of `5e-11` yields
weight vectors that agree with the direct solution to several decimal
places.

## References

De, S. and Huling, J.D. (2025). Data adaptive covariate balancing for
causal effect estimation for high dimensional data. *arXiv preprint
arXiv:2512.18069*.

## Examples

``` r
# \donttest{
library(grf)
n <- 200
p <- 5
X <- matrix(rnorm(n * p), n, p)
A <- rbinom(n, 1, plogis(X[, 1]))
Y <- X[, 1] + rnorm(n)

forest <- multi_regression_forest(X, cbind(A, Y), num.trees = 500)
K <- forest_kernel(forest)
bal <- kernel_balance(A, K)

# Weighted ATE estimate
w <- bal$weights
ate <- weighted.mean(Y[A == 1], w[A == 1]) -
       weighted.mean(Y[A == 0], w[A == 0])
# }
```
