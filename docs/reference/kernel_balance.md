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
  leaf_matrix = NULL,
  num.trees = NULL,
  solver = c("auto", "direct", "cg", "bj"),
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
  if `Z` is provided. Required for `solver = "direct"` (if not provided
  but `Z` is available, the kernel is formed automatically, though this
  is \\O(n^2)\\ and may be slow for large \\n\\).

- Z:

  Optional sparse indicator matrix from
  [`leaf_node_kernel_Z`](https://jaredhuling.github.io/forestBalance/reference/leaf_node_kernel_Z.md)
  such that \\K = Z Z^\top / B\\. When supplied, the iterative solvers
  (`"cg"`, `"bj"`) can perform matrix-free products without forming the
  full kernel. Required for `solver = "cg"` and `solver = "bj"`.

- leaf_matrix:

  Optional integer matrix of leaf node assignments (observations x
  trees), as returned by
  [`get_leaf_node_matrix`](https://jaredhuling.github.io/forestBalance/reference/get_leaf_node_matrix.md).
  Required for `solver = "bj"` (Block Jacobi preconditioner uses tree
  1's leaf partition). If `NULL` and `solver = "bj"`, falls back to
  `"cg"` with a warning.

- num.trees:

  Number of trees \\B\\. Required when `Z` is provided.

- solver:

  Which linear solver to use. `"auto"` (default) selects the best
  available solver based on the inputs: `"cg"` when `Z` is available and
  \\n \> 5000\\, or `"direct"` otherwise. See Details for solver
  requirements.

- tol:

  Convergence tolerance for iterative solvers. Default is `1e-8`.

- maxiter:

  Maximum iterations for iterative solvers. Default is 2000.

## Value

A list with the following elements:

- weights:

  A numeric vector of length \\n\\ containing the balancing weights.
  Treated weights sum to \\n_1\\ and control weights sum to \\n_0\\.

- solver:

  The solver that was used.

## Details

The modified kernel \\K_q\\ used in the optimization is block-diagonal:
the treated–control cross-blocks are zero because \\K_q(i,j) = 0\\
whenever \\A_i \neq A_j\\. All solvers exploit this structure by working
on the treated and control blocks independently.

**Solver requirements:**

|  |  |  |
|----|----|----|
| Solver | Required inputs | Optional inputs |
| `"direct"` | `kern` (or `Z` + `num.trees`) |  |
| `"cg"` | `Z` + `num.trees` |  |
| `"bj"` | `Z` + `num.trees` + `leaf_matrix` | (falls back to `"cg"` if `leaf_matrix` is missing) |

The **direct** solver extracts sub-blocks of the kernel and solves via
sparse Cholesky. If only `Z` is provided, the kernel is formed as \\K =
Z Z^\top / B\\, which requires \\O(n^2)\\ time and memory.

The **CG** solver uses the factored representation \\K = Z Z^\top / B\\
to perform matrix–vector products without forming any kernel matrix.

The **Block Jacobi** solver (`"bj"`) uses the first tree's leaf
partition (from `leaf_matrix`) to define a block-diagonal preconditioner
for CG. Each leaf block is a small dense system that is cheap to factor.

Only 2 linear solves per block are needed (not 3) because the third
right-hand side is a linear combination of the first two.

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
