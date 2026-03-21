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
  if `Z` is provided.

- Z:

  Optional sparse indicator matrix from
  [`leaf_node_kernel_Z`](https://jaredhuling.github.io/forestBalance/reference/leaf_node_kernel_Z.md)
  such that \\K = Z Z^\top / B\\. When supplied, the solver can avoid
  forming the full kernel matrix.

- leaf_matrix:

  Optional integer matrix of leaf node assignments (observations x
  trees), as returned by
  [`get_leaf_node_matrix`](https://jaredhuling.github.io/forestBalance/reference/get_leaf_node_matrix.md).
  Used by the Block Jacobi preconditioner (`solver = "bj"`) to partition
  observations into leaf groups. If `NULL` and `solver = "bj"`, falls
  back to `"cg"`.

- num.trees:

  Number of trees \\B\\. Required when `Z` is provided.

- solver:

  Which linear solver to use. `"auto"` (default) selects the fastest
  available solver: `"bj"` (Block Jacobi preconditioned CG) when
  `leaf_matrix` is available and \\n \> 5000\\, `"cg"` (plain CG) when
  only `Z` is available, or `"direct"` (sparse Cholesky) for small
  problems. See Details.

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
whenever \\A_i \neq A_j\\. Both solvers exploit this structure by
working on the treated and control blocks independently.

The **Block Jacobi** solver (`"bj"`) uses the first tree's leaf
partition to define a block-diagonal preconditioner for CG. Each leaf
block is a small dense system (~`min.node.size` x `min.node.size`) that
is cheap to factor. This typically reduces CG iterations by 5–10x,
giving a ~20x overall speedup at large \\n\\.

The **CG** solver uses the factored representation \\K = Z Z^\top / B\\
to perform matrix–vector products without forming any kernel matrix.

The **direct** solver extracts sub-blocks and solves via sparse
Cholesky.

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
