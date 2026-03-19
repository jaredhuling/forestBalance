# Estimate ATE using forest-based kernel energy balancing

Fits a multivariate random forest that jointly models the relationship
between covariates, treatment, and outcome, computes a random forest
proximity kernel, and then uses kernel energy balancing to produce
weights for estimating the average treatment effect (ATE).

## Usage

``` r
forest_balance(
  X,
  A,
  Y,
  num.trees = 1000,
  min.node.size = 10,
  scale.outcomes = TRUE,
  solver = c("auto", "direct", "cg"),
  tol = 5e-11,
  ...
)
```

## Arguments

- X:

  A numeric matrix or data frame of covariates (\\n \times p\\).

- A:

  A binary (0/1) vector of treatment assignments.

- Y:

  A numeric vector of outcomes.

- num.trees:

  Number of trees to grow in the forest. Default is 1000.

- min.node.size:

  Minimum number of observations per leaf node. Default is 10.

- scale.outcomes:

  If `TRUE` (default), the joint outcome matrix `cbind(A, Y)` is
  column-standardized before fitting the forest. This ensures that
  treatment and outcome contribute equally to the splits.

- solver:

  Which linear solver to use for the balancing weights. `"auto"`
  (default) selects `"direct"` for \\n \le 5000\\ and `"cg"` for \\n \>
  5000\\. See
  [`kernel_balance`](http://jaredhuling.org/forestBalance/reference/kernel_balance.md)
  for details.

- tol:

  Convergence tolerance for the CG solver. Default is `5e-11`.

- ...:

  Additional arguments passed to
  [`multi_regression_forest`](https://rdrr.io/pkg/grf/man/multi_regression_forest.html).

## Value

An object of class `"forest_balance"` (a list) with the following
elements:

- ate:

  The estimated average treatment effect (Hajek estimator).

- weights:

  The balancing weights from
  [`kernel_balance`](http://jaredhuling.org/forestBalance/reference/kernel_balance.md).

- kernel:

  The \\n \times n\\ forest proximity kernel (sparse matrix), or `NULL`
  when the CG solver is used (since the kernel is never formed).

- forest:

  The trained
  [`multi_regression_forest`](https://rdrr.io/pkg/grf/man/multi_regression_forest.html)
  object.

- X, A, Y:

  The input data (stored for use by
  [`summary.forest_balance`](http://jaredhuling.org/forestBalance/reference/summary.forest_balance.md)).

- n, n1, n0:

  Total, treated, and control sample sizes.

- solver:

  The solver that was used (`"direct"` or `"cg"`).

The object has `print` and `summary` methods. Use
[`summary.forest_balance`](http://jaredhuling.org/forestBalance/reference/summary.forest_balance.md)
for covariate balance diagnostics.

## Details

The method proceeds in three steps:

1.  A
    [`multi_regression_forest`](https://rdrr.io/pkg/grf/man/multi_regression_forest.html)
    is fit on covariates `X` with a bivariate response `(A, Y)`. This
    jointly models the relationship between covariates, treatment
    assignment, and outcome.

2.  The forest's leaf co-membership structure defines a proximity
    kernel: \\K(i,j)\\ is the proportion of trees where \\i\\ and \\j\\
    share a leaf. Because the forest splits on both \\A\\ and \\Y\\,
    this kernel captures confounding structure.

3.  [`kernel_balance`](http://jaredhuling.org/forestBalance/reference/kernel_balance.md)
    computes balancing weights via the closed-form kernel energy
    distance solution. The ATE is then estimated using the Hajek (ratio)
    estimator with these weights.

For \\n \le 5000\\ (or when `solver = "direct"`), the kernel matrix is
formed explicitly and the balancing system is solved via sparse
Cholesky. For \\n \> 5000\\ (or when `solver = "cg"`), the kernel is
never formed; instead, a conjugate gradient solver operates on the
factored representation \\K = Z Z^\top / B\\, saving both time and
memory.

## References

De, S. and Huling, J.D. (2025). Data adaptive covariate balancing for
causal effect estimation for high dimensional data. *arXiv preprint
arXiv:2512.18069*.

## Examples

``` r
# \donttest{
n <- 500
p <- 10
X <- matrix(rnorm(n * p), n, p)
A <- rbinom(n, 1, plogis(0.5 * X[, 1]))
Y <- X[, 1] + rnorm(n)  # true ATE = 0

result <- forest_balance(X, A, Y)
result$ate
#> [1] 0.1645587
# }
```
