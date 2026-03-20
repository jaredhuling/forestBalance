# Estimate ATE using forest-based kernel energy balancing

Fits a multivariate random forest that jointly models the relationship
between covariates, treatment, and outcome, computes a random forest
proximity kernel, and then uses kernel energy balancing to produce
weights for estimating the average treatment effect (ATE). By default,
K-fold cross-fitting is used to avoid overfitting bias from estimating
the kernel on the same data used for treatment effect estimation.

## Usage

``` r
forest_balance(
  X,
  A,
  Y,
  num.trees = 1000,
  min.node.size = NULL,
  cross.fitting = TRUE,
  num.folds = 2,
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

  Minimum number of observations per leaf node. If `NULL` (default), an
  adaptive heuristic is used:
  `max(20, min(floor(n/200) + p, floor(n/50)))`. This scales the leaf
  size with both the sample size and the number of covariates, which
  empirically yields low bias. See Details.

- cross.fitting:

  Logical; if `TRUE` (default), use K-fold cross-fitting to construct
  the kernel from held-out data, reducing overfitting bias. If `FALSE`,
  the kernel is estimated on the full sample.

- num.folds:

  Number of cross-fitting folds. Default is 2. Only used when
  `cross.fitting = TRUE`.

- scale.outcomes:

  If `TRUE` (default), the joint outcome matrix `cbind(A, Y)` is
  column-standardized before fitting the forest. This ensures that
  treatment and outcome contribute equally to the splits.

- solver:

  Which linear solver to use for the balancing weights. `"auto"`
  (default) selects `"direct"` for small fold sizes and `"cg"` for large
  fold sizes. See
  [`kernel_balance`](https://jaredhuling.github.io/forestBalance/reference/kernel_balance.md)
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

  The estimated average treatment effect. When cross-fitting is used,
  this is the average of per-fold Hajek estimates (DML1).

- weights:

  The balancing weight vector (length \\n\\). When cross-fitting is
  used, these are the concatenated per-fold weights.

- kernel:

  The \\n \times n\\ forest proximity kernel (sparse matrix), or `NULL`
  when cross-fitting or the CG solver is used.

- forest:

  The trained forest object. When cross-fitting is used, this is the
  last fold's forest.

- X, A, Y:

  The input data.

- n, n1, n0:

  Total, treated, and control sample sizes.

- solver:

  The solver that was used (`"direct"` or `"cg"`).

- crossfit:

  Logical indicating whether cross-fitting was used.

- num.folds:

  Number of folds (if cross-fitting was used).

- fold_ates:

  Per-fold ATE estimates (if cross-fitting was used).

- fold_ids:

  Fold assignments (if cross-fitting was used).

The object has `print` and `summary` methods. Use
[`summary.forest_balance`](https://jaredhuling.github.io/forestBalance/reference/summary.forest_balance.md)
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

3.  [`kernel_balance`](https://jaredhuling.github.io/forestBalance/reference/kernel_balance.md)
    computes balancing weights via the closed-form kernel energy
    distance solution. The ATE is then estimated using the Hajek (ratio)
    estimator with these weights.

**Cross-fitting** (default): For each fold \\k\\, the forest is trained
on all data *except* fold \\k\\, and the kernel for fold \\k\\ is built
from that held-out forest's leaf predictions. This breaks the dependence
between the kernel and the outcomes, reducing overfitting bias. The
final ATE is the average of the per-fold Hajek estimates (DML1).

**Adaptive leaf size**: The default `min.node.size` is set adaptively
via `max(20, min(floor(n/200) + p, floor(n/50)))`. Larger leaves produce
smoother kernels that generalize better, while the cap at `n/50`
prevents kernel degeneracy. This heuristic was calibrated empirically to
minimize RMSE across a range of sample sizes and dimensions.

## References

Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C.,
Newey, W. and Robins, J. (2018). Double/debiased machine learning for
treatment and structural parameters. *The Econometrics Journal*, 21(1),
C1–C68.

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

# Default: cross-fitting with adaptive leaf size
result <- forest_balance(X, A, Y)
result
#> Forest Kernel Energy Balancing (cross-fitted)
#> -------------------------------------------------- 
#>   n = 500  (n_treated = 237, n_control = 263)
#>   Trees: 1000
#>   Cross-fitting: 2 folds
#>   Solver: cg
#>   ATE estimate: 0.0576
#>   Fold ATEs: 0.0499, 0.0652
#>   ESS: treated = 172/237 (73%)   control = 198/263 (75%)
#> -------------------------------------------------- 
#> Use summary() for covariate balance details.

# Without cross-fitting
result_nocf <- forest_balance(X, A, Y, cross.fitting = FALSE)
# }
```
