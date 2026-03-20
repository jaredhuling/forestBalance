# Summarize a forest_balance object

Produces a detailed summary of the forest balance fit, including the ATE
estimate, covariate balance diagnostics (SMD, ESS, energy distance), and
kernel sparsity information.

## Usage

``` r
# S3 method for class 'forest_balance'
summary(object, X.trans = NULL, threshold = 0.1, energy.dist = TRUE, ...)

# S3 method for class 'summary.forest_balance'
print(x, ...)
```

## Arguments

- object:

  A `forest_balance` object returned by
  [`forest_balance`](http://jaredhuling.org/forestBalance/reference/forest_balance.md).

- X.trans:

  An optional matrix of nonlinear or transformed covariates (\\n \times
  q\\) on which to additionally assess balance (e.g., interactions,
  squared terms).

- threshold:

  SMD threshold for flagging imbalanced covariates. Default is 0.1.

- energy.dist:

  Logical; if `TRUE` (default), compute the weighted energy distance
  (skipped if \\n \> 5000\\).

- ...:

  Ignored.

- x:

  A `summary.forest_balance` object.

## Value

Invisibly returns a list of class `"summary.forest_balance"` containing:

- ate:

  The estimated ATE.

- balance_weighted:

  A `forest_balance_diag` object for the weighted sample.

- balance_unweighted:

  A `forest_balance_diag` object for the unweighted sample.

- kernel_sparsity:

  Fraction of nonzero entries in the kernel matrix.

- n:

  Total sample size.

- n1:

  Number of treated.

- n0:

  Number of control.

- num.trees:

  Number of trees in the forest.

## Examples

``` r
# \donttest{
n <- 500; p <- 10
X <- matrix(rnorm(n * p), n, p)
A <- rbinom(n, 1, plogis(0.5 * X[, 1]))
Y <- X[, 1] + rnorm(n)

fit <- forest_balance(X, A, Y)
summary(fit)
#> Forest Kernel Energy Balancing
#> ============================================================ 
#>   n = 500  (n_treated = 267, n_control = 233)
#>   Trees: 1000
#>   Solver: CG (kernel not materialized)
#> 
#>   ATE estimate: 0.0391
#> ============================================================ 
#> 
#> Covariate Balance (|SMD|)
#> ------------------------------------------------------------ 
#>   Covariate     Unweighted      Weighted
#>   ----------  ------------  ------------
#>   X1              0.5147 *        0.0259
#>   X2                0.0254        0.0088
#>   X3                0.0334        0.0136
#>   X4                0.0265        0.0048
#>   X5              0.1943 *        0.0310
#>   X6                0.0621        0.0282
#>   X7                0.0906        0.0088
#>   X8              0.1452 *        0.0322
#>   X9                0.0517        0.0195
#>   X10               0.0038        0.0053
#>   ----------  ------------  ------------
#>   Max |SMD|         0.5147        0.0322
#>   (* indicates |SMD| > 0.10)
#> 
#> Effective Sample Size
#> ------------------------------------------------------------ 
#>   Treated: 207 / 267  (77%)
#>   Control: 178 / 233  (76%)
#> 
#> Energy Distance
#> ------------------------------------------------------------ 
#>   Unweighted: 0.0941
#>   Weighted:   0.0214
#> ============================================================ 

# With nonlinear balance assessment
X.nl <- cbind(X[,1]^2, X[,1]*X[,2])
colnames(X.nl) <- c("X1^2", "X1*X2")
summary(fit, X.trans = X.nl)
#> Forest Kernel Energy Balancing
#> ============================================================ 
#>   n = 500  (n_treated = 267, n_control = 233)
#>   Trees: 1000
#>   Solver: CG (kernel not materialized)
#> 
#>   ATE estimate: 0.0391
#> ============================================================ 
#> 
#> Covariate Balance (|SMD|)
#> ------------------------------------------------------------ 
#>   Covariate     Unweighted      Weighted
#>   ----------  ------------  ------------
#>   X1              0.5147 *        0.0259
#>   X2                0.0254        0.0088
#>   X3                0.0334        0.0136
#>   X4                0.0265        0.0048
#>   X5              0.1943 *        0.0310
#>   X6                0.0621        0.0282
#>   X7                0.0906        0.0088
#>   X8              0.1452 *        0.0322
#>   X9                0.0517        0.0195
#>   X10               0.0038        0.0053
#>   ----------  ------------  ------------
#>   Max |SMD|         0.5147        0.0322
#>   (* indicates |SMD| > 0.10)
#> 
#> Transformed Covariate Balance (|SMD|)
#> ------------------------------------------------------------ 
#>   Transform     Unweighted      Weighted
#>   ----------  ------------  ------------
#>   X1^2              0.0438        0.0478
#>   X1*X2             0.0648        0.0431
#>   ----------  ------------  ------------
#>   Max |SMD|         0.0648        0.0478
#> 
#> Effective Sample Size
#> ------------------------------------------------------------ 
#>   Treated: 207 / 267  (77%)
#>   Control: 178 / 233  (76%)
#> 
#> Energy Distance
#> ------------------------------------------------------------ 
#>   Unweighted: 0.0941
#>   Weighted:   0.0214
#> ============================================================ 
# }
```
