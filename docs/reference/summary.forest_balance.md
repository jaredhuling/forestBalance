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
  [`forest_balance`](https://jaredhuling.github.io/forestBalance/reference/forest_balance.md).

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

- threshold:

  The SMD threshold used for flagging imbalanced covariates.

The input `x`, invisibly. Called for its side effect of printing a
detailed balance summary to the console.

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
#>   n = 500  (n_treated = 266, n_control = 234)
#>   Trees: 1000
#>   Solver: CG (kernel not materialized)
#> 
#>   ATE estimate: 0.0025
#> ============================================================ 
#> 
#> Covariate Balance (|SMD|)
#> ------------------------------------------------------------ 
#>   Covariate     Unweighted      Weighted
#>   ----------  ------------  ------------
#>   X1              0.4659 *        0.0155
#>   X2                0.0757        0.0142
#>   X3                0.0844        0.0083
#>   X4              0.1315 *        0.0309
#>   X5                0.0121        0.0325
#>   X6              0.1600 *        0.0590
#>   X7              0.1450 *        0.0631
#>   X8                0.0711        0.0231
#>   X9                0.0233        0.0442
#>   X10               0.0582        0.0017
#>   ----------  ------------  ------------
#>   Max |SMD|         0.4659        0.0631
#>   (* indicates |SMD| > 0.10)
#> 
#> Effective Sample Size
#> ------------------------------------------------------------ 
#>   Treated: 201 / 266  (75%)
#>   Control: 176 / 234  (75%)
#> 
#> Energy Distance
#> ------------------------------------------------------------ 
#>   Unweighted: 0.0833
#>   Weighted:   0.0220
#> ============================================================ 

# With nonlinear balance assessment
X.nl <- cbind(X[,1]^2, X[,1]*X[,2])
colnames(X.nl) <- c("X1^2", "X1*X2")
summary(fit, X.trans = X.nl)
#> Forest Kernel Energy Balancing
#> ============================================================ 
#>   n = 500  (n_treated = 266, n_control = 234)
#>   Trees: 1000
#>   Solver: CG (kernel not materialized)
#> 
#>   ATE estimate: 0.0025
#> ============================================================ 
#> 
#> Covariate Balance (|SMD|)
#> ------------------------------------------------------------ 
#>   Covariate     Unweighted      Weighted
#>   ----------  ------------  ------------
#>   X1              0.4659 *        0.0155
#>   X2                0.0757        0.0142
#>   X3                0.0844        0.0083
#>   X4              0.1315 *        0.0309
#>   X5                0.0121        0.0325
#>   X6              0.1600 *        0.0590
#>   X7              0.1450 *        0.0631
#>   X8                0.0711        0.0231
#>   X9                0.0233        0.0442
#>   X10               0.0582        0.0017
#>   ----------  ------------  ------------
#>   Max |SMD|         0.4659        0.0631
#>   (* indicates |SMD| > 0.10)
#> 
#> Transformed Covariate Balance (|SMD|)
#> ------------------------------------------------------------ 
#>   Transform     Unweighted      Weighted
#>   ----------  ------------  ------------
#>   X1^2              0.0953        0.0186
#>   X1*X2             0.0305        0.0579
#>   ----------  ------------  ------------
#>   Max |SMD|         0.0953        0.0579
#> 
#> Effective Sample Size
#> ------------------------------------------------------------ 
#>   Treated: 201 / 266  (75%)
#>   Control: 176 / 234  (75%)
#> 
#> Energy Distance
#> ------------------------------------------------------------ 
#>   Unweighted: 0.0833
#>   Weighted:   0.0220
#> ============================================================ 
# }
```
