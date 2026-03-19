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
#>   n = 500  (n_treated = 240, n_control = 260)
#>   Trees: 1000
#>   Kernel density: 42.2% nonzero
#> 
#>   ATE estimate: 0.1775
#> ============================================================ 
#> 
#> Covariate Balance (|SMD|)
#> ------------------------------------------------------------ 
#>   Covariate     Unweighted      Weighted
#>   ----------  ------------  ------------
#>   X1              0.5187 *        0.0446
#>   X2                0.0150        0.0314
#>   X3                0.0680        0.0021
#>   X4                0.0512        0.0108
#>   X5                0.0002        0.0021
#>   X6                0.0557        0.0342
#>   X7                0.0189        0.0382
#>   X8                0.0667        0.0231
#>   X9                0.0604        0.0171
#>   X10               0.0015        0.0047
#>   ----------  ------------  ------------
#>   Max |SMD|         0.5187        0.0446
#>   (* indicates |SMD| > 0.10)
#> 
#> Effective Sample Size
#> ------------------------------------------------------------ 
#>   Treated: 175 / 240  (73%)
#>   Control: 199 / 260  (77%)
#> 
#> Energy Distance
#> ------------------------------------------------------------ 
#>   Unweighted: 0.0774
#>   Weighted:   0.0187
#> ============================================================ 

# With nonlinear balance assessment
X.nl <- cbind(X[,1]^2, X[,1]*X[,2])
colnames(X.nl) <- c("X1^2", "X1*X2")
summary(fit, X.trans = X.nl)
#> Forest Kernel Energy Balancing
#> ============================================================ 
#>   n = 500  (n_treated = 240, n_control = 260)
#>   Trees: 1000
#>   Kernel density: 42.2% nonzero
#> 
#>   ATE estimate: 0.1775
#> ============================================================ 
#> 
#> Covariate Balance (|SMD|)
#> ------------------------------------------------------------ 
#>   Covariate     Unweighted      Weighted
#>   ----------  ------------  ------------
#>   X1              0.5187 *        0.0446
#>   X2                0.0150        0.0314
#>   X3                0.0680        0.0021
#>   X4                0.0512        0.0108
#>   X5                0.0002        0.0021
#>   X6                0.0557        0.0342
#>   X7                0.0189        0.0382
#>   X8                0.0667        0.0231
#>   X9                0.0604        0.0171
#>   X10               0.0015        0.0047
#>   ----------  ------------  ------------
#>   Max |SMD|         0.5187        0.0446
#>   (* indicates |SMD| > 0.10)
#> 
#> Transformed Covariate Balance (|SMD|)
#> ------------------------------------------------------------ 
#>   Transform     Unweighted      Weighted
#>   ----------  ------------  ------------
#>   X1^2            0.1326 *        0.0154
#>   X1*X2             0.0120        0.0103
#>   ----------  ------------  ------------
#>   Max |SMD|         0.1326        0.0154
#> 
#> Effective Sample Size
#> ------------------------------------------------------------ 
#>   Treated: 175 / 240  (73%)
#>   Control: 199 / 260  (77%)
#> 
#> Energy Distance
#> ------------------------------------------------------------ 
#>   Unweighted: 0.0774
#>   Weighted:   0.0187
#> ============================================================ 
# }
```
