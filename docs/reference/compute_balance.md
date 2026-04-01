# Compute covariate balance diagnostics for a set of weights

Computes standardized mean differences (SMD), effective sample sizes
(ESS), and optionally the weighted energy distance for a given set of
balancing weights. Can also assess balance on user-supplied nonlinear
transformations of the covariates.

## Usage

``` r
compute_balance(X, trt, weights, X.trans = NULL, energy.dist = TRUE)

# S3 method for class 'forest_balance_diag'
print(x, threshold = 0.1, ...)
```

## Arguments

- X:

  An \\n \times p\\ numeric covariate matrix.

- trt:

  A binary (0/1) vector of treatment assignments of length \\n\\.

- weights:

  A numeric weight vector of length \\n\\. Treated weights should sum to
  \\n_1\\ and control weights to \\n_0\\ (as returned by
  [`kernel_balance`](https://jaredhuling.github.io/forestBalance/reference/kernel_balance.md)
  or
  [`forest_balance`](https://jaredhuling.github.io/forestBalance/reference/forest_balance.md)).

- X.trans:

  An optional matrix of nonlinear or transformed covariates (\\n \times
  q\\) on which to additionally assess balance (e.g., interactions,
  squared terms). If `NULL` (default), only linear covariate balance is
  reported.

- energy.dist:

  Logical; if `TRUE`, compute the weighted energy distance between the
  treated and control groups. This requires computing an \\n \times n\\
  distance matrix and is only feasible for moderate \\n\\ (automatically
  skipped when \\n \> 5000\\). Default is `TRUE`.

- x:

  A `forest_balance_diag` object.

- threshold:

  SMD threshold for flagging imbalanced covariates. Default is 0.1, a
  standard threshold in the causal inference literature.

- ...:

  Ignored.

## Value

An object of class `"forest_balance_diag"` containing:

- smd:

  Named vector of \|SMD\| for each covariate.

- max_smd:

  Maximum \|SMD\| across covariates.

- smd_trans:

  Named vector of \|SMD\| for transformed covariates (if `X.trans` was
  supplied), otherwise `NULL`.

- max_smd_trans:

  Maximum \|SMD\| for transformed covariates, or `NA`.

- energy_dist:

  Weighted energy distance, or `NA` if not computed.

- ess_treated:

  Effective sample size for the treated group as a fraction of \\n_1\\.

- ess_control:

  Effective sample size for the control group as a fraction of \\n_0\\.

- n:

  Total sample size.

- n1:

  Number of treated units.

- n0:

  Number of control units.

The input `x`, invisibly. Called for its side effect of printing balance
diagnostics to the console.

## Details

The standardized mean difference for covariate \\j\\ is defined as
\$\$\text{SMD}\_j = \frac{\|\bar X\_{j,1}^w - \bar
X\_{j,0}^w\|}{s_j},\$\$ where \\\bar X\_{j,a}^w\\ is the weighted mean
of covariate \\j\\ in group \\a\\ and \\s_j\\ is the pooled (unweighted)
standard deviation.

The effective sample size for a group is \$\$\text{ESS} = \frac{(\sum_i
w_i)^2}{\sum_i w_i^2},\$\$ reported as a fraction of the group size.

The weighted energy distance is \$\$E = 2 \sum\_{i,j} p_i q_j \\X_i -
X_j\\ - \sum\_{i,j} p_i p_j \\X_i - X_j\\ - \sum\_{i,j} q_i q_j \\X_i -
X_j\\,\$\$ where \\p\\ and \\q\\ are the normalized treated and control
weights.

## Examples

``` r
# \donttest{
n <- 500; p <- 10
X <- matrix(rnorm(n * p), n, p)
A <- rbinom(n, 1, plogis(0.5 * X[, 1]))
Y <- X[, 1] + rnorm(n)

fit <- forest_balance(X, A, Y)
bal <- compute_balance(X, A, fit$weights)
bal
#> Covariate Balance Diagnostics
#>   n = 500  (n_treated = 259, n_control = 241)
#> ------------------------------------------------------------ 
#>   ESS (treated):  75.1%
#>   ESS (control):  73.5%
#>   Energy distance: 0.0210
#> 
#>   |SMD| for covariates (10 features):
#>     median = 0.0308   Q75 = 0.0405   Q90 = 0.0515   max = 0.0571
#>     All |SMD| below 0.10
#> ------------------------------------------------------------ 

# With nonlinear features
X.nl <- cbind(X[,1]^2, X[,1]*X[,2])
colnames(X.nl) <- c("X1^2", "X1*X2")
bal2 <- compute_balance(X, A, fit$weights, X.trans = X.nl)
bal2
#> Covariate Balance Diagnostics
#>   n = 500  (n_treated = 259, n_control = 241)
#> ------------------------------------------------------------ 
#>   ESS (treated):  75.1%
#>   ESS (control):  73.5%
#>   Energy distance: 0.0210
#> 
#>   |SMD| for covariates (10 features):
#>     median = 0.0308   Q75 = 0.0405   Q90 = 0.0515   max = 0.0571
#>     All |SMD| below 0.10
#> 
#>   |SMD| for transformed covariates (2 features):
#>     median = 0.0436   Q75 = 0.0533   Q90 = 0.0591   max = 0.0629
#>     All |SMD| below 0.10
#> ------------------------------------------------------------ 
# }
```
