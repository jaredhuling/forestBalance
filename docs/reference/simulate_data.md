# Simulate observational study data with confounding

Generates data from a design with nonlinear confounding, where
covariates jointly influence both treatment assignment and the outcome
through non-trivial functions. The true average treatment effect is
known, allowing evaluation of estimator performance.

## Usage

``` r
simulate_data(n = 500, p = 10, ate = 0, rho = -0.25, sigma = 1)
```

## Arguments

- n:

  Sample size. Default is 500.

- p:

  Number of covariates. Must be at least 5. Default is 10.

- ate:

  True average treatment effect. Default is 0.

- rho:

  Correlation parameter for the AR(1) covariance structure among
  covariates: \\\Sigma\_{jk} = \rho^{\|j-k\|}\\. Default is \\-0.25\\.

- sigma:

  Noise standard deviation for the outcome. Default is 1.

## Value

A list with the following elements:

- X:

  The \\n \times p\\ covariate matrix.

- A:

  Binary (0/1) treatment assignment vector.

- Y:

  Observed outcome vector.

- propensity:

  True propensity scores \\P(A=1 \mid X)\\.

- ate:

  The true ATE used in the simulation.

- n:

  Sample size.

- p:

  Number of covariates.

## Details

The data generating process is:

**Covariates:** \\X \sim N(0, \Sigma)\\ where \\\Sigma\_{jk} =
\rho^{\|j-k\|}\\.

**Propensity score:** \\P(A = 1 \mid X) = 0.25(1 + B(X_1; 2, 4))\\,
where \\B(\cdot; 2, 4)\\ is the Beta(2,4) density. This creates a
nonlinear relationship between \\X_1\\ and treatment assignment.

**Outcome model:** \$\$Y = 2(X_1 - 1) + 2\\B(X_1; 2, 4) + X_2 +
2\\B(X_5; 2, 4) + \tau \cdot A + \varepsilon,\$\$ where \\\varepsilon
\sim N(0, \sigma^2)\\ and \\\tau\\ is the ATE. Confounding arises
because \\X_1\\ affects both the propensity score and the outcome
nonlinearly.

## Examples

``` r
dat <- simulate_data(n = 500, p = 10, ate = 0)
str(dat)
#> List of 7
#>  $ X         : num [1:500, 1:10] -1.17 -1.51 2.69 0.23 -0.17 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:10] "X1" "X2" "X3" "X4" ...
#>  $ A         : int [1:500] 0 0 0 1 0 0 1 0 0 0 ...
#>  $ Y         : num [1:500] -4.703 -5.379 -0.876 1.799 2.629 ...
#>  $ propensity: num [1:500] 0.25 0.25 0.25 0.775 0.25 ...
#>  $ ate       : num 0
#>  $ n         : num 500
#>  $ p         : num 10

# True ATE is 0, naive estimate is biased
mean(dat$Y[dat$A == 1]) - mean(dat$Y[dat$A == 0])
#> [1] 0.8173965
```
