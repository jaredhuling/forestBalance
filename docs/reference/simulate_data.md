# Simulate observational study data with confounding

Generates data from a design with nonlinear confounding, where
covariates jointly influence both treatment assignment and the outcome
through non-trivial functions. The true average treatment effect is
known, allowing evaluation of estimator performance.

## Usage

``` r
simulate_data(n = 500, p = 10, ate = 0, rho = -0.25, sigma = 1, dgp = 1)
```

## Arguments

- n:

  Sample size. Default is 500.

- p:

  Number of covariates. Must be at least 5 for `dgp = 1` and at least 8
  for `dgp = 2`. Default is 10.

- ate:

  True average treatment effect. Default is 0.

- rho:

  Correlation parameter for the AR(1) covariance structure among
  covariates: \\\Sigma\_{jk} = \rho^{\|j-k\|}\\. Default is \\-0.25\\.

- sigma:

  Noise standard deviation for the outcome. Default is 1.

- dgp:

  Integer selecting the data generating process. Default is 1. See
  Details.

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

- dgp:

  The DGP that was used.

## Details

Both DGPs generate covariates from \\X \sim N(0, \Sigma)\\ where
\\\Sigma\_{jk} = \rho^{\|j-k\|}\\.

**DGP 1** (default): Confounding through \\X_1\\ via a Beta density.

- Propensity: \\P(A=1\|X) = 0.25(1 + B(X_1; 2, 4))\\ where \\B\\ is the
  Beta(2,4) density.

- Outcome: \\Y = 2(X_1-1) + 2 B(X_1;2,4) + X_2 + 2 B(X_5;2,4) + \tau A +
  \varepsilon\\.

**DGP 2**: Rich outcome surface with moderate confounding. Designed to
illustrate the benefit of the augmented estimator. Confounding operates
through \\X_1\\ and \\X_2\\, while the outcome depends on \\X_1, \ldots,
X_8\\ with interactions and nonlinearities.

- Propensity: \\P(A=1\|X) = \mathrm{logit}^{-1}(0.6 X_1 - 0.4 X_2 + 0.2
  X_1 X_2)\\.

- Outcome: \\Y = 2 X_1 + X_2^2 - 1.5 X_3 + \sin(2 X_4) + X_5 X_1 + 0.8
  X_6 - \cos(X_7) + 0.5 X_8 + \tau A + \varepsilon\\.

## Examples

``` r
dat1 <- simulate_data(n = 500, p = 10, dgp = 1)
dat2 <- simulate_data(n = 500, p = 20, dgp = 2)
```
