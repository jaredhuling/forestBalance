# Augmented (Doubly-Robust) Estimation

## Overview

The default
[`forest_balance()`](https://jaredhuling.github.io/forestBalance/reference/forest_balance.md)
estimator uses kernel energy balancing weights alone. When the kernel
captures all confounding structure, this produces an unbiased ATE
estimate. In practice, however, residual confounding may remain —
particularly when the confounding operates through complex nonlinear
functions of the covariates.

The **augmented** (doubly-robust) estimator combines the balancing
weights with **group-specific** outcome regression models
$`\hat\mu_1(X) = E[Y \mid X, A{=}1]`$ and
$`\hat\mu_0(X) = E[Y \mid X, A{=}0]`$. The ATE is estimated as:

``` math
\hat\tau_{\mathrm{aug}} =
  \frac{1}{n}\sum_{i=1}^n \bigl[\hat\mu_1(X_i) - \hat\mu_0(X_i)\bigr]
  +
  \frac{\sum_i w_i A_i \bigl(Y_i - \hat\mu_1(X_i)\bigr)}
       {\sum_i w_i A_i}
  -
  \frac{\sum_i w_i (1 - A_i) \bigl(Y_i - \hat\mu_0(X_i)\bigr)}
       {\sum_i w_i (1 - A_i)}.
```

The first term is the outcome-model-based regression estimate of the
ATE. The second and third terms are weighted bias corrections that use
the balancing weights to account for errors in the outcome models.

This has a key robustness property: the estimator is consistent if
**either** the kernel (balancing weights) or the outcome models are
correctly specified. When both are reasonable, the augmented estimator
typically has lower bias and RMSE than either approach alone.

## Basic usage

Enable augmentation with `augmented = TRUE`:

``` r
set.seed(123)
dat <- simulate_data(n = 2000, p = 20, ate = 0, dgp = 2)

# Standard (weighting only)
fit <- forest_balance(dat$X, dat$A, dat$Y, num.trees = 500)

# Augmented (doubly-robust)
fit_aug <- forest_balance(dat$X, dat$A, dat$Y, num.trees = 500,
                           augmented = TRUE)

c("Standard" = fit$ate, "Augmented" = fit_aug$ate, "Truth" = 0)
#>   Standard  Augmented      Truth 
#> 0.01664814 0.06350369 0.00000000
```

The `print` method indicates when augmentation is active:

``` r
fit_aug
#> Forest Kernel Energy Balancing (cross-fitted)
#> -------------------------------------------------- 
#>   n = 2,000  (n_treated = 1042, n_control = 958)
#>   Trees: 500
#>   Cross-fitting: 2 folds
#>   Augmented: yes (doubly-robust)
#>   Solver: direct
#>   ATE estimate: 0.0635
#>   Fold ATEs: 0.0898, 0.0372
#>   ESS: treated = 795/1042 (76%)   control = 712/958 (74%)
#> -------------------------------------------------- 
#> Use summary() for covariate balance details.
```

## How it works with cross-fitting

When both `cross.fitting = TRUE` (the default) and `augmented = TRUE`,
the outcome model is cross-fitted **in lockstep** with the kernel:

For each fold $`k`$:

1.  Train the joint forest (for the kernel) on folds $`-k`$.
2.  Train two outcome regression forests on folds $`-k`$: one on treated
    observations ($`\hat\mu_1`$) and one on control observations
    ($`\hat\mu_0`$).
3.  Predict leaf nodes, $`\hat\mu_1`$, and $`\hat\mu_0`$ for fold $`k`$
    out-of-sample.
4.  Compute the fold-level doubly-robust ATE using the cross-fitted
    weights and outcome predictions.

This ensures that neither the kernel nor the outcome models overfit to
the estimation sample. The user does not need to manage the
cross-fitting of $`\hat\mu_1`$ and $`\hat\mu_0`$ separately — it is
handled automatically.

## Simulation comparison

We compare four configurations using DGP 2, which has a rich outcome
surface with moderate confounding ($`n = 5{,}000`$, $`p = 20`$):

``` r
set.seed(123)
nreps <- 5
n <- 5000; p <- 20

methods <- c("Standard", "CF only", "Aug only", "CF + Aug")
res <- matrix(NA, nreps, length(methods), dimnames = list(NULL, methods))

for (r in seq_len(nreps)) {
  dat <- simulate_data(n = n, p = p, ate = 0, dgp = 2)

  res[r, "Standard"] <- forest_balance(dat$X, dat$A, dat$Y, num.trees = 500,
                                        cross.fitting = FALSE,
                                        min.node.size = 10)$ate
  res[r, "CF only"]  <- forest_balance(dat$X, dat$A, dat$Y,
                                        num.trees = 500)$ate
  res[r, "Aug only"] <- forest_balance(dat$X, dat$A, dat$Y, num.trees = 500,
                                        cross.fitting = FALSE,
                                        augmented = TRUE,
                                        min.node.size = 10)$ate
  res[r, "CF + Aug"] <- forest_balance(dat$X, dat$A, dat$Y, num.trees = 500,
                                        augmented = TRUE)$ate
}
```

| Method   |   Bias |     SD |   RMSE |
|:---------|-------:|-------:|-------:|
| Standard | 0.1204 | 0.0437 | 0.1281 |
| CF only  | 0.0151 | 0.0673 | 0.0690 |
| Aug only | 0.0515 | 0.0535 | 0.0743 |
| CF + Aug | 0.0055 | 0.0609 | 0.0611 |

DGP 2: n = 5,000, p = 20, true ATE = 0, 5 reps.

The combination of cross-fitting and augmentation (CF + Aug) achieves
the lowest bias. Cross-fitting alone handles the kernel overfitting,
while augmentation accounts for outcome variation that the kernel does
not fully capture. DGP 2 has a rich outcome surface depending on 8
covariates, of which only 2 are confounders, making the outcome model
particularly valuable.

## User-supplied outcome predictions

You can supply your own group-specific outcome predictions via `mu.hat`
as a list with components `mu1` (predictions of $`E[Y|X, A{=}1]`$) and
`mu0` (predictions of $`E[Y|X, A{=}0]`$):

``` r
set.seed(123)
dat <- simulate_data(n = 1000, p = 20, ate = 0, dgp = 2)

# Fit your own group-specific outcome models
df_all <- data.frame(Y = dat$Y, dat$X)
df1 <- df_all[dat$A == 1, ]
df0 <- df_all[dat$A == 0, ]
lm1 <- lm(Y ~ ., data = df1)
lm0 <- lm(Y ~ ., data = df0)
mu_custom <- list(
  mu1 = predict(lm1, newdata = df_all),
  mu0 = predict(lm0, newdata = df_all)
)

fit_custom <- forest_balance(dat$X, dat$A, dat$Y, num.trees = 500,
                              mu.hat = mu_custom)
fit_custom
#> Forest Kernel Energy Balancing (cross-fitted)
#> -------------------------------------------------- 
#>   n = 1,000  (n_treated = 505, n_control = 495)
#>   Trees: 500
#>   Cross-fitting: 2 folds
#>   Augmented: yes (doubly-robust)
#>   Solver: direct
#>   ATE estimate: 0.1576
#>   Fold ATEs: 0.2724, 0.0428
#>   ESS: treated = 391/505 (77%)   control = 380/495 (77%)
#> -------------------------------------------------- 
#> Use summary() for covariate balance details.
```

When `mu.hat` is supplied with `cross.fitting = TRUE`, the package uses
the provided predictions as-is and emits a message reminding you to
ensure they were cross-fitted externally.

## When to use augmentation

Augmentation is most beneficial when:

- **The confounding is nonlinear** and the kernel alone cannot fully
  capture it.
- **The outcome model is reasonably accurate**, even if not perfectly
  specified.
- **You want robustness** against partial misspecification of either the
  kernel or the outcome model.

The computational cost of augmentation is modest: fitting two additional
`regression_forest` models (one per treatment group) per fold, which is
typically fast relative to the joint forest and kernel construction.

## References

Robins, J.M., Rotnitzky, A. and Zhao, L.P. (1994). Estimation of
regression coefficients when some regressors are not always observed.
*Journal of the American Statistical Association*, 89(427), 846–866.

De, S. and Huling, J.D. (2025). Data adaptive covariate balancing for
causal effect estimation for high dimensional data. *arXiv preprint
arXiv:2512.18069*.
