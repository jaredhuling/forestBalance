# Cross-Fitting for Debiased Kernel Estimation

## The overfitting problem

The forest balance estimator uses the data twice: once to fit the random
forest that defines the kernel, and again to estimate the treatment
effect using that kernel. This creates a subtle **overfitting bias**
that persists even at large sample sizes.

To see this, we compare the standard (no cross-fitting) estimator with
small leaf size against the cross-fitted estimator with adaptive leaf
size and an oracle that uses the true propensity scores. We use
$`n = 10{,}000`$ and $`p = 50`$ covariates, where the bias is most
pronounced:

``` r
library(forestBalance)

set.seed(123)
nreps <- 5

res <- matrix(NA, nreps, 3,
              dimnames = list(NULL, c("No CF (mns=10)", "CF (default)", "Oracle IPW")))

for (r in seq_len(nreps)) {
  dat <- simulate_data(n = 10000, p = 50, ate = 0)

  # Standard: no cross-fitting, small min.node.size
  fit_nocf <- forest_balance(dat$X, dat$A, dat$Y,
                              cross.fitting = FALSE, min.node.size = 10,
                              num.trees = 500)
  res[r, "No CF (mns=10)"] <- fit_nocf$ate

  # Cross-fitted with adaptive leaf size (package default)
  fit_cf <- forest_balance(dat$X, dat$A, dat$Y, num.trees = 500)
  res[r, "CF (default)"] <- fit_cf$ate

  # Oracle IPW (true propensity scores)
  ps <- dat$propensity
  w_ipw <- ifelse(dat$A == 1, 1 / ps, 1 / (1 - ps))
  res[r, "Oracle IPW"] <- weighted.mean(dat$Y[dat$A == 1], w_ipw[dat$A == 1]) -
                           weighted.mean(dat$Y[dat$A == 0], w_ipw[dat$A == 0])
}
```

| Method         |    Bias |     SD |   RMSE |
|:---------------|--------:|-------:|-------:|
| No CF (mns=10) |  0.1737 | 0.0361 | 0.1774 |
| CF (default)   |  0.0215 | 0.0383 | 0.0439 |
| Oracle IPW     | -0.0058 | 0.1385 | 0.1387 |

n = 10,000, p = 50, true ATE = 0, 5 replications.

The no-cross-fitting estimator with small leaves has substantial bias.
The default cross-fitted estimator with adaptive leaf size achieves much
lower bias and RMSE.

## Cross-fitting details

The idea follows the **double/debiased machine learning** (DML)
framework of Chernozhukov et al. (2018), adapted to kernel energy
balancing.

Let $`K`$ denote the $`n \times n`$ proximity kernel built from a random
forest trained on the full sample $`(X, A, Y)`$. The kernel captures
which observations are “similar” in terms of confounding structure.
However, because $`K`$ was estimated from the same data used to compute
the ATE, the kernel overfits: it encodes information about the specific
realisation of outcomes, not just the confounding structure. This
creates a bias that does not vanish as $`n \to \infty`$.

### K-fold cross-fitting

Given $`K`$ folds, the cross-fitted estimator proceeds as follows:

1.  **Randomly partition** the data into $`K`$ roughly equal folds
    $`F_1, \ldots, F_K`$.

2.  **For each fold** $`k = 1, \ldots, K`$:

    1.  Train a `multi_regression_forest` on the data in folds
        $`F_{-k}`$ (all folds except $`k`$).
    2.  Using this held-out forest, predict leaf node memberships for
        the observations in fold $`F_k`$.
    3.  Build the proximity kernel $`K_k`$ ($`n_k \times n_k`$) from
        these out-of-sample leaf predictions.
    4.  Compute kernel energy balancing weights $`w_k`$ for the
        observations in fold $`F_k`$ using $`K_k`$.
    5.  Estimate the fold-level ATE via the Hajek estimator:
        ``` math
        \hat\tau_k = \frac{\sum_{i \in F_k} w_i A_i Y_i}
          {\sum_{i \in F_k} w_i A_i}
          - \frac{\sum_{i \in F_k} w_i (1-A_i) Y_i}
          {\sum_{i \in F_k} w_i (1-A_i)}.
        ```

3.  **Average** the fold-level estimates (DML1):
    ``` math
    \hat\tau_{\mathrm{CF}} = \frac{1}{K} \sum_{k=1}^{K} \hat\tau_k.
    ```

## The role of leaf size

Cross-fitting alone is not sufficient to eliminate bias. The **minimum
leaf size** (`min.node.size`) plays a crucial role:

- **Small leaves** (e.g., `min.node.size = 5`): the kernel is very
  granular, distinguishing observations at a fine scale. But with
  out-of-sample predictions, small leaves lead to noisy similarity
  estimates—two observations may share a tiny leaf by chance rather than
  true similarity.

- **Large leaves** (e.g., `min.node.size = 50--100`): the kernel
  captures broader confounding structure. Out-of-sample predictions are
  more stable because large leaves better represent population-level
  similarity.

The optimal leaf size scales with both $`n`$ (more data supports finer
leaves) and $`p`$ (more covariates require coarser leaves to avoid the
curse of dimensionality). `forestBalance` uses an adaptive heuristic:

``` math
\mathrm{min.node.size} = \max\!\Big(20,\;\min\!\big(\lfloor n/200 \rfloor + p,\;\lfloor n/50 \rfloor\big)\Big).
```

``` r
set.seed(123)
nreps <- 5
node_sizes <- c(5, 10, 20, 50, 75, 100)
n <- 10000; p <- 50

leaf_res <- do.call(rbind, lapply(node_sizes, function(mns) {
  ates <- replicate(nreps, {
    dat <- simulate_data(n = n, p = p, ate = 0)
    forest_balance(dat$X, dat$A, dat$Y, num.trees = 500,
                   min.node.size = mns)$ate
  })
  data.frame(mns = mns, bias = mean(ates), sd = sd(ates),
             rmse = sqrt(mean(ates)^2 + var(ates)))
}))

heuristic <- max(20, min(floor(n / 200) + p, floor(n / 50)))
```

| min.node.size |   Bias |     SD |   RMSE |     |
|--------------:|-------:|-------:|-------:|:----|
|             5 | 0.1549 | 0.0489 | 0.1624 |     |
|            10 | 0.0940 | 0.0180 | 0.0957 |     |
|            20 | 0.0811 | 0.0144 | 0.0824 |     |
|            50 | 0.0217 | 0.0156 | 0.0267 |     |
|            75 | 0.0109 | 0.0270 | 0.0292 |     |
|           100 | 0.0101 | 0.0294 | 0.0311 | \<– |

Cross-fitted estimator (n = 10,000, p = 50, 5 reps). Arrow marks the
adaptive default (mns = 100).

Bias decreases with larger leaves until variance begins to dominate. The
adaptive default balances bias reduction against variance.

## Practical usage

The default call uses cross-fitting with the adaptive leaf size:

``` r
dat <- simulate_data(n = 2000, p = 10, ate = 0)
fit <- forest_balance(dat$X, dat$A, dat$Y)
fit
#> Forest Kernel Energy Balancing (cross-fitted)
#> -------------------------------------------------- 
#>   n = 2,000  (n_treated = 690, n_control = 1310)
#>   Trees: 1000
#>   Cross-fitting: 2 folds
#>   Solver: direct
#>   ATE estimate: 0.0646
#>   Fold ATEs: 0.0218, 0.1075
#>   ESS: treated = 469/690 (68%)   control = 995/1310 (76%)
#> -------------------------------------------------- 
#> Use summary() for covariate balance details.
```

To disable cross-fitting (e.g., for speed or to inspect the kernel):

``` r
fit_nocf <- forest_balance(dat$X, dat$A, dat$Y, cross.fitting = FALSE)
fit_nocf
#> Forest Kernel Energy Balancing
#> -------------------------------------------------- 
#>   n = 2,000  (n_treated = 690, n_control = 1310)
#>   Trees: 1000
#>   Solver: direct
#>   ATE estimate: 0.0664
#>   ESS: treated = 399/690 (58%)   control = 909/1310 (69%)
#> -------------------------------------------------- 
#> Use summary() for covariate balance details.
```

To manually set the leaf size:

``` r
fit_custom <- forest_balance(dat$X, dat$A, dat$Y, min.node.size = 50)
fit_custom
#> Forest Kernel Energy Balancing (cross-fitted)
#> -------------------------------------------------- 
#>   n = 2,000  (n_treated = 690, n_control = 1310)
#>   Trees: 1000
#>   Cross-fitting: 2 folds
#>   Solver: direct
#>   ATE estimate: 0.0274
#>   Fold ATEs: -0.02, 0.0749
#>   ESS: treated = 400/690 (58%)   control = 851/1310 (65%)
#> -------------------------------------------------- 
#> Use summary() for covariate balance details.
```

## Choosing the number of folds

The default is `num.folds = 2` (sample splitting). With two folds, each
fold retains half the observations, producing high-quality per-fold
kernels. More folds train the forest on more data but produce smaller
per-fold kernels. Our experiments show that the number of folds has a
modest effect compared to the leaf size. Values of 2–5 all work well:

``` r
set.seed(123)
nreps <- 5
n <- 5000

fold_res <- do.call(rbind, lapply(c(2, 3, 5, 10), function(nfolds) {
  ates <- replicate(nreps, {
    dat <- simulate_data(n = n, p = 10, ate = 0)
    forest_balance(dat$X, dat$A, dat$Y, num.trees = 500,
                   num.folds = nfolds)$ate
  })
  data.frame(folds = nfolds, bias = round(mean(ates), 4),
             sd = round(sd(ates), 4),
             rmse = round(sqrt(mean(ates)^2 + var(ates)), 4))
}))
```

| Folds |   Bias |     SD |   RMSE |
|------:|-------:|-------:|-------:|
|     2 | 0.0193 | 0.0405 | 0.0448 |
|     3 | 0.0398 | 0.0357 | 0.0534 |
|     5 | 0.0609 | 0.0446 | 0.0755 |
|    10 | 0.2014 | 0.0456 | 0.2065 |

Effect of number of folds (n = 5,000, adaptive leaf size).

## References

Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C.,
Newey, W. and Robins, J. (2018). Double/debiased machine learning for
treatment and structural parameters. *The Econometrics Journal*, 21(1),
C1–C68.

De, S. and Huling, J.D. (2025). Data adaptive covariate balancing for
causal effect estimation for high dimensional data. *arXiv preprint
arXiv:2512.18069*.
