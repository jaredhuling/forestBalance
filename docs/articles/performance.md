# Performance and Scalability

## Overview

This vignette benchmarks the `forestBalance` pipeline to characterize
how speed and memory scale with sample size ($`n`$) and number of trees
($`B`$).

`forestBalance` adaptively selects between two linear solvers based on
the per-fold sample size and the kernel density (which depends on
`min.node.size`):

- **Direct**: sparse Cholesky on the treated and control sub-blocks of
  the kernel matrix. Exact solution. Preferred for smaller problems.
- **Conjugate gradient (CG)**: iterative solver using the factored $`Z`$
  representation ($`K = Z Z^\top / B`$), so the full kernel matrix is
  never formed. Preferred when $`n`$ is large or the adaptive
  `min.node.size` creates dense kernels (e.g., high $`p`$).

``` r
library(forestBalance)
library(grf)
library(Matrix)
```

## Mathematical background

### The kernel energy balancing system

Given a kernel matrix $`K \in \mathbb{R}^{n \times n}`$ and a binary
treatment vector $`A \in \{0,1\}^n`$ with $`n_1`$ treated and $`n_0`$
control units, the kernel energy balancing weights $`w`$ are obtained by
solving the linear system
``` math
K_q \, w = z,
```
where $`K_q`$ is the *modified kernel* and $`z`$ is a right-hand side
vector determined by a linear constraint.

The modified kernel is defined element-wise as
``` math
K_q(i,j)
  = \frac{A_i \, A_j \, K(i,j)}{n_1^2}
  + \frac{(1-A_i)(1-A_j) \, K(i,j)}{n_0^2}.
```

A critical structural observation is that $`K_q(i,j) = 0`$ whenever
$`A_i \neq A_j`$: the treated–control cross-blocks are identically zero.
Therefore $`K_q`$ is **block-diagonal**:
``` math
K_q = \begin{pmatrix} K_{tt} / n_1^2 & 0 \\ 0 & K_{cc} / n_0^2 \end{pmatrix},
```
where $`K_{tt} = K[A{=}1,\; A{=}1]`$ is the $`n_1 \times n_1`$ treated
sub-block and $`K_{cc} = K[A{=}0,\; A{=}0]`$ is the $`n_0 \times n_0`$
control sub-block.

The right-hand side vector $`b`$ has a similarly separable structure.
Writing $`\mathbf{r} = K \mathbf{1}`$ (the row sums of $`K`$), we have
``` math
b_i = \frac{A_i \, r_i}{n_1 \, n} + \frac{(1-A_i) \, r_i}{n_0 \, n}.
```

The full system decomposes into two **independent** sub-problems:

- **Treated block:** solve $`K_{tt} \, w_t = n_1^2 \, z_t`$ of size
  $`n_1`$,
- **Control block:** solve $`K_{cc} \, w_c = n_0^2 \, z_c`$ of size
  $`n_0`$,

where $`z_t`$ and $`z_c`$ are the constraint-adjusted right-hand sides
for each group (each requiring two preliminary solves to determine the
Lagrange multiplier).

### The kernel factorization

The proximity kernel is $`K = Z Z^\top / B`$, where $`Z`$ is a sparse
$`n \times L`$ indicator matrix ($`L = \sum_{b=1}^B L_b`$, total leaves
across all trees). Each row of $`Z`$ has exactly $`B`$ nonzero entries
(one per tree), so $`Z`$ has $`nB`$ nonzeros total.

### Direct solver (block Cholesky)

For moderate $`n`$, we form the sub-block kernels explicitly:
``` math
K_{tt} = Z_t Z_t^\top / B, \qquad K_{cc} = Z_c Z_c^\top / B,
```
where $`Z_t = Z[A{=}1, \,\cdot\,]`$ and $`Z_c = Z[A{=}0, \,\cdot\,]`$.
Each sub-block is computed via a sparse `tcrossprod`. The linear systems
are then solved by sparse Cholesky factorization.

**Computational Cost:** $`O(nB \cdot \bar{s})`$ for the two sub-block
cross-products (where $`\bar{s}`$ is the average leaf size), plus
$`O(n_1^{3/2} + n_0^{3/2})`$ for sparse Cholesky (in the best case).

### CG solver (matrix-free)

For large $`n`$, forming $`K_{tt}`$ and $`K_{cc}`$ becomes expensive.
The conjugate gradient (CG) solver avoids forming *any* kernel matrix by
operating on the factored representation. To solve
``` math
K_{tt} \, x = r \quad \Longleftrightarrow \quad Z_t Z_t^\top x = B \, r,
```
CG only needs matrix–vector products of the form
``` math
v \;\mapsto\; Z_t \bigl(Z_t^\top v\bigr),
```
each costing $`O(n_1 B)`$ via two sparse matrix–vector multiplies. The
same applies to the control block with $`Z_c`$.

Here, the memory use is $`O(nB)`$ for $`Z`$ alone, versus $`O(n^2)`$ for
the kernel. At $`n = 25{,}000`$ with $`B = 1{,}000`$ this is ~300 MB vs
~4.7 GB. Each CG iteration costs $`O(n_g B)`$ (where $`n_g`$ is the
group size). Convergence typically requires 100–200 iterations,
independent of $`n`$, so the total cost is
$`O(n_g B \cdot T_{\text{iter}})`$. The six required solves (three per
block) are independent and each converges in a similar number of
iterations.

**Computational Cost:**
$`O\bigl((n_1 + n_0) \cdot B \cdot T_{\text{iter}}\bigr)`$ where
$`T_{\text{iter}} \approx 100\text{--}200`$. This scales linearly in
both $`n`$ and $`B`$, making it the preferred solver for large problems.

## End-to-end timing

We benchmark the full
[`forest_balance()`](http://jaredhuling.org/forestBalance/reference/forest_balance.md)
pipeline (forest fitting, leaf extraction, kernel/Z construction, and
weight computation) across a range of sample sizes with $`B = 1{,}000`$
trees. To show the effect of solver choice, we run each $`n`$ with both
`solver = "direct"` and `solver = "cg"` where feasible:

``` r
n_vals <- c(500, 1000, 2500, 5000, 10000, 25000)
B <- 1000
p <- 10

bench <- do.call(rbind, lapply(n_vals, function(nn) {
  set.seed(123)
  dat <- simulate_data(n = nn, p = p)

  # Auto (default)
  t_auto <- system.time(
    fit_auto <- forest_balance(dat$X, dat$A, dat$Y, num.trees = B)
  )["elapsed"]

  # Direct (skip for n > 10000 — too slow)
  if (nn <= 10000) {
    t_dir <- system.time(
      fit_dir <- forest_balance(dat$X, dat$A, dat$Y, num.trees = B,
                                solver = "direct")
    )["elapsed"]
  } else {
    t_dir <- NA
  }

  # CG (run for all n)
  t_cg <- system.time(
    fit_cg <- forest_balance(dat$X, dat$A, dat$Y, num.trees = B,
                             solver = "cg")
  )["elapsed"]

  data.frame(n = nn, trees = B,
             auto = t_auto, direct = t_dir, cg = t_cg,
             auto_solver = fit_auto$solver)
}))
```

|     n | Trees | Direct (s) | CG (s) | Auto picks |
|------:|------:|-----------:|:-------|:-----------|
|   500 |  1000 |       0.09 | 0.29   | cg         |
|  1000 |  1000 |       0.21 | 0.56   | direct     |
|  2500 |  1000 |       0.89 | 1.64   | direct     |
|  5000 |  1000 |       2.87 | 4.51   | direct     |
| 10000 |  1000 |      10.90 | 16.56  | direct     |
| 25000 |  1000 |          – | 77.29  | cg         |

Full pipeline time by solver.

![](performance_files/figure-html/timing-plot-1.png)

The circled points show which solver `"auto"` selects at each $`n`$. The
switchover depends on both the per-fold sample size and the adaptive
`min.node.size` (which creates denser kernels at higher $`p`$).

## Scaling with number of trees

``` r
tree_vals <- c(200, 500, 1000, 2000)
n_test <- c(1000, 5000)

tree_bench <- do.call(rbind, lapply(n_test, function(nn) {
  do.call(rbind, lapply(tree_vals, function(B) {
    set.seed(123)
    dat <- simulate_data(n = nn, p = 10)
    t <- system.time(
      fit <- forest_balance(dat$X, dat$A, dat$Y, num.trees = B)
    )["elapsed"]
    data.frame(n = nn, trees = B, time = t)
  }))
}))
```

|    n | Trees | Time (s) |
|-----:|------:|---------:|
| 1000 |   200 |     0.06 |
| 1000 |   500 |     0.11 |
| 1000 |  1000 |     0.20 |
| 1000 |  2000 |     0.38 |
| 5000 |   200 |     0.93 |
| 5000 |   500 |     1.37 |
| 5000 |  1000 |     2.19 |
| 5000 |  2000 |     4.30 |

Pipeline time across tree counts.

![](performance_files/figure-html/tree-plot-1.png)

The pipeline scales approximately linearly in both $`n`$ and $`B`$.

## Kernel sparsity and memory

When the direct solver is used ($`n \le 5{,}000`$), the kernel is formed
as a sparse matrix. The kernel becomes sparser as $`n`$ grows, because
fewer pairs of observations share a leaf in any given tree.

For $`n > 5{,}000`$, the CG solver avoids forming the kernel entirely;
memory usage is then dominated by the sparse indicator matrix $`Z`$
($`n \times L`$ with $`nB`$ nonzeros). The table and plot below show the
actual memory usage for both regimes:

``` r
mem_data <- do.call(rbind, lapply(c(500, 1000, 2500, 5000, 10000, 25000),
  function(nn) {
    set.seed(123)
    dat <- simulate_data(n = nn, p = 10)
    B_val <- 1000

    # Use adaptive min.node.size (same as forest_balance default)
    mns <- max(20L, min(floor(nn / 200) + ncol(dat$X), floor(nn / 50)))

    fit_forest <- multi_regression_forest(
      dat$X, scale(cbind(dat$A, dat$Y)),
      num.trees = B_val, min.node.size = mns
    )
    leaf_mat <- get_leaf_node_matrix(fit_forest, dat$X)

    # Z matrix (always computed)
    Z <- leaf_node_kernel_Z(leaf_mat)
    z_mb <- as.numeric(object.size(Z)) / 1e6

    # Kernel
    K <- leaf_node_kernel(leaf_mat)
    nnz <- length(K@x)
    pct_nz <- round(100 * nnz / as.numeric(nn)^2, 1)
    k_mb <- round(as.numeric(object.size(K)) / 1e6, 1)

    # Which solver would auto pick?
    # (approximate: uses per-fold size n/2 for CF K=2)
    n_fold <- nn %/% 2
    auto_solver <- if (n_fold > 5000 || mns > n_fold / 20) "cg" else "direct"

    dense_mb <- round(8 * as.numeric(nn)^2 / 1e6, 0)

    data.frame(n = nn, mns = mns, pct_nz = pct_nz,
               sparse_MB = k_mb, Z_MB = round(z_mb, 1),
               dense_MB = dense_mb, solver = auto_solver)
  }
))
```

|     n | mns | Solver | K % nonzero | Stored     | Actual (MB) | Dense K (MB) | Actual / Dense |
|------:|----:|:-------|:------------|:-----------|------------:|-------------:|:---------------|
|   500 |  20 | cg     | 38%         | Z (factor) |         6.0 |            2 | 300%           |
|  1000 |  20 | direct | 29.9%       | K (sparse) |         3.6 |            8 | 45%            |
|  2500 |  22 | direct | 20.2%       | K (sparse) |        15.1 |           50 | 30.2%          |
|  5000 |  35 | direct | 15.4%       | K (sparse) |        46.3 |          200 | 23.1%          |
| 10000 |  60 | direct | 11.9%       | K (sparse) |       143.0 |          800 | 17.9%          |
| 25000 | 135 | cg     | 8.4%        | Z (factor) |       300.4 |         5000 | 6%             |

Memory usage with adaptive leaf size (p = 10).

![](performance_files/figure-html/memory-plot-1.png)

At $`n = 25{,}000`$, a dense kernel would require **4.7 GB**. The CG
solver stores only the $`Z`$ matrix (**~300 MB**), a **16-fold**
reduction.

## Summary

- The full
  [`forest_balance()`](http://jaredhuling.org/forestBalance/reference/forest_balance.md)
  pipeline scales to **$`n = 25{,}000`$** in under a minute with 1,000
  trees.
- The adaptive `min.node.size` heuristic adjusts leaf size to both $`n`$
  and $`p`$, creating denser kernels than the old default of 10.
- The solver switchover between direct and CG depends on both the
  per-fold sample size and the kernel density. The `"auto"` setting
  handles this transparently.
- The pipeline scales approximately **linearly** in both $`n`$ and
  $`B`$.
