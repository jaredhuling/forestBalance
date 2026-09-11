test_that("kernel_balance direct solver satisfies weight constraints", {
  set.seed(123)
  dat <- simulate_data(n = 500, p = 5)
  forest <- grf::multi_regression_forest(dat$X, scale(cbind(dat$A, dat$Y)),
                                          num.trees = 100, min.node.size = 10)
  lm <- get_leaf_node_matrix(forest, dat$X)
  K <- leaf_node_kernel(lm)
  bal <- kernel_balance(dat$A, kern = K, solver = "direct")

  n1 <- sum(dat$A == 1); n0 <- sum(dat$A == 0)
  expect_equal(sum(bal$weights[dat$A == 1]), n1, tolerance = 1e-6)
  expect_equal(sum(bal$weights[dat$A == 0]), n0, tolerance = 1e-6)
  expect_equal(bal$solver, "direct")
})

test_that("kernel_balance CG solver satisfies weight constraints", {
  set.seed(123)
  dat <- simulate_data(n = 500, p = 5)
  forest <- grf::multi_regression_forest(dat$X, scale(cbind(dat$A, dat$Y)),
                                          num.trees = 100, min.node.size = 10)
  lm <- get_leaf_node_matrix(forest, dat$X)
  Z <- leaf_node_kernel_Z(lm)
  bal <- kernel_balance(dat$A, Z = Z, num.trees = 100, solver = "cg")

  n1 <- sum(dat$A == 1); n0 <- sum(dat$A == 0)
  expect_equal(sum(bal$weights[dat$A == 1]), n1, tolerance = 0.1)
  expect_equal(sum(bal$weights[dat$A == 0]), n0, tolerance = 0.1)
  expect_equal(bal$solver, "cg")
})

test_that("kernel_balance errors on invalid input", {
  expect_error(kernel_balance(c(1, 0, 1)), "Either")
  expect_error(kernel_balance(c(1, 1, 1), kern = diag(3)), "both treated")
  expect_error(kernel_balance(c(1, 0), kern = diag(3)), "dimensions")
})

test_that("kernel_balance direct and CG give similar ATEs", {
  set.seed(123)
  dat <- simulate_data(n = 300, p = 5)
  forest <- grf::multi_regression_forest(dat$X, scale(cbind(dat$A, dat$Y)),
                                          num.trees = 100, min.node.size = 20)
  lm <- get_leaf_node_matrix(forest, dat$X)
  K <- leaf_node_kernel(lm)
  Z <- leaf_node_kernel_Z(lm)

  w_dir <- kernel_balance(dat$A, kern = K, solver = "direct")$weights
  w_cg  <- kernel_balance(dat$A, Z = Z, num.trees = 100, solver = "cg")$weights

  ate_dir <- weighted.mean(dat$Y[dat$A == 1], w_dir[dat$A == 1]) -
             weighted.mean(dat$Y[dat$A == 0], w_dir[dat$A == 0])
  ate_cg  <- weighted.mean(dat$Y[dat$A == 1], w_cg[dat$A == 1]) -
             weighted.mean(dat$Y[dat$A == 0], w_cg[dat$A == 0])

  expect_equal(ate_dir, ate_cg, tolerance = 0.05)
})

test_that("lambda = 0 reproduces the unpenalized solution and lambda > 0 keeps constraints", {
  set.seed(123)
  dat <- simulate_data(n = 400, p = 5)
  forest <- grf::multi_regression_forest(dat$X, scale(cbind(dat$A, dat$Y)),
                                          num.trees = 100, min.node.size = 10)
  lm <- get_leaf_node_matrix(forest, dat$X)
  K <- leaf_node_kernel(lm)
  Z <- leaf_node_kernel_Z(lm)
  n1 <- sum(dat$A == 1); n0 <- sum(dat$A == 0)

  w0 <- kernel_balance(dat$A, kern = K, solver = "direct")$weights
  w0b <- kernel_balance(dat$A, kern = K, solver = "direct", lambda = 0)$weights
  expect_equal(w0, w0b)

  for (lam in c(0.01, 1)) {
    bd <- kernel_balance(dat$A, kern = K, solver = "direct", lambda = lam)
    bc <- kernel_balance(dat$A, Z = Z, num.trees = 100, solver = "cg", lambda = lam, tol = 1e-10)
    bb <- kernel_balance(dat$A, Z = Z, leaf_matrix = lm, num.trees = 100, solver = "bj",
                         lambda = lam, tol = 1e-10)
    expect_equal(sum(bd$weights[dat$A == 1]), n1, tolerance = 1e-6)
    expect_equal(sum(bd$weights[dat$A == 0]), n0, tolerance = 1e-6)
    expect_equal(bd$weights, bc$weights, tolerance = 1e-4)
    expect_equal(bd$weights, bb$weights, tolerance = 1e-4)
    expect_true(is.finite(bc$iters) && bc$iters >= 1)
    # KKT check for the direct solution on the treated block:
    # (K_tt + lam I) w_t - n1^2 b_t must be constant across treated units.
    idx_t <- which(dat$A == 1)
    K_tt <- as.matrix(K[idx_t, idx_t])
    rs <- as.numeric(Matrix::rowSums(K))
    b_t <- rs[idx_t] / (n1 * length(dat$A))
    resid <- (K_tt + lam * diag(n1)) %*% bd$weights[idx_t] - n1^2 * b_t
    expect_lt(diff(range(resid)), 1e-6 * max(1, abs(mean(resid))))
  }
  expect_error(kernel_balance(dat$A, kern = K, lambda = -1), "nonnegative")
})
