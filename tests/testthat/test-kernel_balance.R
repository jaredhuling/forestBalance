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
