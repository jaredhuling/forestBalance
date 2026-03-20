test_that("get_leaf_node_matrix returns correct dimensions", {
  dat <- simulate_data(n = 100, p = 5)
  forest <- grf::multi_regression_forest(dat$X, cbind(dat$A, dat$Y),
                                          num.trees = 20, min.node.size = 10)
  lm <- get_leaf_node_matrix(forest, dat$X)
  expect_equal(nrow(lm), 100)
  expect_equal(ncol(lm), 20)
  expect_type(lm, "integer")
})

test_that("get_leaf_node_matrix works with NULL newdata (uses training data)", {
  dat <- simulate_data(n = 80, p = 5)
  forest <- grf::multi_regression_forest(dat$X, cbind(dat$A, dat$Y),
                                          num.trees = 10)
  lm <- get_leaf_node_matrix(forest)
  expect_equal(nrow(lm), 80)
})

test_that("get_leaf_node_matrix works out-of-sample", {
  dat <- simulate_data(n = 100, p = 5)
  forest <- grf::multi_regression_forest(dat$X, cbind(dat$A, dat$Y),
                                          num.trees = 10)
  X_new <- matrix(rnorm(30 * 5), 30, 5)
  lm <- get_leaf_node_matrix(forest, X_new)
  expect_equal(nrow(lm), 30)
  expect_equal(ncol(lm), 10)
})

test_that("leaf_node_kernel returns sparse symmetric matrix", {
  dat <- simulate_data(n = 50, p = 5)
  forest <- grf::multi_regression_forest(dat$X, cbind(dat$A, dat$Y),
                                          num.trees = 20)
  lm <- get_leaf_node_matrix(forest, dat$X)
  K <- leaf_node_kernel(lm)
  expect_s4_class(K, "dsCMatrix")
  expect_equal(nrow(K), 50)
  expect_equal(ncol(K), 50)
  # Diagonal should be 1 (each obs shares leaf with itself in every tree)
  expect_equal(as.numeric(Matrix::diag(K)), rep(1, 50))
  # All entries in [0, 1]
  expect_true(all(K@x >= 0 & K@x <= 1))
})

test_that("leaf_node_kernel_Z has correct dimensions", {
  dat <- simulate_data(n = 50, p = 5)
  forest <- grf::multi_regression_forest(dat$X, cbind(dat$A, dat$Y),
                                          num.trees = 20)
  lm <- get_leaf_node_matrix(forest, dat$X)
  Z <- leaf_node_kernel_Z(lm)
  expect_equal(nrow(Z), 50)
  # Each row has exactly B=20 nonzero entries
  expect_equal(as.numeric(Matrix::rowSums(Z)), rep(20, 50))
})

test_that("leaf_node_kernel matches Z Z^T / B", {
  dat <- simulate_data(n = 50, p = 5)
  forest <- grf::multi_regression_forest(dat$X, cbind(dat$A, dat$Y),
                                          num.trees = 20)
  lm <- get_leaf_node_matrix(forest, dat$X)
  K <- leaf_node_kernel(lm, sparse = FALSE)
  Z <- leaf_node_kernel_Z(lm)
  K_from_Z <- as.matrix(Matrix::tcrossprod(Z) / 20)
  expect_equal(K, K_from_Z, tolerance = 1e-12)
})
