test_that("compute_balance with uniform weights gives ESS = 100%", {
  dat <- simulate_data(n = 100, p = 5)
  bal <- compute_balance(dat$X, dat$A, rep(1, 100), energy.dist = FALSE)
  expect_equal(bal$ess_treated, 1, tolerance = 1e-10)
  expect_equal(bal$ess_control, 1, tolerance = 1e-10)
})

test_that("compute_balance returns correct structure", {
  dat <- simulate_data(n = 100, p = 5)
  bal <- compute_balance(dat$X, dat$A, rep(1, 100), energy.dist = FALSE)
  expect_s3_class(bal, "forest_balance_diag")
  expect_length(bal$smd, 5)
  expect_true(all(bal$smd >= 0))
  expect_equal(bal$n, 100)
  expect_equal(bal$n1 + bal$n0, 100L)
  expect_true(is.na(bal$energy_dist))  # skipped
  expect_null(bal$smd_trans)
})

test_that("compute_balance X.trans works", {
  dat <- simulate_data(n = 100, p = 5)
  X.nl <- cbind(dat$X[, 1]^2, dat$X[, 1] * dat$X[, 2])
  colnames(X.nl) <- c("X1sq", "X1X2")
  bal <- compute_balance(dat$X, dat$A, rep(1, 100), X.trans = X.nl,
                          energy.dist = FALSE)
  expect_length(bal$smd_trans, 2)
  expect_equal(names(bal$smd_trans), c("X1sq", "X1X2"))
})

test_that("compute_balance energy distance is computed for small n", {
  dat <- simulate_data(n = 100, p = 5)
  bal <- compute_balance(dat$X, dat$A, rep(1, 100), energy.dist = TRUE)
  expect_true(!is.na(bal$energy_dist))
  expect_true(bal$energy_dist >= 0)
})

test_that("compute_balance print works", {
  dat <- simulate_data(n = 100, p = 5)
  bal <- compute_balance(dat$X, dat$A, rep(1, 100), energy.dist = FALSE)
  expect_output(print(bal), "Covariate Balance Diagnostics")
})
