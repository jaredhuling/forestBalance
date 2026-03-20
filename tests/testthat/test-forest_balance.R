test_that("forest_balance returns correct class and structure", {
  dat <- simulate_data(n = 200, p = 5)
  fit <- forest_balance(dat$X, dat$A, dat$Y, num.trees = 50)
  expect_s3_class(fit, "forest_balance")
  expect_type(fit$ate, "double")
  expect_length(fit$weights, 200)
  expect_equal(fit$n, 200)
  expect_equal(fit$n1 + fit$n0, 200L)
})

test_that("forest_balance cross.fitting=TRUE is default", {
  dat <- simulate_data(n = 200, p = 5)
  fit <- forest_balance(dat$X, dat$A, dat$Y, num.trees = 50)
  expect_true(fit$crossfit)
  expect_length(fit$fold_ates, 2)  # default num.folds=2
  expect_null(fit$kernel)
})

test_that("forest_balance cross.fitting=FALSE returns kernel", {
  dat <- simulate_data(n = 200, p = 5)
  fit <- forest_balance(dat$X, dat$A, dat$Y, num.trees = 50,
                         cross.fitting = FALSE)
  expect_false(isTRUE(fit$crossfit))
  expect_false(is.null(fit$kernel) && fit$solver == "direct")
})

test_that("forest_balance adaptive min.node.size works", {
  # Heuristic: max(20, min(floor(n/200) + p, floor(n/50)))
  # For n=200, p=5: max(20, min(1+5, 4)) = max(20, 4) = 20
  # For n=10000, p=50: max(20, min(50+50, 200)) = max(20, 100) = 100
  expect_equal(max(20, min(floor(200/200) + 5, floor(200/50))), 20)
  expect_equal(max(20, min(floor(10000/200) + 50, floor(10000/50))), 100)
})

test_that("forest_balance errors on invalid input", {
  dat <- simulate_data(n = 100, p = 5)
  expect_error(forest_balance(dat$X, dat$A[-1], dat$Y), "same number")
  expect_error(forest_balance(dat$X, rep(2, 100), dat$Y), "binary")
})

test_that("forest_balance print and summary work", {
  dat <- simulate_data(n = 200, p = 5)
  fit <- forest_balance(dat$X, dat$A, dat$Y, num.trees = 50)
  expect_output(print(fit), "Forest Kernel Energy Balancing")
  expect_output(print(fit), "cross-fitted")
  s <- summary(fit)
  expect_s3_class(s, "summary.forest_balance")
  expect_output(print(s), "Covariate Balance")
})
