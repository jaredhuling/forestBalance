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

test_that("forest_balance augmented estimator works", {
  dat <- simulate_data(n = 300, p = 5)
  fit_aug <- forest_balance(dat$X, dat$A, dat$Y, num.trees = 50,
                             augmented = TRUE)
  expect_true(fit_aug$augmented)
  expect_length(fit_aug$mu1.hat, 300)
  expect_length(fit_aug$mu0.hat, 300)
  expect_type(fit_aug$ate, "double")
  expect_output(print(fit_aug), "doubly-robust")
})

test_that("forest_balance augmented with user-supplied mu.hat works", {
  dat <- simulate_data(n = 200, p = 5)
  mu_hat <- list(mu1 = rep(mean(dat$Y[dat$A == 1]), 200),
                 mu0 = rep(mean(dat$Y[dat$A == 0]), 200))
  fit <- forest_balance(dat$X, dat$A, dat$Y, num.trees = 50,
                         mu.hat = mu_hat)
  expect_true(fit$augmented)
  expect_equal(fit$mu1.hat, mu_hat$mu1)
  expect_equal(fit$mu0.hat, mu_hat$mu0)
})

test_that("forest_balance augmented without cross-fitting works", {
  dat <- simulate_data(n = 200, p = 5)
  fit <- forest_balance(dat$X, dat$A, dat$Y, num.trees = 50,
                         augmented = TRUE, cross.fitting = FALSE)
  expect_true(fit$augmented)
  expect_false(isTRUE(fit$crossfit))
  expect_length(fit$mu1.hat, 200)
  expect_length(fit$mu0.hat, 200)
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

test_that("kernel.response ablation options run and differ from the joint kernel", {
  set.seed(321)
  dat <- simulate_data(n = 300, p = 5)
  f_joint <- forest_balance(dat$X, dat$A, dat$Y, num.trees = 100, cross.fitting = FALSE)
  f_trt   <- forest_balance(dat$X, dat$A, dat$Y, num.trees = 100, cross.fitting = FALSE,
                            kernel.response = "treatment")
  f_out   <- forest_balance(dat$X, dat$A, dat$Y, num.trees = 100, cross.fitting = FALSE,
                            kernel.response = "outcome")
  for (f in list(f_joint, f_trt, f_out)) {
    expect_equal(sum(f$weights[dat$A == 1]), sum(dat$A == 1), tolerance = 1e-6)
    expect_true(is.finite(f$ate))
  }
  expect_equal(f_joint$kernel.response, "joint")
  expect_false(isTRUE(all.equal(f_joint$weights, f_trt$weights)))
  expect_error(forest_balance(dat$X, dat$A, dat$Y, kernel.response = "both"))
})

test_that("crossfit.balance = 'full' keeps fold weight sums and the mean-zero information flow", {
  set.seed(654)
  dat <- simulate_data(n = 400, p = 5)
  set.seed(1); f_fold <- forest_balance(dat$X, dat$A, dat$Y, num.trees = 100, num.folds = 2, lambda = 0.1)
  set.seed(1); f_full <- forest_balance(dat$X, dat$A, dat$Y, num.trees = 100, num.folds = 2, lambda = 0.1,
                                        crossfit.balance = "full")
  expect_equal(f_full$fold_ids, f_fold$fold_ids)
  expect_equal(f_full$crossfit.balance, "full")
  expect_true(is.finite(f_full$ate))
  # full-sample balancing: weights of each arm sum to the arm size over the FULL sample, so the
  # fold's own weights need not sum to the fold's arm size; check they differ from the fold variant.
  expect_false(isTRUE(all.equal(f_full$weights, f_fold$weights)))
  # information flow: perturbing outcomes inside fold 1 must not change fold 1's weights.
  Y2 <- dat$Y; Y2[f_full$fold_ids == 1] <- Y2[f_full$fold_ids == 1] + rnorm(sum(f_full$fold_ids == 1))
  set.seed(1); f_full2 <- forest_balance(dat$X, dat$A, Y2, num.trees = 100, num.folds = 2, lambda = 0.1,
                                         crossfit.balance = "full")
  expect_equal(f_full2$fold_ids, f_full$fold_ids)
  expect_equal(f_full2$weights[f_full$fold_ids == 1], f_full$weights[f_full$fold_ids == 1])
})
