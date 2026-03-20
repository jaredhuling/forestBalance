test_that("simulate_data returns correct structure", {
  dat <- simulate_data(n = 100, p = 10, ate = 0)
  expect_type(dat, "list")
  expect_equal(nrow(dat$X), 100)
  expect_equal(ncol(dat$X), 10)
  expect_length(dat$A, 100)
  expect_length(dat$Y, 100)
  expect_length(dat$propensity, 100)
  expect_equal(dat$ate, 0)
  expect_equal(dat$n, 100)
  expect_equal(dat$p, 10)
  expect_true(all(dat$A %in% c(0, 1)))
  expect_true(all(dat$propensity > 0 & dat$propensity < 1))
})

test_that("simulate_data respects ate parameter", {
  dat <- simulate_data(n = 200, p = 5, ate = 5)
  expect_equal(dat$ate, 5)
})

test_that("simulate_data errors for p < 5", {
  expect_error(simulate_data(n = 100, p = 3), "p must be at least 5")
})

test_that("simulate_data column names are set", {
  dat <- simulate_data(n = 50, p = 7)
  expect_equal(colnames(dat$X), paste0("X", 1:7))
})
