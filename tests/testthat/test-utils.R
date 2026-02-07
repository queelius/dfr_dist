# =============================================================================
# Tests for get_params()
# =============================================================================

test_that("get_params returns par when both par and default provided", {
  expect_equal(get_params(c(1, 2), c(3, 4)), c(1, 2))
})

test_that("get_params returns default when par is NULL", {
  expect_equal(get_params(NULL, c(3, 4)), c(3, 4))
})

test_that("get_params returns default when par is all NA", {
  expect_equal(get_params(c(NA, NA), c(3, 4)), c(3, 4))
})

test_that("get_params returns NULL when both par and default are NULL", {
  expect_null(get_params(NULL, NULL))
})

test_that("get_params merges NA values from default", {
  expect_equal(get_params(c(1, NA), c(3, 4)), c(1, 4))
})

test_that("get_params errors on length mismatch", {
  expect_error(get_params(c(1, 2, 3), c(4, 5)), "Parameter length mismatch")
})

test_that("get_params returns par when default is NULL", {
  expect_equal(get_params(c(1, 2), NULL), c(1, 2))
})

# =============================================================================
# Tests for require_params()
# =============================================================================

test_that("require_params returns par when available", {
  expect_equal(require_params(c(1, 2), c(3, 4)), c(1, 2))
})

test_that("require_params returns default when par is NULL", {
  expect_equal(require_params(NULL, c(3, 4)), c(3, 4))
})

test_that("require_params errors when both par and default are NULL", {
  expect_error(require_params(NULL, NULL), "Parameters required")
})

# =============================================================================
# Tests for get_delta()
# =============================================================================

test_that("get_delta extracts delta column when present", {
  df <- data.frame(t = c(1, 2, 3), delta = c(1, 0, 1))
  expect_equal(get_delta(df), c(1, 0, 1))
})

test_that("get_delta returns all-ones when delta column absent", {
  df <- data.frame(t = c(1, 2, 3))
  expect_equal(get_delta(df), c(1, 1, 1))
})

test_that("get_delta respects custom delta_col", {
  df <- data.frame(t = c(1, 2, 3), censor = c(0, 1, 0))
  expect_equal(get_delta(df, delta_col = "censor"), c(0, 1, 0))
})

test_that("get_delta returns all-ones when custom delta_col is absent", {
  df <- data.frame(t = c(1, 2, 3), delta = c(1, 0, 1))
  expect_equal(get_delta(df, delta_col = "censor"), c(1, 1, 1))
})
