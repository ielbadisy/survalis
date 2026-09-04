test_that("default_times returns a positive increasing grid", {
  g <- default_times(veteran$time, veteran$status)
  expect_true(is.numeric(g))
  expect_true(all(g > 0))
  expect_false(is.unsorted(g, strictly = TRUE))
  expect_lte(length(g), 50L)
})

test_that("default_times respects n and uses unique event times when few", {
  g10 <- default_times(veteran$time, veteran$status, n = 10)
  expect_lte(length(g10), 10L)

  few <- default_times(c(5, 10, 40, 90), c(1, 0, 1, 1), n = 50)
  expect_equal(few, c(5, 40, 90))
})

test_that("default_times falls back when no events", {
  g <- default_times(c(3, 7, 11), c(0, 0, 0))
  expect_equal(g, c(3, 7, 11))
})

test_that(".resolve_times validates explicit input and defers otherwise", {
  expect_equal(.resolve_times(c(90, 30, 30), NULL), c(30, 90))
  expect_error(.resolve_times(c(-1, 5), NULL), "finite and positive")
  g <- .resolve_times(NULL, veteran$time, veteran$status, n = 8)
  expect_lte(length(g), 8L)
})
