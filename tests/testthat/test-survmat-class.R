test_that("new_survmat sets the contract", {
  m <- matrix(c(0.9, 0.8, 0.7, 0.6), nrow = 2, byrow = TRUE)
  s <- new_survmat(m, times = c(30, 90))
  expect_s3_class(s, "survmat")
  expect_true(is.matrix(s))
  expect_identical(attr(s, "times"), c(30, 90))
  expect_identical(colnames(s), c("t=30", "t=90"))
  expect_identical(attr(s, "quantity"), "survival")
})

test_that("validate_survmat enforces structure", {
  ok <- new_survmat(matrix(c(0.9, 0.5), nrow = 1), times = c(10, 20))
  expect_invisible(validate_survmat(ok))

  bad_na <- ok; bad_na[1, 1] <- NA_real_
  expect_error(validate_survmat(bad_na), "missing values")

  bad_times <- ok; attr(bad_times, "times") <- c(20, 10)
  expect_error(validate_survmat(bad_times), "strictly increasing")

  bad_range <- new_survmat(matrix(c(0.9, 1.5), nrow = 1), times = c(10, 20))
  expect_error(validate_survmat(bad_range), "outside \\[0, 1\\]")

  chf <- new_survmat(matrix(c(0.1, 2.3), nrow = 1), times = c(10, 20), quantity = "chf")
  expect_invisible(validate_survmat(chf))
})

test_that("as_survmat coerces data frames and vectors", {
  df <- data.frame(`t=1` = c(0.9, 0.8), `t=2` = c(0.7, 0.6), check.names = FALSE)
  s <- as_survmat(df, times = c(1, 2))
  expect_s3_class(s, "survmat")
  expect_equal(unclass(as.matrix(s)), as.matrix(df), ignore_attr = TRUE)

  v <- as_survmat(c(0.9, 0.5, 0.2), times = c(1, 2, 3))
  expect_equal(nrow(v), 1L)
})

test_that("survmat_times reads attribute then falls back to colnames", {
  s <- new_survmat(matrix(0.5, 1, 2), times = c(5, 9))
  expect_identical(survmat_times(s), c(5, 9))
  df <- data.frame(`t=5` = 0.5, `t=9` = 0.4, check.names = FALSE)
  expect_identical(survmat_times(df), c(5, 9))
})

test_that("as.data.frame.survmat drops class and attributes", {
  s <- new_survmat(matrix(c(0.9, 0.8), 1, 2), times = c(1, 2))
  d <- as.data.frame(s)
  expect_s3_class(d, "data.frame")
  expect_null(attr(d, "times"))
  expect_identical(names(d), c("t=1", "t=2"))
})
