f <- Surv(time, status) ~ age + karno + celltype

skip_verb <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("pec")
}


test_that("predict() returns a survmat for matrix types and records the grid", {
  skip_verb()
  m <- fit(f, veteran, "coxph")
  S <- predict(m, veteran[1:5, ], times = c(90, 180, 365))
  expect_s3_class(S, "survmat")
  expect_identical(attr(S, "times"), c(90, 180, 365))
  expect_equal(dim(S), c(5L, 3L))
  expect_true(all(S >= 0 & S <= 1))

  R <- predict(m, veteran[1:5, ], times = c(90, 180, 365), type = "risk")
  expect_identical(attr(R, "quantity"), "risk")
  expect_equal(unclass(as.matrix(R)), 1 - unclass(as.matrix(S)), ignore_attr = TRUE)

  H <- predict(m, veteran[1:5, ], times = c(90, 180, 365), type = "chf")
  expect_identical(attr(H, "quantity"), "chf")
})

test_that("predict() default grid is used and recorded when times is NULL", {
  skip_verb()
  m <- fit(f, veteran, "coxph")
  S <- predict(m, veteran[1:3, ])
  expect_true(length(attr(S, "times")) > 1L)
  expect_equal(attr(S, "times"), default_times(veteran$time, veteran$status))
})

test_that("predict() scalar types and their required args", {
  skip_verb()
  m <- fit(f, veteran, "coxph")
  r <- predict(m, veteran[1:5, ], times = seq(30, 360, 30), type = "rmst", tau = 300)
  expect_length(r, 5L)
  expect_true(all(r > 0))

  expect_error(predict(m, veteran[1:5, ], type = "quantile"), "requires `probs`")
  q <- predict(m, veteran[1:5, ], times = seq(30, 600, 30), type = "quantile", probs = 0.5)
  med <- predict(m, veteran[1:5, ], times = seq(30, 600, 30), type = "median")
  expect_equal(q, med)

  qq <- predict(m, veteran[1:5, ], times = seq(30, 600, 30), type = "quantile",
                probs = c(0.25, 0.5))
  expect_equal(dim(qq), c(5L, 2L))
})

test_that("predict() rejects newdata with missing model values", {
  skip_verb()
  m <- fit(f, veteran, "coxph")
  nd <- veteran[1:5, ]
  nd$age[2] <- NA
  expect_error(predict(m, nd, times = 90), "missing values in model variable")
})
