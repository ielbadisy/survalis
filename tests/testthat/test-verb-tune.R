f <- Surv(time, status) ~ age + karno + celltype
skip_verb <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("ranger")
}

test_that("tune() returns a survalis_tune with a usable refit fit", {
  skip_verb()
  g <- expand.grid(num.trees = c(100, 200), mtry = 2, min.node.size = c(5, 15))
  t <- tune(f, veteran, "ranger", grid = g, times = c(90, 180),
            resampling = cv(3, seed = 1), metric = "ibs")
  expect_s3_class(t, "survalis_tune")
  expect_s3_class(t$fit, "survalis_fit")
  expect_equal(nrow(t$results), 4L)
  expect_true(all(c("num.trees", "mtry", "min.node.size") %in% names(t$fit$spec)))
  expect_output(print(t), "survalis_tune.*ranger")

  p <- predict(t, veteran[1:4, ], times = c(90, 180))
  expect_s3_class(p, "survmat")
  expect_equal(dim(p), c(4L, 2L))
})

test_that("tune() rejects non-tunable models and non-cv resampling", {
  skip_verb()
  expect_error(tune(f, veteran, "coxph", times = 90), "not tunable")
  expect_error(
    tune(f, veteran, "ranger", times = 90, resampling = bootstrap(5)),
    "cv\\(\\) only"
  )
})

test_that("evaluate() works on a tune result", {
  skip_verb()
  g <- expand.grid(num.trees = 100, mtry = 2, min.node.size = c(5, 15))
  t <- tune(f, veteran, "ranger", grid = g, times = c(90, 180),
            resampling = cv(2, seed = 1), metric = "ibs")
  e <- evaluate(t, times = c(90, 180), resampling = cv(2, seed = 2), metrics = "cindex")
  expect_s3_class(e, "survalis_eval")
})
