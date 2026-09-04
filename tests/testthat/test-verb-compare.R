f <- Surv(time, status) ~ age + karno + celltype
skip_verb <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("ranger")
  testthat::skip_if_not_installed("pec")
}

test_that("compare() evaluates every model on shared folds", {
  skip_verb()
  cmp <- compare(f, veteran, models = c("coxph", "ranger"),
                 specs = list(ranger = list(num.trees = 100)),
                 times = c(90, 180, 365), resampling = cv(3, seed = 1),
                 metrics = c("cindex", "ibs"))
  expect_s3_class(cmp, "survalis_compare")
  expect_setequal(unique(cmp$per_split$model), c("coxph", "ranger"))
  expect_equal(nrow(cmp$leaderboard), 4L)
  # paired: same fold ids across models
  ids_cox <- sort(unique(cmp$per_split$id[cmp$per_split$model == "coxph"]))
  ids_rgr <- sort(unique(cmp$per_split$id[cmp$per_split$model == "ranger"]))
  expect_equal(ids_cox, ids_rgr)
  expect_output(print(cmp), "survalis_compare.*2 models")
})

test_that("compare() needs >= 2 known models", {
  skip_verb()
  expect_error(compare(f, veteran, models = "coxph", times = 90), "at least two")
  expect_error(compare(f, veteran, models = c("coxph", "nope"), times = 90), "Unknown model")
})

test_that("compare(tune = TRUE) messages about non-nested tuning", {
  skip_verb()
  expect_message(
    compare(f, veteran, models = c("coxph", "ranger"), times = c(90, 180),
            resampling = cv(2, seed = 1), metrics = "cindex", tune = TRUE),
    "non-nested"
  )
})

test_that("summary() returns a wide model-by-metric table and plot() a ggplot", {
  skip_verb()
  cmp <- compare(f, veteran, models = c("coxph", "ranger"),
                 times = c(90, 180), resampling = cv(2, seed = 1), metrics = "cindex")
  s <- summary(cmp)
  expect_true("model" %in% names(s))
  expect_s3_class(plot(cmp), "ggplot")
})
