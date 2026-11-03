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

test_that("tune() needs a grid without a built-in tuner, and cv() for built-in grids", {
  skip_verb()
  expect_error(tune(f, veteran, "coxph", times = 90), "no built-in grid")
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

test_that("tune() with a partial grid gives omitted tuned parameters their defaults", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rpart")
  f <- Surv(time, status) ~ age + karno + celltype
  # rpart's built-in tuner tunes minsplit, cp and maxdepth; the grid names cp only
  t <- tune(f, veteran, "rpart", grid = data.frame(cp = c(0.001, 0.01)),
            times = c(90, 180), resampling = cv(2, seed = 1), metric = "ibs")
  expect_s3_class(t, "survalis_tune")
  expect_named(t$results, c("cp", "ibs"))
  expect_equal(nrow(t$results), 2L)
  expect_named(t$fit$spec, "cp")
  # the refit uses fit_rpart()'s defaults for the omitted parameters
  ref <- fit(f, veteran, "rpart", spec = list(cp = t$fit$spec$cp))
  expect_equal(predict(t, veteran[1:5, ], times = c(90, 180)),
               predict(ref, veteran[1:5, ], times = c(90, 180)))

  # a grid with every tuned column still goes through the built-in tuner
  full <- tune(f, veteran, "rpart",
               grid = expand.grid(minsplit = 10, cp = c(0.001, 0.01), maxdepth = 5),
               times = c(90, 180), resampling = cv(2, seed = 1), metric = "ibs")
  expect_setequal(names(full$fit$spec), c("minsplit", "cp", "maxdepth"))
})

test_that("partial grids work for every learner with a built-in tuner", {
  testthat::skip_on_cran()
  f <- Surv(time, status) ~ age + karno + celltype
  grids <- list(
    bnnsurv = data.frame(k = c(3, 5)),
    xgboost = data.frame(eta = c(0.1, 0.3)),
    coxboost = data.frame(stepno = c(20, 50)),
    ranger = data.frame(num.trees = c(50, 100)),
    rsf = data.frame(ntree = c(50, 100)),
    cforest = data.frame(ntree = c(30, 60)),
    blackboost = data.frame(mstop = c(20, 50)),
    survsvm = data.frame(gamma.mu = c(0.1, 1)),
    orsf = data.frame(n_tree = c(20, 40)),
    selectcox = data.frame(rule = c("aic", "p")),
    flexsurvreg = data.frame(dist = c("weibull", "lognormal"))
  )
  pkgs <- c(bnnsurv = NA, xgboost = "xgboost", coxboost = "CoxBoost",
            ranger = "ranger", rsf = "randomForestSRC", cforest = "partykit",
            blackboost = "mboost", survsvm = NA, orsf = "aorsf",
            selectcox = "pec", flexsurvreg = "flexsurv")
  for (id in names(grids)) {
    if (!is.na(pkgs[[id]]) && !requireNamespace(pkgs[[id]], quietly = TRUE)) next
    t <- tune(f, veteran, id, grid = grids[[id]], times = c(90, 180),
              resampling = cv(2, seed = 1), metric = "ibs")
    expect_s3_class(t, "survalis_tune")
    expect_equal(nrow(t$results), 2L, info = id)
    expect_named(t$fit$spec, names(grids[[id]]), info = id)
  }
})
