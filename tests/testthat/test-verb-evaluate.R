f <- Surv(time, status) ~ age + karno + celltype
skip_verb <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("pec")
}

test_that("evaluate() on a formula returns a survalis_eval with a summary", {
  skip_verb()
  e <- evaluate(f, veteran, "coxph", times = c(90, 180, 365), resampling = cv(3, seed = 1))
  expect_s3_class(e, "survalis_eval")
  expect_setequal(e$summary$metric, c("cindex", "ibs"))
  expect_equal(e$times, c(90, 180, 365))
  expect_equal(nrow(e$per_split), 3L * 2L)
  expect_true(all(e$summary$mean >= 0))
  expect_output(print(e), "survalis_eval.*coxph")
})

test_that("evaluate() on a fit reuses the stored recipe", {
  skip_verb()
  m <- fit(f, veteran, "coxph")
  e1 <- evaluate(m, times = c(90, 180), resampling = cv(3, seed = 42))
  e2 <- evaluate(f, veteran, "coxph", times = c(90, 180), resampling = cv(3, seed = 42))
  expect_equal(e1$per_split$value, e2$per_split$value)
})

test_that("evaluate() honors the resampling spec (holdout, repeats)", {
  skip_verb()
  eh <- evaluate(f, veteran, "coxph", times = 180, resampling = holdout(0.7, seed = 1),
                 metrics = "cindex")
  expect_equal(nrow(eh$per_split), 1L)

  er <- evaluate(f, veteran, "coxph", times = c(90, 180), resampling = cv(3, repeats = 2, seed = 1))
  expect_equal(length(unique(er$per_split$id)), 6L)
})

test_that("evaluate() default grid is recorded", {
  skip_verb()
  e <- evaluate(f, veteran, "coxph", resampling = cv(3, seed = 1))
  expect_equal(e$times, default_times(veteran$time, veteran$status))
})

test_that("bootstrap .632, .632+ and optimism estimators follow their definitions", {
  f <- Surv(time, status) ~ age + karno + celltype
  tm <- c(30, 90, 180)
  ev <- function(est) evaluate(f, veteran, model = "coxph", times = tm, metrics = c("cindex", "ibs", "brier"),
                               resampling = bootstrap(8, seed = 1, estimator = est))
  oob <- ev("oob"); e <- ev(".632+")
  expect_null(oob$bootstrap)
  b <- e$bootstrap
  expect_equal(b$oob, oob$summary$mean[match(b$metric, oob$summary$metric)])
  m <- fit(f, veteran, model = "coxph")
  S <- as.matrix(predict(m, veteran, times = tm)); y <- Surv(veteran$time, veteran$status); n <- nrow(S)
  # apparent score and the exact n^2-pair no-information Brier score
  expect_equal(b$apparent[b$metric == "cindex"], unname(cindex_survmat(y, S, max(tm))))
  ni <- vapply(1:3, function(k) mean(vapply(1:n, function(j) survalis:::brier(y, rep(S[j, k], n), tm[k]), 0)), 0)
  expect_equal(b$noinf[match(paste0("brier@", tm), b$metric)], ni)
  expect_equal(b$noinf[b$metric == "ibs"], sum(ni[1:2] * diff(tm)) / diff(range(tm)))
  expect_lt(abs(b$noinf[b$metric == "cindex"] - 0.5), 0.05)
  # the .632+ formula, on the loss scale for the C-index
  sg <- ifelse(b$metric == "cindex", -1, 1)
  R <- with(b, pmin(pmax((sg * oob - sg * apparent) / (sg * noinf - sg * apparent), 0), 1))
  expect_equal(b$R, R)
  expect_equal(b[[".632+"]], with(b, (1 - w) * apparent + w * oob))
  expect_equal(b[[".632"]], with(b, 0.368 * apparent + 0.632 * oob))
  expect_equal(e$summary$mean, b[[".632+"]][match(e$summary$metric, b$metric)])
  expect_true(all(is.na(e$summary$se)))
  # optimism: apparent plus mean over resamples of (score on all rows - score on own resample)
  p <- e$per_split[e$per_split$metric == "ibs", ]
  expect_equal(b$optimism[b$metric == "ibs"], b$apparent[b$metric == "ibs"] + mean(p$on_full - p$on_inbag))
  expect_output(print(e), "\\.632\\+ estimate")
})

test_that("bootstrap() validates the estimator and defaults to out-of-bag", {
  expect_error(bootstrap(5, estimator = "x"))
  expect_equal(bootstrap(5)$params$estimator, "oob")
})
