f <- Surv(time, status) ~ age + karno + celltype
skip_verb <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")
  testthat::skip_if_not_installed("pec")
}

test_that("estimate() returns a survalis_estimate with a bootstrap CI", {
  skip_verb()
  e <- estimate(f, veteran, "coxph", estimand = "RMST_diff", treatment = "trt",
                tau = 300, R = 30, refit = FALSE, seed = 1)
  expect_s3_class(e, "survalis_estimate")
  expect_equal(e$estimand, "RMST_diff")
  expect_identical(e$levels, sort(unique(veteran$trt)))
  expect_true(e$conf_low <= e$estimate && e$estimate <= e$conf_high)
  expect_equal(nrow(e$curves), length(e$curves$time))
  expect_output(print(e), "RMST difference at tau = 300")
})

test_that("estimate() supports the four estimands", {
  skip_verb()
  for (q in c("RMST_diff", "RMST_ratio", "survival_diff", "survival_ratio")) {
    e <- estimate(f, veteran, "coxph", estimand = q, treatment = "trt",
                  tau = 200, R = 15, refit = FALSE, seed = 1)
    expect_true(is.finite(e$estimate), info = q)
  }
})

test_that("estimate() validates treatment, tau and grid coverage", {
  skip_verb()
  expect_error(estimate(f, veteran, "coxph", treatment = "nope", tau = 100),
               "must name a column")
  vv <- veteran; vv$grp <- vv$age  # continuous -> not binary
  expect_error(estimate(Surv(time, status) ~ karno, vv, "coxph", treatment = "grp", tau = 100),
               "exactly two distinct values")
  expect_error(estimate(f, veteran, "coxph", treatment = "trt", tau = 100,
                        times = c(10, 20, 30)),
               "must reach `tau`")
})

test_that("estimate(refit = TRUE) refits per replicate", {
  skip_verb()
  e <- estimate(f, veteran, "coxph", treatment = "trt", tau = 200, R = 20,
                refit = TRUE, seed = 1)
  expect_true(e$refit)
  expect_length(e$boot, 20L)
  expect_s3_class(plot(e), "ggplot")
})

test_that("as.data.frame() of an estimate has plain column names", {
  skip_verb()
  e <- estimate(f, veteran, "coxph", estimand = "RMST_diff", treatment = "trt",
                tau = 300, R = 5, refit = FALSE, seed = 1)
  d <- as.data.frame(e)
  expect_false(any(grepl("[.=]", names(d))))
  expect_true(all(c("estimate", "conf_low", "conf_high") %in% names(d)))
})

test_that("estimate() works with tree learners (numeric or factor treatment)", {
  skip_on_cran()
  d <- survival::veteran; d$trt <- d$trt - 1
  tt <- c(50, 100, 200)
  e <- estimate(Surv(time, status) ~ age + karno, d, model = "ranger", treatment = "trt",
                tau = 200, times = tt, R = 2, seed = 1, spec = list(num.trees = 50))
  # point estimate = G-computation with the main-effects ranger fit
  lrn <- .resolve_learner("ranger")
  set.seed(1)
  m <- lrn$fit(Surv(time, status) ~ trt + age + karno, d, num.trees = 50)
  S <- function(a) { x <- d; x$trt <- a; colMeans(as.matrix(lrn$predict(m, x, tt))) }
  ref <- survmat_to_rmst(matrix(S(1), 1), tt, 200) - survmat_to_rmst(matrix(S(0), 1), tt, 200)
  expect_equal(e$estimate, ref, tolerance = 1e-10)
  d$trt <- factor(d$trt); d$celltype <- factor(d$celltype)
  for (l in c("ranger", "rpart", "orsf")) {
    e <- estimate(Surv(time, status) ~ age + karno + celltype, d, model = l, treatment = "trt",
                  tau = 200, times = tt, R = 2, seed = 1)
    expect_true(is.finite(e$estimate))
  }
})
