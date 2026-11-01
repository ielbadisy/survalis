test_that("dwscore closed forms match numerical integration", {
  r <- c(0.03, 0.1, 0.3, 0.6)
  num <- function(w) list(
    l0 = sapply(r, function(x) stats::integrate(function(c) w(c) * c, 0, x)$value),
    l1 = sapply(r, function(x) stats::integrate(function(c) w(c) * (1 - c), x, 1)$value)
  )
  chk <- function(wt, w, tol = 1e-6) {
    n <- num(w)
    expect_equal(wt$l0(r), n$l0, tolerance = tol)
    expect_equal(wt$l1(r), n$l1, tolerance = tol)
  }
  chk(.dw_weight("brier"), function(c) 2 + 0 * c)
  chk(.dw_weight("survacc"), function(c) (c^2 + (1 - c)^2)^(-3 / 2))
  chk(.dw_weight("beta", shape = c(2, 18)), function(c) stats::dbeta(c, 2, 18))
  chk(.dw_weight("range", thresholds = c(0.05, 0.2)),
      function(c) (c >= 0.05 & c <= 0.2) / 0.15, tol = 1e-4)
})

test_that("dwscore returns one loss per time", {
  set.seed(1)
  idx <- sample(nrow(veteran), 90)
  m <- fit(Surv(time, status) ~ age + karno, veteran[idx, ], "coxph")
  te <- veteran[-idx, ]
  S <- predict(m, te, times = c(60, 120))
  y <- Surv(te$time, te$status)
  out <- dwscore(S, y, weight = "brier")
  expect_s3_class(out, "data.frame")
  expect_equal(out$time, c(60, 120))
  expect_true(all(out$loss >= 0))
  # the "brier" metric and dwscore(weight = "brier") are one estimator
  b <- vapply(1:2, function(k) unname(brier(y, S[, k], c(60, 120)[k])), 0)
  expect_equal(b, round(out$loss, 6))
})

test_that("decision-weighted loss equals weighted area under the IPCW net benefit curve", {
  set.seed(2)
  idx <- sample(nrow(veteran), 90)
  m <- fit(Surv(time, status) ~ age + karno, veteran[idx, ], "coxph")
  te <- veteran[-idx, ]
  S <- predict(m, te, times = 100)
  y <- Surv(te$time, te$status)
  r <- 1 - S[, 1]
  # IPCW weights as in .dw_ipcw
  km <- survival::survfit(survival::Surv(te$time, 1 - te$status) ~ 1)
  Gl <- function(q) { i <- findInterval(q, km$time, left.open = TRUE)
    o <- rep(1, length(q)); o[i > 0] <- km$surv[i[i > 0]]; o }
  Gt <- summary(km, times = 100, extend = TRUE)$surv
  w0 <- (te$time > 100) / Gt
  ev <- te$time <= 100 & te$status == 1
  w1 <- numeric(nrow(te)); w1[ev] <- 1 / Gl(te$time[ev])
  n <- nrow(te); pibar <- sum(w1) / n
  nb <- function(c) sum(w1[r > c]) / n - sum(w0[r > c]) / n * c / (1 - c)
  # range weight: L = int_{cl}^{cu} h(c) (pibar - NB(c)) dc, h = (1 - c) / (cu - cl)
  cs <- sort(unique(c(0.05, 0.2, r[r > 0.05 & r < 0.2])))
  area <- sum(vapply(seq_len(length(cs) - 1), function(k) {
    stats::integrate(function(c) vapply(c, function(cc) (1 - cc) / 0.15 * (pibar - nb(cc)), 0),
                     cs[k], cs[k + 1])$value
  }, 0))
  expect_equal(dwscore(S, y, weight = "range", thresholds = c(0.05, 0.2))$loss, area,
               tolerance = 1e-6)
})

test_that("survacc and range ids work in evaluate, compare and tune", {
  S <- predict(fit(Surv(time, status) ~ age + karno, veteran, "coxph"), veteran, times = 100)
  y <- Surv(veteran$time, veteran$status)
  ev <- evaluate(S, y = y, metrics = c("survacc", "range:0.05-0.20"))
  d <- dwscore(S, y, weight = "survacc")$loss
  expect_equal(ev$summary$mean[1], round(1 - d, 3))
  expect_true(ev$summary$mean[2] >= 0)
  expect_error(evaluate(S, y = y, metrics = "range:0.3-0.1"), "lower < upper")

  cmp <- compare(Surv(time, status) ~ age + karno, veteran, models = c("coxph", "aalen"),
                 times = 100, resampling = cv(3, seed = 1),
                 metrics = c("survacc", "range:0.05-0.20"))
  expect_setequal(unique(cmp$leaderboard$metric), c("survacc", "range:0.05-0.20"))
  expect_true(survalis:::.metric_maximize("survacc"))
  expect_false(survalis:::.metric_maximize("range:0.05-0.20"))

  tu <- tune(Surv(time, status) ~ age + karno, veteran, "glmnet",
             grid = data.frame(alpha = c(0, 1)), times = 100,
             resampling = cv(3, seed = 1), metric = "range:0.05-0.20")
  expect_s3_class(tu, "survalis_tune")
})
