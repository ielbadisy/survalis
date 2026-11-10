test_that("kmplot() draws survival curves starting at 1", {
  p <- kmplot(Surv(time, status) ~ trt, data = veteran, risk.table = FALSE)
  expect_s3_class(p, "ggplot")
  d <- p$data
  expect_equal(d$est[d$time == 0], c(1, 1))
  fit <- survival::survfit(Surv(time, status) ~ trt, data = veteran)
  expect_equal(d$est[d$time > 0], fit$surv)
  expect_equal(sum(d$cens), sum(fit$n.censor))
})

test_that("type = 'risk' is one minus the Kaplan-Meier estimate", {
  s <- kmplot(Surv(time, status) ~ trt, data = veteran, risk.table = FALSE)$data
  r <- kmplot(Surv(time, status) ~ trt, data = veteran, type = "risk",
              risk.table = FALSE)$data
  expect_equal(r$est, 1 - s$est)
  expect_equal(r$lo, 1 - s$hi)
  expect_equal(r$hi, 1 - s$lo)
  expect_true(all(r$lo <= r$est + 1e-12, na.rm = TRUE))
})

test_that("kmplot() p-value matches survdiff and accepts a string", {
  sd <- survival::survdiff(Surv(time, status) ~ trt, data = veteran)
  p <- stats::pchisq(sd$chisq, 1, lower.tail = FALSE)
  expect_equal(.kmplot_pval(TRUE, Surv(time, status) ~ trt, veteran, 2),
               sprintf("Log-rank p = %.3f", p))
  expect_equal(.kmplot_pval("p = 0.5", NULL, NULL, 2), "p = 0.5")
  expect_null(.kmplot_pval(TRUE, Surv(time, status) ~ 1, veteran, 1))
})

test_that("kmplot() handles one group, several factors and the risk table", {
  p1 <- kmplot(Surv(time, status) ~ 1, data = veteran, risk.table = FALSE)
  expect_equal(levels(p1$data$group), "All")
  expect_equal(.kmplot_labels("trt=1, celltype=large"), "1, large")
  skip_if_not_installed("patchwork")
  p <- kmplot(Surv(time, status) ~ celltype, data = veteran, conf.int = TRUE,
              breaks = seq(0, 1000, 250))
  expect_s3_class(p, "patchwork")
  expect_silent(ggplot2::ggplot_build(p[[1]]))
  expect_error(kmplot(Surv(time, status) ~ celltype, data = veteran,
                      colors = "red"), "at least 4")
})

test_that("type = 'cumhaz' and 'cloglog' transform the Kaplan-Meier estimate", {
  fit <- survival::survfit(Surv(time, status) ~ trt, data = veteran)
  h <- kmplot(Surv(time, status) ~ trt, data = veteran, type = "cumhaz",
              risk.table = FALSE)$data
  ref <- -log(fit$surv)
  ref[!is.finite(ref)] <- NA
  expect_equal(h$est[h$time > 0], ref)
  expect_true(all(h$lo <= h$est + 1e-12, na.rm = TRUE))
  g <- kmplot(Surv(time, status) ~ trt, data = veteran, type = "cloglog")
  expect_s3_class(g, "ggplot")
  expect_false(inherits(g, "patchwork"))
  ok <- fit$surv < 1 & fit$surv > 0
  expect_equal(sort(g$data$est), sort(log(-log(fit$surv[ok]))))
  expect_silent(ggplot2::ggplot_build(g))
})

test_that("competing risks draw the Aalen-Johansen cumulative incidence", {
  d <- survival::mgus2
  d$etime <- ifelse(d$pstat == 1, d$ptime, d$futime)
  d$event <- factor(ifelse(d$pstat == 1, 1, 2 * d$death), 0:2,
                    c("censored", "pcm", "death"))
  aj <- survival::survfit(Surv(etime, event) ~ sex, data = d)
  p <- kmplot(Surv(etime, event) ~ sex, data = d, event = "death",
              risk.table = FALSE)
  x <- p$data
  expect_equal(x$est[x$time > 0], unname(aj$pstate[, "death"]))
  expect_equal(x$est[x$time == 0], c(0, 0))
  expect_equal(p$labels$y, "Cumulative incidence of death")
  expect_error(kmplot(Surv(etime, event) ~ sex, data = d, type = "survival"),
               "Competing-risks")
  expect_error(kmplot(Surv(etime, event) ~ sex, data = d, event = "x"), "must be one of")
  skip_if_not_installed("patchwork")
  expect_s3_class(kmplot(Surv(etime, event) ~ sex, data = d), "patchwork")
})

test_that("model curves average the predictions within or across groups", {
  m <- fit(Surv(time, status) ~ trt + karno, data = veteran, model = "coxph")
  p <- kmplot(Surv(time, status) ~ trt, data = veteran, model = m,
              risk.table = FALSE)
  mc <- p$layers[[2]]$data
  expect_equal(levels(mc$group), c("1", "2"))
  t1 <- mc$time[mc$time > 0][10]
  s1 <- mean(predict(m, veteran[veteran$trt == 1, ], times = t1))
  expect_equal(mc$est[mc$group == "1" & mc$time == t1], s1)

  a <- kmplot(Surv(time, status) ~ trt, data = veteran, model = m,
              model.curve = "adjusted", risk.table = FALSE)
  ac <- a$layers[[2]]$data
  nd <- veteran
  nd$trt <- 2
  expect_equal(ac$est[ac$group == "2" & ac$time == t1], mean(predict(m, nd, times = t1)))
  expect_error(kmplot(Surv(time, status) ~ trt + celltype, data = veteran, model = m,
                      model.curve = "adjusted"), "exactly one")
  expect_error(kmplot(Surv(time, status) ~ trt, data = veteran, model = 1), "fit\\(\\)")
})
