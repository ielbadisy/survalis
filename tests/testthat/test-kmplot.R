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

test_that("conf.type sets the band, log-log by default, as survfit()", {
  for (ct in c("log-log", "log", "plain")) {
    d <- kmplot(Surv(time, status) ~ trt, data = veteran, conf.int = TRUE, conf.type = ct,
                risk.table = FALSE)$data
    f <- survival::survfit(Surv(time, status) ~ trt, data = veteran, conf.type = ct)
    expect_equal(d$lo[d$time > 0], f$lower)
    expect_equal(d$hi[d$time > 0], f$upper)
  }
  d <- kmplot(Surv(time, status) ~ 1, data = veteran, risk.table = FALSE)$data
  f <- survival::survfit(Surv(time, status) ~ 1, data = veteran, conf.type = "log-log")
  expect_equal(d$lo[-1], f$lower)
})

test_that("risk table drops the duplicated time axis and scales its height", {
  skip_if_not_installed("patchwork")
  p <- kmplot(Surv(time, status) ~ 1, data = veteran)
  expect_equal(p[[1]]$labels$x, "Time")
  r <- p[[2]]
  expect_null(r$labels$x)
  expect_s3_class(r$theme$axis.text.x, "element_blank")
  expect_s3_class(r$theme$axis.line, "element_blank")
  h1 <- p$patches$layout$heights
  h4 <- kmplot(Surv(time, status) ~ celltype, data = veteran)$patches$layout$heights
  expect_lt(h1[2], h4[2])
  expect_lt(h1[2] / h1[1], 0.25)
  expect_equal(kmplot(Surv(time, status) ~ 1, data = veteran,
                      table.height = 0.5)$patches$layout$heights, c(1, 0.5))
})

test_that("model and naive curves put the linetype legend above the panel", {
  m <- fit(Surv(time, status) ~ trt + karno, data = veteran, model = "coxph")
  p <- kmplot(Surv(time, status) ~ trt, data = veteran, model = m, model.curve = "adjusted",
              risk.table = FALSE)
  lt <- p$scales$get_scales("linetype")
  expect_equal(lt$guide$params$position, "top")
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("several competing events in one call, with the naive 1 - KM", {
  d <- survival::mgus2
  d$etime <- ifelse(d$pstat == 1, d$ptime, d$futime)
  d$event <- factor(ifelse(d$pstat == 1, 1, 2 * d$death), 0:2, c("censored", "pcm", "death"))
  aj <- survival::survfit(Surv(etime, event) ~ 1, data = d)
  p <- kmplot(Surv(etime, event) ~ 1, data = d, event = c("pcm", "death"), naive = TRUE,
              risk.table = FALSE)
  x <- p$data
  expect_equal(levels(x$group), c("pcm", "death"))
  expect_equal(x$est[x$group == "death" & x$time > 0], unname(aj$pstate[, "death"]))
  expect_equal(p$labels$y, "Cumulative incidence")
  nk <- p$layers[[2]]$data
  km <- survival::survfit(Surv(etime, event == "pcm") ~ 1, data = d)
  expect_equal(nk$est[nk$group == "pcm" & nk$time > 0], 1 - km$surv)
  s <- kmplot(Surv(etime, event) ~ sex, data = d, event = c("pcm", "death"), risk.table = FALSE)
  expect_equal(levels(s$data$group), c("F, pcm", "F, death", "M, pcm", "M, death"))
  expect_error(kmplot(Surv(time, status) ~ trt, data = veteran, naive = TRUE), "competing")
  skip_if_not_installed("patchwork")
  expect_s3_class(kmplot(Surv(etime, event) ~ sex, data = d, event = c("pcm", "death")),
                  "patchwork")
})

test_that("rmst shades the area between two curves, or under one, up to tau", {
  # areas checked against the restricted means of survfit (print rmean)
  area <- function(rb, v) sum(diff(rb$time)[c(TRUE, FALSE)] * v[c(TRUE, FALSE)])
  f <- survival::survfit(Surv(time, status) ~ sex, data = survival::lung)
  s <- split(summary(f, times = 0:500, extend = TRUE)$surv, rep(1:2, each = 501))
  p <- kmplot(Surv(time, status) ~ sex, data = survival::lung, rmst = 500, risk.table = FALSE)
  rb <- p$layers[[1]]$data
  expect_equal(max(rb$time), 500)
  # lung times are integers, so the unit grid is exact; the curves cross, hence |S1 - S2|
  expect_equal(area(rb, rb$hi - rb$lo), sum(abs(s[[1]] - s[[2]])[-501]), tolerance = 1e-8)
  one <- kmplot(Surv(time, status) ~ 1, data = survival::lung, rmst = 500, risk.table = FALSE)
  rm1 <- summary(survival::survfit(Surv(time, status) ~ 1, data = survival::lung), rmean = 500)$table
  expect_equal(area(one$layers[[1]]$data, one$layers[[1]]$data$hi), unname(rm1["rmean"]),
               tolerance = 1e-8)
  expect_silent(ggplot2::ggplot_build(one))
  expect_error(kmplot(Surv(time, status) ~ trt, data = veteran, type = "risk", rmst = 200),
               "survival curves")
})
