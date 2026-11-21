# Pass-4 book issue: compare/evaluate plots squashed metrics onto one shared value axis.
set.seed(1)
ps <- data.frame(model = rep(c("coxph", "ranger"), each = 20), metric = rep(c("cindex", "ibs"), each = 10),
                 value = c(runif(10, 0.66, 0.70), runif(10, 0.15, 0.17), runif(10, 0.67, 0.71),
                           runif(10, 0.16, 0.18)))
lb <- aggregate(value ~ model + metric, ps, mean); names(lb)[3] <- "mean"
cmp <- structure(list(per_split = ps, leaderboard = lb, metrics = c("cindex", "ibs")),
                 class = "survalis_compare")
xr <- function(p) lapply(ggplot2::ggplot_build(p)$layout$panel_params, function(q) q$x.range)

test_that("plot.survalis_compare gives each metric its own value axis", {
  r <- xr(plot(cmp))
  expect_length(r, 2L)
  expect_true(r[[1]][1] > 0.6 && r[[2]][2] < 0.2)
  expect_equal(xr(plot(cmp, scales = "fixed"))[[1]], xr(plot(cmp, scales = "fixed"))[[2]])
})

test_that("plot.survalis_compare labels direction and puts the best model on top", {
  p <- plot(cmp, title = NULL)
  expect_equal(levels(p$data$metric), c("cindex (higher is better)", "ibs (lower is better)"))
  expect_equal(levels(p$data$model), c("coxph", "ranger"))
  expect_null(p$labels$title)
})

test_that("plot.survalis_eval gives each metric its own value axis", {
  r <- xr(cv_plot(ps[ps$model == "coxph", ]))
  expect_length(r, 2L)
  expect_true(r[[1]][1] > 0.6 && r[[2]][2] < 0.2)
})

test_that("evaluate() scores Brier at several times and plots it against time", {
  m <- fit(Surv(time, status) ~ age + karno, veteran, "coxph")
  tt <- c(30, 90, 180)
  ev <- evaluate(m, times = tt, resampling = cv(3, seed = 1), metrics = c("cindex", "brier"))
  expect_equal(ev$summary$metric, c("brier@30", "brier@90", "brier@180", "cindex"))
  S <- predict(m, veteran, times = tt)
  y <- Surv(veteran$time, veteran$status)
  one <- evaluate(S, y = y, metrics = "brier")
  ref <- vapply(seq_along(tt), function(k) unname(brier(y, S[, k], tt[k])), numeric(1))
  expect_equal(one$summary$mean, ref)
  # single time: unchanged metric name and value
  s1 <- evaluate(predict(m, veteran, times = 90), y = y, metrics = "brier")
  expect_equal(s1$summary$metric, "brier")
  expect_equal(s1$summary$mean, ref[2])
  expect_error(tune(Surv(time, status) ~ age, veteran, "glmnet", times = tt, metric = "brier"),
               "single evaluation time")
  p <- plot(ev, type = "time", reference = one, labels = c("CV", "Apparent"), xlab = "Days")
  expect_equal(sort(unique(p$data$time)), tt)
  expect_equal(levels(p$data$curve), c("CV", "Apparent"))
  expect_equal(p$labels$y, "Brier score")
  expect_silent(ggplot2::ggplot_build(p))
  expect_error(plot(evaluate(S, y = y, metrics = "ibs"), type = "time"), "several")
})

test_that("evaluate plot has no default title and spaced value axes", {
  m <- fit(Surv(time, status) ~ age + karno, veteran, "coxph")
  p <- plot(evaluate(m, times = 90, resampling = cv(3, seed = 1),
                     metrics = c("cindex", "auc", "brier")))
  expect_null(p$labels$title)
  expect_equal(p$theme$panel.spacing.x, ggplot2::unit(1.5, "lines"))
  expect_equal(plot(evaluate(m, times = 90, resampling = cv(3, seed = 1)), title = "A")$labels$title, "A")
})

test_that("ROC, DCA and PDP take colours and time labels, and recolour without a message", {
  m <- fit(Surv(time, status) ~ age + karno, veteran, "coxph")
  d <- plot(interpret(m, "dca", eval_time = 90))
  expect_null(d$labels$title)
  rc <- ggplot2::scale_colour_manual(values = c("red", "green", "blue"))
  expect_no_message(b <- ggplot2::ggplot_build(d + rc))
  expect_setequal(unique(b$data[[1]]$colour), c("red", "green", "blue"))
  b <- ggplot2::ggplot_build(plot(interpret(m, "dca", eval_time = 90),
                                 colors = c("black", "grey50", "grey80")))
  expect_setequal(unique(b$data[[1]]$colour), c("black", "grey50", "grey80"))
  r <- plot(interpret(m, "roc", eval_time = 90), colors = "red")
  expect_equal(r$layers[[2]]$aes_params$colour, "red")
  pd <- interpret(m, "pdp", feature = "karno", times = c(90, 180))
  p <- plot(pd, time.labels = function(t) paste(t, "days"))
  expect_equal(p$facet$params$labeller(data.frame(time = c("90", "180")))[[1]], c("90 days", "180 days"))
  expect_equal(.time_lab("days", 90), "t = 90 days")
  expect_equal(.time_lab(NULL, 90), "t = 90")
  expect_silent(ggplot2::ggplot_build(p))
  a <- plot(interpret(m, "ale", feature = "karno", times = c(180, 90)), time.labels = "days")
  expect_equal(levels(a$data$time), c("t = 90 days", "t = 180 days"))
})

test_that("ALE keeps empty bins in place (discrete feature, interpolated grid)", {
  m <- fit(Surv(time, status) ~ age + karno, veteran, "coxph")
  tt <- c(90, 180)
  a <- compute_ale(m, veteran, "karno", times = tt, grid.size = 20)
  g <- a$ale$feature_value
  K <- length(g) - 1
  b <- cut(veteran$karno, g, include.lowest = TRUE, labels = FALSE)
  expect_true(any(tabulate(b, K) == 0))
  pr <- function(v) { d <- veteran; d$karno <- v; as.matrix(predict(m, d, times = tt)) }
  dif <- pr(g[b + 1]) - pr(g[b])
  loc <- t(vapply(seq_len(K), function(k) if (any(b == k)) colMeans(dif[b == k, , drop = FALSE]) else
    c(0, 0), numeric(2)))
  ale <- rbind(0, apply(loc, 2, cumsum))
  nk <- tabulate(b, K)
  mid <- (ale[-1, ] + ale[-(K + 1), ]) / 2
  ref <- sweep(ale, 2, colSums(nk * mid) / sum(nk))
  expect_equal(unname(as.matrix(a$ale[, -1])), unname(ref), tolerance = 1e-10)
})
