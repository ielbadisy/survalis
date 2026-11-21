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
