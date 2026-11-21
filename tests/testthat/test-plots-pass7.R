test_that("PDP plot: sentence-case y label and readable ICE/PDP legend keys", {
  m <- fit(Surv(time, status) ~ age + karno, veteran, "coxph")
  p <- plot(interpret(m, "pdp", feature = "karno", times = c(90, 180)))
  expect_equal(p$labels$y, "Survival probability")
  expect_equal(levels(p$data$key), c("Individual (ICE)", "Average (PDP)"))
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("kmplot moves the log-rank p off curves that cross the default spot", {
  set.seed(3); g <- factor(rep(1:4, each = 100)); t <- rexp(400, c(3, 1.5, 0.3, 0.15)[g])
  d <- data.frame(t = pmin(t, 12), e = as.integer(t < 12), g = g)
  pv <- function(p) { l <- p[[1]]$layers; l[[length(l)]]$data }
  a <- pv(kmplot(Surv(t, e) ~ g, data = d))
  expect_gt(a$x, 6); expect_gt(a$y, 0.08)
  # default spot kept when it is free
  b <- pv(kmplot(Surv(time, status) ~ sex, data = survival::lung))
  expect_equal(b$y, 0.08); expect_lt(b$x, 50)
})
