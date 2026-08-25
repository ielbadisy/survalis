test_that("compute_shap_matrix returns one time-integrated phi per observation/feature", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")

  mod <- fit_coxph(Surv(time, status) ~ age + karno + diagtime, data = veteran)
  sm <- compute_shap_matrix(mod, newdata = veteran[1:6, ], baseline_data = veteran,
                            times = c(50, 100, 200), sample.size = 5)

  expect_true(all(c("observation", "feature", "phi", "raw_value") %in% colnames(sm)))
  expect_equal(nrow(sm), 6 * length(setdiff(colnames(veteran[1:6, ]), c("time", "status", "event"))))
  expect_true(all(is.finite(sm$phi)))
})

test_that("plot_shap_beeswarm returns a ggplot and respects top_n", {
  testthat::skip_on_cran()

  mod <- fit_coxph(Surv(time, status) ~ age + karno + diagtime, data = veteran)
  sm <- compute_shap_matrix(mod, newdata = veteran[1:6, ], baseline_data = veteran,
                            times = c(50, 100, 200), sample.size = 5)

  p_all <- plot_shap_beeswarm(sm)
  expect_s3_class(p_all, "ggplot")

  p_top2 <- plot_shap_beeswarm(sm, top_n = 2)
  expect_equal(nlevels(p_top2$data$feature), 2)
})
