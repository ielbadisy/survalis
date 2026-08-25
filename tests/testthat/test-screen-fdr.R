test_that("screen_fdr selects covariates and controls FDR via BH", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")

  sel <- screen_fdr(Surv(time, status) ~ age + karno + diagtime + prior,
                    data = veteran, alpha = 0.2, minscreen = 1)
  tbl <- attr(sel, "screen_table")

  expect_true(is.character(sel))
  expect_true(is.data.frame(tbl))
  expect_setequal(tbl$feature, c("age", "karno", "diagtime", "prior"))
  expect_true("karno" %in% sel)
  expect_true(all(tbl$q_value[!is.na(tbl$q_value)] >= tbl$p_value[!is.na(tbl$p_value)] - 1e-8))
})

test_that("screen_fdr backfills to minscreen when too few pass the FDR threshold", {
  testthat::skip_on_cran()

  sel <- screen_fdr(Surv(time, status) ~ age + karno + diagtime + prior,
                    data = veteran, alpha = 1e-6, minscreen = 3)
  expect_length(sel, 3)
})

test_that("screen_fdr validates inputs", {
  expect_error(screen_fdr(Surv(time, status) ~ age, data = veteran, alpha = 0), "alpha")
  expect_error(screen_fdr(Surv(time, status) ~ age, data = veteran, minscreen = -1), "minscreen")
})
