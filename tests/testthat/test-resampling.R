test_that("constructors validate and carry parameters", {
  expect_s3_class(cv(), "survalis_resampling")
  expect_equal(cv(10, repeats = 2)$params$v, 10L)
  expect_equal(cv(10, repeats = 2)$params$repeats, 2L)
  expect_error(cv(1), ">= 2")
  expect_error(holdout(prop = 1), "\\(0, 1\\)")
  expect_error(group_cv(), "single column name")
  expect_error(bootstrap(0), ">= 1")
})

test_that("print is compact and informative", {
  expect_output(print(cv(5, seed = 1)), "5-fold CV.*seed = 1")
  expect_output(print(bootstrap(30)), "bootstrap \\(30 resamples")
})

test_that(".make_splits builds disjoint train/test index sets for CV", {
  sp <- .make_splits(cv(4, seed = 1), veteran, status_col = "status")
  expect_length(sp, 4L)
  for (s in sp) {
    expect_true(length(intersect(s$analysis, s$assessment)) == 0L)
    expect_setequal(sort(c(s$analysis, s$assessment)), seq_len(nrow(veteran)))
  }
})

test_that(".make_splits handles repeats, holdout and bootstrap", {
  rep_sp <- .make_splits(cv(3, repeats = 2, seed = 1), veteran, status_col = "status")
  expect_length(rep_sp, 6L)

  ho <- .make_splits(holdout(0.7, seed = 1), veteran, status_col = "status")
  expect_length(ho, 1L)
  expect_gt(length(ho[[1]]$analysis), length(ho[[1]]$assessment))

  bs <- .make_splits(bootstrap(10, seed = 1), veteran, status_col = "status")
  expect_length(bs, 10L)
  # OOB assessment is a strict subset, analysis has replacement duplicates
  expect_true(all(bs[[1]]$assessment %in% seq_len(nrow(veteran))))
})

test_that(".make_splits stratification falls back gracefully on unknown column", {
  sp <- .make_splits(cv(3, strata = "nonexistent", seed = 1), veteran, status_col = "status")
  expect_length(sp, 3L)
})
