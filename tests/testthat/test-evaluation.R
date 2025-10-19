
test_that("cindex_survmat works on simple no-censoring data and handles ties/columns", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")

  time   <- c(1, 2, 3, 4)
  status <- c(1, 1, 1, 1)
  y      <- survival::Surv(time, status)

  sp <- cbind("t=5" = c(0.1, 0.3, 0.6, 0.8))
  c1 <- cindex_survmat(y, predicted = sp, t_star = 5)
  expect_named(c1, "C index")
  expect_equal(unname(c1), 1, tolerance = 1e-12)

  c2 <- cindex_survmat(y, predicted = sp, t_star = NULL)
  expect_equal(unname(c2), unname(c1), tolerance = 1e-12)

  sp_tie <- cbind("t=5" = c(0.1, 0.1, 0.3, 0.4))
  c_tie <- cindex_survmat(y, predicted = sp_tie, t_star = 5)
  expect_gt(unname(c_tie), 0.5)
  expect_lt(unname(c_tie), 1)

  expect_error(
    cindex_survmat(y, predicted = sp, t_star = 10),
    "not found in predicted survival matrix"
  )

  expect_error(
    cindex_survmat(list(), predicted = sp, t_star = 5),
    "survival object"
  )
})

test_that("brier reduces to unweighted score when no censoring; checks length", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")

  time   <- c(1, 2, 4, 6)
  status <- c(1, 1, 1, 1)
  y      <- survival::Surv(time, status)
  t_star <- 3

  pre_sp <- c(0.2, 0.5, 0.7, 0.9)
  exp_val <- mean(c(0.2^2, 0.5^2, (1 - 0.7)^2, (1 - 0.9)^2))

  b <- brier(y, pre_sp = pre_sp, t_star = t_star)
  expect_named(b, "brier")
  expect_equal(unname(b), round(exp_val, 6), tolerance = 1e-12)

  expect_error(
    brier(y, pre_sp = pre_sp[-1], t_star = t_star),
    "Length of predictions"
  )
})

test_that("ibs_survmat integrates Brier over multiple times", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")

  set.seed(1)
  n <- 8
  time   <- c(1,2,3,4,5,6,7,8)
  status <- rep(1, n)
  y      <- survival::Surv(time, status)

  times <- c(2, 4, 6)
  sp <- cbind("t=2" = rep(0.8, n),
              "t=4" = rep(0.6, n),
              "t=6" = rep(0.3, n))

  b_at <- function(t, Scol) {
    early <- which(time < t)
    late  <- which(time >= t)
    mean(c(rep(Scol[1]^2, length(early)), rep((1 - Scol[1])^2, length(late))))
  }
  b2 <- b_at(2, sp[,1]); b4 <- b_at(4, sp[,2]); b6 <- b_at(6, sp[,3])

  dt <- diff(times)
  exp_ibs <- (b2*dt[1] + b4*dt[2]) / (max(times) - min(times))

  ibs <- ibs_survmat(y, sp_matrix = sp, times = times)
  expect_named(ibs, "ibs")
  expect_equal(unname(ibs), round(exp_ibs, 6), tolerance = 1e-12)

  ibs1 <- ibs_survmat(y, sp_matrix = sp[, 1, drop = FALSE], times = times[1])
  expect_equal(unname(ibs1), unname(brier(y, pre_sp = sp[,1], t_star = times[1])))
})

test_that("iae_survmat and ise_survmat are ~0 when mean prediction matches KM", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("survival")

  time   <- c(1,2,3,5,8,9,10,12)
  status <- c(1,1,0,1,0,1,0,1)
  y      <- survival::Surv(time, status)

  km <- survival::survfit(y ~ 1)
  grid <- km$time[km$n.event > 0]
  s_km <- km$surv[km$n.event > 0]

  sp <- matrix(rep(s_km, each = length(time)),
               nrow = length(time), byrow = FALSE)

  iae <- iae_survmat(y, sp_matrix = sp, times = grid)
  ise <- ise_survmat(y, sp_matrix = sp, times = grid)

  expect_named(iae, "iae"); expect_named(ise, "ise")
  expect_equal(unname(iae), 0, tolerance = 1e-8)
  expect_equal(unname(ise), 0, tolerance = 1e-8)
})
