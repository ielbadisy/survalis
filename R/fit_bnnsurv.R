#' Fit a Bagged Nearest-Neighbour Survival Model
#'
#' Fits a bagged ensemble of k-nearest-neighbour Kaplan-Meier learners
#' (Lowsky et al., 2013, bagged as in Wright's \pkg{bnnSurvival}). The engine
#' is implemented natively in \pkg{survalis}.
#'
#' @section Method:
#' Each base learner \eqn{b = 1, \dots, B} draws a bootstrap (or subsample) of
#' the training rows and, optionally, a random subset of the predictors. For a
#' new observation \eqn{x}, the learner finds the \eqn{k} training rows closest
#' to \eqn{x} and returns their weighted Kaplan-Meier curve
#' \deqn{\hat S_b(t \mid x) = \prod_{t_j \le t}
#'   \left(1 - \frac{\sum_i w_i \delta_i 1\{T_i = t_j\}}
#'                  {\sum_i w_i 1\{T_i \ge t_j\}}\right),}
#' with weights \eqn{w_i = w(d(x, x_i))}. The ensemble prediction is the
#' average \eqn{\hat S(t \mid x) = B^{-1} \sum_b \hat S_b(t \mid x)}, a
#' non-increasing step function in \eqn{t}. Training rows tied with the
#' \eqn{k}-th nearest distance share the remaining neighbour slots equally, so
#' predictions do not depend on the row order of the data.
#'
#' The Mahalanobis distance uses the covariance \eqn{\Sigma_b} of the base
#' learner's own sample: \eqn{d(x, y)^2 = (x - y)^\top \Sigma_b^{+} (x - y)},
#' where \eqn{\Sigma_b^{+}} is the Moore-Penrose pseudo-inverse. Constant and
#' collinear columns (for example a dummy that never varies in a bootstrap
#' sample) therefore contribute nothing instead of breaking the inversion.
#' \code{metric = "euclidean"} uses the Euclidean distance on standardized
#' columns.
#'
#' Differences from \pkg{bnnSurvival}: factors are one-hot encoded (not coded
#' as integer labels) and new data reuse the training levels; feature subsets
#' are drawn from all predictors (a factor counts as one predictor); constant
#' columns are detected from the variances; prediction runs in the current
#' process (no forking); and predictions at arbitrary times are read off the
#' step function rather than linearly interpolated. With numeric predictors,
#' all predictors per learner, constant weights and no tied distances, the
#' two engines give the same curves for the same seed.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
#' @param data A data frame containing all variables referenced in \code{formula}.
#' @param k Integer, number of neighbours for each base learner. Default \code{5}.
#' @param num_base_learners Integer, number of base learners in the ensemble.
#'   Default \code{10}.
#' @param num_features_per_base_learner Integer or \code{NULL} (all
#'   predictors). Number of predictors drawn at random for each base learner.
#' @param metric Distance, \code{"mahalanobis"} (default) or
#'   \code{"euclidean"} (on standardized columns).
#' @param weighting_function Function mapping a vector of non-negative
#'   distances to non-negative neighbour weights. Default: constant weights.
#' @param replace Logical; sample rows with replacement for each base learner.
#'   Default \code{TRUE}.
#' @param sample_fraction Optional numeric in \code{(0, 1]}; fraction of rows
#'   drawn for each base learner. Defaults to \code{1} when
#'   \code{replace = TRUE} and \code{0.6321} otherwise.
#'
#' @return An object of class \code{"mlsurv_model"} with elements:
#' \describe{
#'   \item{model}{The fitted ensemble: design information, time grid and the
#'     base learners (row indices, columns and whitening maps).}
#'   \item{learner}{Scalar \code{"bnnsurv"}.}
#'   \item{engine}{Scalar \code{"survalis"}.}
#'   \item{formula, data}{Inputs preserved for downstream use.}
#'   \item{time, status}{Character names of the survival outcome fields.}
#' }
#'
#' @references
#' Lowsky DJ, Ding Y, Lee DKK, McCulloch CE, Ross LF, Thistlethwaite JR,
#' Zenios SA (2013). A K-nearest neighbors survival probability prediction
#' method. \emph{Statistics in Medicine}, 32(12), 2062-2069.
#' \doi{10.1002/sim.5673}
#'
#' @seealso `predict_bnnsurv()`, `tune_bnnsurv()`
#'
#' @examples
#' mod <- fit_bnnsurv(
#'   Surv(time, status) ~ age + karno + diagtime + prior,
#'   data = veteran
#' )
#' head(predict_bnnsurv(mod, newdata = veteran[1:5, ], times = c(50, 100)))
#'
#' @keywords internal
#' @noRd

fit_bnnsurv <- function(formula, data,
                        k = 5,
                        num_base_learners = 10,
                        num_features_per_base_learner = NULL,
                        metric = c("mahalanobis", "euclidean"),
                        weighting_function = function(x) x * 0 + 1,
                        replace = TRUE,
                        sample_fraction = NULL) {
  metric <- match.arg(metric)
  stopifnot(is.function(weighting_function))
  k <- as.integer(k)
  num_base_learners <- as.integer(num_base_learners)
  if (length(k) != 1L || is.na(k) || k < 1L) stop("`k` must be a positive integer.")
  if (length(num_base_learners) != 1L || is.na(num_base_learners) || num_base_learners < 1L) {
    stop("`num_base_learners` must be a positive integer.")
  }
  if (is.null(sample_fraction)) sample_fraction <- if (replace) 1 else 0.6321
  if (!is.numeric(sample_fraction) || length(sample_fraction) != 1L ||
      !(sample_fraction > 0 && sample_fraction <= 1)) {
    stop("`sample_fraction` must be a number in (0, 1].")
  }

  data <- as.data.frame(data)
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  y <- stats::model.response(mf)
  if (!inherits(y, "Surv") || attr(y, "type") != "right") {
    stop("`formula` must have a right-censored Surv(time, status) response.")
  }
  tt <- stats::terms(mf)
  X <- .bnn_design(tt, mf, xlev = NULL)
  if (!ncol(X)) stop("`formula` must contain at least one predictor.")
  time <- as.numeric(y[, 1L])
  status <- as.numeric(y[, 2L])
  n <- nrow(X)

  # predictors are sampled as whole terms: a factor is one predictor
  assign <- attr(X, "assign")
  p <- max(assign)
  m <- if (is.null(num_features_per_base_learner)) p else as.integer(num_features_per_base_learner)
  if (length(m) != 1L || is.na(m) || m < 1L || m > p) {
    stop("`num_features_per_base_learner` must be between 1 and ", p, ".")
  }

  n_draw <- if (!replace && sample_fraction == 1) n else floor(n * sample_fraction)
  if (n_draw < 2L) stop("Each base learner needs at least two rows; increase `sample_fraction`.")

  learners <- lapply(seq_len(num_base_learners), function(b) {
    rows <- if (!replace && sample_fraction == 1) seq_len(n) else
      sample.int(n, n_draw, replace = replace)
    vars <- if (m == p) seq_len(p) else sort(sample.int(p, m))
    cols <- which(assign %in% vars)
    map <- .bnn_whitener(X[rows, cols, drop = FALSE], metric)
    list(rows = rows, cols = cols, center = map$center, W = map$W,
         Z = .bnn_project(X[rows, cols, drop = FALSE], map))
  })

  model <- list(
    terms = stats::delete.response(tt),
    xlevels = stats::.getXlevels(tt, mf),
    columns = colnames(X),
    time = time,
    status = status,
    grid = sort(unique(time)),
    learners = learners,
    k = k,
    metric = metric,
    weighting_function = weighting_function,
    replace = replace,
    sample_fraction = sample_fraction
  )

  structure(list(
    model   = model,
    learner = "bnnsurv",
    engine  = "survalis",
    formula = formula,
    data    = data,
    time    = all.vars(formula)[[2]],
    status  = all.vars(formula)[[3]]
  ), class = "mlsurv_model")
}


#' Predict Survival with a Bagged Nearest-Neighbour Model
#'
#' Generates survival probabilities from an \code{mlsurv_model} fitted by
#' \code{fit_bnnsurv}.
#'
#' @param object A fitted \code{mlsurv_model} returned by \code{fit_bnnsurv()}.
#' @param newdata A data frame of new observations for prediction.
#' @param times Optional numeric vector of evaluation time points. If \code{NULL}
#'   (default), survival is returned at every distinct training time.
#'
#' @return A base \code{data.frame} with one row per observation and columns
#'   named \code{"t=<time>"} containing survival probabilities in \code{[0, 1]},
#'   non-increasing over time.
#'
#' @details
#' The ensemble curve is a step function that jumps only at training event
#' times, so the value at time \eqn{t} is the curve at the last training time
#' not after \eqn{t} (and \code{1} before the first one). Factor levels in
#' \code{newdata} are matched to the training levels; unseen levels and
#' missing predictor values are errors.
#'
#' @seealso `fit_bnnsurv()`, `tune_bnnsurv()`
#'
#' @examples
#' mod <- fit_bnnsurv(Surv(time, status) ~ age + karno + diagtime + prior, data = veteran)
#' pred <- predict_bnnsurv(mod, newdata = veteran[1:3, ], times = c(50, 100, 200))
#' pred
#'
#' @keywords internal
#' @noRd

predict_bnnsurv <- function(object, newdata, times = NULL) {
  stopifnot(object$learner == "bnnsurv")
  fit <- object$model
  newdata <- as.data.frame(newdata)

  mf <- stats::model.frame(fit$terms, data = newdata, xlev = fit$xlevels,
                           na.action = stats::na.pass)
  Xnew <- .bnn_design(fit$terms, mf, xlev = fit$xlevels)
  if (!identical(colnames(Xnew), fit$columns)) {
    stop("Predictors in `newdata` do not match the training design.")
  }
  if (anyNA(Xnew)) stop("`newdata` has missing predictor values.")

  grid <- fit$grid
  S <- matrix(0, length(grid), nrow(Xnew))
  for (bl in fit$learners) {
    S <- S + .bnn_learner_surv(bl, Xnew, fit)
  }
  S <- t(S / length(fit$learners))

  if (is.null(times)) {
    return(.finalize_survmat(S, times = grid))
  }

  times <- as.numeric(times)
  stopifnot(length(times) > 0L, all(is.finite(times)))
  # read the step function: S(t) = S(t_j) for the last grid time t_j <= t
  j <- findInterval(times, grid)
  survmat <- cbind(1, S)[, j + 1L, drop = FALSE]
  .finalize_survmat(survmat, times = times)
}


# Design matrix without intercept; treatment-coded dummies for factors.
.bnn_design <- function(tt, mf, xlev) {
  X <- stats::model.matrix(tt, mf, xlev = xlev)
  keep <- colnames(X) != "(Intercept)"
  assign <- attr(X, "assign")[keep]
  X <- X[, keep, drop = FALSE]
  attr(X, "assign") <- assign
  X
}

# Linear map W such that squared Euclidean distances of (x - center) %*% W
# are the requested squared distances. For Mahalanobis, W = V diag(1/sqrt(l))
# over the positive eigenvalues of the sample covariance, which gives the
# pseudo-inverse metric and silently drops constant or collinear directions.
.bnn_whitener <- function(X, metric) {
  center <- colMeans(X)
  Xc <- sweep(X, 2L, center)
  p <- ncol(X)
  if (metric == "euclidean") {
    s <- sqrt(colSums(Xc^2) / max(nrow(X) - 1L, 1L))
    inv <- ifelse(s > 0, 1 / s, 0)
    return(list(center = center, W = diag(inv, nrow = p, ncol = p)))
  }
  S <- crossprod(Xc) / max(nrow(X) - 1L, 1L)
  e <- eigen(S, symmetric = TRUE)
  tol <- max(dim(S)) * max(e$values, 0) * .Machine$double.eps
  pos <- e$values > tol
  W <- if (any(pos)) {
    e$vectors[, pos, drop = FALSE] %*% diag(1 / sqrt(e$values[pos]), sum(pos))
  } else {
    matrix(0, p, 1L)
  }
  list(center = center, W = W)
}

.bnn_project <- function(X, map) {
  sweep(X, 2L, map$center) %*% map$W
}

# Survival curves of one base learner for every row of Xnew, on fit$grid
# (times x rows).
.bnn_learner_surv <- function(bl, Xnew, fit) {
  Z <- bl$Z
  Znew <- .bnn_project(Xnew[, bl$cols, drop = FALSE], bl)
  # squared distances, newdata rows x learner rows
  D2 <- -2 * tcrossprod(Znew, Z) + rowSums(Znew^2)
  D2 <- D2 + rep(rowSums(Z^2), each = nrow(D2))
  D2[D2 < 0] <- 0
  k <- min(fit$k, ncol(D2))
  nn <- .bnn_neighbours(D2, k)

  w <- fit$weighting_function(sqrt(nn$d2))
  if (length(w) != length(nn$d2) || any(!is.finite(w)) || any(w < 0)) {
    stop("`weighting_function` must return one finite, non-negative weight per distance.")
  }
  # a neighbourhood whose kernel weights all vanish falls back to equal weights
  zero <- rowsum(w * nn$share, nn$row)[, 1L] == 0
  w[zero[nn$row]] <- 1
  w <- w * nn$share

  idx <- bl$rows[nn$col]
  .weighted_km(fit$time[idx], fit$status[idx], w,
               row = nn$row, nrow = nrow(D2), grid = fit$grid)
}

# The k nearest columns of each row of the squared-distance matrix D2.
# Columns tied with the k-th distance share the remaining neighbour slots
# equally (share < 1), so the neighbourhood does not depend on row order.
.bnn_neighbours <- function(D2, k) {
  kth <- if (ncol(D2) == k) apply(D2, 1L, max) else
    apply(D2, 1L, function(d) sort.int(d, partial = k)[k])
  eps <- 1e-10 * pmax(kth, 1e-12)
  lo <- kth - eps
  hi <- kth + eps
  closer <- rowSums(D2 < lo)
  tied <- rowSums(D2 <= hi) - closer
  sel <- which(D2 <= hi, arr.ind = TRUE)
  row <- sel[, 1L]
  d2 <- D2[sel]
  share <- ifelse(d2 < lo[row], 1, (k - closer[row]) / tied[row])
  list(row = row, col = sel[, 2L], d2 = d2, share = share)
}

# Weighted Kaplan-Meier for many small samples at once. Observation i belongs
# to curve row[i]; curves are evaluated on `grid`, which contains every time.
# Only the (curve, time) cells holding observations are touched: the at-risk
# weight of a cell is a reverse cumulative sum within its curve, and
# log S = cumulative sum of log(1 - deaths / at risk) over the event cells.
# After the last observation a curve keeps its last value. The result is
# times x curves.
.weighted_km <- function(time, status, w, row, nrow, grid) {
  G <- length(grid)
  cell <- (row - 1) * G + match(time, grid)
  ucell <- sort(unique(cell))
  g <- match(cell, ucell)
  leave <- rowsum(w, g)[, 1L]
  death <- rowsum(w * (status == 1), g)[, 1L]
  r <- (ucell - 1) %/% G + 1
  j <- (ucell - 1) %% G + 1

  tail_w <- rev(cumsum(rev(leave)))            # weight at or after each cell
  last <- c(which(diff(r) != 0), length(r))   # last cell of each curve
  after <- c(tail_w[-1L], 0)[last]            # weight in later curves
  at_risk <- tail_w - after[match(r, r[last])]

  ev <- death > 0
  L <- matrix(0, G, nrow)
  L[cbind(j[ev], r[ev])] <- log1p(-pmin(death[ev] / at_risk[ev], 1))
  exp(apply(L, 2L, cumsum))
}


#' Tune Bagged Nearest-Neighbour Survival Hyperparameters (Cross-Validation)
#'
#' Cross-validates \code{fit_bnnsurv} over a user-supplied grid and
#' aggregates metrics (for example, \code{"cindex"}, \code{"ibs"}). Optionally
#' refits and returns the best configuration.
#'
#' @param formula A survival formula of the form \code{Surv(time, status) ~ predictors}.
#' @param data Training data frame.
#' @param times Numeric vector of evaluation time points used during CV.
#' @param param_grid A data frame (for example from \code{expand.grid()}) with
#'   columns \code{k}, \code{num_base_learners}, and \code{sample_fraction}.
#'   Defaults cover a small illustrative search.
#' @param metrics Character vector of metric names to compute and summarize.
#'   The first entry is treated as the primary ranking metric.
#' @param folds Integer number of cross-validation folds. Default \code{5}.
#' @param seed Integer random seed for reproducibility.
#' @param ncores Integer number of CPU cores passed to \code{cv_survlearner}
#'   for fold evaluation (default \code{1}).
#' @param refit_best Logical; if \code{TRUE}, refits the best configuration on
#'   the full data and returns that model. Otherwise returns the results table.
#'
#' @return If \code{refit_best = FALSE}, a \code{tibble} with one row per
#'   configuration containing \code{metrics} and the tuning parameters, with a
#'   \code{failed} flag for combinations that errored during CV. If
#'   \code{refit_best = TRUE}, an \code{mlsurv_model} refit at the selected
#'   hyperparameters.
#'
#' @details
#' Evaluation is performed by \code{cv_survlearner()} using
#' \code{fit_bnnsurv()} and \code{predict_bnnsurv()} to match production code
#' paths. Results are ordered by the first entry in \code{metrics}.
#'
#' @seealso `fit_bnnsurv()`, `predict_bnnsurv()`, `cv_survlearner()`
#'
#' @examples
#' \donttest{
#' grid <- expand.grid(
#'   k = c(2, 3),
#'   num_base_learners = c(5, 10),
#'   sample_fraction = c(0.5, 1),
#'   stringsAsFactors = FALSE
#' )
#'
#' res <- tune_bnnsurv(
#'   formula = Surv(time, status) ~ age + karno + diagtime + prior,
#'   data = veteran,
#'   times = c(100, 200),
#'   param_grid = grid,
#'   refit_best = FALSE
#' )
#' res
#'
#' best <- tune_bnnsurv(
#'   formula = Surv(time, status) ~ age + karno + diagtime + prior,
#'   data = veteran,
#'   times = c(100, 200),
#'   param_grid = grid,
#'   refit_best = TRUE
#' )
#' head(predict_bnnsurv(best, newdata = veteran[1:5, ], times = c(50, 100)))
#' }
#' 
#' @keywords internal
#' @noRd

tune_bnnsurv <- function(formula, data, times,
                         param_grid = expand.grid(
                           k = c(2, 3),
                           num_base_learners = c(30, 50),
                           sample_fraction = c(0.5, 1),
                           stringsAsFactors = FALSE
                         ),
                         metrics = c("cindex", "ibs"),
                         folds = 5,
                         seed = 123,
                         ncores = 1,
                         refit_best = FALSE) {
  results <- .pmap_rbind_dt(param_grid, function(k, num_base_learners, sample_fraction) {
    params <- list(k = k, num_base_learners = num_base_learners, sample_fraction = sample_fraction)

    cv_result <- tryCatch({
      cv_survlearner(
        formula = formula,
        data = data,
        fit_fun = fit_bnnsurv,
        pred_fun = predict_bnnsurv,
        times = times,
        metrics = metrics,
        folds = folds,
        seed = seed,
        ncores = ncores,
        k = k,
        num_base_learners = num_base_learners,
        sample_fraction = sample_fraction
      )
    }, error = function(e) NULL)

    if (is.null(cv_result)) {
      return(do.call(data.frame, c(params, list(failed = TRUE), stringsAsFactors = FALSE)))
    }

    summary <- cv_summary(cv_result)
    .wide_metric_row(c(params, list(failed = FALSE)), summary)
  })

  results <- results[!results$failed, , drop = FALSE]

  if (nrow(results) == 0) {
    warning("All tuning combinations failed.")
    return(data.frame())
  }

  maximize <- .metric_maximize(metrics[1])
  results <- .arrange_by_metric_dt(results, metrics[1], maximize)

  if (refit_best) {
    best <- results[1, ]
    model <- fit_bnnsurv(
      formula = formula,
      data = data,
      k = best$k,
      num_base_learners = best$num_base_learners,
      sample_fraction = best$sample_fraction
    )
    return(model)
  }

  return(results)
}
