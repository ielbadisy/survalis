
#' Fit a Stacked Survival Meta‑Learner (Time‑Varying NNLS)
#'
#' Learns nonnegative, time‑specific stacking weights over a set of base survival
#' learners by solving a nonnegative least squares (NNLS) problem at each
#' requested time \eqn{t^{*}}. For each \eqn{t^{*}}, the target is the indicator
#' \eqn{Y = I(T > t^{*})}, restricted to subjects whose status at \eqn{t^{*}} is
#' known and weighted by inverse probability of censoring weights; the features are the base learners' predicted survival
#' probabilities \eqn{S_\ell(t^{*} \mid x)}. Weights are constrained to be
#' nonnegative and are normalized to sum to 1 per time point.
#'
#' @param base_preds A named list of matrices/data frames, one per base learner,
#'   each of dimension \code{n x length(times)}, with columns named \code{"t=<time>"}.
#'   When \code{NULL} (the default), the base learners in \code{learners} are fitted
#'   from \code{formula} and \code{data}, and \code{base_preds} are their
#'   out-of-fold predictions.
#' @param time Numeric vector of observed event/censoring times (length \code{n}).
#'   Taken from \code{data} when \code{base_preds} is \code{NULL}.
#' @param status Numeric/binary vector of event indicators (1=event, 0=censor) (length \code{n}).
#'   Taken from \code{data} when \code{base_preds} is \code{NULL}.
#' @param times Numeric vector of evaluation times at which to learn weights. When
#'   \code{NULL} and \code{base_preds} is \code{NULL}, a grid of 20 times from
#'   \code{default_times()} is used.
#' @param base_models A named list of fitted base learner objects; names must
#'   match \code{names(base_preds)} and the learner names used by \code{predict_*()}.
#'   Fitted automatically when \code{base_preds} is \code{NULL}.
#' @param formula A \code{Surv()} formula that was used for the base models
#'   (stored for metadata and downstream scoring).
#' @param data The training data frame used for the base models (stored for metadata).
#' @param learners Base learners fitted when \code{base_preds} is \code{NULL}
#'   (default \code{c("coxph", "glmnet", "rsf")}). A learner that fails in
#'   cross-fitting is dropped with a warning.
#' @param folds Number of cross-fitting folds used to build the out-of-fold
#'   predictions (default \code{5}).
#' @param seed Seed for the cross-fitting folds.
#'
#' @details
#' For each \code{t} in \code{times}, this function fits
#' \deqn{\min_{w \ge 0} \| Y - S w \|_2^2}
#' where \eqn{Y = I(T > t)} and \eqn{S} is the matrix of base survival
#' probabilities at \eqn{t}. Only subjects with \eqn{T > t} or an observed event
#' by \eqn{t} enter, each weighted by \eqn{1 / G(\min(T, t)-)} with \eqn{G} the
#' reverse Kaplan-Meier estimate of the censoring distribution; the weights enter
#' the NNLS fit by \eqn{\sqrt{w}} row scaling. The NNLS solution from \pkg{nnls} is renormalized
#' to sum to 1 afterwards (if the solution is all zeros, equal weights are used).
#'
#' @return An object of class \code{c("mlsurv_model","survmetalearner")} with elements:
#' \describe{
#'   \item{weights}{Matrix \code{L x T} of nonnegative stacking weights, rows=learners, cols=\code{"t=<time>"}.}
#'   \item{base_models}{Named list of base learner fits.}
#'   \item{base_preds}{Named list of base prediction matrices on training data.}
#'   \item{learners}{Character vector of learner names (from \code{names(base_preds)}).}
#'   \item{formula, data, time, status}{Training metadata for scoring/reporting.}
#'   \item{learner}{The string \code{"survmetalearner"} (for \code{predict_*} dispatch).}
#' }
#'
#' @seealso \code{predict_survmetalearner}
#'
#' @examples
#' form <- Surv(time, status) ~ age + karno + trt
#' times <- c(80, 160)
#' mod_cox <- fit_coxph(form, data = veteran)
#' mod_rpart <- fit_rpart(form, data = veteran)
#' base_models <- list(coxph = mod_cox, rpart = mod_rpart)
#' base_preds <- list(
#'   coxph = predict_coxph(mod_cox, newdata = veteran, times = times),
#'   rpart = predict_rpart(mod_rpart, newdata = veteran, times = times)
#' )
#' meta_model <- fit_survmetalearner(
#'   base_preds = base_preds,
#'   time = veteran$time,
#'   status = veteran$status,
#'   times = times,
#'   base_models = base_models,
#'   formula = form,
#'   data = veteran
#' )
#' meta_model$weights
#' @keywords internal
#' @noRd

fit_survmetalearner <- function(base_preds = NULL, time = NULL, status = NULL, times = NULL,
                                base_models = NULL, formula = NULL, data = NULL,
                                learners = c("coxph", "glmnet", "rsf"),
                                folds = 5, seed = 123L) {
  stopifnot(is.data.frame(data))

  if (is.null(base_preds)) {
    # Formula-and-data interface (as used by evaluate() and compare()): fit the base
    # learners, and build the out-of-fold predictions that the stacking weights are
    # estimated from, so that the weights do not reward over-fitted base learners.
    stopifnot(!is.null(formula))
    parsed <- .parse_surv_formula(formula, data)
    time <- data[[parsed$time_col]]
    status <- .recode_status_vec(parsed, data)
    if (is.null(times)) times <- default_times(time, status, n = 20L)
    fit_one <- function(l, d) get(paste0("fit_", l), envir = asNamespace("survalis"))(formula = formula, data = d)
    pred_one <- function(l, m, d) {
      as.matrix(.finalize_survmat(get(paste0("predict_", l), envir = asNamespace("survalis"))(m, newdata = d, times = times), times = times))
    }
    n <- nrow(data)
    set.seed(seed)
    fold_id <- integer(n)
    for (s_val in unique(status)) {
      idx <- which(status == s_val)
      fold_id[idx] <- sample(rep_len(seq_len(folds), length(idx)))
    }
    base_preds <- list()
    for (l in learners) {
      oof <- matrix(NA_real_, nrow = n, ncol = length(times))
      err <- NULL
      for (k in seq_len(folds)) {
        tr <- data[fold_id != k, , drop = FALSE]
        te <- data[fold_id == k, , drop = FALSE]
        p <- tryCatch(pred_one(l, fit_one(l, tr), te), error = function(e) {
          err <<- conditionMessage(e)
          NULL
        })
        if (is.null(p)) break
        oof[fold_id == k, ] <- p
      }
      if (is.null(err)) base_preds[[l]] <- oof
      else warning("survmetalearner: base learner '", l, "' failed in cross-fitting ",
                   "and was dropped: ", err, call. = FALSE)
    }
    if (!length(base_preds)) stop("survmetalearner: every base learner failed.", call. = FALSE)
    base_models <- lapply(setNames(names(base_preds), names(base_preds)), function(l) fit_one(l, data))
  }
  stopifnot(is.list(base_preds), is.list(base_models))

  base_preds <- lapply(base_preds, function(x) {
    if (!is.matrix(x)) as.matrix(x) else x
  })

  learners <- names(base_preds)
  T <- length(times)
  weight_matrix <- matrix(NA_real_, nrow = length(learners), ncol = T,
                          dimnames = list(learners, paste0("t=", times)))

  # Reverse Kaplan-Meier of censoring, left limit G(u-) = P(C >= u).
  cens_km <- survival::survfit(survival::Surv(time, 1 - status) ~ 1)
  G_left <- function(u) c(1, cens_km$surv)[findInterval(u, cens_km$time, left.open = TRUE) + 1L]

  for (j in seq_along(times)) {
    t_star <- times[j]
    # Keep subjects whose status at t* is known (still at risk, or event by t*), and
    # weight them by IPCW 1 / G(min(T, t*)-).
    G <- G_left(pmin(time, t_star))
    keep <- (time > t_star | status == 1) & G > 0
    X <- sapply(base_preds, function(S) S[, j])
    X <- matrix(X, ncol = length(base_preds))
    keep <- keep & stats::complete.cases(X)
    Y <- as.numeric(time[keep] > t_star)
    sw <- sqrt(1 / G[keep])
    fit <- nnls::nnls(X[keep, , drop = FALSE] * sw, Y * sw)
    w <- coef(fit)
    # all-zero NNLS solution: fall back to equal weights rather than 0/0
    w <- if (sum(w) > 0) w / sum(w) else rep(1 / length(w), length(w))
    weight_matrix[, j] <- w
  }

  structure(list(
    weights = weight_matrix,
    base_models = base_models,
    base_preds = base_preds,
    times = times,
    formula = formula,
    data = data,
    time = time,
    status = status,
    learner = "survmetalearner",         # for predict_* lookup
    learners = learners
  ), class = c("mlsurv_model", "survmetalearner"), engine = "survalis")  # dual class

}

#' Predict with a Stacked Survival Meta‑Learner
#'
#' Produces stacked survival probabilities by combining base learner predictions
#' via the time‑specific weights learned by \code{fit_survmetalearner}.
#'
#' @param model A \code{"survmetalearner"} object returned by \code{fit_survmetalearner()}.
#' @param newdata A data frame of new observations for prediction.
#' @param times Numeric vector of evaluation times. When they differ from the times the
#'   meta-learner was trained at, the weights are linearly interpolated in time.
#'
#' @details
#' For each base learner listed in \code{model$learners}, the corresponding
#' \code{predict_<learner>} function is called to obtain \eqn{S_\ell(t \mid x)}.
#' Stacked predictions are computed as \eqn{\hat{S}(t \mid x) = \sum_\ell w_\ell(t) S_\ell(t \mid x)},
#' where \eqn{w_\ell(t)} are the learned nonnegative weights for time \eqn{t}.
#'
#' @return A data frame with one row per observation and one column per requested
#'   time (columns named \code{"t=<time>"}), containing stacked survival probabilities.
#'
#' @seealso \code{fit_survmetalearner}
#'
#' @examples
#' form <- Surv(time, status) ~ age + karno + trt
#' times <- c(80, 160)
#' mod_cox <- fit_coxph(form, data = veteran)
#' mod_rpart <- fit_rpart(form, data = veteran)
#' base_models <- list(coxph = mod_cox, rpart = mod_rpart)
#' base_preds <- list(
#'   coxph = predict_coxph(mod_cox, newdata = veteran, times = times),
#'   rpart = predict_rpart(mod_rpart, newdata = veteran, times = times)
#' )
#' meta_model <- fit_survmetalearner(
#'   base_preds = base_preds,
#'   time = veteran$time,
#'   status = veteran$status,
#'   times = times,
#'   base_models = base_models,
#'   formula = form,
#'   data = veteran
#' )
#' predict_survmetalearner(meta_model, newdata = veteran[1:3, ], times = times)
#' @keywords internal
#' @noRd

predict_survmetalearner <- function(model, newdata, times) {
  stopifnot(inherits(model, "survmetalearner"))

  learners <- model$learners
  base_models <- model$base_models
  requested_tnames <- paste0("t=", times)
  if (all(requested_tnames %in% colnames(model$weights))) {
    W <- model$weights[, requested_tnames, drop = FALSE]
  } else {
    # The evaluation times differ from the times the weights were fitted at: interpolate.
    fit_times <- as.numeric(sub("^t=", "", colnames(model$weights)))
    W <- .rows_by_times(apply(model$weights, 1L, function(w)
      stats::approx(x = fit_times, y = w, xout = times, rule = 2)$y))
    dimnames(W) <- list(rownames(model$weights), requested_tnames)
  }

  base_preds_test <- lapply(learners, function(learner) {
    pred_fun <- get(paste0("predict_", learner), envir = asNamespace("survalis"), mode = "function")
    pred <- pred_fun(base_models[[learner]], newdata = newdata, times = times)
    pred <- as.matrix(pred)
    colnames(pred) <- requested_tnames
    pred
  })
  names(base_preds_test) <- learners

  n <- nrow(base_preds_test[[1]])
  T <- length(times)
  survmat <- matrix(NA_real_, nrow = n, ncol = T)
  for (j in seq_len(T)) {
    survmat[, j] <- Reduce(`+`, lapply(seq_along(learners), function(k) base_preds_test[[k]][, j] * W[k, j]))
  }

  .finalize_survmat(survmat, times = times)
}
