#' Predict from a fitted survival learner
#'
#' Produces predictions from a [fit()] object. Every matrix-valued `type`
#' returns a [survmat] (an `n x length(times)` numeric matrix carrying its grid
#' in the `"times"` attribute); the scalar-valued types (`"rmst"`, `"quantile"`,
#' `"median"`) return a numeric vector (or an `n x length(probs)` matrix for a
#' vector `probs`).
#'
#' @param object A `survalis_fit` from [fit()].
#' @param newdata A data frame of predictors. Defaults to the training data.
#'   Rows with missing values in model variables are rejected.
#' @param times Numeric evaluation grid. When `NULL`, [default_times()] is used
#'   and recorded on the result.
#' @param type One of:
#'   \describe{
#'     \item{`"survival"`}{predicted survival probability \eqn{S(t)} (default).}
#'     \item{`"risk"`}{cumulative risk \eqn{1 - S(t)}.}
#'     \item{`"chf"`}{cumulative hazard \eqn{-\log S(t)}.}
#'     \item{`"hazard"`}{piecewise-constant hazard on the grid.}
#'     \item{`"rmst"`}{restricted mean survival time; upper limit `tau` (`...`),
#'       default `max(times)`.}
#'     \item{`"quantile"`}{survival-time quantile(s); requires `probs` (`...`).}
#'     \item{`"median"`}{median survival time (`probs = 0.5`).}
#'   }
#' @param ... Type-specific arguments: `tau` for `"rmst"`, `probs` for
#'   `"quantile"`.
#'
#' @return A [survmat] for `"survival"`, `"risk"`, `"chf"`, `"hazard"`; a numeric
#'   vector for `"rmst"`, `"median"` and scalar-`probs` `"quantile"`; an
#'   `n x length(probs)` matrix for vector `probs`.
#'
#' @examplesIf requireNamespace("ranger", quietly = TRUE)
#' m <- fit(Surv(time, status) ~ age + karno, veteran, "ranger",
#'          spec = list(num.trees = 100))
#' predict(m, veteran[1:5, ], times = c(90, 180, 365))
#' predict(m, veteran[1:5, ], times = c(90, 180, 365), type = "risk")
#' predict(m, veteran[1:5, ], type = "rmst", tau = 365)
#' @method predict survalis_fit
#' @export
predict.survalis_fit <- function(object, newdata = NULL, times = NULL,
                                 type = c("survival", "risk", "chf", "hazard",
                                          "rmst", "quantile", "median"),
                                 ...) {
  type <- match.arg(type)
  dots <- list(...)

  status <- .fit_status_vector(object)
  times <- .resolve_times(times, object$data[[object$time]], status)

  if (is.null(newdata)) newdata <- object$data

  mv <- tryCatch(all.vars(stats::delete.response(stats::terms(object$formula))),
                 error = function(e) character(0))
  mv <- intersect(mv, names(newdata))
  if (length(mv)) {
    na_cols <- mv[colSums(is.na(newdata[, mv, drop = FALSE])) > 0L]
    if (length(na_cols)) {
      stop("`newdata` has rows with missing values in model variable(s): ",
           paste(na_cols, collapse = ", "), ".", call. = FALSE)
    }
  }

  pred_fun <- get(paste0("predict_", object$learner), envir = asNamespace("survalis"))
  raw <- pred_fun(object, newdata = newdata, times = times)
  S <- as_survmat(as.matrix(raw), times = times, quantity = "survival")
  Sdf <- as.data.frame(S)

  switch(type,
    survival = S,
    risk  = as_survmat(1 - as.matrix(S), times = times, quantity = "risk"),
    chf   = as_survmat(survmat_to_chf(Sdf), times = times, quantity = "chf"),
    hazard = as_survmat(survmat_to_haz(Sdf, times = times), times = times, quantity = "hazard"),
    rmst = {
      tau <- .null_default(dots$tau, max(times))
      out <- survmat_to_rmst(Sdf, times = times, tau = tau)
      stats::setNames(as.numeric(out), rownames(newdata))
    },
    quantile = {
      if (is.null(dots$probs)) {
        stop("type = \"quantile\" requires `probs` (a probability or vector).", call. = FALSE)
      }
      probs <- as.numeric(dots$probs)
      out <- vapply(probs, function(p) survmat_to_quantile(Sdf, times = times, p = p),
                    numeric(nrow(Sdf)))
      if (length(probs) == 1L) return(as.numeric(out))
      colnames(out) <- paste0("p=", probs)
      out
    },
    median = as.numeric(survmat_to_quantile(Sdf, times = times, p = 0.5))
  )
}
