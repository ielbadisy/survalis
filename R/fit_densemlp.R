#' Fit a densemlp Survival Model
#'
#' Fits a compact feedforward neural network for right-censored survival data
#' with the \pkg{densemlp} package (native C++ backend, no \pkg{torch}
#' dependency) and returns an object compatible with the `mlsurv_model`
#' interface. The IPCW integrated-Brier loss (`loss = "brier"`) is used by
#' default so that a full survival-probability curve is available; predictions
#' are interpolated from the network's discrete-time grid onto the requested
#' evaluation times.
#'
#' @param formula A survival formula \code{Surv(time, status) ~ predictors}.
#' @param data A \code{data.frame} containing the formula variables.
#' @param loss Survival loss: \code{"brier"} (default, discrete-time hazard grid,
#'   required for survival-curve prediction) or \code{"cox"} (a single linear
#'   risk score only).
#' @param n_bins Number of discrete-time bins for \code{loss = "brier"}
#'   (default \code{10}).
#' @param hidden_units Integer vector of hidden-layer widths (default
#'   \code{c(32, 16)}).
#' @param epochs Maximum training epochs (default \code{100}).
#' @param lr Adam learning rate (default \code{1e-3}).
#' @param ... Additional arguments passed to \code{densemlp::densemlp()} (for
#'   example \code{batch_size}, \code{dropout}, \code{batch_norm},
#'   \code{early_stopping}, \code{patience}, \code{ensemble}, \code{seed},
#'   \code{ncores}).
#'
#' @return An object of class \code{"mlsurv_model"} with elements \code{model},
#'   \code{learner} (\code{"densemlp"}), \code{formula}, \code{data},
#'   \code{time}, \code{status}, and \code{loss}.
#'
#' @seealso \code{\link{predict_densemlp}}, \code{\link[densemlp]{densemlp}},
#'   \code{\link{fit}} for the verb interface
#'
#' @examplesIf requireNamespace("densemlp", quietly = TRUE)
#' mod <- fit_densemlp(Surv(time, status) ~ age + karno + celltype, veteran,
#'                     n_bins = 8, epochs = 20, hidden_units = c(16, 8))
#' summary(mod)
#'
#' @keywords internal
#' @export
fit_densemlp <- function(formula, data, loss = c("brier", "cox"),
                         n_bins = 10, hidden_units = c(32, 16),
                         epochs = 100, lr = 1e-3, ...) {
  stopifnot(requireNamespace("densemlp", quietly = TRUE))
  loss <- match.arg(loss)

  model <- densemlp::densemlp(
    formula = formula, data = data, task = "survival", loss = loss,
    n_bins = n_bins, hidden_units = hidden_units, epochs = epochs, lr = lr,
    verbose = FALSE, ...
  )

  time_status <- all.vars(formula[[2]])
  structure(list(
    model = model,
    learner = "densemlp",
    formula = formula,
    data = data,
    time = time_status[1],
    status = time_status[2],
    loss = loss
  ), class = "mlsurv_model", engine = "densemlp")
}

#' Predict Survival Probabilities from a densemlp Model
#'
#' Generates predicted survival probabilities at the requested \code{times} for a
#' model fitted with \code{\link{fit_densemlp}} (\code{loss = "brier"}). The
#' network returns survival on its own discrete-time grid; values are linearly
#' interpolated (and constant-extrapolated) onto \code{times}.
#'
#' @param object An \code{"mlsurv_model"} from \code{\link{fit_densemlp}}.
#' @param newdata A \code{data.frame} of predictors.
#' @param times Numeric vector of evaluation time points.
#'
#' @return A \code{data.frame} of survival probabilities with one row per
#'   observation in \code{newdata} and one column per time point (columns named
#'   \code{"t=<time>"}).
#'
#' @seealso \code{\link{fit_densemlp}}
#'
#' @examplesIf requireNamespace("densemlp", quietly = TRUE)
#' mod <- fit_densemlp(Surv(time, status) ~ age + karno, veteran,
#'                     n_bins = 8, epochs = 20)
#' predict_densemlp(mod, veteran[1:5, ], times = c(100, 200, 300))
#'
#' @keywords internal
#' @export
predict_densemlp <- function(object, newdata, times) {
  if (!is.null(object$learner) && object$learner != "densemlp") {
    warning("Object passed to predict_densemlp() may not come from fit_densemlp().")
  }
  stopifnot(requireNamespace("densemlp", quietly = TRUE))
  if (!is.null(object$loss) && object$loss != "brier") {
    stop("predict_densemlp() needs loss = 'brier'; refit with fit_densemlp(..., loss = 'brier').")
  }

  pred <- stats::predict(object$model, newdata = newdata, type = "survival")
  pred <- as.matrix(pred)
  if (is.null(dim(pred))) pred <- matrix(pred, nrow = 1)

  grid <- suppressWarnings(as.numeric(sub("^t=", "", colnames(pred))))
  if (anyNA(grid)) grid <- as.numeric(object$model$survival_breaks)
  if (length(grid) != ncol(pred) || anyNA(grid)) {
    stop("Could not recover the densemlp discrete-time grid for interpolation.")
  }

  out <- matrix(NA_real_, nrow = nrow(pred), ncol = length(times))
  for (i in seq_len(nrow(pred))) {
    out[i, ] <- stats::approx(x = grid, y = pred[i, ], xout = times,
                              method = "linear", rule = 2, ties = "ordered")$y
  }
  .finalize_survmat(out, times = times)
}
