#' Standardized RMST Estimate or Causal-Style RMST Contrast via G-Computation
#'
#' Estimates, from a fitted \pkg{survalis} model (\code{fit_*()} output),
#' either (i) the standardized marginal restricted mean survival time (RMST)
#' at a fixed horizon, or (ii) a covariate-adjusted contrast in RMST between
#' the two levels of a binary grouping variable via g-computation -- with
#' percentile bootstrap inference adapted from the \pkg{tvrmst} package's
#' dynamic-RMST bootstrap design. The two variants are dispatched by the
#' \code{trt_col} argument (\code{NULL} = marginal estimate, matching the
#' \code{group = NULL} dispatch already used by \code{\link{plot_survmat}}),
#' rather than by separate function names.
#'
#' @param model An \code{mlsurv_model} object from a \code{fit_*()} learner
#'   (must expose \code{model$learner} so \code{predict_<learner>()} can be
#'   dispatched, as in \code{\link{compute_calibration}}).
#' @param data A data.frame of covariates to standardize over.
#' @param tau Positive numeric restriction time for RMST.
#' @param times Numeric vector of evaluation times passed to
#'   \code{predict_<learner>()}; must be dense enough to integrate up to
#'   \code{tau} (see \code{\link{survmat_to_rmst}}).
#' @param trt_col Optional character; name of a binary grouping/exposure
#'   column in \code{data}. If \code{NULL} (default), \code{estimate_rmst()}
#'   returns the standardized marginal RMST at \code{tau} for \code{data} as
#'   given. If supplied, its two observed values are used as counterfactual
#'   levels (the larger value is the "treated"/A=1 arm) and the function
#'   returns the standardized contrast between them.
#' @param R Number of bootstrap replicates. Default \code{200}.
#' @param conf Confidence level in (0, 1). Default \code{0.95}.
#' @param seed Optional RNG seed for reproducibility.
#'
#' @details
#' \strong{Marginal mode (\code{trt_col = NULL}).} Individual RMSTs
#' \eqn{\mathrm{RMST}_i(\tau)} are obtained via \code{\link{survmat_to_rmst}}
#' from \eqn{\hat S(t \mid X_i)}, and the standardized estimate is
#' \eqn{\widehat{\mathrm{RMST}}(\tau) = n^{-1}\sum_i \mathrm{RMST}_i(\tau)}.
#'
#' \strong{Contrast mode (\code{trt_col} supplied).} Counterfactual survival
#' curves \eqn{\hat S(t \mid X_i, A{=}1)} and \eqn{\hat S(t \mid X_i, A{=}0)}
#' are obtained by setting \code{trt_col} to each level for the entire data
#' set and predicting from the fixed, already-fitted \code{model}
#' (g-computation/standardization). The standardized contrast is
#' \deqn{
#'   \widehat{\Delta}_{\mathrm{RMST}}(\tau) =
#'   \frac{1}{n}\sum_{i=1}^n \mathrm{RMST}_i(\tau, A{=}1) -
#'   \frac{1}{n}\sum_{i=1}^n \mathrm{RMST}_i(\tau, A{=}0).
#' }
#'
#' In both modes, uncertainty is quantified by a nonparametric percentile
#' bootstrap over subjects (rows of \code{data}), holding \code{model} fixed
#' (i.e. conditional on the fitted ensemble, as in \pkg{tvrmst}'s
#' \code{boot_rmst_delta()}): each replicate resamples rows with replacement
#' and recomputes the estimate.
#'
#' When \code{trt_col} corresponds to a manipulable intervention and standard
#' causal identification assumptions hold, \eqn{\widehat{\Delta}_{\mathrm{RMST}}(\tau)}
#' admits a causal interpretation; for non-manipulable attributes (e.g. sex,
#' biomarker status) it should be read as a covariate-adjusted marginal
#' contrast rather than a causal effect.
#'
#' @return A list of class \code{"rmst_estimate"} with elements
#'   \code{estimate}, \code{lo}, \code{hi}, \code{conf}, \code{tau},
#'   \code{trt_col} (\code{NULL} in marginal mode), \code{levels} (the two
#'   levels compared in contrast mode, or \code{NULL}), \code{mode}
#'   (\code{"marginal"} or \code{"contrast"}), and \code{boot} (the vector of
#'   bootstrap replicate estimates).
#'
#' @references
#' Royston P, Parmar MKB (2013). "Restricted Mean Survival Time: An
#' Alternative to the Hazard Ratio for the Design and Analysis of Randomized
#' Trials with a Time-to-Event Outcome." \emph{BMC Medical Research
#' Methodology}, 13, 152.
#'
#' @seealso [survmat_to_rmst()], [plot_rmst()]
#'
#' @examples
#' \donttest{
#' mod <- fit_coxph(Surv(time, status) ~ age + karno + trt, data = veteran)
#' times <- seq(10, 300, by = 10)
#'
#' marg <- estimate_rmst(mod, data = veteran, tau = 200, times = times, R = 20, seed = 1)
#' marg$estimate
#'
#' contrast <- estimate_rmst(mod, data = veteran, tau = 200, times = times,
#'                           trt_col = "trt", R = 20, seed = 1)
#' contrast$estimate
#' }
#' @keywords internal
#' @export
estimate_rmst <- function(model, data, tau, times, trt_col = NULL,
                          R = 200, conf = 0.95, seed = NULL) {
  .legacy_notice()
  if (is.null(model$learner) || !exists(paste0("predict_", model$learner))) {
    stop("Could not infer prediction function. Ensure model was created using fit_*() and includes a valid 'learner'.")
  }
  if (!is.finite(tau) || tau <= 0) stop("`tau` must be a positive finite number.")
  if (!is.finite(conf) || conf <= 0 || conf >= 1) stop("`conf` must be in (0,1).")
  if (!is.null(seed)) set.seed(seed)

  predict_function <- get(paste0("predict_", model$learner))
  n <- nrow(data)

  if (is.null(trt_col)) {
    mode <- "marginal"
    levels_trt <- NULL

    stat_fun <- function(idx) {
      d <- data[idx, , drop = FALSE]
      S <- predict_function(model, newdata = d, times = times)
      mean(survmat_to_rmst(S, times = times, tau = tau))
    }
  } else {
    mode <- "contrast"
    if (!trt_col %in% colnames(data)) stop("`trt_col` not found in `data`.")
    levels_trt <- sort(unique(data[[trt_col]]))
    if (length(levels_trt) != 2) stop("`trt_col` must have exactly two distinct values.")
    lo_level <- levels_trt[1]
    hi_level <- levels_trt[2]

    stat_fun <- function(idx) {
      d <- data[idx, , drop = FALSE]
      d1 <- d; d1[[trt_col]] <- hi_level
      d0 <- d; d0[[trt_col]] <- lo_level

      S1 <- predict_function(model, newdata = d1, times = times)
      S0 <- predict_function(model, newdata = d0, times = times)

      rmst1 <- survmat_to_rmst(S1, times = times, tau = tau)
      rmst0 <- survmat_to_rmst(S0, times = times, tau = tau)
      mean(rmst1) - mean(rmst0)
    }
  }

  estimate <- stat_fun(seq_len(n))
  boot <- vapply(seq_len(R), function(r) stat_fun(sample.int(n, n, replace = TRUE)), numeric(1))

  alpha <- (1 - conf) / 2
  lo <- unname(stats::quantile(boot, probs = alpha, na.rm = TRUE))
  hi <- unname(stats::quantile(boot, probs = 1 - alpha, na.rm = TRUE))

  structure(
    list(
      estimate = estimate,
      lo = lo,
      hi = hi,
      conf = conf,
      tau = tau,
      trt_col = trt_col,
      levels = levels_trt,
      mode = mode,
      boot = boot
    ),
    class = "rmst_estimate"
  )
}

#' RMST Estimate or Contrast Curve Over a Grid of Restriction Times
#'
#' Evaluates \code{\link{estimate_rmst}} across a grid of restriction times
#' \eqn{\tau} and plots the resulting curve with percentile bootstrap
#' confidence bands. As in \code{\link{estimate_rmst}}, marginal vs. contrast
#' mode is dispatched by \code{trt_col} rather than by separate function names.
#'
#' @param model An \code{mlsurv_model} object from a \code{fit_*()} learner.
#' @param data A data.frame of covariates.
#' @param tau_seq Numeric vector of restriction times to evaluate.
#' @param times Numeric vector of evaluation times passed to
#'   \code{predict_<learner>()}.
#' @param trt_col Optional character; name of a binary grouping/exposure
#'   column. \code{NULL} (default) plots the marginal RMST(tau) curve;
#'   supplying it plots the RMST contrast curve.
#' @param R Number of bootstrap replicates per \code{tau}. Default \code{100}.
#' @param conf Confidence level in (0, 1). Default \code{0.95}.
#' @param seed Optional RNG seed.
#' @param title Plot title. If missing, an automatically generated title is
#'   used. Pass \code{NULL} to omit the title entirely.
#'
#' @return A \pkg{ggplot2} object with an attribute \code{"curve"} holding the
#'   underlying \code{data.frame} (\code{tau}, \code{estimate}, \code{lo},
#'   \code{hi}).
#'
#' @seealso [estimate_rmst()]
#'
#' @examples
#' \donttest{
#' mod <- fit_coxph(Surv(time, status) ~ age + karno + trt, data = veteran)
#' times <- seq(10, 300, by = 10)
#' plot_rmst(mod, data = veteran, tau_seq = seq(50, 250, by = 50),
#'                 times = times, trt_col = "trt", R = 10, seed = 1)
#' }
#' @export
plot_rmst <- function(model, data, tau_seq, times, trt_col = NULL,
                            R = 100, conf = 0.95, seed = NULL, title) {
  curve <- do.call(rbind, lapply(tau_seq, function(tau) {
    est <- estimate_rmst(model, data = data, tau = tau, times = times,
                         trt_col = trt_col, R = R, conf = conf, seed = seed)
    data.frame(tau = tau, estimate = est$estimate, lo = est$lo, hi = est$hi)
  }))

  if (missing(title)) {
    title <- if (is.null(trt_col)) {
      "Standardized marginal RMST"
    } else {
      paste0("Adjusted marginal RMST contrast (", trt_col, ")")
    }
  }

  ylab <- if (is.null(trt_col)) "RMST" else expression(Delta ~ RMST)

  p <- ggplot2::ggplot(curve, ggplot2::aes(x = tau, y = estimate))
  if (!is.null(trt_col)) {
    p <- p + ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")
  }
  p <- p +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.2, fill = .survalis_palette[1]) +
    ggplot2::geom_line(color = .survalis_palette[1], linewidth = 1) +
    ggplot2::labs(x = "Restriction time (tau)", y = ylab) +
    theme_survalis()
  if (!is.null(title)) p <- p + ggplot2::labs(title = title)

  attr(p, "curve") <- curve
  p
}
