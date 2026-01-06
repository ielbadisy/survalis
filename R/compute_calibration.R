#' Calibration of Survival Predictions at a Single Time Point
#'
#' Computes a nonparametric calibration curve for a survival model at one
#' evaluation time by binning predicted survival probabilities and comparing
#' bin-wise means to Kaplan–Meier-based observed survival, with bootstrap CIs.
#'
#' @param model A fitted model object created by \code{fit_*()} in \pkg{survalis}.
#'   Must contain a \code{learner} field used to dispatch to \code{predict_<learner>()}.
#' @param data A data.frame containing predictors and the survival outcome.
#' @param time Either a numeric vector of observed times (same length as \code{nrow(data)})
#'   or a single string giving the time column name in \code{data}.
#' @param status Either a numeric/integer vector of event indicators (1=event, 0=censor)
#'   or a single string giving the status column name in \code{data}.
#' @param eval_time A single numeric value at which to evaluate calibration.
#' @param n_bins Integer. Number of bins used to group predictions (default 10).
#' @param n_boot Integer. Number of bootstrap resamples for confidence intervals (default 100).
#' @param seed Integer or \code{NULL}. Random seed for bootstrap resampling (default 123).
#'   If \code{NULL}, no seed is set.
#' @param learner_name Optional string label for the learner. If \code{NULL},
#'   uses \code{model$learner}.
#'
#' @return A list with components:
#' \describe{
#'   \item{calibration_table}{Data frame with columns \code{bin}, \code{n_bin},
#'     \code{mean_pred_surv}, \code{observed_surv}, \code{lower_ci}, \code{upper_ci}.}
#'   \item{eval_time}{The evaluation time used.}
#'   \item{n_bins}{Number of bins.}
#'   \item{n_boot}{Number of bootstrap resamples.}
#'   \item{learner}{Learner label.}
#'   \item{ece}{Expected Calibration Error at \code{eval_time}.}
#' }
#' @export
compute_calibration <- function(model, data, time, status,
                                eval_time, n_bins = 10L, n_boot = 100L,
                                seed = 123, learner_name = NULL) {

  if (length(eval_time) != 1L || !is.finite(eval_time)) {
    stop("`eval_time` must be a single finite numeric value.")
  }
  n_bins <- as.integer(n_bins)
  n_boot <- as.integer(n_boot)
  if (n_bins < 1L) stop("`n_bins` must be >= 1.")
  if (n_boot < 0L) stop("`n_boot` must be >= 0.")

  if (!is.null(seed)) set.seed(seed)

  # Extract time/status if passed as column names
  if (is.character(time) && length(time) == 1L)   time   <- data[[time]]
  if (is.character(status) && length(status) == 1L) status <- data[[status]]

  if (length(time) != nrow(data) || length(status) != nrow(data)) {
    stop("`time` and `status` must have length equal to nrow(data).")
  }

  # Infer prediction function
  if (is.null(model$learner)) {
    stop("`model$learner` is NULL. Cannot dispatch to predict_<learner>().")
  }
  pred_name <- paste0("predict_", model$learner)
  if (!exists(pred_name, mode = "function")) {
    stop("Could not find prediction function: ", pred_name)
  }
  predict_function <- get(pred_name, mode = "function")

  # Predict survival at eval_time
  pred_surv_raw <- predict_function(model, newdata = data, times = eval_time)

  pred_surv <- if (is.data.frame(pred_surv_raw) && paste0("t=", eval_time) %in% names(pred_surv_raw)) {
    pred_surv_raw[[paste0("t=", eval_time)]]
  } else if (is.matrix(pred_surv_raw)) {
    as.vector(pred_surv_raw[, 1L])
  } else if (is.numeric(pred_surv_raw)) {
    as.vector(pred_surv_raw)
  } else {
    stop("Unexpected prediction format from ", pred_name,
         ". Expected numeric vector, matrix, or data.frame containing column 't=<eval_time>'.")
  }

  if (length(pred_surv) != nrow(data)) {
    stop("Predictions returned by ", pred_name, " do not match nrow(data).")
  }

  # Bin predictions by quantiles
  qs <- stats::quantile(pred_surv, probs = seq(0, 1, length.out = n_bins + 1L), na.rm = TRUE)
  qs <- unique(qs)

  if (length(qs) <= 1L) {
    bins <- rep.int(1L, length(pred_surv))
    n_bins_eff <- 1L
  } else {
    bins <- cut(pred_surv, breaks = qs, include.lowest = TRUE, labels = FALSE)
    n_bins_eff <- max(bins, na.rm = TRUE)
  }

  df <- data.frame(pred_surv = pred_surv, time = time, status = status, bin = bins)

  km_at_t <- function(tt, ss, t0) {
    sf <- survival::survfit(survival::Surv(tt, ss) ~ 1)
    s  <- summary(sf, times = t0, extend = TRUE)$surv
    if (length(s) == 0L) NA_real_ else as.numeric(s[1L])
  }

  # Per-bin summaries
  mean_pred <- tapply(df$pred_surv, df$bin, mean, na.rm = TRUE)
  n_bin     <- tapply(df$pred_surv, df$bin, length)

  obs_surv <- vapply(seq_len(n_bins_eff), function(b) {
    idx <- which(df$bin == b)
    if (length(idx) == 0L) return(NA_real_)
    km_at_t(df$time[idx], df$status[idx], eval_time)
  }, numeric(1L))

  calibration_table <- data.frame(
    bin = seq_len(n_bins_eff),
    n_bin = as.integer(n_bin[as.character(seq_len(n_bins_eff))]),
    mean_pred_surv = as.numeric(mean_pred[as.character(seq_len(n_bins_eff))]),
    observed_surv = as.numeric(obs_surv),
    lower_ci = NA_real_,
    upper_ci = NA_real_
  )

  # Bootstrap CIs (aligned to bins 1..n_bins_eff)
  if (n_boot > 0L) {
    boot_mat <- matrix(NA_real_, nrow = n_bins_eff, ncol = n_boot)

    for (bb in seq_len(n_boot)) {
      idx <- sample.int(nrow(df), replace = TRUE)
      dfb <- df[idx, , drop = FALSE]

      # compute observed survival per bin; bins may be missing in bootstrap sample
      for (b in seq_len(n_bins_eff)) {
        ii <- which(dfb$bin == b)
        boot_mat[b, bb] <- if (length(ii) == 0L) NA_real_ else km_at_t(dfb$time[ii], dfb$status[ii], eval_time)
      }
    }

    cis <- t(apply(boot_mat, 1L, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE))
    calibration_table$lower_ci <- cis[, 1L]
    calibration_table$upper_ci <- cis[, 2L]
  }

  # ECE (your existing metric function)
  surv_obj <- survival::Surv(time = df$time, event = df$status)
  sp_mat <- matrix(pred_surv, ncol = 1L, dimnames = list(NULL, paste0("t=", eval_time)))

  ece_value <- ece_survmat(
    object    = surv_obj,
    sp_matrix = sp_mat,
    t_star    = eval_time,
    n_bins    = n_bins_eff
  )

  list(
    calibration_table = calibration_table,
    eval_time = eval_time,
    n_bins = n_bins_eff,
    n_boot = n_boot,
    learner = if (is.null(learner_name)) model$learner else learner_name,
    ece = ece_value
  )
}


#' Plot Calibration Curve for Survival Predictions
#'
#' Produces a calibration plot comparing mean predicted survival (x-axis) to
#' observed survival with bootstrap CIs (y-axis) at a single evaluation time.
#'
#' @param calib_output Output list returned by [compute_calibration()].
#' @param smooth Logical; if `TRUE`, overlays a LOESS smooth of
#'   `observed_surv ~ mean_pred_surv`.
#'
#' @return A `ggplot2` object showing bin-wise calibration points, bootstrap
#' error bars, the 45° reference line, and (optionally) a smooth curve.
#'
#' @details
#' Points above the diagonal indicate underprediction (observed survival higher
#' than predicted), while points below indicate overprediction.
#'
#' @examples
#' # calib <- compute_calibration(model, data, time = "time", status = "status", eval_time = 365)
#' # p <- plot_calibration(calib)
#' # p
#'
#' @seealso [compute_calibration()]
#' @export
plot_calibration <- function(calib_output, smooth = TRUE) {
  if (!is.list(calib_output) || is.null(calib_output$calibration_table)) {
    stop("`calib_output` must be the list returned by compute_calibration().")
  }
  ct <- calib_output$calibration_table
  req <- c("mean_pred_surv", "observed_surv", "lower_ci", "upper_ci")
  miss <- setdiff(req, names(ct))
  if (length(miss) > 0L) {
    stop("`calibration_table` is missing required columns: ", paste(miss, collapse = ", "))
  }

  etime <- calib_output$eval_time
  nboot <- calib_output$n_boot
  nbins <- calib_output$n_bins
  ece   <- calib_output$ece

  subtitle_text <- paste0(
    "n_bins = ", nbins,
    ", bootstrap = ", nboot,
    if (!is.null(ece) && is.finite(ece)) paste0(", ECE = ", signif(ece, 3)) else ""
  )

  learner_label <- calib_output$learner
  if (is.null(learner_label) || !nzchar(learner_label)) learner_label <- "model"

  p <- ggplot2::ggplot(ct, ggplot2::aes(x = mean_pred_surv, y = observed_surv)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = lower_ci, ymax = upper_ci),
      width = 0.02
    ) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      x = "Mean Predicted Survival",
      y = "Observed Survival",
      title = paste0("Calibration at t = ", etime, " | ", learner_label),
      subtitle = subtitle_text
    ) +
    ggplot2::theme_minimal()

  if (isTRUE(smooth)) {
    p <- p + ggplot2::geom_smooth(method = "loess", se = FALSE, formula = y ~ x)
  }

  p
}
