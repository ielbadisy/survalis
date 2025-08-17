#' Calibration of Survival Predictions at a Single Time Point
#'
#' Computes a nonparametric calibration curve for a survival model at one
#' evaluation time by binning predicted survival probabilities and comparing
#' bin-wise means to Kaplan-Meier-based observed survival, with bootstrap CIs.
#'
#' @param model An `mlsurv_model` fitted via a `fit_*()` function; must include a
#'   valid `learner` so the corresponding `predict_<learner>()` can be dispatched,
#'   and (ideally) a `$data` field for level alignment.
#' @param data A data frame with predictors and survival outcome columns.
#' @param time Survival time; either a numeric vector of the same length as
#'   `nrow(data)` or a single string giving the column name in `data`.
#' @param status Event indicator; either a numeric/logical vector or a single
#'   string giving the column name in `data`. Events should be coded 1, censoring 0.
#' @param eval_time Single numeric time at which to assess calibration.
#' @param n_bins Integer number of quantile-based bins used to group predictions.
#' @param n_boot Integer number of bootstrap resamples for confidence intervals.
#' @param seed Integer seed for reproducibility of binning/bootstrap.
#' @param learner_name Optional character override for labeling the learner in
#'   downstream plots (defaults to `model$learner`).
#'
#' @details
#' Predicted survival at `eval_time` is obtained from the appropriate
#' `predict_<learner>()`. Predictions are split into `n_bins` quantile bins.
#' For each bin, the function reports:
#' mean predicted survival, observed survival at `eval_time` from a Kaplan-Meier
#' fit on the bin's rows, and bootstrap percentile (2.5%, 97.5%) CIs on the
#' observed survival computed by resampling rows with replacement.
#'
#' @return A list with components:
#' \describe{
#'   \item{calibration_table}{A data frame with columns `bin`, `mean_pred_surv`,
#'         `observed_surv`, `lower_ci`, `upper_ci`.}
#'   \item{eval_time}{The evaluation time used.}
#'   \item{n_bins}{Number of bins.}
#'   \item{n_boot}{Number of bootstrap resamples.}
#'   \item{learner}{The learner label (from `learner_name` or `model$learner`).}
#' }
#'
#' @examples
#' # mod_cox <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = veteran)
#' # calib <- compute_calibration(
#' #   model = mod_cox, data = veteran,
#' #   time = "time", status = "status",
#' #   eval_time = 80, n_bins = 10, n_boot = 30
#' # )
#' # plot_calibration(calib)
#'
#' @seealso [plot_calibration()]
#' @export

compute_calibration <- function(model, data, time, status,
                                eval_time, n_bins = 10, n_boot = 100,
                                seed = 123, learner_name = NULL) {
  set.seed(seed)

  # check eval_time is a single value
  if (length(eval_time) != 1) {
    stop("`eval_time` must be a single numeric value. Calibration at multiple time points is not supported.")
  }

  # infer prediction function
  if (is.null(model$learner) || !exists(paste0("predict_", model$learner))) {
    stop("Could not infer prediction function. Ensure model was created using fit_*() and includes a valid 'learner'.")
  }
  predict_function <- get(paste0("predict_", model$learner))

  # if time/status are passed as strings, extract from data
  if (is.character(time) && length(time) == 1) time <- data[[time]]
  if (is.character(status) && length(status) == 1) status <- data[[status]]

  # predict survival at eval_time
  pred_surv_raw <- predict_function(model, newdata = data, times = eval_time)

  pred_surv <- if (is.data.frame(pred_surv_raw) && paste0("t=", eval_time) %in% colnames(pred_surv_raw)) {
    pred_surv_raw[[paste0("t=", eval_time)]]
  } else if (is.matrix(pred_surv_raw)) {
    as.vector(pred_surv_raw[, 1])
  } else if (is.numeric(pred_surv_raw)) {
    pred_surv_raw
  } else {
    stop("Unexpected prediction format. Ensure predict_<learner>() returns a data.frame or matrix with survival probabilities.")
  }

  # bin predictions
  bins <- cut(pred_surv,
              breaks = quantile(pred_surv, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE),
              include.lowest = TRUE, labels = FALSE)
  df <- data.frame(pred_surv = pred_surv, time = time, status = status, bin = bins)

  # calibration table
  calibration_table <- df |>
    group_by(bin) |>
    summarise(
      mean_pred_surv = mean(pred_surv, na.rm = TRUE),
      observed_surv = {
        surv_fit <- survfit(Surv(time, status) ~ 1, data = pick(everything()))
        surv_summary <- summary(surv_fit, times = eval_time, extend = TRUE)
        if (length(surv_summary$surv) == 0) NA else surv_summary$surv
      },
      .groups = "drop"
    )

  # bootstrap CIs
  boot_results <- replicate(n_boot, {
    idx <- sample(seq_len(nrow(df)), replace = TRUE)
    df_boot <- df[idx, ]
    df_boot |>
      group_by(bin) |>
      summarise(
        observed_surv = {
          surv_fit <- survfit(Surv(time, status) ~ 1, data = pick(everything()))
          surv_summary <- summary(surv_fit, times = eval_time, extend = TRUE)
          if (length(surv_summary$surv) == 0) NA else surv_summary$surv
        }, .groups = "drop"
      ) |>
      pull(observed_surv)
  }, simplify = "matrix")

  boot_ci <- apply(boot_results, 1, function(x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
  calibration_table$lower_ci <- boot_ci[1, ]
  calibration_table$upper_ci <- boot_ci[2, ]

  list(
    calibration_table = calibration_table,
    eval_time = eval_time,
    n_bins = n_bins,
    n_boot = n_boot,
    learner = learner_name %||% model$learner
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
#' # p <- plot_calibration(calib)
#' # p
#'
#' @seealso [compute_calibration()]
#' @export

plot_calibration <- function(calib_output, smooth = TRUE) {
  ct <- calib_output$calibration_table
  etime <- calib_output$eval_time
  nboot <- calib_output$n_boot
  nbins <- calib_output$n_bins

  p <- ggplot(ct, aes(x = mean_pred_surv, y = observed_surv)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.02) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      x = "Mean Predicted Survival",
      y = "Observed Survival",
      title = paste0("Calibration at t = ", etime, " | ", calib_result$learner)#,
      #subtitle = paste0("n_bins = ", nbins, ", Bootstrap = ", nboot, " resamples")
    ) +
    theme_minimal()

  if (smooth) {
    p <- p + geom_smooth(method = "loess", se = FALSE, color = "blue", formula = y ~ x)
  }

  return(p)
}


