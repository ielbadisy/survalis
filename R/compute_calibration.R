#' Calibration of Survival Predictions at a Single Time Point
#'
#' Computes a nonparametric calibration curve for a survival model at one
#' evaluation time by binning predicted survival probabilities and comparing
#' bin-wise means to Kaplan-Meier-based observed survival, with bootstrap CIs.
#'
#' @return A list with components:
#' \describe{
#'   \item{calibration_table}{Data frame with `bin`, `n_bin`, `mean_pred_surv`,
#'         `observed_surv`, `lower_ci`, `upper_ci`.}
#'   \item{eval_time}{The evaluation time used.}
#'   \item{n_bins}{Number of bins.}
#'   \item{n_boot}{Number of bootstrap resamples.}
#'   \item{learner}{Learner label.}
#'   \item{ece}{Expected Calibration Error at \code{eval_time}.}
#' }
#' @export
compute_calibration <- function(model, data, time, status,
  eval_time, n_bins = 10, n_boot = 100,
  seed = 123, learner_name = NULL) {
set.seed(seed)

if (length(eval_time) != 1) {
stop("`eval_time` must be a single numeric value. Calibration at multiple time points is not supported.")
}

# infer prediction function
if (is.null(model$learner) || !exists(paste0("predict_", model$learner))) {
stop("Could not infer prediction function. Ensure model was created using fit_*() and includes a valid 'learner'.")
}
predict_function <- get(paste0("predict_", model$learner))

# if time/status are passed as strings, extract from data
if (is.character(time) && length(time) == 1)  time   <- data[[time]]
if (is.character(status) && length(status) == 1) status <- data[[status]]

# Predict survival at eval_time
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

# Bin predictions
qs <- stats::quantile(
pred_surv,
probs = seq(0, 1, length.out = n_bins + 1L),
na.rm = TRUE
)
qs <- unique(qs)

if (length(qs) <= 1L) {
bins <- rep.int(1L, length(pred_surv))
} else {
bins <- cut(pred_surv, breaks = qs, include.lowest = TRUE, labels = FALSE)
}

df <- data.frame(pred_surv = pred_surv, time = time, status = status, bin = bins)

# Calibration table
calibration_table <- df |>
dplyr::group_by(bin) |>
dplyr::summarise(
n_bin = dplyr::n(),
mean_pred_surv = mean(pred_surv, na.rm = TRUE),
observed_surv = {
surv_fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = dplyr::pick(dplyr::everything()))
surv_summary <- summary(surv_fit, times = eval_time, extend = TRUE)
if (length(surv_summary$surv) == 0) NA_real_ else surv_summary$surv
},
.groups = "drop"
)

# Bootstrap CIs
if (n_boot > 0L) {
boot_results <- replicate(n_boot, {
idx <- sample(seq_len(nrow(df)), replace = TRUE)
df_boot <- df[idx, ]
df_boot |>
dplyr::group_by(bin) |>
dplyr::summarise(
observed_surv = {
surv_fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = dplyr::pick(dplyr::everything()))
surv_summary <- summary(surv_fit, times = eval_time, extend = TRUE)
if (length(surv_summary$surv) == 0) NA_real_ else surv_summary$surv
},
.groups = "drop"
) |>
dplyr::pull(observed_surv)
}, simplify = "matrix")

boot_ci <- apply(boot_results, 1, function(x) stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
calibration_table$lower_ci <- boot_ci[1, ]
calibration_table$upper_ci <- boot_ci[2, ]
} else {
calibration_table$lower_ci <- NA_real_
calibration_table$upper_ci <- NA_real_
}

# ECE using the metric layer
surv_obj <- survival::Surv(time = df$time, event = df$status)
sp_mat   <- matrix(
pred_surv,
ncol = 1L,
dimnames = list(NULL, paste0("t=", eval_time))
)

ece_value <- ece_survmat(
object   = surv_obj,
sp_matrix = sp_mat,
t_star   = eval_time,
n_bins   = n_bins
)

list(
calibration_table = calibration_table,
eval_time = eval_time,
n_bins = n_bins,
n_boot = n_boot,
learner = learner_name %||% model$learner,
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
#' # p <- plot_calibration(calib)
#' # p
#'
#' @seealso [compute_calibration()]
#' @export
plot_calibration <- function(calib_output, smooth = TRUE) {
  ct    <- calib_output$calibration_table
  etime <- calib_output$eval_time
  nboot <- calib_output$n_boot
  nbins <- calib_output$n_bins
  ece   <- calib_output$ece

  subtitle_text <- paste0(
    "n_bins = ", nbins,
    ", bootstrap = ", nboot,
    if (!is.null(ece) && !is.na(ece)) paste0(", ECE = ", signif(ece, 3)) else ""
  )

  p <- ggplot2::ggplot(ct, ggplot2::aes(x = mean_pred_surv, y = observed_surv)) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lower_ci, ymax = upper_ci), width = 0.02) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggplot2::coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Mean Predicted Survival",
      y = "Observed Survival",
      title = paste0("Calibration at t = ", etime, " | ", calib_output$learner),
      subtitle = subtitle_text
    ) +
    ggplot2::theme_minimal()

  if (smooth) {
    p <- p + ggplot2::geom_smooth(method = "loess", se = FALSE, formula = y ~ x)
  }

  p
}
