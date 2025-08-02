
compute_calibration <- function(model, predict_function, data, time, status,
                                eval_time, n_bins = 10, n_boot = 100, seed = 123) {
  library(dplyr)
  library(survival)

  set.seed(seed)

  # if time/status are passed as names, extract from data
  if (is.character(time) && length(time) == 1) {
    time <- data[[time]]
  }
  if (is.character(status) && length(status) == 1) {
    status <- data[[status]]
  }

  # predict survival probabilities
  pred_surv_raw <- predict_function(model, newdata = data, times = eval_time)
  if (is.data.frame(pred_surv_raw) && paste0("t=", eval_time) %in% colnames(pred_surv_raw)) {
    pred_surv <- pred_surv_raw[[paste0("t=", eval_time)]]
  } else if (is.matrix(pred_surv_raw)) {
    pred_surv <- as.vector(pred_surv_raw[, 1])
  } else if (is.numeric(pred_surv_raw)) {
    pred_surv <- pred_surv_raw
  } else {
    stop("Unexpected prediction format. Ensure predict_function returns a data.frame with column named t=eval_time.")
  }

  # binning
  bins <- cut(pred_surv, breaks = quantile(pred_surv, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE),
              include.lowest = TRUE, labels = FALSE)
  df <- data.frame(pred_surv = pred_surv, time = time, status = status, bin = bins)

  # calibration table
  calibration_table <- df %>%
    group_by(bin) %>%
    summarise(
      mean_pred_surv = mean(pred_surv, na.rm = TRUE),
      observed_surv = {
        surv_fit <- survfit(Surv(time, status) ~ 1, data = cur_data())
        surv_summary <- summary(surv_fit, times = eval_time, extend = TRUE)
        if (length(surv_summary$surv) == 0) NA else surv_summary$surv
      },
      .groups = "drop"
    )

  # bootstrap CI
  boot_results <- replicate(n_boot, {
    idx <- sample(seq_len(nrow(df)), replace = TRUE)
    df_boot <- df[idx, ]
    df_boot %>%
      group_by(bin) %>%
      summarise(
        observed_surv = {
          surv_fit <- survfit(Surv(time, status) ~ 1, data = cur_data())
          surv_summary <- summary(surv_fit, times = eval_time, extend = TRUE)
          if (length(surv_summary$surv) == 0) NA else surv_summary$surv
        }, .groups = "drop"
      ) %>%
      pull(observed_surv)
  }, simplify = "matrix")

  boot_ci <- apply(boot_results, 1, function(x) {
    quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  })
  calibration_table$lower_ci <- boot_ci[1, ]
  calibration_table$upper_ci <- boot_ci[2, ]

  list(
    calibration_table = calibration_table,
    eval_time = eval_time,
    n_bins = n_bins,
    n_boot = n_boot
  )
}


plot_calibration <- function(calib_output, smooth = TRUE) {
  library(ggplot2)
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
      title = paste0("Calibration at t = ", etime)#,
      #subtitle = paste0("n_bins = ", nbins, ", Bootstrap = ", nboot, " resamples")
    ) +
    theme_minimal()

  if (smooth) {
    p <- p + geom_smooth(method = "loess", se = FALSE, color = "blue", formula = y ~ x)
  }

  return(p)
}


calib_result <- compute_calibration(
  model = mod_cox,
  predict_function = predict_coxph,
  data = veteran,
  time = "time",
  status = "status",
  eval_time = 80,
  n_bins = 10,
  n_boot = 30
)


plot_calibration(calib_result)

## may be interesting to think about time varying calibration plot!!!
