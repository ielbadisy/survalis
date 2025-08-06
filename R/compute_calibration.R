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
  bins <- cut(pred_surv,
              breaks = quantile(pred_surv, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE),
              include.lowest = TRUE, labels = FALSE)
  df <- data.frame(pred_surv = pred_surv, time = time, status = status, bin = bins)

  # Calibration table
  calibration_table <- df %>%
    group_by(bin) %>%
    summarise(
      mean_pred_surv = mean(pred_surv, na.rm = TRUE),
      observed_surv = {
        surv_fit <- survfit(Surv(time, status) ~ 1, data = pick(everything()))
        surv_summary <- summary(surv_fit, times = eval_time, extend = TRUE)
        if (length(surv_summary$surv) == 0) NA else surv_summary$surv
      },
      .groups = "drop"
    )

  # Bootstrap CIs
  boot_results <- replicate(n_boot, {
    idx <- sample(seq_len(nrow(df)), replace = TRUE)
    df_boot <- df[idx, ]
    df_boot %>%
      group_by(bin) %>%
      summarise(
        observed_surv = {
          surv_fit <- survfit(Surv(time, status) ~ 1, data = pick(everything()))
          surv_summary <- summary(surv_fit, times = eval_time, extend = TRUE)
          if (length(surv_summary$surv) == 0) NA else surv_summary$surv
        }, .groups = "drop"
      ) %>%
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

calib_result <- compute_calibration(
  model = mod_cox,
  data = veteran,
  time = "time",
  status = "status",
  eval_time = c(80),
  n_bins = 10,
  n_boot = 30
)



plot_calibration(calib_result)
