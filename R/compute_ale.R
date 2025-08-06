compute_ale <- function(model, newdata, feature, times, grid.size = 20) {
  if (!feature %in% names(newdata)) stop("Feature not found in data.")

  if (is.factor(newdata[[feature]]) || is.character(newdata[[feature]])) {
    stop("ALE is only defined for numeric continuous features. For categorical features, use PDP instead.")
  }

  # Infer prediction function from model$learner
  if (is.null(model$learner) || !exists(paste0("predict_", model$learner))) {
    stop("Invalid or missing 'learner' in model. Make sure model was created using fit_*() functions.")
  }
  predict_function <- get(paste0("predict_", model$learner))

  # Create quantile grid
  grid <- quantile(newdata[[feature]], probs = seq(0, 1, length.out = grid.size), na.rm = TRUE)
  grid <- unique(grid)
  K <- length(grid) - 1
  n_times <- length(times)
  ale_mat <- matrix(0, nrow = K + 1, ncol = n_times)

  for (i in seq_len(K)) {
    data_lower <- newdata
    data_upper <- newdata
    data_lower[[feature]] <- grid[i]
    data_upper[[feature]] <- grid[i + 1]

    pred_lower <- predict_function(model, data_lower, times = times)
    pred_upper <- predict_function(model, data_upper, times = times)

    delta <- colMeans(pred_upper - pred_lower, na.rm = TRUE)
    ale_mat[i + 1, ] <- ale_mat[i, ] + delta
  }

  ale_mat_centered <- scale(ale_mat, scale = FALSE)
  ale_df <- data.frame(feature_value = grid, ale_mat_centered)
  colnames(ale_df)[-1] <- paste0("t=", times)

  integrated <- if (n_times > 1) {
    int_effect <- rowMeans(ale_mat_centered, na.rm = TRUE)
    data.frame(feature_value = grid, integrated_ale = int_effect)
  } else {
    NULL
  }

  return(list(ale = ale_df, integrated = integrated))
}



library(survival)
data(veteran)

# Fit a model
mod_coxph <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = veteran)

# Compute ALE
ale_result <- compute_ale(
  model = mod_coxph,
  newdata = veteran,
  feature = "karno",
  times = c(100, 200, 300)
)

# Plot per-time
plot_ale(ale_result, feature = "karno", which = "per_time")

# Plot integrated ALE
plot_ale(ale_result, feature = "karno", which = "integrated", smooth = TRUE)


