
compute_ale <- function(model, predict_function, newdata, feature, times, grid.size = 20) {
  if (!feature %in% names(newdata)) stop("Feature not found in data.")

  if (is.factor(newdata[[feature]]) || is.character(newdata[[feature]])) {
    stop("ALE is only defined for numeric continuous features. For categorical features, use PDP instead.")
  }

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


#' Plot ALE result
#'
#' @param ale_result A list returned by `compute_ale()`
#' @param feature Name of the feature used (for labeling)
#' @param which Type of ALE plot: "per_time" or "integrated"
#' @param smooth Whether to add a smoothed line (only for numeric)
#' @export
plot_ale <- function(ale_result, feature, which = c("per_time", "integrated"), smooth = FALSE) {
  which <- match.arg(which)
  library(ggplot2)

  if (which == "per_time") {
    df <- ale_result$ale
    df_long <- reshape(
      df,
      varying = names(df)[-1],
      v.names = "ale",
      timevar = "time",
      times = names(df)[-1],
      idvar = "feature_value",
      direction = "long"
    )
    df_long$time <- gsub("t=", "", df_long$time)

    p <- ggplot(df_long, aes(x = feature_value, y = ale, color = time, group = time)) +
      geom_line() +
      labs(title = paste("ALE curves for", feature), x = feature, y = "ALE effect") +
      theme_minimal()
    if (smooth) p <- p + geom_smooth(se = FALSE, method = "loess")
    return(p)
  }

  if (which == "integrated" && !is.null(ale_result$integrated)) {
    df <- ale_result$integrated
    p <- ggplot(df, aes(x = feature_value, y = integrated_ale)) +
      geom_line() +
      geom_point() +
      labs(title = paste("Integrated ALE curve for", feature), x = feature, y = "Integrated ALE") +
      theme_minimal()
    if (smooth) p <- p + geom_smooth(se = FALSE, method = "loess")
    return(p)
  }

  stop("Invalid 'which' argument or missing integrated data.")
}


library(survival)
data(veteran, package = "survival")

# fit a survival model
mod <- fit_cox(Surv(time, status) ~ age + karno + celltype, data = veteran)

# Compute ALE for continuous feature 'age' at multiple time points
ale_result <- compute_ale(
  model = mod,
  predict_function = predict_cox,
  newdata = veteran,
  feature = "age",
  times = c(100, 200, 300)
)

# plot per-time ALE
plot_ale(ale_result, feature = "age", which = "per_time")

# plot integrated ALE with smoothing
plot_ale(ale_result, feature = "age", which = "integrated", smooth = TRUE)

# attempting ALE on a categorical feature (should error)
try({compute_ale(mod, predict_ranger, veteran, feature = "celltype", times = 100)})

