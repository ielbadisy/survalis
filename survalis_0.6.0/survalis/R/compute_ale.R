#' Accumulated Local Effects (ALE) for Survival Models
#'
#' Computes ALE curves for a numeric (continuous) feature with respect to
#' survival probabilities at one or more evaluation times. ALE summarizes the
#' average *local* effect of changing a feature within small intervals, is
#' robust to correlated features, and is centered to have mean zero.
#'
#' @param model An `mlsurv_model` created by a `fit_*()` function. Must include
#'   a valid `learner` so that `predict_<learner>()` can be dispatched.
#' @param newdata Data frame used to compute ALE (typically the training set or
#'   a representative sample).
#' @param feature Single **numeric/continuous** feature name for which to
#'   compute ALE. Categorical features are not supported here (use PDP/ICE).
#' @param times Numeric vector of time points at which to evaluate survival
#'   probabilities.
#' @param grid.size Integer number of quantile cut points used to build the ALE
#'   grid (default 20). The algorithm uses quantiles of `newdata[[feature]]`.
#'
#' @details
#' For consecutive quantile bins \eqn{[z_k, z_{k+1}]} of the target feature,
#' ALE integrates the *local* change in the model prediction when moving the
#' feature from \eqn{z_k} to \eqn{z_{k+1}} while holding all other features at
#' their observed values, and then accumulates these differences across bins.
#' For survival models, predictions are survival probabilities at `times`.
#' The returned ALE curves are centered (mean zero across the grid) per time.
#'
#' @return A list with:
#' \describe{
#'   \item{ale}{A data frame with columns `feature_value` and one column per
#'   time (`"t=<time>"`) containing centered ALE effects.}
#'   \item{integrated}{If multiple times were provided, a data frame with
#'   columns `feature_value` and `integrated_ale` equal to the mean of per-time
#'   ALE effects across `times`; otherwise `NULL`.}
#' }
#'
#' @examples
#' # mod <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = veteran)
#' # ale_res <- compute_ale(mod, newdata = veteran, feature = "karno",
#' #                        times = c(100, 200, 300))
#' # head(ale_res$ale)
#'
#' @seealso [plot_ale()], [compute_pdp()]
#' @export

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

#' Plot ALE Curves for Survival Models
#'
#' Visualizes ALE results produced by [compute_ale()] either as per-time curves
#' (one curve per evaluation time) or as an integrated curve averaged across
#' times.
#'
#' @param ale_result A list returned by [compute_ale()].
#' @param feature Character name of the feature (for axis labeling only).
#' @param which Either `"per_time"` to plot one ALE curve per time, or
#'   `"integrated"` to plot the time-averaged ALE (requires multiple `times`).
#' @param smooth Logical; if `TRUE`, overlays a LOESS smooth on the plotted
#'   curve(s).
#'
#' @details
#' Per-time plots show how the feature's local effect varies across different
#' evaluation times. The integrated plot summarizes the average effect over the
#' supplied time grid (simple mean across times of the centered ALE values).
#'
#' @return A `ggplot2` object.
#'
#' @examples
#' # p1 <- plot_ale(ale_res, feature = "karno", which = "per_time")
#' # p2 <- plot_ale(ale_res, feature = "karno", which = "integrated", smooth = TRUE)
#'
#' @seealso [compute_ale()], [compute_pdp()], [plot_pdp()]
#' @export

plot_ale <- function(ale_result, feature, which = c("per_time", "integrated"), smooth = FALSE) {
  which <- match.arg(which)

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

