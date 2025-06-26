compute_rmstdiff_iptw <- function(model, data, A, W, times, predict_fun) {
  weights <- compute_iptw_weights(A, W)
  survmat <- predict_fun(model, newdata = data, times = times)

  # Compute individual RMSTs
  rmst_indiv <- compute_rmst(survmat, times)

  # Weighted RMST averages
  treated <- A == 1
  control <- A == 0

  rmst_treated <- sum(rmst_indiv[treated] * weights[treated]) / sum(weights[treated])
  rmst_control <- sum(rmst_indiv[control] * weights[control]) / sum(weights[control])

  list(
    method = "iptw",
    estimand = "rmst_difference",
    tau = max(times),
    rmst_treated = rmst_treated,
    rmst_control = rmst_control,
    ate_rmst = rmst_treated - rmst_control
  )
}


plot_rmstdiff_iptw <- function(res) {
  df <- data.frame(
    Group = c("Treated", "Control"),
    RMST = c(res$rmst_treated, res$rmst_control)
  )

  ggplot(df, aes(x = Group, y = RMST, fill = Group)) +
    geom_col(width = 0.5, alpha = 0.8) +
    labs(
      title = sprintf("IPTW RMST Comparison (τ = %.1f)", res$tau),
      y = "Restricted Mean Survival Time", x = ""
    ) +
    theme_minimal() +
    scale_fill_manual(values = c("blue", "red")) +
    theme(legend.position = "none")
}

# Fit model
model <- fit_xgboost(Surv(time, status) ~ A + X1 + X2, data = df)
times <- seq(0, 5, 0.1)

# Compute IPTW RMST difference
res_rmst <- compute_rmstdiff_iptw(
  model = model,
  data = df,
  A = df$A,
  W = df[, c("X1", "X2")],
  times = times,
  predict_fun = predict_xgboost
)

# Plot
plot_rmstdiff_iptw(res_rmst)
