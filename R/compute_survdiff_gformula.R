compute_survdiff_gformula <- function(model, data, times, predict_fun) {
  data1 <- data; data1$A <- 1
  data0 <- data; data0$A <- 0

  # Predict survival under A = 1 and A = 0
  S1 <- predict_fun(model, newdata = data1, times = times)
  S0 <- predict_fun(model, newdata = data0, times = times)

  S1_avg <- colMeans(S1)
  S0_avg <- colMeans(S0)

  list(
    method = "gformula",
    estimand = "survival_difference",
    time = times,
    surv_treated = S1_avg,
    surv_control = S0_avg,
    ate = S1_avg - S0_avg
  )
}

plot_survdiff_gformula <- function(res) {
  df <- data.frame(time = res$time, ATE = res$ate)

  ggplot(df, aes(x = time, y = ATE)) +
    geom_line(color = "darkgreen", size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = "ATE(t) Curve via G-Formula",
      x = "Time",
      y = "Survival Difference: A=1 - A=0"
    ) +
    theme_minimal()
}

plot_survival_curves_gformula <- function(res) {
  df <- data.frame(
    time = res$time,
    Treated = res$surv_treated,
    Control = res$surv_control
  )

  df_long <- tidyr::pivot_longer(df, cols = c("Treated", "Control"),
                                 names_to = "Group", values_to = "Survival")

  ggplot(df_long, aes(x = time, y = Survival, color = Group)) +
    geom_line(size = 1) +
    labs(title = "Counterfactual Survival Curves via G-Formula",
         x = "Time", y = "Survival Probability") +
    scale_color_manual(values = c("blue", "red")) +
    theme_minimal()
}


# Fit model
model <- fit_xgboost(Surv(time, status) ~ A + X1 + X2, data = df)
times <- seq(0, 5, 0.1)

# Compute survival diff via G-formula
res_gf <- compute_survdiff_gformula(model, data = df, times = times, predict_fun = predict_xgboost)

# Plot ATE curve
plot_survdiff_gformula(res_gf)

# Plot counterfactual survival curves
plot_survival_curves_gformula(res_gf)

