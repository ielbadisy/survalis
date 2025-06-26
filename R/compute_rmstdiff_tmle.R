compute_rmstdiff_tmle <- function(model, data, A, W, times, predict_fun, ps_fun = NULL) {
  data1 <- data; data1$A <- 1
  data0 <- data; data0$A <- 0

  S1 <- predict_fun(model, newdata = data1, times = times)
  S0 <- predict_fun(model, newdata = data0, times = times)

  rmst1 <- compute_rmst(S1, times)
  rmst0 <- compute_rmst(S0, times)

  # Propensity score estimation
  if (is.null(ps_fun)) {
    ps_model <- glm(A ~ ., data = as.data.frame(W), family = binomial)
    ps <- predict(ps_model, type = "response")
  } else {
    ps <- ps_fun(A, W)
  }

  H1 <- A / ps
  H0 <- (1 - A) / (1 - ps)

  eps1 <- coef(glm(rmst1 ~ -1 + offset(rmst1) + H1, family = gaussian()))
  eps0 <- coef(glm(rmst0 ~ -1 + offset(rmst0) + H0, family = gaussian()))

  rmst_star_1 <- rmst1 + eps1 * H1
  rmst_star_0 <- rmst0 + eps0 * H0

  list(
    method = "tmle",
    estimand = "rmst_difference",
    tau = max(times),
    rmst_treated = mean(rmst_star_1),
    rmst_control = mean(rmst_star_0),
    ate_rmst = mean(rmst_star_1 - rmst_star_0)
  )
}



plot_rmstdiff_tmle <- function(res) {
  df <- data.frame(
    Group = c("Treated", "Control"),
    RMST = c(res$rmst_treated, res$rmst_control)
  )

  ggplot(df, aes(x = Group, y = RMST, fill = Group)) +
    geom_col(width = 0.5, alpha = 0.8) +
    labs(
      title = sprintf("TMLE RMST Comparison (τ = %.1f)", res$tau),
      y = "Restricted Mean Survival Time", x = ""
    ) +
    theme_minimal() +
    scale_fill_manual(values = c("blue", "red")) +
    theme(legend.position = "none")
}



model <- fit_xgboost(Surv(time, status) ~ A + X1 + X2, data = df)
times <- seq(0, 5, 0.1)

res_tmle <- compute_rmstdiff_tmle(
  model = model,
  data = df,
  A = df$A,
  W = df[, c("X1", "X2")],
  times = times,
  predict_fun = predict_xgboost
)

plot_rmstdiff_tmle(res_tmle)
