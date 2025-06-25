compute_survdiff_iptw <- function(model, data, A, W, times, predict_fun) {
  weights <- compute_iptw_weights(A, W)
  survmat <- predict_fun(model, newdata = data, times = times)

  treated <- A == 1
  control <- A == 0

  S1 <- colSums(survmat[treated, , drop = FALSE] * weights[treated]) / sum(weights[treated])
  S0 <- colSums(survmat[control, , drop = FALSE] * weights[control]) / sum(weights[control])

  list(
    method = "iptw",
    estimand = "survival_difference",
    time = times,
    surv_treated = S1,
    surv_control = S0,
    ate = S1 - S0
  )
}




plot_survdiff_iptw <- function(res) {
  if (!all(c("time", "ate") %in% names(res))) {
    stop("Result must include 'time' and 'ate'")
  }

  df <- data.frame(time = res$time, ATE = res$ate)

  ggplot(df, aes(x = time, y = ATE)) +
    geom_line(color = "darkblue", size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = "ATE(t) Curve via IPTW",
      x = "Time",
      y = "Survival Difference: A=1 - A=0"
    ) +
    theme_minimal()
}


plot_survival_curves_iptw <- function(res) {
  df <- data.frame(
    time = res$time,
    Treated = res$surv_treated,
    Control = res$surv_control
  )

  df_long <- tidyr::pivot_longer(df, cols = c("Treated", "Control"), names_to = "Group", values_to = "Survival")

  ggplot(df_long, aes(x = time, y = Survival, color = Group)) +
    geom_line(size = 1) +
    labs(title = "IPTW-Weighted Survival Curves",
         x = "Time", y = "Survival Probability") +
    scale_color_manual(values = c("blue", "red")) +
    theme_minimal()
}



# Simulate data
set.seed(42)
n <- 500
X1 <- rnorm(n)
X2 <- rbinom(n, 1, 0.5)
A <- rbinom(n, 1, plogis(0.5 * X1 - 0.5 * X2))
hr_A <- ifelse(X2 == 1, 0.4, 0.9)
lp <- log(hr_A) * A + log(1.3) * X1 + log(1.5) * X2
time <- rexp(n, rate = 0.1 * exp(lp))
cens <- rexp(n, rate = 0.02)
status <- as.integer(time <= cens)
time <- pmin(time, cens)
df <- data.frame(time, status, A, X1, X2)

# Fit model and predict
model <- fit_xgboost(Surv(time, status) ~ A + X1 + X2, data = df)
times <- seq(0, 5, by = 0.1)

# Compute and plot
res <- compute_survdiff_iptw(model, df, df$A, df[, c("X1", "X2")], times, predict_xgboost)
plot_survdiff_iptw(res)


plot_survival_curves_iptw(res)
