compute_rmst <- function(surv_matrix, times) {
  # Ensure numeric matrix
  surv_matrix <- as.matrix(surv_matrix)

  # Trapezoidal rule: midpoint average × time difference
  delta_t <- diff(times)  # length T-1
  mid_surv <- (surv_matrix[, -1, drop = FALSE] + surv_matrix[, -ncol(surv_matrix), drop = FALSE]) / 2
  rmst <- mid_surv %*% delta_t  # matrix mult: (n × T-1) × (T-1 × 1) → (n × 1)
  as.numeric(rmst)
}



estimate_ate_rmst <- function(rmst, A, weights) {
  treated <- A == 1
  control <- A == 0

  rmst_treated <- sum(rmst[treated] * weights[treated]) / sum(weights[treated])
  rmst_control <- sum(rmst[control] * weights[control]) / sum(weights[control])

  list(
    rmst_treated = rmst_treated,
    rmst_control = rmst_control,
    ate_rmst = rmst_treated - rmst_control
  )
}

# 1. Fit any survival model
model <- fit_xgboost(Surv(time, status) ~ A + X1 + X2, data = df)

# 2. Predict survival
times <- seq(0, 5, by = 0.1)
survmat <- predict_xgboost(model, newdata = df, times = times)

# 3. Compute IPTW weights
weights <- compute_iptw_weights(A = df$A, X = df[, c("X1", "X2")])

# 4. Compute RMST for each subject
rmst_i <- compute_rmst(survmat, times)

# 5. Estimate ATE on RMST
rmst_ate <- estimate_ate_rmst(rmst_i, A = df$A, weights = weights)
print(rmst_ate)


library(ggplot2)

df$rmst <- rmst_i
df$weights <- weights

ggplot(df, aes(x = rmst, fill = factor(A))) +
  geom_density(aes(weight = weights), alpha = 0.5, color = NA) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8"),
                    labels = c("Control (A=0)", "Treated (A=1)")) +
  labs(title = "Weighted RMST Distribution by Treatment",
       x = "Restricted Mean Survival Time",
       y = "Weighted Density", fill = "Group") +
  theme_minimal()

