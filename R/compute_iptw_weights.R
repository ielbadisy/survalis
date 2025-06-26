compute_iptw_weights <- function(A, X, method = "glm") {
  if (method == "glm") {
    ps_model <- glm(A ~ ., data = X, family = binomial)
    ps <- predict(ps_model, type = "response")
  } else {
    stop("Only 'glm' method is currently implemented")
  }
  pA <- mean(A)
  weights <- ifelse(A == 1, pA / ps, (1 - pA) / (1 - ps))
  return(weights)
}


estimate_ate_surv_iptw <- function(model, surv_matrix, A, weights, times) {
  # surv_matrix: survival probs from predict_*(), shape n x length(times)
  treated <- A == 1
  control <- A == 0

  # Weighted average survival probabilities for each group
  surv_treated  <- colSums(surv_matrix[treated, ] * weights[treated]) / sum(weights[treated])
  surv_control  <- colSums(surv_matrix[control, ] * weights[control]) / sum(weights[control])

  # ATE = S(t|A=1) - S(t|A=0)
  ate_curve <- surv_treated - surv_control

  return(list(
    time = times,
    surv_treated = surv_treated,
    surv_control = surv_control,
    ate = ate_curve
  ))
}


set.seed(2025)
n <- 1000  # sample size

# Covariates
X1 <- rnorm(n)
X2 <- rbinom(n, 1, 0.5)

# Propensity score (treatment assignment)
logit_ps <- 0.5 * X1 - 0.5 * X2
ps <- plogis(logit_ps)
A <- rbinom(n, 1, ps)  # binary treatment

# Baseline hazard
lambda <- 0.1
hr_A <- 0.6  # hazard ratio for treatment
hr_X1 <- 1.3
hr_X2 <- 1.5

# Linear predictor for hazard
linpred <- log(hr_A) * A + log(hr_X1) * X1 + log(hr_X2) * X2

# Survival times from exponential distribution
T_true <- rexp(n, rate = lambda * exp(linpred))

# Censoring
C <- rexp(n, rate = 0.02)  # administrative censoring
time <- pmin(T_true, C)
status <- as.integer(T_true <= C)

# Final dataset
df <- data.frame(
  time = time,
  status = status,
  A = A,
  X1 = X1,
  X2 = X2
)

head(df)


# This should now work correctly:
model <- fit_xgboost(Surv(time, status) ~ A + X1 + X2, data = df)


# Assume you already have A, X, time, status
model <- fit_xgboost(Surv(time, status) ~ A + X1 + X2, data = df)

# Predict survival probs
times <- seq(0, 5, by = 0.1)
survmat <- predict_xgboost(model, newdata = df, times = times)

# Compute weights and estimate ATE
weights <- compute_iptw_weights(A = df$A, X = df[, c("X1", "X2")])
ate_result <- estimate_ate_surv_iptw(model, survmat, A = df$A, weights, times)

# Plot ATE curve
matplot(ate_result$time, cbind(ate_result$surv_treated, ate_result$surv_control),
        type = "l", col = 1:2, lty = 1, ylab = "Survival", xlab = "Time")
legend("topright", legend = c("A=1", "A=0"), col = 1:2, lty = 1)

