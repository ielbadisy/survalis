
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
