fit_aftgee <- function(formula, data, corstr = "independence", ...) {
  stopifnot(requireNamespace("aftgee", quietly = TRUE))

  # create default id
  #id <- seq_len(nrow(data))

  model <- aftgee::aftgee(
    formula = formula,
    data = data,
    id = NULL, ## explicitly no clustering supported <=> one observation per subject
    corstr = corstr,
    ...
  )

  structure(list(
    model = model,
    learner = "aftgee",
    formula = formula,
    data = data
  ), class = "mlsurv_model", engine = "aftgee")
}


predict_aftgee <- function(object, newdata, times = NULL, ...) {
  stopifnot(object$learner == "aftgee")
  stopifnot(requireNamespace("aftgee", quietly = TRUE))

  beta <- object$model$coef.res
  X <- model.matrix(delete.response(terms(object$formula)), newdata)
  log_time_pred <- as.numeric(X %*% beta)

  if (is.null(times)) {
    # use quantiles or range from newdata to set a reasonable time grid
    t_max <- max(newdata[[all.vars(object$formula[[2]])[1]]], na.rm = TRUE)
    times <- seq(1, t_max, length.out = 20)
  }


  # approximate survival using lognormal-like assumption
  surv_mat <- sapply(times, function(t) {
    1 - pnorm((log(t) - log_time_pred))  # approximate S(t) from log-T model
  })

  colnames(surv_mat) <- paste0("t=", times)
  rownames(surv_mat) <- paste0("ID_", seq_len(nrow(newdata)))
  as.data.frame(surv_mat)
}



library(aftgee)
set.seed(123)

# simulated example
datgen <- function(n = 100, tau = 0.3, dim = 2) {
  x1 <- rbinom(dim * n, 1, 0.5)
  x2 <- rnorm(dim * n)
  e <- c(t(exp(MASS::mvrnorm(n = n, mu = rep(0, dim), Sigma = tau + (1 - tau) * diag(dim)))))
  tt <- exp(2 + x1 + x2 + e)
  cen <- runif(n, 0, 100)
  data.frame(Time = pmin(tt, cen), status = 1 * (tt < cen),
             x1 = x1, x2 = x2, id = rep(1:n, each = dim))
}
dat <- datgen(n = 50, dim = 1)

cv_survlearner(
  formula = Surv(Time, status) ~ x1 + x2,
  data = dat,
  fit_fun = fit_aftgee,
  pred_fun = predict_aftgee,
  times = c(10, 20, 40),
  metrics = c("cindex", "ibs"),
  folds = 5
)
