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


  if (!is.null(object$learner) && object$learner != "aftgee") {
    warning("Object passed to predict_aftgee() may not come from fit_aftgee().")
  }

  stopifnot(requireNamespace("aftgee", quietly = TRUE))

  beta <- object$model$coef.res
  X <- model.matrix(delete.response(terms(object$formula)), newdata)
  log_time_pred <- as.numeric(X %*% beta)

  if (is.null(times)) {
    t_max <- max(newdata[[all.vars(object$formula[[2]])[1]]], na.rm = TRUE)
    times <- seq(1, t_max, length.out = 20)
    }


  # approximate survival using lognormal-like assumption
  survmat <- sapply(times, function(t) {
    1 - pnorm((log(t) - log_time_pred))  # approximate S(t) from log-T model
  })

  colnames(survmat) <- paste0("t=", times)
  as.data.frame(survmat)
}



library(aftgee)

library(survival)

mod <- fit_aftgee(Surv(time, status) ~ trt + karno + age,
                 data = veteran)



pred <- predict_aftgee(mod, newdata = veteran, times = 0:100)
pred



cv_survlearner(
  formula = Surv(time, status) ~ trt + karno + age,
  data = veteran,
  fit_fun = fit_aftgee,
  pred_fun = predict_aftgee,
  times = c(0, 20, 40),
  metrics = c("cindex", "ibs"),
  folds = 5
)
