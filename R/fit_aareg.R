fit_aareg <- function(formula, data, max.time = NULL, n.sim = 0, resample.iid = 1, ...) {
  stopifnot(requireNamespace("timereg", quietly = TRUE))
  
  model <- timereg::aalen(
    formula = formula,
    data = data,
    max.time = max.time,
    n.sim = n.sim,
    resample.iid = resample.iid,
    ...
  )
  
  structure(list(
    model = model,
    learner = "aareg",
    formula = formula,
    data = data
  ), class = "mlsurv")
}


predict_aareg <- function(object, newdata, times = NULL, ...) {
  stopifnot(object$learner == "aareg")
  
  pout <- predict(object$model, newdata = newdata, times = times, ...)
  
  surv_probs <- pout$S0
  if (is.null(surv_probs)) {
    stop("Survival probabilities (S0) not found in predict() output.")
  }
  
  # label columns by  -> should be unifed later through a common helper
  colnames(surv_probs) <- paste0("t=", pout$time)
  rownames(surv_probs) <- paste0("ID_", seq_len(nrow(surv_probs)))
  
  as.data.frame(surv_probs)
}


# TEST

library(survival)
library(timereg)

data(sTRACE)

## test
mod_aareg <- fit_aareg(Surv(time, status == 9) ~ sex + diabetes + chf + vf,
                       data = sTRACE, max.time = 7)

## test predict
pred_aareg <- predict_aareg(mod_aareg, newdata = sTRACE[1:3, ], time = c(0, 1, 10)) 

pred_aareg
