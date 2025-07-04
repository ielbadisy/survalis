#' Fit Parametric Survival Model (AFT)
fit_survreg <- function(formula, data,
                        dist = "weibull",
                        scale = 1,
                        parms = list(),
                        ...) {
  stopifnot(requireNamespace("rms", quietly = TRUE))
  stopifnot(requireNamespace("Hmisc", quietly = TRUE))

  model <- rms::psm(
    formula = formula,
    data = data,
    dist = dist,
    scale = scale,
    parms = parms,
    control = survival::survreg.control(...)
  )

  structure(list(
    model = model,
    learner = "survreg",
    formula = formula,
    data = data,
    time = all.vars(formula)[[2]],
    status = all.vars(formula)[[3]]
  ), class = "mlsurv_model", engine = "survreg")
}


#' Predict Survival Probabilities from Parametric Model
predict_survreg <- function(object, newdata, times) {
  stopifnot(object$learner == "survreg")
  stopifnot(requireNamespace("rms", quietly = TRUE))

  newdata <- as.data.frame(newdata)

  pred <- rms::survest(object$model, newdata = newdata, times = times, conf.int = FALSE)

  mat <- if (is.list(pred)) {
    matrix(pred$surv, nrow = 1)
  } else {
    pred
  }

  colnames(mat) <- paste0("t=", times)
  rownames(mat) <- paste0("ID_", seq_len(nrow(mat)))
  as.data.frame(mat)
}

# Example
library(survival)
data(veteran)

mod_weib <- fit_survreg(Surv(time, status) ~ age + karno + celltype,
                        data = veteran, dist = "weibull")

pred <- predict_survreg(mod_weib, newdata = veteran[1:5, ], times = c(100, 200, 300))
print(round(pred, 3))

## NOT TUNABLE
