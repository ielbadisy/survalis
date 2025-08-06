
fit_stpm2 <- function(formula, data, df = 3, ...) {

  stopifnot(requireNamespace("rstpm2", quietly = TRUE))

  model <- rstpm2::stpm2(formula = formula, data = data, df = df, ...)

  structure(
    list(
      model = model,
      learner = "stpm2", # Parametric and penalised generalised survival models
      formula = formula,
      data = data,
      time = all.vars(formula)[[2]],
      status = all.vars(formula)[[3]]
    ),
    class = "mlsurv_model",
    engine = "rstpm2"
  )
}


predict_stpm2 <- function(object, newdata, times, ...) {

  if (!is.null(object$learner) && object$learner != "stpm2") {
    warning("Object passed to predict_stpm2() may not come from fit_stpm2().")
    }

  stopifnot(requireNamespace("rstpm2", quietly = TRUE))

  n <- nrow(newdata)
  expanded_data <- newdata[rep(seq_len(n), each = length(times)), , drop = FALSE]
  expanded_data$time <- rep(times, times = n)

  surv_probs <- predict(
    object$model,
    newdata = expanded_data,
    type = "surv",
    newtime = expanded_data$time
    )

  survmat <- matrix(surv_probs, nrow = n, byrow = TRUE)
  colnames(survmat) <- paste0("t=", times)

  as.data.frame(survmat)
}


library(survival)
library(rstpm2)

data(veteran)

mod_stpm2 <- fit_stpm2(Surv(time, status) ~ age + karno + celltype, data = veteran, df = 4)
predict_stpm2(mod_stpm2, newdata = veteran[1:5, ], times = c(100, 200, 300))



summary(mod_stpm2)
