#' Fit a flexible parametric survival model using rstpm2
#'
#' @param formula A survival formula (e.g., Surv(time, status) ~ x1 + x2)
#' @param data A data.frame with variables
#' @param df Degrees of freedom for spline basis (default = 3)
#' @param ... Additional arguments passed to rstpm2::stpm2
#'
#' @return An object of class 'mlsurv_model'
#' @export
fit_rstpm2 <- function(formula, data, df = 3, ...) {
  stopifnot(requireNamespace("rstpm2", quietly = TRUE))

  model <- rstpm2::stpm2(formula = formula, data = data, df = df, ...)

  structure(
    list(
      model = model,
      learner = "rstpm2",
      formula = formula,
      data = data,
      time = all.vars(formula)[[2]],
      status = all.vars(formula)[[3]]
    ),
    class = "mlsurv_model",
    engine = "rstpm2"
  )
}


predict_rstpm2 <- function(object, newdata, times, ...) {
  stopifnot(object$learner == "rstpm2")
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

  surv_matrix <- matrix(surv_probs, nrow = n, byrow = TRUE)
  colnames(surv_matrix) <- paste0("t=", times)
  rownames(surv_matrix) <- paste0("ID_", seq_len(n))

  as.data.frame(surv_matrix)
}


library(survival)
library(rstpm2)

data(veteran)

mod_rstpm2 <- fit_rstpm2(Surv(time, status) ~ age + karno + celltype, data = veteran, df = 4)
predict_rstpm2(mod_rstpm2, newdata = veteran[1:5, ], times = c(100, 200, 300))
