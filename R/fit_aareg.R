
#************** more save approach
#' Fit an Aalen Additive Survival Model (Functional Style)
#'
#' @param formula A survival formula, e.g., Surv(time, status) ~ x1 + x2
#' @param data A data frame containing the variables
#' @param max.time Optional maximum follow-up time
#' @param n.sim Number of simulations (for variance estimation)
#' @param resample.iid Logical or numeric for resampling type
#' @param ... Additional arguments passed to timereg::aalen
#' @return A named list representing the fitted model, with attribute "engine" = "aareg"
#' @export
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

  result <- list(
    fit = model,
    formula = formula,
    data = data
  )

  attr(result, "engine") <- "aareg"
  return(result)
}

#' Predict Survival Probabilities from Aalen Model (Functional Style)
#'
#' @param model A fitted model from `fit_aareg`
#' @param newdata New data frame to predict on
#' @param times Optional vector of times
#' @param ... Additional arguments passed to `predict.timereg`
#' @return A data frame of survival probabilities (rows = obs, cols = times)
#' @export
predict_aareg <- function(model, newdata, times = NULL, ...) {
  pout <- predict(model$fit, newdata = newdata, times = times, ...)

  surv_probs <- pout$S0
  if (is.null(surv_probs)) {
    stop("Survival probabilities (S0) not found in predict() output.")
  }

  colnames(surv_probs) <- paste0("t=", pout$time)
  rownames(surv_probs) <- paste0("ID_", seq_len(nrow(surv_probs)))

  as.data.frame(surv_probs)
}


library(survival)
library(timereg)

data(sTRACE)

mod <- fit_aareg(Surv(time, status == 9) ~ sex + diabetes + chf + vf,
                 data = sTRACE, max.time = 7)



pred2 <- predict_aareg(mod, newdata = sTRACE[1:3, ], times = c(0, 1, 3, 5))
pred2
##------------- add tuner

## fix max.time = max(data$time) and remove tuning. Otherwise this will gives inconsistent prediction horizons -> harder to compare within model
