#' Fit a survival model using <LEARNER>
#'
#' @param formula A survival formula (e.g., Surv(time, status) ~ x1 + x2)
#' @param data A data.frame containing the variables
#' @param ... Additional learner-specific hyperparameters
#'
#' @return A list of class 'mlsurv_model' with attributes for unified prediction
#' @export
fit_<learner> <- function(formula, data, ...) {
  stopifnot(requireNamespace("<PACKAGE>", quietly = TRUE))

  # Fit model using the specific survival learner
  model <- <PACKAGE>::<learner_function>(
    formula = formula,
    data = data,
    ...
  )

  structure(list(
    model = model,
    learner = "<learner>",     # short name, e.g., "aareg", "ranger", etc.
    formula = formula,
    data = data
  ), class = "mlsurv_model", engine = "<learner>")
}



#' Predict survival probabilities using a <LEARNER> model
#'
#' @param object A model fitted with fit_<learner>()
#' @param newdata A new data.frame for prediction
#' @param times A numeric vector of time points for prediction
#' @param ... Extra arguments passed to the prediction method
#'
#' @return A data.frame: rows = individuals, columns = t=100, t=200, ...
#' @export
predict_<learner> <- function(object, newdata, times, ...) {
  stopifnot(object$learner == "<learner>")
  stopifnot(requireNamespace("<PACKAGE>", quietly = TRUE))

  # === Core prediction logic ===
  pred_output <- <PACKAGE>::<predict_function>(
    object$model,
    newdata = newdata,
    ...
  )

  # === Extract or interpolate survival probabilities ===
  # Ensure you return a matrix or data.frame of shape [nrow(newdata), length(times)]
  surv_probs <- <custom code for extracting or interpolating predicted survival>

    # === Format output ===
    colnames(surv_probs) <- paste0("t=", times)
  rownames(surv_probs) <- paste0("ID_", seq_len(nrow(newdata)))


  as.data.frame(surv_probs)
}



mod <- fit_<learner>(Surv(time, status) ~ x1 + x2, data = your_data)
pred <- predict_<learner>(mod, newdata = your_data[1:5, ], times = c(100, 200, 300))
print(pred)

