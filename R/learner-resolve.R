# Central learner lookup used by every verb.
#
# `model` is either a built-in learner id (a string from list_survlearners())
# or a custom learner: a list with `fit = function(formula, data, ...)` and
# `predict = function(object, newdata, times)`, the latter returning an
# n x length(times) matrix of survival probabilities. An optional `name`
# labels the learner in printed output.

#' @noRd
.is_custom_learner <- function(model) {
  is.list(model) && !is.null(model$fit) && !is.null(model$predict)
}

#' @noRd
.resolve_learner <- function(model) {
  if (.is_custom_learner(model)) {
    if (!is.function(model$fit) || !is.function(model$predict)) {
      stop("A custom learner must be list(fit = <function>, predict = <function>).",
           call. = FALSE)
    }
    id <- .null_default(model$name, "custom")
    user_fit <- model$fit
    user_predict <- model$predict
    pred_fun <- function(object, newdata, times) {
      out <- user_predict(object$model, newdata = newdata, times = times)
      out <- as.matrix(out)
      if (!is.numeric(out) || nrow(out) != nrow(newdata) || ncol(out) != length(times)) {
        stop("The custom learner's predict() must return a numeric matrix with ",
             nrow(newdata), " rows and ", length(times), " columns (one per time).",
             call. = FALSE)
      }
      colnames(out) <- paste0("t=", times)
      out
    }
    fit_fun <- function(formula, data, ...) {
      raw <- user_fit(formula, data, ...)
      structure(list(model = raw, learner = id, formula = formula, data = data,
                     predict_fun = pred_fun),
                class = "mlsurv_model")
    }
    return(list(id = id, fit = fit_fun, predict = pred_fun, tune = NULL,
                custom = TRUE, def = model))
  }

  known <- list_survlearners()$learner
  if (length(model) != 1L || !is.character(model) || !model %in% known) {
    lab <- if (is.character(model)) paste(model, collapse = ", ") else class(model)[1]
    stop("Unknown model ", sQuote(lab), ". Use an id from list_survlearners() ",
         "or a custom learner list(fit = , predict = ).", call. = FALSE)
  }
  ns <- asNamespace("survalis")
  list(id = model,
       fit = get(paste0("fit_", model), envir = ns),
       predict = get(paste0("predict_", model), envir = ns),
       tune = get0(paste0("tune_", model), envir = ns, inherits = FALSE),
       custom = FALSE, def = model)
}

# Prediction function for a fitted model object (built-in or custom).
#' @noRd
.predict_fun <- function(model) {
  if (is.function(model$predict_fun)) return(model$predict_fun)
  if (is.null(model$learner)) {
    stop("`model` has no learner id; fit it with fit().", call. = FALSE)
  }
  fn <- get0(paste0("predict_", model$learner), envir = asNamespace("survalis"),
             mode = "function")
  if (is.null(fn)) stop("No predict method for learner ", sQuote(model$learner), ".",
                        call. = FALSE)
  fn
}

# Learner label for printing: id for built-ins, `name` (or "custom") otherwise.
#' @noRd
.learner_label <- function(model) {
  if (.is_custom_learner(model)) .null_default(model$name, "custom") else model
}
