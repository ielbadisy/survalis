# Central learner lookup used by every verb.
#
# `model` is either a built-in learner id (a string from listlearners())
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
      # same clamping and monotone repair as the built-in learners, so every
      # verb scores the same values
      as.matrix(.finalize_survmat(out, times = times))
    }
    fit_fun <- function(formula, data, ...) {
      s <- .std_outcome(formula, data)
      raw <- user_fit(formula = s$formula, data = s$data, ...)
      structure(list(model = raw, learner = id, formula = formula, data = data,
                     predict_fun = pred_fun),
                class = "mlsurv_model")
    }
    return(list(id = id, fit = fit_fun, predict = pred_fun, tune = NULL,
                custom = TRUE, def = model))
  }

  known <- listlearners()$learner
  if (length(model) != 1L || !is.character(model) || !model %in% known) {
    lab <- if (is.character(model)) paste(model, collapse = ", ") else class(model)[1]
    stop("Unknown model ", sQuote(lab), ". Use an id from listlearners() ",
         "or a custom learner list(fit = , predict = ).", call. = FALSE)
  }
  ns <- asNamespace("survalis")
  engine_fit <- get(paste0("fit_", model), envir = ns)
  engine_tune <- get0(paste0("tune_", model), envir = ns, inherits = FALSE)
  list(id = model,
       fit = function(formula, data, ...) {
         s <- .std_outcome(formula, data)
         engine_fit(formula = s$formula, data = s$data, ...)
       },
       predict = get(paste0("predict_", model), envir = ns),
       tune = engine_tune,
       custom = FALSE, def = model)
}

# The learners read the outcome columns by name and differ in how they treat
# a 1/2 or logical status, or ignore a `status == k` recoding. Hand every
# learner the same outcome: a 0/1 event column and a plain Surv(time, status).
.std_outcome <- function(formula, data) {
  parsed <- .parse_surv_formula(formula, data)
  data[[parsed$status_col]] <- .recode_status_vec(parsed, data)
  formula[[2]][[3]] <- as.name(parsed$status_col)
  list(formula = formula, data = data)
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
