#' Parse a survival outcome from a model formula
#'
#' Internal helper that extracts the time and status components from a formula
#' whose left-hand side is `Surv(...)` or `pkg::Surv(...)`.
#'
#' @param formula A model formula with a survival outcome on the left-hand side.
#' @param data A data.frame used to resolve status recoding expressions such as
#'   `status == 2`.
#'
#' @return A list with `time_col`, `status_col`, `event_value`, and
#'   `recode_status`.
#' @keywords internal
.parse_surv_formula <- function(formula, data) {
  tf <- terms(formula, data = data)
  outcome <- attr(tf, "variables")[[2]]

  head <- outcome[[1]]
  is_surv_call <-
    inherits(outcome, "call") &&
    (
      identical(head, as.name("Surv")) ||
        (
          is.call(head) &&
          as.character(head[[1]]) %in% c("::", ":::") &&
          identical(head[[3]], as.name("Surv"))
        )
    )

  if (!is_surv_call) {
    stop("Formula must begin with a Surv(...) outcome.")
  }

  if (length(outcome) != 3L) {
    stop("The outcome must be Surv(time, status) with exactly two arguments.",
         call. = FALSE)
  }
  time_expr <- outcome[[2]]
  status_expr <- outcome[[3]]

  if (is.call(status_expr) && identical(status_expr[[1]], as.name("=="))) {
    status_sym <- status_expr[[2]]
    event_value <- eval(status_expr[[3]], envir = data, enclos = parent.frame())
    recode_status <- TRUE
  } else {
    status_sym <- status_expr
    event_value <- 1
    recode_status <- FALSE
  }

  # The learners read the outcome columns from `data` by name, so a computed
  # time or status (e.g. Surv(time / 365.25, status)) would silently use the
  # raw column. Ask for a plain column instead.
  for (s in list(time_expr, status_sym)) {
    if (!is.name(s)) {
      stop("Surv() arguments must be column names, not expressions (got `",
           deparse(s), "`). Compute it as a column first, e.g. ",
           "data$years <- data$time / 365.25.", call. = FALSE)
    }
    if (!as.character(s) %in% names(data)) {
      stop("Column `", as.character(s), "` is not in `data`.", call. = FALSE)
    }
  }
  time_col <- as.character(time_expr)
  status_col <- as.character(status_sym)

  list(
    time_col = time_col,
    status_col = status_col,
    event_value = event_value,
    recode_status = recode_status
  )
}

## Predictors of a fitted model, in formula order, restricted to `cols` (the
## columns of the data being explained). The explanation methods use this so
## that they never perturb or report the outcome or columns the model did not
## use. A `.` on the right-hand side stands for every non-outcome column.
#' @keywords internal
#' @noRd
.model_features <- function(model, cols) {
  f <- model$formula
  parsed <- .parse_surv_formula(f, model$data)
  outcome <- c(parsed$time_col, parsed$status_col)
  v <- all.vars(f[[3]])
  if ("." %in% v) v <- union(setdiff(v, "."), names(model$data))
  v <- setdiff(v, outcome)
  v[v %in% cols]
}

## Engines such as ranger and aorsf only accept plain column names on the
## right-hand side. When the formula has transformed terms (log1p(x), I(x^2),
## factor(x)), these helpers evaluate each predictor expression into its own
## column with a syntactic name, so the engine sees plain columns. The same
## spec is applied to new data at prediction time. Returns NULL when every
## predictor is already a plain column (nothing to do).
#' @keywords internal
#' @noRd
.derive_predictors <- function(formula, data) {
  tt <- stats::terms(formula, data = data)
  vars <- as.list(attr(tt, "variables"))[-1L]
  resp <- attr(tt, "response")
  if (resp > 0L) vars <- vars[-resp]
  if (!length(vars) || all(vapply(vars, is.name, logical(1)))) return(NULL)
  labels <- vapply(vars, function(e) paste(deparse(e, width.cutoff = 500L), collapse = ""),
                   character(1))
  outcome <- all.vars(formula[[2L]])
  nm <- make.names(c(outcome, labels), unique = TRUE)[-seq_along(outcome)]
  spec <- list(exprs = vars, labels = labels, names = nm, outcome = outcome,
               env = environment(formula))
  spec$formula <- stats::reformulate(nm, response = formula[[2L]], env = environment(formula))
  spec
}

#' @keywords internal
#' @noRd
.apply_derived <- function(spec, data, outcome = TRUE) {
  cols <- lapply(seq_along(spec$exprs), function(j) {
    v <- eval(spec$exprs[[j]], as.data.frame(data), spec$env)
    if (!is.null(dim(v)) && NCOL(v) > 1L) {
      stop("Predictor `", spec$labels[j], "` gives a matrix; this engine needs ",
           "one column per term. Create the columns in the data first.", call. = FALSE)
    }
    if (!is.null(dim(v))) v <- as.vector(v)
    v
  })
  names(cols) <- spec$names
  out <- as.data.frame(cols, stringsAsFactors = FALSE, optional = TRUE)
  keep <- intersect(spec$outcome, names(data))
  if (outcome && length(keep)) out <- cbind(as.data.frame(data)[keep], out)
  out
}
