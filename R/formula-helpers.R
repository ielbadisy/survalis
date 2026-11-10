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
