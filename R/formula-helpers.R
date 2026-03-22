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

  time_col <- as.character(outcome[[2]])
  status_expr <- outcome[[3]]

  if (is.call(status_expr) && identical(status_expr[[1]], as.name("=="))) {
    status_col <- as.character(status_expr[[2]])
    event_value <- eval(status_expr[[3]], envir = data, enclos = parent.frame())
    recode_status <- TRUE
  } else {
    status_col <- as.character(status_expr)
    event_value <- 1
    recode_status <- FALSE
  }

  list(
    time_col = time_col,
    status_col = status_col,
    event_value = event_value,
    recode_status = recode_status
  )
}
