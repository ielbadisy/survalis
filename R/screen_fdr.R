#' FDR-Based Variable Screening via Univariate Cox Models
#'
#' Screens candidate covariates for a right-censored outcome by fitting a
#' univariate Cox model \eqn{\lambda(t \mid X_j) = \lambda_0(t)\exp(\beta_j X_j)}
#' per covariate, Wald-testing \eqn{H_0: \beta_j = 0}, and controlling the false
#' discovery rate (FDR) across all tests via the Benjamini-Hochberg (1995)
#' step-up procedure.
#'
#' @param formula A model formula with a \code{Surv(time, status)} outcome on
#'   the left-hand side and candidate screening covariates on the right-hand
#'   side (as used throughout \pkg{survalis}, e.g. \code{\link{fit_coxph}}).
#' @param data A data.frame containing the variables in \code{formula}.
#' @param alpha Target FDR level in (0, 1). Covariates with BH-adjusted
#'   p-value (q-value) \code{<= alpha} are retained. Default \code{0.1}.
#' @param minscreen Minimum number of covariates to retain. If fewer than
#'   \code{minscreen} covariates pass the FDR threshold, the smallest-p-value
#'   covariates are added back until \code{minscreen} is reached (subject to
#'   the number of candidate covariates available). Default \code{5}.
#'
#' @details
#' For each candidate covariate \eqn{X_j}, \pkg{survalis} fits
#' \code{coxph(Surv(time, status) ~ X_j, data = data)} and extracts the Wald
#' p-value for \eqn{\beta_j}. Raw p-values \eqn{p_1, \dots, p_p} are ordered
#' \eqn{p_{(1)} \le \cdots \le p_{(p)}} and the BH-adjusted q-value for the
#' \eqn{k}-th smallest p-value is
#' \deqn{q_{(k)} = \min_{k \le l \le p} \left\{ \frac{p \cdot p_{(l)}}{l} \right\},}
#' guaranteeing \eqn{\mathrm{E}[\text{FDR}] \le \alpha} under the retained set
#' \eqn{\{j : q_j \le \alpha\}} (Benjamini and Hochberg 1995).
#'
#' @return A character vector of selected covariate names, with an attribute
#'   \code{"screen_table"} holding the full \code{data.frame} of
#'   \code{feature}, \code{p_value}, \code{q_value}, and \code{selected}.
#'
#' @references
#' Benjamini Y, Hochberg Y (1995). "Controlling the False Discovery Rate: A
#' Practical and Powerful Approach to Multiple Testing." \emph{Journal of the
#' Royal Statistical Society: Series B}, 57(1), 289-300.
#'
#' @examples
#' selected <- screen_fdr(Surv(time, status) ~ age + karno + diagtime + prior,
#'                        data = veteran, alpha = 0.2, minscreen = 2)
#' selected
#' attr(selected, "screen_table")
#'
#' @export
screen_fdr <- function(formula, data, alpha = 0.1, minscreen = 5) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame.")
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0,1).")
  if (!is.finite(minscreen) || minscreen < 0) stop("`minscreen` must be a non-negative integer.")

  parsed <- .parse_surv_formula(formula, data)
  time <- data[[parsed$time_col]]
  status <- if (parsed$recode_status) {
    as.numeric(data[[parsed$status_col]] == parsed$event_value)
  } else {
    data[[parsed$status_col]]
  }

  candidates <- attr(terms(formula, data = data), "term.labels")
  candidates <- candidates[candidates %in% colnames(data)]
  if (!length(candidates)) stop("No usable candidate covariates found on the right-hand side of `formula`.")

  y <- survival::Surv(time, status)

  p_values <- vapply(candidates, function(feat) {
    fit <- tryCatch(
      survival::coxph(y ~ data[[feat]]),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NA_real_)
    s <- summary(fit)$coefficients
    unname(s[1, "Pr(>|z|)"])
  }, numeric(1))

  q_values <- stats::p.adjust(p_values, method = "BH")

  screen_table <- data.frame(
    feature = candidates,
    p_value = as.numeric(p_values),
    q_value = as.numeric(q_values),
    stringsAsFactors = FALSE
  )
  screen_table <- screen_table[order(screen_table$p_value, na.last = TRUE), , drop = FALSE]
  screen_table$selected <- !is.na(screen_table$q_value) & screen_table$q_value <= alpha

  n_selected <- sum(screen_table$selected)
  if (n_selected < minscreen) {
    fill <- which(!screen_table$selected & !is.na(screen_table$p_value))
    n_fill <- min(length(fill), minscreen - n_selected)
    if (n_fill > 0) screen_table$selected[fill[seq_len(n_fill)]] <- TRUE
  }

  screen_table <- screen_table[order(screen_table$feature), , drop = FALSE]
  rownames(screen_table) <- NULL

  selected <- screen_table$feature[screen_table$selected]
  attr(selected, "screen_table") <- screen_table
  selected
}
