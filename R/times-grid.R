#' Default evaluation-time grid
#'
#' Builds the time grid used by the verb-level functions ([predict()],
#' [evaluate()], [tune()], [compare()], [interpret()]) when the caller does not
#' pass `times` explicitly. The grid spans the observed event times so that
#' integrated metrics (IBS, IAE, ISE) and survival-curve summaries are computed
#' over the region where the data are informative.
#'
#' The chosen grid is recorded on every result object (`res$times`) and printed,
#' so that runs remain reproducible and comparable even though the default is
#' data-dependent.
#'
#' @param time Numeric vector of observed follow-up times.
#' @param status Optional 0/1 event indicator aligned with `time`. When supplied,
#'   the grid is built from event times only; otherwise all positive times are
#'   used.
#' @param n Integer; maximum number of grid points (default `50`).
#' @param range Two-value probability vector giving the lower and upper quantiles
#'   of the event-time distribution spanned by the grid (default
#'   `c(0.05, 0.95)`), which trims unstable tails.
#'
#' @return A numeric vector, strictly increasing and positive, of length at most
#'   `n`.
#'
#' @examples
#' default_times(veteran$time, veteran$status)
#' default_times(veteran$time, veteran$status, n = 10)
#' @export
default_times <- function(time, status = NULL, n = 50L, range = c(0.05, 0.95)) {
  time <- as.numeric(time)
  if (!is.null(status)) {
    status <- as.numeric(status)
    evt <- time[is.finite(time) & time > 0 & status == 1]
  } else {
    evt <- time[is.finite(time) & time > 0]
  }
  if (length(evt) < 2L) {
    evt <- sort(unique(time[is.finite(time) & time > 0]))
    if (length(evt) < 2L) stop("Cannot build a time grid: need >= 2 positive times.")
    return(evt)
  }
  n <- max(2L, as.integer(n))
  u <- sort(unique(evt))
  if (length(u) <= n) return(u)
  probs <- seq(range[1], range[2], length.out = n)
  g <- unique(as.numeric(stats::quantile(evt, probs = probs, type = 7, names = FALSE)))
  g <- g[g > 0]
  if (length(g) < 2L) g <- u[seq(1L, length(u), length.out = n)]
  sort(unique(g))
}

#' Resolve a caller-supplied `times` argument against the default grid
#'
#' @param times `NULL` or a numeric vector.
#' @param time,status Passed to [default_times()] when `times` is `NULL`.
#' @param n Grid size for the default.
#' @return A validated numeric vector (strictly increasing, positive).
#' @keywords internal
.resolve_times <- function(times, time, status = NULL, n = 50L) {
  if (is.null(times)) return(default_times(time, status, n = n))
  times <- sort(unique(as.numeric(times)))
  if (any(!is.finite(times)) || any(times <= 0)) {
    stop("`times` must be finite and positive.", call. = FALSE)
  }
  times
}
