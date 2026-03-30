#' Validate and standardize a survival-probability matrix (survmat)
#'
#' Internal utility to validate the basic structure of a `surv_mat` and coerce it
#' to a numeric matrix. Optionally checks column names against `times`.
#'
#' @param S A `surv_mat`: numeric matrix or data.frame of survival probabilities.
#' @param times Optional numeric vector of time points. If provided, requires
#'   `length(times) == ncol(S)` and enforces the column names `"t=<time>"`.
#' @param strict_colnames Logical; if `TRUE` and `times` is provided, checks that
#'   existing colnames are exactly `"t=<time>"` (after coercion to character).
#'   If `FALSE`, colnames are overwritten to match the contract.
#' @param clamp Logical; if `TRUE`, clamps `S` into `[0,1]`.
#'
#' @return A numeric matrix.
#' @keywords internal
.survmat_as_matrix <- function(S,
                               times = NULL,
                               strict_colnames = FALSE,
                               clamp = FALSE) {
  if (is.data.frame(S)) S <- as.matrix(S)
  if (!is.matrix(S) || !is.numeric(S)) {
    stop("`S` must be a numeric matrix or data.frame (surv_mat).")
  }

  if (clamp) S <- pmin(pmax(S, 0), 1)

  if (!is.null(times)) {
    times <- as.numeric(times)
    if (length(times) != ncol(S)) {
      stop("`times` must have length equal to ncol(S).")
    }
    if (any(!is.finite(times)) || any(times < 0)) {
      stop("`times` must be non-negative and finite.")
    }
    if (is.unsorted(times, strictly = TRUE)) {
      stop("`times` must be strictly increasing.")
    }

    target <- paste0("t=", times)
    if (strict_colnames) {
      cn <- colnames(S)
      if (is.null(cn) || !identical(cn, target)) {
        stop("Column names of `S` must be exactly: ", paste(target, collapse = ", "))
      }
    } else {
      colnames(S) <- target
    }
  }

  S
}

#' Standardize learner prediction output to the survmat contract
#'
#' Internal utility used by `predict_*()` methods and evaluators to coerce raw
#' predictions to a validated survival-probability matrix contract.
#'
#' @param S Raw survival predictions as matrix, data.frame, or vector.
#' @param times Numeric vector of evaluation times aligned with columns of `S`.
#' @param clamp Logical; if `TRUE`, clamp values to `[0, 1]`.
#' @param enforce_monotone Logical; if `TRUE`, enforce non-increasing survival
#'   over increasing `times` via cumulative minima.
#'
#' @return A base data.frame with columns named `t=<time>`.
#' @keywords internal
.finalize_survmat <- function(S,
                              times,
                              clamp = TRUE,
                              enforce_monotone = TRUE) {
  if (is.null(times) || !length(times)) {
    stop("`times` must be a non-empty numeric vector.")
  }

  times <- as.numeric(times)
  if (any(!is.finite(times)) || any(times < 0)) {
    stop("`times` must contain finite, non-negative values.")
  }

  if (is.vector(S) && !is.list(S)) {
    S <- matrix(S, nrow = 1)
  }

  S <- .survmat_as_matrix(S, clamp = FALSE)

  if (ncol(S) != length(times)) {
    stop("Prediction output must have one column per requested time.")
  }

  ord <- order(times)
  S_ord <- S[, ord, drop = FALSE]

  if (clamp) {
    S_ord <- pmin(pmax(S_ord, 0), 1)
  }

  if (enforce_monotone && ncol(S_ord) > 1L) {
    S_ord <- t(apply(S_ord, 1L, cummin))
  }

  inv_ord <- order(ord)
  out <- S_ord[, inv_ord, drop = FALSE]
  colnames(out) <- paste0("t=", times)
  as.data.frame(out)
}

#' Infer time points from survmat column names
#'
#' Internal helper that parses \code{"t=<time>"} column names when an explicit
#' time vector is not supplied.
#'
#' @param S A survival-probability matrix or data frame.
#'
#' @return Numeric vector of inferred times.
#' @keywords internal
.infer_survmat_times <- function(S) {
  cn <- colnames(S)
  if (is.null(cn) || !all(grepl("^t=", cn))) {
    stop("`times` must be supplied when `S` does not have `t=<time>` column names.")
  }

  times <- suppressWarnings(as.numeric(sub("^t=", "", cn)))
  if (anyNA(times)) {
    stop("Could not parse numeric times from `S` column names.")
  }

  times
}

#' Reconstruct Cox-style survival curves from linear predictors
#'
#' Internal helper shared by Cox-family learners to produce comparable survival
#' probabilities from a baseline cumulative hazard and subject-specific linear
#' predictors.
#'
#' @param lp_train Numeric vector of linear predictors on the training data.
#' @param lp_new Numeric vector of linear predictors on `newdata`.
#' @param y_train A `Surv` object for the training data.
#' @param times Numeric evaluation times.
#'
#' @return A standardized survival-probability data.frame.
#' @keywords internal
.predict_cox_from_lp <- function(lp_train, lp_new, y_train, times) {
  base_haz <- survival::basehaz(survival::coxph(y_train ~ lp_train), centered = FALSE)

  survmat <- vapply(times, function(t) {
    H0_t <- stats::approx(base_haz$time, base_haz$hazard, xout = t, rule = 2)$y
    exp(-H0_t * exp(lp_new))
  }, numeric(length(lp_new)))

  if (is.vector(survmat)) {
    survmat <- matrix(survmat, ncol = length(times))
  }

  .finalize_survmat(survmat, times = times)
}

#' Convert a survival-probability matrix (survmat) to cumulative hazard
#'
#' Uses the identity \eqn{\Lambda(t \mid x) = -\log S(t \mid x)} to compute
#' cumulative hazards for each observation and time point.
#'
#' @param S A `surv_mat`: numeric matrix/data.frame of survival probabilities.
#' @param eps Numeric in (0,1). Stabilizer to avoid \code{log(0)}. Values of
#'   \code{S} are clamped to \eqn{[\varepsilon, 1-\varepsilon]} before logging.
#'   Default is \code{1e-12}.
#'
#' @return A numeric matrix with the same dimensions as \code{S}, containing
#' cumulative hazards \eqn{\Lambda(t \mid x)}.
#'
#' @examples
#' S <- data.frame(`t=1` = c(0.9, 0.8), `t=2` = c(0.7, 0.6))
#' survmat_to_chf(S)
#'
#' @export
survmat_to_chf <- function(S, eps = 1e-12) {
  S <- .survmat_as_matrix(S, clamp = FALSE)

  if (!is.finite(eps) || eps <= 0 || eps >= 1) {
    stop("`eps` must be in (0,1).")
  }

  S2 <- pmax(pmin(S, 1 - eps), eps)
  -log(S2)
}

#' Convert a survival-probability matrix (survmat) to hazards on a time grid
#'
#' Computes *individual-specific* (i.e., conditional on covariates) hazards
#' from predicted survival curves \eqn{S(t \mid x_i)} on a discrete time grid.
#'
#' @details
#' This function returns a \emph{piecewise-constant hazard} per interval:
#' \deqn{
#' h_{i,j} \approx \frac{\Lambda_i(t_j) - \Lambda_i(t_{j-1})}{t_j - t_{j-1}}
#' = \frac{\log S_i(t_{j-1}) - \log S_i(t_j)}{t_j - t_{j-1}}.
#' }
#' with \eqn{S_i(0) = 1} and \eqn{\Lambda_i(t) = -\log S_i(t)}.
#'
#' @param S A `surv_mat`: numeric matrix/data.frame of survival probabilities.
#' @param times Numeric vector of time points corresponding to the columns of \code{S}.
#'   Must be strictly increasing and non-negative.
#' @param eps Numeric in (0,1). Stabilizer to avoid \code{log(0)} and division by
#'   zero; survival probabilities are clamped to \eqn{[\varepsilon, 1-\varepsilon]}.
#'   Default is \code{1e-12}.
#' @param t0 Numeric scalar; left boundary for the first interval. Default is \code{0}.
#'   If your grid already starts at 0 and you want the first interval to be (0, t1],
#'   keep \code{t0 = 0}.
#'
#' @return A numeric matrix of hazards with the same dimensions as \code{S}.
#' Each column corresponds to the interval ending at \code{times[j]}.
#'
#' @examples
#' times <- c(1, 2, 3)
#' S <- matrix(c(0.9, 0.8, 0.7,
#'               0.95,0.9,0.85), nrow = 2, byrow = TRUE)
#' colnames(S) <- paste0("t=", times)
#' H <- survmat_to_haz(S, times)
#' H
#'
#' @export
survmat_to_haz <- function(S, times, eps = 1e-12, t0 = 0) {
  S <- .survmat_as_matrix(S, times = times, strict_colnames = FALSE, clamp = FALSE)

  if (!is.finite(eps) || eps <= 0 || eps >= 1) stop("`eps` must be in (0,1).")
  if (!is.finite(t0) || t0 < 0) stop("`t0` must be finite and >= 0.")

  times <- as.numeric(times)
  if (times[1] <= t0) {
    stop("`times[1]` must be > t0 (default t0=0) to define the first interval.")
  }

  # Clamp S for stability, then compute CHF and interval hazard
  S2 <- pmax(pmin(S, 1 - eps), eps)
  logS <- log(S2)

  # dt_j = t_j - t_{j-1}, with t_0 = t0 and log S(t0) = log(1)=0
  dt <- c(times[1] - t0, diff(times))
  if (any(dt <= 0)) stop("`times` must be strictly increasing and > t0.")

  logS_prev <- cbind(rep(0, nrow(S2)), logS[, -ncol(logS), drop = FALSE])
  haz <- (logS_prev - logS) / matrix(dt, nrow = nrow(S2), ncol = length(dt), byrow = TRUE)

  # Hazards should be >= 0; numeric artifacts can create tiny negatives
  pmax(haz, 0)
}

#' Compute restricted mean survival time (RMST) from a survival-probability matrix (survmat)
#'
#' Computes \eqn{\mathrm{RMST}(\tau) = \int_0^\tau S(t)\, dt} for each observation,
#' approximated via the trapezoidal rule on the provided time grid.
#'
#' If the grid does not include \code{0}, the function prepends \eqn{(0, S(0)=1)}
#' to correctly integrate from time 0.
#'
#' @param S A `surv_mat`: numeric matrix/data.frame of survival probabilities.
#' @param times Numeric vector of time points corresponding to columns of \code{S}.
#'   Must be strictly increasing and non-negative.
#' @param tau Positive numeric scalar; upper integration limit. Default is \code{max(times)}.
#'
#' @return Numeric vector of RMST values (length = \code{nrow(S)}).
#'
#' @examples
#' times <- c(1, 2, 3)
#' S <- matrix(c(0.9, 0.8, 0.7,
#'               0.95,0.9,0.85), nrow = 2, byrow = TRUE)
#' colnames(S) <- paste0("t=", times)
#' survmat_to_rmst(S, times, tau = 3)
#'
#' @export
survmat_to_rmst <- function(S, times, tau = max(times)) {
  S <- .survmat_as_matrix(S, times = times, strict_colnames = FALSE, clamp = FALSE)
  times <- as.numeric(times)

  if (!is.finite(tau) || tau <= 0) stop("`tau` must be a positive finite number.")

  # Augment with (0,1) if needed
  if (times[1] > 0) {
    times2 <- c(0, times)
    S2 <- cbind(rep(1, nrow(S)), S)
  } else {
    times2 <- times
    S2 <- S
  }

  keep <- times2 <= tau
  if (sum(keep) < 2) stop("Need at least 2 time points <= tau to compute RMST.")
  t <- times2[keep]
  Sg <- S2[, keep, drop = FALSE]

  dt <- diff(t)
  left  <- Sg[, -ncol(Sg), drop = FALSE]
  right <- Sg[, -1, drop = FALSE]

  rmst <- rowSums(((left + right) / 2) * matrix(dt, nrow = nrow(Sg), ncol = length(dt), byrow = TRUE))
  as.numeric(rmst)
}


#' Compute a survival-time quantile from a survival-probability matrix (survmat) with grid-based approach
#'
#' Returns the smallest grid time \eqn{t} such that \eqn{S(t) \le 1-p}.
#' This is a right-continuous approximation on the provided time grid.
#'
#' For \code{p = 0.5}, this yields a grid-based **median survival time**.
#'
#' @param S A `surv_mat`: numeric matrix/data.frame of survival probabilities.
#' @param times Numeric vector of time points corresponding to columns of \code{S}.
#'   Must be strictly increasing and non-negative.
#' @param p Probability in (0,1). Default is \code{0.5}.
#'
#' @return Numeric vector of quantile times (length = \code{nrow(S)}). Returns
#' \code{NA} if the threshold is not crossed on the grid (e.g., survival stays
#' above \eqn{1-p} over all \code{times}).
#'
#' @examples
#' times <- c(1, 2, 3)
#' S <- matrix(c(0.9, 0.6, 0.4,
#'               0.95,0.9,0.85), nrow = 2, byrow = TRUE)
#' colnames(S) <- paste0("t=", times)
#' survmat_to_quantile(S, times, p = 0.5)  # median on the grid
#'
#' @export
survmat_to_quantile <- function(S, times, p = 0.5) {
  S <- .survmat_as_matrix(S, times = times, strict_colnames = FALSE, clamp = FALSE)
  times <- as.numeric(times)

  if (!is.finite(p) || p <= 0 || p >= 1) stop("`p` must be in (0,1).")

  thr <- 1 - p
  out <- rep(NA_real_, nrow(S))

  for (i in seq_len(nrow(S))) {
    j <- which(S[i, ] <= thr)[1L]
    out[i] <- if (is.na(j)) NA_real_ else times[j]
  }
  out
}

#' Plot Predicted Survival Curves from a survmat
#'
#' Plots one or more predicted survival curves from a survival-probability
#' matrix. By default, each row of \code{S} is shown as an individual curve. If
#' a \code{group} vector is supplied, curves are summarized by group using the
#' requested aggregation function.
#'
#' @param S A `surv_mat`: numeric matrix/data.frame of survival probabilities.
#' @param times Optional numeric vector of time points corresponding to columns
#'   of \code{S}. If omitted, times are inferred from column names of the form
#'   \code{"t=<time>"}.
#' @param group Optional vector of group labels of length \code{nrow(S)}. When
#'   supplied, grouped summary curves are plotted.
#' @param ids Optional integer vector of row ids to subset before plotting.
#' @param summary_fun Aggregation for grouped curves; one of \code{"mean"}
#'   (default) or \code{"median"}.
#' @param show_individual Logical. If \code{NULL} (default), individual curves
#'   are shown only when \code{group} is omitted. If \code{TRUE}, individual
#'   curves are overlaid beneath the grouped summary curves.
#' @param alpha Alpha transparency for individual curves (default \code{0.2}).
#' @param linewidth Line width for plotted curves (default \code{0.7}).
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' S <- data.frame(`t=1` = c(0.95, 0.90, 0.92),
#'                 `t=2` = c(0.80, 0.70, 0.78),
#'                 `t=3` = c(0.60, 0.45, 0.55))
#' plot_survmat(S)
#' plot_survmat(S, group = c("A", "B", "A"))
#'
#' @export
plot_survmat <- function(S,
                         times = NULL,
                         group = NULL,
                         ids = NULL,
                         summary_fun = c("mean", "median"),
                         show_individual = NULL,
                         alpha = 0.2,
                         linewidth = 0.7) {
  summary_fun <- match.arg(summary_fun)

  if (is.null(times)) {
    times <- .infer_survmat_times(S)
  }

  S <- .survmat_as_matrix(S, times = times, strict_colnames = FALSE, clamp = TRUE)
  plot_df <- as.data.frame(S)
  plot_df$.id <- seq_len(nrow(plot_df))

  if (!is.null(ids)) {
    plot_df <- plot_df[plot_df$.id %in% ids, , drop = FALSE]
  }

  if (!nrow(plot_df)) {
    stop("No rows remain to plot.")
  }

  if (is.null(show_individual)) {
    show_individual <- is.null(group)
  }

  if (!is.null(group)) {
    if (length(group) != nrow(S)) {
      stop("`group` must have length equal to nrow(S).")
    }
    plot_df$group <- factor(group[plot_df$.id])
  }

  plot_long <- tidyr::pivot_longer(
    plot_df,
    cols = grep("^t=", names(plot_df), value = TRUE),
    names_to = "time",
    values_to = "surv_prob"
  )
  plot_long$time <- as.numeric(sub("^t=", "", plot_long$time))

  if (is.null(group)) {
    if (isTRUE(show_individual)) {
      return(
        ggplot2::ggplot(
          plot_long,
          ggplot2::aes(x = time, y = surv_prob, group = .id)
        ) +
          ggplot2::geom_line(alpha = alpha, linewidth = linewidth, color = "steelblue") +
          ggplot2::labs(
            title = "Predicted survival curves",
            x = "Time",
            y = "Survival probability"
          ) +
          ggplot2::theme_minimal()
      )
    }

    summary_df <- plot_long |>
      dplyr::group_by(time) |>
      dplyr::summarise(
        surv_prob = if (summary_fun == "mean") mean(surv_prob) else stats::median(surv_prob),
        .groups = "drop"
      )

    return(
      ggplot2::ggplot(summary_df, ggplot2::aes(x = time, y = surv_prob)) +
        ggplot2::geom_line(linewidth = linewidth + 0.2, color = "steelblue") +
        ggplot2::labs(
          title = paste("Predicted survival curve (", summary_fun, " summary)", sep = ""),
          x = "Time",
          y = "Survival probability"
        ) +
        ggplot2::theme_minimal()
    )
  }

  summary_df <- plot_long |>
    dplyr::group_by(group, time) |>
    dplyr::summarise(
      surv_prob = if (summary_fun == "mean") mean(surv_prob) else stats::median(surv_prob),
      .groups = "drop"
    )

  p <- ggplot2::ggplot(summary_df, ggplot2::aes(x = time, y = surv_prob, color = group)) +
    ggplot2::geom_line(linewidth = linewidth + 0.2) +
    ggplot2::labs(
      title = paste("Predicted survival curves by group (", summary_fun, " summary)", sep = ""),
      x = "Time",
      y = "Survival probability",
      color = "Group"
    ) +
    ggplot2::theme_minimal()

  if (isTRUE(show_individual)) {
    p <- p + ggplot2::geom_line(
      data = plot_long,
      ggplot2::aes(group = .id),
      alpha = alpha,
      linewidth = linewidth
    )
  }

  p
}
