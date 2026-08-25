#' Time-Integrated SHAP Values for a Set of Observations
#'
#' Computes \code{\link{compute_shap}} for every row of \code{newdata} and,
#' instead of averaging across observations (as the internal
#' \code{compute_shap_mean()} helper does), keeps one time-integrated SHAP
#' value per (observation, feature) pair together with each observation's raw
#' feature value -- the input needed for \code{\link{plot_shap_beeswarm}}.
#'
#' @param model An \code{mlsurv_model} object from a \code{fit_*()} learner.
#' @param newdata A data.frame of observations to explain (one row per
#'   observation).
#' @param baseline_data A data.frame to sample background values from
#'   (typically the training data used to fit \code{model}).
#' @param times Numeric vector of evaluation times passed to
#'   \code{\link{compute_shap}}.
#' @param sample.size Integer, number of random feature orderings sampled per
#'   observation and time point (passed to \code{\link{compute_shap}}).
#'   Default \code{100}.
#' @param features Optional character vector of covariate names to explain,
#'   passed through to \code{\link{compute_shap}}. If \code{NULL} (default),
#'   inferred as \code{setdiff(colnames(newdata), c(model$time, model$status))}.
#'
#' @details
#' For each observation \eqn{x}, per-time Shapley values
#' \eqn{\varphi_j(x, t)} (from \code{\link{compute_shap}}) are time-integrated
#' via the trapezoidal rule over the requested time grid,
#' \deqn{
#'   \Phi_j(x) = \int_0^{T} \varphi_j(x, t)\, dt
#'   \approx \sum_{k=2}^{m} \frac{\varphi_j(x, t_k) + \varphi_j(x, t_{k-1})}{2}
#'   (t_k - t_{k-1}),
#' }
#' using the same integration already implemented by
#' \code{compute_shap(..., aggregate = TRUE, method = "integral")}. This
#' produces one summary contribution per (observation, feature) that reflects
#' impact over the *entire* prediction horizon in a single value, which is
#' what \code{\link{plot_shap_beeswarm}} visualizes -- as opposed to a
#' per-timepoint SurvSHAP(t) panel.
#'
#' @return A long data.frame with columns \code{observation}, \code{feature},
#'   \code{phi} (time-integrated SHAP value), and \code{raw_value} (the
#'   observation's value for that feature).
#'
#' @seealso [compute_shap()], [plot_shap_beeswarm()]
#'
#' @examples
#' \donttest{
#' mod <- fit_coxph(Surv(time, status) ~ age + karno + diagtime, data = veteran)
#' sm <- compute_shap_matrix(mod, newdata = veteran[1:15, ],
#'                           baseline_data = veteran, times = c(50, 100, 200),
#'                           sample.size = 5)
#' head(sm)
#' }
#' @export
compute_shap_matrix <- function(model, newdata, baseline_data, times, sample.size = 100, features = NULL) {
  stopifnot(nrow(newdata) >= 1)
  feature_names <- if (is.null(features)) {
    setdiff(colnames(newdata), c(model$time, model$status))
  } else {
    features
  }

  rows <- lapply(seq_len(nrow(newdata)), function(i) {
    phi_df <- compute_shap(
      model = model,
      newdata = newdata[i, , drop = FALSE],
      baseline_data = baseline_data,
      times = times,
      sample.size = sample.size,
      features = feature_names,
      aggregate = TRUE,
      method = "integral"
    )
    data.frame(
      observation = i,
      feature = phi_df$feature,
      phi = phi_df$phi,
      raw_value = suppressWarnings(as.numeric(newdata[i, phi_df$feature])),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' SHAP Beeswarm Plot (Time-Integrated)
#'
#' Plots the distribution of time-integrated SHAP contributions
#' (\code{\link{compute_shap_matrix}} output) across observations, one row per
#' feature, with points coloured by each observation's (min-max scaled) raw
#' feature value -- summarizing global feature impact across the entire
#' prediction horizon in a single panel.
#'
#' @param shap_matrix A data.frame returned by \code{\link{compute_shap_matrix}}.
#' @param top_n Optional integer; if supplied, only the \code{top_n} features
#'   with the largest mean absolute SHAP value are shown.
#' @param title Plot title. If missing, an automatically generated title is
#'   used. Pass \code{NULL} to omit the title entirely.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @seealso [compute_shap_matrix()]
#'
#' @examples
#' \donttest{
#' mod <- fit_coxph(Surv(time, status) ~ age + karno + diagtime, data = veteran)
#' sm <- compute_shap_matrix(mod, newdata = veteran[1:15, ],
#'                           baseline_data = veteran, times = c(50, 100, 200),
#'                           sample.size = 5)
#' plot_shap_beeswarm(sm)
#' }
#' @export
plot_shap_beeswarm <- function(shap_matrix, top_n = NULL, title) {
  features <- unique(shap_matrix$feature)
  imp <- vapply(features, function(f) mean(abs(shap_matrix$phi[shap_matrix$feature == f])), numeric(1))
  ord <- features[order(imp)]

  if (!is.null(top_n)) {
    ord <- utils::tail(ord, top_n)
  }

  plot_df <- shap_matrix[shap_matrix$feature %in% ord, , drop = FALSE]
  plot_df$scaled_value <- stats::ave(plot_df$raw_value, plot_df$feature, FUN = function(x) {
    if (all(is.na(x)) || diff(range(x, na.rm = TRUE)) == 0) return(rep(0.5, length(x)))
    (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE))
  })
  plot_df$feature <- factor(plot_df$feature, levels = ord)

  if (missing(title)) title <- "SHAP beeswarm (time-integrated over prediction horizon)"

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = phi, y = feature, colour = scaled_value)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
    ggplot2::geom_jitter(height = 0.3, width = 0, alpha = 0.85, size = 1.8) +
    ggplot2::scale_colour_gradient(
      low = .survalis_palette[3], high = .survalis_palette[2], na.value = "grey60",
      breaks = c(0, 1), labels = c("Low", "High"), name = "Feature value"
    ) +
    ggplot2::labs(x = expression(Phi ~ "(time-integrated SHAP value)"), y = NULL) +
    theme_survalis() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(hjust = 0))
  if (!is.null(title)) p <- p + ggplot2::labs(title = title)
  p
}
