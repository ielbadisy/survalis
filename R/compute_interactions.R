#' Compute Feature Interactions for Survival Predictions
#'
#' Estimates global and time-varying interaction strengths of model predictions,
#' using a Friedman-H style decomposition adapted to survival partial dependence.
#' Works with any `mlsurv_model` that has a matching `predict_*()` method returning
#' survival probabilities.
#'
#' @param model An `mlsurv_model` created via a `fit_*()` function; must include a
#'   valid `learner` so the appropriate `predict_<learner>()` can be dispatched.
#' @param data Data frame used to probe the model (typically training data).
#' @param times Numeric vector of evaluation times used for prediction.
#' @param target_time Single time at which to quantify interactions for
#'   \code{type = "1way"} or \code{type = "heatmap"}. Required for these types.
#' @param features Optional character vector of feature names to evaluate.
#'   Defaults to all predictors in \code{data} except the outcome columns.
#' @param type One of:
#'   \itemize{
#'     \item \code{"1way"} -- per-feature Friedman-H interaction strength vs. all others at \code{target_time};
#'     \item \code{"heatmap"} -- pairwise feature interaction matrix at \code{target_time};
#'     \item \code{"time"} -- per-feature interaction strength across \code{times}.
#'   }
#' @param grid.size Integer; number of observations sampled from \code{data}
#'   at which the partial dependence functions are evaluated, and over which
#'   they are averaged (default 30). The cost grows with its square.
#' @param batch.size Reserved for future batching support (currently unused).
#'
#' @details
#' For a target time \eqn{t^{*}}, let \eqn{f(x)} be the predicted survival
#' probability, and let \eqn{PD_S} be the partial dependence function of the
#' features in \eqn{S}, each centred to mean zero over the sampled observations
#' \eqn{x^{(1)}, \dots, x^{(m)}} (\eqn{m} = \code{grid.size}), which also
#' serve as the marginalizing sample. Following Friedman and Popescu (2008),
#' the one-way statistic of feature \eqn{j} is
#' \deqn{ H_j^2 = \frac{\sum_i \{ f(x^{(i)}) - PD_j(x^{(i)}_j) - PD_{-j}(x^{(i)}_{-j}) \}^2}{\sum_i f(x^{(i)})^2}, }
#' with \eqn{f} centred too, and the pairwise statistic is
#' \deqn{ H_{jk}^2 = \frac{\sum_i \{ PD_{jk}(x^{(i)}_j, x^{(i)}_k) - PD_j(x^{(i)}_j) - PD_k(x^{(i)}_k) \}^2}{\sum_i PD_{jk}(x^{(i)}_j, x^{(i)}_k)^2}. }
#' Both are 0 when the model is additive in the features involved. The
#' reported value is \eqn{H} (not clipped), and is \code{NA} when the
#' prediction does not vary.
#'
#' Larger values indicate stronger non-additivity (interaction) involving the feature(s).
#' The \code{"time"} mode repeats the computation across all \code{times} to show dynamics.
#'
#' @return
#' A data frame whose structure depends on \code{type}:
#' \itemize{
#'   \item \code{"1way"}: columns \code{feature}, \code{interaction}.
#'   \item \code{"heatmap"}: columns \code{feature1}, \code{feature2}, \code{interaction}
#'         (symmetric with zeros on the diagonal).
#'   \item \code{"time"}: columns \code{feature}, \code{time}, \code{interaction}.
#' }
#'
#' @references
#' Friedman, J. H., and Popescu, B. E. (2008). Predictive learning via rule ensembles.
#' \emph{Annals of Applied Statistics}. (Friedman's \eqn{H} interaction measure.)
#'
#' @examples
#' mod <- fit_coxph(Surv(time, status) ~ age + karno + trt, data = veteran)
#' times <- c(80, 160)
#' compute_interactions(
#'   model = mod,
#'   data = veteran,
#'   times = times,
#'   target_time = 80,
#'   features = c("age", "karno"),
#'   type = "1way",
#'   grid.size = 6
#' )
#'
#' @noRd

compute_interactions <- function(model, data, times,
                                 target_time = NULL,
                                 features = NULL,
                                 type = c("1way", "heatmap", "time"),
                                 grid.size = 30,
                                 batch.size = 100) {
  type <- match.arg(type)

  if (type != "time" && is.null(target_time)) stop("target_time must be specified for type = '1way' or 'heatmap'")
  if (is.null(features)) {
    features <- .model_features(model, colnames(data))
  }

  predict_function <- .predict_fun(model)
  eval_times <- if (type == "time") times else target_time
  m <- min(as.integer(grid.size), nrow(data))
  rows <- data[sample.int(nrow(data), m), , drop = FALSE]
  pred <- function(nd) {
    as.matrix(.finalize_survmat(predict_function(model, nd, eval_times), eval_times))
  }
  ctr <- function(M) sweep(M, 2L, colMeans(M))
  # Partial dependence of the features in S at each sampled row (m x times),
  # averaging over the sampled rows for every other feature, centred.
  ii <- rep(seq_len(m), times = m)   # row supplying the S values
  kk <- rep(seq_len(m), each = m)    # row supplying everything else
  pd <- function(S) {
    big <- rows[kk, , drop = FALSE]
    for (s in S) big[[s]] <- rows[[s]][ii]
    ctr(rowsum(pred(big), ii) / m)
  }
  hstat <- function(num, den) {
    out <- sqrt(colSums(num^2) / colSums(den^2))
    out[!is.finite(out)] <- NA_real_
    out
  }
  full <- ctr(pred(rows))

  if (type %in% c("1way", "time")) {
    H <- vapply(features, function(f) {
      hstat(full - pd(f) - pd(setdiff(features, f)), full)
    }, numeric(length(eval_times)))
    H <- matrix(H, nrow = length(eval_times))
    if (type == "1way") {
      return(data.frame(feature = features, interaction = H[1, ]))
    }
    return(data.frame(feature = rep(features, each = length(eval_times)),
                      time = rep(eval_times, length(features)),
                      interaction = as.vector(H)))
  }

  pd1 <- lapply(stats::setNames(features, features), pd)
  pairs <- t(utils::combn(features, 2))
  h <- vapply(seq_len(nrow(pairs)), function(i) {
    f1 <- pairs[i, 1]; f2 <- pairs[i, 2]
    pd12 <- pd(c(f1, f2))
    hstat(pd12 - pd1[[f1]] - pd1[[f2]], pd12)[1]
  }, numeric(1))
  df <- data.frame(feature1 = pairs[, 1], feature2 = pairs[, 2], interaction = h)
  basetable::rbindfill(list(
    df,
    data.frame(feature1 = df$feature2, feature2 = df$feature1, interaction = df$interaction),
    data.frame(feature1 = features, feature2 = features, interaction = 0)
  ), fill = TRUE)
}



plot_interactions <- function(object, type = c("1way", "heatmap", "time"), title) {
  type <- match.arg(type)
  if (type == "1way") {
  if (missing(title)) title <- NULL
  p <- ggplot(object, ggplot2::aes(x = interaction, y = reorder(feature, interaction))) +
      geom_col(fill = .survalis_palette[1]) +
      theme_survalis() +
      labs(x = "Interaction strength (H)", y = "Feature")
  if (!is.null(title)) p <- p + labs(title = title)
  p
  } else if (type == "heatmap") {
    if (missing(title)) title <- "Pairwise Interaction Heatmap"
    p <- ggplot(object, ggplot2::aes(x = feature1, y = feature2, fill = interaction)) +
      geom_tile(color = "white", linewidth = 0.6) +
      ggplot2::scale_fill_viridis_c(
        direction = -1,
        na.value = "grey85",
        guide = ggplot2::guide_colorbar(barwidth = ggplot2::unit(5, "cm"), barheight = ggplot2::unit(0.4, "cm"))
      ) +
      theme_survalis()
    if (!is.null(title)) p <- p + labs(title = title)
    p
  } else {
    if (missing(title)) title <- "Time-Varying Feature Interactions"
    p <- ggplot(object, ggplot2::aes(x = time, y = interaction, color = feature)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 2) +
      scale_color_survalis() +
      theme_survalis() +
      labs(x = "Time", y = "Interaction Strength")
    if (!is.null(title)) p <- p + labs(title = title)
    p
  }
}
