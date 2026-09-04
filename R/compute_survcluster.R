#' Cluster Predicted Survival Curves
#'
#' Clusters the individualized survival curves predicted by a fitted
#' `mlsurv_model` using [unsurv::unsurv()] (Partitioning Around Medoids on a
#' weighted feature representation of the curves).
#'
#' @param model An `mlsurv_model` fitted via a `fit_*()` function; must include a
#'   valid `learner` so the corresponding `predict_<learner>()` can be dispatched.
#' @param newdata A data frame with predictors for the subjects to cluster.
#' @param times Numeric vector of evaluation times (the common time grid on
#'   which survival curves are predicted and clustered).
#' @param K Optional integer number of clusters. If `NULL`, selected by
#'   silhouette width (see [unsurv::unsurv()]).
#' @param K_max Maximum `K` considered when `K` is `NULL`.
#' @param distance Distance type passed to [unsurv::unsurv()]: `"L2"` or `"L1"`.
#' @param seed Optional integer seed forwarded to [unsurv::unsurv()].
#' @param ... Additional arguments forwarded to [unsurv::unsurv()] (e.g.
#'   `weights`, `enforce_monotone`, `smooth_median_width`, `standardize_cols`).
#'
#' @details
#' Predicted survival is obtained from the model's `predict_<learner>()` on
#' `times`, coerced to a numeric matrix, and passed to [unsurv::unsurv()].
#' Requires the \pkg{unsurv} package.
#'
#' @return A list with components:
#' \describe{
#'   \item{fit}{The fitted `unsurv` object (clusters, medoids, silhouette, ...).}
#'   \item{survmat}{The numeric survival-probability matrix that was clustered.}
#'   \item{times}{The evaluation time grid used.}
#'   \item{learner}{The learner label (from `model$learner`).}
#' }
#'
#' @references
#' El Badisy, I. (2026). unsurv: Clustering Individualized Survival Curves.
#' Bioinformatics Advances, vbag218. \doi{10.1093/bioadv/vbag218}
#'
#' @examples
#' if (requireNamespace("unsurv", quietly = TRUE) &&
#'     requireNamespace("cluster", quietly = TRUE)) {
#'   mod <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = veteran)
#'   time_grid <- seq(10, 300, length.out = 20)
#'   cl <- compute_survcluster(mod, newdata = veteran, times = time_grid, K = 2, seed = 1)
#'   table(cl$fit$clusters)
#' }
#'
#' @seealso [plot_survcluster()]
#' @export
compute_survcluster <- function(model, newdata, times, K = NULL, K_max = 10,
                                 distance = c("L2", "L1"), seed = NULL, ...) {
  if (!requireNamespace("unsurv", quietly = TRUE)) {
    stop("Package 'unsurv' is required for compute_survcluster(). Please install it.", call. = FALSE)
  }

  distance <- match.arg(distance)

  if (is.null(model$learner) || !exists(paste0("predict_", model$learner))) {
    stop("Could not infer prediction function. Ensure model was created using fit_*() and includes a valid 'learner'.")
  }
  predict_function <- get(paste0("predict_", model$learner))

  pred_raw <- predict_function(model, newdata = newdata, times = times)
  survmat <- .survmat_as_matrix(pred_raw, times = times)

  fit <- unsurv::unsurv(
    S = survmat, times = times,
    K = K, K_max = K_max,
    distance = distance,
    seed = seed,
    ...
  )

  list(
    fit = fit,
    survmat = survmat,
    times = as.numeric(times),
    learner = model$learner
  )
}

#' Plot Clustered Survival Curves
#'
#' Plots the individualized survival curves from [compute_survcluster()],
#' colored by cluster membership, with medoid curves overlaid.
#'
#' @param survcluster_output Output list returned by [compute_survcluster()].
#' @param title Plot title. If missing, an automatically generated title is
#'   used. Pass `NULL` to omit the title entirely (e.g., for journals
#'   requiring caption-only figures).
#'
#' @return A `ggplot2` object showing per-subject survival curves colored by
#' cluster, with cluster medoid curves overlaid.
#'
#' @examples
#' if (requireNamespace("unsurv", quietly = TRUE) &&
#'     requireNamespace("cluster", quietly = TRUE)) {
#'   mod <- fit_coxph(Surv(time, status) ~ age + karno + celltype, data = veteran)
#'   time_grid <- seq(10, 300, length.out = 20)
#'   cl <- compute_survcluster(mod, newdata = veteran, times = time_grid, K = 2, seed = 1)
#'   plot_survcluster(cl)
#' }
#'
#' @seealso [compute_survcluster()]
#' @export
plot_survcluster <- function(survcluster_output, title) {
  fit <- survcluster_output$fit
  survmat <- survcluster_output$survmat
  times <- survcluster_output$times
  clusters <- factor(fit$clusters)

  if (missing(title)) {
    title <- paste0("Survival Curve Clusters (K = ", fit$K, ") | ", survcluster_output$learner)
  }

  curve_df <- data.frame(
    subject = rep(seq_len(nrow(survmat)), each = length(times)),
    time = rep(times, times = nrow(survmat)),
    survival = as.vector(t(survmat)),
    cluster = rep(clusters, each = length(times))
  )

  medoids <- fit$medoids
  medoid_df <- data.frame(
    time = rep(times, times = nrow(medoids)),
    survival = as.vector(t(medoids)),
    cluster = factor(sort(unique(fit$clusters)), levels = levels(clusters))[
      rep(seq_len(nrow(medoids)), each = length(times))
    ]
  )

  p <- ggplot(curve_df, aes(x = time, y = survival, group = subject, color = cluster)) +
    geom_line(linewidth = 0.25, alpha = 0.20) +
    geom_line(
      data = medoid_df,
      aes(x = time, y = survival, group = cluster, color = cluster),
      linewidth = 1.1
    ) +
    scale_color_survalis() +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Time", y = "Predicted Survival Probability", color = "Cluster") +
    theme_survalis()

  if (!is.null(title)) p <- p + labs(title = title)

  p
}
