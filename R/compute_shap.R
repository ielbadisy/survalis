#' Compute local SHAP-like contributions for survival predictions
#'
#' Estimates per-feature contributions (à la SHAP) to the predicted survival
#' probability for a **single** observation at one or more time points. For each
#' time, the method samples random feature orders, marginalizes future features
#' with values from a baseline dataset, and accumulates the marginal effects of
#' adding each feature back.
#'
#' @param model A fitted survival model produced by \code{fit_*()} that exposes
#'   \code{model$learner} (e.g., "coxph") and \code{model$data} (the training frame).
#' @param newdata A data frame with **exactly one row** (the instance to explain).
#' @param baseline_data A data frame to sample background values from (typically
#'   the training data used to fit \code{model}).
#' @param times Numeric vector of evaluation times (same scale as the outcome).
#' @param sample.size Integer, number of random feature orderings to sample per
#'   time (default \code{100}).
#' @param features Optional character vector of covariate names to explain. If
#'   \code{NULL} (default), the model's predictors (the right-hand side of its
#'   formula) that are columns of \code{newdata}.
#' @param aggregate Logical; if \code{TRUE}, aggregate contributions across
#'   \code{times} into a single value per feature (see \code{method}).
#' @param method Character; aggregation method if \code{aggregate=TRUE}:
#'   \code{"meanabs"} (mean absolute contribution) or \code{"integral"}
#'   (trapezoidal integral over time). Default \code{"meanabs"}.
#'
#' @return
#' If \code{aggregate = FALSE}: a data frame with columns \code{feature}, \code{phi},
#' and \code{time} (one row per feature per time). If \code{aggregate = TRUE}:
#' a data frame with columns \code{feature} and \code{phi} (one row per feature),
#' with attribute \code{"shap_method"} set to the aggregation used.
#'
#' @details
#' The prediction function is inferred from \code{model$learner} as
#' \code{predict_<learner>} and called with signature
#' \code{predict_fun(model, newdata, times = times)}. Factor levels in
#' \code{newdata} are harmonized to those in \code{model$data}.
#'
#' @seealso `plot_shap()`, the various \code{predict_*} methods (e.g. `predict_coxph()`)
#'
#' @examples
#' \donttest{
#' mod <- fit_coxph(survival::Surv(time, status) ~ age + karno + celltype, data = veteran)
#'
#' shap_td <- compute_shap(
#'   model         = mod,
#'   newdata       = veteran[100, , drop = FALSE],
#'   baseline_data = veteran,
#'   times         = c(100, 200),
#'   sample.size   = 5,
#'   aggregate     = FALSE
#' )
#' head(shap_td)
#'
#' shap_meanabs <- compute_shap(
#'   model         = mod,
#'   newdata       = veteran[100, , drop = FALSE],
#'   baseline_data = veteran,
#'   times         = c(100, 200),
#'   sample.size   = 5,
#'   aggregate     = TRUE,
#'   method        = "meanabs"
#' )
#' head(shap_meanabs)
#' }
#'
#' @noRd

compute_shap <- function(model, newdata, baseline_data,
  times, sample.size = 100,
  features = NULL,
  aggregate = FALSE,
  method = c("meanabs", "integral")) {
stopifnot(nrow(newdata) == 1)
method <- match.arg(method)

if (is.null(model$learner) && !is.function(model$predict_fun)) {
stop("`model` must come from fit().", call. = FALSE)
}



predict_function <- .predict_fun(model)
feature_names <- if (is.null(features)) {
  .model_features(model, colnames(newdata))
} else {
  features
}
n_features <- length(feature_names)
times <- sort(times)

# Permutation estimator (Strumbelj and Kononenko 2014). For each draw, take a
# random feature order and a random baseline row z, and walk from z to the
# explained row x by switching the features to their x value one at a time in
# that order. Each switch is credited to the feature switched, so in every draw
# the contributions sum to f(x) - f(z). All rows are predicted in one call.
p1 <- n_features + 1L
zi <- sample.int(nrow(baseline_data), sample.size, replace = TRUE)
walk <- baseline_data[rep(zi, each = p1), , drop = FALSE]
step <- rep(seq_len(p1) - 1L, sample.size)          # 0 = z, ..., p = x
pos <- t(replicate(sample.size, sample.int(n_features)))  # position of each feature
if (n_features == 1L) pos <- matrix(pos, ncol = 1L)
colnames(pos) <- feature_names
draw <- rep(seq_len(sample.size), each = p1)
for (j in feature_names) {
  on <- step >= pos[draw, j]
  walk[[j]][on] <- newdata[[j]]
}
P <- as.matrix(.finalize_survmat(predict_function(model, walk, times), times))
gain <- P[-seq(1L, nrow(P), by = p1), , drop = FALSE] -
  P[-seq(p1, nrow(P), by = p1), , drop = FALSE]     # one row per (draw, step)
credited <- as.vector(apply(pos, 1L, function(r) feature_names[order(r)]))
phi <- rowsum(gain, factor(credited, levels = feature_names)) / sample.size

all_results <- lapply(seq_along(times), function(k) {
  data.frame(feature = feature_names, phi = as.numeric(phi[, k]), time = times[k])
})

final_result <- basetable::rbindfill(unname(all_results), fill = TRUE)

if (aggregate) {
phi_long <- final_result[c("feature", "time", "phi")]
phi_wide <- basetable::towide(
  phi_long, names = "time", values = "phi", idcols = "feature", fun = identity
)
time_cols <- as.character(sort(unique(phi_long$time)))
phi_wide <- phi_wide[, c("feature", time_cols)]

phi_matrix <- as.matrix(phi_wide[, -1, drop = FALSE])
rownames(phi_matrix) <- phi_wide$feature

agg_phi <- switch(
method,
meanabs = apply(phi_matrix, 1, function(row) mean(abs(row))),
integral = apply(phi_matrix, 1, function(row) trapz(times, row))
)

result_agg <- data.frame(
feature = rownames(phi_matrix),
phi = agg_phi
)
attr(result_agg, "shap_method") <- method
return(result_agg)
}

attr(final_result, "shap_method") <- method
return(final_result)
}






#' Plot SHAP-like contributions for survival models
#'
#' Plots time-dependent SHAP estimates (lines over time) or aggregated SHAP
#' (bar chart) returned by `compute_shap()`.
#'
#' @param shapley_result A data frame returned by `compute_shap()`. If it has a
#'   \code{time} column, a time plot is produced; otherwise an aggregated bar plot.
#' @param type One of \code{"auto"} (default), \code{"time"}, or \code{"integrated"}.
#'   With \code{"auto"}, the presence of a \code{time} column determines the plot type.
#' @param title Plot title. If missing, an automatically generated title is
#'   used. Pass \code{NULL} to omit the title entirely (e.g., for journals
#'   requiring caption-only figures).
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' \donttest{
#' mod <- fit_coxph(survival::Surv(time, status) ~ age + karno + celltype, data = veteran)
#' shap_td <- compute_shap(mod, veteran[10, , drop = FALSE],
#' veteran, times = c(50, 100), sample.size = 5)
#' p1 <- plot_shap(shap_td)      # auto -> time plot
#'
#' shap_ag <- compute_shap(mod, veteran[10, , drop = FALSE],
#' veteran, times = c(50, 100),
#' sample.size = 5, aggregate = TRUE, method = "meanabs")
#' p2 <- plot_shap(shap_ag)      # auto -> bar plot
#' }
#' @noRd

plot_shap <- function(shapley_result, type = c("auto"), title) {
type <- match.arg(type)

is_time_dependent <- "time" %in% colnames(shapley_result)

if (type == "auto") {
type <- if (is_time_dependent) "time" else "integrated"
}

if (type == "time") {
if (missing(title)) title <- "Shapley Contributions Over Time"
p <- ggplot(shapley_result, ggplot2::aes(x = time, y = phi, color = feature)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_survalis() +
  labs(
    x = "Time",
    y = "SHAP value"
  ) +
  theme_survalis()
if (!is.null(title)) p <- p + labs(title = title)
p
} else {
shap_method <- attr(shapley_result, "shap_method", exact = TRUE)
shap_method <- if (is.null(shap_method)) "integral" else shap_method

if (missing(title)) {
  title <- switch(
    shap_method,
    meanabs = "Mean Absolute SHAP",
    integral = "Integrated SHAP",
    "SHAP Summary"
  )
}

shapley_result <- basetable::orderrows(shapley_result, by = "phi")
shapley_result$feature <- factor(shapley_result$feature, levels = shapley_result$feature)
shapley_result$direction <- ifelse(shapley_result$phi >= 0, "Positive", "Negative")

p <- ggplot(shapley_result, ggplot2::aes(x = phi, y = feature, fill = direction)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = c("Positive" = .survalis_palette[1], "Negative" = .survalis_palette[2])) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    x = "SHAP value",
    y = "Feature"
  ) +
  theme_survalis()
if (!is.null(title)) p <- p + labs(title = title)
p
}
}

