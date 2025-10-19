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
#' @seealso [plot_shap()], the various \code{predict_*} methods (e.g. [predict_coxph()])
#'
#' @examples
#' \donttest{
#' mod <- fit_coxph(survival::Surv(time, status) ~ age + karno + celltype, data = veteran)
#'
#' # Explain one patient at a few times
#' shap_td <- compute_shap(
#'   model         = mod,
#'   newdata       = veteran[100, , drop = FALSE],
#'   baseline_data = veteran,
#'   times         = c(100, 200, 300),
#'   sample.size   = 50,         # small for example speed
#'   aggregate     = FALSE
#' )
#' head(shap_td)
#'
#' # Aggregate across time (mean absolute contribution)
#' shap_meanabs <- compute_shap(
#'   model         = mod,
#'   newdata       = veteran[100, , drop = FALSE],
#'   baseline_data = veteran,
#'   times         = c(100, 200, 300),
#'   sample.size   = 50,
#'   aggregate     = TRUE,
#'   method        = "meanabs"
#' )
#' head(shap_meanabs)
#' }
#'
#' @export

compute_shap <- function(model, newdata, baseline_data,
  times, sample.size = 100,
  aggregate = FALSE,
  method = c("meanabs", "integral")) {
stopifnot(nrow(newdata) == 1)
method <- match.arg(method)

if (is.null(model$learner) || !exists(paste0("predict_", model$learner))) {
stop("Could not infer prediction function. Ensure model has a valid `learner` and was created using `fit_*()`.")
}



predict_function <- get(paste0("predict_", model$learner))
feature_names <- setdiff(colnames(newdata), c("time", "status", "event"))
n_features <- length(feature_names)
all_results <- list()

for (current_time in times) {
shapley_values <- replicate(sample.size, {
order_features <- sample(feature_names)
baseline_sample <- baseline_data[sample(nrow(baseline_data), 1), , drop = FALSE]
instance_with <- newdata
instance_without <- newdata
contributions <- numeric(n_features)

for (k in seq_len(n_features)) {
feature_k <- order_features[k]
if (k < n_features) {
future_features <- order_features[(k + 1):n_features]
instance_with[future_features] <- baseline_sample[future_features]
}

instance_without <- instance_with
instance_without[feature_k] <- baseline_sample[[feature_k]]

pred_with <- predict_function(model, instance_with, times)
pred_without <- predict_function(model, instance_without, times)

idx_time <- which.min(abs(times - current_time))
pred_with <- if (is.vector(pred_with)) pred_with[idx_time] else pred_with[, idx_time]
pred_without <- if (is.vector(pred_without)) pred_without[idx_time] else pred_without[, idx_time]

idx_feature <- which(feature_names == feature_k)
contributions[idx_feature] <- contributions[idx_feature] + (pred_with - pred_without)
}

contributions
}, simplify = TRUE)

shapley_df <- as.data.frame(t(shapley_values))
colnames(shapley_df) <- feature_names
mean_phi <- colMeans(shapley_df)

result_df <- data.frame(
feature = feature_names,
phi = as.numeric(mean_phi),
time = current_time
)
all_results[[as.character(current_time)]] <- result_df
}

final_result <- do.call(rbind, all_results)

if (aggregate) {
phi_wide <- final_result |>
select(feature, time, phi) |>
pivot_wider(names_from = time, values_from = phi)

phi_matrix <- as.matrix(phi_wide[, -1])
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
#' (bar chart) returned by [compute_shap()].
#'
#' @param shapley_result A data frame returned by [compute_shap()]. If it has a
#'   \code{time} column, a time plot is produced; otherwise an aggregated bar plot.
#' @param type One of \code{"auto"} (default), \code{"time"}, or \code{"integrated"}.
#'   With \code{"auto"}, the presence of a \code{time} column determines the plot type.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' \donttest{
#' # Time-dependent SHAP
#' mod <- fit_coxph(survival::Surv(time, status) ~ age + karno + celltype, data = veteran)
#' shap_td <- compute_shap(mod, veteran[10, , drop = FALSE],
#' veteran, times = c(50, 100, 150), sample.size = 20)
#' p1 <- plot_shap(shap_td)      # auto -> time plot
#'
#' # Aggregated SHAP
#' shap_ag <- compute_shap(mod, veteran[10, , drop = FALSE],
#' veteran, times = c(50, 100, 150),
#' sample.size = 20, aggregate = TRUE, method = "meanabs")
#' p2 <- plot_shap(shap_ag)      # auto -> bar plot
#' }
#' @export

plot_shap <- function(shapley_result, type = c("auto")) {
type <- match.arg(type)

is_time_dependent <- "time" %in% colnames(shapley_result)

if (type == "auto") {
type <- if (is_time_dependent) "time" else "integrated"
}

if (type == "time") {
ggplot2::ggplot(shapley_result, ggplot2::aes(x = time, y = phi, color = feature)) +
ggplot2::geom_line(size = 1) +
ggplot2::geom_point(size = 2) +
ggplot2::labs(
title = "Shapley Contributions Over Time",
x = "Time",
y = "SHAP value"
) +
ggplot2::theme_minimal(base_size = 14)
} else {
shap_method <- attr(shapley_result, "shap_method", exact = TRUE)
shap_method <- if (is.null(shap_method)) "integral" else shap_method

title <- switch(
shap_method,
meanabs = "Mean Absolute SHAP",
integral = "Integrated SHAP",
"SHAP Summary"
)

shapley_result <- shapley_result[order(shapley_result$phi), ]
shapley_result$feature <- factor(shapley_result$feature, levels = shapley_result$feature)
shapley_result$direction <- ifelse(shapley_result$phi >= 0, "Positive", "Negative")

ggplot2::ggplot(shapley_result, ggplot2::aes(x = phi, y = feature, fill = direction)) +
ggplot2::geom_col(show.legend = FALSE) +
ggplot2::scale_fill_manual(values = c("Positive" = "#4CAF50", "Negative" = "#E53935")) +
ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
ggplot2::labs(
title = title,
x = "SHAP value",
y = "Feature"
) +
ggplot2::theme_minimal(base_size = 14)
}
}





