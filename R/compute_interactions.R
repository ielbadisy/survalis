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
#' @param grid.size Integer; number of random grid values / replicates used for Monte Carlo
#'   marginalization (default 30).
#' @param batch.size Reserved for future batching support (currently unused).
#'
#' @details
#' For a target time \eqn{t^{*}}, let \eqn{f(x)} be the predicted survival probability.
#' For feature \eqn{j}, we approximate a decomposition:
#' \deqn{ f(x) \approx f_j(x_j) + f_{-j}(x_{-j}) }
#' by Monte Carlo marginalization over subsets of features using random sampling from the
#' empirical distribution in \code{data}. The reported interaction strength is:
#' \deqn{ H_j(t^{*}) = \sqrt{ \frac{\sum ( f(x) - \{ \tilde f_j(x_j) + \tilde f_{-j}(x_{-j}) \})^2 }{ \sum f(x)^2 } } , }
#' clipped to \[0, 1]. Pairwise heatmaps are computed analogously for \eqn{(j, k)}.
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
#' @export

compute_interactions <- function(model, data, times,
                                 target_time = NULL,
                                 features = NULL,
                                 type = c("1way", "heatmap", "time"),
                                 grid.size = 30,
                                 batch.size = 100) {
  type <- match.arg(type)

  if (type != "time" && is.null(target_time)) stop("target_time must be specified for type = '1way' or 'heatmap'")
  if (is.null(features)) {
    features <- setdiff(colnames(data), c(model$time, model$status))
  }

  predict_function <- get(paste0("predict_", model$learner))
  idx_time <- function(t) which.min(abs(times - t))

  if (type == "1way") {
    res <- lapply(features, function(feature) {
      grid <- unique(data[[feature]])
      grid <- sort(sample(grid, min(length(grid), grid.size)))
      rows <- data[sample(1:nrow(data), grid.size), , drop = FALSE]

      full <- predict_function(model, rows, times)[, idx_time(target_time)]

      marg_j <- rowMeans(replicate(grid.size, {
        rows[[feature]] <- sample(grid, nrow(rows), replace = TRUE)
        predict_function(model, rows, times)[, idx_time(target_time)]
      }))

      marg_rest <- rowMeans(replicate(grid.size, {
        others <- setdiff(features, feature)
        for (f in others) rows[[f]] <- sample(unique(data[[f]]), nrow(rows), replace = TRUE)
        predict_function(model, rows, times)[, idx_time(target_time)]
      }))

      h <- sqrt(sum((full - (marg_j + marg_rest))^2) / sum(full^2))
      h <- pmin(h, 1)
      data.frame(feature = feature, interaction = h)
    })
    return(do.call(rbind, res))
  }

  if (type == "heatmap") {
    pairs <- t(combn(features, 2))
    res <- lapply(seq_len(nrow(pairs)), function(i) {
      f1 <- pairs[i, 1]; f2 <- pairs[i, 2]
      grid1 <- sort(sample(unique(data[[f1]]), min(length(unique(data[[f1]])), grid.size)))
      grid2 <- sort(sample(unique(data[[f2]]), min(length(unique(data[[f2]])), grid.size)))
      rows <- data[sample(1:nrow(data), grid.size), , drop = FALSE]

      full <- predict_function(model, rows, times)[, idx_time(target_time)]

      m1 <- rowMeans(replicate(grid.size, {
        rows[[f1]] <- sample(grid1, nrow(rows), replace = TRUE)
        predict_function(model, rows, times)[, idx_time(target_time)]
      }))
      m2 <- rowMeans(replicate(grid.size, {
        rows[[f2]] <- sample(grid2, nrow(rows), replace = TRUE)
        predict_function(model, rows, times)[, idx_time(target_time)]
      }))

      h <- sqrt(sum((full - (m1 + m2))^2) / sum(full^2))
      data.frame(feature1 = f1, feature2 = f2, interaction = h)
    })

    df <- do.call(rbind, res)
    df_sym <- rbind(
      df,
      data.frame(feature1 = df$feature2, feature2 = df$feature1, interaction = df$interaction),
      data.frame(feature1 = features, feature2 = features, interaction = 0)
    )
    return(df_sym)
  }

  if (type == "time") {
    res <- lapply(times, function(t) {
      lapply(features, function(f) {
        grid <- unique(data[[f]])
        grid <- sort(sample(grid, min(length(grid), grid.size)))
        rows <- data[sample(1:nrow(data), grid.size), , drop = FALSE]

        full <- predict_function(model, rows, times)[, idx_time(t)]

        marg_j <- rowMeans(replicate(grid.size, {
          rows[[f]] <- sample(grid, nrow(rows), replace = TRUE)
          predict_function(model, rows, times)[, idx_time(t)]
        }))

        marg_rest <- rowMeans(replicate(grid.size, {
          others <- setdiff(features, f)
          for (k in others) rows[[k]] <- sample(unique(data[[k]]), nrow(rows), replace = TRUE)
          predict_function(model, rows, times)[, idx_time(t)]
        }))

        h <- sqrt(sum((full - (marg_j + marg_rest))^2) / sum(full^2))
        h <- pmin(h, 1)
        data.frame(feature = f, time = t, interaction = h)
      })
    })

    return(do.call(rbind, unlist(res, recursive = FALSE)))
  }
}


#' Plot Interaction Strengths for Survival Models
#'
#' Visualizes interaction outputs from \code{\link{compute_interactions}} as
#' (i) a ranked bar chart for one-way interactions, (ii) a pairwise heatmap,
#' or (iii) time-varying interaction trajectories.
#'
#' @param object A data frame returned by \code{compute_interactions()} whose
#'   columns match the requested \code{type}.
#' @param type One of \code{"1way"}, \code{"heatmap"}, or \code{"time"}.
#'
#' @details
#' \strong{1way}: Bars rank features by Friedman-H interaction strength at the target time.  
#' \strong{heatmap}: Tiles show pairwise interaction magnitudes (symmetric).  
#' \strong{time}: Lines show interaction strength vs. time for each feature.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' mod <- fit_coxph(Surv(time, status) ~ age + karno + trt, data = veteran)
#' times <- c(80, 160)
#' ia <- compute_interactions(
#'   model = mod,
#'   data = veteran,
#'   times = times,
#'   target_time = 80,
#'   features = c("age", "karno"),
#'   type = "1way",
#'   grid.size = 6
#' )
#' plot_interactions(ia, type = "1way")
#'
#' @export

plot_interactions <- function(object, type = c("1way", "heatmap", "time")) {
  type <- match.arg(type)
  if (type == "1way") {
  ggplot(object, ggplot2::aes(x = interaction, y = reorder(feature, interaction))) +
      geom_col(fill = "steelblue") +
      theme_minimal() +
      labs(x = "Interaction strength (H)", y = "Feature")
  } else if (type == "heatmap") {
    ggplot(object, ggplot2::aes(x = feature1, y = feature2, fill = interaction)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "steelblue") +
      theme_minimal() +
      labs(title = "Pairwise Interaction Heatmap")
  } else {
    ggplot(object, ggplot2::aes(x = time, y = interaction, color = feature)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 2) +
      theme_minimal() +
      labs(title = "Time-Varying Feature Interactions", x = "Time", y = "Interaction Strength")
  }
}
