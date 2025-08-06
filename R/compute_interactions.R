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



plot_interactions <- function(object, type = c("1way", "heatmap", "time")) {
  type <- match.arg(type)
  if (type == "1way") {
    ggplot2::ggplot(object, ggplot2::aes(x = interaction, y = reorder(feature, interaction))) +
      ggplot2::geom_col(fill = "steelblue") +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Interaction strength (H)", y = "Feature")
  } else if (type == "heatmap") {
    ggplot2::ggplot(object, ggplot2::aes(x = feature1, y = feature2, fill = interaction)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Pairwise Interaction Heatmap")
  } else {
    ggplot2::ggplot(object, ggplot2::aes(x = time, y = interaction, color = feature)) +
      ggplot2::geom_line(size = 1.2) +
      ggplot2::geom_point(size = 2) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = "Time-Varying Feature Interactions", x = "Time", y = "Interaction Strength")
  }
}




# Required packages
library(survival)
library(ggplot2)

# Fit a survival model using ranger (must be defined in your framework)
data(veteran, package = "survival")
mod_ranger <- fit_ranger(
  Surv(time, status) ~ age + karno + celltype,
  data = veteran,
  num.trees = 300
)

# Define prediction times and features
times <- c(50, 100, 150, 200, 250, 300, 350)
target_time <- 200
features <- c("age", "karno", "celltype")  # specify a subset or leave NULL to include all

# --- 1. One-way interaction analysis (Friedman's H per feature)
interaction_1way <- compute_interactions(
  model = mod_ranger,
  data = veteran,
  times = times,
  target_time = target_time,
  features = features,
  type = "1way"
)

print(interaction_1way)
plot_interactions(interaction_1way, type = "1way")

# --- 2. Pairwise interaction heatmap
interaction_heatmap <- compute_interactions(
  model = mod_ranger,
  data = veteran,
  times = times,
  target_time = target_time,
  features = features,
  type = "heatmap"
)

print(interaction_heatmap)
plot_interactions(interaction_heatmap, type = "heatmap")

# --- 3. Time-varying interaction strengths
interaction_time <- compute_interactions(
  model = mod_ranger,
  data = veteran,
  times = times,
  features = features,
  type = "time"
)

print(interaction_time)
plot_interactions(interaction_time, type = "time")

