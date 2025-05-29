plot_pdp <- function(pdp_ice_output, feature,
                     method = "pdp+ice", ids = NULL,
                     which = c("per_time", "integrated"),
                     alpha_ice = 0.2, smooth = FALSE) {

  which <- match.arg(which)
  results <- pdp_ice_output$results
  feature_data <- results[[feature]]
  is_categorical <- is.factor(feature_data) || is.character(feature_data)

  if (which == "integrated") {

    if (is.null(pdp_ice_output$pdp_integrated)) stop("No integrated PDP found.")
    if (!"integrated_surv" %in% names(pdp_ice_output$pdp_integrated)) stop("Invalid integrated PDP.")

    if (method %in% c("ice", "pdp+ice")) {
      ice_data <- results[results$type == "ice", ]

      # Compute integrated survival for each individual by group
      ice_integrated <- do.call(rbind, lapply(split(ice_data, list(ice_data$.id, ice_data[[feature]])), function(df) {
        df <- df[order(df$time), ]
        int_surv <- sum((head(df$surv_prob, -1) + tail(df$surv_prob, -1)) / 2 * diff(df$time))
        data.frame(.id = unique(df$.id), integrated_surv = int_surv, feature_value = unique(df[[feature]]))
      }))
      names(ice_integrated)[3] <- feature
    }

    pdp_integrated <- pdp_ice_output$pdp_integrated

    if (is_categorical) {
      p <- ggplot2::ggplot() +
        ggplot2::geom_boxplot(
          data = if (method %in% c("ice", "pdp+ice")) ice_integrated else NULL,
          mapping = ggplot2::aes_string(x = feature, y = "integrated_surv"),
          alpha = alpha_ice, fill = "pink"
        ) +
        ggplot2::geom_point(
          data = if (method %in% c("pdp", "pdp+ice")) pdp_integrated else NULL,
          mapping = ggplot2::aes_string(x = feature, y = "integrated_surv"),
          shape = 21, size = 3, fill = "black"
        ) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::labs(
          title = paste("Integrated PDP with ICE Boxplot for", feature),
          x = feature,
          y = "Integrated Survival"
        )
      return(p)
    } else {
      p <- ggplot2::ggplot(pdp_integrated, ggplot2::aes_string(x = feature, y = "integrated_surv")) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::labs(
          title = paste("Integrated PDP over Survival Time:", feature),
          x = feature,
          y = "Integrated Survival"
        )
      if (smooth) {
        p <- p + ggplot2::geom_smooth(method = "loess", se = FALSE, color = "steelblue", size = 1.2)
      } else {
        p <- p + ggplot2::geom_line(color = "steelblue", size = 1.2)
      }
      return(p)
    }
  }

  # per-time plot
  plot_data <- results[results$type %in% switch(method,
                                                "pdp" = "pdp",
                                                "ice" = "ice",
                                                "pdp+ice" = c("pdp", "ice"),
                                                stop("Invalid method")), ]

  if (!is.null(ids) && "ice" %in% plot_data$type) {
    plot_data <- plot_data[!(plot_data$type == "ice" & !(plot_data$.id %in% ids)), ]
  }

  if (is_categorical) {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = feature, y = "surv_prob", fill = "type")) +
      ggplot2::facet_wrap(~ time, scales = "free_y") +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(
        title = paste("Survival", toupper(method), "for", feature),
        x = feature,
        y = "Survival Probability",
        fill = "Type"
      )

    if ("ice" %in% plot_data$type) {
      p <- p + ggplot2::geom_boxplot(data = plot_data[plot_data$type == "ice", ],
                                     mapping = ggplot2::aes_string(x = feature, y = "surv_prob", fill = "type"),
                                     alpha = alpha_ice, position = "dodge")
    }

    if ("pdp" %in% plot_data$type) {
      p <- p + ggplot2::stat_summary(data = plot_data[plot_data$type == "pdp", ],
                                     mapping = ggplot2::aes_string(x = feature, y = "surv_prob"),
                                     fun = mean, geom = "point",
                                     shape = 21, fill = "black", size = 3,
                                     position = ggplot2::position_dodge(width = 0.75))
    }

  } else {
    p <- ggplot2::ggplot() +
      ggplot2::facet_wrap(~ time, scales = "free_y") +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(
        title = paste("Survival", toupper(method), "for", feature),
        x = feature,
        y = "Survival Probability",
        color = NULL
      ) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1)))

    if ("ice" %in% plot_data$type) {
      p <- p + ggplot2::geom_line(data = plot_data[plot_data$type == "ice", ],
                                  mapping = ggplot2::aes_string(x = feature, y = "surv_prob",
                                                                group = "interaction(.id, time)", color = "type"),
                                  alpha = alpha_ice)
    }

    if ("pdp" %in% plot_data$type) {
      p <- p + ggplot2::geom_line(data = plot_data[plot_data$type == "pdp", ],
                                  mapping = ggplot2::aes_string(x = feature, y = "surv_prob", group = "time"),
                                  size = 1.2)
    }
  }

  return(p)
}


## test
plot_pdp(pdp_ice_result, feature = "age", which = "per_time")
plot_pdp(pdp_ice_result, feature = "age", which = "integrated", smooth = TRUE)

plot_pdp(pdp_cat_result, feature = "celltype", which = "per_time")
plot_pdp(pdp_cat_result, feature = "celltype", which = "integrated")

