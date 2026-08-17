.survcurve_data <- function(fit) {
  s <- summary(fit)
  strata_vec <- if (!is.null(s$strata)) as.character(s$strata) else "All"

  dt <- data.table::data.table(
    time = s$time,
    surv = s$surv,
    lower = s$lower,
    upper = s$upper,
    n.risk = s$n.risk,
    n.event = s$n.event,
    strata = strata_vec
  )

  starts <- dt[, list(
    time = 0, surv = 1, lower = 1, upper = 1,
    n.risk = max(n.risk), n.event = 0
  ), by = strata]

  out <- data.table::rbindlist(list(starts, dt), use.names = TRUE, fill = TRUE)
  data.table::setorderv(out, c("strata", "time"))
  out[]
}

.survcurve_pvalue <- function(formula, data) {
  sd <- survival::survdiff(formula, data = data)
  df <- length(sd$n) - 1
  if (df < 1) return(NA_real_)
  stats::pchisq(sd$chisq, df = df, lower.tail = FALSE)
}

#' Kaplan-Meier Survival Curve with Risk Table (survminer-Inspired)
#'
#' Produces a styled Kaplan-Meier survival curve with an optional confidence
#' band, log-rank p-value annotation, and an aligned number-at-risk table
#' beneath the curve, in the spirit of \code{survminer::ggsurvplot()} but
#' implemented natively with \pkg{ggplot2}/\pkg{data.table} (no dependency on
#' \pkg{survminer}).
#'
#' @param formula A survival formula, e.g. \code{Surv(time, status) ~ 1} or
#'   \code{Surv(time, status) ~ group}. The right-hand side determines
#'   stratification (as in \code{\link[survival]{survfit}}).
#' @param data A data frame containing the variables in \code{formula}.
#' @param conf.int Logical; draw a confidence ribbon around each curve
#'   (default \code{TRUE}).
#' @param risk.table Logical; draw an aligned number-at-risk table beneath
#'   the curve (default \code{TRUE}). Requires the \pkg{patchwork} package.
#' @param pval Logical; annotate the plot with a log-rank test p-value when
#'   there are 2 or more strata (default \code{TRUE}).
#' @param break.time.by Numeric spacing of x-axis breaks (and risk-table
#'   columns). If \code{NULL} (default), a suitable spacing is chosen
#'   automatically via \code{pretty()}.
#' @param xlab,ylab Axis labels (defaults \code{"Time"} and
#'   \code{"Survival probability"}).
#' @param legend.title Legend/strata title. If \code{NULL} (default), uses
#'   the stratifying variable name, or is omitted entirely for an
#'   unstratified curve.
#' @param title Plot title. \code{NULL} (default) omits it, matching the
#'   package's other \code{plot_*()} functions and suiting journals that
#'   require caption-only figures.
#'
#' @return If \code{risk.table = FALSE}, a single \pkg{ggplot2} object. If
#'   \code{risk.table = TRUE}, a combined \pkg{patchwork} object (KM curve
#'   stacked above the risk table).
#'
#' @examples
#' plot_survcurve(Surv(time, status) ~ 1, data = veteran, risk.table = FALSE)
#' \donttest{
#' plot_survcurve(Surv(time, status) ~ trt, data = veteran)
#' }
#'
#' @seealso [theme_survalis()], [scale_color_survalis()]
#' @export
plot_survcurve <- function(formula, data,
                           conf.int = TRUE,
                           risk.table = TRUE,
                           pval = TRUE,
                           break.time.by = NULL,
                           xlab = "Time",
                           ylab = "Survival probability",
                           legend.title = NULL,
                           title = NULL) {
  stopifnot(is.data.frame(data))

  fit <- survival::survfit(formula, data = data)
  dt <- .survcurve_data(fit)
  n_strata <- length(unique(dt$strata))

  if (is.null(legend.title)) {
    rhs_vars <- all.vars(formula[[3]])
    legend.title <- if (n_strata > 1 && length(rhs_vars)) rhs_vars[1] else "Strata"
  }

  if (is.null(break.time.by)) {
    breaks <- pretty(c(0, max(dt$time)), n = 5)
  } else {
    breaks <- seq(0, max(dt$time), by = break.time.by)
  }
  xlim <- c(0, max(breaks, dt$time))

  p <- ggplot2::ggplot(dt, ggplot2::aes(x = time, y = surv, color = strata))

  if (isTRUE(conf.int)) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper, fill = strata),
      alpha = 0.15, color = NA
    )
  }

  p <- p +
    ggplot2::geom_step(linewidth = 0.9) +
    scale_color_survalis() +
    scale_fill_survalis() +
    ggplot2::scale_x_continuous(breaks = breaks, limits = xlim) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = xlab, y = ylab, color = legend.title, fill = legend.title) +
    theme_survalis()

  if (!is.null(title)) p <- p + ggplot2::labs(title = title)

  if (n_strata < 2) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  if (isTRUE(pval) && n_strata >= 2) {
    pv <- .survcurve_pvalue(formula, data)
    if (is.finite(pv)) {
      p <- p + ggplot2::annotate(
        "text", x = xlim[1], y = 0.05,
        label = paste0("Log-rank p ", format.pval(pv, digits = 3, eps = 1e-4)),
        hjust = 0, size = 3.8
      )
    }
  }

  if (!isTRUE(risk.table)) {
    return(p)
  }

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("Package 'patchwork' is required for risk.table = TRUE; returning the curve without a risk table.")
    return(p)
  }

  risk_dt <- summary(fit, times = breaks, extend = TRUE)
  risk_strata <- if (!is.null(risk_dt$strata)) as.character(risk_dt$strata) else "All"
  risk_tbl <- data.table::data.table(
    time = risk_dt$time,
    strata = risk_strata,
    n.risk = risk_dt$n.risk
  )

  rt <- ggplot2::ggplot(risk_tbl, ggplot2::aes(x = time, y = strata, label = n.risk, color = strata)) +
    ggplot2::geom_text(size = 3.6, show.legend = FALSE) +
    scale_color_survalis() +
    ggplot2::scale_x_continuous(breaks = breaks, limits = xlim) +
    ggplot2::labs(x = xlab, y = NULL, title = "Number at risk") +
    theme_survalis() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = ggplot2::rel(0.85), face = "plain"),
      legend.position = "none"
    )

  patchwork::wrap_plots(p, rt, ncol = 1, heights = c(3, 1))
}
