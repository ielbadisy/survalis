#' Kaplan-Meier curves with a number-at-risk table
#'
#' Draws Kaplan-Meier curves in a publication style: step curves with
#' censoring marks, the legend inside the panel, an optional log-rank p-value,
#' and a number-at-risk table coloured by group under the plot. Besides the
#' survival function, it draws the risk function, the cumulative hazard, the
#' complementary log-log diagnostic of proportional hazards, and cumulative
#' incidence under competing risks. A fitted survalis model adds its predicted
#' or covariate-adjusted curves, dashed, over the Kaplan-Meier curves.
#'
#' @param formula A survival formula such as `Surv(time, status) ~ 1` or
#'   `Surv(time, status) ~ group`, as in [survival::survfit()]. When the
#'   status is a factor whose first level is censoring and whose other levels
#'   are competing events, the plot shows the Aalen-Johansen cumulative
#'   incidence of `event`.
#' @param data A data frame holding the variables in `formula`.
#' @param type What to draw: `"survival"` (default) for \eqn{S(t)}, `"risk"`
#'   for \eqn{1 - S(t)}, `"cumhaz"` for the cumulative hazard
#'   \eqn{H(t) = -\log S(t)}, or `"cloglog"` for \eqn{\log\{-\log S(t)\}}
#'   against \eqn{\log t}, where parallel curves support proportional hazards.
#'   Competing-risks data take `"risk"` only, which is also their default.
#' @param event Competing-risks data only: the event(s) whose cumulative
#'   incidence is drawn. `NULL` (default) takes the first event level. Several
#'   events, for example `c("relapse", "death")`, are drawn in one panel, one
#'   colour per event (and group).
#' @param model Optional fitted model from [fit()]. Its curves are drawn dashed
#'   over the Kaplan-Meier curves, so that the two can be compared. Not
#'   available for competing-risks data.
#' @param model.curve `"predicted"` (default) averages the predicted survival
#'   of the subjects in each group, a calibration check of the model.
#'   `"adjusted"` sets the single grouping variable to each of its levels for
#'   every subject and averages the predictions (G-computation), which gives
#'   covariate-adjusted curves.
#' @param title Panel title, for example `"A. Overall survival"`. `NULL`
#'   (default) draws none.
#' @param xlab,ylab Axis labels. `ylab = NULL` uses a label for `type`.
#' @param breaks Time breaks for the x axis and the risk table. `NULL`
#'   (default) uses [pretty()] over the follow-up.
#' @param pval Log-rank p-value when there are two or more groups: `TRUE`
#'   (default) computes it with [survival::survdiff()], `FALSE` omits it, and
#'   a character string is printed as given. It is not computed for
#'   competing-risks data or for `type = "cloglog"`, where only a string is
#'   printed.
#' @param risk.table Logical; draw the number-at-risk table (default `TRUE`).
#'   Needs the 'patchwork' package. Not drawn for `type = "cloglog"`, whose
#'   time axis is on the log scale.
#' @param conf.int Logical; draw pointwise 95% confidence bands (default
#'   `FALSE`).
#' @param conf.type Transformation of the pointwise band, passed to
#'   [survival::survfit()]: `"log-log"` (default), `"log"` (the survfit
#'   default) or `"plain"`.
#' @param colors Group colours. `NULL` (default) uses blue (#0072B2) and
#'   vermillion (#D55E00) for the first two groups, then the Okabe-Ito colours.
#' @param legend.pos Legend position inside the panel, `c(x, y)` in `[0, 1]`
#'   units of the panel. `NULL` (default) chooses a corner away from the
#'   curves.
#' @param pval.pos P-value position, `c(time, y)` in data units. `NULL`
#'   (default) chooses a corner away from the curves.
#' @param base_size Base font size (default `10`).
#' @param rmst Restricted mean survival time horizon \eqn{\tau}, survival
#'   curves only. With two groups the area between the curves up to
#'   \eqn{\tau} (the RMST difference) is shaded; otherwise the area under each
#'   curve (the RMST). A dashed line marks \eqn{\tau}. `NULL` (default) draws
#'   none.
#' @param naive Competing-risks data only: also draw, dashed, the naive
#'   \eqn{1 - } Kaplan-Meier estimate of each event that censors the competing
#'   events (default `FALSE`), to show how far it overstates the cumulative
#'   incidence.
#' @param table.height Height of the risk table relative to the curve panel.
#'   `NULL` (default) scales it with the number of rows.
#' @param subtitle,caption Plot subtitle and caption; `NULL` (default) draws
#'   none.
#'
#' @return A 'patchwork' object (curve over risk table) when the risk table is
#'   drawn, otherwise a 'ggplot2' object. Panels combine with patchwork
#'   operators, for example `kmplot(...) | kmplot(...)`, and `p[[1]]` is the
#'   curve panel of a patchwork `p`, to which layers can be added.
#'
#' @examples
#' kmplot(Surv(time, status) ~ trt, data = veteran, risk.table = FALSE)
#' kmplot(Surv(time, status) ~ trt, data = veteran, type = "cloglog")
#'
#' # Competing risks: progression versus death without progression
#' d <- survival::mgus2
#' d$etime <- ifelse(d$pstat == 1, d$ptime, d$futime)
#' d$event <- factor(ifelse(d$pstat == 1, 1, 2 * d$death), 0:2,
#'                   c("censored", "progression", "death"))
#' kmplot(Surv(etime, event) ~ sex, data = d, risk.table = FALSE)
#' kmplot(Surv(etime, event) ~ 1, data = d, event = c("progression", "death"),
#'        naive = TRUE, risk.table = FALSE)
#' kmplot(Surv(time, status) ~ trt, data = veteran, rmst = 365, risk.table = FALSE)
#' \donttest{
#' if (requireNamespace("patchwork", quietly = TRUE)) {
#'   kmplot(Surv(time, status) ~ trt, data = veteran,
#'          title = "A. Overall survival", xlab = "Days",
#'          breaks = seq(0, 1000, 200))
#'
#'   # Covariate-adjusted curves from a Cox model, over the Kaplan-Meier curves
#'   m <- fit(Surv(time, status) ~ trt + karno + age, data = veteran,
#'            model = "coxph")
#'   kmplot(Surv(time, status) ~ trt, data = veteran, model = m,
#'          model.curve = "adjusted", breaks = seq(0, 1000, 200))
#' }
#' }
#' @export
kmplot <- function(formula, data,
                   type = c("survival", "risk", "cumhaz", "cloglog"),
                   event = NULL, model = NULL,
                   model.curve = c("predicted", "adjusted"),
                   title = NULL, xlab = "Time", ylab = NULL,
                   breaks = NULL, pval = TRUE, risk.table = TRUE,
                   conf.int = FALSE, colors = NULL, legend.pos = NULL,
                   pval.pos = NULL, base_size = 10,
                   conf.type = c("log-log", "log", "plain"), rmst = NULL, naive = FALSE,
                   table.height = NULL, subtitle = NULL, caption = NULL) {
  stopifnot(is.data.frame(data))
  conf.type <- match.arg(conf.type)
  fit <- survival::survfit(formula, data = data, conf.type = conf.type)
  multi <- inherits(fit, "survfitms")
  if (multi && missing(type)) type <- "risk"
  type <- match.arg(type)
  model.curve <- match.arg(model.curve)
  if (multi) {
    if (type != "risk")
      stop("Competing-risks data take `type = \"risk\"` (cumulative incidence) only.",
           call. = FALSE)
    if (!is.null(model))
      stop("`model` curves are not available for competing-risks data.", call. = FALSE)
    events <- fit$states[-1]
    if (is.null(event)) event <- events[1]
    if (!all(event %in% events))
      stop("`event` must be one of: ", paste(events, collapse = ", "), ".", call. = FALSE)
  } else if (isTRUE(naive)) stop("`naive` applies to competing-risks data only.", call. = FALSE)
  if (!is.null(rmst) && (multi || type != "survival"))
    stop("`rmst` shading needs survival curves (`type = \"survival\"`).", call. = FALSE)
  loglog <- type == "cloglog"
  if (loglog) risk.table <- FALSE
  if (is.null(ylab)) {
    ylab <- switch(type, survival = "Survival", cumhaz = "Cumulative hazard",
                   cloglog = "log(-log survival)",
                   risk = if (!multi) "Risk" else if (length(event) == 1)
                     paste("Cumulative incidence of", event) else "Cumulative incidence")
  }

  curve <- .kmplot_events(fit, type, event)
  groups <- levels(.kmplot_curve(fit, type, event[1])$group)
  keys <- levels(curve$group)
  ng <- length(groups)
  pal <- .kmplot_colors(length(keys), colors)

  tmax <- max(fit$time)
  if (is.null(breaks)) breaks <- pretty(c(0, tmax), n = 6)
  xlim <- c(0, max(breaks, tmax) * 1.03)
  if (is.null(legend.pos)) legend.pos <- .kmplot_legend_corner(curve, xlim, type)

  mc <- NULL
  if (!is.null(model)) mc <- .kmplot_model(model, formula, data, fit, groups,
                                           type, model.curve)
  if (isTRUE(naive)) mc <- .kmplot_naive(formula, data, event, conf.type)
  lt <- c(if (multi) "Aalen-Johansen" else "Kaplan-Meier",
          if (multi) "1 - Kaplan-Meier" else if (model.curve == "adjusted") "Adjusted" else "Predicted")

  p <- ggplot2::ggplot(curve, ggplot2::aes(.data$time, .data$est,
                                           colour = .data$group))
  if (!is.null(rmst)) p <- p + .kmplot_rmst(curve, rmst, pal, base_size)
  if (isTRUE(conf.int)) {
    p <- p + ggplot2::geom_ribbon(
      data = .kmplot_stepband(curve),
      ggplot2::aes(ymin = .data$lo, ymax = .data$hi, fill = .data$group),
      alpha = 0.12, colour = NA
    )
  }
  if (isTRUE(conf.int) || (!is.null(rmst) && ng != 2))
    p <- p + ggplot2::scale_fill_manual(values = pal, guide = "none")
  if (is.null(mc)) {
    p <- p + ggplot2::geom_step(linewidth = 0.6)
  } else {
    p <- p +
      ggplot2::geom_step(ggplot2::aes(linetype = lt[1]), linewidth = 0.6) +
      (if (multi) ggplot2::geom_step else ggplot2::geom_line)(
        data = mc, ggplot2::aes(linetype = lt[2]), linewidth = 0.6) +
      ggplot2::scale_linetype_manual(values = stats::setNames(c("solid", "22"), lt),
                                     breaks = lt,
                                     guide = ggplot2::guide_legend(position = "top"))
  }
  if (!loglog) {
    p <- p + ggplot2::geom_point(data = curve[curve$cens > 0, , drop = FALSE],
                                 shape = 3, size = 1.4)
  }
  p <- p + switch(type,
    survival = , risk = ggplot2::scale_y_continuous(
      labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)),
    cumhaz = ggplot2::scale_y_continuous(limits = c(0, NA)),
    cloglog = ggplot2::scale_y_continuous()
  )
  p <- p + if (loglog) ggplot2::scale_x_log10() else
    ggplot2::scale_x_continuous(breaks = breaks, limits = xlim)
  p <- p +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::labs(x = xlab, y = ylab, title = title, subtitle = subtitle, caption = caption,
                  linetype = NULL) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = if (length(keys) > 1 || !is.null(mc)) "inside" else "none",
      legend.position.inside = legend.pos,
      legend.justification = legend.pos,
      legend.title = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.key.width = ggplot2::unit(1.6, "lines"),
      legend.margin = ggplot2::margin(0, 0, 0, 0)
    )
  if (length(keys) < 2) p <- p + ggplot2::guides(colour = "none")

  if ((multi || loglog) && !is.character(pval)) pval <- FALSE
  label <- .kmplot_pval(pval, formula, data, ng)
  if (!is.null(label)) {
    if (is.null(pval.pos)) pval.pos <- .kmplot_pval_pos(curve, xlim, type, legend.pos)
    p <- p + ggplot2::annotate("text", x = pval.pos[1], y = pval.pos[2],
                               hjust = if (length(pval.pos) > 2) pval.pos[3] else 0,
                               size = base_size * 0.3, label = label)
  }

  if (!isTRUE(risk.table)) return(p)
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("Package 'patchwork' is needed for the risk table; returning the curve only.",
            call. = FALSE)
    return(p)
  }

  nr <- summary(fit, times = breaks, extend = TRUE)
  n <- if (is.matrix(nr$n.risk)) nr$n.risk[, 1] else nr$n.risk
  tab <- data.frame(
    time = nr$time, n = n,
    group = if (is.null(nr$strata)) "All" else .kmplot_labels(as.character(nr$strata))
  )
  tab$group <- factor(tab$group, levels = rev(groups))
  # rows take the curve colours when the colours code the groups only
  tcol <- if (identical(keys, groups)) rev(pal) else rep("grey20", ng)
  r <- ggplot2::ggplot(tab, ggplot2::aes(.data$time, .data$group,
                                         label = .data$n, colour = .data$group)) +
    ggplot2::geom_text(size = base_size * 0.3) +
    ggplot2::scale_x_continuous(breaks = breaks, limits = xlim) +
    ggplot2::scale_colour_manual(values = tcol) +
    ggplot2::labs(x = NULL, y = NULL, subtitle = "Number at risk") +
    ggplot2::theme_classic(base_size = base_size - 1) +
    ggplot2::theme(legend.position = "none",
                   plot.subtitle = ggplot2::element_text(size = base_size - 1),
                   axis.line = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank())
  # the curve panel carries the time axis; the table needs only its rows
  if (is.null(table.height)) table.height <- 0.08 + 0.07 * ng
  p / r + patchwork::plot_layout(heights = c(1, table.height))
}

# Curves of one or several competing events. With several, the colour key is the event
# (one group) or "group, event".
.kmplot_events <- function(fit, type, event) {
  if (length(event) < 2) return(.kmplot_curve(fit, type, event))
  d <- lapply(event, function(e) {
    x <- .kmplot_curve(fit, type, e)
    x$event <- e
    x
  })
  d <- do.call(rbind, d)
  ng <- nlevels(d$group)
  key <- if (ng == 1) d$event else paste(d$group, d$event, sep = ", ")
  lev <- if (ng == 1) event else paste(rep(levels(d$group), each = length(event)), event, sep = ", ")
  d$group <- factor(key, lev)
  d$event <- NULL
  d
}

# Naive 1 - Kaplan-Meier of each event, censoring the competing events, keyed as the
# Aalen-Johansen curves.
.kmplot_naive <- function(formula, data, event, conf.type) {
  data$.kmt <- eval(formula[[2]][[2]], data, environment(formula))
  st <- eval(formula[[2]][[3]], data, environment(formula))
  f <- formula
  f[[2]] <- quote(survival::Surv(.kmt, .kme))
  d <- do.call(rbind, lapply(event, function(e) {
    data$.kme <- as.integer(st == e)
    x <- .kmplot_curve(survival::survfit(f, data = data, conf.type = conf.type), "risk")
    x$event <- e
    x
  }))
  ng <- nlevels(d$group)
  if (length(event) > 1)
    d$group <- if (ng == 1) d$event else paste(d$group, d$event, sep = ", ")
  d$event <- NULL
  d
}

# RMST shading up to tau: between the curves for two groups, under each curve otherwise,
# plus a dashed line at tau.
.kmplot_rmst <- function(curve, tau, pal, base_size) {
  sf <- lapply(split(curve, curve$group), function(x) {
    x <- x[x$time <= tau, , drop = FALSE]
    function(t) x$est[findInterval(t, x$time)]
  })
  u <- sort(unique(c(curve$time[curve$time <= tau], tau)))
  # step form: the value at u[k] holds on [u[k], u[k + 1])
  step <- function(f) rep(f(u[-length(u)]), each = 2)
  xs <- rep(u, each = 2)[-c(1, 2 * length(u))]
  vline <- list(ggplot2::geom_vline(xintercept = tau, linetype = 2, colour = "grey50"),
                ggplot2::annotate("text", x = tau, y = 0.02, hjust = -0.1, size = base_size * 0.28,
                                  label = paste("tau ==", tau), parse = TRUE, colour = "grey30"))
  if (length(sf) == 2) {
    a <- step(sf[[1]]); b <- step(sf[[2]])
    d <- data.frame(time = xs, lo = pmin(a, b), hi = pmax(a, b))
    return(c(list(ggplot2::geom_ribbon(data = d, ggplot2::aes(.data$time, ymin = .data$lo,
                                                               ymax = .data$hi),
                                       inherit.aes = FALSE, fill = "grey60", alpha = 0.35)), vline))
  }
  d <- do.call(rbind, lapply(names(sf), function(g)
    data.frame(time = xs, hi = step(sf[[g]]), group = factor(g, levels(curve$group)))))
  c(list(ggplot2::geom_ribbon(data = d, ggplot2::aes(.data$time, ymin = 0, ymax = .data$hi,
                                                     fill = .data$group),
                              inherit.aes = FALSE, alpha = 0.15)), vline)
}

# Step data for each group, with a row at time 0 so every curve starts at
# S = 1 (risk 0). Columns: time, est, lo, hi, cens, group.
.kmplot_curve <- function(fit, type, event = NULL) {
  ns <- if (is.null(fit$strata)) length(fit$time) else fit$strata
  labs <- if (is.null(fit$strata)) "All" else .kmplot_labels(names(fit$strata))
  col <- function(x, j) if (is.null(x)) NA_real_ else if (is.matrix(x)) x[, j] else x
  if (inherits(fit, "survfitms")) {
    j <- match(event, fit$states)
    d <- data.frame(time = fit$time, est = fit$pstate[, j], lo = col(fit$lower, j),
                    hi = col(fit$upper, j), cens = col(fit$n.censor, 1),
                    group = rep(labs, ns))
    start <- 0
  } else {
    d <- data.frame(time = fit$time, est = fit$surv, lo = col(fit$lower, 1),
                    hi = col(fit$upper, 1), cens = fit$n.censor,
                    group = rep(labs, ns))
    start <- 1
  }
  d <- do.call(rbind, lapply(labs, function(g) {
    rbind(data.frame(time = 0, est = start, lo = start, hi = start, cens = 0, group = g),
          d[d$group == g, , drop = FALSE])
  }))
  if (start == 1) d <- .kmplot_transform(d, type)
  d$group <- factor(d$group, levels = labs)
  rownames(d) <- NULL
  if (type == "cloglog") d <- d[is.finite(d$est) & d$time > 0, , drop = FALSE]
  d
}

# Survival scale (est, lo, hi) to the scale of `type`. The transforms are
# decreasing in S, so the band limits swap.
.kmplot_transform <- function(d, type) {
  f <- switch(type, survival = identity, risk = function(s) 1 - s,
              cumhaz = function(s) -log(s), cloglog = function(s) log(-log(s)))
  lo <- d$lo
  d$est <- f(d$est)
  if (type == "survival") return(d)
  d$lo <- f(d$hi)
  d$hi <- f(lo)
  for (v in c("est", "lo", "hi")) d[[v]][!is.finite(d[[v]])] <- NA_real_
  d
}

# Band in step form: each value holds until the next event time.
.kmplot_stepband <- function(curve) {
  do.call(rbind, lapply(split(curve, curve$group), function(x) {
    k <- nrow(x)
    if (k < 2) return(x)
    i <- c(1, rep(2:k, each = 2))
    j <- c(rep(1:(k - 1), each = 2), k)
    data.frame(time = x$time[i], est = x$est[j], lo = x$lo[j], hi = x$hi[j],
               group = x$group[1])
  }))
}

# Legend corner. Survival and risk: top right or bottom right, whichever
# leaves more room between the legend and the curves at the right of the
# panel. The cumulative hazard and the log-log curve increase, so the top
# left and bottom right are free.
.kmplot_legend_corner <- function(curve, xlim, type = "survival") {
  if (type == "cumhaz") return(c(0.02, 0.98))
  if (type == "cloglog") return(c(0.98, 0.02))
  right <- curve$time >= 0.6 * xlim[2]
  last <- !duplicated(curve$group, fromLast = TRUE)
  v <- curve$est[right | last]
  if (1 - max(v) > min(v)) c(0.98, 0.98) else c(0.98, 0.02)
}

# P-value position c(x, y, hjust) in data units.
.kmplot_pval_pos <- function(curve, xlim, type, legend.pos) {
  if (type %in% c("cumhaz", "cloglog")) {
    yr <- range(curve$est[is.finite(curve$est)])
    if (type == "cumhaz") yr[1] <- 0
    return(c(0.98 * xlim[2], yr[1] + 0.05 * diff(yr), 1))
  }
  # Bottom left (top left for risk) unless curves cross the label there: then the first spot with the
  # fewest crossings, going up the free right column (not the legend corner), then up the left one.
  ys <- c(0.08, 0.15, 0.22, 0.29, 0.36); if (type == "risk") ys <- 1 - ys
  lg <- legend.pos[1] > 0.5 && (legend.pos[2] > 0.5) == (type == "risk")
  cand <- c(list(c(0.02 * xlim[2], ys[1], 0)), if (!lg) lapply(ys, function(y) c(0.98 * xlim[2], y, 1)),
            lapply(ys[-1], function(y) c(0.02 * xlim[2], y, 0)))
  w <- 0.35 * xlim[2]
  segs <- lapply(split(curve[c("time", "est")], curve$group), function(g) {
    n <- nrow(g)
    data.frame(a = g$time[-n], b = g$time[-1], lo = pmin(g$est[-n], g$est[-1]), hi = pmax(g$est[-n], g$est[-1]))
  })
  hits <- vapply(cand, function(c) {
    x0 <- if (c[3] == 0) c[1] else c[1] - w
    sum(vapply(segs, function(s) any(s$b >= x0 & s$a <= x0 + w & s$hi >= c[2] - 0.05 & s$lo <= c[2] + 0.05,
                                     na.rm = TRUE), logical(1)))
  }, numeric(1))
  cand[[which.min(hits)]]
}

# Model curves on a grid of 100 times, one per Kaplan-Meier group: the mean
# prediction within each group ("predicted"), or the mean prediction with the
# grouping variable set to each level for everyone ("adjusted").
.kmplot_model <- function(model, formula, data, fit, groups, type, how) {
  if (!inherits(model, "survalis_fit"))
    stop("`model` must be a model fitted by fit().", call. = FALSE)
  vars <- unique(c(all.vars(formula), all.vars(model$formula)))
  vars <- intersect(vars, names(data))
  data <- data[stats::complete.cases(data[vars]), , drop = FALSE]
  grid <- seq(0, max(fit$time), length.out = 101)[-1]
  # The plot time may be a rescaling of the model's time column (for example
  # Surv(time / 365.25, status) for a model fitted in days). Predict on the
  # model's scale and plot on the kmplot scale.
  t_km <- eval(formula[[2]][[2]], data, environment(formula))
  t_model <- data[[model$time]]
  if (is.null(t_model))
    stop("`data` must contain the model's time column `", model$time, "`.", call. = FALSE)
  ok <- is.finite(t_km) & is.finite(t_model) & t_model > 0
  r <- unique(signif(t_km[ok] / t_model[ok], 8))
  if (length(r) != 1L || !is.finite(r) || r <= 0)
    stop("With `model`, the time in `formula` must be the model's time column `",
         model$time, "` or a rescaling of it.", call. = FALSE)
  mean_surv <- function(nd) {
    c(1, colMeans(as.matrix(stats::predict(model, newdata = nd, times = grid / r,
                                           type = "survival"))))
  }
  rhs <- all.vars(stats::delete.response(stats::terms(formula)))
  if (how == "adjusted") {
    rhs_terms <- attr(stats::terms(formula), "term.labels")
    if (length(rhs) != 1L || !identical(rhs_terms, rhs))
      stop("`model.curve = \"adjusted\"` needs exactly one grouping variable in ",
           "`formula`, given as a plain column name.", call. = FALSE)
    if (!rhs %in% .model_features(model, names(data)))
      warning("The model does not use `", rhs, "`, so the adjusted curves are ",
              "identical.", call. = FALSE)
    x <- data[[rhs]]
    lev <- if (is.factor(x)) levels(droplevels(x)) else sort(unique(x))
    est <- lapply(lev, function(l) {
      nd <- data
      nd[[rhs]] <- if (is.factor(x)) factor(l, levels = levels(x)) else l
      mean_surv(nd)
    })
    lab <- as.character(lev)
  } else {
    lab <- if (length(rhs)) .kmplot_group(formula, data) else rep("All", nrow(data))
    est <- lapply(groups, function(g) mean_surv(data[lab == g, , drop = FALSE]))
    lab <- groups
  }
  d <- data.frame(time = rep(c(0, grid), length(lab)), est = unlist(est),
                  lo = NA_real_, hi = NA_real_, group = rep(lab, each = length(grid) + 1))
  d <- .kmplot_transform(d, type)
  if (type == "cloglog") d <- d[is.finite(d$est) & d$time > 0, , drop = FALSE]
  d$group <- factor(d$group, levels = groups)
  d[!is.na(d$group), , drop = FALSE]
}

# Group label of each row of `data`, as survfit() names its strata.
.kmplot_group <- function(formula, data) {
  mf <- stats::model.frame(stats::delete.response(stats::terms(formula)), data,
                           na.action = stats::na.pass)
  as.character(do.call(survival::strata, c(as.list(mf), shortlabel = TRUE, sep = ", ")))
}

# "trt=1" -> "1"; "trt=1, sex=M" -> "1, M"
.kmplot_labels <- function(x) gsub("(^|, )[^=,]+=", "\\1", x)

.kmplot_colors <- function(n, colors) {
  if (!is.null(colors)) {
    if (length(colors) < n) stop("`colors` needs at least ", n, " values.", call. = FALSE)
    return(unname(colors[seq_len(n)]))
  }
  .okabe(n)
}

.kmplot_pval <- function(pval, formula, data, ng) {
  if (is.character(pval)) return(pval[1])
  if (!isTRUE(pval) || ng < 2) return(NULL)
  sd <- survival::survdiff(formula, data = data)
  p <- stats::pchisq(sd$chisq, df = length(sd$n) - 1, lower.tail = FALSE)
  if (!is.finite(p)) return(NULL)
  sprintf("Log-rank p %s", if (p < 0.001) "< 0.001" else sprintf("= %.3f", p))
}
