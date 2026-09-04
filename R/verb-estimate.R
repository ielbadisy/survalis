#' Standardised causal-style contrast for a binary treatment
#'
#' Estimates a covariate-adjusted marginal contrast between the two levels of a
#' binary `treatment` by G-computation (outcome standardisation): a survival
#' model is fitted for `Surv(time, status) ~ treatment * (confounders)`, every
#' subject's outcome is predicted under both treatment levels, and the
#' counterfactual mean survival curves are contrasted.
#'
#' Uncertainty is a nonparametric bootstrap over subjects; by default the model
#' is refitted in each replicate (`refit = TRUE`), which propagates model
#' estimation error. Set `refit = FALSE` for a faster interval that is
#' conditional on the fitted model.
#'
#' @section Interpretation:
#' The estimate has a causal reading only under standard identification
#' assumptions: consistency, no unmeasured confounding given the `formula`
#' covariates, positivity, no interference, and a correctly specified outcome
#' model. For a non-manipulable attribute it is a covariate-adjusted marginal
#' contrast, not a causal effect. This is a plug-in G-formula estimator; it is
#' not doubly robust.
#'
#' @param formula A survival formula whose right-hand side lists the
#'   **confounders** (the treatment is added by the function, with
#'   treatment-by-confounder interactions).
#' @param data A data frame.
#' @param model A learner id (`list_survlearners()$learner`).
#' @param estimand One of `"RMST_diff"` (default), `"RMST_ratio"`,
#'   `"survival_diff"`, `"survival_ratio"`.
#' @param treatment Name of a binary column in `data`. The larger value is the
#'   "treated" (A = 1) level.
#' @param tau Horizon: the RMST restriction time for the RMST estimands, or the
#'   evaluation time for the survival estimands.
#' @param times Prediction grid; `NULL` uses [default_times()]. Must reach `tau`.
#' @param spec Engine arguments for `fit_<model>()` (named list).
#' @param R Bootstrap replicates (default `200`).
#' @param conf Confidence level (default `0.95`).
#' @param refit Logical; refit the model inside each bootstrap replicate
#'   (default `TRUE`).
#' @param seed Optional RNG seed.
#' @param ... Unused.
#'
#' @return A `survalis_estimate` object with `estimand`, `estimate`,
#'   `conf.low`, `conf.high`, `conf`, `treatment`, `levels`, `tau`, `R`,
#'   `refit`, `boot` (replicate estimates) and `curves` (counterfactual mean
#'   survival for each arm on `times`).
#'
#' @examplesIf requireNamespace("survival", quietly = TRUE)
#' est <- estimate(Surv(time, status) ~ age + karno + celltype, veteran,
#'                 model = "coxph", estimand = "RMST_diff",
#'                 treatment = "trt", tau = 300, R = 50, seed = 1)
#' est
#' @export
estimate <- function(formula, data, model, estimand = c("RMST_diff", "RMST_ratio",
                                                        "survival_diff", "survival_ratio"),
                     treatment, tau, times = NULL, spec = list(),
                     R = 200, conf = 0.95, refit = TRUE, seed = NULL, ...) {
  estimand <- match.arg(estimand)
  known <- list_survlearners()$learner
  if (length(model) != 1L || !model %in% known) {
    stop("Unknown model ", sQuote(model), ". See list_survlearners().", call. = FALSE)
  }
  if (missing(treatment) || !is.character(treatment) || length(treatment) != 1L ||
      !treatment %in% names(data)) {
    stop("`treatment` must name a column in `data`.", call. = FALSE)
  }
  if (missing(tau) || !is.finite(tau) || tau <= 0) {
    stop("`tau` must be a positive finite number.", call. = FALSE)
  }
  if (!is.finite(conf) || conf <= 0 || conf >= 1) stop("`conf` must be in (0, 1).", call. = FALSE)
  if (!is.null(seed)) set.seed(as.integer(seed))

  op <- options(survalis.quiet_legacy = TRUE)
  on.exit(options(op), add = TRUE)

  parsed <- .parse_surv_formula(formula, data)
  data <- .complete_cases_df(data, unique(c(all.vars(formula), treatment)))

  lv <- sort(unique(data[[treatment]]))
  if (length(lv) != 2L) stop("`treatment` must have exactly two distinct values.", call. = FALSE)
  lo <- lv[1]; hi <- lv[2]

  confounders <- setdiff(attr(stats::terms(formula), "term.labels"), treatment)
  lhs <- deparse(formula[[2]])
  rhs <- if (length(confounders)) {
    sprintf("%s * (%s)", treatment, paste(confounders, collapse = " + "))
  } else treatment
  work_formula <- stats::as.formula(paste(lhs, "~", rhs))

  status <- if (isTRUE(parsed$recode_status)) {
    as.integer(data[[parsed$status_col]] == parsed$event_value)
  } else as.integer(data[[parsed$status_col]] != 0)
  times <- .resolve_times(times, data[[parsed$time_col]], status)
  if (max(times) < tau) {
    stop("`times` must reach `tau` (max(times) = ", max(times), " < ", tau, ").", call. = FALSE)
  }

  fit_fun  <- get(paste0("fit_", model), envir = asNamespace("survalis"))
  pred_fun <- get(paste0("predict_", model), envir = asNamespace("survalis"))
  spec <- as.list(spec)

  contrast_of <- function(S1bar, S0bar) {
    r1 <- survmat_to_rmst(matrix(S1bar, nrow = 1), times = times, tau = tau)
    r0 <- survmat_to_rmst(matrix(S0bar, nrow = 1), times = times, tau = tau)
    st <- function(v) stats::approx(times, v, xout = tau, rule = 2)$y
    switch(estimand,
      RMST_diff     = r1 - r0,
      RMST_ratio    = r1 / r0,
      survival_diff = st(S1bar) - st(S0bar),
      survival_ratio = st(S1bar) / st(S0bar)
    )
  }

  gcomp <- function(df, do_fit) {
    m <- if (do_fit) {
      do.call(fit_fun, c(list(formula = work_formula, data = df), spec))
    } else fixed_model
    d1 <- df; d1[[treatment]] <- hi
    d0 <- df; d0[[treatment]] <- lo
    S1 <- .finalize_survmat(pred_fun(m, newdata = d1, times = times), times = times)
    S0 <- .finalize_survmat(pred_fun(m, newdata = d0, times = times), times = times)
    list(S1bar = colMeans(as.matrix(S1)), S0bar = colMeans(as.matrix(S0)))
  }

  fixed_model <- do.call(fit_fun, c(list(formula = work_formula, data = data), spec))
  point_curves <- gcomp(data, do_fit = FALSE)
  point <- contrast_of(point_curves$S1bar, point_curves$S0bar)

  n <- nrow(data)
  boot <- vapply(seq_len(R), function(r) {
    idx <- sample.int(n, n, replace = TRUE)
    cc <- tryCatch(gcomp(data[idx, , drop = FALSE], do_fit = isTRUE(refit)),
                   error = function(e) NULL)
    if (is.null(cc)) return(NA_real_)
    contrast_of(cc$S1bar, cc$S0bar)
  }, numeric(1))

  a <- (1 - conf) / 2
  ci <- unname(stats::quantile(boot, c(a, 1 - a), na.rm = TRUE))

  structure(list(
    estimand = estimand, estimate = point, conf.low = ci[1], conf.high = ci[2],
    conf = conf, treatment = treatment, levels = c(lo, hi), tau = tau,
    R = R, refit = isTRUE(refit), model = model, boot = boot,
    curves = data.frame(time = times, surv_treated = point_curves$S1bar,
                        surv_control = point_curves$S0bar)
  ), class = "survalis_estimate")
}

#' @export
print.survalis_estimate <- function(x, ...) {
  lab <- c(RMST_diff = "RMST difference", RMST_ratio = "RMST ratio",
           survival_diff = "survival difference", survival_ratio = "survival ratio")[x$estimand]
  cat(sprintf("<survalis_estimate> %s at tau = %s  (model: %s)\n", lab, format(x$tau), x$model))
  cat(sprintf("  treatment: %s  (%s vs %s)\n", x$treatment, x$levels[2], x$levels[1]))
  cat(sprintf("  estimate: %.4f   %.0f%% CI [%.4f, %.4f]   (%d bootstrap%s)\n",
              x$estimate, 100 * x$conf, x$conf.low, x$conf.high, x$R,
              if (x$refit) ", refit" else ", model fixed"))
  cat("  plug-in G-formula; causal reading needs exchangeability, positivity,\n",
      "  consistency, no interference and a correct outcome model.\n", sep = "")
  invisible(x)
}

#' @export
summary.survalis_estimate <- function(object, ...) {
  data.frame(estimand = object$estimand, tau = object$tau,
             estimate = object$estimate, conf.low = object$conf.low,
             conf.high = object$conf.high, model = object$model)
}

#' @export
plot.survalis_estimate <- function(x, title, ...) {
  df <- x$curves
  long <- rbind(
    data.frame(time = df$time, surv = df$surv_treated, arm = paste0(x$treatment, " = ", x$levels[2])),
    data.frame(time = df$time, surv = df$surv_control, arm = paste0(x$treatment, " = ", x$levels[1]))
  )
  if (missing(title)) {
    title <- sprintf("Counterfactual mean survival (%s)", x$treatment)
  }
  p <- ggplot2::ggplot(long, ggplot2::aes(x = time, y = surv, color = arm)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_vline(xintercept = x$tau, linetype = "dashed", color = "gray50") +
    scale_color_survalis() +
    ggplot2::labs(x = "Time", y = "Survival probability", color = NULL) +
    ggplot2::ylim(0, 1) +
    theme_survalis()
  if (!is.null(title)) p <- p + ggplot2::labs(title = title)
  p
}
