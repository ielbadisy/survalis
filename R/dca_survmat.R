#' Decision Curve Analysis for Right-Censored Survival Outcomes
#'
#' Computes net benefit across a range of risk thresholds at a fixed horizon
#' \code{t_star} for a survival model's predicted risk, alongside the
#' "treat all" and "treat none" reference strategies (Vickers and Elkin 2006),
#' adapted to right-censored outcomes via the same inverse-probability-of-
#' censoring weighting (IPCW) estimator of \eqn{\hat G(t)} used by
#' \code{brier} and \code{auc_survmat}.
#'
#' @param object A \code{Surv} object of length \eqn{n}.
#' @param predicted An \code{n x k} matrix or data frame of survival
#'   probabilities with columns named \code{"t=<time>"}.
#' @param t_star Numeric horizon at which risk is evaluated.
#' @param thresholds Numeric vector of risk thresholds in (0, 1). Default
#'   \code{seq(0.01, 0.99, by = 0.01)}.
#'
#' @details
#' The predicted risk at the horizon is \eqn{\pi_i = 1 - \hat S(t^\star \mid X_i)}.
#' At threshold \code{pt}, subjects with \eqn{\pi_i > pt} are classified
#' high-risk. Using the IPCW identity for the Kaplan-Meier estimator (Graf et
#' al. 1999), the true- and false-positive fractions are estimated as
#' \deqn{
#'   \widehat{\mathrm{TPF}}(pt) = \frac{1}{n}\sum_{i=1}^n
#'     I(\pi_i > pt)\, \frac{\Delta_i\, I(\tilde T_i \le t^\star)}{\hat G(\tilde T_i-)},
#'   \qquad
#'   \widehat{\mathrm{FPF}}(pt) = \frac{1}{n}\sum_{i=1}^n
#'     I(\pi_i > pt)\, \frac{I(\tilde T_i > t^\star)}{\hat G(t^\star)},
#' }
#' and the net benefit is
#' \deqn{
#'   \widehat{\mathrm{NB}}(pt) = \widehat{\mathrm{TPF}}(pt) -
#'     \widehat{\mathrm{FPF}}(pt) \cdot \frac{pt}{1-pt}.
#' }
#' "Treat all" (classify everyone high-risk) uses the same weights with
#' \eqn{I(\pi_i > pt)} replaced by one. "Treat none" has net benefit \eqn{0}
#' by definition. These are the weights of [dwscore()], so a decision-weighted
#' loss equals the weighted area between treat all and the model curve exactly.
#'
#' @return An object of class \code{"dca_survmat"}: a list with \code{curve}
#'   (a data.frame with \code{threshold}, \code{model}, \code{treat_all},
#'   \code{treat_none}) and \code{t_star}.
#'
#' @seealso `plot.dca_survmat()`, `brier()`, `auc_survmat()`
#'
#' @references
#' Vickers AJ, Elkin EB (2006). "Decision Curve Analysis: A Novel Method for
#' Evaluating Prediction Models." \emph{Medical Decision Making}, 26(6), 565-574.
#'
#' Vickers AJ, Van Calster B, Steyerberg EW (2019). "A Simple, Step-by-Step
#' Guide to Interpreting Decision Curve Analysis." \emph{Diagnostic and
#' Prognostic Research}, 3, 18.
#'
#' Graf E, Schmoor C, Sauerbrei W, Schumacher M (1999). "Assessment and
#' Comparison of Prognostic Classification Schemes for Survival Data."
#' \emph{Statistics in Medicine}, 18(17-18), 2529-2545.
#'
#' @examples
#' y <- survival::Surv(time = veteran$time, event = veteran$status)
#' sp <- matrix(stats::plogis(scale(veteran$karno)), ncol = 1,
#'             dimnames = list(NULL, "t=100"))
#' d <- dca_survmat(y, predicted = sp, t_star = 100)
#' head(d$curve)
#' @noRd
dca_survmat <- function(object, predicted, t_star, thresholds = seq(0.01, 0.99, by = 0.01)) {
  if (!inherits(object, "Surv")) stop("object must be a survival object (from Surv())")
  if (any(thresholds <= 0) || any(thresholds >= 1)) {
    stop("`thresholds` must lie strictly between 0 and 1.")
  }

  time <- object[, 1]
  status <- object[, 2]
  n <- length(time)

  t_name <- paste0("t=", t_star)
  if (!(t_name %in% colnames(predicted))) {
    stop("t_star = ", t_star, " not found in predicted survival matrix.")
  }
  surv_prob <- if (is.matrix(predicted)) predicted[, t_name] else predicted[[t_name]]
  risk <- 1 - surv_prob

  w <- .ipcw_weights(object, t_star)
  if (is.null(w)) stop("Net benefit is undefined because G(t_star) is NA or zero.")

  nb <- function(treat, pt) mean(treat * w$w1) - mean(treat * w$w0) * pt / (1 - pt)
  net_benefit_model <- vapply(thresholds, function(pt) nb(risk > pt, pt), numeric(1))
  net_benefit_all <- vapply(thresholds, function(pt) nb(rep(TRUE, n), pt), numeric(1))

  curve <- data.frame(
    threshold = thresholds,
    model = net_benefit_model,
    treat_all = net_benefit_all,
    treat_none = 0
  )

  out <- list(curve = curve, t_star = t_star)
  class(out) <- "dca_survmat"
  out
}

#' @param x A \code{dca_survmat} object.
#' @param title Plot title. If missing, an automatically generated title is
#'   used. Pass \code{NULL} to omit the title entirely.
#' @param ... Additional arguments (unused).
#' @rdname dca_survmat
#' @export
#' @noRd
plot.dca_survmat <- function(x, title, ...) {
  curve <- x$curve
  long <- data.frame(
    threshold = rep(curve$threshold, 3L),
    net_benefit = c(curve$model, curve$treat_all, curve$treat_none),
    strategy = factor(
      rep(c("Model", "Treat all", "Treat none"), each = nrow(curve)),
      levels = c("Model", "Treat all", "Treat none")
    )
  )
  y_floor <- max(min(curve$model, curve$treat_all, curve$treat_none, na.rm = TRUE), -0.1)

  if (missing(title)) title <- sprintf("Decision curve analysis at t = %s", x$t_star)

  p <- ggplot2::ggplot(long, ggplot2::aes(x = threshold, y = net_benefit, colour = strategy)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::coord_cartesian(ylim = c(y_floor, NA)) +
    ggplot2::scale_colour_manual(values = stats::setNames(
      c(.survalis_palette[1], .survalis_palette[2], "grey50"),
      c("Model", "Treat all", "Treat none")
    )) +
    ggplot2::labs(x = "Threshold probability", y = "Net benefit", colour = NULL) +
    theme_survalis()
  if (!is.null(title)) p <- p + ggplot2::labs(title = title)
  p
}

#' @param x A \code{dca_survmat} object.
#' @param ... Additional arguments passed to \code{plot.dca_survmat}.
#' @rdname dca_survmat
#' @noRd
plot_dca <- function(x, ...) {
  plot.dca_survmat(x, ...)
}
