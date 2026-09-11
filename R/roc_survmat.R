#' Time-Dependent ROC Curve from a Survival-Probability Matrix
#'
#' Computes the full cumulative/dynamic time-dependent ROC curve (sensitivity
#' vs. specificity across risk-score thresholds) at a fixed horizon
#' \code{t_star}, using the same inverse-probability-of-censoring weighted
#' (IPCW) case/control definition as \code{\link{auc_survmat}} (Uno et al.
#' 2007; Heagerty and Zheng 2005). Its trapezoidal AUC closely tracks
#' \code{\link{auc_survmat}}'s pairwise-concordance AUC (small numerical
#' differences arise from tie handling between the step-function ROC curve
#' and the exact pairwise comparison).
#'
#' @param object A \code{\link[survival]{Surv}} object of length \eqn{n}.
#' @param predicted An \code{n x k} matrix or data frame of survival
#'   probabilities with columns named \code{"t=<time>"}.
#' @param t_star Numeric time at which to evaluate the ROC curve.
#'
#' @details
#' As in \code{\link{auc_survmat}}, cases are subjects with an observed event
#' strictly before \code{t_star} (IPCW-weighted by \eqn{1/\hat G(\tilde T_i)}),
#' and controls are subjects known to survive beyond \code{t_star}. The risk
#' score is \eqn{1 - \hat S(t^\star \mid x)}. Sensitivity and specificity are
#' swept over every observed risk-score value \eqn{c}:
#' \deqn{
#'   \widehat{\mathrm{Se}}(c) = \frac{\sum_{i \in \text{cases}} w_i\, I(\text{risk}_i \ge c)}
#'   {\sum_{i \in \text{cases}} w_i}, \qquad
#'   \widehat{\mathrm{Sp}}(c) = \frac{\sum_{i \in \text{controls}} I(\text{risk}_i < c)}
#'   {|\text{controls}|}.
#' }
#'
#' @return An object of class \code{"roc_survmat"}: a list with \code{curve}
#'   (a data.frame with \code{threshold}, \code{sensitivity},
#'   \code{specificity}), \code{auc}, and \code{t_star}.
#'
#' @seealso [auc_survmat()] for the scalar AUC at one horizon,
#'   [timeroc_survmat()] for the AUC(t) trajectory across horizons,
#'   [plot.roc_survmat()].
#'
#' @references
#' Uno H, Cai T, Tian L, Wei LJ (2007). Evaluating prediction rules for
#' t-year survivors with censored regression models. \doi{10.1198/016214507000000149}
#'
#' Heagerty PJ, Zheng Y (2005). Survival model predictive accuracy and ROC
#' curves. \doi{10.1111/j.0006-341X.2005.030814.x}
#'
#' @examples
#' y <- survival::Surv(time = veteran$time, event = veteran$status)
#' sp <- matrix(stats::plogis(scale(veteran$karno)), ncol = 1,
#'             dimnames = list(NULL, "t=100"))
#' r <- roc_survmat(y, predicted = sp, t_star = 100)
#' r$auc
#' head(r$curve)
#' @export
roc_survmat <- function(object, predicted, t_star) {
  if (!inherits(object, "Surv")) stop("object must be a survival object (from Surv())")

  time <- object[, 1]
  status <- object[, 2]

  t_name <- paste0("t=", t_star)
  if (!(t_name %in% colnames(predicted))) {
    stop("t_star = ", t_star, " not found in predicted survival matrix.")
  }
  surv_prob <- if (is.matrix(predicted)) predicted[, t_name] else predicted[[t_name]]
  risk_score <- 1 - surv_prob

  cases <- which(time < t_star & status == 1)
  controls <- which(time > t_star)
  if (!length(cases) || !length(controls)) {
    stop("ROC curve is undefined because no comparable case/control pairs are available.")
  }

  km_fit <- survival::survfit(survival::Surv(time, 1 - status) ~ 1)
  eval_times <- sort(unique(c(t_star, time[cases])))
  G_all <- summary(km_fit, times = eval_times, extend = TRUE)$surv
  names(G_all) <- as.character(eval_times)

  G_case <- G_all[as.character(time[cases])]
  valid_cases <- which(is.finite(G_case) & G_case > 0)
  if (!length(valid_cases)) stop("ROC curve is undefined because censoring weights are not available.")

  cases <- cases[valid_cases]
  case_weights <- 1 / G_case[valid_cases]
  w_total <- sum(case_weights)

  risk_cases <- risk_score[cases]
  risk_controls <- risk_score[controls]
  n_controls <- length(controls)

  thresholds <- sort(unique(risk_score))

  curve <- do.call(rbind, lapply(thresholds, function(c_) {
    sens <- sum(case_weights[risk_cases >= c_]) / w_total
    spec <- sum(risk_controls < c_) / n_controls
    data.frame(threshold = c_, sensitivity = sens, specificity = spec)
  }))
  curve <- curve[order(curve$threshold), , drop = FALSE]
  rownames(curve) <- NULL

  fpr <- 1 - curve$specificity
  ord <- order(fpr)
  auc_val <- abs(pracma::trapz(fpr[ord], curve$sensitivity[ord]))

  out <- list(curve = curve, auc = auc_val, t_star = t_star)
  class(out) <- "roc_survmat"
  out
}

#' @param x A \code{roc_survmat} object.
#' @param title Plot title. If missing, an automatically generated title is
#'   used (includes the AUC). Pass \code{NULL} to omit the title entirely.
#' @param ... Additional arguments (unused).
#' @rdname roc_survmat
#' @export
plot.roc_survmat <- function(x, title, ...) {
  if (missing(title)) {
    title <- sprintf("Time-dependent ROC at t = %s (AUC = %.3f)", x$t_star, x$auc)
  }
  p <- ggplot2::ggplot(x$curve, ggplot2::aes(x = 1 - specificity, y = sensitivity)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey70", linetype = "dashed") +
    ggplot2::geom_line(colour = .survalis_palette[1], linewidth = 0.9) +
    ggplot2::labs(x = "1 - Specificity", y = "Sensitivity") +
    theme_survalis()
  if (!is.null(title)) p <- p + ggplot2::labs(title = title)
  p
}
