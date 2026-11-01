#' Fit a Survival Support Vector Machine
#'
#' Fits a kernel survival support vector regression on the log-time scale and
#' returns an \code{mlsurv_model}. The engine is implemented natively in
#' \pkg{survalis}.
#'
#' @section Method:
#' The model is an accelerated failure time model
#' \eqn{\log T = f(x) + \varepsilon} with
#' \eqn{f(x) = \sum_i \beta_i k(x, x_i) + b} in the reproducing kernel Hilbert
#' space of the kernel \eqn{k}. With \eqn{y_i = \log T_i}, \eqn{f} minimizes
#' \deqn{\frac12 \|f\|_k^2 + \frac{\gamma}{2} \sum_{i:\,\delta_i = 1} (y_i - f(x_i))^2
#'   + \frac{\gamma}{2} \sum_{i:\,\delta_i = 0} \max(0,\, y_i - f(x_i))^2 ,}
#' the survival support vector regression of Shivaswamy et al. (2007) and Van
#' Belle et al. (2011) with squared losses: an event time is fitted from both
#' sides, a censored time is only penalized when the prediction falls below
#' it. The problem is strictly convex and piecewise quadratic, and is solved
#' exactly by an active-set Newton method (Chapelle, 2007). Given the active
#' set \eqn{A} (events, and censored points predicted too early), the optimum
#' solves the least-squares SVM system
#' \deqn{\begin{pmatrix} 0 & 1^\top \\ 1 & K_{AA} + I/\gamma \end{pmatrix}
#'   \begin{pmatrix} b \\ \beta_A \end{pmatrix} =
#'   \begin{pmatrix} 0 \\ y_A \end{pmatrix}, \qquad \beta_{A^c} = 0,}
#' by one Cholesky factorization; \eqn{A} is updated, with a line search,
#' until it no longer changes. No quadratic programming solver is involved.
#' Censored points outside \eqn{A} have \eqn{\beta_i = 0} and are not stored.
#'
#' Survival curves follow from the AFT form:
#' \eqn{S(t \mid x) = S_\varepsilon(\log t - f(x))}, where \eqn{S_\varepsilon}
#' is the Kaplan-Meier estimate of the cross-fitted residuals
#' \eqn{y_i - f_{-k(i)}(x_i)}: each residual comes from a model fitted without
#' the observation's fold, so the curves are not made overconfident by the
#' in-sample fit. See \code{predict_survsvm} for alternatives.
#'
#' @section Kernels:
#' \describe{
#'   \item{\code{"add_kernel"}}{Clinical (additive) kernel of Daemen et al.
#'     (2012): the average over predictors of
#'     \eqn{\max(0, 1 - |x_j - z_j| / r_j)} for numeric and ordered predictors
#'     (\eqn{r_j} the training range) and \eqn{1\{x_j = z_j\}} for unordered
#'     factors. Scale free, positive semi-definite, and a factor counts as one
#'     predictor. Constant predictors add a constant.}
#'   \item{\code{"lin_kernel"}}{\eqn{x^\top z / p} on standardized columns
#'     (factors in treatment coding).}
#'   \item{\code{"rbf_kernel"}}{\eqn{\exp(-\|x - z\|^2 / (2 \sigma^2))} on
#'     standardized columns; \code{kernel.pars} is \eqn{\sigma^2}, default
#'     the number of columns \eqn{p}.}
#'   \item{\code{"poly_kernel"}}{\eqn{(x^\top z / p + c)^d} on standardized
#'     columns; \code{kernel.pars} is \eqn{(c, d)}, default \eqn{(1, 2)}.}
#' }
#'
#' @param formula A survival formula of the form `Surv(time, status) ~ predictors`.
#' @param data A `data.frame` containing all variables referenced in `formula`.
#' @param type Only \code{"regression"} is supported: the ranking types of
#'   \pkg{survivalsvm} give a risk score, not a survival time, and therefore no
#'   survival curve.
#' @param gamma.mu Error weight \eqn{\gamma > 0}; larger values follow the
#'   data more closely. Default \code{0.1}.
#' @param opt.meth,diff.meth Ignored. Kept so that older calls keep working;
#'   the native solver needs no optimizer choice.
#' @param kernel Kernel: \code{"add_kernel"} (default), \code{"lin_kernel"},
#'   \code{"rbf_kernel"} or \code{"poly_kernel"}.
#' @param kernel.pars Optional kernel parameters, see the Kernels section.
#' @param resid.folds Number of folds for the cross-fitted residuals (default
#'   \code{5}); \code{1} uses in-sample residuals. Folds interleave the ordered
#'   survival times, so fitting does not touch the random number generator.
#'
#' @return An object of class `mlsurv_model`, a named list with elements
#'   \code{model} (the fitted engine), \code{learner} (\code{"survsvm"}),
#'   \code{formula}, \code{data}, \code{time} and \code{status}; the engine
#'   name \code{"survalis"} is stored in \code{attr(, "engine")}.
#'
#' @seealso `predict_survsvm()`, `tune_survsvm()`
#'
#' @references
#' Shivaswamy PK, Chu W, Jansche M (2007). A support vector approach to
#' censored targets. \emph{Seventh IEEE International Conference on Data
#' Mining}, 655-660. \doi{10.1109/ICDM.2007.93}
#'
#' Van Belle V, Pelckmans K, Van Huffel S, Suykens JAK (2011). Support vector
#' methods for survival analysis: a comparison between ranking and regression
#' approaches. \emph{Artificial Intelligence in Medicine}, 53(2), 107-118.
#' \doi{10.1016/j.artmed.2011.06.006}
#'
#' Daemen A, Timmerman D, Van den Bosch T, et al. (2012). Improved modeling of
#' clinical data with kernel methods. \emph{Artificial Intelligence in
#' Medicine}, 54(2), 103-114. \doi{10.1016/j.artmed.2011.11.001}
#'
#' Chapelle O (2007). Training a support vector machine in the primal.
#' \emph{Neural Computation}, 19(5), 1155-1178.
#' \doi{10.1162/neco.2007.19.5.1155}
#'
#' @examples
#' mod_svm <- fit_survsvm(Surv(time, status) ~ age + celltype + karno,
#'                        data = veteran, kernel = "add_kernel")
#' predict_survsvm(mod_svm, newdata = veteran[1:5, ], times = c(100, 300, 500))
#'
#' @keywords internal
#' @noRd

fit_survsvm <- function(formula, data,
                        type = "regression",
                        gamma.mu = 0.1,
                        opt.meth = NULL,
                        kernel = c("add_kernel", "lin_kernel", "rbf_kernel", "poly_kernel"),
                        diff.meth = NULL,
                        kernel.pars = NULL,
                        resid.folds = 5) {
  if (!identical(type, "regression")) {
    stop("Only type = \"regression\" is supported: the ranking types give a ",
         "risk score, not a survival time, so no survival curve.", call. = FALSE)
  }
  kernel <- match.arg(kernel)
  gamma.mu <- as.numeric(gamma.mu)
  if (length(gamma.mu) != 1L || !is.finite(gamma.mu) || gamma.mu <= 0) {
    stop("`gamma.mu` must be a positive number.", call. = FALSE)
  }
  resid.folds <- as.integer(resid.folds)
  if (length(resid.folds) != 1L || is.na(resid.folds) || resid.folds < 1L) {
    stop("`resid.folds` must be a positive integer.", call. = FALSE)
  }

  data <- as.data.frame(data)
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  surv <- stats::model.response(mf)
  if (!inherits(surv, "Surv") || attr(surv, "type") != "right") {
    stop("`formula` must have a right-censored Surv(time, status) response.", call. = FALSE)
  }
  time <- as.numeric(surv[, 1L])
  delta <- as.numeric(surv[, 2L])
  if (any(time <= 0)) {
    stop("survsvm models log(time), so all survival times must be positive.", call. = FALSE)
  }
  y <- log(time)

  enc <- .ssvm_encoder(stats::terms(mf), mf, kernel, kernel.pars)
  Z <- .ssvm_encode(enc, mf)
  K <- .ssvm_kernel(enc, Z, Z)

  sol <- .ssvm_solve(K, y, delta, gamma.mu)

  # cross-fitted residuals for the error distribution
  n <- length(y)
  if (resid.folds > 1L && n >= 2L * resid.folds) {
    fold <- integer(n)
    fold[order(y, -delta)] <- rep_len(seq_len(resid.folds), n)
    f_out <- numeric(n)
    for (k in seq_len(resid.folds)) {
      tr <- which(fold != k)
      te <- which(fold == k)
      s <- .ssvm_solve(K[tr, tr, drop = FALSE], y[tr], delta[tr], gamma.mu)
      f_out[te] <- drop(K[te, tr[s$active], drop = FALSE] %*% s$beta) + s$b
    }
  } else {
    f_out <- drop(K[, sol$active, drop = FALSE] %*% sol$beta) + sol$b
  }
  resid <- y - f_out
  km <- survival::survfit(survival::Surv(resid, delta) ~ 1)

  model <- list(
    encoder = enc,
    sv = Z[sol$active, , drop = FALSE],
    beta = sol$beta,
    b = sol$b,
    gamma = gamma.mu,
    resid = resid,
    status = delta,
    resid_km = list(time = km$time, surv = km$surv),
    iterations = sol$iterations
  )

  structure(list(
    model = model,
    learner = "survsvm",
    formula = formula,
    data = data,
    time = all.vars(formula)[[2]],
    status = all.vars(formula)[[3]]
  ), class = "mlsurv_model", engine = "survalis")
}


#' Predict Survival Probabilities from a Survival SVM
#'
#' Turns the fitted log-time predictor \eqn{f(x)} of a model from
#' `fit_survsvm()` into survival probabilities at the requested times.
#'
#' @param object A fitted `mlsurv_model` returned by `fit_survsvm()`.
#' @param newdata A `data.frame` of new observations with the predictors used
#'   in the formula.
#' @param times Numeric vector of evaluation time points (same scale as the
#'   training times).
#' @param dist How the distribution of the error \eqn{\log T - f(x)} is
#'   obtained:
#'   \describe{
#'     \item{\code{"km"}}{(default) Kaplan-Meier estimate of the cross-fitted
#'       residuals, \eqn{S(t \mid x) = \hat S_\varepsilon(\log t - f(x))}.}
#'     \item{\code{"lognormal"}}{normal errors, mean and scale estimated by
#'       censored maximum likelihood on the cross-fitted residuals.}
#'     \item{\code{"exp"}, \code{"weibull"}}{the earlier fixed mappings with
#'       \eqn{e^{f(x)}} as the time scale: \eqn{\exp(-t / e^{f(x)})} and
#'       \eqn{\exp(-(t / e^{f(x)})^{shape})}.}
#'   }
#' @param shape Weibull shape parameter (only for \code{dist = "weibull"}).
#'
#' @return A base `data.frame` with one row per observation in `newdata` and
#'   columns `t=<time>` holding survival probabilities, non-increasing in time.
#'
#' @seealso `fit_survsvm()`, `tune_survsvm()`
#'
#' @examples
#' mod <- fit_survsvm(Surv(time, status) ~ age + celltype + karno, data = veteran)
#' predict_survsvm(mod, newdata = veteran[1:5, ], times = c(100, 300, 500))
#' predict_survsvm(mod, newdata = veteran[1:5, ], times = c(100, 300, 500),
#'                 dist = "lognormal")
#'
#' @keywords internal
#' @noRd

predict_survsvm <- function(object, newdata, times,
                            dist = c("km", "lognormal", "exp", "weibull"),
                            shape = 1) {
  if (!is.null(object$learner) && object$learner != "survsvm") {
    warning("Object passed to predict_survsvm() may not come from fit_survsvm().")
  }
  dist <- match.arg(dist)
  fit <- object$model
  times <- as.numeric(times)
  stopifnot(length(times) > 0L, all(is.finite(times)), all(times >= 0))

  f <- .ssvm_predict_f(fit, as.data.frame(newdata))
  u <- matrix(log(times), length(f), length(times), byrow = TRUE) - f

  S <- switch(dist,
    km = {
      km <- fit$resid_km
      c(1, km$surv)[findInterval(u, km$time) + 1L]
    },
    lognormal = {
      p <- .ssvm_lognormal(fit$resid, fit$status)
      stats::pnorm((u - p[1L]) / p[2L], lower.tail = FALSE)
    },
    exp = exp(-exp(u)),
    weibull = exp(-exp(shape * u))
  )
  .finalize_survmat(matrix(S, nrow = length(f)), times = times)
}


#' Tune Survival SVM Hyperparameters (Cross-Validation)
#'
#' Cross-validates Survival SVM models over a user-specified grid and selects
#' the best configuration based on the chosen metric. Optionally refits the
#' best model on the full dataset.
#'
#' @param formula A survival formula of the form `Surv(time, status) ~ predictors`.
#'   The left-hand side must be a `Surv()` object from the `survival` package.
#' @param data A `data.frame` containing all variables referenced in `formula`.
#' @param times Numeric vector of evaluation time points (same scale as the
#'   survival time used at training). Must be non-negative and finite.
#' @param metrics Character vector of metrics to evaluate/optimize
#'   (e.g., `"cindex"`, `"ibs"`). The first entry is the primary selection metric.
#' @param param_grid A named list of candidate hyperparameters with entries
#'   `gamma.mu` and `kernel`.
#' @param folds Number of cross-validation folds.
#' @param seed Integer random seed for reproducibility.
#' @param ncores Integer number of CPU cores passed to \code{cv_survlearner}
#'   for fold evaluation (default \code{1}).
#' @param refit_best Logical; if `TRUE`, refits the best configuration on the
#'   full data and returns it as the result (with tuning attributes).
#' @param dist,shape Passed to `predict_survsvm()` during cross-validation.
#'
#' @return If `refit_best = FALSE`, a `data.frame` (class `"tune_surv"`) of grid
#' results with metric columns and tuning parameters. If `refit_best = TRUE`, a
#' fitted `mlsurv_model` (class augmented with `"tune_surv"`) with the full
#' results attached in `attr(, "tuning_results")`.
#'
#' @details
#' Internally calls `cv_survlearner()` with `fit_survsvm()`/`predict_survsvm()`
#' to keep code paths consistent with production usage. A grid point that
#' fails is skipped with a warning.
#'
#' @seealso `fit_survsvm()`, `predict_survsvm()`
#'
#' @examples
#' \donttest{
#' grid <- list(
#'   gamma.mu = c(0.01, 0.1),
#'   kernel = c("lin_kernel", "add_kernel")
#' )
#'
#' res_svm <- tune_survsvm(
#'   formula = Surv(time, status) ~ age + celltype + karno,
#'   data = veteran,
#'   times = c(100, 300, 500),
#'   metrics = c("cindex", "ibs"),
#'   param_grid = grid,
#'   folds = 3
#' )
#' res_svm
#' }
#'
#' @keywords internal
#' @noRd

tune_survsvm <- function(formula, data, times,
                         metrics = "cindex",
                         param_grid = list(gamma.mu = c(0.01, 0.1, 1),
                                           kernel = c("lin_kernel", "add_kernel")),
                         folds = 5, seed = 42, ncores = 1,
                         refit_best = FALSE,
                         dist = "km", shape = 1) {
  stopifnot(is.list(param_grid))
  param_df <- do.call(expand.grid, c(param_grid, list(stringsAsFactors = FALSE)))

  results <- .pmap_rbind_dt(param_df, function(gamma.mu, kernel) {
    tryCatch({
      res_cv <- cv_survlearner(
        formula = formula,
        data = data,
        fit_fun = fit_survsvm,
        pred_fun = function(...) predict_survsvm(..., dist = dist, shape = shape),
        times = times,
        metrics = metrics,
        folds = folds,
        seed = seed,
        ncores = ncores,
        gamma.mu = gamma.mu,
        kernel = kernel
      )
      .wide_metric_row(list(gamma.mu = gamma.mu, kernel = kernel), cv_summary(res_cv))
    }, error = function(e) {
      warning(glue("Skipping (gamma.mu={gamma.mu}, kernel={kernel}): {conditionMessage(e)}"),
              call. = FALSE)
      NULL
    })
  })

  if (is.null(results) || nrow(results) == 0) {
    warning("All tuning combinations failed.")
    return(NULL)
  }

  maximize <- metrics[1] %in% c("cindex", "auc", "accuracy")
  results <- .arrange_by_metric_dt(results, metrics[1], maximize)

  if (refit_best) {
    model <- fit_survsvm(
      formula = formula,
      data = data,
      gamma.mu = results$gamma.mu[1],
      kernel = results$kernel[1]
    )
    attr(model, "tuning_results") <- results
    class(model) <- c("tune_surv", class(model))
    return(model)
  }

  class(results) <- c("tune_surv", class(results))
  results
}


# ---- engine ---------------------------------------------------------------

# Exact minimizer of
#   J(beta, b) = 1/2 beta' K beta + gamma/2 sum_{i in A(f)} (y_i - f_i)^2,
#   f = K beta + b,  A(f) = {delta_i = 1} U {delta_i = 0, y_i > f_i},
# by active-set Newton steps. At the optimum beta_i = gamma (y_i - f_i) on A,
# zero elsewhere, and sum(beta) = 0. Returns beta on the support only.
.ssvm_solve <- function(K, y, delta, gamma, maxit = 100L) {
  n <- length(y)
  event <- delta == 1

  objective <- function(beta, b) {
    Kb <- drop(K %*% beta)
    r <- y - Kb - b
    0.5 * sum(beta * Kb) + 0.5 * gamma * sum(r[event | r > 0]^2)
  }
  # LS-SVM solution restricted to A
  newton <- function(A) {
    H <- K[A, A, drop = FALSE]
    diag(H) <- diag(H) + 1 / gamma
    R <- chol(H)
    rhs <- backsolve(R, forwardsolve(t(R), cbind(1, y[A])))
    b <- sum(rhs[, 2L]) / sum(rhs[, 1L])
    beta <- numeric(n)
    beta[A] <- rhs[, 2L] - b * rhs[, 1L]
    list(beta = beta, b = b)
  }

  # start from the LS-SVM that treats every time as an event
  A_exact <- seq_len(n)
  cur <- newton(A_exact)
  obj <- objective(cur$beta, cur$b)
  it <- 0L
  repeat {
    f <- drop(K %*% cur$beta) + cur$b
    A <- which(event | y > f)
    if (identical(A, A_exact)) break         # Newton point with a stable active set
    it <- it + 1L
    if (it > maxit) {
      warning("survsvm solver reached the iteration limit.", call. = FALSE)
      break
    }
    cand <- newton(A)
    step <- 1
    repeat {                                   # backtracking keeps J decreasing
      beta <- cur$beta + step * (cand$beta - cur$beta)
      b <- cur$b + step * (cand$b - cur$b)
      o <- objective(beta, b)
      if (o <= obj + 1e-12 * max(1, abs(obj)) || step < 1e-10) break
      step <- step / 2
    }
    cur <- list(beta = beta, b = b)
    obj <- o
    A_exact <- if (step == 1) A else NULL
  }

  active <- which(cur$beta != 0)
  list(beta = cur$beta[active], b = cur$b, active = active, iterations = it)
}

# Everything needed to encode new data exactly as the training data.
.ssvm_encoder <- function(tt, mf, kernel, kernel.pars) {
  tt <- stats::delete.response(tt)
  enc <- list(terms = tt, xlevels = stats::.getXlevels(tt, mf), kernel = kernel)
  if (kernel == "add_kernel") {
    if (any(attr(tt, "order") > 1L)) {
      stop("add_kernel works on single predictors; drop interaction terms.", call. = FALSE)
    }
    vars <- attr(tt, "term.labels")
    X <- mf[vars]
    is_cat <- vapply(X, function(x) (is.factor(x) && !is.ordered(x)) ||
                       is.character(x) || is.logical(x), logical(1))
    enc$vars <- vars
    enc$is_cat <- unname(is_cat)
    enc$levels <- lapply(X, function(x) {
      if (is.factor(x)) levels(x) else if (is.character(x) || is.logical(x))
        sort(unique(as.character(x))) else NULL
    })
    enc$range <- vapply(seq_along(X), function(j) {
      if (is_cat[[j]]) NA_real_ else diff(range(.ssvm_num(X[[j]])))
    }, numeric(1))
  } else {
    X <- .bnn_design(tt, mf, xlev = enc$xlevels)
    enc$columns <- colnames(X)
    enc$center <- colMeans(X)
    s <- sqrt(colSums(sweep(X, 2L, enc$center)^2) / max(nrow(X) - 1L, 1L))
    enc$scale <- ifelse(s > 0, s, 1)
    enc$pars <- if (!is.null(kernel.pars)) as.numeric(kernel.pars) else
      switch(kernel, rbf_kernel = ncol(X), poly_kernel = c(1, 2), NULL)
  }
  enc
}

.ssvm_num <- function(x) if (is.factor(x)) as.numeric(unclass(x)) else as.numeric(x)

.ssvm_encode <- function(enc, mf) {
  if (enc$kernel == "add_kernel") {
    Z <- vapply(seq_along(enc$vars), function(j) {
      x <- mf[[enc$vars[j]]]
      if (!enc$is_cat[[j]]) return(.ssvm_num(x))
      code <- match(as.character(x), enc$levels[[j]])
      if (anyNA(code[!is.na(x)])) {
        stop("`", enc$vars[j], "` has levels not seen in training.", call. = FALSE)
      }
      as.numeric(code)
    }, numeric(nrow(mf)))
    return(matrix(Z, nrow = nrow(mf), dimnames = list(NULL, enc$vars)))
  }
  X <- .bnn_design(enc$terms, mf, xlev = enc$xlevels)
  if (!identical(colnames(X), enc$columns)) {
    stop("Predictors in `newdata` do not match the training design.", call. = FALSE)
  }
  sweep(sweep(X, 2L, enc$center), 2L, enc$scale, "/")
}

.ssvm_kernel <- function(enc, Z1, Z2) {
  p <- max(ncol(Z1), 1L)
  switch(enc$kernel,
    add_kernel = {
      K <- matrix(0, nrow(Z1), nrow(Z2))
      for (j in seq_len(ncol(Z1))) {
        r <- enc$range[[j]]
        K <- K + if (enc$is_cat[[j]]) {
          outer(Z1[, j], Z2[, j], "==")
        } else if (r > 0) {
          pmax(0, 1 - abs(outer(Z1[, j], Z2[, j], "-")) / r)
        } else {
          1
        }
      }
      K / p
    },
    lin_kernel = tcrossprod(Z1, Z2) / p,
    rbf_kernel = {
      D2 <- outer(rowSums(Z1^2), rowSums(Z2^2), "+") - 2 * tcrossprod(Z1, Z2)
      exp(-pmax(D2, 0) / (2 * enc$pars[1L]))
    },
    poly_kernel = (tcrossprod(Z1, Z2) / p + enc$pars[1L])^enc$pars[2L]
  )
}

# f(x) for new data.
.ssvm_predict_f <- function(fit, newdata) {
  enc <- fit$encoder
  mf <- stats::model.frame(enc$terms, data = newdata, xlev = enc$xlevels,
                           na.action = stats::na.pass)
  Z <- .ssvm_encode(enc, mf)
  if (anyNA(Z)) stop("`newdata` has missing predictor values.", call. = FALSE)
  drop(.ssvm_kernel(enc, Z, fit$sv) %*% fit$beta) + fit$b
}

# Censored normal fit (mean, sd) of the residuals.
.ssvm_lognormal <- function(resid, status) {
  rf <- survival::survreg(survival::Surv(exp(resid), status) ~ 1, dist = "lognormal")
  c(unname(stats::coef(rf)), rf$scale)
}
