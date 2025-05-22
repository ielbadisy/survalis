fit_cforest <- function(formula, data,
                        teststat = "quad",
                        testtype = "Univariate",
                        mincriterion = 0,
                        ntree = 500,
                        mtry = 5,
                        replace = TRUE,
                        fraction = 0.632,
                        ...) {
  stopifnot(requireNamespace("party", quietly = TRUE))

  model <- party::cforest(
    formula = formula,
    data = data,
    controls = party::cforest_control(
      teststat = teststat,
      testtype = testtype,
      mincriterion = mincriterion,
      ntree = ntree,
      mtry = mtry,
      replace = replace,
      fraction = fraction,
      ...
    )
  )

  structure(
    list(
      model = model,
      learner = "cforest",
      formula = formula,
      data = data
    ),
    class = "mlsurv_model",
    engine = "cforest"
  )
}


predict_cforest <- function(object, newdata, times, ...) {
  stopifnot(object$learner == "cforest")
  stopifnot(requireNamespace("party", quietly = TRUE))

  newdata <- as.data.frame(newdata)

  # get survival predictions from cforest
  survlist <- predict(object$model, newdata = newdata, type = "prob")

  # interpolate at requested times
  prob_matrix <- t(sapply(survlist, function(sfit) {
    s_time <- sfit$time
    s_surv <- sfit$surv

    if (is.null(s_time) || is.null(s_surv) || length(s_surv) == 0) {
      return(rep(NA, length(times)))
    }

    # Linear interpolation with rule = 2 for extrapolation
    approx(x = s_time, y = s_surv, xout = times, method = "linear", rule = 2)$y
  }))

  # label columns and rows as in other learners
  colnames(prob_matrix) <- paste0("t=", times)
  rownames(prob_matrix) <- paste0("ID_", seq_len(nrow(newdata)))

  return(as.data.frame(prob_matrix))
}


library(survival)
data(veteran)

mod_cforest <- fit_cforest(Surv(time, status) ~ age + celltype + karno, data = veteran)
predict_cforest(mod_cforest, newdata = veteran[1:5, ], times = c(100, 200, 300))



# --- Evaluate with CV ---
cv_results_cforest <- cv_survlearner(
  formula = Surv(time, status) ~ age + celltype + karno,
  data = veteran,
  fit_fun = fit_cforest,
  pred_fun = predict_cforest,
  times = c(100, 200, 300),
  metrics = c("cindex", "ibs"),
  folds = 5,
  seed = 42
)

# --- Print Results ---
print(cv_results_cforest)
