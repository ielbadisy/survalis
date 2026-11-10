## Inputs that used to be accepted silently and give wrong results, or that
## failed only after a full resampling run.

test_that("a computed time or status in Surv() is rejected, not misread", {
  expect_error(fit(Surv(time / 365, status) ~ age, veteran, "coxph"),
               "column names, not expressions")
  expect_error(fit(Surv(time, status > 0) ~ age, veteran, "coxph"),
               "column names, not expressions")
  expect_error(fit(Surv(nope, status) ~ age, veteran, "coxph"), "not in `data`")
  expect_s3_class(fit(Surv(time, status == 1) ~ age, veteran, "coxph"),
                  "survalis_fit")
})

test_that("metric requests are validated before any model is fitted", {
  f <- Surv(time, status) ~ age + karno
  expect_error(evaluate(f, veteran, model = "coxph", metrics = "foo"),
               "Unknown metric")
  expect_error(evaluate(f, veteran, model = "coxph", times = c(50, 100),
                        metrics = c("cindex", "ece")), "single evaluation time")
  expect_error(compare(f, veteran, models = c("coxph", "rpart"),
                       times = c(50, 100), metrics = "brier"),
               "single evaluation time")
  expect_error(tune(f, veteran, model = "rpart", times = c(50, 100),
                    metric = "brier"), "single evaluation time")
  m <- fit(f, veteran, "coxph")
  expect_error(evaluate(m, newdata = veteran, times = c(50, 100),
                        metrics = "ece"), "single evaluation time")
})

test_that("metric values are stored at full precision", {
  f <- Surv(time, status) ~ age + karno
  e <- evaluate(f, veteran, model = "coxph", times = c(50, 100, 200),
                resampling = cv(3, seed = 1), metrics = c("cindex", "ibs"))
  expect_true(any(e$per_split$value != round(e$per_split$value, 3)))
  expect_true(any(e$summary$mean != round(e$summary$mean, 3)))
  s <- cv_summary(data.frame(metric = "ibs", value = c(0.12345, 0.23456)))
  expect_equal(s$mean, mean(c(0.12345, 0.23456)))
})

test_that("explanation methods use only the model's predictors", {
  g <- survival::gbsg
  g$junk <- seq_len(nrow(g))
  m <- fit(Surv(rfstime, status) ~ age + size + nodes, g, "coxph")
  expect_setequal(.model_features(m, names(g)), c("age", "size", "nodes"))
  vi <- as.data.frame(interpret(m, "varimp", times = c(500, 1000)))
  expect_setequal(vi$feature, c("age", "size", "nodes"))
  sh <- as.data.frame(interpret(m, "shap", newdata = g[1, ],
                                times = c(500, 1000)))
  expect_setequal(unique(sh$feature), c("age", "size", "nodes"))
  ts <- as.data.frame(interpret(m, "tree_surrogate", times = c(500, 1000)))
  expect_true(all(ts$feature %in% c("age", "size", "nodes")))
  m2 <- fit(Surv(rfstime, status) ~ ., g[c("rfstime", "status", "age", "nodes")],
            "coxph")
  expect_equal(.model_features(m2, names(g)), c("age", "nodes"))
})
