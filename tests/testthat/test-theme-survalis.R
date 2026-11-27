test_that("theme_survalis() and scale helpers produce valid ggplot components", {
  testthat::skip_if_not_installed("ggplot2")

  th <- theme_survalis()
  expect_s3_class(th, "theme")

  p <- ggplot2::ggplot(veteran, ggplot2::aes(age, karno, color = celltype)) +
    ggplot2::geom_point() +
    scale_color_survalis() +
    theme_survalis()

  built <- ggplot2::ggplot_build(p)
  expect_s3_class(built, "ggplot_built")
})

test_that(".survalis_pal() recycles/interpolates beyond the base palette length", {
  expect_length(.survalis_pal(3), 3)
  expect_length(.survalis_pal(8), 8)
  expect_length(.survalis_pal(12), 12)
  expect_true(all(grepl("^#", .survalis_pal(12))))
})

test_that("single-time ALE and ICE-only plots build (one discrete colour)", {
  m <- fit(Surv(time, status) ~ age + karno, survival::veteran, model = "coxph")
  p1 <- plot(interpret(m, "ale", feature = "karno", times = 100))
  p2 <- plot(interpret(m, "ice", feature = "karno", times = c(50, 100)))
  for (p in list(p1, p2)) expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
})
