test_that("D-optimal design returns valid weights (generic regression)", {
  quad_reg <- function(x) c(1, x, x^2)
  u <- seq(-1, 1, length.out = 21)

  res <- calc_Dopt(u, quad_reg)

  expect_s3_class(res, "gt_design")
  expect_equal(sum(res$weights), 1, tolerance = 1e-6)
  expect_true(all(res$weights >= 0))
})

test_that("homogeneous group-testing Fisher helper is consistent", {
  p <- 0.05
  m <- 10L
  f <- gt_sqrt_fisher_regressor_homogeneous(p)
  fi <- f(m)
  I <- gt_homogeneous_pool_fisher(p, m)
  expect_equal(as.numeric(fi)^2, I, tolerance = 1e-8)
})

test_that("D-opt on pool sizes uses gt_design class", {
  u <- 1:20
  p <- 0.08
  f <- gt_sqrt_fisher_regressor_homogeneous(p)
  res <- calc_Dopt(u, f)
  expect_s3_class(res, "gt_d_design")
  expect_equal(sum(res$weights), 1, tolerance = 1e-5)
})
