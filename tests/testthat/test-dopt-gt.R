test_that("D-optimal design returns valid weights (generic regression)", {
  quad_reg <- function(x) c(1, x, x^2)
  u <- seq(-1, 1, length.out = 21)

  res <- calc_Dopt(u, quad_reg)

  expect_s3_class(res, "gt_design")
  expect_equal(sum(res$weights), 1, tolerance = 1e-6)
  expect_true(all(res$weights >= 0))
})

test_that("D-opt on one-parameter pool-size information uses gt_design class", {
  u <- 1:20
  p <- 0.08
  fisher_sqrt <- function(m) {
    m <- as.integer(m)
    q <- 1 - (1 - p)^m
    dq <- m * (1 - p)^(m - 1)
    sqrt(dq^2 / (q * (1 - q)))
  }
  f <- function(x) fisher_sqrt(as.integer(x))
  res <- calc_Dopt(u, f)
  expect_s3_class(res, "gt_d_design")
  expect_equal(sum(res$weights), 1, tolerance = 1e-5)
})
