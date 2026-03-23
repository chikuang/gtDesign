test_that("Huang et al. (2020) D-opt matches arXiv:2508.08445 Table 1 (q=0, M=61)", {
  theta <- c(p0 = 0.07, p1 = 0.93, p2 = 0.96)
  u <- 1:61
  f <- gt_huang2020_regressor(theta, q = 0)
  res <- calc_Dopt(u, f, drop_tol = 1e-6)

  sup <- sort(res$design$point)
  expect_equal(sup, c(1, 17, 61))
  expect_equal(res$design$weight, rep(1 / 3, 3), tolerance = 1e-4)
})

test_that("gt_huang2020_pi matches explicit formula", {
  theta <- c(p0 = 0.1, p1 = 0.9, p2 = 0.95)
  x <- 5L
  z <- (1 - theta["p0"])^x
  pi_manual <- theta["p1"] - (theta["p1"] + theta["p2"] - 1) * z
  expect_equal(gt_huang2020_pi(x, theta), as.numeric(pi_manual))
})
