test_that("compute_maximin_design accepts D reference loss from calc_Dopt (can be negative)", {
  theta <- c(p0 = 0.07, p1 = 0.93, p2 = 0.96)
  u <- 1:15L
  f <- gt_huang2020_regressor(theta, q = 0)
  res_d <- calc_Dopt(u, f, drop_tol = 1e-6)
  res_a <- calc_Aopt(u, f, drop_tol = 1e-6)
  loss_ref <- list(D = -res_d$value, A = res_a$value)
  res <- compute_maximin_design(u, f, loss_ref, criteria = c("D", "A"))
  expect_true(length(res$weights) == length(u))
  expect_true(is.finite(res$tstar))
})
