test_that("exact_design_efficiency_maximin returns min <= 1 for D and A", {
  theta <- c(p0 = 0.07, p1 = 0.93, p2 = 0.96)
  u <- 1:20L
  q_cost <- 0.2
  f <- gt_huang2020_regressor(theta, q_cost)
  res_d <- calc_Dopt(u, f, drop_tol = 1e-5)
  res_a <- calc_Aopt(u, f, drop_tol = 1e-5)
  loss_ref <- list(D = -res_d$value, A = res_a$value)
  mm <- compute_maximin_design(u, f, loss_ref, criteria = c("D", "A"), drop_tol = 1e-5)
  out <- round_gt_design_budget(
    mm,
    u,
    theta,
    C = 50,
    q_cost = q_cost,
    criterion = "D"
  )
  ed <- exact_design_efficiency_maximin(out$M_exact, loss_ref, c("D", "A"))
  expect_true(length(ed$efficiencies) >= 2L)
  expect_true(is.finite(ed$min_efficiency))
  expect_true(ed$min_efficiency > 0 && ed$min_efficiency <= 1.01)
})

test_that("round_gt_design_budget_maximin reproduces DD-AA Table 4 rows", {
  theta <- c(p0 = 0.07, p1 = 0.93, p2 = 0.96)
  u <- 1:61L
  q_cost <- 0.2
  f <- gt_huang2020_regressor(theta, q_cost)

  res_d <- calc_Dopt(u, f, drop_tol = 1e-6)
  res_a <- calc_Aopt(u, f, drop_tol = 1e-6)
  loss_ref <- list(D = -res_d$value, A = res_a$value)
  mm <- compute_maximin_design(u, f, loss_ref, criteria = c("D", "A"), drop_tol = 1e-6)

  out100 <- round_gt_design_budget_maximin(mm, u, theta, 100, q_cost, loss_ref, c("D", "A"))
  out500 <- round_gt_design_budget_maximin(mm, u, theta, 500, q_cost, loss_ref, c("D", "A"))

  expect_equal(as.integer(out100$design_exact[1, ]), c(1L, 10L, 59L, 61L))
  expect_equal(as.integer(out100$design_exact[2, ]), c(26L, 8L, 1L, 3L))
  expect_equal(out100$min_efficiency, 0.932, tolerance = 0.01)

  expect_equal(as.integer(out500$design_exact[1, ]), c(1L, 10L, 60L, 61L))
  expect_equal(as.integer(out500$design_exact[2, ]), c(130L, 44L, 1L, 18L))
  expect_equal(out500$min_efficiency, 0.948, tolerance = 0.01)
})
