test_that("rounding_run_size_combinations matches small enumeration", {
  r <- rounding_run_size_combinations(2L, 3L)
  expect_true(nrow(r) > 0L)
  expect_equal(ncol(r), 4L)
  expect_true(all(rowSums(r[, 1:2, drop = FALSE]) <= 3))
})

test_that("rounding_budget_combinations returns tight rows", {
  r <- rounding_budget_combinations(c(1, 2), 3)
  expect_true(nrow(r) > 0L)
  expect_equal(ncol(r), 4L)
})

test_that("round_gt_design_budget returns finite efficiency (D)", {
  theta <- c(p0 = 0.07, p1 = 0.93, p2 = 0.96)
  u <- 1:35L
  q_cost <- 0.2
  f <- gt_huang2020_regressor(theta, q_cost)
  res <- calc_Dopt(u, f, drop_tol = 1e-6)
  out <- round_gt_design_budget(
    res,
    u,
    theta,
    C = 75,
    q_cost = q_cost,
    criterion = "D"
  )
  expect_true(is.finite(out$efficiency))
  expect_true(out$efficiency > 0 && out$efficiency <= 1.01)
})

test_that("round_gt_design_budget returns finite efficiency (E)", {
  skip_if_not_installed("CVXR")
  theta <- c(p0 = 0.07, p1 = 0.93, p2 = 0.96)
  u <- 1:35L
  q_cost <- 0.2
  f <- gt_huang2020_regressor(theta, q_cost)
  res <- calc_Eopt(u, f, drop_tol = 1e-6)
  out <- round_gt_design_budget(
    res,
    u,
    theta,
    C = 75,
    q_cost = q_cost,
    criterion = "E"
  )
  expect_true(is.finite(out$efficiency))
  expect_true(out$efficiency > 0 && out$efficiency <= 1.01)
})
