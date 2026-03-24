test_that("maximin_design_workflow runs D + A on small grid", {
  theta <- c(p0 = 0.07, p1 = 0.93, p2 = 0.96)
  u <- 1:15L
  f <- gt_huang2020_regressor(theta, q = 0)
  out <- maximin_design_workflow(
    u = u,
    f = f,
    criteria = c("D", "A"),
    drop_tol = 1e-6,
    tol_eta = 1e-3,
    tol_equiv = 1e-2,
    check_equiv = TRUE
  )
  expect_true(is.list(out$loss_ref))
  expect_true(is.finite(out$maximin$tstar))
  expect_true(length(out$eta) == 2L)
  expect_true(out$equivalence$all_nonpositive %in% c(TRUE, FALSE))
})
