test_that("E-optimal directional derivative matches (v'f)^2 - lambda_min when r*=1", {
  M <- matrix(c(3, 0, 0, 5), 2L, 2L)
  ev <- eigen(M, symmetric = TRUE)
  lambda_min <- min(ev$values)
  v <- ev$vectors[, which.min(ev$values), drop = FALSE]
  f <- function(x) {
    if (x == 1) {
      c(1, 0)
    } else {
      c(0, 1)
    }
  }
  u <- c(1, 2)
  dd <- calc_directional_derivatives(u, M, f, criteria = "E")
  expect_named(dd, "dE")
  fx1 <- f(1)
  fx2 <- f(2)
  d1 <- as.numeric((t(v) %*% fx1)^2 - lambda_min)
  d2 <- as.numeric((t(v) %*% fx2)^2 - lambda_min)
  expect_equal(dd$dE[1], d1, tolerance = 1e-10)
  expect_equal(dd$dE[2], d2, tolerance = 1e-10)
})

test_that("E-only path does not require invertible M", {
  M <- matrix(c(1, 1, 1, 1), 2L, 2L)
  f <- function(x) c(1, x)
  u <- c(0, 1)
  dd <- calc_directional_derivatives(
    u, M, f,
    criteria = "E",
    use_ginv = FALSE
  )
  expect_length(dd$dE, 2L)
  expect_true(all(is.finite(dd$dE)))
})
