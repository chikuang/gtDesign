# Example: E-optimal approximate design (Huang 2020, M = 61, q = 0)
# --------------------------------------------------------------------
# Run:
#   Rscript -e 'devtools::load_all(); source("inst/examples/example_calc_Eopt.R")'

options(scipen = 999)

if (!exists("calc_Eopt", mode = "function")) {
  stop("Load package first: devtools::load_all()", call. = FALSE)
}

theta <- c(p0 = 0.07, p1 = 0.93, p2 = 0.96)
M <- 61L
u <- seq_len(M)
f <- gt_huang2020_regressor(theta, q = 0)

res_E <- calc_Eopt(u, f, drop_tol = 1e-6)

print("E-optimal approximate design (pool size, weight):")
print(res_E$design)

print(paste("Loss value = -lambda_min(M):", format(res_E$value, digits = 8)))
print(paste("lambda_min(M) at optimum:", format(-res_E$value, digits = 8)))

print(paste("Solver status:", res_E$status))




# Equivalence theorem -----------------------------------------------------

library(gtDesign)
theta <- c(p0 = 0.07, p1 = 0.93, p2 = 0.96)
u <- seq_len(150L) #
f <- gt_huang2020_regressor(theta, q = 0.5)

res <- calc_Eopt(u, f, drop_tol = 1e-4)

res$design |> round(3)
eq <- check_equivalence(res, f, tol = 1e-4)

eq$all_nonpositive       # should be TRUE（max d_E ≤ tol）
eq$support_equal_zero    # should be true
eq$max_violation         # must be small

## plot out d_E
plot(eq$candidate_points, eq$directional_derivative, type = "l",
     xlab = "pool size", ylab = expression(d[E]))
abline(h = 0, lty = 2)
points(eq$support_points, eq$support_values, pch = 16, col = "red")
