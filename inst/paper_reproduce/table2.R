# Reproduce Table 2 (maximin approximate designs, M = 150, q = 0)
# -----------------------------------------------------------------
# Paper Table 2: rows include D–A, D–A–D_s, and (here) D–A–c–E with E-optimality.
# Reference losses match single-objective solves; maximin uses the same CVX
# formulation as MATLAB `compute_maximin_design` (E: -lambda_min(M) <= loss_ref.E/tstar).
#
# Run:
#   Rscript -e 'devtools::load_all(); source("inst/paper_reproduce/table2.R")'

options(scipen = 999)

if (!exists("compute_maximin_design", mode = "function")) {
  stop("Load package first: devtools::load_all()", call. = FALSE)
}

theta <- c(p0 = 0.07, p1 = 0.93, p2 = 0.96)
u <- seq_len(61L)
q_val <- 0.8
f_q <- gt_huang2020_regressor(theta, q = q_val)
drop_tol <- 1e-6

# --- Reference designs (single-objective) ---
res_d <- calc_Dopt(u, f_q, drop_tol = drop_tol)
res_a <- calc_Aopt(u, f_q, drop_tol = drop_tol)
res_ds <- calc_copt(u, f_q, cVec = c(1, 0, 0), drop_tol = drop_tol)
res_c <- calc_copt(u, f_q, cVec = c(0, 1, 1), drop_tol = drop_tol)
res_e <- calc_Eopt(u, f_q, drop_tol = drop_tol)

loss_ref_da <- list(
  D = -res_d$value,
  A = res_a$value
)

loss_ref_dads <- list(
  D = -res_d$value,
  A = res_a$value,
  Ds = res_ds$value
)

loss_ref_dace <- list(
  D = -res_d$value,
  A = res_a$value,
  c = res_c$value,
  E = res_e$value
)

opts_ds <- list(cVec_Ds = c(1, 0, 0))
opts_c <- list(cVec_c = c(0, 1, 1))

# =============================================================================
# (1) D–A (paper Table 2, first block, q = 0, M = 150)
# =============================================================================
print("=== Table 2: maximin (D, A), M = 150, q = 0 ===")
print("Reference losses (D, A):")
print(loss_ref_da)

res_da <- compute_maximin_design(
  u = u,
  f = f_q,
  loss_ref = loss_ref_da,
  criteria = c("D", "A")
)

print("Design:")
print(res_da$design)
print("Efficiencies:")
print(res_da$efficiency)
print(paste("tstar =", res_da$tstar))

tol <- 1e-4
eq_eff <- abs(res_da$efficiency["D"] - res_da$efficiency["A"])
eq_t <- abs(min(res_da$efficiency) - 1 / res_da$tstar)
print(paste(
  "Numerical check (eff balance & 1/t*):",
  eq_eff < tol && eq_t < tol
))

# =============================================================================
# (2) D–A–D_s (paper Table 2, triple-objective block, q = 0, M = 150)
# =============================================================================
print("")
print("=== Table 2: maximin (D, A, Ds), M = 150, q = 0 ===")
print("Reference losses (D, A, Ds):")
print(loss_ref_dads)

res_dads <- compute_maximin_design(
  u = u,
  f = f_q,
  loss_ref = loss_ref_dads,
  criteria = c("D", "A", "Ds"),
  opts = opts_ds
)

print("Design:")
print(res_dads$design)
print("Efficiencies:")
print(res_dads$efficiency)
print(paste("tstar =", res_dads$tstar))
print(paste(
  "1/t* vs min efficiency:",
  1 / res_dads$tstar,
  min(res_dads$efficiency)
))

# =============================================================================
# (3) D–A–c–E (multi-objective with E-optimality; same q, M)
# =============================================================================
print("")
print("=== Maximin (D, A, c, E), M = 150, q = 0 ===")
print("Reference losses (D, A, c, E); note E is -lambda_min(M) from calc_Eopt:")
print(loss_ref_dace)

res_dace <- compute_maximin_design(
  u = u,
  f = f_q,
  loss_ref = loss_ref_dace,
  criteria = c("D", "A", "c", "E"),
  opts = opts_c
)

print("Design:")
print(res_dace$design)
print("Efficiencies:")
print(res_dace$efficiency |> round(3))
print(paste("tstar =", res_dace$tstar |> round(3)))
print(paste(
  "1/t* vs min efficiency:",
  (1 / res_dace$tstar) |> round(3),
  min(res_dace$efficiency) |> round(3)
))
