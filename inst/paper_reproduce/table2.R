# Reproduce Table 2 (maximin D–A approximate design, M = 150, q = 0)
# --------------------------------------------------------------------
# From README.qmd chunk `table2`; defines theta explicitly (same as Table 1 intro).
#
# Run:
#   Rscript -e 'devtools::load_all(); source("inst/paper_reproduce/table2.R")'

options(scipen = 999)

if (!exists("compute_maximin_design", mode = "function")) {
  stop("Load package first: devtools::load_all()", call. = FALSE)
}

theta <- c(p0 = 0.07, p1 = 0.93, p2 = 0.96)
u <- seq_len(150L)
f_q <- gt_huang2020_regressor(theta, q = 0.0)

res_d_q <- calc_Dopt(u, f_q, drop_tol = 1e-6)
res_a_q <- calc_Aopt(u, f_q, drop_tol = 1e-6)
loss_ref <- list(
  D = -res_d_q$value,
  A = res_a_q$value
)

print("Reference losses (D, A):")
print(loss_ref)

res_da <- compute_maximin_design(
  u = u,
  f = f_q,
  loss_ref = loss_ref,
  criteria = c("D", "A")
)

print("Maximin (D, A) approximate design:")
print(res_da$design)

print("Efficiencies:")
print(res_da$efficiency)

print(paste("tstar =", res_da$tstar))

tol <- 1e-4
eq_eff <- abs(res_da$efficiency["D"] - res_da$efficiency["A"])
eq_t <- abs(min(res_da$efficiency) - 1 / res_da$tstar)
print(paste(
  "Numerical check (eff balance & 1/t*): eq_eff < tol && eq_t < tol =",
  eq_eff < tol && eq_t < tol
))
