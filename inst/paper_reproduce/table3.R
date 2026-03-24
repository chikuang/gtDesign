# Reproduce Table 3 (paper): single-objective exact designs via budget rounding
# ------------------------------------------------------------------------------
# Title in paper: M = 150, theta* = (0.07, 0.93, 0.96), q = 0.2; criteria D, A,
# D_s, c with c_vec = (0, 1, 1) for c-optimality.
# C in {100, 500, 10000}.
#
# IMPORTANT: Table 3 matches MATLAB / `GT_Single_budget`-style rounding where
# zero-fix is NOT applied after floor (`fix_zero_floor = FALSE`). With the
# default `fix_zero_floor = TRUE`, floor counts (and C_r, Delta) differ.
#
# phi(I^-1) in the paper:  for D use det(M)^{-1}; for A use tr(M^{-1}); for D_s
# and c use c' M^{-1} c with the same contrasts as in the optimizations.
#
# Run:
#   Rscript -e 'devtools::load_all(); source("inst/paper_reproduce/table3.R")'

options(scipen = 999)

if (!exists("round_gt_design_budget", mode = "function")) {
  stop("Load package first: devtools::load_all()", call. = FALSE)
}

theta <- c(p0 = 0.07, p1 = 0.93, p2 = 0.96)
u <- seq_len(150L)
q_cost <- 0.2
C_values <- c(100, 500, 10000)
support_tol <- 1e-6

f <- gt_huang2020_regressor(theta, q_cost)
opts <- list(
  cVec_c = c(0, 1, 1),
  cVec_Ds = c(1, 0, 0)
)

res_D <- calc_Dopt(u, f, drop_tol = support_tol)
res_A <- calc_Aopt(u, f, drop_tol = support_tol)
res_Ds <- calc_copt(u, f, cVec = opts$cVec_Ds, drop_tol = support_tol)
res_c <- calc_copt(u, f, cVec = opts$cVec_c, drop_tol = support_tol)

#' Scalar phi(I^-1) as in Table 3
phi_invinfo <- function(M, criterion, opts) {
  if (criterion == "D") {
    return(1 / as.numeric(det(M)))
  }
  Minv <- solve(M)
  if (criterion == "A") {
    return(sum(diag(Minv)))
  }
  if (criterion == "Ds") {
    cvec <- as.numeric(opts$cVec_Ds)
    return(as.numeric(t(cvec) %*% Minv %*% cvec))
  }
  if (criterion == "c") {
    cvec <- as.numeric(opts$cVec_c)
    return(as.numeric(t(cvec) %*% Minv %*% cvec))
  }
  stop("Unknown criterion", call. = FALSE)
}

design_string <- function(dm) {
  paste(paste0(as.integer(dm[1, ]), ":", as.integer(dm[2, ])), collapse = "; ")
}

delta_string <- function(delta) {
  if (!ncol(delta)) {
    return("(none)")
  }
  paste(paste0(as.integer(delta[1, ]), ":", as.integer(delta[2, ])), collapse = "; ")
}

total_cost <- function(dm, q) {
  x <- as.numeric(dm[1, ])
  n <- as.numeric(dm[2, ])
  sum(n * vapply(x, function(xi) gt_huang2020_cost(xi, q), numeric(1)))
}

row_order <- c("D", "A", "Ds", "c")

run_table3_row <- function(approx, criterion) {
  res_list <- vector("list", length(C_values))
  for (i in seq_along(C_values)) {
    C <- C_values[i]
    out <- round_gt_design_budget(
      approx_design = approx,
      u = u,
      theta = theta,
      C = C,
      q_cost = q_cost,
      criterion = criterion,
      opts = opts,
      fix_zero_floor = FALSE,
      repair_floor_budget = TRUE
    )
    C_used <- total_cost(out$design_exact, q_cost)
    Cp <- max(0, round(as.numeric(C - C_used), 4L))
    phi <- phi_invinfo(out$M_exact, criterion, opts)
    res_list[[i]] <- data.frame(
      Criterion = criterion,
      C = C,
      xi_RAC = design_string(out$design_exact),
      C_r = out$C_remaining,
      Delta = delta_string(out$delta),
      C_prime_r = Cp,
      phi_inv = round(as.numeric(phi), 3),
      eff = round(as.numeric(out$efficiency), 3),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, res_list)
}

tab3 <- rbind(
  run_table3_row(res_D, "D"),
  run_table3_row(res_A, "A"),
  run_table3_row(res_Ds, "Ds"),
  run_table3_row(res_c, "c")
)
rownames(tab3) <- NULL

print("Table 3 (paper: q = 0.2, M = 150, fix_zero_floor = FALSE)")
if (requireNamespace("knitr", quietly = TRUE)) {
  print(knitr::kable(tab3, format = "pipe"))
} else {
  print(tab3, row.names = FALSE)
}
