# Quick validation for Supplementary Table S2
# -------------------------------------------
# Focused criteria from S2:
#   - D-Ds
#   - D-Ds-c   (with c = (0,1,1))
#
# Run:
#   Rscript -e 'devtools::load_all(); source("inst/paper_reproduce/tableS2.R")'

options(scipen = 999)

fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
relabel_x_matrix <- function(mat) {
  out <- mat
  if (!is.null(rownames(out)) && length(rownames(out)) >= 1L) {
    rownames(out)[1L] <- "x"
  }
  out
}

if (!exists("calc_Dopt", mode = "function")) {
  stop("Load package first: devtools::load_all()", call. = FALSE)
}

theta <- c(p0 = 0.07, p1 = 0.93, p2 = 0.96)
u <- 1:61L
cvec_ds <- c(1, 0, 0)
cvec_c <- c(0, 1, 1)

make_refs <- function(q_cost, criteria) {
  f <- gt_huang2020_regressor(theta, q_cost)
  res_d <- calc_Dopt(u, f, drop_tol = 1e-6)

  loss_ref <- list(D = -res_d$value)
  opts <- list()

  if ("Ds" %in% criteria) {
    res_ds <- calc_copt(u, f, cVec = cvec_ds, drop_tol = 1e-6)
    loss_ref$Ds <- res_ds$value
    opts$cVec_Ds <- cvec_ds
  }
  if ("c" %in% criteria) {
    res_c <- calc_copt(u, f, cVec = cvec_c, drop_tol = 1e-6)
    loss_ref$c <- res_c$value
    opts$cVec_c <- cvec_c
  }

  list(f = f, loss_ref = loss_ref, opts = opts)
}

run_case <- function(criteria, q_cost, n_vals = NULL, C_vals = NULL) {
  refs <- make_refs(q_cost, criteria)

  mm <- compute_maximin_design(
    u = u,
    f = refs$f,
    loss_ref = refs$loss_ref,
    criteria = criteria,
    opts = refs$opts,
    drop_tol = 1e-6
  )

  print("=======================================================")
  print(paste("S2 block:", paste(criteria, collapse = "-"), "q =", q_cost))
  print(paste("1/t*:", fmt(1 / mm$tstar)))
  print(mm$design)

  if (!is.null(n_vals)) {
    print("-- Fixed n (Algorithm I + modified Step II) --")
    for (n in n_vals) {
      out <- round_gt_design_n_maximin(
        approx_design = mm,
        u = u,
        theta = theta,
        n = n,
        q_cost = q_cost,
        loss_ref = refs$loss_ref,
        criteria = criteria,
        opts = refs$opts
      )
      print(paste("n =", n, "MinEff =", fmt(out$min_efficiency)))
      print(relabel_x_matrix(out$design_exact))
      print("Delta matrix:")
      print(relabel_x_matrix(out$delta))
      print(paste("Eff by criterion:", paste(names(out$efficiencies), fmt(out$efficiencies), collapse = ", ")))
    }
  }

  if (!is.null(C_vals)) {
    print("-- Fixed C (Algorithm II + modified Step II) --")
    for (C in C_vals) {
      out <- round_gt_design_budget_maximin(
        approx_design = mm,
        u = u,
        theta = theta,
        C = C,
        q_cost = q_cost,
        loss_ref = refs$loss_ref,
        criteria = criteria,
        opts = refs$opts
      )
      print(paste("C =", C, "MinEff =", fmt(out$min_efficiency), "C_remaining =", fmt(out$C_remaining)))
      print(relabel_x_matrix(out$design_exact))
      print("Delta matrix:")
      print(relabel_x_matrix(out$delta))
      print(paste("Eff by criterion:", paste(names(out$efficiencies), fmt(out$efficiencies), collapse = ", ")))
    }
  }
}

# Table S2 rows
run_case(criteria = c("D", "Ds"), q_cost = 0.0, n_vals = c(10L, 25L, 50L))
run_case(criteria = c("D", "Ds", "c"), q_cost = 0.0, n_vals = c(10L, 25L, 50L))
run_case(criteria = c("D", "Ds"), q_cost = 0.2, C_vals = c(100L, 500L))
run_case(criteria = c("D", "Ds", "c"), q_cost = 0.2, C_vals = c(100L, 500L))

print("Done: Supplementary Table S2 quick check.")
