# Reproduce Table 6 (fixed individuals) in a minimal script
# ----------------------------------------------------------
# Run:
#   Rscript -e 'devtools::load_all(); source("inst/paper_reproduce/table6.R")'

options(scipen = 999)
if (!exists("calc_Dopt", mode = "function")) {
  stop("Load package first: devtools::load_all()", call. = FALSE)
}

theta <- c(p0 = 0.022, p1 = 0.92, p2 = 0.965)
u <- seq_len(15L)
n_subject_allow <- 6632
q_cost <- 0
support_tol <- 1e-5
f <- gt_huang2020_regressor(theta, q_cost)

design_to_string <- function(des_mat) {
  paste(paste0(as.integer(des_mat[1, ]), ":", as.integer(des_mat[2, ])), collapse = ", ")
}

# Helper for one Table-6 row.
# Input:
#   - criterion: one of "D", "A", "E", "c", "Ds"
#   - cvec: contrast vector only used for c / Ds criteria
# What it does:
#   1) solve the approximate optimal design under this criterion;
#   2) call round_gt_design_subject_budget() to enforce subject constraint
#      (n_subject_allow) by converting it to test budget internally;
#   3) format one row for final Table-6-like output.
# Output:
#   A one-row data.frame with criterion, exact design string, tests used,
#   individuals used, and efficiency.
run_one <- function(criterion, cvec = NULL) {
  approx <- switch(
    criterion,
    D = calc_Dopt(u, f, drop_tol = support_tol),
    A = calc_Aopt(u, f, drop_tol = support_tol),
    E = calc_Eopt(u, f, drop_tol = support_tol),
    c = calc_copt(u, f, cVec = cvec, drop_tol = support_tol),
    Ds = calc_copt(u, f, cVec = cvec, drop_tol = support_tol)
  )

  rounded <- round_gt_design_subject_budget(
    approx_design = approx,
    u = u,
    theta = theta,
    n_subject_allow = n_subject_allow,
    q_cost = q_cost,
    criterion = criterion,
    opts = list(
      cVec_c = c(0, 1, 1),
      cVec_Ds = c(1, 0, 0)
    )
  )
  de <- rounded$design_exact

  data.frame(
    Criterion = criterion,
    `Exact design (pool:count)` = design_to_string(de),
    `# Tests (C)` = rounded$C_tests,
    `# Individuals` = rounded$n_subjects_used,
    Efficiency = as.numeric(formatC(rounded$efficiency, format = "f", digits = 3)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# Create Table 6
table6 <- rbind(
  run_one("D"),
  run_one("A"),
  run_one("E"),
  run_one("c", c(0, 1, 1)),
  run_one("Ds", c(1, 0, 0))
)

print("Table 6 final output:")
if (requireNamespace("knitr", quietly = TRUE)) {
  print(knitr::kable(table6, format = "pipe"))
} else {
  print(table6, row.names = FALSE)
}
