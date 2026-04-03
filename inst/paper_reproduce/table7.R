# Reproduce Table 7 (fixed assays/tests) in a minimal script
# -----------------------------------------------------------
# Run:
#   Rscript -e 'devtools::load_all(); source("inst/paper_reproduce/table7.R")'

options(scipen = 999)
if (!exists("calc_Dopt", mode = "function")) {
  stop("Load package first: devtools::load_all()", call. = FALSE)
}

theta <- c(p0 = 0.022, p1 = 0.92, p2 = 0.965)
u <- seq_len(15L)
C_tests <- 500
q_cost <- 0
support_tol <- 1e-5
f <- gt_huang2020_regressor(theta, q_cost)

design_to_string <- function(des_mat) {
  paste(paste0(as.integer(des_mat[1, ]), ":", as.integer(des_mat[2, ])), collapse = ", ")
}

# Helper for one Table-7 row under fixed test budget C_tests.
run_one <- function(criterion, cvec = NULL) {
  approx <- switch(
    criterion,
    D = calc_Dopt(u, f, drop_tol = support_tol),
    A = calc_Aopt(u, f, drop_tol = support_tol),
    E = calc_Eopt(u, f, drop_tol = support_tol),
    c = calc_copt(u, f, cVec = cvec, drop_tol = support_tol),
    Ds = calc_copt(u, f, cVec = cvec, drop_tol = support_tol)
  )

  rounded <- round_gt_design_budget(
    approx_design = approx,
    u = u,
    theta = theta,
    C = C_tests,
    q_cost = q_cost,
    criterion = criterion,
    opts = list(
      cVec_c = c(0, 1, 1),
      cVec_Ds = c(1, 0, 0)
    )
  )

  de <- rounded$design_exact
  n_ind <- sum(as.integer(de[1, ]) * as.integer(de[2, ]))

  data.frame(
    Criterion = criterion,
    `Exact design (pool:count)` = design_to_string(de),
    `# Tests (C)` = C_tests,
    `# Individuals` = n_ind,
    Efficiency = as.numeric(formatC(rounded$efficiency, format = "f", digits = 3)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# Match paper row order (E after A).
table7 <- rbind(
  run_one("D"),
  run_one("A"),
  run_one("E"),
  run_one("c", c(0, 1, 1)),
  run_one("Ds", c(1, 0, 0))
)

print("Table 7 final output:")
if (requireNamespace("knitr", quietly = TRUE)) {
  print(knitr::kable(table7, format = "pipe"))
} else {
  print(table7, row.names = FALSE)
}

#
# (1 - q_cost * q_cost*1) * 78
# (1 - q_cost * q_cost*7) * 259
# (1 - q_cost * q_cost*15) * 163
# 78 + 259 * 7+ 163 * 15
