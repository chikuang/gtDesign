# Reproduce Table 7: efficient exact designs for 6,632 individuals.
# ---------------------------------------------------------------------------
# Resource cost c(x) = x (q_cost = 1), with no separate assay-count limit.
# Recompute the approximate design per individual before budget rounding;
# scaling an assay-optimal (q = 0) design to a cohort solves a different problem.
#
# With an installed gtDesign package:
#   Rscript inst/paper_reproduce/table7.R
# From a development checkout:
#   Rscript -e 'devtools::load_all(); source("inst/paper_reproduce/table7.R")'
# Optional output directory for CSV results, diagnostics, and session information:
#   Rscript inst/paper_reproduce/table7.R table7_results

options(scipen = 999)
if (!exists("compute_design_SO", mode = "function")) {
  suppressPackageStartupMessages(library(gtDesign))
}

theta <- c(p0 = 0.022, p1 = 0.92, p2 = 0.965)
u <- seq_len(15L)
n_subject_allow <- 6632L
q_cost <- 1
criteria <- c("D", "A", "E", "c", "Ds")
opts <- list(cVec_c = c(0, 1, 1), cVec_Ds = c(1, 0, 0))

# The package includes the cost scaling in the regressor: h(x)h(x)' = H(x)/x.
h <- gt_huang2020_regressor(theta, q_cost)
F <- t(vapply(u, h, numeric(3)))
eig <- eigen(crossprod(F) / length(u), symmetric = TRUE)

design_to_string <- function(des_mat) {
  paste(sprintf("%d (%d)", des_mat[1, ], des_mat[2, ]), collapse = ", ")
}

physical_loss <- function(M, criterion) {
  switch(criterion,
    D = -as.numeric(determinant(M, logarithm = TRUE)$modulus),
    A = sum(diag(solve(M))),
    E = -min(eigen(M, symmetric = TRUE, only.values = TRUE)$values),
    c = as.numeric(crossprod(opts$cVec_c, solve(M, opts$cVec_c))),
    Ds = as.numeric(crossprod(opts$cVec_Ds, solve(M, opts$cVec_Ds)))
  )
}

run_one <- function(criterion) {
  # Improve numerical conditioning without changing the chosen criterion.
  # Arbitrary congruence preserves D-optimality. For c/Ds, transform the
  # contrast as well. A allows common scaling; E allows scaling and rotation.
  T <- if (criterion == "A") {
    10 * diag(3)
  } else if (criterion == "E") {
    100 * t(eig$vectors)
  } else {
    diag(1 / sqrt(eig$values)) %*% t(eig$vectors)
  }
  f_scaled <- local({
    TT <- T
    hh <- h
    function(x) as.numeric(TT %*% hh(x))
  })
  solver_opts <- opts
  if (criterion %in% c("c", "Ds")) {
    key <- if (criterion == "c") "cVec_c" else "cVec_Ds"
    cc <- as.numeric(T %*% opts[[key]])
    solver_opts[[key]] <- cc / sqrt(sum(cc^2))
  }

  approx <- if (criterion == "A") {
    calc_Aopt(
      u, f_scaled, solver = "CLARABEL", drop_tol = 1e-6,
      tol_gap_abs = 1e-9, tol_gap_rel = 1e-9, tol_feas = 1e-9,
      max_iter = 500
    )
  } else {
    compute_design_SO(
      u, f_scaled, criterion = criterion, opts = solver_opts,
      solver = "CLARABEL", support_tol = 1e-6, drop_tol = 1e-10,
      tol_gap_abs = 1e-11, tol_gap_rel = 1e-11, tol_feas = 1e-11,
      max_iter = 500, static_regularization_constant = 1e-12
    )
  }
  stopifnot(approx$status %in% c("optimal", "optimal_inaccurate"))

  # Verify approximate optimality in the ORIGINAL parameter coordinates.
  # A solver status alone is insufficient for these poorly conditioned models.
  M_approx <- crossprod(F, approx$weights * F)
  stopifnot(min(eigen(M_approx, symmetric = TRUE, only.values = TRUE)$values) > 0)
  contrast <- if (criterion == "c") {
    opts$cVec_c
  } else if (criterion == "Ds") {
    opts$cVec_Ds
  } else {
    NULL
  }
  derivatives <- calc_directional_derivatives(
    u, M_approx, h,
    criteria = if (criterion == "Ds") "c" else criterion,
    cVec = contrast, use_ginv = FALSE
  )
  loss_approx <- physical_loss(M_approx, criterion)
  gap_scale <- if (criterion == "D") 3 else abs(loss_approx)
  relative_gap <- max(derivatives[[1]]) / gap_scale
  stopifnot(is.finite(relative_gap), relative_gap < 1e-5)

  # Both approximate optimization and integer rounding use the PERSON budget.
  # Rounding rebuilds the physical information, so use the ORIGINAL contrasts.
  rounded <- round_gt_design_budget(
    approx_design = approx, u = u, theta = theta,
    C = n_subject_allow, q_cost = q_cost, criterion = criterion,
    opts = opts, n_index = 2L
  )
  de <- rounded$design_exact
  stopifnot(
    all(is.finite(de)), all(de[1, ] %in% u),
    all(de[2, ] >= 0), all(de[2, ] == floor(de[2, ])),
    sum(de[1, ] * de[2, ]) == n_subject_allow
  )

  # Independently verify that the exact matrix is TOTAL information / 6632,
  # not information per assay. q = 0 supplies the physical per-assay regressor.
  g <- gt_huang2020_regressor(theta, 0)
  G <- t(vapply(de[1, ], g, numeric(3)))
  J <- crossprod(G, de[2, ] * G)
  M_exact <- J / n_subject_allow
  stopifnot(
    min(eigen(M_exact, symmetric = TRUE, only.values = TRUE)$values) > 0,
    isTRUE(all.equal(M_exact, rounded$M_exact,
                     tolerance = 1e-9, check.attributes = FALSE))
  )
  loss_exact <- physical_loss(M_exact, criterion)
  efficiency <- if (criterion == "D") {
    exp((loss_approx - loss_exact) / 3)
  } else if (criterion == "E") {
    loss_exact / loss_approx
  } else {
    loss_approx / loss_exact
  }
  stopifnot(
    is.finite(efficiency), efficiency > 0.9999, efficiency <= 1 + 1e-5,
    isTRUE(all.equal(efficiency, as.numeric(rounded$efficiency), tolerance = 1e-6))
  )

  list(
    row = data.frame(
      Criterion = criterion,
      `Exact design (pool size and assay count)` = design_to_string(de),
      `# Assays` = sum(de[2, ]),
      `# Individuals` = sum(de[1, ] * de[2, ]),
      `Fixed-person efficiency` = efficiency,
      check.names = FALSE, stringsAsFactors = FALSE
    ),
    diagnostics = data.frame(
      Criterion = criterion, solver_status = approx$status,
      relative_sensitivity_gap = relative_gap,
      loss_approx = loss_approx, loss_exact = loss_exact
    ),
    weights = data.frame(
      Criterion = criterion, pool_size = u,
      individual_budget_proportion = approx$weights
    ),
    counts = data.frame(
      Criterion = criterion, pool_size = as.integer(de[1, ]),
      assays = as.integer(de[2, ]), individuals = as.integer(de[1, ] * de[2, ])
    )
  )
}

table7_results <- setNames(lapply(criteria, run_one), criteria)
table7 <- do.call(rbind, lapply(table7_results, function(x) x$row))
table7_diagnostics <- do.call(rbind, lapply(table7_results, function(x) x$diagnostics))
rownames(table7) <- rownames(table7_diagnostics) <- NULL

cat("Table 7: fixed 6,632 individuals, q = 1; assay counts may vary.\n")
if (requireNamespace("knitr", quietly = TRUE)) {
  print(knitr::kable(table7, format = "pipe", digits = 6))
} else {
  print(table7, row.names = FALSE, digits = 9)
}
cat("\nApproximate-design verification in the original parameter coordinates:\n")
print(table7_diagnostics, row.names = FALSE, digits = 6)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0L) {
  output_dir <- args[1L]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(table7, file.path(output_dir, "table7.csv"), row.names = FALSE)
  write.csv(table7_diagnostics, file.path(output_dir, "table7_diagnostics.csv"),
            row.names = FALSE)
  write.csv(do.call(rbind, lapply(table7_results, function(x) x$weights)),
            file.path(output_dir, "table7_approximate_weights.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lapply(table7_results, function(x) x$counts)),
            file.path(output_dir, "table7_exact_counts.csv"), row.names = FALSE)
  capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
}
