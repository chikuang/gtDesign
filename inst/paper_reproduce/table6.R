# Reproduce Table 6: optimal approximate designs for q = 0 and q = 1.
# ---------------------------------------------------------------------------
# Weights are resource proportions: assays when q = 0, individuals when q = 1.
# Approximate designs do not depend on the total budget. Tables 7 and 8 handle
# integer allocations under the individual-count and assay-count constraints.
#
# With an installed gtDesign package:
#   Rscript inst/paper_reproduce/table6.R
# From a development checkout:
#   Rscript -e 'devtools::load_all(); source("inst/paper_reproduce/table6.R")'
# Optional output directory for CSV results, diagnostics, and session information:
#   Rscript inst/paper_reproduce/table6.R table6_results

options(scipen = 999)
if (!exists("compute_design_SO", mode = "function")) {
  suppressPackageStartupMessages(library(gtDesign))
}

theta <- c(p0 = 0.022, p1 = 0.92, p2 = 0.965)
u <- seq_len(15L)
q_values <- c(0, 1)
criteria <- c("D", "A", "E", "c", "Ds")
opts <- list(cVec_c = c(0, 1, 1), cVec_Ds = c(1, 0, 0))

physical_loss <- function(M, criterion) {
  switch(criterion,
    D = -as.numeric(determinant(M, logarithm = TRUE)$modulus),
    A = sum(diag(solve(M))),
    E = -min(eigen(M, symmetric = TRUE, only.values = TRUE)$values),
    c = as.numeric(crossprod(opts$cVec_c, solve(M, opts$cVec_c))),
    Ds = as.numeric(crossprod(opts$cVec_Ds, solve(M, opts$cVec_Ds)))
  )
}

run_one <- function(q_cost, criterion) {
  # h(x)h(x)' = H(x)/c(x;q), with H(x) the physical per-assay information.
  h <- gt_huang2020_regressor(theta, q_cost)
  F <- t(vapply(u, h, numeric(3)))
  eig <- eigen(crossprod(F) / length(u), symmetric = TRUE)

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
  stopifnot(
    approx$status %in% c("optimal", "optimal_inaccurate"),
    length(approx$weights) == length(u),
    all(is.finite(approx$weights)), all(approx$weights >= 0),
    abs(sum(approx$weights) - 1) < 1e-8
  )

  # Verify optimality in the ORIGINAL parameter coordinates over all 15 sizes.
  # Solver status alone is insufficient for these poorly conditioned models.
  M_approx <- crossprod(F, approx$weights * F)
  min_eigenvalue <- min(eigen(M_approx, symmetric = TRUE,
                              only.values = TRUE)$values)
  stopifnot(min_eigenvalue > 0)
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

  list(
    weights = data.frame(
      q = q_cost, Criterion = criterion, pool_size = u,
      resource_weight = approx$weights
    ),
    diagnostics = data.frame(
      q = q_cost, Criterion = criterion, solver_status = approx$status,
      relative_sensitivity_gap = relative_gap,
      loss_approx = loss_approx, min_eigenvalue = min_eigenvalue
    )
  )
}

table6_results <- list()
for (q_cost in q_values) {
  for (criterion in criteria) {
    key <- paste(q_cost, criterion, sep = "_")
    table6_results[[key]] <- run_one(q_cost, criterion)
  }
}
table6_weights <- do.call(rbind, lapply(table6_results, function(x) x$weights))
table6_diagnostics <- do.call(rbind, lapply(table6_results, function(x) x$diagnostics))
rownames(table6_weights) <- rownames(table6_diagnostics) <- NULL

# These are the unions of the supports in each panel of Table 6. In particular,
# the q = 1 D-optimal design has FOUR support points, including both 4 and 5.
panel_sizes <- list(`0` = c(1L, 7L, 15L), `1` = c(1L, 4L, 5L, 15L))
table6_full_precision <- data.frame(Criterion = criteria, stringsAsFactors = FALSE)
for (q_cost in q_values) {
  sizes <- panel_sizes[[as.character(q_cost)]]
  panel <- table6_weights[table6_weights$q == q_cost, ]
  stopifnot(all(panel$resource_weight[!panel$pool_size %in% sizes] < 1e-6))
  for (size in sizes) {
    column <- sprintf("q = %d: pool size %d", q_cost, size)
    table6_full_precision[[column]] <- vapply(criteria, function(criterion) {
      fit <- table6_results[[paste(q_cost, criterion, sep = "_")]]
      fit$weights$resource_weight[match(size, u)]
    }, numeric(1))
  }
}

# Merge A and E only for display, and only if every printed weight agrees.
# Their distinct unrounded weights remain available in the full-precision data.
table6 <- table6_full_precision
table6[-1] <- lapply(table6[-1], round, digits = 3)
if (all(table6[table6$Criterion == "A", -1] ==
        table6[table6$Criterion == "E", -1])) {
  table6 <- table6[table6$Criterion != "E", ]
  table6$Criterion[table6$Criterion == "A"] <- "A, E"
}
rownames(table6) <- NULL

cat("Table 6: Chlamydia optimal approximate designs, M = 15.\n")
cat("q = 0: assay-budget proportions; q = 1: individual-budget proportions.\n")
cat("Displayed weights are rounded to three decimals and may not sum to one.\n")
if (requireNamespace("knitr", quietly = TRUE)) {
  print(knitr::kable(table6, format = "pipe", digits = 3))
} else {
  print(table6, row.names = FALSE, digits = 3)
}
cat("\nApproximate-design verification in the original parameter coordinates:\n")
print(table6_diagnostics, row.names = FALSE, digits = 6)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0L) {
  output_dir <- args[1L]
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(table6, file.path(output_dir, "table6.csv"), row.names = FALSE)
  write.csv(table6_full_precision, file.path(output_dir, "table6_full_precision.csv"),
            row.names = FALSE)
  write.csv(table6_weights, file.path(output_dir, "table6_approximate_weights.csv"),
            row.names = FALSE)
  write.csv(table6_diagnostics, file.path(output_dir, "table6_diagnostics.csv"),
            row.names = FALSE)
  capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
}
