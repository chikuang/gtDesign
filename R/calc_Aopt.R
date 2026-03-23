#' A-optimal approximate design for group testing experiments
#'
#' @inheritParams calc_Dopt
#'
#' @return An object of class `"gt_a_design"` and `"gt_design"`.
#' @export
calc_Aopt <- function(u,
                      f,
                      solver = "CLARABEL",
                      ...,
                      drop_tol = 1e-8) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package `CVXR` is required.", call. = FALSE)
  }

  validate_design_inputs(u, f)
  Fmat <- build_model_matrix(u, f)

  n <- nrow(Fmat)

  w <- CVXR::Variable(n, nonneg = TRUE)

  M_expr <- Reduce(
    `+`,
    lapply(seq_len(n), function(i) {
      w[i] * tcrossprod(Fmat[i, ])
    })
  )

  problem <- CVXR::Problem(
    CVXR::Minimize(CVXR::tr_inv(M_expr)),
    constraints = list(sum(w) == 1)
  )

  optval <- CVXR::psolve(problem, solver = solver, ...)
  solver_status <- CVXR::status(problem)

  w_opt <- as.numeric(CVXR::value(w))
  w_opt[abs(w_opt) < drop_tol] <- 0
  w_opt <- w_opt / sum(w_opt)

  M_opt <- info_matrix(w_opt, Fmat)
  crit_val <- sum(diag(solve(M_opt)))

  design <- data.frame(point = u, weight = w_opt)
  design <- design[design$weight > 0, , drop = FALSE]
  rownames(design) <- NULL

  out <- list(
    design = design,
    weights = w_opt,
    candidates = u,
    Fmat = Fmat,
    info_matrix = M_opt,
    criterion = "A",
    value = crit_val,
    optval = optval,
    status = solver_status,
    solver = solver
  )

  class(out) <- c("gt_a_design", "gt_design")
  out
}
