#' D-optimal approximate design for group testing experiments
#'
#' Computes a D-optimal approximate design on a finite candidate set. The
#' regression map `f` typically encodes contributions to the Fisher information
#' matrix (e.g. score vectors or square-root information); see
#' [gt_huang2020_regressor()] for the Huang et al.\ (2020) three-parameter model.
#'
#' @param u Candidate design points (e.g. pool sizes or coded pool layouts).
#' @param f Function returning a numeric vector at one design point.
#' @param solver Solver name passed to CVXR. Default is `"CLARABEL"`.
#' @param ... Additional arguments passed to [CVXR::psolve()].
#' @param drop_tol Threshold for removing near-zero weights.
#' @param ridge Small diagonal ridge for numerical stability in `log_det`.
#' @param use_matrix_form If `TRUE`, build `M` via dense matrix multiply.
#'
#' @return An object of class `"gt_d_design"` and `"gt_design"`.
#' @export
calc_Dopt <- function(u,
                      f,
                      solver = "CLARABEL",
                      ...,
                      drop_tol = 1e-8,
                      ridge = 1e-8,
                      use_matrix_form = FALSE) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package `CVXR` is required.", call. = FALSE)
  }

  validate_design_inputs(u, f)
  Fmat <- build_model_matrix(u, f)

  n <- nrow(Fmat)
  p <- ncol(Fmat)

  w <- CVXR::Variable(n, nonneg = TRUE)

  if (use_matrix_form) {
    F_cvx <- CVXR::Constant(Fmat)
    M_expr <- t(F_cvx) %*% CVXR::DiagVec(w) %*% F_cvx
  } else {
    M_expr <- Reduce(
      `+`,
      lapply(seq_len(n), function(i) {
        fi <- matrix(Fmat[i, ], ncol = 1)
        w[i] * (fi %*% t(fi))
      })
    )
  }

  problem <- CVXR::Problem(
    CVXR::Maximize(CVXR::log_det(M_expr + ridge * diag(p))),
    constraints = list(sum(w) == 1)
  )

  optval <- CVXR::psolve(problem, solver = solver, ...)
  solver_status <- CVXR::status(problem)

  w_opt <- as.numeric(CVXR::value(w))
  w_opt[abs(w_opt) < drop_tol] <- 0
  w_opt <- w_opt / sum(w_opt)

  M_opt <- info_matrix(w_opt, Fmat)
  crit_val <- as.numeric(determinant(M_opt, logarithm = TRUE)$modulus)

  design <- data.frame(point = u, weight = w_opt)
  design <- design[design$weight > 0, , drop = FALSE]
  rownames(design) <- NULL

  out <- list(
    design = design,
    weights = w_opt,
    candidates = u,
    Fmat = Fmat,
    info_matrix = M_opt,
    criterion = "D",
    value = crit_val,
    optval = optval,
    status = solver_status,
    solver = solver,
    ridge = ridge
  )

  class(out) <- c("gt_d_design", "gt_design")
  out
}
