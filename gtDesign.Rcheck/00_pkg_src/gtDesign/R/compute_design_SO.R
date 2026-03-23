#' Single-objective optimal approximate design (group testing)
#'
#' Solves a convex optimization problem over a finite candidate set for one
#' criterion. The regression function `f` encodes rank-one contributions to the
#' (approximate) information matrix, as in nonlinear design or group testing
#' Fisher information assembled from pool outcomes.
#'
#' @param u Candidate design points.
#' @param f A function returning the regression vector at a single design point.
#' @param criterion One of `"D"`, `"A"`, `"Ds"`, `"c"`, or `"E"`.
#' @param opts Named list of extra options. For `"Ds"` use `opts$cVec_Ds`; for
#'   `"c"` use `opts$cVec_c`.
#' @param info_weight Optional function returning a nonnegative scalar multiplier
#'   for each rank-one information contribution.
#' @param solver Solver passed to [CVXR::psolve()].
#' @param ... Additional arguments passed to [CVXR::psolve()].
#' @param support_tol Weights smaller than this are dropped from the reported
#'   support.
#' @param drop_tol Numerical tolerance for tiny solver noise before support
#'   cleanup.
#'
#' @return A list of class `"gt_so_design"` and `"gt_design"`.
#' @export
compute_design_SO <- function(u,
                              f,
                              criterion = c("D", "A", "Ds", "c", "E"),
                              opts = list(),
                              info_weight = NULL,
                              solver = "CLARABEL",
                              ...,
                              support_tol = 1e-4,
                              drop_tol = 1e-8) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package `CVXR` is required.", call. = FALSE)
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package `tibble` is required.", call. = FALSE)
  }

  criterion <- match.arg(criterion)
  cr <- toupper(criterion)

  validate_design_inputs(u, f)

  p <- length(eval_regvec(u[1L], f))
  n <- length(u)

  w <- CVXR::Variable(n, nonneg = TRUE)
  del <- CVXR::Variable(1)

  M_expr <- fim_matrix_expr(u, w, f, info_weight)

  constraints <- list(
    sum(w) == 1,
    w >= 0
  )

  if (cr == "D") {
    constraints <- c(constraints, list(-CVXR::log_det(M_expr) <= del))
  } else if (cr == "A") {
    constraints <- c(constraints, list(CVXR::tr_inv(M_expr) <= del))
  } else if (cr == "E") {
    constraints <- c(constraints, list(-CVXR::lambda_min(M_expr) <= del))
  } else if (cr %in% c("DS", "C")) {
    cvec <- get_contrast_vec(criterion, opts, p)
    constraints <- c(constraints, list(CVXR::matrix_frac(cvec, M_expr) <= del))
  } else {
    stop("Unsupported criterion: ", criterion, call. = FALSE)
  }

  problem <- CVXR::Problem(
    objective = CVXR::Minimize(del),
    constraints = constraints
  )

  optval <- CVXR::psolve(problem, solver = solver, ...)
  status <- CVXR::status(problem)

  if (!status %in% c("optimal", "optimal_inaccurate")) {
    stop("Solver did not return an optimal solution. Status: ", status, call. = FALSE)
  }

  w_opt <- as.numeric(CVXR::value(w))

  if (any(!is.finite(w_opt))) {
    stop("Solver returned non-finite weights.", call. = FALSE)
  }
  if (sum(w_opt) <= 0) {
    stop("Solver returned invalid weights.", call. = FALSE)
  }

  w_opt[abs(w_opt) < drop_tol] <- 0
  w_opt[w_opt < support_tol] <- 0

  if (sum(w_opt) <= 0) {
    stop("All weights were removed by `drop_tol` and `support_tol`.", call. = FALSE)
  }

  w_opt <- w_opt / sum(w_opt)

  M_hat <- fim_matrix(u, w_opt, f, info_weight)
  loss <- scalar_loss_from_M(M_hat, criterion, opts)

  design <- tibble::tibble(
    point = u,
    weight = w_opt
  )

  design <- design[design$weight > 0, ]

  out <- list(
    criterion = criterion,
    design = design,
    weights = w_opt,
    candidates = u,
    info_matrix = M_hat,
    loss = loss,
    value = loss,
    optval = as.numeric(optval),
    status = status,
    solver = solver,
    support_tol = support_tol,
    drop_tol = drop_tol
  )

  class(out) <- c("gt_so_design", "gt_design")
  out
}
