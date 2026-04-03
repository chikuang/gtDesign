#' E-optimal approximate design for group testing experiments
#'
#' Minimizes \eqn{-\lambda_{\min}(\mathbf{M})} over \eqn{\mathbf{w}} on a finite
#' candidate set via \code{CVXR::Minimize(-CVXR::lambda_min(\ldots))}. The
#' reported objective is \eqn{\mathrm{loss}(\mathbf{M}) = -\lambda_{\min}(\mathbf{M})},
#' matching MATLAB `loss = -lambda_min(FIM)` (e.g.\ `calc_loss_E`).
#'
#' This matches [compute_design_SO()] with `criterion = "E"` (same scalar loss).
#'
#' @inheritParams calc_Dopt
#'
#' @return An object of class `"gt_e_design"` and `"gt_design"`. Components
#'   `value` and `optval` are the minimization objective
#'   \eqn{-\lambda_{\min}(\mathbf{M})}: `value` is recomputed at the cleaned
#'   support, `optval` is the solver optimum of the minimize objective (they may
#'   differ slightly after `drop_tol`).
#'
#' @examples
#' \donttest{
#' # Huang et al. (2020) model on {1,...,61}, q = 0 (same setting as D-opt example)
#' theta <- c(p0 = 0.07, p1 = 0.93, p2 = 0.96)
#' u <- seq_len(61L)
#' f <- gt_huang2020_regressor(theta, q = 0)
#' res_E <- calc_Eopt(u, f, drop_tol = 1e-6)
#' res_E$design
#' res_E$value # = -lambda_min(M), same as -min(eigen(M)$values)
#' -res_E$value # lambda_min(M) at the optimum
#' res_E$status
#' }
#'
#' @export
calc_Eopt <- function(u,
                      f,
                      solver = "CLARABEL",
                      ...,
                      drop_tol = 1e-8,
                      ridge = 1e-8) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package `CVXR` is required.", call. = FALSE)
  }

  validate_design_inputs(u, f)
  Fmat <- build_model_matrix(u, f)

  n <- nrow(Fmat)
  p <- ncol(Fmat)

  w <- CVXR::Variable(n, nonneg = TRUE)

  M_expr <- Reduce(
    `+`,
    lapply(seq_len(n), function(i) {
      w[i] * tcrossprod(Fmat[i, ])
    })
  )

  M_psd <- M_expr + ridge * diag(p)
  problem <- CVXR::Problem(
    CVXR::Minimize(-CVXR::lambda_min(M_psd)),
    constraints = list(sum(w) == 1)
  )

  optval <- as.numeric(CVXR::psolve(problem, solver = solver, ...))
  solver_status <- CVXR::status(problem)

  w_opt <- as.numeric(CVXR::value(w))
  w_opt[abs(w_opt) < drop_tol] <- 0
  w_opt <- w_opt / sum(w_opt)

  M_opt <- info_matrix(w_opt, Fmat)
  crit_val <- scalar_loss_from_M(M_opt, "E", list())

  design <- data.frame(point = u, weight = w_opt)
  design <- design[design$weight > 0, , drop = FALSE]
  rownames(design) <- NULL

  out <- list(
    design = design,
    weights = w_opt,
    candidates = u,
    Fmat = Fmat,
    info_matrix = M_opt,
    criterion = "E",
    value = crit_val,
    optval = optval,
    status = solver_status,
    solver = solver,
    ridge = ridge
  )

  class(out) <- c("gt_e_design", "gt_design")
  out
}
