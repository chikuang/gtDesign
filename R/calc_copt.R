#' c-optimal approximate design for group testing experiments
#'
#' @inheritParams calc_Dopt
#' @param cVec Numeric vector specifying the linear combination \eqn{c^\top \theta}.
#'
#' @return An object of class `"gt_c_design"` and `"gt_design"`.
#' @export
calc_copt <- function(u,
                      f,
                      cVec,
                      solver = "CLARABEL",
                      ...,
                      drop_tol = 1e-10) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package `CVXR` is required.", call. = FALSE)
  }

  validate_design_inputs(u, f)
  Fmat <- build_model_matrix(u, f)

  p <- ncol(Fmat)

  if (missing(cVec)) {
    stop("`cVec` must be provided.", call. = FALSE)
  }
  if (!is.numeric(cVec)) {
    stop("`cVec` must be a numeric vector.", call. = FALSE)
  }
  if (length(cVec) != p) {
    stop("Length of `cVec` must match the dimension of the regression vector.", call. = FALSE)
  }

  n <- nrow(Fmat)
  w <- CVXR::Variable(n, nonneg = TRUE)

  M_expr <- Reduce(
    `+`,
    lapply(seq_len(n), function(i) {
      w[i] * tcrossprod(Fmat[i, ])
    })
  )

  cmat <- matrix(as.numeric(cVec), ncol = 1)

  problem <- CVXR::Problem(
    CVXR::Minimize(CVXR::matrix_frac(cmat, M_expr)),
    constraints = list(sum(w) == 1)
  )

  optval <- CVXR::psolve(problem, solver = solver, ...)
  solver_status <- CVXR::status(problem)

  w_opt <- as.numeric(CVXR::value(w))
  w_opt[abs(w_opt) < drop_tol] <- 0

  if (sum(w_opt) <= 0) {
    stop("Solver returned invalid weights.", call. = FALSE)
  }

  w_opt <- w_opt / sum(w_opt)

  M_opt <- info_matrix(w_opt, Fmat)

  crit_val <- as.numeric(optval)

  direct_val <- tryCatch(
    as.numeric(t(cVec) %*% solve(M_opt, cVec)),
    error = function(e) NA_real_
  )

  if (!is.na(direct_val) && is.finite(direct_val)) {
    crit_val <- direct_val
    used_ginv <- FALSE
  } else {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      warning(
        "Matrix inversion failed and package `MASS` is not available; using solver objective value.",
        call. = FALSE
      )
      used_ginv <- NA
    } else {
      crit_val <- as.numeric(t(cVec) %*% MASS::ginv(M_opt) %*% cVec)
      used_ginv <- TRUE
    }
  }

  design <- data.frame(point = u, weight = w_opt)
  design <- design[design$weight > 0, , drop = FALSE]
  rownames(design) <- NULL

  out <- list(
    design = design,
    weights = w_opt,
    candidates = u,
    Fmat = Fmat,
    info_matrix = M_opt,
    cVec = as.numeric(cVec),
    criterion = "c",
    value = crit_val,
    optval = as.numeric(optval),
    status = solver_status,
    solver = solver,
    used_ginv = used_ginv
  )

  class(out) <- c("gt_c_design", "gt_design")
  out
}
