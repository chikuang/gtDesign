#' Maximin multi-criterion approximate design (group testing)
#'
#' Joint convex formulation for maximizing the minimum efficiency across several
#' criteria, using reference losses from single-objective designs. Suitable
#' when multiple experimental goals (e.g. overall prevalence versus contrasts)
#' must be balanced in group testing allocation.
#'
#' @inheritParams compute_design_SO
#' @param loss_ref Named list of reference losses from single-objective
#'   designs, on the same scale as [compute_design_SO()] internal losses.
#' @param criteria Character vector containing any of `"D"`, `"A"`, `"Ds"`,
#'   `"c"`, or `"E"`.
#'
#' @return A list of class `"gt_maximin_design"` and `"gt_design"`.
#' @export
compute_maximin_design <- function(u,
                                   f,
                                   loss_ref,
                                   criteria,
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

  validate_design_inputs(u, f)

  criteria <- standardize_criteria(criteria)
  validate_loss_ref(loss_ref, criteria)

  p <- length(eval_regvec(u[1L], f))
  n <- length(u)

  w <- CVXR::Variable(n, nonneg = TRUE)
  tstar <- CVXR::Variable(1)

  M_expr <- fim_matrix_expr(u, w, f, info_weight)

  constraints <- list(
    tstar >= 1e-8,
    w >= 0,
    sum(w) == 1
  )

  for (cr in criteria) {
    if (cr == "D") {
      constraints <- c(constraints, list(
        -CVXR::log_det(M_expr) <= loss_ref$D + p * log(tstar)
      ))
    } else if (cr == "A") {
      constraints <- c(constraints, list(
        CVXR::tr_inv(M_expr) <= tstar * loss_ref$A
      ))
    } else if (cr == "Ds") {
      c_ds <- get_contrast_vec("Ds", opts, p)
      constraints <- c(constraints, list(
        CVXR::matrix_frac(c_ds, M_expr) <= tstar * loss_ref$Ds
      ))
    } else if (cr == "c") {
      c_c <- get_contrast_vec("c", opts, p)
      constraints <- c(constraints, list(
        CVXR::matrix_frac(c_c, M_expr) <= tstar * loss_ref$c
      ))
    } else if (cr == "E") {
      constraints <- c(constraints, list(
        -CVXR::lambda_min(M_expr) <= loss_ref$E * CVXR::inv_pos(tstar)
      ))
    } else {
      stop("Unsupported criterion: ", cr, call. = FALSE)
    }
  }

  problem <- CVXR::Problem(
    objective = CVXR::Minimize(tstar),
    constraints = constraints
  )

  optval <- CVXR::psolve(problem, solver = solver, ...)
  status <- CVXR::status(problem)

  if (!status %in% c("optimal", "optimal_inaccurate")) {
    stop(
      "Solver did not return an optimal solution. Status: ",
      status,
      call. = FALSE
    )
  }

  w_opt <- as.numeric(CVXR::value(w))
  t_opt <- as.numeric(CVXR::value(tstar))

  if (any(!is.finite(w_opt))) {
    stop("Solver returned non-finite weights.", call. = FALSE)
  }
  if (!is.finite(t_opt) || t_opt <= 0) {
    stop("Solver returned an invalid `tstar` value.", call. = FALSE)
  }
  if (sum(w_opt) <= 0) {
    stop("Solver returned invalid weights.", call. = FALSE)
  }

  w_opt[abs(w_opt) < drop_tol] <- 0
  w_opt[w_opt < support_tol] <- 0

  if (sum(w_opt) <= 0) {
    stop(
      "All weights were removed by `drop_tol` and `support_tol`.",
      call. = FALSE
    )
  }

  w_opt <- w_opt / sum(w_opt)

  M_hat <- fim_matrix(u, w_opt, f, info_weight)
  loss <- compute_losses_at_M(M_hat, criteria, opts)
  efficiency <- compute_efficiencies_maximin(loss_ref, loss, criteria, q = p)

  design <- tibble::tibble(
    point = u,
    weight = w_opt
  )
  design <- design[design$weight > 0, ]

  out <- list(
    criterion = criteria,
    design = design,
    weights = w_opt,
    candidates = u,
    loss = loss,
    efficiency = efficiency,
    value = min(efficiency, na.rm = TRUE),
    tstar = t_opt,
    M = M_hat,
    info_matrix = M_hat,
    optval = as.numeric(optval),
    status = status,
    solver = solver,
    support_tol = support_tol,
    drop_tol = drop_tol
  )

  class(out) <- c("gt_maximin_design", "gt_design")
  out
}
