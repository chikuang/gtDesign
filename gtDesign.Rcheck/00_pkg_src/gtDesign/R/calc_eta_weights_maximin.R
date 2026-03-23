#' Compute active-constraint gaps for a maximin design
#'
#' @param tstar Optimal maximin scalar.
#' @param loss_ref Named list of reference losses.
#' @param loss_model Named list of achieved losses at the maximin design.
#' @param criteria Character vector of criteria.
#' @param q Parameter dimension.
#'
#' @return A named numeric vector of constraint gaps.
#' @noRd
compute_maximin_gaps <- function(tstar, loss_ref, loss_model, criteria, q) {
  criteria <- standardize_criteria(criteria)

  gaps <- numeric(length(criteria))
  names(gaps) <- criteria

  for (cr in criteria) {
    key <- canon_crit_key(cr)

    if (key == "D") {
      gaps[key] <- loss_model$D - (loss_ref$D + q * log(tstar))
    } else if (key == "E") {
      gaps[key] <- loss_model$E - loss_ref$E / tstar
    } else {
      gaps[key] <- loss_model[[key]] - tstar * loss_ref[[key]]
    }
  }

  gaps
}


#' Derivatives of the maximin right-hand sides
#'
#' @param tstar Optimal maximin scalar.
#' @param loss_ref Named list of reference losses.
#' @param criteria Character vector of criteria.
#' @param q Parameter dimension.
#'
#' @return A named numeric vector of derivatives h'_j(tstar).
#' @noRd
compute_maximin_rhs_derivatives <- function(tstar, loss_ref, criteria, q) {
  criteria <- standardize_criteria(criteria)

  derivs <- numeric(length(criteria))
  names(derivs) <- criteria

  for (cr in criteria) {
    key <- canon_crit_key(cr)

    if (key == "D") {
      derivs[key] <- q / tstar
    } else if (key == "E") {
      derivs[key] <- -loss_ref$E / (tstar^2)
    } else {
      derivs[key] <- loss_ref[[key]]
    }
  }

  derivs
}


#' Eta weights for the maximin equivalence theorem
#'
#' Solves for nonnegative weights in the maximin equivalence theorem. The
#' returned vector `eta` satisfies the normalization condition induced by the
#' derivative of the active right-hand sides and enforces a nonpositive weighted
#' combination of directional derivatives across the candidate set.
#'
#' @param tstar Optimal maximin scalar returned by [compute_maximin_design()].
#' @param loss_ref Named list of reference losses from single-objective designs.
#' @param loss_model Named list of achieved losses at the maximin design.
#' @param directional_derivatives Named list of directional derivative vectors,
#'   such as `dD`, `dA`, `dDs`, `dc`, and `dE`.
#' @param criteria Character vector of criteria, e.g. `c("D", "A")`.
#' @param q Parameter dimension.
#' @param tol Numerical tolerance.
#' @param solver Solver passed to [CVXR::psolve()].
#' @param ... Additional arguments passed to [CVXR::psolve()].
#'
#' @return A named numeric vector of eta weights.
#' @export
calc_eta_weights_maximin <- function(tstar,
                                     loss_ref,
                                     loss_model,
                                     directional_derivatives,
                                     criteria,
                                     q,
                                     tol = 1e-6,
                                     solver = "CLARABEL",
                                     ...) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    stop("Package `CVXR` is required.", call. = FALSE)
  }

  criteria <- standardize_criteria(criteria)
  validate_loss_ref(loss_ref, criteria)

  if (!is.list(loss_model)) {
    stop("`loss_model` must be a named list.", call. = FALSE)
  }
  if (!is.list(directional_derivatives)) {
    stop("`directional_derivatives` must be a named list.", call. = FALSE)
  }
  if (!is.numeric(tstar) || length(tstar) != 1L || !is.finite(tstar) || tstar <= 0) {
    stop("`tstar` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(q) || length(q) != 1L || !is.finite(q) || q <= 0) {
    stop("`q` must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a positive finite scalar.", call. = FALSE)
  }

  crit_to_dd <- list(
    D = "dD", A = "dA", Ds = "dDs", c = "dc", E = "dE"
  )

  dd_names <- vapply(criteria, function(cr) {
    key <- as.character(canon_crit_key(cr))
    out <- crit_to_dd[[key]]
    if (is.null(out)) {
      stop("Unsupported criterion: ", cr, call. = FALSE)
    }
    out
  }, character(1))

  missing_dd <- dd_names[!dd_names %in% names(directional_derivatives)]
  if (length(missing_dd) > 0) {
    stop(
      "Missing directional derivative field(s): ",
      paste(unique(missing_dd), collapse = ", "),
      call. = FALSE
    )
  }

  N <- length(directional_derivatives[[dd_names[1]]])
  if (N == 0) {
    stop("Directional derivative vectors must be non-empty.", call. = FALSE)
  }

  for (nm in dd_names) {
    if (length(directional_derivatives[[nm]]) != N) {
      stop("All directional derivative vectors must have the same length.", call. = FALSE)
    }
  }

  gaps <- compute_maximin_gaps(
    tstar = tstar,
    loss_ref = loss_ref,
    loss_model = loss_model,
    criteria = criteria,
    q = q
  )

  rhs_derivs <- compute_maximin_rhs_derivatives(
    tstar = tstar,
    loss_ref = loss_ref,
    criteria = criteria,
    q = q
  )

  K <- length(criteria)
  eta <- CVXR::Variable(K, nonneg = TRUE)

  norm_expr <- 0
  for (k in seq_len(K)) {
    key <- canon_crit_key(criteria[k])
    norm_expr <- norm_expr + eta[k] * rhs_derivs[[key]]
  }

  constraints <- list(
    norm_expr == 1,
    eta >= 0
  )

  for (k in seq_len(K)) {
    key <- canon_crit_key(criteria[k])
    delta <- gaps[[key]]

    constraints <- c(constraints, list(
      eta[k] * delta <= tol,
      -eta[k] * delta <= tol
    ))
  }

  sum_deriv <- 0
  for (k in seq_len(K)) {
    field <- dd_names[k]
    sum_deriv <- sum_deriv + eta[k] * directional_derivatives[[field]]
  }

  constraints <- c(constraints, list(sum_deriv <= tol))

  problem <- CVXR::Problem(
    objective = CVXR::Minimize(sum(eta)),
    constraints = constraints
  )

  optval <- CVXR::psolve(problem, solver = solver, ...)
  status <- CVXR::status(problem)

  if (!status %in% c("optimal", "optimal_inaccurate")) {
    stop("Solver did not return an optimal solution. Status: ", status, call. = FALSE)
  }

  eta_opt <- as.numeric(CVXR::value(eta))
  names(eta_opt) <- criteria

  eta_opt
}
