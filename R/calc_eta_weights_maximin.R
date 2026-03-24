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


#' Build CVXR eta LP (internal)
#'
#' @return List with `problem` and `eta` ([CVXR::Variable()]).
#' @noRd
build_eta_problem <- function(K,
                              criteria,
                              dd_names,
                              directional_derivatives,
                              gaps,
                              rhs_derivs,
                              tol_deriv,
                              tol_comp,
                              use_complementary_slack) {
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

  if (isTRUE(use_complementary_slack)) {
    for (k in seq_len(K)) {
      key <- canon_crit_key(criteria[k])
      delta <- gaps[[key]]

      constraints <- c(constraints, list(
        eta[k] * delta <= tol_comp,
        -eta[k] * delta <= tol_comp
      ))
    }
  }

  sum_deriv <- 0
  for (k in seq_len(K)) {
    field <- dd_names[k]
    sum_deriv <- sum_deriv + eta[k] * directional_derivatives[[field]]
  }

  constraints <- c(constraints, list(sum_deriv <= tol_deriv))

  problem <- CVXR::Problem(
    objective = CVXR::Minimize(sum(eta)),
    constraints = constraints
  )

  list(problem = problem, eta = eta)
}


#' Eta weights for the maximin equivalence theorem
#'
#' Solves for nonnegative weights in the maximin equivalence theorem. The
#' returned vector `eta` satisfies the normalization condition induced by the
#' derivative of the active right-hand sides and enforces a nonpositive weighted
#' combination of directional derivatives across the candidate set.
#'
#' The linear program from Gao et al.\ also imposes approximate **complementary
#' slackness** `|eta_k * gap_k| <= tol`. With tight `tol`, that set can be
#' **empty** together with the normalization constraint (common with three or
#' more criteria). By default this function **retries without** complementary
#' slackness if the strict problem is infeasible, and cycles a few solvers
#' (works well with **CVXR** 1.8.x).
#'
#' @param tstar Optimal maximin scalar returned by [compute_maximin_design()].
#' @param loss_ref Named list of reference losses from single-objective designs.
#' @param loss_model Named list of achieved losses at the maximin design.
#' @param directional_derivatives Named list of directional derivative vectors,
#'   such as `dD`, `dA`, `dDs`, `dc`, and `dE`.
#' @param criteria Character vector of criteria, e.g. `c("D", "A")`.
#' @param q Parameter dimension.
#' @param tol Numerical tolerance for the combined directional derivative and (if
#'   used) complementary slackness.
#' @param complementary_slack Logical; if `TRUE`, first attempt includes
#'   complementary slackness. If that LP is infeasible and `fallback_no_slack` is
#'   `TRUE`, a second attempt omits those constraints.
#' @param fallback_no_slack Logical; default `TRUE`. See `complementary_slack`.
#' @param complementary_tol_mult Positive multiplier applied to `tol` for the
#'   complementary slackness bounds (helps feasibility).
#' @param solver Solver passed to [CVXR::psolve()] for the first attempt.
#' @param solvers_fallback Character vector of solvers to try if the first
#'   attempt fails (default includes SCS and ECOS, which often succeed when
#'   CLARABEL reports infeasible on small LPs).
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
                                     complementary_slack = TRUE,
                                     fallback_no_slack = TRUE,
                                     complementary_tol_mult = 100,
                                     solver = "CLARABEL",
                                     solvers_fallback = c("SCS", "ECOS"),
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
  if (!is.numeric(complementary_tol_mult) ||
        length(complementary_tol_mult) != 1L ||
        !is.finite(complementary_tol_mult) ||
        complementary_tol_mult <= 0) {
    stop("`complementary_tol_mult` must be a positive finite scalar.", call. = FALSE)
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
  tol_comp <- tol * complementary_tol_mult
  solvers <- unique(c(solver, solvers_fallback))
  state <- list(attempts = character())

  run_solvers <- function(use_slack) {
    for (slv in solvers) {
      built <- build_eta_problem(
        K = K,
        criteria = criteria,
        dd_names = dd_names,
        directional_derivatives = directional_derivatives,
        gaps = gaps,
        rhs_derivs = rhs_derivs,
        tol_deriv = tol,
        tol_comp = tol_comp,
        use_complementary_slack = use_slack
      )
      ok <- tryCatch(
        {
          CVXR::psolve(built$problem, solver = slv, ...)
          TRUE
        },
        error = function(e) FALSE
      )
      st <- if (ok) CVXR::status(built$problem) else "solver_error"
      state$attempts <- c(state$attempts, paste0(slv, "/", use_slack, ":", st))
      if (ok && st %in% c("optimal", "optimal_inaccurate")) {
        eta <- as.numeric(CVXR::value(built$eta))
        if (all(is.finite(eta))) {
          return(eta)
        }
      }
    }
    NULL
  }

  eta_opt <- NULL
  if (isTRUE(complementary_slack)) {
    eta_opt <- run_solvers(TRUE)
  }
  if (is.null(eta_opt)) {
    if (isTRUE(complementary_slack) && isTRUE(fallback_no_slack)) {
      warning(
        "Complementary slackness constraints were infeasible; ",
        "solving again without them (normalization + derivative constraints only).",
        call. = FALSE
      )
    }
    if (isTRUE(fallback_no_slack) || !isTRUE(complementary_slack)) {
      eta_opt <- run_solvers(FALSE)
    }
  }

  if (is.null(eta_opt)) {
    stop(
      "Could not solve for eta. Attempts: ",
      paste(state$attempts, collapse = "; "),
      call. = FALSE
    )
  }

  names(eta_opt) <- criteria
  eta_opt
}
