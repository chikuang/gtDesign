#' Check the equivalence theorem on a finite candidate set
#'
#' Uses the same directional derivatives as [calc_directional_derivatives()]
#' with a single criterion. For maximin designs, use
#' [check_equivalence_maximin()].
#'
#' @param design_obj Output from [calc_Dopt()], [calc_Aopt()], or [calc_copt()].
#' @param f Regression function returning a numeric vector.
#' @param u Candidate design points. If omitted, uses `design_obj$candidates`.
#' @param tol Numerical tolerance for checking nonpositivity and equality.
#' @param use_ginv Passed to [calc_directional_derivatives()] when `M` is singular.
#'
#' @return A list with directional derivative values and theorem checks.
#' @export
check_equivalence <- function(design_obj,
                              f,
                              u = NULL,
                              tol = 1e-6,
                              use_ginv = TRUE) {
  if (missing(design_obj) || is.null(design_obj$criterion)) {
    stop("`design_obj` must be a valid design object from gtDesign.", call. = FALSE)
  }

  if (inherits(design_obj, "gt_maximin_design") ||
        identical(design_obj$criterion, "maximin")) {
    stop(
      "Use check_equivalence_maximin() for maximin design objects.",
      call. = FALSE
    )
  }

  if (missing(f) || !is.function(f)) {
    stop("`f` must be a regression function.", call. = FALSE)
  }

  crit <- toupper(design_obj$criterion)
  if (!crit %in% c("D", "A", "C")) {
    stop("Equivalence theorem supports criteria 'D', 'A', and 'c'.", call. = FALSE)
  }

  if (is.null(u)) {
    u <- design_obj$candidates
  }

  M <- design_obj$info_matrix
  c_vec <- if (crit == "C") design_obj$cVec else NULL

  dd_list <- calc_directional_derivatives(
    u, M, f,
    criteria = crit,
    cVec = c_vec,
    use_ginv = use_ginv
  )
  nm <- switch(
    tolower(crit),
    d = "dD",
    a = "dA",
    c = "dc"
  )
  deriv_vals <- dd_list[[nm]]

  support_pts <- design_obj$design$point
  support_idx <- match(support_pts, u)

  max_violation <- max(deriv_vals, na.rm = TRUE)
  support_vals <- deriv_vals[support_idx]

  out <- list(
    candidate_points = u,
    directional_derivative = deriv_vals,
    support_points = support_pts,
    support_values = support_vals,
    max_violation = max_violation,
    all_nonpositive = all(deriv_vals <= tol, na.rm = TRUE),
    support_equal_zero = all(abs(support_vals) <= max(10 * tol, tol), na.rm = TRUE),
    criterion = design_obj$criterion,
    tol = tol
  )

  class(out) <- "gt_equivalence"
  out
}
