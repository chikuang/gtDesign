#' Check the equivalence theorem for a maximin design
#'
#' Computes the weighted combined directional derivative
#' \eqn{\sum_j \eta_j d_j(x)} over the candidate set and checks whether it is
#' nonpositive, with approximate equality at the support points.
#'
#' @param design_obj Output from [compute_maximin_design()].
#' @param directional_derivatives Named list of directional derivative vectors,
#'   such as `dD`, `dA`, `dDs`, `dc`, and `dE`.
#' @param eta Named numeric vector of eta weights, typically returned by
#'   [calc_eta_weights_maximin()].
#' @param tol Numerical tolerance used in the checks.
#'
#' @return A list of class `gt_equivalence_maximin`.
#' @export
check_equivalence_maximin <- function(design_obj,
                                      directional_derivatives,
                                      eta,
                                      tol = 1e-6) {
  if (missing(design_obj) || !is.list(design_obj)) {
    stop("`design_obj` must be a valid maximin design object.", call. = FALSE)
  }
  if (is.null(design_obj$candidates) || is.null(design_obj$design)) {
    stop("`design_obj` must contain `candidates` and `design`.", call. = FALSE)
  }
  if (!is.list(directional_derivatives) || length(directional_derivatives) == 0) {
    stop("`directional_derivatives` must be a non-empty named list.", call. = FALSE)
  }
  if (!is.numeric(eta) || is.null(names(eta))) {
    stop("`eta` must be a named numeric vector.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a positive finite scalar.", call. = FALSE)
  }

  u <- design_obj$candidates
  support_pts <- design_obj$design$point

  d_multi <- calc_multi_directional_derivative(
    dd_list = directional_derivatives,
    eta = eta
  )

  if (length(d_multi) != length(u)) {
    stop("Combined directional derivative has incompatible length.", call. = FALSE)
  }

  support_idx <- match(support_pts, u)
  if (anyNA(support_idx)) {
    warning(
      "Some support points could not be matched to `design_obj$candidates`.",
      call. = FALSE
    )
  }

  support_vals <- d_multi[support_idx]
  max_violation <- max(d_multi, na.rm = TRUE)
  min_value <- min(d_multi, na.rm = TRUE)

  out <- list(
    candidate_points = u,
    support_points = support_pts,
    eta = eta,
    combined_directional_derivative = d_multi,
    support_values = support_vals,
    max_violation = max_violation,
    min_value = min_value,
    all_nonpositive = all(d_multi <= tol, na.rm = TRUE),
    support_equal_zero = all(abs(support_vals) <= max(10 * tol, tol), na.rm = TRUE),
    tol = tol
  )

  class(out) <- "gt_equivalence_maximin"
  out
}
