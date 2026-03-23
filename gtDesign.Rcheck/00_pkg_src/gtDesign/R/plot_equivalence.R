#' Plot equivalence theorem directional derivative
#'
#' @importFrom graphics abline plot points
#'
#' @param eq_obj Output from [check_equivalence()].
#' @param main Plot title.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param pch_support Plotting symbol for support points.
#' @param ... Passed to [graphics::plot()].
#'
#' @return Invisibly returns `eq_obj`.
#' @export
plot_equivalence <- function(eq_obj,
                             main = NULL,
                             xlab = "Design point",
                             ylab = "Directional derivative",
                             pch_support = 19,
                             ...) {
  if (!inherits(eq_obj, "gt_equivalence")) {
    stop("`eq_obj` must be an object returned by check_equivalence().", call. = FALSE)
  }

  if (is.null(main)) {
    main <- paste(eq_obj$criterion, "-optimality equivalence theorem")
  }

  plot(
    eq_obj$candidate_points,
    eq_obj$directional_derivative,
    type = "l",
    xlab = xlab,
    ylab = ylab,
    main = main,
    ...
  )
  abline(h = 0, lty = 2)

  support_idx <- match(eq_obj$support_points, eq_obj$candidate_points)
  points(
    eq_obj$candidate_points[support_idx],
    eq_obj$directional_derivative[support_idx],
    pch = pch_support
  )

  invisible(eq_obj)
}
