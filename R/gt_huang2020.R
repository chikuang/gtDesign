#' Parse prevalence / sensitivity / specificity parameter vector
#'
#' @param theta Named or unnamed length-3 vector \eqn{(p_0, p_1, p_2)}.
#'
#' @return A named numeric vector with names `p0`, `p1`, `p2`.
#' @noRd
parse_theta_p012 <- function(theta) {
  if (missing(theta) || length(theta) != 3L) {
    stop("`theta` must have length 3 (p0, p1, p2).", call. = FALSE)
  }
  if (is.null(names(theta)) || any(!nzchar(names(theta)))) {
    names(theta) <- c("p0", "p1", "p2")
  }
  if (!all(c("p0", "p1", "p2") %in% names(theta))) {
    stop("`theta` must be named p0, p1, p2 (or unnamed in that order).", call. = FALSE)
  }
  p0 <- as.numeric(theta[["p0"]])
  p1 <- as.numeric(theta[["p1"]])
  p2 <- as.numeric(theta[["p2"]])
  if (length(p0) != 1L || p0 <= 0 || p0 >= 1) {
    stop("`p0` must lie strictly between 0 and 1.", call. = FALSE)
  }
  if (length(p1) != 1L || p1 <= 0.5 || p1 > 1) {
    stop("`p1` (sensitivity) must lie in (0.5, 1].", call. = FALSE)
  }
  if (length(p2) != 1L || p2 <= 0.5 || p2 > 1) {
    stop("`p2` (specificity) must lie in (0.5, 1].", call. = FALSE)
  }
  c(p0 = p0, p1 = p1, p2 = p2)
}


#' Positive test probability (Huang et al. 2020 / arXiv:2508.08445 Sec. 2)
#'
#' \deqn{\pi(x \mid \theta) = p_1 - (p_1 + p_2 - 1)(1 - p_0)^x}
#'
#' @param x Pool size (positive integer).
#' @param theta \eqn{(p_0, p_1, p_2)} prevalence, sensitivity, specificity.
#'
#' @return Scalar probability in `(0, 1)` when well-defined.
#' @export
gt_huang2020_pi <- function(x, theta) {
  th <- parse_theta_p012(theta)
  m <- as.integer(x)
  if (length(m) != 1L || !is.finite(m) || m < 1L) {
    stop("`x` must be a positive integer pool size.", call. = FALSE)
  }
  p0 <- th[["p0"]]
  p1 <- th[["p1"]]
  p2 <- th[["p2"]]
  z <- (1 - p0)^m
  p1 - (p1 + p2 - 1) * z
}


#' Standardized cost for pool size (Huang et al. 2020)
#'
#' \deqn{c(x) = 1 - q + q x}
#' with \eqn{q = q_1 / (q_0 + q_1)} for assay cost \eqn{q_0} and enrollment cost
#' \eqn{q_1}. Value \eqn{q = 0} gives \eqn{c(x) \equiv 1}.
#'
#' @param x Pool size (positive integer).
#' @param q Cost ratio in `[0, 1]`.
#'
#' @return Scalar cost.
#' @export
gt_huang2020_cost <- function(x, q) {
  m <- as.integer(x)
  if (length(m) != 1L || !is.finite(m) || m < 1L) {
    stop("`x` must be a positive integer pool size.", call. = FALSE)
  }
  if (length(q) != 1L || !is.finite(q) || q < 0 || q > 1) {
    stop("`q` must be a scalar in [0, 1].", call. = FALSE)
  }
  1 - q + q * m
}


#' Weight \eqn{\lambda(x)} in the information matrix (Huang et al. 2020)
#'
#' \deqn{\lambda(x) = [c(x)\,\pi(x)\{1-\pi(x)\}]^{-1}}
#'
#' @param x Pool size.
#' @param theta Same as for \code{\link{gt_huang2020_pi}}.
#' @param q Same as for \code{\link{gt_huang2020_cost}}.
#'
#' @return Scalar (positive when \eqn{\pi} is interior).
#' @export
gt_huang2020_lambda <- function(x, theta, q) {
  pi_x <- gt_huang2020_pi(x, theta)
  if (pi_x <= 0 || pi_x >= 1) {
    stop(
      "`pi(x)` is on the boundary; Fisher weight `lambda` is not finite.",
      call. = FALSE
    )
  }
  cx <- gt_huang2020_cost(x, q)
  1 / (cx * pi_x * (1 - pi_x))
}


#' Gradient vector \eqn{\mathbf{f}(x)} for the Huang et al. (2020) model
#'
#' Partial derivatives of \eqn{\pi(x\mid\theta)} with respect to
#' \eqn{(p_0,p_1,p_2)} as in arXiv:2508.08445 (Eq. below their Eq. (1)).
#'
#' @param x Pool size.
#' @param theta \eqn{(p_0, p_1, p_2)}.
#'
#' @return Length-3 numeric vector.
#' @export
gt_huang2020_f <- function(x, theta) {
  th <- parse_theta_p012(theta)
  m <- as.integer(x)
  if (length(m) != 1L || !is.finite(m) || m < 1L) {
    stop("`x` must be a positive integer pool size.", call. = FALSE)
  }
  p0 <- th[["p0"]]
  p1 <- th[["p1"]]
  p2 <- th[["p2"]]
  a <- p1 + p2 - 1
  z0 <- (1 - p0)^m
  z1 <- if (m == 1L) 1 else (1 - p0)^(m - 1)
  f1 <- m * a * z1
  f2 <- 1 - z0
  f3 <- -z0
  c(f1, f2, f3)
}


#' Effective regressor \eqn{\sqrt{\lambda(x)}\,\mathbf{f}(x)} for convex design code
#'
#' The information matrix in Huang et al. (2020) is
#' \eqn{\sum_i w_i \lambda(x_i) \mathbf{f}(x_i)\mathbf{f}(x_i)^\top}.
#' The same matrix equals \eqn{\sum_i w_i \mathbf{h}(x_i)\mathbf{h}(x_i)^\top} with
#' \eqn{\mathbf{h}(x)=\sqrt{\lambda(x)}\mathbf{f}(x)}, which matches the
#' regression form used by [calc_Dopt()] and [check_equivalence()] without an
#' `info_weight` argument.
#'
#' @param theta Nominal \eqn{(p_0,p_1,p_2)} for local optimality.
#' @param q Cost ratio in `[0, 1]`; use `0` for constant cost per test.
#'
#' @return A function `function(x)` returning a length-3 vector.
#'
#' @references
#' Yeh, C.-K., Wong, W. K., and Zhou, J. (2025). Single and multi-objective
#' optimal designs for group testing experiments. arXiv:2508.08445.
#'
#' Huang, S.-Y., Chen, Y.-H., and Wang, W. (2020). Optimal group testing designs
#' for estimating prevalence with imperfect tests. Journal of the Royal
#' Statistical Society Series C.
#'
#' @export
gt_huang2020_regressor <- function(theta, q = 0) {
  th <- parse_theta_p012(theta)
  if (length(q) != 1L || !is.finite(q) || q < 0 || q > 1) {
    stop("`q` must be a scalar in [0, 1].", call. = FALSE)
  }

  function(x) {
    lam <- gt_huang2020_lambda(x, th, q)
    fv <- gt_huang2020_f(x, th)
    sqrt(lam) * fv
  }
}
