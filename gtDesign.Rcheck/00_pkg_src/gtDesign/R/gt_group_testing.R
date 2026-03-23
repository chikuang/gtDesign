#' Fisher information for one pool under homogeneous prevalence
#'
#' For independent subjects with infection probability \eqn{p \in (0,1)}, a pool
#' of size \eqn{m \geq 1} tests positive with probability
#' \eqn{q = 1 - (1-p)^m}. Under a Bernoulli pool outcome with perfect assay, the
#' Fisher information for \eqn{p} from a single pool is
#' \eqn{(q')^2 / (q(1-q))} where \eqn{q' = m(1-p)^{m-1}}.
#'
#' @param p Scalar infection probability in `(0, 1)`.
#' @param pool_size Positive integer pool size (number of subjects pooled).
#'
#' @return Scalar Fisher information for one such pool observation.
#' @export
gt_homogeneous_pool_fisher <- function(p, pool_size) {
  if (length(p) != 1L || !is.finite(p) || p <= 0 || p >= 1) {
    stop("`p` must be a scalar in (0, 1).", call. = FALSE)
  }
  m <- as.integer(pool_size)
  if (length(m) != 1L || !is.finite(m) || m < 1L) {
    stop("`pool_size` must be a positive integer.", call. = FALSE)
  }

  q <- 1 - (1 - p)^m
  dq_dp <- m * (1 - p)^(m - 1)
  dq_dp^2 / (q * (1 - q))
}


#' Regression vector factory for sqrt-Fisher (single-parameter prevalence)
#'
#' Returns a function `f(x)` suitable for \code{build_model_matrix()} and optimal
#' design routines, with `f(x)` equal to \eqn{\sqrt{I(p; x)}} where `x` is
#' interpreted as pool size. This matches the rank-one Fisher contribution
#' \eqn{f(x) f(x)^\top} for a scalar parameter.
#'
#' @inheritParams gt_homogeneous_pool_fisher
#'
#' @return A function with signature `function(x)` returning a length-one
#'   numeric vector.
#' @export
gt_sqrt_fisher_regressor_homogeneous <- function(p) {
  if (length(p) != 1L || !is.finite(p) || p <= 0 || p >= 1) {
    stop("`p` must be a scalar in (0, 1).", call. = FALSE)
  }

  function(x) {
    m <- as.integer(x)
    if (length(m) != 1L || !is.finite(m) || m < 1L) {
      stop("Design point must be a positive integer pool size.", call. = FALSE)
    }
    sqrt(gt_homogeneous_pool_fisher(p, m))
  }
}
