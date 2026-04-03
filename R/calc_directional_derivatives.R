#' @keywords internal
#' @noRd
.e_opt_subspace_setup <- function(M, eig_tol = NULL) {
  if (is.null(eig_tol)) {
    eig_tol <- max(1e-10, sqrt(.Machine$double.eps))
  }
  ev <- eigen(M, symmetric = TRUE)
  vals <- ev$values
  vecs <- ev$vectors
  lambda_min <- min(vals)
  mult_idx <- which(abs(vals - lambda_min) <= eig_tol * (1 + abs(lambda_min)))
  if (!length(mult_idx)) {
    mult_idx <- which.min(vals)
  }
  r_star <- length(mult_idx)
  Q <- vecs[, mult_idx, drop = FALSE]
  list(lambda_min = lambda_min, Q = Q, r_star = r_star)
}


#' @keywords internal
#' @noRd
.e_directional_scalar <- function(fx, setup) {
  z <- crossprod(setup$Q, fx)
  as.numeric(sum(z^2) / setup$r_star - setup$lambda_min)
}


#' Directional derivatives on a finite candidate set
#'
#' @param u Candidate design points.
#' @param M Information matrix.
#' @param f Regression function returning a numeric vector.
#' @param criteria Character vector of criteria, e.g. `c("D", "A", "c", "E")`.
#' @param cVec Optional vector for c-optimality.
#' @param use_ginv Logical; if `TRUE`, use [MASS::ginv()] when `M` is singular.
#' @param eig_tol Tolerance for grouping eigenvalues with
#'   \eqn{\lambda_{\min}(\mathbf{M})} when computing **E**-optimal derivatives
#'   (multiplicity \eqn{r^*}). Passed to the eigenvalue comparison.
#'
#' @section E-optimality:
#' Implements \eqn{d_E} from Sec.\ 4 (Eq.\ (9)) of Yeh, Wong, and Zhou
#' (arXiv:2508.08445): with orthonormal \eqn{\mathbf{Q}(\mathbf{w})} spanning the
#' eigenspace of \eqn{\lambda_{\min}(\mathbf{M})} and multiplicity \eqn{r^*},
#' uses \eqn{\mathbf{B} = (1/r^*)\mathbf{I}_{r^*}} (trace one, positive
#' semi-definite), so
#' \deqn{d_E(u_i,\mathbf{w}) = \frac{1}{r^*}\|\mathbf{Q}^\top \mathbf{f}(u_i)\|^2
#'   - \lambda_{\min}(\mathbf{M}).}
#' When \eqn{r^*=1}, this is \eqn{(\mathbf{v}^\top \mathbf{f}(u_i))^2 -
#' \lambda_{\min}(\mathbf{M})} for a unit eigenvector \eqn{\mathbf{v}}.
#'
#' @return A named list of directional derivative vectors.
#' @export
calc_directional_derivatives <- function(u,
                                         M,
                                         f,
                                         criteria = c("D"),
                                         cVec = NULL,
                                         use_ginv = TRUE,
                                         eig_tol = NULL) {
  if (missing(u) || length(u) == 0) {
    stop("`u` must be a non-empty candidate set.", call. = FALSE)
  }
  if (missing(M) || !is.matrix(M)) {
    stop("`M` must be a matrix.", call. = FALSE)
  }
  if (missing(f) || !is.function(f)) {
    stop("`f` must be a function.", call. = FALSE)
  }

  criteria <- unique(criteria)
  criteria_low <- tolower(criteria)

  p <- ncol(M)

  needs_minv <- any(criteria_low %in% c("d", "a", "c"))
  Minv <- NULL
  if (needs_minv) {
    Minv <- tryCatch(
      solve(M),
      error = function(e) {
        if (!use_ginv) {
          stop("`M` is singular and `use_ginv = FALSE`.", call. = FALSE)
        }
        if (!requireNamespace("MASS", quietly = TRUE)) {
          stop("Package `MASS` is required for generalized inverse.", call. = FALSE)
        }
        MASS::ginv(M)
      }
    )
  }

  e_setup <- NULL
  if (any(criteria_low == "e")) {
    e_setup <- .e_opt_subspace_setup(M, eig_tol = eig_tol)
  }

  out <- list()

  for (crit in criteria_low) {
    dd <- vapply(u, function(x) {
      fx <- f(x)
      if (is.matrix(fx) && ncol(fx) == 1) {
        fx <- as.vector(fx)
      }
      fx <- as.numeric(fx)

      if (crit == "d") {
        as.numeric(t(fx) %*% Minv %*% fx - p)

      } else if (crit == "a") {
        as.numeric(t(fx) %*% Minv %*% Minv %*% fx - sum(diag(Minv)))

      } else if (crit == "c") {
        if (is.null(cVec)) {
          stop("`cVec` must be supplied for c-optimality.", call. = FALSE)
        }
        as.numeric((t(fx) %*% Minv %*% cVec)^2 -
                     as.numeric(t(cVec) %*% Minv %*% cVec))

      } else if (crit == "e") {
        .e_directional_scalar(fx, e_setup)

      } else if (crit == "ds") {
        stop("`Ds` directional derivative is not implemented yet.", call. = FALSE)

      } else {
        stop(sprintf("Unknown criterion: %s", crit), call. = FALSE)
      }
    }, numeric(1))

    nm <- switch(
      crit,
      d = "dD",
      a = "dA",
      c = "dc",
      ds = "dDs",
      e = "dE"
    )

    out[[nm]] <- dd
  }

  out
}
