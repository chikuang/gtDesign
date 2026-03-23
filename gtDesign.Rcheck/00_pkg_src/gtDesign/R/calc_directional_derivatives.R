#' Directional derivatives on a finite candidate set
#'
#' @param u Candidate design points.
#' @param M Information matrix.
#' @param f Regression function returning a numeric vector.
#' @param criteria Character vector of criteria, e.g. `c("D", "A", "c")`.
#' @param cVec Optional vector for c-optimality.
#' @param use_ginv Logical; if `TRUE`, use [MASS::ginv()] when `M` is singular.
#'
#' @return A named list of directional derivative vectors.
#' @export
calc_directional_derivatives <- function(u,
                                         M,
                                         f,
                                         criteria = c("D"),
                                         cVec = NULL,
                                         use_ginv = TRUE) {
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

  p <- ncol(M)

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
      ds = "dDs"
    )

    out[[nm]] <- dd
  }

  out
}
