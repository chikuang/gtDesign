#' Build the candidate regression matrix
#'
#' @param u Candidate design points (e.g. pool sizes or design labels).
#' @param f Regression function returning a numeric vector (e.g. score / sqrt
#'   Fisher contribution for group testing).
#'
#' @return A numeric matrix of dimension length(u) x p.
#' @noRd
build_model_matrix <- function(u, f) {
  validate_design_inputs(u, f)

  Fmat <- t(vapply(
    X = u,
    FUN = function(x) {
      out <- f(x)
      if (is.matrix(out) && ncol(out) == 1) {
        out <- as.vector(out)
      }
      as.numeric(out)
    },
    FUN.VALUE = as.numeric(f(u[1]))
  ))

  colnames(Fmat) <- paste0("beta", seq_len(ncol(Fmat)) - 1)
  Fmat
}


#' Validate inputs for optimal design problems
#'
#' @param u Candidate design points.
#' @param f Regression function returning a numeric vector.
#'
#' @return Invisibly returns TRUE if inputs are valid.
#' @noRd
validate_design_inputs <- function(u, f) {
  if (missing(u) || length(u) == 0) {
    stop("`u` must be a non-empty candidate set.", call. = FALSE)
  }

  if (!is.function(f)) {
    stop("`f` must be a function.", call. = FALSE)
  }

  test_val <- f(u[1])

  if (!is.numeric(test_val)) {
    stop("`f(u[i])` must return a numeric vector.", call. = FALSE)
  }

  if (is.matrix(test_val) && ncol(test_val) == 1) {
    test_val <- as.vector(test_val)
  }

  if (!is.vector(test_val)) {
    stop("`f(u[i])` must return a numeric vector.", call. = FALSE)
  }

  invisible(TRUE)
}


#' Compute the information matrix for a finite approximate design
#'
#' @param weights A numeric vector of weights.
#' @param Fmat Candidate regression matrix.
#'
#' @return Information matrix.
#' @noRd
info_matrix <- function(weights, Fmat) {
  if (length(weights) != nrow(Fmat)) {
    stop("Length of `weights` must match the number of candidate points.", call. = FALSE)
  }

  crossprod(Fmat, weights * Fmat)
}


#' Evaluate regression vector at one design point
#'
#' @param x A single design point.
#' @param f Regression function.
#'
#' @return Numeric column vector.
#' @noRd
eval_regvec <- function(x, f) {
  out <- f(x)
  if (is.matrix(out) && ncol(out) == 1) {
    out <- as.vector(out)
  }
  as.numeric(out)
}


#' Compute the weighted information matrix
#'
#' @param u Candidate design points.
#' @param w Numeric vector of design weights.
#' @param f Regression function returning a numeric vector.
#' @param info_weight Optional function returning a nonnegative scalar weight.
#'
#' @return A numeric information matrix.
#' @noRd
fim_matrix <- function(u, w, f, info_weight = NULL) {
  if (is.null(info_weight)) {
    info_weight <- function(x) 1
  }

  terms <- lapply(seq_along(u), function(i) {
    fi <- matrix(eval_regvec(u[i], f), ncol = 1L)
    ai <- as.numeric(info_weight(u[i]))

    if (length(ai) != 1L || !is.finite(ai) || ai < 0) {
      stop("`info_weight(x)` must return a finite nonnegative scalar.", call. = FALSE)
    }

    as.numeric(w[i]) * ai * tcrossprod(fi)
  })

  Reduce(`+`, terms)
}


#' Build the weighted information matrix expression for CVXR
#'
#' @param u Candidate design points.
#' @param w CVXR variable of weights.
#' @param f Regression function returning a numeric vector.
#' @param info_weight Optional function returning a nonnegative scalar weight.
#'
#' @return A CVXR matrix expression.
#' @noRd
fim_matrix_expr <- function(u, w, f, info_weight = NULL) {
  if (is.null(info_weight)) {
    info_weight <- function(x) 1
  }

  terms <- lapply(seq_along(u), function(i) {
    fi <- matrix(eval_regvec(u[i], f), ncol = 1L)
    ai <- as.numeric(info_weight(u[i]))

    if (length(ai) != 1L || !is.finite(ai) || ai < 0) {
      stop("`info_weight(x)` must return a finite nonnegative scalar.", call. = FALSE)
    }

    w[i] * ai * tcrossprod(fi)
  })

  Reduce(`+`, terms)
}


#' Extract the contrast vector needed by a criterion
#'
#' @param criterion Criterion name.
#' @param opts Named list of options.
#' @param p Dimension of the regression vector.
#'
#' @return A numeric column vector or NULL.
#' @noRd
get_contrast_vec <- function(criterion, opts, p) {
  criterion <- toupper(criterion)

  if (criterion == "DS") {
    cvec <- opts$cVec_Ds
    nm <- "cVec_Ds"
  } else if (criterion == "C") {
    cvec <- opts$cVec_c
    nm <- "cVec_c"
  } else {
    return(NULL)
  }

  if (is.null(cvec)) {
    stop("For criterion `", criterion, "`, supply `opts$", nm, "`.", call. = FALSE)
  }

  if (!is.numeric(cvec) || length(cvec) != p) {
    stop("`opts$", nm, "` must be a numeric vector of length ", p, ".", call. = FALSE)
  }

  matrix(as.numeric(cvec), ncol = 1L)
}


#' Compute the scalar loss from an information matrix
#'
#' @param M Information matrix.
#' @param criterion One of "D", "A", "Ds", "c", "E".
#' @param opts Named list of options.
#'
#' @return A numeric scalar loss.
#' @noRd
scalar_loss_from_M <- function(M, criterion, opts = list()) {
  cr <- toupper(as.character(criterion))

  safe_solve <- function(A, b = NULL) {
    tryCatch(
      {
        if (is.null(b)) solve(A) else solve(A, b)
      },
      error = function(e) {
        if (!requireNamespace("MASS", quietly = TRUE)) {
          stop(
            "Matrix inversion failed and package `MASS` is required for `ginv()`.",
            call. = FALSE
          )
        }
        G <- MASS::ginv(A)
        if (is.null(b)) G else G %*% b
      }
    )
  }

  if (cr == "D") {
    return(-as.numeric(determinant(M, logarithm = TRUE)$modulus))
  }

  if (cr == "A") {
    Minv <- safe_solve(M)
    return(sum(diag(Minv)))
  }

  if (cr == "DS") {
    cvec <- get_contrast_vec("Ds", opts, ncol(M))
    return(as.numeric(t(cvec) %*% safe_solve(M, cvec)))
  }

  if (cr == "C") {
    cvec <- get_contrast_vec("c", opts, ncol(M))
    return(as.numeric(t(cvec) %*% safe_solve(M, cvec)))
  }

  if (cr == "E") {
    evals <- eigen((M + t(M)) / 2, symmetric = TRUE, only.values = TRUE)$values
    return(-min(evals))
  }

  stop("Unknown criterion: ", criterion, call. = FALSE)
}


#' Standardize a criterion name
#'
#' @param cr A criterion name.
#'
#' @return A standardized criterion key.
#' @noRd
canon_crit_key <- function(cr) {
  key <- toupper(as.character(cr))

  if (key == "DS") {
    return("Ds")
  }
  if (key == "C") {
    return("c")
  }

  key
}


#' Standardize a vector of criterion names
#'
#' @param criteria A character vector of criterion names.
#'
#' @return A character vector of unique standardized criterion names.
#' @noRd
standardize_criteria <- function(criteria) {
  unique(vapply(criteria, canon_crit_key, character(1)))
}


#' Validate reference losses for maximin design
#'
#' @param loss_ref A named list of reference losses.
#' @param criteria A character vector of criteria.
#'
#' @return Invisibly returns TRUE.
#' @noRd
validate_loss_ref <- function(loss_ref, criteria) {
  if (!is.list(loss_ref)) {
    stop("`loss_ref` must be a named list.", call. = FALSE)
  }

  for (cr in criteria) {
    key <- canon_crit_key(cr)
    val <- loss_ref[[key]]

    if (is.null(val)) {
      stop(
        "`loss_ref$", key, "` is required when including criterion `", cr, "`.",
        call. = FALSE
      )
    }

    if (!is.numeric(val) || length(val) != 1L || !is.finite(val)) {
      stop("`loss_ref$", key, "` must be a finite scalar.", call. = FALSE)
    }

    # D-criterion loss is -log det(M); it can be negative when det(M) > 1.
    if (key != "D" && val <= 0) {
      stop(
        "`loss_ref$", key, "` must be a positive finite scalar for this criterion.",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}


#' Compute criterion-specific losses at a given information matrix
#'
#' @param M An information matrix.
#' @param criteria A character vector of criteria.
#' @param opts A named list of criterion-specific options.
#'
#' @return A named list of losses.
#' @noRd
compute_losses_at_M <- function(M, criteria, opts = list()) {
  criteria <- standardize_criteria(criteria)

  out <- list()
  for (cr in criteria) {
    key <- canon_crit_key(cr)
    out[[key]] <- scalar_loss_from_M(M, cr, opts)
  }

  out
}


#' Compute maximin efficiencies from reference and achieved losses
#'
#' @param loss_ref A named list of reference losses.
#' @param loss A named list of achieved losses.
#' @param criteria A character vector of criteria.
#' @param q Parameter dimension.
#'
#' @return A named numeric vector of efficiencies.
#' @noRd
compute_efficiencies_maximin <- function(loss_ref, loss, criteria, q) {
  criteria <- standardize_criteria(criteria)

  eff <- numeric(0)

  for (cr in criteria) {
    key <- canon_crit_key(cr)
    ref_val <- loss_ref[[key]]
    cur_val <- loss[[key]]

    if (is.null(ref_val) || is.null(cur_val)) {
      next
    }

    if (key == "D") {
      eff["D"] <- exp(((-cur_val) - (-ref_val)) / q)
    } else {
      eff[key] <- ref_val / cur_val
    }
  }

  eff
}
