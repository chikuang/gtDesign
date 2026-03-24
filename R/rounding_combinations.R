#' Enumerate nonnegative integer allocations (run-size rounding)
#'
#' All vectors \eqn{(z_1,\ldots,z_m)} with \eqn{z_j \ge 0}, \eqn{\sum_j z_j \le}
#' `max_run_size`, matching the supplementary **Rounding Algorithm I** search
#' (arXiv:2508.08445, Sec. 5.1).
#'
#' @param n_items Number of locations (columns).
#' @param max_run_size Upper bound on \eqn{\sum_j z_j}.
#'
#' @return Numeric matrix with columns \eqn{z_1,\ldots,z_m}, then total count,
#'   then remaining slots (`max_run_size -` total), sorted lexicographically
#'   with the last coordinate changing fastest (MATLAB `sortrows` order).
#'
#' @noRd
.sortrows_lastcoord_first <- function(mat) {
  if (!nrow(mat)) {
    return(mat)
  }
  nc <- ncol(mat)
  args <- lapply(rev(seq_len(nc)), function(j) mat[, j])
  ord <- do.call(order, args)
  mat[ord, , drop = FALSE]
}


#' @export
rounding_run_size_combinations <- function(n_items, max_run_size) {
  if (!is.numeric(n_items) || length(n_items) != 1L || n_items < 1L) {
    stop("`n_items` must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(max_run_size) || length(max_run_size) != 1L ||
        !is.finite(max_run_size) || max_run_size < 0) {
    stop("`max_run_size` must be a nonnegative finite scalar.", call. = FALSE)
  }

  n_items <- as.integer(n_items)
  max_run_size <- as.numeric(max_run_size)

  max_vals <- rep(max_run_size, n_items)
  ranges <- lapply(max_vals, function(m) 0:m)
  grid <- do.call(expand.grid, ranges)
  combos <- as.matrix(grid)
  total_count <- rowSums(combos)
  valid_idx <- total_count <= max_run_size + 1e-9
  valid_combos <- combos[valid_idx, , drop = FALSE]
  valid_counts <- total_count[valid_idx]
  remaining_slots <- max_run_size - valid_counts

  results <- cbind(valid_combos, valid_counts, remaining_slots)
  .sortrows_lastcoord_first(results)
}


#' Enumerate budget-feasible integer allocations (cost-aware rounding)
#'
#' All nonnegative integer vectors \eqn{\Delta} with \eqn{\Delta^\top c \le C_r},
#' then restricted to **tight** rows where remaining budget is insufficient for
#' any additional unit at minimum cost (Step II in **Rounding Algorithm II**,
#' arXiv:2508.08445, Sec. 5.1).
#'
#' @param cxx Per-unit costs \eqn{c(x_j)} for each extended support location
#'   (same length as the number of columns in the search grid).
#' @param Cr Remaining budget after floor allocation.
#'
#' @return Matrix whose first `length(cxx)` columns are counts, then used cost,
#'   then remaining budget. Sorted like the MATLAB reference implementation.
#'
#' @export
rounding_budget_combinations <- function(cxx, Cr) {
  cxx <- as.numeric(cxx)
  if (!length(cxx) || any(!is.finite(cxx)) || any(cxx <= 0)) {
    stop("`cxx` must be a positive finite cost vector.", call. = FALSE)
  }
  Cr <- max(0, as.numeric(Cr))
  if (!is.finite(Cr) || length(Cr) != 1L) {
    stop("`Cr` must be a finite scalar.", call. = FALSE)
  }

  n <- length(cxx)
  max_vals <- floor(Cr / cxx)
  if (any(max_vals < 0)) {
    stop("Invalid `Cr` relative to `cxx`.", call. = FALSE)
  }

  ranges <- lapply(max_vals, function(m) 0:m)
  grid <- do.call(expand.grid, ranges)
  combos <- as.matrix(grid)
  total_cost <- as.numeric(round(combos %*% cxx, 6L))
  Cr_rounded <- round(as.numeric(Cr), 6L)
  remaining_budget <- round(Cr_rounded - total_cost, 6L)
  valid_idx <- total_cost <= Cr_rounded + 1e-9
  valid_combos <- combos[valid_idx, , drop = FALSE]
  valid_costs <- total_cost[valid_idx]
  remaining_budget <- remaining_budget[valid_idx]

  min_item_cost <- min(cxx)
  tight_idx <- remaining_budget < min_item_cost - 1e-6
  if (!any(tight_idx)) {
    return(matrix(numeric(0), nrow = 0L, ncol = n + 2L))
  }

  results <- cbind(
    valid_combos[tight_idx, , drop = FALSE],
    valid_costs[tight_idx],
    remaining_budget[tight_idx]
  )
  .sortrows_lastcoord_first(results)
}
