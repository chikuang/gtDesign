#' Exact maximin design under fixed run size (Algorithm I + search)
#'
#' Rounds a maximin approximate design to an exact design with total run size
#' `n`, then applies the modified Step II search over nearby support points
#' (\eqn{x_i \pm} `n_index`) and selects the extension with the largest
#' minimum efficiency (MinEff) across `criteria`.
#'
#' @param approx_design A maximin/approximate design object with component
#'   `design` containing columns `point` and `weight` (e.g. output of
#'   [compute_maximin_design()]).
#' @param u Integer candidate pool sizes (same grid as used for
#'   `approx_design`).
#' @param theta Nominal \eqn{(p_0,p_1,p_2)}.
#' @param n Total run size (fixed sample size).
#' @param q_cost Cost ratio \eqn{q} in \eqn{c(x)=1-q+qx}.
#' @param loss_ref Named list of reference losses from single-objective designs
#'   (same scale as [compute_maximin_design()]).
#' @param criteria Character vector of criteria to evaluate MinEff, e.g.
#'   `c("D","A")` or `c("D","A","Ds")`.
#' @param opts Contrast options: `cVec_c` and/or `cVec_Ds` when needed.
#' @param n_index Half-width of extension window around each support point.
#' @param ... Reserved.
#'
#' @return A list with `design_round1`, `design_exact`, `delta`,
#'   `M_round1`, `M_exact`, `efficiencies` (named vector), `min_efficiency`,
#'   and `extension_table` (tibble sorted by decreasing `min_efficiency`).
#'
#' @export
round_gt_design_n_maximin <- function(approx_design,
                                      u,
                                      theta,
                                      n,
                                      q_cost,
                                      loss_ref,
                                      criteria,
                                      opts = list(),
                                      n_index = 2L,
                                      ...) {
  if (length(n) != 1L || !is.finite(n) || n < 1) {
    stop("`n` must be a positive scalar.", call. = FALSE)
  }
  n <- as.integer(round(n))
  if (n < 1L) {
    stop("`n` must be at least 1.", call. = FALSE)
  }

  d <- approx_design$design
  x <- as.integer(d$point)
  ww <- as.numeric(d$weight)
  if (any(ww < 0) || abs(sum(ww) - 1) > 1e-5) {
    stop("`approx_design` weights must be nonnegative and sum to 1.", call. = FALSE)
  }

  n_prime <- floor(n * ww)
  m <- n - sum(n_prime)
  if (m < 0L) {
    stop("Internal error: negative remaining runs after floor.", call. = FALSE)
  }
  design_round1 <- rbind(x, n_prime)
  dr1 <- rounding_prune_and_name_design(design_round1)

  M_round1 <- rounding_fim_from_counts(dr1[1, ], dr1[2, ], n, theta, q_cost)

  if (m == 0L) {
    ed <- exact_design_efficiency_maximin(M_round1, loss_ref, criteria, opts = opts)
    return(list(
      design_round1 = dr1,
      design_exact = dr1,
      delta = rounding_delta_matrix(dr1, dr1),
      M_round1 = M_round1,
      M_exact = M_round1,
      efficiencies = ed$efficiencies,
      min_efficiency = ed$min_efficiency,
      extension_table = tibble::tibble()
    ))
  }

  x_temp <- rounding_extend_pool_sizes(x, n_index, u)
  L <- length(x_temp)
  if (L == 0L) {
    stop("No extended support points for Step II search.", call. = FALSE)
  }

  des2 <- rounding_run_size_combinations(L, m)
  des2 <- des2[des2[, L + 1L] == m, , drop = FALSE]
  if (!nrow(des2)) {
    stop("No feasible Step II combinations for remaining runs.", call. = FALSE)
  }

  n_cases <- nrow(des2)
  min_eff <- numeric(n_cases)
  eff_list <- vector("list", n_cases)

  for (k in seq_len(n_cases)) {
    add_counts <- des2[k, seq_len(L), drop = TRUE]
    merged <- rounding_merge_floor_extension(design_round1, x_temp, add_counts)
    M_k <- rounding_fim_from_counts(merged[1, ], merged[2, ], n, theta, q_cost)
    ed_k <- exact_design_efficiency_maximin(M_k, loss_ref, criteria, opts = opts)
    min_eff[k] <- ed_k$min_efficiency
    eff_list[[k]] <- ed_k$efficiencies
  }

  ext_full <- cbind(des2, min_efficiency = min_eff)
  colnames(ext_full) <- c(
    rounding_extension_colnames(L, ncol(des2), x_temp),
    "min_efficiency"
  )
  ord <- order(min_eff, decreasing = TRUE)
  ext_sorted <- ext_full[ord, , drop = FALSE]

  best_ext <- ext_sorted[1L, seq_len(L), drop = TRUE]
  design_exact <- rounding_merge_floor_extension(design_round1, x_temp, best_ext)
  dex <- rounding_prune_and_name_design(design_exact)
  M_exact <- rounding_fim_from_counts(dex[1, ], dex[2, ], n, theta, q_cost)
  ed_best <- exact_design_efficiency_maximin(M_exact, loss_ref, criteria, opts = opts)

  list(
    design_round1 = dr1,
    design_exact = dex,
    delta = rounding_delta_matrix(dr1, dex),
    M_round1 = M_round1,
    M_exact = M_exact,
    efficiencies = ed_best$efficiencies,
    min_efficiency = ed_best$min_efficiency,
    extension_table = tibble::as_tibble(as.data.frame(ext_sorted, check.names = FALSE))
  )
}


#' Exact maximin design under fixed budget (Algorithm II + search)
#'
#' Applies cost-aware floor rounding, then modified Step II search over nearby
#' support points (\eqn{x_i \pm} `n_index`) and selects the extension with the
#' largest minimum efficiency (MinEff) across `criteria`.
#'
#' @inheritParams round_gt_design_n_maximin
#' @param C Total cost budget (sum of run costs \eqn{n_i c(x_i)}).
#' @param fix_zero_floor If `TRUE` (default), run MATLAB-style zero-fix after
#'   floor allocation.
#' @param repair_floor_budget If `TRUE` (default), cap floor design to satisfy
#'   budget by removing runs from the most expensive active pools.
#'
#' @return A list with `design_round1`, `design_exact`, `delta`, `M_round1`,
#'   `M_exact`, `efficiencies`, `min_efficiency`, `C_remaining` after floor,
#'   and `extension_table` (tibble sorted by decreasing `min_efficiency`).
#'
#' @export
round_gt_design_budget_maximin <- function(approx_design,
                                           u,
                                           theta,
                                           C,
                                           q_cost,
                                           loss_ref,
                                           criteria,
                                           opts = list(),
                                           n_index = 2L,
                                           fix_zero_floor = TRUE,
                                           repair_floor_budget = TRUE,
                                           ...) {
  if (length(C) != 1L || !is.finite(C) || C <= 0) {
    stop("`C` must be a positive scalar.", call. = FALSE)
  }

  d <- approx_design$design
  x <- as.integer(d$point)
  ww <- as.numeric(d$weight)
  if (any(ww < 0) || abs(sum(ww) - 1) > 1e-5) {
    stop("`approx_design` weights must be nonnegative and sum to 1.", call. = FALSE)
  }

  cxx <- vapply(x, function(xi) gt_huang2020_cost(xi, q_cost), numeric(1))
  n_prime <- floor(C * ww / cxx)
  if (isTRUE(fix_zero_floor)) {
    n_prime <- rounding_fix_zero_floor_counts(n_prime)
  }
  if (isTRUE(repair_floor_budget)) {
    n_prime <- rounding_repair_floor_budget(n_prime, cxx, C)
  }

  Cr <- max(0, round(as.numeric(C - sum(n_prime * cxx)), 4L))
  design_round1 <- rbind(x, n_prime)
  dr1 <- rounding_prune_and_name_design(design_round1)
  M_round1 <- rounding_fim_from_counts(dr1[1, ], dr1[2, ], C, theta, q_cost)

  x_temp <- rounding_extend_pool_sizes(x, n_index, u)
  L <- length(x_temp)
  cxx_temp <- if (L) {
    vapply(x_temp, function(xi) gt_huang2020_cost(xi, q_cost), numeric(1))
  } else {
    numeric(0)
  }

  des2 <- if (L) {
    rounding_budget_combinations(cxx_temp, Cr)
  } else {
    matrix(numeric(0), nrow = 0L, ncol = 0L)
  }

  if (!nrow(des2)) {
    ed <- exact_design_efficiency_maximin(M_round1, loss_ref, criteria, opts = opts)
    return(list(
      design_round1 = dr1,
      design_exact = dr1,
      delta = rounding_delta_matrix(dr1, dr1),
      M_round1 = M_round1,
      M_exact = M_round1,
      efficiencies = ed$efficiencies,
      min_efficiency = ed$min_efficiency,
      C_remaining = Cr,
      extension_table = tibble::tibble()
    ))
  }

  n_cases <- nrow(des2)
  min_eff <- numeric(n_cases)

  for (k in seq_len(n_cases)) {
    add_counts <- des2[k, seq_len(L), drop = TRUE]
    merged <- rounding_merge_floor_extension(design_round1, x_temp, add_counts)
    M_k <- rounding_fim_from_counts(merged[1, ], merged[2, ], C, theta, q_cost)
    min_eff[k] <- exact_design_efficiency_maximin(
      M_k,
      loss_ref,
      criteria,
      opts = opts
    )$min_efficiency
  }

  ext_full <- cbind(des2, min_efficiency = min_eff)
  colnames(ext_full) <- c(
    rounding_extension_colnames(L, ncol(des2), x_temp),
    "min_efficiency"
  )
  ord <- order(min_eff, decreasing = TRUE)
  ext_sorted <- ext_full[ord, , drop = FALSE]

  best_ext <- ext_sorted[1L, seq_len(L), drop = TRUE]
  design_exact <- rounding_merge_floor_extension(design_round1, x_temp, best_ext)
  dex <- rounding_prune_and_name_design(design_exact)
  M_exact <- rounding_fim_from_counts(dex[1, ], dex[2, ], C, theta, q_cost)
  ed_best <- exact_design_efficiency_maximin(M_exact, loss_ref, criteria, opts = opts)

  list(
    design_round1 = dr1,
    design_exact = dex,
    delta = rounding_delta_matrix(dr1, dex),
    M_round1 = M_round1,
    M_exact = M_exact,
    efficiencies = ed_best$efficiencies,
    min_efficiency = ed_best$min_efficiency,
    C_remaining = Cr,
    extension_table = tibble::as_tibble(as.data.frame(ext_sorted, check.names = FALSE))
  )
}
