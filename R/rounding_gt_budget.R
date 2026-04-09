# After moving runs to fix zero counts, total cost can exceed `C` because c(x)
# varies with pool size; drop runs from the most expensive active pools until
# the total cost is at most `C`.
#' @noRd
rounding_repair_floor_budget <- function(n_prime, cxx, C) {
  n_prime <- as.integer(n_prime)
  cxx <- as.numeric(cxx)
  spent <- sum(n_prime * cxx)
  eps <- 1e-8
  while (spent > C + eps && any(n_prime > 0L)) {
    idx <- which(n_prime > 0L)
    j <- idx[which.max(cxx[idx])]
    n_prime[j] <- n_prime[j] - 1L
    spent <- sum(n_prime * cxx)
  }
  if (spent > C + 1e-5) {
    stop(
      "Could not cap floor allocation to budget `C` (check `C` and support).",
      call. = FALSE
    )
  }
  n_prime
}


# Matches `GT_MO_rounding_budget.m` / `GT_Single_budget.m` redistribution: fixed
# `deficit_idx` from the initial floor, paired with `excess_idx(1)` each step.
#' @noRd
rounding_fix_zero_floor_counts <- function(n_prime) {
  n_prime <- as.integer(round(n_prime))
  deficit_idx <- which(n_prime == 0L)
  if (length(deficit_idx) == 0L) {
    return(n_prime)
  }
  excess_idx <- which(n_prime >= 2L)
  for (i in seq_along(deficit_idx)) {
    if (length(excess_idx) == 0L) {
      stop("Not enough excess runs to fix zero floor counts.", call. = FALSE)
    }
    from <- excess_idx[1L]
    to <- deficit_idx[i]
    n_prime[from] <- n_prime[from] - 1L
    n_prime[to] <- n_prime[to] + 1L
    if (n_prime[from] < 2L) {
      excess_idx <- excess_idx[-1L]
    }
  }
  n_prime
}


#' @noRd
rounding_extend_pool_sizes <- function(x, n_index, u) {
  x <- as.integer(x)
  u <- as.integer(u)
  n_index <- as.integer(n_index)
  ext <- integer(0)
  for (xi in x) {
    ext <- c(ext, seq.int(xi - n_index, xi + n_index))
  }
  sort(unique(ext[ext %in% u]))
}


#' @noRd
rounding_merge_floor_extension <- function(design_round1, x_temp, w_ext) {
  x <- rbind(as.numeric(design_round1[1, ]), as.numeric(design_round1[2, ]))
  x_temp <- as.numeric(x_temp)
  w_ext <- as.numeric(w_ext)
  all_keys <- sort(unique(c(as.numeric(x[1, ]), x_temp)))
  values <- numeric(length(all_keys))
  ix <- match(x[1, ], all_keys)
  values[ix] <- values[ix] + x[2, ]
  iy <- match(x_temp, all_keys)
  if (anyNA(iy)) {
    stop("Internal error: extended support alignment.", call. = FALSE)
  }
  values[iy] <- values[iy] + w_ext
  nz <- values > 0
  rbind(all_keys[nz], values[nz])
}


#' @noRd
rounding_fim_from_counts <- function(pool_sizes, counts, C, theta, q) {
  f <- gt_huang2020_regressor(theta, q)
  cxx <- vapply(pool_sizes, function(x) gt_huang2020_cost(x, q), numeric(1))
  w <- counts * cxx / C
  fim_matrix(pool_sizes, w, f)
}


# Drop support columns with zero runs; no column names on returned matrix.
#' @noRd
rounding_prune_and_name_design <- function(dm) {
  if (!ncol(dm)) {
    return(dm)
  }
  keep <- dm[2L, ] > 0
  if (!any(keep)) {
    stop("Rounded design has no positive run counts.", call. = FALSE)
  }
  out <- dm[, keep, drop = FALSE]
  colnames(out) <- NULL
  rownames(out) <- c("x", "n_prime")
  out
}


# Per pool size: Delta_n = n_exact - n_round1; only columns with nonzero change.
#' @noRd
rounding_delta_matrix <- function(dm_round1, dm_exact) {
  pools <- sort(unique(c(
    as.numeric(dm_round1[1, ]),
    as.numeric(dm_exact[1, ])
  )))
  n1 <- rep(0, length(pools))
  n2 <- rep(0, length(pools))
  for (i in seq_along(pools)) {
    x <- pools[i]
    i1 <- which(as.numeric(dm_round1[1, ]) == x)
    i2 <- which(as.numeric(dm_exact[1, ]) == x)
    if (length(i1)) {
      n1[i] <- as.numeric(dm_round1[2, i1[1L]])
    }
    if (length(i2)) {
      n2[i] <- as.numeric(dm_exact[2, i2[1L]])
    }
  }
  d <- n2 - n1
  keep <- d != 0
  if (!any(keep)) {
    out <- matrix(numeric(0), nrow = 2L, ncol = 0L)
    rownames(out) <- c("x", "Delta_n")
    return(out)
  }
  out <- rbind(pools[keep], d[keep])
  colnames(out) <- NULL
  rownames(out) <- c("x", "Delta_n")
  out
}


# Extension columns named by pool sizes (x_ext), then used_cost, remaining, loss.
#' @noRd
rounding_extension_colnames <- function(L, ncol_mat, x_ext = NULL) {
  if (L > 0L) {
    idx <- if (length(x_ext) == L) {
      as.character(as.integer(x_ext))
    } else {
      as.character(seq_len(L))
    }
  } else {
    idx <- character(0)
  }
  if (ncol_mat == L + 3L) {
    c(idx, "used_cost", "remaining", "loss")
  } else if (ncol_mat == L + 2L) {
    c(idx, "used_cost", "remaining")
  } else if (ncol_mat == L) {
    idx
  } else {
    as.character(seq_len(ncol_mat))
  }
}


#' @noRd
efficiency_exact_vs_approx <- function(M_app, M_ex, criterion, opts, p) {
  la <- scalar_loss_from_M(M_app, criterion, opts)
  le <- scalar_loss_from_M(M_ex, criterion, opts)
  key <- canon_crit_key(criterion)
  if (key == "D") {
    ## Avoid Inf - Inf => NaN when log-det losses are infinite (singular/near-singular M)
    if (!is.finite(la) || !is.finite(le)) {
      return(NA_real_)
    }
    val <- exp((la - le) / p)
    if (!is.finite(val)) {
      return(NA_real_)
    }
    val
  } else if (key == "E") {
    if (!is.finite(la) || !is.finite(le) || la == 0) {
      return(NA_real_)
    }
    le / la
  } else {
    if (!is.finite(la) || !is.finite(le) || le == 0) {
      return(NA_real_)
    }
    la / le
  }
}


#' Exact group-testing design under a fixed budget (Rounding Algorithm II)
#'
#' Implements the floor allocation, zero-fix, extended-support search, and
#' best-merge steps in Sec. 5.1 of Yeh, Wong, and Zhou (arXiv:2508.08445), as
#' in the reference MATLAB scripts `GT_Single_budget.m` /
#' `calc_budget_round_combinations.m`.
#'
#' @param approx_design Output of [calc_Dopt()], [calc_Aopt()], [calc_copt()],
#'   [calc_Eopt()], or [compute_design_SO()] on the same candidate set `u`.
#' @param u Integer candidate pool sizes (same grid as used for `approx_design`).
#' @param theta Nominal \eqn{(p_0,p_1,p_2)}.
#' @param C Total cost budget (sum of run costs \eqn{n_i c(x_i)}).
#' @param q_cost Cost ratio \eqn{q} in \eqn{c(x)=1-q+qx}.
#' @param criterion One of `"D"`, `"A"`, `"c"`, `"Ds"`, `"E"`.
#' @param opts Contrast options: `cVec_c` and/or `cVec_Ds` when needed (not used for `"E"`).
#' @param n_index Half-width of the extended support window (\eqn{x_i\pm} `n_index`).
#' @param fix_zero_floor If `TRUE` (default), redistribute runs so no support
#'   pool has zero runs after `floor`, as in `GT_MO_rounding_budget.m`. If
#'   `FALSE`, skip this step (as in `GT_Single_budget.m`, which only floors).
#' @param repair_floor_budget If `TRUE` (default), after zero-fix, remove runs
#'   from the most expensive active pools until \eqn{\sum_i n_i c(x_i) \le C},
#'   so `design_round1` never exceeds the budget. Zero-fix can increase total
#'   cost because \eqn{c(x)} varies by pool size; without this step the floor
#'   design can strictly overspend. Set to `FALSE` only to match reference
#'   MATLAB MO scripts that skip this repair (then overspend is possible).
#' @param ... Reserved.
#'
#' @return A list with `design_round1` (2-row matrix: pool sizes, floor counts;
#'   columns with zero runs are removed; no column names),
#'   `design_exact` (same format after merge),
#'   `delta` (2-row matrix: pool sizes `x` and `Delta_n` \eqn{= n_{\mathrm{exact}}
#'   - n_{\mathrm{round1}}}; only pool sizes with nonzero change are columns),
#'   `M_approx`,
#'   `M_exact`, `efficiency`, `loss_approx`, `loss_exact`, `C_remaining` after
#'   the floor step (with default `repair_floor_budget`, total cost is at most
#'   `C`): \eqn{\max\{0, \mathrm{round}(C - \sum_i n_i c(x_i), 4)\}}
#'   (never negative), `extension_table` ([tibble::tibble()] of tight extensions;
#'   first columns named by extended-support pool sizes, then
#'   \code{used_cost}, \code{remaining}, \code{loss} when present), and
#'   `criterion`.
#'
#' @export
round_gt_design_budget <- function(approx_design,
                                   u,
                                   theta,
                                   C,
                                   q_cost,
                                   criterion = c("D", "A", "c", "Ds", "E"),
                                   opts = list(),
                                   n_index = 2L,
                                   fix_zero_floor = TRUE,
                                   repair_floor_budget = TRUE,
                                   ...) {
  criterion <- match.arg(
    criterion,
    choices = c("D", "A", "c", "Ds", "E")
  )
  if (length(C) != 1L || !is.finite(C) || C <= 0) {
    stop("`C` must be a positive scalar.", call. = FALSE)
  }

  f <- gt_huang2020_regressor(theta, q_cost)
  p <- length(eval_regvec(u[1L], f))

  d <- approx_design$design
  x <- as.integer(d$point)
  ww <- as.numeric(d$weight)
  if (any(ww < 0) || abs(sum(ww) - 1) > 1e-5) {
    stop("`approx_design` weights must be nonnegative and sum to 1.", call. = FALSE)
  }

  cxx <- vapply(x, function(xi) gt_huang2020_cost(xi, q_cost), numeric(1))
  nn_i <- C * ww / cxx
  n_prime <- floor(nn_i)
  if (isTRUE(fix_zero_floor)) {
    n_prime <- rounding_fix_zero_floor_counts(n_prime)
  }
  if (isTRUE(repair_floor_budget)) {
    n_prime <- rounding_repair_floor_budget(n_prime, cxx, C)
  }

  spent_floor <- sum(n_prime * cxx)
  Cr <- max(0, round(as.numeric(C - spent_floor), 4L))

  design_round1 <- rbind(x, n_prime)

  x_temp <- rounding_extend_pool_sizes(x, n_index, u)
  L <- length(x_temp)
  cxx_all <- if (length(x_temp)) {
    vapply(x_temp, function(xi) gt_huang2020_cost(xi, q_cost), numeric(1))
  } else {
    numeric(0)
  }

  des2 <- if (length(cxx_all)) {
    rounding_budget_combinations(cxx_all, Cr)
  } else {
    matrix(numeric(0), nrow = 0L, ncol = 0L)
  }
  if (ncol(des2) > 0L) {
    colnames(des2) <- rounding_extension_colnames(L, ncol(des2), x_temp)
  }

  if (criterion == "c" && is.null(opts$cVec_c)) {
    stop("For criterion `c`, supply `opts$cVec_c`.", call. = FALSE)
  }
  if (criterion == "Ds" && is.null(opts$cVec_Ds)) {
    stop("For criterion `Ds`, supply `opts$cVec_Ds`.", call. = FALSE)
  }

  M_app <- fim_matrix(x, ww, f)

  loss_eval <- function(M) {
    scalar_loss_from_M(M, criterion, opts)
  }
  loss_app <- loss_eval(M_app)

  if (nrow(des2) == 0L) {
    warning(
      "No tight budget extensions for remaining Cr after floor = ",
      Cr,
      "; using floor allocation only."
    )
    pool_sizes <- as.numeric(design_round1[1, ])
    counts <- as.numeric(design_round1[2, ])
    M_ex <- rounding_fim_from_counts(pool_sizes, counts, C, theta, q_cost)
    ext_tbl <- if (ncol(des2) > 0L) {
      tibble::as_tibble(
        as.data.frame(des2, check.names = FALSE),
        .name_repair = "minimal"
      )
    } else {
      tibble::tibble()
    }
    dr1 <- rounding_prune_and_name_design(design_round1)
    return(list(
      design_round1 = dr1,
      design_exact = dr1,
      delta = rounding_delta_matrix(dr1, dr1),
      M_approx = M_app,
      M_exact = M_ex,
      efficiency = efficiency_exact_vs_approx(M_app, M_ex, criterion, opts, p),
      loss_approx = loss_app,
      loss_exact = loss_eval(M_ex),
      C_remaining = Cr,
      extension_table = ext_tbl,
      criterion = criterion
    ))
  }

  n_cases <- nrow(des2)
  losses <- numeric(n_cases)

  for (k in seq_len(n_cases)) {
    w_ext <- des2[k, seq_len(L), drop = TRUE]
    merged <- rounding_merge_floor_extension(design_round1, x_temp, w_ext)
    pool_sizes <- merged[1, ]
    counts <- merged[2, ]
    M_k <- rounding_fim_from_counts(pool_sizes, counts, C, theta, q_cost)
    losses[k] <- loss_eval(M_k)
  }

  ext_full <- cbind(des2, loss = losses)
  colnames(ext_full) <- rounding_extension_colnames(L, ncol(ext_full), x_temp)
  ord <- order(losses)
  ext_sorted <- ext_full[ord, , drop = FALSE]
  colnames(ext_sorted) <- colnames(ext_full)

  ext_tbl <- tibble::as_tibble(
    as.data.frame(ext_sorted, check.names = FALSE),
    .name_repair = "minimal"
  )

  best_ext <- ext_sorted[1L, seq_len(L), drop = TRUE]
  design_exact <- rounding_merge_floor_extension(design_round1, x_temp, best_ext)

  pool_sizes <- design_exact[1, ]
  counts <- design_exact[2, ]
  M_ex <- rounding_fim_from_counts(pool_sizes, counts, C, theta, q_cost)

  dr1 <- rounding_prune_and_name_design(design_round1)
  dex <- rounding_prune_and_name_design(design_exact)

  list(
    design_round1 = dr1,
    design_exact = dex,
    delta = rounding_delta_matrix(dr1, dex),
    M_approx = M_app,
    M_exact = M_ex,
    efficiency = efficiency_exact_vs_approx(M_app, M_ex, criterion, opts, p),
    loss_approx = loss_app,
    loss_exact = loss_eval(M_ex),
    C_remaining = Cr,
    extension_table = ext_tbl,
    criterion = criterion
  )
}
