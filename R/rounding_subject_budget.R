#' Exact design under subject-count constraint via budget rounding
#'
#' Converts a subject cap (\code{n_subject_allow}) into a test budget
#' \code{C_tests = floor(n_subject_allow / sum(w_i x_i))} using the approximate
#' design, runs [round_gt_design_budget()], then selects the best rounded
#' extension under the subject-count constraint.
#'
#' Selection rule follows the application script logic:
#' among feasible candidates (\code{subjects <= n_subject_allow}), prefer using
#' more subjects (closest to the cap), then smaller loss. If none are feasible,
#' pick the minimum-excess candidate.
#'
#' @param approx_design Output of [calc_Dopt()], [calc_Aopt()], [calc_copt()],
#'   or [compute_design_SO()] on candidate set \code{u}.
#' @param u Integer candidate pool sizes.
#' @param theta Nominal \eqn{(p_0,p_1,p_2)}.
#' @param n_subject_allow Subject-count cap.
#' @param q_cost Cost ratio \eqn{q} in \eqn{c(x)=1-q+qx}.
#' @param criterion One of \code{"D"}, \code{"A"}, \code{"c"}, \code{"Ds"}, \code{"E"}.
#' @param opts Contrast options: \code{cVec_c} and/or \code{cVec_Ds} when needed.
#' @param n_index Half-width of extension window.
#' @param fix_zero_floor Passed to [round_gt_design_budget()].
#' @param repair_floor_budget Passed to [round_gt_design_budget()].
#' @param ... Reserved.
#'
#' @return A list with \code{C_tests}, \code{n_subject_allow},
#'   \code{n_subjects_used}, \code{design_round1}, \code{design_exact},
#'   \code{delta}, \code{M_approx}, \code{M_exact}, \code{loss_approx},
#'   \code{loss_exact}, \code{efficiency}, \code{criterion}, and
#'   \code{extension_table}.
#'
#' @export
round_gt_design_subject_budget <- function(approx_design,
                                           u,
                                           theta,
                                           n_subject_allow,
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
  if (length(n_subject_allow) != 1L || !is.finite(n_subject_allow) || n_subject_allow <= 0) {
    stop("`n_subject_allow` must be a positive scalar.", call. = FALSE)
  }

  d <- approx_design$design
  x <- as.numeric(d$point)
  ww <- as.numeric(d$weight)
  if (any(ww < 0) || abs(sum(ww) - 1) > 1e-5) {
    stop("`approx_design` weights must be nonnegative and sum to 1.", call. = FALSE)
  }

  mean_subjects_per_test <- sum(ww * x)
  if (mean_subjects_per_test <= 0) {
    stop("Invalid approximate design: non-positive mean pool size.", call. = FALSE)
  }
  C_tests <- floor(n_subject_allow / mean_subjects_per_test)
  if (C_tests < 1) {
    stop("Derived `C_tests < 1`; increase `n_subject_allow`.", call. = FALSE)
  }

  base <- round_gt_design_budget(
    approx_design = approx_design,
    u = u,
    theta = theta,
    C = C_tests,
    q_cost = q_cost,
    criterion = criterion,
    opts = opts,
    n_index = n_index,
    fix_zero_floor = fix_zero_floor,
    repair_floor_budget = repair_floor_budget
  )

  merge_counts <- function(base_design, x_ext, add_counts) {
    keys <- sort(unique(c(as.integer(base_design[1, ]), as.integer(x_ext))))
    vals <- setNames(numeric(length(keys)), as.character(keys))
    vals[as.character(base_design[1, ])] <- vals[as.character(base_design[1, ])] + as.numeric(base_design[2, ])
    vals[as.character(x_ext)] <- vals[as.character(x_ext)] + as.numeric(add_counts)
    keep <- vals > 0
    rbind(as.integer(names(vals)[keep]), as.numeric(vals[keep]))
  }

  chosen <- base$design_exact
  chosen_subjects <- sum(as.numeric(chosen[1, ]) * as.numeric(chosen[2, ]))
  ext_df <- as.data.frame(base$extension_table, check.names = FALSE)
  x_cols <- intersect(as.character(as.integer(u)), colnames(ext_df))

  if (nrow(ext_df) > 0L && length(x_cols) > 0L) {
    n_cases <- nrow(ext_df)
    cand_design <- vector("list", n_cases)
    cand_subjects <- numeric(n_cases)
    cand_loss <- as.numeric(ext_df$loss)

    for (i in seq_len(n_cases)) {
      add_counts <- as.numeric(ext_df[i, x_cols, drop = TRUE])
      cand_design[[i]] <- merge_counts(base$design_round1, as.integer(x_cols), add_counts)
      cand_subjects[i] <- sum(as.numeric(cand_design[[i]][1, ]) * as.numeric(cand_design[[i]][2, ]))
    }

    feasible <- cand_subjects <= n_subject_allow
    if (any(feasible)) {
      idx <- which(feasible)
      subject_gap <- n_subject_allow - cand_subjects[idx]
      pick <- idx[order(subject_gap, cand_loss[idx])[1L]]
    } else {
      pick <- which.min(cand_subjects - n_subject_allow)
    }
    chosen <- cand_design[[pick]]
    chosen_subjects <- cand_subjects[pick]
  }

  p <- length(eval_regvec(u[1L], gt_huang2020_regressor(theta, q_cost)))
  M_exact <- rounding_fim_from_counts(chosen[1, ], chosen[2, ], C_tests, theta, q_cost)
  loss_exact <- scalar_loss_from_M(M_exact, criterion, opts)
  eff <- efficiency_exact_vs_approx(base$M_approx, M_exact, criterion, opts, p)

  dr1 <- rounding_prune_and_name_design(base$design_round1)
  dex <- rounding_prune_and_name_design(chosen)

  list(
    C_tests = C_tests,
    n_subject_allow = as.numeric(n_subject_allow),
    n_subjects_used = as.numeric(chosen_subjects),
    design_round1 = dr1,
    design_exact = dex,
    delta = rounding_delta_matrix(dr1, dex),
    M_approx = base$M_approx,
    M_exact = M_exact,
    loss_approx = base$loss_approx,
    loss_exact = loss_exact,
    efficiency = eff,
    criterion = criterion,
    extension_table = base$extension_table
  )
}
