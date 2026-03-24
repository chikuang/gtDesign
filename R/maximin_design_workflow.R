#' Maximin design with reference losses, equivalence, and eta (workflow)
#'
#' Runs the usual multi-step maximin pipeline: single-objective designs to build
#' `loss_ref`, [compute_maximin_design()], [calc_directional_derivatives()],
#' [calc_eta_weights_maximin()], and optionally [check_equivalence_maximin()].
#'
#' @details
#' Only criteria **D**, **A**, and **c** are supported end-to-end (directional
#' derivatives for **Ds** / **E** are not wired here). For **c**, supply
#' `opts = list(cVec_c = ...)` as in [compute_maximin_design()].
#'
#' @param u Candidate design points.
#' @param f Regression function (same as in single-objective and maximin calls).
#' @param criteria Character vector of criteria, subset of `c("D", "A", "c")`.
#' @param opts Named list passed to contrast-based criteria; use `cVec_c` for **c**.
#' @param drop_tol Passed to [calc_Dopt()], [calc_Aopt()], [calc_copt()].
#' @param q Parameter dimension for [calc_eta_weights_maximin()]; default is the
#'   length of `f(u[1])`.
#' @param tol_eta Tolerance passed to [calc_eta_weights_maximin()] (unless
#'   overridden in `eta_args`).
#' @param tol_equiv Tolerance for [check_equivalence_maximin()] when
#'   `check_equiv` is `TRUE`.
#' @param check_equiv If `TRUE`, run [check_equivalence_maximin()].
#' @param keep_reference_designs If `TRUE`, include full outputs of each
#'   single-objective solve in `reference_designs`.
#' @param info_weight Passed to [compute_maximin_design()].
#' @param solver Passed to all convex solves in this pipeline.
#' @param eta_args Named list of extra arguments for [calc_eta_weights_maximin()]
#'   (e.g. `list(solver = "SCS", complementary_slack = FALSE)`).
#' @param deriv_args Named list of extra arguments for [calc_directional_derivatives()].
#' @param make_figure If `TRUE`, draw the maximin equivalence figure via
#'   [plot_equivalence_maximin()] (opens a graphics device).
#' @param plot_args Named list of extra arguments for [plot_equivalence_maximin()]
#'   (merged over defaults that pass `tol = tol_equiv`).
#' @param ... Additional arguments passed to [compute_maximin_design()] (e.g.
#'   `support_tol`).
#'
#' @return A list with components: `loss_ref`, `maximin` (output of
#'   [compute_maximin_design()]), `directional_derivatives`, `eta`,
#'   `equivalence` (or `NULL` if `check_equiv` is `FALSE`), `q`, and optionally
#'   `reference_designs`. Plotting is a side effect when `make_figure` is `TRUE`.
#'
#' @export
maximin_design_workflow <- function(u,
                                    f,
                                    criteria,
                                    opts = list(),
                                    drop_tol = 1e-6,
                                    q = NULL,
                                    tol_eta = 1e-3,
                                    tol_equiv = 0.002,
                                    check_equiv = TRUE,
                                    keep_reference_designs = FALSE,
                                    info_weight = NULL,
                                    solver = "CLARABEL",
                                    eta_args = list(),
                                    deriv_args = list(),
                                    make_figure = FALSE,
                                    plot_args = list(),
                                    ...) {
  criteria <- standardize_criteria(criteria)
  bad <- setdiff(criteria, c("D", "A", "c"))
  if (length(bad) > 0L) {
    stop(
      "`maximin_design_workflow` supports only criteria D, A, and c. Got: ",
      paste(criteria, collapse = ", "),
      call. = FALSE
    )
  }

  p <- length(eval_regvec(u[1L], f))
  if (is.null(q)) {
    q <- as.numeric(p)
  }

  if ("c" %in% criteria && is.null(opts$cVec_c)) {
    stop("For criterion `c`, supply `opts$cVec_c` (length ", p, ").", call. = FALSE)
  }

  loss_ref <- list()
  reference_designs <- list()

  for (cr in criteria) {
    key <- canon_crit_key(cr)
    if (key == "D") {
      r <- calc_Dopt(u, f, solver = solver, drop_tol = drop_tol)
      loss_ref$D <- -r$value
      if (keep_reference_designs) {
        reference_designs$D <- r
      }
    } else if (key == "A") {
      r <- calc_Aopt(u, f, solver = solver, drop_tol = drop_tol)
      loss_ref$A <- r$value
      if (keep_reference_designs) {
        reference_designs$A <- r
      }
    } else if (key == "c") {
      cvec <- as.numeric(opts$cVec_c)
      r <- calc_copt(u, f, cVec = cvec, solver = solver, drop_tol = drop_tol)
      loss_ref$c <- r$value
      if (keep_reference_designs) {
        reference_designs$c <- r
      }
    }
  }

  maximin <- compute_maximin_design(
    u = u,
    f = f,
    loss_ref = loss_ref,
    criteria = criteria,
    opts = opts,
    info_weight = info_weight,
    solver = solver,
    ...
  )

  dd_base <- list(
    u = u,
    M = maximin$info_matrix,
    f = f,
    criteria = criteria
  )
  if ("c" %in% criteria) {
    dd_base$cVec <- as.numeric(opts$cVec_c)
  }
  directional_derivatives <- do.call(
    calc_directional_derivatives,
    modifyList(dd_base, deriv_args)
  )

  eta <- do.call(
    calc_eta_weights_maximin,
    modifyList(
      list(
        tstar = maximin$tstar,
        loss_ref = loss_ref,
        loss_model = maximin$loss,
        directional_derivatives = directional_derivatives,
        criteria = criteria,
        q = q,
        tol = tol_eta
      ),
      eta_args
    )
  )

  equivalence <- NULL
  if (isTRUE(check_equiv)) {
    equivalence <- check_equivalence_maximin(
      maximin,
      directional_derivatives,
      eta,
      tol = tol_equiv
    )
  }

  out <- list(
    loss_ref = loss_ref,
    maximin = maximin,
    directional_derivatives = directional_derivatives,
    eta = eta,
    equivalence = equivalence,
    q = q
  )
  if (keep_reference_designs) {
    out$reference_designs <- reference_designs
  }

  if (isTRUE(make_figure)) {
    do.call(
      plot_equivalence_maximin,
      modifyList(
        list(
          design_obj = maximin,
          directional_derivatives = directional_derivatives,
          eta = eta,
          criteria = criteria,
          tol = tol_equiv
        ),
        plot_args
      )
    )
  }

  out
}
