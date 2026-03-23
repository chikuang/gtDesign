#' Plot equivalence diagnostics for a maximin design
#'
#' @importFrom graphics abline grid par plot
#'
#' @param design_obj Output from [compute_maximin_design()].
#' @param directional_derivatives Named list of directional derivative vectors.
#' @param eta Named numeric vector of eta weights.
#' @param criteria Character vector of criteria to display.
#' @param main_prefix Optional prefix for panel titles.
#' @param line_width Line width for curves.
#' @param show_support Logical; if `TRUE`, show vertical lines at support points.
#' @param tol Numerical tolerance used for the horizontal reference line.
#' @param cex_lab,cex_axis,cex_main,mar Graphics parameters.
#'
#' @return Invisibly returns the object produced by [check_equivalence_maximin()].
#' @export
plot_equivalence_maximin <- function(design_obj,
                                     directional_derivatives,
                                     eta,
                                     criteria = NULL,
                                     main_prefix = NULL,
                                     line_width = 2,
                                     show_support = TRUE,
                                     tol = 1e-6,
                                     cex_lab = 1.2,
                                     cex_axis = 1,
                                     cex_main = 1.1,
                                     mar = c(5, 6.5, 3, 1)) {
  if (missing(design_obj) || !is.list(design_obj)) {
    stop("`design_obj` must be a valid maximin design object.", call. = FALSE)
  }
  if (!is.list(directional_derivatives) || length(directional_derivatives) == 0) {
    stop("`directional_derivatives` must be a non-empty named list.", call. = FALSE)
  }
  if (!is.numeric(eta) || is.null(names(eta))) {
    stop("`eta` must be a named numeric vector.", call. = FALSE)
  }

  eq_obj <- check_equivalence_maximin(
    design_obj = design_obj,
    directional_derivatives = directional_derivatives,
    eta = eta,
    tol = tol
  )

  u <- design_obj$candidates
  support_pts <- design_obj$design$point

  if (is.null(criteria)) {
    criteria <- names(eta)
  }

  criteria <- standardize_criteria(criteria)
  k <- length(criteria)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  if (k == 1) {
    par(mfrow = c(1, 2), mar = mar)
  } else if (k == 2) {
    par(mfrow = c(1, 3), mar = mar)
  } else if (k == 3) {
    par(mfrow = c(2, 2), mar = mar)
  } else {
    par(mfrow = c(ceiling((k + 1) / 2), 2), mar = mar)
  }

  label_map <- list(
    D = expression(d[D](u[i], w^"**")),
    A = expression(d[A](u[i], w^"**")),
    Ds = expression(d[D[s]](u[i], w^"**")),
    c = expression(d[c](u[i], w^"**")),
    E = expression(d[E](u[i], w^"**"))
  )
  name_map <- list(D = "dD", A = "dA", Ds = "dDs", c = "dc", E = "dE")

  dd_label <- function(cr) {
    key <- as.character(canon_crit_key(cr))
    out <- label_map[[key]]
    if (is.null(out)) {
      stop("Unsupported criterion: ", cr, call. = FALSE)
    }
    out
  }

  dd_name <- function(cr) {
    key <- as.character(canon_crit_key(cr))
    out <- name_map[[key]]
    if (is.null(out)) {
      stop("Unsupported criterion: ", cr, call. = FALSE)
    }
    out
  }

  panel_letters <- letters[seq_len(k + 1)]

  for (j in seq_along(criteria)) {
    cr <- criteria[j]
    nm <- dd_name(cr)

    if (is.null(directional_derivatives[[nm]])) {
      stop("Missing directional derivative field `", nm, "`.", call. = FALSE)
    }

    plot(
      u,
      directional_derivatives[[nm]],
      type = "l",
      lwd = line_width,
      xlab = expression(u[i]),
      ylab = dd_label(cr),
      main = paste0("(", panel_letters[j], ")")
    )
    abline(h = 0, lty = 2)

    if (show_support) {
      abline(v = support_pts, lty = 3)
    }

    grid()
  }

  plot(
    u,
    eq_obj$combined_directional_derivative,
    type = "l",
    lwd = line_width,
    xlab = expression(u[i]),
    ylab = "Combined directional derivative",
    main = paste0("(", panel_letters[k + 1], ")")
  )
  abline(h = 0, lty = 2)

  if (show_support) {
    abline(v = support_pts, lty = 3)
  }

  grid()

  invisible(eq_obj)
}
