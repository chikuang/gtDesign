#' Weighted multi-objective directional derivative
#'
#' @param dd_list A named list of directional derivative vectors.
#' @param eta A named numeric vector of weights.
#'
#' @return A numeric vector of combined directional derivative values.
#' @export
calc_multi_directional_derivative <- function(dd_list, eta) {
  if (!is.list(dd_list) || length(dd_list) == 0) {
    stop("`dd_list` must be a non-empty list.", call. = FALSE)
  }

  if (!is.numeric(eta) || is.null(names(eta))) {
    stop("`eta` must be a named numeric vector.", call. = FALSE)
  }

  out <- NULL

  for (nm in names(eta)) {
    dd_name <- switch(
      tolower(nm),
      "d" = "dD",
      "a" = "dA",
      "c" = "dc",
      "ds" = "dDs",
      "e" = "dE",
      stop(sprintf("Unknown eta name: %s", nm), call. = FALSE)
    )

    if (is.null(dd_list[[dd_name]])) {
      stop(sprintf("Directional derivative `%s` not found in `dd_list`.", dd_name),
           call. = FALSE)
    }

    if (is.null(out)) {
      out <- eta[[nm]] * dd_list[[dd_name]]
    } else {
      out <- out + eta[[nm]] * dd_list[[dd_name]]
    }
  }

  out
}
