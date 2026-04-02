# Reproduce Table 1 (approximate OADs, q = 0.8, M = 61)
# -------------------------------------------------------
# From README.qmd chunk `table1`.
#
# Run (from package root):
#   Rscript -e 'devtools::load_all(); source("inst/paper_reproduce/table1.R")'

options(scipen = 999)

if (!exists("calc_Dopt", mode = "function")) {
  stop("Load package first: devtools::load_all()", call. = FALSE)
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("Install dplyr for this script: install.packages(\"dplyr\")", call. = FALSE)
}


# User input --------------------------------------------------------------
theta <- c(p0 = 0.07, p1 = 0.93, p2 = 0.96)
u <- seq_len(150L)
q_val <- 0.8


f_q0 <- gt_huang2020_regressor(theta, q = q_val)

make_summary <- function(design, criterion) {
  design |>
    dplyr::filter(weight > 1e-4) |>
    dplyr::select(1, 2) |>
    dplyr::rename(group_size = 1, weight = 2) |>
    dplyr::mutate(weight = round(weight, 3)) |>
    dplyr::summarise(
      Criterion = criterion,
      `Support points` = paste(group_size, collapse = ", "),
      Weights = paste(weight, collapse = ", "),
      .groups = "drop"
    )
}

tab_summary <- dplyr::bind_rows(
  make_summary(calc_Dopt(u, f_q0)$design, "D-opt"),
  make_summary(calc_Aopt(u, f_q0)$design, "A-opt"),
  make_summary(calc_copt(u, f_q0, cVec = c(1, 0, 0))$design, "D_s-opt"),
  make_summary(calc_copt(u, f_q0, cVec = c(0, 1, 1))$design, "c-opt"),
  make_summary(calc_Eopt(u, f_q0)$design, "E-opt")
)

caption <- paste0(
  "Support points and weights for approximate optimal designs ",
  "for the Huang (2020) group testing model with q = ", q_val,
  " on the design space {", min(u), ", ..., ", max(u), "}."
)

print("Table 1")
if (requireNamespace("knitr", quietly = TRUE)) {
  print(knitr::kable(tab_summary, format = "pipe", caption = caption))
} else {
  print(tab_summary, row.names = FALSE)
}
