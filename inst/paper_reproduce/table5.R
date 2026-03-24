# Reproduce Table 5 from Research_GroupTesting-1.pdf
# ---------------------------------------------------
# Table 5: Optimal approximate designs for Chlamydia with
# theta = (0.022, 0.92, 0.965), M = 15, q = 0.
#
# Run:
#   Rscript -e 'devtools::load_all(); source("inst/paper_reproduce/table5.R")'

options(scipen = 999)

if (!exists("calc_Dopt", mode = "function")) {
  stop("Load package first: devtools::load_all()", call. = FALSE)
}

fmt3 <- function(x) formatC(x, format = "f", digits = 3)

theta <- c(p0 = 0.022, p1 = 0.92, p2 = 0.965)
q_cost <- 0
u <- 1:15L
f <- gt_huang2020_regressor(theta, q_cost)
support_tol <- 1e-5

# Compute single-objective OADs
res_D  <- calc_Dopt(u, f, drop_tol = support_tol)
res_A  <- calc_Aopt(u, f, drop_tol = support_tol)
res_Ds <- calc_copt(u, f, cVec = c(1, 0, 0), drop_tol = support_tol)
res_c  <- calc_copt(u, f, cVec = c(0, 1, 1), drop_tol = support_tol)

extract_support <- function(des_obj, criterion_name, tol = support_tol) {
  d <- des_obj$design
  keep <- as.numeric(d$weight) > tol
  out <- data.frame(
    Criterion = criterion_name,
    `Pool size` = as.integer(d$point[keep]),
    Weight = as.numeric(d$weight[keep]),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  out[order(out$`Pool size`), , drop = FALSE]
}

supp_c <- extract_support(res_c, "c")
supp_D <- extract_support(res_D, "D")
supp_A <- extract_support(res_A, "A")
supp_Ds <- extract_support(res_Ds, "Ds")

print(paste("Optimal supports with weight >", support_tol))
support_out <- rbind(supp_c, supp_D, supp_A, supp_Ds)
support_out$Weight <- as.numeric(fmt3(support_out$Weight))
if (requireNamespace("knitr", quietly = TRUE)) {
  print(knitr::kable(support_out, format = "pipe", digits = 3))
} else {
  print(support_out, row.names = FALSE)
}

# Optional sanity checks against paper values (tolerance 0.002)
# paper <- data.frame(
#   Criterion = c("c", "D", "A", "Ds"),
#   `Pool size 1` = c(0.155, 0.333, 0.159, 0.173),
#   `Pool size 7` = c(0.519, 0.333, 0.517, 0.526),
#   `Pool size 15` = c(0.3260, 0.333, 0.324, 0.301),
#   check.names = FALSE,
#   stringsAsFactors = FALSE
# )
#
# for (i in seq_len(nrow(out))) {
#   diffs <- abs(as.numeric(out[i, 2:4]) - as.numeric(paper[i, 2:4]))
#   print(paste(
#     out$Criterion[i],
#     "max abs diff vs paper =",
#     formatC(max(diffs), format = "f", digits = 4)
#   ))
# }
