# Re-draw Figure 3 from paper (new version)
# ------------------------------------------
# Figure 3 in Research_GroupTesting-1.pdf:
# (a) MinEff vs n for q = 0, criterion D-A-Ds, M = 150
# (b) MinEff vs C for q = 0.2, criterion D-A-Ds, M = 150
# Dashed line: 1/t* from corresponding maximin approximate design.
#
# Run from project root:
#   Rscript -e 'devtools::load_all(); source("inst/examples/redraw_figure3.R")'

options(scipen = 999)

fmt <- function(x, digits = 4) formatC(x, format = "f", digits = digits)

if (!exists("calc_Dopt", mode = "function")) {
  stop("Load package first: devtools::load_all()", call. = FALSE)
}

theta <- c(p0 = 0.07, p1 = 0.93, p2 = 0.96)
u <- 1:150L
criteria <- c("D", "A", "Ds")
opts_ds <- list(cVec_Ds = c(1, 0, 0))

# Figure ranges requested:
# 10 <= n <= 100, 100 <= C <= 1000
n_grid <- seq(10L, 100L, by = 10L)
C_grid <- seq(100L, 1000L, by = 100L)

build_refs <- function(q_cost) {
  f <- gt_huang2020_regressor(theta, q_cost)
  res_d <- calc_Dopt(u, f, drop_tol = 1e-6)
  res_a <- calc_Aopt(u, f, drop_tol = 1e-6)
  res_ds <- calc_copt(u, f, cVec = c(1, 0, 0), drop_tol = 1e-6)

  loss_ref <- list(
    D = -res_d$value,
    A = res_a$value,
    Ds = res_ds$value
  )

  list(f = f, loss_ref = loss_ref)
}

compute_curve_n <- function(mm, loss_ref, q_cost, n_vals) {
  min_eff <- numeric(length(n_vals))
  for (i in seq_along(n_vals)) {
    out <- round_gt_design_n_maximin(
      approx_design = mm,
      u = u,
      theta = theta,
      n = n_vals[i],
      q_cost = q_cost,
      loss_ref = loss_ref,
      criteria = criteria,
      opts = opts_ds
    )
    min_eff[i] <- out$min_efficiency
  }
  data.frame(n = n_vals, MinEff = min_eff)
}

compute_curve_C <- function(mm, loss_ref, q_cost, C_vals) {
  min_eff <- numeric(length(C_vals))
  for (i in seq_along(C_vals)) {
    out <- round_gt_design_budget_maximin(
      approx_design = mm,
      u = u,
      theta = theta,
      C = C_vals[i],
      q_cost = q_cost,
      loss_ref = loss_ref,
      criteria = criteria,
      opts = opts_ds
    )
    min_eff[i] <- out$min_efficiency
  }
  data.frame(C = C_vals, MinEff = min_eff)
}

# Panel (a): q = 0, MinEff vs n
refs_q0 <- build_refs(q_cost = 0)
mm_q0 <- compute_maximin_design(
  u = u,
  f = refs_q0$f,
  loss_ref = refs_q0$loss_ref,
  criteria = criteria,
  opts = opts_ds,
  drop_tol = 1e-6
)
inv_t_q0 <- 1 / mm_q0$tstar
curve_n <- compute_curve_n(mm_q0, refs_q0$loss_ref, q_cost = 0, n_vals = n_grid)

# Panel (b): q = 0.2, MinEff vs C
refs_q02 <- build_refs(q_cost = 0.2)
mm_q02 <- compute_maximin_design(
  u = u,
  f = refs_q02$f,
  loss_ref = refs_q02$loss_ref,
  criteria = criteria,
  opts = opts_ds,
  drop_tol = 1e-6
)
inv_t_q02 <- 1 / mm_q02$tstar
curve_C <- compute_curve_C(mm_q02, refs_q02$loss_ref, q_cost = 0.2, C_vals = C_grid)

print(paste("1/t* for q=0:", fmt(inv_t_q0)))
print(paste("1/t* for q=0.2:", fmt(inv_t_q02)))
print("Head of MinEff vs n (q=0):")
print(utils::head(curve_n, 8))
print("Head of MinEff vs C (q=0.2):")
print(utils::head(curve_C, 8))

# Draw Figure 3-style plot with ggplot2
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package `ggplot2` is required to draw Figure 3.", call. = FALSE)
}

# Global ggplot theme settings (easy to update in one place)
ggplot2::theme_set(
  ggplot2::theme_minimal(base_size = 16) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "#E7E7E7"),
      panel.grid.major.y = ggplot2::element_line(color = "#E7E7E7"),
      axis.title = ggplot2::element_text(face = "bold", size = 18),
      axis.text = ggplot2::element_text(size = 14, color = "#222222"),
      plot.title = ggplot2::element_text(face = "bold", size = 18),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
)

plot_df_n <- data.frame(
  x = curve_n$n,
  MinEff = curve_n$MinEff,
  inv_tstar = inv_t_q0
)
plot_df_C <- data.frame(
  x = curve_C$C,
  MinEff = curve_C$MinEff,
  inv_tstar = inv_t_q02
)

make_plot <- function(df, x_label, q_label, inv_t, x_breaks) {
  ann_x <- min(df$x) + 0.62 * (max(df$x) - min(df$x))
  ann_y <- inv_t + 0.0008

  ggplot2::ggplot(df, ggplot2::aes(x = x, y = MinEff)) +
    ggplot2::geom_line(color = "#2A6FBB", linewidth = 1.1) +
    ggplot2::geom_point(color = "#1F4E79", fill = "#5FA8D3", shape = 21, size = 2.8, stroke = 0.5) +
    ggplot2::geom_hline(yintercept = inv_t, linetype = "dashed", linewidth = 0.8, color = "#555555") +
    ggplot2::annotate(
      "label",
      x = ann_x,
      y = ann_y,
      label = paste0("1/t* = ", fmt(inv_t, 3)),
      size = 5.0,
      fill = "white",
      color = "#333333"
    ) +
    ggplot2::labs(
      x = x_label,
      y = "MinEff",
      title = paste0("D-A-Ds exact design efficiency (q = ", q_label, ")")
    ) +
    ggplot2::scale_x_continuous(breaks = x_breaks)
}

g_n <- make_plot(
  df = plot_df_n,
  x_label = "n",
  q_label = "0",
  inv_t = inv_t_q0,
  x_breaks = c(10, 30, 50, 70, 90)
)
g_C <- make_plot(
  df = plot_df_C,
  x_label = "C",
  q_label = "0.2",
  inv_t = inv_t_q02,
  x_breaks = c(100, 300, 500, 700, 900)
)


gridExtra::grid.arrange(g_n, g_C, nrow = 1)
# plot_file_n <- "inst/examples/figure3_redraw_q0.png"
# plot_file_C <- "inst/examples/figure3_redraw_q02.png"
# ggplot2::ggsave(filename = plot_file_n, plot = g_n, width = 6.2, height = 4.6, dpi = 320)
# ggplot2::ggsave(filename = plot_file_C, plot = g_C, width = 6.2, height = 4.6, dpi = 320)

# print(paste("Saved Figure 3 redraw (q=0) to", plot_file_n))
# print(paste("Saved Figure 3 redraw (q=0.2) to", plot_file_C))

# Return data invisibly when sourced interactively
# invisible(list(
#   q0 = list(curve = curve_n, inv_tstar = inv_t_q0, plot = g_n),
#   q02 = list(curve = curve_C, inv_tstar = inv_t_q02, plot = g_C)
# ))
