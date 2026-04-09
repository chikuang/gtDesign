u <- seq_len(150L)
cVec <- c(1, 0, 0)
theta = c(0.07, 0.93, 0.96)
f_q <- gt_huang2020_regressor(theta, q = 0.2)
out <- maximin_design_workflow(
  u = u,
  f = f_q,
  criteria = c("D", "A", "c"),
  opts = list(cVec_c = cVec),
  make_figure = TRUE
)
out$maximin$efficiency
out$maximin$value
out$maximin$tstar
out$maximin$design
