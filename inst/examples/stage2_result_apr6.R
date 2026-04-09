library(dplyr)
library(tidyr)
library(tibble)
library(purrr)
library(ggplot2)
library(forcats)

# ============================================================
# Settings
# ============================================================
theta_true <- c(0.07, 0.93, 0.96)
Npop <- 10000
k_grid <- 1:20

theta_scenarios <- tribble(
  ~label,                          ~p0,   ~p1,   ~p2,

  # --------------------------------------------------
  # Truth
  # --------------------------------------------------
  "True theta",                    0.07,  0.93,  0.96,

  # --------------------------------------------------
  # Mild misspecification
  # --------------------------------------------------
  "Bad 1: mild joint",             0.09,  0.91,  0.94,
  "Bad 2: mild p0 only",           0.10,  0.93,  0.96,
  "Bad 3: mild p1 only",           0.07,  0.89,  0.96,
  "Bad 4: mild p2 only",           0.07,  0.93,  0.92,

  # --------------------------------------------------
  # Moderate misspecification
  # --------------------------------------------------
  "Bad 5: moderate joint",         0.12,  0.89,  0.93,
  "Bad 6: moderate p0 high",       0.15,  0.93,  0.96,
  "Bad 7: moderate p1 low",        0.07,  0.85,  0.96,
  "Bad 8: moderate p2 low",        0.07,  0.93,  0.88,

  # --------------------------------------------------
  # Severe / stress cases
  # --------------------------------------------------
  "Bad 9: severe joint",           0.15,  0.87,  0.90,
  "Bad 10: high Sp",               0.15,  0.87,  0.99,
  "Bad 11: low Se/Sp",             0.10,  0.80,  0.85,

  # --------------------------------------------------
  # New cases: prevalence underestimated
  # --------------------------------------------------
  "Bad 12: p0 too low",            0.04,  0.93,  0.96,
  "Bad 13: p0 far too low",        0.02,  0.93,  0.96,

  # --------------------------------------------------
  # New cases: prevalence overestimated more strongly
  # --------------------------------------------------
  "Bad 14: p0 very high",          0.20,  0.93,  0.96,
  "Bad 15: p0 extreme high",       0.25,  0.93,  0.96,

  # --------------------------------------------------
  # New cases: overly optimistic assay assumptions
  # --------------------------------------------------
  "Bad 16: optimistic Se/Sp",      0.07,  0.97,  0.99,
  "Bad 17: optimistic all",        0.04,  0.97,  0.99,

  # --------------------------------------------------
  # New cases: overly pessimistic assay assumptions
  # --------------------------------------------------
  "Bad 18: pessimistic Se/Sp",     0.07,  0.80,  0.85,
  "Bad 19: pessimistic all",       0.15,  0.80,  0.85,

  # --------------------------------------------------
  # New cases: crossed-direction errors
  # --------------------------------------------------
  "Bad 20: low p0 high Se/Sp",     0.03,  0.97,  0.99,
  "Bad 21: high p0 low Se/Sp",     0.18,  0.82,  0.88,
  "Bad 22: low p0 low Se/Sp",      0.03,  0.82,  0.88,
  "Bad 23: high p0 high Se/Sp",    0.18,  0.97,  0.99
)

theta_scenarios
# ============================================================
# Expected tests under imperfect Dorfman testing
# ============================================================
pi_group_imperfect <- function(k, theta) {
  p0 <- theta[1]
  p1 <- theta[2]
  p2 <- theta[3]

  pi_k <- p1 - (p1 + p2 - 1) * (1 - p0)^k
  pmin(pmax(pi_k, 1e-12), 1 - 1e-12)
}

expected_tests_dorfman_imperfect <- function(k, theta, Npop) {
  n_groups <- ceiling(Npop / k)
  pi_k <- pi_group_imperfect(k, theta)
  n_groups + Npop * pi_k
}

choose_best_k <- function(theta, k_grid, Npop) {
  tibble(
    k = k_grid,
    expected_tests = map_dbl(
      k_grid,
      ~ expected_tests_dorfman_imperfect(.x, theta = theta, Npop = Npop)
    )
  ) %>%
    arrange(expected_tests, k) %>%
    slice(1)
}

# ============================================================
# Oracle under true theta
# ============================================================
oracle_tbl <- tibble(
  k = k_grid,
  expected_tests = map_dbl(
    k_grid,
    ~ expected_tests_dorfman_imperfect(.x, theta = theta_true, Npop = Npop)
  )
) %>%
  arrange(expected_tests, k)

oracle_k <- oracle_tbl$k[1]
oracle_tests <- oracle_tbl$expected_tests[1]

print(oracle_tbl)
cat("Oracle k =", oracle_k, "\n")
cat("Oracle expected tests =", oracle_tests, "\n")

# ============================================================
# Build stage2_compare_tbl robustly
# ============================================================
stage2_compare_tbl <- pmap_dfr(
  theta_scenarios,
  function(label, p0, p1, p2) {
    theta_used <- c(p0, p1, p2)
    chosen <- choose_best_k(theta_used, k_grid = k_grid, Npop = Npop)

    actual_under_true <- expected_tests_dorfman_imperfect(
      k = chosen$k,
      theta = theta_true,
      Npop = Npop
    )

    tibble(
      label = label,
      p0 = p0,
      p1 = p1,
      p2 = p2,
      k_chosen = chosen$k,
      expected_tests_under_used_theta = chosen$expected_tests,
      actual_expected_tests_under_true = actual_under_true,
      oracle_k = oracle_k,
      oracle_expected_tests_under_true = oracle_tests,
      regret = actual_under_true - oracle_tests,
      pct_increase_vs_oracle = 100 * (actual_under_true - oracle_tests) / oracle_tests,
      type = ifelse(label == "True theta", "True", "Misspecified")
    )
  }
) %>%
  mutate(
    label = factor(label, levels = theta_scenarios$label)
  )

print(stage2_compare_tbl, n = Inf)

# ============================================================
# Build curve_tbl robustly
# ============================================================
curve_tbl <- pmap_dfr(
  theta_scenarios,
  function(label, p0, p1, p2) {
    theta_used <- c(p0, p1, p2)

    tibble(
      label = label,
      k = k_grid,
      expected_tests = map_dbl(
        k_grid,
        ~ expected_tests_dorfman_imperfect(.x, theta = theta_used, Npop = Npop)
      )
    )
  }
) %>%
  mutate(
    label = factor(label, levels = theta_scenarios$label)
  )

print(curve_tbl)

# ============================================================
# Plot 1: regret
# ============================================================
p_regret <- stage2_compare_tbl %>%
  mutate(label = fct_reorder(label, regret, .desc = FALSE)) %>%
  ggplot(aes(x = regret, y = label, colour = type, shape = type)) +
  geom_point(size = 3.5) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(
    title = "Extra essays needed from using misspecified parameter values in Stage II",
    x = "Extra expected tests under true theta",
    y = NULL
  ) +
  theme_bw(base_size = 13)

print(p_regret)

# ============================================================
# Plot 2: chosen k
# ============================================================
p_k <- stage2_compare_tbl %>%
  mutate(label = fct_reorder(label, k_chosen, .desc = FALSE)) %>%
  ggplot(aes(x = k_chosen, y = label, colour = type, shape = type)) +
  geom_point(size = 3.5) +
  geom_vline(xintercept = oracle_k, linetype = 2) +
  labs(
    title = "Chosen Stage II pool size under each parameter scenario",
    x = "Chosen pool size k",
    y = NULL
  ) +
  theme_bw(base_size = 13)

print(p_k)

# ============================================================
# Plot 3: actual expected tests under true theta
# ============================================================
p_actual <- stage2_compare_tbl %>%
  mutate(label = fct_reorder(label, actual_expected_tests_under_true, .desc = FALSE)) %>%
  ggplot(aes(x = actual_expected_tests_under_true, y = label, colour = type, shape = type)) +
  geom_point(size = 3.5) +
  geom_vline(xintercept = oracle_tests, linetype = 2) +
  labs(
    title = "Actual expected tests under true theta",
    x = "Expected tests evaluated at true theta",
    y = NULL
  ) +
  theme_bw(base_size = 13)

print(p_actual)

# ============================================================
# Plot 4: expected test curves in facets
# ============================================================
keep_labels <- c(
  "True theta",
  "Bad 2: mild p0 only",
  "Bad 5: moderate joint",
  "Bad 9: severe joint",
  "Bad 10: high Sp",
  "Bad 13: p0 far too low",
  "Bad 14: p0 very high",
  "Bad 20: low p0 high Se/Sp",
  "Bad 22: low p0 low Se/Sp",
  "Bad 23: high p0 high Se/Sp"
)


## Plot the regret

stage2_compare_tbl %>%
  filter(label %in% keep_labels) %>%
  ggplot(aes(
    x = regret,
    y = forcats::fct_reorder(label, regret),
    colour = label
  )) +
  geom_point(size = 3.5) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(
    title = "Regret from using misspecified parameter values in Stage II",
    x = "Extra expected tests relative to oracle",
    y = NULL,
    colour = "Scenario"
  ) +
  theme_bw(base_size = 13)
#
curve_tbl_sub <- curve_tbl %>%
  filter(label %in% keep_labels)
#
# p_curve_sub <- ggplot(
#   curve_tbl_sub,
#   aes(x = k, y = expected_tests, colour = label, group = label)
# ) +
#   geom_line(linewidth = 1) +
#   geom_point(size = 2) +
#   labs(
#     title = "Expected test curves for selected parameter scenarios",
#     x = "Pool size k",
#     y = "Expected number of tests",
#     colour = "Scenario"
#   ) +
#   theme_bw(base_size = 12)
#
# print(p_curve_sub)


library(dplyr)
library(ggplot2)
library(ggrepel)

# --------------------------------------------------
# only keep one true-theta curve
# --------------------------------------------------
true_curve_tbl <- tibble(
  k = k_grid,
  expected_tests_true = sapply(
    k_grid,
    expected_tests_dorfman_imperfect,
    theta = theta_true,
    Npop = Npop
  )
)

# --------------------------------------------------
# selected scenarios: use only their chosen k,
# then evaluate under true theta
# --------------------------------------------------
chosen_label_tbl <- stage2_compare_tbl %>%
  filter(label %in% keep_labels) %>%
  transmute(
    label,
    k = k_chosen,
    y = actual_expected_tests_under_true,
    label_text = paste0(label, "\n", "k = ", k_chosen)
  ) %>%
  arrange(k, label) %>%
  group_by(k) %>%
  mutate(
    x_plot = k + seq(-0.18, 0.18, length.out = n())
  ) %>%
  ungroup()

# --------------------------------------------------
# plot: one true-theta curve + labelled chosen k's
# --------------------------------------------------
p_curve_sub <- ggplot(true_curve_tbl, aes(x = k, y = expected_tests_true)) +
  geom_line(linewidth = 1.4, colour = "black") +
  geom_point(size = 2.8, colour = "black") +
  geom_point(
    data = chosen_label_tbl,
    aes(x = x_plot, y = y, colour = label),
    size = 3.8
  ) +
  geom_text_repel(
    data = chosen_label_tbl,
    aes(x = x_plot, y = y, label = label_text, colour = label),
    size = 5.2,
    box.padding = 0.35,
    point.padding = 0.25,
    segment.alpha = 0.7,
    show.legend = FALSE
  ) +
  scale_x_continuous(
    limits = c(1, 10),
    breaks = 1:10,
    labels = 1:10
  ) +
  labs(
    title = "Stage II choices evaluated under the true parameter values",
    subtitle = "Black curve: expected tests under true theta; labels show the chosen pool size k for each scenario",
    x = "Pool size k",
    y = "Expected number of tests",
    colour = "Scenario"
  ) +
  theme_bw(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 15),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "none"
  )

print(p_curve_sub)



p_regret_grouped <- ggplot(
  regret_group_tbl,
  aes(x = regret, y = fct_rev(pattern))
) +
  geom_point(
    aes(colour = point_col, size = point_size),
    show.legend = FALSE
  ) +
  geom_text(
    aes(x = label_x, label = sprintf("%.3f%%", pct_increase_vs_oracle)),
    hjust = 0,
    size = 4.8
  ) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(
    limits = c(-10, 520),
    breaks = c(0, 100, 200, 300, 400, 500)
  ) +
  scale_colour_identity() +
  scale_size_identity() +
  labs(
    title = "Stage II regret under parameter misspecification",
    x = "Extra expected tests relative to the oracle choice",
    y = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

print(p_regret_grouped)
