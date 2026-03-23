# gtDesign


- [Overview](#overview)
- [Statistical model (Huang et al. 2020; Sec. 2 of
  arXiv:2508.08445)](#statistical-model-huang-et-al-2020-sec-2-of-arxiv250808445)
- [Installation](#installation)
- [Package contents (exported)](#package-contents-exported)
- [Example: D-optimal design (Table 1, M = 61, q =
  0)](#example-d-optimal-design-table-1-m--61-q--0)
- [Example: A-optimal design (same theta, q =
  0)](#example-a-optimal-design-same-theta-q--0)
- [Example: Cost depending on pool size (q \>
  0)](#example-cost-depending-on-pool-size-q--0)
- [Example: c-optimality](#example-c-optimality)
- [Example: E-optimality via
  `compute_design_SO`](#example-e-optimality-via-compute_design_so)
- [Example: Equivalence theorem check
  (D-opt)](#example-equivalence-theorem-check-d-opt)
- [Maximin designs (brief)](#maximin-designs-brief)
- [References](#references)
- [License](#license)

[![](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**Authors**

- Chi-Kuang Yeh (Georgia State University) · [ORCID
  0000-0001-7057-2096](https://orcid.org/0000-0001-7057-2096)
- Weng Kee Wong (University of California, Los Angeles) · [ORCID
  0000-0001-5568-3054](https://orcid.org/0000-0001-5568-3054)
- Julie Zhou (University of Victoria)

## Overview

`gtDesign` is an R package for **locally optimal approximate designs**
on a **finite candidate set** for **group testing** (pooled testing).
Designs are found by **convex optimization** via
**[CVXR](https://cran.r-project.org/package=CVXR)** (typically with the
**CLARABEL** solver).

The main application is the **Huang et al. (2020)** model for
prevalence, sensitivity, and specificity with optional **cost that
depends on pool size**, as used in **Yeh, Wong, and Zhou (2025)**; see
[arXiv:2508.08445](https://arxiv.org/abs/2508.08445). The same interface
also supports **generic nonlinear design** whenever the (approximate)
information matrix is a sum of rank-one terms
$\sum_i w_i h(x_i) h(x_i)^\top$.

This README is generated from **`README.qmd`**. Regenerate with:

``` bash
quarto render README.qmd
```

## Statistical model (Huang et al. 2020; Sec. 2 of arXiv:2508.08445)

Let $\theta = (p_0, p_1, p_2)^\top$ denote **prevalence**,
**sensitivity**, and **specificity**, with $p_0 \in (0,1)$ and
$p_1, p_2 \in (0.5, 1]$. For a pool of size $x \in \{1,\ldots,M\}$, the
**positive response probability** is

$$
\pi(x \mid \theta) = p_1 - (p_1 + p_2 - 1)(1 - p_0)^x .
$$

**Standardized cost** (assay + enrollment) is

$$
c(x) = 1 - q + q x, \qquad q = \frac{q_1}{q_0 + q_1} \in [0,1],
$$

so $q = 0$ gives **constant cost per test** ($c(x) \equiv 1$), matching
the “equal cost” setting in Huang et al.

The **Fisher information matrix** for $\theta$ under independent
Bernoulli pool outcomes can be written as

$$
\mathbf{I}(\mathbf{w}, \theta) = \sum_{j=1}^{M} w_j \, \lambda(x_j) \, \mathbf{f}(x_j) \mathbf{f}(x_j)^\top,
$$

where $w_j$ are design weights on candidate pool sizes $x_j = j$,
$\sum_j w_j = 1$,

$$
\lambda(x) = \frac{1}{c(x)\,\pi(x\mid\theta)\{1-\pi(x\mid\theta)\}},
$$

and $\mathbf{f}(x) = \nabla_\theta \pi(x\mid\theta)$ is the gradient of
$\pi$ with respect to $\theta$ (see the paper for the explicit three
components).

### Regressor used in this package

Many functions (`calc_Dopt`, `calc_Aopt`, `check_equivalence`, …) assume
a **regression form**

$$
\mathbf{M}(\mathbf{w}) = \sum_j w_j \, \mathbf{h}(x_j) \mathbf{h}(x_j)^\top
$$

with **no separate `info_weight`**. That matches the Huang information
matrix if we set

$$
\mathbf{h}(x) = \sqrt{\lambda(x)} \, \mathbf{f}(x).
$$

The factory **`gt_huang2020_regressor(theta, q)`** returns the function
`function(x)` that computes $\mathbf{h}(x)$. You can still use
**`compute_design_SO(..., info_weight = ...)`** with raw $\mathbf{f}$
and $\lambda$ if you prefer the factored form.

**Local optimality:** all designs are **local** with respect to a
**nominal** $\theta^\ast$ (and fixed $q$, $M$).

## Installation

``` r
# install.packages("remotes")
remotes::install_github("chikuang/gtDesign")
```

``` r
# Alternative
devtools::install_github("chikuang/gtDesign")
```

**Requirements:** R ≥ 4.0 recommended; packages **CVXR**, **tibble**,
**MASS** (see `DESCRIPTION`). A conic solver reachable by CVXR
(e.g. **CLARABEL**) is required at runtime.

## Package contents (exported)

| Area | Functions |
|----|----|
| Huang (2020) building blocks | `gt_huang2020_pi`, `gt_huang2020_cost`, `gt_huang2020_lambda`, `gt_huang2020_f`, `gt_huang2020_regressor` |
| Classical approximate designs | `calc_Dopt`, `calc_Aopt`, `calc_copt` |
| Single-objective (D, A, Ds, c, E) | `compute_design_SO` |
| Maximin multi-objective | `compute_maximin_design`, `calc_eta_weights_maximin` |
| Equivalence | `check_equivalence`, `check_equivalence_maximin` |
| Plots | `plot_equivalence`, `plot_equivalence_maximin` |
| Directional derivatives | `calc_directional_derivatives`, `calc_multi_directional_derivative` |
| Toy 1-parameter perfect assay | `gt_homogeneous_pool_fisher`, `gt_sqrt_fisher_regressor_homogeneous` |

## Example: D-optimal design (Table 1, M = 61, q = 0)

Nominal $\theta^\ast = (0.07, 0.93, 0.96)^\top$ as in the paper;
D-optimal approximate design on $\{1,\ldots,61\}$ with $q=0$:

``` r
library(gtDesign)

theta <- c(p0 = 0.07, p1 = 0.93, p2 = 0.96)
M <- 61L
u <- seq_len(M)
f <- gt_huang2020_regressor(theta, q = 0)

res_d <- calc_Dopt(u, f, drop_tol = 1e-6)
res_d$design
```

      point    weight
    1     1 0.3333346
    2    17 0.3333309
    3    61 0.3333345

``` r
res_d$status
```

    [1] "optimal"

Support points $\{1, 17, 61\}$ with weights $1/3$ each match **Table 1**
(D-criterion, $q=0$) in arXiv:2508.08445.

## Example: A-optimal design (same theta, q = 0)

``` r
res_a <- calc_Aopt(u, f, drop_tol = 1e-6)
res_a$design
```

      point    weight
    1     1 0.4160718
    2    16 0.2133409
    3    61 0.3705873

## Example: Cost depending on pool size (q \> 0)

Larger $q$ puts more weight on enrollment in $c(x) = 1 - q + q x$. Use
`f <- gt_huang2020_regressor(theta, q)` with your chosen $q \in (0,1]$.

``` r
f_q <- gt_huang2020_regressor(theta, q = 0.2)
res_d_q <- calc_Dopt(u, f_q, drop_tol = 1e-6)
res_d_q$design
```

      point    weight
    1     1 0.3333358
    2    10 0.3333339
    3    61 0.3333304

## Example: c-optimality

Minimize the asymptotic variance of $\mathbf{c}^\top \hat{\theta}$ for a
user-specified $\mathbf{c}$ (length 3). Example
$\mathbf{c} = (0,1,1)^\top$ as in Table 1 of the paper:

``` r
c_vec <- c(0, 1, 1)
res_c <- calc_copt(u, f, cVec = c_vec, drop_tol = 1e-8)
subset(res_c$design, weight > 0.01)
```

      point    weight
    1     1 0.5212666
    4    56 0.1799662
    5    57 0.2987671

## Example: E-optimality via `compute_design_SO`

``` r
res_e <- compute_design_SO(
  u = u,
  f = f,
  criterion = "E",
  solver = "CLARABEL"
)
res_e$design
```

    # A tibble: 3 × 2
      point weight
      <int>  <dbl>
    1     1  0.415
    2    16  0.250
    3    61  0.335

## Example: Equivalence theorem check (D-opt)

``` r
eq_d <- check_equivalence(res_d, f = f, tol = 1e-4)
eq_d$max_violation
```

    [1] 2.21661e-05

``` r
eq_d$all_nonpositive
```

    [1] TRUE

``` r
plot_equivalence(eq_d, main = "D-opt: equivalence derivative")
```

<img src="README_files/figure-commonmark/unnamed-chunk-9-1.png"
data-fig-alt="Directional derivative curve for D-optimality; support points highlighted." />

## Maximin designs (brief)

**`compute_maximin_design`** implements the maximin efficiency
formulation in the paper (Sec. 4). You must supply **reference losses**
`loss_ref` from **single-objective** optimal designs (e.g. D, A, c, E)
at the same $\theta^\ast$, $q$, and candidate grid, then pass `criteria`
and optional contrast vectors in `opts`. See `?compute_maximin_design`
and `?calc_eta_weights_maximin` for details.

## References

1.  Yeh, C.-K., Wong, W. K., Zhou, J. (2025). Single and multi-objective
    optimal designs for group testing experiments. *arXiv* 2508.08445.
    <https://doi.org/10.48550/arXiv.2508.08445>

2.  Huang, S.-Y., Chen, Y.-H., Wang, W. (2020). Optimal group testing
    designs for estimating prevalence with imperfect tests. *Journal of
    the Royal Statistical Society Series C*.

3.  Pukelsheim, F. (2006). *Optimal Design of Experiments*. SIAM.

4.  Fedorov, V. V. (1972). *Theory of Optimal Experiments*. Academic
    Press.

## License

GPL-3 — see the `License` field in `DESCRIPTION`.
