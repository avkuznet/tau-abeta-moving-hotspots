# Moving-hotspot solver

MATLAB implementation of the reaction–diffusion model in *Tau creates spatially
moving hotspots of amyloid-β secondary nucleation: a reaction–diffusion model*
(A. V. Kuznetsov). Base MATLAB only — no toolboxes. Tested on R2023a.

```matlab
cd matlab
run_all            % every number quoted in Sections 4 and 5 (a few minutes)
run_all('quick')   % short horizon, coarse grid; checks that the code runs
make_figures       % Figs. 2–9 of the manuscript, written to ../figures
```

## Files

| file | purpose |
|---|---|
| `hotspot_params.m` | the reference and calibrated parameter sets (Table 1) and the physical scales |
| `hotspot_solve.m`  | method-of-lines solver; `opts.qss` selects quasi-steady or time-dependent monomer |
| `hotspot_post.m`   | diagnostics: the three maxima, Damköhler numbers, screening lengths, Eq. (18) check |
| `run_all.m`        | reproduces the numerical results quoted in the manuscript |
| `make_figures.m`   | generates Figs. 2–9 (Fig. 1 is a schematic) |

## Model

Eqs. (11)–(15) on 0 ≤ x ≤ 1 with no-flux boundaries for T and m:

    T_t = D_T T_xx + r_T T (1 - T)
    m_t = D_m m_xx + lam (1 - m) - 2 k_E m P - nu (R_tau + R_2)
    P_t = R_tau + R_2 - mu_P P
    M_t = 2 k_E m P + nu (R_tau + R_2) - mu_M M
    R_tau = k_tau T m^2,  R_2 = k_2 A(M) m^2,  A(M) = M / (1 + (M/M_c)^(1/3))

Optional extensions, switched on by fields of the parameter struct `p`:

| fields | extension | manuscript |
|---|---|---|
| `D_S`, `k_conv` | diffusing oligomer intermediate, Eqs. (21)–(22) | Section 5.7, Fig. 8 |
| `k_off`, `release_surface` | reversible deposition, release from the whole mass or from the accessible shell | Section 5.8, Fig. 9 |
| `sig_lo` | site density that depends on the monomer level at which the fibril grew | Section 5.8, Fig. 9 |

## Parameter sets

`hotspot_params('reference')` is the dimensionless reference set of Table 1,
used for Sections 5.1–5.4.

`hotspot_params('calibrated', Lmm)` is the dimensionally mapped set used for
Figs. 7–9. The two tau parameters fix the units against Meisl et al. (2021):
one time unit is 0.50 months and, for `Lmm = 7.20`, one length unit is
7.20 mm. `lam` and `D_m` are then set from Cirrito et al. (2003) and Waters
(2010), and `k_E` so that the kinetic Damköhler number in tau-exposed tissue
is of the order of the 1.2–2.3 derived from Hong et al. (2011). For longer
segments `D_T` and `D_m` are rescaled with `Lmm` so that the physical
coefficients are unchanged, and the initial tau pulse keeps its physical
width of 0.288 mm.

## The quasi-steady reduction

With `opts.qss = true` the monomer balance Eq. (12) is replaced by its
quasi-steady form, Eq. (16), solved by damped Newton iteration on a sparse
tridiagonal Jacobian (at most three iterations in practice). The operator is
negative definite for `m > 0`, so the positive solution is unique.

The reduction is valid only where the time scales separate:

| parameter set | `lam/r_T` | max &#124;m_qss − m_full&#124; |
|---|---|---|
| `reference`  | 1   | 1.4e-1 |
| `calibrated` | 422 | 7.3e-4 |

The reference runs therefore use `qss = false`; the calibrated runs use the
quasi-steady form.

## Resolution

The calibrated runs hold the physical grid spacing at 0.05 mm, not the number
of grid points, so that the 0.66 mm tau front stays resolved on long segments.
Fitted over 5–15 months, the outer-maximum speed is 1.650 mm per month on 30,
40 and 50 mm segments, with identical positions; a 20 mm segment gives 1.56.

## Conventions used in the report

Positions "at grid resolution" are the grid point of largest R_2;
`hotspot_post` returns parabolic sub-grid positions. Speeds are least-squares
slopes over t = 10–30 (reference run) and 5–15 months (mapped runs).

## License

MIT; see `LICENSE`.
