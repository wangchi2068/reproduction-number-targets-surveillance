# Results summary: revision experiments (revise-experiments step, 2026-10-20)

All numbers come from `analysis/results.json`. Each record there gives its script, inputs, seed and mode. The file holds 510 numbers: 329 measured, 117 analytic and 64 simulated. Of these, 359 come from the revise-experiments step and 151 from round 3 (section 4); the 359 earlier values are unchanged.

Intervals are 95% percentile intervals from the two-way (state × week) pigeonhole bootstrap unless stated otherwise. Data are the NHSN HRD vintage of 2026-10-20 (51 jurisdictions); the data decision is in `analysis/data_decision.md`.

**Specification names:**
- **RA_rule4 (primary):** renewal-type ratio R_A in the rule-defined early-growth windows, four seasons (S = 4).
- **G_rule4:** the submitted weighted growth ratio G in the same windows.
- **G_cal3 (continuity):** G in the submitted calendar windows, S = 3.

## 1. Headline findings

1. **The submitted point estimates reproduce.** G_cal3 gives a = 0.0991 against the manuscript's 0.09906, with b/a = 22.6 and N = 1,267.
2. **Two-way intervals are much wider than state-clustered ones.**

   | Specification | a [95% CI] | b/a [95% CI] |
   |---|---|---|
   | RA_rule4 | 0.0150 [0.0100, 0.0208] | 37.4 [19.1, 71.0] |
   | G_rule4 | 0.0840 [0.0498, 0.1248] | 23.5 [8.4, 59.7] |
   | G_cal3 | 0.0991 [0.0569, 0.1382] | 22.6 [9.0, 52.2] |

   The bootstrap SD of a is 1.95× (RA_rule4) to 4.68× (G_cal3) its state-only value.
3. **About half of the intercept is a common national-week component.** Removing week fixed effects lowers a by 41% (R_A), 51% (G_rule4) and 52% (G_cal3). State-season effects add almost nothing.

   For RA_rule4:
   - a falls from 0.0150 to 0.0083 [0.0064, 0.0141];
   - the week ICC is 0.25.

   A floor that is shared across states at the same calendar week is not a per-state reporting floor. It is better read as a common, time-varying signal.
4. **The plateau is not identified against a power law.**
   - Parametric-bootstrap p-values: 0.11 (R_A), 0.44 (G_rule4), 0.42 (G_cal3).
   - The test never reaches 80% power, even when the true intercept equals the estimated a.
   - Leave-season-out cross-validation favours the affine model in 2 of 4 folds (R_A), 0 of 4 (G_rule4) and 3 of 3 (G_cal3).
   - Leave-state-out cross-validation favours it in 26, 19 and 28 of 51 folds respectively.
   - Neither form is preferred consistently.
5. **Theory and data disagree once the manuscript's own correction is applied.**
   - With Eq. b-delta, the theoretical b/a is 38.9–83.5 (G_cal3) or 40.7–87.5 (G_rule4).
   - The empirical b/a of about 23 lies below that range.
   - The empirical a is 3.1–6.9 times the theoretical a.
6. **Simulation shows the fitted intercept is not v_R/R̄².** Series were simulated with known R_t, generation interval and admission delay.
   - The fitted a tracks the noise-free limit of each statistic: the ratio of fitted to mean-field a has median 1.01 (G) and 1.02 (R_A).
   - Relative to v_R/R̄², that limit is inflated about 2.7× for G and shrunk about 0.48× for R_A.
   - b/a is recovered roughly 1:1 only on the statistic's own scale after the b-delta correction (median ratio 1.03 for G, 0.85 for R_A).
   - **Consequence:** the empirical a cannot be read directly as transmission variance. R_A's implied v_R is 0.0217, and its "excess over v_R" at 0.04/0.08/0.16 is negative (−0.013/−0.040/−0.095).
7. **Reporting completeness explains a negligible share of a.**

   | Specification | Δa [95% CI] | Equivalent at τ0 = 0.10? | Equivalent at τ0 = 0.05? | Seasons equivalent |
   |---|---|---|---|---|
   | G_rule4 | 0.0004 [−0.0036, 0.0033] | yes | no (zero power) | 3/4 |
   | RA_rule4 | 0.0003 [−0.0003, 0.0010] | yes | — | 4/4 |

   Margins were fixed before fitting.

   Calendar-matched May–October comparison, IQR of log G:

   | Period | IQR of log G |
   |---|---|
   | 2023, mandatory (reporting share ≈ 94%) | 0.48 |
   | 2024, voluntary (reporting share ≈ 46%) | 0.76 [0.48, 1.04] |
   | 2025, mandatory (reporting share ≈ 92%) | 0.48 |

   Low reporting coincides with more volatility. The conclusion must therefore be limited to "facility-reporting completeness in mandatory periods".

   Laboratory positivity proxy (44 states): it shares 0.23 [0.04, 0.42] of the R_A intercept.
8. **Seasons are few and the four-season estimate is uncertain.**

   | Specification | Pooled a [95% CI] | I² | Other |
   |---|---|---|---|
   | G_rule4 | 0.062 [0.007, 0.117] | 0.54 | — |
   | G_cal3 | 0.074 [−0.056, 0.203] | 0.78 | τ² Q-profile [0.0001, 0.115] |
   | RA_rule4 | 0.0125 [0.0022, 0.0227] | 0.47 | τ² ≈ 0 (1.9×10⁻⁵; Q-profile [0, 0.0006]) |

   Corollary 2 must use √(2/(S−1)): that gives 1.00 at S = 3 and 0.82 at S = 4. The manuscript's 0.82 and 0.71 are √(2/S).
9. **The cross-pathogen "replication" does not discriminate.** Rerunning the pure-reporting growth null gives ranges that already contain q1 < 1 < q2. For influenza in calendar windows, the ranges are q1 0.63–0.89 and q2 1.13–1.58. Over 40 tests (pathogen × window × criterion):
   - 21 are significant before correction, 10 after Holm and 18 after BH;
   - every Holm-significant test is an influenza ratio against 1 or an influenza excess > 0;
   - no ratio lies outside its null range after correction.

   COVID-19 and RSV show no significant result after Holm.
10. **Theory numbers reproduce, with three corrections.**
    - Reproduced: m× = 53.08, v_H = 0.0307, m×,H = 30.36, m0* = 72.86, T† = 4.484, B† = 775.0 (655.3–762.8 for m× 21.3–49.2), and 0.505 (at ξ = 33.06).
    - Corrected: the reporting lever term is 7.88, not 7.85; the lever total is 22.57, not 22.5; and the Corollary 2 values are √(2/(S−1)) (see item 8).

## 2. revision_plan items: what was done

Status: **done** = recomputed or new analysis in this workspace; **partial** = part done, gap stated; **rewrite** = text-only and left to the paper step; **deferred** = not done in this step.

| Item | Status | Evidence (results.json keys / files) | Notes |
|---|---|---|---|
| M1 ③ (n, ρ) bounds, Fig 1 curves | done | `s07_theory.V_floor_rho1_n*`, `n_over_rho_lever_ratio`, `mc_*`; Fig 1B–D | Monte Carlo maximum relative deviation of V is 0.040. ①②④ are rewrite (literature, estimator class) |
| M1 ④ T-period-mean MSE comparison | done | `s07_theory.rse_vR_T*`, `s05_seasons.design_*`; Fig 6D | — |
| M2 ① renewal R_t statistic | done (aggregate renewal ratio) | `s02_inference.RA_*`, `s01_data.rule4_*_RA_mid`; GI sensitivity in `s03` forest | R_A uses the log-ratio with a gamma GI. It is not a full EpiEstim posterior, so the paper must call it a "renewal-type ratio" |
| M2 ③ simulation with known R_t / GI / delay | done | `s08_sim.*`; Fig 3; `logs/s08_sim_rows.jsonl` | k is not identifiable and the delay is assumed |
| M3 ② hierarchical decomposition | done | `s03_decomposition.*_a_weekFE`, `*_a_both`, `*_share_week`, `*_icc_week`; Fig 4C | — |
| M3 ③ plateau vs power law with power | done | `*_plateau_contrast_p`, `*_power_exponent`, CV folds; Table S2 | Underpowered; no minimum detectable fraction ≤ 1 |
| M3 ④ b/a across specifications | done | Forest block in `s03`; Fig 4D | R_A: trimmed 24.7, drop-season 31.0–55.5, weekFE 59.2, peak 97.1. G_rule4: 13.3, 16.3–37.2, 42.7, 86.8 |
| M3 ① expectation with CV²(Λ) | rewrite + simulation | `s08_sim` | Derivation belongs to the paper |
| M4 corrected theory b/a | done | `s03_decomposition` theory block | 38.9–83.5 (G_cal3), 40.7–87.5 (G_rule4); a over theory 3.1–6.9× |
| M5 two-way clustering, one interval type | done | `s02_inference.*` (4 schemes + CGM), seeds and B in every record | Share of a > 0 in the bootstrap is 99.95–100%. TOST is redone with paired two-way draws |
| M6 ② pre-set TOST margins and power | done | `s04_reporting.*_tost_*`, power curves; Fig 5B–C, Table S3 | — |
| M6 ③ excess at v_R = 0.04/0.08/0.16 | done | `s03_decomposition` excess; `s04_reporting.RA_*excess*` | G_cal3 0.085/0.070/0.042 |
| M6 ④ positivity proxy | partial | `s04_reporting.lab_*`; Delphi data | Missing for 7 jurisdictions; the proxy is laboratory growth, not detection probability |
| M6 ① scope wording | rewrite | — | — |
| M7 ③ Holm/BH over pathogen × window × criterion | done | `s06_pathogens.*`, block `multiplicity`; Fig 6C, Table S6 | — |
| M7 ④ FluSurv-NET "insufficient", ILINet intervals | done | `s06_pathogens` blocks `flusurv_net`, `ilinet` | ILINet values reused from the vintage rerun (q1 = 0.687, q2 = 1.467) |
| M7 ② five classification variants | done (reused) | vintage `r4_results.json` `C_classification` / `C_posthoc` | Across mean, median, ratio-of-means and alternative schemes: 0–1 "replicated", 3–4 "partial", 3 "direction only", 0 "not replicated". This classification predates the A_all null and is superseded by Table S6 |
| M7 ① analysis-plan revision dates | rewrite | `analysis/data_decision.md` §3 (rule v1 → v2 disclosure) | — |
| M8 reproducible package | partial | `run_all.sh`, `code/README.md`, `results.json`, `data/DATA_MANIFEST.md` (hashes), figure scripts 1–6 | Zenodo DOI, generative-AI statement and the user README fixes belong to the paper and repo |
| M9 figures redrawn in English, Dark2, theme_bw, A–D | done | `figures/{en,zh}/fig1–fig6`, figure notes with 80–150-word captions | Fig 1A schematic is left to image-gen in the paper step |
| m1 0.505 general form | done | `s07_theory.xi`, `fisher_const_0505` | — |
| m2 Corollary 2 | done | `s05_seasons.relSE_bound_S*`, `s07_theory.rse_season_S*`; Fig 6B | — |
| m3 S = 3 vs S = 4, Q-profile | done | `s05_seasons.*`; Table S5, Fig 6A | — |
| m4 Proposition 7 | done | `s07_theory.B_dagger_*`, `T_dagger_phi08`, `m0_star_*` | — |
| m5 lever decomposition, equal budget | done | `s07_theory.lever_*`, `eqbudget_*`; Fig 2D | At 2× budget: volume 19.35 vs periods 11.98 (×10⁻³) |
| m6 Wald/UMP | deferred | — | The sim_theory code was not ported in this step. The paper should use "UMP at the attained level" wording or mid-p |
| m7 shrinkage example | deferred | — | The paper should label it as an ideal bound |
| m8 numeric consistency | done | `s04_reporting` regime N 6,009/420/4,077; ICC in `s03` | Full-vintage counts; 10,306 rows when cut at the original end date |
| m9 calendar-matched voluntary period; IV | partial | `s04_reporting.mayoct*`; Fig 5D, Table S4 | Calendar match done. The FGLS intercept is unstable in 2024 (negative), so SD and IQR are reported. The IV was not rerun; the paper should treat it as a consistency check |
| m12 D_t weights and zeros | partial | Forest variants (single-week 51.4, equal-weights 12.2); `s01_data.*_imputed_zero_rows` = 0 | R_A replaces the ad hoc weights with GI-derived weights. No separate GI × delay-derived G weights |
| m10, m11, m13 | rewrite | — | — |
| m14 | done | as M5 | — |

## 3. Caveats the paper must carry

- **Window rule.** Rule v2 was chosen after v1 failed in 2024–25. It is a forking path and must be disclosed; the calendar S = 3 results are kept for continuity.
- **Rule-window null range.** For influenza in rule windows the null range is taken from the calendar-window simulation, because no rule-window cells were simulated.
- **Fig 4B bin error bars** treat observations within a bin as independent, so they are too narrow. Inference is from the bootstrap.
- **Fig 6D and the design surface** use the apparent-scale R_A intercept. Because of item 6, this is illustrative.
- **The R_A generation interval** (mean 3.6 d, SD 1.6 d) is a literature value. Short and long GI variants change b/a by less than 3%.

## 4. 第三轮补充 (round 3: inversion, multiscale, reporting regime, returns, Parag contrast)

The numbers in this section come from `code/s09_inversion.py` and are merged into `results.json` under the prefixes `inversion.*`, `multiscale.*`, `reporting_regime.*`, `returns.*` and `parag.*`. Each record gives its script, the SHA-256 of its inputs (`inputs_sha256`), its seed and its mode. Seeds: simulation cells use `SeedSequence([20261020, 9, cell, chunk])`, the joint bootstrap 20261020 + 901, and the calibration-line draws 20261020 + 902.

**Reproducibility.** The command `cd experiment/code && python3 s09_inversion.py all && python3 s99_collect.py` was run from scratch with no cached parts or simulation rows. It took 370 s on 8 CPUs (simulation stage 329 s) and exited with code 0. All 462 numeric fields and both draw files matched the first run exactly, and the 359 earlier keys in `results.json` are byte-identical.

**Figures** (hidden working copy; the paper step chooses which to publish):
- `experiment/figures/{en,zh}/fig3_sim_inversion.{pdf,png}`, with notes and captions in `figures/figure_notes_fig3_inversion.md`;
- `experiment/figures/{en,zh}/fig6_design_nomogram.{pdf,png}`, with notes and captions in `figures/figure_notes_fig6_nomogram.md`.

All four PNGs were opened and checked. The decomposition figure was not redrawn, because Fig 6C carries the multiscale and J-pooling result.

### 4.1 Two-step inversion (A)

**Method.** The data-matched simulation (setting M0) reproduces the observed seasons, rule-window lengths, state sizes, national week shocks (share 0.43 against 0.41 measured for R_A) and admission delay. Each cell has 200 panels.
- **Calibration line.** E[â] = c0 + c1·v_R is fitted by WLS over v_R ∈ {0.02, 0.04, 0.06, 0.09}.
- **Inversion.** v_R^impl = (a_emp − c0)/c1.
- **Interval.** Each draw combines a joint two-way bootstrap draw of the empirical intercept (B = 2000) with a Monte Carlo draw of (c0, c1).

**Units.** The scale factor is κ = E[â]/(v_R/R̄²), so v_R = a·R̄²/κ. Here R̄ is the true mean reproduction number of the simulated windows, 1.176. The value 1.206 is the sample mean of R_A, which carries a convolution bias.

| Quantity | R_A | G | Key |
|---|---|---|---|
| v_R^impl, M0 (95% CI) | 0.040 (0.026, 0.056) | 0.043 (0.025, 0.064) | `inversion.M0_{RA,G}_vR_implied` |
| Difference R_A − G | −0.003 (−0.012, 0.005) | | `inversion.M0_diff_vR_RA_minus_G` |
| Ratio R_A / G | 0.93 (0.78, 1.18) | | `inversion.M0_ratio_vR_RA_over_G` |
| κ fitted / noise-free limit | 0.516 / 0.493 | 2.729 / 2.745 | `inversion.M0_{RA,G}_kappa_fit`, `_kappa_inf` |
| κ, s08 design (bridge) | 0.485 | 2.698 | `inversion.bridge_kappa_inf_{RA,G}` |
| Simulated vs empirical ICC | 0.258 vs 0.253 | 0.310 vs 0.310 | `inversion.M0_{RA,G}_sim_icc`, `multiscale.{RA,G}_rule4_icc_week` |

- **Settings sensitivity** (10 data-matched settings, figure 3B/3C). Across settings v_R^impl ranges from 0.013 to 0.063 (`inversion.vR_implied_range_across_settings`).
  - The true generation interval matters most: a shorter GI gives 0.013/0.014 and a longer GI 0.058/0.063 (R_A/G).
  - A synchronous reporting shock is the only setting in which the two statistics disagree (0.019 vs 0.029).
  - The two statistics agree in 9 of 10 settings (`inversion.n_settings_agree`).
  - v_R^impl falls inside 0.04–0.08 in 4 settings for R_A and 7 for G (`inversion.n_settings_{RA,G}_in_band`).
- **Calibrated b/a** over the 8 calibration cells. The pooled E[b]/E[a] ranges are 13.9–31.1 for R_A and 8.0–18.8 for G (`inversion.ba_calibrated_{RA,G}_range`).
  - The empirical values, 37.4 (R_A) and 23.5 (G), both lie above these ranges.
  - They lie inside the 95% per-panel sampling range (R_A 8.8–42.7, G 3.6–28.6; `inversion.ba_calibrated_{RA,G}_replicate_q95`), at its 93rd and 91st percentiles (`inversion.ba_empirical_{RA,G}_percentile_in_sim`).

### 4.2 Multiscale decomposition (B)

| Quantity | R_A | G | Key |
|---|---|---|---|
| National-week share of a (95% CI) | 0.41 (0.12, 0.59) | 0.51 (0.21, 0.67) | `multiscale.{RA,G}_rule4_share_week` |
| ICC of weekly deviations | 0.253 | 0.310 | `multiscale.{RA,G}_rule4_icc_week` |
| SD of â: state / week / two-way bootstrap | 0.0014 / 0.0021 / 0.0028 | 0.0084 / 0.0158 / 0.0193 | `multiscale.{RA,G}_rule4_a_sd_{state,week,twoway}` |
| Floor share s²_week/V1 | 0.23 | 0.29 | `multiscale.{RA,G}_rule4_pool_floor_share` |
| V_J/V1 at J = 5 / 10 / 51 | 0.39 / 0.31 / 0.25 | 0.43 / 0.36 / 0.31 | `multiscale.{RA,G}_rule4_pool_VJ_over_V1_J{5,10,51}` |
| J at which V_J is within 10% of the floor | 33 | 25 | `multiscale.{RA,G}_rule4_pool_J_within10pct_floor` |

Week clustering contributes more uncertainty than state clustering (SD ratio 1.46 for R_A, 1.87 for G; `multiscale.*_a_sd_ratio_week_over_state`). Pooling states removes the state component; the national-week component sets a floor at about a quarter of single-state variance.

### 4.3 Mandatory and voluntary reporting (C)

All contrasts in this subsection are observational. Reporting regime coincides with calendar time and with other changes in the 2024 off-season. The confounding checks are therefore matching on admissions D, matching on calendar month, a 2025/2023 placebo, and two further pathogens. None of them identifies a causal effect.

- **The "17%" figure** is the pre-set TOST equivalence margin expressed as a share of a. At v_R = 0.04 it is 0.17 (0.12, 0.25) for R_A and 0.12 (0.08, 0.20) for G (`reporting_regime.{RA,G}_rule4_margin_share_of_a_vR0.04`). The margin is unchanged at v_R^impl.
- **Measured contribution of reporting to a.** The measured Δa share is 0.020 (−0.024, 0.061) for R_A and 0.004 (−0.044, 0.040) for G (`reporting_regime.{RA,G}_rule4_delta_a_share`), with one-sided upper bounds of 0.053 and 0.035.
- **The "58%" figure** is the IQR ratio of log G, May–Oct 2024 over 2023: 1.58 (0.95, 2.23) (`reporting_regime.iqr_ratio_2024_over_2023`). The interval covers 1.
  - Against 2025 the ratio is 1.57 (0.95, 2.39).
  - The placebo 2025/2023 is 1.01 (0.72, 1.25) (`reporting_regime.iqr_ratio_2025_over_2023_placebo`).
  - D-matched ratios are 1.68 and 1.58 (`reporting_regime.iqr_ratio_2024_over_*_Dmatched`).
  - Month-matched ratios have median 1.33, with 9 of 12 above 1 (`reporting_regime.iqr_ratio_month_matched_*`).
  - COVID-19 gives 1.62 against 2023 and 1.18 against 2025; RSV gives 1.36 against 2025 (`reporting_regime.covid_*`, `rsv_*`).
- **Participation.** In the first mandatory period the median facility reporting share per state-week is 94.4% (mean 94.0%). 46.7% of state-weeks are at or above 95% and 84.6% at or above 90% (`reporting_regime.regime1_*`).
  - The voluntary period has a median of 35.3% (`reporting_regime.regime2_state_week_median`).
  - The second mandatory period has a median of 92.9% (`reporting_regime.regime3_state_week_median`).
  - Within the rule windows the median is 93.7% (`reporting_regime.rule4_windows_coverage_median`).

### 4.4 Diminishing returns (D)

- **Balance scale.** m× = 37.4 (19.1, 71.0) for R_A and 23.5 (8.4, 59.7) for G (`returns.{RA,G}_rule4_m_cross`).
- **Precision share.** The share of attainable precision, m/(m + m×), is 33%, 50%, 67% and 80% at m = 0.5, 1, 2 and 4 m× (`returns.precision_fraction_at_*_mcross`).
- **Gain from one doubling** of m starting from those points: variance cuts of 33%, 25%, 17% and 10% (`returns.variance_cut_doubling_from_*_mcross`).
- **Ninety per cent** of attainable precision needs m = 9 m× (`returns.m_for_90pct_precision_over_mcross`).
- **Optimal budget allocation.** With cost ratio c0/c1 = 10, m* = 19.4 (13.8, 26.6) for R_A and 15.3 for G. With c0/c1 = 100, m* = 61.2 and 48.4 (`returns.{RA,G}_rule4_m_star_c{10,100}`).
- **Rows above m×.** 55% of R_A rows and 65% of G rows lie above m× (`returns.{RA,G}_rule4_share_rows_above_mcross`).

### 4.5 Contrast with Parag et al. 2022 (E)

Under Poisson renewal the relative variance of R is R/m, with no floor; reporting only rescales m. The model here is v_R + A/m. On the latent scale the two variances differ by a factor of 1.5 at m = 13, 2 at m = 28, 5 at m = 116 and 10 at m = 262 (`parag.latent_m_ratio*`).

On the apparent scale of R_A:
- The fitted b is 1.25 times the Poisson value b_P = 0.446 (`parag.RA_fitted_b_over_poisson_b`, `parag.RA_poisson_b`).
- The ratio reaches 2 at m = 22 and 10 at m = 261 (`parag.RA_apparent_m_ratio*`).

The two models therefore agree at small counts and diverge near m×.

### 4.6 User numbers and measured values

| User's number (plan r3) | Measured value | Key | Verdict |
|---|---|---|---|
| G amplified 2.68× | 2.68 (s08 design); 2.73 in the data-matched simulation | `s08_sim.G_ainf_over_nominal_median`, `inversion.M0_G_kappa_fit` | Supported; cite the data-matched value |
| R_A compressed to 0.484× | 0.484 (s08 design); 0.516 in the data-matched simulation | `s08_sim.RA_ainf_over_nominal_median`, `inversion.M0_RA_kappa_fit` | Supported |
| v_R^impl ≈ 0.0451 (both statistics) | R_A 0.040 (0.026, 0.056); G 0.043 (0.025, 0.064) | `inversion.M0_{RA,G}_vR_implied` | Replace with ≈ 0.04. The user's 0.0451 uses the sample mean 1.206 for R̄. The written G formula 0.0840/(2.68 × 1.206²) gives 0.0215 and has the wrong units. The correct form is a·R̄²/κ |
| Theory b/a 13.4–28.7 contains the empirical values | Earlier figure: 13.4–28.7 for G only (`review/checks/revised_checks.json`, `C_derived.G_rule_theory_ba_delta_times_sim_ratio`). Data-matched: R_A 13.9–31.1, G 8.0–18.8 | `inversion.ba_calibrated_{RA,G}_range` | Not supported. The earlier G range contains 23.5 but not 37.4. The data-matched point ranges contain neither empirical point estimate, but both empirical 95% intervals overlap them (R_A 19.1–71.0 against 13.9–31.1), and both points lie inside the per-panel sampling range |
| 23.5 (G), 37.4 (R_A) | 23.46, 37.44 | `s02_inference.{G,RA}_rule4_ba` | Matches |
| National week share ≈ 50% | R_A 0.41 (0.12, 0.59); G 0.51 (0.21, 0.67) | `multiscale.{RA,G}_rule4_share_week` | Write "about 40–50%", with intervals |
| ICC 0.29 | 0.293 (G, calendar windows); 0.31 (G) and 0.25 (R_A) in rule windows | `s03_decomposition.G_cal3_icc_week`, `multiscale.*_icc_week` | Matches for G_cal3 |
| Reporting contribution within 17% | 17% is the equivalence margin; the measured share is 2.0% with a one-sided upper bound of 5.3% (R_A) | `reporting_regime.RA_rule4_margin_share_of_a_vR0.04`, `reporting_regime.RA_rule4_delta_a_share` | Reword: the contribution is below the 17% margin (equivalence shown) |
| 2024 voluntary period amplifies growth variability by 58% | IQR ratio 1.58 (0.95, 2.23); placebo 1.01 | `reporting_regime.iqr_ratio_2024_over_2023`, `reporting_regime.iqr_ratio_2025_over_2023_placebo` | Point estimate matches. The interval covers 1, so state it as an estimate with its interval |
| Participation above 95% in mandatory periods | Median 94.4%; 46.7% of state-weeks at or above 95% | `reporting_regime.regime1_state_week_median`, `reporting_regime.regime1_share_ge95` | Not supported; write "about 94%" |
| 2.1万 state-weeks | 21,641 raw rows (all jurisdictions); 10,506 state-weeks in the analysis panel | `s01_data.raw_rows`, `s01_data.panel_rows_Dge5` | Write "about 10,500 state-weeks" (约 1.05 万州周) |

### 4.7 措辞可支撑度 (how far each claim is supported)

| Claim | Rating | Strongest supportable wording (EN) | Strongest supportable wording (ZH) |
|---|---|---|---|
| The two statistics converge to the same latent variance | 有条件支持 | After simulation-based inversion, G and R_A imply similar latent variances (0.043 and 0.040; ratio 0.93, 95% CI 0.78–1.18), and they agree in 9 of 10 simulation settings | 经仿真反演，G 与 R_A 隐含的潜伏环境方差相近（0.043 与 0.040，比值 0.93，95% CI 0.78–1.18），10 种仿真设定中 9 种一致 |
| The implied v_R lies inside the theory band 0.04–0.08 | 有条件支持 | The implied v_R lies at the lower edge of the 0.04–0.08 band; the intervals extend below it, and the value depends on the assumed generation interval (0.013–0.063 across settings) | 反演所得 v_R 位于 0.04–0.08 区间下沿，区间估计向下越出该范围，且取值依赖于假定的代际间隔（各设定 0.013–0.063） |
| "Perfect alignment" (完美对齐) | 不支持 | The simulation reproduces the observed ICC and national-week share | 仿真复现了实测的 ICC 与全国周份额 |
| Theory b/a range "tightly contains" the empirical values (严密包容) | 不支持 | Empirical b/a lies above the calibrated theory range but inside the simulated sampling range (91st–93rd percentile) | 实证 b/a 高于校准后的理论范围，但位于模拟抽样分布之内（第 91–93 百分位） |
| First quantification (首次) | 不支持 | This study quantifies the distortion that delay and smoothing introduce into apparent statistics | 本文量化了延迟与平滑对表观统计量的畸变 |
| About half the intercept is a national-week effect | 支持 | About 40–50% of the intercept is shared nationally by week (R_A 0.41, G 0.51), and state pooling cannot remove it | 截距中约 40%–50% 为全国同步周分量（R_A 0.41，G 0.51），跨州合并无法消除 |
| Mandatory reporting keeps reporting noise within 17% | 有条件支持 | During mandatory reporting, the reporting contribution to the intercept was equivalent to zero within a pre-set margin of 17% (measured 2%, upper bound 5%) | 强制填报期内，报告波动对截距的贡献在预设 17% 等效界内与零等效（实测 2%，上界 5%） |
| The voluntary period amplified growth variability by 58% | 有条件支持 | In May–October 2024, under voluntary reporting, the dispersion of log growth was 1.58 times that of 2023 (95% CI 0.95–2.23; placebo 1.01). This is an observational association | 2024 年 5–10 月自愿填报期，对数增长的离散度为 2023 年同期的 1.58 倍（95% CI 0.95–2.23；安慰剂 1.01），属观察性关联 |
| Proves the institutional value of mandatory reporting (证明) | 不支持 | Consistent with higher facility participation being associated with lower variability; the design cannot exclude confounding | 结果与较高机构参与率伴随较低波动一致，本设计无法排除混杂 |
| Participation above 95% | 不支持 | Median facility participation was about 94% in mandatory periods and 35% in the voluntary period | 强制期机构参与率中位数约 94%，自愿期约 35% |
| A universal critical scale (普适) | 不支持 | Beyond m×, each doubling of cases buys at most a 25% variance cut, falling to 10% at 4 m×. This holds under the a + b/m model for these data | 超过 m× 后，病例数每翻倍至多降低 25% 的方差，至 4 m× 时降至 10%；该结论限于本数据与 a + b/m 模型 |
| The model departs from Parag et al.'s Poisson variance | 支持 | Poisson renewal variance has no floor; the fitted variance exceeds it twofold at about 22–28 cases per week and tenfold at about 260 | 泊松更新方差无下限；拟合方差在每周约 22–28 例时为其 2 倍，约 260 例时为其 10 倍 |

**Scope and extensions.**
- **What the results rest on.**
  - One country, one surveillance system and four influenza seasons.
  - Rule-defined early-growth windows.
  - A simulated generation interval that the inversion depends on.
- **What the data do not settle.** They cannot separate a plateau from slow power-law decay (section 1); m× is the scale at which the fitted floor equals the count term.
- **What the results do not show.** The reporting-regime contrast does not identify a causal effect.
- **Extensions that would sharpen the claims:**
  - external estimates of the generation interval;
  - a second country or a non-influenza pathogen with mandatory reporting;
  - multi-season shrinkage estimators that reduce the floor term.

## 5. 终修补充 (final polish, s10_polish.py, part `polish`)

Every value below is in `results.json` under `polish.<key>`. The simulated values come from the s09 data-matched panel. Each new setting uses 200 panels per v_R; the r1 cells for the s09 settings use 100 panels at v_R = 0.04.

### 5.1 Independence of weekly shocks: lag-1 autocorrelation (review P3/D5)

| Quantity | R_A | G | Key |
|---|---|---|---|
| Data r1 (state-cluster bootstrap 95% CI) | −0.095 (−0.168, −0.019) | 0.123 (0.039, 0.201) | `r1_emp_RA`, `r1_emp_G` |
| Simulated mean r1, φ = 0 (M0) [data percentile] | −0.074 [31st] | 0.114 [62nd] | `r1_sim_M0_*` |
| Simulated mean r1, φ = 0.25 [data percentile] | 0.007 [3rd] | 0.193 [10th] | `r1_sim_ar025_*` |
| Simulated mean r1, φ = 0.5 [data percentile] | 0.072 [0th] | 0.259 [0th] | `r1_sim_ar05_*` |
| Matched φ (interval from data CI) | 0.00 (0.00, 0.17) | 0.03 (0.00, 0.28) | `phi_matched_*` |

- r1 is computed on x/(block mean) − 1, with blocks defined by state and season. The apparent statistics already have a non-zero r1 under independent shocks, because of their construction. The test is therefore a comparison with simulation, not with zero.
- The data fit independent weekly shocks well. φ = 0.5 lies outside the simulated range for both statistics. φ = 0.25 is at the edge of the range for R_A.
- r1 also discriminates other settings:
  - shorter true GI of 2.2 d: data at the 0th percentile for both statistics;
  - synchronous count shock: data at the 99th percentile for G;
  - admission VMR 1.5: data at the 97.5th (R_A) and 100th (G) percentiles.

### 5.2 Generation interval bound from Cowling et al. 2009 (serial interval 3.6 d, 95% CI 2.9–4.3)

| True GI mean | v_R^impl R_A | v_R^impl G | Key |
|---|---|---|---|
| 2.9 d (new setting `gi_cow_lo`) | 0.026 (0.017, 0.036) | 0.027 (0.016, 0.040) | `vR_impl_*_GI_2.9_to_4.3` |
| 3.6 d (primary) | 0.040 (0.026, 0.056) | 0.043 (0.025, 0.064) | inversion part |
| 4.3 d (s09 `gi_long`) | 0.058 (0.038, 0.081) | 0.063 (0.037, 0.095) | `vR_impl_*_GI_2.9_to_4.3` |

- Across the external 95% interval for the interval mean, the implied variance spans 0.026 to 0.063.
- The analyst-side interval stays fixed at 3.6 d.
- In the 2.9 d setting, the ratio of the two statistics is 0.95 (0.79, 1.19).
- In the 2.9 d setting, the empirical b/a sits at the 52.5th (R_A) and 58th (G) percentiles of the simulated b/a at v_R = 0.02, the grid value nearest its implied variance (pooled b/a 36.1 and 21.2; keys `gi_cow_lo_*_ba_pooled_at_nearest_v`, `emp_ba_percentile_at_nearest` in parts/polish.json). A shorter generation interval therefore also accounts for the b/a excess. The r1 values (2.5th and 9th percentiles) count against it, and the two explanations imply v_R of about 0.026 and 0.040.

### 5.3 Admission overdispersion and the b/a excess (review P4/D8)

| Setting | v_R^impl R_A | v_R^impl G | pooled b/a at v_R = 0.04, R_A / G | Data percentile, R_A / G |
|---|---|---|---|---|
| M0 (Poisson admissions) | 0.040 | 0.043 | 28.0 / 16.2 | 93rd / 91st (s09 replicate distribution) |
| VMR 1.5 (`od05`) | 0.040 (0.026, 0.055) | 0.044 (0.026, 0.066) | 68.0 / 40.2 | 0th / 3.5th |
| VMR 2.0 (`od10`) | 0.039 (0.026, 0.055) | 0.045 (0.027, 0.066) | 81.0 / 48.4 | 0th / 0.5th |
| AR φ = 0.25 (`ar025`) | 0.033 (0.021, 0.047) | 0.033 (0.019, 0.049) | 23.4 / 12.6 | 98.5th / 96.5th |

- Extra-Poisson variation in weekly admissions raises the slope b. It leaves the implied variance almost unchanged.
- With VMR 1.5, b/a already exceeds both empirical values, 37.4 and 23.5.
- An admission overdispersion between 1 and 1.5 can therefore account for the excess. The simulation does not identify its size. The r1 values in 5.1 argue against VMR ≥ 1.5.

### 5.4 Latent balance scale (review P4/D7)

- Proposition 1 at the simulator values gives the latent balance scale. The inputs are R̄ = 1.176, an admission probability of 0.01, k = 0.35 and A = 1.216, applied to the v_R^impl draws.
- m× = 30.3 (95% CI: 21.7, 46.3) admissions per generation for R_A, about 59 (42, 90) per week at 7/3.6 generations per week. Keys: `latent_mcross_RA`, `_per_week`.
- For G the values are 28.3 (19.0, 48.3) per generation and 55 (37, 94) per week.
- The figure of 37.4 is the apparent ratio b/a on the scale of R_A. It is not the latent m×.

### 5.5 Agreement between statistics as equivalence (review D4)

- The 90% interval of the ratio of v_R^impl (R_A over G) is 0.81 to 1.13, which lies inside [0.8, 1.25]. Key: `ratio_vR_RA_over_G_ci90`.
- A TOST at α = 0.05 with a ±20% (log-symmetric) margin therefore declares the two estimates equivalent.

### 5.6 What the completeness bound allows (review D2)

- **The bound.** The one-sided 95% upper bound on the completeness share of a is 0.053. At the theory baseline (R̄ = 1.15, v_R = 0.04), this corresponds to τ = 0.040. Key: `tau_allowed_by_upper_share`.
- **Levers at this τ.** Under the Proposition 6 baseline (T = 10, m = 50, φ = ψ = 0.5), the total MSE is 15.22×10⁻³. Standardising reporting gives 14.65×10⁻³. A tenfold increase in cases gives 11.40×10⁻³.
- **Interpretation.** At the reporting volatility the data permit, more cases per period reduce the MSE more than standardising does. The τ = 0.15 example, where standardising gives 14.65 against 18.71, is therefore a hypothetical scenario and not an empirical estimate.

### 5.7 Checks on numbers from the user plan (final polish)

| Plan number | Status | Value used |
|---|---|---|
| 2.68 / 0.484 | Measured: median κ over the 17 scenarios of the first simulation | Data-matched κ is 2.729 / 0.516 |
| 0.0451, R̄ = 1.206 | R̄ = 1.2058 is measured. The value 0.0451 is not used, because the data-matched inversion gives 0.040 / 0.043 | 0.040 / 0.043 |
| 0.04–0.08 | Scenario band. Both point values lie at its lower edge | — |
| 13.4–28.7 | Not a calibrated b/a range. 13.4 coincides with `parag.latent_m_ratio1.5` | 13.9–31.1 (R_A), 8.0–18.8 (G) |
| 23.5, 37.4 (19.1, 71.0) | Measured | Unchanged |
| ~50% | National-week share | 41% / 51% |
| 17% | TOST margin at τ0 = 0.10 | Unchanged |
| 58% | Voluntary-period IQR ratio 1.58 (0.95, 2.23), observational | Unchanged |
| >95% | Not supported | Median about 93–94% |
| 21k state-weeks | 21,641 raw panel rows; the analysis windows use 1748 / 1756 | 1748 / 1756 |
