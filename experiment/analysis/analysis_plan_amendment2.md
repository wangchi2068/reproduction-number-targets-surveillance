# Analysis plan, Amendment 2 (second revision round), 2026-11-05

Paper: "What more cases buy: inferential targets for reproduction-number estimates and what cross-sectional dispersion diagnostics measure".

This amendment belongs to `analysis_plan.md`, which holds the locked plan and Amendment 1. Its sha256 and UTC time are recorded in `analysis_plan_amendment2.lock`, and the lock was written before any analysis listed below was run. Nothing in the locked plan or in Amendment 1 is edited; Section A corrects a statement of Amendment 1 by addition.

All analyses below are added in response to the second review round (items N2, N3, N6, N8, N9, M4, M5, minor 4, minor 12 and questions 2–6). The confirmatory analysis (E1, Section 2 of the plan) and its decision rule are unchanged. No analysis below is confirmatory: each new test is a pre-specified secondary analysis, and none of them can change the E1 verdict, α, the calibration, or the data. Every result below is reported whatever its direction. If a result goes against the current main line, the main line is rewritten rather than the specification changed.

## A. Correction of the timing statement in Amendment 1 (review N4) and the 06:03:14 edit (question 5)

### A.1 What Amendment 1 said, and what the logs show

Amendment 1 (written 06:02 UTC) says that only `rev_rho_grid.py` and `rev_generators.py` had finished by that time, and that only smoke tests with a few dozen draws had been run of `rev_intercept.py` and `rev_season.py`. That statement is wrong. File modification times and the logs in `experiment/logs/` give the following sequence (all UTC, 2026-11-05):

| Time | Event |
|---|---|
| 05:42:50 | Plan locked |
| 05:48:36 / 05:48:46–05:48:52 | `rev_rho_grid.py` code saved / full run |
| 05:50:54 / 05:51:16–05:52:25 | `rev_generators.py` code saved / full run |
| 05:52:23 / 05:52:54–05:52:57 | `rev_daily_renewal.py` code saved / full 200-panel run (not mentioned in Amendment 1) |
| 05:57:59 / 05:58:10 | `rev_intercept.py` code saved / smoke run |
| 05:58:22–05:59:18 | `rev_intercept.py` full run (wild bootstrap-t p = 0.0001 logged at 05:58:56; `rev_intercept.json` written 05:59:18) |
| 06:00:43 | `rev_season.py` smoke run (B = 40) |
| 06:00:54–06:01:53 | `rev_season.py` full run, B = 9,999 |
| 06:01:39 / 06:01:40 | `rev_identification.py` code saved / run finished |
| 06:02:26 | Amendment 1 written |
| 06:03:14 | `rev_season.py` edited |
| 06:03:15–06:04:15 | `rev_season.py` full rerun; its outputs are the ones in `results.json` |

Correct statement: by the time Amendment 1 was written, all six revision scripts (`rev_rho_grid`, `rev_generators`, `rev_daily_renewal`, `rev_intercept`, `rev_season`, `rev_identification`) had completed full runs. Amendment 1 recorded implementation details that were already in the code, not details fixed before the full runs. For E1 the code with the restricted fit and the 1e-7 clamp was saved at 05:57:59, before the first full E1 run at 05:58:22. No evidence was found that any of these details were changed after a result was seen. The primary E1 inference (two-way sandwich, p = 5.1e-9) does not depend on any of them.

### A.2 The 06:03:14 edit to `rev_season.py`

The text of the edit cannot be recovered: the working copy is not under version control, and the only other copy (the clean-copy run of 06:13 UTC) is byte-identical to the current 17,517-byte file. What can be established is the following comparison of the two full runs, read from `logs/rev_season.log`, which logs both runs:

| Logged quantity | Run 06:00:54 (before the edit) | Run 06:03:15 (after the edit) |
|---|---|---|
| Primary calendar-week permutation p (MCSE) | 0.1022 (0.00303) | 0.1022 (0.00303) |
| Week-cluster wild bootstrap p | 0.0108 | 0.0108 |
| S-a4 / S-a6 / S-a8 permutation p | 0.031 / 0.2607 / 0.0808 | 0.031 / 0.2607 / 0.0808 |
| Circular-shift p | 0.40476 | 0.40476 |
| Minimum detectable ratio | 4.4753 | 4.4753 |

Conclusion: every number that both runs logged is identical, so the edit did not change the primary season test, any pre-specified calibration, the data sensitivities, or the power analysis. Blocks that the log does not print (the effect-size table and the descriptive coverage table) cannot be compared, so an edit affecting only those blocks, or only the output format, cannot be ruled out. The text of the edit is reported as not recoverable. To prevent a repeat, every script run of this round is preceded by recording its sha256 in `revision2/script_hashes.jsonl` (Section H).

## B. Functional form, extrapolation and population size (review N2, questions 2 and 3)

Script `rev_r2_floor.py`. Sample R_A, x and L exactly as in the locked plan, Section 1. Seeds `[20261105, 70, j]`.

### B.1 Observed mean loss by count decile (N2 a), descriptive

- Deciles of D = C_t over the 1,748 rows (`np.quantile` cut points; ties go to the lower decile).
- Per decile: n, mean C, mean L, and the relative RMSE √(mean L), with 95% percentile intervals from the calendar-week block bootstrap (B = 1,999, seed `[20261105, 70, 1]`; x̄ recomputed per draw; decile cut points held at the observed values).
- Model predictions per decile, as the mean of the fitted E(L|C) over the decile's rows, from three fits to all rows: the affine a + b/C (primary model); the power law c·C^(−γ); the three-parameter model a + b·C^(−γ) (B.3). Bands come from the same bootstrap draws (refit per draw).
- Table `tab_r2_deciles.tex` and figure `fig_r2_deciles.{pdf,png}`.

### B.2 Population size (N2 b, question 2)

**Data.** State resident population from the U.S. Census Bureau Population Estimates Program, vintage 2025 (`NST-EST2025-ALLDATA.csv`, https://www2.census.gov/programs-surveys/popest/datasets/2020-2025/state/totals/NST-EST2025-ALLDATA.csv). The file is stored in `experiment/data/external/` with its sha256 recorded in `DATA_MANIFEST.md` and `data_contract.md`. The covariate is the log of each state's mean of POPESTIMATE2022–POPESTIMATE2025, the four years spanned by the sample seasons. z_pop is this log population standardised to mean 0 and SD 1 over the 51 jurisdictions (states weighted equally). The covariate is fixed per state.

**Primary population model (pre-specified secondary test).** E(L | C, state) = a0 + a1·z_pop + b/C, gamma quasi-likelihood (IRLS with weights 1/μ², fitted means required positive). Inference for a1 uses the two-way CGM sandwich of the stacked (x̄, a0, a1, b) equations (`rev_common.cgm_from_loss`), with t on 41 df. Two-sided α = 0.05.

Decision rule:
- "a falls significantly with population" if the 95% CI for a1 lies entirely below 0;
- "a rises significantly with population" if it lies entirely above 0;
- otherwise "no population gradient detected", reported with the CI and the implied ratio of a at the 90th and 10th population percentiles, (a0 + a1·z90)/(a0 + a1·z10).

Secondary inference for a1: the state-cluster bootstrap (states resampled, B = 1,999, seed `[20261105, 70, 2]`) and the calendar-week block bootstrap (B = 1,999, seed `[20261105, 70, 3]`), percentile intervals.

Secondary specifications (reported, never replacing the primary):
- the log-link version a·exp(θ·z_pop) + b/C (quasi-deviance minimisation; θ with the same two bootstraps);
- the primary model with b also varying linearly in z_pop.

**Population tertiles.** States are ranked by population and split 17/17/17. For each tertile there is a separate affine fit a_k + b_k/C, with intervals from the two-way sandwich on t with min(17, weeks) − 1 = 16 df. There is also a pooled model with tertile intercepts and a common slope (`fit_group_affine`), with the ratio of the highest to the lowest tertile intercept and its state-cluster bootstrap CI (B = 1,999, seed `[20261105, 70, 4]`).

**What the result means for the main line (fixed now).**
- If a falls significantly with population, `results_summary_r2.md` must state that:
  - the floor is not common to jurisdictions;
  - the cross-jurisdiction statement "more cases help little beyond a few dozen" must be restricted to comparisons within a jurisdiction of given size;
  - the error at a given count must be reported per population tertile.
- If no gradient is detected, the main line stands, and the CI for a1 is reported as the bound on how much of the floor can be attributed to population size.

### B.3 Count exponent left free (N2 c, question 3)

- Model E(L|C) = a + b·C^(−γ), gamma quasi-likelihood: minimisation of Σ(L/μ + log μ) over (a, b, γ), with γ ∈ [0.05, 3], μ > 0 required, and a allowed to take any sign subject to μ > 0. Starting values come from the affine fit and from a grid over γ.
- Intervals: calendar-week block bootstrap (B = 1,999, seed `[20261105, 70, 5]`) as primary, and state-cluster bootstrap (B = 999, seed `[20261105, 70, 6]`) as secondary; percentile intervals for a, γ and b.
- Also reported:
  - the share of draws with a ≤ 0 and with a < 0.001;
  - the quasi-deviance difference to the affine fit;
  - the share of draws in which the three-parameter fit improves the quasi-deviance by more than 2.
- Rule: "γ = 1 compatible with the data" if the primary 95% interval for γ contains 1. "The floor survives a free exponent" if the primary 95% interval for a excludes 0. Both outcomes are reported whichever way they fall.

### B.4 Model-free floor at the high-count end (N2 d)

For rows with D at or above the 75th, 90th and 95th percentiles of D (the 90th percentile defines the top-decile row set):
- n and mean C;
- mean L and relative RMSE √(mean L);
- the floor estimates mean(L) − mean(b̂/C) (affine b̂) and mean(L) − mean(b_P/C) (Poisson reference b_P, which uses no fitted parameter), with relative SD √(floor).

The 95% percentile intervals come from the calendar-week block bootstrap (B = 1,999, seed `[20261105, 70, 7]`), with the cut points held fixed and b̂ and b_P refitted in each draw. This is the number the abstract may quote for the high-count end in place of the extrapolation to C = 1,000.

## C. Bias of â relative to the transmission-plus-reporting total (review N3, question 4)

Script `rev_r2_total_bias.py`.

- **Estimand.** The relative bias E(â)/E(v_R + v_H) − 1 under each generator with a reporting shock, where v_R and v_H are the per-panel transmission and reporting variances over the retained rows, as already stored by `rev_generators.py`. Its 95% MC interval comes from a bootstrap over panels (B = 2,000, seed `[20261105, 71, 1]`). The MC 2.5–97.5% range of the per-panel ratio â/(v_R + v_H) is also reported.
- **Generators.**
  - From the existing panels: FD-S8, FD-S9-25, FD-S9-50, FD-S9-75 (the S9 laws are post hoc).
  - FD-S1, which has no reporting shock, as a reference.
  - New share sweep: the FD-S8 law with transmission share 1.0, 0.75 and 0.25 (0.5 is FD-S8 itself), calibrated exactly as FD-S8 (calibration stream 910 + j), 200 panels each, seed `[20261105, 41, j]`.
- **Effective weight.** w_H = (E(â) − r_1·E(v_R))/E(v_H), where r_1 is the FD-S1 recovery ratio. It is the weight with which the reporting variance enters the intercept, with an interval from the same panel bootstrap.
- **Wording rule (fixed now).**
  - If the relative bias interval excludes 0 under every generator with a reporting shock, the manuscript must say that the intercept "reflects" transmission and reporting variability, with the reporting part entering with weight w_H. It may not call â an estimate of the total.
  - A correction is reported only as a mechanism-conditional rescaling â/t_g, where t_g = E(â)/E(v_R + v_H), over the generator range. It is not used in any headline number.

## D. Generator realism (M4, M5) and external reporting data (question 6)

### D.1 A state-realistic generator family (script `rev_r2_generator.py`)

**Aim.** To test whether one generator reproduces at once:
- the between-state baseline dispersion;
- the upper tail of the counts;
- the counting slope b/b_P at an assumed realistic ρ.

**Law "SR"** (state-realistic), a fresh field each replicate as in FD-*:

1. **Trend-anchored transmission.**
   - Growth for row i: g_i = g0·(g_trend,i/g0)^τ·exp(s_i·z_i − s_i²/2), with z_i iid standard normal, redrawn each replicate.
   - g_trend,i = S_fol/S_cur, where S is the three-point geometric moving average of the observed counts within each state run (two-point at the run ends; counts floored at 0.5).
   - τ = 1 unless stated.
   - The mean path follows `Design.mean_path` (run level anchored to the observed geometric mean, times λ), so the expected trajectory follows the smoothed observed epidemic curve.
2. **State baseline multiplier m_state** (mean 1 over rows), which scales both the transmission variance (s_i² = σ²·m_state) and the reporting-shock variance.
   - "pop" version: m_state ∝ max(â0 + â1·z_pop, 0.1·â0), using the point estimates of the primary population model of B.2.
   - "fe" version: m_state ∝ max(â_state, 0.1·mean(â_state)), using the state intercepts of `fit_group_affine` on the observed losses.
3. **Offspring and reporting.**
   - Latent branching with NB offspring k = 0.3, started at mean/ρ as in `Design.draw`.
   - Reporting probability ρ·shock, with a lognormal state-week shock carrying half of the variance target (as FD-S8; scaled by √m_state).
   - Clustered admissions: c ~ Binomial(latent, min(1, ρ·shock)/(1 + θ_c)) and C = c + Poisson(θ_c·c), so that E(C) is unchanged and the index of dispersion of reporting rises with θ_c.
4. **Calibration.** Identical to the shipped and FD generators (`ma_audit.calibrate`: λ to median retained D = 45, amp to MC mean â = a_obs, λ refitted). σ is calibrated so that the relative variance of wl(g) over all rows equals the transmission share of amp·a_obs. If the trend alone exceeds the target, so that σ = 0 is not enough, the cell is reported as "calibration failed" and its panels are still summarised.

**Cells, 200 panels each.** Seed `[20261105, 42, i]`; calibration stream 920 + i.

| i | Cell | State multiplier | θ_c | ρ | τ |
|---|---|---|---|---|---|
| 0 | SR-pop-c0 | pop | 0 | 0.01 | 1 |
| 1 | SR-pop-c25 | pop | 0.25 | 0.01 | 1 |
| 2 | SR-pop-c50 | pop | 0.5 | 0.01 | 1 |
| 3 | SR-fe-c25 | fe | 0.25 | 0.01 | 1 |
| 4 | SR-pop-c25-r003 | pop | 0.25 | 0.003 | 1 |
| 5 | SR-pop-c25-t50 | pop | 0.25 | 0.01 | 0.5 |

**Pre-specified tests** (each on the MC 2.5–97.5% range of the per-panel statistic):

| Test | Statistic | Observed value | Passes if |
|---|---|---|---|
| T1 | b/b_P | 1.257 | the range contains it |
| T2 | 90th percentile of retained D | 388.6 | the MC mean lies in [350, 428] (±10%) |
| T3 | pooled/within slope ratio b/b_w | 1.476 | the range contains it |
| T4 | between-state baseline dispersion: coefficient of variation of the state intercepts from `fit_group_affine` | observed value computed in the same script | the range contains it |

- **Success rule.** The generator "reproduces the data's key features" if T1–T4 all pass in at least one cell; otherwise M4/M5 are reported as unresolved, naming the failing tests per cell.
- **Also reported per cell:**
  - the SD ratio (observed SE / MC SD of â);
  - two-way coverage of E(â) and of the true total;
  - a/v_R and a/(v_R + v_H);
  - the 10th percentile and median of D.
- **Labelling.** The θ_c grid and the use of the population and state-FE multipliers are chosen with knowledge of the observed counting slope and the observed between-state dispersion. They are labelled as data-informed generator design. No inferential result of the manuscript uses these generators.

### D.2 External reporting data (question 6, script `rev_r2_reporting.py`)

1. **Hospital reporting coverage.** The NHSN panel's own field: the percentage of hospitals reporting, per state-week. Coverage-adjusted counts C·100/coverage are formed per row, x is recomputed from adjusted counts, and the affine model is refitted (same rows, two-way sandwich). The change a_adj − a, with a calendar-week block bootstrap interval of the difference (B = 1,999, seed `[20261105, 72, 1]`), bounds the part of the intercept that is due to week-to-week changes in hospital participation under the proportional-participation assumption.
2. **Revision history (backfill).** First-issue and latest-issue values of `confirmed_admissions_flu_ew` per state and epiweek come from the Delphi Epidata API (`nhsn` source; https://api.delphi.cmu.edu/epidata/covidcast/), fetched for the sample weeks that it covers. The raw responses are stored in `experiment/data/external/` with sha256. On the sample rows covered by both issues:
   - the share of rows whose latest-issue count equals the vintage count used in the paper;
   - the revision variance v_rev = mean(((x_first − x_latest)/x̄_latest)²);
   - the fitted intercepts on first-issue and latest-issue counts and their difference, with calendar-week block bootstrap intervals (B = 1,999, seed `[20261105, 72, 2]`).

   This is descriptive. It bounds the backfill component (provisional vs. settled counts), not the residual reporting error that remains in the settled vintage, and it is reported with that limitation.
3. **Search record.** The search for other public data on state reporting completeness is recorded with URLs in `revision2/reporting_sources.md`, including sources that were found but give no state-week completeness measure.

## E. Wild bootstrap clamp (review N8)

Script `rev_r2_clamp.py`. The E1 week-cluster wild bootstrap-t, exactly as in Amendment 1 item 2 (same seed `[20261105, 21]`, B = 9,999, so the Webb weights are the same draws), under these variants:
- clamp thresholds 1e-9, 1e-7 (shipped), 1e-5 and 1e-3;
- no clamp: negative bootstrap losses are allowed, and the linear gamma-QL fit requires only positive fitted means;
- "exclude clamped draws": the shipped 1e-7 run restricted to bootstrap samples in which no loss was clamped.

Reported for each variant:
- the share of row-level losses clamped (the shipped quantity, 0.135);
- the share of bootstrap samples with at least one clamped loss;
- the p-value with MCSE and the 2.5/50/97.5% quantiles of t*.

If fewer than 100 draws remain in the exclusion variant, it is reported as not evaluable. The wild bootstrap is a secondary calibration of E1, so no outcome here changes the E1 verdict.

## F. Bootstrap conservativeness (N9, minor 12) and the Poisson normalisation (minor 4)

### F.1 Pigeonhole bootstrap against fresh-draw MC SD (script `rev_r2_bootstrap_sd.py`)

- **Panels.** The first 100 panels of each fresh-draw generator FD-S1, FD-S8 and FD-S9-25/50/75 are regenerated exactly: same calibration, as stored in `rev_generators_calibration.jsonl`, and the same panel stream `[20261105, 40, i]`. The regenerated â must match the stored â; any mismatch is reported.
- **Bootstrap.** Each panel gets B = 199 state × calendar-week pigeonhole bootstrap draws (`revision_stats.bootstrap_weights`, L = 1, seed `[20261105, 73, i, panel]`).
- **Reported per generator:**
  - the mean bootstrap SD of â divided by the MC SD of â over the 200 stored panels;
  - the same ratio for the two-way sandwich SE;
  - the coverage of the percentile interval for E(â).
- **Rule.** "The pigeonhole bootstrap is conservative" only if the bootstrap/MC SD ratio exceeds 1 under every fresh-draw generator. It is compared with the shipped fixed-field ratios of 1.6–2.5 and with √3 ≈ 1.73, the upper limit of the pigeonhole variance inflation for a two-way random-effects mean. If the rule fails, the word "conservative" is withdrawn.

### F.2 Exact Poisson reference (script `rev_r2_poisson_ref.py`)

The shipped b_P is a plug-in: the delta approximation with row R_A/x̄ replaced by 1. Two corrections are reported:
- **Normalised plug-in.** b_P,norm = mean[κ(g)²(1 + 1/g)(x_i/x̄)²] over the sample rows, which replaces the unit normalisation by the observed row ratio.
- **Simulated exact reference.**
  - B = 400 panels of independent Poisson counts with means equal to the observed C_t and C_{t+1} (seed `[20261105, 74]`).
  - Rows are kept if the simulated C_t ≥ 5, and x̄ is recomputed per panel.
  - The counting loss is L_cnt = ((x_sim − x_μ)/x̄_sim)², where x_μ is the transform of the mean ratio. b_P,sim is the slope of the affine gamma-QL fit of L_cnt on simulated D, pooled over panels, with the MC 2.5–97.5% range of the per-panel slope.
  - The fitted intercept of L_cnt is also reported; it should be near 0.

Reported: b/b_P under each of the three references and the relative change against the shipped 1.257. The interpretation of b/b_P (CI includes 1) is re-stated under each.

## G. Interval types (review N6)

Script `rev_r2_interval_types.py` writes a table listing every interval quoted in the manuscript's abstract, results, discussion and tables. For each interval it gives:
- the quantity;
- the interval;
- its type: two-way CGM sandwich (t, df), calendar-week block bootstrap percentile, state-cluster bootstrap, pigeonhole bootstrap, MC 2.5–97.5% range, or propagated envelope;
- the source key in `results.json`.

It is descriptive and has no decision rule.

## H. Reproducibility for this round (N5, M7)

- All new numbers go to `analysis/results.json` under the key `revision_2026_11_r2`, each block naming its source script and seeds.
- The `analysis_plan` field of `results.json` is changed to the published relative paths: `analysis/analysis_plan.md`, `analysis/analysis_plan.lock`, `analysis/analysis_plan_amendment2.md` and `analysis/analysis_plan_amendment2.lock`.
- The plan, its lock, Amendment 1 (inside the plan file) and this amendment with its lock are published in `output/rn-targets-peer-review/experiment/analysis/`.
- `run_all.sh` gains the new scripts. One clean-copy run of the whole pipeline is made, and its key numbers are compared with the working copy.
- Before each full run of a new script, its sha256 and the UTC time are appended to `revision2/script_hashes.jsonl`.
- Any change to a script after its first full run is logged there with the reason. If such a change alters a reported number, both values are reported.
- Tables `tables/tab_r2_*.tex` and figures `figures/fig_r2_*.{pdf,png}` are produced by `rev_r2_render.py` from `results.json`.
