# Pre-specified analysis plan for the revision (locked before any new analysis was run)

Paper: "What more cases buy: inferential targets for reproduction-number estimates and what cross-sectional dispersion diagnostics measure".
Working copy: `.work/rn-targets-peer-review/experiment/`. Lock: `revision/analysis_plan.lock` (sha256 of this file + UTC time).
Changes after locking are appended at the end as dated Amendments with reasons; the text above the Amendments section is never edited.

## 0. Disclosure of what was known when this plan was written

This plan is written after the peer review, not before it. The following results were known to the analyst at the time of writing and are listed so that readers can judge every choice below against them:

- Pooled intercept a = 0.01495 (two-way CGM 95% CI 0.0109–0.0191); counting slope b = 0.560; b/a = 37.4 (24–58); b/b_P = 1.257.
- Season test (W_CGM = 17.80, 3 df) p-values: asymptotic χ² 0.0005; bootstrap-covariance χ² 0.064–0.220; generator calibrations 0.005–0.045; state-cluster wild bootstrap 0.0125; within-state permutation 0.0455; **week-cluster wild bootstrap 0.014; calendar-week permutation 0.101** (review recomputation, B = 4,999).
- Week centring lowers a by 34.7% (paper definition) to 42.6% (review definition).
- Generator MC SD of a is smaller than the observed sandwich SE; S9 and the decision rule were fixed after seeing the data.

Because both calendar-respecting season calibrations were known, the choice of the primary season calibration below is justified only by design arguments (Section 3.2), and the consequence of the opposite choice is reported in the results summary.

## 1. Data and sample (frozen)

- NHSN Hospital Respiratory Data, vintage 2026-10-20 (`data/external/nhsn_hrd_vintage_2026-10-20.csv.gz`, sha256 9685ddd4…), influenza admissions.
- Primary sample R_A: `revision_data.load_sample("RA", 3.6)` with the authors' rule windows (national first-growth phase, threshold 1.05): 1,748 state-weeks, 51 jurisdictions, 42 calendar weeks, 4 seasons. No other window, subset or covariate is used for any primary analysis.
- x = Wallinga–Lipsitch transform of (C_{t+1}+0.5)/(C_t+0.5), gamma generation interval mean 3.6 d, SD 1.6 d; analyst scale D = C_t; floor C_t ≥ 5.

## 2. Confirmatory primary analysis: the count-independent intercept

- Estimand E1: intercept a of E(L | C) = a + b/C, L = ((x − x̄)/x̄)², gamma quasi-likelihood / FGLS fit (`revision_stats.fit_loss`).
- H0: a = 0 against H1: a ≠ 0 (two-sided). α = 0.05.
- Primary inference: Wald t-statistic with the two-way (state × calendar week) CGM sandwich SE of the stacked (x̄, a, b) estimator (`ma_calibration.sandwich`), reference t with 41 df (G_week − 1).
- Robustness rule (pre-specified): the claim "between-period relative dispersion does not vanish with more cases" is called **robust** if (i) the primary test rejects and (ii) the 95% lower limit of a stays above 0 in every sensitivity analysis S-a1 to S-a9 below. If (i) holds but (ii) fails anywhere, the claim is reported as significant but sensitivity-dependent, naming the failing analyses. If (i) fails, the claim is demoted.
- Secondary inference for a (reported regardless): state-only and week-only sandwich; week-cluster wild bootstrap-t with a = 0 imposed (Webb 6-point weights, B = 9,999); calendar-week cluster (block) bootstrap percentile CI (B = 1,999).
- Sensitivity analyses for a, b, b/a (all reported; none replaces the primary):
  - S-a1 generation-interval mean 3.2 and 4.0 d (SD 1.6);
  - S-a2 generation-interval SD 1.2 and 2.0 d (mean 3.6);
  - S-a3 count floors C_t ≥ 10, 20, 50;
  - S-a4 the authors' fixed calendar windows (`revision_data.CALENDAR`, 3 seasons);
  - S-a5 rule-window thresholds 1.02 and 1.10;
  - S-a6 excluding transitions with either week ending between 2025-10-03 and 2025-11-13 (CDC data-channel pause);
  - S-a7 leave-one-season-out (4 fits);
  - S-a8 reporting completeness: rows with hospital reporting coverage ≥ 90% in both weeks;
  - S-a9 leave-one-calendar-week-out influence (42 fits; report min and max of a and of its lower limit).
- Decomposition estimand E3 (secondary, estimand change, not a sensitivity): within-week intercept a_W from week-centred losses. Primary definition: leave-state-out week centre, pooled-mean scale (removes self-inclusion). Also reported: paper definition (week mean incl. self, pooled-mean scale), review definition (week mean centre and scale), season centring. Inference: calendar-week block bootstrap (B = 1,999) percentile CI; H0: a_W = 0 judged by whether the 95% CI excludes 0. Share of a removed by week centring reported with CI.

## 3. Secondary pre-specified hypothesis test: equal season intercepts

### 3.1 Statistic, hypothesis, decision

- Model: season-specific intercepts a_s, common slope b. H0: log a_1 = … = log a_4 (omnibus, 3 df; a Wald statistic is two-sided by construction). α = 0.05. No one-sided version is used.
- Statistic: W_CGM, the Wald statistic with the two-way CGM sandwich covariance (`ma_calibration.season_wald`), as shipped (observed 17.797).
- **Primary reference distribution: calendar-week permutation.** Season labels are permuted across the 42 calendar weeks (each calendar week keeps all its state rows; the number of weeks per season, 9/12/8/13, is preserved); W_CGM is recomputed for each permutation; p = (1 + #{W* ≥ W_obs}) / (B + 1), B = 9,999, seed `[20261105, 31]`. MCSE reported.
- Decision rule: the season difference is called statistically significant only if the primary p < 0.05. Otherwise every season statement is reported as exploratory/descriptive; the abstract and main conclusions may not rest on it. The same rule applies whatever the secondary calibrations show.
- Multiplicity: E1 (Section 2) and the season test are separate claims; Holm-adjusted p-values over the two are reported in addition.

### 3.2 Why calendar-week permutation is primary (design argument)

The season regressor is constant within a calendar week and is assigned at the calendar-week level; with 42 weeks and 4 seasons the effective number of "treated" clusters is small. For few treated clusters, MacKinnon & Webb show that wild cluster bootstrap tests can be badly mis-sized and recommend randomization inference at the level of assignment. Permutation at the level where the regressor varies (calendar week) is therefore the reference distribution that matches the design; it imposes exchangeability of weeks under H0 but no parametric or generator assumptions. Known limitation: it ignores serial dependence between adjacent weeks.

### 3.3 Secondary calibrations (all reported in one table, with MCSE where Monte Carlo)

Asymptotic χ²(3); bootstrap-covariance χ² (as shipped); generator calibrations S1/S6/S8/S9 (as shipped, fixed transmission field) and the fresh-draw generators of Section 5; state-cluster wild bootstrap (as shipped, B = 2,000); within-state permutation (as shipped); week-cluster wild bootstrap (Webb, null imposed, B = 9,999, seed `[20261105, 32]`); contiguous-block calendar permutation (season blocks of consecutive weeks kept intact; all 4! = 24 orderings of the four seasonal blocks cannot change the statistic, so instead circular shifts of the season-label sequence over the 42 ordered weeks are used: 42 shifts, exact p).
Data sensitivities with the primary calibration: S-a6 (pause exclusion), S-a8 (coverage ≥ 90%), S-a4 (fixed calendar windows, 2 df).

### 3.4 Effect sizes, intervals and power (reported regardless of significance)

- Season intercepts a_s with 95% intervals from the season-stratified calendar-week block bootstrap (B = 1,999, percentile on the log scale) and from the CGM sandwich (delta method on log a_s).
- Ratios a_s / a_pooled; highest/lowest ratio with bootstrap CI (noting that a max/min ratio is biased upward). The contrast 2023-24 vs the other three seasons is labelled post hoc (chosen after seeing the data).
- Design-based power and minimum detectable ratio: null-imposed losses L0 = L · (a + b/C)/(a_s + b̂/C); injected pattern L_δ = L0 · (a e^{δ c_s} + b/C)/(a + b/C), c_s = the observed centred log-intercept pattern; season labels assigned by random calendar-week permutation; loss-level W (x̄ fixed) and its permutation critical value under δ = 0 (2,000 permutations); 400 draws per δ ∈ {0, 0.25, …, 2.0}; power curve and MDR (highest/lowest ratio at 80% power, linear interpolation). Seed `[20261105, 33]`.
- Reporting participation: mean hospital-reporting coverage per season over the sample rows, and the dated context (voluntary reporting May–Oct 2024; data-channel pause Oct–Nov 2025). Season intercepts are interpreted only as transmission-plus-reporting totals (Proposition 2).

## 4. Counting slope and crossover (estimation, no test)

b, b_P (independent-Poisson reference, `s03b_thinning.poisson_reference`), b/b_P, m× = b/a with sandwich (log-scale delta) and calendar-week block bootstrap intervals; the same under S-a1 to S-a8. "What more cases buy": predicted relative error components a and b/C at C = 10, 37, 100, 1,000 with intervals.

## 5. Identification and generator analyses (M2, M4, M5)

### 5.1 Implied transmission variability (M2)

For each generator g with recovery ratio r_g = E(â)/v_R (transmission) and t_g = E(â)/(v_R + v_H) (total): implied relative SD = sqrt(a_obs / r_g) and sqrt(a_obs / t_g), with intervals propagating the two-way 95% CI of a and the MC 2.5–97.5% range of the per-replicate ratio. Generators: S1, S3, S4, S5, S6, S7a, S7b, S8, S9 (shipped), the fresh-draw generators FD-* (5.2), the extended (ρ, k) grid (5.3) and the daily-renewal generator (5.4); variants for GI mean 3.2/4.0 d and for the within-week estimand. Report the envelope over all mechanisms (min lower – max upper) and over the shipped share/scale grid. Wording rule: the implied transmission SD is reported as "not identified; mechanism-dependent envelope", never as a single value; sqrt(a) is reported as the relative SD of the transmission-plus-reporting total. The (ρ, k) collection is called a grid.
Declared rule for grid cells: a cell "reproduces the counting slope" if the MC 2.5–97.5% range of b/b_P contains the observed 1.257.

### 5.2 Generator variance mismatch (M4)

Hypothesis fixed in advance: the shipped generators hold the transmission field fixed across replicates (one draw of z, seed `[SEED, 1]`), so their MC SD of â reflects only counting and reporting noise, while the sandwich SE also reflects sampling of the heterogeneity itself. Fresh-draw generators redraw the heterogeneity field each replicate with the calibrated σ:
FD-S1 (iid, S1 law), FD-S8 (S8 law), FD-S9-25 / FD-S9-50 / FD-S9-75 (S9 law with national calendar-week share 0.25 / 0.5 / 0.75). λ and amp are calibrated exactly as in `ma_audit` (median retained D = 45; MC mean of â = a_obs). 200 panels each, seed `[20261105, 40, i]`.
Reported for each: MC SD of â and its ratio to the observed two-way SE (0.00203); coverage of the two-way 95% CI for the pseudo-true value E(â) and for the true transmission+reporting total; coverage of the b/a CI; b/b_P; pooled/within slope ratio; week-centring reduction. Labels: S9 and the 30-cell decision rule are marked post hoc in all outputs.
Pre-specified check SC5: the variance mismatch is "explained by fixed-field conditioning" if the SD ratio for at least one fresh-draw generator lies in [0.67, 1.5]; otherwise it is reported as unexplained.
Pre-specified check SC4: the sandwich CI for a is called adequate if its coverage of E(â) is ≥ 0.90 under every fresh-draw generator.

### 5.3 Extended (ρ, k) grid (M5)

ρ ∈ {0.001, 0.003, 0.01, 0.03, 0.1, 0.25, 0.5, 1} × k ∈ {0.1, 0.3, 1, 3, ∞}, S1 transmission field as in the shipped grid, λ recalibrated per cell, 100 panels per cell, seed `[20261105, 50, i, j]`. ρ = 0.001–0.01 is labelled an assumed range for the probability that an infection becomes a reported laboratory-confirmed admission, based on CDC influenza burden estimates (hospitalisations per symptomatic illness about 1–2% and illness per infection below 1); it is not estimated.
Pre-specified checks: SC2 intercept recovery a/v_R within [0.9, 1.1] in every cell; SC3 at ρ ≤ 0.01 the MC mean of b/b_P lies within [0.9, 1.1] for every k (if so: offspring dispersion cannot explain the observed 1.257 at realistic ρ and k is not identified from the counting slope).

### 5.4 Two generations per week (M5)

Daily-renewal generator: daily infections I_d ~ Poisson or NB(k) with mean R_{w(d)} Σ_s w_s I_{d−s}, w = discretised gamma(3.6, 1.6 d) generation interval; weekly R_j drawn afresh each replicate from the S1 law (iid lognormal heterogeneity around the observed geometric-mean growth, σ as calibrated for S1); runs of the real design; weekly reported counts Binomial(Σ_week I, ρ); λ calibrated to median retained D = 45. Cells: ρ ∈ {0.003, 0.01, 0.1} × offspring {Poisson, NB k = 0.3}; 200 panels each, seed `[20261105, 60, i]`.
Reported: a/v_R (v_R = relative variance of the weekly R_j over retained rows), a/v_R̄ with R̄ the mean of R_j and R_{j+1}, b/b_P, lag-1 correlation of x.
Pre-specified check SC1: the one-generation recovery map carries over if a/v_R lies in [0.8, 1.25] in all cells; otherwise the implied-SD conversion is declared mapping-dependent and its envelope widened by the renewal recovery range.

## 6. Reproducibility (M7)

`run_all.sh` runs every script whose output the manuscript cites, including `ma_wild_bootstrap.py`, `ma_manifold_envelope.py` and all new `rev_*.py`; one clean-copy run with logs in `logs/`. All new numbers go to `analysis/results.json` (key `revision_2026_11`) with intervals, source script and MCSE; tables and figures are produced by scripts in the same run.

## 7. What will be reported regardless of outcome

Every row of Sections 2–5 is reported, significant or not. No calibration, window, subset, covariate, generator or α is chosen or dropped on the basis of a p-value. If the primary season test is not significant, it is demoted to exploratory and the main line rests on Sections 2, 4 and 5.

## Amendments

(none at lock time)

### Amendment 1 — 2026-11-05 06:02 UTC (implementation details; no decision rule changed)

Verification of the locked text: the first 14,195 bytes of this file have sha256 `7d3a2c1e8ca019d10f7d9bc6f6cb07bab63d1276357c24b1da573fba0e7533c9`, as recorded in `analysis_plan.lock`; everything after that is this and later amendments.

The locked plan left the following implementation details open. They were fixed in code before the full runs of the scripts concerned. Two scripts had already finished by this time: `rev_rho_grid.py` and `rev_generators.py`. For `rev_intercept.py` and `rev_season.py`, only a smoke test with a few dozen draws had been run. No detail below was changed after a result was seen. None of them changes a hypothesis, a statistic, a primary calibration, α, or a decision rule.

1. **Seeds not listed in the plan**
   - E1 wild bootstrap-t: `[20261105, 21]`.
   - E1 block bootstrap: `[20261105, 22]`.
   - E3 block bootstrap: `[20261105, 23]`.
   - Season data-sensitivity permutations: `[20261105, 34, j]`, with j = 1 (S-a6), 2 (S-a8), 3 (S-a4).
   - Season-stratified effect-size bootstrap: `[20261105, 35]`.
   - Calibration streams of `ma_audit.fit_lambda` / `calibrate`: 500 + 10i + j (ρ × k grid), 900 + i (fresh-draw generators), 600 + i (daily renewal).
2. **E1 wild bootstrap-t (a = 0 imposed)**
   - Restricted gamma quasi-likelihood fit: μ0 = b0/C with b0 = mean(C·L), the closed-form solution.
   - Bootstrap losses: L* = μ0 + v_week(L − μ0), clamped at 1e-7.
   - The stacked two-way sandwich is recomputed with x̄ held at the observed mean.
   - The p-value is two-sided in |t|.
3. **Season week-cluster wild bootstrap**
   - Null fit: the pooled (equal-intercept) affine fit.
   - Same clamp, and x̄ held fixed in the same way.
4. **Daily-renewal generator (Section 5.4)**
   - 35 seed days of exponential growth at the design's geometric-mean weekly growth g0.
   - The level is anchored to the shipped design's expected first-week count of each state run.
   - The generation interval is discretised to days 1–21.
   - Weekly R is redrawn in every replicate with the S1 law; the panel at each ρ follows the same rule.
   - SC1 is evaluated as written (a/v_R in [0.8, 1.25]). The two-week-mean estimand v_R̄ is reported in addition, as a descriptive ratio only.
5. **MCSE of bootstrap percentile endpoints:** batch method (10 equal batches; SD of the batch quantiles / √10).
6. **Implied-SD grid (Section 5.1)**
   - It also lists generation-interval variants (S-a1/S-a2; r = 1) and the within-week estimand E3 (r = 1).
   - These are labelled as separate estimands and are not mixed into a single number.
   - Envelopes are reported three ways: over all rows, over pre-specified laws only, and over rows whose b/b_P range contains the observed value.
7. **S-a5 threshold 1.02** reproduces the primary windows exactly; this is reported as such, not dropped.
8. **Fresh-draw generator season p-values:** the null is the generator's own MC mean of log a_s per season, exactly as `ma_calibration` does for the shipped generators.
