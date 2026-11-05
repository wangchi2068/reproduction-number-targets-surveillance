# Revision results summary

This summary covers each review comment that needed new analysis. For each one it gives the new result, the numbers that changed against the original `analysis/results.json`, and whether the conclusion changes. All values are read from `analysis/results.json` and the JSON files it cites (`ma_audit.json`, `ma_audit_extra.json`, `ma_coverage.json`, `ma_extensions.json`, `china_covid19_empirical.json`). Intervals are 95%. Monte Carlo (MC) standard errors refer to the simulation, not to the data.

## Overview

- The empirical influenza estimates do not change. The pooled R_A fit is a = 0.01495 (0.01001–0.02092), b = 0.5600, b/a = 37.44 (18.93–69.92) and b/b_P = 1.257. Of the 484 entries in the original `results.json`, 44 differ:
  - 25 are S8 audit values, because the shipped S8 output came from an earlier version of the shock code;
  - 3 are week-centred reductions, now reported as "not defined";
  - 16 differ only by floating-point rounding (relative differences below 1e-12).
- There are 382 new entries (866 in total), and none was removed.
- The conclusions that change are listed below:
  - The S8 contrast quantile is now 6.9% (MC SE 0.8%). It replaces the range 4.5–8.0% (M2).
  - The intercept is not a transmission-variance estimate unless the transmission share is known: across the generators that are compatible with the data, a/v_R runs from 1.04 to 3.85 (M2, M3).
  - The China result no longer supports a "pattern the framework predicts" (M4).

## M5 and m9: regenerability and provenance

- One command, `bash code/run_all.sh`, regenerates every output. From a clean copy of `code/` and `data/` it exited 0 in about 17 minutes and reproduced every JSON output (`logs/cleancheck_run_all.log`). With a self-check made to fail it exited 1 (`logs/faultcheck_run_all.log`).
- Text values without a source are now in `results.json`:
  - the window growth factor g0 = 1.4072;
  - the S8 shock scale s = 0.1330 (shipped output 0.1341; text 0.133);
  - each generator's amplitude, count scale lambda and calibration targets.
- All self-checks pass, including the new S8 check that the current-value error is recovered with shocks (`audit_design.self_checks`).
- Upper-tail reweighting (post-stratification of the retained rows to the observed count deciles) changes the recovery a/v_R by +0.06% (S1), −3.2% (S4), −0.5% (S6) and −2.8% (S7a). The text's claim of "a few percent at most" is therefore supported.
- **Conclusion:** unchanged.

## M2: the audit as a simulation study

- ADEMP description: `audit_design.ademp`. MC SEs of the mean and of both interval limits are attached to every audit summary.
- Observed score-contrast quantile, pooled over five independent streams of 200 panels (MC SE in parentheses). Originals are from a single stream.

| Generator | Original | Five-stream pooled (MC SE) | SD across streams |
|---|---|---|---|
| S0 Poisson | 0.20 | 0.195 (0.013) | 0.012 |
| S2 NB reporting | 0.19 | 0.152 (0.011) | 0.039 |
| S2P branching Poisson | 0.205 | 0.171 (0.012) | 0.021 |
| S1, S3, S4, S5, S6 | ≤ 0.01 | 0.000–0.003 | ≤ 0.005 |
| S7a scale dependence 0.25 | 0.43 | 0.495 (0.016) | 0.027 |
| S7b scale dependence 0.5 | 1.00 | 0.987 (0.004) | 0.007 |
| S8 combination | 0.08 (text: 4.5–8.0% across streams) | 0.069 (0.008) | 0.024 |

- Transmission-share grid around S8, calibrated at each point (`fig_09`, `tab_ma_share_grid`):
  - The observed a, b/b_P, b/a and contrast all lie inside the simulated 2.5–97.5% ranges only at four points: gamma = 0.25 with phi = 0.25, 0.5 or 0.75, and gamma = 0.5 with phi = 0.5.
  - At these four points a/v_R is 3.85, 1.75, 1.04 and 1.55.
  - With gamma = 0, every phi fails: the contrast quantile is ≤ 0.011, and b/b_P is about 1.14–1.15 against the observed 1.257.
  - The uncalibrated variant (amplitude 1) fails at every phi.
- **Conclusion:** changed in emphasis. S8 remains compatible, with a contrast quantile of 6.9% (MC SE 0.8%), so it is not rejected at the 5% level. But compatible generators place between 25% and 75% of the total relative variance in transmission. The transmission share is therefore not identified by these summaries, and the intercept overstates transmission variance by a factor of 1.0 to 3.9 among compatible mechanisms.

## M3: what the intercept estimates

Recovery against transmission (a/v_R) and against the total (a/(v_R + v_H)), with the implied transmission relative SD (observed a divided by the recovery, square-rooted):

| Generator | a/v_R | a/total | Implied transmission SD |
|---|---|---|---|
| S1 | 1.007 | 1.007 | 12.2% |
| S3 | 1.017 | 1.017 | 12.1% |
| S4 | 1.115 | 1.115 | 11.6% |
| S5 | 1.149 | 1.149 | 11.4% |
| S6 (shocks, phi = 0.5) | 2.016 | 1.009 | 8.6% |
| S7a | 0.684 | 0.684 | 14.8% |
| S7b | 0.434 | 0.434 | 18.6% |
| S8 | 1.752 | 0.854 | 9.2% |

- Across generators the implied transmission SD runs from 8.6% to 18.6%, against the 12% that follows from reading a as v_R.
- S8 overstates transmission variance, with a/v_R = 1.75 (1.47–2.04 over panels), and understates the total, with a/total = 0.854.
- S8 values against the shipped output:
  - a: 0.01501 → 0.01478;
  - b/a: 40.85 → 42.15 (30.81–58.24);
  - b/b_P: 1.291 → 1.306 (1.081–1.520);
  - s: 0.1341 → 0.1330.

  These now equal the step-1 regeneration. The observed b/a (37.44) and b/b_P (1.257) remain inside the S8 ranges.
- **Conclusion:** the bias direction for S8 is "overstates transmission variance". The 12% figure must be given with its 9–19% range.

## M4: China COVID-19

| Specification | N rows (units) | a | b/a | b/b_P |
|---|---|---|---|---|
| First version (own-row rule, GI 5.2/1.7, 34 rows) | 896 (34) | 0.210 (0.085–0.500) | 23.2 (−3.0 to 99.9) | 6.16 (−1.84 to 13.37) |
| New primary (national windows, period GI, 31 mainland units, to 2022-12-03) | 711 (31) | 0.133 (0.074–0.223) | 13.2 (1.3–42.3) | 0.54 (0.05–2.78) |
| GI 5.2/1.7 instead of the period GI | 711 (31) | 0.298 (0.100–0.700) | | 4.62 (−0.08 to 9.41) |
| Own-row growth rule instead of national windows | | 0.054 (0.028–0.107) | | 2.18 (−0.40 to 4.61) |
| Full period to 2023-01-21 | | 0.146 (0.082–0.240) | | 0.60 (0.07–3.01) |
| Window threshold 1.02 / 1.10 | | 0.144 / 0.135 | | 0.54 / 0.49 |
| With Hong Kong SAR and Macao SAR | | 0.130 (0.078–0.210) | | 0.58 (0.06–2.57) |
| Pre-2022 / from 2022 | | 0.132 / 0.119 | | 1.29 / 0.47 |

- **Conclusion:** changed. The first-version slope ratio of about 6 came from the 5.2/1.7-day interval and the outcome-dependent row selection. Under the redesigned specification the slope-ratio interval (0.05–2.78) contains both 1 and values well below it, so the China data are uninformative about b/b_P > 1. The intercept, a ≈ 0.13 with interval 0.07–0.22, is about nine times the influenza value. This difference is robust to the specification choices other than the generation interval, which moves it to 0.30.

## M6: interval validity and windows

Coverage of the pigeonhole bootstrap 95% intervals (400 known-truth panels × 200 draws; binomial MC SE at most 0.004 near 1 and 0.016 at 0.89):

| Generator | a vs E(a) | b vs E(b) | b/a vs ratio of means | a vs true total | Bootstrap SD / MC SD (a, b, b/a) |
|---|---|---|---|---|---|
| S1 | 1.000 | 0.995 | 1.000 | 1.000 | 2.52, 1.69, 2.01 |
| S6 | 0.998 | 0.995 | 0.998 | 0.998 | 1.76, 1.61, 1.72 |
| S8 | 1.000 | 0.998 | 1.000 | 0.893 (0.016) | 1.71, 1.76, 1.82 |

- The intervals are conservative: the bootstrap SD is 1.6 to 2.5 times the MC SD of the estimator.
- Under S8 the shortfall in coverage of the true total (0.89) reflects the bias of a for the total (a/total = 0.857), not interval width.
- Window threshold:
  - At 1.02 the windows select the same 1,748 rows as at 1.05, so b/a is 37.44 (18.80–69.98) and b/b_P 1.257. The interval differs only through its own bootstrap stream.
  - At 1.10, a = 0.01406, b/a = 38.63 (20.58–70.33) and b/b_P = 1.255.
- **Conclusion:** unchanged for the influenza estimates; the intervals are valid and conservative. The text's "slightly conservative" should read "conservative, with SDs 1.6–2.5 times too wide".

## m3: Table 1

- The values are unchanged: the observed-denominator ratio is 2.42 at m = 10, 1.59 at m = 50, 1.51 at m = 200 and 1.49 at m = 1,000.
- The MC SE of the m = 10 observed-denominator MSE is 0.0048 (value 0.5145), or about 0.023 on the ratio.
- The share of zero denominators is 0 at every m, so the exclusion rule never applies.
- The +0.5 variant gives 0.4253 (MC SE 0.0028) at m = 10, a ratio of 2.00, and lies within 2.5% of the main value from m = 50 upward.
- The independent check in the review (2.392) lies about 1.4 MC SEs below 2.42.
- **Conclusion:** unchanged.

## m4: branching-Poisson baseline

- The ratio of simulated to analytic slope is 0.929, 0.932 and 0.927 at rho = 0.1, 0.25 and 0.5, and 0.907 at rho = 1. So gamma FGLS recovers about 93% of the analytic slope on this design.
- The S0 recovery of b/b_P is 0.93, the same shortfall.
- The unnormalised comparison at rho = 0.25 is simulated 0.802 (0.734–0.878) against analytic 0.861.
- **Conclusion:** "exact alignment" becomes a 7% finite-sample shortfall of the slope, quantified on this design.

## m5: week-centred reduction

The reduction is reported as "not defined (constant transmission)" for S0, S2 and S2P. These replace −30.74 (−34.19 to 12.89), 0.44 (−20.56 to 35.68) and −2.64 (−65.20 to 22.66). The other generators are unchanged, except S8: −0.018 (−0.061 to 0.023), previously −0.017.

## M7: tables

The long tables (including S8 and S9) are rendered as `longtable` and break cleanly in a test build. The check in the compiled manuscript belongs to the manuscript revision.

## Items without new analysis

M1, m1, m2, m6, m7, m8, m10 and m11 are text changes.

## Gap-closing round: calibrated interval, six-test share grid, held-out reference

All values are in `analysis/results.json` under `gaps/*`, with the full outputs in `ma_calibration.json`, `ma_share_tests.json` and `ma_heldout_sim.json`. No earlier entry changed.

**Interval for a under crossed state × week dependence (B-M2).** The two-way cluster-robust sandwich covers the estimator mean of a in 0.9925 (S1), 0.925 (S6), 0.925 (S8) and 0.9825 (S9) of 400 panels (MC SE 0.004–0.013). The figures for b are 0.9175–0.945 and for b/a 0.915–0.97. Its mean SE is 0.94–1.40 times the MC SD of a, against 1.71–2.52 for the pigeonhole bootstrap. By the fixed rule the sandwich is the applied interval: a = 0.01495 (0.01085–0.01906), b = 0.560 (0.430–0.690), b/a = 37.4 (24.1–58.3). The primary bootstrap gave 0.0101–0.0206 and 18.4–72.4 (`primary/RA/*`). Under S8 the interval covers the true total relative variance in only 0.42 of panels (S9: 0.785), because a is biased for that target; the interval is calibrated for the estimand of a, not for v_R + v_H.

**Season homogeneity (B-M2).** With a common slope, the season intercepts are 0.0159, 0.0071, 0.0154 and 0.0180 (2022-23 to 2025-26). The sandwich Wald statistic is 17.8 (3 df). Its simulated size at the χ² critical value is 0.04–0.22 across generators, and the calibrated p-value is 0.005–0.045 (S1 0.0125, S6 0.045, S8 0.0275, S9 0.005). The statistic with the bootstrap covariance (5.52, nominal p = 0.138) has simulated size 0 and null 95th percentiles of 1.3–2.7, so its calibrated p is below 0.005 under every generator. **Conclusion changes:** the seasons differ, driven by the low 2023-24 intercept. The calibration assumes the generators' dependence; the observed sandwich relative SEs of the season intercepts (0.17–0.25) exceed the simulated MC SDs of log a_s (0.07–0.15), so the data are noisier than the generators and the calibrated p-values are best read as a range.

**Share grid under all six tests (B-M1).** None of the 30 cells passes all six tests. Pass counts across cells: b/b_P 25, b/a 26, thinning 24, pooled/within 0, week-centring reduction 3, score contrast 12. The pooled/within slope ratio (observed 1.476) lies above every cell's upper limit (largest 1.247, MC SE 0.028), as it does for every generator. The best cell (γ = 0.25, φ = 0.75, national share 0.5) passes five tests; its implied transmission relative SD is 11.5% (9.9–13.0% with the interval for a propagated). Across all cells with transmission heterogeneity the implied SD is 6.0–18.7% at the observed a and 5.1–20.5% propagated. The cells without a national component reproduce the earlier grid (γ = 0.25, φ = 0.25–0.75: 6.2–12.0%). **Conclusion changes:** the 6–12% range is conditional on the four tests used before, and no tested configuration of the share reproduces all six summaries; the range is a sensitivity statement, not an estimate. Nineteen of the 180 decisions lie within 2 MC SE of a limit, none of them for the five-test cell or the pooled/within test.

**Held-out contrast reference (B-M5).** The observed held-out contrasts (power minus affine) are −0.0020 (leave one state out) and 0.0046 (leave one season out). For the state folds the observed value lies at quantile 0 of the reference under every count-independent heterogeneity generator (S1, S3, S4, S5, S6; mean contrasts 0.011–0.024, MC SE ≤ 0.0006). It lies at 0.05 under S8 and 0.055 under S9, and inside the central range under S0, S2, S2P and S7a (0.125–0.475). Under S7b it lies at 0.995. The season folds discriminate less (S8 0.22, S9 0.145). **Conclusion changes:** the held-out state contrast is informative; it is incompatible with purely count-independent heterogeneity and at the margin for S8 and S9, rather than "not distinguishing" the forms. Fourteen of 2,400 panels had one season fold without a fit (S0 4, S2 4, S2P 6).

**Share-grid subsets (B-M1).** Nine cells pass the slope ratio, the crossover ratio b/a and the held-out contrast; thinning excludes none of them. Their implied transmission relative SD is 6.2–13.9% at the observed a and 5.3–15.5% with the interval for a propagated. The range is conditional on k = 0.3 and ρ = 0.25.

**Tail reweighting (B-m12).** Reweighting the synthetic count distribution to the observed upper tail changes the recovery ratio a/(v_R + v_H) by at most 3.3%. Because the implied transmission relative SD scales as the inverse square root of that ratio, its relative shift is −0.0003 (S1), 0.017 (S4), 0.003 (S6) and 0.014 (S7a), at most 1.7%, small against the ranges above.

**Reproducibility.** A clean-copy rerun of `code/run_all.sh` (empty directory, code and data only) exited 0 and reproduced all 1,296 estimates and every generated table exactly (`logs/cleancheck_gaps_run_all.log`).
