# Analysis changelog (revision of 2026-11-04)

This file lists the changes to the analysis code, data products and outputs, keyed to the review comments (M1–M7, m1–m11). Items marked *manuscript* are changes to the text and are listed only for completeness. Numerical effects are given in `analysis/revision-results-summary.md`.

## Pipeline and provenance (M5, m9)

- `code/run_all.sh` is now the single driver. It replaces `run_authors_pipeline.sh` and `run_ma_audit.sh`, which were removed. It runs four lanes in parallel on at most seven cores and then the renderer, reads only `data/`, and writes only `analysis/`, `tables/`, `figures/` and `logs/`. It also records the start, end and exit code of every step in `logs/run_timeline.log`. Any failing step, including a failed audit self-check, makes the run exit non-zero. Test: an injected self-check failure in a copy stopped the run with exit code 1 (`logs/faultcheck_run_all.log`). A full run from a clean copy of `code/` and `data/` exited 0 and reproduced every JSON output (`logs/cleancheck_run_all.log`).
- `ma_audit.py` no longer reads the external design file `C03`. The observed design statistics (a, b, b/a, b/b_P, the score contrast, g0, the mean kappa) are recomputed from the frozen data at start-up and stored in `meta.observed`.
- `ma_thinning_exact.py` no longer falls back to an external copy of the authors' thinning output. It stops with an error if the output is missing.
- Values cited in the text that were missing from `analysis/results.json` are now written there: g0 = 1.407, the shock scale s (S6 and S8), each generator's amplitude and count scale lambda, and the calibration targets.
- S8 now has its own self-check (recovery of the current-value error v_H with shocks) and amplitude-calibration check (miss < 2%), like the other generators. All self-checks pass, and their results are in `ma_audit.json` and in `results.json` under `audit_design.self_checks`.
- The claim that reweighting toward the upper count tail changes the recovery by less than a few percent is now backed by code. `ma_audit.py` post-stratifies the retained rows to the observed deciles of the count, and the relative shifts in recovery are written to `results.json` (`audit/*/tail_reweighting/*`).
- `analysis/results.json` keeps its structure. The new top-level block `audit_design` holds the ADEMP description, the self-checks and the generator settings. There are 382 new entries under `estimates` (484 → 866); no entry was removed.
- The step-1 rerun outputs, including the withdrawn calibration, are archived outside the package.

## Known-truth audit as a simulation study (M2)

- ADEMP description (aims, data-generating mechanisms, estimands, methods, performance measures) in `ma_audit.json` `meta.ademp` and in `results.json` `audit_design.ademp`.
- Calibration targets are labelled as such (`scenarios.*.calibration`). The fitted intercept matches the observed a by construction and is no longer listed among the checks.
- Monte Carlo standard errors for every summary: the SE of the mean (SD/sqrt(n)) and the SEs of the 2.5% and 97.5% limits (500 resamples of the replicates, `ma_models.summ_mc`). Simulated proportions have a binomial MC SE (`ma_models.binom_mcse`).
- New `code/ma_audit_extra.py`:
  - The observed score-contrast quantile of every generator over five independent streams of 200 panels, with the per-stream values, their SD and the pooled value with its MC SE. This replaces the single-stream values.
  - A transmission-share grid around S8: scale dependence gamma ∈ {0, 0.25, 0.5} by transmission share phi ∈ {0, 0.25, 0.5, 0.75, 1}, calibrated, plus an uncalibrated variant at gamma = 0.25. Each point records the compatibility of the observed a, b/b_P, b/a and contrast with the simulated 2.5–97.5% range.
- New figure `figures/fig_09_share_grid.pdf` (share grid; Fig. 3 of the manuscript, where it is renamed `fig_03_share_grid.pdf`) and table `tab_ma_share_grid`. The audit table shows the five-stream contrast quantile with its MC SE.

## Identification of transmission variance (M3)

- The audit table has two truth columns: the transmission relative variance v_R, and the total v_R + v_H, where v_H is the current-value error from reporting shocks. For each generator it reports a/v_R and a/(v_R + v_H) and the implied transmission relative SD (the observed a divided by the generator's recovery, square-rooted).
- The S8 values changed slightly against the shipped output and now equal the step-1 regeneration. The shipped `ma_audit.json` had been produced by an earlier version of the shock code.
- *Manuscript:* the statement of the bias direction and the 12% figure (step 3).

## China COVID-19 cross-regime check (M4)

`code/china_covid_empirical.py` was rewritten:

- **Units:** the 31 mainland provincial units. The 'Unknown' row is excluded. Hong Kong SAR and Macao SAR enter only as a sensitivity specification. Taiwan, China is not in the extract.
- **Selection:** rows are selected by national growth windows, using the US window rule applied to each national peak, instead of each row's own growth ratio (which selects on C_{t+1}). Window thresholds are 1.05 (primary), 1.02 and 1.10.
- **Generation interval:** period-specific and cited. Before 2022, Ganyani et al. (2020), 3.95/1.51 d. From 2022, Backer et al. (2022), 3.5/2.4 d. The first-version 5.2/1.7 d is a sensitivity specification.
- **Period:** the primary period ends with the week ending 2022-12-03, before routine mass testing ended. The full extract is a sensitivity specification.
- **Negative increments:** set to zero, or the affected weeks dropped (sensitivity).
- **Other sensitivities:** minimum count 20 and the pre-2022 and 2022 subsets. The first-version own-row rule is kept as a comparison only.
- **Panel and table:** the panel `data/panel_china_covid19_2020_2023.csv` was rebuilt (5,270 rows × 12 columns), and `tab_china_covid19` lists all eleven specifications. `data/DATA_MANIFEST.md` was updated.
- *Manuscript:* the claims, the move to the supplement and the data citation (step 3).

## Interval validity and analysis windows (M6)

- New `code/ma_coverage.py` estimates the coverage of the state × week pigeonhole bootstrap 95% intervals on known-truth panels for S1, S6 and S8 (400 panels × 200 draws each, binomial MC SEs). Coverage is assessed against the true total relative variance and against the estimator's own expectation. The run also reports the ratio of bootstrap SD to Monte Carlo SD. Results are in table `tab_ma_coverage`.
- A window-threshold sensitivity at 1.02 was added to the existing 1.10 (`ma_extensions.py`). The bootstrap streams of the existing specifications are unchanged, so their intervals equal the original ones exactly.
- *Manuscript:* the window rule is stated as a definition, and the "slightly conservative" sentence is revised (step 3).

## Tables and length (M7)

- `ma_render.py` writes the within-state and fixed-effects tables (S8, S9) and the other long tables as `longtable`, so they break across pages. Column widths are fixed so that no box is overfull in a test build.
- *Manuscript:* the final check that Tables S8 and S9 appear in the PDF, and the cuts in length (step 3).

## Minor comments

- **m3:** Table 1 states the zero-denominator rule: draws with a zero previous count are excluded from the ratio estimator, and the share of such draws is reported (0 at every m, so the rule never applies). The table also reports the MC SE of the observed-denominator value and a +0.5 variant, (C_t + 0.5)/(C_prev + 0.5), which keeps every draw.
- **m4:** the branching-Poisson baseline reports the simulated and analytic values and their ratio without normalising by the S0 recovery. The recovery shortfall is reported as a quantity in its own right.
- **m5:** the week-centred reduction is reported as "not defined" for generators with constant transmission (S0, S2, S2P).
- **m1, m2, m6, m7, m8, m10, m11:** *manuscript only* (step 3).
- **m9:** covered by M5.

## Other code fixes

- `ma_render.py`: the boxplot call uses the current Matplotlib `orientation=` argument (the old argument was deprecated). The coverage note is reworded. The Table 1 share column prints exact zeros as "0".
- `ma_extensions.py`: each section merges its results into `ma_extensions.json` instead of overwriting the file, so sections can be rerun one at a time. The calibration section, which was withdrawn, was removed.

## Re-review round

These analyses were added after a second review, once all earlier results were known.

- New `code/ma_tests.py`. It applies every declared test to every generator with one decision rule: the observed point value must lie inside the generator's 2.5–97.5% Monte Carlo range, and the observed score contrast between its 0.025 and 0.975 quantiles. Tests: b/b_P, b/a, thinning c/b_P, pooled/within slope ratio, week-centring reduction and the score contrast. Output: `analysis/ma_tests.json`, `tables/tab_ma_tests.tex`, and `tests/*` entries in `analysis/results.json`. Result: no generator passes every test. S8 fails the pooled/within slope ratio (MC range 0.92–1.15 against the observed 1.48) and the week-centring reduction (−6.1% to 2.3% against 34.7%).
- New `code/ma_audit_national.py`: generator S9, which is S8 with half of the transmission heterogeneity carried by a national calendar-week component. The component is the S4/S5 week factor and the share 0.5 is the S4/S5 design value; neither was tuned. S9 is run through the main audit (200 panels), the contrast (5 × 200 panels) and the coverage study (400 panels × 200 bootstrap draws). Output: `analysis/ma_audit_national.json` and `national/S9/*` entries in `analysis/results.json`. S9 also fails the pooled/within ratio (0.92–1.17) and the week-centring reduction (12.5–25.9%).
- `ma_tests.py` also records the season Wald statistic rescaled by the bootstrap-to-Monte-Carlo variance ratio of a from the coverage study (S1, S6, S8, S9; ratio 2.9–6.4). It is an illustration of how much the season test depends on interval calibration, not a corrected test.
- Sign convention: `revision_design.contrast()` returns power minus affine, so positive values favour the affine curve. The legacy key `contrast_affine_minus_power` keeps its name for compatibility, and the same value is also written as `tests/observed/contrast_power_minus_affine`.
- `ma_tests.py` appends the S9 row to `tables/tab_ma_coverage.tex`, which `ma_render.py` writes for S1, S6 and S8.
- `run_all.sh` runs `ma_audit_national.py` after the coverage study in lane 4, and `ma_tests.py` after the renderer.

## Gap-closing round

These analyses close the open analytic items of the final review (B-M1, B-M2, B-M5, B-m5 and the ADEMP replicates line). The decision rules were fixed before the results were seen, and every new value is written to `analysis/results.json` under `gaps/*` only. No existing entry changed: the diff of all 997 earlier estimates is empty, and the only other changes are the timestamp and the ADEMP replicates line.

- New `code/ma_calibration.py` (B-M2). It adds a two-way cluster-robust (state and week) sandwich for the stacked estimating equations of the pooled mean and the affine fit, with t critical values on min(G) − 1 = 41 df and b/a on the log scale. Its coverage of the estimator mean is checked on the 400 coverage panels of S1, S6, S8 and S9 (streams [SEED, 70, stream], as in `ma_coverage.py` and `ma_audit_national.py`). The applied interval was chosen by a rule fixed in advance: use the sandwich if its coverage for a is at least 0.925 under every generator, otherwise the leave-one-generator-out scaled pigeonhole bootstrap. The sandwich was chosen (coverage of a 0.925–0.9925, MC SE 0.004–0.013). The script also calibrates the season homogeneity test (3 df) against simulated nulls in which each season is centred at its MC mean, using both the sandwich statistic and the statistic with the bootstrap covariance. Output: `analysis/ma_calibration.json`, `tables/tab_ma_calibration.tex`.
- New `code/ma_share_tests.py` (B-M1). It reruns the share grid (γ ∈ {0, 0.25, 0.5} × φ ∈ {0, 0.25, 0.5, 0.75, 1} × national share ∈ {0, 0.5}, 30 cells) under all six tests of `ma_tests.py` with the same decision rule. It propagates the uncertainty in a by recalibrating each cell at the limits of the applied interval for a. Output: `analysis/ma_share_tests.json`, `tables/tab_ma_share_tests.tex`. Result: no cell passes all six tests, and every cell fails the pooled/within slope ratio.
- New `code/ma_heldout_sim.py` (B-M5). It computes the reference distribution of the held-out score contrast (leave one state out, leave one season out) on 200 synthetic panels from each of the 12 generators, with the fold design and score of the data. Output: `analysis/ma_heldout_sim.json`, `tables/tab_ma_heldout_sim.tex`. Fourteen of the 2,400 panels have one season fold without a fit; they are kept, and the count is recorded.
- New `code/ma_gaps.py`. It runs consistency self-checks on the three outputs and exits 1 on any failure: the decision rule is recomputed, panel counts are checked, pass counts must agree with the surviving list, every point estimate must lie inside its propagated range, and the fold counts are checked. It then replaces the `gaps/*` entries of `results.json` and writes the three tables above. It also writes `tables/tab_s03_thinning_summary.tex` (B-m5): six rows by subset and fit set, replacing the two-page per-q listing in the appendix; the full listing remains in `analysis/authors/thinning.json`.
- `run_all.sh` runs `ma_calibration.py`, `ma_share_tests.py`, `ma_heldout_sim.py` (seven workers each) and `ma_gaps.py` after `ma_tests.py`. The results do not depend on the number of workers.
- ADEMP replicates line: the string "grid 100" in `ma_audit.py` now reads "rho x k grid 100 per cell". `ma_render.py` writes a replicates line to `results.json` (`audit_design.ademp.replicates`) that also lists the three gap-closing designs.
- `ma_gaps.py` also records three derived checks. (i) The transfer check for the season calibration compares the observed sandwich standard errors with the generators' Monte Carlo SDs (`gaps/calibration/*`). (ii) The share-grid subset summaries cover the cells passing the slope ratio, crossover and contrast, and all cells with heterogeneity (`gaps/share_tests/*`). (iii) The tail-reweighting bound (B-m12) is the relative shift of the implied transmission relative SD, sqrt(unweighted/reweighted recovery) − 1, from the `audit/*/tail_reweighting` entries of `ma_audit.py`: −0.0003 (S1), 0.0166 (S4), 0.0025 (S6), 0.0143 (S7a), with maximum 0.0166 (`gaps/tail_reweighting/*`, 5 entries). In total 299 `gaps/*` entries.
- The notes of `tab_ma_calibration.tex` and `tab_ma_share_tests.tex` were rewritten as complete sentences, and the calibration note now refers to the manuscript section for the definition of the sandwich instead of repeating it. Numbers are unchanged.
- Clean-copy check: `run_all.sh` was rerun from an empty copy of `code/` and `data/` and exited 0. All 1,296 estimates and all generated tables are identical. The other JSON outputs differ only in timestamps, apart from the corrected ADEMP string in `ma_audit.json`; the clean-run versions were adopted (`logs/cleancheck_gaps_run_all.log`).
