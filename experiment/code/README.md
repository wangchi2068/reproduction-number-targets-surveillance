# Code: inferential targets and case volume in reproduction-number surveillance

Every number in `analysis/`, `tables/` and `figures/` is produced by the scripts below; nothing is typed by hand.

## Environment

Python 3.11 with the pinned packages in `requirements.txt` (`python -m pip install -r requirements.txt`). CPU only; each step runs as one process with at most two threads. The full sequence needs well under 8 GiB of memory.

## Input data

`data/external/nhsn_hrd_vintage_2026-10-20.csv.gz` (NHSN Hospital Respiratory Data, weekly influenza admissions by jurisdiction, frozen vintage). `revision_data.py` builds the analysis samples inside the rule-based season windows (R_A: 1,748 state-weeks with C_t >= 5, statistic from the weekly growth ratio; G: the renewal ratio C_{t+1}/D_t with D_t >= 5) and records their hashes in `analysis/authors/input_audit.json`.

## Run order

| Step | Command (from `code/`) | Output | Seed | Wall time |
|---|---|---|---|---|
| 1 | `bash run_authors_pipeline.sh` | `analysis/authors/` (pooled fits, state fixed effects, block bootstraps, theory Monte Carlo, secondary fits, shared-count thinning) | authors' seeds in each script | about 9 min |
| 2 | `bash run_ma_audit.sh` | `analysis/ma_audit.json`, `analysis/ma_audit_draws.npz` | 2026110351 | about 7 min |
| 3 | `python ma_thinning_exact.py` | `analysis/ma_thinning_exact.json` | 2026110353 | about 1 min |
| 4 | `python ma_extensions.py table1 ratios within sensitivity seasons validation` | `analysis/ma_extensions.json` | 2026110352 | about 14 min |
| 5 | `python ma_fe_decomp.py` | `analysis/ma_fe_decomp.json`, `analysis/fe_rev2/{L1,block2,block3}/` | 2026110352 | about 2 min |
| 6 | `python ma_render.py` | `tables/*.tex`, `tables/manifest.json`, `figures/*.pdf`, `figures/*.png`, `figures/manifest.json`, `analysis/results.json` | none | under 2 min |

Step 2 needs only the data. Steps 3 to 5 read `analysis/authors/` or the data. Step 6 reads everything. The original design-power stage (`revision_design.py`) is not run because its table is withdrawn; the calibration stage is withdrawn as well (see below).

`ma_extensions.py` takes section names (default `table1 ratios within sensitivity seasons validation`). Each section merges its result into `analysis/ma_extensions.json`, so sections must run one after another, never concurrently. `calibration` is still available as a section name but is not in the default list.

## Primary bootstrap

All pooled R_A intervals in the tables come from one bootstrap: the RA frame with the 3.6/1.6-day generation interval, B = 2000, `numpy.random.default_rng([2026110352, 1])`, with a, b and b_P computed on the same state x week weights. `ma_extensions.py` (`ratios`) stores it; the pooled rows of the sensitivity, within-state, generation-interval and fixed-effects (L = 1) analyses use the same weights and agree with it exactly. `ma_render.py` checks the point estimates before it substitutes the primary interval into the authors' Table S1, the week-centring rows and the fixed-effects pooled rows, and records the result as `primary/RA/*` in `analysis/results.json`.

Reproduction record: the authors' own pooled run (their rng stream) gives b/a = 37.44 [18.93, 69.92]. It is kept under `pooled/` in `analysis/results.json` as evidence that their pipeline reproduces, and is not used in any table.

## Scripts

Authors' pipeline, used unchanged except for table wording in `render_revision.py`: `revision_data.py`, `revision_stats.py` (affine and power loss fits, gamma FGLS), `revision_empirical.py`, `s03b_thinning.py`, `s03c_state_fe.py`, `revision_theory.py`, `revision_secondary.py`, `revision_calibration*.py`, `revision_design.py`, `render_revision.py`, `run_revision.py`. The remaining authors' files (`build_papers.py`, `check_manuscript.py`, `prune_bibliography.py`) belong to their manuscript build and are not used here.

Added for this analysis:

- `ma_models.py`: shared helpers. Weighted two-way (state and week) demeaning by alternating projections; a precision-weighted within slope of the loss on 1/D with one-step or iterated (damped) gamma working weights; bootstrap and Monte Carlo summaries.
- `ma_audit.py`: known-truth audit on the real R_A design (state, week, observed starting counts). Eleven data-generating processes: Poisson counts with constant R (S0); independent lognormal heterogeneity in R (S1); negative-binomial offspring with binomial reporting (S2); a branching process with Poisson offspring and the same reporting (S2P); S1 plus S2 (S3); state and week structured heterogeneity (S4); S4 plus the S3 offspring law (S5); S3 plus independent state-week reporting shocks (S6); S3 with heterogeneity that falls with expected count, exponent 0.25 (S7a) and 0.5 (S7b); S3 combining scale-dependent heterogeneity (exponent 0.25) with reporting shocks (S8 Combo). Every simulated growth window is anchored to the geometric mean of the observed counts and scaled by a per-generator factor lambda chosen by bisection so that the median retained denominator equals the observed 45; the heterogeneity scale sigma (and the S6/S8 shock scale s) is then calibrated to the observed intercept. Each replicate applies the analyst's own pipeline (simulated C_t as the count scale, C_t >= 5 filter, affine and power fits, state fixed effects, week centring, shared-count thinning) and compares the fitted intercept with the true relative variance on the retained rows. Also a 4 x 5 grid over reporting probability rho and offspring dispersion k, and the branching-Poisson baseline (analytic on the observed rows and simulated).
- `ma_thinning_exact.py`: exact expected thinned loss by summation over binomial supports (continuity correction +0.5 and positive-part truncation included), per-q increment rates, least-squares leverage of each q, and an independent Monte Carlo check.
- `ma_extensions.py`: observed-denominator Monte Carlo for Table 1; the primary bootstrap and the paired generation-interval comparison (3.6/1.6 vs 3.2/1.6 vs 3.2/2.1 days, common resampling indices) with b relative to the branching-Poisson reference; within-state, within-week and two-way slopes; filter, continuity, window and season sensitivity; a season-stratified model with a common slope and homogeneity tests on the season intercepts (Wald test with bootstrap covariance, Cochran's Q, minimum detectable difference at 80 % power); held-out comparisons with paired standard errors and bootstrap intervals.
- `ma_fe_decomp.py`: decomposition of the pooled versus state fixed-effects crossover count into a slope ratio and an intercept ratio, on the primary weights for L = 1 (B = 2000) and on week-block bootstraps for L = 2 and 3 (B = 1000). It calls the authors' `s03c_state_fe.run_analysis` and writes their output format to `analysis/fe_rev2/`.
- `china_covid_empirical.py`: cross-pathogen and cross-regime validation on China COVID-19 (34 provinces, 2020-01-25 to 2023-01-21, JHU CSSE). Builds the weekly panel, keeps rows with own growth g >= 1.05 and C_t >= 5, applies the Wallinga-Lipsitch transform with a gamma generation interval of mean 5.2 d and SD 1.7 d, fits the same gamma quasi-score loss model, and bootstraps over provinces x weeks. Writes `analysis/china_covid19_empirical.json` and the processed panel.
- `ma_render.py`: renders all tables (long tables as `longtable`) and figures from the JSON outputs and writes `analysis/results.json` (`estimates`: name, point, lower, upper, ci95, interval, source).

## Withdrawn: calibration (2026-11-03)

The calibration section (the map from the fitted intercept to the variance of R, the authors' Table S6 and its diagnostics) is withdrawn and no longer appears in `analysis/results.json`, the tables or the figures. At the nearest grid cell the calibration simulator did not reproduce the data it was calibrated to: fitted slope 0.413 against the observed 0.560, retained sample size 1,568 against 1,748, and national lag-1 autocorrelation 0.244 against 0.494. The inversion therefore has no support from its own checks. The earlier outputs are kept in `analysis/superseded/` for the record.

## Rerun record (second revision, 2026-11-03 to 2026-11-04)

- Audit rerun with count alignment and the generators S2P, S6, S7a and S7b; grid rerun with the same alignment.
- Primary bootstrap introduced and propagated to every pooled row; fixed-effects decomposition added; season homogeneity tests added; descriptive season ratio kept but labelled as descriptive.
- Calibration removed from `analysis/ma_extensions.json` (backup `analysis/superseded/v2_ma_extensions_before_rev2.json`).
- Logs: `logs/rev2_audit.log`, `logs/rev2_extensions.log`, `logs/rev2_fe.log`.
- Rendering only (2026-11-04, no rerun): the pooled-to-state-FE slope row of the within-slope table now takes `fe/L1/slope_ratio` from `analysis/fe_rev2/L1/` (the same primary 2,000-draw weights as the pooled row), so the table and the fixed-effects decomposition report one interval. The paper stage copies the tables through its own staging script, which typesets notation (for example C_{t+1}), adds the note that G rows are descriptive and not audited, and sets the tables that follow a long table in place; numerical cells are not changed.

## Notes

- Bootstrap intervals are percentile intervals over resampled states and weeks (the authors' convention). Ratio intervals report the share of non-positive denominators; tables flag a row when it exceeds 2.5 % or when draws failed to converge.
- The undamped iterated two-way FGLS oscillates with period two; the reported iterated estimate uses a damping factor of 0.5 (1,962 of 2,000 draws converged). The one-step estimate is a linear weighted least-squares estimator and always exists.
- `analysis/superseded/` keeps a first audit run in which the analyst's count scale was the real rather than the simulated C_t (an errors-in-variables design flaw), the pre-alignment audit, and the withdrawn calibration. None of it is used.
