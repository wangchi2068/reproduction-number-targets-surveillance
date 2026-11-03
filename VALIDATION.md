# Validation record

This file records what was actually executed for the corrective revision, what passed, and what remains limited. It is not a peer-review certificate and does not assign a journal decision.

Final release commit: `3a5f732` (corrective analysis `342bd41`, wording `01e515c`, build-report fix `3a5f732`).

## Environment

- Windows 10, Python 3.11.9 (venv at `D:/trae/venv`), numpy 1.26.4, scipy 1.16.3, pandas 2.3.3, matplotlib 3.10.9, threadpoolctl 3.6.0, pypdf 6.14.2.
- Typesetting: Tectonic. Fonts: Times New Roman (Latin); SimSun (Chinese).
- All analyses were run single-threaded (`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`) with bytecode writing disabled.

## Analyses executed at publication budget

| Stage | Budget | Outcome |
|---|---|---|
| Data audit and sample construction | full raw snapshot (21,641 rows, 51 jurisdictions) | RA n=1,748; G n=1,756; 0 incomplete rows in windows |
| Pooled inference and centering | 2,000 bootstrap draws; 1,000 for centering/block | converged; see `empirical.json` |
| State fixed effects | 1,000 bootstrap draws | 999/1,000 valid fits; 1 rejected failure retained |
| State fixed effects, 2- and 3-week blocks | 1,000 each | 1,000/1,000 and 999/1,000 valid |
| Binomial thinning | 1,000 outer × 50 inner, 400 point realisations | 1,049,998/1,050,000 inner fits converged; precision criterion met |
| Conditional calibration | 24 scenarios × 5 variance values × 200 panels = 24,000; 2,000 paired draws | 0 RA/G fit failures; 0 nonpositive implied variances; 27 (RA) and 81 (G) of 48,000 draws outside grid, retained |
| Design illustration | 500 replicates × 3 stages × 12 cells | null size near 0.05; power saturated near 1 (reported as a limitation) |
| Secondary sensitivities | 1,000 bootstrap draws | season and participation summaries |

Production wall-clock: calibration 28.4 min for 24,000 panels plus merge (29.8 min total); thinning initial pass approximately 23 min.

## Checks passed

- `pytest tests` — 53 tests passed (statistics, observation model, data and validation).
- `check_manuscript.py` — 0 missing citations, 0 missing figures, 0 missing tables, 0 undefined numeric macros, 0 duplicate labels in both languages.
- Both `numbers.tex` files are byte-identical.
- English PDF: 25 pages, 0 overfull hbox, 0 undefined references, 0 missing characters.
- Chinese PDF: 26 pages, 0 overfull hbox, 0 undefined references, 0 missing characters.
- Bibliography pruned to the 37 cited entries, byte-identical across languages; DOI verification evidence in `references/`.

## Corrections that change reported results

1. State-intercept bootstrap now weights by sampled state multiplicity; omitted states are not zero. The within-state ratio interval changes from `[3.3, 51.4]` to `[3.3, 47.8]`, and the pooled-to-within paired ratio from `[0.93, 7.75]` to `[1.04, 7.74]`.
2. The `D_t` sensitivity now uses the weighted lagged scale; ratio `17.71 [3.27, 47.81]` instead of a duplicate of the `C_t` fit.
3. Week-only and state-plus-week diagnostics report convergence status; the state-plus-week fit did not converge and is not given a slope or interval.
4. Thinning intervals now reflect the nested average-response estimator: pooled rate `0.405 [0.366, 0.441]`, excluding `q=0.1` `0.485 [0.439, 0.528]`, ratio to the approximate elasticity reference `0.909 [0.810, 0.999]`.
5. Calibration observation noise is generated directly from latent intensity; the previous double-noise construction is gone.
6. Cross-validation centres and scoring weights are estimated within training folds; the old residual-bootstrap platform p-value is not carried forward.
7. Multiscale intercept reductions use one common normalisation: national-week centring reduces the intercept by 34.7% (RA) and 33.3% (G), not the previously reported 41–51%.

## Known limitations

- Windows and the primary statistic were chosen in earlier rounds; not preregistered.
- The calibration simulator reproduces the empirical intercept but not the empirical slope (0.413 vs 0.560 for RA), sample size (1,568 vs 1,748) or national lag-one autocorrelation (0.244 vs 0.494). These are reported in `calibration.json` response checks, not tuned away.
- The design exercise saturates; it cannot estimate the number of seasons needed in reality.
- Two of 1,050,000 thinning inner fits failed; the failure policy makes any affected level average undefined, and this is retained in diagnostics.
- The independent clean-copy rerun covers the empirical, fixed-effect, block, thinning, theory, secondary and design stages, with the production and clean values identical to the last printed digit for the compared quantities. Calibration was executed once at full budget with its own 30 implementation tests and deterministic-seed reproducibility checks; it was not rerun a second time from the clean copy.

## Final build report

`validation/build/build_report.json` records both languages: English 25 pages, Chinese 26 pages, no unresolved question marks, and zero overfull boxes, undefined references, missing characters, LaTeX errors or duplicate destinations in the final pass.
