# Code: revision experiments for "case-count bottleneck"

To rebuild the analysis with one command, run `bash run_all.sh --fresh` from this folder (`code/`).
- Scope: `run_all.sh` rebuilds every analysis number, figure and table **from the downloaded data and the shipped vintage growth-null bands**. It does not regenerate those bands; they come from the earlier vintage round and are committed as `../analysis/vintage_r4_results.json`.
- Inputs: `../data/external/` (see `../data/DATA_MANIFEST.md`), plus the vintage growth-null bands, which `s06_pathogens.py` locates automatically in `../analysis/vintage_r4_results.json` or `../vintage_2026-10-20/analysis/r4_results.json`.
- Outputs:
  - `../analysis/parts/*.json`: per-script results with provenance;
  - `../analysis/results.json`: merged results, one record per number with its script, inputs, seed and mode (measured, simulated or analytic);
  - `../figures/{en,zh}/`: figures;
  - `../tables/{en,zh}/`: tables.
- The growth-null bands are the one intermediate that is reused rather than recomputed. No manual path edit is required; if neither location exists, `s06_pathogens.py` raises a descriptive error naming both expected paths. Tables are in `analysis/tables/`; figures are released with the paper.
- Base seed: 20261020. Each script adds its own offset, and each bootstrap records its seed.

## Scripts

| Script | Purpose | revision_plan items | Runtime (8 CPU) |
|---|---|---|---|
| `cc_common.py` | Paths, loading, statistics (G, G1, R_A), windows (calendar S=3; rule v2 S=4), iterated FGLS `fit_stat`, two-way pigeonhole bootstrap weights, provenance helpers `Part`/`Log` | — | — |
| `s01_data_panel.py` | Builds the state-week panel, national series, windows, point fits | M2, M3④, m8, m12 | 5 s |
| `s02_inference.py [B]` | Bootstrap schemes (two-way, state, week, season-week), CGM two-way sandwich, positive-share | M5, m14 | ~60 s |
| `s03_decomposition.py [B_dec B_for N_pb]` | Week/state-season FE decomposition, ICC, plateau vs power-law test with power, leave-one-out CV, b/a forest, theory b/a with b-delta | M3, M4 | 44 s |
| `s03b_thinning.py [B]` | Binomial thinning experiment under shared-count protocol (q in [0.8..0.1], bootstrap intervals, leverage) | M4, m2 | ~35 s |
| `s03c_state_fe.py [B]` | Within-state fixed effects with bootstrap, mean/median/weighted state intercepts, paired ratio | M2, S4 | ~45 s |
| `s04_reporting.py` | Reporting regimes, Δa with paired bootstrap, TOST with pre-set margins and power, excess at three v_R, calendar-matched May–Oct, laboratory proxy | M6, m8, m9 | 28 s |
| `s05_seasons.py` | Season intercepts, DerSimonian–Laird, Q-profile, sqrt(2/(S−1)), design surface | m2, m3 | ~20 s |
| `s06_pathogens.py [B]` | Cross-pathogen q1/q2 with growth-null bands and Holm/BH correction, FluSurv-NET status, ILINet reuse | M7③④ | 28 s |
| `s07_theory.py` | Closed-form quantities (m×, v_H, m*, B†, T†, lever terms, 0.505, Corollary 2) plus a Monte Carlo check | M1③④, m1, m4, m5 | <5 s |
| `s08_sim_apparent.py` | Apparent-series simulation with known R_t, generation interval and admission delay; maps fitted a, b to v_R, A | M2③ | 12 s |
| `s09_inversion.py [sim\|A\|B\|C\|D\|E\|all]` | Round 3. Simulation matched to the data; inversion of the calibration line to v_R^impl, with bootstrap and Monte Carlo uncertainty combined; sensitivity over 11 settings; calibrated b/a (A). Multiscale share, ICC, SD by clustering, J-pooling (B). Mandatory vs voluntary reporting: TOST margin share, Δa, IQR ratios, participation (C). Diminishing returns (D). Comparison with the Poisson renewal variance of Parag et al. 2022 (E). Writes `inversion`, `multiscale`, `reporting_regime`, `returns` and `parag` parts plus `s09_*.npz` | round 3 | 369 s (`all`; sim stage 329 s) |
| `s10_polish.py [sim\|fit\|all]` | Final polish. Lag-1 autocorrelation of within-block deviations in data vs simulated panels (M0, AR φ=0.25/0.5, every s09 setting); four new data-matched settings (true GI 2.9 d = lower end of Cowling et al. 2009, admission VMR 1.5 and 2.0, AR φ=0.25); latent balance scale with interval; 90% interval of the v_R^impl ratio for a TOST on [0.8, 1.25]; τ permitted by the 5.3% completeness bound and Proposition 6 levers at that τ. Writes the `polish` part | final polish | ~130 s sim + 2 s fit |
| `s99_collect.py` | Merges the parts into `results.json`. For round-3 parts only, it adds `inputs_sha256` | M8 | <1 s |
| `fig_style.py` | Figure style that imitates Eales & Riley (2024): Dark2, theme_bw, bold A–D, Liberation Serif plus Noto Serif CJK SC; font, mixed-script and collision gates | M9 | — |
| `make_fig1.py` … `make_fig6.py --lang en/zh` | One script per figure | M9 | few s each |
| `make_fig3_inversion.py`, `make_fig6_nomogram.py --lang en/zh` | Round-3 Fig 3 (simulation-based inversion) and Fig 6 (design nomogram), from the s09 parts | round 3 | few s each |
| `make_tables.py` | LaTeX tables S1–S7 (en/zh) from `results.json`, plus `tables/manifest.json` | — | <1 s |
| `make_tables_r3.py` | LaTeX tables S8 and S9, updates S1/S2/S3/S6 and `manifest.json` | — | <1 s |
| `make_tables_polish.py` | LaTeX table S10, updates `manifest.json` | — | <1 s |
| `make_tables_s11.py` | LaTeX table S11 (thinning & FE), updates `manifest.json` in paper_en and paper_zh | — | <1 s |

## Reused code (sources and changes)

- **`s06_pathogens.py`** ports `growth_weeks` and `build_rows` from the user's `experiment/original/code/r3_c_pathogens.py` and `r3_common.py`. Changes:
  - the state-cluster BCa bootstrap is replaced by the two-way bootstrap;
  - the null bands come from the vintage growth-null rerun instead of hard-coded constants;
  - Holm and BH corrections are added;
  - a rule-window combination is added.
- **`s05_seasons.py`** copies `dl` and `q_profile` from the user's `r4_e_seasons.py`. Changes:
  - the per-season SE is now two-way;
  - S = 3 and S = 4 are reported separately.
- **`cc_common.fit_stat`** reimplements the iterated FGLS of the user's `common.py`. It reproduces a = 0.099055 against 0.099060.
- **Growth-null bands** come from the user's `r4_a_growth_null.py`, rerun unchanged on the 2026-10-20 vintage (`../vintage_2026-10-20/code/`, `run_all.sh` there; 1,386 s on 8 CPU). Copies of its cell output and log are in `../logs/`.
- **Classification variants (M7②)** are read from the vintage `r4_results.json` (`C_classification`, `C_posthoc`). They are not recomputed.

No script calls an LLM API.
