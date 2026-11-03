#!/usr/bin/env bash
# Rebuild the analysis, figures and tables from the downloaded data and the shipped
# vintage growth-null bands (analysis/vintage_r4_results.json).
#
# Scope note: this is NOT a fully end-to-end rebuild. The growth-null bands produced by
# the earlier vintage round are reused rather than regenerated; regenerating them
# requires the archived vintage workspace. Everything downstream of them is recomputed.
#
# Usage: bash code/run_all.sh [--fresh]
#   --fresh  delete cached analysis parts first (each script otherwise skips blocks already in its part file)
# Wall time on 8 CPUs: about 10 min without the growth-null rerun (see code/README.md).
set -euo pipefail
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [[ "${1:-}" == "--fresh" ]]; then rm -f ../analysis/parts/*.json ../analysis/parts/s02_draws_*.npy ../analysis/parts/s09_*.npz; fi
python3 s01_data_panel.py            # data panel, windows, point estimates
python3 s02_inference.py 2000        # two-way / state / week / season-week bootstrap, CGM
python3 s03_decomposition.py 1000 500 200   # FE decomposition, ICC, plateau test, CV, forest, theory b/a
python3 s03b_thinning.py 1000        # binomial thinning experiment (shared-count protocol)
python3 s03c_state_fe.py 1000        # state fixed effects with bootstrap
python3 s04_reporting.py             # reporting completeness, TOST, calendar-matched, laboratory proxy
python3 s05_seasons.py               # season intercepts, random effects, Q-profile, design surface
python3 s06_pathogens.py 2000        # cross-pathogen q1/q2, growth-null bands, Holm/BH
python3 s07_theory.py                # closed-form theory numbers
python3 s08_sim_apparent.py          # apparent-series simulation
python3 s09_inversion.py all         # round 3: inversion, multiscale, reporting regime, returns, Parag contrast (~6 min)
python3 s10_polish.py all          # final polish: lag-1 autocorrelation, Cowling GI bound, admission overdispersion, latent balance scale (~2.5 min)
python3 s99_collect.py               # -> analysis/results.json
for f in 1 2 3 4 5 6; do for l in en zh; do python3 make_fig$f.py --lang $l; done; done
for l in en zh; do python3 make_fig3_inversion.py --lang $l; python3 make_fig6_nomogram.py --lang $l; done
python3 make_tables.py
python3 make_tables_r3.py ../../paper_en/tables
python3 make_tables_r3.py ../../paper_zh/tables --lang zh
python3 make_tables_polish.py ../../paper_en/tables
python3 make_tables_polish.py ../../paper_zh/tables --lang zh
python3 make_tables_s11.py
