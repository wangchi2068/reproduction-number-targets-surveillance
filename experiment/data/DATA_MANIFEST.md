# DATA_MANIFEST: case-count-bottleneck experiment

Snapshot date: 2026-10-20 (NHSN vintage). All files are kept under `experiment/data/`. Downloaded files are in `data/external/`, derived and legacy files are at the top level. Checksums are SHA-256 hashes of the bytes on disk; `code/rev_verify_data.py` recomputes them and is the first step of `code/run_all.sh` (a mismatch in a file the pipeline reads stops the run with exit code 4).

Corrected on 2026-11-05: the JHU file size and hash were wrong in the earlier version, the China panel had no hash, the derived NHSN panel was attributed to a script that is not in the repository, and two files were listed as inputs although no script reads them.

## Inputs read by the pipeline (enforced)

| File | Source URL | Access | Size | Rows / records | Fields | Time range | SHA-256 |
|---|---|---|---|---|---|---|---|
| `external/nhsn_hrd_vintage_2026-10-20.csv.gz` | https://data.cdc.gov/Public-Health-Surveillance/Weekly-Hospital-Respiratory-Data-HRD-Metrics-by-Ju/ua7e-t2fy (full CSV export, dataset id `ua7e-t2fy`) | HTTP 200, 2026-10-20 | 4,794,721 B (gz) | 21,641 | 322 | week ending 2020-08-08 to 2026-10-10 | `9685ddd4b0035694754334e530d5b356bd3503f57b82821f7ee7883c4f5fb94b` |
| `external/jhu_covid19_china_cases.csv` | https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv, filtered to `Country/Region == China` | HTTP 200, 2026-11-04; China subset only | 184,502 B | 34 rows (31 mainland provincial units, Hong Kong SAR, Macao SAR, 'Unknown'); Taiwan, China is not in the extract | cumulative confirmed cases, 1,143 daily columns (1/22/20 to 3/9/23) | 2020-01-22 to 2023-03-09 | `246d97ea96a105063c6381d7c4e1036fbc68432456ac87f9c347419b81a07d57` |
| `external/census_NST-EST2025-ALLDATA.csv` | https://www2.census.gov/programs-surveys/popest/datasets/2020-2025/state/totals/NST-EST2025-ALLDATA.csv (Census Bureau Vintage 2025 state totals) | HTTP 200, 2026-11-05 | 53,555 B | 66 (nation, regions, divisions, 50 states, DC, Puerto Rico) | 97 | estimates 2020 to 2025 | `92188e29cb0a67dcf95afa7d6c47359409782f086478b70ea4128eb70e223ca9` |
| `external/delphi_nhsn/` (51 files `delphi_nhsn_flu_<state>.json.gz` + `MANIFEST.json`) | https://api.delphi.cmu.edu/epidata/covidcast/ (`data_source=nhsn`, signal `confirmed_admissions_flu_ew`, weekly, state, all issues); per-file URL and SHA-256 in `MANIFEST.json` | HTTP 200, 2026-11-05 | 3,048,182 B in total (gz) | 654,504 issue records | time_value, issue, lag, value | issues from epiweek 2024-47 | `MANIFEST.json`: `fab6e1742518f37453bdf96e523359997d2285ab2a8db73d47a217a4378dbd9a`; each data file is checked against the manifest |

The Census file and the Delphi issue history were added for the second revision (`analysis/analysis_plan_amendment2.md`). Population is the mean of `POPESTIMATE2022`–`POPESTIMATE2025` over the 51 jurisdictions with `SUMLEV == 40` (`population()` in `code/rev_r2_lib.py`). The Delphi files are read only by `code/rev_r2_reporting.py` (backfill comparison); records with a null value are skipped; the script fetches a file only if it is absent and leaves `MANIFEST.json` unchanged otherwise.

NHSN fields used:
- `Geographic aggregation` and `Week Ending Date`;
- `Total Influenza Admissions`, `Total COVID-19 Admissions` and `Total RSV Admissions`;
- the percentage of hospitals reporting each pathogen, which gives the reporting-completeness proxy p_hat.

Fifty-one jurisdictions (50 states and DC) are used; territories and national aggregates are dropped. The state-week analysis panel is built in memory from the raw NHSN file by `build_panel()` in `code/revision_data.py` (window rule, $R_A$ and $G$ statistics); no derived NHSN file is read.

The China subset supports the cross-pathogen and cross-regime validation (`code/china_covid_empirical.py`). Its 1,143 daily columns are aggregated to weekly province-weeks; the analysed window is 2020-01-25 to 2023-01-21, and rows enter when their own week-to-week growth factor is at least 1.05 and $C_t \ge 5$. The series ends on 2023-03-09 because the JHU CSSE global feed stopped updating after China discontinued daily case reporting.

## Derived (rewritten by the pipeline; hash reported, not enforced)

| File | Content | Rows × cols | SHA-256 (2026-11-05) |
|---|---|---|---|
| `panel_china_covid19_2020_2023.csv` | Province-week China COVID-19 panel written by `code/china_covid_empirical.py` from the JHU extract (all 34 rows; the analysis keeps the 31 mainland units and uses the SARs only in a sensitivity specification). Columns: `state`, `week_end`, `C`, `neg`, `days`, `C_t`, `C_next`, `neg_next`, `gap_ok`, `neg_any`, `g_own`, `period` | 5,270 × 12 | `1b518b4ad1b723e7313ee9df6da1b1a5ccf31f1934a665d53c988f8e5ed7797f` |

## Legacy files (not read by any script; kept for provenance)

| File | Content | Size / rows | SHA-256 |
|---|---|---|---|
| `external/delphi_fluview_clinical_raw.json` | Delphi `fluview_clinical` state requests, epiweeks 2022-35 to 2026-20, 44 of 51 jurisdictions (none for ak, dc, de, nh, ny, ri, ut) | 1,460,745 B; 7,675 records | `7c14547acead56bf52333fad18449efc13485d67b0db3f635cb573ae7283651a` |
| `panel_nhsn_vintage_2026-10-20.csv.gz` | Earlier state-week influenza panel ($D_t \ge 5$); its builder is not in the repository and the current code rebuilds the panel with `build_panel()` | 10,506 × 21 | `244a025f9833c54bd07a3fbb1d122e19dbe696e89f8e33aa022dfb6fa9d86c7a` |
| `national_flu.csv` | National weekly influenza admissions and growth; not used by the current window rule | 323 × 4 | `91cb2a0d5e223aeb723f480f842b28bd1c1a7dab6d1d87cf0893d87c5d9d635d` |

## Integrity checks performed

- **Checksums.** `code/rev_verify_data.py` (run 2026-11-05, after the second-revision inputs were added): all eight files above and the 51 Delphi files listed in `external/delphi_nhsn/MANIFEST.json` match.
- **NHSN.**
  - Row count, column count and week range were read back after decompression; the values in the table above are those readings.
  - Primary analysis sample ($R_A$, four early-growth seasons): N = 1,748 state-weeks, 51 jurisdictions, 42 calendar weeks, $C_t \ge 5$. Pooled affine fit: $a = 0.01495$ (two-way cluster-robust 95% CI 0.0109–0.0191), $b = 0.560$ (0.430–0.690); see `analysis/results.json` key `revision_2026_11`.
  - The $G$ renewal sample has N = 1,756 with $D_t \ge 5$.
  - No early-window rows come from imputed zeros: 0 of 1,748 rows for $R_A$ and 0 of 1,756 rows for $G$.
  - The earlier statement about a three-season calendar-window baseline (a = 0.0991, N = 1,267) could not be reproduced by any script in the repository and has been removed; the calendar-window specification reproduced by the current code (sensitivity S-a4) gives N = 1,263 and a = 0.0202.
