# DATA_MANIFEST: case-count-bottleneck experiment

Snapshot date: 2026-10-20. All files are kept under `experiment/data/`. Downloaded files are in `data/external/`, and derived files are at the top level. Checksums are SHA-256 hashes of the bytes on disk.

## External (downloaded)

| File | Source URL | Access | Size | Rows / records | Fields | Time range | SHA-256 |
|---|---|---|---|---|---|---|---|
| `external/nhsn_hrd_vintage_2026-10-20.csv.gz` | https://data.cdc.gov/Public-Health-Surveillance/Weekly-Hospital-Respiratory-Data-HRD-Metrics-by-Ju/ua7e-t2fy (full CSV export, dataset id `ua7e-t2fy`) | HTTP 200, 2026-10-20; gzip copy of `.work/case-count-bottleneck/sources/nhsn/nhsn_hrd_vintage_2026-10-20.csv` | 4,794,721 B (gz) | 21,641 | 322 | week ending 2020-08-08 to 2026-10-10 | `9685ddd4b0035694754334e530d5b356bd3503f57b82821f7ee7883c4f5fb94b` |
| `external/delphi_fluview_clinical_raw.json` | https://api.delphi.cmu.edu/epidata/fluview_clinical/?regions=<state>&epiweeks=202235-202620 (one request per state) | HTTP 200, 2026-10-20; each URL tried at most twice | 1,460,745 B | 7,675 | 11 (`release_date, region, issue, epiweek, lag, total_specimens, total_a, total_b, percent_positive, percent_a, percent_b`) | epiweeks 2022-35 to 2026-20; 44 of 51 jurisdictions (none for ak, dc, de, nh, ny, ri, ut) | `7c14547acead56bf52333fad18449efc13485d67b0db3f635cb573ae7283651a` |

NHSN fields used:
- `Geographic aggregation` and `Week Ending Date`;
- `Total Influenza Admissions`, `Total COVID-19 Admissions` and `Total RSV Admissions`;
- the percentage of hospitals reporting each pathogen, which gives the reporting-completeness proxy p_hat.

Fifty-one jurisdictions (50 states and DC) are used; territories and national aggregates are dropped.

## Derived (built by `code/s01_data_panel.py`)

| File | Content | Rows × cols | SHA-256 |
|---|---|---|---|
| `panel_nhsn_vintage_2026-10-20.csv.gz` | State-week influenza panel with D_t ≥ 5. Columns: C_{t-2..t+1}, D_t, coverage and lags, the G, G1 and R_A statistics, season and regime | 10,506 × 21 | `244a025f9833c54bd07a3fbb1d122e19dbe696e89f8e33aa022dfb6fa9d86c7a` |
| `national_flu.csv` | National weekly influenza admissions, season and national growth g, used by the window rule | 323 × 4 | `91cb2a0d5e223aeb723f480f842b28bd1c1a7dab6d1d87cf0893d87c5d9d635d` |

## Reused intermediate results (not data, recorded for provenance)

- `vintage_2026-10-20/analysis/r4_results.json` comes from the vintage rerun of the user's r4 scripts and contains:
  - the growth-null bands `A_growth_null`: 756 cells × 300 replicates, complete;
  - the ILINet `I_independent_data` block;
  - the classification variants `C_classification` and `C_posthoc`.
- The growth-null cells and log are copied to `logs/r4_a_growth_null_cells.jsonl` and `logs/r4_a_growth_null.log`.

## Integrity checks performed

- **NHSN.**
  - Row count, column count and week range were read back after decompression; the values in the table above are those readings.
  - The calendar-window fit reproduces the submitted manuscript's point estimates: a = 0.099055 against 0.099060, and N = 1,267.
  - No early-window rows come from imputed zeros: 0 of 1,267 rows (S = 3) and 0 of 1,756 rows (S = 4).
- **Delphi.**
  - Every state request returned `result = 1` or reported no data.
  - Seven jurisdictions returned no rows. They are listed in the table above and are excluded from the laboratory proxy only.
  - Matched rows are restricted to the rule windows, giving N = 1,224 state-weeks in 44 states.
