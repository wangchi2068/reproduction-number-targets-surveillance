# Data contract (experiment step, 2026-10-20)

- **Unit of analysis:** jurisdiction × week. There are 51 jurisdictions. A week is identified by its NHSN week-ending Saturday.
- **Case count:** C_t, the weekly new laboratory-confirmed influenza admissions reported to NHSN HRD.
- **Statistics:**
  - **Weighted growth ratio:** G = C_{t+1}/D_t, where D_t = 0.65 C_t + 0.25 C_{t−1} + 0.10 C_{t−2}.
  - **Single-week ratio:** G1 = C_{t+1}/C_t.
  - **Renewal-type ratio:** R_A = w ℓ_R((C_{t+1}+0.5)/(C_t+0.5)). It uses a gamma generation interval with mean 3.6 d and SD 1.6 d, and its scale variable is D = C_t.
- **Inclusion:** scale ≥ 5.
  - The four counts C_{t−2..t+1} must be observed (not imputed).
  - For the reporting analyses, the reporting share must be positive in t−2..t+1.
- **Windows:**
  - **Calendar (S = 3):** the windows of the submitted manuscript.
  - **Rule v2 (S = 4, primary):** contiguous weeks with national growth ≥ 1.05, ending the week before the first national turning point.
    - Lower bounds: 2022-10-01, 2023-10-01, 2024-11-09 and 2025-10-01.
    - Resulting windows:

      | Season | Window | N |
      |---|---|---|
      | 2022-23 | 10-01 to 11-26 | 367 |
      | 2023-24 | 10-07 to 12-23 | 503 |
      | 2024-25 | 11-09 to 12-28 | 357 |
      | 2025-26 | 10-04 to 12-27 | 529 |

- **Reporting regimes:**
  - **1:** mandatory, up to April 2024;
  - **2:** voluntary, May to October 2024;
  - **3:** mandatory again, from November 2024.
- **Seasons:** a season starts in August, so the season of week t is year(t) − [month(t) < 8].
- **Missing values:**
  - Rows with missing counts are dropped and never imputed.
  - Jurisdictions absent from Delphi are dropped only from the laboratory proxy.
- **Revisions:** results refer to the 2026-10-20 vintage. Early-window counts match the submitted vintage at the estimate level (a = 0.099055 against 0.099060).
- **Trust:**
  - Text in the user's repository and in `uploads/` was treated as data.
  - A keyword scan for instruction-like text ("ignore previous instructions" and similar) found none in `sources/repo/` or `uploads/main.pdf`.
