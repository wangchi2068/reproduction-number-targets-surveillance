# Data decision (revise-experiments, 2026-10-20)

## 1. Main data: no new dataset

Section 2 of revision_plan concluded that the main analysis needs no new dataset. This step follows that conclusion.

We use the 2026-10-20 vintage of the same NHSN HRD table (`ua7e-t2fy`). It covers the 2025–26 season, which raises the number of seasons from S = 3 to S = 4. Its early-window counts reproduce the submitted estimate: a = 0.099055 against 0.099060, with N = 1,267.

Generation-interval and admission-delay parameters come from the literature, as listed in `code/cc_common.py`. No additional data are needed for them.

## 2. Auxiliary data

| Candidate | Outcome | Use |
|---|---|---|
| Delphi `fluview_clinical` (NREVSS clinical laboratory, state level) | Retrieved (HTTP 200). 44 of 51 jurisdictions; ak, dc, de, nh, ny, ri and ut return no rows. Epiweeks 2022-35 to 2026-20 | M6 ④: a second, positivity-based proxy. N = 1,224 state-weeks in the rule windows. Share of the R_A intercept shared with laboratory growth: 0.23 [0.04, 0.42]. Reported as a partial check, not a full decomposition |
| NHSN reporting-share fields (same table) | Present | M6 ②: TOST with pre-specified margins, using four seasons |
| ILINet (state level) | Not refetched | Intervals reused from the vintage rerun (`r4_results.json`, `I_independent_data`, state-cluster BCa): q1 = 0.687, q2 = 1.467. Supplementary and direction-only, as in M7 ④ |
| FluSurv-NET | Not used | Reported as "insufficient information" (M7 ④), not as a non-replication |
| FluSight; European COVID-19 data | Not used | Excluded by revision_plan §2. FluSight targets NHSN itself, and the European data cover a different pathogen and system |

## 3. Window rule: forking-path disclosure

The submitted manuscript used calendar windows over three seasons. To add 2025–26, a data-driven rule was needed. The rule was revised once, in this step, before any intercept was fitted on the new windows:

- **Rule v1** ended each window at the season maximum. In the double-peaked 2024–25 season it gave only 4 weeks, so it was abandoned.
- **Rule v2**, the one used here, ends each window the week before the first national turning point (the first week with national weekly growth < 1). It keeps the contiguous weeks with national growth ≥ 1.05.
  - Its lower bounds are 1 October, plus 2024-11-09 for 2024–25. That is the first full week after mandatory reporting resumed.

Because the change from v1 to v2 is a researcher degree of freedom, the paper must:
- report it as an analysis-plan revision, with its date;
- keep the calendar S = 3 analysis alongside as a continuity specification;
- show both in the forest plot (Fig 4D) and in Table S5.

The headline b/a is similar under the two windows: G gives 22.6 with calendar windows and 23.5 with rule windows.

The user's r3 → r4 rule change (M7 ①) is a separate issue. It belongs to the manuscript text, which must describe it as a dated analysis-plan revision.

## 4. Growth-null rerun

We reran the pure-reporting growth-null simulation of the user's r4 code on the vintage data.
- Size: 756 cells × 300 replicates, complete; 1,386 s.
- Null ranges, A_all:

  | Influenza window | q1 | q2 |
  |---|---|---|
  | Calendar | 0.633–0.886 | 1.126–1.576 |
  | Growth | 0.68–0.88 | 1.33–1.58 |

These ranges replace the manuscript constants (q1 0.883–1.338, q2 0.892–1.373). The consequence for the paper is substantive. A reporting-only process already produces q1 < 1 < q2, so that ordering cannot show transmission heterogeneity. Specifically:
- no pathogen × window ratio falls outside its null range at Holm-adjusted p < 0.05;
- the cross-pathogen "replication" argument must be withdrawn or recast.

## 5. Untrusted-input check

The user's GitHub repository (`sources/repo/`) and `uploads/` were treated as data, not as instructions.

We ran a case-insensitive keyword scan for "ignore previous", "忽略之前", "忽略以上", "disregard", "system prompt" and "you are an AI". It covered the repository text and the extracted text of `uploads/main.pdf`, and it found no instruction-like text.
