# Referee report

**Manuscript:** "What more cases buy: inferential targets for reproduction-number estimates and what cross-sectional dispersion diagnostics measure"

**Submission reviewed:** https://github.com/wangchi2068/reproduction-number-targets-surveillance, branch `main`, commit `2f5c4074a11cd5dbc8744ce4bf0b706c19171156` (compiled `paper/paper.pdf`, 38 pages, dated 4 November 2026; analysis code and data under `experiment/`)

**Date of review:** 5 November 2026

**Recommendation:** Major revision. **Score:** 4/10

Locations below use the section, table and proposition numbers of the compiled PDF, followed by the source file and line in the repository (each paragraph of the LaTeX source is one line), for example Section 4.4 (`paper/sections/diagnostics.tex:134`). Where I say a result was reproduced or not, I re-ran the relevant computation from the shipped data and code; the checks are summarised in Section 4.3.

## 1. Summary

The manuscript asks what an additional observed case buys in the precision of a reproduction-number estimate, and what the cross-sectional dispersion of surveillance counts can tell us about transmission heterogeneity.

**Theory (Sections 2–3).** In a one-generation branching model with negative binomial offspring (dispersion k) and binomial reporting (probability ρ), the authors distinguish three inferential targets: the current value R_t, the realised mean over the observed periods, and the process mean. They derive exact mean squared errors for each. Counting error falls as A/m in the expected count m for all three targets, while the process mean also retains the between-period variance v_R, which more cases cannot remove. The count at which counting error equals v_R, m× = A/v_R, is identified as the half-weight point of the corresponding shrinkage estimator. Using the previous observed count as the denominator inflates the counting error about 1.5-fold at large counts and 2.4-fold at ten expected cases under illustrative parameters (Table 3). Proposition 2 shows that, when the reporting probability fluctuates between periods, single-period second moments cannot separate reporting volatility from transmission variation.

**Empirical diagnostics (Sections 4–5).** The authors fit an affine model, E(L | C) = a + b/C, for the relative squared deviation L of an apparent reproduction-number statistic (a Wallinga–Lipsitch transform of the ratio of consecutive weekly counts) on inverse count. The data are 1,748 state-weeks of US influenza hospital admissions from NHSN, in the early-growth windows of four seasons (2022–2025). A known-truth simulation audit runs twelve generators (S0–S9, including variants) on the real design. It reports that the intercept a recovers transmission variance under fixed reporting, the transmission-plus-reporting total under reporting shocks, and less than transmission variance when heterogeneity scales with count. The observed intercept is a = 0.0150 (95% CI 0.0109–0.0191, two-way cluster-robust sandwich), the counting slope is 1.26 times an independent-Poisson reference, and the dispersion-to-counting crossover is 37 weekly admissions.

**Headline empirical claims (abstract).**
- Across a (ρ, k) "parameter manifold", the implied transmission relative standard deviation is "tightly bound" at 11.6–12.3%, within a global envelope of 5.0–20.4% across mechanisms.
- None of the audited mechanisms reproduces every tested summary.
- Non-parametric calibrations of a season test reject equal season intercepts (p = 0.013 and 0.045), "corroborating seasonal transmission divergence".

**Claimed contribution.** An operational estimand linking surveillance dispersion to transmission and reporting, a formalisation of m× as an actionable shrinkage balance point, and a known-truth audit of the diagnostics on the real sampling design (Introduction, `paper/sections/introduction.tex:10–12`).

## 2. Strengths

1. **A useful conceptual clarification.** Separating the current value, the realised mean and the process mean, and showing which error component each target retains (Table 2), addresses a real ambiguity in how reproduction-number precision is reported. The practical consequence — that a case-count rule controls counting error for one target but not the between-period variance relevant to the others — is stated clearly (`paper/sections/introduction.tex:6`; Section 6.5 checklist).
2. **Correct and verifiable theory.** The MSE results and the observed-denominator correction are internally consistent. An independent recomputation (closed form plus a 400,000-draw Monte Carlo) gave A = 2.1232 and m× = 53.08 for the stated parameters, and reproduced the observed-to-expected denominator error ratios (1.59, 1.50 and 1.49 at m = 50, 200 and 1,000 against the paper's 1.588, 1.508 and 1.489; Table 3).
3. **An honest identification result.** Proposition 2 is a clean statement of what single-period count moments cannot identify, and the authors carry it into the definition of the intercept (`paper/sections/diagnostics.tex:37`) and into the practice checklist (`paper/sections/discussion.tex:44`).
4. **Known-truth audit on the real design.** Running the generators on the actual state-week structure, with the same filters and windows, is good simulation practice. Reporting Monte Carlo standard errors, coverage against both the estimator mean and the true total, and pass/fail decisions near the limits (`paper/sections/diagnostics.tex:132`) is more transparent than is usual.
5. **Candour about limitations.** The manuscript states that the decision rule was fixed after the generators were first compared with the data (`paper/sections/diagnostics.tex:68`), that the generators vary less than the data (`paper/sections/diagnostics.tex:144`), and that no mechanism reproduces every summary (Section 4.4).
6. **Open and largely reproducible materials.** A dated data vintage with checksums, analysis code and result files are provided. The primary sample (1,748 state-weeks, 51 jurisdictions, 42 weeks), the affine fit (a = 0.014955, b = 0.5600, b/a = 37.44), the season intercepts, every interval in the abstract and the shipped season calibrations all reproduced exactly from the shipped files, even under newer versions of numpy, scipy and pandas than those pinned.

## 3. Weaknesses

1. **The season finding is attributed to transmission, contrary to the paper's own Proposition 2** (Comment M1). The abstract and Discussion say the season test "corroborates" or "confirms" seasonal transmission divergence, although the intercept contains reporting volatility that the paper proves cannot be separated, and the windows overlap documented changes in reporting.
2. **The season test's significance depends on the choice of reference distribution** (Comment M3). The two non-parametric calibrations in the abstract do not respect the calendar-week structure of the season contrast; a calendar-week permutation gives p = 0.101. Uncalibrated tests that use the bootstrap covariance give p = 0.064–0.220, and the paper's own power analysis puts the observed effect below the minimum detectable ratio.
3. **The "tightly bound 11.6–12.3%" transmission SD is close to a restatement of √a** (Comment M2). It is computed under the S1 mechanism that the paper's own curve-shape tests weigh against, with v_R set equal to the observed intercept, over a 20-point grid described as a "manifold".
4. **The audit generators do not reproduce key features of the data** (Comment M4), so their operating characteristics (recovery, coverage, test size) may not transfer to the observed panel.
5. **The mapping from the branching model to weekly hospital admissions is loose** (Comment M5). A week spans about two influenza generations, and the effective reporting probability for laboratory-confirmed admissions is plausibly one to two orders of magnitude below the grid's smallest ρ = 0.1.
6. **The novelty claim is broader than the contribution** (Comment M6). Several results are classical, and closely related work on realised reproduction numbers, window choice and random-effects over-dispersion is not discussed.
7. **The availability statement does not hold for the submitted version** (Comment M7). The cited commit lacks scripts behind abstract-level results, the documented pipeline does not run them, and the data manifest and authorship information are incomplete.

## 4. Detailed comments

### 4.1 Major comments

#### M1. The season conclusion contradicts Proposition 2

**Location.** Abstract (`paper/sections/abstract.tex:2`, "corroborating seasonal transmission divergence"); Discussion (`paper/sections/discussion.tex:6`, "confirming seasonal transmission divergence across surveillance cycles"; `paper/sections/discussion.tex:35`, "independently corroborate the seasonal divergence"); Proposition 2 (`paper/sections/reporting.tex:9`); definition of a (`paper/sections/diagnostics.tex:37`).

**Problem.** By the paper's own account, a estimates the relative transmission variance only when reporting probabilities are constant between weeks; otherwise it estimates transmission plus the variance induced by reporting shocks, and "single-period cross-sectional moments cannot separate the two components" (`paper/sections/diagnostics.tex:37`). The season contrast is driven by 2023–24, whose intercept is 0.0071 against 0.0154–0.0180 in the other seasons (Table S13; `paper/sections/application.tex:36`). A lower intercept in one season is equally consistent with steadier reporting in that season. Reporting regimes did differ across the windows:
- Hospital reporting was voluntary from May to October 2024; the 2024–25 window starts on 9 November 2024, after the return to mandatory reporting (`paper/sections/diagnostics.tex:7–9`).
- The 2025–26 window (4 October – 27 December 2025) overlaps a lapse in federal funding during which CDC's respiratory-illness data channel was not updated from 3 October to 13 November 2025, and some data were omitted because of processing delays [1]. The manuscript does not mention this.

**Why it matters.** This is the only substantive empirical finding stated in the abstract, and it attributes to transmission a quantity that the paper proves is not identified as transmission.

**Requested revision.**
- Describe the result as a season difference in count-independent relative dispersion, and remove "corroborating", "confirming" and "transmission divergence" from the abstract, Section 5 and Section 6.
- Discuss the reporting regime of each season window, including the 2025 publication pause, and state whether the 20 October 2026 vintage contains back-filled or revised values for the affected weeks; add a sensitivity analysis that excludes them.
- If a transmission interpretation is wanted, bring external information that bounds reporting volatility by season (for example hospital reporting participation by week, which the authors already use as an adjustment in Section 5.1, or a second surveillance stream).

#### M2. The implied transmission SD of 11.6–12.3% overstates what the data identify

**Location.** Abstract; Section 4.4 (`paper/sections/diagnostics.tex:134`); Table 6 (`paper/tables/tab_ma_manifold.tex`); Discussion (`paper/sections/discussion.tex:6`, "tightly bound"; `paper/sections/discussion.tex:35`).

**Problem.**
- The range is computed under S1 (count-independent heterogeneity, fixed reporting). The abstract's last sentence reports that in-sample and held-out curve shapes "weighed against heterogeneity that does not depend on count", i.e. against S1.
- Under S1 the generator's v_R is set to the observed intercept (Table 6 note: v_R = 0.0150), and S1 is exactly the case in which the theory says a estimates v_R. The recovery ratio a/v_R is therefore close to 1 by construction (0.996–1.102 over the "compatible" cells), and the implied SD is essentially √(0.0150/recovery) ≈ 12%. The narrow 11.6–12.3% spread reflects the variation of the recovery ratio across k and ρ, not information in the data about transmission. What Table 6 does show — that intercept recovery under S1 is robust to k and ρ — is useful, but it is a property of the estimator, not a bound on transmission heterogeneity.
- The "entire (ρ, k) parameter manifold" is a 4 × 5 grid (ρ ∈ {0.1, 0.25, 0.5, 1}, k ∈ {0.1, 0.3, 1, 3, ∞}).
- "Compatible" means that a cell's 95% Monte Carlo range overlaps the observed bootstrap CI [0.81, 1.73] for b/b_P. This is looser than the decision rule declared in Section 4.3 (observed point value inside the simulated 2.5–97.5% range). Under the declared rule only 3 of the 14 "compatible" cells contain the observed 1.257: (ρ, k) = (0.1, 0.1), (0.25, 0.3) and (0.5, 0.3).
- Serial-interval uncertainty is not propagated. With a mean of 3.2 days the intercept falls to 0.0111–0.0119 (Table S15), corresponding to an implied SD of about 10.5–10.9% instead of 12.2%. The China COVID-19 check shows the same sensitivity in stronger form (a rises from 0.133 to 0.298 under a different interval; `paper/sections/appendix.tex:146`).

**Why it matters.** Readers will quote "about 12%" as an estimate of between-state variation in influenza transmission. It is a model-conditional transform of one number, under a mechanism the paper's own tests disfavour.

**Requested revision.**
- Lead with the widest defensible range (the mechanism envelope, with serial-interval and centring sensitivity), and present the S1 value as a conditional illustration.
- Call the grid a grid, apply the declared decision rule when defining compatibility, and report the cells that pass it.
- Remove evaluative wording such as "tightly bound" and "remarkably stable" (`paper/sections/diagnostics.tex:134`).

#### M3. Inference for the season comparison is fragile and selectively summarised

**Location.** Abstract; Section 4.4 (`paper/sections/diagnostics.tex:144`); Section 5.3 (`paper/sections/application.tex:36`); Table 11 (`paper/tables/tab_ma_calibration.tex`); Table S13 (`paper/tables/tab_ma_seasons.tex:27–31`); `experiment/code/ma_wild_bootstrap.py`.

**Problem.** The season regressor is constant across states within a calendar week, so the contrast is identified only by variation between weeks and, ultimately, between four seasons. The manuscript's own results already show that the conclusion depends on the reference distribution:

| Reference distribution for the 3-df test | p-value | Source |
|---|---|---|
| χ² (3 df), two-way sandwich covariance (W_CGM = 17.8) | 0.0005 | Table 11 |
| χ² (3 df), bootstrap covariance, log intercepts (W = 5.52) | 0.138 | Table S13 |
| χ² (3 df), bootstrap covariance, raw intercepts | 0.064 | Table S13 |
| χ² (3 df), separate season fits | 0.220 | Table S13 |
| Generator-simulated nulls for W_CGM (S1, S6, S8, S9) | 0.005–0.045 | Table 11 |
| Wild cluster bootstrap, Webb weights, state clusters | 0.013 | Table 11; reproduced |
| Within-state permutation of season labels | 0.045 | Table 11; reproduced |
| Wild cluster bootstrap, Webb weights, calendar-week clusters (reviewer check, B = 4,999) | 0.014 | Section 4.3, check D3 |
| Permutation of whole calendar-week cross-sections between seasons (reviewer check, B = 4,999) | 0.101 | Section 4.3, check D4 |

Specific concerns:
- The asymptotic χ² (3 df) reference for W_CGM is badly oversized (size 0.195 and 0.217 at nominal 0.05 under S6 and S8; Table 11), so the 0.0005 is not informative.
- The generator-calibrated p-values are valid only if the generators reproduce the dependence in the data, which they do not (M4; the observed sandwich SE of a is 1.7–2.9 times the generator Monte Carlo SD, `paper/sections/diagnostics.tex:144`). Under S6 the calibrated p is 0.045, at the threshold.
- The shipped wild bootstrap draws one Webb weight per state (`experiment/code/ma_wild_bootstrap.py:45–59`). Clustering by state treats common weekly shocks — which the week-centring result shows are substantial (Section 5.1) — as independent across clusters; guidance is to cluster at the level at which the regressor of interest and the errors are correlated, here calendar time [2]. With effectively four "treated" groups, wild cluster bootstrap tests can be badly mis-sized [3].
- The permutation reshuffles season labels among a state's own rows (`experiment/code/ma_wild_bootstrap.py:117`). This breaks the alignment of calendar weeks across states and the serial order within a state; it does not "preserve all intra-state dependent correlation" as stated at `paper/sections/diagnostics.tex:144`, and it assumes exchangeability of seasons within state, which common national shocks violate. A permutation that keeps each weekly cross-section intact does not reject (p = 0.101). I do not claim this alternative is definitive either; the point is that the result is not robust to defensible choices.
- The paper's own power analysis gives a minimum detectable highest/lowest intercept ratio of 3.71 against an observed 2.54 (Table S13), and the between-season variance can be estimated only to a relative SE of at least 0.82 with four seasons (Appendix B.3).
- The six-point weights are those of Webb [4], who is not cited. Bootstrap theory for two-dimensional dependence [5] and multiway wild bootstraps [6] are relevant and not discussed.

**Why it matters.** The abstract reports only the two favourable non-parametric calibrations, and Section 5.3 states that "under all calibrated specifications the test rejects". A reader cannot tell from the abstract that the most dependence-respecting reference does not reject.

**Requested revision.**
- Report all reference distributions in one table, as above, and justify the exchangeability unit of each.
- Add a calibration that respects calendar time (week or season-week clusters, or a two-way wild bootstrap by state and week).
- Describe the season comparison as exploratory and descriptive, consistent with the power analysis.
- Consider a model with explicit week effects, or a hierarchical model with season and week random effects, which would make the source of the season contrast explicit.

#### M4. The audit generators do not reproduce key features of the data

**Location.** Section 4.3 (`paper/sections/diagnostics.tex:68`); Section 4.4 (`paper/sections/diagnostics.tex:132`, `:134`, `:140`, `:144`); Tables 9–11; Table S2 (`paper/tables/tab_ma_share_tests.tex`).

**Problem.**
- The pooled/within slope ratio fails in all 30 cells of the share grid; the largest simulated upper limit is 1.25 against the observed 1.48 (`paper/sections/diagnostics.tex:132`).
- The observed sandwich SE of a is 1.7–2.9 times the generator Monte Carlo SD (`paper/sections/diagnostics.tex:144`).
- The generators' 90th percentiles of retained counts are 263–324 against 389 observed (`paper/sections/diagnostics.tex:68`).
- Under S8 the sandwich interval covers the true total in only 0.42 of panels (`paper/sections/diagnostics.tex:140`).
- The decision rule and the week-centring test were fixed after the generators were first compared with the data, and S9 was added after the first eleven generators had been compared with the data (`paper/sections/diagnostics.tex:68`, `:70`).
- At `paper/sections/diagnostics.tex:134` the pooled/within divergence is then described as "statistically compatible with sample variability" because its bootstrap CI is 0.988–3.662. Either this statistic is informative, in which case every generator fails it, or it is too noisy to discriminate, in which case it should not be one of the six tests that eliminate every cell. The manuscript currently uses it both ways.

**Why it matters.** The statements that the intercept "recovers" v_R or the total, and that the sandwich interval is calibrated, are established only under generators that are demonstrably less variable than the data.

**Requested revision.**
- Add at least one generator that reproduces between-state differences in baseline dispersion (the authors attribute the between-state slope to "smaller jurisdictions exhibiting higher baseline dispersion", `paper/sections/diagnostics.tex:134`) and the upper count tail, and report recovery, coverage and season-test size under it.
- Resolve the pooled/within inconsistency, for example by taking the sampling uncertainty of the observed statistic into account in the decision rule.
- Label the post hoc rules and S9 as exploratory throughout, not only in Section 4.3.

#### M5. The branching model is mapped loosely onto weekly hospital admissions

**Location.** Section 2 (`paper/sections/model.tex`); Section 4.1 (`paper/sections/diagnostics.tex:22–24`); Table 6; Section 6.3 (`paper/sections/discussion.tex:28`).

**Problem.**
- The theory is a one-generation branching experiment. Influenza generation and serial intervals are about 3–4 days [7, 8], so a week spans roughly two generations and the previous week's count is not the parent generation of the next week's count. The ratio statistic is mapped to a reproduction number through the growth-rate relation [9], a non-linear step whose output depends on the interval distribution (M2).
- The theory's ρ is the probability that an infection is counted. For laboratory-confirmed influenza admissions this is the product of the infection-to-symptomatic, symptomatic-to-hospitalised and hospitalised-to-confirmed probabilities. CDC estimates for 2010–2025 imply on the order of 1–1.5% of symptomatic illnesses are hospitalised [10]; laboratory-confirmed admissions are a fraction of true hospitalisations (under-detection multipliers 2.1–5.2 by age [11]); and most infections are asymptomatic [12]. My rough combination of these figures gives an effective ρ of about 0.001–0.01, one to two orders of magnitude below the smallest grid value of 0.1. This is a back-of-envelope calculation, not a published estimate, but the authors should either justify the grid range or extend it.
- At such ρ, binomial thinning makes the reported counts close to Poisson whatever the offspring law, so k is not identifiable from these counts; the manuscript already finds all five offspring laws slope-compatible at ρ = 0.1 (Table 6). Work on estimating k typically relies on cluster-size data rather than incidence [13].
- The weights of the denominator D_t = 0.65C_t + 0.25C_{t−1} + 0.10C_{t−2} (`paper/sections/diagnostics.tex:13`) are not derived from the interval distribution.

**Why it matters.** The empirical part borrows the theory's vocabulary (ρ, k, v_R) for quantities that are not the parameters of the data-generating process for weekly admissions, and the implied transmission variances inherit this mismatch.

**Requested revision.**
- Either derive the dispersion reference values for an aggregated renewal process with weekly binning and a hospitalisation-scale ρ, or present the empirical part explicitly as describing a count statistic and drop transmission-variance language.
- Extend the (ρ, k) grid to ρ ≤ 0.01, or justify its lower limit with a source.
- Justify the D_t weights from the interval distribution, or show sensitivity to them.

#### M6. Novelty relative to existing work is overstated

**Location.** Introduction (`paper/sections/introduction.tex:6`, `:10`, `:12`); Discussion (`paper/sections/discussion.tex:4`, `:13`).

**Problem.**
- The MSE decompositions follow from the law of total variance, and m× is the classical shrinkage half-weight point [14]; the authors acknowledge both (`paper/sections/introduction.tex:10`). The affine variance form a + b/C is the additive random-effects over-dispersion model σ²/n + τ² used in institutional comparisons and funnel plots [15, 16], which is not cited.
- The "12 cases" rule attributed to Cori et al. (`paper/sections/introduction.tex:6`) is in their Supplementary Web Appendix 3, where it follows from a target posterior coefficient of variation of 0.3, counts cumulative cases, and is one of three criteria for starting estimation [17]. It is a heuristic for posterior precision, not a claim about target-specific MSE; describing the criterion as "incomplete" for other targets is fair but should be framed accordingly. The same supplement (Web Appendix 8) already shows that constant binomial under-reporting leaves mean R estimates essentially unchanged, and Gostic et al. make the same point while showing that time-varying reporting does bias R_t [18]. These should be credited where the manuscript states that constant ascertainment cancels.
- Closely related work is not discussed: the case reproduction number of Wallinga and Teunis, an estimator of a realised, case-level quantity and the clearest prior example of the realised/process distinction [19]; data-adaptive window selection for renewal estimators, which addresses the same bias–variance trade-off over how many cases a window needs [20]; and renewal models with an explicit observation layer pooled across regions [21].

**Why it matters.** Readers need to know what is new beyond a clear synthesis.

**Requested revision.** State the contribution as a synthesis together with a diagnostic and audit framework; engage directly with the work above; soften the framing of the Cori et al. criterion and cite its supplement precisely.

#### M7. The reproducibility and availability statements do not hold for the submitted version

**Location.** Data and code availability (`paper/sections/availability.tex:4`); `experiment/code/run_all.sh`; `experiment/data/DATA_MANIFEST.md`; `experiment/code/check_manuscript.py`; author block (`paper/main.tex:144`); CRediT statement (`paper/sections/declarations.tex:5`).

**Problem.**
- The availability statement cites commit `6bbb8458a5ba3a9ba644d755320df55dbe06746b` (4 November 2026). That commit does not contain `ma_wild_bootstrap.py` or `ma_manifold_envelope.py`, which produce the season calibrations and the manifold range quoted in the abstract.
- The statement says that running `run_all.sh` "regenerates every reported estimate, table and figure". It does not call either script; `ma_gaps.py` only reads their stored JSON output.
- The data manifest says the derived panel is built by `code/s01_data_panel.py`, which is not in the repository.
- The shipped JHU CSSE China extract does not match the manifest (184,502 bytes against the stated 184,537), and the derived China panel has no recorded checksum.
- `check_manuscript.py` fails as shipped because it expects `paper_en/` and `paper_zh/` folders that do not exist.
- The author block reads "[Author names and affiliations to be added before submission]", and the CRediT statement lists roles without names.
- The repository contains review and response files from an earlier round (`review/`), including a self-assigned assessment. These have no bearing on the scientific evaluation but should not be part of a submission package.

**Positive note.** Everything I re-ran from the shipped data reproduced exactly or within Monte Carlo error (Section 4.3), which is to the authors' credit.

**Requested revision.** Cite the final commit, preferably an archived release with a DOI; include every analysis script in `run_all.sh`; fix the manifest and checksums; make the manuscript checker run; complete the author and CRediT information; remove internal review material from the submission.

### 4.2 Minor comments

1. **Envelope rounding** (abstract; `paper/sections/discussion.tex:6`; `paper/sections/diagnostics.tex:134`). The text gives 5.0–20.4%, whereas Table S2, `paper/sections/diagnostics.tex:132` and the stored result (0.05124–0.20451) give 5.1–20.5%. Harmonise.
2. **Table 11 rows are misaligned** (`paper/tables/tab_ma_calibration.tex`). The values are shifted down by one row: the section header "Non-parametric calibrations on real data" carries the wild-bootstrap values (9.8, 0.013); the "Wild cluster bootstrap" row carries the permutation values (17.1, 0.045); and the "Within-state season permutation" row carries the χ² (3 df) reference values (7.8, 0.0005, 0.138). As printed, the table misattributes p-values to methods.
3. **Centring and the estimand** (Section 5.1, `paper/sections/application.tex:9`; Table S8). Week centring removes 34.7% of a (6.0–51.9). In my recomputation, centring on the cross-state mean of the sampled rows in each week removed 42.6%. A large part of the primary intercept is therefore common weekly (national) variation, not between-unit heterogeneity. State which estimand is primary and report the within-week intercept next to the pooled one in the abstract.
4. **Normalisation of the reference slope** (`paper/sections/diagnostics.tex:40`). b_P is built from the row-level elasticity κ(g), whereas L is normalised by the pooled mean X̄². Please state the approximation this implies (row values of R_A close to X̄) and quantify its effect on b/b_P.
5. **Window selection** (`paper/sections/diagnostics.tex:9`; `paper/sections/discussion.tex:35`). The windows are selected on the national growth trajectory, and the 2024–25 start is set to follow the return to mandatory reporting. Report results for fixed calendar windows alongside, or give the rule in a form that could have been pre-specified.
6. **Count-floor dependence of the crossover** (Table 12). With count floors of 10 and 50 the crossover is 49.7 and 53.1 (bootstrap 25.0–96.6 and −8.4–180.1), against 37 in the primary analysis. Mention this next to the headline value.
7. **Thinning check** (Section 4.2, `paper/sections/diagnostics.tex:61`). The paper notes that c/b_P is a check of the delta approximation rather than of mechanism. Given its low diagnostic value, consider moving it to the supplement.
8. **China COVID-19 appendix** (Appendix B.5, `paper/sections/appendix.tex:140–148`). The authors list imports, batch reporting and changing testing regimes, and the intercept moves from 0.054 to 0.298 across reasonable specifications. The appendix adds little to the influenza analysis and could be shortened or removed.
9. **Abstract length** (`paper/sections/abstract.tex`). At about 340 words (counting markup) and dense with numbers, the abstract mixes theory, audit and application. Shorten it to the main messages and check the journal's limit.
10. **Register** (`paper/sections/diagnostics.tex:134`, `:144`). Sentences beginning "To resolve the conditionality…" and "To resolve reliance on synthetic generators…" read as responses to a previous review; rewrite as neutral methods text, and remove evaluative adverbs.
11. **Pre-specification of the serial interval** (`paper/sections/discussion.tex:28`). The statement that the interval was "fixed before the analyses reported here" cannot be checked from the materials; point to a time-stamped protocol or remove it. Note also that the Cowling et al. estimate comes from a Weibull fit to 14 household pairs [7] and is used here as a gamma with the same moments; say so.
12. **Owen (2007) does not explain the observed bootstrap inflation** (`paper/sections/diagnostics.tex:138`). The pigeonhole bootstrap variance is at most about three times the true variance [22], an SD ratio of at most about 1.73. The reported SD ratios of 1.6–2.5 exceed that at the upper end, so "in line with" Owen is only partly true; explain the excess (for example the non-linear ratio estimator or the unbalanced design). Owen and Eckles extend the method to higher-order arrays [23], which may be relevant to the state × week × season structure.
13. **Bibliography** (`paper/bibliography.bib`). None of the 58 entries has a DOI. The NHSN dataset entry lacks its identifier (ua7e-t2fy) and URL; the CMS memo entry lacks its date (22 October 2024) and URL. The key `frutos2025evaluation` carries the year 2026 (published online December 2025, issue February 2026 [24]).
14. **Notation.** "Wild cluster bootstrap" is capitalised inconsistently in the text and tables.

### 4.3 Reproducibility and internal-consistency checks

All checks were run on the shipped files at the commit reviewed. The full simulation pipeline was not re-run; simulation-based tables were checked against their stored results.

| ID | Check | Paper | Re-computed | Result |
|---|---|---|---|---|
| A3–A4 | Primary sample rebuilt from the shipped panel and national series | 1,748 state-weeks, 51 jurisdictions, 42 weeks | 1,748 / 51 / 42; season rows 365 / 499 / 356 / 528 | Match |
| B1 | Gamma quasi-likelihood affine fit | a = 0.0150, b = 0.560, b/a = 37.4 | a = 0.014955, b = 0.5600, b/a = 37.44 | Match |
| B2 | Season intercepts with common slope | 0.0159 / 0.0071 / 0.0154 / 0.0180; slope 0.607 | 0.01585 / 0.00709 / 0.01537 / 0.01800; slope 0.6072 | Match |
| B3 | Abstract intervals traced to stored results | sandwich 0.0109–0.0191; bootstrap 0.0101–0.0206; b/b_P 1.26 (0.81–1.73) | 0.010852–0.019058; 0.0101–0.0206; 1.2568 (0.8117–1.7341) | Match |
| B4 | Global envelope | 5.0–20.4% (text) | 5.124–20.451% (stored) | Rounding inconsistency |
| B5 | Manifold "compatible" cells containing the observed b/b_P = 1.257 | 14 of 20 compatible | 3 of 14 contain 1.257 | Interpretation issue (M2) |
| B6 | Effect of centring on a | week centring −34.7% | season −3.3%; week −42.6% (cross-state mean of sampled rows) | Consistent in direction |
| C1–C3 | Observed-denominator Monte Carlo (400,000 draws) | ratios 2.423 / 1.588 / 1.508 / 1.489 at m = 10 / 50 / 200 / 1,000 | 2.381 / 1.594 / 1.502 / 1.485 | Match within Monte Carlo error |
| D1–D2 | Shipped season calibrations | W_CGM = 17.8; p = 0.013 and 0.045 | 17.797; p = 0.0125 and 0.0455 | Match |
| D3 | Wild cluster bootstrap with calendar-week clusters, B = 4,999 | — | p = 0.014 (MCSE 0.002) | Rejects at 5% |
| D4 | Calendar-week permutation, B = 4,999 | — | p = 0.101 (MCSE 0.004) | Does not reject |
| D5 | Rendering of Table 11 | one row per method | values shifted by one row | Table defect |
| E1, E5 | Pipeline and cited commit contain all scripts behind the abstract | yes | two scripts missing from both | Not as stated |
| E2 | Manuscript checker | runs | fails as shipped (missing folders) | Not as stated |
| A1–A2 | Data manifest | checksums and builder listed | China extract mismatch; panel builder absent | Not as stated |

### 4.4 Questions to the authors

1. Did reporting participation or completeness differ between the four season windows, and does the 2023–24 season differ in any reporting metric, such as the hospital participation used in the Section 5.1 adjustment?
2. Do any of the audit generators include season-level or week-level national shocks of the size implied by the week-centring result? If S9 is the only one, how does the season test behave under S9 with larger national shares?
3. What is the effective reporting probability ρ you consider plausible for laboratory-confirmed NHSN influenza admissions, and how does Table 6 change for ρ ≤ 0.01?
4. Was the S1 generator's v_R tied to the observed intercept in every cell of Table 6, and if so, what does the table add beyond showing that a is unbiased under S1?
5. Which commit, or which archived release, should readers use to reproduce every number in the submitted manuscript?

## 5. Overall recommendation

**Recommendation: major revision. Score: 4/10.**

The theoretical part is correct, clearly written and useful, although much of it is a careful synthesis of known results (W6). The three-target framing, the observed-denominator correction and Proposition 2 deserve publication in some form, and the known-truth audit on the real design and the largely reproducible code are real strengths.

The empirical messages promoted in the abstract are not supported as stated. The "tightly bound" 12% transmission SD is close to a restatement of the observed intercept under a mechanism the paper's own tests disfavour (M2). The "seasonal transmission divergence" attributes to transmission a quantity that Proposition 2 shows is not identified (M1), and its significance does not survive a reference distribution that respects calendar time (M3). The validation generators fail key features of the data (M4), the model-to-data mapping for hospital admissions is loose (M5), and the availability statement does not hold for the submitted version (M7).

These problems can be addressed by reframing the empirical claims, reporting all season-test calibrations, adding one more realistic generator and tidying the reproducibility package; none requires new data. That is why I recommend major revision rather than rejection. If M1–M4 were resolved and the abstract scaled back to what the audit supports, I would expect to rate the work around 6–7.

## References

**[1]** Centers for Disease Control and Prevention. Respiratory Illnesses Data Channel: Frequently Asked Questions (release notes). CDC website, 2025. https://www.cdc.gov/respiratory-viruses/data/faqs.html

**[2]** MacKinnon JG, Nielsen MO, Webb MD. Cluster-robust inference: A guide to empirical practice. Journal of Econometrics 232(2):272-299, 2023. https://doi.org/10.1016/j.jeconom.2022.04.001

**[3]** MacKinnon JG, Webb MD. The wild bootstrap for few (treated) clusters. Econometrics Journal 21(2):114-135, 2018. https://doi.org/10.1111/ectj.12107

**[4]** Webb MD. Reworking wild bootstrap-based inference for clustered errors. Canadian Journal of Economics 56(3):839-858, 2023. https://doi.org/10.1111/caje.12661

**[5]** Menzel K. Bootstrap With Cluster-Dependence in Two or More Dimensions. Econometrica 89(5):2143-2188, 2021. https://doi.org/10.3982/ECTA15383

**[6]** MacKinnon JG, Nielsen MO, Webb MD. Wild Bootstrap and Asymptotic Inference With Multiway Clustering. Journal of Business & Economic Statistics 39(2):505-519, 2021. https://doi.org/10.1080/07350015.2019.1677473

**[7]** Cowling BJ, Fang VJ, Riley S, Peiris JSM, Leung GM. Estimation of the Serial Interval of Influenza. Epidemiology 20(3):344-347, 2009. https://doi.org/10.1097/EDE.0b013e31819d1092

**[8]** Chan LYH, et al. Estimating the generation time for influenza transmission using household data in the United States. Epidemics 50:100815, 2025. https://doi.org/10.1016/j.epidem.2025.100815

**[9]** Wallinga J, Lipsitch M. How generation intervals shape the relationship between growth rates and reproductive numbers. Proceedings of the Royal Society B 274(1609):599-604, 2007. https://doi.org/10.1098/rspb.2006.3754

**[10]** Centers for Disease Control and Prevention. About Estimated Flu Burden. CDC website, 2025. https://www.cdc.gov/flu-burden/php/about/index.html

**[11]** Reed C, Chaves SS, Daily Kirley P, et al. Estimating Influenza Disease Burden from Population-Based Surveillance Data in the United States. PLOS ONE 10(3):e0118369, 2015. https://doi.org/10.1371/journal.pone.0118369

**[12]** Hayward AC, Fragaszy EB, Bermingham A, et al. Comparative community burden and severity of seasonal and pandemic influenza: results of the Flu Watch cohort study. Lancet Respiratory Medicine 2(6):445-454, 2014. https://doi.org/10.1016/S2213-2600(14)70034-7

**[13]** Endo A, Abbott S, Kucharski AJ, Funk S. Estimating the overdispersion in COVID-19 transmission using outbreak sizes outside China. Wellcome Open Research 5:67, 2020. https://doi.org/10.12688/wellcomeopenres.15842.2

**[14]** Fay RE, Herriot RA. Estimates of Income for Small Places: An Application of James-Stein Procedures to Census Data. Journal of the American Statistical Association 74(366a):269-277, 1979. https://doi.org/10.1080/01621459.1979.10482505

**[15]** Spiegelhalter DJ. Funnel plots for comparing institutional performance. Statistics in Medicine 24(8):1185-1202, 2005. https://doi.org/10.1002/sim.1970

**[16]** Spiegelhalter D, Sherlaw-Johnson C, Bardsley M, et al. Statistical methods for healthcare regulation: rating, screening and surveillance. Journal of the Royal Statistical Society Series A 175(1):1-47, 2012. https://doi.org/10.1111/j.1467-985X.2011.01010.x

**[17]** Cori A, Ferguson NM, Fraser C, Cauchemez S. A New Framework and Software to Estimate Time-Varying Reproduction Numbers During Epidemics. American Journal of Epidemiology 178(9):1505-1512, 2013. https://doi.org/10.1093/aje/kwt133

**[18]** Gostic KM, McGough L, Baskerville EB, Abbott S, et al. Practical considerations for measuring the effective reproductive number, Rt. PLOS Computational Biology 16(12):e1008409, 2020. https://doi.org/10.1371/journal.pcbi.1008409

**[19]** Wallinga J, Teunis P. Different Epidemic Curves for Severe Acute Respiratory Syndrome Reveal Similar Impacts of Control Measures. American Journal of Epidemiology 160(6):509-516, 2004. https://doi.org/10.1093/aje/kwh255

**[20]** Parag KV, Donnelly CA. Adaptive Estimation for Epidemic Renewal and Phylogenetic Skyline Models. Systematic Biology 69(6):1163-1179, 2020. https://doi.org/10.1093/sysbio/syaa035

**[21]** Bhatt S, Ferguson N, Flaxman S, Gandy A, Mishra S, Scott JA. Semi-mechanistic Bayesian modelling of COVID-19 with renewal processes. Journal of the Royal Statistical Society Series A 186(4):601-615, 2023. https://doi.org/10.1093/jrsssa/qnad030

**[22]** Owen AB. The pigeonhole bootstrap. Annals of Applied Statistics 1(2):386-411, 2007. https://doi.org/10.1214/07-AOAS122

**[23]** Owen AB, Eckles D. Bootstrapping data arrays of arbitrary order. Annals of Applied Statistics 6(3):895-927, 2012. https://doi.org/10.1214/12-AOAS547

**[24]** Frutos AM, et al. Evaluation of a Novel Data Source for National Influenza Surveillance: Influenza Hospitalization Data in the National Healthcare Safety Network, United States, September 2021–April 2024. American Journal of Public Health 116(2):251-255, 2026. https://doi.org/10.2105/AJPH.2025.308282
