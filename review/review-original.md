# Peer review: "Inferential targets for reproduction-number estimates: counting error, between-period variation and what cross-sectional surveillance diagnostics measure"

Manuscript reviewed: https://github.com/wangchi2068/reproduction-number-targets-surveillance (commit 6bbb8458a5ba3a9ba644d755320df55dbe06746b; `paper.pdf`, 32 pages, with LaTeX sources, code and data).

Basis of this review: every page and figure of the PDF was read. The propositions were checked symbolically and by simulation. The full analysis pipeline was rerun in the pinned environment and its outputs compared number by number with the shipped results. Each decimal number in the text was traced to the shipped result files. Section and line references below are to `paper/sections/*.tex`; table and figure numbers are those of the PDF.

## 1. Summary

The manuscript asks what extra case counts buy when the reproduction number is estimated, and argues that the answer depends on the inferential target. In a one-generation branching model with negative binomial offspring (dispersion k) and binomial reporting (probability ρ), it derives exact mean squared errors for three targets: the current value R_t, the realised mean over the observed periods, and the process mean R̄. Counting error falls as A/m in the expected count m for every target, but the process mean keeps a between-period variance v_R that more cases cannot remove; the two terms cross at m× = A/v_R. Further results cover:

- the cost of dividing by the previous observed count rather than a known scale (about 1.5 times the counting error at the illustrative parameters);
- the inability of single-period count moments to separate period-level reporting volatility from transmission variation (Proposition 2);
- AR(1) window factors, optimal allocation and a coverage corollary (Appendix B).

The empirical part fits an affine model, E(loss | count) = a + b/C, to 1,748 state-weeks of US influenza hospital admissions (2022–2026). It gives an intercept a = 0.0150 (0.0101–0.0206), a slope 1.26 (0.81–1.73) times the independent-Poisson reference, and a dispersion-to-counting crossover b/a of 37 admissions per week (18–72). A known-truth audit with ten synthetic generators (S0–S8) on the real design asks what these diagnostics recover. A further section replicates the diagnostics on weekly COVID-19 cases in 34 Chinese provincial units (2020–2023). The paper concludes that cross-sectional dispersion diagnostics describe a count-independent relative dispersion plus counting noise; they are not a precision limit or a sample-size rule, and separating under-reporting, superspreading and structured transmission needs temporal or external data. A four-item reporting checklist closes the discussion.

## 2. Strengths

1. **The central distinction is useful and correctly worked out.** Separating the current value, the realised mean and the process mean, and showing that counting error and between-period variation behave differently, is a clear way to discipline how Rt uncertainty is reported. All propositions, corollaries and closed forms that were checked are correct. These include the conditional moments, the Fisher information and attainment of the Cramér–Rao bound, the AR(1) window factor for T = 1 to 8, the gap formula (A.5), the reporting-volatility result, the allocation rule and Corollary B.3. The observed-denominator inflation in Table 2 was reproduced by independent simulation: 1.59, 1.51 and 1.49 at m = 50, 200 and 1,000 against 1.591, 1.509 and 1.487.
2. **Claims are kept descriptive where the data are weak.** The abstract, highlights and conclusion avoid over-reading the crossover count, and the paper says plainly that counting slopes and thinning do not identify overdispersion.
3. **The known-truth audit runs on the real design.** Using the actual 1,748-row structure, shared counts and filters, with count alignment to the observed median, is better practice than auditing on idealised panels. It reveals that the intercept tracks transmission variance only without reporting shocks and without scale-dependent heterogeneity, and the main text reports this.
4. **Sensitivity analysis is extensive.** It covers filters, continuity corrections, window rules, seasons, centring, fixed effects, block lengths, the generation interval (a cited influenza household serial interval [1] as primary, with an intrinsic generation-time estimate as a sensitivity analysis and the distinction between the two acknowledged) and held-out model comparison.
5. **The work is largely reproducible.** Code, data snapshot and pinned environment are public. Rerunning the pipeline gives the same values for every headline number, for the empirical analyses and for generators S0–S7b (exact or within 1e-3 relative).

## 3. Weaknesses

1. The theoretical results are standard consequences of the law of total variance, the negative binomial variance function and the Cramér–Rao bound. The paper does not say clearly which results are new. Nor does it show that the practice it corrects (reading a dispersion crossover as a sample-size rule) exists.
2. The known-truth audit is not designed as a simulation study with pre-specified aims and performance measures. Each generator is tuned to reproduce the observed intercept. The combined generator S8 was added with an extra tuned parameter and is then presented as satisfying "three empirical constraints", one of which it meets by construction.
3. The interpretation of the intercept as transmission variation is not identified. Across audited mechanisms the intercept is 0.43 to 2.02 times the true transmission variance, and one sentence in the application section states the direction of the bias for S8 wrongly.
4. The China COVID-19 replication has design problems that make it uninformative: selection on the outcome, a pseudo-unit for unassigned cases, a generation interval taken from another setting and variant, and missing citations for both the interval and the data. Its intervals are compatible with almost any pattern, yet the text says they show the pattern "the framework predicts".
5. Several numbers in the text cannot be regenerated from the shipped code. The drivers depend on paths outside the repository and hide failures.
6. A supplementary table that the main text relies on is partly missing from the PDF. At 32 pages with 20 tables, the manuscript is longer and denser than its contribution needs.

## 4. Detailed comments

### Major comments

**M1. State what is new, and show that the problem exists in practice.**
- *Location:* introduction.tex (contribution list); highlights 1 and 5; model.tex Propositions 1–2; Appendix B.
- *Evidence:* The MSE decompositions are direct applications of the law of total variance to a negative binomial–binomial model. Efficiency of the estimator C/m is the standard Cramér–Rao result, and the allocation rule is the textbook square-root rule. The introduction cites the main Rt frameworks [2], [3], [4], [5], funnel-plot overdispersion [6] and offspring overdispersion [7], but no study that treats a dispersion crossover, or the counting-error formula, as a sample-size or precision threshold. Highlight 5 ("A crossover near 37 weekly admissions is descriptive, not a sample-size rule") therefore answers a reading nobody is shown to hold.
- *Requested change:*
  - Add a short paragraph that separates known results (cite a standard reference for each) from the paper's contributions: the target taxonomy, the observed-denominator correction, Proposition 2 and the audit.
  - Give at least one documented instance in which surveillance practice or guidance conflates the targets or uses count-based precision as a threshold, for example in how intervals from [2] or [5] are interpreted. Otherwise, rephrase the motivation and highlight 5 as a clarification rather than a correction.

**M2. Redesign and report the known-truth audit as a pre-specified simulation study.**
- *Location:* diagnostics.tex:60–92 (design and results); Tables 3–4; Figure 2; discussion.tex:33.
- *Evidence:*
  - Each heterogeneity amplitude is recalibrated until the Monte Carlo mean of the fitted intercept equals the observed a = 0.014955 (`ma_audit.py`, alignment step). Every generator with heterogeneity therefore reproduces the observed intercept by construction. The text says this for S7a but not when presenting S8.
  - The S8 paragraph (diagnostics.tex:90) presents a fitted intercept of 0.0150 (95% interval 0.0131 to 0.0172, original results file) as one of three empirical constraints S8 "simultaneously satisfies". S8 also has a reporting-shock scale s that is "aligned separately", so it has one more tuned parameter than any other generator.
  - The shipped results file lists S8 in neither its alignment block nor its self-checks, although the shipped code writes both for every generator. S8 was therefore added in a separate run after the other generators.
  - The curve-shape result for S8, "at the 4.5–8.0% quantile across random streams", rests on two runs: 0.080 in the shipped results and 0.045 when the shipped code is rerun. With 200 panels per generator the quantile resolution is 0.005, so a 5% threshold is decided by about ten panels.
  - No Monte Carlo standard errors are given for any quantile or recovery ratio, all generators use one seed, and the performance measures are not stated in advance.
  - These are the conditions under which a post hoc best-fitting generator looks like confirmation [8].
- *Requested change:* Restructure the audit around the ADEMP elements [9]. Specifically:
  - add a table of aims, data-generating mechanisms, estimands, methods and performance measures written before results;
  - mark in Tables 3–4 which statistics are calibration targets and which are tests;
  - report Monte Carlo SEs, or quantiles over at least 5 independent streams, for the recovery ratios and the contrast quantile;
  - recast S8 as one member of a grid over (γ, s) rather than a single tuned point, reporting the region of that grid compatible with the three statistics that are not calibration targets.

**M3. The intercept does not identify transmission variance, and one statement of the bias direction is wrong.**
- *Location:* application.tex:7; Table 3 ("True rel. variance" and "a / truth" columns); abstract sentence on the intercept.
- *Evidence:*
  - application.tex:7 says a = 0.0150 corresponds to a relative SD of about 12%. It calls this an upper bound for transmission SD when heterogeneity does not scale with size, and says that "under scale-dependent mechanisms (S7a, S7b, S8) the intercept understates transmission variance".
  - Under S8 the intercept is 1.69 (1.48–1.94) times the true transmission relative variance; it overstates transmission variance. It understates only the transmission-plus-reporting total (0.86). The confusion comes from Table 3, where "True rel. variance" means transmission only for S1–S5 and S7, but transmission plus reporting for S6 and S8 (0.0147, 0.0175). The header does not say so.
  - Across the audited mechanisms the intercept is 0.43 (S7b) to 2.02 (S6) times the transmission variance. The 12% figure therefore corresponds to a transmission relative SD of about 9% to 19%, and the direction of the error is unknown.
- *Requested change:*
  - Correct the sentence in application.tex:7.
  - Split the Table 3 truth column into transmission-only and total, with a/truth for each.
  - Report the implied transmission SD as a range across mechanisms rather than a single 12% value with a one-sided qualifier. The abstract wording ("transmission variation under fixed reporting, and its sum with reporting shocks otherwise") is accurate and can be kept, but should add "and smaller than transmission variation when heterogeneity scales with count".

**M4. The China COVID-19 replication should be redesigned or removed.**
- *Location:* application.tex:55–75; Table 6; `experiment/code/china_covid_empirical.py`.
- *Evidence:* The section has six problems.
  - *Selection.* Rows enter when their own growth factor (C_{t+1}+0.5)/(C_t+0.5) is at least 1.05. R_A is a monotone transform of that same ratio, so the rule truncates the outcome and compresses its variance from below, which biases a and b in an unknown direction. The text calls this "the direct analogue of the early-growth rule applied to the US sample". It is not: the US windows are national early-growth segments (diagnostics.tex:9), which do not select on each row's own outcome.
  - *Units.* The 34 units include an "Unknown" pseudo-unit for cases not assigned to a province. It contributes 70 of the 896 rows, the largest share of any unit. The units also include Hong Kong SAR and Macao SAR, whose surveillance systems differ from the mainland's. Taiwan, China is not in the extract, and the text does not say why or state the inclusion rule.
  - *Generation interval.* The analysis applies a gamma interval with mean 5.2 days and SD 1.7 days to all of 2020–2023. These values match the Singapore estimate of Ganyani et al. [10], who give 3.95 days for Tianjin, China; Omicron intervals were shorter still [11]. The Wallinga–Lipsitch transform [12] is sensitive to this choice. No source is cited for the interval or for the JHU CSSE data [13].
  - *Data handling.* Negative daily revisions are set to zero with no sensitivity analysis.
  - *Inference.* The intervals are uninformative: b/a = 23.2 (−3.0 to 99.9) and b/b_P = 6.16 (−1.84 to 13.37). Even so, the text says the contrast with the US "is structured in the way the framework predicts rather than uniform".
  - *Reproducibility.* None of these numbers appear in the main results file.
- *Requested change:* Either remove the section or redo it:
  - select rows with a regional or national window rule that does not condition on the row's own outcome;
  - exclude the unassigned-case pseudo-unit and state the unit list and its rationale;
  - use period-specific generation intervals with citations, plus a sensitivity analysis;
  - cite the data source;
  - analyse negative revisions;
  - restrict claims to what the intervals support (for example, "the intercept is larger; the crossover is not estimable with useful precision").

**M5. Make every reported number regenerable from the shipped repository.**
- *Location:* `experiment/code/` drivers; diagnostics.tex:46, 62, 64, 75, 90; discussion.tex:33; availability.tex.
- *Evidence:*
  - Running the shipped code from a clean clone fails or silently departs from the shipped results. `ma_audit.py:67` reads a file from a `revision2/` directory outside the repository and crashes after six minutes. `run_authors_pipeline.sh:20` copies calibration files from `../../review/repro/out/`, so those files cannot be regenerated. Default interpreter paths point outside the repository.
  - `run_ma_audit.sh` prints the Python exit code and then exits 0, so failures are invisible.
  - Once the crash is bypassed, the audit reports a failed self-check (the S8 intercept-error check, which is inappropriate for a generator with reporting shocks and is exempted for S6) and exits with code 3.
  - Rerun from the shipped code, S8 gives s = 0.133, intercept 0.0148 (0.0125–0.0171), slope ratio 1.306 (1.081–1.520) and contrast quantile 0.045. The text mixes s = 0.133 from this run with the intercept 0.0150 (95% interval 0.0131 to 0.0172, original results file) and slope ratio 1.291 (1.100–1.508) of the shipped results.
  - Two claims have no code at all: the upper-tail reweighting result "0.001–0.043 (at most 6.3%)" (diagnostics.tex:64), and g0 = 1.407, which is computed inside the audit but not written to any output.
  - `paper.html` is an older version that lacks S8, the China section and the checklist. The repository README gives 30 pages and four highlights against 32 and five.
- *Requested change:*
  - Remove all paths outside the repository.
  - Provide one command that rebuilds every result file and table from a clean clone, and make drivers exit non-zero on any error or failed self-check.
  - Fix the S8 self-check.
  - Write every number used in the text (including g0, s and the reweighting result) to the results file and regenerate the text tables from it.
  - Rebuild `paper.html` from the current sources, and archive the release with a persistent identifier.

**M6. Check interval validity under the dependence the paper itself identifies, and justify the analysis windows.**
- *Location:* appendix.tex:81; discussion.tex:35; diagnostics.tex:9 and 50.
- *Evidence:*
  - The appendix derives that adjacent rows share a count, so their errors have correlation about −1/2. It then states that ignoring this makes multi-week intervals "slightly conservative". That is true for linear sums of errors. But the affine fit is built on squared relative deviations, and for Gaussian errors with correlation −1/2 the squares have correlation +1/4, so intervals for a and b may be too narrow rather than too wide.
  - The pigeonhole bootstrap [14] and the block sensitivities partly address dependence (see also [15]), but their coverage has not been checked, although the audit machinery could do so with known truth.
  - The windows are described as "inherited from the original analysis pipeline and not tuned here" (diagnostics.tex:9, discussion.tex:35). This does not tell the reader how the rule was chosen, or whether it was chosen with knowledge of these outcomes.
- *Requested change:*
  - Report the coverage of the bootstrap intervals for a, b and b/a under S1, S6 and S8, using the existing audit panels.
  - State the window rule as a definition with its rationale rather than as a pipeline inheritance.
  - Add a sensitivity analysis over the growth threshold (for example 1.02, 1.05 and 1.10) if Table 5 does not already cover it.

**M7. Fix the missing table and shorten the manuscript.**
- *Location:* pages 25–26 (Tables S8–S9); application.tex:18; whole manuscript.
- *Evidence:*
  - Table S8 (within-state and two-way slopes) is cut off at the foot of page 25: only its first row is visible, and five rows and the table note are lost.
  - Table S9 (`tab:fe`) has lost its first page, with the caption, the L = 1 rows and part of the L = 2 rows. Page 26 starts with "(continued)", and the caption "Table S9" appears nowhere in the PDF.
  - The main text cites this table for the cross-count ratios 2.16, 1.48 and 1.47, which readers therefore cannot check. The same defect appears when the sources are rebuilt.
  - The manuscript has 32 pages, 6 main and 14 supplementary tables and 7 figures. Many sentences exceed 60 words (diagnostics.tex:64 has 121), and the audit results paragraph (diagnostics.tex:84–90) lists more than 40 interval estimates in running text.
- *Requested change:*
  - Make Table S8 a longtable (or reduce it) and restore the first page of Table S9.
  - Move numerical detail from running text into tables.
  - Aim for a main text of about 20 pages, with the China section (if kept), the thinning analysis and the centring stages in the supplement.

### Minor comments

**m1. Justification of the τ² bound.** reporting.tex:7: "Reporting remains independent across infections given p_t, so τ² ≤ (1−ρ)/ρ". The bound follows from p_t ∈ [0,1] with mean ρ (Var(p_t) ≤ ρ(1−ρ)), not from independence across infections. Correct the reason.

**m2. The "1.5×" rule in the checklist is parameter-specific.** discussion.tex:43 (checklist item 2) recommends inflating sampling variance "by approximately 1.5×". The delta-method ratio is 1.48 at R̄ = 1.15, v_R = 0.04, k = 0.35, ρ = 0.25 and large m. Under the same parameters the exact ratio is 2.39 at m = 10 (Table 2). Replace the constant with the formula (Eq. observed-denominator) or a short lookup table by m. "Right-shift shrinkage half-weights accordingly" should give the shifted value explicitly.

**m3. Small-denominator treatment in Table 2.** Independent simulation gives 2.392 at m = 10 against 2.42 in Table 2; the other values agree to the third decimal. State how zero and small observed denominators are handled and report the Monte Carlo SE.

**m4. Post hoc normalisation by the S0 recovery factor.** diagnostics.tex:86 divides the simulated branching-Poisson slope ratio by the S0 recovery factor (0.930) to obtain "exact alignment" with the analytic 0.861. This is a post hoc rescaling chosen to make the two agree. The supporting claim of a "common finite-sample downward bias of gamma FGLS" has no citation. Either cite and quantify the bias, or report the unnormalised comparison only.

**m5. Week-centred reduction is undefined when a ≈ 0.** Table 3, last column: for S0, S2 and S2P the relative reduction is −30.74 (−34.19 to 12.89), 0.44 (−20.56 to 35.68) and −2.64 (−65.20 to 22.66), ratios of near-zero quantities. Report "not defined" or give the absolute change.

**m6. Front matter.** main.tex:144–145 has an author placeholder and `\date{\today}`. declarations.tex:5 lists the CRediT roles under the heading "Conceptualization", with no contributor names. Complete these before submission.

**m7. Data and code availability.** availability.tex says that code for the audit and sensitivity analyses "is provided as supplementary material", but it is in the public repository. Give the repository URL, the release or commit, and an archived DOI.

**m8. Process language.** "inherited from the original analysis pipeline" (discussion.tex:35; similar wording at diagnostics.tex:9) describes the production history of the analysis rather than the method. Describe the rule itself.

**m9. Number provenance.** g0 = 1.407 (diagnostics.tex:46, 62) and s = 0.133 (diagnostics.tex:75, 90) are not in the results file. The shipped audit output has s = 0.134. Use one source for both (see M5).

**m10. Repository hygiene.** `paper/figures/` contains stale duplicates of figures no longer included. The README's page and highlight counts are out of date.

**m11. Notation load.** The paper uses R_A, G, D_t, b_P, r(ρ), m×, m×,H, A, U, Q and τ² with definitions spread across sections. A notation table at the start of Section 2 would help.

## 5. Recommendation and score

**Recommendation: major revision. Overall score: 6 / 10.**

The theory is correct and the target distinction is worth publishing, and the empirical work on US influenza data is careful and largely reproducible. The score is held down by the design and reporting of the known-truth audit (M2), an incorrect statement of the bias direction (M3), the China replication (M4), several numbers that the shipped code cannot regenerate (M5), and a table missing from the PDF that the main text relies on (M7). None of these touches the abstract's headline numbers, all of which reproduce exactly. All of them can be fixed with work the authors' own infrastructure already supports.

Scoring used a rubric fixed before scoring:

| Dimension | Weight | Score | Main reasons |
|---|---|---|---|
| D1 Originality and significance | 0.15 | 6 | Useful synthesis and framing; theory is largely textbook; the practice being corrected is not documented (M1) |
| D2 Technical correctness | 0.20 | 7 | All propositions verified; misstated reason for the τ² bound (m1); parameter-specific 1.5× rule presented as general (m2) |
| D3 Methodological rigour | 0.20 | 5 | Audit tuned to the observed intercept, post hoc S8 with an extra tuned parameter, one seed, no Monte Carlo SEs (M2); interval coverage unchecked (M6) |
| D4 Data and empirical validity | 0.15 | 5 | China section selects on the outcome and has unit and generation-interval problems (M4); intercept-to-transmission mapping not identified, with one wrong direction statement (M3) |
| D5 Reproducibility | 0.10 | 6 | Headline numbers and S0–S7b reproduce exactly; external paths, masked failures, S8 numbers and two claims not regenerable (M5) |
| D6 Clarity and presentation | 0.10 | 5 | Table S9 partly missing and Table S8 truncated in the PDF; 32 pages of dense running-text numbers; front-matter placeholders (M7, m6) |
| D7 Positioning and references | 0.10 | 6 | Good coverage of Rt methods; missing sources for China data and its generation interval; known versus new results not separated (M1, M4) |

Weighted mean: 0.15×6 + 0.20×7 + 0.20×5 + 0.15×5 + 0.10×6 + 0.10×5 + 0.10×6 = 5.75. Rounded to the nearest integer, this gives 6. The cap rules were checked:
- **Headline-claim cap (7).** No open Major comment affects a headline claim (abstract, highlights or conclusion), so this cap does not apply. It would not be binding anyway.
- **Theorem-or-headline-number cap (5).** No error was found in a main theorem, and every headline number regenerates, so this cap does not apply either.

The overall score is 6, which corresponds to major revision.

To reach 8 (minor revision), a revision would need to:
- resolve M2–M5 and M7;
- address M1 by separating known from new results and grounding the motivation, or by toning down the corrective framing;
- deal with the minor comments.

## 6. Reproducibility statement

The repository was cloned at commit 6bbb8458a5ba3a9ba644d755320df55dbe06746b. All drivers were run in a Python 3.11 environment with the pinned package versions (numpy 1.26.4, scipy 1.16.3, pandas 2.3.3, matplotlib 3.10.9). The run took about 35 minutes on 4 CPUs. One file reference outside the repository had to be bypassed to complete the audit; the change skips a single self-check and alters no computation.

Comparison of shipped and regenerated outputs:
- The authors' pipeline and the auxiliary analyses reproduce: of 76,755 numbers in the pipeline outputs, 99.9% agree to within 1e-3 relative. The 77 larger differences are bootstrap tail quantiles that move by one draw.
- The main results file reproduces in all 2,792 entries except those for generator S8; the S8 discrepancies are documented in M5.
- The China results reproduce exactly.

The LaTeX sources compile to the shipped 32-page PDF with no errors, including the table defect described in M7.

## References

[1] Cowling BJ, Fang VJ, Riley S, Malik Peiris JS, Leung GM. Estimation of the serial interval of influenza. Epidemiology. 2009;20(3):344–347. https://doi.org/10.1097/EDE.0b013e31819d1092

[2] Cori A, Ferguson NM, Fraser C, Cauchemez S. A new framework and software to estimate time-varying reproduction numbers during epidemics. American Journal of Epidemiology. 2013;178(9):1505–1512. https://doi.org/10.1093/aje/kwt133

[3] Gostic KM, McGough L, Baskerville EB, Abbott S, Joshi K, Tedijanto C, et al. Practical considerations for measuring the effective reproductive number, Rt. PLOS Computational Biology. 2020;16(12):e1008409. https://doi.org/10.1371/journal.pcbi.1008409

[4] Parag KV, Donnelly CA. Fundamental limits on inferring epidemic resurgence in real time using effective reproduction numbers. PLOS Computational Biology. 2022;18(4):e1010004. https://doi.org/10.1371/journal.pcbi.1010004

[5] Abbott S, Hellewell J, Thompson RN, Sherratt K, Gibbs HP, Bosse NI, et al. Estimating the time-varying reproduction number of SARS-CoV-2 using national and subnational case counts. Wellcome Open Research. 2020;5:112. https://doi.org/10.12688/wellcomeopenres.16006.1

[6] Spiegelhalter DJ. Funnel plots for comparing institutional performance. Statistics in Medicine. 2005;24(8):1185–1202. https://doi.org/10.1002/sim.1970

[7] Lloyd-Smith JO, Schreiber SJ, Kopp PE, Getz WM. Superspreading and the effect of individual variation on disease emergence. Nature. 2005;438(7066):355–359. https://doi.org/10.1038/nature04153

[8] Gelman A, Loken E. The statistical crisis in science. American Scientist. 2014;102(6):460. https://doi.org/10.1511/2014.111.460

[9] Morris TP, White IR, Crowther MJ. Using simulation studies to evaluate statistical methods. Statistics in Medicine. 2019;38(11):2074–2102. https://doi.org/10.1002/sim.8086

[10] Ganyani T, Kremer C, Chen D, Torneri A, Faes C, Wallinga J, et al. Estimating the generation interval for coronavirus disease (COVID-19) based on symptom onset data, March 2020. Eurosurveillance. 2020;25(17):2000257. https://doi.org/10.2807/1560-7917.ES.2020.25.17.2000257

[11] Backer JA, Eggink D, Andeweg SP, Veldhuijzen IK, van Maarseveen N, Vermaas K, et al. Shorter serial intervals in SARS-CoV-2 cases with Omicron BA.1 variant compared with Delta variant, the Netherlands, 13 to 26 December 2021. Eurosurveillance. 2022;27(6):2200042. https://doi.org/10.2807/1560-7917.ES.2022.27.6.2200042

[12] Wallinga J, Lipsitch M. How generation intervals shape the relationship between growth rates and reproductive numbers. Proceedings of the Royal Society B: Biological Sciences. 2007;274(1609):599–604. https://doi.org/10.1098/rspb.2006.3754

[13] Dong E, Du H, Gardner L. An interactive web-based dashboard to track COVID-19 in real time. The Lancet Infectious Diseases. 2020;20(5):533–534. https://doi.org/10.1016/S1473-3099(20)30120-1

[14] Owen AB. The pigeonhole bootstrap. The Annals of Applied Statistics. 2007;1(2):386–411. https://doi.org/10.1214/07-AOAS122

[15] Cameron AC, Gelbach JB, Miller DL. Robust inference with multiway clustering. Journal of Business & Economic Statistics. 2011;29(2):238–249. https://doi.org/10.1198/jbes.2010.07136
