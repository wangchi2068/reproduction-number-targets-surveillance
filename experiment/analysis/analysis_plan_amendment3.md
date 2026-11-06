# Analysis plan, Amendment 3 (third revision round), 2026-11-06

Paper: "What more cases buy: inferential targets for reproduction-number estimates and what cross-sectional dispersion diagnostics measure".

This amendment belongs to `analysis_plan.md`, which holds the locked plan and Amendment 1, and to `analysis_plan_amendment2.md`. Its sha256 and UTC time are recorded in `analysis_plan_amendment3.lock`, and the lock was written before any analysis listed below was run. Nothing in the locked plan or in Amendments 1 and 2 is edited.

All analyses below respond to the third review round (items N10, N11, N12, N13, N14 and N15). The confirmatory analysis (E1, Section 2 of the plan) and its decision rule are unchanged. No analysis below is confirmatory: each is a pre-specified secondary analysis, and none can change the E1 verdict, alpha, the calibration, or the data. Every result is reported whatever its direction. If a result goes against the current main line, the main line is rewritten rather than the specification changed.

Three facts were established from the data and the shipped results before this amendment was written, and they fix the scope below: (i) the relative bias of the fitted intercept against the transmission-plus-reporting total ranges from -6.3% to -29.3% across the fresh-draw mechanism cells already in `results.json`, so the abstract's "6--22%" is not the full range; (ii) the sensitivity analysis S-a5 at threshold 1.02 selects exactly the primary sample (1,748 rows, identical `a` and interval), so it is a duplicate check and not an independent sensitivity; (iii) the NHSN vintage carries COVID-19 and RSV admission columns with the same jurisdiction and week structure as influenza, so a cross-series replication is feasible without new data.

## A3.1 Population gradient under a free count exponent (review N10)

Script `rev_r3_population.py`. Sample and `z` exactly as in Amendment 2, Section B.2: `z` is the log of each jurisdiction's mean of `POPESTIMATE2022`--`POPESTIMATE2025`, standardised to mean 0 and SD 1 over the 51 jurisdictions, fixed per state. `L`, `x`, `D` and the sample filters are those of the locked plan, Section 1.

Two models are fitted by the same gamma quasi-likelihood used for E1:

- **M1 (affine with population):** `E(L|C,z) = a0 + a1*z + b/C`.
- **M2 (free exponent with population):** `E(L|C,z) = a0 + a1*z + b*C^(-gamma)`, `gamma > 0`.

Inference for both: two-way (state x calendar week) CGM cluster-robust sandwich on the stacked estimator, Wald t with `G - 1 = 41` degrees of freedom, alpha = 0.05, exactly as for E1. Seeds `[20261105, 80, j]`.

Reported quantities, all with 95% intervals:

1. `a1` in M1 and in M2, its standard error, t and two-sided p.
2. The implied intercept ratio between the 10th and 90th population percentiles, `a(z_P10)/a(z_P90)`, computed from `a0 + a1*z` in each model, with a state-cluster bootstrap interval (B = 1,999).
3. The matched-count comparison already computed for the review, moved into the main text as a descriptive table: within count bands 20--49, 50--99, 100--199 and 200--399, the mean loss net of the fitted counting term by population tertile, with a state-cluster bootstrap interval (B = 999) on the small-minus-large difference. Bands with fewer than 10 states in a tertile are reported but not tested.

**Decision rule, fixed here.** The population gradient is called *robust to the functional form* if `a1 < 0` and its 95% two-way interval excludes zero in M2. If it does not, the abstract and the conclusion drop the magnitude claim and state that the gradient is present under the affine form but does not survive a free count exponent. Either way the abstract reports the ratio together with the model that produces it; the phrase "halved" is removed unless the ratio survives in M2.

## A3.2 Conditioning of the b/b_P interval on the reference (review N12)

Script `rev_r3_reference.py`. No new estimator is introduced. The three reference definitions already reported in Table 12 (row-normalised plug-in, exact Poisson, simulated exact) are recomputed with the calendar-week block bootstrap of `rev_intercept.py` (B = 1,999, seed `[20261105, 22]`), and the following are reported for each: the point value of `b/b_P`, its 95% percentile interval, and the share of draws with `b/b_P <= 1`.

**Decision rule, fixed here.** The text states that each interval is conditional on its reference definition and that the three references are not nested. The claim that "the lower limit lies just above 1" is retained only if the lower limit exceeds 1 under all three references; otherwise the text reports the point range and says that whether `b/b_P` exceeds 1 depends on the reference. No joint resampling of a parametric reference is claimed, because the plug-in reference is a deterministic functional of the observed counts and a parametric reference would require a model the paper does not estimate.

## A3.3 Bias range over all simulated mechanisms (review N11)

Script `rev_r3_bias_range.py`. No new simulation. The script reads `results.json`, key `revision_2026_11_r2.C_total_bias.generators`, and reports for every fresh-draw cell with a successful calibration: the relative bias of the fitted intercept against the transmission-plus-reporting total and its 95% interval. The reported range is the minimum and maximum of the point estimates over those cells. The state-realistic cells of `revision_2026_11_r2.D1_state_realistic_generators` are reported separately, with their calibration status and their `a_over_total` ratio; they do not enter the numeric range because none of them passes the six declared tests.

**Decision rule, fixed here.** The abstract reports the range as "6--30%" only if the computed minimum and maximum round to those bounds; otherwise it reports the computed bounds to the nearest integer. The discussion states explicitly that no shipped generator reproduces `b/b_P`, the upper tail of `C_t` and the pooled-to-within slope ratio at the same time, and that the mechanism is therefore not identified from the data.

## A3.4 Cross-series replication on NHSN COVID-19 and RSV admissions (review N10 and the reviewer's recommended addition)

Script `rev_r3_crosspathogen.py`. Data: the same frozen NHSN vintage, columns `Total COVID-19 Admissions` and `Total RSV Admissions`, same 51 jurisdictions and week grid, same `R_A` transform with the same generation interval (mean 3.6 d, SD 1.6 d).

**Window rule, fixed here and applied without exception.** For each pathogen and season in `RULE_START`, let `C_nat(t)` be the national count and `g(t) = C_nat(t+1)/C_nat(t)`. The season peak is the week of the maximum national count between the season start and 30 April. The analysis window is the maximal contiguous run of weeks ending at that peak in which every `g(t) >= 1.0`. A season is included only if the run has at least 4 weeks and at least 45 of the 51 jurisdictions report a non-missing count in at least 90% of the run's weeks. Seasons that fail either condition are listed as excluded with the reason, and no alternative window is substituted. The rule uses national counts only and never the outcome `L`.

Sample filters after the window: identical to the locked plan, Section 1 (four consecutive weeks, `D_t >= 5`, `C_t >= 5`).

Estimands, each reported with its interval:

1. The affine intercept `a` and its two-way cluster-robust 95% interval, and the count exponent `gamma` from the free-exponent model.
2. The model-free floor: the relative RMSE of `R_A` above the 90th percentile of `C_t`, with a calendar-week block bootstrap interval (B = 1,999).
3. The population coefficient `a1` in model M2 of A3.1.

**Replication criteria, fixed here.** A series *replicates the floor* if `a > 0` with a 95% interval excluding zero and the model-free floor lies inside the influenza 95% interval. A series *replicates the gradient* if `a1 < 0` with a 95% interval excluding zero. A series that does not replicate is reported as not replicating, in the abstract and in the results, and the paper's claim is limited to the series that do. Influenza is run through the same rule as an internal check and reported alongside; the primary influenza result is unchanged.

## A3.5 Sensitivity count and table notes (review N13, N14)

Script `rev_r3_tables.py`. No new estimation.

- The sensitivity table labels S-a5 at threshold 1.02 as a duplicate of the primary sample (identical window set, 1,748 rows), and the text and abstract report **16 independent sensitivity analyses**, not 17.
- The bootstrap-t table note states that the studentised statistic uses the stacked sandwich that includes the estimate of `xbar`, and that this choice is required for stability on clamped draws.
- Tables that only support supplementary material (state-realistic generator details, clamp and interval-type lists) move to the supplementary file; the main text keeps the tables that carry the abstract's numbers.

## A3.6 Reproduction and legacy files (review N15)

- `run_all.sh` gains a fourth lane that runs `rev_r3_population.py`, `rev_r3_reference.py`, `rev_r3_bias_range.py`, `rev_r3_crosspathogen.py` and `rev_r3_tables.py`, in that order, after the 2026-11 revision lane.
- The China COVID-19 code, its two data files and its two logs move to `experiment/legacy/china_covid19/`; `run_all.sh` no longer runs them. `experiment/code/README.md` records the move.
- Every script run of this round records its sha256 in `analysis/revision3/script_hashes.jsonl` before execution.
- All new numbers go to `results.json` under the key `revision_2026_11_r3`, with a `meta` block naming the script, the plan section, the seeds and the UTC time.

## A3.7 What happens if a result is unfavourable

If the population gradient does not survive M2, the abstract loses the magnitude and the conclusion loses the jurisdiction-size sentence. If the bias range is wider than 6--30%, the abstract reports the computed range. If a cross-series fails to replicate, that series is reported as a failure and the abstract says the floor was replicated only where it was. No specification above is changed after a result is seen; a change would require a further amendment, locked before the changed analysis is run.
