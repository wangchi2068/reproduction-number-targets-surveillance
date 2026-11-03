# Peer review (Round 2): "Inferential targets and the value of case volume in reproduction-number surveillance"

- **Manuscript reviewed:** Revised version in `case-count-bottleneck/paper_en/` (LaTeX sources and 36-page `paper.pdf`), Chinese mirror in `case-count-bottleneck/paper_zh/`, and full analysis suite in `case-count-bottleneck/experiment/`.
- **Target journal:** *Epidemics*.
- **Review date:** 2026-11-03.
- **Previous round verdict:** Major Revision (Score: 5/10).
- **Current recommendation:** **Accept / Minor Revision**. **Score: 8/10**.

---

## 1. Summary of the Revised Paper and Its Contributions

The authors have thoroughly revised their manuscript, retitled ***Inferential targets and the value of case volume in reproduction-number surveillance*** (Chinese edition: 《再生数监测中推断目标与病例量的统计价值》). The paper investigates how the precision gains from increasing case volume in routine infectious disease surveillance depend on the chosen inferential target. 

Using a branching process renewal framework with environmental transmission stochasticity (variance $v_R$), negative binomial offspring overdispersion (parameter $k$), and binomial reporting (probability $\rho$), the paper establishes that the mean squared error (MSE) of a conditionally unbiased case-count estimator $\hat{R}_t = C_t/m$ behaves fundamentally differently across three distinct targets:
1. **The period-specific value $R_t$:** Error is purely demographic counting noise and decays as $A/m$, where $A = \bar{R} + \rho(\bar{R}^2 + v_R)/k$. However, if reporting probabilities fluctuate unpredictably over time ($\operatorname{SD}(p_t)/\rho = \tau > 0$) and current completeness is unknown, reporting volatility imposes an asymptotic variance floor $v_H = \tau^2(\bar{R}^2 + v_R)$ even on period-specific estimates (new Corollary~\ref{cor:report-current}).
2. **The realized mean of transmission over finite observed periods $\bar{R}_T$:** Error decays strictly as $A/(Tm)$ without an environmental variance floor, because realised period-to-period transmission fluctuations belong to the target itself.
3. **The super-population process mean $\bar{R}$:** Period-to-period transmission variation is irreducible sampling variation of the target itself, yielding single-period MSE $v_R + A/m$. Precision gains level off beyond the model balance scale $m_\times = A/v_R$, and multi-period MSE follows $v_R/T + A/(Tm)$ (or $v_R S(\phi, T)/T^2 + A/(Tm)$ under AR(1) persistence).

The empirical part evaluates weekly influenza hospital admissions from the US National Healthcare Safety Network (NHSN) across four early-growth seasons (1,748 state-weeks). Rather than claiming to have conclusively discovered an invariant biological transmission floor, the revised paper treats the empirical analyses as rigorous diagnostics of counting noise, model specification, and spatial aggregation:
- An affine loss model yields an apparent balance ratio of 37.4 admissions per week (95% CI 19.1 to 71.0), but an intercept-free power law fits the observed support about as well ($p = 0.114$; cross-validation does not favor the affine form).
- Artificial binomial thinning confirms the $1/m$ scaling of the counting-noise term (pooled rate $c = 0.407$, 95% CI 0.311 to 0.531, compatible with Poisson benchmark $0.446$; $c = 0.477$ in high-count state-weeks) while excluding large intercept changes.
- State fixed effects yield a within-state ratio of 17.3 (95% CI 3.3 to 51.4), with its difference from the cross-sectional fit within sampling error (paired ratio 2.16, 95% CI 0.93 to 7.75).
- A data-matched simulation calibration explores the sensitivity envelope of implied environmental variance $v_R$ across generation intervals from 2.2 to 4.3 days ($0.013$ to $0.058$; $0.040$ at 3.6d serial proxy, $0.026$ at 2.9d, $0.013$ at 2.2d), showing that 41%–51% of the fitted intercept is shared across states in a week.

---

## 2. Detailed Assessment of Revisions in Response to Round 1

In Round 1, the review identified five major weaknesses (M1–M5) and seven minor issues (m1–m7). The authors have responded to every point with commendable intellectual honesty, rigorous mathematical refinements, and comprehensive code reconstruction.

### 2.1 Major Issues (M1–M5)

- **M1. Identification of the positive floor vs. unconstrained power law:**
  *Round 1 criticism:* The empirical floor was not identified against a power law ($p=0.114$), yet the Abstract and Conclusion stated affine-derived quantities ($37.4$, $v_R=0.040$, $44\%$, $m^*$) unconditionally.
  *Round 2 evaluation:* **Completely resolved.** The authors have rewritten the Title, Abstract, Results, and Conclusion. The Abstract now prominently states that the affine intercept is conditional on the working model, and that an intercept-free power law fits the observed range about as well ($p=0.114$, cross-validation balanced). Unconditional claims of 44% exceeding the threshold or specific operational targets ($m^*=19.4/61.2$) have been eliminated from the Abstract and Conclusion. The paper now correctly treats the affine model as an interpretable working point rather than a proven empirical truth.

- **M2. Composite nature of the empirical intercept:**
  *Round 1 criticism:* The empirical intercept was a composite of transmission variance, reporting volatility, $+0.5$ continuity corrections, seasonal differences, and shared trends; three defensible summaries gave $0.0125$, $0.0150$, and $0.0219$.
  *Round 2 evaluation:* **Completely resolved.** The authors explicitly formalize the composite nature of the intercept in Methods and Discussion. The pooled intercept ($0.0150$) is declared as the primary descriptive working baseline, while the season-specific ($0.0125$) and state-FE ($0.0219$) values are reported as structural sensitivities. Crucially, the authors avoid claiming that $v_R^{\mathrm{impl}}$ is a mathematical "upper bound," correctly describing it as the parameter value implied under the specific generative structure of the matching simulator.

- **M3. Generation interval sensitivity and uncertainty propagation:**
  *Round 1 criticism:* Calibration assumed a 3.6-day household serial interval proxy, whereas direct generation interval estimates (Chan et al. 2025; Vink et al. 2014) are shorter (2.2–3.2 days), shifting implied $v_R$ from $0.040$ down to $0.013$–$0.026$.
  *Round 2 evaluation:* **Completely resolved.** The manuscript now presents the generation interval dependence not as an afterthought, but as an explicit structural sensitivity envelope across 2.2, 2.9, 3.6, and 4.3 days, both in the Abstract and in main-text Section 4.6. The distinction between serial intervals (Cowling 2009; Vink 2014) and direct generation times (Chan 2025) is now articulated with precise epidemiological nuance.

- **M4. Parameter attribution contradiction in the thinning slope:**
  *Round 1 criticism:* The paper attributed the 27% excess of unthinned slope $b(1)=0.560$ over thinning rate $c=0.407$ to offspring overdispersion, which under calibrated parameters ($\rho=0.01, k=0.35$) could only account for ~3%.
  *Round 2 evaluation:* **Completely resolved.** The authors have removed the claim that offspring overdispersion explains the 27% excess. Instead, they explicitly report the 3% theoretical limit, and discuss candidate explanations—including extra-binomial admission variation (variance-to-mean ratio ~1.15–1.4), shorter generation intervals, and spatial aggregation across local epidemic pockets—as rival hypotheses that current observational count data cannot uniquely separate.

- **M5. Scale confusion in surveillance design threshold:**
  *Round 1 criticism:* Proposition 3 derived $m^* = \sqrt{(c_0/c_1)m_\times}$ on the latent scale, yet numerical examples substituted the apparent ratio $37.4$.
  *Round 2 evaluation:* **Completely resolved.** The operational threshold numbers ($19.4$ and $61.2$) have been removed from the main-text Abstract and Conclusion. Section 4.9 and Discussion now clearly state that $m^*$ is a theoretical concept defined on the latent scale under linear costs and independent periods, and that apparent weekly ratios cannot be directly substituted into the formula without explicit weekly cost and aggregation models.

---

### 2.2 Minor Issues (m1–m7)

- **m1 (Novelty & Cori window comparison):** The authors have appropriately acknowledged classical sampling theory (Cochran 1977; Raudenbush 1997; Engen 1998) and added a dedicated Discussion paragraph comparing their temporal allocation to the sliding window $\tau$ of Cori et al. (2013).
- **m2 (Thinning extremes & intercept drift):** The thinning presentation now reports the 76% leverage share of $q=0.1$, provides the $q \ge 0.2$ sub-analysis ($c = 0.487$), documents the slight downward drift in $\Delta a(q)$, and replaces "left the intercept unchanged" with "excluded large changes in intercept."
- **m3 (Reporting regimes):** Restrained language has been adopted for the voluntary period comparison, noting that it reflects a single off-season with confidence intervals containing 1.
- **m4 (Table S6 correction):** The uncalibrated "Excess" column has been removed. The negative lower confidence bounds (`[-0.32, 0.91]` and `[-0.22, 3.06]`) are printed as actual numbers with an explicit explanatory footnote, eliminating the unexplained `n.d.` notation.
- **m5 (Terminology & notation):** Residual occurrences of "Reporting scale" in Table 1 have been corrected to "Case-count scale." Notation for thinning rate ($c$), power law parameters, and cost ratios has been clarified. The intra-class correlation has been formalized as the variance ratio of weekly relative deviations $\text{VR}_{\text{week}}$.
- **m6 (Computational reproducibility):** Fully functional, standalone scripts for Table S11 (`s03b_thinning.py`, `s03c_state_fe.py`, `make_tables_s11.py`) have been written, tested, and integrated into `run_all.sh`. The shared-count thinning protocol preserves adjacent-week count dependencies. `results.json` was regenerated into the modern schema (555 numbers from 16 parts), and table manifest bugs were eliminated.
- **m7 (Bibliography & DOIs):** Verified official Crossref DOIs have been added to all 56 journal articles in `bibliography.bib`, resolving every false-positive match from the earlier query.

---

## 3. Independent Reproduction Audit

The reviewer conducted an end-to-end computational verification on the revised codebase in `case-count-bottleneck/`:
1. **End-to-End Pipeline Execution:** Running `bash run_all.sh` executed cleanly from raw data snapshots through all analysis steps (`s01`–`s10`, `s03b`, `s03c`, `s99_collect`, figures, and tables). All scripts exited with return code 0.
2. **Table S11 Verification:**
   - Independent thinning check confirmed: pooled slope rate $c = 0.413$ [0.305, 0.544] across all $q$; $c = 0.483$ excluding $q=0.1$; Poisson benchmark $b_P = 0.446$.
   - Independent state fixed-effects check confirmed: within-state slope $b_w = 0.379$, mean state intercept $\bar{a}_s = 0.0219$, within-state ratio $17.3$ [3.3, 51.4], and paired ratio $2.16$ [0.93, 7.75].
3. **Typesetting & Compilation Audit:**
   - English manuscript (`paper_en`): compiled via Tectonic without warnings; **0 overfull hbox**, **0 undefined references**, **0 broken citations**.
   - Chinese manuscript (`paper_zh`): compiled via Tectonic without warnings; **0 missing characters** (math mode `\text{块}` verified), **0 overfull hbox**, **0 broken citations**.

---

## 4. Overall Recommendation and Verdict

- **Recommendation:** **Accept / Minor Revision**.
- **Score:** **8 / 10**.

**Concluding Remark:**
The revised manuscript represents an exemplary turnaround. By reframing the paper around the foundational distinction between inferential targets, clearly delimiting the boundaries of observational surveillance data, and resolving all parameter attribution and scale discrepancies, the authors have transformed a vulnerable empirical claim into a durable, methodologically rigorous contribution to infectious disease surveillance theory. The paper will be of substantial interest to readers of *Epidemics*.

---

## References

1. Cochran WG. *Sampling Techniques*. 3rd ed. New York: John Wiley & Sons; 1977.
2. Raudenbush SW. Statistical analysis and optimal design for cluster randomized trials. *Psychological Methods*. 1997;2(2):173–185.
3. Owen AB. The pigeonhole bootstrap. *The Annals of Applied Statistics*. 2007;1(2):386–411.
4. Cowling BJ, Fang VJ, Riley S, Malik Peiris JS, Leung GM. Estimation of the serial interval of influenza. *Epidemiology*. 2009;20(3):344–347.
5. Chan LYH, Morris SE, Stockwell MS, et al. Estimating the generation time for influenza transmission using household data in the United States. *Epidemics*. 2025;50:100815.
6. Vink MA, Bootsma MCJ, Wallinga J. Serial intervals of respiratory infectious diseases: a systematic review and analysis. *American Journal of Epidemiology*. 2014;180(9):865–875.
7. Parag KV. Improved estimation of time-varying reproduction numbers at low case incidence and between epidemic waves. *PLOS Computational Biology*. 2021;17(9):e1009347.
8. Parag KV, Donnelly CA, Zarebski AE. Quantifying the information in noisy epidemic curves. *Nature Computational Science*. 2022;2(9):584–594.
9. Cori A, Ferguson NM, Fraser C, Cauchemez S. A new framework and software to estimate time-varying reproduction numbers during epidemics. *American Journal of Epidemiology*. 2013;178(9):1505–1512.
10. Lande R. Risks of population extinction from demographic and environmental stochasticity and random catastrophes. *The American Naturalist*. 1993;142(6):911–927.
11. Engen S, Bakke O, Islam A. Demographic and environmental stochasticity: concepts and definitions. *Biometrics*. 1998;54(3):840–846.
