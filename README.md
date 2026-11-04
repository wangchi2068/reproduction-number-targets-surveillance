# Inferential Targets for Reproduction-Number Estimates: Counting Error, Between-Period Variation and What Cross-Sectional Surveillance Diagnostics Measure

[![LaTeX Compilation](https://img.shields.io/badge/LaTeX-Tectonic%20%2F%20XeLaTeX-brightgreen.svg)]()
[![Python](https://img.shields.io/badge/Python-3.11-blue.svg)]()
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()

This repository contains the replication codebase, empirical pipeline, and complete LaTeX manuscript for the paper:

> **Inferential targets for reproduction-number estimates: counting error, between-period variation and what cross-sectional surveillance diagnostics measure**

The study investigates how case counts determine the statistical precision of effective reproduction number ($R_t$) estimates under explicit branching processes, proves an information-theoretic non-identifiability boundary under fluctuating reporting, and provides a known-truth simulation audit evaluating what commonly used cross-sectional surveillance diagnostics actually measure on real-world epidemic designs.

---

## Repository Structure

```
.
├── paper.pdf                  # Compiled 30-page manuscript PDF (clean, publication-ready)
├── paper.html                 # Self-contained offline interactive HTML manuscript
├── highlights.txt             # Four key scientific highlights for journal submission
├── paper/                     # Complete LaTeX / BibTeX manuscript source package
│   ├── main.tex               # Main LaTeX document entry point
│   ├── compile.sh             # Build script
│   ├── bibliography.bib       # Verified BibTeX database with official Crossref DOIs
│   ├── unsrtnat-noid.bst      # BibTeX citation style with clickable DOI support
│   ├── sections/              # Modular LaTeX section files (11 modules)
│   ├── figures/               # High-resolution publication vector figures (PDF)
│   └── tables/                # Programmatically rendered LaTeX tables (.tex)
└── experiment/                # Reproducible computational pipeline and data
    ├── code/                  # Python simulation, audit, and econometric scripts
    ├── analysis/              # Results bundle (results.json) containing all 461 validated estimates
    └── data/                  # Input datasets, vintage hospital admissions, and DATA_MANIFEST.md
```

---

## Key Methodological Contributions

1. **Taxonomy of Inferential Targets & Risk Structure**:
   - *Current instantaneous value ($R_t$)*: Mean squared error $\operatorname{MSE} = A/m$ decays to zero inversely with expected case count $m$ under constant reporting.
   - *Realised mean over $T$ observed periods ($\bar{R}_T$)*: $\operatorname{MSE} = A/(Tm)$, decaying strictly with total sample size without an environmental variance floor.
   - *Super-population stationary process mean ($\bar{R}$)*: $\operatorname{MSE} = v_R/T + A/(Tm)$ (or $v_R S(\phi, T)/T^2 + A/(Tm)$ under AR(1) persistence), retaining between-period transmission variance $v_R$. Extra case volume per period cannot remove cross-period sampling variation.
   - *Noisy denominator inflation*: When the denominator is itself an observed count $C_t$, sampling variance rises by $\sim 1.5\times$, shifting the empirical Bayes shrinkage half-weight scale $m_\times = A/v_R$ from 53.1 to 78.6.

2. **Reporting Volatility and Non-Identifiability (Proposition 2)**:
   - When period-specific reporting completeness varies with variance $\rho^2 \tau^2$, current-value risk retains an asymptotic reporting floor $v_H = \tau^2(\bar{R}^2 + v_R)$.
   - Single-period count moments depend on transmission variance $v_R$ and reporting volatility $\tau$ only through $U = (1+\tau^2)(\bar{R}^2 + v_R)$, proving that cross-sectional count moments cannot separate biological transmission variation from reporting noise.

3. **Known-Truth Diagnostic Audit on Real Surveillance Design**:
   - Synthetic audits across 9 data-generating mechanisms on the real 1,748 state-week design demonstrate that fitted affine dispersion intercepts ($a \approx 0.0150$) recover count-independent relative dispersion ($v_R$ under constant reporting; $v_R + v_H$ under volatile reporting) rather than an estimation-error floor.
   - Counting slopes, slope ratios to Poisson benchmarks, and binomial thinning responses do not separate offspring overdispersion from under-reporting.
   - Incorporating moderate scale-dependent heterogeneity with reporting shocks (S8 Combo) simultaneously reconciles the observed intercept ($0.0148$), crossover count ($41.7$, covering observed $37.4$), slope ratio ($1.298$, covering $1.257$), and curve contrast ($p=0.080$).

---

## Replication & Build Instructions

### Building the Manuscript PDF

The manuscript can be compiled cleanly using [Tectonic](https://tectonic-typesetting.github.io/) or standard XeLaTeX / BibTeX:

```bash
# Option 1: Compile with Tectonic (recommended, 0 dependencies)
cd paper
tectonic main.tex

# Option 2: Compile with XeLaTeX / BibTeX
cd paper
xelatex -interaction=nonstopmode main.tex
bibtex main
xelatex -interaction=nonstopmode main.tex
xelatex -interaction=nonstopmode main.tex
```

### Reproducing Computational Experiments

Python 3.11 with pinned dependencies (`experiment/code/requirements.txt`):

```bash
cd experiment/code
python -m pip install -r requirements.txt

# Run the known-truth audit across simulation generators (S0 to S8)
python ma_audit.py

# Run state fixed effects, dual Poisson baselines, and formal season tests
python ma_extensions.py

# Run exact binomial thinning summation
python ma_thinning_exact.py

# Render tables and compile results bundle
python ma_render.py
```

All 461 numerical estimates cited in the text and supplementary tables are programmatically verified and exported to `experiment/analysis/results.json`.

---

## Data Availability

Analyses utilize weekly jurisdiction-level influenza hospital admissions from the National Healthcare Safety Network (NHSN) frozen vintage as of October 20, 2026. Data integrity hashes and schemas are documented in `experiment/data/DATA_MANIFEST.md` and `experiment/data/data_contract.md`.

---

## Citation

```bibtex
@article{wang2026inferential,
  title={Inferential targets for reproduction-number estimates: counting error, between-period variation and what cross-sectional surveillance diagnostics measure},
  author={Wang, Chi and collaborators},
  journal={Working Paper / Under Review},
  year={2026}
}
```

---

## License

This repository is distributed under the MIT License. See [LICENSE](LICENSE) for details.
