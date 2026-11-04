"""Render the mainline-A tables, figures and results.json from analysis JSON.

No model is fitted here: every number is read from
  analysis/authors/*.json            authors' pipeline rerun (run_authors_pipeline.sh)
  analysis/ma_audit.json, *_draws.npz known-truth audit (ma_audit.py)
  analysis/ma_thinning_exact.json    exact thinning expectation (ma_thinning_exact.py)
  analysis/ma_extensions.json        extensions and the primary bootstrap (ma_extensions.py)
  analysis/ma_fe_decomp.json, fe_rev2/  fixed-effects decomposition (ma_fe_decomp.py)
Primary bootstrap: every pooled row (Table S1, centring, fixed effects, within,
sensitivity) is taken from ma_extensions.ratios (B = 2000, rng [SEED, 1]); the
authors' own pooled intervals stay in results.json under pooled/ as
reproduction evidence. The authors' supplementary tables S1-S5 and S8 and
figures 1-4 are rebuilt with their own builders (render_revision.py); S3 is a
true longtable. S6 and S7 are not rendered (withdrawn / design stage removed).
Outputs: tables/*.tex, figures/*.pdf + *.png, tables/manifest.json,
figures/manifest.json and analysis/results.json (estimates array).
"""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import render_revision as rr
from ma_models import utc, write

EXP = Path(__file__).resolve().parents[1]
AN, AUTH, FE = EXP / "analysis", EXP / "analysis/authors", EXP / "analysis/fe_rev2"
TAB, FIG = EXP / "tables", EXP / "figures"
SCEN = {"S0_poisson": "S0 Poisson, no heterogeneity",
        "S1_heterogeneity": "S1 iid heterogeneity, Poisson counts",
        "S2_nb_reporting": "S2 NB offspring + binomial reporting, no heterogeneity",
        "S3_heterogeneity_nb": "S3 iid heterogeneity + NB/reporting",
        "S4_structured": "S4 state + week heterogeneity, Poisson counts",
        "S5_structured_nb": "S5 state + week heterogeneity + NB/reporting",
        "S2P_branching_poisson": "S2P Poisson offspring + binomial reporting, no heterogeneity",
        "S6_reporting_shock": "S6 S3 + state-week reporting shock (half of a)",
        "S7a_scale_dep": "S7a S3 with scale-dependent heterogeneity, gamma 0.25",
        "S7b_scale_dep": "S7b S3 with scale-dependent heterogeneity, gamma 0.5",
        "S8_combo": "S8 Combo S6 shock + S7a scale-dep (gamma 0.25)"}
SHORT = {k: k.split("_")[0] for k in SCEN}
PRIMARY = {"RA": "RA_gi3.6", "G": "G_gi3.6", "RA32": "RA_gi3.2", "RA32source": "RA_gi3.2_sd2.1"}
RHOS = ("0.1", "0.25", "0.5", "1")
QS = ("0.8", "0.6", "0.5", "0.4", "0.3", "0.2", "0.1")
COL = rr.COLORS


def rel(p):
    return str(Path(p).resolve().relative_to(EXP.parents[1]))


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load(p):
    return json.loads(Path(p).read_text())


def same(x, y, what, tol=1e-9):
    if not (np.isfinite(x) and np.isfinite(y) and abs(x - y) <= tol * max(1.0, abs(y))):
        raise rr.RenderError(f"{what}: point {x} differs from the primary bootstrap point {y}")


def summary(r):
    """boot_ci dict -> the authors' ci_summary layout read by render_revision."""
    v, t = r["valid_draws"], r["total_draws"]
    return {"value": r["value"], "ci95": r["ci95"], "finite_fraction": v / t, "nonpositive_fraction": r.get("nonpositive_fraction", 0.0),
            "n_draws": t, "n_finite": v, "interval_scope": "finite_draws", "source": "primary bootstrap"}


class MAInputs(rr.Inputs):
    """Authors' inputs from the rerun with every pooled row replaced by the primary bootstrap."""

    def __init__(self, ext, fe):
        self.results_dir, self.repo_root, self.allow_incomplete = AUTH, EXP, False
        self.data, self.sources, self.issues = {}, {}, []
        for name in ("empirical", "thinning", "secondary", "theory"):
            self.load(name, AUTH / f"{name}.json", required=True)
        self.load("fixed_effects", FE / "L1/fixed_effects.json", required=True)
        for block in (2, 3):
            self.load(f"block{block}", FE / f"block{block}/fixed_effects.json", required=True)
        self.validate()
        for v in self.sources.values():
            v["path"] = rel(v["path"])
        self.authors_pooled = copy.deepcopy(self.data["empirical"]["pooled"])
        self.overlay_primary(ext["ratios"], fe)

    def overlay_primary(self, ratios, fe):
        emp = self.data["empirical"]
        for stat, frame in PRIMARY.items():
            src, dst = ratios[frame], emp["pooled"][stat]
            for m, k in (("a", "a"), ("b", "b"), ("b_over_a", "ratio")):
                same(dst[k]["value"], src[m]["value"], f"empirical {stat} {k}")
                dst[k] = dict(src[m], source="primary bootstrap")
            pooled_centre = emp["centering"].get(stat, {}).get("pooled")
            if pooled_centre:
                for m in ("a", "b"):
                    same(pooled_centre[m]["value"], src[m]["value"], f"centering {stat} pooled {m}")
                    pooled_centre[m] = dict(src[m], source="primary bootstrap")
        p = ratios["RA_gi3.6"]
        for name, L in (("fixed_effects", 1), ("block2", 2), ("block3", 3)):
            est = self.data[name]["estimates"]
            for k, m in (("pool_a", "a"), ("pool_b", "b"), ("pool_ratio", "b_over_a")):
                same(est[k]["value"], p[m]["value"], f"{name} {k}")
                est[k] = summary(p[m])
            for k in ("slope_ratio", "intercept_ratio"):
                est[k] = summary(fe[f"L{L}"][k])


# ---------------------------------------------------------------- formatting
def f(v, d=3):
    return "--" if v is None or not np.isfinite(v) else f"{v:.{d}f}"


def ci(b, d=3):
    return "--" if not b or any(x is None or not np.isfinite(x) for x in b) else f"[{b[0]:.{d}f}, {b[1]:.{d}f}]"


def bc(r, d=3):
    """boot_ci dict -> 'value [lo, hi]'."""
    return f"{f(r.get('value'), d)} {ci(r.get('ci95'), d)}"


def ms(s, d=3):
    """Monte Carlo summ dict -> 'mean [2.5%, 97.5%]'."""
    return f"{f(s.get('mean'), d)} {ci([s.get('q025'), s.get('q975')], d)}"


def mp(s, d=3):
    """Monte Carlo summ dict with a 'point' (ratio of means) -> 'point [2.5%, 97.5%]'."""
    return f"{f(s.get('point'), d)} {ci([s.get('q025'), s.get('q975')], d)}"


def longtable_tex(caption, label, headers, rows, notes=()):
    """Longtable with repeated header; same column scheme and escaping as the authors' table_tex."""
    n = len(headers)
    if any(len(r) != n for r in rows):
        raise rr.RenderError(f"{label}: inconsistent column counts")
    widths = ([.24] + [.70 / (n - 1)] * (n - 1)) if n > 3 else ([.35] + [.60 / (n - 1)] * (n - 1))
    spec = "@{}" + "".join(r">{\raggedright\arraybackslash}p{" + f"{w:.4f}" + r"\linewidth}" for w in widths) + "@{}"
    row = lambda r: " & ".join(rr.tex_escape(c) for c in r) + r" \\"
    head = [r"\toprule", row(headers), r"\midrule"]
    lines = [r"\begingroup\scriptsize\setlength{\tabcolsep}{2pt}", r"\begin{longtable}{" + spec + "}",
             r"\caption{" + rr.tex_escape(caption) + r"}\label{" + label + r"}\\", *head, r"\endfirsthead",
             r"\multicolumn{" + str(n) + r"}{@{}l}{\textit{(continued)}}\\", *head, r"\endhead",
             r"\midrule", r"\multicolumn{" + str(n) + r"}{r@{}}{\textit{(continued on next page)}}\\", r"\endfoot",
             r"\bottomrule", r"\endlastfoot"]
    lines += [row(r) for r in rows] + [r"\end{longtable}"]
    if notes:
        lines += [r"\vspace{-.5\baselineskip}", r"\begin{minipage}{\linewidth}\footnotesize\raggedright",
                  rr.tex_escape(" ".join(map(str, notes))), r"\end{minipage}\par"]
    return "\n".join(lines + [r"\endgroup", ""])


# ---------------------------------------------------------------- results.json
class Estimates:
    def __init__(self):
        self.items = []

    def add(self, name, point, ci95=None, source="", kind="bootstrap percentile 95%", **extra):
        point = None if point is None or not np.isfinite(point) else float(point)
        if ci95 is not None and (len(ci95) != 2 or any(x is None or not np.isfinite(x) for x in ci95)):
            ci95 = None
        ci95 = None if ci95 is None else [float(x) for x in ci95]
        # flat point/lower/upper on every item; lower/upper are null when no interval applies
        self.items.append({"name": name, "point": point, "lower": None if ci95 is None else ci95[0],
                           "upper": None if ci95 is None else ci95[1], "ci95": ci95,
                           "interval": kind if ci95 is not None else None, "source": source, **extra})

    def boot(self, name, r, source, **extra):
        self.add(name, r.get("value"), r.get("ci95"), source, valid_draws=r.get("valid_draws", r.get("n_finite")),
                 total_draws=r.get("total_draws", r.get("n_draws")), **extra)

    def mc(self, name, s, source, point_key="mean", **extra):
        if not s:  # undefined for this generator (e.g. recovery ratio when the truth is 0)
            return
        self.add(name, s.get(point_key), [s.get("q025"), s.get("q975")], source,
                 kind=f"Monte Carlo 2.5-97.5% over synthetic panels (point = {'MC mean' if point_key == 'mean' else 'ratio of MC means'})",
                 n_valid=s.get("n_valid"), n_total=s.get("n_total"), **extra)


SRC = {"audit": "code/ma_audit.py -> analysis/ma_audit.json",
       "thin": "code/ma_thinning_exact.py -> analysis/ma_thinning_exact.json",
       "ext": "code/ma_extensions.py -> analysis/ma_extensions.json",
       "primary": "code/ma_extensions.py section ratios (B=2000, rng [2026110352, 1]) -> analysis/ma_extensions.json",
       "fe": "code/ma_fe_decomp.py -> analysis/ma_fe_decomp.json",
       "emp": "code/revision_empirical.py (run_authors_pipeline.sh) -> analysis/authors/empirical.json",
       "thinA": "code/s03b_thinning.py (run_authors_pipeline.sh) -> analysis/authors/thinning.json"}


# ---------------------------------------------------------------- new tables
def dq(s):
    """Monte Carlo mean retained D: 'median [P10, P90]'."""
    return f"{s['D_median']['mean']:.0f} [{s['D_p10']['mean']:.0f}, {s['D_p90']['mean']:.0f}]"


def t_audit(au, ext, fe):
    obs, sc = au["observed"], au["scenarios"]
    p = ext["ratios"]["RA_gi3.6"]
    w = ext["within"]
    h1 = ["Generator", "n retained", "D median [P10, P90]", "True rel. variance (retained rows)", "Intercept a", "a / truth", "True-error intercept", "Week-centred a reduction"]
    rows1 = [["Observed (R_A, GI 3.6)", str(obs["n"]), f"{obs['D_median']:.0f} [{obs['D_p10']:.0f}, {obs['D_p90']:.0f}]", "unknown",
              bc(p["a"], 4), "--", "--", f(1 - obs["a_week"] / obs["a"], 2)]]
    for k, s in sc.items():
        tv = s["true_total_kept_rows"]["mean"]  # S6: transmission plus reporting-shock variance
        at = s["a_over_true_total"]
        rows1.append([SCEN[k], f(s["n_retained"]["mean"], 0), dq(s), f(tv, 4) if tv > 0 else "0", ms(s["a"], 4),
                      mp(at, 2) if at.get("point") is not None else "truth 0", ms(s["err_a"], 5), ms(s["week_reduction"], 2)])
    lam = [s["lambda"] for s in sc.values()]
    p90 = [s["D_p90"]["mean"] for s in sc.values()]
    rs = sc["S6_reporting_shock"]["reporting_shock"]
    notes1 = [f"{au['meta']['reps']} synthetic panels per generator on the real {obs['n']}-row R_A design; analyst scale D = simulated C_t and the filter C_t >= 5 applied to each panel; n retained and the D quantiles are Monte Carlo means over panels.",
              f"Count alignment: each run of consecutive weeks is level-anchored to the observed counts of that run, and a global scale per generator ({min(lam):.3f}-{max(lam):.3f}) sets the mean retained median D to 45 (observed 45). The upper tail is not matched: the generator P90 of D is {min(p90):.0f}-{max(p90):.0f} against the observed {obs['D_p90']:.0f}.",
              "Intervals are Monte Carlo 2.5-97.5% ranges across panels; a / truth is the ratio of panel means with the range of panel ratios. Truth is the relative variance of the target R_i on retained rows; for S6 it also includes the relative variance of the reporting-shock term. True-error intercept: affine fit of ((x - R_i)/xbar)^2. Week-centred reduction = 1 - a_week/a. Observed a: primary bootstrap (2,000 draws).",
              f"S6 (shock log-SD {rs['s']:.3f}): a / v_R = {mp(rs['a_over_vR'], 2)}; true-error intercept / reporting-shock variance = {mp(rs['err_a_over_vH'], 2)}; lag-1 correlation of the shock term between consecutive rows {ms(rs['lag1_corr'], 2)}."]
    h2 = ["Generator", "b/a", "b/b_P", "b_w/b_P (state FE)", "Pooled / within slope", "Thinning c/b_P", "In-sample contrast", "Quantile of observed contrast"]
    athin = load(AUTH / "thinning.json")["estimates"]
    inv = fe["L1"]["slope_ratio"]  # pooled / state-FE slope from the fixed-effects decomposition (primary weights)
    rows2 = [["Observed", bc(p["b_over_a"], 1), bc(p["b_over_bP"]), bc(w["over_bP"]["state_authors"]), bc(inv, 2), bc(athin["full/full/c_over_bP"]),
              f(obs["contrast_affine_minus_power"], 4), "--"]]
    for k, s in sc.items():
        defined = s["a_over_true_total"].get("point") is not None
        rows2.append([SCEN[k], ms(s["b_over_a"], 1) if defined else "---", ms(s["b_over_bP"]), ms(s["bw_over_bP"]),
                      ms(s["pooled_over_within_slope"], 2), ms(s["thinning"]["c_over_bP"]), ms(s["contrast"], 4), f(s["observed_contrast_quantile"], 3)])
    undefined = ", ".join(SHORT[k] for k, s in sc.items() if s["a_over_true_total"].get("point") is None)
    notes2 = ["b_P is the independent-Poisson reference computed on each panel. S2, S3, S5, S6 and S7 use negative-binomial offspring (k = 0.3) with binomial reporting probability 0.25; S2P uses Poisson offspring with the same reporting.",
              f"b/a undefined when a ≈ 0 ({undefined}: true relative variance 0), shown as ---. Pooled / within slope = (b/b_P)/(b_w/b_P) with state fixed effects. Observed row: primary bootstrap (2,000 draws); its pooled / within slope is the slope ratio of the fixed-effects decomposition.",
              f"Thinning uses {au['meta']['thin_outer']} panels x {au['meta']['thin_inner']} shared-count binomial draws per q and the authors' least-squares rate through the origin.",
              "Contrast = in-sample affine minus power mean score (authors' revision_design.contrast); its last column is the share of synthetic panels with a contrast at or below the observed value."]
    return (("tab_ma_audit_loss", "Known-truth audit on the real R_A design: retained counts, intercept, true error and recovery by generator", "tab:ma_audit_loss", h1, rows1, notes1),
            ("tab_ma_audit_slopes", "Known-truth audit: b/a, slope references, thinning and the in-sample contrast by generator", "tab:ma_audit_slopes", h2, rows2, notes2))


def t_baselines(au, ext):
    b, p = au["baseline"], ext["ratios"]["RA_gi3.6"]
    h = ["Reference", "Analytic reference / b_P", "Simulated b/b_P (no heterogeneity)", "Observed b / reference"]
    rows = [[f"Branching-Poisson, rho = {r}", f(b["analytic"][f"rho{r}"]["over_bP"]), ms(b["simulation"][f"rho{r}"]["b_over_bP"]),
             bc(p[f"b_over_bPbr_rho{r}"], 2)] for r in RHOS]
    rows.append(["Independent-Poisson b_P", "1", ms(au["scenarios"]["S0_poisson"]["b_over_bP"]), bc(p["b_over_bP"], 2)])
    notes = [f"The branching-Poisson reference is b_Pbr(rho) = mean over the observed retained rows of kappa^2 [(1 - rho) + (C_t + 0.5)/(C_t+1 + 0.5)]; independent-Poisson reference b_P = {b['independent_poisson_bP']:.4f}.",
             f"Simulated: pure branching chain with Poisson offspring and binomial reporting rho (independent-Poisson counts, S0, in the last row), no transmission heterogeneity, count-aligned as in the audit; {au['meta']['baseline_sim_reps']} panels each, Monte Carlo 2.5-97.5% in brackets.",
             f"Observed: primary bootstrap ({p['b']['total_draws']:,} draws), b_Pbr recomputed with the weights of each draw. The ratio of the observed slope to the reference depends on the reporting probability."]
    return ("tab_ma_baselines", "Branching-Poisson and independent-Poisson slope references: analytic, simulated and observed", "tab:ma_baselines", h, rows, notes)


def t_rho_k(au, ext):
    ks = [0.1, 0.3, 1.0, 3.0, None]
    grid = {(g["rho"], g["k"]): g for g in au["rho_k_grid"]}
    rhos = sorted({g["rho"] for g in au["rho_k_grid"]})
    an = au["baseline"]["analytic"]
    h = ["Reporting probability rho", "Branching r(rho)"] + [f"k = {k:g}" if k else "Poisson offspring" for k in ks]
    rows = [[f"{r:g}", f(an[f"rho{r:g}"]["over_bP"])] + [ms(grid[(r, k)]["b_over_bP"], 2) for k in ks] for r in rhos]
    o = ext["ratios"]["RA_gi3.6"]["b_over_bP"]
    notes = [f"Unthinned b/b_P, with b_P the independent-Poisson reference, under S1 iid heterogeneity (true relative variance = observed a) with negative-binomial offspring of dispersion k and binomial reporting rho; {au['meta']['grid_reps']} panels per cell, Monte Carlo 2.5-97.5% in brackets.",
             "r(rho) = b_Pbr(rho)/b_P is the analytic branching-Poisson reference on the observed rows (no heterogeneity).",
             f"Observed b/b_P = {bc(o, 2)} (primary bootstrap, {o['valid_draws']:,} draws)."]
    return ("tab_ma_rho_k", "Unthinned slope relative to the independent-Poisson reference over reporting probability and offspring dispersion", "tab:ma_rho_k", h, rows, notes)


def t_thinning(th):
    h = ["q", "Exact rate", "Exact rate / b_P", "Authors' MC rate [95% CI]", "Leverage in c", "Positive-part mass (max)"]
    mc = {f"{m['q']:.1f}": m for m in th["authors_mc"]}
    rows = []
    for r in th["per_q"]:
        q = f"{r['q']:.1f}"
        rows.append([q, f(r["rate_exact"], 4), f(r["rate_exact_over_bP"], 3), f"{f(mc[q]['rate_mc'], 4)} {ci(mc[q]['ci95'], 4)}",
                     f(th["lsq_leverage_by_q"][q], 3), f"{r['positive_part_mass_max']:.1e}"])
    athin = load(AUTH / "thinning.json")["estimates"]
    rows.append(["c, all q", f(th["c_exact_all_q"], 4), f(th["c_exact_all_q_over_bP"], 3), bc(athin["full/full/c"], 4), "1.000", "--"])
    rows.append(["c, without q = 0.1", f(th["c_exact_without_q01"], 4), f(th["c_exact_without_q01_over_bP"], 3), "--", "--", "--"])
    chk = "; ".join(f"q = {c['q']}: independent-thinning MC {c['independent_mc_rate']:.4f} (SE {c['mc_se']:.4f}) vs expected-loss {c['expected_loss_rate']:.4f}" for c in th["mc_check"])
    notes = [f"Rate = (b_q - b_1)/(1/q - 1) from the affine FGLS fit of the exact expected thinned loss (binomial supports summed, mass below 1e-13 dropped); b_P = {th['bP_delta']:.4f} is the delta-method prediction.",
             "Leverage: share of the least-squares rate through the origin carried by each q.", "Check: " + chk + "."]
    return ("tab_ma_thinning_exact", "Exact expected thinning increment against the delta prediction b_P", "tab:ma_thinning_exact", h, rows, notes)


def t_varlog(th):
    Cs = sorted({v["C"] for v in th["var_log_small_count"]})
    h = ["C"] + [f"q = {q}" for q in QS]
    look = {(v["C"], f"{v['q']:.1f}"): v["exact_added_over_delta"] for v in th["var_log_small_count"]}
    rows = [[str(C)] + [f(look[(C, q)], 2) for q in QS] for C in Cs]
    notes = ["Exact added variance Var log(Bin(C, q) + 0.5) - Var log(C + 0.5) divided by the delta value (1 - q)/(qC). Values below 1 mean that thinning adds less variance than the delta method assumes."]
    return ("tab_ma_varlog", "Small-count variance of the log thinned count relative to the delta approximation", "tab:ma_varlog", h, rows, notes)


def t_targets(ext):
    t = ext["table1"]
    p = t["parameters"]
    h = ["m", "Known denominator A/m", "MC (SE)", "Observed denominator (A + (1 - rho)(Rbar^2 + v_R))/m", "MC (SE)", "MC observed / A/m"]
    rows = [[f"{r['m']:g}", f(r["known_denominator_analytic"], 5), f"{r['known_denominator_mc']:.5f} ({r['known_denominator_mc_se']:.5f})",
             f(r["observed_denominator_delta"], 5), f"{r['observed_denominator_mc']:.5f} ({r['observed_denominator_mc_se']:.5f})", f(r["observed_over_known"], 3)]
            for r in t["rows"]]
    notes = [f"Parameters Rbar = {p['Rbar']}, v_R = {p['v_R']}, k = {p['k']}, rho = {p['rho']} (A = {p['A']:.4f}); {p['reps']:,} Monte Carlo draws per m.",
             f"Crossing scale m_x = A/v_R = {t['m_cross_known']:.2f} with a known denominator and (A + (1 - rho)(Rbar^2 + v_R))/v_R = {t['m_cross_observed_delta']:.2f} when the denominator is the observed previous count (delta approximation).",
             "At m = 10 the delta approximation understates the Monte Carlo error of the observed-denominator estimator."]
    return ("tab_ma_targets_observed", "Point-target error with a known versus an observed denominator", "tab:ma_targets_observed", h, rows, notes)


def t_gi(ext):
    r = ext["ratios"]
    h = ["Generation interval", "n", "a", "b", "b/a", "b_P", "b/b_P", "Paired change in b/a", "Paired change in b/b_P"]
    rows = []
    for name, lab in (("RA_gi3.6", "Gamma 3.6 / 1.6 (primary)"), ("RA_gi3.2", "Gamma 3.2 / 1.42"), ("RA_gi3.2_sd2.1", "Gamma 3.2 / 2.1")):
        s = r[name]
        d1 = r.get(f"{name}_minus_gi3.6_b_over_a")
        d2 = r.get(f"{name}_minus_gi3.6_b_over_bP")
        rows.append([lab, str(s["n"]), bc(s["a"], 4), bc(s["b"], 3), bc(s["b_over_a"], 1), bc(s["bP"], 3), bc(s["b_over_bP"], 2),
                     bc(d1, 2) if d1 else "reference", bc(d2, 3) if d2 else "reference"])
    notes = [f"State-week bootstrap with {r['RA_gi3.6']['b']['total_draws']} draws; the same weights are reused for every generation interval (paired = {r['paired_weights']}), so the last two columns are paired differences from the primary specification.",
             "Gamma 3.2 / 1.42 keeps the primary coefficient of variation; 3.2 / 2.1 is the intrinsic generation-time mean and SD estimated from US household data (Chan et al., 2025)."]
    return ("tab_ma_gi", "Generation-interval sensitivity of the pooled loss model (bootstrap 95% confidence interval in brackets)", "tab:ma_gi", h, rows, notes)




def t_validation(ext):
    v = ext["validation"]
    h = ["Held-out unit", "Folds", "Mean power - affine score", "Paired fold SE", "z", "Folds favouring affine", "Bootstrap 95% CI", "Bootstrap share favouring affine"]
    rows = [[g.capitalize(), str(v[g]["folds"]), f(v[g]["mean_power_minus_affine"], 4), f(v[g]["paired_fold_se"], 4), f(v[g]["z"], 2),
             f"{v[g]['affine_better_folds']}/{v[g]['folds']}", ci(v[g]["bootstrap"]["ci95"], 4), f(v[g]["bootstrap_share_affine_better"], 2)] for g in ("state", "season")]
    notes = ["Score: weighted mean of L/mu + log mu on held-out rows; positive differences favour the affine model.",
             f"Paired fold SE from the between-fold variation of the paired differences; state-week bootstrap with {v['state']['bootstrap']['total_draws']} draws refits both models in every fold."]
    return ("tab_ma_validation", "Held-out comparison of the affine and power loss models with paired uncertainty", "tab:ma_validation", h, rows, notes)




def t_within(ext):
    w = ext["within"]
    pf = w["point_fits"]
    lab = {"pooled_b": "Pooled affine slope b", "bP": "b_P, independent-Poisson reference", "state_authors": "State FE, authors' iterated FGLS (b_w)",
           "state_onestep": "State FE, one-step WLS", "week_onestep": "Week FE, one-step WLS", "two_way_onestep": "State + week FE, one-step WLS",
           "two_way_iterated": "State + week FE, damped iterated FGLS"}
    h = ["Estimator", "Slope [95% CI]", "Slope / b_P", "Slope / pooled b", "Fit status"]
    rows = []
    for k, L in lab.items():
        st = "--"
        if isinstance(pf.get(k), dict):
            p = pf[k]
            st = f"{'converged' if p['converged'] else 'not converged'}, {p['iterations']} it." if p.get("iterations", 1) > 1 else "single step"
        rows.append([L, bc(w["slopes"][k]), bc(w["over_bP"][k], 2) if k in w["over_bP"] else "--",
                     bc(w["over_pooled_b"][k], 2) if k in w["over_pooled_b"] else "--", st])
    u = pf["two_way_iterated_undamped"]
    notes = [f"Primary bootstrap: state-week weights of the primary stream, {w['slopes']['pooled_b']['total_draws']:,} draws, so the pooled row equals Table S1; draws whose iterated fit does not converge are dropped ({w['slopes']['two_way_iterated']['valid_draws']} valid for the damped two-way slope).",
             f"Undamped two-way iteration: {'converged' if u['converged'] else 'did not converge'} after {u['iterations']} iterations (last values {', '.join(f'{x:.4f}' for x in u['history'][-2:])}); damping 0.5 removes the two-cycle.",
             "One-step WLS uses the pooled affine working precision; the iterated estimators re-weight with the fixed-effect fitted means."]
    return ("tab_ma_within", "Within-state, within-week and two-way slopes against the independent-Poisson reference (bootstrap 95% confidence interval in brackets)", "tab:ma_within", h, rows, notes)


def t_sens(ext):
    s = ext["sensitivity"]["specs"]
    lab = {"primary": "Primary (C_t, D_t >= 5; cc 0.5; growth rule 1.05)", "min_count_10": "Counts >= 10", "min_count_20": "Counts >= 20",
           "min_count_50": "Counts >= 50", "continuity_1.0": "Continuity correction 1.0", "window_threshold_1.10": "Growth-window rule 1.10"}
    h = ["Specification", "n", "a", "b", "b/a", "b/b_P"]
    rows = []
    for k, r in s.items():
        L = lab.get(k, "Season " + k.replace("season_", "") if k.startswith("season_") else k)
        ra = r["b_over_a"]
        flag = ""
        if (ra.get("nonpositive_fraction") or 0) > .025:
            flag = f" (b/a <= 0 in {100 * ra['nonpositive_fraction']:.1f}% of draws)"
        elif ra["valid_draws"] < ra["total_draws"]:
            flag = f" ({ra['total_draws'] - ra['valid_draws']} draws not converged)"
        rows.append([L, str(r["n"]), bc(r["a"], 4), bc(r["b"], 3), bc(r["b_over_a"], 1) + flag, bc(r["b_over_bP"], 2)])
    B = s["primary"]["b"]["total_draws"]
    Bo = next(r["b"]["total_draws"] for k, r in s.items() if k != "primary")
    notes = [f"Each row refits the authors' affine FGLS on the subset or variant with a state-week bootstrap; the primary row is the primary bootstrap ({B:,} draws, identical to Table S1), the other rows use {Bo:,} draws each.",
             "Where b/a <= 0 (b <= 0) in more than 2.5% of draws its percentile interval is unstable; draws whose FGLS fit did not converge are dropped and counted."]
    return ("tab_ma_sensitivity", "Sensitivity of the pooled loss model to filters, continuity correction, window rule and season (bootstrap 95% confidence interval in brackets)", "tab:ma_sensitivity", h, rows, notes)


def t_seasons(ext):
    s = ext["seasons"]
    p, c, t = s["pooled_loss_season_intercepts"], s["season_centred_common_slope"], s["test"]
    h = ["Quantity", "Pooled loss, season intercepts", "Season-centred loss, season intercepts"]
    rows = [["Common slope b", bc(p["b_common"]), bc(c["b_common"])]]
    for se in s["seasons"]:
        rows.append([f"a, {se}", bc(p[f"a_{se}"], 4), bc(c[f"a_{se}"], 4)])
    for se in s["seasons"]:
        rows.append([f"b/a, {se}", bc(p[f"b_over_a_{se}"], 1), bc(c[f"b_over_a_{se}"], 1)])
    wl, q = t["wald_log_a"], t["cochran_q_log_a"]
    wr, ws = t["wald_raw_a"], t["wald_sep_log_ba"]
    rows += [[f"Wald test of equal log a_s (df {wl['df']})", f"W = {wl['stat']:.2f}, p = {wl['p']:.3f}", "--"],
             [f"Cochran Q, log a_s (df {q['df']})", f"Q = {q['Q']:.2f}, p = {q['p']:.3f}, I2 = {100 * q['I2']:.0f}%", "--"],
             [f"Wald test of equal a_s (df {wr['df']})", f"W = {wr['stat']:.2f}, p = {wr['p']:.3f}", "--"],
             [f"Wald test, separate season fits, log b_s/a_s (df {ws['df']})", f"W = {ws['stat']:.2f}, p = {ws['p']:.3f}", "--"],
             ["Minimum detectable highest/lowest a (80% power)", f"{t['mdd_max_over_min_a']:.2f} (observed {t['observed_max_over_min_a']:.2f})", "--"]]
    B = p["b_common"]["total_draws"]
    d = s["descriptive"]["max_over_min_b_over_a"]
    notes = [f"Affine loss with one intercept per season and a common slope (authors' group FGLS); state-week bootstrap, {B:,} draws ({p['b_common']['valid_draws']:,} valid). The season-centred column centres x on season means before forming the loss.",
             f"Tests use the bootstrap covariance of the season parameters ({t['n_draws_used']:,} of {t['n_draws_total']:,} draws with all parameters finite); with a common slope, equal log a_s is equivalent to equal season b/a. Cochran Q uses inverse-variance weights. The separate-fit test refits the loss in each season on its own bootstrap stream.",
             f"Minimum detectable difference: noncentrality {t['noncentrality_power80']:.2f} gives 80% power for a level-0.05 chi-square({wl['df']}) test; scaling the observed log a_s deviations by f = {t['mdd_factor']:.2f} reaches it, which corresponds to the highest/lowest a shown.",
             f"The highest to lowest season b/a ratio ({bc(d, 2)}) is descriptive and is not a homogeneity test."]
    return ("tab_ma_seasons", "Season-stratified loss model with a common slope and tests of equal season intercepts (bootstrap 95% confidence interval in brackets)", "tab:ma_seasons", h, rows, notes)


# ---------------------------------------------------------------- figures
def save(fig, name, manifest, deps, caption):
    for ext in ("pdf", "png"):
        p = FIG / f"{name}.{ext}"
        if ext == "pdf":
            fig.savefig(p, format="pdf", bbox_inches="tight", metadata={"Creator": "ma_render.py", "CreationDate": None, "ModDate": None})
        else:
            fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    manifest.append({"name": name, "files": [f"figures/{name}.pdf", f"figures/{name}.png"], "inputs": deps, "caption_hint": caption,
                     "sha256_pdf": sha(FIG / f"{name}.pdf")})


def fig_audit(au, ext, draws):
    sc = list(au["scenarios"])
    obs = au["observed"]
    prim = ext["ratios"]["RA_gi3.6"]
    athin = load(AUTH / "thinning.json")["estimates"]
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.4))
    y = np.arange(len(sc))[::-1]
    lab = [SHORT[k] for k in sc]

    def panel(ax, key, obs_rec, letter, title, xlabel, sub=None):
        for yy, k in zip(y, sc):
            s = au["scenarios"][k][key] if sub is None else au["scenarios"][k][sub][key]
            ax.plot([s["q025"], s["q975"]], [yy, yy], color=COL[0], lw=1)
            ax.plot(s["mean"], yy, "o", color=COL[0])
        ax.axvline(obs_rec["value"], color=COL[1], lw=1, label="observed")
        ax.axvspan(*obs_rec["ci95"], color=COL[1], alpha=.12, lw=0)
        ax.set_yticks(y, lab)
        ax.set_xlabel(xlabel)
        rr.panel_title(ax, letter, title)

    panel(axs[0, 0], "b_over_bP", prim["b_over_bP"], "A", "Unthinned slope", r"$b/b_{\mathrm{P}}$")
    axs[0, 0].axvline(1, color="#777777", lw=.7, ls=":")
    panel(axs[0, 1], "c_over_bP", athin["full/full/c_over_bP"], "B", "Thinning rate", r"$c/b_{\mathrm{P}}$", sub="thinning")
    axs[0, 1].axvline(1, color="#777777", lw=.7, ls=":")
    ax = axs[1, 0]
    for yy, k in zip(y, sc):
        s = au["scenarios"][k]["a"]
        ax.plot([s["q025"], s["q975"]], [yy, yy], color=COL[0], lw=1)
        ax.plot(s["mean"], yy, "o", color=COL[0], label=r"intercept $a$" if k == sc[0] else None)
        ax.plot(au["scenarios"][k]["true_total_kept_rows"]["mean"], yy, "D", mfc="none", color=COL[2],
                label="true relative variance (retained rows)" if k == sc[0] else None)
    ax.axvline(prim["a"]["value"], color=COL[1], lw=1, label=r"observed $a$")
    ax.axvspan(*prim["a"]["ci95"], color=COL[1], alpha=.12, lw=0)
    ax.set_yticks(y, lab)
    ax.set_xlabel(r"Intercept $a$")
    ax.legend(frameon=False, loc="upper right", fontsize=6)
    rr.panel_title(ax, "C", "Intercept and true variance")
    ax = axs[1, 1]
    data = [draws[f"{k}__contrast"] for k in sc]
    ax.boxplot(data, vert=False, positions=y, widths=.55, showfliers=False, medianprops={"color": COL[0]})
    ax.axvline(obs["contrast_affine_minus_power"], color=COL[1], lw=1)
    for yy, k in zip(y, sc):
        ax.text(1.0, yy, f"{au['scenarios'][k]['observed_contrast_quantile']:.3f}", transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=6.5)
    ax.set_yticks(y, lab)
    ax.set_xlabel(r"In-sample contrast (affine $-$ power)")
    rr.panel_title(ax, "D", "Score contrast; quantile of observed at right")
    fig.tight_layout()
    return fig


def fig_rho_k(au, ext):
    ks = [0.1, 0.3, 1.0, 3.0, None]
    rhos = sorted({g["rho"] for g in au["rho_k_grid"]})
    grid = {(g["rho"], g["k"]): g["b_over_bP"] for g in au["rho_k_grid"]}
    o = ext["ratios"]["RA_gi3.6"]["b_over_bP"]
    b = au["baseline"]
    cmap = plt.get_cmap("viridis")
    fig, ax = plt.subplots(figsize=(5.0, 3.4))
    ax.axhspan(*o["ci95"], color=COL[1], alpha=.12, lw=0)
    ax.axhline(o["value"], color=COL[1], lw=1, label=rf"observed $b/b_{{\mathrm{{P}}}}$ = {o['value']:.2f} [{o['ci95'][0]:.2f}, {o['ci95'][1]:.2f}]")
    off = np.linspace(-.012, .012, len(ks))
    for j, k in enumerate(ks):
        c = cmap(j / (len(ks) - 1) * .9)
        x = np.array(rhos) + off[j]
        m = [grid[(r, k)] for r in rhos]
        ax.plot(x, [s["mean"] for s in m], "-o", color=c, ms=3, lw=1, label=(f"heterogeneity, $k$ = {k:g}" if k else "heterogeneity, Poisson offspring"))
        for xx, s in zip(x, m):
            ax.plot([xx, xx], [s["q025"], s["q975"]], color=c, lw=.8)
    rg = np.array(rhos)
    ax.plot(rg, [b["analytic"][f"rho{r:g}"]["over_bP"] for r in rhos], "--", color="black", lw=1, label=r"branching-Poisson reference $r(\rho)$")
    sim = [b["simulation"][f"rho{r:g}"]["b_over_bP"] for r in rhos]
    ax.errorbar(rg + .02, [s["mean"] for s in sim], yerr=[[s["mean"] - s["q025"] for s in sim], [s["q975"] - s["mean"] for s in sim]],
                fmt="s", mfc="white", color="black", ms=3.5, lw=.8, capsize=0, label="simulated branching chain, no heterogeneity")
    ax.axhline(1, color="#777777", lw=.7, ls=":")
    ax.set_xlabel(r"Reporting probability $\rho$")
    ax.set_ylabel(r"$b/b_{\mathrm{P}}$")
    ax.set_xticks(rhos, [f"{r:g}" for r in rhos])
    ax.legend(frameon=False, fontsize=5.8, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.tight_layout()
    return fig


def fig_thinning(th, au):
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.6))
    ax = axs[0]
    q = np.array([r["q"] for r in th["per_q"]])
    ax.plot(q, [r["rate_exact"] for r in th["per_q"]], "-o", color=COL[0], label="exact expectation")
    for m in th["authors_mc"]:
        ax.plot([m["q"], m["q"]], m["ci95"], color=COL[1], lw=1)
    ax.plot([m["q"] for m in th["authors_mc"]], [m["rate_mc"] for m in th["authors_mc"]], "s", color=COL[1], ms=3, label="Monte Carlo")
    ax.axhline(th["bP_delta"], color="#555555", ls="--", lw=.8, label=r"delta prediction $b_{\mathrm{P}}$")
    ax.set_xlabel(r"Retention probability $q$")
    ax.set_ylabel("Slope increment rate")
    ax.legend(frameon=False, loc="lower right", fontsize=6)
    rr.panel_title(ax, "A", r"Increment rate by $q$")
    ax = axs[1]
    lev = th["lsq_leverage_by_q"]
    ax.bar(range(len(QS)), [lev[k] for k in QS], color=COL[0])
    ax.set_xticks(range(len(QS)), QS)
    ax.set_xlabel(r"$q$")
    ax.set_ylabel(r"Leverage in pooled $c$")
    rr.panel_title(ax, "B", r"Weight of each $q$")
    ax = axs[2]
    vl = th["var_log_small_count"]
    for j, qq in enumerate((0.5, 0.3, 0.1)):
        pts = sorted((v["C"], v["exact_added_over_delta"]) for v in vl if abs(v["q"] - qq) < 1e-9)
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "-o", color=COL[j], label=rf"$q={qq}$")
    ax.axhline(1, color="#555555", ls="--", lw=.8)
    ax.set_xscale("log")
    ax.set_xlabel(r"Count $C$")
    ax.set_ylabel("Exact / delta added variance")
    ax.legend(frameon=False, fontsize=6)
    rr.panel_title(ax, "C", "Small-count variance")
    fig.tight_layout()
    return fig




def fig_sensitivity(ext):
    s = ext["sensitivity"]["specs"]
    w = ext["within"]
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.2))
    ax = axs[0]
    names = list(s)
    y = np.arange(len(names))[::-1]
    lab = {"primary": "primary", "min_count_10": r"$C_t\geq10$", "min_count_20": r"$C_t\geq20$", "min_count_50": r"$C_t\geq50$",
           "continuity_1.0": "correction 1.0", "window_threshold_1.10": "window 1.10"}
    for yy, k in zip(y, names):
        r = s[k]["b_over_a"]
        lo, hi = r["ci95"]
        unstable = (r.get("nonpositive_fraction") or 0) > .025
        c = COL[3] if unstable else COL[0]
        ax.plot([max(lo, 2), hi], [yy, yy], color=c, lw=1)
        if lo < 2:
            ax.plot(2, yy, "<", color=c, ms=4)
        ax.plot(r["value"], yy, "o", color=c)
    ax.set_xlim(left=1.6)
    ax.axvline(s["primary"]["b_over_a"]["value"], color=COL[1], lw=.8, ls="--")
    ax.set_xscale("log")
    ax.set_yticks(y, [lab.get(k, k.replace("season_", "season ")) for k in names])
    ax.set_xlabel(r"$b/a$ (log scale)")
    rr.panel_title(ax, "A", r"$b/a$ across specifications")
    ax = axs[1]
    keys = ["pooled_b", "state_authors", "state_onestep", "week_onestep", "two_way_onestep", "two_way_iterated"]
    labs = ["pooled", "state FE iterated", "state FE one-step", "week FE one-step", "two-way one-step", "two-way iterated"]
    y = np.arange(len(keys))[::-1]
    rr.points_ci(ax, y, [w["over_bP"][k] for k in keys], horizontal=True)
    ax.axvline(1, color="#555555", ls="--", lw=.8)
    ax.set_yticks(y, labs)
    ax.set_xlabel(r"slope / $b_{\mathrm{P}}$")
    rr.panel_title(ax, "B", "Fixed-effect slopes")
    fig.tight_layout()
    return fig




# ---------------------------------------------------------------- main
def main():
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    au, th, ext, fe = (load(AN / n) for n in ("ma_audit.json", "ma_thinning_exact.json", "ma_extensions.json", "ma_fe_decomp.json"))
    if "calibration" in ext:
        raise rr.RenderError("ma_extensions.json still carries the withdrawn calibration section")
    i = MAInputs(ext, fe)
    draws = np.load(AN / "ma_audit_draws.npz")
    inputs = {"authors_and_fe": i.sources, **{n: {"path": rel(AN / n), "sha256": sha(AN / n)} for n in
                                              ("ma_audit.json", "ma_audit_draws.npz", "ma_thinning_exact.json", "ma_extensions.json", "ma_fe_decomp.json")}}
    tmani = []
    # authors' supplementary tables (S6 calibration withdrawn, S7 design stage dropped)
    for index, (name, builder) in enumerate(zip(rr.TABLES, rr.TABLE_BUILDERS), 1):
        if name in ("design", "calibration"):
            continue
        caption, headers, rows, notes, long = builder(i, "en")
        fname = f"tab_s{index:02d}_{name}.tex"
        as_long = name == "thinning" or long
        tex = longtable_tex(caption, "tab:" + name, headers, rows, notes) if as_long else rr.table_tex(caption, "tab:" + name, headers, rows, notes)
        (TAB / fname).write_text(tex, encoding="utf-8")
        tmani.append({"file": f"tables/{fname}", "label": "tab:" + name, "rows": len(rows), "environment": "longtable" if as_long else "tabularx",
                      "builder": f"render_revision.{builder.__name__}", "origin": f"authors' Table S{index}",
                      **({"pooled_rows": "primary bootstrap"} if name in ("pooled", "fe", "centering") else {})})
    for stale in ("tab_s06_calibration.tex", "tab_s07_design.tex", "tab_ma_calibration_diag.tex"):
        if (TAB / stale).exists():
            (TAB / stale).unlink()
    caption, headers, rows, notes, _ = rr.summer_table(i, "en")
    (TAB / "tab_ma_summer.tex").write_text(rr.table_tex(caption, "tab:summer", headers, rows, notes), encoding="utf-8")
    tmani.append({"file": "tables/tab_ma_summer.tex", "label": "tab:summer", "rows": len(rows), "environment": "tabularx",
                  "builder": "render_revision.summer_table", "origin": "authors' Table S8 summer rows (descriptive, separated)"})
    builders = [*t_audit(au, ext, fe), t_baselines(au, ext), t_rho_k(au, ext), t_thinning(th), t_varlog(th), t_targets(ext), t_gi(ext),
                t_within(ext), t_sens(ext), t_seasons(ext), t_validation(ext)]
    for fname, caption, label, h, rows, notes in builders:
        long = len(rows) > 12 or len(h) > 7
        tex = longtable_tex(caption, label, h, rows, notes) if long else rr.table_tex(caption, label, h, rows, notes)
        (TAB / f"{fname}.tex").write_text(tex, encoding="utf-8")
        tmani.append({"file": f"tables/{fname}.tex", "label": label, "rows": len(rows), "environment": "longtable" if long else "tabularx",
                      "builder": "ma_render.py", "origin": "mainline-A analysis"})
    write(TAB / "manifest.json", {"generated_utc": utc(), "generator": "code/ma_render.py", "inputs": inputs,
                                  "latex_packages": ["booktabs", "array", "tabularx", "longtable"],
                                  "dropped": {"tab_s07_design": "design stage removed from the manuscript; observed-contrast quantiles under the audit generators replace it (tab_ma_audit_slopes)",
                                              "tab_s06_calibration": "calibration withdrawn (2026-11-03): its simulator did not reproduce the observed slope, row count or autocorrelation; see code/README.md",
                                              "tab_ma_calibration_diag": "calibration withdrawn (2026-11-03)"},
                                  "tables": tmani})
    # figures
    fmani = []
    rc, fonts = rr.font_config("en")
    with matplotlib.rc_context(rc):
        caps = {"fig_01_targets": "authors' Figure 1 (targets)", "fig_02_models": "authors' Figure 2 (loss models)",
                "fig_03_thinning": "authors' Figure 3 (thinning)", "fig_04_heterogeneity": "authors' Figure 4 (heterogeneity)"}
        for name, builder, deps in zip(rr.FIGURES, rr.FIGURE_BUILDERS, rr.FIGURE_INPUTS):
            save(builder(i, "en", plt, np), name, fmani, [i.sources[d]["path"] for d in deps if d in i.sources], caps[name])
        save(fig_audit(au, ext, draws), "fig_05_audit", fmani, ["analysis/ma_audit.json", "analysis/ma_audit_draws.npz", "analysis/ma_extensions.json"],
             "known-truth audit over ten generators: b/b_P, c/b_P, intercept vs truth, in-sample contrast")
        save(fig_thinning(th, au), "fig_06_thinning_exact", fmani, ["analysis/ma_thinning_exact.json"],
             "exact thinning expectation vs authors' MC and delta b_P; leverage by q; small-count variance")
        save(fig_rho_k(au, ext), "fig_07_rho_k_grid", fmani, ["analysis/ma_audit.json", "analysis/ma_extensions.json"],
             "b/b_P over rho by k with the branching-Poisson reference r(rho) and the observed interval")
        save(fig_sensitivity(ext), "fig_08_sensitivity", fmani, ["analysis/ma_extensions.json"],
             "b/a across specifications; fixed-effect slopes relative to b_P")
    write(FIG / "manifest.json", {"generated_utc": utc(), "generator": "code/ma_render.py", "fonts": fonts, "figures": fmani})
    write(AN / "results.json", {"generated_utc": utc(), "generator": "code/ma_render.py", "inputs": inputs,
                                "primary_bootstrap": ext["ratios"].get("primary"), "estimates": estimates(i, au, th, ext, fe)})
    print("tables", len(tmani), "figures", len(fmani))


def estimates(i, au, th, ext, fe):
    E = Estimates()
    for k, n in (("a", "a"), ("b", "b"), ("ratio", "b_over_a")):
        E.boot(f"pooled/RA/{n}", i.authors_pooled["RA"][k], SRC["emp"], role="authors' pooled interval, reproduction record (not used in tables)")
    r = ext["ratios"]
    for stat, frame in (("RA", "RA_gi3.6"), ("G", "G_gi3.6")):
        for m in ("a", "b", "b_over_a", "bP", "b_over_bP") + (tuple(f"b_over_bPbr_rho{x}" for x in RHOS) if stat == "RA" else ()):
            E.boot(f"primary/{stat}/{m}", r[frame][m], SRC["primary"])
    for name, s in r.items():
        if isinstance(s, dict) and "b_over_bP" in s:
            for m in ("a", "b", "b_over_a", "bP", "b_over_bP"):
                E.boot(f"gi/{name}/{m}", s[m], SRC["ext"])
        elif isinstance(s, dict) and "value" in s:
            E.boot(f"gi/{name}", s, SRC["ext"], interval_note="paired bootstrap difference")
    for L in (1, 2, 3):
        d = fe[f"L{L}"]
        for k, s in d.items():
            if k != "meta":
                E.boot(f"fe/L{L}/{k}", s, SRC["fe"], block_length=L, B=d["meta"]["B"])
    athin = i.data["thinning"]["estimates"]
    for k in ("full/full/c", "full/full/c_over_bP", "full/bP", "full/unthinned_b"):
        E.boot(f"thinning/authors/{k.replace('full/', '', 1)}", athin[k], SRC["thinA"])
    w = ext["within"]
    for k, s in w["slopes"].items():
        E.boot(f"within/{k}/b", s, SRC["ext"])
    for k, s in w["over_bP"].items():
        E.boot(f"within/{k}/over_bP", s, SRC["ext"])
    for k, s in w["pooled"].items():
        E.boot(f"within/pooled/{k}", s, SRC["ext"])
    for k, s in ext["sensitivity"]["specs"].items():
        for m in ("a", "b", "b_over_a", "b_over_bP"):
            E.boot(f"sensitivity/{k}/{m}", s[m], SRC["ext"], n=s["n"])
    se = ext["seasons"]
    for k, s in se["pooled_loss_season_intercepts"].items():
        E.boot(f"seasons/common_slope/{k}", s, SRC["ext"])
    t = se["test"]
    for k in ("wald_log_a", "wald_raw_a", "wald_sep_log_ba"):
        E.add(f"seasons/test/{k}", t[k]["stat"], None, SRC["ext"], kind=None, df=t[k]["df"], p=t[k]["p"], param=t[k].get("param"))
    q = t["cochran_q_log_a"]
    E.add("seasons/test/cochran_q_log_a", q["Q"], None, SRC["ext"], kind=None, df=q["df"], p=q["p"], I2=q["I2"])
    E.add("seasons/test/mdd_factor", t["mdd_factor"], None, SRC["ext"], kind=None, noncentrality_power80=t["noncentrality_power80"],
          mdd_max_over_min_a=t["mdd_max_over_min_a"], observed_max_over_min_a=t["observed_max_over_min_a"], definition=t["mdd_definition"])
    E.add("seasons/test/mdd_max_over_min_a", t["mdd_max_over_min_a"], None, SRC["ext"], kind=None,
          observed_max_over_min_a=t["observed_max_over_min_a"])
    E.add("seasons/test/n_draws_used", t["n_draws_used"], None, SRC["ext"], kind=None, n_draws_total=t["n_draws_total"],
          per_test={k: d["n_draws_used"] for k, d in t["detail"].items()},
          rule="draws with any non-finite season value are dropped: a failed fit, a season with no resampled weight, or a non-positive intercept on the log scale")
    E.boot("seasons/descriptive/max_over_min_b_over_a", se["descriptive"]["max_over_min_b_over_a"], SRC["ext"], note="descriptive; not a homogeneity test")
    for g in ("state", "season"):
        v = ext["validation"][g]
        E.boot(f"validation/{g}/power_minus_affine", v["bootstrap"], SRC["ext"], paired_fold_se=v["paired_fold_se"], z=v["z"],
               affine_better_folds=v["affine_better_folds"], folds=v["folds"])
    t1 = ext["table1"]
    E.add("table1/m_cross_known", t1["m_cross_known"], None, SRC["ext"], kind=None)
    E.add("table1/m_cross_observed_delta", t1["m_cross_observed_delta"], None, SRC["ext"], kind=None)
    for row in t1["rows"]:
        m = int(row["m"])
        sd = row["observed_denominator_mc_se"]
        E.add(f"table1/m{m}/observed_denominator_mc", row["observed_denominator_mc"], [row["observed_denominator_mc"] - 1.96 * sd, row["observed_denominator_mc"] + 1.96 * sd],
              SRC["ext"], kind="Monte Carlo +/- 1.96 SE")
        E.add(f"table1/m{m}/observed_denominator_delta", row["observed_denominator_delta"], None, SRC["ext"], kind=None)
        E.add(f"table1/m{m}/known_denominator_analytic", row["known_denominator_analytic"], None, SRC["ext"], kind=None)
    o = au["observed"]
    for k in ("a", "b", "b_over_a", "bP", "b_over_bP", "bw", "bw_over_bP", "a_week", "contrast_affine_minus_power", "n", "D_p10", "D_p90"):
        E.add(f"audit/observed/{k}", o[k], None, SRC["audit"], kind=None)
    E.add("audit/observed/D_p50", o["D_median"], None, SRC["audit"], kind=None)
    for sk, s in au["scenarios"].items():
        defined = s["a_over_true_total"].get("point") is not None
        for m in ("a", "err_a", "b", "bP", "b_over_bP", "bw_over_bP", "week_reduction", "contrast", "pooled_over_within_slope", "n_retained", "true_total_kept_rows") + (("b_over_a", "a_recovery_kept_rows") if defined else ()):
            E.mc(f"audit/{sk}/{m}", s[m], SRC["audit"])
        for m, k in (("D_p10", "D_p10"), ("D_p50", "D_median"), ("D_p90", "D_p90")):
            E.mc(f"audit/{sk}/{m}", s[k], SRC["audit"])
        if defined:
            E.mc(f"audit/{sk}/a_over_true_total", s["a_over_true_total"], SRC["audit"], point_key="point")
        E.mc(f"audit/{sk}/c_over_bP", s["thinning"]["c_over_bP"], SRC["audit"])
        for m in ("true_relative_variance", "observed_contrast_quantile", "shape_contrast_quantile", "lambda", "sigma"):
            if s.get(m) is not None:
                E.add(f"audit/{sk}/{m}", s[m], None, SRC["audit"], kind=None)
    rs = au["scenarios"]["S6_reporting_shock"]["reporting_shock"]
    for m in ("vR_rel", "vH_rel", "lag1_corr"):
        E.mc(f"audit/S6_reporting_shock/{m}", rs[m], SRC["audit"])
    for m in ("a_over_vR", "err_a_over_vH"):
        E.mc(f"audit/S6_reporting_shock/{m}", rs[m], SRC["audit"], point_key="point")
    E.add("audit/S6_reporting_shock/s", rs["s"], None, SRC["audit"], kind=None)
    E.add("audit/S6_reporting_shock/tau2", rs["tau2"], None, SRC["audit"], kind=None)
    E.add("audit/S6_reporting_shock/vH_analytic", rs["vH_analytic_2kappa2s2"], None, SRC["audit"], kind=None)
    b = au["baseline"]
    E.add("baseline/independent_poisson_bP", b["independent_poisson_bP"], None, SRC["audit"], kind=None)
    for rho in RHOS:
        an = b["analytic"][f"rho{rho}"]
        E.add(f"baseline/branching_analytic/rho{rho}/over_bP", an["over_bP"], None, SRC["audit"], kind=None)
        E.add(f"baseline/branching_analytic/rho{rho}/bPbr", an["bPbr"], None, SRC["audit"], kind=None)
        E.mc(f"baseline/branching_sim/rho{rho}/b_over_bP", b["simulation"][f"rho{rho}"]["b_over_bP"], SRC["audit"])
    for g in au["rho_k_grid"]:
        E.mc(f"grid/rho{g['rho']:g}/k{'poisson' if g['k'] is None else format(g['k'], 'g')}/b_over_bP", g["b_over_bP"], SRC["audit"])
    for row in th["per_q"]:
        E.add(f"thinning_exact/q{row['q']:.1f}/rate", row["rate_exact"], None, SRC["thin"], kind=None, over_bP=row["rate_exact_over_bP"])
    for k in ("c_exact_all_q", "c_exact_all_q_over_bP", "c_exact_without_q01", "c_exact_without_q01_over_bP"):
        E.add(f"thinning_exact/{k}", th[k], None, SRC["thin"], kind=None)
    for qq, l in th["lsq_leverage_by_q"].items():
        E.add(f"thinning_exact/q{qq}/leverage", l, None, SRC["thin"], kind=None)
    for c in th["mc_check"]:
        E.add(f"thinning_exact/mc_check/q{c['q']}", c["independent_mc_rate"], [c["independent_mc_rate"] - 1.96 * c["mc_se"], c["independent_mc_rate"] + 1.96 * c["mc_se"]],
              SRC["thin"], kind="Monte Carlo +/- 1.96 SE", expected_loss_rate=c["expected_loss_rate"])
    names = [e["name"] for e in E.items]
    dup = {n for n in names if names.count(n) > 1}
    if dup:
        raise rr.RenderError(f"duplicate estimate names: {sorted(dup)}")
    return E.items


if __name__ == "__main__":
    main()
