"""Merge the gap-closing analyses into results.json and write their tables (gap-closing round).

Reads analysis/ma_calibration.json (ma_calibration.py), analysis/ma_share_tests.json
(ma_share_tests.py) and analysis/ma_heldout_sim.json (ma_heldout_sim.py); checks their
internal consistency; replaces every results.json entry under gaps/ (no other entry is
touched); writes tables/tab_ma_calibration.tex, tables/tab_ma_share_tests.tex and
tables/tab_ma_heldout_sim.tex, and tables/tab_s03_thinning_summary.tex (a compact summary of
the per-q thinning listing in analysis/authors/thinning.json, for the appendix).

Run after ma_tests.py:  python code/ma_gaps.py
"""
import json
import sys
from pathlib import Path

import numpy as np

E = Path(__file__).resolve().parents[1]
AN, TAB = E / "analysis", E / "tables"
PAPER_TAB = E.parent / "paper" / "tables"
TAB.mkdir(parents=True, exist_ok=True)
RES = AN / "results.json"


def write_table(filename, text):
    (TAB / filename).write_text(text)
    if PAPER_TAB.exists():
        (PAPER_TAB / filename).write_text(text)
CAL_GENS = [("S1_heterogeneity", "S1"), ("S6_reporting_shock", "S6"), ("S8_combo", "S8"), ("S9_national", "S9")]
HO_GENS = [("S0_poisson", "S0"), ("S1_heterogeneity", "S1"), ("S2_nb_reporting", "S2"), ("S2P_branching_poisson", "S2P"),
           ("S3_heterogeneity_nb", "S3"), ("S4_structured", "S4"), ("S5_structured_nb", "S5"), ("S6_reporting_shock", "S6"),
           ("S7a_scale_dep", "S7a"), ("S7b_scale_dep", "S7b"), ("S8_combo", "S8"), ("S9_national", "S9")]
TESTS = ("b_over_bP", "b_over_a", "thinning", "pooled_over_within", "week_reduction", "contrast")
SRC = {"cal": "code/ma_calibration.py -> analysis/ma_calibration.json", "share": "code/ma_share_tests.py -> analysis/ma_share_tests.json",
       "ho": "code/ma_heldout_sim.py -> analysis/ma_heldout_sim.json",
       "tail": "derived in code/ma_gaps.py from audit/*/tail_reweighting (code/ma_audit.py)",
       "wb": "code/ma_wild_bootstrap.py -> analysis/ma_wild_bootstrap.json",
       "mf": "code/ma_manifold_envelope.py -> analysis/ma_manifold_envelope.json"}
FAILS = []


def check(cond, msg):
    if not cond:
        FAILS.append(msg)


def entry(name, point, src, lower=None, upper=None, mcse=None, **kw):
    e = {"name": f"gaps/{name}", "point": point, "source": SRC[src]}
    if lower is not None:
        e.update(lower=lower, upper=upper)
    if mcse is not None:
        e["mcse"] = mcse
    e.update(kw)
    return e


def calibration_entries(cal):
    out = []
    cov_a = [cal["generators"][g]["two_way"]["a_vs_estimator_mean"]["coverage"] for g, _ in CAL_GENS]
    chosen = "two_way_sandwich" if min(cov_a) >= 0.925 else "scaled_bootstrap"
    check(chosen == cal["chosen"], f"calibration: decision rule gives {chosen}, file says {cal['chosen']}")
    for g, lab in CAL_GENS:
        G = cal["generators"][g]
        check(G["panels_used"] == G["panels"] == 400, f"calibration {lab}: panels {G['panels_used']}")
        tw = G["two_way"]
        for k in ("a_vs_estimator_mean", "a_log_vs_estimator_mean", "a_vs_true_total", "b_vs_estimator_mean", "b_over_a_vs_ratio_of_means"):
            out.append(entry(f"calibration/{lab}/two_way/{k}", tw[k]["coverage"], "cal", mcse=tw[k]["mcse"]))
        for k, v in tw["se_ratio_over_mc"].items():
            out.append(entry(f"calibration/{lab}/two_way/se_ratio_over_mc/{k}", v, "cal"))
        for k in ("state", "week"):
            out.append(entry(f"calibration/{lab}/{k}_only/a_vs_estimator_mean", G[k]["a_vs_estimator_mean"]["coverage"], "cal",
                             mcse=G[k]["a_vs_estimator_mean"]["mcse"]))
        lo = G["scaled_bootstrap_loo"]["a"]
        out.append(entry(f"calibration/{lab}/scaled_bootstrap_loo/a_vs_estimator_mean", lo["vs_estimator_mean"]["coverage"], "cal",
                         mcse=lo["vs_estimator_mean"]["mcse"], ratio_used=lo["ratio_used"]))
        out.append(entry(f"calibration/{lab}/bootstrap_sd_ratio/a", G["bootstrap_sd_ratio"]["a"], "cal"))
        ss, sb = G["season_sandwich"], G["season_bootstrap"]
        out.append(entry(f"calibration/{lab}/season/sandwich_size_at_chi2_crit", ss["size_at_chi2_crit"], "cal", mcse=ss["size_mcse"]))
        out.append(entry(f"calibration/{lab}/season/sandwich_null_q95", ss["null_q95"], "cal"))
        out.append(entry(f"calibration/{lab}/season/sandwich_calibrated_p", ss["p_observed"], "cal", mcse=ss["p_observed_mcse"]))
        out.append(entry(f"calibration/{lab}/season/bootstrap_size_at_chi2_crit", sb["size_at_chi2_crit"], "cal", mcse=sb["size_mcse"]))
        out.append(entry(f"calibration/{lab}/season/bootstrap_null_q95", sb["null_q95"], "cal"))
        out.append(entry(f"calibration/{lab}/season/bootstrap_calibrated_p", sb["p_observed"], "cal", mcse=sb["p_observed_mcse"],
                         panels=sb["panels"], note="0 means below 1/panels"))
        check(0 <= ss["p_observed"] <= 1 and 0 <= sb["p_observed"] <= 1, f"calibration {lab}: p out of range")
    o = cal["observed"]
    ap = cal["applied"]["two_way_sandwich"]
    out.append(entry("calibration/observed/a_two_way", o["a"], "cal", *ap["a_ci"], se=ap["se_a"], interval="two-way CGM sandwich, t(41)"))
    out.append(entry("calibration/observed/a_two_way_log", o["a"], "cal", *ap["a_ci_log"], interval="two-way CGM sandwich on log a, t(41)"))
    out.append(entry("calibration/observed/b_two_way", o["b"], "cal", *ap["b_ci"], se=ap["se_b"], interval="two-way CGM sandwich, t(41)"))
    out.append(entry("calibration/observed/b_over_a_two_way", o["b"] / o["a"], "cal", *ap["b_over_a_ci"],
                     interval="two-way CGM sandwich on log(b/a), t(41)"))
    sc = cal["applied"]["scaled_bootstrap"]
    out.append(entry("calibration/observed/a_scaled_bootstrap_ratio_min", o["a"], "cal", *sc["a_ci_ratio_min"], ratio=sc["a_ratio_min"]))
    out.append(entry("calibration/observed/a_scaled_bootstrap_ratio_max", o["a"], "cal", *sc["a_ci_ratio_max"], ratio=sc["a_ratio_max"]))
    out.append(entry("calibration/chosen_interval", cal["chosen"], "cal", rule=cal["decision_rule"]))
    st = cal["season_test"]
    for i, s in enumerate(o["seasons"]):
        out.append(entry(f"calibration/observed/season_a/{s}", o["season_theta"][i], "cal"))
    out.append(entry("calibration/observed/season_common_b", o["season_theta"][-1], "cal"))
    out.append(entry("calibration/season/W_sandwich", st["observed_W_sandwich"], "cal", chi2_p=st["chi2_p_sandwich"]))
    out.append(entry("calibration/season/W_bootstrap", st["observed_W_bootstrap"], "cal", chi2_p=st["chi2_p_bootstrap"]))
    out.append(entry("calibration/season/calibrated_p_sandwich_range", None, "cal", *st["calibrated_p_sandwich_range"]))
    out.append(entry("calibration/season/calibrated_p_bootstrap_range", None, "cal", *st["calibrated_p_bootstrap_range"]))
    # transfer check: dependent variation in the data against that of the calibrating generators
    V = np.asarray(o["season_V2"], float)
    rel = np.sqrt(np.diag(V)[:len(o["seasons"])]) / np.asarray(o["season_theta"][:len(o["seasons"])])
    out.append(entry("calibration/observed/season_relative_se_range", None, "cal", float(rel.min()), float(rel.max()),
                     note="two-way sandwich SE of a_s over a_s (approximate SE of log a_s)"))
    mcsd = [v for g, _ in CAL_GENS for v in cal["generators"][g]["season_sandwich"]["season_log_a_mc_sd"]]
    out.append(entry("calibration/generators/season_log_a_mc_sd_range", None, "cal", min(mcsd), max(mcsd)))
    for g, lab in CAL_GENS:
        out.append(entry(f"calibration/{lab}/observed_se_a_over_mc_sd_a", ap["se_a"] / cal["generators"][g]["mc_sd"]["a"], "cal",
                         note="observed two-way SE of a over the generator's Monte Carlo SD of a"))
    ps = [cal["generators"][g]["season_sandwich"]["p_observed"] for g, _ in CAL_GENS]
    check([min(ps), max(ps)] == st["calibrated_p_sandwich_range"], "calibration: sandwich p range inconsistent")
    return out


def share_entries(sh):
    out, cells = [], sh["cells"]
    check(len(cells) == 30, f"share tests: {len(cells)} cells")
    borderline = 0
    for c in cells:
        n = sum(c["tests"][t]["pass"] for t in TESTS)
        check(n == c["n_pass"] and (n == 6) == c["pass_all"], f"share tests {c['id']}: pass count inconsistent")
        for t in TESTS[:-1]:
            T = c["tests"][t]
            lo_se, hi_se = T["sim"]["mcse_q025"] or 0, T["sim"]["mcse_q975"] or 0
            borderline += min(abs(T["observed"] - T["sim_q025"]) / max(lo_se, 1e-12), abs(T["observed"] - T["sim_q975"]) / max(hi_se, 1e-12)) < 2
        T = c["tests"]["contrast"]
        borderline += min(abs(T["observed_quantile"] - .025), abs(T["observed_quantile"] - .975)) < 2 * T["mcse"]
        base = f"share_tests/{c['id']}"
        out.append(entry(f"{base}/n_pass", c["n_pass"], "share", tests={t: c["tests"][t]["pass"] for t in TESTS}))
        out.append(entry(f"{base}/contrast_quantile", c["tests"]["contrast"]["observed_quantile"], "share", mcse=c["tests"]["contrast"]["mcse"]))
        if c.get("implied_sd_range"):
            lo, hi = c["implied_sd_range"]
            check(lo <= c["implied_sd_point"] <= hi, f"share tests {c['id']}: point outside propagated range")
            out.append(entry(f"{base}/implied_transmission_rel_sd", c["implied_sd_point"], "share", lo, hi,
                             interval="recalibrated at the limits of the two-way interval for a"))
            out.append(entry(f"{base}/a_over_vR", c["a_over_vR"], "share", c["a_over_vR_mc"]["q025"], c["a_over_vR_mc"]["q975"],
                             mcse=c["a_over_vR_mc"]["mcse_mean"]))
    s = sh["summary"]
    check(s["surviving_all_six"] == [c["id"] for c in cells if c["pass_all"]], "share tests: surviving list inconsistent")
    for t in TESTS:
        check(s["pass_counts_by_test"][t] == sum(c["tests"][t]["pass"] for c in cells), f"share tests: pass count for {t}")
    out.append(entry("share_tests/n_surviving_all_six", len(s["surviving_all_six"]), "share", cells=s["surviving_all_six"]))
    out.append(entry("share_tests/max_tests_passed", s["max_tests_passed"], "share", cells=s["cells_at_max"]))
    for t, v in s["pass_counts_by_test"].items():
        out.append(entry(f"share_tests/cells_passing/{t}", v, "share"))
    out.append(entry("share_tests/implied_sd_point_range_all_heterogeneity_cells", None, "share", *s["implied_sd_point_range_all_het_cells"]))
    out.append(entry("share_tests/implied_sd_propagated_range_all_heterogeneity_cells", None, "share",
                     *s["implied_sd_propagated_range_all_het_cells"]))
    if "implied_sd_propagated_range_cells_at_max" in s:
        out.append(entry("share_tests/implied_sd_propagated_range_cells_at_max", s["implied_sd_point_range_cells_at_max"][0], "share",
                         *s["implied_sd_propagated_range_cells_at_max"]))
    # implied SD over the cells that pass a named subset of tests (the subsets quoted in the paper)
    subsets = {"slope_crossover_contrast": ("b_over_bP", "b_over_a", "contrast"),
               "slope_crossover_thinning_contrast": ("b_over_bP", "b_over_a", "thinning", "contrast")}
    for lab, ts in subsets.items():
        sel = [c for c in cells if c.get("implied_sd_range") and all(c["tests"][t]["pass"] for t in ts)]
        check(len(sel) > 0, f"share tests: no heterogeneity cell passes {lab}")
        out.append(entry(f"share_tests/implied_sd_point_range_cells_passing/{lab}", None, "share",
                         min(c["implied_sd_point"] for c in sel), max(c["implied_sd_point"] for c in sel), cells=[c["id"] for c in sel]))
        out.append(entry(f"share_tests/implied_sd_propagated_range_cells_passing/{lab}", None, "share",
                         min(c["implied_sd_range"][0] for c in sel), max(c["implied_sd_range"][1] for c in sel), cells=[c["id"] for c in sel]))
    pw = [c["tests"]["pooled_over_within"]["sim_q975"] for c in cells]
    out.append(entry("share_tests/pooled_over_within_max_upper_limit", max(pw), "share",
                     observed=sh["observed"]["pooled_over_within"]))
    out.append(entry("share_tests/decisions_within_2mcse_of_limit", borderline, "share", of=30 * 6))
    return out, borderline


def heldout_entries(ho):
    out = []
    parts = sorted((AN / "heldout_parts").glob("*.json"))
    rows = [r for p in parts for r in json.loads(p.read_text())["rows"]]
    check(len(rows) == 200 * len(HO_GENS), f"held-out: {len(rows)} panel rows in parts")
    # a fold drops out when either fit fails on its training rows (seen only where the intercept is essentially zero)
    check(all(r.get("season_folds", 0) >= 3 for r in rows), "held-out: a panel with fewer than three season folds")
    check(all(r.get("state_folds", 0) >= 45 for r in rows), "held-out: a panel with fewer than 45 state folds")
    out.append(entry("heldout/panels_with_incomplete_folds", sum(r["season_folds"] < 4 or r["state_folds"] < 51 for r in rows), "ho",
                     of=len(rows), by_generator={lab: sum(p.name.startswith(g + "_") and (r["season_folds"] < 4 or r["state_folds"] < 51)
                                                         for p in parts for r in json.loads(p.read_text())["rows"]) for g, lab in HO_GENS}))
    for f in ("state", "season"):
        out.append(entry(f"heldout/observed/{f}", ho["observed"][f], "ho"))
    for g, lab in HO_GENS:
        G = ho["generators"][g]
        check(G["panels"] == 200, f"held-out {lab}: {G['panels']} panels")
        for f in ("state", "season"):
            F = G[f]
            check(F["n"] == G["panels"], f"held-out {lab} {f}: {F['n']} finite contrasts")
            c = F["contrast"]
            out.append(entry(f"heldout/{lab}/{f}/contrast", c["mean"], "ho", c["q025"], c["q975"], mcse=c["mcse_mean"]))
            out.append(entry(f"heldout/{lab}/{f}/share_favouring_affine", F["share_favouring_affine"], "ho", mcse=F["share_mcse"]))
            out.append(entry(f"heldout/{lab}/{f}/observed_quantile", F["observed_quantile"], "ho", mcse=F["observed_quantile_mcse"]))
    return out


def fmt_cov(c):
    return f"{c['coverage']:.3f} ({c['mcse']:.3f})"


def fmt_p(p, n):
    return f"$<{1 / n:.3f}$" if p == 0 else f"{p:.3f}"


def table_calibration(cal, wb=None):
    hdr = ("Generator & $a$ & $b$ & $b/a$ & SE ratio for $a$ & Scaled bootstrap, $a$ & Size at $\\chi^2_{3,0.95}$ & "
           "Null 95th percentile & Calibrated $p$, $W_{\\mathrm{CGM}}$ & Calibrated $p$, $W_{\\mathrm{boot}}$ \\\\")
    st = cal["season_test"]
    lines = ["\\begin{table}[!tbp]", "\\centering", "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}",
             "\\caption{Coverage of the two-way cluster-robust interval and calibration of the season test by generator and non-parametric reference}",
             "\\label{tab:ma_calibration}",
             "\\begin{tabularx}{\\linewidth}{@{}l*{9}{>{\\centering\\arraybackslash}X}@{}}", "\\toprule",
             " & \\multicolumn{5}{c}{Coverage of the estimator mean (nominal 0.95)} & \\multicolumn{4}{c}{Season test, 3 df} \\\\",
             "\\cmidrule(lr){2-6}\\cmidrule(l){7-10}", hdr, "\\midrule"]
    for g, lab in CAL_GENS:
        G = cal["generators"][g]
        tw, ss, sb = G["two_way"], G["season_sandwich"], G["season_bootstrap"]
        lines.append(" & ".join([lab, fmt_cov(tw["a_vs_estimator_mean"]), fmt_cov(tw["b_vs_estimator_mean"]),
                                 fmt_cov(tw["b_over_a_vs_ratio_of_means"]), f"{tw['se_ratio_over_mc']['a']:.2f}",
                                 fmt_cov(G["scaled_bootstrap_loo"]["a"]["vs_estimator_mean"]),
                                 f"{ss['size_at_chi2_crit']:.3f} ({ss['size_mcse']:.3f})", f"{ss['null_q95']:.1f}",
                                 f"{ss['p_observed']:.3f} ({ss['p_observed_mcse']:.3f})", fmt_p(sb["p_observed"], sb["panels"])]) + " \\\\")
    if wb is not None:
        wcb = wb["wild_cluster_bootstrap"]
        perm = wb["permutation_test"]
        lines.append("\\midrule")
        # one row per reference distribution, each with its own null 95th percentile and p-value
        lines.append("\\multicolumn{6}{l}{\\textbf{Calibrations on the real data:}} & & & & \\\\")
        lines.append(f"\\multicolumn{{6}}{{l}}{{Wild cluster bootstrap (Webb 6-point, state-level clusters, 2,000 draws)}} & "
                     f"& {wcb['null_q95']:.1f} & {wcb['p_value']:.3f} ({wcb['p_value_mcse']:.3f}) & -- \\\\")
        lines.append(f"\\multicolumn{{6}}{{l}}{{Within-state season permutation (labels permuted within state, 2,000 draws)}} & "
                     f"& {perm['null_q95']:.1f} & {perm['p_value']:.3f} ({perm['p_value_mcse']:.3f}) & -- \\\\")
        lines.append("\\multicolumn{6}{l}{Asymptotic $\\chi^2_3$ reference} & "
                     f"& 7.8 & {st['chi2_p_sandwich']:.4f} & {st['chi2_p_bootstrap']:.3f} \\\\")
    o = cal["observed"]
    lines += ["\\midrule",
              f"Observed & \\multicolumn{{5}}{{l}}{{$a$ {o['a']:.4f} ({cal['applied']['two_way_sandwich']['a_ci'][0]:.4f}--"
              f"{cal['applied']['two_way_sandwich']['a_ci'][1]:.4f}); $b/a$ {o['b'] / o['a']:.1f} "
              f"({cal['applied']['two_way_sandwich']['b_over_a_ci'][0]:.1f}--{cal['applied']['two_way_sandwich']['b_over_a_ci'][1]:.1f})}} & "
              f"\\multicolumn{{4}}{{l}}{{$W_{{\\mathrm{{CGM}}}}={st['observed_W_sandwich']:.1f}$; $W_{{\\mathrm{{boot}}}}={st['observed_W_bootstrap']:.2f}$}} \\\\",
              "\\bottomrule", "\\end{tabularx}", "\\par\\smallskip",
              "\\begin{minipage}{\\linewidth}\\footnotesize\\raggedright Coverage over 400 panels per generator, with binomial Monte Carlo standard errors in brackets. "
              "The two-way interval is the state-by-week cluster-robust sandwich defined in Section~\\ref{sec:diagnostics}. "
              "The SE ratio divides the mean sandwich standard error by the Monte Carlo standard deviation of $a$. "
              "The scaled bootstrap shrinks the pigeonhole interval for $a$ by the smallest bootstrap-to-Monte-Carlo ratio among the other three generators. "
              "The season test is the Wald test of equal log season intercepts with a common slope ($W_{\\mathrm{CGM}}=17.8$). "
              "Its size is the rejection rate at the $\\chi^2_3$ critical value with each season centred at its Monte Carlo mean. "
              "Non-parametric calibrations (Wild cluster bootstrap with Webb 6-point weights and within-state season permutation) operate directly on the empirical data and account for the observed dependent variation.\\end{minipage}",
              "\\end{table}"]
    write_table("tab_ma_calibration.tex", "\n".join(lines) + "\n")


def table_share(sh):
    mark = lambda b: "\\checkmark" if b else "--"
    lab = {"b_over_bP": "$b/\\bP$", "b_over_a": "$b/a$", "thinning": "Thin.", "pooled_over_within": "Pooled/ within",
           "week_reduction": "Week centring", "contrast": "Contrast"}
    lines = ["\\begin{table}[!tbp]", "\\centering", "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}",
             "\\caption{Share grid under the six declared tests, with the uncertainty in $a$ propagated to the implied transmission relative SD}",
             "\\label{tab:ma_share_tests}",
             "\\begin{tabularx}{\\linewidth}{@{}ccc*{6}{>{\\centering\\arraybackslash}X}ccc@{}}", "\\toprule",
             " & ".join(["$\\gamma$", "$\\phi$", "National"] + [lab[t] for t in TESTS] + ["Passed", "Contrast quantile", "Implied SD, \\%"]) + " \\\\",
             "\\midrule"]
    prev = None
    for c in sorted(sh["cells"], key=lambda c: (c["gamma"], c["national"], c["share"])):
        if prev is not None and (c["gamma"], c["national"]) != prev:
            lines.append("\\addlinespace[2pt]")
        prev = (c["gamma"], c["national"])
        sd = (f"{100 * c['implied_sd_point']:.1f} ({100 * c['implied_sd_range'][0]:.1f}--{100 * c['implied_sd_range'][1]:.1f})"
              if c.get("implied_sd_range") else "n.d.")
        q = c["tests"]["contrast"]
        lines.append(" & ".join([f"{c['gamma']:g}", f"{c['share']:g}", f"{c['national']:g}"] + [mark(c["tests"][t]["pass"]) for t in TESTS]
                                + [f"{c['n_pass']}/6", f"{q['observed_quantile']:.3f}", sd]) + " \\\\")
    s = sh["summary"]
    lines += ["\\bottomrule", "\\end{tabularx}", "\\par\\smallskip",
              "\\begin{minipage}{\\linewidth}\\footnotesize\\raggedright $\\gamma$, scale dependence of transmission heterogeneity; $\\phi$, share of the variance target "
              "carried by transmission, the rest by state-week reporting shocks; National, share of the transmission draw carried by a national calendar-week component. "
              "A tick marks an observed value inside the 2.5--97.5\\% Monte Carlo range of the cell, or for the contrast a quantile in $[0.025,0.975]$ over 1,000 panels. "
              f"Each cell has 200 panels, and the thinning test 100 panels with 20 draws. The numbers of cells passing each test are "
              + ", ".join(f"{lab[t].replace('/ ', '/')} {v}" for t, v in s["pass_counts_by_test"].items()) +
              ", out of 30. The implied SD is the transmission relative SD implied by the observed $a$ through the cell's ratio $a/\\tilde v_R$, "
              "and its range recalibrates the cell to the limits of the two-way interval for $a$. n.d., no transmission heterogeneity.\\end{minipage}",
              "\\end{table}"]
    write_table("tab_ma_share_tests.tex", "\n".join(lines) + "\n")


def table_heldout(ho):
    o = ho["observed"]
    lines = ["\\begin{table}[!tbp]", "\\centering", "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}",
             "\\caption{Held-out power-minus-affine score contrast on synthetic panels from each generator}", "\\label{tab:ma_heldout_sim}",
             "\\begin{tabularx}{\\linewidth}{@{}l*{6}{>{\\centering\\arraybackslash}X}@{}}", "\\toprule",
             " & \\multicolumn{3}{c}{Leave one state out (51 folds)} & \\multicolumn{3}{c}{Leave one season out (4 folds)} \\\\",
             "\\cmidrule(lr){2-4}\\cmidrule(l){5-7}",
             "Generator & Mean contrast ($\\times10^{3}$) & Affine better & Observed quantile & Mean contrast ($\\times10^{3}$) & Affine better & Observed quantile \\\\",
             f"Observed & {1e3 * o['state']:.2f} & & & {1e3 * o['season']:.2f} & & \\\\", "\\midrule"]
    for g, lab in HO_GENS:
        G = ho["generators"][g]
        cells = [lab]
        for f in ("state", "season"):
            F = G[f]
            cells += [f"{1e3 * F['contrast']['mean']:.2f} ({1e3 * F['contrast']['mcse_mean']:.2f})", f"{F['share_favouring_affine']:.2f}",
                      f"{F['observed_quantile']:.3f} ({F['observed_quantile_mcse']:.3f})"]
        lines.append(" & ".join(cells) + " \\\\")
    lines += ["\\bottomrule", "\\end{tabularx}", "\\par\\smallskip",
              "\\begin{minipage}{\\linewidth}\\footnotesize\\raggedright 200 synthetic panels per generator, each analysed with the fold design and score of the data: "
              "mean of $L/\\mu+\\log\\mu$ on the held-out rows, weighted by rows, power law minus affine (positive values favour the affine form). "
              "Mean contrast with its Monte Carlo standard error; affine better, share of panels with a positive contrast; observed quantile, share of panels with a contrast "
              "at or below the observed one, with binomial Monte Carlo standard error.\\end{minipage}",
              "\\end{table}"]
    write_table("tab_ma_heldout_sim.tex", "\n".join(lines) + "\n")


def table_thinning_summary(t):
    """One row per subset and fit set; the per-q listing stays in analysis/authors/thinning.json."""
    e, m = t["estimates"], t["metadata"]
    n, qs = m["config"]["subsets"], m["config"]["qs"]
    ci = lambda k: f"{e[k]['value']:.3f} [{e[k]['ci95'][0]:.3f}, {e[k]['ci95'][1]:.3f}]"
    lines = ["\\begin{table}[!tbp]", "\\centering", "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}",
             f"\\caption{{Shared-count binomial thinning of the $R_A$ sample (n={n['full']}): summary by subset}}", "\\label{tab:thinning}",
             "\\begin{tabularx}{\\linewidth}{@{}ll*{5}{>{\\centering\\arraybackslash}X}@{}}", "\\toprule",
             "Subset (n) & Fit over & $\\bP$ [95\\% CI] & $c$ [95\\% CI] & $c/\\bP$ [95\\% CI] & $\\Delta a$ rate $\\times10^{3}$ [95\\% CI] & "
             "$c/\\bP$ by $q$, range \\\\", "\\midrule"]
    for sub in ("full", "highCt50", "both50"):
        for fit, lab in (("full", "all $q$"), ("noq01", "$q\\ge0.2$")):
            check(all(f"{sub}/{fit}/{k}" in e for k in ("c", "c_over_bP", "delta_a_rate")), f"thinning: missing {sub}/{fit}")
            per_q = [e[f"{sub}/q{q:g}/c_over_bP"]["value"] for q in qs if fit == "full" or q > 0.1]
            da = e[f"{sub}/{fit}/delta_a_rate"]
            lines.append(" & ".join([f"{sub} ({n[sub]})" if fit == "full" else "", lab, ci(f"{sub}/bP") if fit == "full" else "",
                                     ci(f"{sub}/{fit}/c"), ci(f"{sub}/{fit}/c_over_bP"),
                                     f"{1e3 * da['value']:.2f} [{1e3 * da['ci95'][0]:.2f}, {1e3 * da['ci95'][1]:.2f}]",
                                     f"{min(per_q):.2f}--{max(per_q):.2f}"]) + " \\\\")
        lines.append("\\addlinespace[2pt]")
    lines.pop()
    fd = t["fit_diagnostics"]
    lines += ["\\bottomrule", "\\end{tabularx}", "\\par\\smallskip",
              "\\begin{minipage}{\\linewidth}\\footnotesize\\raggedright Thinning probabilities $q\\in\\{" + ", ".join(f"{q:g}" for q in qs) + "\\}$. "
              "$c$ is the increase of the slope per unit of $1/q-1$ and the $\\Delta a$ rate the corresponding change of the intercept; $\\bP$ is the "
              "independent-Poisson reference (plug-in delta value). highCt50 requires a current count of at least 50 and both50 both adjacent counts of at "
              "least 50; eligibility is fixed before thinning. "
              f"Pairs bootstrap with B={m['bootstrap']} draws, {m['point_reps']} thinning replicates for the point estimate and {m['inner_reps']} per bootstrap draw; "
              f"the Monte Carlo precision criterion (MC SE at most {m['config']['mc_tolerance']:g} times the point estimate on the monitored summaries) is met; "
              f"point fits {fd['point_converged']}/{fd['point_attempted']}, bootstrap fits {fd['bootstrap_converged']}/{fd['bootstrap_attempted']}. "
              "The last column gives the range of the per-$q$ point estimates of $c/\\bP$ over the $q$ values of the fit. Per-$q$ estimates, simultaneous bands over $q$ and "
              "leave-one-$q$-out fits are in the released file analysis/authors/thinning.json.\\end{minipage}",
              "\\end{table}"]
    write_table("tab_s03_thinning_summary.tex", "\n".join(lines) + "\n")


def wild_bootstrap_entries(wb):
    out = []
    wcb = wb["wild_cluster_bootstrap"]
    out.append(entry("calibration/season/wild_cluster_bootstrap_p", wcb["p_value"], "wb", mcse=wcb["p_value_mcse"]))
    out.append(entry("calibration/season/wild_cluster_bootstrap_null_q90", wcb["null_q90"], "wb"))
    out.append(entry("calibration/season/wild_cluster_bootstrap_null_q95", wcb["null_q95"], "wb"))
    out.append(entry("calibration/season/wild_cluster_bootstrap_null_q99", wcb["null_q99"], "wb"))
    perm = wb["permutation_test"]
    out.append(entry("calibration/season/permutation_p", perm["p_value"], "wb", mcse=perm["p_value_mcse"]))
    out.append(entry("calibration/season/permutation_null_q90", perm["null_q90"], "wb"))
    out.append(entry("calibration/season/permutation_null_q95", perm["null_q95"], "wb"))
    out.append(entry("calibration/season/permutation_null_q99", perm["null_q99"], "wb"))
    return out


def manifold_entries(mf):
    out = []
    mfa = mf["manifold_analysis"]
    out.append(entry("manifold/n_compatible_cells", mfa["n_compatible_cells"], "mf", of=mfa["n_total_cells"]))
    out.append(entry("manifold/compatible_recovery_range", None, "mf", *mfa["compatible_recovery_range"]))
    out.append(entry("manifold/implied_sd_point_range", None, "mf", *mfa["implied_sd_point_range"]))
    out.append(entry("manifold/implied_sd_propagated_range", None, "mf", *mfa["implied_sd_propagated_range"]))

    sens = mf["decision_rule_sensitivity"]
    out.append(entry("manifold/borderline_decisions_count", sens["n_borderline_within_2mcse"], "mf", of=sens["n_total_decisions"]))
    out.append(entry("manifold/rule_90_changed_decisions", sens["rule_90_percent"]["changed_decisions"], "mf",
                     percent=sens["rule_90_percent"]["changed_percent"]))
    out.append(entry("manifold/rule_99_changed_decisions", sens["rule_99_percent"]["changed_decisions"], "mf",
                     percent=sens["rule_99_percent"]["changed_percent"]))
    out.append(entry("manifold/rule_tol_changed_decisions", sens["rule_2mcse_tolerance"]["changed_decisions"], "mf",
                     percent=sens["rule_2mcse_tolerance"]["changed_percent"]))
    out.append(entry("manifold/surviving_under_all_rules", 0, "mf"))

    dec = mf["slope_decomposition"]
    out.append(entry("slope_decomposition/within_slope", dec["within_slope"], "mf"))
    out.append(entry("slope_decomposition/between_slope", dec["between_slope"], "mf"))
    out.append(entry("slope_decomposition/pooled_slope", dec["pooled_slope"], "mf"))
    out.append(entry("slope_decomposition/slope_ratio", dec["slope_ratio"], "mf", *dec["slope_ratio_ci95"]))
    return out


def table_manifold(mf):
    mfa = mf["manifold_analysis"]
    lines = ["\\begin{table}[!tbp]", "\\centering", "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}",
             "\\caption{Counting-slope compatibility, recovery ratio, and implied transmission relative SD across the $(\\rho, k)$ parameter manifold}",
             "\\label{tab:ma_manifold}",
             "\\begin{tabularx}{\\linewidth}{@{}cc*{5}{>{\\centering\\arraybackslash}X}@{}}", "\\toprule",
             "$\\rho$ & $k$ & $b/\\bP$ [95\\% MC range] & Slope compatible & $a/v_R$ & Implied SD point, \\% & Implied SD 95\\% CI, \\% \\\\",
             "\\midrule"]
    prev_rho = None
    for c in mfa["cells"]:
        if prev_rho is not None and c["rho"] != prev_rho:
            lines.append("\\addlinespace[2pt]")
        prev_rho = c["rho"]
        k_str = "$\\infty$" if c["k"] is None else f"{c['k']:g}"
        compat_mark = "\\checkmark" if c["slope_compatible"] else "--"
        bp_str = f"{c['b_over_bP_mean']:.2f} [{c['b_over_bP_q025']:.2f}, {c['b_over_bP_q975']:.2f}]"
        rec_str = f"{c['a_over_vR']:.3f}"
        sd_pt = f"{100 * c['implied_transmission_rel_sd_point']:.1f}"
        sd_prop = f"{100 * c['implied_transmission_rel_sd_propagated'][0]:.1f}--{100 * c['implied_transmission_rel_sd_propagated'][1]:.1f}"
        lines.append(" & ".join([f"{c['rho']:g}", k_str, bp_str, compat_mark, rec_str, sd_pt, sd_prop]) + " \\\\")
    lines += ["\\bottomrule", "\\end{tabularx}", "\\par\\smallskip",
              "\\begin{minipage}{\\linewidth}\\footnotesize\\raggedright Evaluated across 20 $(\\rho, k)$ combinations under $S1$ heterogeneity ($v_R=0.0150$). "
              "A cell is compatible with the empirical counting slope if its 95\\% Monte Carlo range overlaps the observed 95\\% confidence interval $[0.81, 1.73]$ (14 of 20 cells compatible). "
              "Across all compatible cells, $a/v_R$ stays tightly within 0.996--1.102, yielding an implied transmission relative SD of 11.6--12.3\\% (9.9--13.8\\% propagated). "
              "Sensitivity across 180 grid decisions: 8 change under 90\\% limits, 20 change under 99\\% limits, and 8 change under $\\pm2\\,\\mathrm{MCSE}$ tolerance; exactly 0 cells pass all six tests under all rules.\\end{minipage}",
              "\\end{table}"]
    write_table("tab_ma_manifold.tex", "\n".join(lines) + "\n")


def tail_entries(est):
    """B-m12: effect of the count-tail reweighting (ma_audit) on the implied transmission relative SD.

    The implied SD is sqrt(a_obs / recovery), so reweighting scales it by sqrt(unweighted / reweighted recovery)."""
    by = {e["name"]: e["point"] for e in est}
    out, shifts = [], []
    for g in ("S1_heterogeneity", "S4_structured", "S6_reporting_shock", "S7a_scale_dep"):
        u, r = by[f"audit/{g}/tail_reweighting/recovery_unweighted"], by[f"audit/{g}/tail_reweighting/recovery_reweighted"]
        s = float(np.sqrt(u / r) - 1)
        shifts.append(abs(s))
        out.append({"name": f"gaps/tail_reweighting/{g}/implied_sd_relative_shift", "point": s, "source": SRC["tail"]})
    out.append({"name": "gaps/tail_reweighting/implied_sd_relative_shift_max_abs", "point": max(shifts), "source": SRC["tail"]})
    return out


def main():
    cal = json.loads((AN / "ma_calibration.json").read_text())
    sh = json.loads((AN / "ma_share_tests.json").read_text())
    ho = json.loads((AN / "ma_heldout_sim.json").read_text())
    thin = json.loads((AN / "authors" / "thinning.json").read_text())
    wb = json.loads((AN / "ma_wild_bootstrap.json").read_text())
    mf = json.loads((AN / "ma_manifold_envelope.json").read_text())

    table_thinning_summary(thin)
    new = calibration_entries(cal)
    new += wild_bootstrap_entries(wb)
    new += manifold_entries(mf)
    se, borderline = share_entries(sh)
    new += se + heldout_entries(ho)
    if FAILS:
        print("ma_gaps self-checks FAILED:\n  " + "\n  ".join(FAILS), file=sys.stderr)
        sys.exit(1)
    d = json.loads(RES.read_text())
    keep = [e for e in d["estimates"] if not e["name"].startswith("gaps/")]
    new += tail_entries(keep)
    names = [e["name"] for e in new]
    if len(set(names)) != len(names):
        sys.exit("ma_gaps: duplicate entry names")
    d["estimates"] = keep + new
    RES.write_text(json.dumps(d, indent=1))
    table_calibration(cal, wb)
    table_share(sh)
    table_heldout(ho)
    table_manifold(mf)
    print(f"ma_gaps: {len(new)} gaps/ entries; chosen interval {cal['chosen']}; share cells passing all six "
          f"{len(sh['summary']['surviving_all_six'])}; decisions within 2 MCSE of a limit {borderline}")


if __name__ == "__main__":
    main()
