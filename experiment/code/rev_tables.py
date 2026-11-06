"""LaTeX tables for the 2026-11 revision analyses and the merge into analysis/results.json (key revision_2026_11).

Reads analysis/rev_intercept.json, rev_season.json, rev_identification.json, rev_generators.json, rev_rho_grid.json,
rev_daily_renewal.json; writes tables/tab_rev_*.tex and updates analysis/results.json in place (other keys untouched).
Usage: python code/rev_tables.py
"""
from __future__ import annotations

import json
import time

import numpy as np

import rev_common as R
from ma_models import Log


def load(name):
    return json.loads((R.AN / name).read_text())


def f3(v, d=4):
    return "--" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:.{d}f}"


def ci(v, d=4):
    if v is None or v[0] is None:
        return "--"
    if v[0] < 0:
        return f"$-${abs(v[0]):.{d}f} to {v[1]:.{d}f}"
    return f"{v[0]:.{d}f}--{v[1]:.{d}f}"


def pfmt(p, mcse=None):
    if p is None:
        return "--"
    s = "$<10^{-4}$" if p < 1e-4 else f"{p:.4f}" if p < .001 else f"{p:.3f}"
    return s + (f" ({mcse:.3f})" if mcse else "")


def pct(v):
    return f"$-${abs(100 * v):.1f}" if 100 * v <= -0.05 else f"{abs(100 * v) if abs(100 * v) < 0.05 else 100 * v:.1f}"


def esc(s):
    return str(s).replace("_", r"\_").replace("%", r"\%")


# Presentation-only relabelling of JSON labels for the manuscript tables (2026-11-05 manuscript revision; numbers unchanged)
TEX_LABELS = [("(PRIMARY, locked)", "(primary, locked plan)"), ("chi2(3)", r"$\chi^2_3$"), ("sandwich W", r"sandwich $W$"),
              ("bootstrap-covariance W", r"bootstrap-covariance $W$"), ("(fixed transmission field (shipped))", "(fixed transmission field)"),
              ("(fixed transmission field (shipped); post hoc law)", "(fixed transmission field; post hoc law)"),
              ("S1_heterogeneity", "S1"), ("S6_reporting_shock", "S6"), ("S8_combo", "S8"), ("S9_national", "S9"),
              ("rho x k grid (S1 field)", r"$\rho\times k$ grid (S1 field)"), ("generation-interval variant (r = 1)", "generation-interval variant ($r=1$)"),
              ("within-week estimand (E3, r = 1)", "within-week estimand (E3, $r=1$)"), ("coverage ge 90", r"coverage $\ge$ 90\%"), ("gi mean", "GI mean"), ("gi sd", "GI SD"), ("floor", "count floor"), (" (shipped)", ""), ("shipped (fixed field)", "original generators (fixed field)")]


def tex_label(s):
    for a, b in TEX_LABELS:
        s = s.replace(a, b)
    return s


def table(path, caption, label, cols, header, rows, notes=None, size=r"\scriptsize"):
    out = [r"\begin{table}[!tbp]", r"\centering", size, r"\setlength{\tabcolsep}{3pt}", rf"\caption{{{caption}}}", rf"\label{{{label}}}",
           rf"\begin{{tabularx}}{{\linewidth}}{{@{{}}{cols}@{{}}}}", r"\toprule", header + r" \\", r"\midrule"]
    for r in rows:
        out.append(r if r.endswith("\\\\") or r == r"\midrule" else r + r" \\")
    out += [r"\bottomrule", r"\end{tabularx}"]
    if notes:
        out.append(rf"\par\smallskip\parbox{{\linewidth}}{{\scriptsize\linespread{{1.15}}\selectfont {notes}}}")
    out.append(r"\end{table}")
    (R.TAB / path).write_text("\n".join(out) + "\n", encoding="utf-8")


def t_intercept(ri):
    e1 = ri["E1_primary"]
    inf = e1["inference"]
    rows = [rf"\multicolumn{{7}}{{l}}{{\emph{{E1, primary: two-way CGM sandwich, $t_{{{e1['t_df']}}}$}}}} \\",
            f"All rows ({e1['n_rows']}) & {f3(e1['a'])} & {ci(inf['two_way']['a_ci'])} & {f3(inf['two_way']['t'], 2)} & {pfmt(inf['two_way']['p_two_sided'])} & "
            f"{f3(e1['b_over_a'], 1)} & {f3(e1['b_over_bP'], 3)}"]
    for lab, k in (("State-cluster sandwich", "state"), ("Week-cluster sandwich", "week")):
        rows.append(f"{lab} & {f3(e1['a'])} & {ci(inf[k]['a_ci'])} & {f3(inf[k]['t'], 2)} & {pfmt(inf[k]['p_two_sided'])} & -- & --")
    w = ri["E1_wild_bootstrap_t"]
    rows.append(f"Week wild bootstrap-$t$, $a=0$ imposed ($B$={w['valid']:,}) & -- & -- & -- & {pfmt(w['p_two_sided'])} & -- & --")
    bb = ri["block_bootstrap"]
    rows.append(f"Calendar-week block bootstrap ($B$={bb['meta']['B']:,}) & -- & {ci(bb['a']['ci'])} & -- & -- & {ci(bb['b_over_a']['ci'], 1)} & {ci(bb['b_over_bP']['ci'], 2)}")
    rows.append(r"\midrule")
    rows.append(r"\multicolumn{7}{l}{\emph{Sensitivity analyses (two-way sandwich)}} \\")
    for k, s in ri["sensitivities"].items():
        if k.startswith("S-a9"):
            rows.append(f"S-a9 leave one week out (42 fits) & {ci(s['a_range'])} & min lower {f3(s['lower_limit_min'])} & -- & max {pfmt(s['max_p'])} & -- & {ci(s['b_over_bP_range'], 3)}")
            continue
        tw = s["inference"]["two_way"]
        lab = tex_label(esc(k.replace("_", " "))) + f" ({s['n_rows']})"
        rows.append(f"{lab} & {f3(s['a'])} & {ci(tw['a_ci'])} & {f3(tw['t'], 2)} & {pfmt(tw['p_two_sided'])} & {f3(s['b_over_a'], 1)} & {f3(s['b_over_bP'], 3)}")
    rr = ri["robustness_rule"]
    table("tab_rev_intercept.tex", "Count-independent intercept $a$: primary test of the analysis plan, secondary calibrations and sensitivity analyses",
          "tab:rev_intercept", r"l*{6}{>{\centering\arraybackslash}X}",
          r"Analysis (rows) & $a$ & 95\% CI & $t$ & $p$ (two-sided) & $b/a$ & $b/b_P$", rows,
          notes=f"Robustness rule (locked analysis plan): verdict \\emph{{{rr['verdict']}}}; analyses with a lower limit $\\le 0$: "
                f"{', '.join(rr['failing_sensitivities']) or 'none'}. S-a5 with threshold 1.02 reproduces the primary windows. "
                "Block-bootstrap columns give percentile intervals.")


def t_error(ri):
    rows = []
    for r in ri["error_components"]:
        rows.append(f"{r['C']:,} & {f3(r['a'])} & {f3(r['b_over_C'])} & {f3(r['total'])} & {ci(r['total_ci_sandwich'])} & {ci(r['total_ci_block_bootstrap'])} & "
                    f"{pct(r['share_from_a'])} & {pct(r['relative_rmse'])}")
    co = ri["crossover_b_over_a"]
    e3 = ri["E3_within_week"]
    rows.append(r"\midrule")
    rows.append(rf"\multicolumn{{8}}{{l}}{{Crossover $C^\times=b/a$ = {co['point']:.1f}; sandwich (log-delta) {ci(co['ci_sandwich_log_delta'], 1)}; block bootstrap {ci(co['ci_block_bootstrap'], 1)}}} \\")
    rows.append(r"\midrule")
    rows.append(r"\multicolumn{8}{l}{\emph{Within-week intercept (E3, estimand change; block bootstrap 95\% CI)}} \\")
    labs = {"week_leave_state_out": "Leave-state-out week centre (primary)", "week": "Week centre incl. own row", "week_review": "Week centre and scale",
            "season": "Season centre", "pooled": "Pooled centre (E1 loss)"}
    for m, lab in labs.items():
        v = e3[m]
        red = v.get("reduction_vs_pooled")
        rtxt = f"share of $a$ removed {pct(red['point'])}\\% ({ci([100 * red['ci'][0], 100 * red['ci'][1]], 1)})" if red else "reference"
        rows.append(rf"\multicolumn{{3}}{{l}}{{{lab}}} & {f3(v['a'])} & {ci(v['a_boot']['ci'])} & \multicolumn{{3}}{{c}}{{{rtxt}}} \\")
    table("tab_rev_error_components.tex", r"What more cases buy: predicted relative squared error $a+b/C$ by count, crossover, and the within-week intercept",
          "tab:rev_error_components", r"r*{7}{>{\centering\arraybackslash}X}",
          r"$C$ & $a$ & $b/C$ & $a+b/C$ & 95\% CI (sandwich) & 95\% CI (block bootstrap) & Share from $a$ (\%) & Relative RMSE (\%)", rows)


def t_season(rs):
    rows = []
    for r in rs["calibration_table"]:
        draws = f"{r['draws']:,}" if r["draws"] else "--"
        rows.append(f"{tex_label(r['calibration'])} & {pfmt(r['p'], r['mcse'])} & {draws}")
    rows.append(r"\midrule")
    rows.append(r"\multicolumn{3}{l}{\emph{Data sensitivities, primary calibration (calendar-week permutation, $B$=9,999)}} \\")
    for k, v in rs["data_sensitivities"].items():
        rows.append(f"{tex_label(k.replace('_', ' '))}: $W$={v['W_obs']:.2f}, {v['df']} df, {v['n_weeks']} weeks & {pfmt(v['p'], v['mcse'])} & {v['valid']:,}")
    h = rs["holm"]["adjusted"]
    table("tab_rev_season_calibration.tex", rf"Equal-season-intercept test ($W_{{\mathrm{{CGM}}}}$ = {rs['W_obs']:.2f}): primary calibration of the locked plan and every secondary calibration",
          "tab:rev_season_calibration", r"X c c", r"Calibration & $p$ (MCSE) & Draws", rows,
          notes=f"Decision rule (locked analysis plan): significant only if the primary $p<0.05$. Verdict: {esc(rs['decision']['verdict'])}. "
                f"Holm-adjusted over the two claims: E1 {pfmt(h['E1_intercept'])}, season {pfmt(h['season_primary'])}. "
                "Shipped generator calibrations hold the transmission field fixed; S9 and FD-S9 laws were specified after the data were seen.")


def t_season_effects(rs):
    e, pw = rs["effects"], rs["power"]
    rows = []
    for s, v in e["seasons"].items():
        rows.append(f"{s} & {v['weeks']} & {v['rows']} & {f3(v['a_s'])} & {ci(v['ci_sandwich_log_delta'])} & {ci(v['ci_bootstrap'])} & {f3(v['ratio_to_pooled'], 2)} "
                    f"({ci(v['ratio_to_pooled_ci_bootstrap'], 2)}) & {v['mean_reporting_coverage']:.1f}")
    mm, ph = e["max_over_min"], e["post_hoc_contrast_2023_24_vs_mean_of_others_log"]
    rows.append(r"\midrule")
    rows.append(rf"\multicolumn{{8}}{{l}}{{Highest/lowest intercept ratio {mm['point']:.2f} (bootstrap {ci(mm['ci_bootstrap'], 2)}; biased upward)}} \\")
    rows.append(rf"\multicolumn{{8}}{{l}}{{Post hoc: 2023-24 vs geometric mean of the other seasons, ratio {ph['ratio']:.2f} (sandwich {ci(ph['ci_sandwich'], 2)}; bootstrap {ci(ph['ci_bootstrap'], 2)}); no test}} \\")
    mdr = pw["mdr_max_over_min_at_80pct_power"]
    rows.append(rf"\multicolumn{{8}}{{l}}{{Design-based power at the observed pattern {pw['power_at_observed_pattern_delta_1']:.2f}; minimum detectable highest/lowest ratio at 80\% power {f3(mdr, 2)}}} \\")
    table("tab_rev_season_effects.tex", "Season intercepts (reflecting transmission and reporting variation) with sandwich and season-stratified block-bootstrap intervals, and design-based power",
          "tab:rev_season_effects", r"l*{7}{>{\centering\arraybackslash}X}",
          r"Season & Weeks & Rows & $a_s$ & 95\% CI (sandwich) & 95\% CI (bootstrap) & $a_s/a$ (95\% CI) & Coverage (\%)", rows,
          notes="Coverage: mean hospital-reporting coverage of the sample rows. Voluntary reporting May--Oct 2024; CDC data-channel pause 3 October to 13 November 2025.")


def t_implied(ri_id):
    rows = []
    by = {}
    for r in ri_id["rows"]:
        by.setdefault(r["source"], []).append(r)
    for src, rs in by.items():
        pts = [r["implied_rel_sd"] for r in rs]
        lo = min(r["implied_rel_sd_interval"][0] for r in rs)
        hi = max(r["implied_rel_sd_interval"][1] for r in rs)
        rep = [r for r in rs if r["reproduces_counting_slope"]]
        rows.append(f"{tex_label(src)} & {len(rs)} & {pct(min(pts))}--{pct(max(pts))} & {pct(lo)}--{pct(hi)} & "
                    f"{len(rep) if rs[0]['b_over_bP_range'] is not None else '--'} & {sum(r['post_hoc_law'] for r in rs)}")
    rows.append(r"\midrule")
    for lab, k in (("All rows", "envelope_all"), ("Pre-specified laws", "envelope_pre_specified_laws"),
                   ("Rows reproducing $b/b_P$", "envelope_reproducing_counting_slope")):
        v = ri_id[k]
        rows.append(f"\\emph{{{lab}}} & {v['n']} & {pct(v['point_range'][0])}--{pct(v['point_range'][1])} & "
                    f"{pct(v['interval_envelope'][0])}--{pct(v['interval_envelope'][1])} & -- & --")
    table("tab_rev_implied_sd.tex", r"Implied relative SD of the reproduction number across simulated mechanisms (a grid, not an estimate)",
          "tab:rev_implied_sd", r"X c c c c c",
          r"Mechanism family & Rows & Point range (\%) & Interval envelope (\%) & Reproduce $b/b_P$ & Post hoc laws", rows,
          notes="Implied SD $=\\sqrt{a_{\\mathrm{obs}}/r}$, $r=E(\\hat a)/v_R$ of each mechanism; intervals combine the two-way CI of $a$ with the MC 2.5--97.5\\% range of $r$. "
                "The transmission share of $a$ is not identified. Renewal rows appear twice (following-week $R$ and two-week mean $R$).")


def t_generators(g):
    rows = []
    for nm, v in g["generators"].items():
        dag = r"$^\dagger$" if v["law"]["post_hoc_law"] else ""
        rows.append(f"{nm}{dag} & {f3(v['mc_sd_a'], 5)} & {f3(v['sd_ratio_observed_se_over_mc_sd'], 2)} & "
                    f"{f3(v['mean_sandwich_se_over_mc_sd'], 2)} & {f3(v['coverage_a_ci_of_estimator_mean']['coverage'], 3)} & {f3(v['coverage_a_ci_of_true_total']['coverage'], 3)} & "
                    f"{f3(v['a_over_transmission']['point'], 3)} & {f3(v['b_over_bP']['mean'], 3)} & {pct(v['week_reduction']['mean'])} & "
                    f"{pct(v['implied_transmission_rel_sd']['point'])} ({pct(v['implied_transmission_rel_sd']['propagated'][0])}--{pct(v['implied_transmission_rel_sd']['propagated'][1])})")
    s5, s4 = g["SC5_variance_mismatch_explained"], g["SC4_sandwich_coverage_ge_0.90"]
    table("tab_rev_generators.tex", "Fresh-draw generators (transmission field redrawn in every replicate, 200 panels each)", "tab:rev_generators",
          r"l*{9}{>{\centering\arraybackslash}X}",
          r"Generator & MC SD of $\hat a$ & Observed SE / MC SD & Sandwich SE / MC SD & Coverage of $E(\hat a)$ & Coverage of true total & $E(\hat a)/v_R$ & $b/b_P$ & Week-centring reduction (\%) & Implied SD of $R$ (\%)",
          rows, notes=f"$^\\dagger$ Post hoc law (national calendar-week share). SC5 (observed SE / MC SD in [0.67, 1.5] for some generator): "
                      f"{'met' if s5['pass'] else 'not met'}, only by a post hoc law; SC4 (coverage $\\ge$ 0.90): {'met' if s4['pass'] else 'not met'}.")


def t_rho(gr, rn):
    ks = gr["design"]["k"]
    cell = {(c["rho"], c["k"]): c for c in gr["cells"]}
    rows = []
    for rho in gr["design"]["rho"]:
        ast = "*" if rho <= 0.01 else ""
        line = [f"{rho:g}{ast}"]
        for k in ks:
            c = cell[(rho, k)]
            mark = r"$^\checkmark$" if c["reproduces_counting_slope"] else ""
            line.append(f"{c['a_over_vR']['ratio_of_means']:.2f} / {c['b_over_bP']['mean']:.2f}{mark}")
        rows.append(" & ".join(line))
    rows.append(r"\midrule")
    rows.append(r"\multicolumn{6}{l}{\emph{Daily renewal, two generations per week} ($a/v_R$ / $a/v_{\bar R}$ / $b/b_P$)} \\")
    for c in rn["cells"]:
        rows.append(rf"\multicolumn{{6}}{{l}}{{$\rho$={c['rho']:g}, {'Poisson' if c['k'] is None else 'NB $k$=' + str(c['k'])}: "
                    rf"{c['a_over_vR']['ratio_of_means']:.2f} / {c['a_over_vRbar']['ratio_of_means']:.2f} / {c['b_over_bP']['mean']:.2f}}} \\")
    head = r"$\rho$ & " + " & ".join("Poisson" if k is None else f"$k$={k:g}" for k in ks)
    table("tab_rev_rho_grid.tex", r"Reporting probability $\rho$ $\times$ offspring dispersion $k$: intercept recovery $E(\hat a)/v_R$ and counting slope $b/b_P$ (100 panels per cell)",
          "tab:rev_rho_grid", r"l*{5}{>{\centering\arraybackslash}X}", head, rows,
          notes=f"Cells give $E(\\hat a)/v_R$ / mean $b/b_P$. * Assumed range for $\\rho$ (0.001--0.01; CDC burden estimates, not estimated). "
                f"$^\\checkmark$ MC 2.5--97.5\\% range of $b/b_P$ contains the observed {gr['observed_b_over_bP']:.3f}. "
                f"SC2 (recovery within 0.9--1.1 in every cell): {'met' if gr['SC2_intercept_recovery_0.9_1.1_all_cells']['pass'] else 'not met'}; "
                f"SC3 ($b/b_P$ within 0.9--1.1 at $\\rho\\le0.01$): {'met' if gr['SC3_b_over_bP_0.9_1.1_at_rho_le_0.01']['pass'] else 'not met'}; "
                f"SC1 (renewal $a/v_R$ within 0.8--1.25): {'met' if rn['SC1_one_generation_map_carries_over']['pass'] else 'not met'}.")


def merge(ri, rs, rid, g, gr, rn):
    path = R.AN / "results.json"
    res = json.loads(path.read_text())
    e1 = ri["E1_primary"]
    strip = lambda d, drop: {k: v for k, v in d.items() if k not in drop}
    res["revision_2026_11"] = {
        "generated_utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),
        "analysis_plan": {"plan": "analysis/analysis_plan.md", "lock": "analysis/analysis_plan.lock",
                          "plan_locked_utc": "2026-11-05 05:42:50",
                          "lock_scope": "sha256 of the first 14,195 bytes (the plan up to '(none at lock time)'; Amendment 1 was appended after the lock)",
                          "amendment_1": "analysis/analysis_plan.md, section Amendment 1 (implementation details)",
                          "amendment_2": "analysis/analysis_plan_amendment2.md", "amendment_2_lock": "analysis/analysis_plan_amendment2.lock",
                          "amendment_2_locked_utc": "2026-11-05 10:16:29"},
        "data": "frozen NHSN 2026-10-20 vintage; R_A sample, 1,748 rows, 51 states, 42 calendar weeks",
        "E1_intercept": {"source": "code/rev_intercept.py", "primary": strip(e1, {"V2_two_way_a_b"}), "wild_bootstrap_t": ri["E1_wild_bootstrap_t"],
                         "block_bootstrap": ri["block_bootstrap"], "sensitivity_lower_limits": ri["sensitivity_lower_limits"],
                         "sensitivities": {k: strip(v, {"V2", "fits"}) for k, v in ri["sensitivities"].items()},
                         "robustness_rule": ri["robustness_rule"]},
        "E3_within_week": {"source": "code/rev_intercept.py", **ri["E3_within_week"]},
        "counting_slope_and_crossover": {"source": "code/rev_intercept.py", "b": e1["b"], "b_P": e1["b_P"], "b_over_bP": e1["b_over_bP"],
                                         "b_over_bP_ci_sandwich_bP_fixed": e1["b_over_bP_ci_two_way"], "b_over_bP_ci_block_bootstrap": ri["block_bootstrap"]["b_over_bP"],
                                         "crossover": ri["crossover_b_over_a"], "error_components": ri["error_components"]},
        "season_test": {"source": "code/rev_season.py", "W_obs": rs["W_obs"], "primary": rs["primary"], "decision": rs["decision"], "holm": rs["holm"],
                        "calibration_table": rs["calibration_table"], "data_sensitivities": rs["data_sensitivities"],
                        "circular_shift": strip(rs["circular_shift"], {"W_by_shift"}), "effects": rs["effects"], "power": rs["power"],
                        "reporting_context": rs["reporting_context"]},
        "implied_sd_grid": {"source": "code/rev_identification.py", **strip(rid, {"rows"}), "n_rows": len(rid["rows"])},
        "fresh_draw_generators": {"source": "code/rev_generators.py", **g},
        "rho_k_grid": {"source": "code/rev_rho_grid.py", **strip(gr, {"cells"}),
                       "cells": [{k: c[k] for k in ("rho", "k", "a_over_vR", "b_over_bP", "reproduces_counting_slope", "lambda", "lambda_at_bound")} for c in gr["cells"]]},
        "daily_renewal": {"source": "code/rev_daily_renewal.py", **rn},
    }
    R.write_json(path, res)


def main():
    log = Log(R.LOGS / "rev_tables.log")
    ri, rs, rid = load("rev_intercept.json"), load("rev_season.json"), load("rev_identification.json")
    g, gr, rn = load("rev_generators.json"), load("rev_rho_grid.json"), load("rev_daily_renewal.json")
    t_intercept(ri)
    t_error(ri)
    t_season(rs)
    t_season_effects(rs)
    t_implied(rid)
    t_generators(g)
    t_rho(gr, rn)
    merge(ri, rs, rid, g, gr, rn)
    log("rev_tables done", sorted(p.name for p in R.TAB.glob("tab_rev_*.tex")), "results.json key revision_2026_11 written")


if __name__ == "__main__":
    main()
