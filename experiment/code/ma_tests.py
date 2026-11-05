"""Apply every declared audit test to every generator with one decision rule.

Reads analysis/results.json (written by ma_render.py) and analysis/ma_audit_national.json
(written by ma_audit_national.py), writes analysis/ma_tests.json and tables/tab_ma_tests.tex,
and merges the decisions and the S9 summaries into analysis/results.json under tests/ and
national/ (replaced on every run).

Decision rule: a generator passes a test when the observed point value lies inside the
generator's 2.5-97.5% Monte Carlo range; for the score contrast, when the pooled quantile of
the observed contrast among simulated contrasts lies in [0.025, 0.975].

Also records the season Wald statistic rescaled by the bootstrap-to-Monte-Carlo variance
ratio of a seen in the coverage study (an illustration of how far the season test depends on
interval calibration, not a corrected test).

Run after ma_render.py and ma_audit_national.py:  python code/ma_tests.py
"""
import json
from pathlib import Path

from scipy.stats import chi2

E = Path(__file__).resolve().parents[1]
RES = E / "analysis" / "results.json"
NAT = E / "analysis" / "ma_audit_national.json"
GENS = [("S0_poisson", "S0"), ("S1_heterogeneity", "S1"), ("S2_nb_reporting", "S2"),
        ("S2P_branching_poisson", "S2P"), ("S3_heterogeneity_nb", "S3"), ("S4_structured", "S4"),
        ("S5_structured_nb", "S5"), ("S6_reporting_shock", "S6"), ("S7a_scale_dep", "S7a"),
        ("S7b_scale_dep", "S7b"), ("S8_combo", "S8")]
UNDEFINED_A = {"S0", "S2", "S2P"}  # intercept essentially zero: b/a and the centring reduction are undefined
TESTS = [("b_over_bP", "b_over_bP", "$b/\\bP$"), ("b_over_a", "b_over_a", "$b/a$"),
         ("thinning", "c_over_bP", "Thinning $c/\\bP$"), ("pooled_over_within", "pooled_over_within_slope", "Pooled / within slope"),
         ("week_reduction", "week_reduction", "Week-centring reduction"), ("contrast", None, "Score contrast")]
SRC_T, SRC_N = "code/ma_tests.py -> analysis/ma_tests.json", "code/ma_audit_national.py -> analysis/ma_audit_national.json"


def ranges_from_results(est, key):
    r = {}
    for tid, gkey, _ in TESTS:
        if tid == "contrast":
            r[tid] = est[f"audit/{key}/observed_contrast_quantile_pooled"]["point"]
        elif f"audit/{key}/{gkey}" in est:
            g = est[f"audit/{key}/{gkey}"]
            r[tid] = (g["lower"], g["upper"])
    return r


def ranges_from_national(n):
    a = n["audit"]
    rng = lambda s: (s["q025"], s["q975"])
    return {"b_over_bP": rng(a["b_over_bP"]), "b_over_a": rng(a["b_over_a"]), "thinning": rng(a["thinning"]["c_over_bP"]),
            "pooled_over_within": rng(a["pooled_over_within_slope"]), "week_reduction": rng(a["week_reduction"]),
            "contrast": n["contrast"]["pooled"]}


def national_estimates(n):
    a, c, cv = n["audit"], n["contrast"], n["coverage"]
    out = []

    def add(name, s, point=None):
        out.append({"name": f"national/S9/{name}", "point": s.get("mean") if point is None else point,
                    "lower": s.get("q025"), "upper": s.get("q975"), "source": SRC_N})
    add("a_over_transmission", a["a_over_transmission"], a["a_over_transmission"]["point"])
    add("a_over_true_total", a["a_over_true_total"], a["a_over_true_total"]["point"])
    add("err_a_over_vH", a["err_a_over_vH"], a["err_a_over_vH"]["point"])
    for k in ("b_over_bP", "b_over_a", "bw_over_bP", "week_reduction", "n_retained", "D_median", "D_p90"):
        add(k, a[k])
    add("pooled_over_within_slope", a["pooled_over_within_slope"])
    add("thinning_c_over_bP", a["thinning"]["c_over_bP"])
    isd = a["implied_transmission_rel_sd"]
    out.append({"name": "national/S9/implied_transmission_rel_sd", "point": isd["point"], "lower": isd["range_from_ratio_limits"][0],
                "upper": isd["range_from_ratio_limits"][1], "source": SRC_N})
    for k in ("amp", "shock_s", "national_share"):
        out.append({"name": f"national/S9/{k}", "point": a[k], "source": SRC_N})
    out.append({"name": "national/S9/observed_contrast_quantile_pooled", "point": c["pooled"], "mcse": c["pooled_mcse"], "source": SRC_N})
    out.append({"name": "national/S9/observed_contrast_quantile_sd_streams", "point": c["sd_across_streams"], "source": SRC_N})
    for k in ("a_vs_true_total", "a_vs_estimator_mean", "b_vs_estimator_mean", "b_over_a_vs_ratio_of_means"):
        out.append({"name": f"national/S9/coverage/{k}", "point": cv[k]["coverage"], "mcse": cv[k]["mcse"], "source": SRC_N})
    for k, v in cv["se_ratio_bootstrap_over_mc"].items():
        out.append({"name": f"national/S9/coverage/se_ratio_bootstrap_over_mc/{k}", "point": v, "source": SRC_N})
    out.append({"name": "national/S9/coverage/a_over_true_total", "point": cv["a_over_true_total"], "source": SRC_N})
    return out


def main():
    d = json.loads(RES.read_text())
    keep = [e for e in d["estimates"] if not e["name"].startswith(("tests/", "national/"))]
    est = {e["name"]: e for e in keep}
    o = lambda k: est[k]["point"]
    obs = {"b_over_bP": o("audit/observed/b_over_bP"), "b_over_a": o("audit/observed/b_over_a"),
           "thinning": o("thinning/authors/full/c_over_bP"), "pooled_over_within": o("audit/observed/b") / o("audit/observed/bw"),
           "week_reduction": 1 - o("audit/observed/a_week") / o("audit/observed/a")}
    gens = [(short, ranges_from_results(est, key)) for key, short in GENS]
    nat = json.loads(NAT.read_text()) if NAT.exists() else None
    if nat:
        gens.append(("S9", ranges_from_national(nat)))
    rows, new = [], []
    for short, rg in gens:
        row = {"generator": short, "tests": {}}
        for tid, _, _ in TESTS:
            if tid == "contrast":
                q = rg[tid]
                row["tests"][tid] = {"observed": o("audit/observed/contrast_affine_minus_power"), "quantile": q,
                                     "pass": bool(0.025 <= q <= 0.975)}
            elif tid in ("b_over_a", "week_reduction") and short in UNDEFINED_A:
                row["tests"][tid] = {"observed": obs[tid], "pass": None, "note": "undefined: intercept essentially zero"}
            else:
                lo, hi = rg[tid]
                row["tests"][tid] = {"observed": obs[tid], "mc_lower": lo, "mc_upper": hi, "pass": bool(lo <= obs[tid] <= hi)}
        dec = [t["pass"] for t in row["tests"].values() if t["pass"] is not None]
        row["n_defined"], row["n_pass"] = len(dec), sum(dec)
        row["passes_all"] = row["n_pass"] == row["n_defined"]
        rows.append(row)
        for tid, t in row["tests"].items():
            if t["pass"] is not None:
                new.append({"name": f"tests/{short}/{tid}/pass", "point": int(t["pass"]), "source": SRC_T})
        new.append({"name": f"tests/{short}/n_pass", "point": row["n_pass"], "source": SRC_T})
        new.append({"name": f"tests/{short}/n_defined", "point": row["n_defined"], "source": SRC_T})
    for tid, v in obs.items():
        new.append({"name": f"tests/observed/{tid}", "point": v, "source": "code/ma_tests.py from analysis/results.json"})
    # the legacy key contrast_affine_minus_power holds power minus affine (revision_design.contrast); positive favours affine
    new.append({"name": "tests/observed/contrast_power_minus_affine", "point": o("audit/observed/contrast_affine_minus_power"),
                "source": "code/ma_tests.py; same value as audit/observed/contrast_affine_minus_power, which is power minus affine"})
    n_all = sum(r["passes_all"] for r in rows)
    new.append({"name": "tests/n_generators_passing_all", "point": n_all, "source": SRC_T})
    new.append({"name": "tests/n_generators", "point": len(rows), "source": SRC_T})

    # season test under the variance inflation of the coverage study (illustration)
    w = est["seasons/test/wald_log_a"]["point"]
    ratios = [est[f"coverage/{g}/se_ratio_bootstrap_over_mc/a"]["point"] for g in ("S1_heterogeneity", "S6_reporting_shock", "S8_combo")]
    if nat:
        ratios.append(nat["coverage"]["se_ratio_bootstrap_over_mc"]["a"])
    lo, hi = min(ratios) ** 2, max(ratios) ** 2
    season = {"wald_log_a": w, "df": 3, "variance_ratio_min": lo, "variance_ratio_max": hi,
              "rescaled_min": w * lo, "rescaled_max": w * hi,
              "p_rescaled_min": float(chi2.sf(w * lo, 3)), "p_rescaled_max": float(chi2.sf(w * hi, 3))}
    for k, v in season.items():
        if k != "df":
            new.append({"name": f"tests/season_rescaled/{k}", "point": v, "source": "code/ma_tests.py (coverage se ratios x seasons/test/wald_log_a)"})
    if nat:
        new += national_estimates(nat)

    out = {"rule": "pass = observed point inside the generator's 2.5-97.5% Monte Carlo range; contrast: pooled quantile in [0.025, 0.975]",
           "observed": obs, "generators": rows, "n_generators_passing_all": n_all, "season_rescaled": season}
    (E / "analysis" / "ma_tests.json").write_text(json.dumps(out, indent=1))
    d["estimates"] = keep + new
    RES.write_text(json.dumps(d, indent=1))
    write_table(rows, obs)
    if nat:
        add_coverage_row(nat["coverage"])
    print(f"ma_tests: {n_all} of {len(rows)} generators pass every defined test; season W x [{lo:.2f}, {hi:.2f}] = "
          f"[{w * lo:.1f}, {w * hi:.1f}]")


def cell(t, tid):
    if t["pass"] is None:
        return "n.d."
    if tid == "contrast":
        return f"{'pass' if t['pass'] else 'fail'} ({t['quantile']:.3f})"
    f = (lambda x: f"{100 * x:.1f}") if tid == "week_reduction" else (lambda x: f"{x:.2f}")
    return f"{'pass' if t['pass'] else 'fail'} ({f(t['mc_lower'])} to {f(t['mc_upper'])})"


def coverage_row(cv):
    f = lambda k: f"{cv[k]['coverage']:.3f} ({cv[k]['mcse']:.3f})"
    sr = cv["se_ratio_bootstrap_over_mc"]
    return (f"S9 S8 with half of transmission heterogeneity national-week & {f('a_vs_true_total')} & {f('a_vs_estimator_mean')} & "
            f"{f('b_vs_estimator_mean')} & {f('b_over_a_vs_ratio_of_means')} & {cv['a_over_true_total']:.3f} & "
            f"{sr['a']:.2f}, {sr['b']:.2f}, {sr['b_over_a']:.2f} \\\\")


def add_coverage_row(cv):
    """Append the S9 row (ma_audit_national.py) to tables/tab_ma_coverage.tex, written by ma_render.py for S1, S6 and S8."""
    f = E / "tables" / "tab_ma_coverage.tex"
    lines = [ln for ln in f.read_text().split("\n") if not ln.startswith("S9 ")]
    i = max(j for j, ln in enumerate(lines) if ln.startswith("S8 "))
    lines.insert(i + 1, coverage_row(cv))
    f.write_text("\n".join(lines))


def write_table(rows, obs):
    hdr = " & ".join(["Generator"] + [lab for _, _, lab in TESTS] + ["Passed"])
    fo = {"b_over_bP": f"{obs['b_over_bP']:.2f}", "b_over_a": f"{obs['b_over_a']:.1f}", "thinning": f"{obs['thinning']:.2f}",
          "pooled_over_within": f"{obs['pooled_over_within']:.2f}", "week_reduction": f"{100 * obs['week_reduction']:.1f}",
          "contrast": "$-0.0014$"}
    lines = ["\\begin{table}[!tbp]", "\\centering", "\\scriptsize", "\\setlength{\\tabcolsep}{3pt}",
             "\\caption{Declared known-truth tests by generator under one decision rule}", "\\label{tab:ma_tests}",
             "\\begin{tabularx}{\\linewidth}{@{}l*{6}{>{\\centering\\arraybackslash}X}c@{}}", "\\toprule", hdr + " \\\\",
             " & ".join(["Observed"] + [fo[t] for t, _, _ in TESTS] + [""]) + " \\\\", "\\midrule"]
    for r in rows:
        lines.append(" & ".join([r["generator"]] + [cell(r["tests"][t], t) for t, _, _ in TESTS] + [f"{r['n_pass']}/{r['n_defined']}"]) + " \\\\")
    lines += ["\\bottomrule", "\\end{tabularx}", "\\par\\smallskip",
              "\\begin{minipage}{\\linewidth}\\footnotesize\\raggedright Pass or fail under the decision rule of Section~\\ref{sec:diagnostics}. "
              "Brackets give the Monte Carlo 2.5--97.5\\% range or, for the power-minus-affine score contrast, the quantile of the observed value pooled over five streams. "
              "The week-centring reduction is in per cent. n.d., not defined because the intercept is essentially zero. "
              "The pooled/within slope ratio is $b/b_w$, with $b_w$ the state fixed-effects slope. S9 is S8 with half of the transmission "
              "heterogeneity carried by a national calendar-week component.\\end{minipage}",
              "\\end{table}"]
    (E / "tables" / "tab_ma_tests.tex").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
