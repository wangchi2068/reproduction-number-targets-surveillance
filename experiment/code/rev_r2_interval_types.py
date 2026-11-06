"""Interval types quoted in the manuscript (analysis_plan_amendment2.md, Section G; review N6). Descriptive.

Every interval quoted in the abstract, results and discussion is listed with its quantity, the quoted limits, the
interval type stated by the producing script and the location of the limits in the analysis outputs. The source is
found mechanically, unless a hint names it: results.json, authors/secondary.json and ma_extensions.json are flattened and every pair of numbers stored as a 2-element
list or as sibling lower/upper fields is compared with the quoted limits at the quoted precision (also after x100
for percentages). Entries whose limits are a range of point estimates or of Monte Carlo summaries are typed as such.

Output: analysis/rev_r2_interval_types.json. Usage: python code/rev_r2_interval_types.py
"""
from __future__ import annotations

import json
from decimal import Decimal

import rev_common as R
import rev_r2_lib as L2

OUT = R.AN / "rev_r2_interval_types.json"
SANDWICH = "two-way CGM sandwich (state x calendar week), t with G-1 = 41 df"
WEEKBOOT = "calendar-week block bootstrap percentile (1,999 draws)"
PIGEON = "state x calendar-week pigeonhole bootstrap percentile"
STATEBOOT = "state bootstrap percentile (states resampled)"
PANELBOOT = "panel bootstrap percentile over simulated panels (2,000 draws)"
RANGE = "range of point values (not an interval)"
MCRANGE = "range of Monte Carlo summaries over generators or cells (not an interval)"
F, FORM, BUY, POP, CONT, BIAS, REP, SLOPE, SEAS, VAL, AB, APP = (
    "Results: count-independent floor", "Results: fit of the affine form", "Results: what more cases buy",
    "Results: jurisdiction size", "Results: what the floor contains", "Results: intercept and its sources",
    "Results: external reporting data", "Results: counting slope", "Results: season intercepts",
    "Results: validity of the intervals", "Abstract", "Appendix: audit intervals")
# (location, quantity, quoted lower, quoted upper, type), in the order of the final manuscript
ENTRIES = [
    (AB + "; " + F, "a, pooled affine intercept (primary test)", "0.0109", "0.0191", SANDWICH),
    (F, "b, counting slope", "0.430", "0.690", SANDWICH),
    (F, "a, clustering by state only", "0.0125", "0.0174", "one-way CGM sandwich by state, t with 50 df"),
    (F, "a, clustering by calendar week only", "0.0112", "0.0188", "one-way CGM sandwich by calendar week, t with 41 df"),
    (F, "a, block bootstrap", "0.0107", "0.0192", WEEKBOOT),
    (F, "a, point estimates over 17 specifications", "0.0117", "0.0202", RANGE),
    (F, "a at C_t >= 50", "0.0088", "0.0172", SANDWICH),
    (F, "a, fixed calendar windows", "0.0092", "0.0311", SANDWICH),
    (F, "a, leave-one-week-out point estimates", "0.0136", "0.0156", RANGE),
    (FORM, "mean L, lowest count decile", "0.0538", "0.0854", WEEKBOOT),
    (FORM, "mean L, highest count decile", "0.0079", "0.0137", WEEKBOOT),
    (FORM, "mean L, fifth count decile", "0.0310", "0.0501", WEEKBOOT),
    (FORM, "gamma, three-parameter fit, week bootstrap", "0.355", "0.870", WEEKBOOT),
    (FORM, "gamma, three-parameter fit, state bootstrap", "0.438", "0.862", STATEBOOT + " (999 draws)"),
    (FORM, "a, three-parameter fit, week bootstrap", "-0.0014", "0.0137", WEEKBOOT),
    (FORM, "a, three-parameter fit, state bootstrap", "0.0034", "0.0144", STATEBOOT + " (999 draws)"),
    (FORM, "mean L, rows at or above P90", "0.0078", "0.0139", WEEKBOOT),
    (AB + "; " + FORM, "relative RMSE, rows at or above P90 (%)", "8.8", "11.8", WEEKBOOT),
    (FORM, "floor L - b_hat/C, P90", "0.0069", "0.0132", WEEKBOOT),
    (FORM, "floor relative SD, P90 (%)", "8.3", "11.5", WEEKBOOT),
    (FORM, "floor L - b_P/C, P90", "0.0071", "0.0132", WEEKBOOT),
    (FORM, "floor relative SD, P75 (%)", "9.8", "13.2", WEEKBOOT),
    (FORM, "floor relative SD, P95 (%)", "8.0", "11.1", WEEKBOOT),
    (BUY, "crossover b/a, two-way", "24.1", "58.3", SANDWICH + ", delta method on the log scale"),
    (BUY, "crossover b/a, block bootstrap", "23.4", "60.5", WEEKBOOT),
    (BUY, "crossover b/a over other sensitivity analyses", "25.6", "55.5", RANGE),
    (POP, "a_1, population coefficient", "-0.0085", "-0.0023", SANDWICH),
    (AB + "; " + POP, "ratio a(P90)/a(P10), linear model", "0.408", "0.609", STATEBOOT + " (1,999 draws)"),
    (POP, "ratio a(P90)/a(P10), log link", "0.324", "0.595", STATEBOOT + " (1,999 draws)"),
    (POP, "b_1, slope coefficient on population", "-0.124", "0.095", SANDWICH),
    (POP, "a, smallest-population tertile", "0.0145", "0.0233", SANDWICH + " (t with 16 df)"),
    (POP, "a, middle tertile", "0.0125", "0.0264", SANDWICH + " (t with 16 df)"),
    (POP, "a, largest-population tertile", "0.0090", "0.0180", SANDWICH + " (t with 16 df)"),
    (POP, "ratio tertile 3 / tertile 1, common slope", "0.421", "0.818", STATEBOOT + " (999 draws)"),
    (CONT, "a after week centring (leave-state-out)", "0.0078", "0.0126", WEEKBOOT),
    (CONT, "reduction of a, leave-state-out week centre (%)", "10.2", "47.8", WEEKBOOT),
    (CONT, "reduction of a, week centre incl. own row (%)", "13.8", "49.8", WEEKBOOT),
    (CONT, "reduction of a, week centre and scale (%)", "22.0", "57.0", WEEKBOOT),
    (CONT, "reduction of a, season centring (%)", "-6.8", "19.1", WEEKBOOT),
    (CONT, "implied transmission SD over 71 mechanisms (%)", "8.6", "18.6", RANGE),
    (CONT, "implied transmission SD, propagated (%)", "6.9", "22.7", "propagated envelope: two-way interval of a combined with MC 2.5-97.5% range of r"),
    (CONT, "implied SD, 13 slope-reproducing mechanisms (%)", "8.6", "13.9", RANGE),
    (CONT, "intercept recovery, daily renewal (%)", "75", "78", MCRANGE),
    (BIAS, "relative bias of a, FD-S1 (%)", "0.5", "2.3", PANELBOOT),
    (BIAS, "relative bias of a, FD-S8 (%)", "-15.1", "-13.3", PANELBOOT),
    (BIAS, "relative bias of a, FD-S9 variants (%)", "-14.9", "-13.8", RANGE),
    (AB + "; " + BIAS, "shortfall against the total with reporting shocks (%)", "6.3", "22.1", RANGE),
    (BIAS, "w_H over laws", "0.10", "0.91", RANGE),
    (BIAS, "mechanism-conditional rescaled a", "0.0160", "0.0192", RANGE),
    (REP, "change in a after participation adjustment", "-0.00086", "0.00027", WEEKBOOT),
    (REP, "a, first-issue counts", "0.0133", "0.0298", WEEKBOOT),
    (REP, "a, latest-issue counts", "0.0096", "0.0183", WEEKBOOT),
    (REP, "a, first minus latest", "0.0016", "0.0146", WEEKBOOT),
    (SLOPE, "b/b_P, plug-in reference", "0.97", "1.59", WEEKBOOT),
    (SLOPE, "b/b_P, normalised plug-in reference", "1.02", "1.64", SANDWICH + ", reference fixed"),
    (SLOPE, "b/b_P, simulated exact reference", "1.08", "1.73", SANDWICH + ", reference fixed"),
    (SLOPE, "b/b_P, state fixed effects", "0.24", "1.47", PIGEON + " (2,000 draws)"),
    (SLOPE, "b/b_P, week fixed effects", "0.96", "1.89", PIGEON + " (2,000 draws)"),
    (SLOPE, "expected b/b_P for rho <= 0.01", "0.85", "0.91", MCRANGE),
    (SEAS, "2023-24 intercept", "0.0050", "0.0100", SANDWICH + ", delta method on the log scale"),
    (SEAS, "other-season intercepts", "0.0154", "0.0180", RANGE),
    (SEAS, "highest / lowest season intercept", "1.75", "4.35", "calendar-week block bootstrap percentile within seasons (1,999 draws)"),
    (SEAS, "2023-24 versus other seasons, ratio", "0.29", "0.65", SANDWICH + ", log scale"),
    (VAL, "coverage of E(a_hat), fresh-draw generators", "0.920", "0.975", MCRANGE),
    (VAL, "SE / MC SD, fresh-draw generators", "0.92", "1.07", MCRANGE),
    (VAL, "coverage of the true total, shock generators", "0.44", "0.59", MCRANGE),
    (VAL, "observed SE / MC SD, fresh-draw laws of the plan", "1.53", "2.08", MCRANGE),
    (VAL + "; " + APP, "pigeonhole bootstrap SD / MC SD, fresh-draw generators", "1.55", "1.78", MCRANGE),
    (APP, "pigeonhole bootstrap SD / MC SD, range over panels", "1.19", "2.22", "range over simulated panels (not an interval)"),
    (VAL, "coverage of E(a_hat), state-realistic cells", "0.96", "1.00", MCRANGE),
    (VAL, "observed SE / MC SD, state-realistic cells", "1.63", "2.31", MCRANGE),
    (VAL, "E(a_hat) / simulated total, state-realistic cells (%)", "70", "78", MCRANGE),
    (VAL, "coverage of the true total, state-realistic cells", "0.015", "0.16", MCRANGE),
]
SEARCH = ("results.json", "authors/secondary.json", "ma_extensions.json")
SR = "results.json.revision_2026_11_r2.D1_state_realistic_generators.cells.*"
HINTS = {  # ranges assembled from several stored values; these take priority over mechanical matches
    "a, point estimates over 17 specifications": "results.json.revision_2026_11.E1_intercept.sensitivities.*.a (min and max)",
    "a, leave-one-week-out point estimates": "results.json.revision_2026_11.E1_intercept.sensitivities.S-a9_leave_one_week_out (min and max)",
    "crossover b/a over other sensitivity analyses": "results.json.revision_2026_11.E1_intercept.sensitivities.*.b_over_a (min and max)",
    "other-season intercepts": "results.json.revision_2026_11.season_test.effects.seasons.*.a_s (min and max over 2022-23, 2024-25, 2025-26)",
    "relative bias of a, FD-S9 variants (%)": "results.json.revision_2026_11_r2.C_total_bias.generators.FD-S9-*.relative_bias (min and max)",
    "shortfall against the total with reporting shocks (%)": "results.json.revision_2026_11_r2.C_total_bias.generators.*.relative_bias (min and max over generators_with_reporting_shock)",
    "ratio a(P90)/a(P10), log link": "results.json.revision_2026_11_r2.B_floor.B2_population.secondary_log_link.theta_ci_state_bootstrap, mapped by exp(theta (z_P90 - z_P10))",
    "w_H over laws": "results.json.revision_2026_11_r2.C_total_bias.generators.*.w_H (min and max)",
    "mechanism-conditional rescaled a": "results.json.revision_2026_11_r2.C_total_bias.mechanism_conditional_rescaling.a_obs_over_t_g_range",
    "coverage of E(a_hat), fresh-draw generators": "results.json.revision_2026_11.fresh_draw_generators.generators.* (min and max over generators)",
    "SE / MC SD, fresh-draw generators": "results.json.revision_2026_11.fresh_draw_generators.generators.* (min and max over generators)",
    "coverage of the true total, shock generators": "results.json.revision_2026_11.fresh_draw_generators.generators.* (min and max over shock generators)",
    "observed SE / MC SD, fresh-draw laws of the plan": "results.json.revision_2026_11.fresh_draw_generators.mc_sd_over_observed_se (reciprocals, FD-S1 to FD-S9-50)",
    "pigeonhole bootstrap SD / MC SD, fresh-draw generators": "results.json.revision_2026_11_r2.F1_bootstrap_sd.rule.ratios (min and max)",
    "pigeonhole bootstrap SD / MC SD, range over panels": "results.json.revision_2026_11_r2.F1_bootstrap_sd.generators.* (min and max of the per-panel ratio ranges)",
    "coverage of E(a_hat), state-realistic cells": SR + ".coverage_of_E_a.coverage (min and max)",
    "observed SE / MC SD, state-realistic cells": SR + ".sd_ratio_observed_se_over_mc_sd (min and max)",
    "E(a_hat) / simulated total, state-realistic cells (%)": SR + ".a_over_total (min and max)",
    "coverage of the true total, state-realistic cells": SR + ".coverage_of_true_total.coverage (min and max)",
}

def flatten(obj, path, out):
    if isinstance(obj, dict):
        lo = obj.get("lower", obj.get("lo"))
        hi = obj.get("upper", obj.get("hi"))
        if isinstance(lo, (int, float)) and isinstance(hi, (int, float)):
            out.append((path, float(lo), float(hi), obj.get("interval") or obj.get("method") or obj.get("type")))
        for k, v in obj.items():
            flatten(v, f"{path}.{k}" if path else str(k), out)
    elif isinstance(obj, list):
        if len(obj) == 2 and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in obj):
            out.append((path, float(obj[0]), float(obj[1]), None))
        elif len(obj) < 400:
            for i, v in enumerate(obj):
                name = v.get("name") if isinstance(v, dict) else None
                flatten(v, f"{path}[{name or i}]", out)


def rounds_to(value, quoted):
    q = Decimal(quoted)
    return Decimal(repr(round(value, -q.as_tuple().exponent))).quantize(q) == q if abs(value) < 1e6 else False


def main():
    L2.register_script(__file__)
    pairs = []
    for name in SEARCH:
        p = R.AN / name
        try:
            js = json.loads(p.read_text())
        except Exception:
            continue
        flatten(js, str(p.relative_to(R.AN)), pairs)
    rows = []
    for loc, qty, lo, hi, typ in ENTRIES:
        hits = []
        for scale in (1, 100):
            for path, a, b, lab in pairs:
                if rounds_to(a * scale, lo) and rounds_to(b * scale, hi):
                    hits.append({"path": path, "scale": scale, "label": lab})
        hits = sorted(hits, key=lambda h: (not h["path"].startswith("results.json"), "meta" in h["path"], "E1_intercept" not in h["path"], len(h["path"])))
        rows.append({"location": loc, "quantity": qty, "quoted": [lo, hi], "type": typ,
                     "source_found": bool(hits) or qty in HINTS,
                     "source": HINTS.get(qty) or (hits[0]["path"] if hits else None),
                     "source_label": None if qty in HINTS else next((h["label"] for h in hits if h["label"]), None),
                     "n_matches": 0 if qty in HINTS else len(hits),
                     "other_sources": [] if qty in HINTS else [h["path"] for h in hits[1:6]]})
    res = {"meta": {"script": "code/rev_r2_interval_types.py", "plan": "analysis/analysis_plan_amendment2.md, Section G",
                    "note": "descriptive; no decision rule"},
           "entries": rows, "n": len(rows), "n_source_found": int(sum(r["source_found"] for r in rows))}
    R.write_json(OUT, res)
    for r in rows:
        print(f"{r['quantity'][:48]:48s} {r['quoted']} {r['source']} | {r['source_label']}")


if __name__ == "__main__":
    main()
