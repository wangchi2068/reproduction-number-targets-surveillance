"""Implied transmission-dispersion grid across mechanisms (analysis plan Section 5.1, review M2).

The intercept a is a total of transmission and reporting dispersion; the transmission share is not identified
from a alone. For every mechanism that was simulated, r = E(a_hat) / v_R (estimator mean over true transmission
relative variance) maps the observed intercept to an implied relative SD of R, sqrt(a_obs / r). The point uses the
ratio of means; the propagated interval combines the observed two-way 95% CI of a with the 2.5-97.5% MC range of r.
The result is a grid of mechanism-dependent values, not an estimate: no mechanism is preferred by the data except
through the declared counting-slope check (MC range of b/b_P containing the observed 1.257).

Sources: shipped scenarios (analysis/ma_audit.json, ma_audit_national.json; fixed field), fresh-draw generators
(rev_generators.json), the rho x k grid (rev_rho_grid.json), the daily renewal cells (rev_daily_renewal.json),
generation-interval variants (rev_intercept.json; r = 1) and the within-week estimand (E3; r = 1).
Output: analysis/rev_identification.json.
"""
from __future__ import annotations

import json

import numpy as np

import rev_common as R
from ma_models import Log


def load(name):
    return json.loads((R.AN / name).read_text())


def row(source, mechanism, a_obs, a_ci, ratio, q025, q975, bbp=None, post_hoc=False, estimand="v_R (following-week R)", note=None):
    imp = R.implied_sd(a_obs, a_ci, ratio, q025, q975)
    return {"source": source, "mechanism": mechanism, "estimand": estimand, "ratio_a_over_v": ratio, "ratio_q025": q025, "ratio_q975": q975,
            "implied_rel_sd": imp["point"], "implied_rel_sd_interval": imp["propagated"],
            "b_over_bP_range": bbp, "reproduces_counting_slope": None if bbp is None else bool(bbp[0] <= OBS_BBP <= bbp[1]),
            "post_hoc_law": post_hoc, "note": note}


OBS_BBP = None


def main():
    global OBS_BBP
    log = Log(R.LOGS / "rev_identification.log")
    ri = load("rev_intercept.json")
    e1 = ri["E1_primary"]
    a_obs, a_ci = e1["a"], e1["inference"]["two_way"]["a_ci"]
    OBS_BBP = e1["b_over_bP"]
    rows = []
    aud = load("ma_audit.json")
    for nm, s in aud["scenarios"].items():
        r = s.get("a_over_transmission")
        if not r or r.get("point") is None or not np.isfinite(r["point"]) or s.get("sigma", 0) == 0:
            continue
        rows.append(row("shipped (fixed field)", nm, a_obs, a_ci, r["point"], r["q025"], r["q975"],
                        [s["b_over_bP"]["q025"], s["b_over_bP"]["q975"]]))
    nat = load("ma_audit_national.json")["audit"]
    r = nat["a_over_transmission"]
    rows.append(row("shipped (fixed field)", "S9_national", a_obs, a_ci, r["point"], r["q025"], r["q975"],
                    [nat["b_over_bP"]["q025"], nat["b_over_bP"]["q975"]], post_hoc=True))
    for nm, g in load("rev_generators.json")["generators"].items():
        r = g["a_over_transmission"]
        rows.append(row("fresh-draw generator", nm, a_obs, a_ci, r["point"], r["q025"], r["q975"],
                        [g["b_over_bP"]["q025"], g["b_over_bP"]["q975"]], post_hoc=g["law"]["post_hoc_law"]))
    for c in load("rev_rho_grid.json")["cells"]:
        r = c["a_over_vR"]
        rows.append(row("rho x k grid (S1 field)", f"rho={c['rho']}, k={c['k'] if c['k'] is not None else 'Poisson'}", a_obs, a_ci,
                        r["ratio_of_means"], r["q025"], r["q975"], [c["b_over_bP"]["q025"], c["b_over_bP"]["q975"]],
                        note="rho in the assumed range 0.001-0.01" if c["assumed_rho_range"] else "rho outside the assumed range"))
    for c in load("rev_daily_renewal.json")["cells"]:
        lab = f"daily renewal, rho={c['rho']}, {'Poisson' if c['k'] is None else 'k=' + str(c['k'])}"
        r = c["a_over_vR"]
        rows.append(row("daily renewal (two generations per week)", lab, a_obs, a_ci, r["ratio_of_means"], r["q025"], r["q975"],
                        [c["b_over_bP"]["q025"], c["b_over_bP"]["q975"]]))
        r = c["a_over_vRbar"]
        rows.append(row("daily renewal (two generations per week)", lab, a_obs, a_ci, r["ratio_of_means"], r["q025"], r["q975"],
                        [c["b_over_bP"]["q025"], c["b_over_bP"]["q975"]], estimand="v of the two-week mean R"))
    for k, s in ri["sensitivities"].items():
        if k.startswith(("S-a1", "S-a2")):
            ci = s["inference"]["two_way"]["a_ci"]
            rows.append({"source": "generation-interval variant (r = 1)", "mechanism": k, "estimand": "v_R under this GI",
                         "ratio_a_over_v": 1.0, "implied_rel_sd": float(np.sqrt(s["a"])),
                         "implied_rel_sd_interval": [float(np.sqrt(max(ci[0], 0))), float(np.sqrt(ci[1]))],
                         "b_over_bP_range": None, "reproduces_counting_slope": None, "post_hoc_law": False, "note": f"GI {s['gi']}"})
    e3 = ri["E3_within_week"]["week_leave_state_out"]
    rows.append({"source": "within-week estimand (E3, r = 1)", "mechanism": "leave-state-out week centring",
                 "estimand": "state-specific dispersion around the calendar-week mean", "ratio_a_over_v": 1.0,
                 "implied_rel_sd": float(np.sqrt(e3["a"])), "implied_rel_sd_interval": [float(np.sqrt(v)) for v in e3["a_boot"]["ci"]],
                 "b_over_bP_range": None, "reproduces_counting_slope": None, "post_hoc_law": False, "note": "block-bootstrap CI"})
    pts = np.array([r["implied_rel_sd"] for r in rows])
    los = np.array([r["implied_rel_sd_interval"][0] for r in rows])
    his = np.array([r["implied_rel_sd_interval"][1] for r in rows])
    pre = [r for r in rows if not r["post_hoc_law"]]
    repro = [r for r in rows if r["reproduces_counting_slope"]]
    env = lambda rs: {"n": len(rs), "point_range": [min(r["implied_rel_sd"] for r in rs), max(r["implied_rel_sd"] for r in rs)],
                      "interval_envelope": [min(r["implied_rel_sd_interval"][0] for r in rs), max(r["implied_rel_sd_interval"][1] for r in rs)]}
    out = {"meta": {"script": "code/rev_identification.py", "plan_section": "5.1", "a_obs": a_obs, "a_ci_two_way": a_ci, "observed_b_over_bP": OBS_BBP,
                    "wording": "the transmission share of a is not identified; the implied relative SD of R is a mechanism-dependent grid, not an estimate"},
           "rows": rows,
           "envelope_all": {"n": len(rows), "point_range": [float(pts.min()), float(pts.max())], "interval_envelope": [float(los.min()), float(his.max())]},
           "envelope_pre_specified_laws": env(pre),
           "envelope_reproducing_counting_slope": env(repro) if repro else None,
           "shipped_headline_range_not_used": "11.6-12.3% (S1-like single mechanism); superseded by the grid"}
    R.write_json(R.AN / "rev_identification.json", out)
    log("rev_identification done", "rows", len(rows), "all", out["envelope_all"], "pre", out["envelope_pre_specified_laws"],
        "repro", out["envelope_reproducing_counting_slope"])


if __name__ == "__main__":
    main()
