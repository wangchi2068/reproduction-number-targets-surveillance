"""Bias of the intercept relative to the transmission-plus-reporting total (analysis_plan_amendment2.md, Section C; review N3).

For each fresh-draw generator: relative bias E(a_hat) / E(v_R + v_H) - 1 with a panel bootstrap interval
(B = 2,000, seed [20261105, 71, 1]), the MC 2.5-97.5% range of the per-panel ratio, the effective reporting
weight w_H = (E(a_hat) - r_1 E(v_R)) / E(v_H) with r_1 the FD-S1 recovery ratio, and the mechanism-conditional
rescaling factor t_g = E(a_hat) / E(v_R + v_H).

Existing panels: analysis/rev_generators_<name>.jsonl (rev_generators.py). New share sweep: the FD-S8 law with
transmission share 1.0, 0.75, 0.25, calibrated as FD-S8 (calibration stream 910 + j), 200 panels, seed
[20261105, 41, j]; panels in analysis/rev_r2_share_<share>.jsonl.

Output: analysis/rev_r2_total_bias.json. Usage: python code/rev_r2_total_bias.py [--reps 200] [--workers 3]
"""
from __future__ import annotations

import argparse
import time
from multiprocessing import Pool

import numpy as np
from scipy import stats

import ma_audit as A
import rev_common as R
import rev_r2_lib as L2
from ma_models import Log
from rev_generators import FreshDraw, SPECS, panel_row

OUT = R.AN / "rev_r2_total_bias.json"
SHARES = {0: 1.0, 1: .75, 2: .25}
OBS_A = 0.0149546


def run_share(job):
    j, reps = job
    share = SHARES[j]
    name = f"FD-S8-share{int(share * 100)}"
    log = Log(R.LOGS / f"rev_r2_share_{int(share * 100)}.log")
    t0 = time.time()
    e, d, obs, _, _ = R.setup_design()
    _, k, _, gamma, nat, _ = SPECS["FD-S8"]
    z0 = np.random.default_rng([A.SEED, 1]).standard_normal(d.n)
    gen = FreshDraw(d, obs["g0"], z0, k, share, gamma, obs["kappa_mean"], nat)
    cal = A.calibrate(gen, 910 + j, log, name)
    log(name, "calibrated", "lambda", gen.lam, "amp", gen.amp, "sigma", gen.sigma, "shock_s", gen.shock_s)
    rng = np.random.default_rng([R.SEED, 41, j])
    season_all = d.e.season.astype(str).to_numpy()
    path = R.AN / f"rev_r2_share_{int(share * 100)}.jsonl"
    if path.exists():
        path.unlink()
    for r in range(reps):
        counts = gen.draw(rng, 1)[0]
        row = panel_row(gen, counts, gen.last_g[0], gen.last_shocks[0], season_all, lambda df: stats.t.ppf(.975, df))
        row.update(generator=name, r=r)
        R.append_jsonl(path, row)
        if r % 50 == 49:
            log(name, "panels", r + 1)
    info = {"generator": name, "share": share, "lambda": gen.lam, "amp": gen.amp, "sigma": gen.sigma,
            "shock_s": gen.shock_s, "calibration": cal, "seconds": round(time.time() - t0, 1)}
    R.append_jsonl(R.AN / "rev_r2_share_calibration.jsonl", info)
    return name


def arrays(rows):
    ok = [r for r in rows if r.get("a") is not None]
    return (np.array([r["a"] for r in ok]), np.array([r["truth_relvar_kept"] for r in ok]),
            np.array([r["vH_kept"] for r in ok]))


def summarise(rows, r1_draws, rng, B=2000):
    a, vR, vH = arrays(rows)
    n = len(a)
    tot = vR + vH
    idx = rng.integers(0, n, (B, n))
    bias_b = a[idx].mean(1) / tot[idx].mean(1) - 1
    out = {"panels": n, "mean_a": float(a.mean()), "mean_vR": float(vR.mean()), "mean_vH": float(vH.mean()),
           "share_vH_of_total": float(vH.mean() / tot.mean()),
           "relative_bias": float(a.mean() / tot.mean() - 1), "relative_bias_ci": R.pct_ci(bias_b),
           "per_panel_ratio_range": [float(v) for v in np.quantile(a / tot, [.025, .975])],
           "t_g": float(a.mean() / tot.mean()), "a_over_vR": float(a.mean() / vR.mean())}
    if vH.mean() > 0 and r1_draws is not None:
        r1 = float(np.mean(r1_draws["point"]))
        out["w_H"] = float((a.mean() - r1 * vR.mean()) / vH.mean())
        wb = (a[idx].mean(1) - r1_draws["boot"] * vR[idx].mean(1)) / vH[idx].mean(1)
        out["w_H_ci"] = R.pct_ci(wb)
        out["excludes_zero"] = bool(out["relative_bias_ci"][0] > 0 or out["relative_bias_ci"][1] < 0)
        out["rescaled_a_obs"] = OBS_A / out["t_g"]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--workers", type=int, default=3)
    args = ap.parse_args()
    L2.register_script(__file__)
    t0 = time.time()
    log = Log(R.LOGS / "rev_r2_total_bias.log")
    cal_path = R.AN / "rev_r2_share_calibration.jsonl"
    if cal_path.exists():
        cal_path.unlink()
    with Pool(args.workers) as pool:
        for nm in pool.imap_unordered(run_share, [(j, args.reps) for j in SHARES]):
            log("finished", nm)
    rng = np.random.default_rng([R.SEED, 71, 1])
    s1 = R.read_jsonl(R.AN / "rev_generators_FD-S1.jsonl")
    a1, v1, _ = arrays(s1)
    idx1 = rng.integers(0, len(a1), (2000, len(a1)))
    r1 = {"point": a1.mean() / v1.mean(), "boot": a1[idx1].mean(1) / v1[idx1].mean(1)}
    gens = {"FD-S1": summarise(s1, None, rng)}
    gens["FD-S1"]["recovery_ratio_r1"] = float(r1["point"])
    gens["FD-S1"]["recovery_ratio_r1_ci"] = R.pct_ci(r1["boot"])
    for nm in ("FD-S8", "FD-S9-25", "FD-S9-50", "FD-S9-75"):
        gens[nm] = summarise(R.read_jsonl(R.AN / f"rev_generators_{nm}.jsonl"), r1, rng)
        gens[nm]["post_hoc_law"] = SPECS[nm][5]
    infos = {r["generator"]: r for r in R.read_jsonl(cal_path)}
    for j, sh in SHARES.items():
        nm = f"FD-S8-share{int(sh * 100)}"
        gens[nm] = summarise(R.read_jsonl(R.AN / f"rev_r2_share_{int(sh * 100)}.jsonl"), r1 if sh < 1 else None, rng)
        gens[nm]["calibration"] = {k: infos[nm][k] for k in ("share", "lambda", "amp", "sigma", "shock_s")}
        if sh == 1.0:
            gens[nm]["note"] = "no reporting shock (share 1.0): v_H = 0, reference for the NB(k = 0.3), rho = 0.25 branching law"
    shock = [nm for nm, g in gens.items() if g["mean_vH"] > 0]
    rule_reflects = all(gens[nm]["excludes_zero"] for nm in shock)
    tg = [gens[nm]["t_g"] for nm in shock]
    res = {"meta": {"script": "code/rev_r2_total_bias.py", "plan": "analysis/analysis_plan_amendment2.md, Section C",
                    "panel_bootstrap": {"B": 2000, "seed": [R.SEED, 71, 1]},
                    "share_sweep": {"law": "FD-S8 (NB k = 0.3 branching, rho = 0.25, scale exponent 0.25)",
                                    "calibration_stream": "910 + j", "panel_seed": [R.SEED, 41, "j"], "reps": args.reps,
                                    "shares": list(SHARES.values())},
                    "seconds": None},
           "generators": gens,
           "generators_with_reporting_shock": shock,
           "wording_rule": {"rule": "if the relative-bias interval excludes 0 under every generator with a reporting shock, "
                                    "the intercept 'reflects' transmission and reporting variability (reporting with weight w_H) "
                                    "and is not called an estimate of the total",
                            "all_exclude_zero": rule_reflects,
                            "verdict": "reflects (not an estimate of the total)" if rule_reflects else "total-estimate wording not ruled out by this rule"},
           "mechanism_conditional_rescaling": {"t_g_range": [float(min(tg)), float(max(tg))],
                                               "a_obs_over_t_g_range": [float(OBS_A / max(tg)), float(OBS_A / min(tg))],
                                               "note": "not used in any headline number"}}
    res["meta"]["seconds"] = round(time.time() - t0, 1)
    R.write_json(OUT, res)
    log("done", res["wording_rule"]["verdict"], "t_g", res["mechanism_conditional_rescaling"]["t_g_range"])


if __name__ == "__main__":
    main()
