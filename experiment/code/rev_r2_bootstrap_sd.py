"""Pigeonhole bootstrap SD against the fresh-draw MC SD (analysis_plan_amendment2.md, Section F.1; review N9, minor 12).

The first 100 panels of each fresh-draw generator (FD-S1, FD-S8, FD-S9-25/50/75) are regenerated from the stored
calibration (analysis/rev_generators_calibration.jsonl) and panel stream [20261105, 40, i]; the regenerated a_hat
is checked against analysis/rev_generators_<name>.jsonl. Each panel gets B = 199 state x calendar-week pigeonhole
bootstrap draws (revision_stats.bootstrap_weights, L = 1, seed [20261105, 73, i, panel]).

Output: analysis/rev_r2_bootstrap_sd.json (+ per-panel rows analysis/rev_r2_bootstrap_sd_panels.jsonl).
Usage: python code/rev_r2_bootstrap_sd.py [--panels 100] [--B 199] [--workers 4]
"""
from __future__ import annotations

import argparse
import time
from multiprocessing import Pool

import numpy as np

import ma_audit as A
import ma_calibration as MC
import rev_common as R
import rev_r2_lib as L2
from ma_models import Log
from rev_generators import SPECS, make
from revision_stats import bootstrap_weights, fit_affine, relative_loss

OUT = R.AN / "rev_r2_bootstrap_sd.json"
PANELS = R.AN / "rev_r2_bootstrap_sd_panels.jsonl"
SHIPPED_FIXED_FIELD_RATIO = [1.6, 2.5]


def run(job):
    name, npan, B = job
    log = Log(R.LOGS / f"rev_r2_bootstrap_sd_{name}.log")
    e, d, obs, _, _ = R.setup_design()
    gen = make(name, d, obs)
    i = SPECS[name][0]
    info = {r["generator"]: r for r in R.read_jsonl(R.AN / "rev_generators_calibration.jsonl")}[name]
    gen.set(amp=info["amp"], lam=info["lambda"])
    sig_ok = bool(np.isclose(gen.sigma, info["sigma"], rtol=1e-8) and np.isclose(gen.shock_s, info["shock_s"], rtol=1e-8, atol=1e-12))
    stored = R.read_jsonl(R.AN / f"rev_generators_{name}.jsonl")
    rng = np.random.default_rng([R.SEED, 40, i])
    rows = []
    for p in range(npan):
        counts = gen.draw(rng, 1)[0]
        ct, cn = d.rows(counts)
        keep = ct >= A.MIN_COUNT
        D = ct[keep]
        x = R.wl((cn[keep] + .5) / (ct[keep] + .5))
        st, wk = d.state[keep], d.week[keep]
        s = MC.sandwich(x, D, st, wk)
        a_hat = float(s["beta"][0]) if s is not None else np.nan
        a_st = stored[p].get("a")
        fr = d.e.loc[keep, ["state", "week_end"]].reset_index(drop=True)
        brng = np.random.default_rng([R.SEED, 73, i, p])
        ab = np.full(B, np.nan)
        for r in range(B):
            w = bootstrap_weights(fr, brng, 1)[0]
            loss, _ = relative_loss(x, w)
            fa = fit_affine(loss, D, w)
            if fa["converged"]:
                ab[r] = fa["a"]
        row = {"generator": name, "panel": p, "a_regenerated": a_hat, "a_stored": a_st,
               "match": bool(a_st is not None and np.isclose(a_hat, a_st, rtol=1e-9)),
               "boot_sd": float(np.nanstd(ab, ddof=1)), "boot_ci": R.pct_ci(ab), "boot_valid": int(np.isfinite(ab).sum())}
        R.append_jsonl(PANELS, row)
        rows.append(row)
        if p % 20 == 19:
            log(name, "panels", p + 1)
    return name, rows, sig_ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", type=int, default=100)
    ap.add_argument("--B", type=int, default=199)
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    L2.register_script(__file__)
    t0 = time.time()
    log = Log(R.LOGS / "rev_r2_bootstrap_sd.log")
    if PANELS.exists():
        PANELS.unlink()
    names = list(SPECS)
    with Pool(args.workers) as pool:
        res_rows = pool.map(run, [(n, args.panels, args.B) for n in names])
    gens = {}
    for name, rows, sig_ok in res_rows:
        stored = R.read_jsonl(R.AN / f"rev_generators_{name}.jsonl")
        a_all = np.array([np.nan if r.get("a") is None else r["a"] for r in stored], float)
        se_all = np.array([np.nan if r.get("se_a") is None else r["se_a"] for r in stored], float)
        mc_sd = float(np.nanstd(a_all, ddof=1))
        ea = float(np.nanmean(a_all))
        bsd = np.array([r["boot_sd"] for r in rows])
        ci = np.array([r["boot_ci"] for r in rows], float)
        cov = float(np.mean((ci[:, 0] <= ea) & (ea <= ci[:, 1])))
        gens[name] = {"panels_bootstrapped": len(rows), "panels_matching_stored_a": int(sum(r["match"] for r in rows)),
                      "calibration_restored_exactly": sig_ok,
                      "max_abs_diff_regenerated_vs_stored": float(np.nanmax([abs(r["a_regenerated"] - r["a_stored"]) for r in rows if r["a_stored"] is not None])),
                      "mc_sd_a_200": mc_sd, "mean_boot_sd": float(bsd.mean()),
                      "boot_sd_over_mc_sd": float(bsd.mean() / mc_sd),
                      "boot_sd_over_mc_sd_range_panels": [float(v) for v in np.quantile(bsd / mc_sd, [.025, .975])],
                      "sandwich_se_over_mc_sd": float(np.nanmean(se_all) / mc_sd),
                      "percentile_coverage_of_E_a": cov, "coverage_mcse": float(np.sqrt(cov * (1 - cov) / len(rows))),
                      "post_hoc_law": SPECS[name][5]}
    ratios = {n: g["boot_sd_over_mc_sd"] for n, g in gens.items()}
    conservative = all(v > 1 for v in ratios.values())
    res = {"meta": {"script": "code/rev_r2_bootstrap_sd.py", "plan": "analysis/analysis_plan_amendment2.md, Section F.1",
                    "panels": args.panels, "B": args.B, "seed": [R.SEED, 73, "i", "panel"],
                    "seconds": round(time.time() - t0, 1)},
           "generators": gens,
           "rule": {"text": "the pigeonhole bootstrap is 'conservative' only if bootstrap SD / MC SD > 1 under every fresh-draw generator",
                    "ratios": ratios, "conservative": conservative,
                    "verdict": "conservative" if conservative else "'conservative' withdrawn",
                    "shipped_fixed_field_ratio": SHIPPED_FIXED_FIELD_RATIO, "sqrt3": float(np.sqrt(3))}}
    R.write_json(OUT, res)
    log("done", res["rule"]["verdict"], ratios)


if __name__ == "__main__":
    main()
