"""Sensitivity of the E1 wild bootstrap-t to the loss clamp (analysis_plan_amendment2.md, Section E; review N8).

Same week-cluster wild bootstrap-t as rev_intercept.wild_bootstrap_t (a = 0 imposed, Webb weights per calendar
week, seed [20261105, 21], B = 9,999, so every variant sees the same weight draws), under clamp thresholds
1e-9, 1e-7 (shipped), 1e-5, 1e-3; no clamp (negative bootstrap losses allowed, linear gamma QL with positive
fitted means only); and the shipped 1e-7 run restricted to bootstrap samples with no clamped loss.

Output: analysis/rev_r2_clamp.json. Usage: python code/rev_r2_clamp.py [--B 9999] [--workers 4]
"""
from __future__ import annotations

import argparse
import time
from multiprocessing import Pool

import numpy as np

import rev_common as R
import rev_r2_lib as L2
from ma_models import Log
from revision_data import load_sample
from revision_stats import fit_affine, relative_loss
from rev_intercept import fit as fit_sample

OUT = R.AN / "rev_r2_clamp.json"
VARIANTS = ("1e-09", "1e-07", "1e-05", "1e-03", "none")


def run_variant(job):
    var, B = job
    f = load_sample("RA", 3.6)
    x, D = f.x.to_numpy(float), f.D.to_numpy(float)
    st, wk = f.state.to_numpy(), f.week_end.to_numpy()
    loss, xbar = relative_loss(x)
    b0 = float(np.mean(D * loss))
    mu0 = b0 / D
    _, wi = np.unique(wk, return_inverse=True)
    X = np.column_stack([np.ones(len(x)), 1 / D])
    rng = np.random.default_rng([R.SEED, 21])
    log = Log(R.LOGS / f"rev_r2_clamp_{var}.log")
    ts = np.full(B, np.nan)
    nclamp = np.zeros(B, int)
    for r in range(B):
        v = R.webb_weights(rng, wi.max() + 1)[wi]
        Lb = mu0 + v * (loss - mu0)
        if var == "none":
            fl = L2.fit_linear_ql(Lb, X)
            beta = fl["beta"] if fl["converged"] else None
            nclamp[r] = int(np.sum(Lb < 0))
        else:
            c = float(var)
            nclamp[r] = int(np.sum(Lb < c))
            Lb = np.maximum(Lb, c)
            fa = fit_affine(Lb, D)
            beta = np.array([fa["a"], fa["b"]]) if fa["converged"] else None
        if beta is None:
            continue
        V, _ = R.cgm_from_loss(x, xbar, Lb, D, X, beta, st, wk)
        if V is None or V[0, 0] <= 0:
            continue
        ts[r] = beta[0] / np.sqrt(V[0, 0])
        if r % 2000 == 1999:
            log(var, r + 1)
    return var, ts, nclamp


def pvalue(ts, t_obs):
    ok = np.isfinite(ts)
    if ok.sum() < 100:
        return {"valid": int(ok.sum()), "evaluable": False}
    p = float((1 + np.sum(np.abs(ts[ok]) >= abs(t_obs))) / (ok.sum() + 1))
    return {"valid": int(ok.sum()), "evaluable": True, "p_two_sided": p, "mcse": float(np.sqrt(p * (1 - p) / ok.sum())),
            "p_min_attainable": float(1 / (ok.sum() + 1)),
            "t_star_q": [float(v) for v in np.quantile(ts[ok], [.025, .5, .975])],
            "share_abs_t_star_ge_t_obs": float(np.mean(np.abs(ts[ok]) >= abs(t_obs)))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=9999)
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    L2.register_script(__file__)
    t0 = time.time()
    log = Log(R.LOGS / "rev_r2_clamp.log")
    f = load_sample("RA", 3.6)
    t_obs = fit_sample(f)["inference"]["two_way"]["t"]
    n = len(f)
    with Pool(args.workers) as pool:
        out = {v: (ts, nc) for v, ts, nc in pool.map(run_variant, [(v, args.B) for v in VARIANTS])}
    rows = {}
    for v in VARIANTS:
        ts, nc = out[v]
        key = "no_clamp" if v == "none" else f"clamp_{v}"
        rows[key] = {**pvalue(ts, t_obs), "row_share_clamped": float(nc.sum() / (args.B * n)),
                     "sample_share_with_clamp": float(np.mean(nc > 0))}
        if v == "none":
            rows[key]["note"] = "no clamp: row and sample shares count negative bootstrap losses (kept as they are)"
    ts7, nc7 = out["1e-07"]
    excl = ts7.copy()
    excl[nc7 > 0] = np.nan
    rows["exclude_clamped_draws"] = {**pvalue(excl, t_obs), "draws_without_clamp": int(np.sum(nc7 == 0))}
    ps = [r["p_two_sided"] for r in rows.values() if r.get("evaluable")]
    res = {"meta": {"script": "code/rev_r2_clamp.py", "plan": "analysis/analysis_plan_amendment2.md, Section E",
                    "B": args.B, "seed": [R.SEED, 21], "t_obs": t_obs, "seconds": round(time.time() - t0, 1),
                    "note": "secondary calibration of E1; no outcome changes the E1 verdict"},
           "variants": rows, "p_range_over_evaluable_variants": [float(min(ps)), float(max(ps))]}
    R.write_json(OUT, res)
    log("done", res["p_range_over_evaluable_variants"])


if __name__ == "__main__":
    main()
