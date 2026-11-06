"""Normalisation of the Poisson counting reference b_P (analysis_plan_amendment2.md, Section F.2; review minor 4).

Shipped plug-in: b_P = mean[kappa(g)^2 (1 + 1/g)] (s03b_thinning.poisson_reference; row R_A / xbar set to 1).
Normalised plug-in: b_P,norm = mean[kappa(g)^2 (1 + 1/g) (x_i / xbar)^2].
Simulated exact reference: B = 400 panels of independent Poisson counts with means equal to the observed C_t and
C_{t+1} (seed [20261105, 74]); rows with simulated C_t >= 5; L_cnt = ((x_sim - x_mu) / xbar_sim)^2 with x_mu the
transform of the observed (mean) ratio; b_P,sim = slope of the affine gamma-QL fit of L_cnt on simulated C_t,
pooled over panels, with the MC 2.5-97.5% range of the per-panel slope and the pooled intercept.

Output: analysis/rev_r2_poisson_ref.json. Usage: python code/rev_r2_poisson_ref.py [--B 400]
"""
from __future__ import annotations

import argparse
import time

import numpy as np

import ma_audit as A
import rev_common as R
import rev_r2_lib as L2
from ma_models import Log
from revision_data import load_sample
from revision_stats import fit_affine, relative_loss
from s03b_thinning import poisson_reference

OUT = R.AN / "rev_r2_poisson_ref.json"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=400)
    args = ap.parse_args()
    L2.register_script(__file__)
    t0 = time.time()
    log = Log(R.LOGS / "rev_r2_poisson_ref.log")
    f = load_sample("RA", 3.6)
    x, D = f.x.to_numpy(float), f.D.to_numpy(float)
    ct, cn = f.C_t.to_numpy(float), f.C_next.to_numpy(float)
    loss, xbar = relative_loss(x)
    fa = fit_affine(loss, D)
    b = fa["b"]
    bP = poisson_reference(f)
    g = (cn + .5) / (ct + .5)
    k2 = A.kappa(g) ** 2
    bP_norm = float(np.mean(k2 * (1 + 1 / g) * (x / xbar) ** 2))
    rng = np.random.default_rng([R.SEED, 74])
    x_mu = R.wl(g)
    Ls, Ds, slopes, inters, kept = [], [], [], [], []
    for r in range(args.B):
        cs, ns = rng.poisson(ct), rng.poisson(cn)
        keep = cs >= A.MIN_COUNT
        xs = R.wl((ns[keep] + .5) / (cs[keep] + .5))
        L = ((xs - x_mu[keep]) / xs.mean()) ** 2
        Ds.append(cs[keep].astype(float))
        Ls.append(L)
        fp = fit_affine(L, cs[keep].astype(float))
        if fp["converged"]:
            slopes.append(fp["b"])
            inters.append(fp["a"])
        kept.append(int(keep.sum()))
        if r % 100 == 99:
            log("panels", r + 1)
    Lp, Dp = np.concatenate(Ls), np.concatenate(Ds)
    fp = fit_affine(Lp, Dp)
    bP_sim = float(fp["b"])
    ref = {"shipped_plugin": bP, "normalised_plugin": bP_norm, "simulated_exact": bP_sim}
    ratio = {k: b / v for k, v in ref.items()}
    shipped_ratio = 1.25682
    res = {"meta": {"script": "code/rev_r2_poisson_ref.py", "plan": "analysis/analysis_plan_amendment2.md, Section F.2",
                    "B": args.B, "seed": [R.SEED, 74], "min_count": A.MIN_COUNT, "seconds": None},
           "b_hat": b, "references": ref, "b_over_bP": ratio,
           "relative_change_vs_shipped_ratio": {k: v / shipped_ratio - 1 for k, v in ratio.items()},
           "simulated": {"pooled_slope": bP_sim, "pooled_intercept": float(fp["a"]),
                         "per_panel_slope_range": [float(v) for v in np.quantile(slopes, [.025, .975])],
                         "per_panel_intercept_range": [float(v) for v in np.quantile(inters, [.025, .975])],
                         "mean_rows_kept": float(np.mean(kept)), "panels_converged": len(slopes)},
           "b_ci_scaled_note": "b/b_P intervals under each reference are the shipped b intervals divided by the reference "
                               "(sandwich, b_P held fixed); rendered in rev_r2_render.py from E1 / counting_slope_and_crossover"}
    res["meta"]["seconds"] = round(time.time() - t0, 1)
    R.write_json(OUT, res)
    log("done", ref, ratio)


if __name__ == "__main__":
    main()
