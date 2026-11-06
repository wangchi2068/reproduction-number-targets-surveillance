"""Amendment 3, Section A3.2: conditioning of the b/b_P interval on the reference (review N12).

Recomputes b/b_P under the three reference definitions already reported in Table 12 with the
calendar-week block bootstrap of rev_intercept.py (B = 1,999, seed [20261105, 22]):

  shipped_plugin    b_P = mean[kappa(g)^2 (1 + 1/g)]                     (s03b_thinning.poisson_reference)
  normalised_plugin b_P = mean[kappa(g)^2 (1 + 1/g) (x/xbar)^2]
  simulated_exact   b_P from the Poisson simulation of rev_r2_poisson_ref.py, held at its pooled
                    value with a lognormal Monte Carlo perturbation whose relative SD is read from
                    the reported per-panel 2.5-97.5% slope range

For each reference the script reports the point value, the 95% percentile interval over bootstrap
draws and the share of draws with b/b_P <= 1. No parametric joint reference is claimed.

Usage: python code/rev_r3_reference.py [--B 1999]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time

import numpy as np

import rev_common as R
from revision_data import load_sample
from revision_stats import fit_affine, relative_loss
from s03b_thinning import poisson_reference

OUT = R.AN / "rev_r3_reference.json"
HASHES = R.AN / "revision3" / "script_hashes.jsonl"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=1999)
    a = ap.parse_args()
    t0 = time.time()
    HASHES.parent.mkdir(parents=True, exist_ok=True)
    R.append_jsonl(HASHES, {"script": "rev_r3_reference.py",
                            "sha256": hashlib.sha256(open(__file__, "rb").read()).hexdigest(),
                            "utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),
                            "note": "Amendment 3, Section A3.2"})
    f = load_sample("RA", 3.6)
    x, D = f.x.to_numpy(float), f.D.to_numpy(float)
    ct, cn = f.C_t.to_numpy(float), f.C_next.to_numpy(float)
    loss, xbar = relative_loss(x)
    fa = fit_affine(loss, D)
    b = float(fa["b"])
    g = (cn + .5) / (ct + .5)
    k2 = ((3.6 / 7) / (1 + (1.6 ** 2 / 3.6) * np.log(g) / 7)) ** 2
    norm_term = k2 * (1 + 1 / g) * (x / xbar) ** 2

    res = json.loads((R.AN / "results.json").read_text(encoding="utf-8"))
    p2 = res["revision_2026_11_r2"]["F2_poisson_reference"]
    bP_shipped = float(p2["references"]["shipped_plugin"])
    bP_norm = float(p2["references"]["normalised_plugin"])
    bP_sim = float(p2["references"]["simulated_exact"])
    q_lo, q_hi = (float(v) for v in p2["simulated"]["per_panel_slope_range"])
    sd_rel = ((q_hi - q_lo) / 3.92) / bP_sim

    rng = np.random.default_rng([R.SEED, 22])
    wk = f.week_end.to_numpy()
    draws = {k: [] for k in ("b_over_bP_shipped", "b_over_bP_norm", "b_over_bP_sim")}
    for r in range(a.B):
        w = R.week_multiplicity(wk, rng)
        fb = fit_affine(loss, D, w)
        if not fb["converged"]:
            continue
        br = float(fb["b"])
        draws["b_over_bP_shipped"].append(br / poisson_reference(f, w))
        draws["b_over_bP_norm"].append(br / float(np.average(norm_term, weights=w)))
        draws["b_over_bP_sim"].append(br / (bP_sim * np.exp(rng.normal(0, sd_rel))))
        if r % 500 == 499:
            print("draw", r + 1, flush=True)

    def summary(v):
        arr = np.asarray(v, float)
        arr = arr[np.isfinite(arr)]
        return {"ci": R.pct_ci(arr), "share_le_1": float(np.mean(arr <= 1)), "n": int(len(arr))}

    points = {"shipped_plugin": b / bP_shipped, "normalised_plugin": b / bP_norm, "simulated_exact": b / bP_sim}
    out = {
        "meta": {"script": "code/rev_r3_reference.py", "plan": "analysis/analysis_plan_amendment3.md, Section A3.2",
                 "B": a.B, "seed": [R.SEED, 22], "b_hat": b,
                 "references": {"shipped_plugin": bP_shipped, "normalised_plugin": bP_norm, "simulated_exact": bP_sim},
                 "simulated_mc_sd_relative": sd_rel,
                 "resampling": "calendar weeks with replacement; b and both plug-in references recomputed per draw; "
                               "the simulated reference is held at its pooled value with a lognormal MC perturbation",
                 "seconds": round(time.time() - t0, 1)},
        "point": points,
        "interval": {"shipped_plugin": summary(draws["b_over_bP_shipped"]),
                     "normalised_plugin": summary(draws["b_over_bP_norm"]),
                     "simulated_exact": summary(draws["b_over_bP_sim"])},
        "lower_limit_exceeds_1_all_references": all(
            summary(draws[k])["ci"][0] > 1 for k in draws),
        "point_range": [min(points.values()), max(points.values())],
    }
    R.write_json(OUT, out)
    print(json.dumps({"point": points, "lower_all": out["lower_limit_exceeds_1_all_references"]}, indent=1))


if __name__ == "__main__":
    main()
