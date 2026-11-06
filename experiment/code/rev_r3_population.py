"""Amendment 3, Section A3.1: population gradient under a free count exponent (review N10).

M1  E(L | C, z) = a0 + a1 z + b / C
M2  E(L | C, z) = a0 + a1 z + b C^-gamma,  gamma > 0

z is the standardised log mean population of Amendment 2, Section B.2.
Both models are fitted by gamma quasi-likelihood. Inference is the two-way (state x calendar week)
CGM sandwich of E1, t with G - 1 = 41 df, conditional on the profiled gamma in M2. A state-cluster
bootstrap (B = 1,999, seed [.., 80, 1]) refits both models and gives intervals for a1 and for the
P90/P10 intercept ratio that do not condition on gamma.

The matched-count table of review N10 is recomputed here with state-level population tertiles
(17/17/17) and net loss L - b/C, with a state-cluster bootstrap (B = 999, seed [.., 80, 2]) on the
small-minus-large difference.

Usage: python code/rev_r3_population.py [--B 1999] [--B-matched 999]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time

import numpy as np
import pandas as pd
from scipy import stats

import ma_calibration as MC
import rev_common as R
import rev_r2_lib as L2
from revision_data import load_sample
from revision_stats import relative_loss

OUT = R.AN / "rev_r3_population.json"
HASHES = R.AN / "revision3" / "script_hashes.jsonl"
BANDS = ((20, 49), (50, 99), (100, 199), (200, 399))
LOG = None


def log(*a):
    line = " ".join(str(x) for x in a)
    print(line, flush=True)
    if LOG is not None:
        with open(LOG, "a", encoding="utf-8") as fh:
            fh.write(line + "\n")


def design(g, z, D, free):
    return np.column_stack([np.ones(len(D)), z, D ** (-g) if free else 1 / D])


def profile_three(y, z, D, w=None):
    """Profile gamma for E(y) = a0 + a1 z + b D^-gamma; returns the fit at the profiled gamma."""
    def prof(g):
        X = design(g, z, D, True)
        f = L2.fit_linear_ql(y, X, w)
        if not f["converged"]:
            return np.inf, None
        return L2.qdev(y, X @ f["beta"], w), f
    grid = np.concatenate([np.linspace(.05, 1.5, 30), np.linspace(1.6, 3.0, 8)])
    vals = [prof(g)[0] for g in grid]
    if not np.any(np.isfinite(vals)):
        return None
    j = int(np.nanargmin(vals))
    lo, hi = grid[max(j - 1, 0)], grid[min(j + 1, len(grid) - 1)]
    gr = (np.sqrt(5) - 1) / 2
    c1, c2 = hi - gr * (hi - lo), lo + gr * (hi - lo)
    f1, f2 = prof(c1)[0], prof(c2)[0]
    for _ in range(40):
        if f1 < f2:
            hi, c2, f2 = c2, c1, f1
            c1 = hi - gr * (hi - lo)
            f1 = prof(c1)[0]
        else:
            lo, c1, f1 = c1, c2, f2
            c2 = lo + gr * (hi - lo)
            f2 = prof(c2)[0]
    g = (lo + hi) / 2
    q, f = prof(g)
    if f is None:
        return None
    return {"gamma": float(g), "beta": f["beta"], "qdev": float(q), "X": design(g, z, D, True)}


def fit_model(y, z, D, free, w=None):
    if free:
        return profile_three(y, z, D, w)
    f = L2.fit_linear_ql(y, design(1.0, z, D, False), w)
    if not f["converged"]:
        return None
    return {"gamma": 1.0, "beta": f["beta"], "qdev": float(L2.qdev(y, design(1.0, z, D, False) @ f["beta"], w)),
            "X": design(1.0, z, D, False)}


def inference(fit, x, xbar, loss, st, wk, z10, z90):
    """Two-way CGM sandwich conditional on the fitted gamma, and the implied intercept ratio."""
    beta = fit["beta"]
    V, mu = R.cgm_from_loss(x, xbar, loss, None, fit["X"], beta, st, wk)
    if V is None:
        return None
    se = np.sqrt(np.diag(V))
    G = min(len(np.unique(st)), len(np.unique(wk)))
    df = G - 1
    tq = stats.t.ppf(.975, df)
    out = {}
    for i, lab in enumerate(("a0", "a1", "b")):
        t = float(beta[i] / se[i]) if se[i] > 0 else np.nan
        out[lab] = {"estimate": float(beta[i]), "se": float(se[i]), "t": t,
                    "p_two_sided": float(2 * stats.t.sf(abs(t), df)),
                    "ci": [float(beta[i] - tq * se[i]), float(beta[i] + tq * se[i])]}
    a0, a1 = float(beta[0]), float(beta[1])
    out["ratio_p90_over_p10"] = {"point": float((a0 + a1 * z90) / (a0 + a1 * z10)),
                                 "a_at_p10": float(a0 + a1 * z10), "a_at_p90": float(a0 + a1 * z90),
                                 "z_p10": float(z10), "z_p90": float(z90)}
    out["gamma"] = float(fit["gamma"])
    out["df"] = int(df)
    return out


def state_multiplicity(states, rng):
    u, inv = np.unique(states, return_inverse=True)
    m = np.zeros(len(u))
    np.add.at(m, rng.integers(0, len(u), len(u)), 1)
    return m[inv]


def matched_table(e, loss, b, z, B, rng):
    """Mean net loss L - b/C by state population tertile and count band (review N10)."""
    order = np.argsort(z)
    tert = np.empty(len(z), dtype=object)
    tert[order[:17]] = "small"
    tert[order[17:34]] = "mid"
    tert[order[34:]] = "large"
    rows = []
    for lo, hi in BANDS:
        mask = (e.D.to_numpy() >= lo) & (e.D.to_numpy() <= hi)
        rec = {"band": f"{lo}-{hi}", "n": int(mask.sum()), "n_states": {}}
        means = {}
        for k in ("small", "mid", "large"):
            mk = mask & (tert == k)
            rec["n_states"][k] = int(len(np.unique(e.state.to_numpy()[mk])))
            means[k] = float(np.mean((loss - b / e.D.to_numpy())[mk])) if mk.any() else None
            rec[k] = means[k]
        # state-cluster bootstrap on small - large
        st_all = e.state.to_numpy()
        net = loss - b / e.D.to_numpy()
        draws = []
        states = np.unique(st_all)
        for _ in range(B):
            pick = states[rng.integers(0, len(states), len(states))]
            idx = np.concatenate([np.flatnonzero(st_all == s) for s in pick])
            sub_mask = mask[idx]
            sub_t = tert[idx]
            sm = sub_mask & (sub_t == "small")
            lg = sub_mask & (sub_t == "large")
            if sm.sum() < 3 or lg.sum() < 3:
                continue
            draws.append(float(net[idx][sm].mean() - net[idx][lg].mean()))
        rec["small_minus_large"] = float(means["small"] - means["large"]) if means["small"] is not None and means["large"] is not None else None
        rec["small_minus_large_ci"] = R.pct_ci(draws) if draws else [None, None]
        rec["n_boot"] = len(draws)
        rows.append(rec)
    return rows


def main():
    global LOG
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=1999)
    ap.add_argument("--B-matched", type=int, default=999)
    a = ap.parse_args()
    LOG = R.LOGS / "rev_r3_population.log"
    t0 = time.time()
    HASHES.parent.mkdir(parents=True, exist_ok=True)
    R.append_jsonl(HASHES, {"script": "rev_r3_population.py",
                            "sha256": hashlib.sha256(open(__file__, "rb").read()).hexdigest(),
                            "utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),
                            "note": "Amendment 3, Section A3.1"})

    e = load_sample("RA")
    pop = L2.population()
    e = e.merge(pop[["state", "z_pop"]], on="state", how="left")
    assert e.z_pop.notna().all() and len(e) == 1748
    z = e.z_pop.to_numpy(float)
    D = e.D.to_numpy(float)
    st, wk = e.state.to_numpy(), e.week_end.to_numpy()
    x = e.x.to_numpy(float)
    loss, xbar = relative_loss(x)
    zs = e.groupby("state").z_pop.first().to_numpy()
    z10, z90 = (float(v) for v in np.quantile(zs, [.1, .9]))

    m1 = fit_model(loss, z, D, False)
    m2 = fit_model(loss, z, D, True)
    r1 = inference(m1, x, xbar, loss, st, wk, z10, z90)
    r2 = inference(m2, x, xbar, loss, st, wk, z10, z90)
    log("M1 a1", r1["a1"], "ratio", r1["ratio_p90_over_p10"]["point"])
    log("M2 gamma", r2["gamma"], "a1", r2["a1"], "ratio", r2["ratio_p90_over_p10"]["point"])

    # state-cluster bootstrap: refit both models
    rng = np.random.default_rng([R.SEED, 80, 1])
    draws = {k: [] for k in ("m1_a1", "m2_a1", "m1_ratio", "m2_ratio", "m2_gamma")}
    for r in range(a.B):
        w = state_multiplicity(st, rng)
        for tag, free in (("m1", False), ("m2", True)):
            f = fit_model(loss, z, D, free, w)
            if f is None:
                continue
            b0, b1 = float(f["beta"][0]), float(f["beta"][1])
            draws[f"{tag}_a1"].append(b1)
            draws[f"{tag}_ratio"].append((b0 + b1 * z90) / (b0 + b1 * z10))
            if free:
                draws[f"{tag}_gamma"].append(f["gamma"])
        if r % 200 == 199:
            log("bootstrap", r + 1, "of", a.B)
    boot = {k: {"ci": R.pct_ci(v), "n": len(v)} for k, v in draws.items()}

    matched = matched_table(e, loss, float(m1["beta"][2]), z, a.B_matched,
                            np.random.default_rng([R.SEED, 80, 2]))

    verdict = ("robust" if (r2["a1"]["estimate"] < 0 and r2["a1"]["ci"][1] < 0) else
               "not robust to the free exponent")
    out = {
        "meta": {"script": "code/rev_r3_population.py", "plan": "analysis/analysis_plan_amendment3.md, Section A3.1",
                 "data": "frozen NHSN 2026-10-20 vintage, R_A sample, 1,748 rows",
                 "seeds": {"bootstrap": [R.SEED, 80, 1], "matched": [R.SEED, 80, 2]},
                 "B": a.B, "B_matched": a.B_matched, "seconds": round(time.time() - t0, 1)},
        "M1_affine_population": r1, "M2_free_exponent_population": r2,
        "bootstrap": boot,
        "matched_count_table": matched,
        "decision_rule": {"rule": "gradient robust iff a1 < 0 and its 95% interval excludes 0 in M2",
                          "verdict": verdict},
    }
    R.write_json(OUT, out)
    log("done", verdict, "seconds", round(time.time() - t0, 1))


if __name__ == "__main__":
    main()
