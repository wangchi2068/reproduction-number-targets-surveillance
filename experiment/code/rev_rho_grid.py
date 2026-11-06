"""Extended reporting-probability x offspring-dispersion grid (analysis plan Section 5.3, review M5).

rho in {0.001, 0.003, 0.01, 0.03, 0.1, 0.25, 0.5, 1} x k in {0.1, 0.3, 1, 3, inf}; S1 transmission field held
fixed as in the shipped grid (ma_audit.py); lambda recalibrated per cell to the observed median retained count;
100 panels per cell, seed [20261105, 50, i, j]. rho = 0.001-0.01 is an assumed range (CDC burden estimates), not
estimated. Each finished cell is appended to analysis/rev_rho_grid.jsonl; the summary goes to rev_rho_grid.json.

Usage: python code/rev_rho_grid.py [--reps 100] [--workers 6]
"""
from __future__ import annotations

import argparse
import time
from multiprocessing import Pool

import numpy as np

import ma_audit as A
import rev_common as R
from ma_models import Log

RHOS = (0.001, 0.003, 0.01, 0.03, 0.1, 0.25, 0.5, 1.0)
KS = (0.1, 0.3, 1.0, 3.0, np.inf)
ASSUMED_RHO = (0.001, 0.003, 0.01)
OUT_JSONL = R.AN / "rev_rho_grid.jsonl"
OUT_JSON = R.AN / "rev_rho_grid.json"
_W = {}


def _init():
    e, d, obs, z, gens = R.setup_design()
    s1 = R.restore(gens["S1_heterogeneity"], "S1_heterogeneity")
    _W.update(d=d, obs=obs, z=z, g_het=s1.g.copy(), sigma=s1.sigma)


def run_cell(task):
    i, j, rho, k, reps = task
    t0 = time.time()
    d, obs, z = _W["d"], _W["obs"], _W["z"]
    gen = A.Generator(d, obs["g0"], z, "none", k=k, rho=rho, g_fixed=_W["g_het"]).set(lam=1.0)
    al = A.fit_lambda(gen, 500 + 10 * i + j)
    o, _ = A.scenario_panels(gen, np.random.default_rng([R.SEED, 50, i, j]), reps, obs["contrast_affine_minus_power"], heavy=False)
    rec = o["a"] / o["truth_relvar_kept"]
    bb = o["b_over_bP"]
    q = R.mc(bb)
    return {"rho": rho, "k": None if np.isinf(k) else k, "offspring": "poisson" if np.isinf(k) else "negative_binomial",
            "assumed_rho_range": rho in ASSUMED_RHO, "lambda": gen.lam, "lambda_at_bound": al["at_bound"],
            "median_D_50panels": al["median_D_50panels"], "sigma_S1_field": _W["sigma"],
            "a": R.mc(o["a"]), "b": R.mc(o["b"]), "b_over_bP": q,
            "a_over_vR": {"ratio_of_means": R.ratio_of_means(o["a"], o["truth_relvar_kept"]), **R.mc(rec)},
            "truth_relvar_kept": R.mc(o["truth_relvar_kept"]), "D_median": R.mc(o["D_median"]), "n_rows": R.mc(o["n"]),
            "reproduces_counting_slope": bool(q["q025"] <= obs["b_over_bP"] <= q["q975"]),
            "reps": reps, "seed": [R.SEED, 50, i, j], "seconds": round(time.time() - t0, 1)}


def summarise(cells, obs):
    sc2 = [c for c in cells if not (.9 <= c["a_over_vR"]["ratio_of_means"] <= 1.1)]
    low = [c for c in cells if c["rho"] <= 0.01]
    sc3 = [c for c in low if not (.9 <= c["b_over_bP"]["mean"] <= 1.1)]
    return {
        "observed_b_over_bP": obs["b_over_bP"], "observed_a": obs["a"],
        "design": {"rho": list(RHOS), "k": [None if np.isinf(k) else k for k in KS], "n_cells": len(cells),
                   "transmission": "S1 heterogeneity field fixed across replicates (as in the shipped grid)",
                   "assumed_rho_range": {"values": list(ASSUMED_RHO),
                                         "label": "assumed range for P(infection -> reported laboratory-confirmed admission), from CDC burden "
                                                  "estimates; not estimated"}},
        "declared_rule": "a cell reproduces the counting slope if the MC 2.5-97.5% range of b/b_P contains the observed value",
        "cells_reproducing_counting_slope": [{"rho": c["rho"], "k": c["k"]} for c in cells if c["reproduces_counting_slope"]],
        "n_cells_reproducing": sum(c["reproduces_counting_slope"] for c in cells),
        "SC2_intercept_recovery_0.9_1.1_all_cells": {"pass": not sc2, "failing": [{"rho": c["rho"], "k": c["k"],
                                                                                  "a_over_vR": c["a_over_vR"]["ratio_of_means"]} for c in sc2],
                                                    "range": [min(c["a_over_vR"]["ratio_of_means"] for c in cells),
                                                              max(c["a_over_vR"]["ratio_of_means"] for c in cells)]},
        "SC3_b_over_bP_0.9_1.1_at_rho_le_0.01": {"pass": not sc3, "failing": [{"rho": c["rho"], "k": c["k"],
                                                                              "b_over_bP": c["b_over_bP"]["mean"]} for c in sc3],
                                                "range": [min(c["b_over_bP"]["mean"] for c in low), max(c["b_over_bP"]["mean"] for c in low)]},
        "cells": cells}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--fresh", action="store_true", help="ignore cells already in the JSONL")
    a = ap.parse_args()
    log = Log(R.LOGS / "rev_rho_grid.log")
    if a.fresh and OUT_JSONL.exists():
        OUT_JSONL.unlink()
    done = {(c["rho"], c["k"]) for c in R.read_jsonl(OUT_JSONL) if c["reps"] == a.reps}
    tasks = [(i, j, rho, k, a.reps) for i, rho in enumerate(RHOS) for j, k in enumerate(KS)
             if (rho, None if np.isinf(k) else k) not in done]
    # slowest (smallest rho, latent counts largest) first
    tasks.sort(key=lambda t: t[2])
    log("rev_rho_grid start", "cells to run", len(tasks), "done", len(done))
    with Pool(a.workers, initializer=_init) as pool:
        for cell in pool.imap_unordered(run_cell, tasks):
            R.append_jsonl(OUT_JSONL, cell)
            log("cell rho", cell["rho"], "k", cell["k"], "lambda", cell["lambda"], "a/vR", cell["a_over_vR"]["ratio_of_means"],
                "b/bP", cell["b_over_bP"]["mean"], "repro", cell["reproduces_counting_slope"], "s", cell["seconds"])
    cells = sorted((c for c in R.read_jsonl(OUT_JSONL) if c["reps"] == a.reps),
                   key=lambda c: (c["rho"], np.inf if c["k"] is None else c["k"]))
    _, _, obs, _, _ = R.setup_design()
    out = summarise(cells, obs)
    out["meta"] = {"script": "code/rev_rho_grid.py", "plan_section": "5.3", "reps_per_cell": a.reps, "seed": [R.SEED, 50, "i", "j"]}
    R.write_json(OUT_JSON, out)
    log("rev_rho_grid done", "SC2", out["SC2_intercept_recovery_0.9_1.1_all_cells"]["pass"],
        "SC3", out["SC3_b_over_bP_0.9_1.1_at_rho_le_0.01"]["pass"], "reproducing", out["n_cells_reproducing"])


if __name__ == "__main__":
    main()
