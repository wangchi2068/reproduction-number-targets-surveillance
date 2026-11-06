"""Daily-renewal generator: two generations per week (analysis plan Section 5.4, review M5).

On the real R_A design (state runs of consecutive weeks), daily infections follow a renewal process
I_d ~ Poisson(R_w(d) * P_d) or, with offspring dispersion k, Poisson(Gamma(k P_d, R_w(d) / k)), where
P_d = sum_s w_s I_(d-s) and w is the gamma(mean 3.6 d, SD 1.6 d) generation interval discretised to days
1..21. Weekly reproduction numbers R_j = wl_R(g0 exp(sigma z_j - sigma^2 / 2)) are drawn afresh in each
replicate with the S1 law (iid z, sigma as calibrated for S1). Each run is seeded with 35 days of
exponential growth at the geometric-mean weekly growth g0 and a level anchored like the shipped design;
weekly reported counts are Binomial(sum of the week's infections, rho). lambda (global level) is
calibrated so that the median retained count is 45 (ma_audit.fit_lambda).

Estimands per scored row (key c -> key f): v_R = relative variance of R_f over retained rows; v_Rbar = relative
variance of (R_c + R_f) / 2. Reported: a / v_R, a / v_Rbar, b / b_P, lag-1 correlation of x.
Cells: rho in {0.003, 0.01, 0.1} x offspring {Poisson, NB k = 0.3}; 200 panels each; seed [20261105, 60, i].

Usage: python code/rev_daily_renewal.py [--reps 200] [--workers 6]
"""
from __future__ import annotations

import argparse
import json
import time
from multiprocessing import Pool

import numpy as np
from scipy import stats

import ma_audit as A
import rev_common as R
from ma_models import Log
from revision_stats import fit_loss

CELLS = [(rho, k) for rho in (0.003, 0.01, 0.1) for k in (np.inf, 0.3)]
GI_MAX, SEED_DAYS, CHUNK = 21, 35, 25
OUT_JSONL = R.AN / "rev_daily_renewal.jsonl"
OUT_JSON = R.AN / "rev_daily_renewal.json"


def gi_weights():
    shape, scale = (A.GI_MEAN / A.GI_SD) ** 2, A.GI_SD ** 2 / A.GI_MEAN
    s = np.arange(1, GI_MAX + 1)
    w = stats.gamma.cdf(s + .5, shape, scale=scale) - stats.gamma.cdf(s - .5, shape, scale=scale)
    return w / w.sum()


class Renewal:
    def __init__(self, d, g0, sigma, rho, k):
        self.d, self.g0, self.sigma, self.rho, self.k = d, g0, sigma, rho, k
        self.lam, self.amp = 1.0, 1.0
        self.w = gi_weights()
        self.run_of_key = d.run
        self.pos = np.zeros(len(d.keys), int)
        for j in range(1, len(d.keys)):
            self.pos[j] = 0 if d.start[j] else self.pos[j - 1] + 1
        self.n_runs, self.L = int(d.run.max() + 1), int(self.pos.max() + 1)
        self.level0 = d.mean_path(np.full(d.n, g0), g0, 1.0)[0][d.start]  # expected count of each run's first week
        r0 = np.log(g0) / 7
        self.seed_shape = np.exp(r0 * np.arange(-SEED_DAYS, 0))
        self.week0_shape = np.exp(r0 * np.arange(7)).sum()
        self.last_R = None

    def set(self, amp=None, lam=None):
        self.lam = self.lam if lam is None else lam
        return self

    def draw(self, rng, reps):
        out, Rs = [], []
        for lo in range(0, reps, CHUNK):
            c, Rk = self._draw(rng, min(CHUNK, reps - lo))
            out.append(c)
            Rs.append(Rk)
        self.last_R = np.concatenate(Rs)
        return np.concatenate(out)

    def _draw(self, rng, reps):
        d, nr, L, w = self.d, self.n_runs, self.L, self.w
        nk = len(d.keys)
        z = rng.standard_normal((reps, nk))
        Rk = R.wl(self.g0 * np.exp(self.sigma * z - self.sigma ** 2 / 2))  # weekly R per key
        Rgrid = np.ones((reps, nr, L))
        Rgrid[:, self.run_of_key, self.pos] = Rk
        days = SEED_DAYS + 7 * L
        I = np.zeros((reps, nr, days))
        I0 = self.lam * self.level0 / self.rho / self.week0_shape
        I[:, :, :SEED_DAYS] = (I0[:, None] * self.seed_shape[None, :])[None]
        ww = w[::-1]
        for t in range(SEED_DAYS, days):
            lo = max(0, t - GI_MAX)
            P = I[:, :, lo:t] @ ww[GI_MAX - (t - lo):]
            mean = Rgrid[:, :, (t - SEED_DAYS) // 7] * P
            if np.isinf(self.k):
                I[:, :, t] = rng.poisson(mean)
            else:
                shape = np.maximum(self.k * P, 1e-300)
                I[:, :, t] = np.where(P > 0, rng.poisson(rng.gamma(shape, np.where(P > 0, mean / shape, 0.0))), 0)
        weekly = I[:, :, SEED_DAYS:].reshape(reps, nr, L, 7).sum(-1)
        lat = weekly[:, self.run_of_key, self.pos]
        counts = lat if self.rho >= 1 else rng.binomial(lat.astype(np.int64), self.rho).astype(float)
        return counts, Rk


def run_cell(task):
    i, rho, k, reps = task
    t0 = time.time()
    e, d, obs, _, gens = R.setup_design()
    s1 = R.restore(gens["S1_heterogeneity"], "S1_heterogeneity")
    gen = Renewal(d, obs["g0"], s1.sigma, rho, k)
    al = A.fit_lambda(gen, 600 + i)
    rng = np.random.default_rng([R.SEED, 60, i])
    counts = gen.draw(rng, reps)
    Rk = gen.last_R
    ct, cn = d.rows(counts)
    rows = {kk: np.full(reps, np.nan) for kk in ("a", "b", "bP", "vR", "vRbar", "lag1", "n", "D_median")}
    for r in range(reps):
        keep = ct[r] >= A.MIN_COUNT
        D = ct[r][keep]
        x = R.wl((cn[r][keep] + .5) / (ct[r][keep] + .5))
        f = fit_loss(x, D)
        if f["converged"]:
            rows["a"][r], rows["b"][r] = f["a"], f["b"]
        rows["bP"][r] = A.bp_of(D, cn[r][keep])
        Rf, Rc = Rk[r][d.fol][keep], Rk[r][d.cur][keep]
        rows["vR"][r] = R.rel_var(Rf)
        rows["vRbar"][r] = R.rel_var((Rc + Rf) / 2)
        xa = np.full(d.n, np.nan)
        xa[keep] = x
        m = keep & (d.next_row >= 0)
        m[m] &= keep[d.next_row[m]]
        rows["lag1"][r] = float(np.corrcoef(xa[m], xa[d.next_row[m]])[0, 1])
        rows["n"][r], rows["D_median"][r] = keep.sum(), np.median(D)
    rec = {"rho": rho, "k": None if np.isinf(k) else k, "offspring": "poisson" if np.isinf(k) else "negative_binomial",
           "lambda": gen.lam, "lambda_at_bound": al["at_bound"], "median_D_50panels": al["median_D_50panels"], "sigma_S1": s1.sigma,
           "a": R.mc(rows["a"]), "b": R.mc(rows["b"]), "b_over_bP": R.mc(rows["b"] / rows["bP"]),
           "a_over_vR": {"ratio_of_means": R.ratio_of_means(rows["a"], rows["vR"]), **R.mc(rows["a"] / rows["vR"])},
           "a_over_vRbar": {"ratio_of_means": R.ratio_of_means(rows["a"], rows["vRbar"]), **R.mc(rows["a"] / rows["vRbar"])},
           "vR": R.mc(rows["vR"]), "vRbar": R.mc(rows["vRbar"]), "lag1_corr_x": R.mc(rows["lag1"]),
           "n_rows": R.mc(rows["n"]), "D_median": R.mc(rows["D_median"]),
           "reps": reps, "seed": [R.SEED, 60, i], "seconds": round(time.time() - t0, 1)}
    R.append_jsonl(OUT_JSONL, rec)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    log = Log(R.LOGS / "rev_daily_renewal.log")
    if OUT_JSONL.exists():
        OUT_JSONL.unlink()
    tasks = [(i, rho, k, a.reps) for i, (rho, k) in enumerate(CELLS)]
    log("rev_daily_renewal start", len(tasks), "cells", "reps", a.reps, "GI weights", np.round(gi_weights(), 4).tolist())
    with Pool(min(a.workers, len(tasks))) as pool:
        for c in pool.imap_unordered(run_cell, tasks):
            log("cell rho", c["rho"], "k", c["k"], "lambda", c["lambda"], "a/vR", c["a_over_vR"]["ratio_of_means"],
                "a/vRbar", c["a_over_vRbar"]["ratio_of_means"], "b/bP", c["b_over_bP"]["mean"], "lag1", c["lag1_corr_x"]["mean"],
                "s", c["seconds"])
    cells = sorted(R.read_jsonl(OUT_JSONL), key=lambda c: (c["rho"], c["k"] is not None))
    rec = [c["a_over_vR"]["ratio_of_means"] for c in cells]
    out = {"meta": {"script": "code/rev_daily_renewal.py", "plan_section": "5.4", "reps": a.reps,
                    "gi": {"mean_days": A.GI_MEAN, "sd_days": A.GI_SD, "discretised_days": [1, GI_MAX]}, "seed": [R.SEED, 60, "i"],
                    "transmission": "weekly R_j redrawn each replicate, S1 law (iid lognormal), constant within the week"},
           "cells": cells,
           "SC1_one_generation_map_carries_over": {"rule": "a / v_R in [0.8, 1.25] in every cell", "range": [min(rec), max(rec)],
                                                   "pass": all(.8 <= v <= 1.25 for v in rec)},
           "a_over_vRbar_range": [min(c["a_over_vRbar"]["ratio_of_means"] for c in cells),
                                  max(c["a_over_vRbar"]["ratio_of_means"] for c in cells)],
           "b_over_bP_range": [min(c["b_over_bP"]["mean"] for c in cells), max(c["b_over_bP"]["mean"] for c in cells)]}
    R.write_json(OUT_JSON, out)
    log("rev_daily_renewal done", "SC1", out["SC1_one_generation_map_carries_over"])


if __name__ == "__main__":
    main()
