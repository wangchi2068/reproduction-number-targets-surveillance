"""State-realistic (SR) fresh-draw generators (analysis_plan_amendment2.md, Section D.1; review M4, M5).

Law SR, field redrawn every replicate:
  g_i = g0 (g_trend,i / g0)^tau exp(s_i z_i - s_i^2 / 2), s_i^2 = sigma^2 m_state, z_i iid N(0, 1);
  g_trend = S_fol / S_cur, S = three-point geometric moving average of the observed counts within each state run
  (two-point at run ends, counts floored at 0.5); mean path from Design.mean_path (run-anchored, times lambda).
  m_state (mean 1 over rows): "pop" proportional to max(a0 + a1 z_pop, 0.1 a0) with the B.2 point estimates,
  "fe" proportional to max(a_state, 0.1 mean a_state) from fit_group_affine on the observed losses.
  Latent NB(k = 0.3) branching started at mean / rho; reporting probability rho * shock with a lognormal
  state-week shock carrying half of the variance target, scaled by sqrt(m_state); clustered admissions
  c ~ Bin(latent, min(1, rho shock) / (1 + theta_c)), C = c + Poisson(theta_c c).
Calibration as the shipped and FD generators (ma_audit.calibrate, stream 920 + i); sigma such that the relative
variance of wl(g) over all rows on the reference field equals the transmission share of amp * a_obs; if the trend
alone exceeds the target the cell is flagged "calibration failed" (sigma = 0) and still summarised.
Cells i = 0..5, 200 panels, seed [20261105, 42, i]. Tests T1-T4 of Amendment 2, D.1.
Data-informed generator design: theta_c grid and state multipliers chosen knowing the observed counting slope and
between-state dispersion; no inferential result uses these generators.

Output: analysis/rev_r2_generator.json, panels analysis/rev_r2_generator_<cell>.jsonl.
Usage: python code/rev_r2_generator.py [--reps 200] [--workers 4]
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
import rev_r2_lib as L2
from ma_models import Log
from rev_generators import SE_OBS, cover, panel_row
from revision_data import load_sample
from revision_stats import relative_loss
from s03c_state_fe import fit_group_affine

OUT = R.AN / "rev_r2_generator.json"
CELLS = {  # name: (i, multiplier, theta_c, rho, tau)
    "SR-pop-c0": (0, "pop", 0.0, .01, 1.0),
    "SR-pop-c25": (1, "pop", .25, .01, 1.0),
    "SR-pop-c50": (2, "pop", .5, .01, 1.0),
    "SR-fe-c25": (3, "fe", .25, .01, 1.0),
    "SR-pop-c25-r003": (4, "pop", .25, .003, 1.0),
    "SR-pop-c25-t50": (5, "pop", .25, .01, .5),
}
SHARE, K = .5, A.K
OBS = {"b_over_bP": 1.25682, "p90": 388.6, "pooled_over_within": 1.476}


def smoothed_trend(d):
    """Per-row g_trend = S_fol / S_cur from 3-point geometric moving averages of observed counts within runs."""
    lc = np.log(np.maximum(d.obs_counts, .5))
    S = np.empty(len(lc))
    for r in np.unique(d.run):
        idx = np.flatnonzero(d.run == r)
        v = lc[idx]
        if len(v) == 1:
            S[idx] = v
            continue
        s = np.empty(len(v))
        s[0], s[-1] = v[:2].mean(), v[-2:].mean()
        if len(v) > 2:
            s[1:-1] = (v[:-2] + v[1:-1] + v[2:]) / 3
        S[idx] = s
    return np.exp(S[d.fol] - S[d.cur])


def state_multiplier(kind, d):
    st = d.state
    if kind == "pop":
        pop = L2.population().set_index("state")
        fl = json.loads((R.AN / "rev_r2_floor.json").read_text())["B2_population"]["primary_linear"]
        z = pop.loc[st, "z_pop"].to_numpy(float)
        m = np.maximum(fl["a0"] + fl["a1"] * z, .1 * fl["a0"])
    else:
        f = load_sample("RA", 3.6)
        loss, _ = relative_loss(f.x.to_numpy(float))
        fg = fit_group_affine(loss, f.D.to_numpy(float), f.state.to_numpy())
        ic = dict(zip(fg["groups"], np.asarray(fg["intercepts"], float)))
        mean_ic = float(np.nanmean(list(ic.values())))
        m = np.array([max(ic[s], .1 * mean_ic) for s in st])
    return m / m.mean()


class SR(A.Generator):
    def __init__(self, d, g0, kap, mult, theta_c, rho, tau):
        z0 = np.random.default_rng([A.SEED, 1]).standard_normal(d.n)
        super().__init__(d, g0, z0, "iid", k=K, rho=rho, share=SHARE, kap=kap)
        self.m_row = state_multiplier(mult, d)
        key_state = np.array([k[0] for k in d.keys])
        row_m = dict(zip(d.state, self.m_row))
        self.m_key = np.array([row_m.get(s, 1.0) for s in key_state])
        self.theta_c, self.tau = theta_c, tau
        self.base = g0 * (smoothed_trend(d) / g0) ** tau
        self.calibration_failed = False
        self.last_g, self.last_shocks = None, None

    def g_of(self, z, sigma=None):
        s = (self.sigma if sigma is None else sigma) * np.sqrt(self.m_row)
        return self.base * np.exp(s * z - s ** 2 / 2)

    def set(self, amp=None, lam=None):
        self.amp = self.amp if amp is None else amp
        self.lam = self.lam if lam is None else lam
        total = self.amp * A.GLOBAL["a_obs"]
        target = total * self.share
        f0 = A.rel_var(R.wl(self.base))
        if f0 >= target:
            self.sigma, self.calibration_failed = 0.0, True
        else:
            self.sigma = A.bisect(lambda s: A.rel_var(R.wl(self.g_of(self.z, s))) - target, 0.0, 3.0)
            self.calibration_failed = False
        self.shock_s = float(np.sqrt(total * (1 - self.share) / (2 * self.kap ** 2)))
        self.g = self.g_of(self.z)
        return self

    def draw(self, rng, reps):
        d = self.d
        out, gs, shocks = [], [], []
        for _ in range(reps):
            g = self.g_of(rng.standard_normal(d.n))
            mu, gt = d.mean_path(g, self.g0, self.lam)
            latent = np.empty(len(mu))
            for j in range(len(mu)):
                if d.start[j]:
                    latent[j] = rng.poisson(mu[j] / self.rho)
                else:
                    parent = latent[j - 1]
                    latent[j] = rng.poisson(rng.gamma(max(K * parent, 1e-300), gt[j - 1] / K)) if parent > 0 else 0
            s = self.shock_s * np.sqrt(self.m_key)
            zs = np.exp(s * rng.standard_normal(len(mu)) - s ** 2 / 2)
            c = rng.binomial(latent.astype(np.int64), np.minimum(1.0, self.rho * zs) / (1 + self.theta_c)).astype(float)
            if self.theta_c > 0:
                c = c + rng.poisson(self.theta_c * c)
            out.append(c)
            gs.append(g)
            shocks.append(zs)
        self.last_g, self.last_shocks = np.array(gs), shocks
        self.d.last_shock = None
        return np.array(out)


def state_cv(loss, D, st, w=None):
    fg = fit_group_affine(loss, D, st, w)
    if not fg["converged"]:
        return np.nan
    ic = np.asarray(fg["intercepts"], float)
    ic = ic[np.isfinite(ic)]
    return float(np.std(ic) / np.mean(ic))


def run_cell(job):
    name, reps = job
    i, mult, theta_c, rho, tau = CELLS[name]
    t0 = time.time()
    log = Log(R.LOGS / f"rev_r2_generator_{name}.log")
    e, d, obs, _, _ = R.setup_design()
    gen = SR(d, obs["g0"], obs["kappa_mean"], mult, theta_c, rho, tau)
    cal = A.calibrate(gen, 920 + i, log, name)
    log(name, "calibrated", "lambda", gen.lam, "amp", gen.amp, "sigma", gen.sigma, "shock_s", gen.shock_s,
        "failed", gen.calibration_failed)
    rng = np.random.default_rng([R.SEED, 42, i])
    season_all = d.e.season.astype(str).to_numpy()
    path = R.AN / f"rev_r2_generator_{name}.jsonl"
    if path.exists():
        path.unlink()
    for r in range(reps):
        counts = gen.draw(rng, 1)[0]
        row = panel_row(gen, counts, gen.last_g[0], gen.last_shocks[0], season_all, lambda df: stats.t.ppf(.975, df))
        ct, cn = d.rows(counts)
        keep = ct >= A.MIN_COUNT
        D = ct[keep]
        x = R.wl((cn[keep] + .5) / (ct[keep] + .5))
        loss, _ = relative_loss(x)
        row.update(generator=name, r=r, p10=float(np.quantile(D, .1)), median=float(np.median(D)),
                   p90=float(np.quantile(D, .9)), state_cv=state_cv(loss, D, d.state[keep]))
        R.append_jsonl(path, row)
        if r % 50 == 49:
            log(name, "panels", r + 1)
    return {"generator": name, "lambda": gen.lam, "amp": gen.amp, "sigma": gen.sigma, "shock_s": gen.shock_s,
            "calibration_failed": gen.calibration_failed, "trend_relvar": A.rel_var(R.wl(gen.base)),
            "calibration": cal, "seconds": round(time.time() - t0, 1)}


def rng_of(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return [float(np.quantile(v, .025)), float(np.quantile(v, .975))]


def summarise(name, rows, info, obs_cv):
    f = lambda k: np.array([np.nan if r.get(k) is None else r[k] for r in rows], float)
    a, b, bP, bw = f("a"), f("b"), f("bP"), f("bw")
    vR, vH = f("truth_relvar_kept"), f("vH_kept")
    sd_a = float(np.nanstd(a, ddof=1))
    t1, t3, t4 = rng_of(b / bP), rng_of(b / bw), rng_of(f("state_cv"))
    p90 = float(np.nanmean(f("p90")))
    tests = {"T1_b_over_bP": {"range": t1, "observed": OBS["b_over_bP"], "pass": t1[0] <= OBS["b_over_bP"] <= t1[1]},
             "T2_p90": {"mc_mean": p90, "band": [350, 428], "observed": OBS["p90"], "pass": 350 <= p90 <= 428},
             "T3_pooled_over_within": {"range": t3, "observed": OBS["pooled_over_within"], "pass": t3[0] <= OBS["pooled_over_within"] <= t3[1]},
             "T4_state_cv": {"range": t4, "observed": obs_cv, "pass": t4[0] <= obs_cv <= t4[1]}}
    i, mult, theta_c, rho, tau = CELLS[name]
    return {"law": {"multiplier": mult, "theta_c": theta_c, "rho": rho, "tau": tau, "offspring_k": K, "transmission_share": SHARE},
            "calibration": {k: info[k] for k in ("lambda", "amp", "sigma", "shock_s", "calibration_failed", "trend_relvar")},
            "panels": len(rows), "tests": tests, "all_pass": all(t["pass"] for t in tests.values()),
            "failing_tests": [k for k, t in tests.items() if not t["pass"]],
            "b_over_bP_mean": float(np.nanmean(b / bP)), "pooled_over_within_mean": float(np.nanmean(b / bw)),
            "state_cv_mean": float(np.nanmean(f("state_cv"))),
            "p10_mean": float(np.nanmean(f("p10"))), "median_mean": float(np.nanmean(f("median"))), "p90_mean": p90,
            "mc_sd_a": sd_a, "sd_ratio_observed_se_over_mc_sd": SE_OBS / sd_a,
            "coverage_of_E_a": cover([r.get("a_ci") for r in rows], np.nanmean(a)),
            "coverage_of_true_total": cover([r.get("a_ci") for r in rows], vR + vH),
            "a_over_vR": R.ratio_of_means(a, vR), "a_over_total": R.ratio_of_means(a, vR + vH), "mean_a": float(np.nanmean(a))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    L2.register_script(__file__)
    t0 = time.time()
    log = Log(R.LOGS / "rev_r2_generator.log")
    f = load_sample("RA", 3.6)
    loss, _ = relative_loss(f.x.to_numpy(float))
    obs_cv = state_cv(loss, f.D.to_numpy(float), f.state.to_numpy())
    with Pool(args.workers) as pool:
        infos = {r["generator"]: r for r in pool.imap_unordered(run_cell, [(n, args.reps) for n in CELLS])}
    cells = {n: summarise(n, R.read_jsonl(R.AN / f"rev_r2_generator_{n}.jsonl"), infos[n], obs_cv) for n in CELLS}
    ok = [n for n, c in cells.items() if c["all_pass"]]
    res = {"meta": {"script": "code/rev_r2_generator.py", "plan": "analysis/analysis_plan_amendment2.md, Section D.1",
                    "reps": args.reps, "seed": [R.SEED, 42, "i"], "calibration_stream": "920 + i",
                    "labelling": "data-informed generator design; no inferential result uses these generators",
                    "seconds": round(time.time() - t0, 1)},
           "observed": {**OBS, "state_cv": obs_cv}, "cells": cells,
           "success_rule": {"text": "T1-T4 all pass in at least one cell", "cells_passing": ok,
                            "verdict": "reproduces the key features" if ok else "M4/M5 unresolved"}}
    R.write_json(OUT, res)
    log("done", res["success_rule"]["verdict"], {n: c["failing_tests"] for n, c in cells.items()})


if __name__ == "__main__":
    main()
