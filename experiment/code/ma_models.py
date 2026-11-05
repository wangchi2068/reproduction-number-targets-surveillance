"""Shared helpers for the mainline-A analyses; import performs no IO.

Adds to the authors' estimators (revision_stats, s03c_state_fe) without
changing them:
  * weighted two-way (state + week) demeaning by alternating projections and a
    convergent precision-weighted within slope for the affine loss model;
  * Monte Carlo / bootstrap summaries and a progress logger.
"""
from __future__ import annotations

import json
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from revision_stats import fit_affine, json_safe


def utc():
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S")


class Log:
    """Line-buffered progress log written incrementally to logs/<name>.log."""

    def __init__(self, path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.t0 = time.time()

    def __call__(self, *parts):
        line = f"{utc()} +{time.time()-self.t0:7.1f}s " + " ".join(str(p) for p in parts)
        with self.path.open("a", encoding="utf-8") as fh:
            fh.write(line + "\n")
        print(line, flush=True)


def summ(values, point=None):
    """Distribution summary: mean, sd, 2.5/97.5 percentiles over finite values."""
    v = np.asarray(values, float)
    good = v[np.isfinite(v)]
    out = {"n_valid": int(len(good)), "n_total": int(len(v))}
    if len(good):
        out.update(mean=float(good.mean()), sd=float(good.std(ddof=1)) if len(good) > 1 else None,
                   q025=float(np.quantile(good, .025)), q975=float(np.quantile(good, .975)))
    if point is not None:
        out["point"] = None if point is None or not np.isfinite(point) else float(point)
    return out


MC_BOOT = 500  # bootstrap resamples of the replicate values for quantile MC SEs


def summ_mc(values, point=None):
    """summ plus Monte Carlo standard errors (Morris et al. 2019): mean SE = sd/sqrt(n);
    2.5/97.5 percentile SEs from MC_BOOT resamples of the replicates (fixed stream)."""
    out = summ(values, point)
    v = np.asarray(values, float)
    good = v[np.isfinite(v)]
    n = len(good)
    if n > 1:
        out["mcse_mean"] = float(good.std(ddof=1) / np.sqrt(n))
        idx = np.random.default_rng([20261104, n]).integers(n, size=(MC_BOOT, n))
        qs = np.quantile(good[idx], (.025, .975), axis=1)
        out["mcse_q025"], out["mcse_q975"] = float(qs[0].std(ddof=1)), float(qs[1].std(ddof=1))
    return out


def binom_mcse(p, n):
    """Monte Carlo SE of an estimated proportion or quantile level from n replicates."""
    return float(np.sqrt(p * (1 - p) / n)) if n else None


def boot_ci(point, draws, ratio=False):
    """Percentile bootstrap interval in the authors' convention (finite_summary)."""
    a = np.asarray(draws, float)
    good = np.isfinite(a)
    out = {"value": None if point is None or not np.isfinite(point) else float(point),
           "ci95": np.percentile(a[good], [2.5, 97.5]).tolist() if good.any() else None,
           "valid_draws": int(good.sum()), "total_draws": int(len(a))}
    if ratio:
        nonpos = float(np.mean(a[good] <= 0)) if good.any() else None
        out["nonpositive_fraction"] = nonpos
        if not good.all() or (nonpos is not None and nonpos > .025):
            out["interval_status"] = "unstable denominator; percentile interval shown but flagged"
    return out


def write(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_safe(obj), indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")


def codes_of(values):
    _, codes = np.unique(np.asarray(values), return_inverse=True)
    return codes


def demean_two_way(v, s, t, w, tol=1e-13, maxit=10000):
    """Weighted projection of v off state and week indicators (alternating projections).

    Returns (residual, sweeps, converged). For a connected state-week design the
    method of alternating projections converges to the exact two-way within
    transformation (the weighted LS residual on both sets of dummies).
    """
    ns, nt = s.max() + 1, t.max() + 1
    ws = np.bincount(s, weights=w, minlength=ns)
    wt = np.bincount(t, weights=w, minlength=nt)
    r = v.astype(float).copy()
    scale = max(float(np.sqrt(w @ (v * v) / w.sum())), 1e-300)
    for sweep in range(1, maxit + 1):
        ms = np.divide(np.bincount(s, weights=w * r, minlength=ns), ws, out=np.zeros(ns), where=ws > 0)
        r -= ms[s]
        mt = np.divide(np.bincount(t, weights=w * r, minlength=nt), wt, out=np.zeros(nt), where=wt > 0)
        r -= mt[t]
        if max(np.max(np.abs(ms[ws > 0])), np.max(np.abs(mt[wt > 0]))) < tol * scale:
            return r, sweep, True
    return r, maxit, False


def demean_one_way(v, g, w):
    ng = g.max() + 1
    tot = np.bincount(g, weights=w, minlength=ng)
    m = np.divide(np.bincount(g, weights=w * v, minlength=ng), tot, out=np.zeros(ng), where=tot > 0)
    return v - m[g]


def within_slope(loss, D, groups, w=None, precision="pooled", iterate=0, floor_frac=.1, damp=1.0):
    """Precision-weighted within slope of loss on 1/D with group fixed effects.

    groups: ("state",), ("week",) or ("state", "week") arrays passed as a tuple
    of integer code vectors. precision="pooled" uses the gamma working
    precision 1/mu^2 from the authors' pooled affine fit (one-step FGLS, a
    linear WLS estimator that always exists); iterate>0 re-weights with the
    fitted fixed-effect means floored at floor_frac*mean(loss) and reports
    whether the slope converged (|change| < 1e-8 relative). damp<1 averages the
    new working mean with the previous one (mu <- damp*new + (1-damp)*old),
    which removes the period-2 oscillation of the undamped two-way iteration.
    """
    y = np.asarray(loss, float)
    z = 1.0 / np.asarray(D, float)
    w0 = np.ones(len(y)) if w is None else np.asarray(w, float)
    keep = w0 > 0
    y, z, w0 = y[keep], z[keep], w0[keep]
    gs = [np.unique(g[keep], return_inverse=True)[1] for g in groups]
    base = fit_affine(y, 1 / z, w0)
    if not base["converged"]:
        return {"b": None, "converged": False, "failure": "pooled_fit_failed"}
    mu = base["a"] + base["b"] * z
    if np.min(mu) <= 0:
        mu = np.maximum(mu, floor_frac * float(np.average(y, weights=w0)))
    history = []
    converged, sweeps_ok = True, True
    b = None
    for it in range(iterate + 1):
        p = w0 / mu ** 2
        if len(gs) == 1:
            ry, rz = demean_one_way(y, gs[0], p), demean_one_way(z, gs[0], p)
        else:
            ry, k1, c1 = demean_two_way(y, gs[0], gs[1], p)
            rz, k2, c2 = demean_two_way(z, gs[0], gs[1], p)
            sweeps_ok = sweeps_ok and c1 and c2
        den = float(p @ (rz * rz))
        if den <= 0:
            return {"b": None, "converged": False, "failure": "no_within_variation"}
        b_new = float(p @ (rz * ry) / den)
        history.append(b_new)
        if b is not None and abs(b_new - b) <= 1e-8 * max(1.0, abs(b)):
            b = b_new
            break
        b = b_new
        if it < iterate:
            fitted = y - (ry - b * rz)          # fixed effects + slope part of the fit
            mu = damp * np.maximum(fitted, floor_frac * float(np.average(y, weights=w0))) + (1 - damp) * mu
    else:
        converged = iterate == 0
    return {"b": b, "converged": bool(converged and sweeps_ok), "projection_converged": bool(sweeps_ok),
            "iterations": len(history), "history": history[-5:], "damp": damp, "pooled_a": base["a"], "pooled_b": base["b"]}
