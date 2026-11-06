"""Shared helpers for the 2026-11 revision analyses (rev_*.py).

The revision analyses follow the locked plan in revision/analysis_plan.md. They reuse the
shipped design (ma_audit.setup) and the shipped calibrated generator parameters
(analysis/ma_audit.json, analysis/ma_audit_national.json) so that every new generator is
a documented variant of a shipped one.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import ma_audit as A
from revision_data import wl_R

EXP = Path(__file__).resolve().parents[1]
AN, LOGS, TAB, FIG = EXP / "analysis", EXP / "logs", EXP / "tables", EXP / "figures"
SEED = 20261105
CHI2_3_95 = 7.814727903251178


def null_log(*_):
    return None


def setup_design():
    """Real R_A design and shipped generators (calibrations are restored from JSON, not refitted)."""
    e, d, obs, z, gens = A.setup(null_log)
    return e, d, obs, z, gens


def restore(gen, name):
    """Set a shipped generator to its calibrated (amp, lambda) as stored in ma_audit*.json."""
    if name == "S9_national":
        r = json.loads((AN / "ma_audit_national.json").read_text())["audit"]
        amp, lam = r["amp"], r["lambda"]
    else:
        al = json.loads((AN / "ma_audit.json").read_text())["alignment"][name]
        amp, lam = (al["amp"]["amp"] if al["amp"] else 1.0), al["lambda"]["lambda"]
    gen.set(amp=amp, lam=lam)
    return gen


def append_jsonl(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(obj, default=_default) + "\n")


def read_jsonl(path):
    path = Path(path)
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def _default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    raise TypeError(type(o))


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, default=_default), encoding="utf-8")


def mc(values):
    """Monte Carlo summary: mean, sd, 2.5/97.5% quantiles, MCSE of the mean."""
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if not len(v):
        return {"n": 0}
    return {"n": int(len(v)), "mean": float(v.mean()), "sd": float(v.std(ddof=1)) if len(v) > 1 else 0.0,
            "q025": float(np.quantile(v, .025)), "q975": float(np.quantile(v, .975)),
            "mcse_mean": float(v.std(ddof=1) / np.sqrt(len(v))) if len(v) > 1 else None}


def ratio_of_means(num, den):
    num, den = np.asarray(num, float), np.asarray(den, float)
    ok = np.isfinite(num) & np.isfinite(den)
    return float(num[ok].mean() / den[ok].mean()) if ok.any() and den[ok].mean() != 0 else None


def rel_var(r):
    r = np.asarray(r, float)
    return float(np.var(r) / np.mean(r) ** 2)


def implied_sd(a_obs, a_ci, ratio_point, ratio_q025=None, ratio_q975=None):
    """sqrt(a / r): point; lower uses the lower CI of a and the upper ratio quantile; upper the reverse."""
    out = {"point": float(np.sqrt(a_obs / ratio_point))}
    lo_r = ratio_q975 if ratio_q975 is not None else ratio_point
    hi_r = ratio_q025 if ratio_q025 is not None else ratio_point
    out["propagated"] = [float(np.sqrt(a_ci[0] / lo_r)), float(np.sqrt(a_ci[1] / hi_r))]
    return out


def webb_weights(rng, size):
    w = np.array([-np.sqrt(1.5), -1.0, -np.sqrt(.5), np.sqrt(.5), 1.0, np.sqrt(1.5)])
    return w[rng.integers(0, 6, size)]


def wl(g):
    return wl_R(g)


def cgm_from_loss(x, xbar, loss, D, X, beta, state, week):
    """Two-way CGM covariance of beta for the stacked (xbar, beta) equations when the loss is given
    (bootstrap losses with xbar held at the observed mean, as in ma_calibration.sandwich).
    Returns (V2 of beta, mu) or (None, None) when a fitted mean is not positive."""
    n, k = len(x), X.shape[1] + 1
    mu = X @ beta
    if np.any(mu <= 0):
        return None, None
    r = (loss - mu) / mu ** 2
    psi = np.column_stack([x - xbar, X * r[:, None]])
    J = np.zeros((k, k))
    J[0, 0] = -n
    u = (x - xbar) / xbar
    J[1:, 0] = X.T @ ((2 * u * (-x / xbar ** 2)) / mu ** 2)
    J[1:, 1:] = -(X.T * (1 / mu ** 2)) @ X
    Ji = np.linalg.inv(J)

    def meat(lab):
        _, ci = np.unique(lab, return_inverse=True)
        S = np.zeros((ci.max() + 1, k))
        np.add.at(S, ci, psi)
        return S.T @ S, ci.max() + 1
    MS, GS = meat(state)
    MW, GW = meat(week)
    V = Ji @ (GS / (GS - 1) * MS + GW / (GW - 1) * MW - psi.T @ psi) @ Ji.T
    w, U = np.linalg.eigh((V + V.T) / 2)
    return ((U * np.maximum(w, 0)) @ U.T)[1:, 1:], mu


def week_multiplicity(week, rng, strata=None):
    """Calendar-week cluster bootstrap: each row gets the number of times its week is drawn.
    strata (per-row labels constant within week): resample weeks within each stratum."""
    weeks, wi = np.unique(week, return_inverse=True)
    m = np.zeros(len(weeks))
    if strata is None:
        np.add.at(m, rng.integers(0, len(weeks), len(weeks)), 1)
    else:
        s_of_week = np.empty(len(weeks), dtype=object)
        s_of_week[wi] = strata
        for s in np.unique(s_of_week):
            idx = np.flatnonzero(s_of_week == s)
            np.add.at(m, idx[rng.integers(0, len(idx), len(idx))], 1)
    return m[wi]


def pct_ci(draws, level=.95):
    v = np.asarray(draws, float)
    v = v[np.isfinite(v)]
    return [float(np.quantile(v, (1 - level) / 2)), float(np.quantile(v, (1 + level) / 2))] if len(v) else [None, None]
