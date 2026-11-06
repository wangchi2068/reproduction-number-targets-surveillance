"""Shared helpers for the second-round revision analyses (rev_r2_*.py; analysis_plan_amendment2.md).

Generic linear gamma quasi-likelihood fits (E(L) = X beta, weights 1/mu^2, positive fitted means),
nonlinear power-law fits by quasi-deviance minimisation, state population covariates from the Census
Bureau vintage-2025 state totals, and the script-hash register of Amendment 2, Section H.
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

import rev_common as R

POP_FILE = R.EXP / "data" / "external" / "census_NST-EST2025-ALLDATA.csv"
POP_URL = "https://www2.census.gov/programs-surveys/popest/datasets/2020-2025/state/totals/NST-EST2025-ALLDATA.csv"
HASHES = R.EXP.parent / "revision2" / "script_hashes.jsonl"
if not HASHES.parent.is_dir():  # published or clean copies: keep the record next to the logs
    HASHES = R.LOGS / "rev_r2_script_hashes.jsonl"
ABBR = {"Alabama": "AL", "Alaska": "AK", "Arizona": "AZ", "Arkansas": "AR", "California": "CA", "Colorado": "CO",
        "Connecticut": "CT", "Delaware": "DE", "District of Columbia": "DC", "Florida": "FL", "Georgia": "GA",
        "Hawaii": "HI", "Idaho": "ID", "Illinois": "IL", "Indiana": "IN", "Iowa": "IA", "Kansas": "KS",
        "Kentucky": "KY", "Louisiana": "LA", "Maine": "ME", "Maryland": "MD", "Massachusetts": "MA",
        "Michigan": "MI", "Minnesota": "MN", "Mississippi": "MS", "Missouri": "MO", "Montana": "MT",
        "Nebraska": "NE", "Nevada": "NV", "New Hampshire": "NH", "New Jersey": "NJ", "New Mexico": "NM",
        "New York": "NY", "North Carolina": "NC", "North Dakota": "ND", "Ohio": "OH", "Oklahoma": "OK",
        "Oregon": "OR", "Pennsylvania": "PA", "Rhode Island": "RI", "South Carolina": "SC", "South Dakota": "SD",
        "Tennessee": "TN", "Texas": "TX", "Utah": "UT", "Vermont": "VT", "Virginia": "VA", "Washington": "WA",
        "West Virginia": "WV", "Wisconsin": "WI", "Wyoming": "WY"}


def register_script(path, note="full run"):
    """Amendment 2, Section H: sha256 and UTC time of a script before its full run."""
    p = Path(path)
    R.append_jsonl(HASHES, {"script": p.name, "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
                            "utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), "argv": sys.argv[1:], "note": note})


def population():
    """State population: mean of POPESTIMATE2022-2025 (vintage 2025), its log and z-score over 51 jurisdictions."""
    raw = pd.read_csv(POP_FILE, encoding="latin1")
    raw = raw[raw.SUMLEV.astype(int) == 40]
    raw = raw[raw.NAME.isin(ABBR)]
    pop = raw[[f"POPESTIMATE{y}" for y in (2022, 2023, 2024, 2025)]].mean(1)
    out = pd.DataFrame({"state": raw.NAME.map(ABBR).to_numpy(), "pop": pop.to_numpy()})
    out["log_pop"] = np.log(out["pop"])
    out["z_pop"] = (out.log_pop - out.log_pop.mean()) / out.log_pop.std(ddof=0)
    return out.sort_values("state").reset_index(drop=True)


def fit_linear_ql(y, X, w=None, start=None, maxit=500, tol=1e-10):
    """Gamma QL for E(y) = X beta: IRLS with weights w / mu^2 and step halving to keep mu > 0.
    y may contain negative values (unclamped wild-bootstrap losses); only mu > 0 is required."""
    y, X = np.asarray(y, float), np.asarray(X, float)
    w = np.ones(len(y)) if w is None else np.asarray(w, float)
    keep = w > 0
    y, X, w = y[keep], X[keep], w[keep]
    scale = float(np.average(np.abs(y), weights=w))
    if not np.isfinite(scale) or scale <= 0:
        return {"converged": False, "beta": None}
    ys = y / scale
    if start is None:
        beta = np.zeros(X.shape[1])
        beta[0] = float(np.average(ys, weights=w)) if np.average(ys, weights=w) > 0 else 1.0
    else:
        beta = np.asarray(start, float) / scale
    mu = X @ beta
    if np.any(mu <= 0):
        beta = np.zeros(X.shape[1])
        beta[0] = 1.0
        mu = X @ beta
    for it in range(maxit):
        W = w / mu ** 2
        try:
            prop = np.linalg.solve((X.T * W) @ X, (X.T * W) @ ys)
        except np.linalg.LinAlgError:
            return {"converged": False, "beta": None}
        step = 1.0
        for _ in range(60):
            cand = beta + step * (prop - beta)
            mc = X @ cand
            if np.all(mc > 0):
                break
            step /= 2
        else:
            return {"converged": False, "beta": None}
        change = float(np.max(np.abs(mc - mu) / np.maximum(np.abs(mu), 1e-12)))
        beta, mu = cand, mc
        if change < tol and step == 1.0:
            return {"converged": True, "beta": beta * scale, "iterations": it + 1}
    return {"converged": False, "beta": beta * scale}


def qdev(y, mu, w=None):
    w = np.ones(len(y)) if w is None else w
    if np.any(mu <= 0) or not np.all(np.isfinite(mu)):
        return np.inf
    return float(np.sum(w * (y / mu + np.log(mu))))


def fit_power(y, D, w=None):
    """E(L) = c D^-gamma by gamma QL (quasi-deviance minimisation over log c, gamma)."""
    y, D = np.asarray(y, float), np.asarray(D, float)
    w = np.ones(len(y)) if w is None else np.asarray(w, float)
    lD = np.log(D)
    m = float(np.average(y, weights=w))

    def obj(p):
        return qdev(y, np.exp(p[0] - p[1] * lD), w)
    best = None
    for g0 in (0.2, 0.4, 0.7, 1.0):
        st = [np.log(m) + g0 * float(np.average(lD, weights=w)), g0]
        o = minimize(obj, st, method="Nelder-Mead", options={"xatol": 1e-9, "fatol": 1e-10, "maxiter": 4000})
        if best is None or o.fun < best.fun:
            best = o
    return {"converged": bool(np.isfinite(best.fun)), "c": float(np.exp(best.x[0])), "gamma": float(best.x[1]),
            "qdev": float(best.fun)}


def fit_three(y, D, w=None, affine=None):
    """E(L) = a + b D^-gamma, gamma in [0.05, 3], signed a with mu > 0 (Amendment 2, B.3).
    Profile over gamma: for fixed gamma the model is linear in (a, b) and fitted by fit_linear_ql."""
    y, D = np.asarray(y, float), np.asarray(D, float)
    w = np.ones(len(y)) if w is None else np.asarray(w, float)

    def prof(g):
        X = np.column_stack([np.ones(len(D)), D ** (-g)])
        f = fit_linear_ql(y, X, w)
        if not f["converged"]:
            return np.inf, None
        return qdev(y, X @ f["beta"], w), f["beta"]
    grid = np.concatenate([np.linspace(.05, 1.5, 30), np.linspace(1.6, 3.0, 8)])
    vals = [prof(g)[0] for g in grid]
    j = int(np.nanargmin(vals))
    lo, hi = grid[max(j - 1, 0)], grid[min(j + 1, len(grid) - 1)]
    # golden-section refinement on the bracket
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
    q, beta = prof(g)
    if beta is None:
        return {"converged": False}
    return {"converged": True, "a": float(beta[0]), "b": float(beta[1]), "gamma": float(g), "qdev": float(q),
            "at_bound": bool(g <= .051 or g >= 2.999)}


def cgm_general(x, xbar, loss, X, beta, state, week):
    """rev_common.cgm_from_loss for an arbitrary linear design X (no D argument needed)."""
    return R.cgm_from_loss(x, xbar, loss, None, X, beta, state, week)


def results_block(key, block):
    """Write one block under revision_2026_11_r2 in analysis/results.json (read-modify-write)."""
    path = R.AN / "results.json"
    res = json.loads(path.read_text(encoding="utf-8"))
    res.setdefault("revision_2026_11_r2", {})[key] = block
    R.write_json(path, res)
