"""Thinning as a test of the delta approximation (mainline A, P1).

For every R_A row and every thinning probability q the thinned counts are
X ~ Bin(C_t, q), Y ~ Bin(C_next, q) and the thinned statistic is
x' = wl_R((Y + .5)/(X + .5)) = max(1 + s log((Y+.5)/(X+.5))/7, 0)^shape.
E[x'] and E[x'^2] are computed exactly by summation over the binomial supports
(mass below 1e-13 dropped and reported), including the positive-part
truncation. The expected thinned loss
    E L'_i ~= (Var x'_i + (E x'_i - xbar'_E)^2) / xbar'_E^2
is fitted with the authors' affine FGLS on the unthinned C_t, giving an exact
(second-order complete) expected slope increment per q. The delta method
predicts a constant increment rate b_P = mean kappa^2 (1 + 1/g) per unit
(1/q - 1). Outputs: rate per q (exact vs delta vs the authors' Monte Carlo),
the least-squares rate through the origin with and without q = 0.1, the exact
Var log(Bin(C,q)+.5) relative to the delta value (1-q)/(qC), the positive-part
mass, and an independent-thinning Monte Carlo check of the expected-loss fit.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.stats import binom

from ma_models import Log, write
from revision_data import load_sample, wl_R
from revision_stats import fit_affine, fit_loss, fit_loss_batch
from s03b_thinning import QS, poisson_reference

EXP = Path(__file__).resolve().parents[1]
SEED = 2026110353
MEAN, SD = 3.6, 1.6
SHAPE, SCALE = (MEAN / SD) ** 2, SD ** 2 / MEAN
TAIL = 1e-13


def support(n, q):
    k = np.arange(int(n) + 1)
    p = binom.pmf(k, int(n), q)
    keep = p > TAIL
    return k[keep], p[keep], float(1 - p[keep].sum())


def row_moments(ct, cn, q):
    xs, px, tx = support(ct, q)
    ys, py, ty = support(cn, q)
    lg = np.log(ys + .5)[:, None] - np.log(xs + .5)[None, :]
    base = 1 + SCALE * lg / 7
    v = np.power(np.maximum(base, 0), SHAPE)
    P = py[:, None] * px[None, :]
    tot = P.sum()
    m1 = float((P * v).sum() / tot)
    m2 = float((P * v * v).sum() / tot)
    bind = float(P[base <= 0].sum() / tot)
    return m1, max(m2 - m1 * m1, 0.0), bind, tx + ty


def exact_var_log(C, q):
    k, p, _ = support(C, q)
    l = np.log(k + .5)
    m = p @ l / p.sum()
    return float(p @ (l - m) ** 2 / p.sum())


def main():
    log = Log(EXP / "logs/ma_thinning_exact.log")
    e = load_sample("RA")
    ct, cn, D = e.C_t.to_numpy(int), e.C_next.to_numpy(int), e.C_t.to_numpy(float)
    x = e.x.to_numpy(float)
    base = fit_loss(x, D)
    bp = poisson_reference(e)
    inv = 1 / QS - 1
    per_q = []
    for q in QS:
        mom = np.array([row_moments(a, b, q) for a, b in zip(ct, cn)])
        m1, var, bind, tail = mom.T
        xb = m1.mean()
        EL = (var + (m1 - xb) ** 2) / xb ** 2
        f = fit_affine(EL, D)
        # same expected-loss construction at q = 1 is the observed loss itself
        rate = (f["b"] - base["b"]) / (1 / q - 1)
        per_q.append({"q": float(q), "expected_b": f["b"], "expected_a": f["a"], "delta_b": f["b"] - base["b"],
                      "rate_exact": rate, "rate_exact_over_bP": rate / bp,
                      "mean_E_x": float(xb), "bias_mean_E_x_minus_xbar": float(xb - x.mean()),
                      "positive_part_mass_max": float(bind.max()), "rows_positive_part_mass_gt_1e-12": int((bind > 1e-12).sum()),
                      "dropped_tail_mass_max": float(tail.max())})
        log("q", q, "exact rate", rate, "rate/bP", rate / bp, "bind max", bind.max())
    rates = np.array([r["rate_exact"] for r in per_q])
    db = rates * inv
    c_exact = float(inv @ db / (inv @ inv))
    c_noq01 = float(inv[:-1] @ db[:-1] / (inv[:-1] @ inv[:-1]))
    # authors' Monte Carlo rates (shared-count thinning, analysis/authors/thinning.json)
    src = EXP / "analysis/authors/thinning.json"
    if not src.exists():  # run_all.sh runs the authors' pipeline first; never read outputs from outside the package
        raise SystemExit(f"missing {src}: run the authors' pipeline (run_all.sh lane 1) first")
    auth = json.loads(src.read_text())["estimates"]
    mc = [{"q": float(q), "rate_mc": auth[f"full/q{q:.1f}/c"]["value"], "ci95": auth[f"full/q{q:.1f}/c"]["ci95"]} for q in QS]
    # leverage of each q in the least-squares rate through the origin
    lev = (inv ** 2 / (inv @ inv)).tolist()
    # small-count variance of log(Bin(C,q)+.5) against the delta value (1-q)/(qC)
    vl = []
    for C in (5, 10, 20, 45, 100, 500):
        for q in QS:
            ex = exact_var_log(C, q) - 0.0
            vl.append({"C": C, "q": float(q), "exact_var_log": ex, "delta_var_log_added": (1 - q) / (q * C),
                       "exact_added_over_delta": (ex - exact_var_log(C, 1.0)) / ((1 - q) / (q * C))})
    # independent-thinning Monte Carlo check of the expected-loss approximation
    rng = np.random.default_rng(SEED)
    check = []
    for q in (0.5, 0.1):
        X = rng.binomial(np.broadcast_to(ct, (200, len(ct))), q)
        Y = rng.binomial(np.broadcast_to(cn, (200, len(cn))), q)
        fb = fit_loss_batch(wl_R((Y + .5) / (X + .5)), D)
        mc_rate = (float(np.nanmean(fb["b"])) - base["b"]) / (1 / q - 1)
        ex = next(r for r in per_q if abs(r["q"] - q) < 1e-9)["rate_exact"]
        check.append({"q": q, "independent_mc_rate": mc_rate, "mc_se": float(np.nanstd(fb["b"], ddof=1) / np.sqrt(200) / (1 / q - 1)),
                      "expected_loss_rate": ex})
        log("mc check q", q, mc_rate, "expected-loss", ex)
    out = {"authors_mc_source": str(src.resolve().relative_to(EXP)), "bP_delta": bp, "observed_b": base["b"], "per_q": per_q, "authors_mc": mc, "lsq_leverage_by_q": dict(zip([f"{q:.1f}" for q in QS], lev)),
           "c_exact_all_q": c_exact, "c_exact_all_q_over_bP": c_exact / bp,
           "c_exact_without_q01": c_noq01, "c_exact_without_q01_over_bP": c_noq01 / bp,
           "var_log_small_count": vl, "mc_check": check,
           "approximation": "expected thinned loss (Var x' + (E x' - mean E x')^2)/mean E x'^2, independent X,Y within a row, cross-row covariance ignored"}
    write(EXP / "analysis/ma_thinning_exact.json", out)
    log("done c_exact", c_exact, "over bP", c_exact / bp, "noq01", c_noq01 / bp)


if __name__ == "__main__":
    main()
