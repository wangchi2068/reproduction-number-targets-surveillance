"""Confirmatory intercept test, sensitivities, within-week decomposition and error components
(analysis plan Sections 2 and 4; review M1-M3, M6).

E1  primary: two-way (state x calendar week) CGM sandwich Wald t of the stacked (xbar, a, b) estimator,
    t with G - 1 = 41 df, two-sided, alpha = 0.05 (ma_calibration.sandwich / ab_intervals).
    secondary: state-only and week-only sandwich; week-cluster wild bootstrap-t with a = 0 imposed
    (Webb weights, B = 9,999, seed [20261105, 21]); calendar-week block bootstrap percentile CIs
    (B = 1,999, seed [20261105, 22]).
S-a1..S-a9 sensitivities (all reported) and the pre-specified robustness rule.
E3  within-week intercept under four centring definitions, block bootstrap (B = 1,999, seed [20261105, 23]).
Sec. 4  b, b_P, b/b_P, b/a and predicted error components a, b/C at C = 10, 37, 100, 1,000.

Every finished block is appended to analysis/rev_intercept.jsonl; the summary goes to analysis/rev_intercept.json.
Usage: python code/rev_intercept.py [--B-wild 9999] [--B-block 1999]
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from scipy import stats

import ma_calibration as MC
import rev_common as R
from ma_models import Log
from revision_data import CALENDAR, build_panel, load_raw, load_sample, rule_windows
from revision_empirical import centered_loss, group_center
from revision_stats import fit_affine, fit_loss, relative_loss
from s03b_thinning import poisson_reference

OUT_JSONL = R.AN / "rev_intercept.jsonl"
OUT_JSON = R.AN / "rev_intercept.json"
PAUSE = (pd.Timestamp("2025-10-03"), pd.Timestamp("2025-11-13"))
C_GRID = (10, 37, 100, 1000)


def fit(frame, gi_mean=3.6, gi_sd=None):
    """Point estimates and sandwich inference for one sample (x already computed for its GI)."""
    x, D = frame.x.to_numpy(float), frame.D.to_numpy(float)
    st, wk = frame.state.to_numpy(), frame.week_end.to_numpy()
    s = MC.sandwich(x, D, st, wk)
    out = {"n_rows": int(len(frame)), "n_states": int(len(np.unique(st))), "n_weeks": int(len(np.unique(wk)))}
    if s is None:
        return {**out, "converged": False}
    G = min(out["n_states"], out["n_weeks"])
    df = G - 1
    iv = MC.ab_intervals(s, stats.t.ppf(.975, df))
    a, b = (float(v) for v in s["beta"])
    bP = poisson_reference(frame, mean=gi_mean, sd=gi_sd)
    res = {}
    for lab in ("two_way", "state", "week"):
        t = a / iv[lab]["se_a"]
        res[lab] = {"se_a": float(iv[lab]["se_a"]), "a_ci": [float(v) for v in iv[lab]["a_ci"]],
                    "t": float(t), "p_two_sided": float(2 * stats.t.sf(abs(t), df)),
                    "se_b": float(iv[lab]["se_b"]), "b_ci": [float(v) for v in iv[lab]["b_ci"]],
                    "b_over_a_ci": [float(v) for v in iv[lab]["b_over_a_ci"]]}
    return {**out, "converged": True, "a": a, "b": b, "b_over_a": b / a if a > 0 else None, "b_P": bP,
            "b_over_bP": b / bP, "b_over_bP_ci_two_way": [v / bP for v in res["two_way"]["b_ci"]],
            "G": G, "t_df": df, "neg_eig": int(s["neg_eig"]), "inference": res, "V2": s["V2"].tolist()}


def batch_mcse_quantile(draws, q, nb=10):
    v = np.asarray(draws, float)
    v = v[np.isfinite(v)]
    parts = np.array_split(v, nb)
    return float(np.std([np.quantile(p, q) for p in parts], ddof=1) / np.sqrt(nb))


def ci_with_mcse(draws):
    lo, hi = R.pct_ci(draws)
    return {"ci": [lo, hi], "mcse": [batch_mcse_quantile(draws, .025), batch_mcse_quantile(draws, .975)],
            "n_finite": int(np.isfinite(np.asarray(draws, float)).sum())}


def wild_bootstrap_t(frame, t_obs, B, log):
    """Week-cluster wild bootstrap-t, a = 0 imposed (restricted gamma QL fit mu0 = b0 / D, b0 = mean(D L))."""
    x, D = frame.x.to_numpy(float), frame.D.to_numpy(float)
    st, wk = frame.state.to_numpy(), frame.week_end.to_numpy()
    loss, xbar = relative_loss(x)
    b0 = float(np.mean(D * loss))
    mu0 = b0 / D
    _, wi = np.unique(wk, return_inverse=True)
    X = np.column_stack([np.ones(len(x)), 1 / D])
    rng = np.random.default_rng([R.SEED, 21])
    ts, clamped = np.full(B, np.nan), 0
    for r in range(B):
        v = R.webb_weights(rng, wi.max() + 1)[wi]
        Lb = mu0 + v * (loss - mu0)
        clamped += int(np.sum(Lb < 1e-7))
        Lb = np.maximum(Lb, 1e-7)
        f = fit_affine(Lb, D)
        if not f["converged"]:
            continue
        beta = np.array([f["a"], f["b"]])
        V, _ = R.cgm_from_loss(x, xbar, Lb, D, X, beta, st, wk)
        if V is None or V[0, 0] <= 0:
            continue
        ts[r] = beta[0] / np.sqrt(V[0, 0])
        if r % 2000 == 1999:
            log("wild bootstrap-t", r + 1)
    ok = np.isfinite(ts)
    p = float((1 + np.sum(np.abs(ts[ok]) >= abs(t_obs))) / (ok.sum() + 1))
    return {"B": B, "valid": int(ok.sum()), "b0_restricted": b0, "t_obs": t_obs, "p_two_sided": p,
            "mcse": float(np.sqrt(p * (1 - p) / ok.sum())), "t_star_q": [float(v) for v in np.quantile(ts[ok], [.025, .5, .975])],
            "clamped_share": clamped / (B * len(x)), "seed": [R.SEED, 21],
            "weights": "Webb 6-point per calendar week; losses clamped at 1e-7; xbar held at the observed mean"}


def block_bootstrap(frame, B, log):
    """Calendar-week cluster bootstrap (row weights = week multiplicities)."""
    x, D = frame.x.to_numpy(float), frame.D.to_numpy(float)
    wk = frame.week_end.to_numpy()
    rng = np.random.default_rng([R.SEED, 22])
    draws = {k: np.full(B, np.nan) for k in ("a", "b", "b_P", "b_over_bP", "b_over_a")}
    comp = {C: np.full(B, np.nan) for C in C_GRID}
    for r in range(B):
        w = R.week_multiplicity(wk, rng)
        f = fit_loss(x, D, w)
        if not f["converged"]:
            continue
        bP = poisson_reference(frame, w)
        draws["a"][r], draws["b"][r], draws["b_P"][r] = f["a"], f["b"], bP
        draws["b_over_bP"][r] = f["b"] / bP
        draws["b_over_a"][r] = f["b"] / f["a"] if f["a"] > 0 else np.nan
        for C in C_GRID:
            comp[C][r] = f["a"] + f["b"] / C
        if r % 500 == 499:
            log("block bootstrap", r + 1)
    out = {k: ci_with_mcse(v) for k, v in draws.items()}
    out["share_a_le_0"] = float(np.mean(draws["a"][np.isfinite(draws["a"])] <= 0))
    out["error_total"] = {str(C): ci_with_mcse(v) for C, v in comp.items()}
    out["meta"] = {"B": B, "seed": [R.SEED, 22], "resampling": "calendar weeks with replacement; xbar and b_P recomputed per draw"}
    return out, draws


def review_week_loss(frame, w):
    """Review definition: week mean as both centre and scale."""
    x = frame.x.to_numpy(float)
    c = group_center(x, frame.week_end, w)
    with np.errstate(divide="ignore", invalid="ignore"):
        L = ((x - c) / c) ** 2
    return L, np.isfinite(L) & (w > 0)


E3_MODES = ("week_leave_state_out", "week", "week_review", "season", "pooled")


def e3_fit(frame, mode, w):
    D = frame.D.to_numpy(float)
    if mode == "week_review":
        L, valid = review_week_loss(frame, w)
    else:
        r = centered_loss(frame, mode, w)
        if r is None:
            return None
        L, valid, _ = r
    f = fit_affine(L[valid], D[valid], w[valid])
    return (float(f["a"]), float(f["b"])) if f["converged"] else None


def within_week(frame, B, log):
    n = len(frame)
    point = {m: e3_fit(frame, m, np.ones(n)) for m in E3_MODES}
    rng = np.random.default_rng([R.SEED, 23])
    wk = frame.week_end.to_numpy()
    dr = {m: np.full(B, np.nan) for m in E3_MODES}
    for r in range(B):
        w = R.week_multiplicity(wk, rng)
        for m in E3_MODES:
            f = e3_fit(frame, m, w)
            if f is not None:
                dr[m][r] = f[0]
        if r % 500 == 499:
            log("E3 bootstrap", r + 1)
    out = {}
    for m in E3_MODES:
        red = 1 - dr[m] / dr["pooled"]
        out[m] = {"a": point[m][0], "b": point[m][1], "a_boot": ci_with_mcse(dr[m]),
                  "ci_excludes_0": bool(R.pct_ci(dr[m])[0] > 0),
                  "reduction_vs_pooled": {"point": 1 - point[m][0] / point["pooled"][0], **ci_with_mcse(red)}}
    out["primary"] = "week_leave_state_out"
    out["meta"] = {"B": B, "seed": [R.SEED, 23],
                   "definitions": {"week_leave_state_out": "centre = week mean excluding the row's own state; scale = pooled mean (primary)",
                                   "week": "centre = week mean including the row (paper definition); scale = pooled mean",
                                   "week_review": "centre and scale = week mean (review definition)",
                                   "season": "centre = season mean; scale = pooled mean",
                                   "pooled": "centre and scale = pooled mean (E1 loss)"}}
    return out


def error_components(f):
    a, b = f["a"], f["b"]
    V = np.asarray(f["V2"])
    tq = stats.t.ppf(.975, f["t_df"])
    rows = []
    for C in C_GRID:
        g = np.array([1.0, 1.0 / C])
        tot = a + b / C
        se = float(np.sqrt(g @ V @ g))
        rows.append({"C": C, "a": a, "b_over_C": b / C, "total": tot, "total_ci_sandwich": [tot - tq * se, tot + tq * se],
                     "share_from_a": a / tot, "relative_rmse": float(np.sqrt(tot)),
                     "relative_rmse_ci_sandwich": [float(np.sqrt(max(tot - tq * se, 0))), float(np.sqrt(tot + tq * se))]})
    return rows


def sensitivities(base, log):
    """S-a1..S-a9. Returns {id: fit} with the label of each variant."""
    panel = build_panel(load_raw())
    win, _ = rule_windows()
    out = {}

    def run(key, frame, gi_mean=3.6, gi_sd=None, **meta):
        f = fit(frame, gi_mean, gi_sd)
        f.update(meta)
        out[key] = f
        R.append_jsonl(OUT_JSONL, {"block": "sensitivity", "id": key, **{k: v for k, v in f.items() if k != "V2"}})
        tw = f.get("inference", {}).get("two_way", {})
        log(key, "n", f["n_rows"], "a", f.get("a"), "ci", tw.get("a_ci"), "p", tw.get("p_two_sided"))

    for m in (3.2, 4.0):
        run(f"S-a1_gi_mean_{m}", load_sample("RA", gi_mean=m, gi_sd=1.6, panel=panel), m, 1.6, group="S-a1", gi=[m, 1.6])
    for sd in (1.2, 2.0):
        run(f"S-a2_gi_sd_{sd}", load_sample("RA", gi_mean=3.6, gi_sd=sd, panel=panel), 3.6, sd, group="S-a2", gi=[3.6, sd])
    for fl in (10, 20, 50):
        run(f"S-a3_floor_{fl}", base[base.C_t >= fl].reset_index(drop=True), group="S-a3", floor=fl)
    run("S-a4_calendar_windows", load_sample("RA", window="calendar", panel=panel), group="S-a4",
        windows=[[str(a_), str(b_)] for a_, b_ in CALENDAR])
    for th in (1.02, 1.10):
        w_th, _ = rule_windows(threshold=th)
        run(f"S-a5_threshold_{th}", load_sample("RA", windows=w_th, panel=panel), group="S-a5",
            windows={k: [str(v[0].date()), str(v[1].date())] for k, v in w_th.items()},
            identical_to_primary=bool(w_th == win))
    nxt = base.week_end + pd.Timedelta(days=7)
    in_pause = base.week_end.between(*PAUSE) | nxt.between(*PAUSE)
    run("S-a6_exclude_pause", base[~in_pause].reset_index(drop=True), group="S-a6", rows_removed=int(in_pause.sum()))
    for s in sorted(base.season.unique()):
        run(f"S-a7_drop_{s}", base[base.season != s].reset_index(drop=True), group="S-a7", dropped_season=s)
    cov = (base.coverage >= 90) & (base.coverage_next >= 90)
    run("S-a8_coverage_ge_90", base[cov].reset_index(drop=True), group="S-a8", rows_removed=int((~cov).sum()))
    weeks = np.sort(base.week_end.unique())
    loo = []
    for w in weeks:
        f = fit(base[base.week_end != w].reset_index(drop=True))
        loo.append({"week": str(pd.Timestamp(w).date()), "a": f.get("a"), "a_ci": f["inference"]["two_way"]["a_ci"],
                    "p": f["inference"]["two_way"]["p_two_sided"], "b_over_bP": f.get("b_over_bP")})
    a_l = np.array([r["a"] for r in loo])
    lo_l = np.array([r["a_ci"][0] for r in loo])
    s9 = {"group": "S-a9", "fits": loo, "a_range": [float(a_l.min()), float(a_l.max())],
          "lower_limit_min": float(lo_l.min()), "week_of_min_lower_limit": loo[int(lo_l.argmin())]["week"],
          "max_p": float(max(r["p"] for r in loo)), "b_over_bP_range": [min(r["b_over_bP"] for r in loo), max(r["b_over_bP"] for r in loo)]}
    R.append_jsonl(OUT_JSONL, {"block": "sensitivity", "id": "S-a9_leave_one_week_out", **{k: v for k, v in s9.items() if k != "fits"}})
    log("S-a9", "a range", s9["a_range"], "min lower limit", s9["lower_limit_min"], "max p", s9["max_p"])
    out["S-a9_leave_one_week_out"] = s9
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--B-wild", type=int, default=9999)
    ap.add_argument("--B-block", type=int, default=1999)
    a = ap.parse_args()
    log = Log(R.LOGS / "rev_intercept.log")
    t0 = time.time()
    if OUT_JSONL.exists():
        OUT_JSONL.unlink()
    base = load_sample("RA")
    prim = fit(base)
    tw = prim["inference"]["two_way"]
    R.append_jsonl(OUT_JSONL, {"block": "E1_primary", **{k: v for k, v in prim.items() if k != "V2"}})
    log("E1 primary", "a", prim["a"], "se", tw["se_a"], "t", tw["t"], "df", prim["t_df"], "p", tw["p_two_sided"], "ci", tw["a_ci"])
    sens = sensitivities(base, log)
    wild = wild_bootstrap_t(base, tw["t"], a.B_wild, log)
    R.append_jsonl(OUT_JSONL, {"block": "E1_wild_bootstrap_t", **wild})
    log("wild bootstrap-t p", wild["p_two_sided"], "mcse", wild["mcse"])
    blk, _ = block_bootstrap(base, a.B_block, log)
    R.append_jsonl(OUT_JSONL, {"block": "block_bootstrap", **blk})
    log("block bootstrap a CI", blk["a"]["ci"], "b/a", blk["b_over_a"]["ci"], "b/bP", blk["b_over_bP"]["ci"])
    e3 = within_week(base, a.B_block, log)
    R.append_jsonl(OUT_JSONL, {"block": "E3", **e3})
    log("E3", {m: (e3[m]["a"], e3[m]["a_boot"]["ci"]) for m in E3_MODES})

    lower = {k: (v["lower_limit_min"] if k.startswith("S-a9") else v["inference"]["two_way"]["a_ci"][0]) for k, v in sens.items()}
    failing = [k for k, v in lower.items() if not v > 0]
    rejects = tw["p_two_sided"] < .05
    verdict = ("demoted" if not rejects else "robust" if not failing else "significant but sensitivity-dependent")
    comps = error_components(prim)
    for row in comps:
        row["total_ci_block_bootstrap"] = blk["error_total"][str(row["C"])]["ci"]
    out = {"meta": {"script": "code/rev_intercept.py", "plan_sections": ["2", "4"], "data": "frozen NHSN 2026-10-20 vintage, R_A sample",
                    "seconds": round(time.time() - t0, 1)},
           "E1_primary": {k: v for k, v in prim.items() if k != "V2"} | {"V2_two_way_a_b": prim["V2"]},
           "E1_wild_bootstrap_t": wild, "block_bootstrap": blk,
           "sensitivities": sens, "sensitivity_lower_limits": lower,
           "robustness_rule": {"rule": "robust iff the primary test rejects and the 95% lower limit of a > 0 in every S-a1..S-a9 analysis",
                               "primary_rejects": rejects, "failing_sensitivities": failing, "verdict": verdict},
           "E3_within_week": e3, "error_components": comps,
           "crossover_b_over_a": {"point": prim["b_over_a"], "ci_sandwich_log_delta": tw["b_over_a_ci"],
                                  "ci_block_bootstrap": blk["b_over_a"]["ci"]}}
    R.write_json(OUT_JSON, out)
    log("rev_intercept done", "verdict", verdict, "failing", failing, "s", round(time.time() - t0, 1))


if __name__ == "__main__":
    main()
