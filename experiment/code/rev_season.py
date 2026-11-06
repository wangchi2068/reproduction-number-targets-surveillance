"""Season-intercept test with the locked primary calibration, all secondary calibrations, effect sizes and power
(analysis plan Section 3; review M3).

Statistic: W_CGM, two-way CGM sandwich Wald test of equal log season intercepts (common slope), as shipped.
Primary reference: calendar-week permutation of season labels (week counts 9/12/8/13 preserved),
B = 9,999, seed [20261105, 31], p = (1 + #{W* >= W_obs}) / (B + 1). Decision: significant only if p < 0.05.
Secondary: asymptotic chi2(3); shipped calibrations (generators with a fixed field, bootstrap covariance,
state-cluster wild bootstrap, within-state permutation); fresh-draw generators (rev_generators.json);
week-cluster wild bootstrap (null imposed, B = 9,999, seed [20261105, 32]); circular shifts of the season-label
sequence (42, exact). Data sensitivities S-a6, S-a8, S-a4 with the primary calibration (seed [20261105, 34, j]).
Effect sizes: season-stratified calendar-week block bootstrap (B = 1,999, seed [20261105, 35]) and sandwich delta CIs.
Power / minimum detectable ratio: plasmode on the observed losses (seed [20261105, 33]).

Each finished task is appended to analysis/rev_season.jsonl; the summary goes to analysis/rev_season.json.
Usage: python code/rev_season.py [--B 9999] [--B-boot 1999] [--workers 6]
"""
from __future__ import annotations

import argparse
import json
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy import stats

import ma_calibration as MC
import rev_common as R
from ma_models import Log
from revision_data import build_panel, load_raw, load_sample
from revision_stats import fit_loss, relative_loss
from s03c_state_fe import fit_group_affine

OUT_JSONL = R.AN / "rev_season.jsonl"
OUT_JSON = R.AN / "rev_season.json"
PAUSE = (pd.Timestamp("2025-10-03"), pd.Timestamp("2025-11-13"))
DELTAS = np.round(np.arange(0, 2.01, .25), 2)


class Design:
    def __init__(self, frame):
        self.f = frame
        self.x, self.D = frame.x.to_numpy(float), frame.D.to_numpy(float)
        self.st, self.wk = frame.state.to_numpy(), frame.week_end.to_numpy()
        self.se = frame.season.astype(str).to_numpy()
        self.loss, self.xbar = relative_loss(self.x)
        self.weeks, self.wi = np.unique(self.wk, return_inverse=True)
        self.week_season = np.empty(len(self.weeks), dtype=object)
        self.week_season[self.wi] = self.se
        self.n = len(self.x)

    def W(self, labels):
        g = MC.sandwich(self.x, self.D, self.st, self.wk, groups=labels)
        return np.nan if g is None else MC.season_wald(g["beta"], g["V2"])

    def W_loss(self, loss_b, labels):
        """W from a given loss vector, xbar held at the observed mean (as in the review's wald_from_loss)."""
        names, codes = np.unique(labels, return_inverse=True)
        fb = fit_group_affine(loss_b, self.D, labels)
        if not fb["converged"] or not np.all(fb["active"]) or not np.all(np.asarray(fb["intercepts"]) > 0):
            return np.nan
        beta = np.r_[np.asarray(fb["intercepts"], float), fb["b"]]
        X = np.column_stack([np.eye(len(names))[codes], 1 / self.D])
        V, _ = R.cgm_from_loss(self.x, self.xbar, loss_b, self.D, X, beta, self.st, self.wk)
        return np.nan if V is None else MC.season_wald(beta, V)

    def perm_labels(self, rng):
        return rng.permutation(self.week_season)[self.wi]


def frame_for(key):
    base = load_sample("RA")
    if key == "primary":
        return base
    if key == "S-a6_exclude_pause":
        nxt = base.week_end + pd.Timedelta(days=7)
        return base[~(base.week_end.between(*PAUSE) | nxt.between(*PAUSE))].reset_index(drop=True)
    if key == "S-a8_coverage_ge_90":
        return base[(base.coverage >= 90) & (base.coverage_next >= 90)].reset_index(drop=True)
    if key == "S-a4_calendar_windows":
        return load_sample("RA", window="calendar", panel=build_panel(load_raw()))
    raise ValueError(key)


def p_upper(null, obs):
    v = np.asarray(null, float)
    v = v[np.isfinite(v)]
    p = float((1 + np.sum(v >= obs)) / (len(v) + 1))
    return {"p": p, "mcse": float(np.sqrt(p * (1 - p) / len(v))), "valid": int(len(v)),
            "null_q95": float(np.quantile(v, .95)), "null_q99": float(np.quantile(v, .99))}


def task_permutation(args):
    key, seed, B = args
    t0 = time.time()
    d = Design(frame_for(key))
    W_obs = d.W(d.se)
    rng = np.random.default_rng(seed)
    Wp = np.array([d.W(d.perm_labels(rng)) for _ in range(B)])
    np.save(R.AN / f"rev_season_null_{key}.npy", Wp)
    out = {"task": "permutation", "data": key, "W_obs": W_obs, "df": len(np.unique(d.se)) - 1, "n_rows": d.n,
           "n_weeks": len(d.weeks), "weeks_per_season": {s: int(np.sum(d.week_season == s)) for s in np.unique(d.se)},
           "asymptotic_p": float(stats.chi2.sf(W_obs, len(np.unique(d.se)) - 1)),
           "B": B, "seed": list(seed), **p_upper(Wp, W_obs), "seconds": round(time.time() - t0, 1)}
    R.append_jsonl(OUT_JSONL, out)
    return out


def task_wild(args):
    _, seed, B = args
    t0 = time.time()
    d = Design(frame_for("primary"))
    W_obs = d.W(d.se)
    f0 = fit_loss(d.x, d.D)
    mu0 = f0["a"] + f0["b"] / d.D
    res0 = d.loss - mu0
    rng = np.random.default_rng(seed)
    Wb = np.full(B, np.nan)
    for r in range(B):
        v = R.webb_weights(rng, len(d.weeks))[d.wi]
        try:
            Wb[r] = d.W_loss(np.maximum(mu0 + v * res0, 1e-7), d.se)
        except np.linalg.LinAlgError:
            pass
    np.save(R.AN / "rev_season_null_wild_week.npy", Wb)
    out = {"task": "wild_week", "W_obs": W_obs, "B": B, "seed": list(seed), **p_upper(Wb, W_obs),
           "construction": "Webb weights per calendar week; null fit = pooled (equal-intercept) affine fit; losses clamped at 1e-7; xbar fixed",
           "seconds": round(time.time() - t0, 1)}
    R.append_jsonl(OUT_JSONL, out)
    return out


def task_circular(args):
    t0 = time.time()
    d = Design(frame_for("primary"))
    W_obs = d.W(d.se)
    seq = d.week_season  # ordered by calendar date
    Ws = np.array([d.W(np.roll(seq, k)[d.wi]) for k in range(len(seq))])
    p = float(np.mean(Ws >= W_obs - 1e-9))
    out = {"task": "circular_shift", "W_obs": W_obs, "shifts": len(seq), "W_by_shift": Ws.tolist(), "p_exact": p,
           "min_attainable_p": 1 / len(seq), "rank_of_observed": int(np.sum(Ws >= W_obs - 1e-9)),
           "seconds": round(time.time() - t0, 1)}
    R.append_jsonl(OUT_JSONL, out)
    return out


def task_effects(args):
    _, seed, B = args
    t0 = time.time()
    d = Design(frame_for("primary"))
    g = MC.sandwich(d.x, d.D, d.st, d.wk, groups=d.se)
    names = np.unique(d.se)
    S = len(names)
    a_s, V = g["beta"][:S], g["V2"]
    a_p = fit_loss(d.x, d.D)["a"]
    df = len(d.weeks) - 1
    tq = stats.t.ppf(.975, df)
    se_log = np.sqrt(np.diag(V)[:S]) / a_s
    post = names == "2023-24"
    c = np.where(post, 1.0, -1.0 / (S - 1))
    gd = c / a_s
    con = float(c @ np.log(a_s))
    se_con = float(np.sqrt(gd @ V[:S, :S] @ gd))
    rng = np.random.default_rng(seed)
    la = np.full((B, S), np.nan)
    ap = np.full(B, np.nan)
    for r in range(B):
        w = R.week_multiplicity(d.wk, rng, strata=d.se)
        loss, _ = relative_loss(d.x, w)
        fb = fit_group_affine(loss, d.D, d.se, w)
        a = np.asarray(fb["intercepts"], float)
        if fb["converged"] and np.all(a > 0):
            la[r] = np.log(a)
        fp = fit_loss(d.x, d.D, w)
        if fp["converged"] and fp["a"] > 0:
            ap[r] = fp["a"]
    rat = np.exp(la) / ap[:, None]
    mm = np.exp(la.max(1) - la.min(1))
    cb = la @ c
    seasons = {}
    for j, s in enumerate(names):
        seasons[s] = {"a_s": float(a_s[j]), "ci_sandwich_log_delta": [float(a_s[j] * np.exp(-tq * se_log[j])), float(a_s[j] * np.exp(tq * se_log[j]))],
                      "ci_bootstrap": [float(np.exp(q)) for q in R.pct_ci(la[:, j])],
                      "ratio_to_pooled": float(a_s[j] / a_p), "ratio_to_pooled_ci_bootstrap": R.pct_ci(rat[:, j]),
                      "weeks": int(np.sum(d.week_season == s)), "rows": int(np.sum(d.se == s)),
                      "mean_reporting_coverage": float(d.f.loc[d.se == s, "coverage"].mean())}
    out = {"task": "effects", "seasons": seasons, "a_pooled": a_p, "common_slope_b": float(g["beta"][S]),
           "max_over_min": {"point": float(a_s.max() / a_s.min()), "ci_bootstrap": R.pct_ci(mm),
                            "note": "max/min of estimated intercepts is biased upward"},
           "post_hoc_contrast_2023_24_vs_mean_of_others_log": {"point": con, "ratio": float(np.exp(con)),
                                                               "ci_sandwich": [float(np.exp(con - tq * se_con)), float(np.exp(con + tq * se_con))],
                                                               "ci_bootstrap": [float(np.exp(q)) for q in R.pct_ci(cb)],
                                                               "label": "post hoc: chosen after seeing the data; no test"},
           "B": B, "seed": list(seed), "valid": int(np.isfinite(la).all(1).sum()), "t_df": df,
           "seconds": round(time.time() - t0, 1)}
    R.append_jsonl(OUT_JSONL, out)
    return out


def task_power(args):
    _, seed, n_crit, n_draw = args
    t0 = time.time()
    d = Design(frame_for("primary"))
    names, codes = np.unique(d.se, return_inverse=True)
    g = MC.sandwich(d.x, d.D, d.st, d.wk, groups=d.se)
    S = len(names)
    a_s, b_s = g["beta"][:S], g["beta"][S]
    f0 = fit_loss(d.x, d.D)
    a, b = f0["a"], f0["b"]
    L0 = d.loss * (a + b / d.D) / (a_s[codes] + b_s / d.D)
    cpat = np.log(a_s) - np.log(a_s).mean()
    rng = np.random.default_rng(seed)
    crit_draws = np.array([d.W_loss(L0, d.perm_labels(rng)) for _ in range(n_crit)])
    crit = float(np.nanquantile(crit_draws, .95))
    curve = []
    for delta in DELTAS:
        Ws = np.full(n_draw, np.nan)
        for r in range(n_draw):
            lab = d.perm_labels(rng)
            _, pc = np.unique(lab, return_inverse=True)
            Ld = L0 * (a * np.exp(delta * cpat[pc]) + b / d.D) / (a + b / d.D)
            Ws[r] = d.W_loss(Ld, lab)
        ok = np.isfinite(Ws)
        pw = float(np.mean(Ws[ok] >= crit))
        curve.append({"delta": float(delta), "max_over_min_ratio": float(np.exp(delta * (cpat.max() - cpat.min()))), "power": pw,
                      "mcse": float(np.sqrt(pw * (1 - pw) / ok.sum())), "valid": int(ok.sum())})
        R.append_jsonl(OUT_JSONL, {"task": "power_point", **curve[-1]})
    pw = np.array([c["power"] for c in curve])
    rr = np.array([c["max_over_min_ratio"] for c in curve])
    mdr = None
    hit = np.flatnonzero(pw >= .8)
    if len(hit):
        j = hit[0]
        mdr = float(rr[j]) if j == 0 else float(rr[j - 1] + (.8 - pw[j - 1]) * (rr[j] - rr[j - 1]) / (pw[j] - pw[j - 1]))
    out = {"task": "power", "crit_W_perm_95": crit, "n_crit": n_crit, "draws_per_delta": n_draw, "curve": curve,
           "observed_pattern_centred_log": dict(zip(names.tolist(), cpat.tolist())),
           "observed_max_over_min": float(a_s.max() / a_s.min()),
           "mdr_max_over_min_at_80pct_power": mdr, "power_at_observed_pattern_delta_1": float(pw[DELTAS == 1.0][0]),
           "seed": list(seed), "seconds": round(time.time() - t0, 1),
           "construction": "null-imposed losses L0 = L (a + b/C) / (a_s + b_hat/C); L_delta = L0 (a exp(delta c_s) + b/C) / (a + b/C); "
                           "season labels by random calendar-week permutation; loss-level W with xbar fixed"}
    R.append_jsonl(OUT_JSONL, out)
    return out


def run_task(t):
    return {"perm": task_permutation, "wild": task_wild, "circ": task_circular, "eff": task_effects, "pow": task_power}[t[0]](t[1])


def holm(ps):
    keys = sorted(ps, key=ps.get)
    out, run = {}, 0.0
    for i, k in enumerate(keys):
        run = max(run, min(1.0, (len(keys) - i) * ps[k]))
        out[k] = run
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=9999)
    ap.add_argument("--B-boot", type=int, default=1999)
    ap.add_argument("--power-crit", type=int, default=2000)
    ap.add_argument("--power-draws", type=int, default=400)
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    log = Log(R.LOGS / "rev_season.log")
    t0 = time.time()
    if OUT_JSONL.exists():
        OUT_JSONL.unlink()
    tasks = [("perm", ("primary", (R.SEED, 31), a.B)),
             ("wild", (None, (R.SEED, 32), a.B)),
             ("pow", (None, (R.SEED, 33), a.power_crit, a.power_draws)),
             ("perm", ("S-a6_exclude_pause", (R.SEED, 34, 1), a.B)),
             ("perm", ("S-a8_coverage_ge_90", (R.SEED, 34, 2), a.B)),
             ("perm", ("S-a4_calendar_windows", (R.SEED, 34, 3), a.B)),
             ("eff", (None, (R.SEED, 35), a.B_boot)),
             ("circ", None)]
    log("rev_season start", len(tasks), "tasks", "B", a.B)
    res = {}
    with Pool(min(a.workers, len(tasks))) as pool:
        for r in pool.imap_unordered(run_task, tasks):
            key = r["task"] + ("_" + r["data"] if r["task"] == "permutation" else "")
            res[key] = r
            log("finished", key, {k: r.get(k) for k in ("p", "mcse", "p_exact", "W_obs", "mdr_max_over_min_at_80pct_power", "seconds")})
    cal = json.loads((R.AN / "ma_calibration.json").read_text())
    wb = json.loads((R.AN / "ma_wild_bootstrap.json").read_text())
    gen = json.loads((R.AN / "rev_generators.json").read_text())
    prim = res["permutation_primary"]
    shipped = {nm: {"p": g["season_sandwich"]["p_observed"], "mcse": g["season_sandwich"]["p_observed_mcse"],
                    "panels": g["season_sandwich"]["panels"], "field": "fixed transmission field (shipped)",
                    "post_hoc_law": nm == "S9_national"} for nm, g in cal["generators"].items()}
    fresh = {nm: {**g["season_test_generator_p"], "field": "redrawn per replicate", "post_hoc_law": g["law"]["post_hoc_law"]}
             for nm, g in gen["generators"].items()}
    table = [
        {"calibration": "Calendar-week permutation (PRIMARY, locked)", "p": prim["p"], "mcse": prim["mcse"], "draws": prim["valid"], "source": "rev_season.py"},
        {"calibration": "Asymptotic chi2(3), sandwich W", "p": prim["asymptotic_p"], "mcse": None, "draws": None, "source": "rev_season.py"},
        {"calibration": "Asymptotic chi2(3), bootstrap-covariance W (shipped)", "p": cal["season_test"]["chi2_p_bootstrap"], "mcse": None, "draws": None,
         "source": "ma_calibration.py"},
        {"calibration": "Week-cluster wild bootstrap (Webb, null imposed)", "p": res["wild_week"]["p"], "mcse": res["wild_week"]["mcse"],
         "draws": res["wild_week"]["valid"], "source": "rev_season.py"},
        {"calibration": "Circular shifts of the season sequence (exact, 42)", "p": res["circular_shift"]["p_exact"], "mcse": None, "draws": 42,
         "source": "rev_season.py"},
        {"calibration": "State-cluster wild bootstrap (shipped)", "p": wb["wild_cluster_bootstrap"]["p_value"], "mcse": wb["wild_cluster_bootstrap"]["p_value_mcse"],
         "draws": wb["wild_cluster_bootstrap"]["n_draws"], "source": "ma_wild_bootstrap.py"},
        {"calibration": "Within-state permutation (shipped)", "p": wb["permutation_test"]["p_value"], "mcse": wb["permutation_test"]["p_value_mcse"],
         "draws": wb["permutation_test"]["n_draws"], "source": "ma_wild_bootstrap.py"},
    ]
    for nm, v in shipped.items():
        table.append({"calibration": f"Generator {nm} ({v['field']}{'; post hoc law' if v['post_hoc_law'] else ''})", "p": v["p"], "mcse": v["mcse"],
                      "draws": v["panels"], "source": "ma_calibration.py"})
    for nm, v in fresh.items():
        table.append({"calibration": f"Generator {nm} ({v['field']}{'; post hoc law' if v['post_hoc_law'] else ''})", "p": v["p"], "mcse": v["mcse"],
                      "draws": v["n"], "source": "rev_generators.py"})
    ri = json.loads((R.AN / "rev_intercept.json").read_text())
    p_e1 = ri["E1_primary"]["inference"]["two_way"]["p_two_sided"]
    hp = holm({"E1_intercept": p_e1, "season_primary": prim["p"]})
    sig = prim["p"] < .05
    out = {"meta": {"script": "code/rev_season.py", "plan_section": "3", "B": a.B, "seconds": round(time.time() - t0, 1)},
           "W_obs": prim["W_obs"], "primary": prim,
           "decision": {"rule": "season difference significant only if the primary calendar-week permutation p < 0.05",
                        "significant": sig, "verdict": "significant" if sig else "not significant: season statements are exploratory/descriptive"},
           "holm": {"raw": {"E1_intercept": p_e1, "season_primary": prim["p"]}, "adjusted": hp},
           "calibration_table": table,
           "data_sensitivities": {k: res[f"permutation_{k}"] for k in ("S-a6_exclude_pause", "S-a8_coverage_ge_90", "S-a4_calendar_windows")},
           "wild_week": res["wild_week"], "circular_shift": res["circular_shift"], "effects": res["effects"], "power": res["power"],
           "reporting_context": {"voluntary_reporting": "May-Oct 2024 (between the 2023-24 and 2024-25 windows)",
                                 "data_channel_pause": "Oct-Nov 2025 (inside the 2025-26 window)",
                                 "interpretation": "season intercepts are transmission-plus-reporting totals (Proposition 2)"}}
    R.write_json(OUT_JSON, out)
    log("rev_season done", "primary p", prim["p"], "significant", sig, "holm", hp)


if __name__ == "__main__":
    main()
