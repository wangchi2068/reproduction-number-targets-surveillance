"""Functional form, high-count floor and population size (analysis_plan_amendment2.md, Section B; review N2).

B.1 observed mean loss by count decile against affine, power-law and three-parameter fits
    (calendar-week block bootstrap, B = 1,999, seed [20261105, 70, 1]).
B.2 population: a0 + a1 z_pop + b/C (two-way CGM sandwich, t 41 df; state-cluster and week bootstraps,
    seeds [.., 70, 2] and [.., 70, 3]); log-link and b-varying secondary models; population tertiles
    (separate fits, t 16 df; tertile intercepts with a common slope, state-cluster bootstrap seed [.., 70, 4]).
B.3 a + b C^-gamma (week bootstrap B = 1,999 seed [.., 70, 5]; state bootstrap B = 999 seed [.., 70, 6]).
B.4 model-free floor at D >= P75 / P90 / P95 (week bootstrap B = 1,999, seed [.., 70, 7]).

Output: analysis/rev_r2_floor.json (+ rev_r2_floor_draws.npz). Usage: python code/rev_r2_floor.py [--B 1999] [--workers 4]
"""
from __future__ import annotations

import argparse
import time
from multiprocessing import Pool

import numpy as np
from scipy import stats
from scipy.optimize import minimize

import ma_calibration as MC
import rev_common as R
import rev_r2_lib as L2
from ma_models import Log
from revision_data import load_sample
from revision_stats import fit_affine, relative_loss
from rev_intercept import fit as fit_sample
from s03b_thinning import poisson_reference
from s03c_state_fe import fit_group_affine

OUT = R.AN / "rev_r2_floor.json"
QS_FLOOR = (.75, .90, .95)


def data():
    f = load_sample("RA", 3.6)
    pop = L2.population()
    f = f.merge(pop, on="state", how="left")
    assert f.z_pop.notna().all() and len(f) == 1748
    return f


def deciles(D):
    cuts = np.quantile(D, np.arange(1, 10) / 10)
    return np.searchsorted(cuts, D, side="left"), cuts


def state_multiplicity(state, rng):
    names, si = np.unique(state, return_inverse=True)
    m = np.bincount(rng.integers(0, len(names), len(names)), minlength=len(names))
    return m[si].astype(float)


def fit_pop_linear(loss, D, z, w=None, b_varies=False):
    X = np.column_stack([np.ones(len(D)), z, 1 / D] + ([z / D] if b_varies else []))
    return L2.fit_linear_ql(loss, X, w), X


def fit_pop_log(loss, D, z, w=None, start=None):
    w = np.ones(len(D)) if w is None else w
    def obj(p):
        return L2.qdev(loss, np.exp(p[0] + p[1] * z) + p[2] / D, w)
    st = start if start is not None else [np.log(.015), 0.0, .56]
    o = minimize(obj, st, method="Nelder-Mead", options={"xatol": 1e-10, "fatol": 1e-12, "maxiter": 6000})
    return {"converged": bool(o.success or np.isfinite(o.fun)), "alpha0": float(o.x[0]), "theta": float(o.x[1]),
            "b": float(o.x[2]), "qdev": float(o.fun)}


def job(args):
    kind, B, seed = args
    f = data()
    x, D, st, wk, z = (f.x.to_numpy(float), f.D.to_numpy(float), f.state.to_numpy(), f.week_end.to_numpy(),
                       f.z_pop.to_numpy(float))
    rng = np.random.default_rng([R.SEED, 70, seed])
    dec, _ = deciles(D)
    thr = np.quantile(D, QS_FLOOR)
    tert_lab = tertile(f)
    log = Log(R.LOGS / f"rev_r2_floor_{kind}.log")
    out = {}
    for r in range(B):
        w = state_multiplicity(st, rng) if kind in ("pop_state", "tertile_state", "three_state") else R.week_multiplicity(wk, rng)
        loss, _ = relative_loss(x, w)
        rec = {}
        if kind in ("deciles", "floor", "three_week", "three_state"):
            fa = fit_affine(loss, D, w)
            rec["a"], rec["b"] = (fa["a"], fa["b"]) if fa["converged"] else (np.nan, np.nan)
        if kind == "deciles":
            fp = L2.fit_power(loss, D, w)
            ft = L2.fit_three(loss, D, w)
            for k in range(10):
                m = (dec == k) & (w > 0)
                rec[f"meanL_{k}"] = float(np.average(loss[m], weights=w[m])) if m.any() else np.nan
                mk = dec == k
                rec[f"aff_{k}"] = float(np.mean(rec["a"] + rec["b"] / D[mk]))
                rec[f"pow_{k}"] = float(np.mean(fp["c"] * D[mk] ** (-fp["gamma"])))
                rec[f"three_{k}"] = float(np.mean(ft["a"] + ft["b"] * D[mk] ** (-ft["gamma"]))) if ft["converged"] else np.nan
        if kind in ("three_week", "three_state"):
            ft = L2.fit_three(loss, D, w)
            if ft["converged"]:
                rec.update(t_a=ft["a"], t_b=ft["b"], t_gamma=ft["gamma"],
                           dq=L2.qdev(loss[w > 0], (rec["a"] + rec["b"] / D)[w > 0], w[w > 0]) - ft["qdev"]
                           if np.isfinite(rec["a"]) else np.nan)
        if kind in ("pop_state", "pop_week"):
            fl, _ = fit_pop_linear(loss, D, z, w)
            if fl["converged"]:
                rec["a0"], rec["a1"], rec["b"] = fl["beta"]
            fg = fit_pop_log(loss, D, z, w)
            rec["theta"] = fg["theta"]
        if kind == "tertile_state":
            fg = fit_group_affine(loss, D, tert_lab, w)
            if fg["converged"]:
                ic = np.asarray(fg["intercepts"], float)
                rec.update({f"tert_{k}": ic[k] for k in range(3)})
                rec["tert_ratio"] = ic[2] / ic[0] if ic[0] > 0 else np.nan
        if kind == "floor":
            bP = poisson_reference(f, w)
            for q, t in zip(QS_FLOOR, thr):
                m = (D >= t) & (w > 0)
                mL = float(np.average(loss[m], weights=w[m]))
                rec[f"meanL_q{int(q * 100)}"] = mL
                rec[f"rrmse_q{int(q * 100)}"] = float(np.sqrt(mL))
                rec[f"floor_b_q{int(q * 100)}"] = mL - float(np.average(rec["b"] / D[m], weights=w[m]))
                rec[f"floor_bP_q{int(q * 100)}"] = mL - float(np.average(bP / D[m], weights=w[m]))
        for k, v in rec.items():
            out.setdefault(k, np.full(B, np.nan))[r] = np.nan if v is None else v
        if r % 200 == 199:
            log(kind, r + 1)
    return kind, out


def tertile(f):
    """Population tertile per row: states ranked by population, split 17/17/17 (0 = smallest)."""
    ps = f.groupby("state")["pop"].first().sort_values()
    lab = {s: i // 17 for i, s in enumerate(ps.index)}
    return f.state.map(lab).to_numpy()


def summ(v, point=None):
    v = np.asarray(v, float)
    o = {"ci": R.pct_ci(v), "n_finite": int(np.isfinite(v).sum())}
    if point is not None:
        o["point"] = point
    return o


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=1999)
    ap.add_argument("--B-state-three", type=int, default=999)
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    L2.register_script(__file__)
    t0 = time.time()
    log = Log(R.LOGS / "rev_r2_floor.log")
    f = data()
    f["tertile_label"] = tertile(f)
    x, D, st, wk, z = (f.x.to_numpy(float), f.D.to_numpy(float), f.state.to_numpy(), f.week_end.to_numpy(),
                       f.z_pop.to_numpy(float))
    loss, xbar = relative_loss(x)
    fa = fit_affine(loss, D)
    fp = L2.fit_power(loss, D)
    ft = L2.fit_three(loss, D)
    q_aff = L2.qdev(loss, fa["a"] + fa["b"] / D)
    log("point fits", fa["a"], fa["b"], fp, ft)
    jobs = [("deciles", args.B, 1), ("pop_state", args.B, 2), ("pop_week", args.B, 3), ("tertile_state", args.B, 4),
            ("three_week", args.B, 5), ("three_state", args.B_state_three, 6), ("floor", args.B, 7)]
    with Pool(args.workers) as pool:
        draws = dict(pool.map(job, jobs))
    res = {"meta": {"script": "code/rev_r2_floor.py", "plan": "analysis/analysis_plan_amendment2.md, Section B",
                    "sample": "R_A, 1,748 rows", "B": args.B, "B_state_three": args.B_state_three,
                    "seeds": {"deciles": [R.SEED, 70, 1], "pop_state": [R.SEED, 70, 2], "pop_week": [R.SEED, 70, 3],
                              "tertile_state": [R.SEED, 70, 4], "three_week": [R.SEED, 70, 5],
                              "three_state": [R.SEED, 70, 6], "floor": [R.SEED, 70, 7]},
                    "population_source": {"url": L2.POP_URL, "file": "data/external/census_NST-EST2025-ALLDATA.csv",
                                          "vintage": 2025, "covariate": "log of mean POPESTIMATE2022-2025, z-scored over 51 jurisdictions"}}}
    # ---- B.1 deciles
    dec, cuts = deciles(D)
    dd = draws["deciles"]
    rows = []
    for k in range(10):
        m = dec == k
        mL = float(loss[m].mean())
        rows.append({"decile": k + 1, "n": int(m.sum()), "C_min": float(D[m].min()), "C_max": float(D[m].max()),
                     "mean_C": float(D[m].mean()), "mean_L": mL, "mean_L_ci": R.pct_ci(dd[f"meanL_{k}"]),
                     "rel_rmse": float(np.sqrt(mL)), "rel_rmse_ci": [float(np.sqrt(max(v, 0))) for v in R.pct_ci(dd[f"meanL_{k}"])],
                     "affine": float(np.mean(fa["a"] + fa["b"] / D[m])), "affine_band": R.pct_ci(dd[f"aff_{k}"]),
                     "power": float(np.mean(fp["c"] * D[m] ** (-fp["gamma"]))), "power_band": R.pct_ci(dd[f"pow_{k}"]),
                     "three": float(np.mean(ft["a"] + ft["b"] * D[m] ** (-ft["gamma"]))), "three_band": R.pct_ci(dd[f"three_{k}"])})
        rows[-1]["inside_affine_band"] = bool(rows[-1]["affine_band"][0] <= mL <= rows[-1]["affine_band"][1])
        rows[-1]["affine_inside_meanL_ci"] = bool(rows[-1]["mean_L_ci"][0] <= rows[-1]["affine"] <= rows[-1]["mean_L_ci"][1])
    res["B1_deciles"] = {"cut_points": cuts.tolist(), "rows": rows,
                         "fits": {"affine": {"a": fa["a"], "b": fa["b"], "qdev": q_aff},
                                  "power": fp, "three_parameter": ft},
                         "deciles_where_affine_outside_observed_ci": [r["decile"] for r in rows if not r["affine_inside_meanL_ci"]]}
    # ---- B.2 population
    fl, X = fit_pop_linear(loss, D, z)
    beta = fl["beta"]
    V, _ = R.cgm_from_loss(x, xbar, loss, D, X, beta, st, wk)
    se = np.sqrt(np.diag(V))
    tc = stats.t.ppf(.975, 41)
    a1_ci = [float(beta[1] - tc * se[1]), float(beta[1] + tc * se[1])]
    zs = f.groupby("state").z_pop.first().to_numpy()
    z10, z90 = np.quantile(zs, [.1, .9])
    ratio9010 = (beta[0] + beta[1] * z90) / (beta[0] + beta[1] * z10)
    if a1_ci[1] < 0:
        verdict = "a falls significantly with population"
    elif a1_ci[0] > 0:
        verdict = "a rises significantly with population"
    else:
        verdict = "no population gradient detected"
    flb, Xb = fit_pop_linear(loss, D, z, b_varies=True)
    Vb, _ = R.cgm_from_loss(x, xbar, loss, D, Xb, flb["beta"], st, wk)
    seb = np.sqrt(np.diag(Vb))
    fg = fit_pop_log(loss, D, z)
    ps, pw = draws["pop_state"], draws["pop_week"]
    lp = f.groupby("state").log_pop.first()
    res["B2_population"] = {
        "covariate": {"log_pop_mean": float(lp.mean()), "log_pop_sd": float(lp.std(ddof=0)),
                      "z_p10": float(z10), "z_p90": float(z90), "pop_p10": float(np.exp(lp.quantile(.1))),
                      "pop_p90": float(np.exp(lp.quantile(.9)))},
        "primary_linear": {"a0": float(beta[0]), "a1": float(beta[1]), "b": float(beta[2]), "se": se.tolist(),
                           "a1_ci_two_way_t41": a1_ci, "t_a1": float(beta[1] / se[1]),
                           "p_a1_two_sided": float(2 * stats.t.sf(abs(beta[1] / se[1]), 41)),
                           "a0_ci_two_way_t41": [float(beta[0] - tc * se[0]), float(beta[0] + tc * se[0])],
                           "a_at_pop_p10": float(beta[0] + beta[1] * z10), "a_at_pop_p90": float(beta[0] + beta[1] * z90),
                           "ratio_a_p90_over_p10": float(ratio9010),
                           "a1_ci_state_bootstrap": R.pct_ci(ps["a1"]), "a1_ci_week_bootstrap": R.pct_ci(pw["a1"]),
                           "ratio_ci_state_bootstrap": R.pct_ci((ps["a0"] + ps["a1"] * z90) / (ps["a0"] + ps["a1"] * z10)),
                           "verdict": verdict},
        "secondary_log_link": {**fg, "theta_ci_state_bootstrap": R.pct_ci(ps["theta"]), "theta_ci_week_bootstrap": R.pct_ci(pw["theta"]),
                               "ratio_a_p90_over_p10": float(np.exp(fg["theta"] * (z90 - z10)))},
        "secondary_b_varies": {"beta": flb["beta"].tolist(), "se": seb.tolist(), "names": ["a0", "a1", "b0", "b1"],
                               "a1_ci_two_way_t41": [float(flb["beta"][1] - tc * seb[1]), float(flb["beta"][1] + tc * seb[1])],
                               "b1_ci_two_way_t41": [float(flb["beta"][3] - tc * seb[3]), float(flb["beta"][3] + tc * seb[3])]},
    }
    tert = []
    for k in range(3):
        sub = f[f.tertile_label == k]
        fs = fit_sample(sub)
        tert.append({"tertile": k + 1, "states": sorted(sub.state.unique().tolist()), "n_rows": fs["n_rows"],
                     "pop_range": [float(sub["pop"].min()), float(sub["pop"].max())], "median_C": float(sub.D.median()),
                     "a": fs.get("a"), "a_ci_two_way": fs["inference"]["two_way"]["a_ci"] if fs["converged"] else None,
                     "b": fs.get("b"), "b_ci_two_way": fs["inference"]["two_way"]["b_ci"] if fs["converged"] else None,
                     "b_over_a": fs.get("b_over_a"), "b_over_a_ci": fs["inference"]["two_way"]["b_over_a_ci"] if fs["converged"] else None,
                     "t_df": fs.get("t_df"), "rel_rmse_at_C": {str(C): float(np.sqrt(fs["a"] + fs["b"] / C)) for C in (37, 100, 390)}
                     if fs["converged"] and fs["a"] + fs["b"] / 390 > 0 else None})
    sg = MC.sandwich(x, D, st, wk, groups=f.tertile_label.to_numpy())
    td = draws["tertile_state"]
    res["B2_population"]["tertiles_separate"] = tert
    res["B2_population"]["tertiles_common_slope"] = {
        "intercepts": sg["beta"][:3].tolist(), "b": float(sg["beta"][3]),
        "se_two_way": np.sqrt(np.diag(sg["V2"]))[:4].tolist(),
        "ratio_high_over_low": float(sg["beta"][2] / sg["beta"][0]),
        "ratio_ci_state_bootstrap": R.pct_ci(td["tert_ratio"]),
        "intercept_ci_state_bootstrap": [R.pct_ci(td[f"tert_{k}"]) for k in range(3)]}
    # ---- B.3 three-parameter
    tw, ts = draws["three_week"], draws["three_state"]
    a_w = tw["t_a"]
    res["B3_three_parameter"] = {
        "point": ft, "affine_qdev_minus_three": float(q_aff - ft["qdev"]),
        "week_bootstrap": {"a": summ(tw["t_a"]), "b": summ(tw["t_b"]), "gamma": summ(tw["t_gamma"]),
                           "share_a_le_0": float(np.nanmean(a_w <= 0)), "share_a_lt_0.001": float(np.nanmean(a_w < .001)),
                           "share_qdev_gain_gt_2": float(np.nanmean(tw["dq"] > 2)), "share_gamma_at_upper_bound": float(np.nanmean(tw["t_gamma"] >= 2.999))},
        "state_bootstrap": {"a": summ(ts["t_a"]), "gamma": summ(ts["t_gamma"]), "share_a_le_0": float(np.nanmean(ts["t_a"] <= 0))},
    }
    gci = R.pct_ci(tw["t_gamma"])
    aci = R.pct_ci(tw["t_a"])
    res["B3_three_parameter"]["gamma_1_compatible"] = bool(gci[0] <= 1 <= gci[1])
    res["B3_three_parameter"]["floor_survives_free_exponent"] = bool(aci[0] > 0)
    # ---- B.4 floor
    fd = draws["floor"]
    bP = poisson_reference(f)
    thr = np.quantile(D, QS_FLOOR)
    fl_rows = []
    for q, t in zip(QS_FLOOR, thr):
        m = D >= t
        mL = float(loss[m].mean())
        tag = f"q{int(q * 100)}"
        fb = mL - float(np.mean(fa["b"] / D[m]))
        fP = mL - float(np.mean(bP / D[m]))
        fl_rows.append({"threshold_quantile": q, "threshold_C": float(t), "n": int(m.sum()), "mean_C": float(D[m].mean()),
                        "mean_L": mL, "mean_L_ci": R.pct_ci(fd[f"meanL_{tag}"]),
                        "rel_rmse": float(np.sqrt(mL)), "rel_rmse_ci": R.pct_ci(fd[f"rrmse_{tag}"]),
                        "floor_minus_bhat": fb, "floor_minus_bhat_ci": R.pct_ci(fd[f"floor_b_{tag}"]),
                        "floor_minus_bP": fP, "floor_minus_bP_ci": R.pct_ci(fd[f"floor_bP_{tag}"]),
                        "floor_rel_sd_bhat": float(np.sqrt(max(fb, 0))), "floor_rel_sd_bhat_ci": [float(np.sqrt(max(v, 0))) for v in R.pct_ci(fd[f"floor_b_{tag}"])],
                        "counting_term_mean_bhat_over_C": float(np.mean(fa["b"] / D[m]))})
    res["B4_model_free_floor"] = {"rows": fl_rows, "b_hat": fa["b"], "b_P": bP}
    res["meta"]["seconds"] = round(time.time() - t0, 1)
    R.write_json(OUT, res)
    np.savez_compressed(R.AN / "rev_r2_floor_draws.npz", **{f"{k}__{kk}": vv for k, d in draws.items() for kk, vv in d.items()})
    log("done", verdict, "s", res["meta"]["seconds"])


if __name__ == "__main__":
    main()
