"""Intervals calibrated under crossed state x week dependence (gap-closing round, B-M2).

The pigeonhole bootstrap (state x calendar week, L = 1) over-covers: on the
coverage generators its SD is 1.6-2.5 times the Monte Carlo SD. This script
calibrates two alternatives against the Monte Carlo truth of exactly the
coverage panels (same calibration, panel streams [SEED, 70, stream] and
retained rows as ma_coverage / ma_audit_national; 400 panels per generator)
for S1, S6, S8 and S9:

  1. Two-way cluster-robust sandwich (Cameron, Gelbach and Miller 2011) for the
     stacked estimating equations of the mean xbar and the gamma-FGLS affine fit
     (psi_i = [x_i - xbar, X_i (y_i - mu_i) / mu_i^2], X_i = (1, 1/D_i));
     V = J^-1 (c_S M_S + c_W M_W - M_SW) J^-T with c_G = G / (G - 1), negative
     eigenvalues set to zero. Intervals: a +- t_(G-1) se with G = min(states,
     weeks) (primary), normal; b likewise; b/a on the log scale (delta method).
     One-way state and one-way week sandwiches are reported for comparison.
  2. Variance-ratio-scaled pigeonhole bootstrap: the stored percentile interval
     shrunk towards the point estimate by the bootstrap/MC SD ratio,
     a + (ci - a) / r, with r taken leave-one-generator-out (smallest ratio
     among the other three generators, which widens the interval most).

Decision rule, fixed before the run: the two-way sandwich is the calibrated
interval if its coverage of the estimator mean is at least 0.925 for every
generator; otherwise the leave-one-out scaled bootstrap is used.

Season comparison: season intercepts with a common slope on the pooled loss
(as in ma_extensions), sandwich covariance from the same stacked equations,
Wald test of equal log a_s (adjacent-difference contrasts, 3 df). Its null
distribution is simulated by centring log a_s at the generator's MC mean per
season; the same is done for the bootstrap-covariance Wald statistic used in
the paper (B pigeonhole draws per panel on the first panels; stream
[SEED, 72, stream, panel]). The observed statistics are referred to these
simulated null distributions (calibrated p-values).

Coverage rows carry binomial MC SEs; quantiles carry MC SEs (summ_mc).
Output: analysis/ma_calibration.json, parts in analysis/calibration_parts/.
Run: python code/ma_calibration.py [--workers 8]
"""
from __future__ import annotations

import argparse
import json
import os
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy import stats

import ma_audit as A
import ma_audit_national as N
import ma_coverage as C
from ma_models import Log, binom_mcse, summ_mc as summ, write
from revision_data import load_sample, wl_R
from revision_stats import bootstrap_weights, fit_loss, relative_loss
from s03c_state_fe import fit_group_affine

GENS = {"S1_heterogeneity": 1, "S6_reporting_shock": 7, "S8_combo": 10, "S9_national": 11}
Z = stats.norm.ppf(.975)
ST = {}


def init(log_path):
    N.init(log_path)  # C.STATE holds d, obs, gens (S9 registered), log
    ST["log"] = C.STATE["log"]


def sandwich(x, D, state, week, groups=None):
    """Stacked (xbar, beta) M-estimator: point, two-way / one-way covariance of beta.
    groups=None: beta = (a, b); else beta = (a_1..a_S, b) with season intercepts."""
    n = len(x)
    loss, xbar = relative_loss(x)
    if groups is None:
        f = fit_loss(x, D)
        if not f["converged"]:
            return None
        beta = np.array([f["a"], f["b"]])
        X = np.column_stack([np.ones(n), 1 / D])
    else:
        names, codes = np.unique(groups, return_inverse=True)
        f = fit_group_affine(loss, D, groups)
        if not np.all(f["active"]) or not np.all(np.isfinite(f["intercepts"])):
            return None
        beta = np.r_[np.asarray(f["intercepts"], float), f["b"]]
        X = np.column_stack([np.eye(len(names))[codes], 1 / D])
    mu = X @ beta
    if np.any(mu <= 0):
        return None
    u = (x - xbar) / xbar
    r = (loss - mu) / mu ** 2
    psi = np.column_stack([x - xbar, X * r[:, None]])
    k = psi.shape[1]
    J = np.zeros((k, k))
    J[0, 0] = -n
    J[1:, 0] = X.T @ ((2 * u * (-x / xbar ** 2)) / mu ** 2)
    J[1:, 1:] = -(X.T * (1 / mu ** 2)) @ X
    Ji = np.linalg.inv(J)

    def meat(lab):
        _, ci = np.unique(lab, return_inverse=True)
        G = ci.max() + 1
        S = np.zeros((G, k))
        np.add.at(S, ci, psi)
        return S.T @ S, G

    MS, GS = meat(state)
    MW, GW = meat(week)
    MR = psi.T @ psi
    def cov(M):
        V = Ji @ M @ Ji.T
        w, U = np.linalg.eigh((V + V.T) / 2)
        return (U * np.maximum(w, 0)) @ U.T, int(np.sum(w < 0))
    V2, neg = cov(GS / (GS - 1) * MS + GW / (GW - 1) * MW - MR)
    VS, _ = cov(GS / (GS - 1) * MS)
    VW, _ = cov(GW / (GW - 1) * MW)
    score = float(np.abs(psi[:, 1:].sum(0)).max() / np.abs(psi[:, 1:]).sum(0).max())
    return {"beta": beta, "V2": V2[1:, 1:], "VS": VS[1:, 1:], "VW": VW[1:, 1:], "neg_eig": neg,
            "G": (int(GS), int(GW)), "score_rel": score}


def ab_intervals(s, tcrit):
    """95% intervals for a, b (linear) and b/a (log scale) from each covariance."""
    a, b = s["beta"]
    out = {}
    for lab, V in (("two_way", s["V2"]), ("state", s["VS"]), ("week", s["VW"])):
        sa, sb = np.sqrt(V[0, 0]), np.sqrt(V[1, 1])
        g = np.array([-1 / a, 1 / b]) if a > 0 and b > 0 else None
        slr = float(np.sqrt(g @ V @ g)) if g is not None else np.nan
        lr = np.log(b / a) if g is not None else np.nan
        out[lab] = {"se_a": sa, "se_b": sb, "se_log_b_over_a": slr,
                    "a_ci": [a - tcrit * sa, a + tcrit * sa], "a_ci_z": [a - Z * sa, a + Z * sa],
                    "a_ci_log": [a * np.exp(-tcrit * sa / a), a * np.exp(tcrit * sa / a)] if a > 0 else [np.nan, np.nan],
                    "b_ci": [b - tcrit * sb, b + tcrit * sb],
                    "b_over_a_ci": [float(np.exp(lr - tcrit * slr)), float(np.exp(lr + tcrit * slr))]}
    return out


def season_wald(theta, V, centre=None):
    S = len(theta) - 1  # theta = (a_1..a_S, b)
    a = np.asarray(theta[:S], float)
    if np.any(a <= 0):
        return np.nan
    lt = np.log(a) - (0 if centre is None else centre)
    Gd = np.diag(1 / a)
    Vl = Gd @ V[:S, :S] @ Gd
    Cm = np.eye(S)[:-1] - np.eye(S)[1:]
    try:
        return float((Cm @ lt) @ np.linalg.solve(Cm @ Vl @ Cm.T, Cm @ lt))
    except np.linalg.LinAlgError:
        return np.nan


def boot_season(x, D, frame, season, rng, B):
    """Bootstrap draws of log a_s (common-slope season model), pigeonhole weights."""
    out = np.full((B, len(np.unique(season))), np.nan)
    for j in range(B):
        w = bootstrap_weights(frame, rng)[0]
        if not all((w[season == s] > 0).any() for s in np.unique(season)):
            continue
        loss, _ = relative_loss(x, w)
        f = fit_group_affine(loss, D, season, w)
        a = np.asarray(f["intercepts"], float)
        if np.all(np.isfinite(a)) and np.all(a > 0):
            out[j] = np.log(a)
    return out


def run(job):
    name, lo, hi, reps, Bs, nbs = job["name"], job["lo"], job["hi"], job["reps"], job["boot_season"], job["boot_panels"]
    gen, o, counts = C.panels(name, reps)
    d = gen.d
    ct, cn = d.rows(counts)
    season_all = d.e.season.astype(str).to_numpy()
    rows = []
    for r in range(lo, hi):
        keep = ct[r] >= A.MIN_COUNT
        D, xr = ct[r][keep], wl_R((cn[r][keep] + .5) / (ct[r][keep] + .5))
        st, wk, se = d.state[keep], d.week[keep], season_all[keep]
        G = min(len(np.unique(st)), len(np.unique(wk)))
        row = {"r": r, "a_cov": o["a"][r], "true_total": o["truth_relvar_kept"][r] + o["vH_kept"][r]}
        s = sandwich(xr, D, st, wk)
        if s is not None:
            row.update(a=s["beta"][0], b=s["beta"][1], neg_eig=s["neg_eig"], score_rel=s["score_rel"],
                       iv=ab_intervals(s, stats.t.ppf(.975, G - 1)))
        g = sandwich(xr, D, st, wk, groups=se)
        if g is not None:
            row.update(season_theta=g["beta"].tolist(), season_V2=g["V2"].tolist())
        if r < nbs:
            rng = np.random.default_rng([A.SEED, 72, GENS[name], r])
            bd = boot_season(xr, D, d.e.loc[keep], se, rng, Bs)
            ok = np.all(np.isfinite(bd), 1)
            row.update(season_boot_V=np.cov(bd[ok].T).tolist(), season_boot_valid=int(ok.sum()))
        rows.append(row)
    write(Path(job["outdir"]) / "calibration_parts" / f"{name}_{lo:04d}.json", {"rows": rows})
    ST["log"](name, "panels", lo, hi, "done")
    return name, rows


def cover(ci, t):
    return C.cover(np.asarray(ci, float), t)


def observed(B_note):
    e = load_sample("RA")
    x, D = e.x.to_numpy(float), e.D.to_numpy(float)
    assert np.allclose(D, e.C_t.to_numpy(float)), "D must equal C_t"
    st, wk, se = e.state.to_numpy(), e.week_end.to_numpy(), e.season.astype(str).to_numpy()
    G = min(len(np.unique(st)), len(np.unique(wk)))
    s = sandwich(x, D, st, wk)
    g = sandwich(x, D, st, wk, groups=se)
    return {"G_states_weeks": [len(np.unique(st)), len(np.unique(wk))], "t_df": G - 1,
            "a": float(s["beta"][0]), "b": float(s["beta"][1]), "neg_eig": s["neg_eig"],
            "score_rel": s["score_rel"], "intervals": ab_intervals(s, stats.t.ppf(.975, G - 1)),
            "season_theta": g["beta"].tolist(), "season_V2": g["V2"].tolist(),
            "season_W_sandwich": season_wald(g["beta"], g["V2"]), "seasons": sorted(np.unique(se).tolist())}


def summarise(rows_by, cov_rows, obs, ext, log):
    res, ratios = {}, {}
    for n, rows in rows_by.items():
        ok = [r for r in rows if "iv" in r]
        a = np.array([r["a"] for r in ok])
        b = np.array([r["b"] for r in ok])
        ea, eb = float(a.mean()), float(b.mean())
        tt = np.array([r["true_total"] for r in ok])
        rat = b / np.where(a > 0, a, np.nan)
        sd = {"a": float(np.std(a, ddof=1)), "b": float(np.std(b, ddof=1)),
              "log_b_over_a": float(np.nanstd(np.log(rat), ddof=1))}
        g = {"panels": len(rows), "panels_used": len(ok), "estimator_mean": {"a": ea, "b": eb, "b_over_a_ratio_of_means": eb / ea},
             "mc_sd": sd, "neg_eig_share": float(np.mean([r["neg_eig"] > 0 for r in ok])),
             "max_score_rel": float(max(r["score_rel"] for r in ok))}
        for lab in ("two_way", "state", "week"):
            iv = [r["iv"][lab] for r in ok]
            g[lab] = {"a_vs_estimator_mean": cover([v["a_ci"] for v in iv], ea),
                      "a_vs_estimator_mean_normal": cover([v["a_ci_z"] for v in iv], ea),
                      "a_log_vs_estimator_mean": cover([v["a_ci_log"] for v in iv], ea),
                      "a_vs_true_total": cover([v["a_ci"] for v in iv], tt),
                      "b_vs_estimator_mean": cover([v["b_ci"] for v in iv], eb),
                      "b_over_a_vs_ratio_of_means": cover([v["b_over_a_ci"] for v in iv], eb / ea),
                      "se_ratio_over_mc": {"a": float(np.mean([v["se_a"] for v in iv]) / sd["a"]),
                                           "b": float(np.mean([v["se_b"] for v in iv]) / sd["b"]),
                                           "log_b_over_a": float(np.nanmean([v["se_log_b_over_a"] for v in iv]) / sd["log_b_over_a"])}}
        # pigeonhole bootstrap (stored coverage parts): SD ratio on the panels' own estimates
        cr = cov_rows[n]
        acov = np.array([c["a"] for c in cr])
        ratios[n] = {"a": float(np.nanmean([c["a_bsd"] for c in cr]) / np.std(acov, ddof=1)),
                     "b": float(np.nanmean([c["b_bsd"] for c in cr]) / np.std([c["b"] for c in cr], ddof=1))}
        g["bootstrap_sd_ratio"] = ratios[n]
        # season sandwich Wald: raw and centred (null) distribution
        th = np.array([r["season_theta"] for r in rows if "season_theta" in r])
        V = [np.array(r["season_V2"]) for r in rows if "season_theta" in r]
        S = th.shape[1] - 1
        la = np.log(np.where(th[:, :S] > 0, th[:, :S], np.nan))
        centre = np.nanmean(la, 0)
        Wc = np.array([season_wald(t, v, centre) for t, v in zip(th, V)])
        Wr = np.array([season_wald(t, v) for t, v in zip(th, V)])
        crit = stats.chi2.ppf(.95, S - 1)
        sz = float(np.nanmean(Wc > crit))
        g["season_sandwich"] = {"panels": int(np.isfinite(Wc).sum()), "size_at_chi2_crit": sz,
                                "size_mcse": binom_mcse(sz, int(np.isfinite(Wc).sum())),
                                "raw_rejection": float(np.nanmean(Wr > crit)), "null_W": summ(Wc),
                                "null_q95": float(np.nanpercentile(Wc, 95)),
                                "season_log_a_mc_mean": centre.tolist(),
                                "season_log_a_mc_sd": np.nanstd(la, 0, ddof=1).tolist(),
                                "p_observed": float(np.nanmean(Wc >= obs["season_W_sandwich"]))}
        g["season_sandwich"]["p_observed_mcse"] = binom_mcse(g["season_sandwich"]["p_observed"], int(np.isfinite(Wc).sum()))
        bth = [(r["season_theta"], np.array(r["season_boot_V"])) for r in rows if "season_boot_V" in r and "season_theta" in r]
        if bth:
            Wb = []
            for t, vb in bth:
                lt = np.log(np.asarray(t[:S], float)) - centre
                Cm = np.eye(S)[:-1] - np.eye(S)[1:]
                Wb.append(float((Cm @ lt) @ np.linalg.solve(Cm @ vb @ Cm.T, Cm @ lt)) if np.all(np.isfinite(lt)) else np.nan)
            Wb = np.array(Wb)
            m = int(np.isfinite(Wb).sum())
            szb = float(np.nanmean(Wb > crit))
            pb = float(np.nanmean(Wb >= ext["wald_log_a"]["stat"]))
            g["season_bootstrap"] = {"panels": m, "boot_draws": int(np.median([r["season_boot_valid"] for r in rows if "season_boot_valid" in r])),
                                     "size_at_chi2_crit": szb, "size_mcse": binom_mcse(szb, m), "null_W": summ(Wb),
                                     "null_q95": float(np.nanpercentile(Wb, 95)), "p_observed": pb, "p_observed_mcse": binom_mcse(pb, m)}
        res[n] = g
        log(n, "two-way a/E", g["two_way"]["a_vs_estimator_mean"]["coverage"], "a/true", g["two_way"]["a_vs_true_total"]["coverage"],
            "b/E", g["two_way"]["b_vs_estimator_mean"]["coverage"], "b/a", g["two_way"]["b_over_a_vs_ratio_of_means"]["coverage"],
            "se ratio", g["two_way"]["se_ratio_over_mc"], "boot ratio", ratios[n],
            "season size", sz, "p_obs", g["season_sandwich"]["p_observed"])
    # leave-one-generator-out scaled bootstrap
    for n in rows_by:
        others = [m for m in rows_by if m != n]
        cr = cov_rows[n]
        a = np.array([c["a"] for c in cr])
        tt = np.array([c["true_total"] for c in cr])
        ea = float(np.mean(a))
        out = {}
        for k in ("a", "b"):
            rr = min(ratios[m][k] for m in others)
            pt = np.array([c[k] for c in cr])
            ci = np.array([c[f"{k}_ci"] for c in cr])
            sc = pt[:, None] + (ci - pt[:, None]) / rr
            tgt = ea if k == "a" else float(np.mean(pt))
            out[k] = {"ratio_used": rr, "vs_estimator_mean": cover(sc, tgt)}
            if k == "a":
                out[k]["vs_true_total"] = cover(sc, tt)
                out[k]["unscaled_vs_estimator_mean"] = cover(ci, tgt)
        res[n]["scaled_bootstrap_loo"] = out
    rule_ok = all(res[n]["two_way"]["a_vs_estimator_mean"]["coverage"] >= .925 for n in res)
    chosen = "two_way_sandwich" if rule_ok else "scaled_bootstrap_loo"
    # applied to the observed influenza estimates
    prim = ext["primary"]
    rmin = min(ratios[n]["a"] for n in ratios)
    rmax = max(ratios[n]["a"] for n in ratios)
    ap, (alo, ahi) = prim["a"]["value"], prim["a"]["ci95"]
    bp_, (blo, bhi) = prim["b"]["value"], prim["b"]["ci95"]
    rbmin = min(ratios[n]["b"] for n in ratios)
    applied = {"two_way_sandwich": obs["intervals"]["two_way"],
               "scaled_bootstrap": {"a_ratio_min": rmin, "a_ratio_max": rmax,
                                    "a_ci_ratio_min": [ap + (alo - ap) / rmin, ap + (ahi - ap) / rmin],
                                    "a_ci_ratio_max": [ap + (alo - ap) / rmax, ap + (ahi - ap) / rmax],
                                    "b_ratio_min": rbmin, "b_ci_ratio_min": [bp_ + (blo - bp_) / rbmin, bp_ + (bhi - bp_) / rbmin]},
               "primary_bootstrap": {"a": [ap, alo, ahi], "b": [bp_, blo, bhi]}}
    ps = [res[n]["season_sandwich"]["p_observed"] for n in res]
    pb = [res[n]["season_bootstrap"]["p_observed"] for n in res if "season_bootstrap" in res[n]]
    season = {"observed_W_sandwich": obs["season_W_sandwich"], "observed_W_bootstrap": ext["wald_log_a"]["stat"],
              "chi2_p_sandwich": float(stats.chi2.sf(obs["season_W_sandwich"], len(obs["seasons"]) - 1)),
              "chi2_p_bootstrap": ext["wald_log_a"]["p"],
              "calibrated_p_sandwich_range": [min(ps), max(ps)], "calibrated_p_bootstrap_range": [min(pb), max(pb)] if pb else None}
    return {"generators": res, "decision_rule": "two-way sandwich if its a coverage of the estimator mean >= 0.925 for all generators, else LOO scaled bootstrap",
            "chosen": chosen, "observed": obs, "applied": applied, "season_test": season}


def main(reps, chunk, workers, boot_season, boot_panels, outdir, logdir, log_name):
    log_path = Path(logdir) / log_name
    log = Log(log_path)
    obs = observed(None)
    log("observed two-way", obs["intervals"]["two_way"], "season W sandwich", obs["season_W_sandwich"])
    ext = json.loads((Path(outdir) / "ma_extensions.json").read_text())
    ext = {"primary": ext["ratios"]["RA_gi3.6"], "wald_log_a": ext["seasons"]["test"]["wald_log_a"]}
    upto = int(os.environ.get("SMOKE_PANELS", reps))  # smoke: full panel set generated, first panels analysed
    jobs = [{"name": n, "lo": lo, "hi": min(lo + chunk, upto), "reps": reps, "boot_season": boot_season,
             "boot_panels": boot_panels, "outdir": str(outdir)} for n in GENS for lo in range(0, upto, chunk)]
    log("calibration jobs", len(jobs), "panels", reps, "season bootstrap", boot_season, "x", boot_panels, "workers", workers)
    with Pool(workers, initializer=init, initargs=(log_path,)) as pool:
        out = pool.map(run, jobs, chunksize=1)
    rows_by = {n: sorted((r for nm, rr in out if nm == n for r in rr), key=lambda r: r["r"]) for n in GENS}
    cov_rows = {}
    for n in GENS:
        parts = sorted((Path(outdir) / "coverage_parts").glob(f"{n}_*.json"))
        cr = sorted((r for p in parts for r in json.loads(p.read_text())["rows"]), key=lambda r: r["r"])
        if not os.environ.get("SMOKE"):
            assert len(cr) == reps, (n, len(cr))
        cr = cr[:reps]
        # self-check: regenerated panels reproduce the stored coverage estimates exactly
        mine = {r["r"]: r["a_cov"] for r in rows_by[n]}
        dev = max(abs(c["a"] - mine[c["r"]]) for c in cr if c["r"] in mine)
        assert dev < 1e-12, (n, "panel mismatch", dev)
        fa = max(abs(r["a"] - r["a_cov"]) for r in rows_by[n] if "a" in r)
        assert fa < 1e-8, (n, "sandwich fit differs from coverage fit", fa)
        cov_rows[n] = [c for c in cr if c["r"] in mine]
    res = summarise(rows_by, cov_rows, obs, ext, log)
    res["meta"] = {"seed": A.SEED, "panels": reps, "panel_rng": "[SEED, 70, stream] (as ma_coverage)",
                   "season_bootstrap": {"draws": boot_season, "panels": boot_panels, "rng": "[SEED, 72, stream, panel]"},
                   "generators": GENS, "t_crit_df": "min(states, weeks) - 1 per panel"}
    write(Path(outdir) / "ma_calibration.json", res)
    log("calibration complete; chosen", res["chosen"])


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--reps", type=int, default=400)
    p.add_argument("--chunk", type=int, default=50)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--boot-season", type=int, default=200)
    p.add_argument("--boot-panels", type=int, default=200)
    p.add_argument("--output-dir", type=Path, default=A.EXP / "analysis")
    p.add_argument("--log-dir", type=Path, default=A.EXP / "logs")
    p.add_argument("--log-name", default="gap_calibration.log")
    a = p.parse_args()
    if os.environ.get("SMOKE"):
        a.chunk, a.boot_season, a.boot_panels, a.log_name = 4, 20, 4, "gap_calibration_smoke.log"
        os.environ["SMOKE_PANELS"] = "8"
        a.output_dir = a.output_dir / "smoke"
        a.output_dir.mkdir(parents=True, exist_ok=True)
        src = A.EXP / "analysis"
        (a.output_dir / "ma_extensions.json").write_text((src / "ma_extensions.json").read_text())
        (a.output_dir / "coverage_parts").mkdir(exist_ok=True)
        for n in GENS:
            for f in sorted((src / "coverage_parts").glob(f"{n}_*.json")):
                (a.output_dir / "coverage_parts" / f.name).write_text(f.read_text())
    main(a.reps, a.chunk, a.workers, a.boot_season, a.boot_panels, a.output_dir, a.log_dir, a.log_name)
