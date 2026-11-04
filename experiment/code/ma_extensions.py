"""P1 empirical extensions for mainline A, on the authors' estimators and bootstrap.

Sections (each written to analysis/ma_extensions.json as soon as it finishes):
  ratios        b/a, b_P and the unthinned b(1)/b_P for R_A under GI 3.6/1.6
                (Cowling) and GI 3.2 (sd scaled, and sd 2.1 as in Chan), paired
                state x week bootstrap (revision_stats.bootstrap_weights)
  within        within slopes on the pooled relative loss: state FE (authors'
                fit_group_affine), week FE and state + week two-way (one-step
                precision-weighted FGLS, ma_models.within_slope) with bootstrap;
                iterated two-way variant reported with its convergence
  sensitivity   b/a and b/b_P under count filter >= 10/20/50, continuity +1.0,
                window threshold 1.10, and per season
  seasons       season-stratified affine model: season intercepts + common
                slope on the pooled loss; fully stratified (season-centred) fits
  validation    held-out affine vs power score by state and by season folds,
                paired fold standard error and state x week bootstrap interval
  table1        observed-denominator (ratio estimator) companion row of Table 1:
                delta-method MSE, Monte Carlo MSE and its m-cross
  calibration   WITHDRAWN 2026-11-03 (kept as a function, not in the default run):
                inversion v and slope / n / autocorrelation diagnostics of the
                authors' calibration (reproduced outputs, not rerun)

Primary bootstrap (revision 2, spec A3): the RA_gi3.6 frame of `ratios`, B=2000,
rng default_rng([SEED, 1]); the sensitivity primary row and the within-table
pooled row reuse exactly this stream on the same rows, so their pooled entries
coincide with the primary interval draw by draw.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from ma_models import Log, boot_ci, codes_of, summ, utc, within_slope, write
from revision_data import build_panel, load_raw, load_sample, national_series, rule_windows, wl_R
from revision_empirical import centered_loss, fit_power
from revision_stats import bootstrap_weights, fit_affine, fit_loss, relative_loss
from scipy import stats

from s03b_thinning import poisson_reference
from s03c_state_fe import fit_group_affine

SEED = 2026110352
EXP = Path(__file__).resolve().parents[1]
SMOKE = bool(os.environ.get("SMOKE"))


BRANCH_RHOS = (0.1, 0.25, 0.5, 1.0)


def branching_reference(e, w=None, gi=(3.6, None), rho=.25):
    """Plug-in reference for a Poisson branching chain with binomial reporting rho:
    weighted mean kappa^2 [(1 - rho) + (C_t+.5)/(C_next+.5)] (independent Poisson: 1 + 1/g)."""
    mean = gi[0]
    sd = mean * 1.6 / 3.6 if gi[1] is None else gi[1]
    ct, cn = e.C_t.to_numpy(float), e.C_next.to_numpy(float)
    weights = np.ones(len(e)) if w is None else np.asarray(w, float)
    use = weights > 0
    g = (cn[use] + .5) / (ct[use] + .5)
    kappa = (mean / 7) / (1 + (sd ** 2 / mean) * np.log(g) / 7)
    return float(np.average(kappa ** 2 * ((1 - rho) + 1 / g), weights=weights[use]))


def ab_bp(e, w=None, gi=(3.6, None)):
    x, D = e.x.to_numpy(float), e.D.to_numpy(float)
    f = fit_loss(x, D, w)
    if not f["converged"]:
        return [np.nan] * len(AB_NAMES)
    bp = poisson_reference(e, w, mean=gi[0], sd=gi[1])
    ratio = f["b"] / f["a"] if f["a"] > 0 else np.nan
    br = [f["b"] / branching_reference(e, w, gi, rho) for rho in BRANCH_RHOS]
    return [f["a"], f["b"], ratio, bp, f["b"] / bp, *br]


AB_NAMES = ("a", "b", "b_over_a", "bP", "b_over_bP", *(f"b_over_bPbr_rho{r:g}" for r in BRANCH_RHOS))


def boot_ab(e, B, rng, gi=(3.6, None)):
    point = ab_bp(e, None, gi)
    draws = np.array([ab_bp(e, bootstrap_weights(e, rng)[0], gi) for _ in range(B)])
    return {"n": len(e), "states": int(e.state.nunique()), "weeks": int(e.week_end.nunique()),
            **{nm: boot_ci(point[j], draws[:, j], ratio=nm == "b_over_a") for j, nm in enumerate(AB_NAMES)}}, draws


def section_ratios(B, log):
    frames = {"RA_gi3.6": (load_sample("RA"), (3.6, None)),
              "RA_gi3.2": (load_sample("RA", 3.2), (3.2, None)),
              "RA_gi3.2_sd2.1": (load_sample("RA", 3.2, gi_sd=2.1), (3.2, 2.1)),
              "G_gi3.6": (load_sample("G"), (3.6, None))}
    out, store = {}, {}
    for name, (e, gi) in frames.items():
        # same seed for every GI frame: the frames share rows, so draw i reuses the
        # same state/week weights and the GI contrasts below are paired
        rng = np.random.default_rng([SEED, 1])
        out[name], store[name] = boot_ab(e, B, rng, gi)
        out[name]["gi"] = {"mean": gi[0], "sd": gi[1] if gi[1] is not None else gi[0] * 1.6 / 3.6}
        out[name]["bPbr"] = {f"rho{r:g}": branching_reference(e, None, gi, r) for r in BRANCH_RHOS}
        log("ratios", name, {k: out[name][k]["value"] for k in AB_NAMES}, "b/bP ci", out[name]["b_over_bP"]["ci95"])
    base = store["RA_gi3.6"]
    same_rows = all(frames[n][0][["state", "week_end"]].equals(frames["RA_gi3.6"][0][["state", "week_end"]]) for n in frames if n.startswith("RA"))
    out["paired_weights"] = bool(same_rows)
    out["primary"] = {"frame": "RA_gi3.6", "B": B, "rng": [SEED, 1], "G_frame": "G_gi3.6 (same rng, own rows; supplement only)"}
    for name in ("RA_gi3.2", "RA_gi3.2_sd2.1"):
        for j, metric in ((2, "b_over_a"), (4, "b_over_bP")):
            pt = out[name][metric]["value"] - out["RA_gi3.6"][metric]["value"]
            out[f"{name}_minus_gi3.6_{metric}"] = boot_ci(pt, store[name][:, j] - base[:, j])
    return out


def iterated_b(loss, D, groups, w=None):
    """Damped iterated FGLS within slope; NaN when it does not converge."""
    r = within_slope(loss, D, groups, w, iterate=500, damp=.5)
    return r["b"] if r["converged"] else np.nan


def two_way_point_and_boot(e, B, log, rng):
    D = e.D.to_numpy(float)
    s, t = codes_of(e.state), codes_of(e.week_end)
    loss, _ = relative_loss(e.x.to_numpy(float))
    point = {}
    point["state_authors"] = fit_group_affine(loss, D, e.state)["b"]
    point["state_onestep"] = within_slope(loss, D, (s,))
    point["week_onestep"] = within_slope(loss, D, (t,))
    point["two_way_onestep"] = within_slope(loss, D, (s, t))
    point["two_way_iterated_undamped"] = within_slope(loss, D, (s, t), iterate=50)
    point["two_way_iterated"] = within_slope(loss, D, (s, t), iterate=500, damp=.5)
    bp = poisson_reference(e)
    names = ("pooled_b", "bP", "state_authors", "state_onestep", "week_onestep", "two_way_onestep", "two_way_iterated", "pooled_a")
    draws = np.full((B, len(names)), np.nan)
    for i in range(B):
        w = bootstrap_weights(e, rng)[0]
        try:
            lw, _ = relative_loss(e.x.to_numpy(float), w)
        except ValueError:
            continue
        f = fit_loss(e.x.to_numpy(float), D, w)
        fs = fit_group_affine(lw, D, e.state, w)
        draws[i] = [f["b"] if f["converged"] else np.nan, poisson_reference(e, w),
                    fs["b"] if fs["converged"] else np.nan,
                    within_slope(lw, D, (s,), w)["b"], within_slope(lw, D, (t,), w)["b"],
                    within_slope(lw, D, (s, t), w)["b"], iterated_b(lw, D, (s, t), w),
                    f["a"] if f["converged"] else np.nan]
        if (i + 1) % 200 == 0:
            log("within bootstrap", i + 1)
    pooled = fit_loss(e.x.to_numpy(float), D)
    pv = [pooled["b"], bp, point["state_authors"], point["state_onestep"]["b"], point["week_onestep"]["b"],
          point["two_way_onestep"]["b"], point["two_way_iterated"]["b"], pooled["a"]]
    slopes = names[:-1]
    res = {"point_fits": point, "slopes": {nm: boot_ci(pv[j], draws[:, j]) for j, nm in enumerate(slopes)}}
    res["over_bP"] = {nm: boot_ci(pv[j] / bp, draws[:, j] / draws[:, 1]) for j, nm in enumerate(slopes) if nm != "bP"}
    res["over_pooled_b"] = {nm: boot_ci(pv[j] / pv[0], draws[:, j] / draws[:, 0]) for j, nm in enumerate(slopes) if nm not in ("bP", "pooled_b")}
    pa = draws[:, -1]
    res["pooled"] = {"a": boot_ci(pooled["a"], pa), "b": boot_ci(pooled["b"], draws[:, 0]),
                     "b_over_a": boot_ci(pooled["b"] / pooled["a"], np.where(pa > 0, draws[:, 0] / pa, np.nan), ratio=True),
                     "bP": boot_ci(bp, draws[:, 1]), "b_over_bP": boot_ci(pooled["b"] / bp, draws[:, 0] / draws[:, 1])}
    res["bootstrap"] = {"B": B, "rng": "primary stream default_rng([SEED, 1])"}
    return res


def section_within(B, log):
    e = load_sample("RA")
    res = two_way_point_and_boot(e, B, log, np.random.default_rng([SEED, 1]))
    log("within", {k: v["value"] for k, v in res["slopes"].items()}, "iterated", res["point_fits"]["two_way_iterated"]["b"],
        res["point_fits"]["two_way_iterated"]["converged"])
    return res


def sens_frames():
    panel = build_panel(load_raw())
    nat = national_series()
    w105, w110 = rule_windows(nat, 1.05)[0], rule_windows(nat, 1.10)[0]
    base = load_sample("RA", panel=panel, windows=w105)
    specs = {"primary": base}
    for m in (10, 20, 50):
        specs[f"min_count_{m}"] = base[(base.C_t >= m) & (base.D_t >= m)].reset_index(drop=True)
    cc = base.copy()
    cc["x"] = wl_R((cc.C_next + 1.0) / (cc.C_t + 1.0))
    specs["continuity_1.0"] = cc
    specs["window_threshold_1.10"] = load_sample("RA", panel=panel, windows=w110)
    for season in sorted(base.season.unique()):
        specs[f"season_{season}"] = base[base.season == season].reset_index(drop=True)
    windows = {"1.05": {s: [a.date().isoformat(), b.date().isoformat()] for s, (a, b) in w105.items()},
               "1.10": {s: [a.date().isoformat(), b.date().isoformat()] for s, (a, b) in w110.items()}}
    return specs, windows


def section_sensitivity(B, log):
    specs, windows = sens_frames()
    out = {"windows": windows, "specs": {}}
    for j, (name, e) in enumerate(specs.items()):
        # the primary specification reuses the primary stream (same rows as ratios/RA_gi3.6)
        r, _ = boot_ab(e, B, np.random.default_rng([SEED, 1] if name == "primary" else [SEED, 3, j]))
        out["specs"][name] = r
        log("sensitivity", name, r["n"], "b/a", r["b_over_a"]["value"], r["b_over_a"]["ci95"], "b/bP", r["b_over_bP"]["value"])
    return out


def season_models(e, w=None):
    """Season intercepts + common slope (pooled loss) and season-centred common slope."""
    x, D = e.x.to_numpy(float), e.D.to_numpy(float)
    loss, _ = relative_loss(x, w)
    f = fit_group_affine(loss, D, e.season, w)
    z = centered_loss(e, "season", w)
    Lc, valid, _ = z
    wc = np.ones(len(e)) if w is None else np.asarray(w, float)
    fc = fit_group_affine(np.where(valid, Lc, 0.0), D, e.season, np.where(valid, wc, 0.0))
    return f, fc


def homogeneity_tests(theta, draws):
    """Wald test of equal season parameters (bootstrap covariance, adjacent-difference
    contrasts; invariant to the contrast basis) and Cochran's Q with inverse-variance
    weights. Draws with any non-finite component are dropped."""
    S = len(theta)
    ok = np.all(np.isfinite(draws), 1)
    d = draws[ok]
    V = np.cov(d.T)
    C = np.eye(S)[:-1] - np.eye(S)[1:]
    W = float((C @ theta) @ np.linalg.solve(C @ V @ C.T, C @ theta))
    wts = 1 / np.diag(V)
    tbar = float(wts @ theta / wts.sum())
    Q = float(wts @ (theta - tbar) ** 2)
    return {"theta": np.asarray(theta, float).tolist(), "se": np.sqrt(np.diag(V)).tolist(), "weighted_mean": tbar,
            "n_draws_used": int(ok.sum()), "n_draws_total": int(len(draws)),
            "wald": {"stat": W, "df": S - 1, "p": float(stats.chi2.sf(W, S - 1))},
            "cochran_q": {"Q": Q, "df": S - 1, "p": float(stats.chi2.sf(Q, S - 1)), "I2": max(0.0, (Q - (S - 1)) / Q) if Q > 0 else 0.0}}


def noncentrality_for_power(df, alpha=.05, power=.8):
    """Noncentrality at which a chi-square(df) test of level alpha has the given power."""
    crit = stats.chi2.ppf(1 - alpha, df)
    lo, hi = 1e-6, 100.0
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if stats.ncx2.sf(crit, df, mid) < power else (lo, mid)
    return (lo + hi) / 2


def section_seasons(B, log):
    e = load_sample("RA")
    f, fc = season_models(e)
    seasons = [str(s) for s in f["groups"]]
    S = len(seasons)
    lab = e.season.astype(str).to_numpy()
    x, D = e.x.to_numpy(float), e.D.to_numpy(float)
    def vec(fit):
        if not fit["converged"]:
            return [np.nan] * (2 * len(seasons) + 1)
        a = np.asarray(fit["intercepts"], float)
        return [fit["b"], *a, *(fit["b"] / np.where(a > 0, a, np.nan))]
    def separate(w=None):
        res = []
        for sname in seasons:
            m = lab == sname
            ww = None if w is None else w[m]
            fs = fit_loss(x[m], D[m], ww) if ww is None or ww.sum() > 0 else {"converged": False}
            res.append((fs["a"], fs["b"]) if fs.get("converged") else (np.nan, np.nan))
        return np.array(res)
    p1, p2 = vec(f), vec(fc)
    d1, d2 = np.full((B, len(p1)), np.nan), np.full((B, len(p2)), np.nan)
    rng = np.random.default_rng([SEED, 4])
    for i in range(B):
        w = bootstrap_weights(e, rng)[0]
        g1, g2 = season_models(e, w)
        if not all((w[e.season.to_numpy() == s] > 0).any() for s in seasons):
            continue
        d1[i], d2[i] = vec(g1), vec(g2)
    # separate season fits (sensitivity) on their own stream, as in revision2/checks/c02
    sep0 = separate()
    SEP = np.full((B, S, 2), np.nan)
    rng2 = np.random.default_rng([SEED, 40])
    for i in range(B):
        w2 = bootstrap_weights(e, rng2)[0]
        if all((w2[lab == s] > 0).any() for s in seasons):
            SEP[i] = separate(w2)
    labels = ["b_common"] + [f"a_{s}" for s in seasons] + [f"b_over_a_{s}" for s in seasons]
    out = {"seasons": seasons,
           "pooled_loss_season_intercepts": {lab_: boot_ci(p1[j], d1[:, j], ratio=lab_.startswith("b_over")) for j, lab_ in enumerate(labels)},
           "season_centred_common_slope": {lab_: boot_ci(p2[j], d2[:, j], ratio=lab_.startswith("b_over")) for j, lab_ in enumerate(labels)},
           "fit_methods": [f.get("method"), fc.get("method")]}
    a_pt, a_dr = np.asarray(p1[1:1 + S], float), d1[:, 1:1 + S]
    # descriptive range only; not a homogeneity test
    ratios = d1[:, 1 + S:]
    out["descriptive"] = {"max_over_min_b_over_a": {**boot_ci(
        np.nanmax(p1[1 + S:]) / np.nanmin(p1[1 + S:]), np.nanmax(ratios, 1) / np.nanmin(ratios, 1)),
        "note": "descriptive; not a homogeneity test"}}
    t_log = homogeneity_tests(np.log(a_pt), np.log(np.where(a_dr > 0, a_dr, np.nan)))
    t_raw = homogeneity_tests(a_pt, a_dr)
    r0 = np.log(sep0[:, 1] / sep0[:, 0])
    R = np.log(np.where((SEP[..., 0] > 0) & (SEP[..., 1] > 0), SEP[..., 1] / SEP[..., 0], np.nan))
    t_sep = homogeneity_tests(r0, R)
    lam_nc = noncentrality_for_power(S - 1)
    W = t_log["wald"]["stat"]
    fac = float(np.sqrt(lam_nc / W))
    dev = np.log(a_pt) - t_log["weighted_mean"]
    out["test"] = {
        "wald_log_a": {**t_log["wald"], "param": "log a_s, season intercepts with common slope (equivalently log b/a_s)"},
        "cochran_q_log_a": t_log["cochran_q"],
        "wald_raw_a": {**t_raw["wald"], "param": "a_s, common-slope model (sensitivity)", "n_draws_used": t_raw["n_draws_used"]},
        "wald_sep_log_ba": {**t_sep["wald"], "param": "log b_s/a_s, separate season fits, rng [SEED,40] (sensitivity)",
                            "n_draws_used": t_sep["n_draws_used"]},
        "noncentrality_power80": lam_nc,
        "mdd_factor": fac,
        "mdd_max_over_min_a": float(np.exp(fac * (dev.max() - dev.min()))),
        "observed_max_over_min_a": float(np.exp(dev.max() - dev.min())),
        "n_draws_used": t_log["n_draws_used"], "n_draws_total": B,
        "detail": {"log_a": t_log, "raw_a": t_raw, "separate_log_b_over_a": t_sep},
        "separate_point": {"a": sep0[:, 0].tolist(), "b": sep0[:, 1].tolist()},
        "contrast": "3 adjacent season differences; the Wald statistic is invariant to the contrast basis",
        "mdd_definition": ("f = sqrt(lambda/W), lambda the noncentrality giving power 0.8 for a level-0.05 chi-square(3) test; "
                           "mdd_max_over_min_a = exp(f * range of log a_s deviations from their inverse-variance weighted mean)")}
    log("seasons common b", p1[0], "season b/a", p1[1 + S:], "Wald", t_log["wald"], "Q", t_log["cochran_q"],
        "raw", t_raw["wald"], "sep", t_sep["wald"], "lambda", lam_nc, "f", fac)
    return out


def fold_scores(e, group, w):
    x, D = e.x.to_numpy(float), e.D.to_numpy(float)
    labels = e[group].astype(str).to_numpy()
    rows = []
    for lab in np.unique(labels):
        test = (labels == lab) & (w > 0)
        train = (labels != lab) & (w > 0)
        if not test.any() or not train.any():
            continue
        loss, center = relative_loss(x[train], w[train])
        af = fit_affine(loss, D[train], w[train])
        pw = fit_power(loss, D[train], w[train])
        held = ((x[test] - center) / center) ** 2
        ma = af["a"] + af["b"] / D[test]
        mp = pw["coefficient"] * D[test] ** (-pw["exponent"])
        if not (af["converged"] and pw["converged"] and np.all(ma > 0) and np.all(mp > 0)):
            continue
        wt = w[test]
        sa = float(np.average(held / ma + np.log(ma), weights=wt))
        sp = float(np.average(held / mp + np.log(mp), weights=wt))
        rows.append((lab, float(wt.sum()), sp - sa))
    return rows


def paired_summary(rows):
    n = np.array([r[1] for r in rows])
    d = np.array([r[2] for r in rows])
    F = len(rows)
    mean = float(n @ d / n.sum())
    se = float(np.sqrt(F / (F - 1) * np.sum((n / n.sum()) ** 2 * (d - mean) ** 2))) if F > 1 else np.nan
    return {"folds": F, "mean_power_minus_affine": mean, "paired_fold_se": se,
            "z": mean / se if se and se > 0 else None, "affine_better_folds": int(np.sum(d > 0)),
            "fold_values": [{"fold": r[0], "weight": r[1], "power_minus_affine": r[2]} for r in rows]}


def section_validation(B, log):
    e = load_sample("RA")
    out = {}
    for j, group in enumerate(("state", "season")):
        point = paired_summary(fold_scores(e, group, np.ones(len(e))))
        rng = np.random.default_rng([SEED, 5, j])
        draws = np.full(B, np.nan)
        for i in range(B):
            rows = fold_scores(e, group, bootstrap_weights(e, rng)[0])
            if len(rows) > 1:
                draws[i] = paired_summary(rows)["mean_power_minus_affine"]
        point["bootstrap"] = boot_ci(point["mean_power_minus_affine"], draws)
        point["bootstrap_share_affine_better"] = float(np.nanmean(draws > 0))
        out[group] = point
        log("validation", group, point["mean_power_minus_affine"], "se", point["paired_fold_se"], point["bootstrap"]["ci95"])
    out["score"] = "mean L/mu + log mu on held-out rows; positive power-minus-affine favours affine"
    return out


def section_table1(reps, log):
    rb, vr, k, rho = 1.15, .04, .35, .25           # revision_theory.py parameters
    A = rb + rho * (rb ** 2 + vr) / k
    extra = (1 - rho) * (rb ** 2 + vr)
    rows = []
    for j, m in enumerate((10., 50., 200., 1000.)):
        rng = np.random.default_rng([SEED, 6, j])
        n = int(round(m / rho))
        R = rng.gamma(rb ** 2 / vr, vr / rb, size=reps)
        lam = rng.gamma(n * k, R / k)
        C = rng.poisson(rho * lam)
        cprev = rng.binomial(n, rho, size=reps)
        ok = cprev > 0
        err = (C[ok] / cprev[ok] - R[ok]) ** 2
        known = (C / m - R) ** 2
        rows.append({"m": m, "known_denominator_analytic": A / m, "known_denominator_mc": float(known.mean()),
                     "known_denominator_mc_se": float(known.std(ddof=1) / np.sqrt(reps)),
                     "observed_denominator_delta": (A + extra) / m, "observed_denominator_mc": float(err.mean()),
                     "observed_denominator_mc_se": float(err.std(ddof=1) / np.sqrt(ok.sum())),
                     "observed_over_known": float(err.mean() / (A / m)), "share_zero_denominator": float(1 - ok.mean())})
        log("table1 m", m, rows[-1]["observed_denominator_mc"], "delta", rows[-1]["observed_denominator_delta"])
    return {"parameters": {"Rbar": rb, "v_R": vr, "k": k, "rho": rho, "A": A, "reps": reps},
            "m_cross_known": A / vr, "m_cross_observed_delta": (A + extra) / vr,
            "estimator": "C_t / C_prev with C_prev ~ Binomial(m/rho, rho) observed previous count (exposure not known)",
            "rows": rows}


def section_calibration(log):
    src = EXP / "analysis/authors/calibration_reproduced"
    cal = json.loads((src / "calibration.json").read_text())
    emp = json.loads((src / "calibration_empirical.json").read_text())
    obs = emp["3.6"]["RA"]
    out = {"source": "analysis/authors/calibration_reproduced/calibration.json (reproduced, bit-identical to authors)",
           "observed": {k: obs[k] for k in ("n", "a", "b", "ratio", "r1_state", "r1_national", "VR_week", "D_median")},
           "scenarios": {}}
    cells = {}
    for c in cal["cells"]:
        cells.setdefault(c["scenario"], []).append(c)
    for name, inf in cal["inference"].items():
        r = inf["RA"]
        vd = r["v_draws"]
        cs = sorted(cells.get(name, []), key=lambda c: c["v"])
        v = r["v_implied"]
        near = min(cs, key=lambda c: abs(c["v"] - v)) if cs else None
        diag = None
        if near is not None:
            ra = near["RA"]
            diag = {"cell_v": near["v"], **{k: ra[k]["mean_finite"] for k in ("b", "n", "r1_state", "r1_national", "VR_week")}}
        out["scenarios"][name] = {"v_implied": v, "v_ci95": [vd["q025_50_975_finite"][0], vd["q025_50_975_finite"][2]],
                                  "line_status": r["line_status"], "nearest_cell": diag}
    m = out["scenarios"]["matched36_binomial"]["nearest_cell"]
    out["primary_mismatch"] = {k: {"simulated": m[k], "observed": out["observed"][k]} for k in ("b", "n", "r1_national", "r1_state", "VR_week")}
    vs = [s["v_implied"] for s in out["scenarios"].values()]
    out["v_range_across_scenarios"] = [float(np.min(vs)), float(np.max(vs))]
    out["decision"] = ("withdraw from the manuscript: the inversion exists but the calibrated simulator does not reproduce "
                       "the pooled slope, the retained sample size or the national lag-1 autocorrelation at the implied v, "
                       "so the map from a to v is not supported by its own checks")
    log("calibration primary mismatch", out["primary_mismatch"], "v range", out["v_range_across_scenarios"])
    return out


def main(sections, outdir, logdir):
    log = Log(Path(logdir) / os.environ.get("MA_EXT_LOG", "ma_extensions.log"))
    path = Path(outdir) / "ma_extensions.json"
    result = json.loads(path.read_text()) if path.exists() else {}
    B = 40 if SMOKE else 1000
    plan = {"ratios": lambda: section_ratios(2000 if not SMOKE else 40, log),
            "within": lambda: section_within(2000 if not SMOKE else 30, log),
            "sensitivity": lambda: section_sensitivity(2000 if not SMOKE else 30, log),
            "seasons": lambda: section_seasons(2000 if not SMOKE else 30, log),
            "validation": lambda: section_validation(1000 if not SMOKE else 6, log),
            "table1": lambda: section_table1(20000 if SMOKE else 400000, log),
            "calibration": lambda: section_calibration(log)}
    for name in sections:
        log("start", name)
        result[name] = plan[name]()
        result.setdefault("meta", {})[name] = {"seed": SEED, "completed_utc": utc(),
                                              "smoke": SMOKE}
        write(path, result)
        log("done", name)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("sections", nargs="*", default=["table1", "ratios", "within", "sensitivity", "seasons", "validation"])
    p.add_argument("--output-dir", type=Path, default=EXP / "analysis")
    p.add_argument("--log-dir", type=Path, default=EXP / "logs")
    a = p.parse_args()
    main(a.sections, a.output_dir, a.log_dir)
