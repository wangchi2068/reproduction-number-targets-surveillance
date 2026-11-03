"""s03: what the intercept a contains and how stable b/a is (M3/D6-D8, M4/D9, m8/D27).

1. Hierarchical decomposition of a (two-way bootstrap intervals, B = 1000):
     national-week component = a_total - a_weekFE
     state-season component  = a_weekFE - a_both
     residual                = a_both
   a_weekFE: deviations from the cross-state mean of the same week (reviewer check
   review/checks/empirical_checks.py, extended with a degrees-of-freedom factor n_w/(n_w-1));
   a_both: deviations from an additive week + state-season fit (alternating projections),
   loss scaled by N/(N - n_week - n_stateseason + 1).
2. Plateau (a + b/D) versus power-law decay (c D^-gamma): leave-one-state-out and
   leave-one-season-out cross-validated prediction error with model-neutral weights, an in-sample
   weighted-SSR contrast calibrated by a residual bootstrap under the power-law null, and the
   minimum detectable plateau (power 0.8).
3. Forest of b/a across specifications (two-way bootstrap, B = 500 each).
4. Theory-data comparison for G with the b-delta correction (M4), scenarios of
   review/checks/theory_checks.py re-evaluated at the window mean of this workspace.
"""
import sys

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from cc_common import (SEED, Log, Part, boot_weights, fgls_affine, fit_stat, load_panel, mask_windows, pct,
                       rel_loss, spec_frame, wl_R)

log = Log("s03_decomposition")
P = Part("s03_decomposition", "s03_decomposition.py", ["data/panel_nhsn_vintage_2026-10-20.csv.gz"], "measured", SEED + 3)
panel = load_panel()
B_DEC = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
B_FOR = int(sys.argv[2]) if len(sys.argv) > 2 else 500
N_PB = int(sys.argv[3]) if len(sys.argv) > 3 else 200
SPECS_MAIN = ["RA_rule4", "G_rule4", "G_cal3"]


# ----------------------------------------------------------------------------- 1. decomposition
def wmean_by(codes, v, w, n):
    num = np.bincount(codes, weights=v * w, minlength=n)
    den = np.bincount(codes, weights=w, minlength=n)
    return np.divide(num, den, out=np.zeros(n), where=den > 0)


def decompose(e, w=None):
    x, D = e.x.to_numpy(), e.D.to_numpy()
    w = np.ones(len(e)) if w is None else w
    k = w > 0
    x, D, w = x[k], D[k], w[k]
    wk = pd.factorize(e.week_end.to_numpy()[k])[0]
    ss = pd.factorize((e.state + "|" + e.season).to_numpy()[k])[0]
    nw, ns = wk.max() + 1, ss.max() + 1
    L0, _ = rel_loss(x, w)
    a0, b0 = fgls_affine(L0, D, w)
    # week FE
    mw = wmean_by(wk, x, w, nw)
    cnt_w = np.bincount(wk, minlength=nw)            # distinct rows per week (dof factor)
    f_w = np.where(cnt_w[wk] > 1, cnt_w[wk] / np.maximum(cnt_w[wk] - 1, 1), np.nan)
    Lw = ((x - mw[wk]) / mw[wk]) ** 2 * f_w
    ok = np.isfinite(Lw)
    a1, b1 = fgls_affine(Lw[ok], D[ok], w[ok])
    # week + state-season additive FE
    mu = np.average(x, weights=w)
    al, be = np.zeros(nw), np.zeros(ns)
    for _ in range(50):
        al = wmean_by(wk, x - mu - be[ss], w, nw)
        be = wmean_by(ss, x - mu - al[wk], w, ns)
    fit = mu + al[wk] + be[ss]
    dof = len(np.unique(np.column_stack([wk, ss]), axis=0))
    f_b = dof / max(dof - nw - ns + 1, 1)
    Lb = ((x - fit) / fit) ** 2 * f_b
    a2, b2 = fgls_affine(Lb, D, w)
    return {"a_total": a0, "b_total": b0, "a_weekFE": a1, "b_weekFE": b1, "a_both": a2, "b_both": b2,
            "comp_week": a0 - a1, "comp_stateseason": a1 - a2, "comp_resid": a2}


for si, spec in enumerate(SPECS_MAIN):
    if f"decomp|{spec}" in P.d["blocks"]:
        continue
    e = spec_frame(panel, spec)
    d0 = decompose(e)
    dev = e.x / e.x.mean() - 1
    wkm = dev.groupby(e.week_end).mean()
    icc = float(wkm.var() / dev.var())
    rng = np.random.default_rng(SEED + 30 + si)
    draws = []
    for r in range(B_DEC):
        draws.append(decompose(e, boot_weights(e, "twoway", rng)))
    dr = pd.DataFrame(draws)
    blk = {k: {"value": v, "ci95": pct(dr[k])} for k, v in d0.items()}
    blk["icc_week"] = icc
    blk["B"] = B_DEC
    P.put(f"decomp|{spec}", blk)
    for k in ("a_total", "a_weekFE", "a_both", "comp_week", "comp_stateseason", "comp_resid"):
        P.num(f"{spec}_{k}", d0[k], ci=blk[k]["ci95"], seed=SEED + 30 + si, B=B_DEC)
    P.num(f"{spec}_share_week", d0["comp_week"] / d0["a_total"])
    P.num(f"{spec}_icc_week", icc)
    log("decomp", spec, {k: round(v, 5) for k, v in d0.items()}, "icc", round(icc, 3))


# ----------------------------------------------------------------------------- 2. plateau vs power law
def fit_power(L, D, w=None, it=30):
    w = np.ones_like(L) if w is None else w
    lD = np.log(D)
    p = np.polyfit(lD, np.log(np.clip(L, 1e-8, None)), 1)
    th = np.array([np.log(max(np.mean(L), 1e-6)) + 0.0, max(-p[0], 0.05)])
    th[0] = np.log(np.mean(L * D ** th[1]))
    for _ in range(it):
        mu = np.exp(th[0] - th[1] * lD)
        sw = np.sqrt(w) / mu
        res = least_squares(lambda t: (L - np.exp(t[0] - t[1] * lD)) * sw, th, method="lm")
        if np.max(np.abs(res.x - th)) < 1e-9:
            th = res.x
            break
        th = res.x
    return float(np.exp(th[0])), float(th[1])


def neutral_w(L, D, nb=10):
    q = pd.qcut(D, nb, labels=False, duplicates="drop")
    m = pd.Series(L).groupby(q).transform("mean").to_numpy()
    return 1 / m ** 2


def wssr(L, mu, wn):
    return float(np.sum(wn * (L - mu) ** 2))


def contrast(L, D, wn):
    a, b = fgls_affine(L, D)
    c, g = fit_power(L, D)
    return wssr(L, c * D ** -g, wn) - wssr(L, a + b / D, wn), (a, b), (c, g)


def cv(e, L, key, wn):
    out = []
    for _, idx in e.groupby(key).indices.items():
        tr = np.ones(len(e), bool); tr[idx] = False
        a, b = fgls_affine(L[tr], e.D.to_numpy()[tr])
        c, g = fit_power(L[tr], e.D.to_numpy()[tr])
        Dh = e.D.to_numpy()[idx]
        out.append((wssr(L[idx], c * Dh ** -g, wn[idx]) - wssr(L[idx], a + b / Dh, wn[idx])) / len(e))
    out = np.array(out)
    return {"sum": float(out.sum()), "folds": len(out), "folds_affine_better": int((out > 0).sum()),
            "se": float(out.std(ddof=1) * np.sqrt(len(out)))}


for si, spec in enumerate(["RA_rule4", "G_rule4", "G_cal3"]):
    if f"plateau|{spec}" in P.d["blocks"]:
        continue
    e = spec_frame(panel, spec)
    L, _ = rel_loss(e.x)
    D = e.D.to_numpy()
    wn = neutral_w(L, D)
    obs, (a, b), (c, g) = contrast(L, D, wn)
    cv_state = cv(e, L, "state", wn)
    cv_season = cv(e, L, "season", wn)
    # residual bootstrap: ratios L/mu resampled within D deciles
    rng = np.random.default_rng(SEED + 40 + si)
    dec = pd.qcut(D, 10, labels=False, duplicates="drop")

    def sim(mu, eps):
        out = np.empty_like(mu)
        for q in np.unique(dec):
            idx = np.where(dec == q)[0]
            out[idx] = mu[idx] * rng.choice(eps[idx], len(idx))
        return out
    mu_p = c * D ** -g
    eps_p = L / mu_p
    null = np.array([contrast(Ls, D, neutral_w(Ls, D))[0] for Ls in (sim(mu_p, eps_p) for _ in range(N_PB))])
    crit = float(np.quantile(null, 0.95))
    pval = float((1 + np.sum(null >= obs)) / (1 + len(null)))
    mu_a = a + b / D
    eps_a = L / mu_a
    power = {}
    for f in (0.25, 0.5, 0.75, 1.0):
        mu_f = f * a + b / D
        st = np.array([contrast(Ls, D, neutral_w(Ls, D))[0] for Ls in (sim(mu_f, eps_a) for _ in range(N_PB // 2))])
        power[str(f)] = float(np.mean(st > crit))
        log("power", spec, f, power[str(f)])
    mdd = next((float(f) * a for f in (0.25, 0.5, 0.75, 1.0) if power[str(f)] >= 0.8), None)
    blk = {"affine": [a, b], "power": [c, g], "contrast_obs": obs, "null_crit95": crit, "p_value": pval,
           "cv_leave_state": cv_state, "cv_leave_season": cv_season, "power_by_fraction_of_a": power,
           "min_detectable_a": mdd, "N_null": N_PB, "N_power": N_PB // 2,
           "note": "contrast > 0 favours the plateau; residuals resampled iid within D deciles (ignores clustering, so power is optimistic)"}
    # binned means for Fig 4B
    qb = pd.qcut(D, 12, labels=False, duplicates="drop")
    bins = pd.DataFrame({"D": D, "L": L, "q": qb}).groupby("q").agg(D_mean=("D", "mean"), D_min=("D", "min"),
                                                                    D_max=("D", "max"), L_mean=("L", "mean"),
                                                                    L_sd=("L", "std"), n=("L", "size"))
    blk["bins"] = bins.reset_index(drop=True).to_dict("records")
    P.put(f"plateau|{spec}", blk)
    P.num(f"{spec}_power_exponent", g, note="power-law fit E(L)=c D^-gamma")
    P.num(f"{spec}_plateau_contrast_p", pval, seed=SEED + 40 + si, note="residual bootstrap under power-law null")
    P.num(f"{spec}_cv_state_folds_affine_better", cv_state["folds_affine_better"], note=f"of {cv_state['folds']}")
    P.num(f"{spec}_cv_season_folds_affine_better", cv_season["folds_affine_better"], note=f"of {cv_season['folds']}")
    P.num(f"{spec}_min_detectable_a", mdd, seed=SEED + 40 + si)
    log("plateau", spec, {k: v for k, v in blk.items() if k not in ("bins",)})


# ----------------------------------------------------------------------------- 3. forest of b/a
def huber_affine(L, D, w=None, c=1.345, it=50):
    w = np.ones_like(L) if w is None else w
    a, b = fgls_affine(L, D, w)
    for _ in range(it):
        mu = np.clip(a + b / D, 1e-6, None)
        r = (L - mu) / mu
        s = np.median(np.abs(r - np.median(r))) / 0.6745
        u = np.abs(r) / (c * s)
        hw = np.where(u <= 1, 1, 1 / u)
        na, nb = fgls_affine(L, D, w * hw)
        if abs(na - a) + abs(nb - b) < 1e-10:
            break
        a, b = na, nb
    return a, b


PEAKS4 = ["2022-12-03", "2023-12-30", "2025-02-08", "2026-01-03"]   # season maxima (original + 2025-26)


def peak_frame(stat, dcol):
    win = [(d, str((pd.Timestamp(d) + pd.Timedelta(weeks=7)).date())) for d in PEAKS4]
    q = panel[mask_windows(panel, win)].copy()
    q = q[q[dcol] >= 5]
    q["x"], q["D"] = q[stat], q[dcol]
    return q[np.isfinite(q.x)].reset_index(drop=True)


def forest_specs():
    out = []
    for base in ("RA_rule4", "G_rule4"):
        e = spec_frame(panel, base)
        out.append((base, "main", e, None))
        L, _ = rel_loss(e.x)
        out.append((base, "trim_top1pct", e[L <= np.quantile(L, 0.99)].reset_index(drop=True), None))
        for s in sorted(e.season.unique()):
            out.append((base, f"drop_{s}", e[e.season != s].reset_index(drop=True), None))
        out.append((base, "huber", e, "huber"))
        out.append((base, "weekFE", e, "weekFE"))
        stat, dcol = ("RA_mid", "C_t") if base.startswith("RA") else ("G", "D_t")
        out.append((base, "peak_8w", peak_frame(stat, dcol), None))
    e = spec_frame(panel, "RA_rule4")
    for gi in ("short", "long"):
        e2 = e.copy(); e2["x"] = e2[f"RA_{gi}"]
        out.append(("RA_rule4", f"GI_{gi}", e2, None))
    out.append(("RA_rule4", "S3_calendar", spec_frame(panel, "RA_cal3"), None))
    out.append(("G_rule4", "S3_calendar", spec_frame(panel, "G_cal3"), None))
    e = spec_frame(panel, "G_rule4")
    for wname, wts in (("single_1w", (1, 0, 0)), ("equal_3w", (1 / 3, 1 / 3, 1 / 3))):
        e2 = e.copy()
        e2["D"] = wts[0] * e2.C_t + wts[1] * e2.C_prev1 + wts[2] * e2.C_prev2
        e2 = e2[e2.D >= 5].copy()
        e2["x"] = e2.C_next / e2.D
        out.append(("G_rule4", f"D_weights_{wname}", e2.reset_index(drop=True), None))
    return out


def est(e, how, w=None):
    if how == "huber":
        L, _ = rel_loss(e.x, w)
        a, b = huber_affine(L, e.D.to_numpy(), w)
        return b / a
    if how == "weekFE":
        d = decompose(e, w)
        return d["b_weekFE"] / d["a_weekFE"]
    return fit_stat(e.x.to_numpy(), e.D.to_numpy(), w)["ba"]


rows = P.d["blocks"].get("forest", [])
done = {(r["base"], r["spec"]) for r in rows}
for fi, (base, name, e, how) in enumerate(forest_specs()):
    if (base, name) in done:
        continue
    pt = est(e, how)
    rng = np.random.default_rng(SEED + 500 + fi)
    dr = [est(e, how, boot_weights(e, "twoway", rng)) for _ in range(B_FOR)]
    rows.append({"base": base, "spec": name, "N": int(len(e)), "ba": float(pt), "ci95": pct(dr),
                 "seed": SEED + 500 + fi, "B": B_FOR})
    P.put("forest", rows)
    log("forest", base, name, len(e), round(pt, 2), [round(v, 1) for v in pct(dr)])
for r in rows:
    P.num(f"forest_{r['base']}_{r['spec']}_ba", r["ba"], ci=r["ci95"], seed=r["seed"], B=r["B"], N=r["N"])

# ----------------------------------------------------------------------------- 4. theory vs data (G, b-delta)
SW2 = sum(w ** 2 for w in (0.65, 0.25, 0.10))
for spec in ("G_cal3", "G_rule4"):
    e = spec_frame(panel, spec)
    Rb = float(e.x.mean())
    sc = []
    for vr in (0.04, 0.08):
        for kk in (0.2, 0.5):
            for rr in (0.0056, 0.021):
                A = Rb + rr * (Rb ** 2 + vr) / kk
                a_th = vr / Rb ** 2
                b_naive = A / Rb ** 2
                b_delta = (A + (1 - rr) * (Rb ** 2 + vr) * SW2) / Rb ** 2
                sc.append({"vR": vr, "k": kk, "rho": rr, "a_th": a_th, "b_naive": b_naive, "b_delta": b_delta,
                           "ba_naive": b_naive / a_th, "ba_delta": b_delta / a_th})
    s = pd.DataFrame(sc)
    f = fit_stat(e.x, e.D)
    P.put(f"theory_scenarios|{spec}", {"Rbar": Rb, "sum_w2": SW2, "rows": sc, "a_emp": f["a"], "b_emp": f["b"]})
    P.num(f"{spec}_theory_ba_delta_range", [float(s.ba_delta.min()), float(s.ba_delta.max())], mode="analytic",
          note="scenarios vR {0.04,0.08} x k {0.2,0.5} x rho {0.0056,0.021}, b-delta corrected")
    P.num(f"{spec}_theory_ba_naive_range", [float(s.ba_naive.min()), float(s.ba_naive.max())], mode="analytic")
    P.num(f"{spec}_a_over_theory_range", [float(f["a"] / s.a_th.max()), float(f["a"] / s.a_th.min())], mode="analytic")
    P.num(f"{spec}_b_over_theory_delta_range", [float(f["b"] / s.b_delta.max()), float(f["b"] / s.b_delta.min())], mode="analytic")
    P.num(f"{spec}_b_over_theory_naive_range", [float(f["b"] / s.b_naive.max()), float(f["b"] / s.b_naive.min())], mode="analytic")
    for vr in (0.04, 0.08, 0.16):
        P.num(f"{spec}_excess_vR{vr}", f["a"] - vr / Rb ** 2, mode="analytic",
              note="a minus vR/Rbar^2 (M6 item 3)")
    log("theory", spec, Rb, s[["ba_naive", "ba_delta"]].agg(["min", "max"]).to_dict())
log.close()
