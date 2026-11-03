"""s04: facility-reporting completeness and a second, laboratory proxy (M6/D12-D13, m8/D27, m9/D28).

A. Reporting-completeness proxy p_hat (share of hospitals reporting) by reporting regime.
B. Delta a = a(statistic) - a(coverage-adjusted statistic), paired two-way bootstrap (B = 2000),
   pooled and per season; TOST with substantive margins fixed before fitting:
     G   : delta = tau0^2 (1 + v_R/Rbar^2)          (theory: a_H = v_H / Rbar^2, v_H = tau^2 (Rbar^2 + v_R))
     R_A : delta = kappa^2 tau0^2 (1 + v_R/Rbar^2)  (same perturbation of the weekly growth ratio, propagated
                                                     through the Wallinga-Lipsitch map; kappa = dlogR/dlog g)
   with tau0 = 0.10 (primary) and 0.05, v_R = 0.04; plus a normal-approximation power curve.
   The adjusted statistic divides each week's count by its p_hat (original r3_common.build_rows logic,
   experiment/original/code/r3_common.py lines 93-118).
C. Excess a - v_R/Rbar^2 for v_R in {0.04, 0.08, 0.16} for R_A (G versions are in s03).
D. Calendar-matched comparison of the voluntary-reporting period (May-Oct 2024) with May-Oct 2023/2025.
E. Second proxy: state laboratory influenza positives (Delphi fluview_clinical, 44 states). Share of the
   hospital intercept shared with laboratory growth (state-week), two-way bootstrap (B = 500).
"""
import datetime as dt
import json

import numpy as np
import pandas as pd

from cc_common import (DATA, GI_FLU, SEED, W_ORIG, Log, Part, boot_weights, fgls_affine, fit_stat, load_panel, pct,
                       rel_loss, spec_frame, wl_R)

log = Log("s04_reporting")
P = Part("s04_reporting", "s04_reporting.py",
         ["data/panel_nhsn_vintage_2026-10-20.csv.gz", "data/external/delphi_fluview_clinical_raw.json"], "measured", SEED + 4)
panel = load_panel()
B, B_E = 2000, 500
ORIG_END = pd.Timestamp("2026-08-29")       # last D_t week of the submitted manuscript's panel (10,303 rows)

# ----------------------------------------------------------------------------- A. p_hat by regime
reg = {}
for lab, sub in (("full_vintage", panel), ("cut_original_end", panel[panel.week_end <= ORIG_END])):
    for k in (1, 2, 3):
        s = sub[sub.period == k]
        P.num(f"regime{k}_N_{lab}", int(len(s)), note="state-weeks with D_t >= 5, all calendar months")
    reg[lab] = int(len(sub))
P.num("panel_N_cut_original_end", reg["cut_original_end"])
P.num("panel_N_full_vintage", reg["full_vintage"])
hist = {}
edges = np.linspace(0, 100, 41)
for k in (1, 2, 3):
    c = panel.loc[panel.period == k, "coverage"].dropna()
    hist[str(k)] = {"counts": np.histogram(c.clip(0, 100), edges)[0].tolist(), "N": int(len(c)),
                    "q": [float(v) for v in np.percentile(c, [5, 25, 50, 75, 95])], "mean": float(c.mean())}
    P.num(f"regime{k}_phat_median", float(np.median(c)), unit="% hospitals reporting")
    P.num(f"regime{k}_phat_share_lt90", float(np.mean(c < 90)))
P.put("phat_hist", {"edges": edges.tolist(), "by_regime": hist})
log("regime N", {k: v["N"] for k, v in hist.items()})


# ----------------------------------------------------------------------------- B. Delta a and TOST
def adjusted(e, stat):
    h = {j: e[c] / 100.0 for j, c in (("n", "coverage_next"), ("0", "coverage"), ("1", "coverage_prev1"), ("2", "coverage_prev2"))}
    ok = np.all([np.isfinite(v) & (v > 0) for v in h.values()], axis=0)
    if stat == "G":
        ds = W_ORIG[0] * e.C_t / h["0"] + W_ORIG[1] * e.C_prev1 / h["1"] + W_ORIG[2] * e.C_prev2 / h["2"]
        xs = (e.C_next / h["n"]) / ds
    else:
        xs = pd.Series(wl_R((e.C_next / h["n"] + 0.5) / (e.C_t / h["0"] + 0.5)), index=e.index)
    return xs.to_numpy(), ok.to_numpy() if hasattr(ok, "to_numpy") else ok


def kappa2(e):
    m, s = GI_FLU["mid"]
    shape, scale = (m / s) ** 2, s ** 2 / m
    lg = np.log((e.C_next + 0.5) / (e.C_t + 0.5))
    k = shape * (scale / 7) / (1 + scale * lg / 7)
    return float(np.mean(k[np.isfinite(k) & (1 + scale * lg / 7 > 0)] ** 2))


def delta_a(x, xs, D, w=None):
    return fit_stat(x, D, w)["a"] - fit_stat(xs, D, w)["a"]


def tost_power(se, margin, true=0.0, z=1.6448536):
    from scipy.stats import norm
    return float(max(0.0, norm.cdf((margin - true) / se - z) + norm.cdf((margin + true) / se - z) - 1))


VR = 0.04
for si, spec in enumerate(["G_rule4", "RA_rule4", "G_cal3"]):
    if f"delta|{spec}" in P.d["blocks"]:
        continue
    stat = "G" if spec.startswith("G") else "RA"
    e0 = spec_frame(panel, spec)
    xs, ok = adjusted(e0, stat)
    e = e0[ok].reset_index(drop=True); xs = xs[ok]
    k2 = kappa2(e) if stat == "RA" else 1.0
    Rb = float(e.x.mean())
    margins = {f"tau0_{t}": k2 * t ** 2 * (1 + VR / Rb ** 2) for t in (0.10, 0.05)}
    out = {"N": int(len(e)), "N_dropped_invalid_coverage": int((~ok).sum()), "kappa2": k2, "margins": margins, "by": {}}
    groups = [("pooled", np.ones(len(e), bool))] + [(s, (e.season == s).to_numpy()) for s in sorted(e.season.unique())]
    for gi, (g, m) in enumerate(groups):
        eg = e[m].reset_index(drop=True); xg = xs[m]
        pt = delta_a(eg.x.to_numpy(), xg, eg.D.to_numpy())
        rng = np.random.default_rng(SEED + 400 + 10 * si + gi)
        dr = np.array([delta_a(eg.x.to_numpy(), xg, eg.D.to_numpy(), boot_weights(eg, "twoway", rng)) for _ in range(B)])
        ci90, ci95 = pct(dr, (5, 95)), pct(dr, (2.5, 97.5))
        se = (ci90[1] - ci90[0]) / (2 * 1.6448536)      # percentile-based SE (draws are heavy-tailed)
        rec = {"N": int(len(eg)), "delta_a": pt, "ci90": ci90, "ci95": ci95, "se": se, "seed": SEED + 400 + 10 * si + gi,
               "a_unadj": fit_stat(eg.x, eg.D)["a"], "a_adj": fit_stat(xg, eg.D.to_numpy())["a"]}
        for mk, mv in margins.items():
            rec[f"equivalent_{mk}"] = bool(-mv < ci90[0] and ci90[1] < mv)
            rec[f"power0_{mk}"] = tost_power(se, mv)
        out["by"][g] = rec
        log("delta", spec, g, {k: (round(v, 5) if isinstance(v, float) else v) for k, v in rec.items()})
    se0 = out["by"]["pooled"]["se"]
    out["power_curve"] = {mk: {"true_delta": list(np.linspace(-mv * 1.5, mv * 1.5, 61)),
                               "power": [tost_power(se0, mv, t) for t in np.linspace(-mv * 1.5, mv * 1.5, 61)]}
                          for mk, mv in margins.items()}
    # power as a function of the number of seasons (SE scaled by sqrt(S_obs/S)), margin tau0 = 0.05
    S_obs = e.season.nunique()
    out["power_vs_S"] = {mk: {str(S): tost_power(se0 * np.sqrt(S_obs / S), mv) for S in (2, 3, 4, 6, 8, 12)}
                         for mk, mv in margins.items()}
    P.put(f"delta|{spec}", out)
    pr = out["by"]["pooled"]
    P.num(f"{spec}_delta_a", pr["delta_a"], ci=pr["ci95"], seed=pr["seed"], B=B, ci90=pr["ci90"])
    for mk, mv in margins.items():
        P.num(f"{spec}_tost_margin_{mk}", mv, mode="analytic")
        P.num(f"{spec}_tost_equivalent_{mk}", pr[f"equivalent_{mk}"])
        P.num(f"{spec}_tost_power0_{mk}", pr[f"power0_{mk}"])
    P.num(f"{spec}_delta_a_seasons_equivalent_tau0_0.1",
          int(sum(v["equivalent_tau0_0.1"] for k, v in out["by"].items() if k != "pooled")),
          note=f"of {len(out['by']) - 1} seasons")

# ----------------------------------------------------------------------------- C. excess for R_A
e = spec_frame(panel, "RA_rule4")
f = fit_stat(e.x, e.D)
for vr in (0.04, 0.08, 0.16):
    P.num(f"RA_rule4_excess_vR{vr}", f["a"] - vr / f["xbar"] ** 2, mode="analytic",
          note="generation-scale v_R compared with the weekly R_A intercept; see s08 simulation for the mapping")
P.num("RA_rule4_implied_vR", f["a"] * f["xbar"] ** 2, mode="analytic")

# ----------------------------------------------------------------------------- D. calendar-matched voluntary period
if "calendar_matched" not in P.d["blocks"]:
    cm = {}
    for yi, y in enumerate((2023, 2024, 2025)):
        s = panel[(panel.week_end >= f"{y}-05-01") & (panel.week_end <= f"{y}-10-31")].reset_index(drop=True)
        s = s.assign(x=s.G, D=s.D_t)
        fs = fit_stat(s.x, s.D)
        rng = np.random.default_rng(SEED + 450 + yi)
        lg = np.log((s.C_next + 0.5) / (s.D + 0.5)).to_numpy()

        def wiqr(v, w):
            o = np.argsort(v); cw = np.cumsum(w[o]) / w.sum()
            return float(v[o][np.searchsorted(cw, 0.75)] - v[o][np.searchsorted(cw, 0.25)])
        dr = np.array([[np.sqrt(np.cov(s.x, aweights=w)), wiqr(lg, w)]
                       for w in (boot_weights(s, "twoway", rng) for _ in range(1000))])
        cm[str(y)] = {"N": int(len(s)), "regime": int(s.period.mode()[0]), "sd_G": float(s.x.std(ddof=1)),
                      "sd_G_ci95": pct(dr[:, 0]), "iqr_logG": wiqr(lg, np.ones(len(lg))), "iqr_logG_ci95": pct(dr[:, 1]),
                      "a_fgls_unstable": fs["a"], "Gbar": fs["xbar"],
                      "D_median": float(s.D.median()), "phat_median": float(s.coverage.median())}
        log("calendar-matched", y, cm[str(y)])
    P.put("calendar_matched", cm)
    for y, v in cm.items():
        P.num(f"mayoct{y}_N", v["N"]); P.num(f"mayoct{y}_sd_G", v["sd_G"], ci=v["sd_G_ci95"], B=1000)
        P.num(f"mayoct{y}_iqr_logG", v["iqr_logG"], ci=v["iqr_logG_ci95"], B=1000,
              note="robust spread; the FGLS intercept is not stable in the voluntary period (negative point estimate)")


# ----------------------------------------------------------------------------- E. laboratory proxy
def mmwr_end(ew: int) -> pd.Timestamp:
    y, w = divmod(int(ew), 100)
    jan4 = dt.date(y, 1, 4)
    start = jan4 - dt.timedelta(days=(jan4.weekday() + 1) % 7)
    return pd.Timestamp(start + dt.timedelta(days=7 * (w - 1) + 6))


if "lab_proxy" not in P.d["blocks"]:
    lab = pd.DataFrame(json.loads((DATA / "external/delphi_fluview_clinical_raw.json").read_text())["epidata"])
    lab["state"] = lab.region.str.upper()
    lab["week_end"] = lab.epiweek.map(mmwr_end)
    lab["pos"] = lab.total_a + lab.total_b
    lab = lab.sort_values(["state", "week_end"])
    lab["pos_next"] = lab.groupby("state").pos.shift(-1)
    lab["wk_next"] = lab.groupby("state").week_end.shift(-1)
    lab = lab[(lab.wk_next - lab.week_end) == pd.Timedelta(days=7)]
    e = spec_frame(panel, "RA_rule4")
    m = e.merge(lab[["state", "week_end", "pos", "pos_next", "total_specimens"]], on=["state", "week_end"], how="inner")
    m = m[m.pos >= 5].reset_index(drop=True)
    m["x_lab"] = wl_R((m.pos_next + 0.5) / (m.pos + 0.5))

    def lab_share(w=None):
        ww = np.ones(len(m)) if w is None else w
        xb = np.average(m.x, weights=ww); lb = np.average(m.x_lab, weights=ww)
        beta = np.sum(ww * (m.x - xb) * (m.x_lab - lb)) / np.sum(ww * (m.x_lab - lb) ** 2)
        xr = m.x - beta * (m.x_lab - lb)
        a_h = fit_stat(m.x.to_numpy(), m.D.to_numpy(), ww)["a"]
        a_r = fit_stat(xr.to_numpy(), m.D.to_numpy(), ww)["a"]
        r = np.sum(ww * (m.x - xb) * (m.x_lab - lb)) / np.sqrt(np.sum(ww * (m.x - xb) ** 2) * np.sum(ww * (m.x_lab - lb) ** 2))
        return a_h, a_r, 1 - a_r / a_h, beta, r
    pt = lab_share()
    rng = np.random.default_rng(SEED + 470)
    dr = np.array([lab_share(boot_weights(m, "twoway", rng)) for _ in range(B_E)])
    blk = {"N": int(len(m)), "n_states": int(m.state.nunique()), "n_states_lab": int(lab.state.nunique()),
           "a_hosp": pt[0], "a_resid": pt[1], "share_shared_with_lab": pt[2], "beta": pt[3], "corr": pt[4],
           "ci95": {k: pct(dr[:, i]) for i, k in enumerate(["a_hosp", "a_resid", "share", "beta", "corr"])},
           "B": B_E, "seed": SEED + 470}
    P.put("lab_proxy", blk)
    P.num("lab_N", blk["N"]); P.num("lab_n_states", blk["n_states"])
    P.num("lab_share_of_a_shared", pt[2], ci=blk["ci95"]["share"], seed=SEED + 470, B=B_E)
    P.num("lab_corr_RA", pt[4], ci=blk["ci95"]["corr"], seed=SEED + 470, B=B_E)
    log("lab proxy", blk)
log.close()
