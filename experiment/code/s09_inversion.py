"""s09: round-3 additions -- two-step inversion, multiscale decomposition, reporting regimes,
diminishing returns and the contrast with Parag et al. (2022, eq. 5).

All new numbers go to five part files (keys inversion.*, multiscale.*, reporting_regime.*, returns.*,
parag.* in results.json); existing parts are read, never rewritten.

A. Two-step inversion (reviewer N1/D2). Step 1: a data-matched version of the s08 renewal simulator maps
   a latent environmental variance v_R to the expected fitted intercept E[a_hat](v_R) of each apparent
   statistic (R_A and G) under the empirical design. Step 2: the empirical intercept is inverted on that
   calibration line. Data-matched settings (M0):
   - the four rule-4 windows with their actual lengths (9, 12, 8, 13 observation weeks) and 51 states per
     season, i.e. one simulated panel has the cross-section and length of the analysis panel;
   - season mean R_s = wl_R(g_s) with g_s the geometric-mean national weekly growth over the window weeks
     (national_flu.csv), so the simulated growth matches the observed growth season by season;
   - initial weekly admissions of state i in season s equal its observed admissions in the first window
     week (raw NHSN, floored at 0.5);
   - R_{i,s,w} = R_s + eta_{s,w} + eps_{i,s,w}: a national synchronous week shock eta (variance rho v_R,
     shared by all states of a panel replicate) plus a state-specific shock (variance (1-rho) v_R),
     Gaussian, clipped at 0.05; rho = 0.45 in M0; optional AR(1) persistence phi across weeks;
   - daily NB renewal (k = 0.35), binomial admission (rho_h = 0.01), gamma infection-to-admission delay
     (5.5, 3.0 d), GI gamma (3.6, 1.6 d) as in s08; statistics computed exactly as in the data
     (R_A with the analyst's GI 3.6/1.6 d; G with W_ORIG), rows kept when D_t >= 5 (both) and C_t >= 5 (R_A).
   Sensitivity settings change one ingredient each (rho 0 / 0.65, AR phi 0.5, true GI short/long with
   the analyst keeping the mid GI, delay 3.5/2 or 8/4 d, k 0.1, a synchronous multiplicative reporting
   shock on counts instead of a transmission shock, and the s08 design as a bridge).
   v_R grid 0.02/0.04/0.06/0.09, NREP_A = 200 panel replicates per cell. The calibration line
   E[a_hat] = c0 + c1 v_R is fitted by weighted least squares (weights 1/MC variance of the cell mean);
   kappa = E[a_hat]/(v_R/Rbar^2) with Rbar the row-weighted true mean R, and kappa_inf from a noise-free
   (mean-field) run of the same design.
   Uncertainty of v_R^implied: B = 2000 joint two-way bootstrap draws of (a_RA, a_G) on shared
   state/week multiplicities, each combined with an independent draw of (c0, c1) from its Monte Carlo
   sampling distribution (normal, WLS covariance); 95% percentile interval.
   Calibrated theoretical b/a: the same M0 design at v_R in {0.04, 0.08} x k in {0.2, 0.5} x
   rho_h in {0.0056, 0.021} (the s03 theory scenarios); pooled b/a = E[b_hat]/E[a_hat] per cell and the
   2.5-97.5% range of per-replicate b/a (what one panel of this size would show).
B. Multiscale decomposition: share of the intercept removed by week fixed effects and the week ICC,
   with joint two-way bootstrap intervals; SD of a_hat under state / week / two-way clustering (s02
   draws); variance of a J-state pooled growth deviation (subsampling, plus the one-way random-effects
   model V_J = s2_week + (V_1 - s2_week)/J).
C. Reporting regime: equivalence margin as a share of a (sample-mean normalisation, as s04), observed
   Delta a / a with joint bootstrap interval; May-October 2024 (voluntary) versus 2023/2025 (mandatory)
   growth-ratio dispersion with a bootstrap interval for the ratio, scale-matched (D reweighting),
   month-by-month and other-pathogen (COVID-19, RSV) checks; participation under the mandate.
D. Diminishing returns from b/a (fraction of the attainable precision, gain per doubling, m*).
E. Poisson Fisher information (Parag et al. 2022, eqs. 3-5) versus v_R + A/m.

Usage: python3 s09_inversion.py [sim|A|B|C|D|E|all]   (default all; sim rows are resume-safe)
Outputs: analysis/parts/{inversion,multiscale,reporting_regime,returns,parag}.json,
logs/s09_inversion.log, logs/s09_sim_rows.jsonl, analysis/parts/s09_joint_draws.npz.
"""
import json
import sys
import time

import numpy as np
import pandas as pd

from cc_common import (DATA, GI_FLU, LOGS, PANEL, PARTS, SEED, W_ORIG, Log, Part, build_panel, fgls_affine,
                       fit_stat, load_panel, load_raw, pct, rel_loss, rule_windows, spec_frame, wl_R)
from s08_sim_apparent import disc_gamma, growth_rate

INPUTS = ["data/panel_nhsn_vintage_2026-10-20.csv.gz", "data/national_flu.csv",
          "data/external/nhsn_hrd_vintage_2026-10-20.csv.gz"]
ROWS = LOGS / "s09_sim_rows.jsonl"
DRAWS = PARTS / "s09_joint_draws.npz"
NST = 51
PRE_W = 4
NREP_A = 200
CHUNK = 25
B_INV = 2000
VGRID = (0.02, 0.04, 0.06, 0.09)
BAND = (0.04, 0.08)                       # theory reference band of v_R used in the manuscript
SW2 = {"G": float(np.sum(np.square(W_ORIG))), "RA": 1.0}
SEED_SIM = "SeedSequence([SEED, 9, cell_index, chunk_index])"

BASE = dict(rho=0.45, phi=0.0, gi="mid", delay=(5.5, 3.0), k=0.35, rho_h=0.01, sig_c=0.0, design="matched")
SETTINGS = {
    "M0": {},
    "rho0": {"rho": 0.0},
    "rho65": {"rho": 0.65},
    "ar05": {"phi": 0.5},
    "gi_short": {"gi": "short"},
    "gi_long": {"gi": "long"},
    "delay_short": {"delay": (3.5, 2.0)},
    "delay_long": {"delay": (8.0, 4.0)},
    "k01": {"k": 0.1},
    "count_shock": {"rho": 0.0, "sig_c": 0.12},
    "bridge_s08": {"design": "s08", "rho": 0.0},
}
LABEL = {"M0": "data-matched (M0)", "rho0": "no national shock", "rho65": "national share 0.65",
         "ar05": "AR(1) phi = 0.5", "gi_short": "true GI short (2.2 d)", "gi_long": "true GI long (4.3 d)",
         "delay_short": "delay 3.5 d", "delay_long": "delay 8 d", "k01": "k = 0.1",
         "count_shock": "synchronous reporting shock", "bridge_s08": "s08 design (bridge)"}
CELLS = [(f"{s}|v{v}", s, dict(BASE, **SETTINGS[s]), v) for s in SETTINGS for v in VGRID]
for k_ in (0.2, 0.5):
    for rh in (0.0056, 0.021):
        for v in (0.04, 0.08):
            CELLS.append((f"ba_k{k_}_r{rh}|v{v}", f"ba_k{k_}_r{rh}", dict(BASE, k=k_, rho_h=rh), v))

log = Log("s09_inversion")


# ============================================================================ design inputs
def design_inputs():
    nat = pd.read_csv(DATA / "national_flu.csv", parse_dates=["week_end"])
    win, _ = rule_windows(nat)
    raw = load_raw("flu")
    states = sorted(raw["Geographic aggregation"].unique())
    out = []
    for s, (a, b) in win.items():
        sub = nat[(nat.week_end >= a) & (nat.week_end <= b)]
        L = len(sub)
        g = float(np.exp(np.log(sub.g).mean()))
        c0 = raw[raw.week_end == a].set_index("Geographic aggregation")["C"].reindex(states).fillna(0).to_numpy()
        out.append({"season": s, "start": a.strftime("%Y-%m-%d"), "end": b.strftime("%Y-%m-%d"), "L": L,
                    "g_nat": g, "d0": np.maximum(c0, 0.5)})
    return out, states


DESIGN, STATES = design_inputs()


# ============================================================================ simulator
def ar_noise(rng, shape, phi):
    z = rng.standard_normal(shape)
    if phi == 0:
        return z
    out = np.empty(shape)
    out[..., 0] = z[..., 0]
    c = np.sqrt(1 - phi ** 2)
    for w in range(1, shape[-1]):
        out[..., w] = phi * out[..., w - 1] + c * z[..., w]
    return out


def sim_block(n, L, rbar, d0, vr, cfg, rng, meanfield=False):
    """n panel replicates x len(d0) states of one season block. Returns weekly admissions (S, W)."""
    S = n * len(d0)
    W = PRE_W + L + 2
    gk = disc_gamma(*GI_FLU[cfg["gi"]], s0=1)
    dk = disc_gamma(*cfg["delay"], s0=0)
    ndays, smax = 7 * W, len(gk)
    hist = smax + len(dk)
    r = growth_rate(rbar, gk)
    d0s = np.tile(d0, n)
    mgf = np.sum(dk * np.exp(-r * np.arange(len(dk))))
    i0 = d0s / (cfg["rho_h"] * mgf * np.exp(7 * r) ** PRE_W * np.sum(np.exp(r * np.arange(7))))
    if cfg["design"] == "s08":                                   # gamma iid weeks, as in s08
        shp = rbar ** 2 / vr
        Rw = rng.gamma(shp, rbar / shp, (S, W))
    else:
        eta = np.repeat(ar_noise(rng, (n, W), cfg["phi"]) * np.sqrt(cfg["rho"] * vr), len(d0), axis=0)
        eps = ar_noise(rng, (S, W), cfg["phi"]) * np.sqrt((1 - cfg["rho"]) * vr)
        Rw = np.clip(rbar + eta + eps, 0.05, None)
    I = np.zeros((S, hist + ndays))
    I[:, :hist] = i0[:, None] * np.exp(r * np.arange(-hist, 0))[None, :]
    if not meanfield:
        I[:, :hist] = np.round(I[:, :hist])
    grev = gk[::-1]
    k = cfg["k"]
    for d in range(ndays):
        j = hist + d
        lam = I[:, j - smax:j] @ grev
        Rd = Rw[:, d // 7]
        if meanfield:
            I[:, j] = Rd * lam
        else:
            x = np.zeros(S)
            pos = lam > 0
            if pos.any():
                x[pos] = rng.gamma(k * lam[pos], Rd[pos] / k)
            I[:, j] = rng.poisson(x)
    nd = hist + ndays
    out = np.zeros((S, nd + len(dk)))
    if meanfield:
        adm = cfg["rho_h"] * I
        for j, p in enumerate(dk):
            out[:, j:j + nd] += p * adm
    else:
        adm = rng.binomial(I.astype(np.int64), cfg["rho_h"])
        M = rng.multinomial(adm, dk)
        for j in range(len(dk)):
            out[:, j:j + nd] += M[:, :, j]
    C = out[:, hist:hist + ndays].reshape(S, W, 7).sum(axis=2)
    if cfg["sig_c"] > 0:                                          # synchronous multiplicative reporting shock
        sc = cfg["sig_c"]
        mult = np.repeat(np.exp(sc * rng.standard_normal((n, W)) - sc ** 2 / 2), len(d0), axis=0)
        C = C * mult if meanfield else np.round(C * mult)
    return C, Rw


def block_rows(C, Rw, n, nst, L, season_idx):
    """Long arrays of the observation weeks of one block: D_t, C_t, G, R_A, true R, rep, state, week id."""
    w = np.asarray(W_ORIG)
    t = np.arange(PRE_W, PRE_W + L)
    D = w[0] * C[:, t] + w[1] * C[:, t - 1] + w[2] * C[:, t - 2]
    Ct, Cn = C[:, t], C[:, t + 1]
    with np.errstate(divide="ignore", invalid="ignore"):
        G = Cn / np.where(D > 0, D, np.nan)
    RA = wl_R((Cn + 0.5) / (Ct + 0.5), *GI_FLU["mid"])
    S = C.shape[0]
    rep = np.repeat(np.arange(n), nst)[:, None] * np.ones((1, L), int)
    st = np.tile(np.arange(nst), n)[:, None] * np.ones((1, L), int)
    wk = season_idx * 100 + np.arange(L)[None, :] * np.ones((S, 1), int)
    return {"D": D.ravel(), "Ct": Ct.ravel(), "G": G.ravel(), "RA": RA.ravel(), "R": Rw[:, t].ravel(),
            "rep": rep.ravel(), "state": st.ravel() + 1000 * season_idx, "week": wk.ravel()}


def wmean_by(codes, v, w, n):
    num = np.bincount(codes, weights=v * w, minlength=n)
    den = np.bincount(codes, weights=w, minlength=n)
    return np.divide(num, den, out=np.zeros(n), where=den > 0)


def a_weekfe(x, D, week, w=None):
    """Intercept after removing week means (s03 decompose, week-FE step, with its dof factor)."""
    w = np.ones(len(x)) if w is None else w
    k = w > 0
    x, D, w, week = x[k], D[k], w[k], week[k]
    wk = pd.factorize(week)[0]
    nw = wk.max() + 1
    mw = wmean_by(wk, x, w, nw)
    cnt = np.bincount(wk, minlength=nw)
    f = np.where(cnt[wk] > 1, cnt[wk] / np.maximum(cnt[wk] - 1, 1), np.nan)
    Lw = ((x - mw[wk]) / mw[wk]) ** 2 * f
    ok = np.isfinite(Lw)
    return fgls_affine(Lw[ok], D[ok], w[ok])[0]


def icc_week(x, week, w=None, w_week=None):
    """s03 definition: var(week means of dev)/var(dev), dev = x/xbar - 1. Under bootstrap weights the
    week means are row-weighted and the across-week variance is weighted by the week multiplicity."""
    w = np.ones(len(x)) if w is None else w
    k = w > 0
    x, w, week = x[k], w[k], week[k]
    xb = np.average(x, weights=w)
    dev = x / xb - 1
    codes, uniq = pd.factorize(week)
    m = wmean_by(codes, dev, w, len(uniq))
    ww = np.ones(len(uniq)) if w_week is None else np.asarray(w_week, float)
    vb = np.sum(ww * (m - np.average(m, weights=ww)) ** 2) / (ww.sum() - 1)
    vt = np.sum(w * (dev - np.average(dev, weights=w)) ** 2) / (w.sum() - 1)
    return float(vb / vt)


def masks(rows):
    okD = rows["D"] >= 5
    return {"RA": okD & (rows["Ct"] >= 5) & np.isfinite(rows["RA"]), "G": okD & np.isfinite(rows["G"])}


def fit_rep(rows, sel):
    out = {"N": {}, "Rbar_true": {}}
    mk = masks(rows)
    for nm in ("RA", "G"):
        m = mk[nm] & sel
        x, D, wk = rows[nm][m], rows["D"][m] if nm == "G" else rows["Ct"][m], rows["week"][m]
        f = fit_stat(x, D)
        out[nm] = {"a": f["a"], "b": f["b"], "xbar": f["xbar"], "a_weekFE": a_weekfe(x, D, wk),
                   "icc": icc_week(x, wk), "D_med": float(np.median(D))}
        out["N"][nm] = int(m.sum())
        out["Rbar_true"][nm] = float(rows["R"][m].mean())
    return out


def sim_panel_rows(cfg, vr, n, rng, meanfield=False):
    blocks = []
    if cfg["design"] == "s08":                                   # bridge: one block, 500 states, 6 weeks
        d0 = np.exp(rng.uniform(np.log(2.0), np.log(150.0), 500)) if not meanfield else np.full(500, 1e6)
        C, Rw = sim_block(n, 6, 1.2, d0, vr, cfg, rng, meanfield)
        blocks.append(block_rows(C, Rw, n, 500, 6, 0))
    else:
        for si, ds in enumerate(DESIGN):
            rbar = float(wl_R(ds["g_nat"], *GI_FLU[cfg["gi"]]))
            C, Rw = sim_block(n, ds["L"], rbar, ds["d0"], vr, cfg, rng, meanfield)
            blocks.append(block_rows(C, Rw, n, NST, ds["L"], si))
    return {k: np.concatenate([b[k] for b in blocks]) for k in blocks[0]}


def done_keys():
    if not ROWS.exists():
        return set()
    out = set()
    for line in ROWS.read_text().splitlines():
        try:
            r = json.loads(line)
            out.add((r["cell"], r["rep"]))
        except Exception:
            pass
    return out


def run_sims():
    done = done_keys()
    for ci, (cid, sname, cfg, vr) in enumerate(CELLS):
        for ch in range(NREP_A // CHUNK):
            reps = range(ch * CHUNK, (ch + 1) * CHUNK)
            if all((cid, r) in done for r in reps):
                continue
            t0 = time.time()
            rng = np.random.default_rng(np.random.SeedSequence([SEED, 9, ci, ch]))
            rows = sim_panel_rows(cfg, vr, CHUNK, rng)
            lines = []
            for j, r in enumerate(reps):
                fr = fit_rep(rows, rows["rep"] == j)
                lines.append(json.dumps({"cell": cid, "setting": sname, "v_R": vr, "rep": r, "chunk": ch,
                                         "seed": [SEED, 9, ci, ch], **fr}))
            with open(ROWS, "a") as f:
                f.write("\n".join(lines) + "\n")
            last = json.loads(lines[-1])
            log(f"{cid} chunk {ch} RA a={last['RA']['a']:.4f} G a={last['G']['a']:.4f} "
                f"N={last['N']} ({time.time() - t0:.1f} s)")


def meanfield(cfg, vr, ci):
    rng = np.random.default_rng(np.random.SeedSequence([SEED, 99, ci]))
    rows = sim_panel_rows(cfg, vr, 100, rng, meanfield=True)
    mk = masks(rows)
    out = {}
    for nm in ("RA", "G"):
        x = rows[nm][mk[nm]]
        out[nm] = {"a_inf": float(np.mean(((x - x.mean()) / x.mean()) ** 2)), "xbar_inf": float(x.mean()),
                   "Rbar_true": float(rows["R"][mk[nm]].mean())}
    return out


# ============================================================================ empirical joint bootstrap
def coverage_adjusted(e, stat):
    """Copy of s04_reporting.adjusted (coverage-proxy adjusted statistic)."""
    h = {j: e[c] / 100.0 for j, c in (("n", "coverage_next"), ("0", "coverage"), ("1", "coverage_prev1"),
                                      ("2", "coverage_prev2"))}
    ok = np.all([np.isfinite(v) & (v > 0) for v in h.values()], axis=0)
    if stat == "G":
        ds = W_ORIG[0] * e.C_t / h["0"] + W_ORIG[1] * e.C_prev1 / h["1"] + W_ORIG[2] * e.C_prev2 / h["2"]
        xs = (e.C_next / h["n"]) / ds
    else:
        xs = pd.Series(wl_R((e.C_next / h["n"] + 0.5) / (e.C_t / h["0"] + 0.5)), index=e.index)
    return xs.to_numpy(), np.asarray(ok)


def joint_bootstrap(panel):
    """B_INV joint two-way draws for RA_rule4 and G_rule4 (shared state and week multiplicities)."""
    frames = {"RA": spec_frame(panel, "RA_rule4"), "G": spec_frame(panel, "G_rule4")}
    states = sorted(set(frames["RA"].state) | set(frames["G"].state))
    weeks = sorted(set(frames["RA"].week_end) | set(frames["G"].week_end))
    sidx = {s: i for i, s in enumerate(states)}
    widx = {w: i for i, w in enumerate(weeks)}
    codes = {nm: (e.state.map(sidx).to_numpy(), e.week_end.map(widx).to_numpy()) for nm, e in frames.items()}
    adj = {nm: coverage_adjusted(e, nm) for nm, e in frames.items()}

    def stats(w, wv=None):
        rec = {}
        for nm, e in frames.items():
            x, D, wk = e.x.to_numpy(), e.D.to_numpy(), e.week_end.to_numpy()
            ww = np.ones(len(e)) if w is None else w[nm]
            f = fit_stat(x, D, ww)
            rec[f"{nm}_a"], rec[f"{nm}_b"], rec[f"{nm}_xbar"] = f["a"], f["b"], f["xbar"]
            rec[f"{nm}_ba"] = f["ba"]
            rec[f"{nm}_a_weekFE"] = a_weekfe(x, D, wk, ww)
            if wv is None:
                rec[f"{nm}_icc"] = icc_week(x, wk)
            else:                                   # week multiplicities in order of first appearance
                kk = ww > 0
                order = pd.factorize(codes[nm][1][kk])[1]
                rec[f"{nm}_icc"] = icc_week(x, codes[nm][1], ww, wv[order])
            xs, ok = adj[nm]
            wo = ww[ok]
            rec[f"{nm}_a_cov"] = fit_stat(x[ok], D[ok], wo)["a"]
            rec[f"{nm}_a_adj"] = fit_stat(xs[ok], D[ok], wo)["a"]
        return rec

    point = stats(None)
    if DRAWS.exists():
        z = np.load(DRAWS, allow_pickle=False)
        dr = pd.DataFrame({k: z[k] for k in z.files})
    else:
        rng = np.random.default_rng(SEED + 901)
        out = []
        t0 = time.time()
        for b in range(B_INV):
            ws = np.bincount(rng.integers(0, len(states), len(states)), minlength=len(states)).astype(float)
            wv = np.bincount(rng.integers(0, len(weeks), len(weeks)), minlength=len(weeks)).astype(float)
            out.append(stats({nm: ws[c[0]] * wv[c[1]] for nm, c in codes.items()}, wv))
            if (b + 1) % 500 == 0:
                log(f"joint bootstrap {b + 1}/{B_INV} ({time.time() - t0:.0f} s)")
        dr = pd.DataFrame(out)
        np.savez_compressed(DRAWS, **{k: dr[k].to_numpy() for k in dr.columns})
    return frames, point, dr


# ============================================================================ A
def part_A(panel):
    P = Part("inversion", "s09_inversion.py", INPUTS, "simulated",
             {"simulation": SEED_SIM, "meanfield": "SeedSequence([SEED, 99, cell_index])",
              "joint_bootstrap": SEED + 901, "calibration_draws": SEED + 902})
    rows = [json.loads(l) for l in ROWS.read_text().splitlines() if l.strip()]
    df = pd.json_normalize(rows)
    P.put("design", {"seasons": [{k: (v if k != "d0" else [float(np.min(v)), float(np.median(v)), float(np.max(v))])
                                  for k, v in d.items()} for d in DESIGN],
                     "d0_note": "d0 summary = [min, median, max] of observed first-window-week admissions by state",
                     "Rbar_season_GImid": [float(wl_R(d["g_nat"], *GI_FLU["mid"])) for d in DESIGN],
                     "base": BASE, "settings": SETTINGS, "labels": LABEL, "v_grid": VGRID, "nrep": NREP_A,
                     "pre_weeks": PRE_W, "states": NST, "cells": [c[0] for c in CELLS]})
    frames, point, dr = joint_bootstrap(panel)
    emp = {nm: {"a": point[f"{nm}_a"], "b": point[f"{nm}_b"], "xbar": point[f"{nm}_xbar"],
                "a_ci": pct(dr[f"{nm}_a"]), "a_sd": float(dr[f"{nm}_a"].std(ddof=1)),
                "icc": point[f"{nm}_icc"], "share_week": 1 - point[f"{nm}_a_weekFE"] / point[f"{nm}_a"],
                "N": int(len(frames[nm]))} for nm in ("RA", "G")}
    P.put("empirical", emp)

    cells = {}
    for ci, (cid, sname, cfg, vr) in enumerate(CELLS):
        sub = df[df.cell == cid]
        if len(sub) < NREP_A:
            log("incomplete cell", cid, len(sub))
            continue
        mf = meanfield(cfg, vr, ci)
        rec = {"setting": sname, "v_R": vr, "n": int(len(sub))}
        for nm in ("RA", "G"):
            a, b = sub[f"{nm}.a"].to_numpy(), sub[f"{nm}.b"].to_numpy()
            rb = float(sub[f"Rbar_true.{nm}"].mean())
            a_nom = vr / rb ** 2
            rec[nm] = {"a_mean": float(a.mean()), "a_se": float(a.std(ddof=1) / np.sqrt(len(a))),
                       "a_sd": float(a.std(ddof=1)), "b_mean": float(b.mean()),
                       "ba_pooled": float(b.mean() / a.mean()), "ba_rep_q": pct(b / a, (2.5, 50, 97.5)),
                       "xbar": float(sub[f"{nm}.xbar"].mean()), "Rbar_true": rb,
                       "kappa_fit": float(a.mean() / a_nom), "kappa_fit_se": float(a.std(ddof=1) / np.sqrt(len(a)) / a_nom),
                       "a_inf": mf[nm]["a_inf"], "kappa_inf": mf[nm]["a_inf"] / (vr / mf[nm]["Rbar_true"] ** 2),
                       "icc": float(sub[f"{nm}.icc"].mean()),
                       "share_week": float(np.mean(1 - sub[f"{nm}.a_weekFE"] / sub[f"{nm}.a"])),
                       "N": float(sub[f"N.{nm}"].mean()), "D_med": float(sub[f"{nm}.D_med"].mean())}
        cells[cid] = rec
        log(f"{cid}: RA a={rec['RA']['a_mean']:.4f}+-{rec['RA']['a_se']:.4f} kfit={rec['RA']['kappa_fit']:.3f} "
            f"kinf={rec['RA']['kappa_inf']:.3f} icc={rec['RA']['icc']:.2f} N={rec['RA']['N']:.0f} | "
            f"G a={rec['G']['a_mean']:.4f} kfit={rec['G']['kappa_fit']:.3f} kinf={rec['G']['kappa_inf']:.3f} "
            f"icc={rec['G']['icc']:.2f} ba RA={rec['RA']['ba_pooled']:.1f} G={rec['G']['ba_pooled']:.1f}")
    P.put("cells", cells)

    # ---------------------------------------------------------- inversion per setting
    rng = np.random.default_rng(SEED + 902)
    sens, vdraws = {}, {}
    for sname in SETTINGS:
        cs = [cells[f"{sname}|v{v}"] for v in VGRID if f"{sname}|v{v}" in cells]
        if len(cs) < len(VGRID):
            continue
        v = np.array([c["v_R"] for c in cs])
        rec = {"label": LABEL[sname]}
        for nm in ("RA", "G"):
            y = np.array([c[nm]["a_mean"] for c in cs]); se = np.array([c[nm]["a_se"] for c in cs])
            X = np.column_stack([np.ones_like(v), v]); Wt = 1 / se ** 2
            cov = np.linalg.inv(X.T @ (X * Wt[:, None]))
            coef = cov @ (X.T @ (Wt * y))
            resid = y - X @ coef
            lin_dev = float(np.max(np.abs(resid) / se))
            a_e = emp[nm]["a"]
            v_hat = (a_e - coef[0]) / coef[1]
            cd = rng.multivariate_normal(coef, cov, len(dr))
            v_b = (dr[f"{nm}_a"].to_numpy() - cd[:, 0]) / cd[:, 1]
            v_mc = (a_e - cd[:, 0]) / cd[:, 1]
            # linear-kappa version, kappa from the cell nearest to the implied v_R
            near = cs[int(np.argmin(np.abs(v - v_hat)))]
            rb = near[nm]["Rbar_true"]
            v_lin = a_e * rb ** 2 / near[nm]["kappa_fit"]
            v_lin_inf = a_e * rb ** 2 / near[nm]["kappa_inf"]
            kap = [c[nm]["kappa_fit"] for c in cs]
            rec[nm] = {"c0": float(coef[0]), "c1": float(coef[1]), "cov": cov.tolist(), "max_abs_resid_z": lin_dev,
                       "v_implied": float(v_hat), "v_ci95": pct(v_b), "v_ci95_mc_only": pct(v_mc),
                       "v_ci95_boot_only": pct((dr[f"{nm}_a"].to_numpy() - coef[0]) / coef[1]),
                       "kappa_fit_grid": kap, "kappa_fit_mean": float(np.mean(kap)),
                       "kappa_inf_grid": [c[nm]["kappa_inf"] for c in cs],
                       "kappa_at_nearest": near[nm]["kappa_fit"], "kappa_inf_at_nearest": near[nm]["kappa_inf"],
                       "Rbar_true": rb, "v_linear_kappa_fit": float(v_lin), "v_linear_kappa_inf": float(v_lin_inf),
                       "sim_icc_at_nearest": near[nm]["icc"], "sim_share_week_at_nearest": near[nm]["share_week"],
                       "sim_a_sd_at_nearest": near[nm]["a_sd"], "sim_xbar_at_nearest": near[nm]["xbar"],
                       "sim_N_at_nearest": near[nm]["N"]}
            vdraws[f"{sname}|{nm}"] = v_b
        d = vdraws[f"{sname}|RA"] - vdraws[f"{sname}|G"]
        rat = vdraws[f"{sname}|RA"] / vdraws[f"{sname}|G"]
        rec["diff_RA_minus_G"] = {"value": rec["RA"]["v_implied"] - rec["G"]["v_implied"], "ci95": pct(d)}
        rec["ratio_RA_over_G"] = {"value": rec["RA"]["v_implied"] / rec["G"]["v_implied"], "ci95": pct(rat)}
        rec["agree"] = bool(pct(d)[0] < 0 < pct(d)[1])
        sens[sname] = rec
        log(f"inversion {sname}: v_RA={rec['RA']['v_implied']:.4f} {np.round(rec['RA']['v_ci95'], 4)} "
            f"v_G={rec['G']['v_implied']:.4f} {np.round(rec['G']['v_ci95'], 4)} "
            f"diff CI {np.round(rec['diff_RA_minus_G']['ci95'], 4)} kRA={rec['RA']['kappa_at_nearest']:.3f} "
            f"kG={rec['G']['kappa_at_nearest']:.3f}")
    P.put("sensitivity", sens)
    np.savez_compressed(PARTS / "s09_vR_draws.npz", **{k.replace("|", "__"): v for k, v in vdraws.items()})

    m0 = sens["M0"]
    for nm in ("RA", "G"):
        r = m0[nm]
        P.num(f"M0_{nm}_vR_implied", r["v_implied"], ci=r["v_ci95"], mode="simulated",
              note="calibration-line inversion, data-matched M0; interval combines joint two-way bootstrap of "
                   "a (B=2000) with Monte Carlo uncertainty of the calibration line",
              ci95_mc_only=r["v_ci95_mc_only"], ci95_bootstrap_only=r["v_ci95_boot_only"])
        P.num(f"M0_{nm}_kappa_fit", r["kappa_at_nearest"], mode="simulated",
              note="E[a_hat]/(v_R/Rbar_true^2) at the grid cell nearest the implied v_R",
              grid=r["kappa_fit_grid"], Rbar_true=r["Rbar_true"])
        P.num(f"M0_{nm}_kappa_inf", r["kappa_inf_at_nearest"], mode="simulated",
              note="noise-free (mean-field) limit a_inf/(v_R/Rbar_true^2), same design")
        P.num(f"M0_{nm}_vR_linear_kappa", r["v_linear_kappa_fit"], mode="simulated",
              note="a_emp * Rbar_true^2 / kappa_fit (one-line formula; ignores the calibration-line offset c0)")
        P.num(f"M0_{nm}_calib_c0", r["c0"], mode="simulated"); P.num(f"M0_{nm}_calib_c1", r["c1"], mode="simulated")
        P.num(f"M0_{nm}_sim_icc", r["sim_icc_at_nearest"], mode="simulated", empirical=emp[nm]["icc"])
        P.num(f"M0_{nm}_sim_share_week", r["sim_share_week_at_nearest"], mode="simulated",
              empirical=emp[nm]["share_week"])
        P.num(f"M0_{nm}_sim_a_sd", r["sim_a_sd_at_nearest"], mode="simulated",
              note="replicate SD of a_hat across simulated panels", empirical_twoway_boot_sd=emp[nm]["a_sd"])
        P.num(f"M0_{nm}_sim_xbar", r["sim_xbar_at_nearest"], mode="simulated", empirical=emp[nm]["xbar"])
        P.num(f"M0_{nm}_sim_N", r["sim_N_at_nearest"], mode="simulated", empirical=emp[nm]["N"])
        P.num(f"emp_{nm}_a_joint", emp[nm]["a"], ci=emp[nm]["a_ci"], mode="measured", seed=SEED + 901, B=B_INV)
    P.num("M0_diff_vR_RA_minus_G", m0["diff_RA_minus_G"]["value"], ci=m0["diff_RA_minus_G"]["ci95"], mode="simulated")
    P.num("M0_ratio_vR_RA_over_G", m0["ratio_RA_over_G"]["value"], ci=m0["ratio_RA_over_G"]["ci95"], mode="simulated")
    vals = [sens[s][nm]["v_implied"] for s in sens if s != "bridge_s08" for nm in ("RA", "G")]
    P.num("vR_implied_range_across_settings", [float(min(vals)), float(max(vals))], mode="simulated",
          note="min-max of point v_R^implied over both statistics and all data-matched settings (bridge excluded)")
    P.num("n_settings_agree", int(sum(sens[s]["agree"] for s in sens if s != "bridge_s08")), mode="simulated",
          note=f"settings (of {len([s for s in sens if s != 'bridge_s08'])}) where the 95% interval of "
               "v_RA - v_G covers 0")
    for nm in ("RA", "G"):
        bands = [s for s in sens if s != "bridge_s08" and BAND[0] <= sens[s][nm]["v_implied"] <= BAND[1]]
        P.num(f"n_settings_{nm}_in_band", len(bands), mode="simulated", settings=bands,
              note="settings whose point v_R^implied lies in [0.04, 0.08]")
    br = sens.get("bridge_s08")
    if br:
        c = [cells[f"bridge_s08|v{v}"] for v in VGRID]
        P.num("bridge_kappa_inf_G", float(np.median([x["G"]["a_inf"] / (x["v_R"] / 1.2 ** 2) for x in c])),
              mode="simulated", note="s08 design re-run with the s09 simulator; s08 reported 2.6812")
        P.num("bridge_kappa_inf_RA", float(np.median([x["RA"]["a_inf"] / (x["v_R"] / 1.2 ** 2) for x in c])),
              mode="simulated", note="s08 reported 0.48448")

    # ---------------------------------------------------------- calibrated theoretical b/a
    ba, ba_reps = {}, {}
    for nm in ("RA", "G"):
        pool, reps = [], []
        for cid, c in cells.items():
            if c["setting"].startswith("ba_"):
                pool.append(c[nm]["ba_pooled"])
                sub = df[df.cell == cid]
                reps.append((sub[f"{nm}.b"] / sub[f"{nm}.a"]).to_numpy())
        reps = np.concatenate(reps)
        ba_reps[nm] = reps
        q = pct(reps, (2.5, 97.5))
        emp_ba = point[f"{nm}_ba"]
        ba[nm] = {"pooled_range": [float(min(pool)), float(max(pool))], "replicate_q95": q,
                  "M0_k0.35_r0.01_v0.04": cells["M0|v0.04"][nm]["ba_pooled"], "empirical": emp_ba,
                  "empirical_in_pooled_range": bool(min(pool) <= emp_ba <= max(pool)),
                  "empirical_in_replicate_range": bool(q[0] <= emp_ba <= q[1]),
                  "empirical_rep_percentile": float(np.mean(reps <= emp_ba))}
        P.num(f"ba_calibrated_{nm}_range", ba[nm]["pooled_range"], mode="simulated",
              note="E[b_hat]/E[a_hat] of the data-matched simulation over v_R {0.04,0.08} x k {0.2,0.5} x "
                   "rho_h {0.0056,0.021} (s03 theory scenarios)")
        P.num(f"ba_calibrated_{nm}_replicate_q95", q, mode="simulated",
              note="2.5-97.5% of per-panel b/a over the same 8 cells (sampling range of one panel of this size)")
        P.num(f"ba_empirical_{nm}_percentile_in_sim", ba[nm]["empirical_rep_percentile"], mode="simulated")
    P.put("ba_calibrated", ba)
    np.savez_compressed(PARTS / "s09_ba_reps.npz", **ba_reps)
    log("b/a calibrated", json.dumps(ba))


# ============================================================================ B
def part_B(panel):
    P = Part("multiscale", "s09_inversion.py", INPUTS, "measured", {"joint_bootstrap": SEED + 901,
                                                                     "subsampling": SEED + 911})
    frames, point, dr = joint_bootstrap(panel)
    for nm, spec in (("RA", "RA_rule4"), ("G", "G_rule4")):
        a, aw = point[f"{nm}_a"], point[f"{nm}_a_weekFE"]
        sh = 1 - dr[f"{nm}_a_weekFE"] / dr[f"{nm}_a"]
        P.num(f"{spec}_a_total", a, ci=pct(dr[f"{nm}_a"]), B=B_INV)
        P.num(f"{spec}_a_weekFE", aw, ci=pct(dr[f"{nm}_a_weekFE"]), B=B_INV)
        P.num(f"{spec}_comp_week", a - aw, ci=pct(dr[f"{nm}_a"] - dr[f"{nm}_a_weekFE"]), B=B_INV)
        P.num(f"{spec}_share_week", 1 - aw / a, ci=pct(sh), B=B_INV,
              note="share of the intercept removed by week fixed effects (national synchronous component)")
        P.num(f"{spec}_icc_week", point[f"{nm}_icc"], ci=pct(dr[f"{nm}_icc"]), B=B_INV)
        sds = {}
        for sch in ("state", "week", "twoway"):
            z = np.load(PARTS / f"s02_draws_{spec}_{sch}.npy")
            sds[sch] = float(np.nanstd(z[:, 0], ddof=1))
        P.num(f"{spec}_a_sd_state", sds["state"], note="SD of a over s02 bootstrap draws (B=2000)")
        P.num(f"{spec}_a_sd_week", sds["week"])
        P.num(f"{spec}_a_sd_twoway", sds["twoway"])
        P.num(f"{spec}_a_sd_ratio_week_over_state", sds["week"] / sds["state"])
        log(nm, "share_week", round(1 - aw / a, 4), pct(sh), "icc", point[f"{nm}_icc"], sds)

        # J-state pooling of the weekly growth deviation (observed scale)
        e = frames[nm]
        dev = (e.x / e.x.mean() - 1).to_numpy()
        wk = pd.factorize(e.week_end)[0]
        st = pd.factorize(e.state)[0]
        nw = wk.max() + 1
        n_w = np.bincount(wk)
        m_w = np.bincount(wk, weights=dev) / n_w
        V1 = float(np.mean((dev - dev.mean()) ** 2))
        within = float(np.sum((dev - m_w[wk]) ** 2) / (len(dev) - nw))
        s2w = float(max(np.var(m_w, ddof=1) - np.mean(within / n_w), 0.0))
        rng = np.random.default_rng(SEED + 911 + (nm == "G"))
        Js = [1, 2, 3, 5, 8, 12, 20, 30, 40, 51]
        emp_curve = {}
        M = np.full((st.max() + 1, nw), np.nan)
        M[st, wk] = dev
        for J in Js:
            vals = []
            for _ in range(300):
                pick = rng.choice(M.shape[0], J, replace=False)
                mm = np.nanmean(M[pick], axis=0)
                mm = mm[np.isfinite(mm)]
                vals.append(np.mean(mm ** 2))
            emp_curve[J] = {"mean": float(np.mean(vals)), "q": pct(vals, (2.5, 97.5))}
        model = {J: s2w + (V1 - s2w) / J for J in Js}
        J90 = int(np.ceil((V1 - s2w) / (0.1 * s2w)))
        P.put(f"Jpool|{spec}", {"J": Js, "empirical": emp_curve, "model": model, "V1": V1, "s2_week": s2w,
                                "within": within, "note": "dev = x/xbar - 1; empirical = mean over weeks of the "
                                "squared J-state mean deviation, 300 random state subsets per J (band = 2.5-97.5% "
                                "over subsets); model = one-way random effects s2_week + (V1 - s2_week)/J"})
        P.num(f"{spec}_pool_V1", V1, note="variance of a single-state weekly deviation (observed scale)")
        P.num(f"{spec}_pool_floor_s2week", s2w, note="between-week variance component (ANOVA), the limit of J-state pooling")
        P.num(f"{spec}_pool_floor_share", s2w / V1)
        P.num(f"{spec}_pool_J_within10pct_floor", J90, mode="analytic",
              note="smallest J with (V1 - s2_week)/J <= 0.1 s2_week")
        for J in (5, 10, 51):
            P.num(f"{spec}_pool_VJ_over_V1_J{J}", (s2w + (V1 - s2w) / J) / V1, mode="analytic")
        log(nm, "pool V1", V1, "s2w", s2w, "share", s2w / V1, "J90", J90)


# ============================================================================ C
def wiqr(v, w):
    o = np.argsort(v)
    cw = np.cumsum(w[o]) / w.sum()
    return float(v[o][np.searchsorted(cw, 0.75)] - v[o][np.searchsorted(cw, 0.25)])


def mayoct(p, y):
    s = p[(p.week_end >= f"{y}-05-01") & (p.week_end <= f"{y}-10-31")].reset_index(drop=True)
    s = s.assign(lg=np.log((s.C_next + 0.5) / (s.D_t + 0.5)))
    return s


def part_C(panel):
    P = Part("reporting_regime", "s09_inversion.py", INPUTS, "measured",
             {"joint_bootstrap": SEED + 901, "season_ratio_bootstrap": SEED + 921})
    frames, point, dr = joint_bootstrap(panel)
    inv = json.loads((PARTS / "inversion.json").read_text())["numbers"] if (PARTS / "inversion.json").exists() else {}
    s04 = json.loads((PARTS / "s04_reporting.json").read_text())
    for nm, spec in (("RA", "RA_rule4"), ("G", "G_rule4")):
        blk = s04["blocks"][f"delta|{spec}"]
        k2 = blk["kappa2"]
        xb = point[f"{nm}_xbar"]
        a = point[f"{nm}_a"]
        for vr_lab, vr in (("vR0.04", 0.04), ("vR_implied", inv.get(f"M0_{nm}_vR_implied", {}).get("value"))):
            if vr is None:
                continue
            margin = k2 * 0.1 ** 2 * (1 + vr / xb ** 2)
            P.num(f"{spec}_margin_share_of_a_{vr_lab}", margin / a, ci=pct(margin / dr[f"{nm}_a"]), mode="analytic",
                  margin=margin, note="TOST margin k2 tau0^2 (1 + v_R/xbar^2), tau0 = 0.1, xbar = sample mean of the "
                                      "statistic (s04 convention), divided by a; interval from the joint bootstrap of a")
        da = point[f"{nm}_a_cov"] - point[f"{nm}_a_adj"]
        sh = (dr[f"{nm}_a_cov"] - dr[f"{nm}_a_adj"]) / dr[f"{nm}_a_cov"]
        q = pct(sh, (2.5, 95, 97.5))
        P.num(f"{spec}_delta_a_share", da / point[f"{nm}_a_cov"], ci=[q[0], q[2]], upper95_one_sided=q[1], B=B_INV,
              note="(a_raw - a_coverage_adjusted)/a_raw on rows with valid coverage; captures only what the "
                   "coverage proxy explains")
        log(spec, "delta share", da / point[f"{nm}_a_cov"], q)

    # ---------------------------------------------------- May-Oct 2024 (voluntary) vs 2023 / 2025 (mandatory)
    def year_sets(p):
        return {y: mayoct(p, y) for y in (2023, 2024, 2025)}

    ys = year_sets(panel)
    states = sorted(panel.state.unique())
    sidx = {s: i for i, s in enumerate(states)}
    rng = np.random.default_rng(SEED + 921)
    draws = []
    for _ in range(2000):
        ws = np.bincount(rng.integers(0, len(states), len(states)), minlength=len(states)).astype(float)
        rec = {}
        for y, s in ys.items():
            wcode = pd.factorize(s.week_end)[0]
            nwk = wcode.max() + 1
            wv = np.bincount(rng.integers(0, nwk, nwk), minlength=nwk).astype(float)
            w = ws[s.state.map(sidx).to_numpy()] * wv[wcode]
            rec[f"iqr{y}"] = wiqr(s.lg.to_numpy(), w)
            rec[f"sd{y}"] = float(np.sqrt(np.cov(s.G, aweights=w)))
        draws.append(rec)
    dd = pd.DataFrame(draws)
    pt = {y: wiqr(s.lg.to_numpy(), np.ones(len(s))) for y, s in ys.items()}
    sdp = {y: float(s.G.std(ddof=1)) for y, s in ys.items()}
    pooled_ref = wiqr(pd.concat([ys[2023], ys[2025]]).lg.to_numpy(), np.ones(len(ys[2023]) + len(ys[2025])))
    for ref in (2023, 2025):
        r = pt[2024] / pt[ref]
        P.num(f"iqr_ratio_2024_over_{ref}", r, ci=pct(dd["iqr2024"] / dd[f"iqr{ref}"]), B=2000,
              note="IQR of log((C_next+0.5)/(D_t+0.5)), May-Oct; state draw shared across years, weeks within year")
        P.num(f"sd_ratio_2024_over_{ref}", sdp[2024] / sdp[ref], ci=pct(dd["sd2024"] / dd[f"sd{ref}"]), B=2000)
    P.num("iqr_ratio_2025_over_2023_placebo", pt[2025] / pt[2023], ci=pct(dd["iqr2025"] / dd["iqr2023"]), B=2000)
    P.num("iqr_ratio_2024_over_pooled_2023_2025", pt[2024] / pooled_ref, mode="measured")
    for y in ys:
        P.num(f"mayoct{y}_D_median", float(ys[y].D_t.median()))
        P.num(f"mayoct{y}_phat_median", float(ys[y].coverage.median()))

    # scale-matched: reweight each year to the pooled 2023+2025 distribution of log D (10 quantile bins)
    ref = np.log(pd.concat([ys[2023], ys[2025]]).D_t.to_numpy())
    edges = np.quantile(ref, np.linspace(0, 1, 11)); edges[0], edges[-1] = -np.inf, np.inf
    pref = np.histogram(ref, edges)[0] / len(ref)
    dm = {}
    for y, s in ys.items():
        b = np.clip(np.searchsorted(edges, np.log(s.D_t.to_numpy()), side="right") - 1, 0, 9)
        py = np.bincount(b, minlength=10) / len(b)
        w = pref[b] / np.where(py[b] > 0, py[b], np.nan)
        w = np.nan_to_num(w)
        dm[y] = wiqr(s.lg.to_numpy(), w)
    for refy in (2023, 2025):
        P.num(f"iqr_ratio_2024_over_{refy}_Dmatched", dm[2024] / dm[refy], mode="measured",
              note="each year reweighted to the pooled 2023+2025 decile distribution of log D_t")
    # month-matched
    mm = {}
    for mth in range(5, 11):
        v = {y: wiqr(s[s.week_end.dt.month == mth].lg.to_numpy(), np.ones(int((s.week_end.dt.month == mth).sum())))
             for y, s in ys.items()}
        mm[mth] = {"2024/2023": v[2024] / v[2023], "2024/2025": v[2024] / v[2025], "2025/2023": v[2025] / v[2023]}
    P.put("month_matched", mm)
    r24 = [mm[m]["2024/2023"] for m in mm] + [mm[m]["2024/2025"] for m in mm]
    P.num("iqr_ratio_month_matched_median", float(np.median(r24)), range=[float(min(r24)), float(max(r24))],
          note="median over May-Oct months of the 2024/2023 and 2024/2025 monthly IQR ratios (12 ratios)")
    P.num("iqr_ratio_month_matched_n_gt1", int(sum(x > 1 for x in r24)), note="of 12 monthly ratios")
    # other pathogens, same months
    for pg in ("covid", "rsv"):
        pp = build_panel(load_raw(pg))
        yp = year_sets(pp)
        v = {y: (wiqr(s.lg.to_numpy(), np.ones(len(s))) if len(s) >= 30 else float("nan")) for y, s in yp.items()}
        nn = [int(len(yp[y])) for y in (2023, 2024, 2025)]
        P.num(f"{pg}_iqr_ratio_2024_over_2023", v[2024] / v[2023], N_2023_2024_2025=nn,
              note=f"{pg.upper()} admissions, May-Oct, same statistic; NaN if a year has < 30 rows with D_t >= 5")
        P.num(f"{pg}_iqr_ratio_2024_over_2025", v[2024] / v[2025])
        P.num(f"{pg}_mayoct_D_median_2023_2024_2025", [float(yp[y].D_t.median()) if len(yp[y]) else float("nan")
                                                         for y in (2023, 2024, 2025)])
        log(pg, v, nn)
    log("iqr", pt, "Dmatched", dm, "month", mm)

    # ---------------------------------------------------- participation under the mandate
    cols = ["Week Ending Date", "Geographic aggregation", "Percent Hospitals Reporting Influenza Admissions",
            "Number Hospitals Reporting Influenza Admissions"]
    from cc_common import RAW_CSV, EXCLUDE, P1_END, P3_START
    raw = pd.read_csv(RAW_CSV, usecols=cols)
    g = raw["Geographic aggregation"].astype(str)
    raw = raw[(g.str.len() == 2) & ~g.isin(EXCLUDE)].copy()
    raw["week_end"] = pd.to_datetime(raw["Week Ending Date"])
    raw = raw[raw.week_end >= "2022-10-01"]
    pc, nh = raw[cols[2]], raw[cols[3]]
    raw["n_total"] = np.where(pc > 0, nh / (pc / 100.0), np.nan)
    raw["regime"] = np.where(raw.week_end <= P1_END, 1, np.where(raw.week_end < P3_START, 2, 3))
    part = {}
    for k in (1, 2, 3):
        s = raw[raw.regime == k]
        c = s[cols[2]].dropna()
        natw = s.groupby("week_end").apply(lambda q: q[cols[3]].sum() / q.n_total.sum() * 100, include_groups=False)
        part[k] = {"state_week_mean": float(c.mean()), "state_week_median": float(c.median()),
                   "share_ge95": float(np.mean(c >= 95)), "share_ge90": float(np.mean(c >= 90)),
                   "national_hosp_weighted_median": float(natw.median()), "national_q": pct(natw, (5, 95)),
                   "N_state_weeks": int(len(c))}
        for kk, vv in part[k].items():
            P.num(f"regime{k}_{kk}", vv, note="from 2022-10-01; raw NHSN state-weeks (not filtered on D_t)")
    # inside the analysis windows
    e = frames["RA"]
    P.num("rule4_windows_coverage_median", float(e.coverage.median()))
    P.num("rule4_windows_coverage_share_ge95", float(np.mean(e.coverage >= 95)))
    P.num("rule4_windows_coverage_share_ge90", float(np.mean(e.coverage >= 90)))
    P.put("participation", part)
    log("participation", part)


# ============================================================================ D
def part_D(panel):
    P = Part("returns", "s09_inversion.py", INPUTS, "analytic", {"s02_draws": "s02 two-way draws (seed 20261024/20261124)"})
    R = json.loads((PARTS / "s02_inference.json").read_text())["numbers"]
    for nm, spec, col in (("RA", "RA_rule4", "C_t"), ("G", "G_rule4", "D_t")):
        mx = R[f"{spec}_ba"]["value"]
        z = np.load(PARTS / f"s02_draws_{spec}_twoway.npy")[:, 2]
        z = z[np.isfinite(z) & (z > 0)]
        P.num(f"{spec}_m_cross", mx, ci=R[f"{spec}_ba"]["ci95"], mode="measured", source=f"s02_inference.{spec}_ba")
        for f_, lab in ((0.5, "half"), (1, "1x"), (2, "2x"), (4, "4x")):
            m = f_ * mx
            P.num(f"precision_fraction_at_{lab}_mcross", m / (m + mx),
                  note="(1/V(m))/(1/a) = m/(m + m_x); identical for any m_x")
            P.num(f"variance_cut_doubling_from_{lab}_mcross", mx / (2 * (m + mx)),
                  note="[V(m) - V(2m)]/V(m) with V(m) = a(1 + m_x/m)")
        for c in (10, 100):
            ms = np.sqrt(c * mx)
            P.num(f"{spec}_m_star_c{c}", ms, ci=pct(np.sqrt(c * z)), note=f"m* = sqrt((c0/c1) m_x), c0/c1 = {c}; "
                  "interval from the two-way bootstrap draws of b/a")
        e = spec_frame(panel, spec)
        P.num(f"{spec}_share_rows_above_mcross", float(np.mean(e.D >= mx)), mode="measured",
              note=f"share of window state-weeks with {col} >= m_x")
        P.num(f"{spec}_share_rows_above_2mcross", float(np.mean(e.D >= 2 * mx)), mode="measured")
    P.num("m_for_90pct_precision_over_mcross", 9.0, note="m/(m+m_x) = 0.9 at m = 9 m_x")


# ============================================================================ E
def part_E(panel):
    P = Part("parag", "s09_inversion.py", INPUTS + ["Parag et al. 2022 (deepread), eqs. 3-5"], "analytic", None)
    inv = json.loads((PARTS / "inversion.json").read_text())["numbers"]
    R = json.loads((PARTS / "s02_inference.json").read_text())["numbers"]
    vr = inv["M0_RA_vR_implied"]["value"]
    rb = inv["M0_RA_kappa_fit"]["Rbar_true"]
    rows = []
    for k in (0.2, 0.5):
        for rh in (0.0056, 0.021):
            for v in (vr, 0.04, 0.08):
                A = rb + rh * (rb ** 2 + v) / k
                ms = {f"m_ratio{r}": (r * rb - A) / v for r in (1.5, 2, 5, 10)}
                rows.append({"k": k, "rho_h": rh, "v_R": v, "A": A, **ms})
    P.put("latent_divergence", {"rows": rows, "Rbar": rb,
                                "formula": "ratio(m) = (v_R + A/m)/(Rbar/m) = (v_R m + A)/Rbar; Poisson renewal FI "
                                           "of R_t (Parag et al. 2022 eq. 3, Lambda/R) gives Var = R/m with m "
                                           "the expected count of the generation; reporting (eq. 5, theta) "
                                           "rescales m to theta m and adds no floor"})
    at = [r for r in rows if r["v_R"] == vr]
    for rr in (1.5, 2, 5, 10):
        vals = [r[f"m_ratio{rr}"] for r in at]
        P.num(f"latent_m_ratio{rr}", float(np.median(vals)), range=[float(min(vals)), float(max(vals))],
              note=f"m at which (v_R + A/m) is {rr}x the Poisson-only variance Rbar/m; v_R = implied (M0, R_A), "
                   "range over k {0.2,0.5} x rho_h {0.0056,0.021}")
    # apparent R_A scale: Poisson delta-method coefficient versus the fitted a + b/m
    e = spec_frame(panel, "RA_rule4")
    m, s = GI_FLU["mid"]
    shape, scale = (m / s) ** 2, s ** 2 / m
    lg = np.log((e.C_next + 0.5) / (e.C_t + 0.5))
    kap = shape * (scale / 7) / (1 + scale * lg / 7)
    bP = float(np.mean(kap ** 2 * (1 + e.C_t / np.maximum(e.C_next, 0.5))))
    a, b = R["RA_rule4_a"]["value"], R["RA_rule4_b"]["value"]
    P.num("RA_poisson_b", bP, note="mean over window rows of (dlnR_A/dlng)^2 (1 + C_t/C_next): Poisson-only "
                                    "relative variance of R_A is b_P/C_t")
    P.num("RA_fitted_b_over_poisson_b", b / bP, note="overdispersion of the 1/m term relative to Poisson")
    for rr in (1.5, 2, 5, 10):
        P.num(f"RA_apparent_m_ratio{rr}", max((rr * bP - b) / a, 0.0),
              note=f"C_t at which a + b/m equals {rr} x b_P/m (0 = already exceeded at every m)")
    curve_m = np.geomspace(2, 5000, 120)
    P.put("curves", {"m": curve_m, "poisson": bP / curve_m, "model": a + b / curve_m, "a": a, "b": b, "bP": bP,
                     "latent_v": vr, "Rbar": rb})
    log("parag", {"bP": bP, "b/bP": b / bP, "rows": at})


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    panel = load_panel() if what != "sim" else None
    if what in ("sim", "all"):
        run_sims()
        log(f"sim stage {time.time() - t0:.0f} s")
    if panel is None:
        panel = load_panel()
    if what in ("A", "all"):
        part_A(panel)
    if what in ("B", "all"):
        part_B(panel)
    if what in ("C", "all"):
        part_C(panel)
    if what in ("D", "all"):
        part_D(panel)
    if what in ("E", "all"):
        part_E(panel)
    log(f"{what} done in {time.time() - t0:.0f} s")
    log.close()
