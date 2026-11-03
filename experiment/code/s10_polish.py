"""s10: final-polish checks requested by the round-3 review (P3, P4, D2, D4). New keys polish.* only;
every earlier part is read, never rewritten. Reuses the s09 data-matched simulator unchanged.

A. Within-window autocorrelation (P3/D5 i). Lag-1 autocorrelation of the relative deviation of each
   statistic from its state-season mean, dev = x / mean_block(x) - 1, pooled over consecutive valid
   weeks: r1 = sum(dev_t dev_{t+1}) / sum(dev_t^2). Computed on the rule-4 windows of the data and on
   simulated panels of the same design (100 replicates at v_R = 0.04) for every s09 setting and for
   AR(1) phi = 0.25. Counting noise makes r1 negative even without persistence, so the reference is
   the simulated distribution, not zero. A row-weighted version restricted to D >= 20 is also given.
B. External generation-interval bound (P3/D5 ii). Cowling et al. (2009) serial interval 3.6 d
   (95% CI 2.9, 4.3). New setting gi_cow_lo: true GI gamma(2.9, 2.9 x 1.6/3.6) with the analyst
   keeping 3.6 d; the upper end (4.3 d) is the s09 gi_long setting. Calibration line and inversion as in
   s09 part A (200 replicates x 4 v_R; same bootstrap draws of a).
C. Overdispersed admissions (P4/D8). Weekly state counts are replaced by a gamma-Poisson mixture with
   the same mean and variance-to-mean ratio 1 + delta (delta = 0.5, 1.0), independent across states and
   weeks (batching of reports). Calibrated b/a and v_R^impl per delta.
D. Arithmetic on existing draws: latent balance scale m_x = A/v_R (Proposition 1 at the data-matched
   simulator values) with its interval and per-week conversion; 90% interval of the ratio v_RA/v_G for a
   TOST-type equivalence statement on [0.8, 1.25]; reporting volatility tau allowed by the 5.3% upper
   bound and the lever comparison of Proposition 6 at that tau.

Usage: python3 s10_polish.py [sim|fit|all]   (sim rows resume-safe in logs/s10_sim_rows.jsonl)
"""
import json
import sys
import time

import numpy as np
import pandas as pd

import s09_inversion as s9
from cc_common import GI_FLU, LOGS, PARTS, SEED, Log, Part, load_panel, pct, spec_frame

ROWS = LOGS / "s10_sim_rows.jsonl"
NREP_R1 = 100
SEED_SIM = "SeedSequence([SEED, 10, cell_index, chunk_index])"
GI_COW_LO = (2.9, 2.9 * 1.6 / 3.6)
s9.GI_FLU["cow_lo"] = GI_COW_LO                     # sim_block looks GI up in this dict by name

NEW = {"gi_cow_lo": {"gi": "cow_lo"}, "od05": {"od": 0.5}, "od10": {"od": 1.0}, "ar025": {"phi": 0.25}}
LABEL = {"gi_cow_lo": "true GI 2.9 d (Cowling 2009 lower bound)", "od05": "admissions VMR 1.5",
         "od10": "admissions VMR 2.0", "ar025": "AR(1) phi = 0.25"}
CELLS = [(f"{s}|v{v}", s, dict(s9.BASE, **NEW[s]), v, s9.NREP_A) for s in NEW for v in s9.VGRID]
CELLS += [(f"r1_{s}|v0.04", s, dict(s9.BASE, **s9.SETTINGS[s]), 0.04, NREP_R1)
          for s in s9.SETTINGS if s != "bridge_s08"]
log = Log("s10_polish")

_sim_block = s9.sim_block


def sim_block_od(n, L, rbar, d0, vr, cfg, rng, meanfield=False):
    C, Rw = _sim_block(n, L, rbar, d0, vr, cfg, rng, meanfield)
    od = cfg.get("od", 0.0)
    if od > 0 and not meanfield:
        lam = np.zeros_like(C, dtype=float)
        pos = C > 0
        lam[pos] = rng.gamma(C[pos] / od, od)
        C = rng.poisson(lam).astype(float)
    return C, Rw


s9.sim_block = sim_block_od


def S_ar_cont(phi, T):
    """Copy of s07_theory.S_ar_cont (importing s07 would rewrite its part file)."""
    return T if phi == 0 else T * (1 + phi) / (1 - phi) - 2 * phi * (1 - phi ** T) / (1 - phi) ** 2


def levers(T, m, tau, phi=0.5, psi=0.5, rb=1.15, vr=0.04, k=0.35, rho=0.25):
    """Copy of s07_theory.levers (Proposition 6 terms)."""
    env = vr * S_ar_cont(phi, T) / T ** 2
    rep = tau ** 2 * (rb ** 2 * S_ar_cont(psi, T) + vr * S_ar_cont(phi * psi, T)) / T ** 2
    samp = (rb + rho * (1 + tau ** 2) * (rb ** 2 + vr) / k) / (T * m)
    return {"env": env, "rep": rep, "samp": samp, "total": env + rep + samp}


def r1_blocks(dev, block, t):
    """Pooled lag-1 autocorrelation over consecutive observations within blocks (dev already centred)."""
    o = np.lexsort((t, block))
    dv, bl, tt = dev[o], block[o], t[o]
    nxt = (bl[1:] == bl[:-1]) & (tt[1:] - tt[:-1] == 1)
    return float(np.sum(dv[:-1][nxt] * dv[1:][nxt]) / np.sum(dv ** 2)), int(nxt.sum())


def centred(x, block):
    codes = pd.factorize(block)[0]
    m = np.bincount(codes, weights=x) / np.bincount(codes)
    cnt = np.bincount(codes)[codes]
    dev = x / m[codes] - 1
    return dev, cnt


def r1_of(x, D, block, t, dmin=5):
    k = np.isfinite(x) & (D >= dmin)
    x, D, block, t = x[k], D[k], block[k], t[k]
    dev, cnt = centred(x, block)
    keep = cnt >= 3
    return r1_blocks(dev[keep], block[keep], t[keep])


def sim_r1(rows, sel):
    out = {}
    mk = s9.masks(rows)
    for nm in ("RA", "G"):
        m = mk[nm] & sel
        D = rows["D"][m] if nm == "G" else rows["Ct"][m]
        block = rows["state"][m].astype(np.int64)
        t = rows["week"][m] % 100
        out[nm] = r1_of(rows[nm][m], D, block, t)[0]
        out[nm + "_D20"] = r1_of(rows[nm][m], D, block, t, dmin=20)[0]
    return out


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
    for ci, (cid, sname, cfg, vr, nrep) in enumerate(CELLS):
        for ch in range(nrep // s9.CHUNK):
            reps = range(ch * s9.CHUNK, (ch + 1) * s9.CHUNK)
            if all((cid, r) in done for r in reps):
                continue
            t0 = time.time()
            rng = np.random.default_rng(np.random.SeedSequence([SEED, 10, ci, ch]))
            rows = s9.sim_panel_rows(cfg, vr, s9.CHUNK, rng)
            lines = []
            for j, r in enumerate(reps):
                sel = rows["rep"] == j
                fr = s9.fit_rep(rows, sel)
                fr["r1"] = sim_r1(rows, sel)
                lines.append(json.dumps({"cell": cid, "setting": sname, "v_R": vr, "rep": r, "chunk": ch,
                                         "seed": [SEED, 10, ci, ch], **fr}))
            with open(ROWS, "a") as f:
                f.write("\n".join(lines) + "\n")
            last = json.loads(lines[-1])
            log(f"{cid} chunk {ch} RA a={last['RA']['a']:.4f} b/a={last['RA']['b'] / last['RA']['a']:.1f} "
                f"r1={last['r1']['RA']:.3f} ({time.time() - t0:.1f} s)")


def empirical_r1(panel):
    out = {}
    for nm, spec in (("RA", "RA_rule4"), ("G", "G_rule4")):
        e = spec_frame(panel, spec)
        season = (e.week_end - pd.DateOffset(months=6)).dt.year
        block = pd.factorize(e.state.astype(str) + "_" + season.astype(str))[0]
        t = ((e.week_end - e.week_end.min()).dt.days // 7).to_numpy()
        x, D = e.x.to_numpy(float), e.D.to_numpy(float)
        r, npair = r1_of(x, D, block, t)
        r20, n20 = r1_of(x, D, block, t, dmin=20)
        # state-cluster bootstrap of r1 (blocks resampled whole)
        rng = np.random.default_rng(SEED + 1001 + (nm == "G"))
        st = pd.factorize(e.state)[0]
        ns = st.max() + 1
        bs = []
        for _ in range(1000):
            pick = rng.integers(0, ns, ns)
            idx = np.concatenate([np.flatnonzero(st == s) for s in pick])
            blk = block[idx] + 100000 * np.repeat(np.arange(ns), [np.sum(st == s) for s in pick])
            bs.append(r1_of(x[idx], D[idx], blk, t[idx])[0])
        out[nm] = {"r1": r, "r1_ci95": pct(bs), "pairs": npair, "r1_D20": r20, "pairs_D20": n20}
        log(f"empirical r1 {nm}: {r:.3f} {np.round(pct(bs), 3)} pairs={npair}; D>=20 {r20:.3f} ({n20})")
    return out


def calib(sub_by_v, nm):
    v = np.array(sorted(sub_by_v))
    y = np.array([sub_by_v[x][f"{nm}.a"].mean() for x in v])
    se = np.array([sub_by_v[x][f"{nm}.a"].std(ddof=1) / np.sqrt(len(sub_by_v[x])) for x in v])
    X = np.column_stack([np.ones_like(v), v])
    Wt = 1 / se ** 2
    cov = np.linalg.inv(X.T @ (X * Wt[:, None]))
    return cov @ (X.T @ (Wt * y)), cov


def fit_all(panel):
    P = Part("polish", "s10_polish.py", s9.INPUTS, "simulated",
             {"simulation": SEED_SIM, "r1_bootstrap": "SEED + 1001/1002", "calibration_draws": SEED + 1003})
    inv = json.loads((PARTS / "inversion.json").read_text())
    s07 = json.loads((PARTS / "s07_theory.json").read_text())["numbers"]
    rr = json.loads((PARTS / "reporting_regime.json").read_text())["numbers"]
    sens = inv["blocks"]["sensitivity"]
    vdr = np.load(PARTS / "s09_vR_draws.npz")
    jd = np.load(PARTS / "s09_joint_draws.npz")
    rows = [json.loads(l) for l in ROWS.read_text().splitlines() if l.strip()]
    df = pd.json_normalize(rows)

    # ---------------------------------------------------------------- A autocorrelation
    emp = empirical_r1(panel)
    P.put("r1_empirical", emp)
    r1sim = {}
    for cid in sorted(set(df.cell)):
        sub = df[df.cell == cid]
        if not (cid.startswith("r1_") or cid.endswith("|v0.04")):
            continue
        s = cid.split("|")[0].replace("r1_", "")
        r1sim[s] = {nm: {"mean": float(sub[f"r1.{nm}"].mean()), "q95": pct(sub[f"r1.{nm}"]),
                         "emp_percentile": float(np.mean(sub[f"r1.{nm}"] <= emp[nm]["r1"]))}
                    for nm in ("RA", "G")}
        for nm in ("RA", "G"):
            r1sim[s][nm]["mean_D20"] = float(sub[f"r1.{nm}_D20"].mean())
        r1sim[s]["n"] = int(len(sub))
    P.put("r1_simulated_v0.04", r1sim)
    for nm in ("RA", "G"):
        P.num(f"r1_emp_{nm}", emp[nm]["r1"], ci=emp[nm]["r1_ci95"], mode="measured", seed=SEED + 1001 + (nm == "G"),
              note="lag-1 autocorrelation of x/mean_block - 1 within state-season rule windows; state-cluster "
                   "bootstrap B=1000", pairs=emp[nm]["pairs"])
        P.num(f"r1_emp_{nm}_D20", emp[nm]["r1_D20"], mode="measured", pairs=emp[nm]["pairs_D20"])
        for s in ("M0", "ar025", "ar05"):
            P.num(f"r1_sim_{s}_{nm}", r1sim[s][nm]["mean"], ci=r1sim[s][nm]["q95"],
                  note="mean and 2.5-97.5% over simulated panels at v_R = 0.04",
                  emp_percentile=r1sim[s][nm]["emp_percentile"], mean_D20=r1sim[s][nm]["mean_D20"])
    # phi consistent with the data: linear interpolation of mean r1 over phi in {0, 0.25, 0.5}
    phis = np.array([0.0, 0.25, 0.5])
    for nm in ("RA", "G"):
        mr = np.array([r1sim[s][nm]["mean"] for s in ("M0", "ar025", "ar05")])
        e_ = emp[nm]["r1"]
        phi_hat = float(np.interp(e_, mr, phis)) if mr[0] <= mr[-1] else float(np.interp(e_, mr[::-1], phis[::-1]))
        lo, hi = emp[nm]["r1_ci95"]
        phi_ci = sorted(float(np.interp(z, mr, phis)) for z in (lo, hi))
        P.num(f"phi_matched_{nm}", phi_hat, ci=phi_ci,
              note="phi at which the simulated mean r1 equals the empirical r1 (interpolation over phi 0/0.25/0.5, "
                   "clamped to [0, 0.5]); interval maps the bootstrap interval of the empirical r1",
              sim_mean_r1=mr.tolist())
    log("r1 sim", json.dumps({s: {nm: round(r1sim[s][nm]["mean"], 3) for nm in ("RA", "G")} for s in r1sim}))

    # ---------------------------------------------------------------- B, C inversion in new settings
    rng = np.random.default_rng(SEED + 1003)
    a_draws = {nm: jd[f"{nm}_a"] for nm in ("RA", "G")}
    emp_a = {nm: inv["blocks"]["empirical"][nm]["a"] for nm in ("RA", "G")}
    emp_ba = {nm: inv["blocks"]["ba_calibrated"][nm]["empirical"] for nm in ("RA", "G")}
    newinv = {}
    for s in NEW:
        sub_by_v = {v: df[df.cell == f"{s}|v{v}"] for v in s9.VGRID}
        if any(len(x) < s9.NREP_A for x in sub_by_v.values()):
            log("incomplete", s)
            continue
        rec = {"label": LABEL[s]}
        vb = {}
        for nm in ("RA", "G"):
            coef, cov = calib(sub_by_v, nm)
            cd = rng.multivariate_normal(coef, cov, len(a_draws[nm]))
            v_hat = (emp_a[nm] - coef[0]) / coef[1]
            vb[nm] = (a_draws[nm] - cd[:, 0]) / cd[:, 1]
            ba_v = {v: float(sub_by_v[v][f"{nm}.b"].mean() / sub_by_v[v][f"{nm}.a"].mean()) for v in s9.VGRID}
            near = min(s9.VGRID, key=lambda v: abs(v - v_hat))
            reps = (sub_by_v[near][f"{nm}.b"] / sub_by_v[near][f"{nm}.a"]).to_numpy()
            rec[nm] = {"gamma0": float(coef[0]), "gamma1": float(coef[1]), "v_implied": float(v_hat),
                       "v_ci95": pct(vb[nm]), "ba_pooled_by_v": ba_v, "ba_at_nearest_v": ba_v[near],
                       "nearest_v": near, "ba_rep_q95_at_nearest": pct(reps),
                       "emp_ba_percentile_at_nearest": float(np.mean(reps <= emp_ba[nm])),
                       "b_mean_at_nearest": float(sub_by_v[near][f"{nm}.b"].mean())}
        rec["ratio_RA_over_G"] = {"value": rec["RA"]["v_implied"] / rec["G"]["v_implied"],
                                  "ci95": pct(vb["RA"] / vb["G"])}
        d = vb["RA"] - vb["G"]
        rec["agree"] = bool(pct(d)[0] < 0 < pct(d)[1])
        newinv[s] = rec
        for nm in ("RA", "G"):
            P.num(f"{s}_{nm}_vR_implied", rec[nm]["v_implied"], ci=rec[nm]["v_ci95"], note=LABEL[s])
            P.num(f"{s}_{nm}_ba_pooled_at_nearest_v", rec[nm]["ba_at_nearest_v"], note=LABEL[s],
                  nearest_v=rec[nm]["nearest_v"], replicate_q95=rec[nm]["ba_rep_q95_at_nearest"],
                  emp_percentile=rec[nm]["emp_ba_percentile_at_nearest"], empirical=emp_ba[nm])
        log(f"{s}: vRA={rec['RA']['v_implied']:.4f} {np.round(rec['RA']['v_ci95'], 4)} "
            f"vG={rec['G']['v_implied']:.4f} {np.round(rec['G']['v_ci95'], 4)} "
            f"ba RA={rec['RA']['ba_at_nearest_v']:.1f} G={rec['G']['ba_at_nearest_v']:.1f} agree={rec['agree']}")
    P.put("new_settings", newinv)
    # M0 reference b/a at the nearest grid value (same definition)
    m0c = inv["blocks"]["cells"]
    for nm in ("RA", "G"):
        P.num(f"M0_{nm}_ba_pooled_v0.04", m0c["M0|v0.04"][nm]["ba_pooled"], note="data-matched M0, v_R = 0.04 (s09)")
    if "gi_cow_lo" in newinv:
        for nm in ("RA", "G"):
            lo = newinv["gi_cow_lo"][nm]["v_implied"]
            hi = sens["gi_long"][nm]["v_implied"]
            P.num(f"vR_impl_{nm}_GI_2.9_to_4.3", [lo, hi],
                  note="point v_R^impl when the true GI mean is at the lower (2.9 d, gi_cow_lo) and upper (4.3 d, "
                       "s09 gi_long) end of the Cowling et al. (2009) 95% interval; analyst GI fixed at 3.6 d",
                  ci_lo=newinv["gi_cow_lo"][nm]["v_ci95"], ci_hi=sens["gi_long"][nm]["v_ci95"])

    # ---------------------------------------------------------------- D arithmetic on existing draws
    rhoh, kk = s9.BASE["rho_h"], s9.BASE["k"]
    gpw = 7 / GI_FLU["mid"][0]
    for nm in ("RA", "G"):
        rb = sens["M0"][nm]["Rbar_true"]
        v = sens["M0"][nm]["v_implied"]
        vd = vdr[f"M0__{nm}"]
        vd = vd[vd > 0]
        A = rb + rhoh * (rb ** 2 + v) / kk
        mx = A / v
        mxd = (rb + rhoh * (rb ** 2 + vd) / kk) / vd
        P.num(f"latent_mcross_{nm}", mx, ci=pct(mxd), unit="admissions per generation", mode="simulated",
              note="Proposition 1 at data-matched simulator values (Rbar_true, rho = admission probability 0.01, "
                   "k = 0.35) and v_R^impl; interval from the v_R^impl draws (positive draws)",
              A=A, Rbar_true=rb, n_nonpositive=int(np.sum(vdr[f"M0__{nm}"] <= 0)))
        P.num(f"latent_mcross_{nm}_per_week", mx * gpw, ci=[x * gpw for x in pct(mxd)],
              unit="admissions per week", mode="simulated", note="per-generation value x 7/3.6 generations per week")
    rat = vdr["M0__RA"] / vdr["M0__G"]
    q90 = pct(rat, (5, 95))
    P.num("ratio_vR_RA_over_G_ci90", float(sens["M0"]["ratio_RA_over_G"]["value"]), ci=q90, mode="simulated",
          note="5-95% of the ratio draws (TOST at alpha = 0.05 on [0.8, 1.25] holds if the 90% interval is inside)",
          within_0_8_1_25=bool(0.8 <= q90[0] and q90[1] <= 1.25), ci95=pct(rat))
    # tau allowed by the 5.3% upper bound on the reporting share of the intercept (theory baseline of s07)
    up = rr["RA_rule4_delta_a_share"]["upper95_one_sided"]
    RB, VR = 1.15, 0.04
    tau = float(np.sqrt(up / (1 - up) * VR / (RB ** 2 + VR)))
    base, x10, std = levers(10, 50, tau), levers(10, 500, tau), levers(10, 50, 0.0)
    P.num("tau_allowed_by_upper_share", tau, mode="analytic",
          note="tau with v_H/(v_R + v_H) equal to the one-sided 95% upper bound of the completeness share of a "
               "(reporting_regime.RA_rule4_delta_a_share), at the s07 theory baseline Rbar 1.15, v_R 0.04",
          upper_share=up)
    P.num("lever_at_tau_allowed_total", base["total"] * 1e3, unit="1e-3", mode="analytic",
          note="Proposition 6 baseline T=10, m=50, phi=psi=0.5 at tau_allowed",
          rep_share=base["rep"] / base["total"], standardised=std["total"] * 1e3, cases_x10=x10["total"] * 1e3)
    log("tau_allowed", tau, "base", base, "std", std, "x10", x10, "ratio CI90", q90)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    if what in ("sim", "all"):
        run_sims()
        log(f"sim stage {time.time() - t0:.0f} s")
    if what in ("fit", "all"):
        fit_all(load_panel())
    log(f"{what} done in {time.time() - t0:.0f} s")
    log.close()
