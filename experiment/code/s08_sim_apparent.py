"""s08: simulation validation of the apparent (hospital-admission) series [reviewer M2/D33].

Question: if transmission follows a daily renewal process with known R_t, and only weekly hospital
admissions are observed, what do the weekly statistics used in the paper (growth factor G with the
original D_t weights, and the apparent reproduction number R_A) recover when the affine loss
L = a + b/D is fitted, and how do a and b/a relate to v_R/Rbar^2 and A/v_R?

Design (adapted from the user's original simulator experiment/original/code/r4_sim_core.py, which
uses the same gamma-Poisson renewal step, discretised gamma kernels and weekly aggregation):
- Independent "states" (trajectories). R is constant within a calendar week and independent across
  weeks: R_w ~ Gamma(shape Rbar^2/v_R, scale v_R/Rbar) (mean Rbar, variance v_R; gamma chosen because
  it is positive and matches the manuscript's Monte Carlo and Prop 5 assumption). Rbar = 1.2 per
  generation (growth phase).
- Daily infections: I_d ~ Poisson(Gamma(k Lambda_d, R_d/k)), Lambda_d = sum_s g_s I_{d-s}, i.e. NB
  offspring with dispersion k (mean R Lambda, variance R Lambda + R^2 Lambda/k); g = discretised gamma
  generation interval, mean 3.6 d, sd 1.6 d (GI_FLU['mid'], Cowling et al. 2009).
- Admissions: each infection is admitted with probability rho_h (binomial), after a gamma
  infection-to-admission delay with mean 5.5 d and sd 3.0 d, distributed multinomially over days.
  The delay is an ASSUMPTION (incubation about 1.4-2 d, Lessler et al. 2009 Lancet Infect Dis
  9:291-300, plus roughly 3-4 d from symptom onset to admission); no sensitivity analysis on it.
- Weekly admissions C_w. Initial weekly admissions at the first observation week are drawn
  log-uniform in [2, 150]; with weekly growth ~1.43 over 4 warm-up + 6 observation weeks this gives
  D_t from below 5 to a few thousand (median ~ 40, close to NHSN early windows).
- Observation weeks: the 6 weeks after a 4-week warm-up (all in growth phase), D_t >= 5.
  G = C_{t+1}/D_t with D_t = 0.65 C_t + 0.25 C_{t-1} + 0.10 C_{t-2} (W_ORIG);
  R_A = wl_R((C_{t+1}+0.5)/(C_t+0.5)) (cc_common.wl_R, Wallinga-Lipsitch, gamma GI 3.6/1.6 d).
  Both are fitted with cc_common.fit_stat (iterated FGLS of L = ((X - Xbar)/Xbar)^2 on 1/D_t).
- Grid: v_R in {0.01, 0.02, 0.04, 0.08, 0.16} x k in {0.1, 0.35, 1.0} at rho_h = 0.01, plus
  rho_h in {0.005, 0.02} at v_R = 0.04, k = 0.35 (17 scenarios). NREP = 20 replicate seeds x NS = 500
  states per scenario (10000 trajectories). Seed of (scenario i, replicate r) = SEED + 800 + 1000 i + r.
  Point estimates: pooled fit over all replicates; MC 95% interval: pooled +- 1.96 sd(replicate
  estimates)/sqrt(NREP).

Mapping to theory (stated honestly, two levels):
- Nominal (as in the manuscript, week treated as one generation): a_nom = v_R/Rbar^2,
  (b/a)_nom = A/v_R with A = Rbar + rho_h (Rbar^2+v_R)/k; with the manuscript's denominator-noise
  (b-delta, Eq. 25) correction (b/a)_nom,delta = [A + (1-rho_h)(Rbar^2+v_R) sum w^2]/v_R, sum w^2 = 0.495 for
  G and 1 for R_A (single-week denominator C_t). Same formula as review/checks/theory_checks.py
  (empirical_scenarios block).
- Statistic scale: a noise-free (mean-field, infinite population) run of the same renewal model with
  the same R_w distribution gives the D -> infinity limit a_inf = Var(X)/E(X)^2 of each statistic X;
  this is the quantity the fitted intercept a should recover. The b/a prediction on this scale uses
  the count-ratio Y underlying X (Y = G for G; Y = C_{t+1}/C_t for R_A): Rbar_Y = E(Y), v_Y = Var(Y)
  in the noise-free run, A_Y = Rbar_Y + rho_h (Rbar_Y^2 + v_Y)/k, (b/a)_stat = A_Y/v_Y and the b-delta
  version [A_Y + (1-rho_h)(Rbar_Y^2 + v_Y) sum w^2]/v_Y.

Outputs: analysis/parts/s08_sim.json, logs/s08_sim_apparent.log, logs/s08_sim_rows.jsonl (one line
per (scenario, replicate), resume-safe), logs/s08_obs/*.npz (per-replicate observations).
Usage: python3 s08_sim_apparent.py [--only-summary]
"""
import json
import sys
import time

import numpy as np
from scipy import stats as sps

from cc_common import LOGS, SEED, W_ORIG, Log, Part, fit_stat, wl_R

RBAR = 1.2
GI = (3.6, 1.6)
DELAY = (5.5, 3.0)
PRE_W, NW = 4, 6
NS, NREP = 500, 20
D0_RANGE = (2.0, 150.0)
SW2 = {"G": float(np.sum(np.square(W_ORIG))), "RA": 1.0}

SCEN = []
for vr in (0.01, 0.02, 0.04, 0.08, 0.16):
    for k in (0.1, 0.35, 1.0):
        SCEN.append({"v_R": vr, "k": k, "rho_h": 0.01})
for rh in (0.005, 0.02):
    SCEN.append({"v_R": 0.04, "k": 0.35, "rho_h": rh})
for i, s in enumerate(SCEN):
    s["id"] = f"v{s['v_R']}_k{s['k']}_r{s['rho_h']}"
    s["idx"] = i
BASE_ID = "v0.04_k0.35_r0.01"

ROWS = LOGS / "s08_sim_rows.jsonl"
OBS = LOGS / "s08_obs"
OBS.mkdir(exist_ok=True)


def disc_gamma(mean, sd, s0=1, q=0.999):
    """Discretised gamma on s >= s0 (mass below s0 folded into s0), truncated at quantile q
    (copied from r4_sim_core.disc_gamma)."""
    shape, scale = (mean / sd) ** 2, sd ** 2 / mean
    smax = int(np.ceil(sps.gamma.ppf(q, shape, scale=scale))) + 1
    s = np.arange(s0, smax + 1)
    F = sps.gamma.cdf(s + 0.5, shape, scale=scale) - sps.gamma.cdf(s - 0.5, shape, scale=scale)
    F[0] += sps.gamma.cdf(s0 - 0.5, shape, scale=scale)
    return F / F.sum()


G_KER = disc_gamma(*GI, s0=1)
D_KER = disc_gamma(*DELAY, s0=0)


def growth_rate(R, g=G_KER):
    """Daily growth rate r solving the discrete Euler-Lotka equation R sum g_s e^{-r s} = 1."""
    s = np.arange(1, len(g) + 1)
    lo, hi = -0.5, 0.5
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if R * np.sum(g * np.exp(-mid * s)) > 1:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def simulate(S, vr, k, rho_h, rng, meanfield=False, d0=None):
    """Return weekly admissions C (S x W) and the weekly R path (S x W)."""
    W = PRE_W + NW + 2
    ndays = 7 * W
    smax = len(G_KER)
    hist = smax + len(D_KER)
    r = growth_rate(RBAR)
    lam = np.exp(7 * r)
    if d0 is None:
        d0 = np.exp(rng.uniform(np.log(D0_RANGE[0]), np.log(D0_RANGE[1]), S))
    # daily infections at day 0 so that admissions in week PRE_W are about d0
    # weekly admissions in week j = rho_h i0 E[e^{-r delay}] e^{7 r j} sum_{u=0..6} e^{r u}
    mgf = np.sum(D_KER * np.exp(-r * np.arange(len(D_KER))))
    i0 = d0 / (rho_h * mgf * lam ** PRE_W * np.sum(np.exp(r * np.arange(7))))
    I = np.zeros((S, hist + ndays))
    hd = np.arange(-hist, 0)
    I[:, :hist] = i0[:, None] * np.exp(r * hd)[None, :]
    if not meanfield:
        I[:, :hist] = np.round(I[:, :hist])
    shp = RBAR ** 2 / vr
    Rw = rng.gamma(shp, RBAR / shp, (S, W))
    grev = G_KER[::-1]
    for d in range(ndays):
        j = hist + d
        L = I[:, j - smax:j] @ grev
        Rd = Rw[:, d // 7]
        if meanfield:
            I[:, j] = Rd * L
        else:
            x = np.zeros(S)
            pos = L > 0
            if pos.any():
                x[pos] = rng.gamma(k * L[pos], Rd[pos] / k)
            I[:, j] = rng.poisson(x)
    # admissions with delay
    nd = hist + ndays
    if meanfield:
        adm_inf = rho_h * I
        out = np.zeros((S, nd + len(D_KER)))
        for j, p in enumerate(D_KER):
            out[:, j:j + nd] += p * adm_inf
    else:
        adm_inf = rng.binomial(I.astype(np.int64), rho_h)
        M = rng.multinomial(adm_inf, D_KER)        # (S, nd, len(D_KER))
        out = np.zeros((S, nd + len(D_KER)))
        for j in range(len(D_KER)):
            out[:, j:j + nd] += M[:, :, j]
    Cd = out[:, hist:hist + ndays]
    C = Cd.reshape(S, W, 7).sum(axis=2)
    return C, Rw


def weekly_stats(C):
    """G, R_A, D_t for observation weeks t = PRE_W .. PRE_W+NW-1 (all states), plus g = C_{t+1}/C_t."""
    w = np.asarray(W_ORIG)
    t = np.arange(PRE_W, PRE_W + NW)
    D = w[0] * C[:, t] + w[1] * C[:, t - 1] + w[2] * C[:, t - 2]
    G = C[:, t + 1] / np.where(D > 0, D, np.nan)
    gcc = (C[:, t + 1] + 0.5) / (C[:, t] + 0.5)
    RA = wl_R(gcc, *GI)
    g1 = C[:, t + 1] / np.where(C[:, t] > 0, C[:, t], np.nan)
    return D, G, RA, g1


def fit_rows(D, G, RA):
    keep = np.isfinite(G) & (D >= 5)
    out = {"N": int(keep.sum()), "D_med": float(np.median(D[keep])),
           "D_q": [float(x) for x in np.percentile(D[keep], [1, 5, 50, 95, 99])]}
    for nm, X in (("G", G), ("RA", RA)):
        f = fit_stat(X[keep], D[keep])
        out[nm] = {k_: float(v) for k_, v in f.items()}
    return out


def meanfield_limits(vr, rng, S=20000):
    C, _ = simulate(S, vr, 1.0, 0.01, rng, meanfield=True, d0=np.full(S, 1e6))
    D, G, RA, g1 = weekly_stats(C)
    res = {}
    for nm, X in (("G", G), ("RA", RA), ("g", g1)):
        X = X.ravel()
        res[nm] = {"mean": float(X.mean()), "var": float(X.var()), "a_inf": float(X.var() / X.mean() ** 2)}
    return res


def done_keys():
    if not ROWS.exists():
        return set()
    out = set()
    for line in ROWS.read_text().splitlines():
        try:
            r = json.loads(line)
            out.add((r["id"], r["rep"]))
        except Exception:
            pass
    return out


def run_sims(log):
    done = done_keys()
    for s in SCEN:
        for rep in range(NREP):
            if (s["id"], rep) in done:
                continue
            seed = SEED + 800 + 1000 * s["idx"] + rep
            rng = np.random.default_rng(seed)
            t0 = time.time()
            C, Rw = simulate(NS, s["v_R"], s["k"], s["rho_h"], rng)
            D, G, RA, _ = weekly_stats(C)
            np.savez_compressed(OBS / f"{s['id']}_r{rep}.npz", D=D, G=G, RA=RA,
                                **({"C": C, "Rw": Rw} if (s["id"] == BASE_ID and rep == 0) else {}))
            row = {"id": s["id"], "rep": rep, "seed": seed, **{k_: s[k_] for k_ in ("v_R", "k", "rho_h")},
                   **fit_rows(D, G, RA), "sec": round(time.time() - t0, 2)}
            with open(ROWS, "a") as f:
                f.write(json.dumps(row) + "\n")
            log(f"{s['id']} rep {rep} N={row['N']} Dmed={row['D_med']:.1f} "
                f"G a={row['G']['a']:.4f} b/a={row['G']['ba']:.1f} | RA a={row['RA']['a']:.4f} "
                f"b/a={row['RA']['ba']:.1f} ({row['sec']} s)")


def summarise(log):
    P = Part("s08_sim", "s08_sim_apparent.py", ["simulation only (no data)"], "simulated",
             {"base": SEED + 800, "rule": "SEED + 800 + 1000*scenario_index + replicate", "meanfield": SEED + 899})
    rows = [json.loads(l) for l in ROWS.read_text().splitlines() if l.strip()]
    P.put("design", {"Rbar": RBAR, "GI_mean_sd_days": GI, "delay_mean_sd_days": DELAY,
                     "R_distribution": "Gamma(mean Rbar, var v_R), constant within a calendar week, iid across weeks",
                     "offspring": "NB(k) via Poisson(Gamma(k Lambda, R/k))", "warmup_weeks": PRE_W,
                     "obs_weeks": NW, "states_per_rep": NS, "reps": NREP, "D0_range": D0_RANGE,
                     "D_min": 5, "weights": W_ORIG, "sum_w2": SW2, "scenarios": SCEN})
    mf = {}
    rng = np.random.default_rng(SEED + 899)
    for vr in sorted({s["v_R"] for s in SCEN}):
        mf[str(vr)] = meanfield_limits(vr, rng)
        log("meanfield", vr, json.dumps(mf[str(vr)]))
    P.put("meanfield", mf)

    table = []
    for s in SCEN:
        rs = [r for r in rows if r["id"] == s["id"]]
        if len(rs) < NREP:
            log("incomplete", s["id"], len(rs)); continue
        Ds, Gs, RAs = [], [], []
        for r in range(NREP):
            z = np.load(OBS / f"{s['id']}_r{r}.npz")
            Ds.append(z["D"]); Gs.append(z["G"]); RAs.append(z["RA"])
        pooled = fit_rows(np.concatenate(Ds), np.concatenate(Gs), np.concatenate(RAs))
        vr, k, rh = s["v_R"], s["k"], s["rho_h"]
        A = RBAR + rh * (RBAR ** 2 + vr) / k
        rec = {"id": s["id"], "v_R": vr, "k": k, "rho_h": rh, "N": pooled["N"], "D_med": pooled["D_med"],
               "D_q": pooled["D_q"], "a_nom": vr / RBAR ** 2, "A": A, "ba_nom": A / vr}
        m = mf[str(vr)]
        for nm, Y in (("G", "G"), ("RA", "g")):
            est = pooled[nm]
            reps_a = np.array([r[nm]["a"] for r in rs]); reps_b = np.array([r[nm]["b"] for r in rs])
            reps_ba = np.array([r[nm]["ba"] for r in rs])

            def ci(x0, xs):
                se = np.nanstd(xs, ddof=1) / np.sqrt(np.isfinite(xs).sum())
                return [float(x0 - 1.96 * se), float(x0 + 1.96 * se)]
            rb, vy = m[Y]["mean"], m[Y]["var"]
            Ay = rb + rh * (rb ** 2 + vy) / k
            rec[nm] = {"a": est["a"], "a_ci": ci(est["a"], reps_a), "b": est["b"], "b_ci": ci(est["b"], reps_b),
                       "ba": est["ba"], "ba_ci": ci(est["ba"], reps_ba), "xbar": est["xbar"],
                       "a_inf": m[nm]["a_inf"],
                       "ba_nom_delta": (A + (1 - rh) * (RBAR ** 2 + vr) * SW2[nm]) / vr,
                       "ba_stat": Ay / vy,
                       "ba_stat_delta": (Ay + (1 - rh) * (rb ** 2 + vy) * SW2[nm]) / vy}
        table.append(rec)
        log(f"{s['id']}: G a={rec['G']['a']:.4f} (nom {rec['a_nom']:.4f}, inf {rec['G']['a_inf']:.4f}) "
            f"b/a={rec['G']['ba']:.1f} (nom {A / vr:.1f}, nom+d {rec['G']['ba_nom_delta']:.1f}, "
            f"stat {rec['G']['ba_stat']:.1f}, stat+d {rec['G']['ba_stat_delta']:.1f}) | "
            f"RA a={rec['RA']['a']:.4f} (inf {rec['RA']['a_inf']:.4f}) b/a={rec['RA']['ba']:.1f} "
            f"(stat+d {rec['RA']['ba_stat_delta']:.1f})")
    P.put("scenarios", table)

    # ---------------------------------------------------------------- summary numbers
    grid = [r for r in table if r["rho_h"] == 0.01]
    for nm in ("G", "RA"):
        x = np.array([r["a_nom"] for r in grid]); y = np.array([r[nm]["a"] for r in grid])
        xi = np.array([r[nm]["a_inf"] for r in grid])
        slope = float(np.sum(x * y) / np.sum(x * x))
        P.num(f"{nm}_a_slope_on_nominal", slope,
              note="least-squares slope through the origin of recovered a on v_R/Rbar^2 (15 grid scenarios, rho_h=0.01)")
        P.num(f"{nm}_a_ratio_to_meanfield_median", float(np.median(y / xi)),
              note="median over grid of recovered a / noise-free limit a_inf of the same statistic")
        P.num(f"{nm}_ainf_over_nominal_median", float(np.median(xi / x)),
              note="noise-free attenuation/amplification factor a_inf/(v_R/Rbar^2)")
        mono = all(np.all(np.diff([r[nm]["a"] for r in sorted([g for g in grid if g["k"] == kk], key=lambda z: z["v_R"])]) > 0)
                   for kk in (0.1, 0.35, 1.0))
        P.num(f"{nm}_a_monotone_in_vR", bool(mono), note="recovered a strictly increasing in v_R within every k")
        for lab, key in (("nom", None), ("nom_delta", "ba_nom_delta"), ("stat", "ba_stat"), ("stat_delta", "ba_stat_delta")):
            th = np.array([r["ba_nom"] if key is None else r[nm][key] for r in grid])
            ba = np.array([r[nm]["ba"] for r in grid])
            P.num(f"{nm}_ba_ratio_{lab}_median", float(np.median(ba / th)),
                  note=f"median over grid of recovered b/a divided by predicted b/a ({lab})",
                  range=[float(np.min(ba / th)), float(np.max(ba / th))])
        base = next(r for r in table if r["id"] == BASE_ID)
        P.num(f"{nm}_base_a", base[nm]["a"], ci=base[nm]["a_ci"], note="baseline v_R=0.04, k=0.35, rho_h=0.01")
        P.num(f"{nm}_base_b", base[nm]["b"], ci=base[nm]["b_ci"])
        P.num(f"{nm}_base_ba", base[nm]["ba"], ci=base[nm]["ba_ci"],
              note=f"nominal A/v_R={base['ba_nom']:.1f}, nominal+delta={base[nm]['ba_nom_delta']:.1f}")
    P.num("base_D_median", next(r for r in table if r["id"] == BASE_ID)["D_med"])
    P.num("base_D_q01_q99", next(r for r in table if r["id"] == BASE_ID)["D_q"],
          note="quantiles 1, 5, 50, 95, 99 % of D_t (D_t >= 5), pooled baseline")
    P.num("sim_runtime_sec", float(sum(r["sec"] for r in rows)), note="sum of per-replicate simulation times")
    P.num("n_trajectories_per_scenario", NS * NREP)

    # ---------------------------------------------------------------- figure blocks
    z = np.load(OBS / f"{BASE_ID}_r0.npz")
    C, Rw = z["C"], z["Rw"]
    D0 = 0.65 * C[:, PRE_W] + 0.25 * C[:, PRE_W - 1] + 0.10 * C[:, PRE_W - 2]
    j = int(np.argsort(np.abs(np.log(np.maximum(D0, 0.5)) - np.log(40)))[0])
    w = np.asarray(W_ORIG)
    t = np.arange(2, C.shape[1] - 1)
    Dj = w[0] * C[j, t] + w[1] * C[j, t - 1] + w[2] * C[j, t - 2]
    P.put("fig3_traj", {"state": j, "week": t, "R_true": Rw[j, t], "G": C[j, t + 1] / Dj,
                        "RA": wl_R((C[j, t + 1] + 0.5) / (C[j, t] + 0.5), *GI), "C": C[j], "D": Dj,
                        "obs_weeks": [PRE_W, PRE_W + NW - 1], "Rbar": RBAR})
    Ds = np.concatenate([np.load(OBS / f"{BASE_ID}_r{r}.npz")["D"] for r in range(NREP)]).ravel()
    bins = {}
    for nm in ("G", "RA"):
        X = np.concatenate([np.load(OBS / f"{BASE_ID}_r{r}.npz")[nm] for r in range(NREP)]).ravel()
        keep = np.isfinite(X) & (Ds >= 5)
        Xk, Dk = X[keep], Ds[keep]
        L = ((Xk - Xk.mean()) / Xk.mean()) ** 2
        qs = np.quantile(1 / Dk, np.linspace(0, 1, 11))
        idx = np.clip(np.searchsorted(qs, 1 / Dk, side="right") - 1, 0, 9)
        bx, bm, blo, bhi = [], [], [], []
        brng = np.random.default_rng(SEED + 898)
        for b in range(10):
            Lb = L[idx == b]
            bx.append(float(np.mean(1 / Dk[idx == b]))); bm.append(float(Lb.mean()))
            boots = [Lb[brng.integers(0, len(Lb), len(Lb))].mean() for _ in range(500)]
            blo.append(float(np.percentile(boots, 2.5))); bhi.append(float(np.percentile(boots, 97.5)))
        base = next(r for r in table if r["id"] == BASE_ID)
        bins[nm] = {"invD": bx, "L_mean": bm, "L_lo": blo, "L_hi": bhi, "a": base[nm]["a"], "b": base[nm]["b"],
                    "a_nom": base["a_nom"], "a_inf": base[nm]["a_inf"],
                    "interval": "bin-mean percentile bootstrap (500 resamples), 10 quantile bins of 1/D"}
    P.put("fig3_bins", bins)
    log("summary written")


if __name__ == "__main__":
    log = Log("s08_sim_apparent")
    if "--only-summary" not in sys.argv:
        run_sims(log)
    summarise(log)
    log.close()
