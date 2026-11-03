"""s05: between-season uncertainty (m2/D22, m3/D22) and the surveillance-design surface (Fig. 6D).

Season-specific intercepts use the season's own mean (as in the submitted manuscript) and two-way
(state x week) bootstrap standard errors (B = 1000). DerSimonian-Laird random effects and the Q-profile
interval for tau^2 are copied from experiment/original/code/r4_e_seasons.py (functions dl, q_profile);
the only change is that the per-season SE is now two-way instead of state-clustered, and S = 3
(calendar windows) and S = 4 (rule windows) are reported separately instead of being spliced.
"""
import numpy as np
from scipy import optimize, stats

from cc_common import SEED, Log, Part, boot_weights, fit_stat, load_panel, pct, spec_frame

log = Log("s05_seasons")
P = Part("s05_seasons", "s05_seasons.py", ["data/panel_nhsn_vintage_2026-10-20.csv.gz"], "measured", SEED + 5)
panel = load_panel()
B = 1000


def dl(a, se):
    w = 1 / se ** 2
    mu_f = np.sum(w * a) / w.sum()
    Q = float(np.sum(w * (a - mu_f) ** 2))
    S = len(a)
    tau2 = max(0.0, (Q - (S - 1)) / (w.sum() - np.sum(w ** 2) / w.sum()))
    ws = 1 / (se ** 2 + tau2)
    mu = float(np.sum(ws * a) / ws.sum())
    se_mu = float(np.sqrt(1 / ws.sum()))
    tq = stats.t.ppf(0.975, S - 1)
    return {"S": S, "Q": Q, "p_Q": float(stats.chi2.sf(Q, S - 1)), "tau2": tau2, "a_RE": mu, "se_RE": se_mu,
            "ci95_t": [mu - tq * se_mu, mu + tq * se_mu], "I2": max(0.0, (Q - (S - 1)) / Q)}


def q_profile(a, se, alpha=0.05):
    S = len(a)

    def Qgen(t2):
        w = 1 / (se ** 2 + t2)
        mu = np.sum(w * a) / w.sum()
        return np.sum(w * (a - mu) ** 2)
    lo_c, hi_c = stats.chi2.ppf([1 - alpha / 2, alpha / 2], S - 1)
    big = 10.0
    lo = 0.0 if Qgen(0) - lo_c < 0 else optimize.brentq(lambda t: Qgen(t) - lo_c, 0, big)
    hi = 0.0 if Qgen(0) - hi_c < 0 else optimize.brentq(lambda t: Qgen(t) - hi_c, 0, big)
    return [float(lo), float(hi)]


for si, spec in enumerate(["RA_rule4", "G_rule4", "G_cal3", "RA_cal3"]):
    if f"seasons|{spec}" in P.d["blocks"]:
        continue
    e = spec_frame(panel, spec)
    rows = []
    for ki, s in enumerate(sorted(e.season.unique())):
        es = e[e.season == s].reset_index(drop=True)
        f = fit_stat(es.x, es.D)
        rng = np.random.default_rng(SEED + 500 + 10 * si + ki)
        dr = np.array([fit_stat(es.x.to_numpy(), es.D.to_numpy(), boot_weights(es, "twoway", rng))["a"] for _ in range(B)])
        rows.append({"season": s, "N": int(len(es)), "n_weeks": int(es.week_end.nunique()), "a": f["a"], "b": f["b"],
                     "xbar": f["xbar"], "se_twoway": float(np.std(dr, ddof=1)), "ci95": pct(dr), "seed": SEED + 500 + 10 * si + ki})
        log(spec, rows[-1])
    a = np.array([r["a"] for r in rows]); se = np.array([r["se_twoway"] for r in rows])
    re = dl(a, se)
    re["tau2_qprofile95"] = q_profile(a, se)
    re["relSE_bound_sqrt_2_over_S_minus_1"] = float(np.sqrt(2 / (len(a) - 1)))
    re["cv_between"] = float(np.std(a, ddof=1) / np.mean(a))
    P.put(f"seasons|{spec}", {"rows": rows, "random_effects": re, "B": B})
    for r in rows:
        P.num(f"{spec}_a_{r['season']}", r["a"], ci=r["ci95"], seed=r["seed"], B=B)
    P.num(f"{spec}_a_RE", re["a_RE"], ci=re["ci95_t"], note=f"DerSimonian-Laird, t interval with S-1 = {len(a) - 1} df")
    P.num(f"{spec}_Q", re["Q"], p=re["p_Q"])
    P.num(f"{spec}_tau2", re["tau2"], ci=re["tau2_qprofile95"], note="Q-profile 95% interval")
    P.num(f"{spec}_I2", re["I2"])
    log(spec, "RE", re)

# relative-SE bound (Corollary 2) for Fig. 6B
S_grid = list(range(2, 21))
P.put("relse_bound", {"S": S_grid, "sqrt_2_over_S_minus_1": [float(np.sqrt(2 / (S - 1))) for S in S_grid]})
for S in (3, 4, 10):
    P.num(f"relSE_bound_S{S}", float(np.sqrt(2 / (S - 1))), mode="analytic")

# design surface: MSE of a T-period mean of R with per-period reported volume m,
# MSE(T, m) = (v_R + A/m)/T with v_R = a Rbar^2 and A = b Rbar^2 from the primary fit
e = spec_frame(panel, "RA_rule4")
f = fit_stat(e.x, e.D)
vR, A = f["a"] * f["xbar"] ** 2, f["b"] * f["xbar"] ** 2
P.num("design_vR_RA", vR, mode="analytic", note="v_R implied by the R_A intercept (weekly apparent scale)")
P.num("design_A_RA", A, mode="analytic")
P.num("design_mcross_RA", A / vR, mode="analytic", note="b/a: volume at which per-period noise equals the floor")
m = np.logspace(0, 4, 120)
targets = [vR / T for T in (4, 8, 16, 32)]
P.put("design", {"vR": vR, "A": A, "m": m.tolist(),
                 "iso": {f"{t:.6g}": ((vR + A / m) / t).tolist() for t in targets},
                 "targets_T_at_infinite_m": [4, 8, 16, 32]})
log.close()
