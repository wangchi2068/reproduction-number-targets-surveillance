"""s02: inference for the intercept a, slope b and b/a under four cluster schemes (M5/D10-D11, m14).

Schemes (B = 2000 each, percentile 95% intervals, xbar re-estimated in every replicate):
  state       resample states (the submitted manuscript's scheme)
  week        resample calendar weeks
  twoway      pigeonhole bootstrap: independent state and week multiplicities, product weights
              (Owen 2007; McCullagh 2000) -- the headline interval of the revision
  season_week resample weeks within season, states fixed
Analytic check: Cameron-Gelbach-Miller (2011) two-way cluster-robust sandwich for the final FGLS step.
Also reports the share of replicates with a > 0 and with b/a above the theory threshold used in the paper.
"""
import sys

import numpy as np

from cc_common import B_DEFAULT, SEED, SPECS, Log, Part, fgls_affine, fit_stat, load_panel, pct, rel_loss, spec_frame, boot_weights

SCHEMES = ["state", "week", "twoway", "season_week"]
B = int(sys.argv[1]) if len(sys.argv) > 1 else B_DEFAULT

log = Log("s02_inference")
P = Part("s02_inference", "s02_inference.py", ["data/panel_nhsn_vintage_2026-10-20.csv.gz"], "measured", SEED + 2)
panel = load_panel()


def cgm(e, a, b):
    """CGM two-way cluster-robust SE of (a, b) for the last weighted LS step (weights 1/mu^2)."""
    L, _ = rel_loss(e.x)
    X = np.column_stack([np.ones(len(e)), 1 / e.D.to_numpy()])
    mu = X @ np.array([a, b])
    W = 1 / mu ** 2
    u = L - mu
    Xw = X * W[:, None]
    bread = np.linalg.inv(X.T @ Xw)

    def meat(groups):
        S = np.zeros((2, 2))
        for _, idx in e.groupby(groups).indices.items():
            s = (Xw[idx] * u[idx, None]).sum(0)
            S += np.outer(s, s)
        return S
    G_s, G_w = e.state.nunique(), e.week_end.nunique()
    V = (bread @ meat("state") @ bread * G_s / (G_s - 1) + bread @ meat("week_end") @ bread * G_w / (G_w - 1)
         - bread @ meat(np.arange(len(e))) @ bread)
    return [float(np.sqrt(max(V[0, 0], 0))), float(np.sqrt(max(V[1, 1], 0)))]


for si, spec in enumerate(SPECS):
    e = spec_frame(panel, spec)
    f = fit_stat(e.x, e.D)
    P.num(f"{spec}_N", int(len(e)))
    P.num(f"{spec}_n_states", int(e.state.nunique()))
    P.num(f"{spec}_n_weeks", int(e.week_end.nunique()))
    P.num(f"{spec}_xbar", f["xbar"])
    se = cgm(e, f["a"], f["b"])
    P.num(f"{spec}_a_cgm_se", se[0]); P.num(f"{spec}_b_cgm_se", se[1])
    res = {}
    for ki, sch in enumerate(SCHEMES):
        if f"{spec}|{sch}" in P.d["blocks"]:
            res[sch] = P.d["blocks"][f"{spec}|{sch}"]
            continue
        rng = np.random.default_rng(SEED + 2 + 100 * si + ki)
        draws = np.empty((B, 3))
        for r in range(B):
            w = boot_weights(e, sch, rng)
            g = fit_stat(e.x.to_numpy(), e.D.to_numpy(), w)
            draws[r] = g["a"], g["b"], g["ba"]
        blk = {"B": B, "seed": SEED + 2 + 100 * si + ki,
               "a_ci": pct(draws[:, 0]), "b_ci": pct(draws[:, 1]), "ba_ci": pct(draws[:, 2]),
               "a_sd": float(np.std(draws[:, 0], ddof=1)), "b_sd": float(np.std(draws[:, 1], ddof=1)),
               "share_a_pos": float(np.mean(draws[:, 0] > 0)), "share_b_pos": float(np.mean(draws[:, 1] > 0))}
        P.put(f"{spec}|{sch}", blk)
        np.save(P.path.parent / f"s02_draws_{spec}_{sch}.npy", draws.astype(np.float32))
        res[sch] = blk
        log(spec, sch, {k: v for k, v in blk.items() if k != "B"})
    tw = res["twoway"]
    P.num(f"{spec}_a", f["a"], ci=tw["a_ci"], note="FGLS point estimate; 95% two-way (state x week) pigeonhole bootstrap",
          seed=tw["seed"], B=B)
    P.num(f"{spec}_b", f["b"], ci=tw["b_ci"], seed=tw["seed"], B=B)
    P.num(f"{spec}_ba", f["ba"], ci=tw["ba_ci"], seed=tw["seed"], B=B)
    P.num(f"{spec}_twoway_share_a_pos", tw["share_a_pos"])
    P.num(f"{spec}_a_sd_ratio_twoway_vs_state", tw["a_sd"] / res["state"]["a_sd"])
log.close()
