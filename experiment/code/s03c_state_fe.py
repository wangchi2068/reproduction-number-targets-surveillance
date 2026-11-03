"""s03c_state_fe.py: Within-state fixed-effects analysis and spatial heterogeneity tests.
Fits L_st = a_s + b_w / C_st via iterated FGLS on RA_rule4 (n = 1748, 51 states).
Compares pooled cross-sectional fit with state-specific intercepts:
- Unweighted mean, median, row-weighted and case-weighted state intercepts.
- Scale D_t sensitivity.
- Week intercepts only and State + Week intercepts.
- Paired ratio of cross-sectional b/a to within-state ratio.
- 1000 two-way bootstrap replicates for paired confidence intervals.
Writes analysis/parts/s11_fe.json.
Usage: python3 s03c_state_fe.py [B=1000]
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from cc_common import (
    ANALYSIS,
    PARTS,
    SEED,
    Log,
    Part,
    boot_weights,
    fgls_affine,
    load_panel,
    pct,
    rel_loss,
    spec_frame,
)

B = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
log = Log("s03c_state_fe")
P = Part("s11_fe", "s03c_state_fe.py", ["data/panel_nhsn_vintage_2026-10-20.csv.gz"], "measured", SEED + 32)

p = load_panel()
e = spec_frame(p, "RA_rule4")
n = len(e)
L, _ = rel_loss(e.x.to_numpy())
D = e.D.to_numpy()
Ct = e.C_t.to_numpy()
states = e.state.astype(str).to_numpy()
weeks = e.week_end.astype(str).to_numpy()

uniq_states, state_codes = np.unique(states, return_inverse=True)
uniq_weeks, week_codes = np.unique(weeks, return_inverse=True)
state_codes = np.asarray(state_codes, int).ravel()
week_codes = np.asarray(week_codes, int).ravel()

Ns = len(uniq_states)
Nw = len(uniq_weeks)
n_per_state = np.bincount(state_codes)
cases_per_state = np.bincount(state_codes, weights=Ct)

# Design matrices
X_state = np.column_stack([np.eye(Ns)[state_codes], 1.0 / Ct])
X_state_D = np.column_stack([np.eye(Ns)[state_codes], 1.0 / D])
X_week = np.column_stack([np.eye(Nw)[week_codes], 1.0 / Ct])
X_state_week = np.column_stack([np.eye(Ns)[state_codes], np.eye(Nw)[:, 1:][week_codes], 1.0 / Ct])


def fit_fgls(L_vec, X_mat, w_vec=None, maxit=100, tol=1e-9):
    sw = np.sqrt(w_vec) if w_vec is not None else np.ones(len(L_vec))
    keep = sw > 0
    Lk, Xk, swk = L_vec[keep], X_mat[keep], sw[keep]
    try:
        beta = np.linalg.lstsq(Xk * swk[:, None], Lk * swk, rcond=None)[0]
    except Exception:
        return None
    for _ in range(maxit):
        mu = np.clip(Xk @ beta, 1e-6, None)
        ww = swk / mu
        try:
            new = np.linalg.lstsq(Xk * ww[:, None], Lk * ww, rcond=None)[0]
        except Exception:
            break
        if np.max(np.abs(new - beta) / (np.abs(beta) + 1e-12)) < tol:
            beta = new
            break
        beta = new
    return beta


# 1. Point estimates
# A. Pooled cross-section
a_pool, b_pool = fgls_affine(L, Ct)
ba_pool = b_pool / a_pool

# B. State intercepts
beta_state = fit_fgls(L, X_state)
as_state = beta_state[:-1]
bw_state = beta_state[-1]
mean_as = float(np.mean(as_state))
median_as = float(np.median(as_state))
row_as = float(np.average(as_state, weights=n_per_state))
case_as = float(np.average(as_state, weights=cases_per_state))

r_mean = bw_state / mean_as
r_median = bw_state / median_as
r_row = bw_state / row_as
r_case = bw_state / case_as
paired_r = ba_pool / r_mean

# C. Scale D_t
beta_state_D = fit_fgls(L, X_state_D)
as_state_D = beta_state_D[:-1]
bw_state_D = beta_state_D[-1]
r_state_D = bw_state_D / float(np.mean(as_state_D))

# D. Week intercepts only
beta_week = fit_fgls(L, X_week)
bw_week = float(beta_week[-1])

# E. State and week intercepts
beta_sw = fit_fgls(L, X_state_week)
bw_sw = float(beta_sw[-1])

log("point estimates", bw_state, mean_as, r_mean, paired_r)

# 2. Bootstrap
log(f"running two-way bootstrap (B={B})...")
boot_rng = np.random.default_rng(20261104)

draws = {
    "pool_b": [], "pool_ba": [],
    "state_b": [], "mean_as": [],
    "r_mean": [], "r_median": [], "r_row": [], "r_case": [],
    "r_state_D": [], "bw_state_D": [],
    "bw_week": [], "bw_sw": [],
    "paired_ratio": []
}

for _ in range(B):
    w = boot_weights(e, "twoway", boot_rng)
    L_b, _ = rel_loss(e.x.to_numpy(), w=w)
    
    # pooled
    ap_b, bp_b = fgls_affine(L_b, Ct, w=w)
    bap_b = bp_b / ap_b if ap_b > 0 else np.nan
    draws["pool_b"].append(bp_b)
    draws["pool_ba"].append(bap_b)
    
    # state FE
    bs_b = fit_fgls(L_b, X_state, w_vec=w)
    if bs_b is not None:
        as_b = bs_b[:-1]
        bw_b = bs_b[-1]
        m_as_b = float(np.mean(as_b))
        med_as_b = float(np.median(as_b))
        row_as_b = float(np.average(as_b, weights=n_per_state))
        case_as_b = float(np.average(as_b, weights=cases_per_state))
        
        rm_b = bw_b / m_as_b if m_as_b > 0 else np.nan
        draws["state_b"].append(bw_b)
        draws["mean_as"].append(m_as_b)
        draws["r_mean"].append(rm_b)
        draws["r_median"].append(bw_b / med_as_b if med_as_b > 0 else np.nan)
        draws["r_row"].append(bw_b / row_as_b if row_as_b > 0 else np.nan)
        draws["r_case"].append(bw_b / case_as_b if case_as_b > 0 else np.nan)
        draws["paired_ratio"].append(bap_b / rm_b if (rm_b > 0 and bap_b > 0) else np.nan)
        
    # state D
    bs_D_b = fit_fgls(L_b, X_state_D, w_vec=w)
    if bs_D_b is not None:
        draws["bw_state_D"].append(bs_D_b[-1])
        m_as_D_b = float(np.mean(bs_D_b[:-1]))
        draws["r_state_D"].append(bs_D_b[-1] / m_as_D_b if m_as_D_b > 0 else np.nan)
        
    # week only
    bw_b = fit_fgls(L_b, X_week, w_vec=w)
    if bw_b is not None:
        draws["bw_week"].append(bw_b[-1])
        
    # state + week
    bsw_b = fit_fgls(L_b, X_state_week, w_vec=w)
    if bsw_b is not None:
        draws["bw_sw"].append(bsw_b[-1])

res = {
    "pooled": {"b": {"value": b_pool, "ci95": pct(draws["pool_b"])},
               "ba": {"value": ba_pool, "ci95": pct(draws["pool_ba"])}},
    "state_fe": {
        "slope": {"value": bw_state, "ci95": pct(draws["state_b"])},
        "mean_as": {"value": mean_as, "ci95": pct(draws["mean_as"])},
        "r_mean": {"value": r_mean, "ci95": pct(draws["r_mean"])},
        "r_median": {"value": r_median, "ci95": pct(draws["r_median"])},
        "r_row": {"value": r_row, "ci95": pct(draws["r_row"])},
        "r_case": {"value": r_case, "ci95": pct(draws["r_case"])},
        "paired_ratio": {"value": paired_r, "ci95": pct(draws["paired_ratio"])},
    },
    "state_fe_D": {
        "slope": {"value": bw_state_D, "ci95": pct(draws["bw_state_D"])},
        "ratio": {"value": r_state_D, "ci95": pct(draws["r_state_D"])},
    },
    "week_only": {"slope": {"value": bw_week, "ci95": pct(draws["bw_week"])}},
    "state_and_week": {"slope": {"value": bw_sw, "ci95": pct(draws["bw_sw"])}},
    "B": B
}

P.put("fe_table_s11b", res)
P.num("fe_bw_state", bw_state, ci=pct(draws["state_b"]))
P.num("fe_r_mean", r_mean, ci=pct(draws["r_mean"]))
P.num("fe_paired_ratio", paired_r, ci=pct(draws["paired_ratio"]))

log("completed state FE analysis", r_mean, paired_r)
