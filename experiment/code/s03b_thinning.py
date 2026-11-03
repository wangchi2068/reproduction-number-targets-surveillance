"""s03b_thinning.py: Binomial thinning experiment on the primary R_A sample (n = 1748).
Tests the 1/m scaling of the counting error term by thinning observed counts: C^(q) ~ Binomial(C, q).
Adheres strictly to the protocol:
- Fixed sample of 1748 state-weeks of RA_rule4.
- Single draw per unique (state, week_end) shared between C_t and C_next across adjacent weeks.
- Retention probabilities q in [0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1].
- Regression on original case volume D (which equals C_t for R_A).
- Two-way (pigeonhole) bootstrap (B = 1000) for paired confidence intervals.
- Sub-analyses for q >= 0.2 and C_t >= 50.
Writes analysis/parts/s11_thinning.json.
Usage: python3 s03b_thinning.py [B=1000]
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
    load_raw,
    pct,
    rel_loss,
    spec_frame,
    wl_R,
)

B = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
log = Log("s03b_thinning")
P = Part("s11_thinning", "s03b_thinning.py", ["data/panel_nhsn_vintage_2026-10-20.csv.gz"], "measured", SEED + 31)

# 1. Load data and setup sample
p = load_panel()
e = spec_frame(p, "RA_rule4")
raw = load_raw("flu")
raw["week_end"] = pd.to_datetime(raw["Week Ending Date"])
raw_dict = dict(zip(zip(raw["Geographic aggregation"], raw["week_end"]), raw["C"].fillna(0).astype(int)))

keys_t = list(zip(e.state, e.week_end))
keys_n = list(zip(e.state, e.week_end + pd.Timedelta(days=7)))
D = e.D.to_numpy()
Ct = e.C_t.to_numpy().astype(int)
Cn = e.C_next.to_numpy().astype(int)

# Unthinned fit
L1, _ = rel_loss(e.x.to_numpy())
a1, b1 = fgls_affine(L1, D)
ba1 = b1 / a1
log("unthinned", a1, b1, ba1)

# Poisson reference b_P: mean of kappa_A^2 * (1 + C_t / C_{t+1})
# where kappa_A = 1/2 approx for RA; empirical formula from methods
r_ratio = (Ct + 0.5) / (Cn + 0.5)
bP_all = float(np.mean(0.25 * (1.0 + r_ratio)))
bP_c50 = float(np.mean(0.25 * (1.0 + r_ratio[Ct >= 50])))

qs = [0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
R_INNER = 100  # inner replicates for point estimates
rng = np.random.default_rng(20261102)

log(f"running thinning point estimates across {len(qs)} levels (R={R_INNER})...")
point_levels = {}
for q in qs:
    bs, as_ = [], []
    for _ in range(R_INNER):
        thinned = {k: rng.binomial(v, q) for k, v in raw_dict.items()}
        ct = np.array([thinned[k] for k in keys_t])
        cn = np.array([thinned[k] for k in keys_n])
        xq = wl_R((cn + 0.5) / (ct + 0.5))
        Lq, _ = rel_loss(xq)
        aq, bq = fgls_affine(Lq, D)
        bs.append(bq)
        as_.append(aq)
    mean_bq = float(np.mean(bs))
    mean_aq = float(np.mean(as_))
    inv_q = 1.0 / q - 1.0
    cq = (mean_bq - b1) / inv_q
    da = (mean_aq - a1) * 1e3
    point_levels[str(q)] = {
        "q": q, "inv_q": inv_q, "b": mean_bq, "a": mean_aq,
        "cq": cq, "da_e3": da
    }

# Pooled slope across all q
xs = np.array([1.0 / q - 1.0 for q in qs])
dbs = np.array([point_levels[str(q)]["b"] - b1 for q in qs])
pooled_c = float(np.sum(xs * dbs) / np.sum(xs ** 2))

# Pooled slope excluding q=0.1
m_no01 = np.array(qs) >= 0.2
pooled_c_no01 = float(np.sum(xs[m_no01] * dbs[m_no01]) / np.sum(xs[m_no01] ** 2))
q01_leverage = float(xs[~m_no01][0] ** 2 / np.sum(xs ** 2))

log(f"point pooled_c={pooled_c:.3f}, excl_q0.1={pooled_c_no01:.3f}, q0.1 leverage={q01_leverage:.2%}")

# 2. Sub-analysis for C_t >= 50
sub50 = Ct >= 50
D50 = D[sub50]
keys_t50 = [keys_t[i] for i in range(len(keys_t)) if sub50[i]]
keys_n50 = [keys_n[i] for i in range(len(keys_n)) if sub50[i]]
L1_50, _ = rel_loss(e.x.to_numpy()[sub50])
a1_50, b1_50 = fgls_affine(L1_50, D50)

point_c50 = {}
for q in qs:
    bs = []
    for _ in range(R_INNER):
        thinned = {k: rng.binomial(v, q) for k, v in raw_dict.items()}
        ct = np.array([thinned[k] for k in keys_t50])
        cn = np.array([thinned[k] for k in keys_n50])
        xq = wl_R((cn + 0.5) / (ct + 0.5))
        Lq, _ = rel_loss(xq)
        aq, bq = fgls_affine(Lq, D50)
        bs.append(bq)
    mean_bq = float(np.mean(bs))
    point_c50[str(q)] = float((mean_bq - b1_50) / (1.0 / q - 1.0))

dbs_50 = np.array([point_c50[str(q)] * xs[i] for i, q in enumerate(qs)])
pooled_c_50 = float(np.sum(xs * dbs_50) / np.sum(xs ** 2))
log(f"point pooled_c_50={pooled_c_50:.3f}")

# 3. Two-way bootstrap for intervals
log(f"running two-way bootstrap (B={B})...")
boot_rng = np.random.default_rng(20261103)
boot_bq = {str(q): [] for q in qs}
boot_cq = {str(q): [] for q in qs}
boot_da = {str(q): [] for q in qs}
boot_c50_q = {str(q): [] for q in qs}
boot_pooled_c = []
boot_pooled_c_50 = []

for b_idx in range(B):
    w = boot_weights(e, "twoway", boot_rng)
    # unthinned under bootstrap weights
    L1_b, _ = rel_loss(e.x.to_numpy(), w=w)
    a1_b, b1_b = fgls_affine(L1_b, D, w=w)
    
    # C_t >= 50 under bootstrap weights
    w50 = w[sub50]
    L1_50_b, _ = rel_loss(e.x.to_numpy()[sub50], w=w50)
    a1_50_b, b1_50_b = fgls_affine(L1_50_b, D50, w=w50)
    
    # Thin one realization per q
    thinned = {k: boot_rng.binomial(v, 0.5) for k, v in raw_dict.items()}  # base draw
    dbs_b = []
    dbs_50_b = []
    
    for q in qs:
        # draw thinning for this level
        th = {k: boot_rng.binomial(v, q) for k, v in raw_dict.items()}
        ct = np.array([th[k] for k in keys_t])
        cn = np.array([th[k] for k in keys_n])
        xq = wl_R((cn + 0.5) / (ct + 0.5))
        Lq, _ = rel_loss(xq, w=w)
        aq, bq = fgls_affine(Lq, D, w=w)
        
        inv_q = 1.0 / q - 1.0
        cq_b = (bq - b1_b) / inv_q
        da_b = (aq - a1_b) * 1e3
        boot_bq[str(q)].append(bq)
        boot_cq[str(q)].append(cq_b)
        boot_da[str(q)].append(da_b)
        dbs_b.append(bq - b1_b)
        
        # sub 50
        ct50 = ct[sub50]
        cn50 = cn[sub50]
        xq50 = wl_R((cn50 + 0.5) / (ct50 + 0.5))
        Lq50, _ = rel_loss(xq50, w=w50)
        aq50, bq50 = fgls_affine(Lq50, D50, w=w50)
        c50_b = (bq50 - b1_50_b) / inv_q
        boot_c50_q[str(q)].append(c50_b)
        dbs_50_b.append(bq50 - b1_50_b)
        
    c_all_b = float(np.sum(xs * np.array(dbs_b)) / np.sum(xs ** 2))
    boot_pooled_c.append(c_all_b)
    c_50_b = float(np.sum(xs * np.array(dbs_50_b)) / np.sum(xs ** 2))
    boot_pooled_c_50.append(c_50_b)

# Summarize results
res_levels = {}
for q in qs:
    sq = str(q)
    res_levels[sq] = {
        "q": q,
        "b": {"value": point_levels[sq]["b"], "ci95": pct(boot_bq[sq])},
        "cq": {"value": point_levels[sq]["cq"], "ci95": pct(boot_cq[sq])},
        "da_e3": {"value": point_levels[sq]["da_e3"], "ci95": pct(boot_da[sq])},
        "c50": {"value": point_c50[sq], "ci95": pct(boot_c50_q[sq])}
    }

out_data = {
    "unthinned": {"a": a1, "b": b1, "ba": ba1},
    "levels": res_levels,
    "pooled_c": {"value": pooled_c, "ci95": pct(boot_pooled_c)},
    "pooled_c_excl_q01": pooled_c_no01,
    "q01_leverage_share": q01_leverage,
    "pooled_c_50": {"value": pooled_c_50, "ci95": pct(boot_pooled_c_50)},
    "bP_all": bP_all,
    "bP_c50": bP_c50,
    "c_over_bP": {"value": pooled_c / bP_all, "ci95": [v / bP_all for v in pct(boot_pooled_c)]},
    "B": B
}

# Put into part file
P.put("thinning_table_s11a", out_data)
P.num("thinning_pooled_c", pooled_c, ci=pct(boot_pooled_c), note="pooled slope rate c")
P.num("thinning_pooled_c_50", pooled_c_50, ci=pct(boot_pooled_c_50), note="pooled slope rate for Ct>=50")
P.num("thinning_c_over_bP", pooled_c / bP_all, ci=[v / bP_all for v in pct(boot_pooled_c)])

log("completed thinning analysis", pooled_c, pct(boot_pooled_c))
