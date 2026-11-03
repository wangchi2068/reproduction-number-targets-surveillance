"""s06: cross-pathogen weight-ratio criteria with a multiplicity correction (M7/D14-D15).

Ported from the original package experiment/original/code/r3_c_pathogens.py and r3_common.py
(functions load, growth_weeks, build_rows; q1 = a_single_1w / a_original, q2 = a_equal_3w / a_original,
each specification on its own rows with D >= 5 under its own weights; rows with a non-positive or missing
reporting share in t-2..t+1 dropped, as in r3_common.build_rows valid_h). Changes:
  * the state-cluster BCa bootstrap is replaced by the two-way (state x week) pigeonhole bootstrap,
    shared across the three weight specifications within a combination (B = 2000);
  * pure-reporting null ranges come from the vintage rerun of the growth-null simulation
    (vintage_2026-10-20/code/r4_a_growth_null.py, band "A_all": three generation intervals x three
    admission-delay settings x four reporting-noise levels); the constants of the submitted manuscript
    (RNG_Q1, RNG_Q2) are the fallback and are reported for comparison;
  * every pathogen x window x criterion is a one-sided bootstrap test; Holm (FWER) and
    Benjamini-Hochberg (FDR) adjustments are applied over the whole family;
  * a flu W-rule4 combination (the primary windows of this revision) is added.
Also records the status of FluSurv-NET (insufficient information) and the ILINet intervals of the vintage
rerun (state-cluster BCa; reused, not recomputed here).
"""
import json

import numpy as np
import pandas as pd

from cc_common import (EARLY_CAL, EXP, P1_END, P3_START, SEED, WEIGHTS, Log, Part, fit_stat, load_raw, pct,
                       windows_for, mask_windows)

log = Log("s06_pathogens")
P = Part("s06_pathogens", "s06_pathogens.py",
         ["data/external/nhsn_hrd_vintage_2026-10-20.csv.gz",
          "analysis/vintage_r4_results.json | vintage_2026-10-20/analysis/r4_results.json"],
         "measured", SEED + 6)
import sys
B = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
VR = 0.04
GROWTH_START, RSV_START, G_NAT_MIN = pd.Timestamp("2022-08-01"), pd.Timestamp("2024-11-01"), 1.05
RNG_Q = {"q1": (0.883, 1.338), "q2": (0.892, 1.373)}       # submitted manuscript, fallback only


def growth_weeks(df, pathogen):
    nat = df.groupby("week_end").C.sum().sort_index()
    c, wk = nat.to_numpy(), nat.index
    w = WEIGHTS["original"]
    t = np.arange(2, len(c) - 1)
    d = w[0] * c[t] + w[1] * c[t - 1] + w[2] * c[t - 2]
    with np.errstate(divide="ignore", invalid="ignore"):
        gn = np.where(d > 0, c[t + 1] / d, np.nan)
    wt = wk[t]
    ok = (gn >= G_NAT_MIN) & (wt >= GROWTH_START) & ~((wt > P1_END) & (wt < P3_START))
    if pathogen == "rsv":
        ok &= wt >= RSV_START
    return set(wt[ok])


def build_rows(df, w):
    w = np.asarray(w, float)
    out = []
    for s, z in df.groupby("Geographic aggregation", sort=True):
        if len(z) < 10:
            continue
        c, h, wk = z.C.to_numpy(), (z.coverage / 100.0).to_numpy(float), z.week_end.to_numpy()
        t = np.arange(2, len(c) - 1)
        d = w[0] * c[t] + w[1] * c[t - 1] + w[2] * c[t - 2]
        hw = np.column_stack([h[t - 2], h[t - 1], h[t], h[t + 1]])
        valid = np.all(np.isfinite(hw) & (hw > 0), axis=1)
        out.append(pd.DataFrame({"state": s, "week_end": wk[t], "D": d, "x": c[t + 1] / np.where(d > 0, d, np.nan),
                                 "valid": valid})[(d >= 5)])
    r = pd.concat(out, ignore_index=True)
    r["week_end"] = pd.to_datetime(r.week_end)
    y = r.week_end.dt.year - (r.week_end.dt.month < 8).astype(int)
    r["season_start"] = y
    return r[r.valid].reset_index(drop=True)


# ----------------------------------------------------------------------------- null bands
# The vintage growth-null bands are shipped under two possible layouts: the released
# package keeps them at analysis/vintage_r4_results.json, while a full working tree
# keeps them under the vintage folder. Probe both so a fresh clone runs unchanged.
R4_CANDIDATES = [
    EXP / "analysis/vintage_r4_results.json",
    EXP / "vintage_2026-10-20/analysis/r4_results.json",
]
R4_PATH = next((p for p in R4_CANDIDATES if p.exists()), None)
if R4_PATH is None:
    raise FileNotFoundError(
        "Growth-null bands not found. Expected one of:\n  "
        + "\n  ".join(str(p) for p in R4_CANDIDATES)
        + "\nRe-run the vintage r4 stage, or copy r4_results.json to analysis/vintage_r4_results.json."
    )
log(f"growth-null bands from {R4_PATH}")
r4 = json.loads(R4_PATH.read_text(encoding="utf-8"))
gn = r4.get("A_growth_null", {})
bands = gn.get("bands", {})
P.put("null_band_source", {"complete": gn.get("complete"), "reps_per_cell": gn.get("reps_per_cell"),
                           "seed": gn.get("seed"), "band": "A_all", "fallback": RNG_Q})

# ----------------------------------------------------------------------------- combinations
combos = []
for pg in ("flu", "covid", "rsv"):
    raw = load_raw(pg)
    gw = growth_weeks(raw, pg)
    rows = {k: build_rows(raw, w) for k, w in WEIGHTS.items()}
    cal = EARLY_CAL[2:] if pg == "rsv" else EARLY_CAL
    sel = {"W-cal": lambda r, cal=cal: mask_windows(r, cal),
           "W-growth": lambda r, gw=gw: r.week_end.isin(gw).to_numpy()}
    if pg == "flu":
        sel["W-growth-2025"] = lambda r, gw=gw: (r.week_end.isin(gw) & (r.season_start == 2025)).to_numpy()
        sel["W-rule4"] = lambda r: mask_windows(r, windows_for("rule4"))
    for wn, fn in sel.items():
        combos.append((f"{pg}|{wn}", {k: v[fn(v)].reset_index(drop=True) for k, v in rows.items()}))


def stats_of(R, wts=None):
    f = {k: fit_stat(R[k].x.to_numpy(), R[k].D.to_numpy(), None if wts is None else wts[k]) for k in R}
    a0 = f["original"]["a"]
    return np.array([f["single_1w"]["a"] / a0, f["equal_3w"]["a"] / a0, a0 - VR / f["original"]["xbar"] ** 2,
                     a0, f["single_1w"]["a"], f["equal_3w"]["a"]])


NAMES = ["q1", "q2", "excess_vR0.04", "a_original", "a_single_1w", "a_equal_3w"]
tests = []
for ci, (name, R) in enumerate(combos):
    if f"combo|{name}" in P.d["blocks"]:
        blk = P.d["blocks"][f"combo|{name}"]
    else:
        st = pd.Index(sorted(set().union(*[set(R[k].state) for k in R])))
        wk = pd.Index(sorted(set().union(*[set(R[k].week_end) for k in R])))
        sc = {k: st.get_indexer(R[k].state) for k in R}
        wc = {k: wk.get_indexer(R[k].week_end) for k in R}
        rng = np.random.default_rng(SEED + 600 + ci)
        pt = stats_of(R)
        dr = np.empty((B, len(NAMES)))
        for r in range(B):
            cs = np.bincount(rng.integers(0, len(st), len(st)), minlength=len(st)).astype(float)
            cw = np.bincount(rng.integers(0, len(wk), len(wk)), minlength=len(wk)).astype(float)
            dr[r] = stats_of(R, {k: cs[sc[k]] * cw[wc[k]] for k in R})
        bkey = "flu|W-cal" if name == "flu|W-rule4" else name     # rule windows were not simulated; flu W-cal band
        band = (bands.get(bkey) or {}).get("A_all")
        blk = {"N": {k: int(len(v)) for k, v in R.items()}, "n_states": int(len(st)), "n_weeks": int(len(wk)),
               "seasons": {str(k): int(v) for k, v in R["original"].season_start.value_counts().sort_index().items()},
               "point": dict(zip(NAMES, pt.tolist())), "ci95": {n: pct(dr[:, i]) for i, n in enumerate(NAMES)},
               "B": B, "seed": SEED + 600 + ci,
               "null_band": band and {"q1": band["q1"], "q2": band["q2"], "n_cells": band["n_cells"], "source": f"growth-null A_all ({bkey})"},
               "null_band_fallback": RNG_Q}
        bnd = blk["null_band"] or {"q1": list(RNG_Q["q1"]), "q2": list(RNG_Q["q2"])}
        p = lambda m: float((1 + np.sum(m)) / (B + 1))
        q1, q2, ex = dr[:, 0], dr[:, 1], dr[:, 2]
        blk["p"] = {
            "q1_lt_1": p(q1 >= 1), "q2_gt_1": p(q2 <= 1), "excess_gt_0": p(ex <= 0),
            "q1_outside_null": min(1.0, 2 * min(p(q1 <= bnd["q1"][1]), p(q1 >= bnd["q1"][0]))),
            "q2_outside_null": min(1.0, 2 * min(p(q2 <= bnd["q2"][1]), p(q2 >= bnd["q2"][0])))}
        blk["direction"] = {"q1_vs_null": "above" if pt[0] > bnd["q1"][1] else ("below" if pt[0] < bnd["q1"][0] else "inside"),
                            "q2_vs_null": "above" if pt[1] > bnd["q2"][1] else ("below" if pt[1] < bnd["q2"][0] else "inside")}
        P.put(f"combo|{name}", blk)
        log(name, blk["N"], {k: round(v, 4) for k, v in blk["point"].items()}, blk["p"], blk["direction"])
    for k, v in blk["p"].items():
        tests.append({"combo": name, "criterion": k, "p": v})

# ----------------------------------------------------------------------------- multiplicity
t = pd.DataFrame(tests)
m = len(t)
o = np.argsort(t.p.to_numpy())
ps = t.p.to_numpy()[o]
holm = np.minimum(1, np.maximum.accumulate((m - np.arange(m)) * ps))
bh = np.minimum(1, np.minimum.accumulate((m / np.arange(1, m + 1) * ps)[::-1])[::-1])
t.loc[t.index[o], "p_holm"] = holm
t.loc[t.index[o], "p_bh"] = bh
t["sig_holm_0.05"] = t.p_holm < 0.05
t["sig_bh_0.05"] = t.p_bh < 0.05
P.put("multiplicity", {"family_size": m, "tests": t.to_dict("records"),
                       "note": "one-sided bootstrap p = (1 + #draws on the null side)/(B + 1); outside-null = 2 x min one-sided"})
P.num("pathogen_family_size", m, mode="analytic")
P.num("pathogen_n_sig_raw_0.05", int((t.p < 0.05).sum()))
P.num("pathogen_n_sig_holm_0.05", int(t["sig_holm_0.05"].sum()))
P.num("pathogen_n_sig_bh_0.05", int(t["sig_bh_0.05"].sum()))
for name, _ in combos:
    blk = P.d["blocks"][f"combo|{name}"]
    for i, n in enumerate(("q1", "q2")):
        P.num(f"{name}_{n}", blk["point"][n], ci=blk["ci95"][n], seed=blk["seed"], B=B)
    P.num(f"{name}_excess_vR0.04", blk["point"]["excess_vR0.04"], ci=blk["ci95"]["excess_vR0.04"], seed=blk["seed"], B=B)
log(t.to_string())

# ----------------------------------------------------------------------------- other surveillance systems
P.put("flusurv_net", {"status": "insufficient information",
                      "reason": "FluSurv-NET covers about 14 states (catchment areas of ~9% of the US population) "
                                "and is not a census of admissions; the original comparison had too few state-weeks "
                                "with >= 5 cases in the early windows to estimate an intercept. Reported as "
                                "insufficient information rather than as a non-replication (revision_plan M7-4)."})
il = r4.get("I_independent_data", {}).get("ILINet", {})
if il:
    keep = {k: il["W-cal"]["stats"].get(k) for k in ("a_original", "rbar_original", "q1", "q2", "excess_scn")
            if k in il["W-cal"]["stats"]}
    P.put("ilinet", {"source": il.get("source"), "file": il.get("file"), "sha256": il.get("sha256"),
                     "N": il["W-cal"].get("N"), "n_states": il["W-cal"].get("n_units"), "stats": keep,
                     "interval": "state-cluster BCa, B = 2000 (vintage rerun r4 I; reused, not recomputed)",
                     "role": "supplementary, direction only"})
    for k, v in keep.items():
        if v:
            P.num(f"ilinet_{k}", v["point"], ci=v.get("bca95"), mode="measured",
                  note="vintage rerun r4_results.json I_independent_data; state-cluster BCa")
log.close()
