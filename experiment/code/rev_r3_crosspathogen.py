"""Amendment 3, Section A3.4: cross-series replication on NHSN COVID-19 and RSV admissions.

Same frozen NHSN vintage, same 51 jurisdictions and week grid, same R_A transform with the same
generation interval (mean 3.6 d, SD 1.6 d). Window rule, fixed in the plan and applied without
exception: for each pathogen and season in RULE_START, the peak is the week of the maximum national
count between the season start and 30 April; the window is the maximal contiguous run of weeks ending
at that peak in which every national growth ratio g(t) = C_nat(t+1)/C_nat(t) is >= 1.0. A season is
included only if the run has >= 4 weeks and >= 45 of the 51 jurisdictions report a non-missing count
in >= 90% of the run's weeks. Excluded seasons are listed with the reason; no alternative window is
substituted. Sample filters after the window are those of the locked plan, Section 1.

Influenza is run through the same rule as an internal check and reported alongside.

Estimands per series (each with an interval):
  1. affine intercept a with the two-way cluster-robust 95% interval, and the free-exponent gamma;
  2. model-free floor: relative RMSE of R_A above the 90th percentile of C_t, calendar-week block
     bootstrap (B = 1,999);
  3. population coefficient a1 in model M2 of A3.1.

Usage: python code/rev_r3_crosspathogen.py [--B 1999]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))

import rev_common as R
import rev_r2_lib as L2
import revision_data as RD
from revision_data import load_sample
from revision_stats import fit_affine, relative_loss

OUT = R.AN / "rev_r3_crosspathogen.json"
HASHES = R.AN / "revision3" / "script_hashes.jsonl"
PATHOGENS = ("flu", "covid", "rsv")
MIN_WEEKS = 4
MIN_STATES = 45
MIN_REPORT_SHARE = 0.90


def window_rule(raw, season, start):
    """Peak week and the maximal contiguous pre-peak run with g(t) >= 1.0, plus the inclusion test."""
    year = int(season[:4])
    nat = raw.groupby("week_end").C.sum(min_count=1).sort_index()
    sub = nat[nat.index.to_series().between(pd.Timestamp(start), pd.Timestamp(f"{year+1}-04-30"))]
    if not len(sub) or sub.isna().all():
        return None, {"season": season, "included": False,
                      "reason": "no non-missing national counts in the season window"}
    g = sub.shift(-1) / sub.replace(0, np.nan)
    peak = sub.idxmax()
    j = int(np.flatnonzero(sub.index == peak)[0])
    while j > 0 and np.isfinite(g.iloc[j - 1]) and g.iloc[j - 1] >= 1.0:
        j -= 1
    run = sub.iloc[j:int(np.flatnonzero(sub.index == peak)[0]) + 1]
    first, last = run.index[0], run.index[-1]
    weeks = len(run)
    block = raw[raw.week_end.between(first, last)]
    share = block.groupby("state").C.apply(lambda s: float(s.notna().mean()))
    n_rep = int((share >= MIN_REPORT_SHARE).sum())
    if weeks < MIN_WEEKS:
        reason = f"run has {weeks} weeks (< {MIN_WEEKS})"
    elif n_rep < MIN_STATES:
        reason = (f"{n_rep} of {raw.state.nunique()} jurisdictions report in >= "
                  f"{MIN_REPORT_SHARE:.0%} of the run weeks (< {MIN_STATES})")
    else:
        reason = None
    meta = {"season": season, "included": reason is None, "reason": reason,
            "turning_week": peak.date().isoformat(),
            "window_start": first.date().isoformat(), "window_end": last.date().isoformat(),
            "weeks": weeks, "reporting_jurisdictions": n_rep}
    return (first, last), meta


def floor_rmse(e, B, rng):
    """Relative RMSE of R_A above the 90th percentile of C_t, calendar-week block bootstrap."""
    c = e.C_t.to_numpy(float); x = e.x.to_numpy(float); w = e.week_end.to_numpy()

    def stat(mask, weight):
        wt = weight[mask]; v = x[mask]
        if wt.sum() <= 0:
            return np.nan
        mean = float(np.average(v, weights=wt))
        return float(np.sqrt(np.average((v - mean) ** 2, weights=wt)) / mean)

    thr = float(np.quantile(c, .90))
    base_mask = c >= thr
    point = stat(base_mask, np.ones(len(c)))
    draws = []
    for _ in range(B):
        mult = R.week_multiplicity(w, rng)
        draws.append(stat(base_mask, mult))
    arr = np.asarray(draws, float)
    arr = arr[np.isfinite(arr)]
    return {"threshold_c_t": thr, "rows_above": int(base_mask.sum()), "relative_rmse": point,
            "ci": R.pct_ci(arr), "n": int(len(arr))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=1999)
    a = ap.parse_args()
    t0 = time.time()
    HASHES.parent.mkdir(parents=True, exist_ok=True)
    R.append_jsonl(HASHES, {"script": "rev_r3_crosspathogen.py",
                            "sha256": hashlib.sha256(open(__file__, "rb").read()).hexdigest(),
                            "utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),
                            "note": "Amendment 3, Section A3.4"})

    series = {}
    for pathogen in PATHOGENS:
        raw = RD.load_raw(pathogen)
        windows, ledger = {}, []
        for season, start in RD.RULE_START.items():
            win, meta = window_rule(raw, season, start)
            ledger.append(meta)
            if win is not None:
                windows[season] = win
        rec = {"window_ledger": ledger, "seasons_included": sorted(windows)}
        if not windows:
            rec["status"] = "no season passes the window rule"
            series[pathogen] = rec
            print(pathogen, "no season passes the window rule", flush=True)
            continue

        e = load_sample("RA", 3.6, pathogen=pathogen, raw=raw, windows=windows)
        e = e.merge(L2.population()[["state", "z_pop"]], on="state", how="left")
        assert e.z_pop.notna().all(), "population missing for some jurisdictions"
        rec["sample"] = {"rows": len(e), "states": int(e.state.nunique()),
                         "weeks": int(e.week_end.nunique()),
                         "count_range": [float(e.D.min()), float(e.D.max())]}

        from rev_r3_population import fit_model, inference, profile_three, state_multiplicity

        x, D = e.x.to_numpy(float), e.D.to_numpy(float)
        loss, xbar = relative_loss(x)
        st, wk = e.state.to_numpy(), e.week_end.to_numpy()
        z = e.z_pop.to_numpy(float)
        zs = e.groupby("state").z_pop.first().to_numpy()
        z10, z90 = (float(v) for v in np.quantile(zs, [.1, .9]))

        fa = fit_affine(loss, D)
        Xa = np.column_stack([np.ones(len(D)), 1 / D])
        beta_a = np.array([fa["a"], fa["b"]])
        V, _ = R.cgm_from_loss(x, xbar, loss, D, Xa, beta_a, st, wk)
        G = min(len(np.unique(st)), len(np.unique(wk)))
        tq = stats.t.ppf(.975, G - 1)
        se_a = float(np.sqrt(V[0, 0])) if V is not None else float("nan")
        a_affine = {"estimate": float(fa["a"]), "se": se_a,
                    "ci": [float(fa["a"] - tq * se_a), float(fa["a"] + tq * se_a)],
                    "b": float(fa["b"]), "b_ci": [float(fa["b"] - tq * np.sqrt(V[1, 1])),
                                                  float(fa["b"] + tq * np.sqrt(V[1, 1]))] if V is not None else None,
                    "df": int(G - 1)}
        m2 = fit_model(loss, z, D, True)
        r2 = inference(m2, x, xbar, loss, st, wk, z10, z90) if m2 is not None else None

        rng = np.random.default_rng([R.SEED, 91, PATHOGENS.index(pathogen)])
        fr = floor_rmse(e, a.B, rng)

        rec["affine_intercept"] = a_affine
        rec["free_exponent"] = ({"gamma": r2["gamma"], "a0": r2["a0"], "a1": r2["a1"]}
                                if r2 is not None else None)
        rec["model_free_floor"] = fr
        series[pathogen] = rec
        print(pathogen, "seasons", sorted(windows), "rows", len(e),
              "a", round(a_affine["estimate"], 6), "floor", round(fr["relative_rmse"], 4), flush=True)

    # replication criteria are judged against the influenza intervals
    ref = series.get("flu", {})
    ref_floor = (ref.get("model_free_floor") or {}).get("ci")
    verdicts = {}
    for pathogen, rec in series.items():
        if "affine_intercept" not in rec:
            verdicts[pathogen] = {"floor": "not assessed", "gradient": "not assessed"}
            continue
        a0 = rec["affine_intercept"]
        floor_ok = bool(a0["estimate"] > 0 and a0["ci"][0] > 0
                        and ref_floor is not None
                        and ref_floor[0] <= rec["model_free_floor"]["relative_rmse"] <= ref_floor[1])
        fe = rec.get("free_exponent")
        grad_ok = bool(fe is not None and fe["a1"]["estimate"] < 0 and fe["a1"]["ci"][1] < 0)
        verdicts[pathogen] = {"floor": "replicates" if floor_ok else "does not replicate",
                              "gradient": "replicates" if grad_ok else "does not replicate"}

    out = {
        "meta": {"script": "code/rev_r3_crosspathogen.py",
                 "plan": "analysis/analysis_plan_amendment3.md, Section A3.4",
                 "data": "frozen NHSN 2026-10-20 vintage; columns Total Influenza / COVID-19 / RSV Admissions",
                 "window_rule": {"min_weeks": MIN_WEEKS, "min_states": MIN_STATES,
                                 "min_report_share": MIN_REPORT_SHARE,
                                 "growth_rule": "national g(t) = C_nat(t+1)/C_nat(t) >= 1.0"},
                 "B": a.B, "seconds": round(time.time() - t0, 1)},
        "series": series,
        "replication": {"reference": "influenza floor 95% interval", "reference_floor_ci": ref_floor,
                        "reference_note": "influenza is run through the A3.4 window rule as an internal check, so its "
                                          "sample (1,799 rows) differs from the primary analysis (1,748 rows) and its "
                                          "rows here are not the primary result",
                        "verdicts": verdicts},
    }
    R.write_json(OUT, out)
    print("done", json.dumps(verdicts, ensure_ascii=False), "seconds", round(time.time() - t0, 1))


if __name__ == "__main__":
    main()
