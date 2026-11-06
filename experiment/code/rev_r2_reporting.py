"""External reporting data: hospital coverage and backfill (analysis_plan_amendment2.md, Section D.2; review question 6).

1. Coverage adjustment: counts C * 100 / coverage (NHSN percentage of hospitals reporting, per state-week, for
   both weeks of a transition); x recomputed from the adjusted counts; affine gamma-QL refit on the same rows with
   the observed C_t as the count (counting noise is set by the reported count); a_adj - a with a calendar-week block
   bootstrap interval of the difference (B = 1,999, seed [20261105, 72, 1]).
2. Backfill: first-issue and latest-issue values of confirmed_admissions_flu_ew from the Delphi Epidata API (nhsn
   source), one request per state, raw responses stored in data/external/delphi_nhsn/ with sha256. The Delphi
   archive starts at issue 2024w47, so a week counts as covered only if its first archived issue has lag <= 2 weeks
   (a provisional value); a row needs both of its weeks covered. Reported: share of covered rows whose latest-issue
   counts equal the vintage counts; v_rev = mean(((x_first - x_latest) / xbar_latest)^2); intercepts on first-issue
   and latest-issue counts and their difference with a calendar-week block bootstrap (B = 1,999, seed [.., 72, 2]).
   Descriptive: it bounds the backfill component, not the residual error in the settled vintage.

Output: analysis/rev_r2_reporting.json. Usage: python code/rev_r2_reporting.py [--B 1999] [--refetch]
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import time
import urllib.request

import numpy as np
import pandas as pd

import ma_calibration as MC
import rev_common as R
import rev_r2_lib as L2
from ma_models import Log
from revision_data import load_sample
from revision_stats import fit_affine, relative_loss

OUT = R.AN / "rev_r2_reporting.json"
RAW = R.EXP / "data" / "external" / "delphi_nhsn"
API = ("https://api.delphi.cmu.edu/epidata/covidcast/?data_source=nhsn&signals=confirmed_admissions_flu_ew"
       "&time_type=week&geo_type=state&geo_values={geo}&time_values=202235-202601&issues=202001-202652")
MAX_LAG = 2


def epiweek(sat):
    wed = pd.Timestamp(sat) - pd.Timedelta(days=3)
    return wed.year * 100 + (wed.dayofyear - 1) // 7 + 1


def fetch(states, refetch, log):
    RAW.mkdir(parents=True, exist_ok=True)
    manifest, fetched = {}, 0
    for s in states:
        p = RAW / f"delphi_nhsn_flu_{s.lower()}.json.gz"
        if refetch or not p.exists():
            fetched += 1
            url = API.format(geo=s.lower())
            for attempt in range(2):
                try:
                    with urllib.request.urlopen(url, timeout=120) as r:
                        body = r.read()
                    break
                except Exception as exc:  # recorded, second attempt only
                    log("fetch failed", s, attempt, exc)
                    body = None
            if body is None:
                continue
            p.write_bytes(gzip.compress(body, mtime=0))
        manifest[p.name] = {"state": s, "sha256": hashlib.sha256(p.read_bytes()).hexdigest(), "url": API.format(geo=s.lower())}
    if not fetched and (RAW / "MANIFEST.json").exists():  # frozen inputs: keep the manifest byte-identical
        return manifest
    (RAW / "MANIFEST.json").write_text(json.dumps({"fetched_utc_of_this_manifest": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),
                                                   "files": manifest}, indent=1), encoding="utf-8")
    return manifest


def issues_table(states):
    recs = []
    for s in states:
        p = RAW / f"delphi_nhsn_flu_{s.lower()}.json.gz"
        if not p.exists():
            continue
        js = json.loads(gzip.decompress(p.read_bytes()))
        for e in js.get("epidata") or []:
            if e.get("value") is None:
                continue
            recs.append((s, int(e["time_value"]), int(e["issue"]), int(e["lag"]), float(e["value"])))
    t = pd.DataFrame(recs, columns=["state", "ew", "issue", "lag", "value"]).sort_values(["state", "ew", "issue"])
    first = t.groupby(["state", "ew"]).first().rename(columns={"issue": "issue_first", "lag": "lag_first", "value": "first"})
    last = t.groupby(["state", "ew"]).last()[["issue", "value"]].rename(columns={"issue": "issue_last", "value": "latest"})
    n_iss = t.groupby(["state", "ew"]).size().rename("n_issues")
    return first.join(last).join(n_iss).reset_index(), len(t)


def two_fit(x1, x2, D1, D2, w):
    l1, _ = relative_loss(x1, w)
    l2, _ = relative_loss(x2, w)
    f1, f2 = fit_affine(l1, D1, w), fit_affine(l2, D2, w)
    if not (f1["converged"] and f2["converged"]):
        return np.nan, np.nan, np.nan
    return f1["a"], f2["a"], f2["a"] - f1["a"]


def boot_diff(x1, x2, D1, D2, week, B, seed):
    rng = np.random.default_rng([R.SEED, 72, seed])
    d = np.full((B, 3), np.nan)
    for r in range(B):
        d[r] = two_fit(x1, x2, D1, D2, R.week_multiplicity(week, rng))
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, default=1999)
    ap.add_argument("--refetch", action="store_true")
    args = ap.parse_args()
    L2.register_script(__file__)
    t0 = time.time()
    log = Log(R.LOGS / "rev_r2_reporting.log")
    f = load_sample("RA", 3.6)
    x, D = f.x.to_numpy(float), f.D.to_numpy(float)
    ct, cn = f.C_t.to_numpy(float), f.C_next.to_numpy(float)
    wk, st = f.week_end.to_numpy(), f.state.to_numpy()
    x_chk = R.wl((cn + .5) / (ct + .5))
    assert np.allclose(x_chk, x), "R_A transform mismatch"
    res = {"meta": {"script": "code/rev_r2_reporting.py", "plan": "analysis/analysis_plan_amendment2.md, Section D.2",
                    "B": args.B, "seeds": {"coverage": [R.SEED, 72, 1], "backfill": [R.SEED, 72, 2]}}}
    # ---- 1. coverage adjustment
    cov, covn = f.coverage.to_numpy(float), f.coverage_next.to_numpy(float)
    x_adj = R.wl((cn * 100 / covn + .5) / (ct * 100 / cov + .5))
    base = MC.sandwich(x, D, st, wk)
    adj = MC.sandwich(x_adj, D, st, wk)
    se_adj = float(np.sqrt(adj["V2"][0, 0]))
    bd = boot_diff(x, x_adj, D, D, wk, args.B, 1)
    res["coverage_adjustment"] = {
        "coverage_pct": {"mean": float(cov.mean()), "min": float(cov.min()), "p10": float(np.quantile(cov, .1)),
                         "share_rows_100": float(np.mean((cov >= 99.95) & (covn >= 99.95))),
                         "mean_abs_week_to_week_change_pct_points": float(np.mean(np.abs(covn - cov)))},
        "a": float(base["beta"][0]), "a_adjusted": float(adj["beta"][0]), "b_adjusted": float(adj["beta"][1]),
        "a_adjusted_se_two_way": se_adj, "difference": float(adj["beta"][0] - base["beta"][0]),
        "difference_ci_week_bootstrap": R.pct_ci(bd[:, 2]), "relative_difference": float(adj["beta"][0] / base["beta"][0] - 1),
        "count_used": "observed C_t (counting noise set by the reported count)",
        "assumption": "proportional participation: unreported hospitals have the same admissions per hospital as reporting ones"}
    # ---- 2. backfill
    states = sorted(np.unique(st))
    man = fetch(states, args.refetch, log)
    it, nrec = issues_table(states)
    it = it.set_index(["state", "ew"])
    ew_t = np.array([epiweek(s) for s in wk])
    ew_n = np.array([epiweek(pd.Timestamp(s) + pd.Timedelta(days=7)) for s in wk])

    def look(ews, col):
        idx = pd.MultiIndex.from_arrays([st, ews])
        return it.reindex(idx)[col].to_numpy()
    lag_t, lag_n = look(ew_t, "lag_first"), look(ew_n, "lag_first")
    covd = (np.nan_to_num(lag_t, nan=99) <= MAX_LAG) & (np.nan_to_num(lag_n, nan=99) <= MAX_LAG)
    in_archive = np.isfinite(look(ew_t, "latest")) & np.isfinite(look(ew_n, "latest"))
    ft, fn, lt, ln = (look(ew_t, "first")[covd], look(ew_n, "first")[covd], look(ew_t, "latest")[covd], look(ew_n, "latest")[covd])
    bf = {"api": API, "raw_files": len(man), "records": nrec, "max_lag_first_issue_weeks": MAX_LAG,
          "rows_in_archive": int(in_archive.sum()),
          "latest_equals_vintage_share_all_archive_rows": float(np.mean((look(ew_t, "latest")[in_archive] == ct[in_archive]) &
                                                                       (look(ew_n, "latest")[in_archive] == cn[in_archive]))),
          "rows_covered": int(covd.sum()), "weeks_covered": int(len(np.unique(wk[covd]))),
          "covered_week_range": [str(pd.Timestamp(wk[covd].min()).date()), str(pd.Timestamp(wk[covd].max()).date())] if covd.any() else None}
    if covd.sum() > 50:
        ct_c, cn_c = ct[covd], cn[covd]
        bf["latest_equals_vintage_share"] = float(np.mean((lt == ct_c) & (ln == cn_c)))
        bf["first_equals_latest_share_current_week"] = float(np.mean(ft == lt))
        bf["mean_relative_revision_current_week"] = float(np.mean((lt - ft) / np.maximum(ft, 1)))
        bf["share_revised_up"] = float(np.mean(lt > ft))
        bf["share_revised_down"] = float(np.mean(lt < ft))
        keep = (ft >= 5) & (lt >= 5)
        xf = R.wl((fn[keep] + .5) / (ft[keep] + .5))
        xl = R.wl((ln[keep] + .5) / (lt[keep] + .5))
        bf["rows_fitted"] = int(keep.sum())
        bf["v_rev"] = float(np.mean(((xf - xl) / xl.mean()) ** 2))
        wkc = wk[covd][keep]
        stc = st[covd][keep]
        sf = MC.sandwich(xf, ft[keep], stc, wkc)
        sl = MC.sandwich(xl, lt[keep], stc, wkc)
        bb = boot_diff(xl, xf, lt[keep], ft[keep], wkc, args.B, 2)
        bf.update({"a_first": float(sf["beta"][0]) if sf else None, "a_latest": float(sl["beta"][0]) if sl else None,
                   "a_first_minus_latest": float(sf["beta"][0] - sl["beta"][0]) if sf and sl else None,
                   "a_first_minus_latest_ci_week_bootstrap": R.pct_ci(bb[:, 2]),
                   "a_first_ci_week_bootstrap": R.pct_ci(bb[:, 1]), "a_latest_ci_week_bootstrap": R.pct_ci(bb[:, 0]),
                   "v_rev_over_a_latest": float(np.mean(((xf - xl) / xl.mean()) ** 2) / sl["beta"][0]) if sl else None,
                   "limitation": "bounds the backfill (provisional vs settled) component only, not the residual error of the settled vintage"})
    res["backfill"] = bf
    res["meta"]["seconds"] = round(time.time() - t0, 1)
    R.write_json(OUT, res)
    log("done", res["coverage_adjustment"]["difference"], bf.get("a_first_minus_latest"))


if __name__ == "__main__":
    main()
