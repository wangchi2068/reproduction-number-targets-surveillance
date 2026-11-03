"""Weekly-count inputs and fixed, disclosed analysis windows; importing writes nothing."""
from __future__ import annotations

from functools import lru_cache
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

EXP = Path(__file__).resolve().parents[1]
RAW = EXP / "data/external/nhsn_hrd_vintage_2026-10-20.csv.gz"
STATES = tuple("AL AK AZ AR CA CO CT DE DC FL GA HI ID IL IN IA KS KY LA ME MD MA MI MN MS MO MT NE NV NH NJ NM NY NC ND OH OK OR PA RI SC SD TN TX UT VT VA WA WV WI WY".split())
COUNT = {"flu": "Total Influenza Admissions", "covid": "Total COVID-19 Admissions", "rsv": "Total RSV Admissions"}
COVERAGE = {k: f"Percent Hospitals Reporting {v} Admissions" for k, v in (("flu", "Influenza"), ("covid", "COVID-19"), ("rsv", "RSV"))}
RULE_START = {"2022-23": "2022-10-01", "2023-24": "2023-10-01", "2024-25": "2024-11-09", "2025-26": "2025-10-01"}
WEIGHTS = np.array([.65, .25, .10])
CALENDAR = [("2022-10-08", "2022-11-26"), ("2023-10-14", "2023-12-16"), ("2024-11-09", "2025-01-18")]


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def season_of(dates):
    year = dates.dt.year - (dates.dt.month < 8).astype(int)
    return year.astype(str) + "-" + ((year + 1) % 100).astype(str).str.zfill(2)


def wl_R(g, mean=3.6, sd=None):
    sd = mean * 1.6 / 3.6 if sd is None else sd
    shape, scale = (mean / sd) ** 2, sd ** 2 / mean
    g = np.asarray(g, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        base = 1 + scale * np.log(g) / 7
        result = np.power(np.maximum(base, 0), shape)
    return np.where((g > 0) & np.isfinite(g), result, np.nan)


def kappa_A(g, mean=3.6, sd=None):
    sd = mean * 1.6 / 3.6 if sd is None else sd
    base = 1 + (sd ** 2 / mean) * np.log(np.asarray(g, dtype=float)) / 7
    return np.divide(mean / 7, base, out=np.full_like(base, np.nan), where=base > 0)


@lru_cache(maxsize=3)
def _raw(pathogen):
    col, cov = COUNT[pathogen], COVERAGE[pathogen]
    original = pd.read_csv(RAW, usecols=["Week Ending Date", "Geographic aggregation", col, cov])
    original = original[original["Geographic aggregation"].isin(STATES)].copy()
    original = original.rename(columns={"Week Ending Date": "week_end", "Geographic aggregation": "state", col: "C", cov: "coverage"})
    original["week_end"] = pd.to_datetime(original.week_end)
    if original.duplicated(["state", "week_end"]).any():
        raise ValueError("Duplicate jurisdiction-week keys in raw snapshot")
    original["C"] = pd.to_numeric(original.C, errors="coerce")
    original.loc[original.C < 0, "C"] = np.nan
    valid = original.C.dropna().to_numpy()
    if not np.all(valid == np.floor(valid)):
        raise ValueError("Admissions must be integer counts")
    original["coverage"] = pd.to_numeric(original.coverage, errors="coerce")
    original.loc[~original.coverage.between(0, 100), "coverage"] = np.nan
    records = []
    for state, group in original.groupby("state", sort=True):
        group = group.sort_values("week_end").set_index("week_end")
        expected = pd.date_range(group.index.min(), group.index.max(), freq="7D")
        if not group.index.isin(expected).all():
            raise ValueError(f"Irregular reporting weekday in {state}")
        aligned = group.reindex(expected)
        aligned["state"] = state
        aligned["missing_week"] = ~expected.isin(group.index)
        aligned["missing_count"] = aligned.C.isna()
        aligned.index.name = "week_end"
        records.append(aligned.reset_index())
    result = pd.concat(records, ignore_index=True).sort_values(["state", "week_end"]).reset_index(drop=True)
    result["season"] = season_of(result.week_end)
    return result


def load_raw(pathogen="flu"):
    return _raw(pathogen).copy()


def national_series(raw=None):
    raw = load_raw() if raw is None else raw
    grouped = raw.groupby("week_end").C
    nat = pd.DataFrame({"C": grouped.sum(min_count=1), "states_observed": grouped.count()}).reset_index()
    nat = nat.sort_values("week_end").reset_index(drop=True)
    nat["season"] = season_of(nat.week_end)
    nat["g"] = nat.C.shift(-1) / nat.C
    return nat


def rule_windows(nat=None, threshold=1.05):
    nat = national_series() if nat is None else nat
    windows, peaks = {}, {}
    for season, start in RULE_START.items():
        year = int(season[:4])
        sub = nat[nat.week_end.between(pd.Timestamp(start), pd.Timestamp(f"{year+1}-04-30"))].reset_index(drop=True)
        turns = np.flatnonzero(sub.g.to_numpy() < 1)
        if not len(turns) or turns[0] < 1:
            raise ValueError(f"No eligible first growth phase for {season}")
        ip = int(turns[0]); j = ip - 1
        while j > 0 and sub.g.iloc[j-1] >= threshold:
            j -= 1
        windows[season] = (sub.week_end.iloc[j], sub.week_end.iloc[ip-1])
        peaks[season] = sub.week_end.iloc[ip]
    return windows, peaks


def build_panel(raw=None):
    raw = load_raw() if raw is None else raw
    g = raw.sort_values(["state", "week_end"]).reset_index(drop=True)
    groups = g.groupby("state", sort=False)
    lag1, lag2, next_count = groups.C.shift(1), groups.C.shift(2), groups.C.shift(-1)
    dates1, dates2, dates_next = groups.week_end.shift(1), groups.week_end.shift(2), groups.week_end.shift(-1)
    z = pd.DataFrame({"state": g.state, "week_end": g.week_end, "C_t": g.C,
                      "C_prev1": lag1, "C_prev2": lag2, "C_next": next_count,
                      "coverage": g.coverage, "coverage_next": groups.coverage.shift(-1)})
    consecutive = ((g.week_end-dates1).dt.days==7) & ((g.week_end-dates2).dt.days==14) & ((dates_next-g.week_end).dt.days==7)
    z["consecutive"] = consecutive
    z["any_missing"] = z[["C_prev2", "C_prev1", "C_t", "C_next"]].isna().any(axis=1) | ~consecutive
    z["D_t"] = .65*z.C_t + .25*z.C_prev1 + .10*z.C_prev2
    z["season"] = season_of(z.week_end)
    return z


def load_sample(stat="RA", gi_mean=3.6, window="rule", pathogen="flu", raw=None, windows=None, gi_sd=None, panel=None):
    panel = build_panel(load_raw(pathogen) if raw is None else raw) if panel is None else panel
    bounds = (windows.values() if isinstance(windows, dict) else windows) if windows is not None else (rule_windows()[0].values() if window == "rule" else CALENDAR)
    included = np.zeros(len(panel), dtype=bool)
    for first, last in bounds:
        included |= panel.week_end.between(pd.Timestamp(first), pd.Timestamp(last)).to_numpy()
    e = panel[included & ~panel.any_missing & (panel.D_t >= 5)].copy()
    if stat == "RA":
        e = e[e.C_t >= 5].copy()
        e["x"] = wl_R((e.C_next+.5)/(e.C_t+.5), mean=gi_mean, sd=gi_sd)
        e["D"] = e.C_t
    elif stat == "G":
        e["x"] = e.C_next/e.D_t
        e["D"] = e.D_t
    else:
        raise ValueError(stat)
    e = e[np.isfinite(e.x)].sort_values(["state", "week_end"]).reset_index(drop=True)
    return e


def audit_inputs(output_dir):
    out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
    raw, panel, nat = load_raw(), build_panel(), national_series()
    windows, peaks = rule_windows(nat)
    meta = {"snapshot_sha256": sha256(RAW), "source": "https://data.cdc.gov/\u0075a7e-t2fy", "snapshot_date": "2026-10-20",
            "raw_state_rows_after_grid": len(raw), "states": int(raw.state.nunique()),
            "missing_count_rows": int(raw.missing_count.sum()), "inserted_missing_weeks": int(raw.missing_week.sum()),
            "windows": {s: {"start": a.date().isoformat(), "end": b.date().isoformat(), "weeks": (b-a).days//7+1,
                              "turning_week": peaks[s].date().isoformat()} for s,(a,b) in windows.items()}, "samples": {}}
    for stat,gi,name in [("RA",3.6,"RA"),("G",3.6,"G"),("RA",3.2,"RA32")]:
        e=load_sample(stat,gi)
        included=np.zeros(len(panel),bool)
        for a,b in windows.values(): included |= panel.week_end.between(a,b).to_numpy()
        meta["samples"][name]={"rows":len(e),"states":int(e.state.nunique()),"weeks":int(e.week_end.nunique()),
                               "candidate_rows":int(included.sum()),"incomplete_rows_in_windows":int((included & panel.any_missing).sum()),
                               "zero_next_counts":int((e.C_next==0).sum()),"count_range":[float(e.D.min()),float(e.D.max())]}
        e.to_csv(out/f"sample_{name}.csv",index=False,float_format="%.12g")
        meta["samples"][name]["file_sha256"]=sha256(out/f"sample_{name}.csv")
    (out/'input_audit.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False),encoding='utf-8')
    return meta
