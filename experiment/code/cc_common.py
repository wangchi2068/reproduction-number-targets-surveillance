"""Shared definitions for the revised case-count-bottleneck analyses (step revise-experiments).

Provenance
----------
- Panel construction (jurisdiction filter, zero-filling of missing counts, D_t with weights
  (0.65, 0.25, 0.10), D_t >= 5, regime and season dates, calendar early windows) is ported from the
  original package `experiment/code/common.py` (wangchi2068/case-count-bottleneck, commit c2c859a),
  so that the continuity specification reproduces the original panel exactly.
- Changes relative to the original package (reviewer comments in brackets):
  * all paths are relative to this file; no hard-coded pipeline layout [M8/D16];
  * a renewal-equation statistic R_A is added next to the growth factor G [M2/D5, m12/D31];
  * a uniform, rule-based early-growth window is defined for all four seasons 2022-23 to 2025-26 [M3, m3];
  * FGLS accepts observation weights so that the state x week pigeonhole bootstrap can be run
    without duplicating rows [M5/D10-D11];
  * one interval type is used throughout: 95% percentile intervals from the two-way bootstrap [M5/D11].
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

CODE = Path(__file__).resolve().parent
EXP = CODE.parent
DATA = EXP / "data"
LOGS = EXP / "logs"
ANALYSIS = EXP / "analysis"
PARTS = ANALYSIS / "parts"
for _d in (DATA, LOGS, PARTS):
    _d.mkdir(parents=True, exist_ok=True)

RAW_CSV = Path(os.environ.get("NHSN_CSV", DATA / "external/nhsn_hrd_vintage_2026-10-20.csv.gz"))
PANEL = DATA / "panel_nhsn_vintage_2026-10-20.csv.gz"

SEED = 20261020                      # base seed of this step; each script offsets it
B_DEFAULT = 2000                     # bootstrap replicates for headline intervals

EXCLUDE = {"PR", "VI", "GU", "MP", "AS", "US"}
PATHOGEN_COLS = {"flu": "Total Influenza Admissions",
                 "covid": "Total COVID-19 Admissions",
                 "rsv": "Total RSV Admissions"}
COVERAGE_COLS = {"flu": "Percent Hospitals Reporting Influenza Admissions",
                 "covid": "Percent Hospitals Reporting COVID-19 Admissions",
                 "rsv": "Percent Hospitals Reporting RSV Admissions"}
W_ORIG = (0.65, 0.25, 0.10)          # original D_t weights (mean lag 1.45 weeks)
WEIGHTS = {"original": W_ORIG, "single_1w": (1.0, 0.0, 0.0), "equal_3w": (1 / 3, 1 / 3, 1 / 3)}
EARLY_CAL = [("2022-10-08", "2022-11-26"), ("2023-10-14", "2023-12-16"), ("2024-11-09", "2025-01-18")]
P1_END, P3_START = pd.Timestamp("2024-04-30"), pd.Timestamp("2024-11-01")
SEASON_MONTHS = {10, 11, 12, 1, 2, 3}
SEASONS = ["2022-23", "2023-24", "2024-25", "2025-26"]

# Influenza generation interval, gamma (mean, sd) in days. "mid" = serial interval of
# Cowling et al. (2009); short/long bracket the literature range used by the original package.
GI_FLU = {"short": (2.2, 0.98), "mid": (3.6, 1.6), "long": (4.3, 1.91)}
GROWTH_RULE = 1.05                   # national weekly growth threshold of the rule-based window


# ----------------------------------------------------------------------------- statistics
def wl_R(g, mean=GI_FLU["mid"][0], sd=GI_FLU["mid"][1]):
    """Euler-Lotka (Wallinga & Lipsitch 2007) map from a weekly growth factor g to R for a gamma
    generation interval: r = log(g)/7, R = (1 + r*scale)**shape. Delays between infection and
    admission shift but do not change r in exponential phases, so applied to admissions this is the
    apparent reproduction number R_A of Eales & Riley (2024)."""
    g = np.asarray(g, float)
    shape, scale = (mean / sd) ** 2, sd ** 2 / mean
    with np.errstate(divide="ignore", invalid="ignore"):
        base = 1.0 + scale * np.log(g) / 7.0
    return np.where(base > 0, np.power(np.clip(base, 0, None), shape), 0.0)


def fgls_affine(L, D, w=None, tol=1e-10, maxit=200):
    """Iterated FGLS for E(L|D) = a + b/D with Var(L|D) proportional to mu^2 (original package,
    common.fgls_affine), extended with optional non-negative observation weights w."""
    L = np.asarray(L, float); D = np.asarray(D, float)
    w = np.ones_like(L) if w is None else np.asarray(w, float)
    keep = w > 0
    L, D, w = L[keep], D[keep], w[keep]
    X = np.column_stack([np.ones_like(D), 1.0 / D])
    sw = np.sqrt(w)
    beta = np.linalg.lstsq(X * sw[:, None], L * sw, rcond=None)[0]
    for _ in range(maxit):
        mu = np.clip(X @ beta, 1e-6, None)
        ww = sw / mu
        new = np.linalg.lstsq(X * ww[:, None], L * ww, rcond=None)[0]
        done = np.max(np.abs(new - beta) / (np.abs(beta) + 1e-12)) < tol
        beta = new
        if done:
            break
    return float(beta[0]), float(beta[1])


def rel_loss(x, w=None):
    """L = ((x - xbar)/xbar)^2 with xbar the (weighted) window mean, re-estimated in every resample."""
    x = np.asarray(x, float)
    xb = np.average(x, weights=w) if w is not None else x.mean()
    return ((x - xb) / xb) ** 2, float(xb)


def fit_stat(x, D, w=None):
    L, xb = rel_loss(x, w)
    a, b = fgls_affine(L, D, w)
    return {"a": a, "b": b, "ba": b / a if a > 0 else np.nan, "xbar": xb}


# ----------------------------------------------------------------------------- bootstrap
def cluster_weights(groups, rng):
    """Multiplicity weights from resampling the unique values of `groups` with replacement."""
    codes, uniq = pd.factorize(groups)
    draw = rng.integers(0, len(uniq), len(uniq))
    cnt = np.bincount(draw, minlength=len(uniq)).astype(float)
    return cnt[codes]


def boot_weights(df, scheme, rng):
    """scheme: 'state' | 'week' | 'twoway' (pigeonhole: state count x week count, Owen 2007)
    | 'season_week' (weeks resampled within season, states fixed)."""
    if scheme == "state":
        return cluster_weights(df["state"].to_numpy(), rng)
    if scheme == "week":
        return cluster_weights(df["week_end"].to_numpy(), rng)
    if scheme == "twoway":
        return cluster_weights(df["state"].to_numpy(), rng) * cluster_weights(df["week_end"].to_numpy(), rng)
    if scheme == "season_week":
        w = np.zeros(len(df))
        for _, idx in df.groupby("season").indices.items():
            sub = df.iloc[idx]
            w[idx] = cluster_weights(sub["week_end"].to_numpy(), rng)
        return w
    raise ValueError(scheme)


def pct(x, q=(2.5, 97.5)):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return [float(v) for v in np.percentile(x, q)] if len(x) else [float("nan")] * len(q)


# ----------------------------------------------------------------------------- data
def load_raw(pathogen="flu") -> pd.DataFrame:
    col, cov = PATHOGEN_COLS[pathogen], COVERAGE_COLS[pathogen]
    df = pd.read_csv(RAW_CSV, usecols=["Week Ending Date", "Geographic aggregation", col, cov])
    g = df["Geographic aggregation"].astype(str)
    df = df[(g.str.len() == 2) & ~g.isin(EXCLUDE)].copy()
    df["week_end"] = pd.to_datetime(df["Week Ending Date"])
    df["missing"] = df[col].isna()
    df["C"] = df[col].fillna(0).astype(float)
    df["coverage"] = df[cov]
    return df.sort_values(["Geographic aggregation", "week_end"]).reset_index(drop=True)


def season_of(ts: pd.Series) -> pd.Series:
    y = ts.dt.year - (ts.dt.month < 8).astype(int)
    return y.astype(str) + "-" + ((y + 1) % 100).astype(str).str.zfill(2)


def build_panel(df: pd.DataFrame, dmin: float = 5.0) -> pd.DataFrame:
    """One row per state-week t with D_t (original weights) >= dmin, as in the original package,
    plus the single-week scale C_t, the growth factors and the renewal statistic R_A."""
    w = np.asarray(W_ORIG)
    out = []
    for state, sub in df.groupby("Geographic aggregation", sort=True):
        c = sub["C"].to_numpy(); miss = sub["missing"].to_numpy()
        d = sub["week_end"].to_numpy(); cov = sub["coverage"].to_numpy()
        n = len(c)
        if n < 10:
            continue
        t = np.arange(2, n - 1)
        dt = w[0] * c[t] + w[1] * c[t - 1] + w[2] * c[t - 2]
        rec = pd.DataFrame({
            "state": state, "week_end": d[t], "C_prev2": c[t - 2], "C_prev1": c[t - 1], "C_t": c[t],
            "C_next": c[t + 1], "D_t": dt, "coverage": cov[t],
            "coverage_next": cov[t + 1], "coverage_prev1": cov[t - 1], "coverage_prev2": cov[t - 2],
            "any_missing": miss[t - 2] | miss[t - 1] | miss[t] | miss[t + 1]})
        out.append(rec[dt >= dmin])
    p = pd.concat(out, ignore_index=True)
    p["week_end"] = pd.to_datetime(p["week_end"])
    p["G"] = p["C_next"] / p["D_t"]
    with np.errstate(divide="ignore", invalid="ignore"):
        p["G1"] = np.where(p["C_t"] > 0, p["C_next"] / p["C_t"], np.nan)
    gcc = (p["C_next"] + 0.5) / (p["C_t"] + 0.5)          # continuity-corrected single-week growth
    for k, (m, s) in GI_FLU.items():
        p[f"RA_{k}"] = wl_R(gcc, m, s)
    p["period"] = np.where(p.week_end <= P1_END, 1, np.where(p.week_end < P3_START, 2, 3))
    p["month"] = p.week_end.dt.month
    p["in_season"] = p.month.isin(SEASON_MONTHS)
    p["season"] = season_of(p.week_end)
    return p


def national_series(raw: pd.DataFrame) -> pd.DataFrame:
    nat = raw.groupby("week_end", as_index=False)["C"].sum().sort_values("week_end").reset_index(drop=True)
    nat["season"] = season_of(nat.week_end)
    nat["g"] = nat["C"].shift(-1) / nat["C"]
    return nat


# Earliest admissible week per season. 2024-25: first full week under the reinstated reporting
# mandate (week 3-9 Nov 2024); earlier weeks are voluntary-period reporting artefacts (regime 2).
RULE_START = {"2022-23": "2022-10-01", "2023-24": "2023-10-01", "2024-25": "2024-11-09", "2025-26": "2025-10-01"}


def rule_windows(nat: pd.DataFrame, seasons=SEASONS, thr=GROWTH_RULE):
    """Uniform early-growth window per season (rule v2, see analysis/data_decision.md):
    the first growth phase ends at the week before the first week (on/after RULE_START) whose
    national growth C(t+1)/C(t) falls below 1 (first national turning point P); the window is the
    contiguous run of weeks t <= P-1 with national growth >= thr, bounded below by RULE_START.
    Rows are selected by the week of D_t, so C(t+1) never exceeds the first turning point.
    Rule v1 (window ending at the season maximum) was dropped after it returned a 4-week window
    in 2024-25, whose national curve has two peaks; the change is disclosed in data_decision.md."""
    win, peaks = {}, {}
    for s in seasons:
        y = int(s[:4])
        sub = nat[(nat.week_end >= RULE_START[s]) & (nat.week_end <= f"{y + 1}-04-30")].reset_index(drop=True)
        ip = int(np.argmax(sub.g.to_numpy() < 1.0))        # first week whose next-week growth < 1
        peaks[s] = sub.week_end[ip]
        j = ip - 1
        while j - 1 >= 0 and sub.g[j - 1] >= thr:
            j -= 1
        win[s] = (sub.week_end[j], sub.week_end[ip - 1])
    return win, peaks


def mask_windows(p: pd.DataFrame, windows) -> np.ndarray:
    m = np.zeros(len(p), bool)
    for a, b in windows:
        m |= ((p.week_end >= pd.Timestamp(a)) & (p.week_end <= pd.Timestamp(b))).to_numpy()
    return m


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ----------------------------------------------------------------------------- results
class Part:
    """Incrementally written results block analysis/parts/<name>.json.

    num(): one headline number with provenance (script, input, seed, mode).
    put(): a detailed block (tables, curves) used by figures and tables."""

    def __init__(self, name: str, script: str, inputs, mode: str, seed=None):
        self.path = PARTS / f"{name}.json"
        self.meta = {"script": f"code/{script}", "inputs": inputs, "mode": mode, "seed": seed}
        self.d = json.loads(self.path.read_text()) if self.path.exists() else {}
        self.d["meta"] = self.meta
        self.d.setdefault("numbers", {})
        self.d.setdefault("blocks", {})

    def _save(self):
        self.d["meta"]["updated"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        tmp = self.path.with_suffix(".tmp")
        tmp.write_text(json.dumps(self.d, indent=1, default=_json_default, ensure_ascii=False))
        tmp.replace(self.path)

    def num(self, key, value, ci=None, unit=None, note=None, mode=None, seed="meta", **extra):
        rec = {"value": value}
        if ci is not None:
            rec["ci95"] = ci
        if unit:
            rec["unit"] = unit
        if note:
            rec["note"] = note
        if mode:
            rec["mode"] = mode
        if seed != "meta":
            rec["seed"] = seed
        rec.update(extra)
        self.d["numbers"][key] = rec
        self._save()

    def put(self, key, obj):
        self.d["blocks"][key] = obj
        self._save()


def _json_default(o):
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (pd.Timestamp,)):
        return o.strftime("%Y-%m-%d")
    raise TypeError(type(o))


class Log:
    def __init__(self, name: str):
        self.f = open(LOGS / f"{name}.log", "a", encoding="utf-8")
        self.t0 = time.time()
        self(f"==== {name} start {time.strftime('%Y-%m-%d %H:%M:%S')}")

    def __call__(self, *a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        self.f.write(s + "\n"); self.f.flush()

    def close(self):
        self(f"elapsed {time.time() - self.t0:.1f} s")
        self.f.close()


def load_panel() -> pd.DataFrame:
    p = pd.read_csv(PANEL, parse_dates=["week_end"])
    return p


# ----------------------------------------------------------------------------- specifications
SPECS = {
    # name: (window set, statistic column, scale column)
    "RA_rule4": ("rule4", "RA_mid", "C_t"),     # primary: renewal R_A, uniform rule windows, S = 4
    "G_rule4": ("rule4", "G", "D_t"),           # original statistic on the primary windows
    "G_cal3": ("cal3", "G", "D_t"),             # continuity with the submitted manuscript (S = 3)
    "RA_cal3": ("cal3", "RA_mid", "C_t"),
}


def windows_for(kind: str, nat: pd.DataFrame | None = None):
    if kind == "cal3":
        return EARLY_CAL
    if kind == "rule4":
        if nat is None:
            nat = pd.read_csv(DATA / "national_flu.csv", parse_dates=["week_end"])
        return list(rule_windows(nat)[0].values())
    raise ValueError(kind)


def spec_frame(panel: pd.DataFrame, name: str, nat=None) -> pd.DataFrame:
    """Rows of a specification with columns x (statistic) and D (scale >= 5)."""
    kind, xcol, dcol = SPECS[name]
    e = panel[mask_windows(panel, windows_for(kind, nat))].copy()
    e = e[e[dcol] >= 5]
    e["x"] = e[xcol].astype(float)
    e["D"] = e[dcol].astype(float)
    e = e[np.isfinite(e.x)].reset_index(drop=True)
    return e
