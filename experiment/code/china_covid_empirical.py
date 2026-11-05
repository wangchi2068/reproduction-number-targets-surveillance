"""China COVID-19 cross-regime check (2020-2023), review M4.

The dispersion diagnostics of the main analysis applied to SARS-CoV-2
confirmed cases in mainland China, with the selection rule, generation
interval and units defined as in the main analysis rather than on the outcome:

  * Units: the 31 mainland provincial units of the JHU CSSE China extract. The
    "Unknown" row is not a geographic unit and is excluded. Hong Kong SAR and
    Macao SAR (separate surveillance systems) enter only as a sensitivity
    specification. Taiwan, China is not in the extract.
  * Windows: the national series is the sum of the mainland units. As in the
    US rule (revision_data.rule_windows), each national peak (g_nat >= 1 into
    the peak week, g_nat < 1 after it) with a peak weekly count >= 100 defines a
    growth window: the weeks before the peak, walking back while the national
    growth factor stays >= 1.05 (1.10 and 1.02 as sensitivities); windows
    shorter than two weeks are dropped. Rows are province-weeks inside the
    windows with C_t >= 5. The first-version rule (own-row growth g >= 1.05,
    which selects on C_next) is reported only as a comparison.
  * Generation interval by period: before 2022, the gamma interval of
    Ganyani et al. (2020) (mean 3.95 d, SD 1.51 d); from 2022 (Omicron), the
    within-household serial interval of Backer et al. (2022) for BA.1 (mean
    3.5 d, SD 2.4 d) as a proxy. The interval 5.2/1.7 d of the first version is a
    sensitivity specification. b_P uses each row's own interval.
  * Period: primary through the week ending 2022-12-03, before the end of
    routine mass testing on 7 December 2022; the full extract (to 2023-01-21) is
    a sensitivity specification.
  * Negative daily increments (downward revisions of the cumulative series)
    are set to zero (primary); the sensitivity drops every province-week that
    contains one.
  * Gamma quasi-score FGLS fit of the relative loss on 1/C_t; state x week
    pigeonhole bootstrap (all converged draws kept; the share with a <= 0 is
    reported and ratio intervals are flagged when it exceeds 2.5%).
"""
from __future__ import annotations

import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ma_models import boot_ci, write
from revision_data import wl_R
from revision_stats import bootstrap_weights, fit_loss

SEED = 2026110355
EXP = Path(__file__).resolve().parents[1]
RAW = EXP / "data/external/jhu_covid19_china_cases.csv"

GI_PERIOD = {"pre2022": (3.95, 1.51), "from2022": (3.5, 2.4)}
GI_FIRST_VERSION = (5.2, 1.7)
SARS = ("Hong Kong", "Macau")
MIN_COUNT, PEAK_MIN, MIN_WEEKS = 5, 100, 2
START, END_PRIMARY, END_FULL = "2020-01-25", "2022-12-03", "2023-01-21"
B_PRIMARY, B_SENS = 1000, 500


def weekly():
    df = pd.read_csv(RAW)
    id_vars = ["Province/State", "Country/Region", "Lat", "Long"]
    m = df.melt(id_vars=id_vars, var_name="date_str", value_name="cum")
    m["date"] = pd.to_datetime(m["date_str"], format="%m/%d/%y")
    m = m.sort_values(["Province/State", "date"]).reset_index(drop=True)
    m["daily_raw"] = m.groupby("Province/State")["cum"].diff().fillna(0)
    m["neg"] = m["daily_raw"] < 0
    m["daily"] = np.maximum(m["daily_raw"], 0)
    m["week_end"] = m["date"].dt.to_period("W-SAT").dt.end_time.dt.normalize()
    w = m.groupby(["Province/State", "week_end"], as_index=False).agg(C=("daily", "sum"), neg=("neg", "any"), days=("daily", "size"))
    w = w.rename(columns={"Province/State": "state"})
    return w[(w.week_end >= START) & (w.week_end <= END_FULL) & (w.days == 7)].reset_index(drop=True)


def panel(w):
    w = w.sort_values(["state", "week_end"]).reset_index(drop=True)
    grp = w.groupby("state")
    w["C_t"] = w.C.astype(float)
    w["C_next"] = grp.C.shift(-1)
    w["neg_next"] = grp.neg.shift(-1)
    w["gap_ok"] = (grp.week_end.shift(-1) - w.week_end).dt.days == 7
    w = w[w.C_next.notna() & w.gap_ok].copy()
    w["neg_any"] = w.neg | w.neg_next.astype(bool)
    w["g_own"] = (w.C_next + .5) / (w.C_t + .5)
    w["period"] = np.where(w.week_end < "2022-01-01", "pre2022", "from2022")
    return w.reset_index(drop=True)


def national(w):
    nat = w[~w.state.isin(SARS + ("Unknown",))].groupby("week_end", as_index=False).C.sum().sort_values("week_end")
    nat["g"] = nat.C.shift(-1) / nat.C
    return nat.reset_index(drop=True)


def growth_windows(nat, threshold, end):
    """US-rule analogue applied to every national peak (see module docstring)."""
    sub = nat[nat.week_end <= end].reset_index(drop=True)
    g, C = sub.g.to_numpy(), sub.C.to_numpy(float)
    out = []
    for ip in range(1, len(sub)):
        if not (g[ip - 1] >= 1 and (g[ip] < 1 or not np.isfinite(g[ip]))) or C[ip] < PEAK_MIN:
            continue
        j = ip - 1
        while j > 0 and g[j - 1] >= threshold:
            j -= 1
        if ip - j >= MIN_WEEKS:
            out.append({"first": sub.week_end[j].date().isoformat(), "last": sub.week_end[ip - 1].date().isoformat(),
                        "peak": sub.week_end[ip].date().isoformat(), "peak_count": int(C[ip]), "weeks": int(ip - j)})
    return out


def in_windows(p, windows):
    m = np.zeros(len(p), dtype=bool)
    for wd in windows:
        m |= p.week_end.between(pd.Timestamp(wd["first"]), pd.Timestamp(wd["last"])).to_numpy()
    return m


def gi_columns(sub, gi):
    if gi == "period":
        mean = sub.period.map({k: v[0] for k, v in GI_PERIOD.items()}).to_numpy(float)
        sd = sub.period.map({k: v[1] for k, v in GI_PERIOD.items()}).to_numpy(float)
    else:
        mean, sd = np.full(len(sub), gi[0]), np.full(len(sub), gi[1])
    return mean, sd


def bp_rows(ct, cn, mean, sd, w=None):
    """b_P as in s03b_thinning.poisson_reference, with each row's own generation interval."""
    g = (cn + .5) / (ct + .5)
    kap = (mean / 7) / (1 + (sd ** 2 / mean) * np.log(g) / 7)
    v = kap ** 2 * (1 + 1 / g)
    return float(np.average(v, weights=None if w is None else w))


def summarise(sub, gi, B, tag):
    sub = sub.reset_index(drop=True)
    mean, sd = gi_columns(sub, gi)
    ct, cn = sub.C_t.to_numpy(float), sub.C_next.to_numpy(float)
    x = np.array([float(wl_R(g, m, s)) for g, m, s in zip(sub.g_own, mean, sd)])
    D = ct
    point = fit_loss(x, D)
    bP = bp_rows(ct, cn, mean, sd)
    rng = np.random.default_rng([SEED, zlib.crc32(tag.encode()) % 100000])
    dr = np.full((B, 4), np.nan)
    for i in range(B):
        w = bootstrap_weights(sub, rng)[0]
        f = fit_loss(x, D, w=w)
        if f["converged"]:
            dr[i] = f["a"], f["b"], (f["b"] / f["a"] if f["a"] > 0 else np.nan), f["b"] / bp_rows(ct, cn, mean, sd, w)
    pa = point["a"] if point["converged"] else np.nan
    pb = point["b"] if point["converged"] else np.nan
    return {"N": int(len(sub)), "provinces": int(sub.state.nunique()), "weeks": int(sub.week_end.nunique()),
            "C_t_median": float(np.median(ct)), "converged": bool(point["converged"]),
            "a": boot_ci(pa, dr[:, 0], ratio=True), "b": boot_ci(pb, dr[:, 1]),
            "b_over_a": boot_ci(pb / pa if pa > 0 else np.nan, dr[:, 2]), "bP": bP,
            "b_over_bP": boot_ci(pb / bP, dr[:, 3]), "bootstrap_draws": B,
            "a_nonpositive_draw_share": float(np.mean(dr[np.isfinite(dr[:, 0]), 0] <= 0))}


def main():
    w = weekly()
    p = panel(w)
    nat = national(w)
    mainland = ~p.state.isin(SARS + ("Unknown",))
    base = (p.C_t >= MIN_COUNT).to_numpy()
    win = {thr: {end: growth_windows(nat, thr, end) for end in (END_PRIMARY, END_FULL)} for thr in (1.02, 1.05, 1.10)}
    upto = lambda end: (p.week_end <= end).to_numpy()
    W = lambda thr, end: in_windows(p, win[thr][end])
    prim = mainland.to_numpy() & base & upto(END_PRIMARY) & W(1.05, END_PRIMARY)
    specs = {
        "primary": (prim, "period", B_PRIMARY),
        "full_period_to_2023-01-21": (mainland.to_numpy() & base & upto(END_FULL) & W(1.05, END_FULL), "period", B_SENS),
        "window_threshold_1.10": (mainland.to_numpy() & base & upto(END_PRIMARY) & W(1.10, END_PRIMARY), "period", B_SENS),
        "window_threshold_1.02": (mainland.to_numpy() & base & upto(END_PRIMARY) & W(1.02, END_PRIMARY), "period", B_SENS),
        "with_Hong_Kong_SAR_and_Macao_SAR": (base & upto(END_PRIMARY) & W(1.05, END_PRIMARY) & (p.state != "Unknown").to_numpy(), "period", B_SENS),
        "gi_5.2_1.7_first_version": (prim, GI_FIRST_VERSION, B_SENS),
        "drop_negative_revision_weeks": (prim & ~p.neg_any.to_numpy(), "period", B_SENS),
        "min_count_20": (prim & (p.C_t >= 20).to_numpy(), "period", B_SENS),
        "pre_2022": (prim & (p.period == "pre2022").to_numpy(), "period", B_SENS),
        "from_2022": (prim & (p.period == "from2022").to_numpy(), "period", B_SENS),
        "own_row_growth_rule_first_version": (mainland.to_numpy() & base & upto(END_PRIMARY) & (p.g_own >= 1.05).to_numpy(), "period", B_SENS),
    }
    res = {"settings": {
        "source": "JHU CSSE COVID-19 time_series_covid19_confirmed_global, China rows (Dong, Du and Gardner 2020)",
        "units_primary": "31 mainland provincial units; 'Unknown' excluded; Hong Kong SAR and Macao SAR only in a sensitivity row; "
                         "Taiwan, China is not in the extract",
        "period_primary": [START, END_PRIMARY], "period_full": [START, END_FULL],
        "generation_interval": {"pre2022": {"mean": 3.95, "sd": 1.51, "source": "Ganyani et al. 2020, Euro Surveill 25(17):2000257"},
                                "from2022": {"mean": 3.5, "sd": 2.4, "source": "Backer et al. 2022, Euro Surveill 27(6):2200042 "
                                                                               "(BA.1 within-household serial interval, proxy)"},
                                "sensitivity": {"mean": GI_FIRST_VERSION[0], "sd": GI_FIRST_VERSION[1]}},
        "window_rule": f"national peaks with weekly count >= {PEAK_MIN}; growth weeks before each peak walking back while national g >= threshold; "
                       f">= {MIN_WEEKS} weeks", "min_count": MIN_COUNT, "negative_increments": "set to zero (primary)",
        "seed": SEED},
        "panel": {"province_weeks": int(len(p)), "units": int(p.state.nunique()),
                  "negative_revision_province_weeks_mainland": int((p.neg_any & mainland).sum())},
        "windows": {f"{thr:g}": {"to_" + end: v for end, v in d.items()} for thr, d in win.items()},
        "specs": {}}
    for name, (mask, gi, B) in specs.items():
        r = summarise(p[mask], gi, B, name)
        res["specs"][name] = r
        print(f"{name}: N={r['N']} prov={r['provinces']} a={r['a']['value']} {r['a']['ci95']} b/a={r['b_over_a']['value']} "
              f"b/bP={r['b_over_bP']['value']} {r['b_over_bP']['ci95']}", flush=True)
    p.to_csv(EXP / "data/panel_china_covid19_2020_2023.csv", index=False)
    write(EXP / "analysis/china_covid19_empirical.json", res)
    print("written analysis/china_covid19_empirical.json")


if __name__ == "__main__":
    main()
