"""China COVID-19 cross-regime validation (2020-2023).

Independent replication of the dispersion-diagnostic analysis on a second
pathogen (SARS-CoV-2) and a second country (China), under a regime of intensive
containment rather than the near-stationary US influenza baseline.

Design choices that mirror the main analysis:
  * Unit of analysis: province x week, 34 provinces.
  * Growth-phase restriction: rows whose own week-to-week growth factor is
    g >= 1.05, the analogue of the main analysis' early-growth windows.
  * Wallinga-Lipsitch transform with the COVID-19 gamma generation interval
    (mean 5.2 d, SD 1.7 d).
  * Gamma quasi-score FGLS loss fit on the inverse count.
  * Pigeonhole bootstrap over provinces x weeks.
"""
from __future__ import annotations

import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from revision_stats import bootstrap_weights, fit_loss
from s03b_thinning import poisson_reference

SEED = 2026110355
EXP = Path(__file__).resolve().parents[1]

GI_MEAN, GI_SD = 5.2, 1.7
THETA = GI_SD ** 2 / GI_MEAN
ALPHA = GI_MEAN / THETA
MIN_COUNT = 5
GROWTH_CUT = 1.05


def build_panel() -> pd.DataFrame:
    raw = EXP / "data/external/jhu_covid19_china_cases.csv"
    df = pd.read_csv(raw)
    id_vars = ["Province/State", "Country/Region", "Lat", "Long"]
    date_cols = [c for c in df.columns if c not in id_vars]
    m = df.melt(id_vars=id_vars, value_vars=date_cols, var_name="date_str", value_name="cum_cases")
    m["date"] = pd.to_datetime(m["date_str"], format="%m/%d/%y")
    m = m.sort_values(["Province/State", "date"]).reset_index(drop=True)
    m["daily"] = m.groupby("Province/State")["cum_cases"].diff().fillna(0)
    m["daily"] = np.maximum(m["daily"], 0)
    m["week_end"] = m["date"].dt.to_period("W-SAT").dt.end_time.dt.normalize()

    w = m.groupby(["Province/State", "week_end"], as_index=False)["daily"].sum()
    w = w[(w["week_end"] >= "2020-01-25") & (w["week_end"] <= "2023-01-21")].copy()
    w = w.sort_values(["Province/State", "week_end"]).reset_index(drop=True)
    w["C_t"] = w["daily"]
    w["C_next"] = w.groupby("Province/State")["C_t"].shift(-1)
    w = w.dropna(subset=["C_next"]).copy()

    w["g_own"] = (w["C_next"] + 0.5) / (w["C_t"] + 0.5)
    w["R_A"] = (1.0 + THETA * np.log(w["g_own"]) / 7.0) ** ALPHA
    w["year"] = w["week_end"].dt.year
    w["state"] = w["Province/State"]
    return w


def summarise(sub: pd.DataFrame, B: int, tag: str) -> dict:
    x, D = sub["R_A"].to_numpy(), sub["C_t"].to_numpy()
    point = fit_loss(x, D)
    bP = float(poisson_reference(sub, mean=GI_MEAN, sd=GI_SD))

    rng = np.random.default_rng([SEED, zlib.crc32(tag.encode()) % 100000])
    a_d, b_d, r_d, bb_d = [], [], [], []
    for _ in range(B):
        w = bootstrap_weights(sub, rng)[0]
        f = fit_loss(x, D, w=w)
        if f["converged"] and f["a"] > 0:
            a_d.append(f["a"]); b_d.append(f["b"]); r_d.append(f["ratio"])
            bb_d.append(f["b"] / poisson_reference(sub, w=w, mean=GI_MEAN, sd=GI_SD))

    def ci(v):
        return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))] if v else None

    return {
        "N": int(len(sub)),
        "provinces": int(sub["state"].nunique()),
        "a": {"point": float(point["a"]), "ci95": ci(a_d)},
        "b": {"point": float(point["b"]), "ci95": ci(b_d)},
        "b_over_a": {"point": float(point["ratio"]), "ci95": ci(r_d)},
        "bP": bP,
        "b_over_bP": {"point": float(point["b"] / bP), "ci95": ci(bb_d)},
        "converged": bool(point["converged"]),
    }


def main() -> None:
    panel = build_panel()
    panel.to_csv(EXP / "data/panel_china_covid19_2020_2023.csv", index=False)
    print(f"panel rows: {len(panel)}; provinces: {panel['state'].nunique()}")

    growth = panel[(panel["g_own"] >= GROWTH_CUT) & (panel["C_t"] >= MIN_COUNT)].reset_index(drop=True)
    print(f"growth-phase sample: {len(growth)}")

    res = {
        "settings": {
            "source": "JHU CSSE time_series_covid19_confirmed_global (China subset)",
            "period": ["2020-01-25", "2023-01-21"],
            "generation_interval_days": {"mean": GI_MEAN, "sd": GI_SD},
            "growth_cut": GROWTH_CUT,
            "min_count": MIN_COUNT,
            "seed": SEED,
        },
        "primary": summarise(growth, 1000, "primary"),
    }
    print("primary:", json.dumps(res["primary"]["a"]), json.dumps(res["primary"]["b_over_a"]))

    res["by_year"] = {str(y): summarise(growth[growth["year"] == y].reset_index(drop=True), 500, f"y{y}")
                      for y in [2020, 2021, 2022]}
    for y, r in res["by_year"].items():
        print(f"  {y}: N={r['N']} a={r['a']['point']:.4f} b/a={r['b_over_a']['point']:.1f} b/bP={r['b_over_bP']['point']:.2f}")

    res["by_threshold"] = {str(c): summarise(growth[growth["C_t"] >= c].reset_index(drop=True), 500, f"c{c}")
                           for c in [10, 20, 50]}
    for c, r in res["by_threshold"].items():
        print(f"  C>={c}: N={r['N']} a={r['a']['point']:.4f} b/a={r['b_over_a']['point']:.1f} b/bP={r['b_over_bP']['point']:.2f}")

    out = EXP / "analysis/china_covid19_empirical.json"
    out.write_text(json.dumps(res, indent=2), encoding="utf-8")
    print("written:", out)


if __name__ == "__main__":
    main()
