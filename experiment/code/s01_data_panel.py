"""s01: verify the NHSN snapshot, build the state-week panel, fix the early-growth windows.

Outputs: data/panel_nhsn_vintage_2026-10-20.csv.gz, data/national_flu.csv, analysis/parts/s01_data.json.
Continuity check: the S=3 calendar-window fit must reproduce a = 0.099055 (vintage) / 0.099060 (original).
"""
import numpy as np
import pandas as pd

from cc_common import (DATA, EARLY_CAL, PANEL, RAW_CSV, SEASONS, Log, Part, build_panel, fit_stat,
                       load_raw, mask_windows, national_series, rule_windows, sha256)

log = Log("s01_data_panel")
P = Part("s01_data", "s01_data_panel.py", [str(RAW_CSV.name)], "measured")

raw_all = pd.read_csv(RAW_CSV)
P.num("raw_rows", int(len(raw_all)))
P.num("raw_columns", int(raw_all.shape[1]))
P.num("raw_sha256_gz", sha256(RAW_CSV))
wk = pd.to_datetime(raw_all["Week Ending Date"])
P.num("raw_week_min", wk.min().strftime("%Y-%m-%d"))
P.num("raw_week_max", wk.max().strftime("%Y-%m-%d"))
log("raw", raw_all.shape, wk.min().date(), wk.max().date())
del raw_all

raw = load_raw("flu")
P.num("n_jurisdictions", int(raw["Geographic aggregation"].nunique()), note="50 states + DC; PR, VI, GU, MP, AS, US excluded")
panel = build_panel(raw)
panel.to_csv(PANEL, index=False, compression="gzip")
P.num("panel_rows_Dge5", int(len(panel)))
log("panel rows", len(panel))

nat = national_series(raw)
nat.to_csv(DATA / "national_flu.csv", index=False)

# regime counts (in-season, D>=5) — m8/D27 one definition
ins = panel[panel.in_season]
for k in (1, 2, 3):
    P.num(f"regime{k}_N_inseason", int((ins.period == k).sum()))

# continuity spec: S=3 calendar windows
m3 = mask_windows(panel, EARLY_CAL)
e3 = panel[m3]
f3 = fit_stat(e3.G, e3.D_t)
P.num("cal3_N", int(len(e3)))
P.num("cal3_imputed_zero_rows", int(e3.any_missing.sum()), note="m12/D31: early window contains no imputed zeros")
P.num("cal3_a_G", f3["a"]); P.num("cal3_b_G", f3["b"]); P.num("cal3_ba_G", f3["ba"]); P.num("cal3_Gbar", f3["xbar"])
log("calendar S=3", len(e3), f3)

# primary spec: rule-based windows, S=4
win, peaks = rule_windows(nat)
P.put("rule_windows", {s: {"start": win[s][0], "end": win[s][1], "national_peak": peaks[s]} for s in SEASONS})
log("rule windows", {s: (str(win[s][0].date()), str(win[s][1].date()), str(peaks[s].date())) for s in SEASONS})
m4 = mask_windows(panel, win.values())
e4 = panel[m4].copy()
P.num("rule4_N", int(len(e4)))
P.num("rule4_imputed_zero_rows", int(e4.any_missing.sum()))
P.num("rule4_n_weeks", int(e4.week_end.nunique()))
P.num("rule4_n_states", int(e4.state.nunique()))
per = e4.groupby("season").agg(N=("G", "size"), weeks=("week_end", "nunique"), Dmed=("D_t", "median"))
P.put("rule4_by_season", per.reset_index().to_dict("records"))
log(per)
for stat in ["G", "RA_mid", "G1"]:
    sub = e4[np.isfinite(e4[stat])] if stat == "G" else e4[e4.C_t >= 5]   # single-week scale C_t >= 5
    f = fit_stat(sub[stat], sub.D_t if stat == "G" else sub.C_t)
    P.num(f"rule4_N_{stat}", int(len(sub)))
    P.num(f"rule4_a_{stat}", f["a"]); P.num(f"rule4_b_{stat}", f["b"]); P.num(f"rule4_ba_{stat}", f["ba"])
    P.num(f"rule4_mean_{stat}", f["xbar"])
    log("rule S=4", stat, len(sub), f)
# national growth inside rule windows (sanity)
gw = [nat[(nat.week_end >= win[s][0]) & (nat.week_end <= win[s][1])].g.min() for s in SEASONS]
P.num("rule4_min_national_growth", float(min(gw)))
log.close()
