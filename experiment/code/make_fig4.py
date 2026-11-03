"""Fig 4 (2x2): empirical floor of the relative loss. Reads only stored results; no recomputation.

A: national weekly influenza admissions (data/national_flu.csv, 2022-08 .. 2026-06, log y); shaded
   rule-based early-growth windows (s01_data.json blocks.rule_windows: start/end, primary, S = 4),
   hatched calendar windows of the submitted manuscript (cc_common.EARLY_CAL, S = 3), and the first
   national turning point per season (rule_windows[*].national_peak).
B: binned relative loss L against scale D (s03_decomposition.json blocks["plateau|RA_rule4"],
   blocks["plateau|G_rule4"]: bins D_mean/L_mean/L_sd/n; affine [a, b]; power [c, gamma]).
   Points = bin means, bars = +-1.96 L_sd/sqrt(n) (iid within bin, ignores clustering);
   solid = a + b/D, dashed = c D^-gamma, dotted = floor a.
C: intercept a after successive fixed-effect adjustments (blocks["decomp|<spec>"]: a_total,
   a_weekFE, a_both with two-way bootstrap 95% CI, B = 1000) for RA_rule4, G_rule4, G_cal3.
D: forest of b/a (blocks["forest"], two-way bootstrap 95% CI, B = 500) for RA_rule4 and G_rule4,
   log x; grey band = theory b/a range for G with the b-delta correction
   (numbers["G_rule4_theory_ba_delta_range"]).
Usage: python3 make_fig4.py --lang en|zh
"""
import json
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.font_manager import FontProperties

import fig_style as fs
from cc_common import DATA, EARLY_CAL, PARTS

lang = fs.cli()
fs.setup(lang)
T = lambda en, zh: fs.T(lang, en, zh)  # noqa: E731
s01 = json.loads((PARTS / "s01_data.json").read_text())["blocks"]
P3 = json.loads((PARTS / "s03_decomposition.json").read_text())
s03, n03 = P3["blocks"], P3["numbers"]
COL = {"RA_rule4": fs.DARK2[0], "G_rule4": fs.DARK2[1], "G_cal3": fs.DARK2[2]}
LAB = {"RA_rule4": "$R_A$", "G_rule4": "$G$"}


def ticklabels_font(labels):
    for t in labels:
        if fs.CJK_RE.search(t.get_text()):
            t.set_fontproperties(FontProperties(family=fs.CJK, size=7.5))


fig, axs = plt.subplots(2, 2, figsize=(fs.FULL_W, fs.FULL_W * 0.86),
                        gridspec_kw={"height_ratios": [1, 1.3]})

# ---------------------------------------------------------------- A
ax = axs[0, 0]
nat = pd.read_csv(DATA / "national_flu.csv", parse_dates=["week_end"])
nat = nat[(nat.week_end >= "2022-08-01") & (nat.week_end <= "2026-06-30")]
ax.plot(nat.week_end, nat.C, color="black", lw=0.8)
half = pd.Timedelta(days=3.5)
for s, w in s01["rule_windows"].items():
    ax.axvspan(pd.Timestamp(w["start"]) - half, pd.Timestamp(w["end"]) + half,
               color=fs.DARK2[0], alpha=0.30, lw=0)
for a, b in EARLY_CAL:
    ax.axvspan(pd.Timestamp(a) - half, pd.Timestamp(b) + half, fill=False, hatch="////",
               edgecolor=fs.DARK2[1], lw=0.6)
pk = [pd.Timestamp(w["national_peak"]) for w in s01["rule_windows"].values()]
pkC = [float(nat.loc[nat.week_end == p, "C"].iloc[0]) for p in pk]
ax.plot(pk, pkC, "v", color=fs.DARK2[2], ms=4, mec="black", mew=0.4, zorder=5)
ax.set_yscale("log")
ax.set_ylim(100, 4e6)
ax.set_xlim(pd.Timestamp("2022-08-01"), pd.Timestamp("2026-06-30"))
ax.set_xticks([pd.Timestamp(f"{y}-01-01") for y in range(2023, 2027)])
ax.set_xticklabels([str(y) for y in range(2023, 2027)])
fs.xlabel(ax, T("Week ending", "周末日期"))
fs.ylabel(ax, T("National admissions per week", "全国每周入院数"))
h = [Patch(facecolor=fs.DARK2[0], alpha=0.30, lw=0, label=T("Rule (primary)", "规则（主分析）")),
     Patch(facecolor="white", edgecolor=fs.DARK2[1], hatch="////", lw=0.6,
           label=T("Calendar (submitted)", "日历（原稿）")),
     Line2D([], [], marker="v", ls="none", color=fs.DARK2[2], mec="black", mew=0.4, ms=4,
            label=T("First turning point", "首个拐点"))]
fs.legend(ax, T("Early-growth window", "早期增长窗口"), handles=h, loc="upper left", ncol=2,
          columnspacing=0.8, handlelength=1.4, borderaxespad=0.4)
fs.panel_letter(ax, "A", dx=-0.17)

# ---------------------------------------------------------------- B
ax = axs[0, 1]
Dg = np.geomspace(5, 4000, 200)
for k, spec in enumerate(["RA_rule4", "G_rule4"]):
    p = s03[f"plateau|{spec}"]
    bins = pd.DataFrame(p["bins"])
    se = bins.L_sd / np.sqrt(bins.n)
    lo = np.maximum(bins.L_mean - 1.96 * se, 1e-4)
    ax.errorbar(bins.D_mean, bins.L_mean, yerr=[bins.L_mean - lo, 1.96 * se], fmt="o", ms=3,
                color=COL[spec], elinewidth=0.7, capsize=1.5, zorder=4)
    a, b = p["affine"]
    c, g = p["power"]
    ax.plot(Dg, a + b / Dg, color=COL[spec], ls="-", lw=0.9)
    ax.plot(Dg, c * Dg ** (-g), color=COL[spec], ls="--", lw=0.9)
    ax.plot([Dg[0], Dg[-1]], [a, a], color=COL[spec], ls=":", lw=0.9)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(4, 5000); ax.set_ylim(0.004, 1.0)
fs.xlabel(ax, T("Scale, $D$ (cases per week)", "规模 $D$（例/周）"))
fs.ylabel(ax, T("Mean relative loss, $L$", "平均相对损失 $L$"))
h = [Line2D([], [], color=COL[s], marker="o", ms=3, ls="none", label=LAB[s]) for s in ["RA_rule4", "G_rule4"]]
h += [Line2D([], [], color="black", ls="-", label=T("$a+b/D$", "仿射 $a+b/D$")),
      Line2D([], [], color="black", ls=(0, (4, 1.5)), label=T("Power law", "幂律")),
      Line2D([], [], color="black", ls=(0, (1, 1.2)), label=T("Intercept $a$", "截距 $a$"))]
fs.legend(ax, T("Statistic; fit", "统计量与拟合"), handles=h, loc="upper right", ncol=2,
          columnspacing=0.8, handlelength=1.6)
fs.panel_letter(ax, "B", dx=-0.17)

# ---------------------------------------------------------------- C
ax = axs[1, 0]
specs = ["RA_rule4", "G_rule4", "G_cal3"]
steps = [("a_total", T("None", "无调整")), ("a_weekFE", T("Week FE", "周固定效应")),
         ("a_both", T("Week + state-season FE", "周与州季固定效应"))]
wbar = 0.26
for j, (key, lab) in enumerate(steps):
    x = np.arange(len(specs)) + (j - 1) * wbar
    v = np.array([s03[f"decomp|{s}"][key]["value"] for s in specs])
    ci = np.array([s03[f"decomp|{s}"][key]["ci95"] for s in specs])
    ax.bar(x, v, width=wbar * 0.92, color=fs.DARK2[3 + j], alpha=0.85, label=lab, lw=0)
    ax.errorbar(x, v, yerr=[v - ci[:, 0], ci[:, 1] - v], fmt="none", ecolor="black",
                elinewidth=0.7, capsize=1.5)
ax.set_xticks(np.arange(len(specs)))
ax.set_xticklabels(["$R_A$\nrule windows", "$G$\nrule windows", "$G$\ncalendar windows"] if lang == "en"
                   else ["$R_A$\n规则窗口", "$G$\n规则窗口", "$G$\n日历窗口"])
ticklabels_font(ax.get_xticklabels())
for t in ax.get_xticklabels():
    t.set_linespacing(1.5)
ax.set_ylim(0, 0.19)
ax.grid(axis="x", visible=False)
fs.ylabel(ax, T("Intercept $a$", "截距 $a$"))
fs.legend(ax, T("Fixed effects removed", "剔除的固定效应"), loc="upper left")
fs.panel_letter(ax, "C", dx=-0.17)

# ---------------------------------------------------------------- D
ax = axs[1, 1]
rows = [("main", "Main", "主分析"), ("trim_top1pct", "Trim top 1%", "剔除最高$1\\%$"),
        ("drop_2022-23", "Drop 2022-23", "剔除$2022$–$23$季"),
        ("drop_2023-24", "Drop 2023-24", "剔除$2023$–$24$季"),
        ("drop_2024-25", "Drop 2024-25", "剔除$2024$–$25$季"),
        ("drop_2025-26", "Drop 2025-26", "剔除$2025$–$26$季"),
        ("huber", "Huber loss", "稳健回归"), ("weekFE", "Week FE", "周固定效应"),
        ("peak_8w", "Pre-peak 8 weeks", "峰前$8$周"), ("GI_short", "Short GI", "代际间隔偏短"),
        ("GI_long", "Long GI", "代际间隔偏长"), ("S3_calendar", "Calendar, $S=3$", "日历窗口 $S=3$"),
        ("D_weights_single_1w", "$D_t$ one week", "$D_t$ 单周"),
        ("D_weights_equal_3w", "$D_t$ equal 3 weeks", "$D_t$ 三周等权")]
fr = {(r["base"], r["spec"]): r for r in s03["forest"]}
th = n03["G_rule4_theory_ba_delta_range"]["value"]
ax.axvspan(th[0], th[1], color=fs.DARK2[1], alpha=0.15, lw=0)
yy = np.arange(len(rows))[::-1]
for k, spec in enumerate(["RA_rule4", "G_rule4"]):
    off = 0.18 if k == 0 else -0.18
    for (key, _, _), y in zip(rows, yy):
        r = fr.get((spec, key))
        if r is None:
            continue
        ax.plot(r["ci95"], [y + off] * 2, color=COL[spec], lw=0.8)
        ax.plot([r["ba"]], [y + off], "o" if k == 0 else "s", color=COL[spec], ms=3)
ax.set_xscale("log")
ax.set_xlim(2, 1000)
ax.set_yticks(yy)
ax.set_yticklabels([T(e, z) for _, e, z in rows])
ticklabels_font(ax.get_yticklabels())
ax.tick_params(axis="y", labelsize=7)
ax.set_ylim(-0.7, len(rows) - 0.3)
ax.grid(axis="y", visible=False)
fs.xlabel(ax, T("Balance scale, $b/a$ (cases per week)", "平衡规模 $b/a$（例/周）"))
h = [Line2D([], [], color=COL["RA_rule4"], marker="o", ms=3, label=LAB["RA_rule4"]),
     Line2D([], [], color=COL["G_rule4"], marker="s", ms=3, label=LAB["G_rule4"]),
     Patch(facecolor=fs.DARK2[1], alpha=0.15, lw=0, label=T("Theory, $G$", "理论范围 $G$"))]
fs.legend(ax, T("Statistic", "统计量"), handles=h, loc="lower right")
fs.panel_letter(ax, "D", dx=-0.42)

fig.tight_layout(h_pad=1.0, w_pad=1.2)
fs.save(fig, "fig4_empirical_floor", lang)
