"""Fig 6 (2x2): between-season uncertainty, cross-pathogen weight ratios and surveillance design.
Reads analysis/parts/s05_seasons.json and s06_pathogens.json only. No recomputation.

A: season-specific intercepts a_s (G statistic) with two-way bootstrap 95% intervals, calendar windows
   (S = 3, continuity) and rule windows (S = 4, primary windows), plus the DerSimonian-Laird a_RE with its
   t interval (blocks seasons|G_cal3, seasons|G_rule4).
B: relative standard error of the between-season variance, sqrt(2/(S-1)) (Corollary 2, corrected), and the
   sqrt(2/S) values printed in the submitted manuscript (block relse_bound); S = 3 and 4 marked.
C: q1 = a_single/a_original and q2 = a_equal/a_original per pathogen x window with two-way 95% intervals;
   grey bars = growth-null ranges (pure reporting, A_all band); filled markers = the one-sided test against 1
   survives Holm correction over the 40-test family (blocks combo|*, multiplicity).
D: periods T needed to reach the precision of T_inf noise-free periods, T(m) = T_inf (1 + m_x/m), with
   v_R and A from the primary R_A fit (block design); dotted line at m_x.
Usage: python3 make_fig6.py --lang en|zh
"""
import json

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import fig_style as fs
from cc_common import PARTS

lang = fs.cli()
fs.setup(lang)
T = lambda en, zh: fs.T(lang, en, zh)  # noqa: E731
s5 = json.loads((PARTS / "s05_seasons.json").read_text())["blocks"]
s6 = json.loads((PARTS / "s06_pathogens.json").read_text())["blocks"]

fig, axs = plt.subplots(2, 2, figsize=(fs.FULL_W, fs.FULL_W * 0.66))

# ---------------------------------------------------------------- A
ax = axs[0, 0]
seasons = ["2022-23", "2023-24", "2024-25", "2025-26"]
xpos = {s: i for i, s in enumerate(seasons + ["RE"])}
for j, (spec, lab, mk) in enumerate([("G_cal3", T("Calendar, S = 3", "日历窗口 $S=3$"), "o"),
                                     ("G_rule4", T("Rule, S = 4", "规则窗口 $S=4$"), "s")]):
    blk = s5[f"seasons|{spec}"]
    off = -0.12 + 0.24 * j
    col = fs.DARK2[j]
    for r in blk["rows"]:
        x = xpos[r["season"]] + off
        ax.errorbar(x, r["a"], yerr=[[r["a"] - r["ci95"][0]], [r["ci95"][1] - r["a"]]], fmt=mk, color=col,
                    ms=3.5, lw=0.8, capsize=1.5)
    re = blk["random_effects"]
    ax.errorbar(xpos["RE"] + off, re["a_RE"], yerr=[[re["a_RE"] - re["ci95_t"][0]], [re["ci95_t"][1] - re["a_RE"]]],
                fmt="D", color=col, mfc="white", ms=4, lw=0.8, capsize=1.5)
    ax.plot([], [], mk, color=col, label=lab, ms=3.5)
ax.axhline(0, color=fs.GREY, lw=0.6, ls="--")
ax.axvline(3.5, color=fs.GREY, lw=0.5)
ax.set_xticks(range(5))
ax.set_xticklabels(seasons + [T("Pooled", "合并")])
for t in ax.get_xticklabels():
    t.set_fontproperties(fs.fp_for(t.get_text()))
    t.set_fontsize(7)
fs.xlabel(ax, T("Season", "流行季"))
fs.ylabel(ax, T("Intercept, $a$", "截距 $a$"))
fs.legend(ax, T("Windows", "窗口"), loc="upper left")
fs.panel_letter(ax, "A", dx=-0.17)

# ---------------------------------------------------------------- B
ax = axs[0, 1]
rb = s5["relse_bound"]
S = np.array(rb["S"], float)
ax.plot(S, rb["sqrt_2_over_S_minus_1"], color=fs.DARK2[0], label=r"$\sqrt{2/(S-1)}$")
ax.plot(S, np.sqrt(2 / S), color=fs.DARK2[1], ls=(0, (4, 1.5)), label=r"$\sqrt{2/S}$")
for s, mk in ((3, "o"), (4, "s")):
    ax.plot(s, np.sqrt(2 / (s - 1)), mk, color=fs.DARK2[0], ms=4)
    ax.plot(s, np.sqrt(2 / s), mk, color=fs.DARK2[1], mfc="white", ms=4)
ax.set_xlim(1.5, 20.5)
ax.set_ylim(0, 1.5)
ax.set_xticks([2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
fs.xlabel(ax, T("Number of seasons, $S$", "流行季数 $S$"))
fs.ylabel(ax, T("Relative SE of season variance", "季间方差相对标准误"))
fs.legend(ax, T("Formula", "公式"), loc="upper right")
fs.panel_letter(ax, "B", dx=-0.17)

# ---------------------------------------------------------------- C
ax = axs[1, 0]
order = ["flu|W-rule4", "flu|W-cal", "flu|W-growth", "flu|W-growth-2025",
         "covid|W-cal", "covid|W-growth", "rsv|W-cal", "rsv|W-growth"]
names = {"flu": T("Influenza", "流感"), "covid": T("COVID-19", "新冠"), "rsv": T("RSV", "合胞病毒")}
wn = {"W-rule4": T("rule", "规则"), "W-cal": T("calendar", "日历"), "W-growth": T("growth", "增长期"),
      "W-growth-2025": T("growth 2025", "增长期·最新季")}
tests = {(t["combo"], t["criterion"]): t for t in s6["multiplicity"]["tests"]}
yl = []
for i, c in enumerate(order):
    b = s6[f"combo|{c}"]
    y = len(order) - 1 - i
    nb = b["null_band"]
    for k, dy in (("q1", 0.15), ("q2", -0.15)):
        ax.add_patch(plt.Rectangle((nb[k][0], y + dy - 0.11), nb[k][1] - nb[k][0], 0.22, color="#cccccc", lw=0))
    for k, dy, col, crit in (("q1", 0.15, fs.DARK2[0], "q1_lt_1"), ("q2", -0.15, fs.DARK2[1], "q2_gt_1")):
        lo, hi = b["ci95"][k]
        lo = max(lo, 0.0)
        sig = tests[(c, crit)]["sig_holm_0.05"]
        ax.plot([lo, hi], [y + dy] * 2, color=col, lw=0.9)
        ax.plot(b["point"][k], y + dy, "o", color=col, mfc=col if sig else "white", ms=3.5, mew=0.9)
    pg, w = c.split("|")
    yl.append(f"{names[pg]} {wn[w]}" if lang == "en" else f"{names[pg]}·{wn[w]}")
ax.plot([1, 1], [-0.6, len(order) - 0.4], color="black", lw=0.6, ls="--")
ax.set_yticks(range(len(order)))
ax.set_yticklabels(yl[::-1])
for t in ax.get_yticklabels():
    t.set_fontproperties(fs.fp_for(t.get_text()))
    t.set_fontsize(7)
ax.set_xlim(0, 3.2)
ax.set_ylim(-0.6, len(order) + 2.9)
ax.grid(axis="y", visible=False)
fs.xlabel(ax, T("Intercept ratio", "截距比"))
h = [Line2D([], [], color=fs.DARK2[0], marker="o", ms=3.5, label=r"$q_1$"),
     Line2D([], [], color=fs.DARK2[1], marker="o", ms=3.5, label=r"$q_2$"),
     Patch(color="#cccccc", label=T("Reporting-only range", "纯报告零假设范围")),
     Line2D([], [], color="black", marker="o", mfc="white", ls="", ms=3.5, label=T("Not significant (Holm)", "校正后不显著"))]
fs.legend(ax, T("Ratio", "比值"), handles=h, loc="upper right", ncol=2, columnspacing=1.0, handlelength=1.5)
fs.panel_letter(ax, "C", dx=-0.30)

# ---------------------------------------------------------------- D
ax = axs[1, 1]
dz = s5["design"]
m = np.array(dz["m"])
mx = dz["A"] / dz["vR"]
for i, Tinf in enumerate(dz["targets_T_at_infinite_m"]):
    ax.plot(m, Tinf * (1 + mx / m), color=fs.DARK2[i], ls=["-", (0, (4, 1.5)), (0, (1, 1.2)), (0, (5, 1.5, 1, 1.5))][i],
            label=f"{Tinf}")
ax.axvline(mx, color=fs.GREY, lw=0.6, ls=":")
fs.text(ax, mx * 1.15, 1.6e3, r"$m_\times$", color=fs.GREY)
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(1, 1e4); ax.set_ylim(2, 4e3)
fs.xlabel(ax, T("Reported admissions per period, $m$", "单期报告入院数 $m$"))
fs.ylabel(ax, T("Periods needed, $T$", "所需期数 $T$"))
fs.legend(ax, T("Noise-free periods", "无噪声等效期数"), loc="upper right")
fs.panel_letter(ax, "D", dx=-0.17)

fig.tight_layout(h_pad=1.0, w_pad=1.5)
fs.save(fig, "fig6_seasons_pathogens_design", lang)
