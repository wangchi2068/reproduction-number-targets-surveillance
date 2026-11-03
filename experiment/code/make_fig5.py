"""Fig 5 (2x2): reporting completeness. Reads analysis/parts/s04_reporting.json only; no recomputation.

A: distribution of the reporting-completeness proxy p_hat (% hospitals reporting influenza admissions)
   by reporting regime (blocks.phat_hist: edges, by_regime[1|2|3].counts / N); density per
   percentage point, step outline with light fill.
B: Delta a = a(statistic) - a(coverage-adjusted statistic), pooled and per season, 90% paired two-way
   bootstrap CI (blocks["delta|G_rule4"|"delta|RA_rule4"].by[*].delta_a / ci90), each divided by its own
   TOST margin at tau0 = 0.1 (blocks[...].margins.tau0_0.1) so both statistics share one axis.
   Bands: |Delta a / margin| <= 1 (tau0 = 0.1) and <= 0.25 (tau0 = 0.05; margins scale with tau0^2).
C: TOST power (normal approximation) against true Delta a, expressed in units of the tau0 = 0.1
   margin (blocks[...].power_curve.tau0_0.1 / tau0_0.05: true_delta, power); solid tau0 = 0.1,
   dashed tau0 = 0.05; dotted reference at power 0.8. power_vs_S is reported in the notes only.
D: calendar-matched May-October comparison (blocks.calendar_matched[2023|2024|2025]): IQR of log G
   with two-way bootstrap 95% CI; 2024 = voluntary reporting (regime 2); labels give median p_hat.
Usage: python3 make_fig5.py --lang en|zh
"""
import json

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import fig_style as fs
from cc_common import PARTS

lang = fs.cli()
fs.setup(lang)
T = lambda en, zh: fs.T(lang, en, zh)  # noqa: E731
d = json.loads((PARTS / "s04_reporting.json").read_text())["blocks"]
REG = {"1": (fs.DARK2[0], T("Mandatory, to Apr 2024", "强制报告（前期）")),
       "2": (fs.DARK2[1], T("Voluntary, May–Oct 2024", "自愿报告")),
       "3": (fs.DARK2[2], T("Mandatory, from Nov 2024", "强制报告（后期）"))}
COL = {"RA_rule4": fs.DARK2[0], "G_rule4": fs.DARK2[1]}
MK = {"RA_rule4": "o", "G_rule4": "s"}
LAB = {"RA_rule4": "$R_A$", "G_rule4": "$G$"}


def ticklabels_font(labels):
    for t in labels:
        if fs.CJK_RE.search(t.get_text()):
            t.set_fontproperties(FontProperties(family=fs.CJK, size=7.5))


fig, axs = plt.subplots(2, 2, figsize=(fs.FULL_W, fs.FULL_W * 0.70))

# ---------------------------------------------------------------- A
ax = axs[0, 0]
h = d["phat_hist"]
edges = np.array(h["edges"])
bw = np.diff(edges)
for r in ["1", "2", "3"]:
    b = h["by_regime"][r]
    dens = np.array(b["counts"]) / b["N"] / bw
    col, lab = REG[r]
    ax.stairs(dens, edges, color=col, lw=1.0, label=lab)
    ax.stairs(dens, edges, color=col, alpha=0.15, fill=True)
ax.set_xlim(0, 100); ax.set_ylim(0, 0.105)
fs.xlabel(ax, T("Hospitals reporting, $\\hat{p}$ (%)", "报告医院比例 $\\hat{p}$（%）"))
fs.ylabel(ax, T("Density (per percentage point)", "密度（每百分点）"))
fs.legend(ax, T("Reporting regime", "报告制度"), loc="upper left")
fs.panel_letter(ax, "A", dx=-0.17)

# ---------------------------------------------------------------- B
ax = axs[0, 1]
rows = ["pooled", "2022-23", "2023-24", "2024-25", "2025-26"]
yy = np.arange(len(rows))[::-1]
ax.axvspan(-1, 1, color="#BDBDBD", alpha=0.25, lw=0)
ax.axvspan(-0.25, 0.25, color="#BDBDBD", alpha=0.45, lw=0)
ax.axvline(0, color="black", lw=0.5)
for k, spec in enumerate(["RA_rule4", "G_rule4"]):
    blk = d[f"delta|{spec}"]
    m = blk["margins"]["tau0_0.1"]
    off = 0.15 if k == 0 else -0.15
    for r, y in zip(rows, yy):
        v = blk["by"][r]
        lo, hi = np.array(v["ci90"]) / m
        ax.plot([lo, hi], [y + off] * 2, color=COL[spec], lw=0.9)
        ax.plot([v["delta_a"] / m], [y + off], MK[spec], color=COL[spec], ms=3.2)
ax.set_xlim(-2.2, 2.2)
ax.set_ylim(-0.6, len(rows) - 0.4)
ax.set_yticks(yy)
ax.set_yticklabels([T("Pooled", "合并") if r == "pooled" else r for r in rows])
ticklabels_font(ax.get_yticklabels())
ax.grid(axis="y", visible=False)
fs.ylabel(ax, T("Season", "流行季"))
fs.xlabel(ax, T("$\\Delta a$ / margin at $\\tau_0=0.1$", "$\\Delta a$ 与 $\\tau_0=0.1$ 界值之比"))
hh = [Line2D([], [], color=COL[s], marker=MK[s], ms=3.2, label=LAB[s]) for s in ["RA_rule4", "G_rule4"]]
hh += [Patch(facecolor="#BDBDBD", alpha=0.25, lw=0, label="$\\tau_0=0.1$"),
       Patch(facecolor="#BDBDBD", alpha=0.6, lw=0, label="$\\tau_0=0.05$")]
fs.legend(ax, T("Statistic; margin", "统计量与界值"), handles=hh, loc="upper right", ncol=1)
fs.panel_letter(ax, "B", dx=-0.17)

# ---------------------------------------------------------------- C
ax = axs[1, 0]
LS = {"tau0_0.1": "-", "tau0_0.05": (0, (4, 1.5))}
for spec in ["RA_rule4", "G_rule4"]:
    blk = d[f"delta|{spec}"]
    m = blk["margins"]["tau0_0.1"]
    for mk in ["tau0_0.1", "tau0_0.05"]:
        pc = blk["power_curve"][mk]
        ax.plot(np.array(pc["true_delta"]) / m, pc["power"], color=COL[spec], ls=LS[mk], lw=1.0)
ax.axhline(0.8, color=fs.GREY, lw=0.6, ls=":")
ax.set_xlim(-1.5, 1.5); ax.set_ylim(-0.03, 1.5)
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
fs.xlabel(ax, T("True $\\Delta a$ / margin at $\\tau_0=0.1$", "真实 $\\Delta a$ 与 $\\tau_0=0.1$ 界值之比"))
fs.ylabel(ax, T("Power of equivalence test", "等效检验功效"))
hh = [Line2D([], [], color=COL[s], ls="-", label=LAB[s]) for s in ["RA_rule4", "G_rule4"]]
hh += [Line2D([], [], color="black", ls="-", label="$\\tau_0=0.1$"),
       Line2D([], [], color="black", ls=LS["tau0_0.05"], label="$\\tau_0=0.05$")]
fs.legend(ax, T("Statistic; margin", "统计量与界值"), handles=hh, loc="upper center", ncol=2,
          columnspacing=1.0)
fs.panel_letter(ax, "C", dx=-0.17)

# ---------------------------------------------------------------- D
ax = axs[1, 1]
cm = d["calendar_matched"]
years = ["2023", "2024", "2025"]
for i, y in enumerate(years):
    v = cm[y]
    col = REG[str(v["regime"])][0]
    lo, hi = v["iqr_logG_ci95"]
    ax.plot([i, i], [lo, hi], color=col, lw=1.0)
    ax.plot([i], [v["iqr_logG"]], "o" if v["regime"] != 2 else "D", color=col, ms=4.5)
    ax.text(i + 0.08, v["iqr_logG"], f"$\\hat{{p}}={v['phat_median']:.0f}\\%$", fontsize=7,
            va="center", ha="left", color="black")
ax.set_xlim(-0.5, 2.7); ax.set_ylim(0.3, 1.35)
ax.set_xticks(range(3)); ax.set_xticklabels(years)
ax.grid(axis="x", visible=False)
fs.xlabel(ax, T("May to October of year", "年份（五月至十月）"))
fs.ylabel(ax, T("IQR of $\\log G$", "$\\log G$ 的四分位距"))
hh = [Line2D([], [], color=REG["1"][0], marker="o", ms=4, label=T("Mandatory (earlier)", "强制报告（前期）")),
      Line2D([], [], color=REG["2"][0], marker="D", ms=4, label=T("Voluntary", "自愿报告")),
      Line2D([], [], color=REG["3"][0], marker="o", ms=4, label=T("Mandatory (later)", "强制报告（后期）"))]
fs.legend(ax, T("Reporting regime", "报告制度"), handles=hh, loc="upper right")
fs.panel_letter(ax, "D", dx=-0.17)

fig.tight_layout(h_pad=1.0, w_pad=1.5)
fs.save(fig, "fig5_reporting", lang)
