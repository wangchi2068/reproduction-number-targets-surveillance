"""Fig 6 (round 3, 2x2): surveillance design nomogram.
Reads analysis/parts/{s02_inference,inversion,returns,multiscale,parag}.json (written by s02/s09). No recomputation
beyond the closed-form curves below.

A: nomogram for R_A (rule-4 windows): weeks T needed so that the relative standard error of a T-week mean of R_A
   in one state equals 2, 3, 5 or 8%: T = (a + b/m)/(xbar^2 relSE^2), m = weekly cases; assumes weekly
   deviations independent within a state (s02 a, b; xbar = sample mean, inversion.empirical). Vertical line =
   m_x = b/a with its two-way bootstrap interval (shaded). Zone labels left/right of m_x.
B: precision fraction m/(m + m_x) against m/m_x (returns.precision_fraction_at_*_mcross), with the variance cut
   bought by one doubling of m (returns.variance_cut_doubling_from_*).
C: J-state pooling of the weekly deviation for R_A and G (multiscale blocks Jpool|*): mean over weeks of the
   squared J-state mean deviation, relative to one state, 300 random subsets per J (band = 2.5-97.5%), and the
   one-way random-effects curve s2_week + (V1 - s2_week)/J with its floor.
D: relative variance of R_A per state-week: Poisson-only renewal variance b_P/m (Parag et al. 2022, eq. 3/5;
   reporting only rescales m) versus the fitted a + b/m (parag.curves); vertical line = m_x.
Usage: python3 make_fig6_nomogram.py --lang en|zh
"""
import json

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter

import fig_style as fs
from cc_common import PARTS

lang = fs.cli()
fs.setup(lang)
T = lambda en, zh: fs.T(lang, en, zh)  # noqa: E731
ZH_LS = {"linespacing": 1.6} if lang == "zh" else {}  # CJK glyph boxes overlap at default spacing
J_ = lambda n: json.loads((PARTS / f"{n}.json").read_text())  # noqa: E731
s02, inv, ret, ms, par = (J_(n) for n in ("s02_inference", "inversion", "returns", "multiscale", "parag"))
N2 = s02["numbers"]
a, b = N2["RA_rule4_a"]["value"], N2["RA_rule4_b"]["value"]
mx, mx_ci = N2["RA_rule4_ba"]["value"], N2["RA_rule4_ba"]["ci95"]
xbar = inv["blocks"]["empirical"]["RA"]["xbar"]
COL = {"RA": fs.DARK2[0], "G": fs.DARK2[1]}
LAB = {"RA": "$R_A$", "G": "$G$"}

fig, axs = plt.subplots(2, 2, figsize=(fs.FULL_W, fs.FULL_W * 0.70))

# ---------------------------------------------------------------- A nomogram
ax = axs[0, 0]
m = np.geomspace(2, 2000, 300)
Tw = np.geomspace(1, 200, 300)
MM, TT = np.meshgrid(m, Tw)
rel = np.sqrt((a + b / MM) / TT) / xbar
lev = [0.02, 0.03, 0.05, 0.08]
cs = ax.contour(MM, TT, rel, levels=lev, colors="black", linewidths=0.8)
mlab = 300.0
ax.clabel(cs, fmt={l: f"{100 * l:g}%" for l in lev}, fontsize=7, inline_spacing=2,
          manual=[(mlab, (a + b / mlab) / (xbar * l) ** 2) for l in lev])
ax.axvspan(*mx_ci, color="#EEEEEE", lw=0, zorder=0)
ax.axvline(mx, color=COL["RA"], lw=1.0, ls="--")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(2, 2000); ax.set_ylim(1, 200)
ax.set_xticks([2, 10, 100, 1000]); ax.set_xticklabels(["2", "10", "100", "1000"])
ax.set_yticks([1, 4, 13, 52, 200]); ax.set_yticklabels(["1", "4", "13", "52", "200"])
fs.text(ax, mx * 1.08, 1.1, T(f"$m_\\times = {mx:.0f}$", f"平衡尺度 ${mx:.0f}$"), color=COL["RA"], ha="left",
        va="bottom")
fs.text(ax, 2.3, 1.1, T("Scale-up\nhigh-return zone", "扩样\n高收益区"), ha="left", va="bottom", **ZH_LS)
fs.text(ax, 1800, 180, T("More weeks / regions\nhigh-return zone", "延长周期与跨区协同\n高收益区"), ha="right", va="top", **ZH_LS)
fs.xlabel(ax, T("Weekly cases per state, $m$", "单州每周病例数"))
fs.ylabel(ax, T("Observation weeks, $T$", "观测周数"))
fs.panel_letter(ax, "A", dx=-0.17)

# ---------------------------------------------------------------- B precision fraction
ax = axs[0, 1]
r = np.geomspace(0.05, 20, 200)
ax.plot(r, r / (1 + r), color="black", lw=1.0, label=T("Precision share", "可达精度占比"))
ax.plot(r, 1 / (2 * (1 + r)), color=fs.GREY, lw=1.0, ls="--",
        label=T("Cut per doubling", "翻倍方差降幅"))
RN = ret["numbers"]
for lab, f_ in (("half", 0.5), ("1x", 1), ("2x", 2), ("4x", 4)):
    pf = RN[f"precision_fraction_at_{lab}_mcross"]["value"]
    vc = RN[f"variance_cut_doubling_from_{lab}_mcross"]["value"]
    ax.plot(f_, pf, "o", color="black", ms=3.5)
    ax.plot(f_, vc, "o", color=fs.GREY, ms=3.5)
    fs.text(ax, f_ / 1.12, pf + 0.03, f"${100 * pf:.0f}\\%$", ha="right", va="bottom")
    if lab != "half":  # both curves equal 1/3 at m = m_x/2; one label serves both
        fs.text(ax, f_ * 1.10, vc + 0.02, f"${100 * vc:.0f}\\%$", color=fs.GREY, ha="left", va="bottom")
ax.axvline(1, color=COL["RA"], lw=0.8)
ax.set_xscale("log"); ax.set_xlim(0.05, 20); ax.set_ylim(0, 1.05)
ax.set_xticks([0.1, 0.5, 1, 2, 4, 10]); ax.set_xticklabels(["0.1", "0.5", "1", "2", "4", "10"])
fs.xlabel(ax, T("Case count relative to balance scale, $m/m_\\times$", "病例数与平衡尺度之比"))
fs.ylabel(ax, T("Fraction", "比例"))
fs.legend(ax, T("Quantity", "指标"), loc="center right", bbox_to_anchor=(1.0, 0.42))
fs.panel_letter(ax, "B", dx=-0.17)

# ---------------------------------------------------------------- C J pooling
ax = axs[1, 0]
for nm, spec in (("RA", "RA_rule4"), ("G", "G_rule4")):
    blk = ms["blocks"][f"Jpool|{spec}"]
    Js = np.array(blk["J"]); V1 = blk["V1"]
    e = blk["empirical"]
    mean = np.array([e[str(j)]["mean"] for j in Js]) / V1
    lo = np.array([e[str(j)]["q"][0] for j in Js]) / V1
    hi = np.array([e[str(j)]["q"][1] for j in Js]) / V1
    ax.fill_between(Js, lo, hi, color=COL[nm], alpha=0.18, lw=0)
    ax.plot(Js, mean, "o", color=COL[nm], ms=3, label=LAB[nm])
    jj = np.geomspace(1, 51, 100)
    ax.plot(jj, (blk["s2_week"] + (V1 - blk["s2_week"]) / jj) / V1, color=COL[nm], lw=0.9)
    ax.axhline(blk["s2_week"] / V1, color=COL[nm], lw=0.7, ls=":")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(0.9, 60); ax.set_ylim(0.12, 1.3)
ax.set_xticks([1, 2, 5, 10, 20, 51]); ax.set_xticklabels(["1", "2", "5", "10", "20", "51"])
ax.set_yticks([0.2, 0.3, 0.5, 1.0]); ax.set_yticklabels(["0.2", "0.3", "0.5", "1"])
fs.xlabel(ax, T("States pooled, $J$", "合并州数"))
fs.ylabel(ax, T("Variance relative to one state", "相对单州的方差"))
fs.legend(ax, T("Statistic", "统计量"), loc="upper right")
fs.text(ax, 1.05, 0.13, T("dotted: national-week component", "点线：全国周同步分量"), color=fs.GREY, va="bottom")
fs.panel_letter(ax, "C", dx=-0.17)

# ---------------------------------------------------------------- D Parag contrast
ax = axs[1, 1]
cv = par["blocks"]["curves"]
cm = np.array(cv["m"])
ax.plot(cm, np.array(cv["poisson"]), color=fs.GREY, lw=1.0, ls="--",
        label=T("Poisson renewal only", "仅泊松更新"))
ax.plot(cm, np.array(cv["model"]), color=COL["RA"], lw=1.2, label=T("Fitted $a+b/m$", "拟合曲线"))
ax.axhline(a, color=COL["RA"], lw=0.7, ls=":")
ax.axvline(mx, color=COL["RA"], lw=0.8, ls="--")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(2, 5000); ax.set_ylim(1e-4, 0.5)
ax.set_xticks([2, 10, 100, 1000]); ax.set_xticklabels(["2", "10", "100", "1000"])
fs.xlabel(ax, T("Weekly cases per state, $m$", "单州每周病例数"))
fs.ylabel(ax, T("Relative variance of $R_A$", "更新比的相对方差"))
fs.legend(ax, T("Variance model", "方差模型"), loc="upper right")
fs.panel_letter(ax, "D", dx=-0.17)

for ax in axs.flat:
    ax.xaxis.set_minor_formatter(NullFormatter()); ax.yaxis.set_minor_formatter(NullFormatter())
fig.tight_layout(h_pad=1.0, w_pad=1.2)
fs.save(fig, "fig6_design_nomogram", lang)
