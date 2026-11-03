"""Fig 3 (round 3, 2x2): simulation-based inversion from apparent statistics to latent dynamics.
Reads analysis/parts/inversion.json, s09_vR_draws.npz and s09_ba_reps.npz (written by s09_inversion.py).
No recomputation.

A: calibration lines of the data-matched simulation (M0): mean fitted intercept E[a_hat] over 200 simulated
   panels, divided by the empirical intercept of the same statistic, against v_R; error bars = 95% Monte
   Carlo intervals of the cell mean; the horizontal line at 1 is the empirical intercept and the dotted drop
   lines mark where each calibration line reaches it (v_R^implied); grey band = v_R in [0.04, 0.08]
   (blocks cells, sensitivity.M0).
B: scale factor kappa = E[a_hat]/(v_R/Rbar^2) at the cell nearest v_R^implied for every simulation setting
   (filled) and the noise-free limit kappa_inf (open); log axis (block sensitivity).
C: distribution of v_R^implied under M0: joint two-way bootstrap of the empirical intercepts (B = 2000)
   combined with Monte Carlo draws of the calibration line; ticks on top = point values of all non-bridge
   sensitivity settings; grey band = [0.04, 0.08] (s09_vR_draws.npz, block sensitivity).
D: balance ratio b/a: per-panel b/a of the data-matched simulation over the 8 calibration cells
   (box = 25-75%, whiskers = 2.5-97.5%), the range of pooled E[b]/E[a] (thick bar), and the empirical
   estimate with its two-way bootstrap interval (s02_inference) (block ba_calibrated).
Usage: python3 make_fig3_inversion.py --lang en|zh
"""
import json

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter
from matplotlib.lines import Line2D

import fig_style as fs
from cc_common import PARTS

lang = fs.cli()
fs.setup(lang)
T = lambda en, zh: fs.T(lang, en, zh)  # noqa: E731
inv = json.loads((PARTS / "inversion.json").read_text())
blk = inv["blocks"]
cells, sens, emp = blk["cells"], blk["sensitivity"], blk["empirical"]
vd = np.load(PARTS / "s09_vR_draws.npz")
bar = np.load(PARTS / "s09_ba_reps.npz")
s02 = json.loads((PARTS / "s02_inference.json").read_text())["numbers"]
VG = blk["design"]["v_grid"]
COL = {"RA": fs.DARK2[0], "G": fs.DARK2[1]}
MK = {"RA": "o", "G": "s"}
LAB = {"RA": "$R_A$", "G": "$G$"}
BAND = (0.04, 0.08)

fig, axs = plt.subplots(2, 2, figsize=(fs.FULL_W, fs.FULL_W * 0.80))

# ---------------------------------------------------------------- A
ax = axs[0, 0]
ax.axvspan(*BAND, color="#EEEEEE", lw=0, zorder=0)
ax.axhline(1.0, color="black", lw=0.8)
xx = np.linspace(0.0, 0.1, 50)
for j, nm in enumerate(("RA", "G")):
    ae = emp[nm]["a"]
    v = np.array([cells[f"M0|v{g}"]["v_R"] for g in VG])
    y = np.array([cells[f"M0|v{g}"][nm]["a_mean"] for g in VG]) / ae
    se = np.array([cells[f"M0|v{g}"][nm]["a_se"] for g in VG]) / ae
    off = (-0.0012, 0.0012)[j]
    ax.errorbar(v + off, y, yerr=1.96 * se, fmt=MK[nm], ms=3.5, color=COL[nm], mfc="white" if nm == "RA" else COL[nm],
                elinewidth=0.6, label=LAB[nm])
    s = sens["M0"][nm]
    ax.plot(xx, (s["c0"] + s["c1"] * xx) / ae, color=COL[nm], lw=0.9)
    vi = s["v_implied"]
    ax.annotate("", xy=(vi, 0.05), xytext=(vi, 1.0),
                arrowprops=dict(arrowstyle="-|>", color=COL[nm], lw=0.8, ls=(0, (2, 1.5)), mutation_scale=7))
    fs.text(ax, vi + (-0.0015 if nm == "RA" else 0.0015), 0.12, f"${vi:.3f}$", color=COL[nm],
            ha="right" if nm == "RA" else "left", va="bottom")
ax.set_xlim(0, 0.1); ax.set_ylim(0, 2.6)
ax.set_xticks([0.02, 0.04, 0.06, 0.08, 0.10]); ax.set_yticks([0.5, 1.0, 1.5, 2.0, 2.5])
fs.text(ax, 0.06, 2.42, T("theory band", "理论区间"), color=fs.GREY, ha="center", va="top")
fs.xlabel(ax, T("Latent environmental variance, $v_R$", "潜伏环境方差"))
fs.ylabel(ax, T("Simulated / empirical intercept", "模拟截距与实证截距之比"))
fs.legend(ax, T("Statistic", "统计量"), loc="upper left")
fs.panel_letter(ax, "A", dx=-0.22)

# ---------------------------------------------------------------- B
ax = axs[0, 1]
order = ["M0", "rho0", "rho65", "ar05", "gi_short", "gi_long", "delay_short", "delay_long", "k01", "count_shock",
         "bridge_s08"]
ZH = {"M0": "数据匹配基准", "rho0": "无全国冲击", "rho65": "全国份额 $0.65$", "ar05": "周间自相关 $0.5$",
      "gi_short": "真实代际间隔偏短", "gi_long": "真实代际间隔偏长", "delay_short": "入院延迟 $3.5$ 天",
      "delay_long": "入院延迟 $8$ 天", "k01": "离散参数 $0.1$", "count_shock": "同步报告冲击",
      "bridge_s08": "原仿真设计"}
EN = {"M0": "Data-matched (M0)", "rho0": "No national shock", "rho65": "National share 0.65",
      "ar05": "AR(1), 0.5", "gi_short": "True GI short", "gi_long": "True GI long", "delay_short": "Delay 3.5 d",
      "delay_long": "Delay 8 d", "k01": "Dispersion 0.1", "count_shock": "Sync. reporting shock",
      "bridge_s08": "s08 design (bridge)"}
order = [o for o in order if o in sens]
yy = np.arange(len(order))[::-1]
for nm in ("RA", "G"):
    kf = [sens[o][nm]["kappa_at_nearest"] for o in order]
    ki = [sens[o][nm]["kappa_inf_at_nearest"] for o in order]
    ax.plot(kf, yy + (0.15 if nm == "RA" else -0.15), MK[nm], color=COL[nm], ms=3.5, ls="none")
    ax.plot(ki, yy + (0.15 if nm == "RA" else -0.15), MK[nm], color=COL[nm], mfc="white", ms=3.5, ls="none")
ax.axvline(1.0, color=fs.GREY, lw=0.5)
ax.set_xscale("log"); ax.set_xlim(0.25, 10)
ax.set_xticks([0.3, 0.5, 1, 2, 4, 8]); ax.set_xticklabels(["0.3", "0.5", "1", "2", "4", "8"])
ax.set_yticks(yy)
labs = [T(EN[o], ZH[o]) for o in order]
ax.set_yticklabels(labs)
for t in ax.get_yticklabels():
    if lang == "zh":
        t.set_fontproperties(fs.fp_for(t.get_text()))
    t.set_fontsize(7)
ax.set_ylim(-0.7, len(order) + 1.6)
fs.xlabel(ax, T("Scale factor, $\\kappa$ (fitted intercept per unit $v_R/\\bar{R}^2$)", "尺度因子（单位名义截距对应的拟合截距）"))
h = [Line2D([], [], color=COL["RA"], marker="o", ls="none", ms=3.5, label=LAB["RA"]),
     Line2D([], [], color=COL["G"], marker="s", ls="none", ms=3.5, label=LAB["G"]),
     Line2D([], [], color="black", marker="o", ls="none", ms=3.5, label=T("Fitted", "拟合")),
     Line2D([], [], color="black", marker="o", mfc="white", ls="none", ms=3.5, label=T("Noise-free limit", "无噪声极限"))]
fs.legend(ax, None, handles=h, loc="upper center", ncol=4, columnspacing=0.8, handletextpad=0.3)
fs.panel_letter(ax, "B", dx=-0.52)

# ---------------------------------------------------------------- C
ax = axs[1, 0]
ax.axvspan(*BAND, color="#EEEEEE", lw=0, zorder=0)
bins = np.linspace(0.0, 0.12, 49)
top = 0
for nm in ("RA", "G"):
    d = vd[f"M0__{nm}"]
    h_, e_ = np.histogram(d[np.isfinite(d)], bins, density=True)
    ax.stairs(h_, e_, color=COL[nm], lw=1.0, label=LAB[nm])
    top = max(top, h_.max())
for nm in ("RA", "G"):
    s = sens["M0"][nm]
    ax.axvline(s["v_implied"], color=COL[nm], lw=0.8, ls="--")
ymax = top * 1.32
for j, nm in enumerate(("RA", "G")):
    pts = [sens[o][nm]["v_implied"] for o in sens if o != "bridge_s08"]
    yb = ymax * (0.95 if nm == "RA" else 0.89)
    ax.plot(pts, [yb] * len(pts), "|", color=COL[nm], ms=6, mew=0.9)
ax.set_xlim(0, 0.12); ax.set_ylim(0, ymax)
ax.set_xticks([0.02, 0.04, 0.06, 0.08, 0.10])
fs.xlabel(ax, T("Implied latent variance, $v_R^{\\mathrm{impl}}$", "反演得到的潜伏环境方差"))
fs.ylabel(ax, T("Density", "密度"))
fs.legend(ax, T("Statistic", "统计量"), loc="center right")
fs.text(ax, 0.118, ymax * 0.80, T("ticks: each setting", "刻线：各设定点估计"), color=fs.GREY, ha="right", va="top")
fs.panel_letter(ax, "C", dx=-0.22)

# ---------------------------------------------------------------- D
ax = axs[1, 1]
bc = blk["ba_calibrated"]
for i, nm in enumerate(("RA", "G")):
    x0 = i * 1.0
    r = bar[nm]
    r = r[np.isfinite(r) & (r > 0)]
    q = np.percentile(r, [2.5, 25, 50, 75, 97.5])
    ax.add_patch(plt.Rectangle((x0 - 0.17, q[1]), 0.34, q[3] - q[1], fc="white", ec=COL[nm], lw=0.8))
    ax.plot([x0 - 0.17, x0 + 0.17], [q[2]] * 2, color=COL[nm], lw=0.9)
    ax.plot([x0, x0], [q[0], q[1]], color=COL[nm], lw=0.7)
    ax.plot([x0, x0], [q[3], q[4]], color=COL[nm], lw=0.7)
    lo, hi = bc[nm]["pooled_range"]
    ax.plot([x0 - 0.30, x0 - 0.30], [lo, hi], color=COL[nm], lw=3.0, solid_capstyle="butt")
    spec = "RA_rule4" if nm == "RA" else "G_rule4"
    e = s02[f"{spec}_ba"]
    ax.errorbar([x0 + 0.32], [e["value"]], yerr=[[e["value"] - e["ci95"][0]], [e["ci95"][1] - e["value"]]],
                fmt="D", ms=4, color="black", elinewidth=0.8, capsize=2)
ax.set_yscale("log"); ax.set_ylim(3, 500); ax.set_xlim(-0.6, 1.6)
ax.set_yticks([5, 10, 20, 50, 100, 200]); ax.set_yticklabels(["5", "10", "20", "50", "100", "200"])
ax.set_xticks([0, 1]); ax.set_xticklabels(["$R_A$", "$G$"])
fs.ylabel(ax, T("Balance ratio, $b/a$ (cases per week)", "平衡比（每周病例数）"))
h = [Line2D([], [], color=fs.GREY, lw=3.0, label=T("Calibrated theory range", "校准后理论范围")),
     Line2D([], [], color=fs.GREY, marker="s", mfc="white", ls="none", ms=5, label=T("Simulated panels", "模拟面板")),
     Line2D([], [], color="black", marker="D", ls="none", ms=4, label=T("Empirical, 95% CI", "实证估计及区间"))]
fs.legend(ax, T("Source", "来源"), handles=h, loc="upper right")
fs.panel_letter(ax, "D", dx=-0.17)

for ax in axs.flat:
    ax.xaxis.set_minor_formatter(NullFormatter()); ax.yaxis.set_minor_formatter(NullFormatter())
fig.tight_layout(h_pad=1.0, w_pad=1.2)
fs.save(fig, "fig3_sim_inversion", lang)
