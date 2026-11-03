"""Fig 3 (2x2): simulation validation of the apparent (weekly hospital-admission) series.
Reads analysis/parts/s08_sim.json only (blocks fig3_traj, fig3_bins, scenarios, meanfield).

A: one simulated state (baseline v_R = 0.04, k = 0.35, rho_h = 0.01): true weekly R (step function),
   apparent R_A = wl_R((C_{t+1}+0.5)/(C_t+0.5)) and growth factor G = C_{t+1}/D_t; grey band = observation
   weeks used in the fit; dashed grey = Rbar.
B: baseline scenario, pooled over 20 replicates: loss L = ((X - Xbar)/Xbar)^2 averaged in 10 quantile bins
   of 1/D (95% bootstrap intervals of the bin mean), with the fitted FGLS line a + b/D, for X = G and R_A.
C: recovered intercept a against the nominal v_R/Rbar^2 for 15 grid scenarios (rho_h = 0.01), colour = k,
   marker = statistic; solid black = 1:1; dashed grey = noise-free (D -> infinity) limit a_inf of each
   statistic from the mean-field run (trivial interpolation of stored values).
D: recovered b/a of G against predicted b/a, log-log, solid black = 1:1. Colour = mapping
   (nominal: week treated as one generation, A/v_R; statistic scale: A_Y/v_Y from the noise-free moments
   of G), marker = without / with the manuscript's b-delta denominator correction (Eq. 25, sum w^2 = 0.495).
Usage: python3 make_fig3.py --lang en|zh
"""
import json

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import fig_style as fs
from cc_common import PARTS

lang = fs.cli()
fs.setup(lang)
T = lambda en, zh: fs.T(lang, en, zh)  # noqa: E731
J = json.loads((PARTS / "s08_sim.json").read_text())
d = J["blocks"]

fig, axs = plt.subplots(2, 2, figsize=(fs.FULL_W, fs.FULL_W * 0.66))

# ---------------------------------------------------------------- A
ax = axs[0, 0]
b = d["fig3_traj"]
wk = np.array(b["week"])
ax.axvspan(b["obs_weeks"][0] - 0.5, b["obs_weeks"][1] + 0.5, color="#EEEEEE", lw=0, zorder=0)
ax.step(wk, b["R_true"], where="mid", color="black", lw=1.0, label=T("True $R_t$", "真实值"))
ax.plot(wk, b["RA"], "o-", color=fs.DARK2[0], ms=3, lw=0.9, label="$R_A$")
ax.plot(wk, b["G"], "s", ls=(0, (4, 1.5)), color=fs.DARK2[1], ms=3, lw=0.9, label="$G$")
ax.axhline(b["Rbar"], color=fs.GREY, lw=0.7, ls="--")
ax.set_xlim(wk[0] - 0.6, wk[-1] + 0.6); ax.set_ylim(0.6, 3.6)
ax.set_xticks(wk)
ax.set_xticklabels([str(int(w - wk[0] + 1)) for w in wk])
fs.xlabel(ax, T("Simulated week", "模拟周次"))
fs.ylabel(ax, T("$R$ or growth factor", "再生数或增长因子"))
fs.legend(ax, T("Series", "序列"), loc="upper left", ncol=3)
fs.panel_letter(ax, "A", dx=-0.17)

# ---------------------------------------------------------------- B
ax = axs[0, 1]
for i, (nm, lab, mk) in enumerate((("G", "$G$", "s"), ("RA", "$R_A$", "o"))):
    bb = d["fig3_bins"][nm]
    x = np.array(bb["invD"]); y = np.array(bb["L_mean"])
    lo, hi = np.array(bb["L_lo"]), np.array(bb["L_hi"])
    col = fs.DARK2[1] if nm == "G" else fs.DARK2[0]
    ax.errorbar(x, y, yerr=[y - lo, hi - y], fmt=mk, ms=3, color=col, mfc="white", elinewidth=0.6, label=lab)
    xx = np.linspace(0, 0.21, 100)
    ax.plot(xx, bb["a"] + bb["b"] * xx, color=col, lw=0.9, ls="-")
ax.set_yscale("log"); ax.set_xlim(0, 0.21); ax.set_ylim(0.008, 0.6)
fs.xlabel(ax, T("Inverse scale, $1/D_t$ (per case)", "规模倒数（每例）"))
fs.ylabel(ax, T("Mean relative squared error, $L$", "平均相对平方误差"))
fs.legend(ax, T("Statistic", "统计量"), loc="lower right")
fs.panel_letter(ax, "B", dx=-0.17)

# ---------------------------------------------------------------- C
ax = axs[1, 0]
sc = [r for r in d["scenarios"] if r["rho_h"] == 0.01]
kcol = {0.1: fs.DARK2[0], 0.35: fs.DARK2[1], 1.0: fs.DARK2[2]}
mk = {"G": "s", "RA": "o"}
for r in sc:
    for nm in ("G", "RA"):
        lo, hi = r[nm]["a_ci"]
        ax.errorbar([r["a_nom"]], [r[nm]["a"]], yerr=[[r[nm]["a"] - lo], [hi - r[nm]["a"]]], fmt=mk[nm], ms=3.5,
                    color=kcol[r["k"]], mfc="white" if nm == "RA" else kcol[r["k"]], elinewidth=0.6, alpha=0.9)
xs = np.array([0.004, 0.2])
ax.plot(xs, xs, color="black", lw=0.8)
mf = d["meanfield"]
vr = sorted(float(v) for v in mf)
Rb = d["design"]["Rbar"]
for nm in ("G", "RA"):
    ax.plot([v / Rb ** 2 for v in vr], [mf[f"{v:g}"][nm]["a_inf"] for v in vr], color=fs.GREY, lw=0.7, ls="--")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.005, 0.15); ax.set_ylim(0.002, 0.6)
fs.text(ax, 0.035, 0.045, "1:1", color="black", ha="right", va="bottom")
fs.xlabel(ax, T("Nominal intercept, $v_R/\\bar{R}^2$", "名义截距"))
fs.ylabel(ax, T("Recovered intercept, $a$", "回收截距"))
hk = [Line2D([], [], color=c, marker="s", ls="none", ms=3.5, label=f"$k={k:g}$") for k, c in kcol.items()]
lg = fs.legend(ax, T("Dispersion", "离散参数"), handles=hk, loc="upper left")
ax.add_artist(lg)
hs = [Line2D([], [], color="black", marker="s", ls="none", ms=3.5, label="$G$"),
      Line2D([], [], color="black", marker="o", mfc="white", ls="none", ms=3.5, label="$R_A$")]
fs.legend(ax, T("Statistic", "统计量"), handles=hs, loc="lower right")
fs.panel_letter(ax, "C", dx=-0.17)

# ---------------------------------------------------------------- D
ax = axs[1, 1]
allsc = d["scenarios"]
yy = np.array([r["G"]["ba"] for r in allsc])
lo = np.array([r["G"]["ba_ci"][0] for r in allsc]); hi = np.array([r["G"]["ba_ci"][1] for r in allsc])
series = [("nom", np.array([r["ba_nom"] for r in allsc]), fs.DARK2[0], "o", "white"),
          ("nom_d", np.array([r["G"]["ba_nom_delta"] for r in allsc]), fs.DARK2[0], "^", fs.DARK2[0]),
          ("stat", np.array([r["G"]["ba_stat"] for r in allsc]), fs.DARK2[1], "o", "white"),
          ("stat_d", np.array([r["G"]["ba_stat_delta"] for r in allsc]), fs.DARK2[1], "^", fs.DARK2[1])]
for _, x, col, m, fc in series:
    ax.errorbar(x, yy, yerr=[yy - lo, hi - yy], fmt=m, ms=3.5, color=col, mfc=fc, elinewidth=0.6)
xs = np.array([1, 400])
ax.plot(xs, xs, color="black", lw=0.8)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1.5, 300); ax.set_ylim(1, 220)
fs.text(ax, 115, 122, "1:1", color="black", ha="right", va="bottom")
fs.xlabel(ax, T("Predicted balance ratio", "理论平衡比"))
fs.ylabel(ax, T("Recovered balance ratio, $b/a$", "回收平衡比"))
hm = [Line2D([], [], color=fs.DARK2[0], marker="s", ls="none", ms=3.5, label=T("Nominal", "名义映射")),
      Line2D([], [], color=fs.DARK2[1], marker="s", ls="none", ms=3.5, label=T("Statistic scale", "统计量尺度"))]
lg = fs.legend(ax, T("Mapping", "映射"), handles=hm, loc="upper left")
ax.add_artist(lg)
hc = [Line2D([], [], color="black", marker="o", mfc="white", ls="none", ms=3.5, label=T("Uncorrected", "未修正")),
      Line2D([], [], color="black", marker="^", ls="none", ms=3.5, label=T("Denominator corrected", "分母噪声修正"))]
fs.legend(ax, T("Prediction", "理论值"), handles=hc, loc="lower right")
fs.panel_letter(ax, "D", dx=-0.17)

fig.tight_layout(h_pad=1.0, w_pad=1.5)
fs.save(fig, "fig3_sim_validation", lang)
