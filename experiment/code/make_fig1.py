"""Fig 1 (panels B-D; panel A is reserved for an image-gen concept schematic): theory of the
case-count bottleneck. Reads analysis/parts/s07_theory.json only (blocks fig1_V_m, fig1_sim_V_m,
fig1_V_rho, fig1_sim_V_rho, fig1_mx_contour). No recomputation.

B: V(m) = v_R + A/m (Target 2, marginal mean) and A/m (Target 1, conditional R_t) against m, log-log,
   with seeded Monte Carlo points (20,000 replications per m, 95% intervals); floor v_R, line m_x = A/v_R.
C: V(n, rho) = v_R + Rbar/(rho n) + (Rbar^2+v_R)/(k n) against rho at n = 40, 200, 1000, with MC points;
   the value at rho = 1 is the residual demographic floor (Remark 1).
D: contours of m_x = [Rbar + rho(Rbar^2+v_R)/k]/v_R in the (v_R, k) plane (Rbar = 1.15, rho = 0.25),
   baseline point marked.
Usage: python3 make_fig1.py --lang en|zh
"""
import json

import matplotlib.pyplot as plt
import numpy as np

import fig_style as fs
from cc_common import PARTS

lang = fs.cli()
fs.setup(lang)
T = lambda en, zh: fs.T(lang, en, zh)  # noqa: E731
d = json.loads((PARTS / "s07_theory.json").read_text())["blocks"]

fig, axs = plt.subplots(1, 3, figsize=(fs.FULL_W, fs.FULL_W * 0.36))
c1, c2, c3 = fs.DARK2[0], fs.DARK2[1], fs.DARK2[2]

# ---------------------------------------------------------------- B
ax = axs[0]
b, s = d["fig1_V_m"], d["fig1_sim_V_m"]
m = np.array(b["m"])
ax.plot(m, b["V"], color=c1, ls="-", label=T("Target 2: $\\bar{R}$", "跨期均值"))
ax.plot(m, b["target1"], color=c2, ls=(0, (4, 1.5)), label=T("Target 1: $R_t$", "当期值"))
ms = np.array(s["m"])
for key, ci, col in (("V_sim", "V_ci", c1), ("T1_sim", "T1_ci", c2)):
    y = np.array(s[key]); lo, hi = np.array(s[ci]).T
    ax.errorbar(ms, y, yerr=[y - lo, hi - y], fmt="o", ms=3, mfc="white", mec=col, ecolor=col, elinewidth=0.6)
ax.plot([], [], "o", ms=3, mfc="white", mec="black", label=T("Simulation", "模拟"))
ax.axhline(b["v_R"], color="black", lw=0.7, ls="--")
ax.axvline(b["m_cross"], color=fs.GREY, lw=0.7, ls="--")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(0.8, 1300); ax.set_ylim(1.2e-3, 4)
fs.text(ax, 1.0, b["v_R"] * 1.25, "$v_R$", color="black")
fs.text(ax, b["m_cross"] * 1.15, 1.6, "$m_\\times$", color=fs.GREY)
fs.xlabel(ax, T("Expected reported cases per period, $m$ (cases)", "单期期望报告数（例）"))
fs.ylabel(ax, T("Error variance", "误差方差"))
fs.legend(ax, T("Estimation target", "估计目标"), loc="lower left")
fs.panel_letter(ax, "B", dx=-0.25)

# ---------------------------------------------------------------- C
ax = axs[1]
b, s = d["fig1_V_rho"], d["fig1_sim_V_rho"]
rho = np.array(b["rho"])
for col, ls, n in zip((c1, c2, c3), ("-", (0, (4, 1.5)), (0, (1, 1.2))), b["n"]):
    ax.plot(rho, b["V"][str(n)], color=col, ls=ls, label=f"$n={n}$")
    y = np.array(s["V_sim"][str(n)]); lo, hi = np.array(s["V_ci"][str(n)]).T
    ax.errorbar(s["rho"], y, yerr=[y - lo, hi - y], fmt="o", ms=3, mfc="white", mec=col, ecolor=col,
                elinewidth=0.6)
ax.axhline(b["v_R"], color="black", lw=0.7, ls="--")
fs.text(ax, 0.04, b["v_R"] * 0.82, "$v_R$", color="black", va="top")
ax.set_yscale("log"); ax.set_xlim(0, 1.03); ax.set_ylim(0.025, 2.0)
fs.xlabel(ax, T("Reporting probability, $\\rho$", "报告概率"))
fs.ylabel(ax, T("Error variance", "误差方差"))
fs.legend(ax, T("Ancestor scale", "祖代规模"), loc="upper right")
fs.panel_letter(ax, "C", dx=-0.25)

# ---------------------------------------------------------------- D
ax = axs[2]
b = d["fig1_mx_contour"]
vg, kg, Z = np.array(b["v_R"]), np.array(b["k"]), np.array(b["m_cross"])
levels = [10, 20, 50, 100, 200]
ax.contour(vg, kg, Z, levels=levels, colors=[fs.DARK2[2]], linewidths=0.8)
# level labels placed where each contour crosses k = 1.5 (closed form: v = [Rbar + rho Rbar^2/k]/(L - rho/k))
kl, rb, rh = 1.5, b["Rbar"], b["rho"]
for L in levels:
    v = (rb + rh * rb ** 2 / kl) / (L - rh / kl)
    fs.text(ax, v, kl, f"{L}", color=fs.DARK2[2], ha="center", va="center",
            bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none"))
# what the contour numbers mean (single-script text artist per language)
fs.text(ax, 0.5, 0.97, T("Contours: balance scale,\n$m_\\times$ (cases)", "等值线：平衡规模（例）"),
        transform=ax.transAxes, color="black", ha="center", va="top", ma="center",
        bbox=dict(boxstyle="square,pad=0.25", fc="white", ec="0.6", lw=0.5))
bl = b["baseline"]
ax.plot([bl["v_R"]], [bl["k"]], "o", color="black", ms=4)
fs.text(ax, bl["v_R"] * 1.15, bl["k"] * 1.12, T("Baseline", "基线"), color="black", va="bottom")
ax.set_xscale("log"); ax.set_yscale("log")
fs.xlabel(ax, T("Environmental variance, $v_R$", "环境方差"))
fs.ylabel(ax, T("Offspring dispersion, $k$", "离散参数"))
fs.panel_letter(ax, "D", dx=-0.25)

fig.tight_layout(w_pad=1.2)
fs.save(fig, "fig1_theory", lang)
