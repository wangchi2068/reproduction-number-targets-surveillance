"""Fig 2 (2x2): reporting heterogeneity and budget design. Reads analysis/parts/s07_theory.json only
(blocks fig2_VH_m, fig2_mxH_tau, fig2_mstar, fig2_equal_budget). No recomputation.

A: V_H(m) = v_R + v_H + A_H/m for tau in {0, 0.1, 0.15, 0.2}; dotted lines are the floors v_R + v_H,
   v_H = tau^2 (Rbar^2 + v_R), A_H = Rbar + rho (1+tau^2)(Rbar^2+v_R)/k (Prop 3).
B: m_x,H = A_H/(v_R + v_H) against tau; dashed reference m_x = A/v_R (tau = 0); marker at tau = 0.15.
C: m* = sqrt(c0/c1 * m_x) against c0/c1 (log-log) for m_x = 21.3 and 49.2 (Prop 4).
D: equal-budget lever comparison (Prop 6 MSE with AR(1) environment phi = 0.5 and AR(1) reporting
   psi = 0.5): budget B = T (c0 + c1 m), c0/c1 = 100, baseline T = 10, m = 50 (B0 = 1500 c1).
   'More volume': T = 10 fixed, m = B/(10 c1) - c0/c1. 'More periods': m = 50 fixed, T = B/(c0 + 50 c1).
   Solid: tau = 0.15; dashed: tau = 0 (reporting standardised). MSE = env + rep + samp with
   env = v_R S(phi,T)/T^2, rep = tau^2 [Rbar^2 S(psi,T) + v_R S(phi psi,T)]/T^2, samp = A_H/(T m).
   At any budget multiple the vertical distance between curves is the equal-budget comparison.
Usage: python3 make_fig2.py --lang en|zh
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
d = json.loads((PARTS / "s07_theory.json").read_text())["blocks"]
LS = ["-", (0, (4, 1.5)), (0, (1, 1.2)), (0, (5, 1.5, 1, 1.5))]

fig, axs = plt.subplots(2, 2, figsize=(fs.FULL_W, fs.FULL_W * 0.64))

# ---------------------------------------------------------------- A
ax = axs[0, 0]
b = d["fig2_VH_m"]
m = np.array(b["m"])
for i, t in enumerate(b["tau"]):
    col = fs.DARK2[i]
    ax.plot(m, b["V"][str(t)], color=col, ls=LS[i], label=f"$\\tau={t:g}$")
    ax.plot([m[0], m[-1]], [b["floor"][str(t)]] * 2, color=col, lw=0.6, ls=":")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1, 1000); ax.set_ylim(0.03, 4)
fs.xlabel(ax, T("Expected reported cases per period, $m$ (cases)", "单期期望报告数（例）"))
fs.ylabel(ax, T("Error variance", "误差方差"))
fs.legend(ax, T("Reporting volatility", "报告波动"), loc="upper right")
fs.panel_letter(ax, "A", dx=-0.17)

# ---------------------------------------------------------------- B
ax = axs[0, 1]
b = d["fig2_mxH_tau"]
tau = np.array(b["tau"]); mx = np.array(b["m_cross_H"])
ax.plot(tau, mx, color=fs.DARK2[0], ls="-")
ax.axhline(b["m_cross"], color=fs.GREY, lw=0.7, ls="--")
fs.text(ax, 0.30, b["m_cross"] * 1.03, "$m_\\times$", color=fs.GREY, va="bottom")
i15 = int(np.argmin(np.abs(tau - 0.15)))
ax.plot([tau[i15]], [mx[i15]], "o", color="black", ms=4)
fs.text(ax, tau[i15] + 0.012, mx[i15] * 1.04, "$\\tau=0.15$", color="black", va="bottom")
ax.set_xlim(0, 0.4); ax.set_ylim(5, 60); ax.set_yticks([10, 20, 30, 40, 50, 60])
fs.xlabel(ax, T("Relative reporting volatility, $\\tau$", "相对报告波动"))
fs.ylabel(ax, T("Balance scale, $m_{\\times,H}$ (cases)", "广义平衡规模（例）"))
fs.panel_letter(ax, "B", dx=-0.17)

# ---------------------------------------------------------------- C
ax = axs[1, 0]
b = d["fig2_mstar"]
c = np.array(b["c0c1"])
for i, mxv in enumerate(b["m_x"]):
    ax.plot(c, b["m_star"][str(mxv)], color=fs.DARK2[i], ls=LS[i], label=f"$m_\\times={mxv:g}$")
ax.fill_between(c, b["m_star"]["21.3"], b["m_star"]["49.2"], color=fs.DARK2[2], alpha=0.12, lw=0)
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(1, 1000)
fs.xlabel(ax, T("Cost ratio per period to per case, $c_0/c_1$", "单期固定成本与单例成本之比"))
fs.ylabel(ax, T("Optimal volume, $m^*$ (cases)", "最优单期规模（例）"))
fs.legend(ax, T("Balance scale", "平衡规模"), loc="upper left")
fs.panel_letter(ax, "C", dx=-0.17)

# ---------------------------------------------------------------- D
ax = axs[1, 1]
b = d["fig2_equal_budget"]
x = np.array(b["budget_multiple"])
strat = [("volume", fs.DARK2[0], T("More volume per period", "扩大单期规模")),
         ("periods", fs.DARK2[1], T("More periods", "增加观测期数"))]
for key, col, _ in strat:
    ax.plot(x, np.array(b["series"][f"{key}_tau0.15"]) * 1e-3, color=col, ls="-")
    ax.plot(x, np.array(b["series"][f"{key}_tau0.0"]) * 1e-3, color=col, ls=(0, (4, 1.5)))
ax.set_xlim(1, 8); ax.set_ylim(0, 0.035)
fs.xlabel(ax, T("Budget (multiple of baseline)", "预算（基线倍数）"))
fs.ylabel(ax, T("Mean squared error", "均方误差"))
hs = [Line2D([], [], color=col, ls="-", label=lab) for _, col, lab in strat]
ht = [Line2D([], [], color="black", ls="-", label="$\\tau=0.15$"),
      Line2D([], [], color="black", ls=(0, (4, 1.5)), label="$\\tau=0$")]
lg = fs.legend(ax, T("Extra budget spent on", "追加预算用途"), handles=hs, loc="upper right")
ax.add_artist(lg)
fs.legend(ax, T("Reporting volatility", "报告波动"), handles=ht, loc="upper left")
fs.panel_letter(ax, "D", dx=-0.17)

fig.tight_layout(h_pad=1.0, w_pad=1.5)
fs.save(fig, "fig2_heterogeneity_budget", lang)
