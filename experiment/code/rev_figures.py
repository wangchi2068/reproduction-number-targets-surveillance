"""Figures for the 2026-11 revision analyses (matplotlib; data-computed only).

fig_rev_intercept_sens.pdf/png  forest plot of a under the primary fit and the S-a1..S-a9 sensitivities, plus E3
fig_rev_season.pdf/png          season intercepts with intervals; calendar-week permutation null of W; power curve
fig_rev_implied_sd.pdf/png      implied relative SD of R by mechanism (grid, not estimate)
fig_rev_rho_grid.pdf/png        rho x k heatmaps of E(a_hat)/v_R and b/b_P
Usage: python code/rev_figures.py
"""
from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import rev_common as R
from ma_models import Log

plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "legend.fontsize": 7, "figure.dpi": 150,
                     "axes.spines.top": False, "axes.spines.right": False, "font.family": "serif",
                     "font.serif": ["Liberation Serif", "Times New Roman", "DejaVu Serif"], "mathtext.fontset": "stix",
                     "pdf.fonttype": 42, "ps.fonttype": 42})
BLUE, RED, GREY, ORANGE = "#1f4e79", "#c0392b", "#7f7f7f", "#d68910"


def load(name):
    return json.loads((R.AN / name).read_text())


def save(fig, stem):
    R.FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(R.FIG / f"{stem}.{ext}", bbox_inches="tight", dpi=300 if ext == "png" else None)
    plt.close(fig)


SENS_LABELS = {"gi": "GI", "ge 90": "$\\geq$ 90%", "floor": "count floor", "drop": "drop"}
FAMILY_LABELS = {"shipped (fixed field)": "Original generators\n(fixed field)", "fresh-draw generator": "Fresh-draw\ngenerators",
                 "rho x k grid (S1 field)": "$\\rho\\times k$ grid\n(S1 field)", "daily renewal (two generations per week)": "Daily renewal\n(two generations/week)",
                 "generation-interval variant (r = 1)": "Generation-interval\nvariants ($r=1$)", "within-week estimand (E3, r = 1)": "Within-week\nestimand (E3, $r=1$)"}


def sens_label(k):
    s = k.replace("_", " ")
    for a, b in SENS_LABELS.items():
        s = s.replace(a, b)
    return s


def fig_intercept(ri):
    e1 = ri["E1_primary"]
    items = [("Primary (all rows)", e1["a"], e1["inference"]["two_way"]["a_ci"], BLUE)]
    for k, s in ri["sensitivities"].items():
        if k.startswith("S-a9"):
            continue
        items.append((sens_label(k), s["a"], s["inference"]["two_way"]["a_ci"], GREY))
    s9 = ri["sensitivities"]["S-a9_leave_one_week_out"]
    items.append(("S-a9 leave one week out (range; min lower limit)", np.mean(s9["a_range"]), [s9["lower_limit_min"], s9["a_range"][1]], GREY))
    e3 = ri["E3_within_week"]["week_leave_state_out"]
    items.append(("E3 within-week (estimand change)", e3["a"], e3["a_boot"]["ci"], ORANGE))
    fig, ax = plt.subplots(figsize=(6.3, 5.2))
    y = np.arange(len(items))[::-1]
    for yi, (lab, a, ci, c) in zip(y, items):
        if lab.startswith("S-a9"):
            ax.plot(s9["a_range"], [yi, yi], color=c, lw=4, alpha=.5, solid_capstyle="butt")
            ax.plot([s9["lower_limit_min"]], [yi], "|", color=c, ms=8)
            continue
        ax.errorbar(a, yi, xerr=[[a - ci[0]], [ci[1] - a]], fmt="o", color=c, ms=3.5, capsize=2, lw=1)
    ax.axvline(e1["a"], color=BLUE, ls=":", lw=.8)
    ax.axvline(0, color="k", lw=.8)
    ax.set_yticks(y)
    ax.set_yticklabels([i[0] for i in items])
    ax.set_xlabel("Count-independent intercept $a$ (relative variance), 95% CI")
    ax.set_xlim(left=-0.001)
    save(fig, "fig_rev_intercept_sens")


def fig_season(rs):
    e, pw = rs["effects"], rs["power"]
    null = np.load(R.AN / "rev_season_null_primary.npy")
    fig, axs = plt.subplots(1, 3, figsize=(9.6, 3.0), gridspec_kw={"width_ratios": [1.1, 1, 1]})
    ax = axs[0]
    ss = list(e["seasons"])
    for i, s in enumerate(ss):
        v = e["seasons"][s]
        ax.errorbar(i - .1, v["a_s"], yerr=[[v["a_s"] - v["ci_sandwich_log_delta"][0]], [v["ci_sandwich_log_delta"][1] - v["a_s"]]],
                    fmt="o", color=BLUE, ms=4, capsize=2, label="Two-way sandwich (log scale)" if i == 0 else None)
        ax.errorbar(i + .1, v["a_s"], yerr=[[v["a_s"] - v["ci_bootstrap"][0]], [v["ci_bootstrap"][1] - v["a_s"]]],
                    fmt="s", color=ORANGE, ms=3.5, capsize=2, label="Season-stratified block bootstrap" if i == 0 else None)
    ax.axhline(e["a_pooled"], color=GREY, ls=":", lw=.8, label="Pooled $a$")
    ax.set_xticks(range(len(ss)))
    ax.set_xticklabels(ss)
    ax.set_ylabel("Season intercept $a_s$, 95% CI")
    ax.set_ylim(0, .041)
    ax.set_title("A. Season intercepts (descriptive)")
    ax.legend(loc="upper left", frameon=False)
    ax = axs[1]
    w = rs["W_obs"]
    hi = np.quantile(null, .995)
    ax.hist(null, bins=60, range=(0, hi), color=GREY, alpha=.7, density=True)
    ax.axvline(w, color=RED, lw=1.2)
    ax.axvline(np.quantile(null, .95), color=BLUE, ls="--", lw=.9)
    ax.text(w, ax.get_ylim()[1] * .9, f"observed W = {w:.1f} \np = {rs['primary']['p']:.3f} ", color=RED, va="top", ha="right")
    ax.text(np.quantile(null, .95), ax.get_ylim()[1] * .55, f" 95th pct = {np.quantile(null, .95):.1f}", color=BLUE, va="top")
    ax.set_xlabel("Cluster-robust Wald statistic W (3 df)")
    ax.set_ylabel("Density")
    ax.set_title("B. Calendar-week permutation null")
    ax = axs[2]
    c = pw["curve"]
    xs = [r["max_over_min_ratio"] for r in c]
    ys = [r["power"] for r in c]
    se = [r["mcse"] for r in c]
    ax.errorbar(xs, ys, yerr=[1.96 * s for s in se], fmt="o-", color=BLUE, ms=3, capsize=1.5, lw=1)
    ax.axhline(.8, color=GREY, ls=":", lw=.8)
    ax.axhline(.05, color=GREY, ls=":", lw=.8)
    ax.axvline(pw["observed_max_over_min"], color=RED, lw=1)
    mdr = pw["mdr_max_over_min_at_80pct_power"]
    if mdr:
        ax.axvline(mdr, color=BLUE, ls="--", lw=.9)
        ax.text(mdr, .15, f" MDR {mdr:.1f}", color=BLUE)
    ax.text(pw["observed_max_over_min"], .93, f"observed {pw['observed_max_over_min']:.2f} \npower {pw['power_at_observed_pattern_delta_1']:.2f} ", color=RED, va="top", ha="right")
    ax.set_xscale("log")
    ax.set_xticks([1, 2, 3, 4, 6])
    ax.set_xticklabels(["1", "2", "3", "4", "6"])
    ax.set_xlabel("Highest/lowest season intercept ratio")
    ax.set_ylabel("Power of the primary test")
    ax.set_ylim(0, 1)
    ax.set_title("C. Design-based power")
    fig.tight_layout()
    save(fig, "fig_rev_season")


def fig_implied(rid):
    rows = rid["rows"]
    fams = list(dict.fromkeys(r["source"] for r in rows))
    fig, ax = plt.subplots(figsize=(6.3, 3.6))
    rng = np.random.default_rng(R.SEED)
    for i, f in enumerate(fams):
        rs = [r for r in rows if r["source"] == f]
        for r in rs:
            yi = i + rng.uniform(-.28, .28)
            c = RED if r["reproduces_counting_slope"] else (GREY if r["reproduces_counting_slope"] is False else BLUE)
            mk = "^" if r["post_hoc_law"] else "o"
            ax.plot(np.array(r["implied_rel_sd_interval"]) * 100, [yi, yi], color=c, lw=.6, alpha=.6)
            ax.plot(r["implied_rel_sd"] * 100, yi, mk, color=c, ms=3)
    env = rid["envelope_all"]
    ax.axvspan(env["point_range"][0] * 100, env["point_range"][1] * 100, color=BLUE, alpha=.06)
    ax.set_yticks(range(len(fams)))
    ax.set_yticklabels([FAMILY_LABELS.get(f, f) for f in fams])
    for t in ax.get_yticklabels():
        t.set_linespacing(1.6)
    ax.invert_yaxis()
    ax.set_xlabel("Implied relative SD of the reproduction number (%), with propagated interval")
    from matplotlib.lines import Line2D
    h = [Line2D([], [], color=RED, marker="o", ls="", ms=3, label="Mechanism reproduces observed $b/b_P$"),
         Line2D([], [], color=GREY, marker="o", ls="", ms=3, label="Does not reproduce $b/b_P$"),
         Line2D([], [], color=BLUE, marker="o", ls="", ms=3, label="No counting-slope check ($r=1$)"),
         Line2D([], [], color="k", marker="^", ls="", ms=3, label="Post hoc law")]
    ax.legend(handles=h, loc="lower right", frameon=False)
    save(fig, "fig_rev_implied_sd")


def fig_rho(gr):
    rhos, ks = gr["design"]["rho"], gr["design"]["k"]
    cell = {(c["rho"], c["k"]): c for c in gr["cells"]}
    A = np.array([[cell[(r, k)]["a_over_vR"]["ratio_of_means"] for k in ks] for r in rhos])
    B = np.array([[cell[(r, k)]["b_over_bP"]["mean"] for k in ks] for r in rhos])
    rep = np.array([[cell[(r, k)]["reproduces_counting_slope"] for k in ks] for r in rhos])
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.3))
    for ax, M, title, cmap, vmin, vmax in ((axs[0], A, r"A. Intercept recovery $E(\hat a)/v_R$", "RdBu_r", .5, 1.5),
                                          (axs[1], B, r"B. Counting slope $b/b_P$ (observed 1.257)", "RdBu_r", .5, 1.5)):
        im = ax.imshow(M, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
        for i in range(len(rhos)):
            for j in range(len(ks)):
                ax.text(j, i, f"{M[i, j]:.2f}" + ("*" if ax is axs[1] and rep[i, j] else ""), ha="center", va="center", fontsize=6.5,
                        color="white" if abs(M[i, j] - 1) > .3 else "black")
        ax.set_xticks(range(len(ks)))
        ax.set_xticklabels(["Poisson" if k is None else f"k={k:g}" for k in ks])
        ax.set_yticks(range(len(rhos)))
        ax.set_yticklabels([f"{r:g}" for r in rhos])
        ax.axhspan(-.5, 2.5, fill=False, ec="k", lw=1.2)
        ax.set_xlabel("Offspring dispersion")
        ax.set_ylabel(r"Reporting probability $\rho$")
        ax.set_title(title)
        fig.colorbar(im, ax=ax, shrink=.8)
    fig.tight_layout()
    save(fig, "fig_rev_rho_grid")


def main():
    log = Log(R.LOGS / "rev_figures.log")
    fig_intercept(load("rev_intercept.json"))
    fig_season(load("rev_season.json"))
    fig_implied(load("rev_identification.json"))
    fig_rho(load("rev_rho_grid.json"))
    log("rev_figures done", sorted(p.name for p in R.FIG.glob("fig_rev_*")))


if __name__ == "__main__":
    main()
