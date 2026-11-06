"""Merge, tables and figures for Amendment 2 (analysis_plan_amendment2.md, Section H). No new computation.

Reads analysis/rev_r2_{floor,total_bias,generator,reporting,clamp,bootstrap_sd,poisson_ref,interval_types}.json and
analysis/rev_intercept.json; writes analysis/results.json key revision_2026_11_r2 (other keys untouched; the
analysis_plan field of revision_2026_11 is set to the published plan paths), tables/tab_r2_*.tex and
figures/fig_r2_{deciles,population,floor}.{pdf,png}. Derived numbers: b/b_P intervals under the alternative Poisson
references (the two-way interval of b divided by each reference, b_P treated as fixed, as for the shipped reference)
and, in fig_r2_population, the log-link P90/P10 ratio interval exp(theta_ci (z_p90 - z_p10)) mapped from the stored
state-bootstrap interval of theta (monotone map, so the percentile interval maps exactly).
Usage: python code/rev_r2_render.py
"""
from __future__ import annotations

import json
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import rev_common as R
import rev_r2_lib as L2
from ma_models import Log
from rev_tables import esc, pfmt, table as _table


def table(*args, size=r"\scriptsize\renewcommand{\arraystretch}{1.15}", **kw):
    """Shared table writer with extra row spacing: at \\scriptsize the minus signs of adjacent rows touch."""
    return _table(*args, size=size, **kw)

plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "legend.fontsize": 7, "figure.dpi": 150,
                     "axes.spines.top": False, "axes.spines.right": False, "font.family": "serif",
                     "font.serif": ["Liberation Serif", "Times New Roman", "DejaVu Serif"], "mathtext.fontset": "stix",
                     "pdf.fonttype": 42, "ps.fonttype": 42})
FIG = R.EXP / "figures"
DAG = r"$^\dagger$"
PLAN = {"plan": "analysis/analysis_plan.md", "lock": "analysis/analysis_plan.lock", "plan_locked_utc": "2026-11-05 05:42:50",
        "lock_scope": "sha256 of the first 14,195 bytes (the plan up to '(none at lock time)'; Amendment 1 was appended after the lock)",
        "amendment_1": "analysis/analysis_plan.md, section Amendment 1 (implementation details)",
        "amendment_2": "analysis/analysis_plan_amendment2.md", "amendment_2_lock": "analysis/analysis_plan_amendment2.lock",
        "amendment_2_locked_utc": "2026-11-05 10:16:29"}
BLOCKS = [("B_floor", "rev_r2_floor.json"), ("C_total_bias", "rev_r2_total_bias.json"), ("D1_state_realistic_generators", "rev_r2_generator.json"),
          ("D2_reporting_data", "rev_r2_reporting.json"), ("E_clamp", "rev_r2_clamp.json"), ("F1_bootstrap_sd", "rev_r2_bootstrap_sd.json"),
          ("F2_poisson_reference", "rev_r2_poisson_ref.json"), ("G_interval_types", "rev_r2_interval_types.json")]


def load(name):
    return json.loads((R.AN / name).read_text())


def num(v, d=4):
    if v is None or not np.isfinite(v):
        return "--"
    return f"$-${abs(v):.{d}f}" if v < 0 else f"{v:.{d}f}"


def ci(v, d=4):
    if v is None or v[0] is None:
        return "--"
    if min(v) < 0:
        return f"{num(v[0], d)} to {num(v[1], d)}"
    return f"{v[0]:.{d}f}--{v[1]:.{d}f}"


def pc(v, d=1):
    return num(100 * v, d)


def pci(v, d=1):
    return ci([100 * v[0], 100 * v[1]], d) if v else "--"


def yes(b):
    return "yes" if b else "no"


# ------------------------------------------------------------------ merge
def merge(js, b_ci):
    path = R.AN / "results.json"
    res = json.loads(path.read_text())
    if "revision_2026_11" in res:
        res["revision_2026_11"]["analysis_plan"] = PLAN
    out = {"generated_utc": time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()), "analysis_plan": PLAN,
           "data": "frozen NHSN 2026-10-20 vintage; R_A sample, 1,748 rows, 51 states, 42 calendar weeks",
           "seed_root": R.SEED}
    for key, name in BLOCKS:
        blk = js[name]
        out[key] = {"source": blk["meta"]["script"], **blk}
    pr = out["F2_poisson_reference"]
    pr["b_ci_two_way_t41"] = b_ci
    pr["b_over_bP_ci_two_way_bP_fixed"] = {k: [b_ci[0] / v, b_ci[1] / v] for k, v in pr["references"].items()}
    res["revision_2026_11_r2"] = out
    R.write_json(path, res)
    return out


# ------------------------------------------------------------------ tables
def t_deciles(fl):
    d = fl["B1_deciles"]
    rows = []
    for r in d["rows"]:
        rows.append(f"{r['decile']} & {r['C_min']:.0f}--{r['C_max']:.0f} & {r['n']} & {r['mean_C']:.1f} & {r['mean_L']:.4f} & {ci(r['mean_L_ci'])} & "
                    f"{r['affine']:.4f} & {r['power']:.4f} & {r['three']:.4f} & {yes(r['affine_inside_meanL_ci'])}")
    ft = d["fits"]
    notes = (r"Rows of $R_A$ (1,748) split at the deciles of the current-week count $C_t$. Mean $L$: observed mean relative squared deviation; "
             r"95\% interval from a calendar-week block bootstrap (1,999 draws). Model columns give the mean prediction in the decile under "
             rf"the affine fit $a+b/C$ ($a={ft['affine']['a']:.4f}$, $b={ft['affine']['b']:.3f}$), the power law $cC^{{-\gamma}}$ "
             rf"($\gamma={ft['power']['gamma']:.3f}$) and the three-parameter fit $a+bC^{{-\gamma}}$ ($a={ft['three_parameter']['a']:.4f}$, "
             rf"$\gamma={ft['three_parameter']['gamma']:.3f}$). Quasi-deviance: affine {ft['affine']['qdev']:.1f}, power {ft['power']['qdev']:.1f}, "
             rf"three-parameter {ft['three_parameter']['qdev']:.1f}. Last column: whether the affine prediction lies inside the interval of the observed mean.")
    table("tab_r2_deciles.tex", r"Observed and fitted loss by decile of the current-week count", "tab:r2-deciles", r"rrrrrlrrrr",
          r"Decile & $C_t$ & $n$ & Mean $C_t$ & Mean $L$ & 95\% interval & Affine & Power & Three-par. & Affine inside", rows, notes)


def t_population(fl):
    p = fl["B2_population"]
    li, lg, bv = p["primary_linear"], p["secondary_log_link"], p["secondary_b_varies"]
    cv = p["covariate"]
    rows = [r"\multicolumn{5}{@{}l}{\textit{Linear intercept $a_0+a_1z$, common slope $b$ (primary)}}",
            rf"$a_0$ & {li['a0']:.4f} & {ci(li['a0_ci_two_way_t41'])} & -- & --",
            rf"$a_1$ (per SD of log population) & {num(li['a1'])} & {ci(li['a1_ci_two_way_t41'])} & {ci(li['a1_ci_state_bootstrap'])} & {ci(li['a1_ci_week_bootstrap'])}",
            rf"$a$ at population P10 / P90 & {li['a_at_pop_p10']:.4f} / {li['a_at_pop_p90']:.4f} & -- & -- & --",
            rf"Ratio $a$(P90) / $a$(P10) & {li['ratio_a_p90_over_p10']:.3f} & -- & {ci(li['ratio_ci_state_bootstrap'], 3)} & --",
            rf"$b$ & {li['b']:.3f} & -- & -- & --",
            r"\midrule",
            r"\multicolumn{5}{@{}l}{\textit{Log-link intercept $\exp(\alpha_0+\theta z)$ (secondary)}}",
            rf"$\theta$ & {num(lg['theta'], 3)} & -- & {ci(lg['theta_ci_state_bootstrap'], 3)} & {ci(lg['theta_ci_week_bootstrap'], 3)}",
            rf"Ratio $a$(P90) / $a$(P10) & {lg['ratio_a_p90_over_p10']:.3f} & -- & {ci(sorted(float(np.exp(t * (cv['z_p90'] - cv['z_p10']))) for t in lg['theta_ci_state_bootstrap']), 3)} & --",
            r"\midrule",
            r"\multicolumn{5}{@{}l}{\textit{Intercept and slope both linear in $z$ (secondary)}}",
            rf"$a_1$ & {num(bv['beta'][1])} & {ci(bv['a1_ci_two_way_t41'])} & -- & --",
            rf"$b_1$ & {num(bv['beta'][3], 3)} & {ci(bv['b1_ci_two_way_t41'], 3)} & -- & --",
            r"\midrule",
            r"\multicolumn{5}{@{}l}{\textit{Population tertiles, separate affine fits (two-way interval, $t$ with $G-1$ df)}}"]
    for t in p["tertiles_separate"]:
        rng = f"{t['pop_range'][0] / 1e6:.2f}--{t['pop_range'][1] / 1e6:.1f} M"
        rows.append(rf"Tertile {t['tertile']} ({len(t['states'])} states, {rng}, {t['n_rows']} rows): $a$ & {t['a']:.4f} & {ci(t['a_ci_two_way'])} & -- & --")
        rows.append(rf"\quad $b$ (median $C_t$ = {t['median_C']:.0f}) & {t['b']:.3f} & {ci(t['b_ci_two_way'], 3)} & -- & --")
    tc = p["tertiles_common_slope"]
    rows += [r"\midrule", r"\multicolumn{5}{@{}l}{\textit{Population tertiles, common slope}}"]
    for i, (a, c) in enumerate(zip(tc["intercepts"], tc["intercept_ci_state_bootstrap"]), 1):
        rows.append(rf"Tertile {i} intercept & {a:.4f} & -- & {ci(c)} & --")
    rows.append(rf"Ratio tertile 3 / tertile 1 & {tc['ratio_high_over_low']:.3f} & -- & {ci(tc['ratio_ci_state_bootstrap'], 3)} & --")
    notes = (rf"$z$: log of the mean 2022--2025 resident population (Census Vintage 2025), standardised over the 51 jurisdictions "
             rf"(P10 = {cv['pop_p10'] / 1e6:.2f} M, P90 = {cv['pop_p90'] / 1e6:.1f} M). Primary test of $a_1=0$: two-way sandwich, $t$ with 41 df, "
             rf"$t={li['t_a1']:.2f}$, $p={li['p_a1_two_sided']:.4f}$. State bootstrap: states resampled (1,999 draws; 999 for the tertile fits); "
             r"week bootstrap: calendar-week blocks (1,999 draws). All intervals 95\%.")
    table("tab_r2_population.tex", r"Dependence of the count-independent intercept on state population", "tab:r2-population", r"Xrlll",
          r"Quantity & Estimate & Two-way & State bootstrap & Week bootstrap", rows, notes)


def t_floor(fl, ri):
    b3, b4 = fl["B3_three_parameter"], fl["B4_model_free_floor"]
    pt, wb, sb = b3["point"], b3["week_bootstrap"], b3["state_bootstrap"]
    e1 = ri["E1_primary"]
    rows = [r"\multicolumn{4}{@{}l}{\textit{Affine fit $a+b/C$ (primary, for reference)}}",
            rf"$a$ (two-way interval) & {e1['a']:.4f} & {ci(e1['inference']['two_way']['a_ci'])} & --",
            r"\midrule", r"\multicolumn{4}{@{}l}{\textit{Three-parameter fit $a+bC^{-\gamma}$}}",
            rf"$a$ & {pt['a']:.4f} & {ci(wb['a']['ci'])} & {ci(sb['a']['ci'])}",
            rf"$b$ & {pt['b']:.3f} & {ci(wb['b']['ci'], 3)} & --",
            rf"$\gamma$ & {pt['gamma']:.3f} & {ci(wb['gamma']['ci'], 3)} & {ci(sb['gamma']['ci'], 3)}",
            rf"Share of draws with $a\le 0$ & -- & {wb['share_a_le_0']:.3f} & {sb['share_a_le_0']:.3f}",
            rf"Quasi-deviance gain over affine & {b3['affine_qdev_minus_three']:.2f} & share $>2$: {wb['share_qdev_gain_gt_2']:.3f} & --",
            r"\midrule", r"\multicolumn{4}{@{}l}{\textit{Model-free floor among rows with $C_t$ at or above a quantile}}"]
    for r in b4["rows"]:
        q = int(round(100 * r["threshold_quantile"]))
        rows += [rf"P{q} ($C_t\ge {r['threshold_C']:.0f}$, $n={r['n']}$, mean $C_t$ {r['mean_C']:.0f}): mean $L$ & {r['mean_L']:.4f} & {ci(r['mean_L_ci'])} & --",
                 rf"\quad mean $L-\hat b/C_t$ & {r['floor_minus_bhat']:.4f} & {ci(r['floor_minus_bhat_ci'])} & --",
                 rf"\quad mean $L-b_P/C_t$ & {r['floor_minus_bP']:.4f} & {ci(r['floor_minus_bP_ci'])} & --",
                 rf"\quad relative SD $(L-\hat b/C_t)^{{1/2}}$ (\%) & {pc(r['floor_rel_sd_bhat'])} & {pci(r['floor_rel_sd_bhat_ci'])} & --"]
    notes = (r"Week bootstrap: calendar-week blocks (1,999 draws); state bootstrap: states resampled (999 draws). Model-free floor: mean loss "
             r"of the rows above the threshold minus the mean counting term, with $\hat b$ the fitted slope or $b_P$ the Poisson reference; "
             r"intervals from the calendar-week block bootstrap with $\hat b$ refitted in each draw. "
             rf"Reading rule of Amendment 2: the floor survives a free exponent if the week-bootstrap interval of $a$ excludes 0 "
             rf"(here: {yes(b3['floor_survives_free_exponent'])}).")
    table("tab_r2_floor.tex", r"The floor under a free exponent and without a model", "tab:r2-floor", r"Xrll",
          r"Quantity & Estimate & Week bootstrap & State bootstrap", rows, notes)


def t_bias(tb):
    rows = []
    for name, g in tb["generators"].items():
        lab = name.replace("FD-S8-share", r"FD-S8, share ").replace("FD-S9-", "FD-S9, ")
        wh = f"{g['w_H']:.3f} ({ci(g['w_H_ci'], 3)})" if "w_H" in g else "--"
        rows.append(rf"{lab}{DAG if g.get('post_hoc_law') else ''} & {pc(g['share_vH_of_total'])} & {g['mean_a']:.4f} & {g['mean_vR']:.4f} & {g['mean_vH']:.4f} & "
                    rf"{pc(g['relative_bias'])} ({pci(g['relative_bias_ci'])}) & {g['a_over_vR']:.3f} & {wh}")
    m = tb["mechanism_conditional_rescaling"]
    notes = (r"Fresh-draw generators of the first revision (200 panels each) and an FD-S8 sweep over the transmission share of the calibrated "
             r"variance (shares 1.0, 0.75, 0.25; 0.5 is FD-S8 itself). $v_R$, $v_H$: true relative variance of transmission and of the reporting "
             r"shock on the kept rows; relative bias $=E(\hat a)/(v_R+v_H)-1$; $w_H$: weight on $v_H$ solving $E(\hat a)=v_R+w_Hv_H$. "
             r"Intervals: panel bootstrap (2,000 draws). $^\dagger$ post hoc law. Wording rule of Amendment 2: if every interval excludes 0, "
             rf"the manuscript says $\hat a$ reflects, not estimates, the total (outcome: {esc(tb['wording_rule']['verdict'])}). "
             rf"Mechanism-conditional rescaling $a_{{\mathrm{{obs}}}}/t_g$ would range {m['a_obs_over_t_g_range'][0]:.4f}--{m['a_obs_over_t_g_range'][1]:.4f}; not used.")
    table("tab_r2_bias.tex", r"Bias of the intercept relative to the total non-counting variance", "tab:r2-bias", r"Xrrrrlrl",
          r"Generator & $v_H$ share (\%) & $E(\hat a)$ & $v_R$ & $v_H$ & Rel. bias (\%) & $E(\hat a)/v_R$ & $w_H$", rows, notes)


def t_generator(gj):
    rows = []
    o = gj["observed"]
    for name, c in gj["cells"].items():
        t, cal = c["tests"], c["calibration"]
        mark = lambda k: "" if t[k]["pass"] else "*"
        sig = "failed" if cal["calibration_failed"] else f"{cal['sigma']:.3f}"
        rows.append(rf"{name} & {sig} & "
                    rf"{ci(t['T1_b_over_bP']['range'], 2)}{mark('T1_b_over_bP')} & {t['T2_p90']['mc_mean']:.0f}{mark('T2_p90')} & "
                    rf"{ci(t['T3_pooled_over_within']['range'], 2)}{mark('T3_pooled_over_within')} & {ci(t['T4_state_cv']['range'], 2)}{mark('T4_state_cv')} & "
                    rf"{c['a_over_vR']:.2f} & {c['sd_ratio_observed_se_over_mc_sd']:.2f} & {c['coverage_of_E_a']['coverage']:.3f}")
    rows.append(rf"Observed & -- & {o['b_over_bP']:.2f} & {o['p90']:.0f} & {o['pooled_over_within']:.2f} & {o['state_cv']:.2f} & -- & -- & --")
    notes = (r"State-realistic generators (200 panels per cell): transmission field around a smoothed observed trend with state-specific variance "
             r"(multiplier proportional to the population-dependent intercept, or to the state fixed-effect intercepts for SR-fe), clustered "
             r"admissions with extra-Poisson parameter $\theta_c$ (c0, c25, c50 $=0, 0.25, 0.5$), reporting probability 0.01 (0.003 for r003), "
             r"trend exponent 1 (0.5 for t50). $\sigma$: transmission SD scale; ``failed'' when the smoothed trend alone exceeds the transmission "
             r"target ($\sigma=0$). T1--T4: Monte Carlo 2.5--97.5\% range of $b/b_P$, mean P90 of $C_t$ (band 350--428), range of pooled/within "
             r"slope ratio and of the between-state CV of intercepts; * the observed value lies outside. Success rule of Amendment 2: all four pass in one "
             rf"cell (outcome: {esc(gj['success_rule']['verdict'])}). Last columns: $E(\hat a)/v_R$, observed SE over Monte Carlo SD of $\hat a$, "
             r"coverage of $E(\hat a)$ by the two-way interval. Data-informed generator design; no inferential result uses these generators.")
    table("tab_r2_generator.tex", r"State-realistic fresh-draw generators against four features of the observed panel", "tab:r2-generator",
          r"Xlllllrrr", r"Cell & $\sigma$ & T1 $b/b_P$ & T2 P90 & T3 pooled/within & T4 state CV & $E(\hat a)/v_R$ & SE/SD & Coverage", rows, notes)


def t_clamp(cl):
    rows = []
    lab = {"clamp_1e-09": r"clamp at $10^{-9}$", "clamp_1e-07": r"clamp at $10^{-7}$ (primary)", "clamp_1e-05": r"clamp at $10^{-5}$",
           "clamp_1e-03": r"clamp at $10^{-3}$", "no_clamp": "no clamp", "exclude_clamped_draws": "exclude draws with any clamped row"}
    for k, v in cl["variants"].items():
        if not v["evaluable"]:
            rows.append(rf"{lab[k]} & {v['valid']} & not evaluable & -- & -- & --")
            continue
        q = v["t_star_q"]
        rows.append(rf"{lab[k]} & {v['valid']:,} & {pfmt(v['p_two_sided'])} & {v['p_min_attainable']:.4f} & {pc(v['row_share_clamped'])} & "
                    rf"{num(q[0], 2)} / {q[1]:.2f} / {q[2]:.2f}")
    notes = (rf"Wild cluster bootstrap-$t$ of the intercept (9,999 draws, calendar-week Webb weights, seed {cl['meta']['seed'][1]}); $t_{{\mathrm{{obs}}}}"
             rf"={cl['meta']['t_obs']:.2f}$. Valid: draws with a converged fit. Minimum attainable $p$: $1/(\text{{valid}}+1)$. Clamped rows: share of "
             r"bootstrap rows whose loss is at or below the clamp (for no clamp: negative losses, kept). $t^*$: 2.5\%, 50\%, 97.5\% quantiles. "
             r"Every draw contains at least one clamped row, so the exclusion variant has no draws. Secondary calibration; it does not change the E1 verdict.")
    table("tab_r2_clamp.tex", r"Wild bootstrap of the intercept under alternative loss clamps", "tab:r2-clamp", r"Xrllll",
          r"Variant & Valid draws & $p$ (two-sided) & Min.\ $p$ & Clamped rows (\%) & $t^*$ quantiles", rows, notes)


def t_poisson(pr, b_ci):
    lab = {"shipped_plugin": "plug-in, as shipped", "normalised_plugin": "plug-in, mean-normalised", "simulated_exact": "simulated exact (Poisson-only panels)"}
    rows = []
    for k, v in pr["references"].items():
        rows.append(rf"{lab[k]} & {v:.4f} & {pr['b_over_bP'][k]:.3f} & {ci([b_ci[0] / v, b_ci[1] / v], 2)} & {pc(pr['relative_change_vs_shipped_ratio'][k])}")
    s = pr["simulated"]
    notes = (rf"$\hat b={pr['b_hat']:.3f}$, two-way 95\% interval {ci(b_ci, 3)}; intervals of $b/b_P$ divide it by the reference ($b_P$ fixed). "
             rf"Simulated exact: {pr['meta']['B']} panels with Poisson counts around the observed means, same transform, row filter ($C_t\ge 5$) and "
             rf"affine fit; pooled slope {s['pooled_slope']:.4f} (per-panel {s['per_panel_slope_range'][0]:.3f}--{s['per_panel_slope_range'][1]:.3f}), "
             rf"pooled intercept {s['pooled_intercept']:.1e}.")
    table("tab_r2_poisson.tex", r"Counting slope against three Poisson references", "tab:r2-poisson", r"Xrrlr",
          r"Reference & $b_P$ & $b/b_P$ & 95\% interval & Change vs.\ shipped (\%)", rows, notes)


def t_bootstrap_sd(bs):
    rows = []
    for name, g in bs["generators"].items():
        r = g["boot_sd_over_mc_sd_range_panels"]
        rows.append(rf"{name}{DAG if g['post_hoc_law'] else ''} & {g['panels_bootstrapped']} & {1e3 * g['mc_sd_a_200']:.3f} & {1e3 * g['mean_boot_sd']:.3f} & "
                    rf"{g['boot_sd_over_mc_sd']:.2f} ({r[0]:.2f}--{r[1]:.2f}) & {g['sandwich_se_over_mc_sd']:.2f} & {g['percentile_coverage_of_E_a']:.2f}")
    notes = (r"Pigeonhole (state $\times$ calendar-week) bootstrap, 199 draws per panel, on the first 100 stored panels of each fresh-draw generator "
             r"(panels regenerated exactly from the stored seeds). MC SD: Monte Carlo SD of $\hat a$ over 200 panels. Ratio range over panels in brackets. "
             r"Sandwich: two-way SE over MC SD. Coverage: share of panels whose percentile interval covers $E(\hat a)$. $^\dagger$ post hoc law. "
             rf"Reading rule of Amendment 2: conservative if every ratio exceeds 1 (outcome: {esc(bs['rule']['verdict'])}).")
    table("tab_r2_bootstrap_sd.tex", r"Pigeonhole bootstrap SD of the intercept against the Monte Carlo SD", "tab:r2-bootstrap-sd", r"Xrrrlrr",
          r"Generator & Panels & MC SD ($\times10^3$) & Boot.\ SD ($\times10^3$) & Boot./MC & Sandwich/MC & Coverage", rows, notes)


def t_reporting(rp):
    c, b = rp["coverage_adjustment"], rp["backfill"]
    cp = c["coverage_pct"]
    rows = [r"\multicolumn{3}{@{}l}{\textit{Hospital participation (NHSN, \% of hospitals reporting)}}",
            rf"Mean (min; P10) participation & {cp['mean']:.1f}\% ({cp['min']:.1f}; {cp['p10']:.1f}) & --",
            rf"Rows with full participation in both weeks & {pc(cp['share_rows_100'])}\% & --",
            rf"Mean week-to-week change & {cp['mean_abs_week_to_week_change_pct_points']:.2f} pp & --",
            rf"$a$, observed counts & {c['a']:.4f} & --",
            rf"$a$, participation-adjusted counts & {c['a_adjusted']:.4f} & two-way SE {c['a_adjusted_se_two_way']:.4f}",
            rf"Difference (relative) & {num(c['difference'], 5)} ({pc(c['relative_difference'])}\%) & {ci(c['difference_ci_week_bootstrap'], 5)}",
            r"\midrule", r"\multicolumn{3}{@{}l}{\textit{Backfill (Delphi Epidata issue history, archive from issue 2024w47)}}",
            rf"Rows of $R_A$ in the archive; latest issue equals the vintage & {b['rows_in_archive']:,}; {pc(b['latest_equals_vintage_share_all_archive_rows'])}\% & --",
            rf"Rows with a provisional first issue (lag $\le$ {b['max_lag_first_issue_weeks']} weeks) & {b['rows_covered']} ({b['weeks_covered']} weeks, {b['covered_week_range'][0]} to {b['covered_week_range'][1]}) & --",
            rf"First issue equals latest (current week) & {pc(b['first_equals_latest_share_current_week'])}\% & --",
            rf"Revised up / down & {pc(b['share_revised_up'])}\% / {pc(b['share_revised_down'])}\% & --",
            rf"Mean relative revision (latest $-$ first) / first & {pc(b['mean_relative_revision_current_week'])}\% & --",
            rf"$a$, latest issue ({b['rows_fitted']} rows) & {b['a_latest']:.4f} & {ci(b['a_latest_ci_week_bootstrap'])}",
            rf"$a$, first issue & {b['a_first']:.4f} & {ci(b['a_first_ci_week_bootstrap'])}",
            rf"Difference first $-$ latest & {b['a_first_minus_latest']:.4f} & {ci(b['a_first_minus_latest_ci_week_bootstrap'])}",
            rf"$v_{{\mathrm{{rev}}}}$ (revision variance of $x$) & {b['v_rev']:.4f} & --"]
    notes = (r"Participation adjustment: both counts scaled by 100/participation, $x$ recomputed, the affine fit repeated with the observed $C_t$ "
             r"as the count; assumes unreported hospitals admit as many patients per hospital as reporting ones. Backfill: rows whose two weeks "
             r"both have a provisional first issue; fits restricted to rows with first and latest counts $\ge 5$. Intervals: calendar-week block "
             r"bootstrap, 1,999 draws. The backfill comparison bounds the provisional-versus-settled component only; the residual reporting error "
             r"of the settled vintage is not identified from public data.")
    table("tab_r2_reporting.tex", r"External reporting data: participation and backfill", "tab:r2-reporting", r"Xll",
          r"Quantity & Value & 95\% interval", rows, notes)


def t_interval_types(it):
    rows = []
    for e in it["entries"]:
        lo, hi = e["quoted"]
        neg = lambda s: f"$-${s[1:]}" if s.startswith("-") else s
        q = f"{neg(lo)} to {neg(hi)}" if lo.startswith("-") else f"{lo}--{hi}"
        src = esc(e["source"] or "--").replace(".", r".\allowbreak{}").replace("[", r"\allowbreak{}[")
        rows.append(rf"{esc(e['location'])} & {esc(e['quantity'])} & {q} & {esc(e['type'])} & {{\ttfamily\tiny {src}}}")
    notes = (rf"{it['n']} intervals and ranges quoted in the abstract, Results and Discussion; source found mechanically in the analysis outputs for "
             rf"{it['n_source_found']} (stored pair matching the quoted limits at the quoted precision); the others are ranges assembled from "
             r"several stored values (path pattern given). CGM: Cameron--Gelbach--Miller. All intervals 95\%.")
    table("tab_r2_interval_types.tex", r"Type and source of every interval quoted in the manuscript", "tab:r2-interval-types",
          r"p{0.11\linewidth}p{0.2\linewidth}lp{0.25\linewidth}X", r"Location & Quantity & Quoted & Type & Source", rows, notes, size=r"\tiny")


# ------------------------------------------------------------------ figures
def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(FIG / f"{name}.{ext}", bbox_inches="tight", dpi=300 if ext == "png" else None)
    plt.close(fig)


def f_deciles(fl):
    d = fl["B1_deciles"]
    ft = d["fits"]
    rows = d["rows"]
    mc = np.array([r["mean_C"] for r in rows])
    ml = np.array([r["mean_L"] for r in rows])
    lo = np.array([r["mean_L_ci"][0] for r in rows])
    hi = np.array([r["mean_L_ci"][1] for r in rows])
    C = np.geomspace(5, 3000, 300)
    aff = ft["affine"]["a"] + ft["affine"]["b"] / C
    pw = ft["power"]["c"] * C ** -ft["power"]["gamma"]
    th = ft["three_parameter"]["a"] + ft["three_parameter"]["b"] * C ** -ft["three_parameter"]["gamma"]
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.8))
    ax[0].errorbar(mc, ml, yerr=[ml - lo, hi - ml], fmt="o", color="k", ms=3.5, capsize=2, lw=0.8, label="observed (95% interval)", zorder=5)
    ax[0].plot(C, aff, color="#c0392b", lw=1.2, label=r"affine $a+b/C$")
    ax[0].plot(C, pw, color="#2471a3", lw=1.2, ls="--", label=r"power $cC^{-\gamma}$")
    ax[0].plot(C, th, color="#239b56", lw=1.2, ls=":", label=r"three-parameter $a+bC^{-\gamma}$")
    ax[0].axhline(ft["affine"]["a"], color="#c0392b", lw=0.6, alpha=0.5)
    ax[0].set(xscale="log", yscale="log", xlabel=r"mean current-week count $C_t$ in decile", ylabel=r"mean loss $L$", title="(a) Loss by count decile")
    ax[0].legend(frameon=False, loc="lower left")
    for key, col, mk, lab in (("affine", "#c0392b", "s", "affine"), ("power", "#2471a3", "^", "power"), ("three", "#239b56", "D", "three-parameter")):
        pred = np.array([r[key] for r in rows])
        x = np.arange(1, 11) + {"affine": -0.18, "power": 0, "three": 0.18}[key]
        ax[1].errorbar(x, ml / pred, yerr=[(ml - lo) / pred, (hi - ml) / pred], fmt=mk, color=col, ms=3, capsize=1.5, lw=0.7, label=lab)
    ax[1].axhline(1, color="0.4", lw=0.7)
    ax[1].set(xticks=range(1, 11), xlabel="decile of $C_t$", ylabel="observed / predicted mean loss", title="(b) Observed over predicted")
    ax[1].legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3)
    fig.tight_layout()
    save(fig, "fig_r2_deciles")


def f_population(fl):
    p = fl["B2_population"]
    cv, li, lg = p["covariate"], p["primary_linear"], p["secondary_log_link"]
    pop = np.geomspace(5e5, 4e7, 200)
    z = (np.log(pop) - cv["log_pop_mean"]) / cv["log_pop_sd"]
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.8))
    ax[0].plot(pop, li["a0"] + li["a1"] * z, color="#c0392b", lw=1.2, label=r"linear $a_0+a_1z$ (primary)")
    ax[0].plot(pop, np.exp(lg["alpha0"] + lg["theta"] * z), color="#2471a3", lw=1.2, ls="--", label=r"log-link $\exp(\alpha_0+\theta z)$")
    for t, a_cs, c_cs in zip(p["tertiles_separate"], p["tertiles_common_slope"]["intercepts"], p["tertiles_common_slope"]["intercept_ci_state_bootstrap"]):
        lo, hi = t["pop_range"]
        xm = np.sqrt(lo * hi)
        ax[0].errorbar(xm, t["a"], yerr=[[t["a"] - t["a_ci_two_way"][0]], [t["a_ci_two_way"][1] - t["a"]]], xerr=[[xm - lo], [hi - xm]],
                       fmt="o", color="k", ms=4, capsize=2, lw=0.8, label="tertile, separate fit" if t["tertile"] == 1 else None)
        ax[0].errorbar(xm * 1.12, a_cs, yerr=[[a_cs - c_cs[0]], [c_cs[1] - a_cs]], fmt="o", mfc="white", color="0.35", ms=4, capsize=2, lw=0.8,
                       label="tertile, common slope" if t["tertile"] == 1 else None)
    for q in ("pop_p10", "pop_p90"):
        ax[0].axvline(cv[q], color="0.7", lw=0.6, ls=":")
    ax[0].set(xscale="log", xlabel="state population", ylabel=r"intercept $a$", title="(a) Intercept against population", ylim=(0, None))
    ax[0].legend(frameon=False, loc="upper right", fontsize=6.5)
    labs = ["linear, P90/P10", "log-link, P90/P10", "tertile 3/1, common slope"]
    vals = [li["ratio_a_p90_over_p10"], lg["ratio_a_p90_over_p10"], p["tertiles_common_slope"]["ratio_high_over_low"]]
    dz = cv["z_p90"] - cv["z_p10"]  # ratio exp(theta dz) is monotone in theta: its percentile interval maps from that of theta
    cis = [li["ratio_ci_state_bootstrap"], sorted(float(np.exp(t * dz)) for t in lg["theta_ci_state_bootstrap"]),
           p["tertiles_common_slope"]["ratio_ci_state_bootstrap"]]
    for i, (v, c) in enumerate(zip(vals, cis)):
        ax[1].errorbar(v, i, xerr=None if c is None else [[v - c[0]], [c[1] - v]], fmt="o", color="k", ms=4, capsize=2, lw=0.8)
    ax[1].axvline(1, color="0.4", lw=0.7)
    ax[1].set(yticks=range(3), yticklabels=labs, xlim=(0, 1.2), ylim=(-0.6, 2.6), xlabel="ratio of intercepts (95% state bootstrap)",
              title="(b) Large-state over small-state intercept")
    fig.tight_layout()
    save(fig, "fig_r2_population")


def f_floor(fl, ri):
    b3, b4 = fl["B3_three_parameter"], fl["B4_model_free_floor"]
    e1 = ri["E1_primary"]
    t3 = fl["B2_population"]["tertiles_separate"][2]
    items = [("affine $a$ (two-way)", e1["a"], e1["inference"]["two_way"]["a_ci"], "#c0392b"),
             ("affine $a$, largest-state tertile (two-way)", t3["a"], t3["a_ci_two_way"], "#c0392b"),
             (r"three-parameter $a$ (week bootstrap)", b3["point"]["a"], b3["week_bootstrap"]["a"]["ci"], "#239b56"),
             (r"three-parameter $a$ (state bootstrap)", b3["point"]["a"], b3["state_bootstrap"]["a"]["ci"], "#239b56")]
    for r in b4["rows"]:
        q = int(round(100 * r["threshold_quantile"]))
        items.append((rf"model-free, $C_t\geq${r['threshold_C']:.0f} (P{q}): $L-\hat b/C_t$", r["floor_minus_bhat"], r["floor_minus_bhat_ci"], "#2471a3"))
    fig, ax = plt.subplots(figsize=(5.6, 2.6))
    for i, (lab, v, c, col) in enumerate(items[::-1]):
        ax.errorbar(v, i, xerr=[[v - c[0]], [c[1] - v]], fmt="o", color=col, ms=4, capsize=2, lw=0.9)
    ax.axvline(0, color="0.4", lw=0.7)
    ax.set(yticks=range(len(items)), yticklabels=[it[0] for it in items[::-1]], xlabel="floor estimate (relative variance, 95% interval)")
    # top axis: relative SD = 100*sqrt(floor); ticks placed at the exact squared positions (a secondary_xaxis with a
    # clipped sqrt is not invertible below 0 and misplaced the 5% tick)
    sds = [5, 8, 10, 12, 14]
    sec = ax.twiny()
    sec.set_xlim(ax.get_xlim())
    sec.set_xticks([(s / 100) ** 2 for s in sds])
    sec.set_xticklabels([str(s) for s in sds])
    sec.set_xlabel("relative SD (%)")
    fig.tight_layout()
    save(fig, "fig_r2_floor")


def main():
    L2.register_script(__file__)
    log = Log(R.LOGS / "rev_r2_render.log")
    FIG.mkdir(exist_ok=True)
    js = {n: load(n) for _, n in BLOCKS}
    ri = load("rev_intercept.json")
    b_ci = ri["E1_primary"]["inference"]["two_way"]["b_ci"]
    merge(js, b_ci)
    fl = js["rev_r2_floor.json"]
    t_deciles(fl)
    t_population(fl)
    t_floor(fl, ri)
    t_bias(js["rev_r2_total_bias.json"])
    t_generator(js["rev_r2_generator.json"])
    t_clamp(js["rev_r2_clamp.json"])
    t_poisson(js["rev_r2_poisson_ref.json"], b_ci)
    t_bootstrap_sd(js["rev_r2_bootstrap_sd.json"])
    t_reporting(js["rev_r2_reporting.json"])
    t_interval_types(js["rev_r2_interval_types.json"])
    f_deciles(fl)
    f_population(fl)
    f_floor(fl, ri)
    log("rev_r2_render done", sorted(p.name for p in R.TAB.glob("tab_r2_*.tex")), sorted(p.name for p in FIG.glob("fig_r2_*")))


if __name__ == "__main__":
    main()
