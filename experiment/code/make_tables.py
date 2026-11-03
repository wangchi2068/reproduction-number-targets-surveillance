"""LaTeX tables (booktabs) in English and Chinese from analysis/results.json and the part blocks.
Every number is read from results.json / parts; nothing is recomputed. Writes tables/{en,zh}/*.tex and
tables/manifest.json (table -> results.json keys used).
Usage: python3 make_tables.py
"""
import json

from cc_common import ANALYSIS, EXP, PARTS

R = json.loads((ANALYSIS / "results.json").read_text())["numbers"]
BL = {p.stem: json.loads(p.read_text())["blocks"] for p in PARTS.glob("s*.json")}
OUT = EXP / "tables"
manifest = {}
SPEC_LAB = {"RA_rule4": ("$R_A$, rule windows ($S=4$)", "$R_A$,规则窗口($S=4$)"),
            "G_rule4": ("$G$, rule windows ($S=4$)", "$G$,规则窗口($S=4$)"),
            "G_cal3": ("$G$, calendar windows ($S=3$)", "$G$,日历窗口($S=3$)"),
            "RA_cal3": ("$R_A$, calendar windows ($S=3$)", "$R_A$,日历窗口($S=3$)")}


def g(key, used):
    used.append(key)
    return R[key]


def fmt(v, d=3):
    if v is None:
        return "--"
    if isinstance(v, bool):
        return v
    if isinstance(v, int):
        return f"{v:,}".replace(",", "{,}")
    if v < 0:
        return "$-$" + fmt(-v, d)
    a = abs(v)
    if a != 0 and (a < 10 ** -d * 10):
        return f"{v:.2g}" if a >= 1e-3 else f"{v:.1e}"
    return f"{v:.{d}f}"


def ci(rec, d=3, s=1.0):
    c = rec.get("ci95")
    v = rec["value"] * s
    return f"{fmt(v, d)} [{fmt(c[0] * s, d)}, {fmt(c[1] * s, d)}]" if c else fmt(v, d)


def yn(b, lang):
    return ("yes" if b else "no") if lang == "en" else ("是" if b else "否")


def write(name, lang, caption, label, colspec, header, rows, notes):
    lines = [r"\begin{table}[htbp]", r"\centering", r"\footnotesize", rf"\caption{{{caption}}}", rf"\label{{{label}}}",
             rf"\begin{{tabular}}{{{colspec}}}", r"\toprule", " & ".join(header) + r" \\", r"\midrule"]
    lines += [" & ".join(str(c) for c in r) + r" \\" for r in rows]
    lines += [r"\bottomrule", r"\end{tabular}"]
    if notes:
        lines += [r"\par\smallskip", r"\begin{minipage}{0.95\linewidth}\scriptsize " + notes + r"\end{minipage}"]
    lines += [r"\end{table}", ""]
    (OUT / lang).mkdir(parents=True, exist_ok=True)
    (OUT / lang / f"{name}.tex").write_text("\n".join(lines), encoding="utf-8")


def T(lang, en, zh):
    return zh if lang == "zh" else en


# ----------------------------------------------------------------------------- 1 headline estimates
for lang in ("en", "zh"):
    used = []
    rows = []
    for s in ("RA_rule4", "G_rule4", "G_cal3", "RA_cal3"):
        p = "s02_inference."
        rows.append([SPEC_LAB[s][lang == "zh"], fmt(g(p + s + "_N", used)["value"]), fmt(g(p + s + "_n_weeks", used)["value"]),
                     fmt(g(p + s + "_xbar", used)["value"]), ci(g(p + s + "_a", used), 4), ci(g(p + s + "_b", used), 3),
                     ci(g(p + s + "_ba", used), 1), fmt(g(p + s + "_a_sd_ratio_twoway_vs_state", used)["value"], 2)])
    write("tab_s1_headline", lang,
          T(lang, "Intercept $a$, slope $b$ and balance ratio $b/a$ of the relative loss by specification.",
            "各设定下相对损失的截距 $a$、斜率 $b$ 与平衡比 $b/a$。"),
          "tab:s1", "lrrrlllr",
          [T(lang, "Specification", "设定"), "$N$", T(lang, "Weeks", "周数"), r"$\bar{x}$", "$a$ [95\\% CI]", "$b$ [95\\% CI]",
           "$b/a$ [95\\% CI]", T(lang, "SD ratio", "标准差比")],
          rows,
          T(lang, "Iterated FGLS fit of $L=((x-\\bar{x})/\\bar{x})^2$ on $a+b/D$; 95\\% percentile intervals from the two-way "
                  "(state $\\times$ week) pigeonhole bootstrap, $B=2000$, $\\bar{x}$ re-estimated in every replicate. SD ratio: "
                  "bootstrap SD of $a$ under two-way versus state-only resampling. 51 jurisdictions throughout.",
            "以迭代 FGLS 将 $L=((x-\\bar{x})/\\bar{x})^2$ 拟合为 $a+b/D$;区间为州 $\\times$ 周双向鸽笼自助法 95\\% 百分位区间,"
            "$B=2000$,每次重抽重新估计 $\\bar{x}$。标准差比:双向重抽与仅按州重抽下 $a$ 的自助标准差之比。均为 51 个辖区。"))
    manifest["tab_s1_headline"] = sorted(set(used))

# ----------------------------------------------------------------------------- 2 decomposition and plateau
for lang in ("en", "zh"):
    used = []
    rows = []
    for s in ("RA_rule4", "G_rule4", "G_cal3"):
        p = "s03_decomposition." + s
        rows.append([SPEC_LAB[s][lang == "zh"], ci(g(p + "_a_total", used), 4), ci(g(p + "_a_weekFE", used), 4),
                     ci(g(p + "_a_both", used), 4), fmt(g(p + "_share_week", used)["value"], 2),
                     fmt(g(p + "_icc_week", used)["value"], 2), fmt(g(p + "_power_exponent", used)["value"], 2),
                     fmt(g(p + "_plateau_contrast_p", used)["value"], 2),
                     f"{g(p + '_cv_state_folds_affine_better', used)['value']}/51"])
    write("tab_s2_decomposition", lang,
          T(lang, "Decomposition of the intercept and the plateau versus power-law comparison.",
            "截距分解及平台与幂律形式的比较。"),
          "tab:s2", "lllrrrrrr",
          [T(lang, "Specification", "设定"), T(lang, "$a$ total", "$a$ 总量"), T(lang, "$a$ | week FE", "$a$ | 周固定效应"),
           T(lang, "$a$ | week + state-season", "$a$ | 周 + 州季"), T(lang, "Week share", "周占比"), "ICC",
           T(lang, "Power exp.", "幂指数"), "$p$", T(lang, "CV folds", "交叉验证折")],
          rows,
          T(lang, "Intervals: two-way bootstrap with $B=1000$ (separate draws from Table S1, hence small differences). Week FE: national week fixed effects removed from the statistic before the fit; state-season effects removed "
                  "in addition in the next column. Week share = 1 - $a$|week FE / $a$ total. ICC: share of the statistic's variance "
                  "between weeks. Power exp.: exponent $c$ of $L = \\alpha D^{-c}$. $p$: parametric-bootstrap test of the affine "
                  "plateau against the power law (200 null draws). CV folds: leave-one-state-out folds in which the affine model "
                  "predicts better. A power analysis found no fraction of the estimated $a$ ($\\le 1$) detectable with 80\\% power.",
            "区间为双向自助法 $B=1000$(与表 S1 抽样独立,故略有差异)。周固定效应:拟合前从统计量中去除全国周效应;下一列再去除州-季效应。周占比 = 1 - ($a$|周固定效应)/($a$ 总量)。"
            "ICC:统计量方差中周间部分的占比。幂指数:$L=\\alpha D^{-c}$ 的 $c$。$p$:仿射平台对幂律的参数自助检验(200 次零分布抽样)。"
            "交叉验证折:逐州留一时仿射模型预测更优的折数。功效分析显示,真实截距为估计值的任意比例($\\le 1$)时均达不到 80\\% 功效。"))
    manifest["tab_s2_decomposition"] = sorted(set(used))

# ----------------------------------------------------------------------------- 3 reporting completeness
for lang in ("en", "zh"):
    used = []
    rows = []
    for s in ("RA_rule4", "G_rule4", "G_cal3"):
        p = "s04_reporting." + s
        rec = g(p + "_delta_a", used)
        c90 = rec["ci90"]
        nS = len(BL["s04_reporting"][f"delta|{s}"]["by"]) - 1
        rows.append([SPEC_LAB[s][lang == "zh"], f"{fmt(rec['value'] * 1e3, 2)} [{fmt(c90[0] * 1e3, 2)}, {fmt(c90[1] * 1e3, 2)}]",
                     fmt(g(p + "_tost_margin_tau0_0.1", used)["value"] * 1e3, 2), yn(g(p + "_tost_equivalent_tau0_0.1", used)["value"], lang),
                     fmt(g(p + "_tost_margin_tau0_0.05", used)["value"] * 1e3, 2), yn(g(p + "_tost_equivalent_tau0_0.05", used)["value"], lang),
                     f"{g(p + '_delta_a_seasons_equivalent_tau0_0.1', used)['value']}/{nS}"])
    write("tab_s3_reporting", lang,
          T(lang, "Contribution of facility-reporting completeness to the intercept: equivalence tests.",
            "医院报告完整性对截距的贡献:等效性检验。"),
          "tab:s3", "llrlrlr",
          [T(lang, "Specification", "设定"), "$10^3\\Delta a$ [90\\% CI]", "$10^3\\delta_{0.10}$", T(lang, "Equiv.", "等效"),
           "$10^3\\delta_{0.05}$", T(lang, "Equiv.", "等效"), T(lang, "Seasons equiv.", "等效季数")],
          rows,
          T(lang, "$\\Delta a$: intercept of the statistic minus that of the statistic recomputed from counts divided by the share "
                  "of hospitals reporting; paired two-way bootstrap, $B=2000$. Margins $\\delta_{\\tau_0}=\\kappa^2\\tau_0^2(1+v_R/\\bar{R}^2)$ "
                  "fixed before fitting ($\\kappa=1$ for $G$), $v_R=0.04$. Equivalence: 90\\% interval inside $(-\\delta,\\delta)$ (TOST, "
                  "$\\alpha=0.05$). Seasons equiv.: seasons equivalent at $\\tau_0=0.10$.",
            "$\\Delta a$:统计量截距减去以报告医院占比校正计数后的统计量截距;配对双向自助法,$B=2000$。等效界 "
            "$\\delta_{\\tau_0}=\\kappa^2\\tau_0^2(1+v_R/\\bar{R}^2)$ 在拟合前确定($G$ 取 $\\kappa=1$),$v_R=0.04$。"
            "等效:90\\% 区间落在 $(-\\delta,\\delta)$ 内(TOST,$\\alpha=0.05$)。等效季数:$\\tau_0=0.10$ 时等效的流行季数。"))
    manifest["tab_s3_reporting"] = sorted(set(used))

# ----------------------------------------------------------------------------- 4 calendar-matched and laboratory proxy
for lang in ("en", "zh"):
    used = []
    rows = []
    for y, reg in (("2023", T(lang, "mandatory", "强制")), ("2024", T(lang, "voluntary", "自愿")), ("2025", T(lang, "mandatory", "强制"))):
        p = f"s04_reporting.mayoct{y}"
        rows.append([f"{y}", reg, fmt(g(p + "_N", used)["value"]), ci(g(p + "_sd_G", used), 2), ci(g(p + "_iqr_logG", used), 2)])
    lab = g("s04_reporting.lab_share_of_a_shared", used)
    cor = g("s04_reporting.lab_corr_RA", used)
    nlab = g("s04_reporting.lab_N", used)["value"]
    write("tab_s4_calendar_lab", lang,
          T(lang, "Calendar-matched comparison of the voluntary-reporting period (May--October).",
            "自愿报告期与同历月(5--10 月)的对照。"),
          "tab:s4", "llrll",
          [T(lang, "Year", "年份"), T(lang, "Reporting", "报告制度"), "$N$", T(lang, "SD of $G$", "$G$ 标准差"),
           T(lang, "IQR of $\\log G$", "$\\log G$ 四分位距")],
          rows,
          T(lang, f"State-weeks with $D_t\\ge 5$ whose week ends between 1 May and 31 October; 95\\% two-way bootstrap intervals, "
                  f"$B=1000$. The FGLS intercept is not stable in the voluntary period (negative point estimate), so robust spread "
                  f"measures are reported. Laboratory proxy (Delphi fluview\\_clinical, {nlab} state-weeks): share of the $R_A$ "
                  f"intercept shared with laboratory growth {ci(lab, 2)}; correlation {ci(cor, 2)}.",
            f"周末日期在 5 月 1 日至 10 月 31 日且 $D_t\\ge 5$ 的州-周;95\\% 双向自助区间,$B=1000$。自愿期 FGLS 截距不稳定"
            f"(点估计为负),故报告稳健离散度。实验室代理(Delphi fluview\\_clinical,{nlab} 个州-周):$R_A$ 截距中与实验室增长"
            f"共有的份额 {ci(lab, 2)};相关系数 {ci(cor, 2)}。"))
    manifest["tab_s4_calendar_lab"] = sorted(set(used))

# ----------------------------------------------------------------------------- 5 seasons
for lang in ("en", "zh"):
    used = []
    specs = ("RA_rule4", "G_rule4", "G_cal3")
    rows = []
    for s_ in ("2022-23", "2023-24", "2024-25", "2025-26"):
        row = [s_]
        for s in specs:
            k = f"s05_seasons.{s}_a_{s_}"
            row.append(ci(g(k, used), 4) if k in R else "--")
        rows.append(row)
    rows.append([T(lang, "$a_{\\mathrm{RE}}$ [95\\% CI]", "$a_{\\mathrm{RE}}$ [95\\% CI]")] + [ci(g(f"s05_seasons.{s}_a_RE", used), 4) for s in specs])
    rows.append(["$Q$ ($p$)"] + [f"{fmt(g(f's05_seasons.{s}_Q', used)['value'], 2)} ({fmt(R[f's05_seasons.{s}_Q']['p'], 3)})" for s in specs])
    rows.append(["$10^3\\tau^2$ [Q-profile]"] + [ci(g(f"s05_seasons.{s}_tau2", used), 3, 1e3) for s in specs])
    rows.append(["$I^2$"] + [fmt(g(f"s05_seasons.{s}_I2", used)["value"], 2) for s in specs])
    write("tab_s5_seasons", lang,
          T(lang, "Season-specific intercepts and between-season heterogeneity, reported separately for $S=4$ and $S=3$.",
            "分季截距与季间异质性($S=4$ 与 $S=3$ 分别报告)。"),
          "tab:s5", "llll",
          [T(lang, "Season", "流行季")] + [SPEC_LAB[s][lang == "zh"] for s in specs],
          rows,
          T(lang, "Season intercepts use the season's own mean; 95\\% two-way bootstrap intervals, $B=1000$. DerSimonian--Laird "
                  "random effects with a $t_{S-1}$ interval; Q-profile 95\\% interval for $\\tau^2$. The relative standard error of a "
                  "between-season variance estimated from $S$ seasons is at least $\\sqrt{2/(S-1)}$: 1.00 for $S=3$ and 0.82 for $S=4$.",
            "分季截距以该季自身均值计算;95\\% 双向自助区间,$B=1000$。随机效应采用 DerSimonian--Laird 估计及 $t_{S-1}$ 区间;"
            "$\\tau^2$ 给出 Q-profile 95\\% 区间。由 $S$ 季估计的季间方差,其相对标准误至少为 $\\sqrt{2/(S-1)}$:$S=3$ 时 1.00,$S=4$ 时 0.82。"))
    manifest["tab_s5_seasons"] = sorted(set(used))

# ----------------------------------------------------------------------------- 6 pathogens
mult = {(t["combo"], t["criterion"]): t for t in BL["s06_pathogens"]["multiplicity"]["tests"]}
PG = {"flu": ("Influenza", "流感"), "covid": ("COVID-19", "新冠"), "rsv": ("RSV", "合胞病毒")}
WN = {"W-rule4": ("rule", "规则"), "W-cal": ("calendar", "日历"), "W-growth": ("growth", "增长期"),
      "W-growth-2025": ("growth, 2025-26", "增长期(2025-26 季)")}
for lang in ("en", "zh"):
    used = []
    rows = []
    for c in ("flu|W-rule4", "flu|W-cal", "flu|W-growth", "flu|W-growth-2025", "covid|W-cal", "covid|W-growth", "rsv|W-cal", "rsv|W-growth"):
        b = BL["s06_pathogens"][f"combo|{c}"]
        pg, w = c.split("|")
        nb = b["null_band"]
        ph = lambda k: fmt(mult[(c, k)]["p_holm"], 3)  # noqa: E731
        rows.append([f"{PG[pg][lang == 'zh']}, {WN[w][lang == 'zh']}" if lang == "en" else f"{PG[pg][1]},{WN[w][1]}",
                     fmt(b["N"]["original"]), ci(g(f"s06_pathogens.{c}_q1", used), 2), ci(g(f"s06_pathogens.{c}_q2", used), 2),
                     f"{fmt(nb['q1'][0], 2)}--{fmt(nb['q1'][1], 2)}", f"{fmt(nb['q2'][0], 2)}--{fmt(nb['q2'][1], 2)}",
                     ph("q1_lt_1"), ph("q2_gt_1"), ph("q1_outside_null"), ph("q2_outside_null")])
    used += ["s06_pathogens.pathogen_family_size", "s06_pathogens.pathogen_n_sig_holm_0.05", "s06_pathogens.pathogen_n_sig_bh_0.05"]
    write("tab_s6_pathogens", lang,
          T(lang, "Cross-pathogen intercept ratios against the pure-reporting growth null, with Holm-adjusted $p$-values.",
            "跨病原截距比与纯报告增长零假设的比较(Holm 校正 $p$ 值)。"),
          "tab:s6", "lrllllrrrr",
          [T(lang, "Pathogen, window", "病原,窗口"), "$N$", "$q_1$ [95\\% CI]", "$q_2$ [95\\% CI]", T(lang, "Null $q_1$", "零假设 $q_1$"),
           T(lang, "Null $q_2$", "零假设 $q_2$"), "$q_1<1$", "$q_2>1$", T(lang, "$q_1$ out", "$q_1$ 越界"),
           T(lang, "$q_2$ out", "$q_2$ 越界")],
          rows,
          T(lang, "$q_1=a_{\\text{1 week}}/a_{\\text{original}}$, $q_2=a_{\\text{equal 3 weeks}}/a_{\\text{original}}$; shared two-way "
                  "bootstrap, $B=2000$. Null ranges: growth-null simulation without transmission heterogeneity (three generation "
                  "intervals, three admission delays, four reporting-noise levels; the rule-window row uses the influenza calendar "
                  "range). The last four columns give one-sided bootstrap $p$-values, Holm-adjusted over the 32 ratio tests; \"out\" "
                  "tests whether the ratio lies outside the null range. Negative lower confidence bounds (such as for COVID-19 calendar "
                  "$q_1$ and RSV growth $q_2$) arise because in a small fraction of bootstrap resamples the reference denominator intercept "
                  "is estimated below zero.",
            "$q_1=a_{\\text{单周}}/a_{\\text{原权重}}$, $q_2=a_{\\text{三周等权}}/a_{\\text{原权重}}$; 共享双向自助法, $B=2000$。零假设范围来自"
            "无传播异质性的增长期零模拟（三种世代间隔、三种入院延迟、四档报告噪声；规则窗口行沿用流感日历窗口范围）。后四列为单侧自助 "
            "$p$ 值，对全部 32 个比值检验做 Holm 校正；越界检验比值是否落在零假设范围之外。部分下置信界为负值（如新冠日历窗口 $q_1$ "
            "及合胞病毒增长期 $q_2$），系因少数自助抽样中基准分母截距出现负值所致。"))
    manifest["tab_s6_pathogens"] = sorted(set(used))

# ----------------------------------------------------------------------------- 7 theory reproduction
THEORY = [("$m_\\times$", "53.1", "s07_theory.m_cross", 2), ("$v_H$ ($\\tau=0.15$)", "0.0307", "s07_theory.v_H_tau015", 4),
          ("$m_{\\times,H}$ ($\\tau=0.15$)", "30.4", "s07_theory.m_cross_H_tau015", 2),
          ("$m_0^*$ ($c_0/c_1=100$)", "72.9", "s07_theory.m0_star_baseline_c100", 2),
          ("$T^\\dagger$ ($\\phi=0.8$)", "4.48", "s07_theory.T_dagger_phi08", 3),
          ("$B^\\dagger$", "775", "s07_theory.B_dagger_baseline", 1),
          ("$B^\\dagger$ ($m_\\times=21.3$)", "655", "s07_theory.B_dagger_mx21.3", 1),
          ("$B^\\dagger$ ($m_\\times=49.2$)", "763", "s07_theory.B_dagger_mx49.2", 1),
          (("Reporting lever term ($\\times10^{-3}$)", "报告杠杆项($\\times10^{-3}$)"), "7.85", "s07_theory.lever_base_rep", 2),
          (("Lever total ($\\times10^{-3}$)", "杠杆合计($\\times10^{-3}$)"), "22.5", "s07_theory.lever_base_total", 2),
          (("Fisher constant", "Fisher 常数"), "0.505", "s07_theory.fisher_const_0505", 3),
          (("Relative SE, $S=3$", "相对标准误,$S=3$"), "0.82", "s07_theory.rse_season_S3", 3),
          (("Relative SE, $S=4$", "相对标准误,$S=4$"), "0.71", "s07_theory.rse_season_S4", 3)]
for lang in ("en", "zh"):
    used = []
    rows = []
    for lab, ms, key, d in THEORY:
        v = g(key, used)["value"]
        ok = abs(float(ms) - v) <= 0.51 * 10 ** -(len(ms.split(".")[1]) if "." in ms else 0)
        lab = lab[lang == "zh"] if isinstance(lab, tuple) else lab
        rows.append([lab, ms, fmt(v, d), T(lang, "reproduced", "复现") if ok else T(lang, "corrected", "更正")])
    write("tab_s7_theory", lang,
          T(lang, "Theory quantities of the submitted manuscript recomputed from the closed forms.",
            "以闭式解重算投稿稿中的理论数值。"),
          "tab:s7", "lrrl",
          [T(lang, "Quantity", "数量"), T(lang, "Manuscript", "原稿"), T(lang, "Recomputed", "重算"), T(lang, "Status", "状态")],
          rows,
          T(lang, "Baseline $\\bar{R}=1.15$, $v_R=0.04$, $k=0.35$, $\\rho=0.25$. The manuscript's relative standard errors for "
                  "$S=3,4$ equal $\\sqrt{2/S}$; the stated bound $\\sqrt{2/(S-1)}$ gives the recomputed values. The lever terms differ "
                  "by rounding.",
            "基线 $\\bar{R}=1.15$,$v_R=0.04$,$k=0.35$,$\\rho=0.25$。原稿 $S=3,4$ 的相对标准误等于 $\\sqrt{2/S}$;按文中给出的界 "
            "$\\sqrt{2/(S-1)}$ 得到重算值。杠杆项差异来自舍入。"))
    manifest["tab_s7_theory"] = sorted(set(used))

TABLE_MAP = {
    "tab_s1_headline": ("tab:s1", "tab_s1_headline.tex"),
    "tab_s2_decomposition": ("tab:s2", "tab_s2_decomposition.tex"),
    "tab_s3_reporting": ("tab:s3", "tab_s3_reporting.tex"),
    "tab_s4_calendar_lab": ("tab:s4", "tab_s4_calendar_lab.tex"),
    "tab_s5_seasons": ("tab:s5", "tab_s5_seasons.tex"),
    "tab_s6_pathogens": ("tab:s6", "tab_s6_pathogens.tex"),
    "tab_s7_theory": ("tab:s7", "tab_s7_theory.tex"),
}
table_list = []
for name, used_keys in manifest.items():
    tid, tpath = TABLE_MAP.get(name, (name, f"{name}.tex"))
    table_list.append({
        "id": tid,
        "path": tpath,
        "sha256": "",
        "generator": "experiment/code/make_tables.py",
        "data_sources": ["experiment/analysis/results.json"],
        "display": "raw",
        "n_seeds": 1
    })

(OUT / "manifest.json").write_text(json.dumps({
    "schema_version": 1,
    "source": "analysis/results.json",
    "script": "code/make_tables.py",
    "tables": table_list
}, indent=1, ensure_ascii=False) + "\n")
print("tables:", list(manifest))
