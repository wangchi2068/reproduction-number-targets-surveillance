"""Round-3 paper tables: S8 (inversion settings) and S9 (analysis history), plus the round-3 edits of
S1, S2, S3 and S6 (no thousands separators, season CV folds, margin wording, undefined ratio bounds).
Every number is read from analysis/results.json or analysis/parts/inversion.json; nothing is recomputed.
Usage: python3 make_tables_r3.py <paper_run_dir> [--lang zh]   (edits <paper_run_dir>/tables in place, updates manifest)
"""
import hashlib
import json
import re
import sys
from pathlib import Path

EXP = Path(__file__).resolve().parents[1]
ANALYSIS = EXP / "analysis"
if not (ANALYSIS / "results.json").exists():  # published copy lives in output/<slug>/experiment
    ANALYSIS = EXP.parents[1] / "output" / EXP.parent.name / "experiment" / "analysis"
R = json.loads((ANALYSIS / "results.json").read_text())["numbers"]
INV = json.loads((ANALYSIS / "parts" / "inversion.json").read_text())["blocks"]["sensitivity"]
TAB = Path(sys.argv[1])
if (TAB / "tables").exists():
    TAB = TAB / "tables"
ZH = "--lang" in sys.argv and sys.argv[sys.argv.index("--lang") + 1] == "zh"
GEN = "experiment/code/make_tables_r3.py"


def f(v, d=3):
    return ("$-$" if v < 0 else "") + f"{abs(v):.{d}f}"


def ci(lo, hi, d=3):
    return f"[{f(lo, d)}, {f(hi, d)}]"


def num(key, d):
    r = R[key]
    return f"{f(r['value'], d)} {ci(*r['ci95'], d)}"


# --- S8: implied environmental variance under each inversion setting ---------------------------------
LAB = {"M0": "Data-matched design (M0)", "rho0": "No national shock", "rho65": "National share 0.65",
       "ar05": "AR(1), $\\phi=0.5$", "gi_short": "True GI mean 2.2 d", "gi_long": "True GI mean 4.3 d",
       "delay_short": "Delay mean 3.5 d", "delay_long": "Delay mean 8 d", "k01": "$k=0.1$",
       "count_shock": "Synchronous count shock", "bridge_s08": "Generic design (bridge)"}
LAB_ZH = {"M0": "数据匹配设计（M0）", "rho0": "无全国冲击", "rho65": "全国份额 0.65",
          "ar05": "AR(1)，$\\phi=0.5$", "gi_short": "真实代际间隔均值 2.2 天", "gi_long": "真实代际间隔均值 4.3 天",
          "delay_short": "延迟均值 3.5 天", "delay_long": "延迟均值 8 天", "k01": "$k=0.1$",
          "count_shock": "报告计数同步冲击", "bridge_s08": "通用设计（衔接）"}
if ZH:
    LAB = LAB_ZH
rows = []
for key, lab in LAB.items():
    s = INV[key]
    ra, g = s["RA"], s["G"]
    rat = s["ratio_RA_over_G"]
    agree = ("是" if s["agree"] else "否") if ZH else ("yes" if s["agree"] else "no")
    if key == "bridge_s08":
        agree = "--"
    rows.append(f"{lab} & {f(ra['v_implied'])} {ci(*ra['v_ci95'])} & {f(g['v_implied'])} {ci(*g['v_ci95'])} & "
                f"{f(ra['kappa_at_nearest'])} & {f(g['kappa_at_nearest'])} & "
                f"{f(rat['value'], 2)} {ci(*rat['ci95'], 2)} & {agree} \\\\")
s8 = r"""\begin{table}[htbp]
\centering
\footnotesize
\setlength{\tabcolsep}{4pt}
\caption{Implied environmental variance $v_R^{\mathrm{impl}}$ from the two-step inversion under each simulation setting, with the 95\% confidence interval of each estimate.}
\label{tab:s8}
\begin{tabular}{llllll}
\toprule
Setting & $v_R^{\mathrm{impl}}$ from $R_A$ & $v_R^{\mathrm{impl}}$ from $G$ & $\kappa_{R_A}$ & $\kappa_G$ & Ratio $R_A/G$ & Agree \\
\midrule
""".replace("{llllll}", "{lllrrll}") + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\par\smallskip
\begin{minipage}{0.95\linewidth}\scriptsize Each setting changes one element of the data-matched design (Section~\ref{sec:meth-inversion}); the analyst's generation interval stays at a mean of 3.6 days throughout. Intervals combine 2000 joint two-way bootstrap draws of the two empirical intercepts with draws of the calibration line from its Monte Carlo sampling distribution. $\kappa=\E(\hat{a})/(v_R/\bar{R}^2)$ at the grid value of $v_R$ nearest to the implied variance. Agree marks settings in which the 95\% interval of the difference between the two implied variances contains zero; the generic design of the first simulation is shown as a bridge and is not counted among the ten settings.\end{minipage}
\end{table}
"""
if ZH:
    s8 = s8.replace(r"\caption{Implied environmental variance $v_R^{\mathrm{impl}}$ from the two-step inversion under each simulation setting, with the 95\% confidence interval of each estimate.}",
                    r"\caption{两步反演在各模拟设定下给出的隐含环境方差 $v_R^{\mathrm{impl}}$，以及各估计值的 95\% 置信区间。}")
    s8 = s8.replace(r"Setting & $v_R^{\mathrm{impl}}$ from $R_A$ & $v_R^{\mathrm{impl}}$ from $G$ & $\kappa_{R_A}$ & $\kappa_G$ & Ratio $R_A/G$ & Agree \\",
                    r"设定 & 由 $R_A$ 得到的 $v_R^{\mathrm{impl}}$ & 由 $G$ 得到的 $v_R^{\mathrm{impl}}$ & $\kappa_{R_A}$ & $\kappa_G$ & 比值 $R_A/G$ & 一致 \\")
    s8 = re.sub(r"\\scriptsize .*?\\end\{minipage\}", lambda m: r"\scriptsize 各设定只改变数据匹配设计中的一个要素（第~\ref{sec:meth-inversion}~节），分析者所用代际间隔的均值始终为 3.6 天。区间合并了两个实证截距的 2000 次联合双向自助抽样，以及校准线按其蒙特卡罗抽样分布的抽取。$\kappa=\E(\hat{a})/(v_R/\bar{R}^2)$，取最接近隐含方差的 $v_R$ 格点值。“一致”表示两个隐含方差之差的 95\% 区间包含零；第一套模拟的通用设计作为衔接设定列出，不计入十种设定。\end{minipage}", s8, flags=re.S)
(TAB / "tab_s8_inversion_settings.tex").write_text(s8)

# --- S9: analysis history ------------------------------------------------------------------------------
P = "s02_inference."
s9 = r"""\begin{table}[htbp]
\centering
\footnotesize
\setlength{\tabcolsep}{4pt}
\caption{Specifications used across the rounds of this analysis, with the 95\% confidence interval of each estimate on the data version of 2026-10-20.}
\label{tab:s9}
\begin{tabular}{p{0.24\linewidth}p{0.16\linewidth}lll}
\toprule
Round and specification & Role & $a$ [95\% CI] & $b/a$ [95\% CI] \\
\midrule
""".replace("{p{0.24\\linewidth}p{0.16\\linewidth}lll}", "{p{0.27\\linewidth}p{0.22\\linewidth}ll}") + "\n".join([
    f"Round 1: $G$, three calendar windows ($S=3$) & Primary in round 1; secondary from round 2 & {num(P + 'G_cal3_a', 4)} & {num(P + 'G_cal3_ba', 1)} \\\\",
    r"Round 2, rule version 1: windows ending at the season maximum & Abandoned (four weeks in 2024--25); no estimates retained & -- & -- \\",
    f"Round 2, rule version 2 (2026-10-20): $G$, rule windows ($S=4$) & Secondary & {num(P + 'G_rule4_a', 4)} & {num(P + 'G_rule4_ba', 1)} \\\\",
    f"Round 2, rule version 2: $R_A$, rule windows ($S=4$) & Primary from round 2 & {num(P + 'RA_rule4_a', 4)} & {num(P + 'RA_rule4_ba', 1)} \\\\",
    f"Round 2: $R_A$, calendar windows ($S=3$) & Sensitivity & {num(P + 'RA_cal3_a', 4)} & {num(P + 'RA_cal3_ba', 1)} \\\\",
    r"Round 3: two-step inversion & Calibration of the round-2 fits; no refit & -- & -- \\",
]) + r"""
\bottomrule
\end{tabular}
\par\smallskip
\begin{minipage}{0.95\linewidth}\scriptsize All estimates were refitted on the same data version with two-way bootstrap intervals ($B=2000$). The renewal-type ratio became the primary statistic in round 2 because its relationship to $R_t$ is explicit through the generation interval; the weighted growth ratio is retained for continuity. The point estimate of each specification lies inside the interval of the other specification of the same statistic.\end{minipage}
\end{table}
"""
if ZH:
    for a, b in [
        (r"\caption{Specifications used across the rounds of this analysis, with the 95\% confidence interval of each estimate on the data version of 2026-10-20.}",
         r"\caption{本分析各轮所用的设定，以及各估计值在 2026-10-20 数据版本上的 95\% 置信区间。}"),
        ("Round and specification & Role &", "轮次与设定 & 角色 &"),
        ("Round 1: $G$, three calendar windows ($S=3$) & Primary in round 1; secondary from round 2",
         "第 1 轮：$G$，三个日历窗口（$S=3$） & 第 1 轮主分析；第 2 轮起为次要分析"),
        ("Round 2, rule version 1: windows ending at the season maximum & Abandoned (four weeks in 2024--25); no estimates retained",
         "第 2 轮，规则第 1 版：窗口止于流行季最大值 & 已弃用（2024--25 季仅四周）；未保留估计"),
        ("Round 2, rule version 2 (2026-10-20): $G$, rule windows ($S=4$) & Secondary",
         "第 2 轮，规则第 2 版（2026-10-20）：$G$，规则窗口（$S=4$） & 次要分析"),
        ("Round 2, rule version 2: $R_A$, rule windows ($S=4$) & Primary from round 2",
         "第 2 轮，规则第 2 版：$R_A$，规则窗口（$S=4$） & 第 2 轮起为主分析"),
        ("Round 2: $R_A$, calendar windows ($S=3$) & Sensitivity", "第 2 轮：$R_A$，日历窗口（$S=3$） & 敏感性分析"),
        ("Round 3: two-step inversion & Calibration of the round-2 fits; no refit", "第 3 轮：两步反演 & 校准第 2 轮拟合；未重新拟合"),
    ]:
        assert a in s9, a
        s9 = s9.replace(a, b)
    s9 = re.sub(r"\\scriptsize .*?\\end\{minipage\}", lambda m: r"\scriptsize 全部估计均在同一数据版本上重新拟合，区间来自双向自助法（$B=2000$）。第 2 轮起以更新型比值为主统计量，因为它通过代际间隔与 $R_t$ 有明确的联系；加权增长比予以保留，以保持连续性。同一统计量的各设定，其点估计均落在另一设定的区间之内。\end{minipage}", s9, flags=re.S)
(TAB / "tab_s9_analysis_history.tex").write_text(s9)

# --- round-3 edits of existing tables -----------------------------------------------------------------
def edit(name, fn):
    p = TAB / name
    p.write_text(fn(p.read_text()))


def nosep(t):
    return re.sub(r"(\d)\{,\}(\d)", r"\1\2", t)


edit("tab_s1_headline.tex", nosep)

FOLDS = {"RA_rule4": 4, "G_rule4": 4, "G_cal3": 3}


def s2(t):
    t = t.replace("Decomposition of the intercept and the plateau versus power-law comparison.",
                  "Decomposition of the intercept and comparison of the affine model with a power law, with the 95\\% confidence interval of each intercept.")
    t = t.replace(r"\begin{tabular}[t]{@{}c@{}}CV\\ folds\end{tabular}",
                  r"\begin{tabular}[t]{@{}c@{}}CV folds\\ (state; season)\end{tabular}")
    for spec, n in FOLDS.items():
        st = R[f"s03_decomposition.{spec}_cv_state_folds_affine_better"]["value"]
        se = R[f"s03_decomposition.{spec}_cv_season_folds_affine_better"]["value"]
        t = re.sub(rf"& {st}/51 \\\\", f"& {st}/51; {se}/{n} \\\\\\\\", t, count=1) if f"{st}/51;" not in t else t
    t = t.replace("test of the affine plateau against this power law", "test of the affine model against this power law")
    t = t.replace("CV folds counts the leave-one-state-out folds in which the affine model predicts better.",
                  "CV folds counts the leave-one-state-out and leave-one-season-out folds in which the affine model predicts better.")
    return t


def s2_zh(t):
    t = t.replace("截距分解及平台形式与幂律形式的比较。", "截距分解及仿射模型与幂律的比较，附各截距的 95\\% 置信区间。")
    t = t.replace(r"\begin{tabular}[t]{@{}c@{}}交叉验证\\ 折数\end{tabular}",
                  r"\begin{tabular}[t]{@{}c@{}}交叉验证折数\\（州；季）\end{tabular}")
    for spec, n in FOLDS.items():
        st = R[f"s03_decomposition.{spec}_cv_state_folds_affine_better"]["value"]
        se = R[f"s03_decomposition.{spec}_cv_season_folds_affine_better"]["value"]
        t = re.sub(rf"& {st}/51 \\\\", f"& {st}/51；{se}/{n} \\\\\\\\", t, count=1) if f"{st}/51；" not in t else t
    t = t.replace("$p$ 来自仿射平台对该幂律的参数自助检验", "$p$ 来自仿射模型对该幂律的参数自助检验")
    t = t.replace("交叉验证折数为逐州留一时仿射模型预测更优的折数。", "交叉验证折数为逐州留一与逐季留一时仿射模型预测更优的折数。")
    return t


def s6_zh(t):
    t = nosep(t)
    return t


if ZH:
    edit("tab_s2_decomposition.tex", s2_zh)
    edit("tab_s3_reporting.tex", lambda t: t.replace(
        "在拟合前确定（$G$ 取 $\\kappa=1$，$v_R=0.04$）",
        "由报告波动阈值与各统计量的样本均值 $\\bar{R}$ 确定（$G$ 取 $\\kappa=1$，$v_R=0.04$）"))
    edit("tab_s6_pathogens.tex", s6_zh)
else:
  edit("tab_s2_decomposition.tex", s2)
  edit("tab_s3_reporting.tex", lambda t: t.replace(
    "were fixed before fitting ($\\kappa=1$ for $G$, $v_R=0.04$).",
    "were set from a reporting-volatility threshold and the sample mean $\\bar{R}$ of each statistic ($\\kappa=1$ for $G$, $v_R=0.04$)."))


def s6(t):
    t = nosep(t)
    return t


if not ZH:
    edit("tab_s6_pathogens.tex", s6)

# --- manifest ------------------------------------------------------------------------------------------
mp = TAB / "manifest.json"
man = json.loads(mp.read_text())
if isinstance(man.get("tables"), list):
    ids = {e["id"]: e for e in man["tables"]}
    for tid, path, src in [("tab:s8", "tab_s8_inversion_settings.tex", ["experiment/analysis/parts/inversion.json"]),
                           ("tab:s9", "tab_s9_analysis_history.tex", ["experiment/analysis/results.json"])]:
        if tid not in ids:
            e = {"id": tid, "path": path, "sha256": "", "generator": GEN, "data_sources": src, "display": "raw",
                 "n_seeds": 1}
            man["tables"].append(e)
            ids[tid] = e
    for tid in ("tab:s1", "tab:s2", "tab:s3", "tab:s6"):
        if tid in ids:
            ids[tid]["generator"] = "experiment/code/make_tables.py; " + GEN
    for e in man["tables"]:
        if (TAB / e["path"]).exists():
            e["sha256"] = hashlib.sha256((TAB / e["path"]).read_bytes()).hexdigest()
    mp.write_text(json.dumps(man, indent=1) + "\n")
    print("tables written:", ", ".join(e["path"] for e in man["tables"]))
