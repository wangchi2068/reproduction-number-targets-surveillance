"""Final-polish paper table S10: lag-1 autocorrelation in data and simulation, and the four added inversion
settings (true GI 2.9 d, admission VMR 1.5 and 2.0, AR phi = 0.25). Every number is read from
analysis/parts/polish.json or analysis/parts/inversion.json; nothing is recomputed.
Usage: python3 make_tables_polish.py <paper_run_dir> [--lang zh]   (writes <paper_run_dir>/tables, updates manifest)
"""
import hashlib
import json
import sys
from pathlib import Path

EXP = Path(__file__).resolve().parents[1]
ANALYSIS = EXP / "analysis"
if not (ANALYSIS / "parts" / "polish.json").exists():  # published copy lives in output/<slug>/experiment
    ANALYSIS = EXP.parents[1] / "output" / EXP.parent.name / "experiment" / "analysis"
POL = json.loads((ANALYSIS / "parts" / "polish.json").read_text())
INV = json.loads((ANALYSIS / "parts" / "inversion.json").read_text())["blocks"]["sensitivity"]
N = POL["numbers"]
R1S = POL["blocks"]["r1_simulated_v0.04"]
R1E = POL["blocks"]["r1_empirical"]
NEW = POL["blocks"]["new_settings"]
TAB = Path(sys.argv[1])
if (TAB / "tables").exists():
    TAB = TAB / "tables"
ZH = "--lang" in sys.argv and sys.argv[sys.argv.index("--lang") + 1] == "zh"
GEN = "experiment/code/make_tables_polish.py"


def f(v, d=3):
    return ("$-$" if v < 0 else "") + f"{abs(v):.{d}f}"


def ci(lo, hi, d=3):
    return f"[{f(lo, d)}, {f(hi, d)}]"


LAB = {"data": "Data (rule windows)", "M0": "Data-matched design (M0)", "ar025": "AR(1), $\\phi=0.25$",
       "ar05": "AR(1), $\\phi=0.5$", "gi_short": "True GI mean 2.2 d", "gi_cow_lo": "True GI mean 2.9 d",
       "gi_long": "True GI mean 4.3 d", "od05": "Admissions VMR 1.5", "od10": "Admissions VMR 2.0",
       "count_shock": "Synchronous count shock"}
LAB_ZH = {"data": "数据（规则窗口）", "M0": "数据匹配设计（M0）", "ar025": "AR(1)，$\\phi=0.25$",
          "ar05": "AR(1)，$\\phi=0.5$", "gi_short": "真实代际间隔均值 2.2 天", "gi_cow_lo": "真实代际间隔均值 2.9 天",
          "gi_long": "真实代际间隔均值 4.3 天", "od05": "住院数方差均值比 1.5", "od10": "住院数方差均值比 2.0",
          "count_shock": "报告计数同步冲击"}
if ZH:
    LAB = LAB_ZH

# Panel A: lag-1 autocorrelation
ra = []
e = R1E
ra.append(f"{LAB['data']} & {f(e['RA']['r1'])} {ci(*e['RA']['r1_ci95'])} & -- & "
          f"{f(e['G']['r1'])} {ci(*e['G']['r1_ci95'])} & -- \\\\")
for key in ("M0", "ar025", "ar05", "gi_short", "gi_cow_lo", "gi_long", "od05", "od10", "count_shock"):
    s = R1S[key]
    ra.append(f"{LAB[key]} & {f(s['RA']['mean'])} {ci(*s['RA']['q95'])} & {100 * s['RA']['emp_percentile']:.1f} & "
              f"{f(s['G']['mean'])} {ci(*s['G']['q95'])} & {100 * s['G']['emp_percentile']:.1f} \\\\")

# Panel B: implied variance and b/a at v_R = 0.04 in the added settings, with M0 for reference
rb = []
m0 = INV["M0"]
rb.append(f"{LAB['M0']} & {f(m0['RA']['v_implied'])} {ci(*m0['RA']['v_ci95'])} & "
          f"{f(m0['G']['v_implied'])} {ci(*m0['G']['v_ci95'])} & "
          f"{N['M0_RA_ba_pooled_v0.04']['value']:.1f} & {N['M0_G_ba_pooled_v0.04']['value']:.1f} \\\\")
for key in ("gi_cow_lo", "od05", "od10", "ar025"):
    s = NEW[key]
    rb.append(f"{LAB[key]} & {f(s['RA']['v_implied'])} {ci(*s['RA']['v_ci95'])} & "
              f"{f(s['G']['v_implied'])} {ci(*s['G']['v_ci95'])} & "
              f"{s['RA']['ba_pooled_by_v']['0.04']:.1f} & {s['G']['ba_pooled_by_v']['0.04']:.1f} \\\\")

if ZH:
    cap = (r"\caption{滞后一阶自相关以及新增反演设定下的隐含环境方差与 $b/a$。}")
    ha = r"设定 & $R_A$ 的 $r_1$ & 数据分位 & $G$ 的 $r_1$ & 数据分位 \\"
    hb = r"设定 & 由 $R_A$ 得到的 $v_R^{\mathrm{impl}}$ & 由 $G$ 得到的 $v_R^{\mathrm{impl}}$ & $R_A$ 的 $b/a$ & $G$ 的 $b/a$ \\"
    pa, pb = "A. 滞后一阶自相关", r"B. 新增设定（$b/a$ 取 $v_R=0.04$ 处的汇总值）"
    note = (r"\scriptsize $r_1$ 为州--季节块内相对偏差 $x/\bar{x}_{\text{块}}-1$ 的汇总滞后一阶自相关。"
            r"数据行的区间来自按州整群的自助抽样（1000 次）；模拟行给出 $v_R=0.04$ 时各面板 $r_1$ 的均值及 2.5\%--97.5\% 分位范围"
            r"（新设定每个 $v_R$ 取 200 个面板，其余设定取 100 个），“数据分位”为数据值在模拟分布中的百分位。"
            r"分析者所用代际间隔均值始终为 3.6 天；2.9 天与 4.3 天对应 \citet{cowling2009estimation} 所估流感序列间隔的 95\% 区间端点。"
            r"住院数方差均值比通过对每周住院数施加伽马--泊松混合实现。")
else:
    cap = r"\caption{Lag-1 autocorrelation, and implied environmental variance and $b/a$ in the added inversion settings.}"
    ha = r"Setting & $r_1$ of $R_A$ & Data pct. & $r_1$ of $G$ & Data pct. \\"
    hb = r"Setting & $v_R^{\mathrm{impl}}$ from $R_A$ & $v_R^{\mathrm{impl}}$ from $G$ & $b/a$ of $R_A$ & $b/a$ of $G$ \\"
    pa, pb = "A. Lag-1 autocorrelation", r"B. Added settings ($b/a$ pooled at $v_R=0.04$)"
    note = (r"\scriptsize $r_1$ is the pooled lag-1 autocorrelation of the relative deviation $x/\bar{x}_{\mathrm{block}}-1$ "
            r"within state-season blocks. The interval for the data comes from a state-cluster bootstrap with 1000 draws; "
            r"simulated rows give the mean and the 2.5\% to 97.5\% range of $r_1$ over panels at $v_R=0.04$ "
            r"(200 panels per value of $v_R$ for the added settings, 100 for the others), and Data pct.\ is the percentile "
            r"of the empirical value in that distribution. The analyst's generation interval stays at a mean of 3.6 days; "
            r"2.9 and 4.3 days are the ends of the 95\% interval for the serial interval of influenza in \citet{cowling2009estimation}. "
            r"Overdispersed admissions use a gamma--Poisson mixture on weekly counts.")

tex = (r"""\begin{table}[htbp]
\centering
\footnotesize
\setlength{\tabcolsep}{4pt}
""" + cap + r"""
\label{tab:s10}
\begin{tabular}{lllll}
\toprule
\multicolumn{5}{l}{\textit{""" + pa + r"""}}\\
""" + ha + r"""
\midrule
""" + "\n".join(ra) + r"""
\midrule
\multicolumn{5}{l}{\textit{""" + pb + r"""}}\\
""" + hb + r"""
\midrule
""" + "\n".join(rb) + r"""
\bottomrule
\end{tabular}
\par\smallskip
\begin{minipage}{0.95\linewidth}""" + note + r"""\end{minipage}
\end{table}
""")
(TAB / "tab_s10_polish.tex").write_text(tex)

mp = TAB / "manifest.json"
man = json.loads(mp.read_text())
if isinstance(man.get("tables"), list):
    ids = {x["id"]: x for x in man["tables"]}
    if "tab:s10" not in ids:
        man["tables"].append({"id": "tab:s10", "path": "tab_s10_polish.tex", "sha256": "", "generator": GEN,
                              "data_sources": ["experiment/analysis/parts/polish.json",
                                               "experiment/analysis/parts/inversion.json"],
                              "display": "raw", "n_seeds": 1})
    for x in man["tables"]:
        if (TAB / x["path"]).exists():
            x["sha256"] = hashlib.sha256((TAB / x["path"]).read_bytes()).hexdigest()
    mp.write_text(json.dumps(man, indent=1) + "\n")
print("written", TAB / "tab_s10_polish.tex")
