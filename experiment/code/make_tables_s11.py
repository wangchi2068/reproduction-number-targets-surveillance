"""make_tables_s11.py: Generate Table S11 (Binomial Thinning and State Fixed Effects)
Reads analysis/parts/s11_thinning.json and analysis/parts/s11_fe.json.
Generates:
- paper_en/tables/tab_s11_thinning_fe.tex
- paper_zh/tables/tab_s11_thinning_fe.tex
Updates manifest.json in both directories.
Usage: python3 make_tables_s11.py
"""
import hashlib
import json
from pathlib import Path

EXP = Path(__file__).resolve().parents[1]
ANALYSIS = EXP / "analysis"
PARTS = ANALYSIS / "parts"

thin_file = PARTS / "s11_thinning.json"
fe_file = PARTS / "s11_fe.json"

if not thin_file.exists() or not fe_file.exists():
    raise FileNotFoundError("Missing s11_thinning.json or s11_fe.json. Run s03b_thinning.py and s03c_state_fe.py first.")

THIN = json.loads(thin_file.read_text())["blocks"]["thinning_table_s11a"]
FE = json.loads(fe_file.read_text())["blocks"]["fe_table_s11b"]

REPO = EXP.parent
TARGETS = [
    ("en", REPO / "paper_en" / "tables"),
    ("zh", REPO / "paper_zh" / "tables"),
    ("en_exp", EXP / "tables" / "en"),
    ("zh_exp", EXP / "tables" / "zh")
]


def f(v, d=3):
    if v is None:
        return "--"
    if v < 0:
        return f"$-${abs(v):.{d}f}"
    return f"{v:.{d}f}"


def ci_str(rec, d=3):
    v = rec["value"]
    c = rec.get("ci95")
    if not c:
        return f(v, d)
    return f"{f(v, d)} [{f(c[0], d)}, {f(c[1], d)}]"


def generate_table(lang="en"):
    is_zh = lang.startswith("zh")
    
    # Panel A rows
    rows_a = []
    qs = [0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
    for q in qs:
        d = THIN["levels"][str(q)]
        rows_a.append(f"{q:.1f} & {ci_str(d['b'])} & {ci_str(d['cq'])} & {ci_str(d['da_e3'], 2)} & {ci_str(d['c50'])} \\\\")
        
    p_c = ci_str(THIN["pooled_c"])
    p_c50 = ci_str(THIN["pooled_c_50"])
    bP_all = f"{THIN['bP_all']:.3f}"
    bP_c50 = f"{THIN['bP_c50']:.3f}"
    
    # Panel B rows
    p_b = ci_str(FE["pooled"]["b"])
    p_ba = ci_str(FE["pooled"]["ba"], 1)
    
    fe_slope = ci_str(FE["state_fe"]["slope"])
    r_mean = ci_str(FE["state_fe"]["r_mean"], 1)
    r_median = ci_str(FE["state_fe"]["r_median"], 1)
    r_row = ci_str(FE["state_fe"]["r_row"], 1)
    r_case = ci_str(FE["state_fe"]["r_case"], 1)
    
    r_D_slope = ci_str(FE["state_fe_D"]["slope"])
    r_D = ci_str(FE["state_fe_D"]["ratio"], 1)
    
    w_slope = ci_str(FE["week_only"]["slope"])
    sw_slope = ci_str(FE["state_and_week"]["slope"], 2)
    paired = ci_str(FE["state_fe"]["paired_ratio"], 2)

    if not is_zh:
        caption = (
            "Binomial thinning of the renewal-type ratio sample and fits with state intercepts. "
            "Panel A gives, for each retention probability $q$, the slope on the original scale $C_t$, "
            "the implied thinning slope $c_q=[b(q)-b(1)]/(1/q-1)$, the change in intercept $a(q)-a(1)$ and $c_q$ "
            "in the subset with $C_t\\ge50$ ($n=828$). Panel B gives the slope on $1/C_t$ and the ratio of slope "
            "to intercept under each specification. Values are point estimates with 95\\% two-way bootstrap intervals in brackets."
        )
        panel_a_head = r"\multicolumn{5}{l}{\textit{A. Binomial thinning ($n=1748$)}}\\"
        col_head_a = r"$q$ & $b(q)$ & $c_q$ & $\Delta a$ ($\times10^{-3}$) & $c_q$, $C_t\ge50$ \\"
        pooled_label = r"Pooled slope $c$"
        poisson_label = r"Poisson coefficient $b_P$"
        
        panel_b_head = r"\multicolumn{3}{l}{\textit{B. Within-state fixed effects ($n=1748$, 51 states)}}\\"
        col_head_b = r"Specification & Slope & Ratio \\"
        spec_labels = [
            ("Pooled cross-section", p_b, p_ba),
            ("State intercepts, mean $a_s$", fe_slope, r_mean),
            ("State intercepts, median $a_s$", "--", r_median),
            ("State intercepts, row-weighted $a_s$", "--", r_row),
            ("State intercepts, case-weighted $a_s$", "--", r_case),
            ("State intercepts, scale $D_t$", r_D_slope, r_D),
            ("Week intercepts only", w_slope, "--"),
            ("State and week intercepts", sw_slope, "--"),
            ("Paired ratio, cross-section to mean $a_s$", "--", paired)
        ]
    else:
        caption = (
            "更新型比值样本的二项稀释与州固定效应拟合。"
            "表 A 列出每个保留概率 $q$ 下原始尺度 $C_t$ 的斜率、隐含稀释斜率 $c_q=[b(q)-b(1)]/(1/q-1)$、"
            "截距变化 $a(q)-a(1)$ 以及 $C_t\\ge50$ 子集 ($n=828$) 的 $c_q$。"
            "表 B 列出各模型下 $1/C_t$ 的斜率与平衡比。括号内为 95\\% 双向自助法置信区间。"
        )
        panel_a_head = r"\multicolumn{5}{l}{\textit{A. 受控二项稀释实验 ($n=1748$)}}\\"
        col_head_a = r"$q$ & $b(q)$ & $c_q$ & $\Delta a$ ($\times10^{-3}$) & $c_q$, $C_t\ge50$ \\"
        pooled_label = r"合并稀释斜率 $c$"
        poisson_label = r"泊松基准系数 $b_P$"
        
        panel_b_head = r"\multicolumn{3}{l}{\textit{B. 州固定效应模型 ($n=1748$, 51 个辖区)}}\\"
        col_head_b = r"模型设定 & 斜率 $b$ & 平衡比 $b/a$ \\"
        spec_labels = [
            ("合并横截面回归 (Pooled)", p_b, p_ba),
            ("州固定效应，算术均值 $a_s$", fe_slope, r_mean),
            ("州固定效应，中位数 $a_s$", "--", r_median),
            ("州固定效应，行加权均值 $a_s$", "--", r_row),
            ("州固定效应，病例加权均值 $a_s$", "--", r_case),
            ("州固定效应，分母规模 $D_t$", r_D_slope, r_D),
            ("仅周固定效应 (Week intercepts)", w_slope, "--"),
            ("州与周双向固定效应", sw_slope, "--"),
            ("横截面与州内均值配对比值", "--", paired)
        ]

    rows_b = [f"{lab} & {sl} & {rt} \\\\" for lab, sl, rt in spec_labels]

    tex = rf"""\begin{{table}}[htbp]
\centering
\footnotesize
\setlength{{\tabcolsep}}{{4pt}}
\caption{{{caption}}}
\label{{tab:s11}}
\begin{{tabular}}{{lllll}}
\toprule
{panel_a_head}
{col_head_a}
\midrule
""" + "\n".join(rows_a) + rf"""
\midrule
{pooled_label} & -- & {p_c} & -- & {p_c50} \\
{poisson_label} & -- & {bP_all} & -- & {bP_c50} \\
\bottomrule
\end{{tabular}}

\bigskip
\begin{{tabular}}{{lll}}
\toprule
{panel_b_head}
{col_head_b}
\midrule
""" + "\n".join(rows_b) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
    return tex


for lang, target_dir in TARGETS:
    target_dir.mkdir(parents=True, exist_ok=True)
    tex_content = generate_table(lang)
    out_tex = target_dir / "tab_s11_thinning_fe.tex"
    out_tex.write_text(tex_content, encoding="utf-8")
    print(f"Generated {out_tex}")
    
    # Update manifest.json if exists
    mp = target_dir / "manifest.json"
    if mp.exists():
        try:
            man = json.loads(mp.read_text(encoding="utf-8"))
            if isinstance(man.get("tables"), list):
                tables = man["tables"]
                ids = {x.get("id"): x for x in tables}
                entry = {
                    "id": "tab:s11",
                    "path": "tab_s11_thinning_fe.tex",
                    "sha256": hashlib.sha256(out_tex.read_bytes()).hexdigest(),
                    "generator": "experiment/code/make_tables_s11.py",
                    "data_sources": [
                        "experiment/analysis/parts/s11_thinning.json",
                        "experiment/analysis/parts/s11_fe.json"
                    ],
                    "display": "raw",
                    "n_seeds": 1
                }
                if "tab:s11" in ids:
                    for i, t in enumerate(tables):
                        if t.get("id") == "tab:s11":
                            tables[i] = entry
                else:
                    tables.append(entry)
                # recalculate hashes for all files
                for t in tables:
                    fp = target_dir / t["path"]
                    if fp.exists():
                        t["sha256"] = hashlib.sha256(fp.read_bytes()).hexdigest()
                mp.write_text(json.dumps(man, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
                print(f"Updated manifest {mp}")
        except Exception as err:
            print(f"Warning: could not update manifest {mp}: {err}")
