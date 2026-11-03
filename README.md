# Inferential Targets and the Value of Case Volume in Reproduction-Number Surveillance
### 再生数监测中推断目标与病例量的统计价值

本仓库包含该研究的中英文论文定稿、全套 LaTeX 源码、插图附表，以及可复现的代码与实证数据。

---

## 📌 研究问题

监测实践中通常默认「报告病例越多，再生数估计越准」。但同一个估计量 $\Rhat_t=C_t/m$ 可以服务于两个不同的推断目标，而两者的误差性质并不相同：

- **当期值 $R_t$**（态势感知）：在报告概率固定且已知的基准模型下，误差是计数误差，随单期病例量累积而趋于零。一旦报告完整度本身波动，这一项还会叠加一个不随规模下降的方差（见下）。
- **跨期平均传播水平 $\bar R$**（季节典型值、政策评估基准）：误差中还包含 $R_t$ 在时期之间的**真实变异**。这不是增加单期病例能消除的计数噪声，而是该推断目标自身的抽样变异。

同一个误差分解给出
$$V(m)=v_R+\frac{A}{m},\qquad A=\bar R+\frac{\rho(\bar R^2+v_R)}{k},$$
其中 $v_R$ 为跨期传播变异，$m=\rho n$ 是固定报告概率 $\rho$ 下的设计病例量。**本文的核心主张因此是条件性的：扩大单期覆盖规模时，对跨期平均水平的估计存在边际收益递减；而不是「病例再多也无法准确估计当期 $R_t$」。**

由此定义两个不同性质的量，全文严格区分：

| 量 | 含义 | 性质 |
|---|---|---|
| $m_\times=A/v_R$ | 抽样贡献与环境贡献相等的位置 | 模型内的方差平衡尺度：越过它，单期病例翻倍的方差降幅不超过 25%。它不是纯生物学常数，而是依赖 $\Rbar$、$v_R$、$k$ 与报告概率 $\rho$ 的模型量 |
| $m^*=\sqrt{(c_0/c_1)m_\times}$ | 固定预算下的最优单期病例量 | **经济**切换门槛：在时期独立、成本线性且 $T$ 连续的前提下，只有 $m>m^*$ 时追加预算投向更多时期才优于投向更多病例 |

$m_\times$ 不是操作阈值，$m^*$ 依赖成本比 $c_0/c_1$ 与上述三项假设，两者不可混用。

---

## 🔬 理论与实证结构

### 理论：三类多期目标与一步预测

除过程均值 $\bar R$ 外，论文区分了另外两个常被混用的目标（推论 1）：

| 目标 | 均方误差 |
|---|---|
| 过程均值 $\bar R$，独立时期 | $v_R/T+A/(Tm)$ |
| 过程均值 $\bar R$，一阶自相关 $\phi$ | $v_R S(\phi,T)/T^2+A/(Tm)$ |
| 有限时期实现均值 $\bar R_T$ | $A/(Tm)$（环境项不出现） |
| 以 $\Rhat_t$ 朴素预测 $R_{t+1}$ | $A/m+2(1-\phi)v_R$ |
| 已知 $\bar R$ 时的最优线性预测 | $v_R[1-\phi^2 m/(m+m_\times)]$，极限 $(1-\phi^2)v_R$ |

报告完整度波动会额外引入一项不随 $m$ 下降的方差（命题 2）：
$$V_H(m)=v_R+v_H+\frac{A_H}{m},\qquad m_{\times,H}=A_H/(v_R+v_H)\le m_\times,$$
当且仅当相对报告波动 $\tau=0$ 时取等号。仅凭计数的均值与方差无法分离 $v_R$ 与 $v_H$。

### 实证：三层证据

| 层次 | 内容 | 能回答什么 |
|---|---|---|
| **模拟真值** | 已知 $R_t$ 时，$(\Rhat_t-R_t)^2\to 0$ 而 $(\Rhat_t-\bar R)^2\to v_R$ | 两类误差是否按理论分离 |
| **受控二项稀释** | 在固定 1,748 个州周上对观测计数做 $C^{(q)}\sim\mathrm{Bin}(C,q)$ | 随病例量下降的计数项是否按 $1/m$ 缩放 |
| **州固定效应** | 含 51 个州截距的拟合，用州内时间变化识别斜率 | 横截面斜率中混入多少空间聚合效应 |

**理论、数据与模拟三者承担的任务必须分开读：**

- 理论证明的是不同推断目标下的误差分解；
- 数据拟合的是**表观统计量的相对损失** $L_t=[(X_t-\bar X)/\bar X]^2$，其中 $X_t$ 是周增长统计量、$\bar X$ 是合并均值，因此**不是**对真实 $R_t$ 的估计误差；
- 模拟在指定机制下把表观量映射到潜在参数。

住院数据提供的是**与该结构相容**的观测尺度关系，不构成对两类真实估计误差的直接验证。

---

## 📊 主要结果与适用边界

**样本口径。** 清洗后的全期面板为 10,506 个州周（2020-08-08 至 2026-10-10，$D_t\ge5$）。**主分析使用的是四个流感季早期增长窗口内的 1,748 个州周（$R_A$）/ 1,756 个州周（$G$），共 42 个日历周、51 个辖区。** 两个数字含义不同，不能互相替代。

| 结果 | 数值 | 解释边界 |
|---|---|---|
| 表观平衡比（$R_A$，仿射模型） | 37.4 例住院/周 [19.1, 71.0] | 无截距幂律在观测范围内拟合**同样良好**；正截距及其派生量均以仿射形式为工作点 |
| 二项稀释速率 | $c=0.407$ [0.311, 0.531]，泊松基准 $b_P=0.446$ | 检验的是计数项的缩放，**不能**识别 $v_R$，也不能分离传播与报告 |
| 州内表观平衡比 | 17.3 [3.3, 51.4]，横截面 37.4 | 配对比值 2.16 [0.93, 7.75]，区间含 1，**差异不能与抽样误差区分** |
| 潜在环境方差 $v_R$ | 0.040 [0.026, 0.056]（代际间隔 3.6 天） | 条件于 3.6 天序列间隔代理与独立周冲击；2.2–2.9 天下为 0.013–0.026 |
| 潜在周平衡规模 | 58.9 [42.3, 90.0] 例/周 | 44% [35, 52] 的州周位于其上；该比例随代际间隔强烈变动 |
| 截距的全国周共有份额 | $R_A$ 41%、$G$ 51% | 是**截距**中的份额，不是总方差份额；可能来自共同趋势、报告或收治实践 |

**三项必须保留的限定：**

1. 环境变异、报告波动、季内共同趋势与统计量构造都能产生正截距。现有数据无法完全分离这些来源。
2. 跨州横截面的 $1/D$ 斜率同时包含病例量效应与空间聚合效应。若后者成立，$b/a$ 描述的是州级汇总网络的表观尺度，而非单一报告单元内部的因果规模收益。
3. 自愿报告期与强制报告期的对比仅来自一个非流行季，比值 1.58 [0.95, 2.23] 区间含 1，只作描述性线索。

---

## 📂 仓库结构

```text
├── paper_en/                     # 英文论文（目标期刊：Epidemics）
│   ├── paper.pdf                 # 36 页定稿
│   ├── paper.html                # 自包含阅读版
│   ├── main.tex, compile.sh      # 编译入口与脚本
│   ├── bibliography.bib, *.bst   # 66 条参考文献与样式表
│   ├── sections/                 # 9 个章节源码（含附录全部证明）
│   ├── tables/                   # 附表 S1–S11 源码
│   └── figures/                  # 10 张插图（9 矢量 PDF + 1 位图）
│
├── paper_zh/                     # 中文论文（数字与英文定稿一致）
│   ├── paper.pdf                 # 35 页定稿
│   ├── paper.html
│   ├── main.tex, compile.sh
│   ├── bibliography.bib, *.bst
│   ├── sections/                 # 9 个章节源码
│   ├── tables/                   # 附表 S1–S11 源码（中文）
│   └── figures/                  # 10 张插图
│
├── derivation/                   # 独立推导手册（29 页）及源码
│
└── experiment/                   # 实验与数据复现工程
    ├── code/                     # 清洗、分析、反演与制图脚本（s01–s10, s99）
    │   ├── run_all.sh            # 全流程脚本
    │   ├── cc_common.py          # 公共数据接口
    │   ├── s01_data_panel.py ~ s10_polish.py
    │   └── make_fig*.py, make_tables*.py
    ├── data/
    │   ├── DATA_MANIFEST.md      # 数据来源、快照日期与 SHA-256
    │   ├── panel_nhsn_vintage_2026-10-20.csv.gz  # 清洗后全期面板（10,506 州周）
    │   ├── national_flu.csv
    │   └── external/             # NHSN 与 Delphi 原始快照
    └── analysis/                 # 计算结果
        ├── results.json          # 汇总结果
        ├── results_summary.md
        ├── vintage_r4_results.json
        └── parts/                # 分模块输出
```

---

## 🛠️ 编译指南

**两篇论文都要求 XeLaTeX**（`fontspec` 与 `unicode-math`），`pdflatex` 会直接报错退出。

### 推荐：Tectonic（自包含，无需本地 TeX 发行版）

```bash
cd paper_en && bash compile.sh     # 输出 paper.pdf
cd paper_zh && bash compile.sh     # 输出 paper.pdf
```

`compile.sh` 优先探测 `tectonic`，缺失时回退到本机 `xelatex` + `bibtex`。

### 手动编译（需已安装 TeX Live / MacTeX）

```bash
cd paper_en
xelatex main.tex && bibtex main && xelatex main.tex && xelatex main.tex
```

中文版同理（`paper_zh/` 内），字体需支持 `SimSun` / `Noto Serif CJK SC`。

---

## 🔁 复现指南

环境要求：`Python >= 3.10`，依赖 `numpy pandas scipy matplotlib`。

```bash
cd experiment/code
bash run_all.sh          # 复用已有分块结果
bash run_all.sh --fresh  # 清空缓存后全量重跑
```

脚本依次执行：

1. 数据面板构建与质量审计（`s01_data_panel.py`）
2. 双向聚类推断与损失模型拟合（`s02_inference.py`）
3. 多尺度方差分解（`s03_decomposition.py`）
4. 强制 / 自愿报告制度分析（`s04_reporting.py`）
5. 流行季合并与分母权重检验（`s05_seasons.py`、`s06_pathogens.py`）
6. 理论闭式量验证与蒙特卡罗（`s07_theory.py`）
7. 表观统计量模拟与两步反演（`s08`–`s10`）
8. 结果汇总、制图与附表生成（`s99_collect.py`、`make_*.py`）

**注意**：`s06_pathogens.py` 需要先前反演轮次的 `r4_results.json`。脚本会依次探测 `analysis/vintage_r4_results.json` 与 `vintage_2026-10-20/analysis/r4_results.json`，两者都缺失时会明确报错说明。

---

## 📄 引用

若使用本仓库的代码或数据，请引用论文 *Two inferential targets, different returns to case volume: theory and evidence for epidemic reproduction numbers*。数据来源与版本见 `experiment/data/DATA_MANIFEST.md`。
