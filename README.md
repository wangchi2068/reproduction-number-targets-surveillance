# 《Inferential targets for reproduction-number estimates》终版论文交付包

本目录（`D:\trae\114\论文`）仅包含学术论文最终版的有效工作文件，已剥离所有中间审稿记录与答辩过程文件（历轮审稿意见与修改说明已归档至上级目录 `D:\trae\114\审稿档案`）。

---

## 目录结构

```
D:\trae\114\论文\
│
├── paper.pdf                  # 【最终论文手稿】30页完整出版级 PDF（包含带 DOI 参考文献与完整附录，0 错误编译通过）
├── paper.html                 # 【离线交互阅读版】自包含网页版手稿（离线内嵌高清图表与数学公式，浏览器双击即开）
├── highlights.txt             # 【期刊投稿亮点】针对期刊投稿的 4 条核心学术贡献 Highlights
├── README.md                  # 【本文档】论文结构与复现指南
│
├── paper/                     # 【LaTeX 排版源码包】可直接本地编译或打包上传 Overleaf / 期刊投稿系统
│   ├── main.tex               # LaTeX 主入口文件
│   ├── compile.sh             # 编译脚本
│   ├── bibliography.bib       # 包含 38 条官方 Crossref DOI 的完整 BibTeX 文献数据库
│   ├── unsrtnat-noid.bst      # BibTeX 参考文献排版样式文件
│   ├── sections/              # 正文分章节源码（abstract, introduction, model, diagnostics 等 11 个模块）
│   ├── figures/               # 出版级矢量图件（fig_01 ~ fig_07.pdf）
│   └── tables/                # 自动渲染的数据表格源码（18 张 .tex 表格与 manifest.json）
│
└── experiment/                # 【计算实验与数据包】全量支持代码、输入数据与全套结果
    ├── code/                  # Python 核心计算脚本（ma_audit.py, ma_extensions.py, s03*.py 等）
    ├── analysis/              # 461 项估计值汇总 results.json 及子模型输出
    └── data/                  # 真实流感面板数据与 SHA-256 校验清单（DATA_MANIFEST.md）
```

---

## 核心科学结论概要

1. **推断目标决定病例量价值**：
   - 当期瞬时值 $R_t$ 的误差随病例数 $m$ 趋零；
   - 有限期实现均值 $\bar{R}_T$ 的均方误差为 $A/(Tm)$，无环境方差底线；
   - 过程均值 $\bar{R}$ 受跨期抽样限制，保留不可消除的环境变异底线 $v_R/T$；
   - 观测分母带噪使得相对方差上升约 1.5 倍，收缩半权尺度由 53.1 升至 78.6。
2. **报告波动的单期不可识别性（命题 2）**：
   - 报告率波动留下不可消除的误差下限 $v_H$，单期横截面矩无法分离生物学传播变异与统计报告波动。
3. **已知真值审计揭示诊断本质**：
   - 截距 $a$ 回收的是与计数无关的相对离散度（实测 $0.0150$），而非估计误差底线；
   - 构造 **S8 Combo**（尺度异质性 + 报告波动）成功同时再现了截距（0.0148）、交叉计数（41.7 覆盖实测 37.4）与斜率比（1.298 覆盖实测 1.257），且曲线对比值分位数提升至 8.0%；
   - 上尾加权敏感性证实 P90 差异对截距回收影响 $<4\%$；四季截距在对数尺度 Wald 检验下未拒绝同质性（$W=5.52, p=0.138$）。

---

## 本地编译指南

- **使用 Tectonic 编译**：
  ```bash
  cd "D:\trae\114\论文\paper"
  tectonic main.tex
  ```
- **使用 XeLaTeX / BibTeX 编译**：
  ```bash
  cd "D:\trae\114\论文\paper"
  xelatex -interaction=nonstopmode main.tex
  bibtex main
  xelatex -interaction=nonstopmode main.tex
  xelatex -interaction=nonstopmode main.tex
  ```
