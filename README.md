# 论文修订与三轮审稿：重复数目标与离散诊断（"What more cases buy"）

审阅对象是 GitHub 仓库 `wangchi2068/reproduction-number-targets-surveillance`（`main` 分支，提交 `2f5c407`）中的稿件。本文件夹包含三轮审稿意见、两轮修订稿与逐条答复：

- 第一轮审稿意见 `peer_review.md`（大修，4/10）及逐条答复 `revision_response.md`；
- 第二轮复审意见 `peer_review_revised.md`（大修但范围缩小，6/10）及逐条答复 `revision_response_r2.md`；
- 第三轮复审意见 `peer_review_revised_r2.md`（中等程度修改，**7/10**），审阅对象是当前发布的修订稿 `paper/paper.pdf`（42 页）。

用户对本轮的质量目标是 8 分及以上。本轮评分为 7/10，**这一目标未达成**，原因与达到 8 分所需的修改见下方关键发现第 1 条。修订遵守一个原则：显著性只能来自锁定的分析计划，所有校准和敏感性分析都要报告，不利结果如实写入。

上一级 output 文件夹中的总结报告 final-report（PDF、DOCX、HTML、Markdown，同一来源）已在最终报告步骤按三轮的最终状态重写，汇总 4→6→7 的评分、全部结果与达到 8 分的路径；数字与本文件夹一致。

## 关键发现

1. **第三轮复审结论：中等程度修改，7/10，比第二轮高 1 分，未达到用户要求的 8 分。** 上升的依据：第二轮的 N1（“预设分析”标签过强）与 N2（截距依赖函数形式）在报告层面已妥善处理，N3–N9 全部解决，独立重算的关键数字全部与稿件一致。阻碍 8 分的是四个研究层面的问题：人口梯度的幅度依赖仿射形式（N10）；“低于合计 6–22%”没有一个通过检验的生成器支撑（N11、M4）；截距中的传播份额不可识别（M5）；b/b_P 区间以泊松参照固定为条件（N12）。达到 8 分需要：写入并锁定修正案 3 后在三参数模型中加入人口协变量并报告同计数比较、据此改写摘要中的“halved”；把偏差幅度改为覆盖全部模拟机制的范围；对 b/b_P 联合重抽参照。若再在第二条 NHSN 序列上复现，预期可达 8–9 分。证据：`peer_review_revised_r2.md` 第 3、5 节，第 4.3 节 N10–N15。
2. **主检验 E1 在本轮再次独立复现。** a = 0.0149546（95% CI 0.010852–0.019058），t（41 自由度）= 7.3608，p = 5.12×10^-9，Holm 调整后 1.02×10^-8，与正文逐位一致。自写实现（x̄ 视为已知）给出相同的 a，SE 大约 10%，p = 4.4×10^-8，结论不变。证据：`.work/rn-targets-peer-review/rereview2/repro_checks_r2.md` 第 1 节；同目录 `repro_checks_r2.json`。
3. **第二轮新增结果全部复现，修正案 2 的锁定早于全部新结果。** 截距随州人口下降（a1 = −0.00538，P90/P10 截距比 0.50）；人口三分位分别拟合 0.0189/0.0195/0.0135（三分位内重算 x̄ 时与作者完全一致）；三参数模型 a = 0.0089、γ = 0.617，周自助区间下端低于 0（3.6% 的抽样 ≤ 0）；计数 P90 处无模型下限相对误差 10.5%，扣除计数项后 10.1%；钳位比例 13.5%。修正案 2 锁定于 10:16:29 UTC（SHA-256 `8f00e540…`，22,263 字节），最早的新分析日志在 10:20:41。证据：`repro_checks_r2.md` 第 1 节；`repro_tertile_check.json`；`repro_clamp_diag.json`。
4. **同计数比较只部分支持人口梯度。** 在相同计数段内比较扣除计数项后的平均损失：大辖区在四个计数段中都不高于中等辖区，100 例以上低 37–50%；但小辖区与大辖区之差只在 20–49 例段可分辨（0.0111，州自助 95% 区间 0.0005–0.0251），50–199 例段区间含 0，100–199 例段点值反向。稿件“小辖区截距约为大辖区两倍”的幅度依赖仿射外推，这是 N10 的依据。证据：`repro_checks_r2.md` 第 4 节；`repro_matched_count.json`；`peer_review_revised_r2.md` N10。
5. **bootstrap-t 的结果对夹心估计的写法敏感，但作者的结果成立。** 同一组 Webb 野自助抽样（B = 999）下，作者计入 x̄ 估计的堆叠夹心得到 p = 0.001（B = 999 可达到的最小值），与作者的 p = 0.0001（B = 9,999）一致；x̄ 视为已知的夹心在大量钳位的抽样上 SE 退化（约 5×10^-5），得到 p = 0.020。审稿意见把这一敏感性列为次要问题（N13），不构成对主结果的反驳。证据：`repro_checks_r2.md` 第 1 节“钳位下的 bootstrap-t”；`repro_clamp_diag.json`。
6. **核验发现并更正了一处附录事实错误。** 附录 C.3 讨论 pigeonhole 自助 SD 比 1.78 超过 Owen 上界约 1.73 [50] 时，交付步骤第一次改稿写成“100 个面板、相对标准误约 7%”；`results.json` 中实为 200 个面板，相对标准误约 5%。已更正并重新编译、发布，超出的 3% 仍在 Monte Carlo 误差之内，结论不变。证据：`repro_checks_r2.md` 第 2 节；`paper/paper.pdf` 附录 C.3。
7. **第二轮答复留下的发布与措辞问题已关闭，作者方投稿前事项单列、不计入评分。** 发布的 `experiment/analysis/` 现在包含分析计划、修正案 2、两个锁文件与含第二轮结果的 `results.json`；表格中的 “pre-specified” 与 “plus reporting total” 已改，表 15 脚注日期、p 值精度与零差异拒绝率措辞已更正。仍保留的是中国附录的 legacy 代码与数据（N15，已在 `experiment/code/README.md` 标为不被当前稿件引用）。另外，作者署名与单位仍是占位文字，CRediT 未对应到人，计划与锁文件没有带时间戳的公共存档 DOI，文献表没有 DOI。这些需要真实的作者名单与作者账户，本任务不能生成，也没有编造。证据：`repro_checks_r2.md` 第 3 节；`revision_response_r2.md` 第 7 节；`peer_review_revised_r2.md` 第 5 节“作者方投稿前事项”。

## 文件导航

| 路径 | 内容 |
|---|---|
| `review/peer_review_revised_r2.md` | 第三轮复审意见（中文；原意见与新问题的状态总表、独立核验表、新问题 N10–N15、给作者的问题、7/10 及达到 8 分所需修改） |
| `paper/paper.pdf` | 当前修订稿 PDF（42 页，xelatex 编译，SHA-256 `261c01da…a117`） |
| `paper/paper.html` | 修订稿自包含网页版 |
| `paper/main.tex` | LaTeX 主文件，可直接编译（`paper/compile.sh`） |
| `paper/sections/` | 各节源文件 |
| `paper/tables/` | 正文与附录引用的脚本生成表格 |
| `paper/figures/` | 正文与附录引用的脚本生成图（fig_01–fig_10） |
| `paper/bibliography.bib` | 参考文献（69 条，全部在正文中引用） |
| `review/revision_response_r2.md` | 第二轮逐条答复（N1–N9、M4/M5/M7 余项、次要意见、问题 1–7；第 7 节为交付步骤补充，记录本轮关闭的事项） |
| `review/peer_review_revised.md` | 第二轮复审意见（6/10，未改动） |
| `review/revision_response.md` | 第一轮逐条答复（M1–M7，次要 1–14，未改动） |
| `review/peer_review.md` | 第一轮审稿意见（英文，4/10，未改动） |
| `experiment/code/run_all.sh` | 全流程重现入口（含第一轮 8 个与第二轮 9 个修订分析步骤） |
| `experiment/code/` | 全部分析脚本；修订新增的以 rev_ 开头；含 china 的文件为 legacy，不被当前稿件引用 |
| `experiment/analysis/analysis_plan.md` | 分析计划（2026-11-05 05:42:50 UTC 锁定，写于探索之后，含修正案 1） |
| `experiment/analysis/analysis_plan_amendment2.md` | 修正案 2（10:16:29 UTC 锁定）；两个 .lock 文件记录哈希与时间 |
| `experiment/analysis/results.json` | 全部数值结果；第一轮在 `revision_2026_11` 键，第二轮在 `revision_2026_11_r2` 键 |
| `experiment/data/` | 数据与数据清单（含 SHA-256） |
| `experiment/logs/` | 干净副本全流程运行日志 |
| `.work/rn-targets-peer-review/rereview2/repro_checks_r2.md` | 第三轮独立核验记录（结果表、已更正错误、剩余事项检查、同计数比较、未重跑部分） |
| `.work/rn-targets-peer-review/rereview2/` | 第三轮核验脚本与输出（E1 与第二轮结果重算、两种夹心下的 bootstrap-t、三分位拟合、同计数比较、剩余事项检查，各脚本与同名 JSON 输出） |
| `.work/rn-targets-peer-review/rereview/repro_checks.md` | 第二轮复审的核验记录 |
| `.work/rn-targets-peer-review/revision2/number_changes.md` | 第三轮修订稿的新旧数字对照 |
| `.work/rn-targets-peer-review/paper/qa/claim-evidence-ledger.json` | 论断与证据对照表 |
| `.work/rn-targets-peer-review/paper/qa/scholarly-review.json` | 学术文风检查的逐条处理说明 |

## 结论溯源

修订稿的数字来源是 `experiment/analysis/results.json`；审稿方核验数字来源是 `.work/rn-targets-peer-review/rereview2/` 下的 JSON 文件。

| 结论 | 数字来源 | 生成脚本 | 位置 |
|---|---|---|---|
| E1：a = 0.0150，p = 5.1×10^-9 | `results.json` → `revision_2026_11.E1_intercept` | `rev_intercept.py` | 稿件第 5.1 节，表 5 |
| 截距随人口下降，P90/P10 比 0.50 | `revision_2026_11_r2.B_floor` | `rev_r2_floor.py` | 稿件第 5 节 |
| 三参数模型 a = 0.0089，γ = 0.617 | `revision_2026_11_r2.B_floor` | `rev_r2_floor.py` | 稿件第 5 节，附录 |
| 无模型下限 10.5%（净 10.1%） | `revision_2026_11_r2.B_floor` | `rev_r2_floor.py` | 稿件第 5 节，表 S10 |
| 生成器 0/6 通过、偏差 6–22% | `revision_2026_11_r2.D1_state_realistic_generators`、`C_total_bias` | `rev_r2_generator.py`、`rev_r2_total_bias.py` | 稿件第 4.5、5.6 节，表 4 |
| 钳位比例 13.5%、自助 SD 比、泊松参照 | `revision_2026_11_r2.E_clamp`、`F1_bootstrap_sd`、`F2_poisson_reference` | `rev_r2_clamp.py`、`rev_r2_bootstrap_sd.py`、`rev_r2_poisson_ref.py` | 稿件第 5 节，附录 C |
| 季节检验 p = 0.102（探索性） | `revision_2026_11.season_test` | `rev_season.py` | 稿件第 5.5 节 |
| 审稿核验：E1、人口、三参数、十分位、下限、钳位、锁定时序 | `repro_checks_r2.json` | `repro_checks_r2.py` | `peer_review_revised_r2.md` 第 4.2 节 |
| 审稿核验：bootstrap-t 与夹心 | `repro_clamp_diag.json` | `repro_clamp_diag.py` | `peer_review_revised_r2.md` N13 |
| 审稿核验：三分位拟合 | `repro_tertile_check.json` | `repro_tertile_check.py` | `repro_checks_r2.md` 第 1 节 |
| 审稿补充：同计数人口梯度 | `repro_matched_count.json` | `repro_matched_count.py` | `peer_review_revised_r2.md` N10 |
| 剩余事项的机械检查 | `check_residuals_r2.json` | `check_residuals_r2.py` | `repro_checks_r2.md` 第 3 节 |

## 过程时间线

时间均为 UTC。

| 时间 | 步骤 | 状态 |
|---|---|---|
| 2026-11-05 03:57:15 | 审稿阶段开始：读取需求与计划 | 已完成 |
| 2026-11-05 03:57:35 | 克隆仓库（`main` @ `2f5c407`），通读稿件与关键脚本 | 已完成 |
| 2026-11-05 04:05:40 | 由随附数据独立重建主样本和仿射拟合 | 已验证（精确重现） |
| 2026-11-05 04:10:19 | 季节检验的替代校准（周野自助、日历周置换，B = 4,999） | 已验证（重抽样检验） |
| 2026-11-05 04:16:06 | 流感住院的有效报告概率 ρ ≈ 0.001–0.01 | 未验证的假设（依据已发表负担估计的粗算） |
| 2026-11-05 04:31:28 | 写成 `peer_review.md`（大修，4/10） | 已完成 |
| 2026-11-05 05:42:50 | 写定并锁定分析计划（写于探索之后，第 0 节披露已知结果） | 已完成 |
| 2026-11-05 05:59:18 | E1 完整运行写出结果（早于 06:02:26 写入的修正案 1） | 已验证（比对日志时间与文件修改时间） |
| 2026-11-05 06:05:00 | 依计划运行 8 个修订分析步骤，E1 判定“稳健”，季节检验不显著并降级 | 已完成 |
| 2026-11-05 06:38:04 | 干净副本全流程 `run_all.sh` 完成，无失败 | 已验证（结果与工作副本逐字节一致） |
| 2026-11-05 08:26:27 | 35 条缺少标识的文献通过 Crossref DOI 记录核实 | 已验证（标题相似度 ≥ 0.9 且年份相差不超过 1） |
| 2026-11-05 08:29:44 | 发布第一轮修订稿与实验目录 | 已验证（哈希比对） |
| 2026-11-05 08:30:52 | 写成 `revision_response.md` | 已完成 |
| 2026-11-05 08:57:38 | 第二轮复审：用随稿函数与独立实现重算 E1、计数分布、替代函数形式；季节置换抽查（B = 999） | 已验证（与 `results.json` 逐项比对） |
| 2026-11-05 09:04:31 | 写成 `peer_review_revised.md`（6/10） | 已完成 |
| 2026-11-05 09:13:41 | 第二轮复审交付总检查（`--lang zh`）退出 0 | 已验证（交付总检查） |
| 2026-11-05 10:16:29 | 写定并锁定修正案 2（SHA-256 `8f00e540…ae8d`，22,263 字节）；锁定前未运行任何新分析 | 已验证（锁文件哈希重算；最早的新分析日志在 10:20:41） |
| 2026-11-05 10:24:25 | N2：截距随州人口下降；三参数形式下截距的周整块自助区间含 0；P90 处无模型下限 10.1% | 已完成 |
| 2026-11-05 10:27:39 | 结构化生成器：6 个设定中 0 个通过全部四项检验（结果不利于原主线，如实报告并改写主线） | 已完成 |
| 2026-11-05 10:52:02 | 外部报告数据：医院覆盖率调整与 Delphi 档案回填（661 行） | 已完成 |
| 2026-11-05 11:14:59 | 干净副本全流程 `run_all.sh`（含 9 个第二轮步骤）完成，无失败步骤 | 已验证（`tab_r2_*` 表逐字节一致） |
| 2026-11-05 11:21:46 | 第三轮改稿开始，按补充分析改写主线（该步在时限内未完成全部条目） | 已完成（部分） |
| 2026-11-05 12:12:07 | 第三轮修订稿编译：42 页，无未定义引用，69 条文献全部被引用 | 已验证（编译日志与 `check_manuscript.py` 退出 0） |
| 2026-11-05 12:12:45 | 发布第三轮修订稿与实验目录 | 已验证（sha256 比对） |
| 2026-11-05 12:28:23 | 写成第二轮逐条答复 `revision_response_r2.md`，第 6 节列出 18 项未完成事项 | 已完成 |
| 2026-11-05 12:36:15 | 第三轮复审与交付步骤开始：读取复审意见、答复与第二轮结果摘要 | 已完成 |
| 2026-11-05 12:38:48 | 表格生成脚本只改措辞（pre-specified、合计、表 15 日期）后重新生成，数值未变 | 已完成 |
| 2026-11-05 12:44:31 | 补发实验目录：分析计划、修正案 2、两个锁文件与含第二轮结果的 `results.json` | 已验证（`check_residuals_r2.json` 目录清单） |
| 2026-11-05 12:49:02 | 独立核验：E1、人口、三参数（周自助 B = 499）、十分位、P90 下限、偏差与 pigeonhole 读回、钳位、锁定时序 | 已验证（与稿件和 `results.json` 逐项比对） |
| 2026-11-05 12:49:29 | 核验发现附录 C.3 面板数与相对标准误写错，更正为 200 个面板、约 5% | 已完成 |
| 2026-11-05 12:51:39 | 同一组野自助抽样比较两种夹心下的 bootstrap-t，作者的 p 值复现 | 已验证（`repro_clamp_diag.json`） |
| 2026-11-05 12:52:13 | 重新编译（42 页）、严格稿件检查与 `check_manuscript.py` 退出 0，重新发布论文 | 已验证（发布 PDF 与编译稿 sha256 一致） |
| 2026-11-05 12:52:36 | 对发布 PDF 与目录机械检查第二轮答复第 6 节剩余事项 | 已验证（`check_residuals_r2.json`） |
| 2026-11-05 12:54:03 | 人口三分位分别拟合（三分位内重算 x̄），与作者完全一致 | 已验证（`repro_tertile_check.json`） |
| 2026-11-05 12:56:06 | 审稿方补充检查：同计数段内的人口梯度（州自助 B = 999） | 已完成 |
| 2026-11-05 12:58:40 | 写成第三轮复审意见 `peer_review_revised_r2.md`（7/10） | 已完成 |
| 2026-11-05 12:58:56 | `revision_response_r2.md` 追加第 7 节交付步骤补充 | 已完成 |
| 2026-11-05 12:59:00 | 州结构生成器、MC SD、偏差模拟、泊松参照与回填分析未从头重跑 | 未验证的假设（以作者日志与 `results.json` 为准） |
| 2026-11-05 13:08:30 | 复审意见、修订答复与本 README 以 `--lang zh` 运行交付总检查，修正上标字符、旧式 DOI 链接与文件导航后退出 0 | 已验证（10 项检查全部通过） |

## 产物缺陷

- 作者署名仍是占位文字 "[Author names and affiliations to be added before submission]"；CRediT 只列贡献类别。原仓库没有作者信息，本任务不编造作者。
- 参考文献表不输出 DOI（文献样式 `unsrtnat-noid`）；标识统一列在下方出处表，投稿时按目标期刊样式补齐。
- 计划、修正案与锁文件没有存入带时间戳的公共仓库，锁定时间只有本地文件时间与校验和自证。
- 中国附录的 legacy 代码、数据与日志仍在 `experiment/` 中（N15）；它们不被当前稿件引用，删除会破坏 `run_all.sh` 对旧审计输出的重现，因此保留并标注。
- 稿件 42 页、16 张正文表、摘要 272 词，偏长（N14）；本步只修措辞与事实错误，没有压缩篇幅。
- 上一级 output 文件夹的总结报告 final-report 已在最终报告步骤更新到第三轮（7/10）的状态；它的四张图由已发布的结果文件重新绘制，没有新增分析。
- 中文交付物中的数学符号用字体可显示的写法（a1、5.12×10^-9 等），原因是交付字体缺少上下标与带波浪号字母的字形。

## 环境与工具问题

- **字体替代。** 运行环境没有 Times New Roman，PDF 正文用度量兼容的 Liberation Serif，数学字体是 STIX；版面与字号不受影响。网页版内嵌 Latin Modern Math，避免公式在无数学字体的设备上显示为方框。
- **Crossref 检索限流。** 检索接口连续返回 HTTP 429，改为按候选 DOI 直接查询 works 记录，35 条全部核实。
- **英文论文的交付检查有两项结构性不通过。** 对 `paper/paper.pdf` 与 `paper/paper.html` 以 `--lang en` 运行交付总检查时，字形检查用中文字体核对英文网页的数学符号，且会顺带以英文核对本中文 README；两项都不是稿件缺陷，未为通过检查而改动内容。
- **总结报告的渲染检查在字体预检阶段退出。** 对 `final-report.docx` 运行交付总检查时，其余各项（样式来源、结构审计、语言、字形、引用、图片）都通过；渲染检查因运行环境把 Songti SC 固定映射到 Noto Serif CJK SC、把 Times New Roman 映射到 Liberation Serif 而在预检时退出，没有进行分页检查。为弥补这一点，对 18 页 PDF 逐页目视检查：无缺字、表格完整、图注与图同页、目录页码已解析，文末与 Markdown 结尾一致。
- **一条旧式 DOI 无法被自动核验。** Gelman 1999 的 DOI 是含尖括号的 Wiley SICI 格式，出处表中改给作者主页上的 PDF 链接。
- **第二轮补充分析的数据问题。** Delphi Epidata 部分记录数值为空，改为跳过空值记录；Delphi 的 NHSN 档案最早只到 2024 年第 47 周，回填分析只覆盖 661 行；PubMed 网页返回 Cookie 验证页，改用 NCBI E-utilities。
- **命令过滤。** 运行环境拦截含区间量词的 shell 命令（包括嵌在命令行里的 Python 正则），本轮的剩余事项检查因此写成脚本文件 `check_residuals_r2.py`，不影响结果。
- **图件检查的字形误报。** 附录 C.3 第一次改稿写入的平方根符号与图件检查的规则冲突，改用文字“square root of 3”。
- **出处表脚本缺记录。** 生成出处表时，文献元数据文件缺少 3 条第二轮新增文献（CDC 数据渠道说明、Census 人口估计、Delphi 开放数据库）的记录，补入后 69 条全部生成。
- **稿件检查路径。** `check_manuscript.py` 需要在实验代码目录下以论文目录为参数运行，第一次从错误目录运行报找不到文件，改正后退出 0。

## 研究局限与取舍

- **修订与审稿由同一系统完成，不是独立同行评议。** 三轮审稿意见与两轮修订都在同一任务中由同一系统完成。审稿尽量当作新投稿处理：独立重写拟合与夹心估计、使用不同随机种子、核对锁定时序、对照已发布 PDF 逐项检查；但它不能替代外部审稿人，评分可作参考，投稿后以期刊审稿为准。
- **7 分低于用户要求的 8 分，没有通过措辞或改分弥补。** 本轮按要求只用真实修改提高质量：关闭了全部不需要新分析的剩余事项，更正了核验发现的事实错误。阻碍 8 分的四个问题（N10、N11/M4、M5、N12）都需要新分析：在三参数模型中加入人口协变量、对 b/b_P 联合重抽参照，都必须先写入并锁定修正案 3 再运行。本步是复审与交付步骤，时间预算和工作单都不包括新一轮锁定分析，因此留作后续修改，路径写在 `peer_review_revised_r2.md` 第 5 节。
- **显著性约束。** 没有为求显著而改动 α、更换校准、增删生成器或删改数据；不利结果（生成器 0/6 通过、三参数下截距区间含 0、同计数比较只部分支持人口梯度）都如实报告。
- **作者方事项不计入评分。** 作者署名、CRediT、公共存档 DOI 与文献 DOI 需要真实作者与账户，不能生成，也未编造；它们不反映研究质量，所以单列在评分之外，但完成前稿件不具备投稿条件。
- **未重跑的长时计算。** 州结构生成器、重抽场生成器的 MC SD、偏差模拟、三种泊松参照与回填分析只做了读回核对，以作者日志与 `results.json` 为准。
- **锁定只有本地自证。** 修正案 2 的锁定时间来自本地文件时间与校验和，没有第三方时间戳。
- **期刊与报告规范按专业判断选定。** 未指定目标期刊，按统计方法类流行病学期刊的通用结构撰写；研究是方法加监测数据应用，不套用 CONSORT/PRISMA。投稿前请按目标期刊作者须知核对篇幅与格式（N14）。

## AI 使用声明草稿

以下段落可直接粘贴到投稿系统；如目标期刊对 AI 使用有具体政策，请先按该政策核对措辞。

> 本稿件的修订过程使用了 AI 辅助工具，用于审阅原稿、编写与运行分析脚本、核对参考文献标识以及润色英文表述。分析依照事先锁定的分析计划及其修正案执行（计划写于探索性分析之后，已在文中披露），所有数字均由随附脚本从公开数据重新计算，并经全流程重现核对。作者对稿件内容、分析选择与结论负全部责任。

## 出处表

编号与正文 [n] 一致。链接优先给 DOI;原文位置是精读时记下的节号或图表号。

| 编号 | 文献 | 链接 | 支撑的结论 | 原文位置 |
|---|---|---|---|---|
| [1] | Fraser, Christophe. 2007. Estimating Individual and Household Reproduction Numbers in an Emerging Epidemic. PLoS ONE | https://doi.org/10.1371/journal.pone.0000758 | Reproduction numbers summarise transmission under specified conditions and are routinely estimated from surveillance counts. | 未定位 |
| [2] | Gostic, Katelyn M. 等. 2020. Practical considerations for measuring the effective reproductive number, Rt. PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1008409 | A reporting level that is constant over consecutive weeks cancels in the ratio, as it does in renewal estimators. | 未定位 |
| [3] | Cori, Anne 等. 2013. A New Framework and Software to Estimate Time-Varying Reproduction Numbers During Epidemics. American Journal of Epidemiology | https://doi.org/10.1093/aje/kwt133 | It is not a true reproduction number or a posterior estimate from a renewal model. | 未定位 |
| [4] | Thompson, R.N. 等. 2019. Improved inference of time-varying reproduction numbers during infectious disease outbreaks. Epidemics | https://doi.org/10.1016/j.epidem.2019.100356 | Renewal methods relate incidence to past infectiousness through a generation-interval distribution, and the width of their credible intervals shrinks as the number of cases grows. | 未定位 |
| [5] | Parag, Kris V.. 2021. Improved estimation of time-varying reproduction numbers at low case incidence and between epidemic waves. PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1009347 | Work on low incidence, reporting noise and timely detection has made the limits of period-specific inference explicit, and temporal aggregation and reporting delays separate the reproduction number of the transmission ... | 未定位 |
| [6] | Parag, Kris V. 等. 2022. Quantifying the information in noisy epidemic curves. Nature Computational Science | https://doi.org/10.1038/s43588-022-00313-1 | Information from other features of the epidemic curve, independent reporting measurements or prior structure changes what can be learned, and estimators that accept bias in exchange for lower variance lie outside that ... | 未定位 |
| [7] | Parag, Kris V. 等. 2025. Asymmetric limits on timely interventions from noisy epidemic data. Communications Physics | https://doi.org/10.1038/s42005-025-02358-w | Work on low incidence, reporting noise and timely detection has made the limits of period-specific inference explicit, and temporal aggregation and reporting delays separate the reproduction number of the transmission ... | 未定位 |
| [8] | Nash, Rebecca K. 等. 2023. Estimating the epidemic reproduction number from temporally aggregated incidence data: A statistical modelling approach and software tool. PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1011439 | Work on low incidence, reporting noise and timely detection has made the limits of period-specific inference explicit, and temporal aggregation and reporting delays separate the reproduction number of the transmission ... | 未定位 |
| [9] | Ogi-Gittins, I. 等. 2024. A simulation-based approach for estimating the time-dependent reproduction number from temporally aggregated disease incidence time series data. Epidemics | https://doi.org/10.1016/j.epidem.2024.100773 | Work on low incidence, reporting noise and timely detection has made the limits of period-specific inference explicit, and temporal aggregation and reporting delays separate the reproduction number of the transmission ... | 未定位 |
| [10] | Eales, Oliver 等. 2024. Differences between the true reproduction number and the apparent reproduction number of an epidemic time series. Epidemics | https://doi.org/10.1016/j.epidem.2024.100742 | It is not a true reproduction number or a posterior estimate from a renewal model. | 未定位 |
| [11] | Searle, Shayle R. 等. 1992. Variance Components. John Wiley & Sons | https://openlibrary.org/works/OL16956832W | Observing the effects themselves gives at least as much information as the counts, and their sample variance is the minimum-variance unbiased estimator. | 未定位 |
| [12] | Cochran, William G.. 1977. Sampling Techniques. John Wiley & Sons | https://openlibrary.org/works/OL1351802W | The optimum in Eq. has the form of cluster-size allocation in survey sampling and multilevel design: cases per period should grow only with the square root of and of the relative cost of a period. | 未定位 |
| [13] | DerSimonian, Rebecca 等. 1986. Meta-analysis in clinical trials. Controlled Clinical Trials | https://doi.org/10.1016/0197-2456(86)90046-2 | Random-effects meta-analysis estimates the between-study component on top of within-study error, and small-area estimation and empirical-Bayes disease mapping weight each direct estimate against a common mean according ... | 未定位 |
| [14] | Fay, Robert E. 等. 1979. Estimates of Income for Small Places: An Application of James-Stein Procedures to Census Data. Journal of the American Statistical Association | https://doi.org/10.1080/01621459.1979.10482505 | The count at which counting error equals between-period variance is the half-weight point of the corresponding shrinkage estimator, as in small-area estimation and disease mapping. | 未定位 |
| [15] | Clayton, David 等. 1987. Empirical Bayes Estimates of Age-Standardized Relative Risks for Use in Disease Mapping. Biometrics | https://doi.org/10.2307/2532003 | The count at which counting error equals between-period variance is the half-weight point of the corresponding shrinkage estimator, as in small-area estimation and disease mapping. | 未定位 |
| [16] | Gelman, Andrew 等. 1999. All maps of parameter estimates are misleading. Statistics in Medicine | https://sites.stat.columbia.edu/gelman/research/published/allmaps.pdf | Random-effects meta-analysis estimates the between-study component on top of within-study error, and small-area estimation and empirical-Bayes disease mapping weight each direct estimate against a common mean according ... | 未定位 |
| [17] | Spiegelhalter, David J.. 2005. Funnel plots for comparing institutional performance. Statistics in Medicine | https://doi.org/10.1002/sim.1970 | It has the form of the additive overdispersion model of funnel plots. | 未定位 |
| [18] | Spiegelhalter, David 等. 2012. Statistical Methods for Healthcare Regulation: Rating, Screening and Surveillance. Journal of the Royal Statistical Society: Series A (Statistics in Society) | https://doi.org/10.1111/j.1467-985X.2011.01010.x | It has the form of the additive overdispersion model of funnel plots. | 未定位 |
| [19] | Engen, Steinar 等. 1998. Demographic and Environmental Stochasticity-Concepts and Definitions. Biometrics | https://doi.org/10.2307/2533838 | Population ecology makes the same split between demographic and environmental stochasticity, as do stochastic epidemic models. | 未定位 |
| [20] | Lande, Russell 等. 2003. Stochastic Population Dynamics in Ecology and Conservation. Oxford University Press | https://openlibrary.org/works/OL8013358W | Population ecology makes the same split between demographic and environmental stochasticity, as do stochastic epidemic models. | 未定位 |
| [21] | Bretó, Carles 等. 2009. Time series analysis via mechanistic models. The Annals of Applied Statistics | https://doi.org/10.1214/08-aoas201 | Population ecology makes the same split between demographic and environmental stochasticity, as do stochastic epidemic models. | 未定位 |
| [22] | He, Daihai 等. 2010. Plug-and-play inference for disease dynamics: measles in large and small populations as a case study. Journal of the Royal Society Interface | https://doi.org/10.1098/rsif.2009.0151 | Population ecology makes the same split between demographic and environmental stochasticity, as do stochastic epidemic models. | 未定位 |
| [23] | Casella, George 等. 2002. Statistical Inference. Duxbury | https://openlibrary.org/works/OL4791041W | The mean squared error decompositions follow from the law of total variance. | 未定位 |
| [24] | Lehmann, Erich L. 等. 1998. Theory of Point Estimation. Springer | https://openlibrary.org/works/OL1869108W | The support does not depend on, so the regular one-parameter Cram'er-Rao inequality applies to estimators unbiased for every and gives the attained bound stated in Section. | 未定位 |
| [25] | Lloyd-Smith, J. O. 等. 2005. Superspreading and the effect of individual variation on disease emergence. Nature | https://doi.org/10.1038/nature04153 | Offspring overdispersion follows, and its estimation from small samples. | 未定位 |
| [26] | Lloyd-Smith, James O.. 2007. Maximum likelihood estimation of the negative binomial dispersion parameter for highly overdispersed data, with applications to infectious diseases. PLoS ONE | https://doi.org/10.1371/journal.pone.0000180 | Offspring overdispersion follows, and its estimation from small samples. | 未定位 |
| [27] | Taylor, L. R.. 1961. Aggregation, variance and the mean. Nature | https://doi.org/10.1038/189732a0 | A variance that scales as a power of the mean is Taylor's law. | 未定位 |
| [28] | Graubard, Barry I. 等. 2002. Inference for superpopulation parameters using sample surveys. Statistical Science | https://doi.org/10.1214/ss/1023798999 | The distinction between the realised quantity of a finite population and the parameter of the process that generated it follows; the case reproduction number of is a realised quantity of this kind. | 未定位 |
| [29] | Wallinga, Jacco 等. 2004. Different Epidemic Curves for Severe Acute Respiratory Syndrome Reveal Similar Impacts of Control Measures. American Journal of Epidemiology | https://doi.org/10.1093/aje/kwh255 | The distinction between the realised quantity of a finite population and the parameter of the process that generated it follows; the case reproduction number of is a realised quantity of this kind. | 未定位 |
| [30] | Parag, Kris V 等. 2020. Adaptive Estimation for Epidemic Renewal and Phylogenetic Skyline Models. Systematic Biology | https://doi.org/10.1093/sysbio/syaa035 | Choosing the estimation window to balance counting noise against changes in transmission is the aim of adaptive windowing, and hierarchical renewal models pool an observation layer across regions. | 未定位 |
| [31] | Bhatt, Samir 等. 2023. Semi-mechanistic Bayesian modelling of COVID-19 with renewal processes. Journal of the Royal Statistical Society: Series A (Statistics in Society) | https://doi.org/10.1093/jrsssa/qnad030 | Choosing the estimation window to balance counting noise against changes in transmission is the aim of adaptive windowing, and hierarchical renewal models pool an observation layer across regions. | 未定位 |
| [32] | Azmon, Amin 等. 2014. On the estimation of the reproduction number based on misreported epidemic data. Statistics in Medicine | https://doi.org/10.1002/sim.6015 | Reporting probabilities may change for all infections in the same period, and separating process from observation variation in time series is a known difficulty. | 未定位 |
| [33] | Stoner, Oliver 等. 2019. A Hierarchical Framework for Correcting Under-Reporting in Count Data. Journal of the American Statistical Association | https://doi.org/10.1080/01621459.2019.1573732 | Reporting probabilities may change for all infections in the same period, and separating process from observation variation in time series is a known difficulty. | 未定位 |
| [34] | Knape, Jonas. 2008. Estimability of density dependence in models of time series data. Ecology | https://doi.org/10.1890/08-0071.1 | Temporal structure can separate persistent transmission variation from independent counting and reporting noise under explicit state-space or renewal models, and external measurements of reporting can constrain and its ... | 未定位 |
| [35] | Auger-M'eth'e, Marie 等. 2016. State-space models' dirty little secrets: even simple linear Gaussian models can have estimation problems. Scientific Reports | https://doi.org/10.1038/srep26677 | Temporal structure can separate persistent transmission variation from independent counting and reporting noise under explicit state-space or renewal models, and external measurements of reporting can constrain and its ... | 未定位 |
| [36] | Burghardt, Keith 等. 2022. Unequal impact and spatial aggregation distort COVID-19 growth rates. Philosophical Transactions of the Royal Society A | https://doi.org/10.1098/rsta.2021.0122 | The cross-sectional spread of growth-rate or reproduction-number estimates has been read as variation in transmission or behaviour between places and periods without subtracting the part of that spread that depends on ... | 未定位 |
| [37] | Biggerstaff, Matthew 等. 2014. Estimates of the reproduction number for seasonal, pandemic, and zoonotic influenza: a systematic review of the literature. BMC Infectious Diseases | https://doi.org/10.1186/1471-2334-14-480 | The cross-sectional spread of growth-rate or reproduction-number estimates has been read as variation in transmission or behaviour between places and periods without subtracting the part of that spread that depends on ... | 未定位 |
| [38] | Ives, Anthony R. 等. 2021. Estimating and explaining the spread of COVID-19 at the county level in the USA. Communications Biology | https://doi.org/10.1038/s42003-020-01609-6 | The cross-sectional spread of growth-rate or reproduction-number estimates has been read as variation in transmission or behaviour between places and periods without subtracting the part of that spread that depends on ... | 未定位 |
| [39] | Bhatia, Rajendra 等. 2000. A better bound on the variance. The American Mathematical Monthly | https://doi.org/10.1080/00029890.2000.12005203 | Because lies in with mean,, so. | 未定位 |
| [40] | Centers for Disease Control 等. 2026. Weekly Hospital Respiratory Data (HRD) Metrics by Jurisdiction, National Healthcare Safety Network (NHSN) | https://data.cdc.gov/resource/ua7e-t2fy | Centers for Disease Control and Prevention under the identifier ua7e-t2fy, in the snapshot retrieved on 20 October 2026. | 未定位 |
| [41] | Frutos, Aaron M. 等. 2026. Evaluation of a Novel Data Source for National Influenza Surveillance: Influenza Hospitalization Data in the National Healthcare Safety Network, United States, September 2021–April 2024. American Journal of Public Health | https://doi.org/10.2105/AJPH.2025.308282 | We used the National Healthcare Safety Network (NHSN) weekly hospital respiratory-data snapshot retrieved on 20 October 2026 and restricted influenza admissions to the 50 US states and the District of Columbia. | 未定位 |
| [42] | Centers for Medicare & Medicaid Services. 2024. QSO-25-05-Hospitals/CAHs: Updates to the Condition of Participation (CoP) Requirements for Hospitals and Critical Access Hospitals (CAHs) To Report Acute Respiratory Illnesses | https://www.cms.gov/files/document/qso-25-05-hospitals-cahs.pdf | Reporting was voluntary from May to October 2024 before it became mandatory again, and the national data channel was interrupted between early October and mid-November 2025. | 未定位 |
| [43] | Centers for Disease Control and Prevention. 2026. Respiratory Illnesses Data Channel: Release Notes and Frequently Asked Questions. CDC web page | https://www.cdc.gov/respiratory-viruses/data/faqs.html | The national Respiratory Illnesses Data Channel was not updated from 3 October to 13 November 2025. | 未定位 |
| [44] | U.S. Census Bureau. 2026. Annual Estimates of the Resident Population for the United States, Regions, States, District of Columbia, and Puerto Rico: April 1, 2020 to July 1, 2025 (NST-EST2025-ALLDATA), Vintage 2025. Census Bureau Population Estimates Program | https://www2.census.gov/programs-surveys/popest/datasets/2020-2025/state/totals/NST-EST2025-ALLDATA.csv | State resident populations from the Census Bureau Vintage 2025 estimates define the population covariate of the intercept. | 未定位 |
| [45] | Reinhart, Alex 等. 2021. An open repository of real-time COVID-19 indicators. Proceedings of the National Academy of Sciences | https://doi.org/10.1073/pnas.2111452118 | The Delphi Epidata archive records the first published and every revised value of each state-week NHSN count from November 2024 on; 661 rows had a provisional first issue. | 未定位 |
| [46] | Wallinga, J 等. 2007. How generation intervals shape the relationship between growth rates and reproductive numbers. Proceedings of the Royal Society B: Biological Sciences | https://doi.org/10.1098/rspb.2006.3754 | With the admissions in week, the apparent reproduction number of a state-week is the growth-to-reproduction transform: equation R_A,t=[0,,1+(/7)(C_t+1+0.5)/(C_t+0.5)]^,=(g_/g_)^2,=g_^2/g_. eq:ra equation The transform ... | 未定位 |
| [47] | Cowling, Benjamin J. 等. 2009. Estimation of the Serial Interval of Influenza. Epidemiology | https://doi.org/10.1097/EDE.0b013e31819d1092 | The primary interval has mean days and standard deviation days, the household serial-interval estimate of. | 未定位 |
| [48] | Vink, Margaretha Annelie 等. 2014. Serial Intervals of Respiratory Infectious Diseases: A Systematic Review and Analysis. American Journal of Epidemiology | https://doi.org/10.1093/aje/kwu209 | A serial interval is not a generation interval. | 未定位 |
| [49] | Chan, Louis Yat Hin 等. 2025. Estimating the generation time for influenza transmission using household data in the United States. Epidemics | https://doi.org/10.1016/j.epidem.2025.100815 | Sensitivity analyses set the mean to 3.2 and 4.0 days and the standard deviation to 1.2 and 2.0 days; the lower mean is close to the intrinsic generation-time estimate from US household data. | 未定位 |
| [50] | Cameron, A. Colin 等. 2011. Robust Inference With Multiway Clustering. Journal of Business & Economic Statistics | https://doi.org/10.1198/jbes.2010.07136 | The primary interval is the two-way cluster-robust sandwich of for the stacked estimating equations of and the affine fit, clustered by state and calendar week, with small-sample factors and critical values from the ... | 未定位 |
| [51] | Cameron, A. Colin 等. 2008. Bootstrap-Based Improvements for Inference with Clustered Errors. Review of Economics and Statistics | https://doi.org/10.1162/rest.90.3.414 | Two resampling checks accompany it: a calendar-week block bootstrap (1,999 draws) that resamples whole weeks with all their state rows, and a wild cluster bootstrap- over calendar weeks with six-point weights and the ... | 未定位 |
| [52] | Webb, Matthew D.. 2023. Reworking wild bootstrap-based inference for clustered errors. Canadian Journal of Economics | https://doi.org/10.1111/caje.12661 | Two resampling checks accompany it: a calendar-week block bootstrap (1,999 draws) that resamples whole weeks with all their state rows, and a wild cluster bootstrap- over calendar weeks with six-point weights and the ... | 未定位 |
| [53] | Owen, Art B.. 2007. The pigeonhole bootstrap. The Annals of Applied Statistics | https://doi.org/10.1214/07-AOAS122 | The pigeonhole bootstrap overestimates the variance of a mean in a crossed design by a bounded factor. | 未定位 |
| [54] | MacKinnon, James G. 等. 2018. The wild bootstrap for few (treated) clusters. The Econometrics Journal | https://doi.org/10.1111/ectj.12107 | With few clusters at the level of assignment, wild cluster bootstrap tests can be badly mis-sized and randomisation inference at that level is recommended. | 未定位 |
| [55] | MacKinnon, James G. 等. 2023. Cluster-robust inference: A guide to empirical practice. Journal of Econometrics | https://doi.org/10.1016/j.jeconom.2022.04.001 | With few clusters at the level of assignment, wild cluster bootstrap tests can be badly mis-sized and randomisation inference at that level is recommended. | 未定位 |
| [56] | Reed, Carrie 等. 2015. Estimating Influenza Disease Burden from Population-Based Surveillance Data in the United States. PLOS ONE | https://doi.org/10.1371/journal.pone.0118369 | The range -0.01 for the probability that an infection becomes a reported admission is assumed from burden estimates, not estimated. SC4 and SC5, interval validity: these use generators that redraw the heterogeneity ... | 未定位 |
| [57] | Hayward, Andrew C 等. 2014. Comparative community burden and severity of seasonal and pandemic influenza: results of the Flu Watch cohort study. The Lancet Respiratory Medicine | https://doi.org/10.1016/S2213-2600(14)70034-7 | The range -0.01 for the probability that an infection becomes a reported admission is assumed from burden estimates, not estimated. SC4 and SC5, interval validity: these use generators that redraw the heterogeneity ... | 未定位 |
| [58] | Parag, Kris V. 等. 2022. Fundamental limits on inferring epidemic resurgence in real time using effective reproduction numbers. PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1010004 | Information from other features of the epidemic curve, independent reporting measurements or prior structure changes what can be learned, and estimators that accept bias in exchange for lower variance lie outside that ... | 未定位 |
| [59] | Höhle, Michael 等. 2014. Bayesian nowcasting during the STEC O104:H4 outbreak in Germany, 2011. Biometrics | https://doi.org/10.1111/biom.12194 | In real time, reporting delays leave the most recent weeks incomplete, and right-censoring biases the latest counts downwards in proportion to the delay distribution. | 未定位 |
| [60] | McGough, Sarah F. 等. 2020. Nowcasting by Bayesian Smoothing: A flexible, generalizable model for real-time epidemic tracking. PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1007735 | In real time, reporting delays leave the most recent weeks incomplete, and right-censoring biases the latest counts downwards in proportion to the delay distribution. | 未定位 |
| [61] | Abbott, Sam 等. 2020. Estimating the time-varying reproduction number of SARS-CoV-2 using national and subnational case counts. Wellcome Open Research | https://doi.org/10.12688/wellcomeopenres.16006.1 | Real-time renewal methods address this by modelling delays and the infection process jointly. | 未定位 |
| [62] | Lison, Adrian 等. 2024. Generative Bayesian modeling to nowcast the effective reproduction number from line list data with missing symptom onset dates. PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1012021 | Real-time renewal methods address this by modelling delays and the infection process jointly. | 未定位 |
| [63] | Pitzer, Virginia E 等. 2021. The Impact of Changes in Diagnostic Testing Practices on Estimates of COVID-19 Transmission in the United States. American Journal of Epidemiology | https://doi.org/10.1093/aje/kwab089 | In practice testing, care seeking and admission thresholds can change with epidemic activity. | 未定位 |
| [64] | Sherratt, Katharine 等. 2021. Exploring surveillance data biases when estimating the reproduction number: with insights into subpopulation transmission of COVID-19 in England. Philosophical Transactions of the Royal Society B | https://doi.org/10.1098/rstb.2020.0283 | In practice testing, care seeking and admission thresholds can change with epidemic activity. | 未定位 |
| [65] | Koyama, Shinsuke 等. 2021. Estimating the time-varying reproduction number of COVID-19 with a state-space method. PLOS Computational Biology | https://doi.org/10.1371/journal.pcbi.1008679 | Temporal structure can separate persistent transmission variation from independent counting and reporting noise under explicit state-space or renewal models, and external measurements of reporting can constrain and its ... | 未定位 |
| [66] | Hamilton, James D.. 1994. Time Series Analysis. Princeton University Press | https://openlibrary.org/works/OL3906201W | Under the stated stationary autocovariance, summing all covariance pairs gives with the first expression in Eq., and the finite geometric sum yields the second. | 未定位 |
| [67] | Raudenbush, Stephen W.. 1997. Statistical analysis and optimal design for cluster randomized trials. Psychological Methods | https://doi.org/10.1037/1082-989x.2.2.173 | The optimum in Eq. has the form of cluster-size allocation in survey sampling and multilevel design: cases per period should grow only with the square root of and of the relative cost of a period. | 未定位 |
| [68] | Morris, Tim P. 等. 2019. Using simulation studies to evaluate statistical methods. Statistics in Medicine | https://doi.org/10.1002/sim.8086 | The audit is a simulation study with known truth, reported under the aims, data-generating mechanisms, estimands, methods and performance measures structure of. | 未定位 |
| [69] | Owen, Art B. 等. 2012. Bootstrapping data arrays of arbitrary order. The Annals of Applied Statistics | https://doi.org/10.1214/12-AOAS547 | The pigeonhole bootstrap overestimates the variance of a mean in a crossed design by a bounded factor. | 未定位 |
