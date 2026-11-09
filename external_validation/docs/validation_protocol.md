# PULSE 复现与外部验证实验协议（准备阶段）

## 任务与输入边界

原论文：Wu et al., *Nature Computational Science* **6**, 954–970 (2026), DOI `10.1038/s43588-026-01026-5`。

**原论文同数据复现**需要 UK Biobank 的原始受试者数据以及论文相同的训练、测试、模态、超参数和评价任务。公开 Su2020、HMP2、iPOP、COVAAR 和 CSIRO 均是**独立队列外部验证**，不能称作 UKB 原实验同数据复现。原始 `pulse/` 与原论文仓库既有 benchmark 程序全部保留，后续取得相应 UKB 权限时应建立另一个隔离的复现入口；不能把外部验证的权重、目标函数和性能替换成原论文数值。

当前可在本项目内准备和正式运行的是以下外部验证：

1. 预测当前 visit 的缺失模态，而非下一次尚未观测的 visit；`protein → metab` 和 `metab → protein` 各自单独训练。
2. 按**患者**划分训练、验证和测试（默认 70/15/15），同一患者全部访视只出现在其中一个划分中。
3. 两个任务分别报告：所有有目标观测的 `current_completion`，和具有至少一次既往可配对访视的 `historical_completion`。由于样本选择可能不同，不以两任务指标直接相减作为历史效益。
4. 纵向任务固定训练、验证、测试患者及当前目标，比较 `pulse_full`、`no_history_input`、`no_alignment`、`ridge_current` 和 `last_observation`。`no_history_input` 为重新训练的模型；`ridge_current` 只使用当前源模态，固定岭回归正则化强度 10，不利用测试目标；`last_observation` 直接使用目标模态之前的最近一次真实测量，不训练。
5. 所有特征选取、均值、方差均仅由训练患者数据确定。按特征分别计算测试患者上的 Pearson r、MAE、MSE，另报告有效特征数、有效单元格总数及相对上一次目标值的变化量表现。少于三个有效样本或零方差的相关性定义为缺失。
6. PULSE 外部验证当前使用 `target MSE + 0.25 context reconstruction + 0.001 KL + 0.01 Jeffrey alignment`。历史对齐消融将最后一项设为零。对抗损失关闭。这是固定的**外部验证训练协议**，不冒称论文参数完全一致。
7. 正式结论至少使用多随机种子、逐特征配对差异与患者层面不确定性分析。若对比特征集合因不同划分变化，必须在相同特征和同一测试患者上重算；不得直接比较互不对应的平均值。

## 数据集适用性和单独约束

| 来源 | 目前确定的事实 | 后续处理 |
|---|---|---|
| Su2020 | 现有 252 个同次双组学样本，既有真实实验已跑过 | 保留输出；正式统一协议运行时另外建新结果目录 |
| HMP2 | 308 个双组学精确样本、82 名患者，74 名患者至少有两次配对访视 | 已生成 canonical CSV；直接使用这三份真实文件 |
| iPOP seasonal | RData 内含组学矩阵与注释，但此前发现矩阵行索引与注释行索引不一致 | 先分析真实 RData 的索引语义，拿到可验证的 specimen crosswalk 才导出，不按行号配对 |
| COVAAR | 已取得 Figshare 代谢组与 Sample_Meta_Data.csv、GEO SomaScan 与 series matrix；两种组学样本列分别采用 Sample_... 与 VA... 标识 | 运行 `cohort_join_audit` 检查跨来源原生 ID 和 GEO 样本注释；仅建立可证明的患者访视映射 |
| CSIRO | 当前已获得 005 蛋白质组编号交叉表、006 CCM 样本与受试者时间信息及 002–004 蛋白质组矩阵 | 使用 `prepare_csiro --inspect-only` 按精确编号检验唯一配对，成功后生成 canonical CSV |
| P100 | 公开补充附件不等于获批访问的个体级测量数据 | 保留受限状态，禁止从群体统计重建假的个体测量 |

参考原始发布页：
- COVAAR GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325178
- COVAAR metabolomics: https://doi.org/10.6084/m9.figshare.31685023.v1
- CSIRO vaccine cohort: https://data.csiro.au/collection/csiro:77992
- iPOP seasonal: https://doi.org/10.6084/m9.figshare.12376508

## 数据适配程序

`python -m external_validation.source_probe` 在本地真实原始文件上生成 `source_profile.json`，解析本地 CSV/TSV/Excel/ZIP 文件的结构；对 iPOP RData，检查各对象的尺寸、索引与注释交集。它**不生成伪造的样本配对**。

`python -m external_validation.prepare_tabular --mapping <verified_mapping.json> --inspect-only` 使用已经证实的 cohort-specific crosswalk 核查样本配对与患者/时间，再移除 `--inspect-only` 导出三个 canonical CSV。mapping 必须写明所用的发布源路径、每个模态的 sample ID 字段、明确的 measurement feature 列以及真实跨模态测量交叉表。以样本 ID 严格对应；`patient_id`/`visit_order` 仅定义患者和时间，不能替代对两组学是否同一 specimen 的证明。工具没有设置任何样本 ID 字符串推测或按行数自动匹配。

JSON mapping 的结构参见 `external_validation/tests/test_prepare_tabular.py` 中的 `setup_data()`。数据文件位置均相对 `/mnt/f/datasets/pulse-omics/`，不进入仓库；真实 crosswalk 来源须记录在正式研究报告。

## 并入 main 的条件

- 既有与新增的外部验证单元测试全部通过。
- `main` 与 `omics` 差异限定在 `external_validation/`，不触动原始 `pulse/`、`scripts/` 和项目配置。
- Su2020/HMP2 的 canonical 输入与元数据已完成核查；iPOP/COVAAR/CSIRO 各有实际 schema 证据和明确的**可纳入 / 不能纳入 / 需要额外来源**裁决。
- 正式训练和统计评价命令确定，数据和实验输出继续保存在 F 盘。

在这些条件完成前不合并、不把模拟测试相关性当作科学结果、不启动完整训练。

## 真实数据新发现后的专用准备程序

- `python -m external_validation.cohort_join_audit --datasets covaar csiro`：读取实际 CSV、GEO series matrix 和 CSI​RO Excel 的样本对应字段，分别生成 `join_audit.json`。报告只包含匹配数和原始字段类别，不输出完整患者标识。
- `python -m external_validation.prepare_csiro --inspect-only`：核对 005 的 `Sample IDs / proteomic#` 与 006 的 `Sample IDS (Barcode)`，将 `proteomic#` 与 002–004 的真实 `Vac_Plasma_N` 列号匹配，检查受试者和访视唯一性。严格匹配成功后，移除 `--inspect-only` 生成模型输入。
- iPOP RData 解析保留非致命序列化警告，若仍失败将显示具体解码异常；其矩阵索引与样本注释关联没有证明前不能导出。
- HMP2 BIOM 稀疏值导出为零后的矩阵需在 `preflight` 中检查全零、稀疏非零特征的数量，正式建模之前确定训练集特征筛选策略。
