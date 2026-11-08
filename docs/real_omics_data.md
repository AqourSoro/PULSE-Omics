# PULSE-Omics：真实纵向多组学数据接入与数据契约

研究阶段：先获取真实数据，完成内容、样本及可配对性研究，再决定合成数据体系。
日期：2026-11-08。仓库分支：omics。大数据根目录：`/mnt/f/datasets/pulse-omics/`。

## 1. 研究目标决定数据的选择

核心任务：对于同一受试者 i 在已知过去观测和当前部分观测的条件下，预测缺失组学及之后的真实观测：

- 同次访问的跨组学补全：已观测蛋白质组 P(i,t)，重建代谢组 M(i,t)，反向亦然。
- 有历史信息的条件补全：P(i,≤t), M(i,<t) → M(i,t)；检验历史状态编码器是否提高精度。
- 未来访问预测：P/M(i,≤t) → P/M(i,t+Δt)，需要真实采样时间间隔和对应的后续观测。
- 跨人群、疾病、时间分布变化：按人分组保持测试独立；不能把重复访问随机分入不同数据集。

数据基础单元：**受试者、一次实际采集的生物样本、一个模态的测量、一个被测特征**。模态不等同于一个 Excel 文件，表中多次测定或重复上机需要分别处理。

## 2. 真实公开资源与数据含义

| 数据集 | 生物组织/模态 | 纵向结构 | 数据地址与格式 | 用途/状态 |
|---|---|---|---|---|
| Su2020 | 血浆蛋白质组、代谢组，临床状态 | 139 人，265 次临床采样，126 人 T1/T2；259 个蛋白质组患者样本、254 个代谢组患者样本，252 个同次双组学配对样本（项目实测） | Mendeley v5 的 Table S1 XLSX，已解析到 `raw/su2020/` 和 `processed/su2020/`；蛋白质组 464 特征、代谢组 1050 特征 | 已真实获取且配对；两次都具有双组学的受试者数尚需另算 |
| Stanford iPOP / T2D | **宿主**血浆蛋白质组，以及血清/尿液代谢组；另有细胞因子等 | 106 名参与者、最长约 4 年、访视间隔不规则且事件期加密 | 公开 HMPDACC S3 的 `ihmp/t2d/host/proteome/`、`proteome_nonpride/`、`metabolome/`；具体文件清单由 S3 匿名枚举产生 | 血液纵向生理状态的首选拓展队列；需逐文件核实模态、体液、特征单位与受试者/时间标识 |
| HMP2 / IBDMDB | **肠道微生物**蛋白质组及粪便代谢组，炎症性肠病随访 | 132 人、约 1 年多次访问；不同组学的采样子集并不完全一致 | 官方元数据 CSV、蛋白质组 TSV.GZ、代谢组 BIOM.GZ（数据矩阵及元数据） | 结构丰富的独立纵向多组学任务，**不可视为血浆宿主组学** |
| NIH COVAAR | SomaScan 血液蛋白质组，疫苗干预的代谢组 | 16 名受试者，接种相关基线、35min、2h、6h、24h 和后续加强针 | GEO `GSE325178` 包括 SomaScan RFU 矩阵及 series matrix，Figshare article 31685023 发布代谢组文件；需按受试者、访视、采样相对时间匹配 | 短时干预引起的动态变化，适合应激/非平稳模拟实验 |
| P100 | 临床指标、血液蛋白质组、代谢组；补充微生物组等 | 108 人、约 9 个月，每 3 个月采样 | 论文公开 supplements（Europe PMC PMC5568837 supplementaryFiles），**完整个人级原始数据在 dbGaP phs001363 受控** | 先获取公开处理数据，检查是否包含可训练的个人级矩阵；不能预设 supplements 有完整个人访问记录 |
| CSIRO 2026 疫苗研究 | 同一血浆样本中蛋白质组、脂质组、代谢组 | 53 人、137 份样本；三种疫苗平台。Pfizer 组有额外 6 个月和加强针时间点 | CSIRO collection `csiro:77992`，提供集合下载 API。文件列表端点实际开放性仍待验证 | 第二个宿主血浆干预队列，适合跨干预验证 |

引用与数据入口：

- Su Y. et al. (Cell 2020): https://data.mendeley.com/datasets/tzydswhhb5/5
- Zhou W. et al. (Nature 2019): https://www.nature.com/articles/s41586-019-1236-x
- HMPDACC 2026 公共 S3 目录索引: https://www.ihmpdcc.org/data.php
- HMP2 下载总目录: https://ibdmdb.org/results
- HMP2 蛋白质组处理矩阵: https://ibdmdb.org/downloads/html/products_HMP2_Proteomics.html
- HMP2 代谢组数据: https://ibdmdb.org/downloads/html/products_MBX_2017-03-01.html
- NIH COVAAR SomaScan GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE325178
- NIH COVAAR 代谢组: https://doi.org/10.6084/m9.figshare.31685023.v1
- P100 论文（PMC）: https://pmc.ncbi.nlm.nih.gov/articles/PMC5568837/
- P100 受控数据: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001363.v1.p1
- CSIRO 研究: https://pmc.ncbi.nlm.nih.gov/articles/PMC13672565/
- CSIRO 集合: https://data.csiro.au/collection/csiro:77992

### 可复制的真实下载流程

代码仓库在 Linux 项目目录，大型数据直接写入 Windows F: 盘的 WSL 挂载目录。

```bash
cd ~/projects/PULSE-Omics
conda activate pulse-omics
git switch omics
git pull --ff-only
python -m pip install -r requirements-data.txt

# 自动下载已定位的公开对象（脚本会尝试下载而非只列清单）
python scripts/download_real_omics.py --datasets ipop hmp2 covaar p100

# CSIRO 文件清单 API 尚待独立核实，单独尝试并记录实际返回
python scripts/download_real_omics.py --datasets csiro

# 逐文件探查格式和头信息
python scripts/profile_real_omics.py
```

脚本特点：

- iPOP：Boto3 **匿名 S3** 枚举并下载文件，不需要 AWS 账户；按当前目录中真实出现的文件名与对象大小决定，限定处理数据格式，排除大量原始质谱与测序输出。
- HMP2：从官方资源链接直取元数据、一个蛋白质组汇总矩阵及带元数据的代谢组 BIOM 文件。
- COVAAR：从 GEO 下载 SomaScan 处理矩阵和样本矩阵；通过 Figshare 公共 API 枚举文章真实附件并下载。
- P100：从 Europe PMC 获取论文公开补充文件包。
- CSIRO：尝试官方集合 `files/list` API，只有确认返回实际的公开数据文件 URL 才下载，不把门户 HTML 误当组学数据。
- 统一记录原始来源、文件大小、SHA-256、下载状态；已存在文件跳过，可断点续传。默认下载限制单文件 2 GiB、一次运行 12 GiB，必要时可使用 `--max-file-gb` / `--max-total-gb` 调整，不修改默认存储位置。不会将私有身份认证信息写入仓库。
- 下载错误写入每个数据集的 `raw/<dataset>/download_manifest.json`，部分成功不会被伪报为全部成功。
- `scripts/profile_real_omics.py` 生成 `raw/<dataset>/data_inventory.json`；列出来源文件实际格式、压缩包成员、工作表结构及列预览。

**当前执行边界：** GitHub 分支的下载代码已提交；当前会话无法直接连接或运行用户 WSL 中的 F: 盘，也没有已执行成功的上述新数据下载记录。以实际生成的下载清单为准，不能把已核实公开下载 URL 等同于文件已经保存到用户机器。

## 3. 训练所需的标准数据结构

建议先保留原生文件，建立一层可追溯的**标准表示**，分为如下逻辑表。不要强行把不同队列的特征按名称拼成一张大表：

1. `participants`：`dataset_id,participant_id,cohort,study_group`；一名参与者多次访问共用一个 ID。
2. `visits`：`dataset_id,participant_id,visit_id,time_value,time_unit,event,intervention`；时间可以为相对基线天数/小时，必须有单位与锚点；若只有先后顺序，保存 `visit_order` 而非伪造时间间隔。
3. `specimens`：`dataset_id,specimen_id,participant_id,visit_id,biospecimen_type,collection_time`；区分血浆、血清、尿液、粪便。若一次访问采集了不同体液，要分别保存。
4. `assays`：`assay_id,specimen_id,modality,platform,batch,raw_file,processing_version`；明确同一样本可以上多次检测平台。
5. `features`：`dataset_id,assay_type,feature_id,feature_name,unit,measurement_scale,annotation_confidence`；蛋白质、代谢物的命名体系可能不同，未鉴定峰需要保留原始特征 ID。
6. `measurements`：`assay_id,feature_id,value,observed,missing_reason`，或等价的宽矩阵加 feature/样本掩码；原始 NA、低于检测限、QC 剔除应分开记录。
7. `splits`：`dataset_id,participant_id,split,seed`，保证测试和训练参与者严格互斥。

核心索引可表示为

```text
participant_id ──< visit_id ──< specimen_id ──< assay_id
                                              |
                                              +── proteomics values [features_P]
                                              +── metabolomics values [features_M]
                                              +── measured-feature masks [features_P/M]
```

注意观察数据的两层缺失：
- **模态级缺失**：某次访问没有蛋白质组或代谢组检测，对应 PULSE 现有 0/1/2 掩码（0 已观测，1 训练主动掩盖，2 自然缺失）。
- **特征级缺失**：一个已检测模态中的部分特征 NaN、低于检测限或因 QC 排除。必须另行保存特征掩码，不能将它们直接喂给当前 MSE 重建损失，也不能无条件把 NA 当 0。

### 模态数据的单位意义

蛋白质组可能是 LC-MS 相对丰度、标准化峰强度或 SomaScan RFU；代谢组可能是质谱峰强度、已归一化比值或特定浓度单位。所谓 `protein_001` 与 `metabolite_001` 都是**测量特征**，并不自动等于一个严格定义的生理状态变量。各来源原始值的变换（是否 log、归一化批次和平台）要在每套数据中独立核查。

原始组学矩阵通常是 `samples × features` 或 `features × samples`，BIOM 的形状需从矩阵元数据判别，不能仅按 CSV 行列数确定方向。所有访问链接必须由原生 ID 和真实 metadata 建立，不允许按行号配对。

## 4. 当前 PULSE 实现与标准数据的接口差距

从现有 `pulse/data/dataset.py`、`pulse/model.py` 和 `pulse/training/trainer.py` 可直接确定：

- 模型能处理每名患者的不同长度访问列表，支持整模态的自然缺失及训练中的主动遮盖。
- 现有时间模块根据访问顺序迭代 GRU/Transformer 状态，**没有显式读入真实时间间隔 Δt**；iPOP 的 3 个月间隔和应激期间 1 天间隔当前会表现为同一次状态更新。为真实长时间序列研究需要补充时间编码或时差动态。
- 原始 Su2020 的蛋白质组缺失约 1.73%，代谢组缺失约 18.77%；现有输入编码器与重建损失不能直接容纳 NaN。需要在保持原始缺失掩码的条件下构建训练集内的特征筛选、处理与观测掩码重建。
- 现有 `create_synthetic_data` 对每个患者/时间/模态独立生成高斯向量，这只能检测训练程序能否运行，不能验证生物纵向状态预测或历史信息的价值。

## 5. 数据接入以后必须回答的统计量

逐数据集应报告真实数值，而非按论文总受试者数推断：

- 独立患者数、每位患者的访问数分布与实际时间间隔分布。
- 两种模态各自观测次数、同一次访问真正配对的次数。
- 拥有至少两次、至少三次**双组学配对访问**的患者数量。
- 模态缺失率、特征缺失率、技术批次/采样体液的差异、无标签特征比例。
- 训练/测试患者独立划分后能利用的有效访问和可评估特征数量。
- 是否有真实干预/治疗变化与相应时间标签，是否满足预测场景。

Su2020 目前配对访问数为 252，但 `linked_patients_with_two_draws=124` 是至少一种组学可用的 T1/T2 患者数，不等于两次都**双组学**配对的患者数；接下来应从 `assay_mapping.csv` 严格统计后者。

## 6. 合成数据研究的启动条件（暂不实施）

先获得各公开队列的文件、确认原始数据结构与有效配对，然后再确定生成数据究竟用于哪一种科学任务：

- **算法可识别性实验**：构建已知潜在生理动力学与双组学观测机制，模型能否恢复状态、历史是否改善预测；需要可解释的真值与独立变化参数。
- **测量与缺失机制实验**：加入真实组学类型特有的测量噪声、检测限、整模态缺失、就诊机制与时差，逐一控制。
- **数据增强/真实分布模拟**：以训练组患者估计分布及相关结构，绝不让生成器接触测试患者信息，分别评估统计保真度和真实测试集收益。

模型的候选机制、参数范围、生成数据格式与判别实验应在**真实队列数据结构完成实证盘点之后**制定，而非现在凭空拟合潜在状态。
