# iPOP 与 P100：已恢复来源的真实文件调查

调查日期：2026-11-08。代码核对基准：`omics` 分支 `148aaccb6b846feebee1404331972137d9cea810`。

本报告记录本次实际取得并打开的 iPOP 季节性研究处理文件及 P100 出版方附件。文件字节数、SHA-256、R 对象、全部 Excel 工作表、原始字段和特征目录见 `docs/data_structure/restored/`。这些文件在远程研究环境中读取，尚未写入 WSL 的 F 盘。Su2020、HMP2、COVAAR、CSIRO 的 F 盘文件在当前环境不可访问，本报告不替代这四个队列的实证结构调查。合成数据研究尚未启动。

## iPOP：已恢复的处理数据版本

原入口为 `s3://hmpdcc/ihmp/t2d/host/`。替代入口来自 Sailani et al. 2020《Deep longitudinal multiomics profiling reveals two biological seasonal patterns in California》的数据声明：

- 论文：<https://pmc.ncbi.nlm.nih.gov/articles/PMC7529769/>
- 作者指定的处理数据：<https://doi.org/10.6084/m9.figshare.12376508.v1>
- 已取得文件：`Multi_Omics_Seasonal.RData`，53,439,956 bytes。
- 文件直链：<https://ndownloader.figshare.com/files/22807673>
- SHA-256：`2763c0db3d98fed964435a190439a50257259f560e502601714bf49124b504e0`。
- MD5：`d1c276cd48130916ef3f68b37bf80e74`，与 Figshare 公布值一致。
- 版本 1，CC BY 4.0。该版本是作者处理后的季节性分析数据，不等于已恢复原 S3 上全部实验文件。

实际反序列化得到 14 个具名 R 数据框，共七组测量矩阵与对应的元数据框。测量矩阵均为观测行、特征列，数值为 R double。完整列名与逐列缺失数保存为 `ipop_native_fields.csv.gz`，共 11,693 条字段记录，其中 11,658 条为测量特征、35 条为元数据字段。

| 原生 R 对象 | 行×列 | 缺失单元格 |
|---|---:|---:|
| `Clinical_labs_annotation_colData` | 969×5 | 250 |
| `clinical_labs_Data` | 969×51 | 3480 |
| `Cytokines_annotation_colData` | 964×5 | 252 |
| `Cytokines_df_Data` | 964×66 | 0 |
| `Gut_annotation_colData` | 855×5 | 1712 |
| `gut_df_Data` | 855×96 | 0 |
| `Metabolomics_annotation_colData` | 981×5 | 240 |
| `Metabolomics_df_Data` | 981×724 | 0 |
| `Nasal_annotation_colData` | 837×5 | 228 |
| `Nasal_df_Data` | 837×80 | 0 |
| `Proteomics_annotation_colData` | 950×5 | 253 |
| `Proteomics_df_Data` | 950×302 | 0 |
| `RNA_annotation_colData` | 882×5 | 241 |
| `RNA_df_Data` | 882×10339 | 0 |

各元数据框均包含五个原始字段：

| 原字段 | 原始类型与实际内容 | 科学含义及处理要求 |
|---|---|---|
| `SubjectID` | R factor，字符型受试者标识；不同模态有少量 NA | 个体连接键；不能用数据框行名代替 |
| `Date` | R POSIXct 的底层秒数，已出现异常公元年份 | 保留原始秒数和类信息；不能直接赋予现代日历日期 |
| `IRIS` | R factor，`IR`、`IS`、NA | 胰岛素抵抗/敏感分组；NA 仍保留 |
| `BMI` | 原文件存储为 R factor，其 levels 为数值字符串 | 体质指数；需要显式数值转换，不能把 factor 内码当 BMI |
| `Time` | R double；蛋白质组和代谢组范围均为 2–355 | 一年中的日序，不是访视序号或距基线天数 |

`Date` 的底层数值按 UTC 秒数解码后，蛋白质组范围为 `0013-05-08 07:52:58` 至 `0017-01-09 07:52:58`，代谢组范围为 `0010-03-29 07:52:58` 至 `0017-01-09 07:52:58`。这与研究开展年代不符。两个模态内，解码日期的年内日序都与 `Time` 一致，说明 `Time` 并未提供独立的现代年份证据。原始日期解析、时区及年份恢复仍需源级核对；本次没有自动加 2000 年。

蛋白质组来自血浆 SWATH-MS，原列名为 `IGLL5`、`MASP2`、`APOL1`、`ATRN`、`ATRN.1` 等。带 `.1` 后缀的列保留为独立原始特征。其 286,900 个单元格均为有限值，没有 NA；最小值 −11.38575，中位数 −1.1877675，最大值 12.71131。它是具有负值的处理后量，不可解释为绝对浓度，也不能再次盲目取对数。R 对象未携带单位、变换公式和处理前缺失掩码；仅凭数值范围不能确定其对数底数或中心化方法。

代谢组来自血浆 LC-MS，实际 724 列。列名前缀分别为 `pHILIC` 144 列、`nHILIC` 173 列、`pRPLC` 127 列、`nRPLC` 280 列，例如 `pHILIC_142.1225_1.3`。命名包含色谱/离子模式及数值标记，完整原名保留，不能直接改写为已确证化合物名称。710,244 个测量单元格均为有限值；最小值 0.195553，中位数 205,793.25，最大值 226,000,000,000。该量级与绝对浓度单位无直接对应关系。

论文方法说明：代谢特征先按检出比例筛选，缺失值由近邻法补齐，批内和批间变异以重复 QC 的 LOESS 归一化校正。因此当前代谢矩阵的零 NA 不能解释为原始实验零缺失。文件未提供原始缺失位置，后续评估只能另行主动遮盖已发布的处理值。论文正文写 722 个代谢物，而这一具体版本有 724 列；本报告采用文件实测数量，不将两者强制统一。

### 原生身份连接及纵向统计

元数据框的行名如 `V2`、`V3`，矩阵行名则为 `1`、`2`。这些是行标签，不能作为跨模态生物样本 ID。文件没有独立的采血样本号、分装样本号、检测批次和疫苗/感染事件字段。下列统计只在两份元数据上使用非空的 `SubjectID + Date` 精确连接；两边该键均唯一，连接采用 `validate='one_to_one'`，没有按行顺序配对测量矩阵。

| 实测统计 | 数值 |
|---|---:|
| 蛋白质组元数据行 | 950 |
| 蛋白质组缺失 SubjectID 行 | 1 |
| 蛋白质组非空独立 SubjectID | 106 |
| 代谢组元数据行、非空独立 SubjectID | 981、102 |
| 两模态共同的 SubjectID + Date 记录 | 912 |
| 拥有共同时间记录的受试者 | 102 |
| 至少两次共同时间记录的受试者 | 91 |
| 至少三次共同时间记录的受试者 | 85 |
| 仅蛋白质组有效键记录 | 37 |
| 仅代谢组有效键记录 | 69 |
| 相邻共同时间记录间隔数 | 810 |
| 按原始秒数相减得到的间隔天数：最小/Q1/中位/Q3/最大 | 1 / 14 / 75 / 104 / 551 |

912 是元数据交集数量。将其提升为可训练的同次测量配对，还需要核实作者原生的矩阵行与各自元数据行对应方式；确认同一生物样本则另需样本身份来源。解析器没有导出按位置拼接的测量矩阵。上述日期间隔只描述发布值相减结果，日历异常仍保留。

逐受试者连接键、未匹配键及缺失 SubjectID 记录可由解析器写入显式指定的 `--private-output-dir`。它们不进入仓库。完整访问次数分布和各对象缺失统计在 `ipop_structure.json`。

作者仓库 <https://github.com/aametwally/ipop_seasonal> 的 `code/stat.R` 还引用 `ls`、`sc` 两个对象，而此公开 RData 仅含上表 14 个对象，没有这两个名称。该脚本不能未经核对就作为当前文件的样本连接依据。

## P100：已恢复的出版方附件

来源：Price et al. 2017《A wellness study of 108 individuals using personal, dense, dynamic data clouds》，DOI <https://doi.org/10.1038/nbt.3870>。Europe PMC 打包接口此前未返回真实附件，本次从出版方页面的实际链接恢复了 11 个文件，逐个核对文件签名并打开 Excel/ZIP。精确直链、字节数和 SHA-256 见 `retrieval_manifest.json`。

| 出版方精确文件名 | 字节数 | 实际内容 |
|---|---:|---|
| `41587_2017_BFnbt3870_MOESM26_ESM.pdf` | 2,352,615 | 补充图和文字 PDF |
| `41587_2017_BFnbt3870_MOESM27_ESM.zip` | 220,066 | 横断面相关网络 Cytoscape 会话 |
| `41587_2017_BFnbt3870_MOESM28_ESM.zip` | 174,054 | 变化量相关网络 Cytoscape 会话 |
| `41587_2017_BFnbt3870_MOESM29_ESM.zip` | 74,404 | 分析代码 ZIP |
| `41587_2017_BFnbt3870_MOESM30_ESM.xlsx` | 325,567 | Table 1：参与者信息和特征目录 |
| `41587_2017_BFnbt3870_MOESM31_ESM.xlsx` | 253,924 | Table 2：跨组学横断面相关 |
| `41587_2017_BFnbt3870_MOESM32_ESM.xlsx` | 1,270,450 | Table 3：全部横断面相关 |
| `41587_2017_BFnbt3870_MOESM33_ESM.xlsx` | 195,941 | Table 4：跨组学变化量相关 |
| `41587_2017_BFnbt3870_MOESM34_ESM.xlsx` | 1,878,234 | Table 5：全部变化量相关 |
| `41587_2017_BFnbt3870_MOESM35_ESM.xlsx` | 52,662 | Table 8：遗传性状与临床比较 |
| `41587_2017_BFnbt3870_MOESM36_ESM.xlsx` | 83,873 | Table 14：年龄/性别校正 |

`MOESM27`、`MOESM28` 各包含一个 `.cys`，继续打开后得到 Cytoscape 的 XGMML 网络及属性表。它们分别保存横断面和变化量相关网络，没有逐人逐次的组学测量表。`MOESM29` 是分析代码，其 README 明确将分析输入指向 dbGaP。

七个 Excel 共 19 张实际工作表，全部读取如下。所有非 Legend 工作表的表头在第 1 行。Legend 是对表的说明，不按患者记录统计。

| 文件缩写 | 原始工作表 | 工作表行×列 | 非空行 | 实际数据行 |
|---|---|---:|---:|---:|
| MOESM30 | `Legend` | 14×2 | 9 | 说明 |
| MOESM30 | `Demographics` | 109×32 | 109 | 108 |
| MOESM30 | `Genetic Traits` | 131×12 | 131 | 130 |
| MOESM30 | `Clinical Lab Tests` | 219×6 | 219 | 218 |
| MOESM30 | `Metabolomics` | 644×12 | 644 | 643 |
| MOESM30 | `Proteomics` | 385×7 | 385 | 384 |
| MOESM30 | `Microbiome` | 4617×8 | 4617 | 4616 |
| MOESM31 | `Legend` | 11×2 | 9 | 说明 |
| MOESM31 | `Table` | 3471×6 | 3471 | 3470 |
| MOESM32 | `Legend` | 18×2 | 9 | 说明 |
| MOESM32 | `Table` | 21126×6 | 21126 | 21125 |
| MOESM33 | `Legend` | 18×2 | 9 | 说明 |
| MOESM33 | `Table` | 2407×6 | 2407 | 2406 |
| MOESM34 | `Legend` | 18×2 | 9 | 说明 |
| MOESM34 | `Table` | 31223×6 | 31223 | 31222 |
| MOESM35 | `Legend` | 18×2 | 11 | 说明 |
| MOESM35 | `Table` | 24×8 | 24 | 23 |
| MOESM36 | `Legend` | 18×2 | 8 | 说明 |
| MOESM36 | `Table` | 1261×5 | 1261 | 1260 |

### 字段语义、特征及可用性

`MOESM30 / Demographics` 一行一名参与者，108 人、32 个原字段。除 `participant_id`、`sex`、`age`、`ancestry`、吸烟/饮酒问卷、`genome_id`、`genome_tech` 外，还包含按平台和轮次组织的观测数量列：

- `clinical.Genova.1` 至 `.3`、`clinical.Quest.1` 至 `.3`。
- `metabolomics.metabolon.1` 至 `.3`。
- `proteomics.olink.1` 至 `.3`、`proteomics.srm_mass_spec.1` 至 `.3`。
- `microbiome.secondgenome.1` 至 `.3`、`average.activity.fitbit`。

这些组学列为该参与者该轮被观测的特征数量。例如代谢组第 1 轮非空值范围 639–642，Olink 各轮非空值均为 184，SRM 第 1 轮为 63–79。它们不是单个蛋白质/代谢物丰度。空单元格、问卷中的 `-` 和 `No response` 分别保留；逐列类型、空单元格数和字面缺失标记见 `p100_native_fields.csv`。

由原始观测数量大于零计算的同轮双模态可用性如下。这些统计描述公开人口学表记录的检测覆盖，不代表测量矩阵已经公开取得，也不证明同一分装样本。

| 蛋白质检测平台 | 第1轮与代谢组共同可用 | 第2轮 | 第3轮 | 至少两轮共同可用的参与者 | 三轮均共同可用 |
|---|---:|---:|---:|---:|---:|
| Olink | 80 | 79 | 75 | 79 | 75 |
| SRM | 96 | 87 | 89 | 94 | 71 |
| 任一蛋白平台 | 104 | 99 | 91 | 100 | 86 |

`MOESM30 / Metabolomics` 为 643 条特征注释，12 字段：`biochemical, super_pathway, sub_pathway, comp_id, platform, chemical_id, ri, mass, cas, pubchem, kegg, hmdb`。平台数量为 `lc/ms-neg` 347、`lc/ms-pos` 159、`lc/ms-polar` 74、`gc-fid` 34、`gc/ms` 29。名称、通路、平台与数据库编号均保存在 `p100_metabolomics_annotations.csv`；原表没有每个参与者的丰度列。

`MOESM30 / Proteomics` 为 384 条检测特征注释，7 字段：`source, abbreviation, name, description, category, uniprot, peptide`，其中 SRM 200 行、Olink 184 行。SRM 的肽段与 Olink 的检测目标不能只按名称合并。384 条注释也不等于每名参与者测得 384 个唯一蛋白质。

`Clinical Lab Tests` 为 218 条检测定义，字段 `name, unit, vendor, vendor_id, hmdb_id, source`；`source` 包括不同材料，单位逐项保留，例如 `micromol/dL`、`mmol/mol creatinine`。`Genetic Traits` 为 130 条遗传性状研究注释；`Microbiome` 为 4,616 条 OTU 分类注释，首列无原始表头，其余为 `Kingdom, Phylum, Class, Order, Family, Genus, Species`。输出保留无表头列，不替原作者虚构字段。

`MOESM31–34 / Table` 均为相关网络的边表，原字段 `coefficient, pval, pval_adj, test, annotations_1, annotations_2`。一行是一对特征的相关结果，不能当作一个患者样本。`MOESM35 / Table` 为 23 条遗传性状与临床指标的比较；`MOESM36 / Table` 为 1,260 条年龄/性别校正记录，字段 `analyte, age_pvalue, age_transformed, sex_pvalue, sex_transformed`。

原始个人级数据申请入口为 <https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001363.v1.p1>，需要 Authorized Access。此次公开附件可支持特征字典、检测覆盖和网络结果研究；无法直接提供 PULSE 所需的逐人逐次双组学输入与持出测量真值。

## 复现此次已完成的文件调查

下载程序仍默认写入 F 盘，已存在有效文件跳过。iPOP 的 S3 发现不可用时回退到上述作者发布的处理文件；已有该文件时先校验本地副本，不再访问 S3 目录。P100 改为这 11 个已经实际验证的出版方附件。

```bash
python -m pip install -r requirements-structure.txt
python scripts/download_real_omics.py --datasets ipop p100
python scripts/inspect_restored_omics.py \
  --private-output-dir /mnt/f/datasets/pulse-omics/processed/restored_sources/private_links
```

解析程序原生读取 RData 及全部 Excel 工作表，不生产合成数据、不训练模型，也不导出个体测量矩阵。仓库中的目录和统计来自本次实际取得的上述版本，不能据此推断 F 盘其他四个队列的结构调查已完成。
