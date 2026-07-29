# Seurat V5 模块化分析流程

本版本把正式分析整理为四个清晰的主入口：

1. `运行_01_入口_至大亚群注释.R`
2. `运行_02_入口_从大亚群注释至小亚群完成.R`
3. `运行_03_入口_小亚群人工注释与汇总.R`
4. `运行_04_入口_单大亚群重新聚类与注释.R`

入口1和入口2支持断点续跑。入口1结束后，入口2直接读取已经完成大亚群注释的
Seurat对象，不会重复QC、合并、PCA、Harmony或大亚群聚类。入口3读取入口2
的全部结果，生成可编辑的XLSX注释总表并输出最终注释对象。入口4用于从任意
大亚群注释对象中挑选一个大亚群，反复尝试PC/res并重新完成聚类和注释。

## 入口1：运行至大亚群注释完成

先修改：

```text
config/01_流程参数.R
config/02_大亚群注释参数.R
```

然后在RStudio中Source：

```r
source(
  "C:/你的路径/SeuratV5模块化流程/运行_01_入口_至大亚群注释.R",
  encoding = "UTF-8"
)
```

入口1依次完成：

```text
原始数据读取
→ metadata与Seurat对象
→ 逐样本QC
→ 合并样本
→ 标准化、细胞周期和PCA
→ Harmony与正式聚类
→ cluster marker
→ 大亚群注释
```

最终结果：

```text
result/5.大亚群注释完成_sce.rds
```

如果第一次运行到大亚群注释时尚未填写cluster映射，请修改
`config/02_大亚群注释参数.R`，再重新Source入口1。断点恢复会直接从大亚群注释
模块继续，不会重复前面的耗时计算。

## 入口2：从大亚群注释运行至小亚群完成

先修改：

```text
config/04_小亚群分析参数.R
config/07_外部大亚群RDS输入参数.R
```

然后在RStudio中Source：

```r
source(
  "C:/你的路径/SeuratV5模块化流程/运行_02_入口_从大亚群注释至小亚群完成.R",
  encoding = "UTF-8"
)
```

入口2读取 `5.大亚群注释完成_sce.rds`。如果大亚群注释配置刚被修改，
入口2会先从正式聚类对象快速刷新大亚群注释；配置没有变化时直接使用现有对象。
随后按照 `sce$celltype` 逐个拆分大亚群。
每个大亚群分别执行：

```text
RNA counts重建
→ NormalizeData
→ FindVariableFeatures
→ ScaleData
→ PCA
→ 可配置Harmony
→ resolution比较
→ 正式小亚群聚类
→ UMAP
→ FindAllMarkers
→ 内置先验marker DotPlot
→ 结果汇总
```

大亚群名称支持关键词兼容匹配。名称中包含以下关键词即可匹配对应面板，
且不区分大小写：

| 名称关键词 | 匹配大亚群 |
| --- | --- |
| `B` | B cell |
| `Epi` | Epithelial |
| `Mye` | Myeloid |
| `T` | T/NK |
| `Fib` | Fibroblast |
| `Endo` | Endothelial |
| `Plasma` | Plasma |

长关键词优先于单字母关键词。如果一个关键词同时命中对象中的多个真实名称，
程序会停止并要求使用完整名称，避免选错。

## 使用任意大亚群注释RDS进入下游

修改：

```text
config/07_外部大亚群RDS输入参数.R
```

例如：

```r
downstream_use_external_annotated_rds <- TRUE
downstream_external_annotated_rds <- "D:/singlecell/annotated_object.rds"
downstream_external_celltype_column <- "major_celltype"
downstream_external_rna_assay <- "RNA"
```

入口2和入口4都会读取这套配置。外部RDS必须是Seurat对象，并包含大亚群注释
metadata列和可用于重建对象的RNA counts。程序会把配置的注释列统一复制为
`celltype`供下游使用。

## 小亚群PC与resolution选择

每个大亚群都会独立计算：

- `co1`：累计解释率达到90%且单个PC贡献率低于5%的首个PC；
- `co2`：解释率下降幅度大于0.1%的最后一个拐点；
- 推荐PC：`min(co1, co2)`，同时受最小PC和实际PCA维数约束。

正式PC自动使用该推荐值，不再要求用户选择PC。聚类算法默认固定为
Louvain（`subcluster_algorithm <- 1`）。

对每个大亚群，程序会先输出：

```text
不同resolution_UMAP.pdf
不同resolution_clustree.pdf
聚类候选参数与自动评分.csv
聚类候选评分图.pdf
PCA自动推荐结果.csv
```

随后在控制台逐个大亚群提问：

```text
请输入正式resolution（候选值：0.2、0.4、0.6、0.8、1）
```

- 输入候选值：采用该resolution；
- 直接回车：采用自动评分最高的resolution；
- 输入非候选值：原地重新询问，不会终止流程，也不会写入错误断点。

如果使用非交互式 `Rscript`，可在
`config/04_小亚群分析参数.R` 中预先填写：

```r
subcluster_resolution_overrides <- list(
  "B cell" = 0.4,
  "Myeloid" = 0.6
)
```

## 小亚群结果

结果保存在：

```text
result/小亚群自动分析/
```

主要汇总文件包括：

- `6.所有大亚群小亚群聚类完成_sce.rds`
- `全部大亚群正式聚类参数汇总.csv`
- `全部小亚群细胞数量汇总.csv`
- `全部细胞与小亚群对应关系.csv`
- `按大亚群命名的独立对象/`

每个大亚群有独立文件夹，保存PCA推荐、resolution比较、正式参数、UMAP、
marker、DotPlot、热图、先验marker点图、细胞映射表、子对象RDS及阶段检查点。
每完成一个大亚群，入口2还会把独立对象复制到
`按大亚群命名的独立对象/`，例如：

```text
T_NK_小亚群聚类完成_sce.rds
Myeloid_小亚群聚类完成_sce.rds
```

## 内置小亚群先验marker

代码包在 `data/先验marker原始表/` 中保留了6个原始工作簿：

- `T_marker_31(1).xlsx`
- `B_markers(1).xlsx`
- `PC_markers(1).xlsx`
- `Mye_markers(1).xlsx`
- `Fib_markers(1).xlsx`
- `EC_markers.xlsx`

运行时不依赖Excel读取包；6个工作簿已经规范化并内置到
`config/05_小亚群marker点图参数.R`。该配置同时加入结直肠癌和肿瘤微环境
单细胞研究中常用的补充marker，并记录marker来自“用户Excel”“文献补充”
或两者共同支持。

入口2在每个大亚群完成小亚群聚类后，会自动匹配 T/NK、B、Plasma、Myeloid、
Fibroblast、Endothelial 等面板，并输出多页DotPlot、基因存在性检查、
各小亚群marker覆盖情况及DotPlot底层数值。先验marker仅用于辅助判读cluster，
不会自动强制覆盖小亚群名称。

## 入口3：小亚群人工注释与汇总

首次使用前安装：

```r
install.packages("openxlsx")
```

可按需修改：

```text
config/06_小亚群人工注释参数.R
```

然后Source：

```r
source(
  "C:/你的路径/SeuratV5模块化流程/运行_03_入口_小亚群人工注释与汇总.R",
  encoding = "UTF-8"
)
```

入口3会生成：

```text
result/小亚群自动分析/小亚群人工注释总表.xlsx
```

工作簿包含：

- `注释填写表`：每个大亚群中的每个cluster占一行；
- `先验亚群名称库`：内置marker面板中的候选小亚群名称，可自由增加或删除；
- `使用说明`：填写规则和输出说明。

程序根据先验marker DotPlot的平均标准化表达评分预填第一、第二候选名称。
黄色的`正式小亚群名称`列默认采用第一候选，但可以任意修改为中文、英文或
自定义名称，也可以删除不需要的候选名称。先验得分只是辅助建议，不会替代
人工检查UMAP、marker和样本构成。

保存并关闭Excel后，回到R控制台按回车，或再次Source入口3。最终文件统一保存到：

```text
result/小亚群注释完成文件/
```

其中包括：

- `所有大亚群小亚群注释完成_sce.rds`；
- 按大亚群名称保存的独立注释完成对象；
- 每个大亚群的注释UMAP；
- `小亚群人工注释对应关系.csv`；
- 填写后用于本次注释的XLSX副本；
- `小亚群注释完成文件清单.csv`。

## 入口4：单大亚群重新聚类与注释

入口4适合在整体流程完成后，对某一个大亚群单独重新探索。先修改：

```text
config/08_单大亚群重分析参数.R
```

然后Source：

```r
source(
  "C:/你的路径/SeuratV5模块化流程/运行_04_入口_单大亚群重新聚类与注释.R",
  encoding = "UTF-8"
)
```

运行过程：

1. 显示对象中的全部大亚群，可输入编号、完整名称或兼容关键词；
2. 只对选中的大亚群重新进行RNA重建、PCA和Harmony；
3. 每轮允许输入任意有效PC数量和任意大于0的resolution；
4. 保存并尝试自动打开本轮UMAP，同时显示每个cluster的细胞数；
5. 询问是否采用本轮方案：
   - 输入`y`：作为最终方案；
   - 输入`n`或直接回车：继续尝试下一组PC/res；
   - 输入`q`：退出，已有试算结果仍保留；
6. 确认后自动运行正式UMAP、FindAllMarkers、热图和先验marker DotPlot；
7. 生成该大亚群的XLSX人工注释表，填写后输出最终注释对象。

结果保存在：

```text
result/单大亚群重新聚类与注释/大亚群名称/
```

其中`PC_res试算记录/`保存每次试算的UMAP、参数和cluster数量，
`PC_res试算历史.csv`记录所有尝试以及最终采用的组合。

## 断点续跑

- 入口1使用 `config/01_流程参数.R` 中的 `resume_mode`。
- 入口2使用 `config/04_小亚群分析参数.R` 中的
  `subcluster_resume_mode`。
- 默认均为 `"auto"`。
- 入口2每完成一个大亚群就保存结果；后续大亚群报错时，已经完成的亚群不会重算。
- 每个大亚群内部还保存PCA/Harmony和正式聚类检查点；marker阶段报错时不会重复
  PCA、Harmony或resolution选择。
- 入口4确认最终PC/res后会保存待人工注释对象；注释阶段中断时重新Source入口4
  可直接恢复，不需要再次试算。
- 需要全部重算时，把相应的resume参数改为 `"restart"`。

## 辅助脚本

- `辅助_仅重新进行大亚群注释.R`：修改大亚群映射后单独重新注释。
- `辅助_小亚群marker点图.R`：入口2已自动绘图；该脚本用于修改marker配置后单独重画。

## 文件结构

```text
SeuratV5模块化流程/
├─ 运行_01_入口_至大亚群注释.R
├─ 运行_02_入口_从大亚群注释至小亚群完成.R
├─ 运行_03_入口_小亚群人工注释与汇总.R
├─ 运行_04_入口_单大亚群重新聚类与注释.R
├─ 辅助_仅重新进行大亚群注释.R
├─ 辅助_小亚群marker点图.R
├─ data/
│  └─ 先验marker原始表/
│     ├─ T_marker_31(1).xlsx
│     ├─ B_markers(1).xlsx
│     ├─ PC_markers(1).xlsx
│     ├─ Mye_markers(1).xlsx
│     ├─ Fib_markers(1).xlsx
│     └─ EC_markers.xlsx
├─ config/
│  ├─ 01_流程参数.R
│  ├─ 02_大亚群注释参数.R
│  ├─ 03_处理前参数覆盖.R
│  ├─ 04_小亚群分析参数.R
│  ├─ 05_小亚群marker点图参数.R
│  ├─ 06_小亚群人工注释参数.R
│  ├─ 07_外部大亚群RDS输入参数.R
│  └─ 08_单大亚群重分析参数.R
├─ R/
│  ├─ 00_载入依赖包.R
│  ├─ 01_辅助函数.R
│  ├─ 01_流程恢复函数.R
│  ├─ 02_识别并读取输入数据.R
│  ├─ 03_metadata与Seurat对象.R
│  ├─ 04_逐样本质量控制.R
│  ├─ 05_合并样本.R
│  ├─ 06_标准化细胞周期与PCA.R
│  ├─ 07_Harmony与正式聚类.R
│  ├─ 08_寻找cluster_marker.R
│  ├─ 09_大亚群注释.R
│  ├─ 10_小亚群自动分析函数.R
│  ├─ 11_小亚群marker点图函数.R
│  ├─ 12_小亚群人工注释函数.R
│  └─ 13_单大亚群重分析函数.R
└─ archive/
   └─ 原始单文件版本.R
```

## 注意事项

- 保持整个文件夹结构不变，不要单独移动入口脚本。
- 请打开入口脚本后运行整个文件的Source，不要只选择部分代码运行。
- Windows路径建议使用 `/`，例如 `C:/singlecell/project`。
- `clustree` 是可选依赖；未安装或绘图失败时不阻断UMAP和正式聚类。
- 入口3和入口4需要`openxlsx`；入口1和入口2不依赖该包。
- 自动评分用于辅助选择，不替代marker、生物学知识和样本构成复核。
