# scRNA-fast2
A script that can accelerate the procedures of scRNA-seq analysis. You can download it and source the main R script to set up
# Seurat V5模块化分析流程

本项目把原来接近两千行的单文件脚本拆分为配置文件、辅助函数和按顺序执行的功能模块。正常分析时不需要逐个打开模块，也不要单独复制某个模块运行；入口脚本会使用 `source()` 自动按正确顺序载入。

## 最简单的运行方法

1. 可以先双击 `SeuratV5模块化流程.Rproj` 打开项目，再打开 `config/01_流程参数.R`，修改 `project_dir`、`raw_data_dir`、`result_dir` 和分析参数。
2. 在RStudio中打开 `运行_01_完整流程.R`，点击 **Source**。
3. 第一次运行时，如果不存在样本metadata，脚本会生成 `样本metadata.xlsx`。在Excel中填写并保存后，回到R控制台按回车继续。
4. 查看PC梯度图并在控制台选择正式PC数量，再查看resolution比较图并选择正式resolution。
5. 流程会生成cluster marker和经典marker图。如果尚未填写注释编号，会在大亚群正式赋值前停止。
6. 打开 `config/02_大亚群注释参数.R`，按照marker结果填写“自定义亚群名称 = cluster编号”。名称可以修改，也可以自由增加或删除大亚群。
7. Source `运行_02_仅重新进行大亚群注释.R`。它会读取正式聚类对象并直接完成大亚群注释，不会重复前面的耗时计算。

## 报错后自动继续

主入口现在会在每个功能块成功后，把下一步需要的对象保存到 `result/流程检查点`，并把当前状态写入 `result/流程运行状态.rds`。

- 默认 `resume_mode <- "auto"`：某个功能块报错后，修正问题并再次 Source `运行_01_完整流程.R`，流程会恢复上一功能块的检查点，从报错功能块重新执行。
- `07_Harmony与正式聚类.R` 还包含更细的 PC 梯度和 resolution 计算检查点；后续绘图报错时不会重复这些耗时计算。
- `clustree` 是可选辅助图。与 `ggplot2 4.x` 不兼容时只产生警告，不再阻断 resolution UMAP 和正式聚类。
- 修改原始数据、metadata、QC阈值或主要聚类参数后，把 `resume_mode` 改成 `"restart"`，避免沿用旧检查点。

也可以在R控制台中执行：

```r
source("C:/你的路径/SeuratV5模块化流程/运行_01_完整流程.R", encoding = "UTF-8")
```

## 文件结构

```text
SeuratV5模块化流程/
├─ SeuratV5模块化流程.Rproj
├─ 运行_01_完整流程.R
├─ 运行_02_仅重新进行大亚群注释.R
├─ config/
│  ├─ 01_流程参数.R
│  ├─ 02_大亚群注释参数.R
│  └─ 03_处理前参数覆盖.R
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
│  └─ 09_大亚群注释.R
└─ archive/
   └─ 原始单文件版本.R
```

## 各配置文件的用途

- `01_流程参数.R`：日常最常修改的文件，包括路径、输入格式、QC阈值、PC梯度和resolution梯度。
- `02_大亚群注释参数.R`：在 `celltype_cluster_map` 中填写不同大亚群对应的cluster编号。列表名称就是最终的 `sce$celltype` 名称，可以使用中文或英文并自由新增、删除、重命名；流程仍保留原脚本的cluster向量直接注释风格，不读入人工注释表。
- `03_处理前参数覆盖.R`：第二处参数覆盖区。只有把 `use_parameter_override` 设置为 `TRUE` 时才生效。

大亚群配置示例：

```r
celltype_cluster_map <- list(
  "T细胞和NK细胞" = c(0, 3, 7),
  "B细胞" = c(2, 8),
  "髓系细胞" = c(1, 5),
  "自定义新增亚群" = c(4)
)
```

引号中的文字会原样写入最终Seurat对象的 `celltype` 列。一个cluster只能属于一个名称；没有填写的cluster会保留为 `unknown`。修改配置后只需重新 Source `运行_02_仅重新进行大亚群注释.R`。

## 注意事项

- 请保持整个文件夹结构不变，不要只移动入口脚本。
- 请打开入口脚本后运行整个文件的 **Source**，不要只选中定位代码点击Run。新定位逻辑在无法取得脚本路径时会停止，不再使用可能被 `setwd()` 改变的工作目录进行猜测。
- Windows路径建议使用 `/`，例如 `C:/singlecell/project`。
- 输入数据仍支持自动识别H5、10X和CSV；10X会自动判断features文件中的基因名称列。
- 如果找不到 `MT-` 或 `mt-` 线粒体基因，流程会停止并显示基因行名示例，不会把线粒体比例错误设置为0。
- 样本metadata模板和“本次实际读入的样本metadata”均使用XLSX。模板带筛选、冻结首行和自动列宽，不再受CSV的UTF-8/GBK编码影响。
- XLSX生成与读取需要 `openxlsx` 包；如未安装，请先运行 `install.packages("openxlsx")`。
- 如果metadata中缺少已经读取到的样本，脚本会自动把缺失样本行补回模板，并提示用户填写后继续。
- 合并样本后会输出 `合并后metadata列检查.csv`；如果某些用户metadata列没有进入 `sce@meta.data`，流程会警告并跳过预览，不会因为预览列缺失而中断。
- 逐样本QC会同时输出单个样本PDF和所有样本合并展示PDF，包括过滤前、过滤后以及过滤前后对比。
- 非交互式 `Rscript` 无法现场等待用户选择PC和resolution，需要提前在参数文件中填写 `final_pc_number` 和 `final_resolution`。
- 在RStudio交互运行时，如果输入了非候选resolution，程序会当场提示并重新询问；非法值不会被写入报错缓存。
- 交互式选择PC时必须输入PC梯度图中已经绘制的候选值；如果误输入非候选值，流程会重新提示，不会把错误PC写入断点缓存。
- 默认PC梯度已收窄为 `seq(10, 30, by = 5)`，默认resolution梯度已收窄为 `seq(0.2, 0.8, by = 0.2)`；需要更精细比较时可在配置文件中自行扩大。
