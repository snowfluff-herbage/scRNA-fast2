# ==============================================================================
# 流程参数：通常只需要修改这一个文件
# 所有模块都会使用这里设置的路径和参数
# ==============================================================================

####一、设置工作目录、载入R包和修改参数####

# Windows路径建议统一使用正斜杠“/”，避免反斜杠被R识别为转义符。
# 请把下面三个路径改成自己的真实路径。
project_dir <- "C:/singlecell/project"
raw_data_dir <- "C:/singlecell/project/rawdata"
result_dir <- "C:/singlecell/project/result"

# 自动建立项目目录和结果目录；raw_data_dir需要由用户提前准备。
dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
setwd(project_dir)

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 30 * 1024^3)  #允许future使用最多30GB内存
set.seed(1234)


# ----需要根据数据修改的参数----

# 输入数据类型。默认使用"auto"自动识别，也可以手动指定为：
#   "10x"：标准10X样本文件夹，或带共同样本前缀并平铺保存的三个矩阵文件；
#   "h5" ：一个或多个10x Genomics生成的.h5/.hdf5文件；
#   "csv"：一个或多个“基因×细胞”的.csv或.csv.gz表达矩阵。
# 如果raw_data_dir中同时放了多种格式，自动模式无法判断应读取哪一种，
# 此时请手动填写"10x"、"h5"或"csv"。
input_data_type <- "auto"

# 仅在10X目录第一次Read10X()读取失败时使用。
#   "ask"：在R控制台询问是否开始整理；
#   "yes"：不询问，直接整理；
#   "no" ：不修改任何文件，直接停止并显示读取错误。
# 整理时不再把样本命名为P1、P2、P3……，而是提取barcode、feature/gene、
# matrix前面的共同文件名前缀，以该前缀建立样本文件夹。例如：
#   GSM123_sampleA_barcodes.tsv.gz
#   GSM123_sampleA_features.tsv.gz
#   GSM123_sampleA_matrix.mtx.gz
# 会建立GSM123_sampleA文件夹，并把三个文件放入其中、改成10X标准文件名。
organize_10x_after_failure <- "ask"

# 用户填写的“样本级metadata”文件路径。
# 本流程统一使用XLSX格式，并且至少包含一列名为sample_id。
# sample_id需要与本脚本读取后得到的sample_name完全一致。
# 其他列可以由用户自由增加，例如group、batch、patient、tissue、condition等；
# 这些列都会自动添加到最终sce@meta.data中。
metadata_file <- file.path(project_dir, "样本metadata.xlsx")

# 建立Seurat对象时的最低要求。
create_min_cells <- 3
create_min_features <- 100

# 正式QC阈值。请结合QC小提琴图修改，不要机械照搬。
min_features <- 200
max_features <- 8000
min_counts <- 500
max_counts <- 100000
max_percent_mito <- 20
max_percent_hb <- 5
min_percent_ribo <- 0

# 标准化、高变基因、PCA参数。
n_variable_features <- 3000
n_pcs_calculate <- 50

# 是否使用HBC培训流程中的co1/co2规则自动选择聚类所用PC数。
# TRUE：PCA完成后计算推荐PC数，作为PC梯度比较的参考；
# FALSE：不计算co1/co2，使用manual_pc_number作为参考值。
auto_choose_pc <- TRUE
manual_pc_number <- 20

# PC数量梯度。脚本会逐一使用这些PC数量完成邻居图、测试聚类和UMAP并绘图。
# 超过实际PCA/Harmony维度的数值会自动删除，co1/co2推荐值会自动加入梯度。
pc_number_list <- seq(10, 30, by = 5)

# 比较不同PC数量时暂时固定使用同一个resolution，避免同时改变两个参数。
pc_test_resolution <- 0.5

# 正式PC数量。默认NA表示绘图后由用户选择：
#   在RStudio交互运行时，控制台会要求输入候选PC数；
#   使用Rscript非交互运行时会在生成PC比较图和中间对象后停止，
#   用户查看图片并在这里填入候选值，重新运行后才会继续比较resolution。
final_pc_number <- NA_integer_

# Harmony使用哪一列作为批次信息。
# 默认每个10X捕获样本是一个批次，因此使用sample_id。
harmony_group <- "sample_id"

# 聚类时先计算多个分辨率，绘图后再选择其中一个作为正式cluster。
resolution_list <- seq(0.2, 0.8, by = 0.2)
# 默认NA表示等待用户查看resolution比较图后输入；也可以预先填写0.5等候选值。
final_resolution <- NA_real_

# 是否启用“正式处理前最终参数覆盖区”。
# 默认FALSE，使用上面的初始值；改成TRUE后，下方第二处参数值将覆盖这里的设置。
use_parameter_override <- FALSE

# 流程报错后的恢复方式：
#   "auto"    ：默认。读取result目录中的“流程运行状态.rds”，自动跳过已经完成的
#               功能块，从上次报错的功能块重新开始；
#   "restart" ：忽略以前的状态和检查点，从输入读取重新开始完整分析。
#
# 修改原始数据、metadata、QC阈值、PCA或聚类参数后，建议先使用"restart"完整重跑。
# 只是修复包兼容问题、填写metadata、选择PC/resolution或填写注释编号时使用"auto"。
resume_mode <- "auto"
