# ==============================================================================
# 正式处理前的第二处参数覆盖区
# ==============================================================================
# 本文件由主脚本在“原始矩阵读取完成”和“CreateSeuratObject/QC开始”之间载入。
# 默认use_parameter_override=FALSE，因此仍使用01_流程参数.R中的初始值。
# 如需在看过输入文件信息后统一覆盖参数：
#   1. 在01_流程参数.R中设置use_parameter_override <- TRUE；
#   2. 修改下面的数值。

if (use_parameter_override) {
  # 建立Seurat对象时的最低要求。
  create_min_cells <- 3
  create_min_features <- 100

  # QC阈值。
  min_features <- 200
  max_features <- 8000
  min_counts <- 500
  max_counts <- 100000
  max_percent_mito <- 20
  max_percent_hb <- 5
  min_percent_ribo <- 0

  # 标准化、高变基因和PCA参数。
  n_variable_features <- 3000
  n_pcs_calculate <- 50
  auto_choose_pc <- TRUE
  manual_pc_number <- 20
  pc_number_list <- seq(10, 50, by = 5)
  pc_test_resolution <- 0.5
  final_pc_number <- NA_integer_

  # Harmony和聚类参数。
  harmony_group <- "sample_id"
  resolution_list <- seq(0.1, 1.0, by = 0.1)
  final_resolution <- NA_real_
}

