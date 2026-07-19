# ==============================================================================
# 功能块05：合并全部样本并核对metadata
# 输入：scRNAlist
# 输出：sce
# ==============================================================================

####五、合并样本和检查metadata####

# 使用merge合并多个Seurat对象。
# Seurat V5会把不同样本保留为不同RNA layers，后面可以直接进行标准化和Harmony。
if (length(scRNAlist) == 1) {
  sce <- scRNAlist[[1]]
} else {
  sce <- merge(x = scRNAlist[[1]], y = scRNAlist[-1], project = "scRNA_project")
}

# 用户填写的每一列已经复制到每个细胞的sce@meta.data中。
# 同时在misc中保存原始的“一行一个样本”表，便于以后直接查看样本设计信息。
sce@misc$sample_metadata <- sample_info

rm(scRNAlist)

dim(sce)
if (!"sample_id" %in% colnames(sce@meta.data)) {
  stop(
    "合并后的sce@meta.data中没有sample_id列，无法继续进行样本分组检查。\n",
    "请检查功能块03中metadata是否成功写入每个Seurat对象。"
  )
}
table(sce$sample_id)

# 显示用户metadata中哪些列已经进入合并后的sce。
user_metadata_columns <- colnames(sample_info)
user_metadata_columns

metadata_columns_found <- intersect(user_metadata_columns, colnames(sce@meta.data))
metadata_columns_missing <- setdiff(user_metadata_columns, colnames(sce@meta.data))

metadata_column_check <- data.frame(
  metadata_column = user_metadata_columns,
  in_sce_meta_data = user_metadata_columns %in% colnames(sce@meta.data)
)
data.table::fwrite(
  metadata_column_check,
  file.path(result_dir, "合并后metadata列检查.csv"),
  bom = TRUE
)

if (length(metadata_columns_missing) > 0) {
  warning(
    "以下用户metadata列没有在合并后的sce@meta.data中找到，",
    "本次预览时已跳过这些列：",
    paste(metadata_columns_missing, collapse = "、"),
    "\n已经输出检查表：",
    file.path(result_dir, "合并后metadata列检查.csv")
  )
}

if (length(metadata_columns_found) > 0) {
  head(sce@meta.data[, metadata_columns_found, drop = FALSE])
}

# group是常用但不是强制列；用户提供时再显示各组细胞数量。
if ("group" %in% colnames(sce@meta.data)) {
  table(sce$group)
}

# nCount与nFeature相关性图可辅助观察异常高UMI细胞。
p.scatter <- FeatureScatter(
  sce,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  group.by = "sample_id"
)
ggsave(
  filename = file.path(result_dir, "QC_nCount与nFeature相关性.pdf"),
  plot = p.scatter,
  width = 8,
  height = 6
)

saveRDS(sce, file.path(result_dir, "2.合并与QC后_sce.rds"))



