# ==============================================================================
# 功能块06：NormalizeData、高变基因、细胞周期、ScaleData和PCA
# 输入：sce
# 输出：完成PCA的sce以及PC候选梯度
# ==============================================================================

####六、标准化、高变基因、细胞周期评分和PCA####

DefaultAssay(sce) <- "RNA"

# merge后的Seurat V5对象通常为每个样本保留独立counts layer。
# 当前脚本使用RunHarmony直接校正PCA，不使用IntegrateLayers，因此这里先把RNA layers
# 合并为统一的counts层。这样CellCycleScoring、AddModuleScore和后续marker分析都能
# 使用明确的RNA data层，避免GetAssayData面对多个data layers时产生歧义。
sce[["RNA"]] <- JoinLayers(sce[["RNA"]])

# NormalizeData把原始UMI counts进行文库大小校正和log1p转换。
# JoinLayers后在统一的RNA counts层上进行标准化；批次信息仍保留在metadata中，
# 后面由Harmony根据harmony_group完成低维空间校正。
sce <- NormalizeData(
  sce,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


####细胞周期评分####

# 细胞周期评分不是必须删除细胞周期效应。
# 这里先计算S.Score、G2M.Score和Phase，供后续判断是否需要回归。
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
s.genes <- intersect(s.genes, rownames(sce))
g2m.genes <- intersect(g2m.genes, rownames(sce))

sce <- CellCycleScoring(
  object = sce,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE
)
sce$CC.Difference <- sce$S.Score - sce$G2M.Score

p.cycle <- VlnPlot(
  sce,
  features = c("S.Score", "G2M.Score"),
  group.by = "sample_id",
  ncol = 2,
  pt.size = 0
)
ggsave(
  filename = file.path(result_dir, "细胞周期评分.pdf"),
  plot = p.cycle,
  width = 10,
  height = 5
)


####高变基因、ScaleData和PCA####

sce <- FindVariableFeatures(
  sce,
  selection.method = "vst",
  nfeatures = n_variable_features
)

head(VariableFeatures(sce), 20)

p.variable <- VariableFeaturePlot(sce)
ggsave(
  filename = file.path(result_dir, "高变基因.pdf"),
  plot = p.variable,
  width = 8,
  height = 6
)

# 按照原脚本风格，对高变基因进行ScaleData。
# 这里不默认回归percent_mito或细胞周期，避免在未检查数据前过度校正生物学信号。
# 如果后续确认线粒体比例或细胞周期明显主导PCA，再自行添加vars.to.regress。
sce <- ScaleData(sce, features = VariableFeatures(sce))

sce <- RunPCA(
  sce,
  features = VariableFeatures(sce),
  npcs = n_pcs_calculate,
  verbose = FALSE
)

p.pca <- DimPlot(sce, reduction = "pca", group.by = "sample_id")
ggsave(
  filename = file.path(result_dir, "PCA_按样本.pdf"),
  plot = p.pca,
  width = 8,
  height = 6
)

p.elbow <- ElbowPlot(sce, ndims = n_pcs_calculate)
ggsave(
  filename = file.path(result_dir, "PCA_ElbowPlot.pdf"),
  plot = p.elbow,
  width = 7,
  height = 5
)


######按照co1/co2规则自动选择聚类使用的PC数######

if (auto_choose_pc) {
  recommended_pc_number <- choose_pc(sce)
} else {
  # 关闭co1/co2推荐时，使用manual_pc_number作为PC梯度的参考值。
  recommended_pc_number <- min(
    manual_pc_number,
    length(sce[["pca"]]@stdev)
  )
  message("已关闭co1/co2推荐，参考PC数为：", recommended_pc_number)
}

# 整理PC梯度：删除小于1或超过实际PCA维数的值，并自动加入co1/co2推荐值。
# 如果用户已经预先填写final_pc_number，也把它加入梯度并绘图，便于核对。
actual_pca_number <- length(sce[["pca"]]@stdev)
if (!is.na(final_pc_number) &&
    (final_pc_number < 1 || final_pc_number > actual_pca_number)) {
  stop("final_pc_number必须位于1到", actual_pca_number, "之间。")
}

pc_number_candidates <- c(pc_number_list, recommended_pc_number)
if (!is.na(final_pc_number)) {
  pc_number_candidates <- c(pc_number_candidates, final_pc_number)
}
pc_number_candidates <- sort(unique(as.integer(pc_number_candidates)))
pc_number_candidates <- pc_number_candidates[
  pc_number_candidates >= 1 & pc_number_candidates <= actual_pca_number
]

if (length(pc_number_candidates) == 0) {
  stop("PC梯度中没有可用数值，请检查pc_number_list和n_pcs_calculate。")
}

message(
  "本次将比较以下PC数量：",
  paste(pc_number_candidates, collapse = "、")
)

# 保存每个PC的贡献比例、累计比例、是否为推荐位置以及是否进入比较梯度。
pc_pct <- sce[["pca"]]@stdev / sum(sce[["pca"]]@stdev) * 100
pc_selection_summary <- data.frame(
  PC = seq_along(pc_pct),
  percent = pc_pct,
  cumulative_percent = cumsum(pc_pct),
  recommended_by_co1_co2 = seq_along(pc_pct) == recommended_pc_number,
  included_in_PC_gradient = seq_along(pc_pct) %in% pc_number_candidates
)
write.csv(
  pc_selection_summary,
  file.path(result_dir, "PC自动选择结果.csv"),
  row.names = FALSE
)




