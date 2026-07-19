# ==============================================================================
# 功能块08：合并Seurat V5表达层并寻找每个cluster的marker
# 输入：完成正式聚类的sce
# 输出：markers、marker表格和热图
# ==============================================================================

####九、寻找每个cluster的marker基因####

# Seurat V5在FindAllMarkers前需要把各样本RNA layers重新合并。
# JoinLayers只合并表达矩阵层，不会删除Harmony降维结果和UMAP坐标。
sce[["RNA"]] <- JoinLayers(sce[["RNA"]])
DefaultAssay(sce) <- "RNA"
Idents(sce) <- "seurat_clusters"

markers <- FindAllMarkers(
  object = sce,
  assay = "RNA",
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)

write.csv(
  markers,
  file.path(result_dir, "所有cluster_marker.csv"),
  row.names = FALSE
)

# 每个cluster按照avg_log2FC选取前10个marker，用于快速查看cluster特征。
top10_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup()

write.csv(
  top10_markers,
  file.path(result_dir, "每个cluster_top10_marker.csv"),
  row.names = FALSE
)

# 热图使用RNA assay的scale.data，因此先对需要展示的marker进行ScaleData。
heatmap_genes <- unique(top10_markers$gene)
sce <- ScaleData(sce, assay = "RNA", features = heatmap_genes, verbose = FALSE)

p.heatmap <- DoHeatmap(
  sce,
  features = heatmap_genes,
  group.by = "seurat_clusters",
  raster = TRUE
) + NoLegend()

ggsave(
  filename = file.path(result_dir, "cluster_top10_marker热图.pdf"),
  plot = p.heatmap,
  width = 12,
  height = max(8, length(heatmap_genes) * 0.12)
)




