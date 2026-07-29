# ==============================================================================
# 功能块09：按照原脚本的cluster向量方法完成大亚群注释
# 注释编号从config/02_大亚群注释参数.R读取
# ==============================================================================

####十、根据经典marker完成大亚群注释####

#####1.使用经典marker检查每个cluster####

# 以下marker只用于大亚群级别判断，不直接决定最终注释。
# 正式注释时需要同时结合FindAllMarkers结果、DotPlot、FeaturePlot和生物学背景。
genes_to_check <- c(
  # Epithelial
  "EPCAM", "KRT8", "KRT18", "KRT19",

  # Endothelial
  "PECAM1", "VWF", "CLDN5", "EMCN",

  # Fibroblast / stromal
  "COL1A1", "COL1A2", "DCN", "COL3A1",

  # T cells
  "CD3D", "CD3E", "TRAC", "CD4", "CD8A", "CD8B",

  # NK cells
  "NKG7", "GNLY", "KLRD1", "NCR1",

  # B cells
  "CD79A", "CD19", "MS4A1", "CD37",

  # Plasma cells
  "MZB1", "JCHAIN", "SDC1", "IGHG1",

  # Myeloid cells
  "LST1", "TYROBP", "FCER1G", "CTSS", "C1QA", "S100A8", "S100A9",

  # Dendritic cells
  "CD1C", "CLEC9A", "LILRA4", "LAMP3",

  # Mast cells
  "TPSAB1", "TPSB2", "CPA3", "MS4A2"
)

# 删除数据中不存在的基因，避免DotPlot因为个别基因缺失而警告。
genes_to_check <- intersect(genes_to_check, rownames(sce))

p.marker.dot <- DotPlot(
  sce,
  features = genes_to_check,
  assay = "RNA",
  group.by = "seurat_clusters"
) +
  coord_flip() +
  theme_bw() +
  ggtitle("大亚群经典marker")

ggsave(
  filename = file.path(result_dir, "大亚群经典marker_DotPlot.pdf"),
  plot = p.marker.dot,
  width = 12,
  height = max(8, length(genes_to_check) * 0.22)
)


#####2.按照原脚本方法直接进行大亚群注释####

# 先查看cluster编号和每个cluster的细胞数量。
# 然后结合上面的DotPlot、top10 marker热图和所有cluster_marker.csv，
# 把相应的cluster编号填入注释配置文件的celltype_cluster_map中。
table(sce@meta.data$seurat_clusters)

# 例如某一类细胞对应cluster 0、3、7，可以在配置文件中写成：
# "自定义亚群名称" = c(0, 3, 7)
# 亚群名称和cluster编号统一从config/02_大亚群注释参数.R读取。

# 如果尚未填写任何cluster，停在这里，避免生成全部为unknown的最终对象。
# 此时marker、热图和经典marker DotPlot已经生成，可以据此完成注释参数。
if (!exists("celltype_cluster_map", inherits = TRUE)) {
  stop(
    "没有找到celltype_cluster_map。请使用新版",
    "config/02_大亚群注释参数.R填写大亚群名称和cluster编号。"
  )
}

# celltype_cluster_map必须是有名称的list：
# 列表名称就是最终写入sce$celltype的大亚群名称，列表内容是对应cluster编号。
if (!is.list(celltype_cluster_map) || is.null(names(celltype_cluster_map))) {
  stop(
    "celltype_cluster_map必须是带名称的list，例如：",
    "list(\"T/NK\" = c(0, 3), \"Myeloid\" = c(1, 5))。"
  )
}

# 去除名称两侧可能误输入的空格；亚群名称不能留空或重复。
annotation_names <- trimws(names(celltype_cluster_map))
if (any(!nzchar(annotation_names))) {
  stop("celltype_cluster_map中存在空白的大亚群名称，请补充名称或删除该项。")
}
if (anyDuplicated(annotation_names)) {
  duplicated_names <- unique(annotation_names[duplicated(annotation_names)])
  stop(
    "celltype_cluster_map中存在重复的大亚群名称：",
    paste(duplicated_names, collapse = "、")
  )
}
names(celltype_cluster_map) <- annotation_names

# cluster编号统一转换为字符，确保数值型0和因子/字符型"0"均可正常匹配。
# 空的c()可以保留在配置文件中；正式赋值时会自动忽略这些空项。
annotation_cluster_list <- lapply(
  celltype_cluster_map,
  function(x) unique(as.character(x[!is.na(x)]))
)
annotation_cluster_list <- annotation_cluster_list[
  lengths(annotation_cluster_list) > 0
]

if (length(annotation_cluster_list) == 0) {
  stop(
    "marker和经典marker图已经生成，但尚未填写大亚群cluster编号。\n",
    "请修改config/02_大亚群注释参数.R，然后运行",
    "“辅助_仅重新进行大亚群注释.R”，不需要重复前面的分析。"
  )
}

# 同一个cluster如果被填入两个大亚群，后写入的类型会覆盖前面的类型。
# 为避免这种不易察觉的错误，在正式赋值前直接检查重复编号。
all_annotated_clusters <- unlist(annotation_cluster_list, use.names = FALSE)
duplicated_clusters <- unique(
  all_annotated_clusters[duplicated(all_annotated_clusters)]
)
if (length(duplicated_clusters) > 0) {
  stop(
    "以下cluster被重复填写到多个大亚群：",
    paste(duplicated_clusters, collapse = "、"),
    "。请修改config/02_大亚群注释参数.R。"
  )
}

existing_clusters <- unique(as.character(sce$seurat_clusters))
invalid_clusters <- setdiff(as.character(all_annotated_clusters), existing_clusters)
if (length(invalid_clusters) > 0) {
  stop(
    "注释参数中包含当前对象不存在的cluster：",
    paste(invalid_clusters, collapse = "、")
  )
}

# 保存配置文件中本次实际填写的关系，方便检查自定义名称和cluster编号。
# 该CSV仅用于记录，不参与注释，也不需要用户人工填写后再读回。
configured_annotation <- do.call(
  rbind,
  lapply(names(annotation_cluster_list), function(current_celltype) {
    data.frame(
      celltype = current_celltype,
      ClusterID = annotation_cluster_list[[current_celltype]],
      stringsAsFactors = FALSE
    )
  })
)
rownames(configured_annotation) <- NULL
write.csv(
  configured_annotation,
  file.path(result_dir, "用户配置的大亚群与cluster关系.csv"),
  row.names = FALSE
)

# 建立cluster与celltype的对应表；所有cluster初始均为unknown。
# 仍然沿用原脚本的直接赋值思路，只是改为自动遍历用户配置的全部名称。
celltype <- data.frame(
  ClusterID = sort(unique(as.character(sce$seurat_clusters))),
  celltype = "unknown",
  stringsAsFactors = FALSE
)

for (current_celltype in names(annotation_cluster_list)) {
  current_clusters <- annotation_cluster_list[[current_celltype]]
  celltype[
    celltype$ClusterID %in% current_clusters,
    "celltype"
  ] <- current_celltype
}

celltype
table(celltype$celltype)


#####3.将大亚群注释写入Seurat对象####

# 先在metadata中增加celltype列，所有细胞初始标记为unknown。
sce@meta.data$celltype = "unknown"

# 按照原脚本方法逐行读取celltype对应表，根据seurat_clusters给每个细胞赋值。
for(i in 1:nrow(celltype)){
  sce@meta.data[
    which(sce@meta.data$seurat_clusters == celltype$ClusterID[i]),
    'celltype'
  ] <- celltype$celltype[i]
}

# 保存本次代码中实际使用的cluster与大亚群对应关系，方便复核；
# 该文件只是结果记录，不参与注释，也不需要人工填表后再读入。
write.csv(
  celltype,
  file.path(result_dir, "大亚群注释对应关系.csv"),
  row.names = FALSE
)

table(sce$celltype)
table(sce$celltype, sce$seurat_clusters)
Idents(sce) <- "celltype"


#####4.绘制大亚群注释结果####

p.celltype <- DimPlot(
  sce,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE,
  shuffle = TRUE,
  raster = FALSE
) + NoLegend()

p.celltype.sample <- DimPlot(
  sce,
  reduction = "umap",
  group.by = "celltype",
  split.by = "sample_id",
  label = FALSE,
  raster = FALSE
)

ggsave(
  filename = file.path(result_dir, "大亚群注释_UMAP.pdf"),
  plot = p.celltype,
  width = 10,
  height = 8
)

ggsave(
  filename = file.path(result_dir, "各样本大亚群注释_UMAP.pdf"),
  plot = p.celltype.sample,
  width = max(12, 4 * length(unique(sce$sample_id))),
  height = 6
)

p.marker.celltype <- DotPlot(
  sce,
  features = genes_to_check,
  assay = "RNA",
  group.by = "celltype"
) +
  coord_flip() +
  theme_bw() +
  ggtitle("大亚群注释后marker复核")

ggsave(
  filename = file.path(result_dir, "大亚群注释后_marker复核.pdf"),
  plot = p.marker.celltype,
  width = 12,
  height = max(8, length(genes_to_check) * 0.22)
)


#####5.保存大亚群注释后的对象####

saveRDS(sce, file.path(result_dir, "5.大亚群注释完成_sce.rds"))
write.csv(sce@meta.data, file.path(result_dir, "大亚群注释后_metadata.csv"))
writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))

cat(
  "\n流程已运行到大亚群注释。\n",
  "Seurat对象：", file.path(result_dir, "5.大亚群注释完成_sce.rds"), "\n",
  "细胞数量：", ncol(sce), "\n",
  "基因数量：", nrow(sce), "\n",
  "正式PC范围：1:", final_pc_number, "\n",
  "正式resolution：", final_resolution, "\n",
  sep = ""
)


# ==============================================================================
# 脚本到此结束：暂不进行细胞亚群重聚类和其他下游分析
# ==============================================================================
