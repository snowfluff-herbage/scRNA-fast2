# ==============================================================================
# 功能块07：Harmony、PC梯度、resolution梯度、UMAP/tSNE和正式聚类
# 输入：完成PCA的sce、pc_number_candidates
# 输出：完成正式聚类的sce
# ==============================================================================

####七、Harmony去除样本或批次效应####

# 功能块07内部还包含Harmony、PC梯度、resolution计算和正式聚类等多个耗时步骤。
# 除了主入口的“功能块检查点”外，这里再保存更细的阶段检查点。
# 如果上次在clustree或后续绘图时报错，重新运行时可直接读取已经计算好的sce，
# 不需要重复Harmony、全部PC梯度和resolution聚类。
module_resume_enabled <- exists("resume_mode", inherits = TRUE) &&
  identical(tolower(resume_mode), "auto")

pc_gradient_checkpoint_file <- file.path(
  result_dir,
  "3.1_PC梯度比较完成_待选择PC_sce.rds"
)
pc_selection_checkpoint_file <- file.path(
  result_dir,
  "3.1.1_正式PC选择结果.rds"
)
resolution_calculation_checkpoint_file <- file.path(
  result_dir,
  "3.2_resolution计算完成_待绘图与选择_sce.rds"
)

resume_from_resolution_checkpoint <- module_resume_enabled &&
  file.exists(resolution_calculation_checkpoint_file)
resume_from_pc_checkpoint <- module_resume_enabled &&
  !resume_from_resolution_checkpoint &&
  file.exists(pc_gradient_checkpoint_file)

if (resume_from_resolution_checkpoint) {
  sce <- readRDS(resolution_calculation_checkpoint_file)
  message(
    "已读取resolution计算检查点，将跳过Harmony、PC梯度和resolution重新计算。"
  )
} else if (resume_from_pc_checkpoint) {
  sce <- readRDS(pc_gradient_checkpoint_file)
  message("已读取PC梯度检查点，将跳过Harmony和全部PC梯度重新计算。")
}

if (module_resume_enabled && file.exists(pc_selection_checkpoint_file)) {
  saved_pc_selection <- readRDS(pc_selection_checkpoint_file)
  saved_pc_number <- suppressWarnings(as.integer(saved_pc_selection$final_pc_number))

  # 旧检查点中的PC值必须仍然属于本次已经绘图的候选值。
  # 如果用户上次误输入了非候选PC，或后来修改了pc_number_list，就忽略旧值并重新选择。
  if (length(saved_pc_number) == 1L &&
      !is.na(saved_pc_number) &&
      saved_pc_number %in% pc_number_candidates) {
    final_pc_number <- saved_pc_number
    message("已恢复上次正式选择的PC数量：", final_pc_number)
  } else {
    final_pc_number <- NA_integer_
    unlink(pc_selection_checkpoint_file)
    message(
      "上次保存的正式PC数量不是当前候选值，已清除旧PC选择，请重新选择。"
    )
  }
}

# Harmony根据PCA低维空间校正批次，不会直接改写RNA counts和data。
# harmony_group指定的metadata列必须存在，并且每个细胞都有对应值。
# 为了能够比较完整PC梯度，Harmony先使用候选梯度中的最大PC数；
# 后面的每次PC测试和正式聚类再分别使用1:当前候选PC数。
if (!resume_from_pc_checkpoint && !resume_from_resolution_checkpoint) {

harmony_pc_num <- seq_len(max(pc_number_candidates))
sce <- RunHarmony(
  object = sce,
  group.by.vars = harmony_group,
  reduction.use = "pca",
  dims.use = harmony_pc_num,
  reduction.save = "harmony",
  verbose = TRUE
)

p.harmony <- DimPlot(sce, reduction = "harmony", group.by = harmony_group)
ggsave(
  filename = file.path(result_dir, "Harmony_按批次.pdf"),
  plot = p.harmony,
  width = 8,
  height = 6
)

saveRDS(sce, file.path(result_dir, "3.标准化_PCA_Harmony后_sce.rds"))



####八、选择PC、聚类分辨率和UMAP/tSNE####

######第一步：按照PC数量梯度分别聚类和绘图######

# PC比较阶段只改变PC数量，其他参数保持一致：
#   reduction固定为harmony；
#   resolution固定为pc_test_resolution；
#   UMAP的n.neighbors、min.dist和随机种子保持一致。
# 这样不同图片之间的主要区别来自纳入的PC数量，而不是resolution同时发生变化。

pc_gradient_plots <- list()
pc_gradient_summary <- data.frame()

for (i in seq_along(pc_number_candidates)) {
  current_pc_number <- pc_number_candidates[i]
  current_pc_dims <- seq_len(current_pc_number)

  message(
    "正在进行PC梯度比较：使用前", current_pc_number,
    "个Harmony维度，测试resolution = ", pc_test_resolution
  )

  # 每次使用临时对象完成邻居图、测试聚类和UMAP。
  # 这样不会在正式sce中同时保存大量测试图结构和UMAP结果，减少对象体积。
  sce_pc_test <- sce
  sce_pc_test <- FindNeighbors(
    sce_pc_test,
    reduction = "harmony",
    dims = current_pc_dims,
    verbose = FALSE
  )
  sce_pc_test <- FindClusters(
    sce_pc_test,
    resolution = pc_test_resolution,
    algorithm = 1,
    random.seed = 1234,
    verbose = FALSE
  )
  sce_pc_test <- RunUMAP(
    sce_pc_test,
    reduction = "harmony",
    dims = current_pc_dims,
    n.neighbors = 30,
    min.dist = 0.3,
    seed.use = 1234,
    verbose = FALSE
  )

  current_cluster_number <- length(unique(sce_pc_test$seurat_clusters))
  current_title <- paste0(
    "PC = 1:", current_pc_number,
    " | resolution = ", pc_test_resolution,
    " | clusters = ", current_cluster_number
  )
  if (current_pc_number == recommended_pc_number) {
    current_title <- paste0(current_title, " | co1/co2推荐")
  }

  pc_gradient_plots[[i]] <- DimPlot(
    sce_pc_test,
    reduction = "umap",
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE,
    shuffle = TRUE
  ) + ggtitle(current_title) + NoLegend()

  # 每个候选PC单独保存一张图，方便放大查看cluster边界和小群。
  ggsave(
    filename = file.path(
      result_dir,
      paste0("PC梯度_PC1-", current_pc_number, "_UMAP.pdf")
    ),
    plot = pc_gradient_plots[[i]],
    width = 8,
    height = 7
  )

  pc_gradient_summary <- rbind(
    pc_gradient_summary,
    data.frame(
      PC_number = current_pc_number,
      test_resolution = pc_test_resolution,
      cluster_number = current_cluster_number,
      recommended_by_co1_co2 = current_pc_number == recommended_pc_number
    )
  )

  rm(sce_pc_test)
  invisible(gc())
}

# 把所有PC梯度图合并到同一个PDF中，便于横向比较。
p.pc.gradient <- wrap_plots(pc_gradient_plots, ncol = 3)
ggsave(
  filename = file.path(result_dir, "PC梯度_全部UMAP比较.pdf"),
  plot = p.pc.gradient,
  width = 15,
  height = 5 * ceiling(length(pc_gradient_plots) / 3)
)

write.csv(
  pc_gradient_summary,
  file.path(result_dir, "PC梯度_聚类数量汇总.csv"),
  row.names = FALSE
)

# 保存绘图完成、尚未正式选择PC的中间对象。
# 该对象已经包含PCA和最大候选维度的Harmony结果，但不包含测试聚类产生的大量临时图。
saveRDS(
  sce,
  file.path(result_dir, "3.1_PC梯度比较完成_待选择PC_sce.rds")
)

} else {
  message("Harmony和PC梯度已经有可用检查点，本次不重复运行。")
}


######用户确定正式PC数量######

if (is.na(final_pc_number)) {
  if (interactive()) {
    repeat {
      pc_answer <- trimws(readline(paste0(
        "请查看result目录中的PC梯度图，并输入正式PC数量（候选值：",
        paste(pc_number_candidates, collapse = "、"),
        "；直接回车使用co1/co2推荐值", recommended_pc_number, "）："
      )))

      if (pc_answer == "") {
        selected_pc_number <- recommended_pc_number
      } else {
        selected_pc_number <- suppressWarnings(as.integer(pc_answer))
      }

      if (!is.na(selected_pc_number) &&
          selected_pc_number %in% pc_number_candidates) {
        final_pc_number <- selected_pc_number
        break
      }

      message(
        "输入的PC数量不是候选值，请重新输入以下候选值之一：",
        paste(pc_number_candidates, collapse = "、")
      )
    }
  } else {
    stop(
      "PC梯度图已经全部生成。当前为非交互式运行，无法等待用户输入。",
      "请查看PC梯度_全部UMAP比较.pdf，将final_pc_number填写为以下候选值之一后重新运行：",
      paste(pc_number_candidates, collapse = "、")
    )
  }
}

if (is.na(final_pc_number) || !final_pc_number %in% pc_number_candidates) {
  final_pc_number <- NA_integer_
  unlink(pc_selection_checkpoint_file)
  stop(
    "final_pc_number必须是本次已经绘图的候选值之一：",
    paste(pc_number_candidates, collapse = "、")
  )
}

pc.num <- seq_len(final_pc_number)
message("正式聚类使用PC：1:", final_pc_number)

# 立即保存用户选择的PC数量。后面即使resolution计算或绘图报错，重新运行时也不再询问。
saveRDS(
  list(
    final_pc_number = final_pc_number,
    pc.num = pc.num
  ),
  pc_selection_checkpoint_file
)


######第二步：固定正式PC数量后比较resolution######

# FindNeighbors只需要运行一次；不同resolution共用同一个邻居图。
# 一次输入多个resolution后，metadata会生成RNA_snn_res.0.1、0.2……等列。
if (!resume_from_resolution_checkpoint) {

sce <- FindNeighbors(
  sce,
  reduction = "harmony",
  dims = pc.num
)

sce <- FindClusters(
  sce,
  resolution = resolution_list,
  algorithm = 1,
  random.seed = 1234
)

# UMAP只依赖选择的reduction和PC，不依赖cluster分辨率，因此只需要计算一次。
sce <- RunUMAP(
  sce,
  reduction = "harmony",
  dims = pc.num,
  n.neighbors = 30,
  min.dist = 0.3,
  seed.use = 1234
)

sce <- RunTSNE(
  sce,
  reduction = "harmony",
  dims = pc.num,
  seed.use = 1234
)

# resolution列、UMAP和tSNE一旦计算完成就立即保存。
# 该检查点必须放在clustree之前，因为clustree只是可选图，不能阻断聚类结果保存。
saveRDS(sce, resolution_calculation_checkpoint_file)

} else {
  message("resolution、UMAP和tSNE已经从检查点恢复，本次不重复计算。")
}


####比较不同聚类分辨率####

# 如果安装了clustree，可以直接观察不同resolution之间cluster的分裂关系。
# clustree/ggraph的部分版本暂时不兼容ggplot2 4.x，可能出现：
# Unknown guide: edge_colourbar。
# clustree只是辅助图，因此使用tryCatch；绘图失败只警告，不再终止整个流程。
if (requireNamespace("clustree", quietly = TRUE)) {
  tryCatch({
    p.clustree <- clustree::clustree(
      sce@meta.data,
      prefix = "RNA_snn_res."
    )
    ggsave(
      filename = file.path(result_dir, "不同分辨率_clustree.pdf"),
      plot = p.clustree,
      width = 12,
      height = 9
    )
  }, error = function(e) {
    warning(
      "clustree绘图失败，但不影响已经完成的聚类、UMAP和tSNE。\n",
      "具体原因：", conditionMessage(e), "\n",
      "请继续查看“不同聚类分辨率_UMAP.pdf”选择resolution。"
    )
  })
} else {
  message("没有安装clustree，已跳过clustree；resolution UMAP仍会正常生成。")
}

# 分别画出不同resolution的UMAP，便于人工选择正式分辨率。
resolution_plots <- list()
for (i in seq_along(resolution_list)) {
  resolution_column <- paste0("RNA_snn_res.", resolution_list[i])
  resolution_plots[[i]] <- DimPlot(
    sce,
    reduction = "umap",
    group.by = resolution_column,
    label = TRUE,
    repel = TRUE,
    shuffle = TRUE
  ) + ggtitle(paste0("resolution = ", resolution_list[i])) + NoLegend()
}

p.resolution <- wrap_plots(resolution_plots, ncol = 3)
ggsave(
  filename = file.path(result_dir, "不同聚类分辨率_UMAP.pdf"),
  plot = p.resolution,
  width = 15,
  height = 5 * ceiling(length(resolution_plots) / 3)
)

# 保存已经固定正式PC并完成resolution梯度计算的中间对象。
saveRDS(
  sce,
  file.path(result_dir, "3.2_resolution比较完成_待选择resolution_sce.rds")
)


#####正式选择聚类分辨率####

# PC数量确定后，才让用户根据UMAP和clustree选择正式resolution。
# 交互运行时，如果输入了非候选值、空值或文字，会立即重新询问，
# 不会让错误值进入断点缓存，也不需要重新Source整个流程。
resolution_value <- suppressWarnings(as.numeric(final_resolution))
resolution_match <- integer()
if (length(resolution_value) == 1L &&
    !is.na(resolution_value) &&
    is.finite(resolution_value)) {
  # 使用数值差值匹配，避免小数在计算机内部的表示误差。
  resolution_match <- which(
    abs(resolution_list - resolution_value) < 1e-10
  )
}
resolution_is_valid <- length(resolution_match) == 1L

if (!resolution_is_valid) {
  # 先清空参数文件或旧错误缓存中的非法值，再进入交互选择。
  final_resolution <- NA_real_

  if (interactive()) {
    repeat {
      resolution_answer <- trimws(readline(paste0(
        "请查看不同聚类分辨率_UMAP.pdf和clustree，并输入正式resolution（候选值：",
        paste(resolution_list, collapse = "、"), "）："
      )))
      resolution_value <- suppressWarnings(as.numeric(resolution_answer))
      resolution_match <- integer()

      if (length(resolution_value) == 1L &&
          !is.na(resolution_value) &&
          is.finite(resolution_value)) {
        resolution_match <- which(
          abs(resolution_list - resolution_value) < 1e-10
        )
      }

      if (length(resolution_match) == 1L) {
        final_resolution <- resolution_list[resolution_match]
        break
      }

      message(
        "输入的resolution不在候选值中，请重新输入：",
        paste(resolution_list, collapse = "、")
      )
    }
  } else {
    stop(
      "resolution梯度图已经生成。当前为非交互式运行，无法等待用户输入。",
      "请查看不同聚类分辨率_UMAP.pdf和clustree，将final_resolution填写为以下候选值之一后重新运行：",
      paste(resolution_list, collapse = "、")
    )
  }
} else {
  # 对来自参数文件或合法断点缓存的值做标准化。
  final_resolution <- resolution_list[resolution_match]
}

# 保存正式resolution选择，便于错误恢复和结果追溯。
saveRDS(
  list(
    final_pc_number = final_pc_number,
    final_resolution = final_resolution
  ),
  file.path(result_dir, "3.2.1_正式resolution选择结果.rds")
)

final_cluster_column <- paste0("RNA_snn_res.", final_resolution)
sce$seurat_clusters <- as.character(sce@meta.data[[final_cluster_column]])
Idents(sce) <- "seurat_clusters"

table(sce$seurat_clusters)

write.csv(
  data.frame(
    final_PC_number = final_pc_number,
    final_PC_range = paste0("1:", final_pc_number),
    final_resolution = final_resolution,
    clustering_algorithm = "Louvain_algorithm_1"
  ),
  file.path(result_dir, "正式聚类参数.csv"),
  row.names = FALSE
)

p.cluster <- DimPlot(
  sce,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  shuffle = TRUE,
  raster = FALSE
) + NoLegend()

p.sample <- DimPlot(
  sce,
  reduction = "umap",
  group.by = "sample_id",
  raster = FALSE
)

p.cluster.sample <- p.cluster + p.sample
ggsave(
  filename = file.path(result_dir, "正式聚类_UMAP.pdf"),
  plot = p.cluster.sample,
  width = 15,
  height = 7
)

p.tsne <- DimPlot(
  sce,
  reduction = "tsne",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + NoLegend()
ggsave(
  filename = file.path(result_dir, "正式聚类_tSNE.pdf"),
  plot = p.tsne,
  width = 8,
  height = 7
)

saveRDS(sce, file.path(result_dir, "4.正式聚类后_sce.rds"))

