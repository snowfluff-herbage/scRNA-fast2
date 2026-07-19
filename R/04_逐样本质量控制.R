# ==============================================================================
# 功能块04：逐样本计算线粒体、核糖体和血红蛋白比例并完成QC
# 输入：scRNAlist及QC参数
# 输出：过滤后的scRNAlist、QC图和汇总表
# ==============================================================================

####四、逐样本计算QC指标和过滤低质量细胞####

# 推荐逐样本QC，而不是先合并再使用同一个阈值。
# 如果不同样本质量差异很大，可以把下面固定阈值改成每个样本独立阈值。

qc_summary <- data.frame()
mito_summary <- data.frame()
qc_before_plot_list <- list()
qc_after_plot_list <- list()
qc_compare_plot_list <- list()

for (i in seq_along(scRNAlist)) {
  current_sample <- names(scRNAlist)[i]

  # 分别统计人类常见的“MT-”和小鼠常见的“mt-”线粒体基因命名。
  # 每个样本独立判断，哪一种匹配到的基因更多，就使用哪一种计算percent_mito。
  mito_genes_upper <- grep("^MT-", rownames(scRNAlist[[i]]), value = TRUE)
  mito_genes_lower <- grep("^mt-", rownames(scRNAlist[[i]]), value = TRUE)

  if (length(mito_genes_upper) >= length(mito_genes_lower) &&
      length(mito_genes_upper) > 0) {
    mito_genes <- mito_genes_upper
    mito_gene_style <- "MT-（人类常见命名）"
  } else if (length(mito_genes_lower) > 0) {
    mito_genes <- mito_genes_lower
    mito_gene_style <- "mt-（小鼠常见命名）"
  } else {
    mito_genes <- character()
    mito_gene_style <- "未识别到MT-/mt-线粒体基因"
  }

  if (length(mito_genes) > 0) {
    scRNAlist[[i]] <- PercentageFeatureSet(
      scRNAlist[[i]],
      features = mito_genes,
      col.name = "percent_mito"
    )
  } else {
    # 不再把percent_mito强行写成0。
    # 如果gene.column读错，或表达矩阵使用Ensembl ID，把它写成0会让线粒体QC失效，
    # 而且后续结果表面上仍能继续运行，因此这里直接停止并展示实际基因行名。
    stop(
      "样本 ", current_sample,
      " 没有识别到MT-/mt-线粒体基因，已停止QC，避免把percent_mito错误设为0。\n",
      "当前对象前20个基因名：",
      paste(head(rownames(scRNAlist[[i]]), 20), collapse = "、"),
      "\n请先查看结果目录中的“10X基因名称列自动判断结果.csv”。",
      "如果这些行名是Gene Expression，说明features文件的基因列选择仍不正确；",
      "如果是ENSG开头的Ensembl ID，则需要先把Ensembl ID转换为gene symbol。"
    )
  }

  # 核糖体和血红蛋白基因比例暂时沿用人类基因符号写法。
  scRNAlist[[i]] <- PercentageFeatureSet(
    scRNAlist[[i]], pattern = "^RP[SL]", col.name = "percent_ribo"
  )
  scRNAlist[[i]] <- PercentageFeatureSet(
    scRNAlist[[i]], pattern = "^HB[ABDEGMQZ]", col.name = "percent_hb"
  )

  # 在过滤细胞之前汇总每个样本的线粒体比例分布。
  # percent_cells_over_threshold表示有多少比例的细胞超过前面设置的max_percent_mito。
  percent_cells_over_threshold <- mean(
    scRNAlist[[i]]$percent_mito > max_percent_mito
  ) * 100

  if (median(scRNAlist[[i]]$percent_mito) > max_percent_mito) {
    mito_judgement <- "样本中位线粒体比例高于QC阈值，建议重点检查"
  } else if (percent_cells_over_threshold > 0) {
    mito_judgement <- "样本中位数未超阈值，但部分细胞超过QC阈值"
  } else {
    mito_judgement <- "当前样本没有细胞超过线粒体QC阈值"
  }

  mito_summary <- rbind(
    mito_summary,
    data.frame(
      sample_id = current_sample,
      mitochondrial_gene_style = mito_gene_style,
      mitochondrial_gene_number = length(mito_genes),
      mean_percent_mito = mean(scRNAlist[[i]]$percent_mito),
      median_percent_mito = median(scRNAlist[[i]]$percent_mito),
      Q1_percent_mito = as.numeric(quantile(scRNAlist[[i]]$percent_mito, 0.25)),
      Q3_percent_mito = as.numeric(quantile(scRNAlist[[i]]$percent_mito, 0.75)),
      max_percent_mito = max(scRNAlist[[i]]$percent_mito),
      QC_mito_threshold = max_percent_mito,
      percent_cells_over_threshold = percent_cells_over_threshold,
      judgement = mito_judgement
    )
  )

  # 先画过滤前QC图，再根据分布调整前面设置的阈值。
  p.qc.before <- VlnPlot(
    scRNAlist[[i]],
    features = c(
      "nFeature_RNA", "nCount_RNA",
      "percent_mito", "percent_ribo", "percent_hb"
    ),
    # 此时尚未运行NormalizeData()，RNA assay只有原始counts层。
    # 明确指定layer="counts"，可避免Seurat反复寻找不存在的data层并产生警告。
    assay = "RNA",
    layer = "counts",
    pt.size = 0,
    ncol = 5
  ) + plot_annotation(title = paste0(current_sample, "_QC过滤前"))

  ggsave(
    filename = file.path(result_dir, paste0(current_sample, "_QC过滤前.pdf")),
    plot = p.qc.before,
    width = 16,
    height = 4
  )
  qc_before_plot_list[[current_sample]] <- p.qc.before

  cells_before <- ncol(scRNAlist[[i]])

  # subset中的条件必须同时满足。
  scRNAlist[[i]] <- subset(
    scRNAlist[[i]],
    subset = nFeature_RNA >= min_features &
      nFeature_RNA <= max_features &
      nCount_RNA >= min_counts &
      nCount_RNA <= max_counts &
      percent_mito <= max_percent_mito &
      percent_hb <= max_percent_hb &
      percent_ribo >= min_percent_ribo
  )

  cells_after <- ncol(scRNAlist[[i]])

  p.qc.after <- VlnPlot(
    scRNAlist[[i]],
    features = c(
      "nFeature_RNA", "nCount_RNA",
      "percent_mito", "percent_ribo", "percent_hb"
    ),
    assay = "RNA",
    layer = "counts",
    pt.size = 0,
    ncol = 5
  ) + plot_annotation(title = paste0(current_sample, "_QC过滤后"))

  ggsave(
    filename = file.path(result_dir, paste0(current_sample, "_QC过滤后.pdf")),
    plot = p.qc.after,
    width = 16,
    height = 4
  )
  qc_after_plot_list[[current_sample]] <- p.qc.after
  qc_compare_plot_list[[current_sample]] <- p.qc.before / p.qc.after

  qc_summary <- rbind(
    qc_summary,
    data.frame(
      sample_id = current_sample,
      cells_before_QC = cells_before,
      cells_after_QC = cells_after,
      median_nFeature = median(scRNAlist[[i]]$nFeature_RNA),
      median_nCount = median(scRNAlist[[i]]$nCount_RNA),
      median_percent_mito = median(scRNAlist[[i]]$percent_mito)
    )
  )
}

# 每个样本的QC图已经分别保存。下面再把所有样本合并到总览PDF中，
# 便于一次性横向查看样本质量差异，而不用逐个打开单独文件。
if (length(qc_before_plot_list) > 0) {
  p.qc.before.all <- wrap_plots(qc_before_plot_list, ncol = 1)
  ggsave(
    filename = file.path(result_dir, "所有样本_QC过滤前_合并展示.pdf"),
    plot = p.qc.before.all,
    width = 16,
    height = max(4, 4 * length(qc_before_plot_list)),
    limitsize = FALSE
  )

  p.qc.after.all <- wrap_plots(qc_after_plot_list, ncol = 1)
  ggsave(
    filename = file.path(result_dir, "所有样本_QC过滤后_合并展示.pdf"),
    plot = p.qc.after.all,
    width = 16,
    height = max(4, 4 * length(qc_after_plot_list)),
    limitsize = FALSE
  )

  p.qc.compare.all <- wrap_plots(qc_compare_plot_list, ncol = 1)
  ggsave(
    filename = file.path(result_dir, "所有样本_QC过滤前后对比_合并展示.pdf"),
    plot = p.qc.compare.all,
    width = 16,
    height = max(8, 8 * length(qc_compare_plot_list)),
    limitsize = FALSE
  )
}

mito_summary
write.csv(
  mito_summary,
  file.path(result_dir, "各样本线粒体基因与比例判断_QC前.csv"),
  row.names = FALSE
)

qc_summary
write.csv(qc_summary, file.path(result_dir, "QC汇总表.csv"), row.names = FALSE)
saveRDS(scRNAlist, file.path(result_dir, "1.逐样本QC后对象列表.rds"))



