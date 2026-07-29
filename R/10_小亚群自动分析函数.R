# ==============================================================================
# 小亚群自动重聚类函数
# 输入：已经完成大亚群注释的Seurat对象
# 输出：每个大亚群独立的Seurat对象、候选评分、正式cluster、marker和图片
# 本文件只定义函数，由“运行_02_入口_从大亚群注释至小亚群完成.R”统一调度。
# ==============================================================================


# 将大亚群名称转换成Windows/Linux均可使用的文件夹名。
safe_subcluster_name <- function(x) {
  y <- trimws(as.character(x))
  y <- gsub("[\\\\/:*?\"<>|]", "_", y)
  y <- gsub("[[:cntrl:]]", "_", y)
  y <- gsub("[. ]+$", "", y)
  if (!nzchar(y)) y <- "未命名大亚群"
  y
}


algorithm_label <- function(algorithm) {
  labels <- c(
    `1` = "Louvain",
    `2` = "Louvain_multilevel_refinement",
    `3` = "SLM",
    `4` = "Leiden"
  )
  answer <- unname(labels[as.character(algorithm)])
  if (length(answer) == 0 || is.na(answer)) {
    answer <- paste0("algorithm_", algorithm)
  }
  answer
}


# 不依赖mclust，直接计算Adjusted Rand Index，用于评估改变随机种子后的聚类稳定性。
adjusted_rand_index <- function(labels1, labels2) {
  if (length(labels1) != length(labels2) || length(labels1) < 2) return(NA_real_)
  tab <- table(labels1, labels2)
  choose2 <- function(x) x * (x - 1) / 2
  sum_cells <- sum(choose2(tab))
  sum_rows <- sum(choose2(rowSums(tab)))
  sum_cols <- sum(choose2(colSums(tab)))
  total_pairs <- choose2(sum(tab))
  if (total_pairs == 0) return(NA_real_)
  expected <- sum_rows * sum_cols / total_pairs
  maximum <- (sum_rows + sum_cols) / 2
  denominator <- maximum - expected
  if (abs(denominator) < .Machine$double.eps) {
    return(if (identical(as.character(labels1), as.character(labels2))) 1 else 0)
  }
  (sum_cells - expected) / denominator
}


mean_silhouette_score <- function(embedding, labels, max_cells = 2000, seed = 1234) {
  labels <- as.character(labels)
  names(labels) <- rownames(embedding)
  common_cells <- intersect(rownames(embedding), names(labels))
  embedding <- embedding[common_cells, , drop = FALSE]
  labels <- labels[common_cells]

  if (length(unique(labels)) < 2 || nrow(embedding) < 3) return(NA_real_)

  if (nrow(embedding) > max_cells) {
    set.seed(seed)
    # 按cluster近似分层抽样，避免稀有cluster在随机抽样中完全消失。
    split_cells <- split(seq_along(labels), labels)
    sampled <- unique(unlist(lapply(split_cells, function(index) {
      quota <- max(2L, ceiling(max_cells * length(index) / length(labels)))
      sample(index, min(length(index), quota))
    }), use.names = FALSE))
    if (length(sampled) > max_cells) sampled <- sample(sampled, max_cells)
    embedding <- embedding[sampled, , drop = FALSE]
    labels <- labels[sampled]
  }

  numeric_labels <- as.integer(factor(labels))
  if (length(unique(numeric_labels)) < 2) return(NA_real_)
  sil <- cluster::silhouette(numeric_labels, stats::dist(embedding))
  mean(sil[, "sil_width"], na.rm = TRUE)
}


choose_subcluster_pc <- function(object, min_pc = 5L) {
  stdev <- object[["pca"]]@stdev
  if (length(stdev) == 0) stop("PCA没有可用维度，无法自动选择PC。")

  # PCA解释率应使用方差（标准差平方），而不是直接使用标准差。
  pct <- stdev^2 / sum(stdev^2) * 100
  cumulative <- cumsum(pct)

  co1_candidates <- which(cumulative >= 90 & pct < 5)
  co1 <- if (length(co1_candidates) > 0) co1_candidates[1] else length(pct)

  drops <- pct[-length(pct)] - pct[-1]
  co2_candidates <- which(drops > 0.1) + 1L
  co2 <- if (length(co2_candidates) > 0) max(co2_candidates) else co1

  recommended <- min(co1, co2, length(pct))
  recommended <- max(min_pc, recommended)
  recommended <- min(recommended, length(pct))

  list(
    recommended = as.integer(recommended),
    co1 = as.integer(co1),
    co2 = as.integer(co2),
    table = data.frame(
      PC = seq_along(pct),
      variance_percent = pct,
      cumulative_percent = cumulative
    )
  )
}


prepare_subcluster_object <- function(parent_object, celltype) {
  selected_cells <- rownames(parent_object@meta.data)[
    as.character(parent_object$celltype) == celltype
  ]
  if (length(selected_cells) == 0) stop("大亚群中没有找到细胞：", celltype)

  selected <- subset(parent_object, cells = selected_cells)
  DefaultAssay(selected) <- "RNA"

  # Seurat V5的多counts层必须先合并，才能安全提取完整原始counts重建对象。
  if (inherits(selected[["RNA"]], "Assay5")) {
    count_layers <- Layers(selected[["RNA"]], search = "^counts")
    if (length(count_layers) == 0) {
      stop("RNA assay中没有counts层，无法从原始计数重建小亚群对象。")
    }
    if (length(count_layers) > 1 || !"counts" %in% count_layers) {
      selected[["RNA"]] <- JoinLayers(selected[["RNA"]])
    }
  }

  counts <- GetAssayData(selected, assay = "RNA", layer = "counts")
  if (nrow(counts) == 0 || ncol(counts) == 0) {
    stop("无法从RNA assay提取counts，请确认大亚群对象保留了RNA counts层。")
  }

  metadata <- selected@meta.data[colnames(counts), , drop = FALSE]
  object <- CreateSeuratObject(
    counts = counts,
    meta.data = metadata,
    project = safe_subcluster_name(celltype),
    min.cells = 0,
    min.features = 0
  )
  object$parent_celltype <- celltype
  object
}


normalize_and_reduce_subcluster <- function(object, celltype, output_dir) {
  object <- NormalizeData(object, assay = "RNA", verbose = FALSE)
  object <- FindVariableFeatures(
    object,
    assay = "RNA",
    selection.method = "vst",
    nfeatures = min(subcluster_n_variable_features, nrow(object)),
    verbose = FALSE
  )

  variable_genes <- VariableFeatures(object)
  excluded <- rep(FALSE, length(variable_genes))
  for (pattern in subcluster_exclude_variable_patterns) {
    excluded <- excluded | grepl(pattern, variable_genes)
  }
  VariableFeatures(object) <- variable_genes[!excluded]
  if (length(VariableFeatures(object)) < 100) {
    stop("排除线粒体/核糖体等基因后高变基因不足100个，请检查基因命名或参数。")
  }

  regress_variables <- intersect(subcluster_vars_to_regress, colnames(object@meta.data))
  missing_regress <- setdiff(subcluster_vars_to_regress, regress_variables)
  if (length(missing_regress) > 0) {
    warning("以下回归变量不在metadata中，已跳过：", paste(missing_regress, collapse = "、"))
  }

  object <- ScaleData(
    object,
    assay = "RNA",
    features = VariableFeatures(object),
    vars.to.regress = if (length(regress_variables) > 0) regress_variables else NULL,
    verbose = FALSE
  )

  maximum_pcs <- min(
    as.integer(subcluster_n_pcs_calculate),
    ncol(object) - 1L,
    length(VariableFeatures(object)) - 1L
  )
  if (maximum_pcs < 2) stop("细胞数或高变基因数过少，无法进行PCA。")

  object <- RunPCA(
    object,
    assay = "RNA",
    features = VariableFeatures(object),
    npcs = maximum_pcs,
    seed.use = 1234,
    verbose = FALSE
  )

  pc_choice <- choose_subcluster_pc(object, min_pc = subcluster_min_pc)
  write.csv(
    pc_choice$table,
    file.path(output_dir, "PCA解释率与累计贡献率.csv"),
    row.names = FALSE
  )
  write.csv(
    data.frame(
      co1 = pc_choice$co1,
      co2 = pc_choice$co2,
      recommended_PC = pc_choice$recommended,
      actual_PCA_dimensions = maximum_pcs
    ),
    file.path(output_dir, "PCA自动推荐结果.csv"),
    row.names = FALSE
  )

  p.elbow <- ElbowPlot(object, ndims = maximum_pcs) +
    geom_vline(xintercept = pc_choice$recommended, linetype = 2, color = "red") +
    ggtitle(paste0(celltype, " | 自动推荐PC = ", pc_choice$recommended))
  ggsave(file.path(output_dir, "PCA_ElbowPlot.pdf"), p.elbow, width = 8, height = 6)

  list(object = object, pc_choice = pc_choice, maximum_pcs = maximum_pcs)
}


run_subcluster_harmony <- function(object, celltype, maximum_pcs) {
  mode <- subcluster_harmony_mode
  override <- subcluster_harmony_overrides[[celltype]]
  if (!is.null(override)) mode <- override
  mode <- tolower(as.character(mode)[1])
  if (!mode %in% c("always", "auto", "none")) {
    stop("Harmony模式只能是'always'、'auto'或'none'。当前大亚群：", celltype)
  }

  should_run <- mode != "none"
  reason <- "按照配置运行Harmony"
  if (!subcluster_harmony_group %in% colnames(object@meta.data)) {
    if (mode == "always") {
      stop("Harmony批次列不存在：", subcluster_harmony_group)
    }
    should_run <- FALSE
    reason <- "批次列不存在，auto模式跳过Harmony"
  } else {
    batch_values <- as.character(object@meta.data[[subcluster_harmony_group]])
    invalid_batch <- is.na(batch_values) | !nzchar(batch_values)
    if (any(invalid_batch) && mode != "none") {
      stop(
        "Harmony批次列“", subcluster_harmony_group, "”中有",
        sum(invalid_batch), "个细胞缺少批次值。请先补全metadata。"
      )
    }
    valid_batches <- unique(batch_values[!is.na(batch_values) & nzchar(batch_values)])
    if (length(valid_batches) < 2) {
      should_run <- FALSE
      reason <- "只有一个有效批次水平，Harmony没有可校正的批次"
    }
  }

  if (should_run) {
    object <- RunHarmony(
      object = object,
      group.by.vars = subcluster_harmony_group,
      reduction.use = "pca",
      dims.use = seq_len(maximum_pcs),
      reduction.save = "harmony",
      verbose = TRUE
    )
    reduction <- "harmony"
  } else {
    reduction <- "pca"
  }

  list(object = object, reduction = reduction, harmony_mode = mode, reason = reason)
}


candidate_size_metrics <- function(labels) {
  cluster_sizes <- table(labels)
  threshold <- max(
    as.integer(subcluster_min_cluster_cells),
    ceiling(length(labels) * subcluster_min_cluster_fraction)
  )
  small_clusters <- names(cluster_sizes)[cluster_sizes < threshold]
  tiny_cell_fraction <- if (length(small_clusters) == 0) {
    0
  } else {
    sum(labels %in% small_clusters) / length(labels)
  }
  list(
    cluster_number = length(cluster_sizes),
    smallest_cluster = min(cluster_sizes),
    tiny_cluster_number = length(small_clusters),
    tiny_cell_fraction = tiny_cell_fraction,
    cluster_size_score = 1 - tiny_cell_fraction
  )
}


evaluate_one_clustering <- function(
    neighbor_object,
    embedding,
    pc,
    resolution,
    algorithm,
    candidate_id) {

  tryCatch({
    first <- FindClusters(
      neighbor_object,
      resolution = resolution,
      algorithm = algorithm,
      random.seed = 1234,
      verbose = FALSE
    )
    labels1 <- as.character(first$seurat_clusters)
    names(labels1) <- colnames(first)

    second <- FindClusters(
      neighbor_object,
      resolution = resolution,
      algorithm = algorithm,
      random.seed = 4321,
      verbose = FALSE
    )
    labels2 <- as.character(second$seurat_clusters)
    names(labels2) <- colnames(second)

    stability <- adjusted_rand_index(labels1, labels2)
    silhouette <- mean_silhouette_score(
      embedding = embedding,
      labels = labels1,
      max_cells = subcluster_evaluation_max_cells,
      seed = 1234
    )
    size_metrics <- candidate_size_metrics(labels1)

    # 将轮廓系数从[-1,1]线性转换为[0,1]；稳定性负值按0处理。
    silhouette_scaled <- if (is.na(silhouette)) 0 else (silhouette + 1) / 2
    stability_scaled <- if (is.na(stability)) 0 else max(0, min(1, stability))
    required_weights <- c("silhouette", "stability", "cluster_size")
    if (!all(required_weights %in% names(subcluster_score_weights))) {
      stop(
        "subcluster_score_weights必须包含以下名称：",
        paste(required_weights, collapse = "、")
      )
    }
    if (any(!is.finite(subcluster_score_weights[required_weights])) ||
        sum(subcluster_score_weights[required_weights]) <= 0) {
      stop("subcluster_score_weights必须是有限数值，且权重总和大于0。")
    }
    weights <- subcluster_score_weights[required_weights] /
      sum(subcluster_score_weights[required_weights])
    score <- if (size_metrics$cluster_number < 2) {
      NA_real_
    } else {
      weights[["silhouette"]] * silhouette_scaled +
        weights[["stability"]] * stability_scaled +
        weights[["cluster_size"]] * size_metrics$cluster_size_score
    }

    data.frame(
      candidate_id = candidate_id,
      PC = pc,
      resolution = resolution,
      algorithm = algorithm,
      algorithm_name = algorithm_label(algorithm),
      cluster_number = size_metrics$cluster_number,
      smallest_cluster = size_metrics$smallest_cluster,
      tiny_cluster_number = size_metrics$tiny_cluster_number,
      tiny_cell_fraction = size_metrics$tiny_cell_fraction,
      mean_silhouette = silhouette,
      stability_ARI = stability,
      recommendation_score = score,
      calculation_status = "success",
      error_message = NA_character_,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      candidate_id = candidate_id,
      PC = pc,
      resolution = resolution,
      algorithm = algorithm,
      algorithm_name = algorithm_label(algorithm),
      cluster_number = NA_integer_,
      smallest_cluster = NA_integer_,
      tiny_cluster_number = NA_integer_,
      tiny_cell_fraction = NA_real_,
      mean_silhouette = NA_real_,
      stability_ARI = NA_real_,
      recommendation_score = NA_real_,
      calculation_status = "failed",
      error_message = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })
}


evaluate_clustering_grid <- function(object, reduction, pcs, resolutions, algorithms) {
  results <- list()
  result_index <- 0L
  candidate_id <- 0L

  for (pc in pcs) {
    message("  构建邻居图：PC 1:", pc)
    neighbor_object <- FindNeighbors(
      object,
      reduction = reduction,
      dims = seq_len(pc),
      verbose = FALSE
    )
    embedding <- Embeddings(object, reduction = reduction)[, seq_len(pc), drop = FALSE]

    for (algorithm in algorithms) {
      for (resolution in resolutions) {
        candidate_id <- candidate_id + 1L
        message(
          "    候选", candidate_id, "：", algorithm_label(algorithm),
          " | resolution=", resolution
        )
        result_index <- result_index + 1L
        results[[result_index]] <- evaluate_one_clustering(
          neighbor_object = neighbor_object,
          embedding = embedding,
          pc = pc,
          resolution = resolution,
          algorithm = algorithm,
          candidate_id = candidate_id
        )
      }
    }
    rm(neighbor_object)
    invisible(gc())
  }

  do.call(rbind, results)
}


save_subcluster_resolution_diagnostics <- function(
    object,
    reduction,
    pc,
    resolutions,
    algorithm,
    celltype,
    output_dir) {
  diagnostic_object <- FindNeighbors(
    object,
    reduction = reduction,
    dims = seq_len(pc),
    verbose = FALSE
  )
  diagnostic_object <- RunUMAP(
    diagnostic_object,
    reduction = reduction,
    dims = seq_len(pc),
    n.neighbors = subcluster_umap_n_neighbors,
    min.dist = subcluster_umap_min_dist,
    seed.use = 1234,
    verbose = FALSE
  )

  resolution_columns <- character(length(resolutions))
  resolution_plots <- vector("list", length(resolutions))
  for (resolution_index in seq_along(resolutions)) {
    current_resolution <- resolutions[resolution_index]
    diagnostic_object <- FindClusters(
      diagnostic_object,
      resolution = current_resolution,
      algorithm = algorithm,
      random.seed = 1234,
      verbose = FALSE
    )
    resolution_column <- paste0(
      "subcluster_res.",
      format(current_resolution, scientific = FALSE, trim = TRUE)
    )
    diagnostic_object[[resolution_column]] <- as.character(
      diagnostic_object$seurat_clusters
    )
    resolution_columns[resolution_index] <- resolution_column
    resolution_plots[[resolution_index]] <- DimPlot(
      diagnostic_object,
      reduction = "umap",
      group.by = resolution_column,
      label = TRUE,
      repel = TRUE,
      shuffle = TRUE,
      raster = FALSE
    ) +
      ggtitle(paste0(
        celltype,
        " | PC 1:", pc,
        " | resolution = ", current_resolution
      )) +
      NoLegend()
  }

  p.resolution <- wrap_plots(resolution_plots, ncol = 3)
  ggsave(
    filename = file.path(output_dir, "不同resolution_UMAP.pdf"),
    plot = p.resolution,
    width = 15,
    height = 5 * ceiling(length(resolution_plots) / 3)
  )

  if (requireNamespace("clustree", quietly = TRUE)) {
    tryCatch({
      p.clustree <- clustree::clustree(
        diagnostic_object@meta.data,
        prefix = "subcluster_res."
      )
      ggsave(
        filename = file.path(output_dir, "不同resolution_clustree.pdf"),
        plot = p.clustree,
        width = 12,
        height = 9
      )
    }, error = function(e) {
      warning(
        "大亚群“", celltype, "”的clustree绘图失败，但不影响resolution选择。\n",
        "具体原因：", conditionMessage(e), "\n",
        "请查看“不同resolution_UMAP.pdf”。"
      )
    })
  } else {
    message(
      "没有安装clustree，已跳过大亚群“", celltype,
      "”的clustree；resolution UMAP仍已生成。"
    )
  }

  invisible(resolution_columns)
}


select_subcluster_resolution <- function(candidate_table, celltype) {
  successful <- candidate_table[
    candidate_table$calculation_status == "success" &
      is.finite(candidate_table$recommendation_score),
    ,
    drop = FALSE
  ]
  if (nrow(successful) == 0) {
    failed_messages <- unique(na.omit(candidate_table$error_message))
    stop(
      "所有聚类候选均计算失败。首个错误：",
      if (length(failed_messages) > 0) failed_messages[1] else "未知错误"
    )
  }

  ordered <- successful[order(
    -successful$recommendation_score,
    -successful$mean_silhouette,
    successful$tiny_cell_fraction,
    successful$cluster_number,
    successful$candidate_id
  ), , drop = FALSE]

  candidate_resolutions <- sort(unique(successful$resolution))
  default_resolution <- ordered$resolution[1]
  configured_resolution <- subcluster_resolution_overrides[[celltype]]

  if (!is.null(configured_resolution) && length(configured_resolution) > 0) {
    configured_resolution <- suppressWarnings(as.numeric(configured_resolution[1]))
    resolution_match <- which(
      abs(candidate_resolutions - configured_resolution) < 1e-10
    )
    if (is.na(configured_resolution) || length(resolution_match) != 1) {
      stop(
        "大亚群“", celltype, "”预设的resolution无效。必须是以下候选值之一：",
        paste(candidate_resolutions, collapse = "、")
      )
    }
    selected_resolution <- candidate_resolutions[resolution_match]
    selected <- successful[
      abs(successful$resolution - selected_resolution) < 1e-10,
      ,
      drop = FALSE
    ][1, , drop = FALSE]
    selected$selection_source <- "config_override"
    message(
      "大亚群“", celltype, "”使用预设resolution：", selected_resolution
    )
    return(selected)
  }

  if (!interactive()) {
    stop(
      "大亚群“", celltype, "”的resolution比较图已经生成，但当前不是交互式R会话。",
      "请在config/04_小亚群分析参数.R的subcluster_resolution_overrides中",
      "为该大亚群填写以下候选值之一：",
      paste(candidate_resolutions, collapse = "、"),
      "，然后重新运行入口2。"
    )
  }

  display_columns <- intersect(
    c(
      "resolution", "cluster_number", "smallest_cluster",
      "mean_silhouette", "stability_ARI", "recommendation_score"
    ),
    colnames(ordered)
  )
  print(ordered[, display_columns, drop = FALSE], row.names = FALSE)

  repeat {
    answer <- trimws(readline(paste0(
      "请查看大亚群“", celltype,
      "”目录中的不同resolution_UMAP.pdf和clustree，",
      "输入正式resolution（候选值：",
      paste(candidate_resolutions, collapse = "、"),
      "；直接回车使用评分推荐值", default_resolution, "）："
    )))

    if (!nzchar(answer)) {
      selected_resolution <- default_resolution
    } else {
      selected_resolution <- suppressWarnings(as.numeric(answer))
    }
    resolution_match <- if (length(selected_resolution) == 1 &&
                            is.finite(selected_resolution)) {
      which(abs(candidate_resolutions - selected_resolution) < 1e-10)
    } else {
      integer()
    }

    if (length(resolution_match) == 1) {
      selected_resolution <- candidate_resolutions[resolution_match]
      selected <- successful[
        abs(successful$resolution - selected_resolution) < 1e-10,
        ,
        drop = FALSE
      ][1, , drop = FALSE]
      selected$selection_source <- if (nzchar(answer)) {
        "interactive"
      } else {
        "interactive_default"
      }
      return(selected)
    }

    message(
      "输入无效：", if (nzchar(answer)) answer else "<空值>",
      "。请输入候选值之一：", paste(candidate_resolutions, collapse = "、"),
      "。本次不会写入错误缓存，请重新输入。"
    )
  }
}


run_final_subcluster <- function(object, reduction, selected) {
  pc <- as.integer(selected$PC[1])
  resolution <- as.numeric(selected$resolution[1])
  algorithm <- as.integer(selected$algorithm[1])

  object <- FindNeighbors(
    object,
    reduction = reduction,
    dims = seq_len(pc),
    verbose = FALSE
  )
  object <- FindClusters(
    object,
    resolution = resolution,
    algorithm = algorithm,
    random.seed = 1234,
    verbose = FALSE
  )
  object <- RunUMAP(
    object,
    reduction = reduction,
    dims = seq_len(pc),
    n.neighbors = subcluster_umap_n_neighbors,
    min.dist = subcluster_umap_min_dist,
    seed.use = 1234,
    verbose = FALSE
  )
  object$subcluster <- as.character(object$seurat_clusters)
  Idents(object) <- "subcluster"
  object
}


save_subcluster_plots <- function(object, celltype, reduction, output_dir) {
  p.cluster <- DimPlot(
    object,
    reduction = "umap",
    group.by = "subcluster",
    label = TRUE,
    repel = TRUE,
    shuffle = TRUE,
    raster = FALSE
  ) + ggtitle(paste0(celltype, " | 小亚群cluster")) + NoLegend()

  group_column <- if (subcluster_harmony_group %in% colnames(object@meta.data)) {
    subcluster_harmony_group
  } else {
    "orig.ident"
  }
  p.batch <- DimPlot(
    object,
    reduction = "umap",
    group.by = group_column,
    shuffle = TRUE,
    raster = FALSE
  ) + ggtitle(paste0("按", group_column, "着色"))

  ggsave(
    file.path(output_dir, "正式小亚群聚类_UMAP.pdf"),
    p.cluster + p.batch,
    width = 15,
    height = 7
  )

  if ("sample_id" %in% colnames(object@meta.data)) {
    sample_number <- length(unique(object$sample_id))
    p.split <- DimPlot(
      object,
      reduction = "umap",
      group.by = "subcluster",
      split.by = "sample_id",
      label = FALSE,
      raster = FALSE
    )
    ggsave(
      file.path(output_dir, "各样本小亚群聚类_UMAP.pdf"),
      p.split,
      width = max(12, 4 * sample_number),
      height = 6
    )
  }

  invisible(NULL)
}


find_and_save_subcluster_markers <- function(object, output_dir) {
  if (inherits(object[["RNA"]], "Assay5")) object[["RNA"]] <- JoinLayers(object[["RNA"]])
  DefaultAssay(object) <- "RNA"
  Idents(object) <- "subcluster"

  markers <- FindAllMarkers(
    object = object,
    assay = "RNA",
    test.use = subcluster_marker_test,
    only.pos = subcluster_marker_only_pos,
    min.pct = subcluster_marker_min_pct,
    logfc.threshold = subcluster_marker_logfc_threshold
  )
  write.csv(markers, file.path(output_dir, "所有小亚群_marker.csv"), row.names = FALSE)

  if (nrow(markers) == 0) {
    warning("FindAllMarkers没有返回marker，已跳过marker筛选表、DotPlot和热图。")
    return(list(object = object, markers = markers, filtered = markers, top_markers = markers))
  }

  filtered <- markers %>%
    filter(!is.na(p_val_adj), p_val_adj < subcluster_marker_padj_threshold)
  write.csv(
    filtered,
    file.path(output_dir, "显著小亚群_marker.csv"),
    row.names = FALSE
  )

  positive <- filtered %>% filter(avg_log2FC > 0)
  top_markers <- positive %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = subcluster_top_marker_number, with_ties = FALSE) %>%
    ungroup()
  write.csv(
    top_markers,
    file.path(output_dir, paste0("每个小亚群_top", subcluster_top_marker_number, "_marker.csv")),
    row.names = FALSE
  )

  if (nrow(top_markers) > 0) {
    dot_genes <- unique(top_markers$gene)
    p.dot <- DotPlot(object, features = dot_genes, assay = "RNA", group.by = "subcluster") +
      coord_flip() + theme_bw() + ggtitle("小亚群top marker")
    ggsave(
      file.path(output_dir, "小亚群_top_marker_DotPlot.pdf"),
      p.dot,
      width = 12,
      height = max(8, length(dot_genes) * 0.18)
    )

    heatmap_markers <- positive %>%
      group_by(cluster) %>%
      slice_max(order_by = avg_log2FC, n = subcluster_heatmap_top_number, with_ties = FALSE) %>%
      ungroup()
    heatmap_genes <- unique(heatmap_markers$gene)
    if (length(heatmap_genes) > 0) {
      object <- ScaleData(object, assay = "RNA", features = heatmap_genes, verbose = FALSE)
      p.heatmap <- DoHeatmap(
        object,
        features = heatmap_genes,
        group.by = "subcluster",
        raster = TRUE
      ) + NoLegend()
      ggsave(
        file.path(output_dir, "小亚群_top_marker热图.pdf"),
        p.heatmap,
        width = 12,
        height = max(8, length(heatmap_genes) * 0.12)
      )
    }
  }

  list(object = object, markers = markers, filtered = filtered, top_markers = top_markers)
}


analyze_one_parent_celltype <- function(parent_object, celltype, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  harmony_checkpoint_file <- file.path(
    output_dir,
    "01_PCA_Harmony完成_中间检查点.rds"
  )
  final_clustering_checkpoint_file <- file.path(
    output_dir,
    "02_正式聚类完成_待marker_中间检查点.rds"
  )
  internal_resume <- exists("subcluster_internal_resume_allowed", inherits = TRUE) &&
    isTRUE(subcluster_internal_resume_allowed)

  resume_from_final_clustering <- internal_resume &&
    file.exists(final_clustering_checkpoint_file)
  resume_from_harmony <- internal_resume &&
    !resume_from_final_clustering &&
    file.exists(harmony_checkpoint_file)

  if (resume_from_final_clustering) {
    saved_stage <- readRDS(final_clustering_checkpoint_file)
    object <- saved_stage$object
    selected <- saved_stage$selected
    reduction <- saved_stage$reduction
    cell_number <- ncol(object)
    message("已恢复正式聚类检查点，将从UMAP/marker与结果保存阶段继续：", celltype)
  } else {
    if (resume_from_harmony) {
      saved_stage <- readRDS(harmony_checkpoint_file)
      object <- saved_stage$object
      pc_choice <- saved_stage$pc_choice
      maximum_pcs <- saved_stage$maximum_pcs
      reduction <- saved_stage$reduction
      cell_number <- ncol(object)
      message("已恢复PCA/Harmony检查点，将从聚类候选评估继续：", celltype)
    } else {
      message("提取并重建大亚群对象：", celltype)
      object <- prepare_subcluster_object(parent_object, celltype)
      cell_number <- ncol(object)

      if (cell_number < subcluster_min_cells) {
        return(list(
          status = "skipped_too_few_cells",
          cell_number = cell_number,
          celltype = celltype
        ))
      }

      normalized <- normalize_and_reduce_subcluster(object, celltype, output_dir)
      object <- normalized$object
      pc_choice <- normalized$pc_choice
      maximum_pcs <- normalized$maximum_pcs

      harmony_result <- run_subcluster_harmony(object, celltype, maximum_pcs)
      object <- harmony_result$object
      reduction <- harmony_result$reduction
      writeLines(
        c(
          paste0("configured_mode=", harmony_result$harmony_mode),
          paste0("reduction_used=", reduction),
          paste0("reason=", harmony_result$reason),
          paste0("group_column=", subcluster_harmony_group)
        ),
        file.path(output_dir, "Harmony运行记录.txt")
      )

      saveRDS(
        list(
          object = object,
          pc_choice = pc_choice,
          maximum_pcs = maximum_pcs,
          reduction = reduction
        ),
        harmony_checkpoint_file
      )
    }

    # 每个大亚群固定使用自身co1/co2规则得出的推荐PC。
    # 用户不再选择PC或聚类算法，只逐个大亚群选择resolution。
    pcs <- as.integer(pc_choice$recommended)
    resolutions <- sort(unique(as.numeric(subcluster_resolution_candidates)))
    algorithms <- as.integer(subcluster_algorithm)[1]

    if (length(pcs) == 0) stop("没有可用的PC候选值。")
    if (length(resolutions) == 0 || any(!is.finite(resolutions) | resolutions <= 0)) {
      stop("resolution候选必须是大于0的有限数值。")
    }
    if (length(algorithms) == 0 || any(!algorithms %in% 1:4)) {
      stop("subcluster_algorithm只能填写1、2、3或4。")
    }

    message(
      "大亚群“", celltype, "”自动使用co1/co2推荐PC：1:", pcs,
      "；仅比较并询问resolution。"
    )
    message("开始评估resolution候选：", celltype)
    candidate_table <- evaluate_clustering_grid(
      object = object,
      reduction = reduction,
      pcs = pcs,
      resolutions = resolutions,
      algorithms = algorithms
    )
    candidate_table$celltype <- celltype
    candidate_table$reduction <- reduction
    candidate_table$co1 <- pc_choice$co1
    candidate_table$co2 <- pc_choice$co2
    candidate_table$recommended_PC_by_co1_co2 <- pc_choice$recommended
    candidate_table <- candidate_table[, c(
      "celltype", "candidate_id", "reduction", "PC", "resolution", "algorithm",
      "algorithm_name", "cluster_number", "smallest_cluster", "tiny_cluster_number",
      "tiny_cell_fraction", "mean_silhouette", "stability_ARI",
      "recommendation_score", "co1", "co2", "recommended_PC_by_co1_co2",
      "calculation_status", "error_message"
    )]

    # 先保存尚未选择的完整候选表；等待resolution输入时表格不会丢失。
    candidate_table$is_selected <- FALSE
    write.csv(
      candidate_table,
      file.path(output_dir, "聚类候选参数与自动评分.csv"),
      row.names = FALSE
    )

    successful_for_plot <- candidate_table[
      candidate_table$calculation_status == "success" &
        is.finite(candidate_table$recommendation_score),
      ,
      drop = FALSE
    ]
    if (nrow(successful_for_plot) > 0) {
      p.score <- ggplot(
        successful_for_plot,
        aes(x = resolution, y = cluster_number, color = recommendation_score)
      ) +
        geom_point(size = 3) +
        geom_line(aes(group = interaction(PC, algorithm)), alpha = 0.35) +
        facet_grid(algorithm_name ~ PC, labeller = label_both) +
        scale_color_viridis_c() +
        theme_bw() +
        labs(title = paste0(celltype, " | 聚类候选评分"), y = "cluster数量")
      ggsave(
        file.path(output_dir, "聚类候选评分图.pdf"),
        p.score,
        width = max(10, 3.2 * length(unique(successful_for_plot$PC))),
        height = max(6, 3.2 * length(unique(successful_for_plot$algorithm)))
      )
    }

    diagnostic_resolutions <- sort(unique(successful_for_plot$resolution))
    if (length(diagnostic_resolutions) > 0) save_subcluster_resolution_diagnostics(
      object = object,
      reduction = reduction,
      pc = pcs,
      resolutions = diagnostic_resolutions,
      algorithm = algorithms,
      celltype = celltype,
      output_dir = output_dir
    )

    selected <- select_subcluster_resolution(candidate_table, celltype)
    selected$is_selected <- TRUE
    candidate_table$is_selected <- candidate_table$candidate_id == selected$candidate_id[1]
    candidate_table <- candidate_table[order(
      candidate_table$calculation_status != "success",
      -candidate_table$recommendation_score
    ), , drop = FALSE]
    write.csv(
      candidate_table,
      file.path(output_dir, "聚类候选参数与自动评分.csv"),
      row.names = FALSE
    )

    selected$selection_mode <- "PC_auto_resolution_per_celltype"
    selected$co1 <- pc_choice$co1
    selected$co2 <- pc_choice$co2
    selected$recommended_PC_by_co1_co2 <- pc_choice$recommended
    selected$reduction <- reduction
    selected$harmony_group <- if (reduction == "harmony") {
      subcluster_harmony_group
    } else {
      NA_character_
    }
    write.csv(selected, file.path(output_dir, "正式小亚群聚类参数.csv"), row.names = FALSE)

    message(
      "正式聚类：PC 1:", selected$PC[1],
      " | resolution=", selected$resolution[1],
      " | ", selected$algorithm_name[1]
    )
    object <- run_final_subcluster(object, reduction, selected)
    saveRDS(
      list(object = object, selected = selected, reduction = reduction),
      final_clustering_checkpoint_file
    )
  }

  save_subcluster_plots(object, celltype, reduction, output_dir)

  marker_result <- find_and_save_subcluster_markers(object, output_dir)
  object <- marker_result$object

  mapping <- data.frame(
    cell = colnames(object),
    celltype = celltype,
    subcluster = as.character(object$subcluster),
    subcluster_label = paste0(celltype, "_", as.character(object$subcluster)),
    stringsAsFactors = FALSE
  )
  write.csv(mapping, file.path(output_dir, "细胞与小亚群对应关系.csv"), row.names = FALSE)

  cluster_summary <- as.data.frame(table(object$subcluster), stringsAsFactors = FALSE)
  colnames(cluster_summary) <- c("subcluster", "cell_number")
  cluster_summary$celltype <- celltype
  cluster_summary$cell_fraction <- cluster_summary$cell_number / ncol(object)
  write.csv(cluster_summary, file.path(output_dir, "小亚群细胞数量汇总.csv"), row.names = FALSE)

  saveRDS(object, file.path(output_dir, "小亚群聚类与marker完成_sce.rds"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

  list(
    status = "completed",
    cell_number = ncol(object),
    celltype = celltype,
    object = object,
    mapping = mapping,
    cluster_summary = cluster_summary,
    selected = selected
  )
}
