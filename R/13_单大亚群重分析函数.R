# ==============================================================================
# 入口4使用的单大亚群反复试算函数
# ==============================================================================


single_recluster_number_stub <- function(x) {
  answer <- format(as.numeric(x), scientific = FALSE, trim = TRUE)
  gsub("[.]", "p", answer)
}


read_single_recluster_pc <- function(default_pc, maximum_pc) {
  repeat {
    answer <- trimws(readline(paste0(
      "请输入本轮使用的PC数量（1-", maximum_pc,
      "；直接回车使用", default_pc, "）："
    )))
    selected <- if (!nzchar(answer)) {
      as.integer(default_pc)
    } else {
      suppressWarnings(as.integer(answer))
    }
    if (length(selected) == 1 &&
        is.finite(selected) &&
        selected >= 1 &&
        selected <= maximum_pc) {
      return(selected)
    }
    message("PC必须是1-", maximum_pc, "之间的整数，请重新输入。")
  }
}


read_single_recluster_resolution <- function(default_resolution) {
  repeat {
    answer <- trimws(readline(paste0(
      "请输入本轮resolution（任意大于0的数；直接回车使用",
      default_resolution, "）："
    )))
    selected <- if (!nzchar(answer)) {
      as.numeric(default_resolution)
    } else {
      suppressWarnings(as.numeric(answer))
    }
    if (length(selected) == 1 &&
        is.finite(selected) &&
        selected > 0) {
      return(selected)
    }
    message("resolution必须是大于0的有限数值，请重新输入。")
  }
}


ask_use_single_recluster_trial <- function() {
  repeat {
    answer <- tolower(trimws(readline(
      "是否把本轮PC/res作为最终聚类方案？输入y采用，n继续试算，q退出："
    )))
    if (answer %in% c("y", "yes")) return(TRUE)
    if (answer %in% c("n", "no", "")) return(FALSE)
    if (answer %in% c("q", "quit", "stop")) {
      stop("用户结束单大亚群重分析；已经生成的试算图和参数均已保留。")
    }
    message("请输入y、n或q。")
  }
}


run_one_single_recluster_trial <- function(
    base_object,
    reduction,
    celltype,
    pc,
    resolution,
    algorithm,
    trial_number,
    group_output_dir) {
  trial_stub <- paste0(
    "试算_", sprintf("%03d", trial_number),
    "_PC", pc,
    "_res", single_recluster_number_stub(resolution)
  )
  trial_dir <- file.path(group_output_dir, "PC_res试算记录", trial_stub)
  dir.create(trial_dir, recursive = TRUE, showWarnings = FALSE)

  selected <- data.frame(
    PC = as.integer(pc),
    resolution = as.numeric(resolution),
    algorithm = as.integer(algorithm),
    algorithm_name = algorithm_label(algorithm),
    stringsAsFactors = FALSE
  )
  trial_object <- run_final_subcluster(
    object = base_object,
    reduction = reduction,
    selected = selected
  )

  cluster_summary <- as.data.frame(
    table(trial_object$subcluster),
    stringsAsFactors = FALSE
  )
  colnames(cluster_summary) <- c("subcluster", "cell_number")
  cluster_summary$cell_fraction <- cluster_summary$cell_number /
    ncol(trial_object)
  write.csv(
    cluster_summary,
    file.path(trial_dir, "试算cluster细胞数量.csv"),
    row.names = FALSE
  )
  write.csv(
    data.frame(
      celltype = celltype,
      trial_number = trial_number,
      reduction = reduction,
      PC = pc,
      resolution = resolution,
      algorithm = algorithm,
      algorithm_name = algorithm_label(algorithm),
      cluster_number = nrow(cluster_summary),
      smallest_cluster = min(cluster_summary$cell_number),
      stringsAsFactors = FALSE
    ),
    file.path(trial_dir, "本轮试算参数.csv"),
    row.names = FALSE
  )

  trial_plot <- DimPlot(
    trial_object,
    reduction = "umap",
    group.by = "subcluster",
    label = TRUE,
    repel = TRUE,
    shuffle = TRUE,
    raster = FALSE
  ) +
    ggtitle(paste0(
      celltype,
      " | PC=", pc,
      " | resolution=", resolution,
      " | cluster=", nrow(cluster_summary)
    )) +
    NoLegend()
  umap_file <- file.path(trial_dir, "本轮试算_UMAP.pdf")
  ggsave(
    filename = umap_file,
    plot = trial_plot,
    width = 10,
    height = 8
  )
  if (isTRUE(single_recluster_save_trial_object)) {
    saveRDS(
      trial_object,
      file.path(trial_dir, "本轮试算_sce.rds")
    )
  }

  list(
    object = trial_object,
    selected = selected,
    cluster_summary = cluster_summary,
    trial_dir = trial_dir,
    umap_file = umap_file
  )
}


save_accepted_single_recluster_results <- function(
    object,
    selected,
    reduction,
    celltype,
    output_dir) {
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
  write.csv(
    mapping,
    file.path(output_dir, "细胞与小亚群对应关系.csv"),
    row.names = FALSE
  )

  cluster_summary <- as.data.frame(
    table(object$subcluster),
    stringsAsFactors = FALSE
  )
  colnames(cluster_summary) <- c("subcluster", "cell_number")
  cluster_summary$celltype <- celltype
  cluster_summary$cell_fraction <- cluster_summary$cell_number / ncol(object)
  write.csv(
    cluster_summary,
    file.path(output_dir, "小亚群细胞数量汇总.csv"),
    row.names = FALSE
  )

  final_parameters <- selected
  final_parameters$celltype <- celltype
  final_parameters$reduction <- reduction
  final_parameters$selection_mode <- "interactive_repeat_until_confirmed"
  write.csv(
    final_parameters,
    file.path(output_dir, "最终确认的PC_res参数.csv"),
    row.names = FALSE
  )

  saveRDS(
    object,
    file.path(output_dir, "最终方案_待人工注释_sce.rds")
  )
  object
}
