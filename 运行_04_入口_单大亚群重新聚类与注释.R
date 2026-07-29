# ==============================================================================
# 入口4：任选一个大亚群，反复试算PC/res，确认后完成聚类与人工注释
# ==============================================================================
# 交互流程：
#   选择大亚群 → PCA/Harmony → 输入任意PC/res → 查看UMAP
#   → 询问是否采用 → 不采用则继续输入下一组 → 采用后运行marker和先验DotPlot
#   → 生成单群注释XLSX → 写回注释并保存最终Seurat对象。
# ==============================================================================


# 自动取得入口脚本所在目录。
command_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(command_file) > 0) {
  entry_script <- sub("^--file=", "", command_file[1])
} else {
  source_files <- vapply(
    sys.frames(),
    function(x) {
      current_file <- x$ofile
      if (is.null(current_file)) NA_character_ else as.character(current_file)[1]
    },
    character(1)
  )
  source_files <- source_files[!is.na(source_files) & nzchar(source_files)]
  if (length(source_files) == 0) {
    stop(
      "无法取得入口脚本位置。请打开",
      "“运行_04_入口_单大亚群重新聚类与注释.R”并点击Source。"
    )
  }
  entry_script <- source_files[length(source_files)]
}
entry_script <- normalizePath(entry_script, mustWork = TRUE)
pipeline_root <- dirname(entry_script)
message("模块化流程目录：", pipeline_root)

required_pipeline_files <- file.path(
  pipeline_root,
  c(
    "config/01_流程参数.R",
    "config/04_小亚群分析参数.R",
    "config/05_小亚群marker点图参数.R",
    "config/06_小亚群人工注释参数.R",
    "config/07_外部大亚群RDS输入参数.R",
    "config/08_单大亚群重分析参数.R",
    "R/00_载入依赖包.R",
    "R/10_小亚群自动分析函数.R",
    "R/11_小亚群marker点图函数.R",
    "R/12_小亚群人工注释函数.R",
    "R/13_单大亚群重分析函数.R"
  )
)
missing_pipeline_files <- required_pipeline_files[!file.exists(required_pipeline_files)]
if (length(missing_pipeline_files) > 0) {
  stop(
    "入口4文件不完整。缺少：\n",
    paste(missing_pipeline_files, collapse = "\n")
  )
}

source(file.path(pipeline_root, "config", "01_流程参数.R"), encoding = "UTF-8")
source(file.path(pipeline_root, "R", "00_载入依赖包.R"), encoding = "UTF-8")
source(file.path(pipeline_root, "config", "04_小亚群分析参数.R"), encoding = "UTF-8")
source(file.path(pipeline_root, "config", "05_小亚群marker点图参数.R"), encoding = "UTF-8")
source(file.path(pipeline_root, "config", "06_小亚群人工注释参数.R"), encoding = "UTF-8")
source(file.path(pipeline_root, "config", "07_外部大亚群RDS输入参数.R"), encoding = "UTF-8")
source(file.path(pipeline_root, "config", "08_单大亚群重分析参数.R"), encoding = "UTF-8")

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop(
    "入口4需要R包openxlsx生成和读取人工注释表。请先运行：\n",
    "install.packages(\"openxlsx\")"
  )
}
source(
  file.path(pipeline_root, "R", "10_小亚群自动分析函数.R"),
  encoding = "UTF-8",
  local = .GlobalEnv
)
source(
  file.path(pipeline_root, "R", "11_小亚群marker点图函数.R"),
  encoding = "UTF-8",
  local = .GlobalEnv
)
source(
  file.path(pipeline_root, "R", "12_小亚群人工注释函数.R"),
  encoding = "UTF-8",
  local = .GlobalEnv
)
source(
  file.path(pipeline_root, "R", "13_单大亚群重分析函数.R"),
  encoding = "UTF-8",
  local = .GlobalEnv
)


# ----------------------------------------------------------------------------
# 读取大亚群注释对象
# ----------------------------------------------------------------------------

uses_external_annotated_rds <- isTRUE(downstream_use_external_annotated_rds)
if (uses_external_annotated_rds) {
  input_file <- trimws(as.character(downstream_external_annotated_rds)[1])
  celltype_column <- trimws(
    as.character(downstream_external_celltype_column)[1]
  )
  configured_rna_assay <- trimws(
    as.character(downstream_external_rna_assay)[1]
  )
} else {
  input_file <- if (!is.null(subcluster_input_file) &&
                    nzchar(as.character(subcluster_input_file)[1])) {
    as.character(subcluster_input_file)[1]
  } else {
    file.path(result_dir, "5.大亚群注释完成_sce.rds")
  }
  celltype_column <- "celltype"
  configured_rna_assay <- "RNA"
}

if (!nzchar(input_file) || !file.exists(input_file)) {
  stop(
    "没有找到大亚群注释完成Seurat对象：", input_file,
    "\n可在config/07中指定任意外部RDS文件。"
  )
}
message("读取大亚群注释对象：", input_file)
parent_object <- readRDS(input_file)
if (!inherits(parent_object, "Seurat")) {
  stop("指定RDS读取后不是Seurat对象。")
}
if (!celltype_column %in% colnames(parent_object@meta.data)) {
  stop("对象metadata中没有配置的大亚群注释列：", celltype_column)
}
parent_object$celltype <- as.character(
  parent_object@meta.data[[celltype_column]]
)
blank_celltype <- is.na(parent_object$celltype) |
  !nzchar(trimws(parent_object$celltype))
if (any(blank_celltype) && !isTRUE(downstream_allow_blank_celltype)) {
  stop(
    "大亚群注释列中有", sum(blank_celltype),
    "个细胞名称为空；请先补充或在config/07中允许空值。"
  )
}
parent_object$celltype[blank_celltype] <- "unknown"

if (!"RNA" %in% Assays(parent_object)) {
  if (!configured_rna_assay %in% Assays(parent_object)) {
    stop(
      "对象中没有RNA assay，也没有配置的assay：",
      configured_rna_assay
    )
  }
  message("把assay“", configured_rna_assay, "”复制为RNA供下游使用。")
  parent_object[["RNA"]] <- parent_object[[configured_rna_assay]]
}


# ----------------------------------------------------------------------------
# 选择一个真实大亚群；支持完整名称、编号或兼容关键词
# ----------------------------------------------------------------------------

available_celltypes <- unique(as.character(parent_object$celltype))
available_celltypes <- available_celltypes[
  !is.na(available_celltypes) &
    nzchar(available_celltypes) &
    !available_celltypes %in% subcluster_exclude_celltypes
]
if (length(available_celltypes) == 0) stop("对象中没有可重分析的大亚群。")

display_table <- data.frame(
  编号 = seq_along(available_celltypes),
  对象中的大亚群名称 = available_celltypes,
  自动匹配面板 = vapply(
    available_celltypes,
    function(x) {
      answer <- tryCatch(resolve_marker_panel_key(x), error = function(e) NA_character_)
      if (is.na(answer)) "未匹配" else answer
    },
    character(1)
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
print(display_table, row.names = FALSE)

if (is.null(single_recluster_target_celltype) ||
    !nzchar(as.character(single_recluster_target_celltype)[1])) {
  if (!interactive()) {
    stop(
      "非交互运行时请在config/08中填写single_recluster_target_celltype。"
    )
  }
  target_answer <- trimws(readline(
    paste0(
      "请输入要重新分析的大亚群编号、完整名称或关键词",
      "（B/Epi/Mye/T/Fib/Endo/Plasma）："
    )
  ))
} else {
  target_answer <- trimws(as.character(single_recluster_target_celltype)[1])
}

target_number <- suppressWarnings(as.integer(target_answer))
if (length(target_number) == 1 &&
    is.finite(target_number) &&
    target_number >= 1 &&
    target_number <= length(available_celltypes) &&
    identical(as.character(target_number), target_answer)) {
  selected_celltype <- available_celltypes[target_number]
} else {
  selected_celltype <- resolve_requested_celltypes(
    target_answer,
    available_celltypes
  )
  if (length(selected_celltype) != 1) {
    stop("入口4每次只能选择一个大亚群。")
  }
}
message("本次重新分析的大亚群：", selected_celltype)

if (is.null(single_recluster_output_root) ||
    !nzchar(as.character(single_recluster_output_root)[1])) {
  single_recluster_output_root <- file.path(
    result_dir,
    "单大亚群重新聚类与注释"
  )
}
group_stub <- safe_subcluster_name(selected_celltype)
group_output_dir <- file.path(single_recluster_output_root, group_stub)
dir.create(group_output_dir, recursive = TRUE, showWarnings = FALSE)
input_info <- file.info(input_file)
single_run_signature <- paste(
  normalizePath(input_file, winslash = "/", mustWork = TRUE),
  input_info$size,
  as.numeric(input_info$mtime),
  celltype_column,
  selected_celltype,
  unname(tools::md5sum(file.path(
    pipeline_root,
    "config",
    "04_小亚群分析参数.R"
  ))),
  sep = "|"
)
writeLines(
  c(
    paste0("input_file=", normalizePath(input_file, winslash = "/", mustWork = TRUE)),
    paste0("celltype_column=", celltype_column),
    paste0("selected_celltype=", selected_celltype)
  ),
  file.path(group_output_dir, "本次输入与大亚群记录.txt")
)


# ----------------------------------------------------------------------------
# 恢复已确认方案，或反复试算任意PC/res直到用户确认
# ----------------------------------------------------------------------------

pending_object_file <- file.path(
  group_output_dir,
  "最终方案_待人工注释_sce.rds"
)
single_state_file <- file.path(group_output_dir, "单大亚群重分析状态.rds")
previous_single_state <- if (file.exists(single_state_file)) {
  tryCatch(readRDS(single_state_file), error = function(e) NULL)
} else {
  NULL
}
resume_pending <- isTRUE(single_recluster_resume_pending_annotation) &&
  file.exists(pending_object_file) &&
  !is.null(previous_single_state) &&
  identical(previous_single_state$run_signature, single_run_signature)

if (resume_pending) {
  message("检测到已经确认的最终聚类对象，将直接恢复人工注释步骤。")
  accepted_object <- readRDS(pending_object_file)
} else {
  preparation_dir <- file.path(group_output_dir, "00_PCA与Harmony准备")
  dir.create(preparation_dir, recursive = TRUE, showWarnings = FALSE)

  message("提取并重建大亚群对象：", selected_celltype)
  base_object <- prepare_subcluster_object(
    parent_object,
    selected_celltype
  )
  normalized <- normalize_and_reduce_subcluster(
    base_object,
    selected_celltype,
    preparation_dir
  )
  harmony_result <- run_subcluster_harmony(
    normalized$object,
    selected_celltype,
    normalized$maximum_pcs
  )
  base_object <- harmony_result$object
  reduction <- harmony_result$reduction
  maximum_pc <- ncol(Embeddings(base_object, reduction = reduction))
  recommended_pc <- normalized$pc_choice$recommended

  default_pc <- suppressWarnings(as.integer(single_recluster_default_pc)[1])
  if (!is.finite(default_pc) || default_pc < 1 || default_pc > maximum_pc) {
    default_pc <- recommended_pc
  }
  default_resolution <- suppressWarnings(
    as.numeric(single_recluster_default_resolution)[1]
  )
  if (!is.finite(default_resolution) || default_resolution <= 0) {
    default_resolution <- 0.4
  }
  algorithm <- suppressWarnings(as.integer(single_recluster_algorithm)[1])
  if (!algorithm %in% 1:4) {
    stop("single_recluster_algorithm只能填写1、2、3或4。")
  }

  trial_root <- file.path(group_output_dir, "PC_res试算记录")
  existing_trials <- if (dir.exists(trial_root)) {
    list.dirs(trial_root, full.names = FALSE, recursive = FALSE)
  } else {
    character()
  }
  existing_numbers <- suppressWarnings(as.integer(
    sub("^试算_([0-9]+).*", "\\1", existing_trials)
  ))
  trial_number <- if (any(is.finite(existing_numbers))) {
    max(existing_numbers[is.finite(existing_numbers)]) + 1L
  } else {
    1L
  }
  history_file <- file.path(group_output_dir, "PC_res试算历史.csv")
  trial_history <- if (file.exists(history_file)) {
    list(read.csv(
      history_file,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ))
  } else {
    list()
  }

  repeat {
    if (interactive()) {
      selected_pc <- read_single_recluster_pc(default_pc, maximum_pc)
      selected_resolution <- read_single_recluster_resolution(default_resolution)
    } else {
      selected_pc <- suppressWarnings(
        as.integer(single_recluster_noninteractive_pc)[1]
      )
      selected_resolution <- suppressWarnings(
        as.numeric(single_recluster_noninteractive_resolution)[1]
      )
      if (!is.finite(selected_pc) ||
          selected_pc < 1 ||
          selected_pc > maximum_pc ||
          !is.finite(selected_resolution) ||
          selected_resolution <= 0) {
        stop(
          "非交互运行时请在config/08中填写有效的",
          "single_recluster_noninteractive_pc和resolution。"
        )
      }
    }

    message(
      "开始试算：PC 1:", selected_pc,
      " | resolution=", selected_resolution,
      " | ", algorithm_label(algorithm)
    )
    trial_result <- tryCatch(
      run_one_single_recluster_trial(
        base_object = base_object,
        reduction = reduction,
        celltype = selected_celltype,
        pc = selected_pc,
        resolution = selected_resolution,
        algorithm = algorithm,
        trial_number = trial_number,
        group_output_dir = group_output_dir
      ),
      error = function(e) e
    )
    if (inherits(trial_result, "error")) {
      if (!interactive()) stop(conditionMessage(trial_result))
      message("本轮试算失败：", conditionMessage(trial_result))
      message("请重新输入另一组PC/res。")
      trial_number <- trial_number + 1L
      next
    }

    print(trial_result$cluster_summary, row.names = FALSE)
    message("本轮UMAP：", trial_result$umap_file)
    if (isTRUE(single_recluster_open_trial_pdf) && interactive()) {
      try(
        utils::browseURL(normalizePath(
          trial_result$umap_file,
          winslash = "/",
          mustWork = TRUE
        )),
        silent = TRUE
      )
    }

    use_as_final <- if (interactive()) {
      ask_use_single_recluster_trial()
    } else {
      TRUE
    }
    trial_history[[length(trial_history) + 1L]] <- data.frame(
      trial_number = trial_number,
      PC = selected_pc,
      resolution = selected_resolution,
      algorithm = algorithm,
      cluster_number = nrow(trial_result$cluster_summary),
      accepted_as_final = use_as_final,
      trial_dir = normalizePath(
        trial_result$trial_dir,
        winslash = "/",
        mustWork = TRUE
      ),
      stringsAsFactors = FALSE
    )
    history_table <- do.call(rbind, trial_history)
    write.csv(
      history_table,
      history_file,
      row.names = FALSE
    )

    if (use_as_final) {
      accepted_object <- save_accepted_single_recluster_results(
        object = trial_result$object,
        selected = trial_result$selected,
        reduction = reduction,
        celltype = selected_celltype,
        output_dir = group_output_dir
      )
      saveRDS(
        list(
          status = "pending_annotation",
          run_signature = single_run_signature,
          selected_celltype = selected_celltype,
          updated_at = Sys.time()
        ),
        single_state_file
      )
      break
    }
    default_pc <- selected_pc
    default_resolution <- selected_resolution
    trial_number <- trial_number + 1L
    rm(trial_result)
    invisible(gc())
  }
}


# ----------------------------------------------------------------------------
# 最终方案：先验marker点图和单大亚群人工注释
# ----------------------------------------------------------------------------

panel_key <- tryCatch(
  resolve_marker_panel_key(selected_celltype),
  error = function(e) NA_character_
)
if (!is.na(panel_key)) {
  prior_plot_dir <- file.path(group_output_dir, "先验marker点图")
  prior_plot_result <- tryCatch(
    plot_prior_markers_for_celltype(
      object = accepted_object,
      actual_celltype = selected_celltype,
      panel_key = panel_key,
      output_dir = prior_plot_dir
    ),
    error = function(e) e
  )
  if (inherits(prior_plot_result, "error")) {
    warning(
      "先验marker点图生成失败，但不影响人工注释：",
      conditionMessage(prior_plot_result)
    )
  }
} else {
  message("该大亚群没有匹配到先验marker面板，将生成空白人工注释名称。")
}

single_folder_map <- data.frame(
  celltype = selected_celltype,
  output_folder = basename(group_output_dir),
  stringsAsFactors = FALSE
)
expected_annotation_table <- build_subcluster_annotation_table(
  dirname(group_output_dir),
  single_folder_map
)
annotation_workbook <- file.path(
  group_output_dir,
  paste0(group_stub, "_小亚群人工注释表.xlsx")
)
workbook_created <- FALSE
if (!file.exists(annotation_workbook) || !resume_pending) {
  create_subcluster_annotation_workbook(
    expected_annotation_table,
    annotation_workbook
  )
  workbook_created <- TRUE
  message("已经生成单大亚群人工注释表：", annotation_workbook)
}

if (isTRUE(subcluster_annotation_wait_for_edit) && interactive()) {
  readline(paste0(
    "请填写、保存并关闭“", annotation_workbook,
    "”，完成后回到R控制台按回车继续："
  ))
} else if (workbook_created && !interactive()) {
  stop(
    "已经生成单大亚群人工注释表：", annotation_workbook, "\n",
    "请填写并保存后重新运行入口4；程序会恢复待注释对象。",
    call. = FALSE
  )
}

annotation_table <- read_and_validate_subcluster_annotations(
  annotation_workbook,
  expected_annotation_table
)
annotated_object <- apply_subcluster_annotations_to_object(
  accepted_object,
  annotation_table,
  allow_unmapped = FALSE
)
Idents(annotated_object) <- "subcluster_annotation"

final_object_file <- file.path(
  group_output_dir,
  paste0(group_stub, "_重新聚类注释完成_sce.rds")
)
saveRDS(annotated_object, final_object_file)
saveRDS(
  list(
    status = "completed",
    run_signature = single_run_signature,
    selected_celltype = selected_celltype,
    final_object_file = final_object_file,
    updated_at = Sys.time()
  ),
  single_state_file
)
write.csv(
  annotation_table,
  file.path(group_output_dir, "重新聚类后_小亚群注释对应关系.csv"),
  row.names = FALSE
)
cell_mapping <- data.frame(
  cell = colnames(annotated_object),
  celltype = as.character(annotated_object$celltype),
  subcluster = as.character(annotated_object$subcluster),
  subcluster_annotation = as.character(annotated_object$subcluster_annotation),
  subcluster_annotation_label = as.character(
    annotated_object$subcluster_annotation_label
  ),
  stringsAsFactors = FALSE
)
write.csv(
  cell_mapping,
  file.path(group_output_dir, "细胞与重新聚类注释对应关系.csv"),
  row.names = FALSE
)

if ("umap" %in% Reductions(annotated_object)) {
  final_plot <- DimPlot(
    annotated_object,
    reduction = "umap",
    group.by = "subcluster_annotation",
    label = TRUE,
    repel = TRUE,
    shuffle = TRUE,
    raster = FALSE
  ) +
    ggtitle(paste0(selected_celltype, " | 重新聚类与人工注释")) +
    NoLegend()
  ggsave(
    filename = file.path(
      group_output_dir,
      paste0(group_stub, "_重新聚类注释完成_UMAP.pdf")
    ),
    plot = final_plot,
    width = 10,
    height = 8
  )
}
writeLines(
  capture.output(sessionInfo()),
  file.path(group_output_dir, "sessionInfo.txt")
)

cat(
  "\n单大亚群重新聚类与注释已经完成。\n",
  "大亚群：", selected_celltype, "\n",
  "最终对象：", final_object_file, "\n",
  "全部试算和结果：", group_output_dir, "\n",
  sep = ""
)
