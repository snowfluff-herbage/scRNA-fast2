# ==============================================================================
# 仅绘制各大亚群内小亚群先验marker的cluster水平DotPlot
# ==============================================================================
# 前提：已经运行“运行_02_入口_从大亚群注释至小亚群完成.R”。
# 本入口只读取每个大亚群已经聚类完成的子对象，不重复PCA、Harmony、聚类或marker分析。
# 修改config/05_小亚群marker点图参数.R后可随时重新Source本辅助脚本。
# ==============================================================================


# 自动取得本入口脚本所在目录。
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
      "无法取得入口脚本的位置。请不要选中部分代码逐行运行；",
      "请打开“辅助_小亚群marker点图.R”并点击Source，",
      "或者使用Rscript运行该文件。"
    )
  }
  entry_script <- source_files[length(source_files)]
}

entry_script <- normalizePath(entry_script, mustWork = TRUE)
pipeline_root <- dirname(entry_script)

required_pipeline_files <- file.path(
  pipeline_root,
  c(
    "config/01_流程参数.R",
    "config/05_小亚群marker点图参数.R",
    "R/00_载入依赖包.R",
    "R/11_小亚群marker点图函数.R"
  )
)
missing_pipeline_files <- required_pipeline_files[!file.exists(required_pipeline_files)]
if (length(missing_pipeline_files) > 0) {
  stop(
    "小亚群marker点图模块文件不完整。当前定位目录：", pipeline_root,
    "\n缺少文件：\n", paste(missing_pipeline_files, collapse = "\n")
  )
}

message("模块化流程目录：", pipeline_root)
source(file.path(pipeline_root, "config", "01_流程参数.R"), encoding = "UTF-8")
source(file.path(pipeline_root, "R", "00_载入依赖包.R"), encoding = "UTF-8")
source(file.path(pipeline_root, "config", "05_小亚群marker点图参数.R"), encoding = "UTF-8")
source(
  file.path(pipeline_root, "R", "11_小亚群marker点图函数.R"),
  encoding = "UTF-8",
  local = .GlobalEnv
)


# ----------------------------------------------------------------------------
# 定位入口2的结果目录和大亚群文件夹
# ----------------------------------------------------------------------------

if (is.null(subcluster_marker_input_root) ||
    !nzchar(as.character(subcluster_marker_input_root)[1])) {
  subcluster_marker_input_root <- file.path(result_dir, "小亚群自动分析")
}
if (!dir.exists(subcluster_marker_input_root)) {
  stop(
    "没有找到小亚群分析结果目录：", subcluster_marker_input_root,
    "\n请先运行“运行_02_入口_从大亚群注释至小亚群完成.R”，",
    "或在config/05中指定subcluster_marker_input_root。"
  )
}

folder_map_file <- file.path(
  subcluster_marker_input_root,
  "大亚群与输出目录对应关系.csv"
)
if (!file.exists(folder_map_file)) {
  stop("缺少大亚群与输出目录对应关系：", folder_map_file)
}
folder_map <- read.csv(folder_map_file, check.names = FALSE, stringsAsFactors = FALSE)
required_map_columns <- c("celltype", "output_folder")
if (!all(required_map_columns %in% colnames(folder_map))) {
  stop("大亚群与输出目录对应关系.csv必须包含celltype和output_folder列。")
}

if (identical(tolower(as.character(subcluster_marker_target_celltypes)[1]), "all")) {
  target_celltypes <- folder_map$celltype
} else {
  target_celltypes <- unique(as.character(subcluster_marker_target_celltypes))
  missing_targets <- setdiff(target_celltypes, folder_map$celltype)
  if (length(missing_targets) > 0) {
    stop("以下目标大亚群不在入口2结果中：", paste(missing_targets, collapse = "、"))
  }
}


# ----------------------------------------------------------------------------
# 输出完整marker目录、来源和原始表中排除条目，便于审计
# ----------------------------------------------------------------------------

global_output_dir <- file.path(subcluster_marker_input_root, "先验marker点图汇总")
dir.create(global_output_dir, recursive = TRUE, showWarnings = FALSE)

all_catalogs <- lapply(
  names(subcluster_marker_panels),
  build_subcluster_marker_catalog
)
all_marker_catalog <- do.call(rbind, all_catalogs)
rownames(all_marker_catalog) <- NULL
write.csv(
  all_marker_catalog,
  file.path(global_output_dir, "全部大亚群_小亚群先验marker目录.csv"),
  row.names = FALSE
)
write.csv(
  subcluster_marker_reference_sources,
  file.path(global_output_dir, "marker文献来源.csv"),
  row.names = FALSE
)
write.csv(
  subcluster_marker_excluded_entries,
  file.path(global_output_dir, "原始Excel中未进入正式点图的marker.csv"),
  row.names = FALSE
)


# ----------------------------------------------------------------------------
# 逐个大亚群读取子对象并绘制DotPlot
# ----------------------------------------------------------------------------

run_results <- list()
result_index <- 0L

for (current_celltype in target_celltypes) {
  map_row <- folder_map[folder_map$celltype == current_celltype, , drop = FALSE]
  current_folder <- file.path(
    subcluster_marker_input_root,
    map_row$output_folder[1]
  )
  object_file <- file.path(current_folder, "小亚群聚类与marker完成_sce.rds")
  panel_key <- tryCatch(
    resolve_marker_panel_key(current_celltype),
    error = function(e) e
  )

  if (inherits(panel_key, "error")) {
    result_index <- result_index + 1L
    run_results[[result_index]] <- data.frame(
      actual_celltype = current_celltype,
      marker_panel = NA_character_,
      status = "failed_name_mapping",
      message = conditionMessage(panel_key),
      stringsAsFactors = FALSE
    )
    warning(conditionMessage(panel_key))
    next
  }

  if (is.na(panel_key)) {
    result_index <- result_index + 1L
    run_results[[result_index]] <- data.frame(
      actual_celltype = current_celltype,
      marker_panel = NA_character_,
      status = "skipped_no_marker_panel",
      message = "未匹配到marker面板；可在config/05的manual_group_map中设置",
      stringsAsFactors = FALSE
    )
    message("跳过没有marker面板的大亚群：", current_celltype)
    next
  }

  if (!file.exists(object_file)) {
    result_index <- result_index + 1L
    run_results[[result_index]] <- data.frame(
      actual_celltype = current_celltype,
      marker_panel = panel_key,
      status = "skipped_no_subcluster_object",
      message = paste0("未找到子对象：", object_file),
      stringsAsFactors = FALSE
    )
    message("跳过没有独立子对象的大亚群：", current_celltype)
    next
  }

  message("\n", strrep("=", 72))
  message("绘制先验marker DotPlot：", current_celltype, " -> ", panel_key)
  message(strrep("=", 72))
  object <- readRDS(object_file)
  output_dir <- file.path(current_folder, "先验marker点图")

  current_result <- tryCatch(
    plot_prior_markers_for_celltype(
      object = object,
      actual_celltype = current_celltype,
      panel_key = panel_key,
      output_dir = output_dir
    ),
    error = function(e) e
  )
  rm(object)
  invisible(gc())

  result_index <- result_index + 1L
  if (inherits(current_result, "error")) {
    run_results[[result_index]] <- data.frame(
      actual_celltype = current_celltype,
      marker_panel = panel_key,
      status = "failed_plotting",
      message = conditionMessage(current_result),
      stringsAsFactors = FALSE
    )
    warning("大亚群“", current_celltype, "”绘图失败：", conditionMessage(current_result))
  } else {
    run_results[[result_index]] <- data.frame(
      actual_celltype = current_result$actual_celltype,
      marker_panel = current_result$marker_panel,
      status = current_result$status,
      message = paste0(
        "cluster=", current_result$cluster_number,
        "; marker覆盖=", current_result$present_marker_number,
        "/", current_result$configured_marker_number,
        "; pages=", current_result$page_number
      ),
      stringsAsFactors = FALSE
    )
    message("点图已经保存：", current_result$output_pdf)
  }
}

run_summary <- do.call(rbind, run_results)
rownames(run_summary) <- NULL
write.csv(
  run_summary,
  file.path(global_output_dir, "先验marker点图运行汇总.csv"),
  row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(global_output_dir, "sessionInfo.txt"))

completed_number <- sum(run_summary$status == "completed")
failed_number <- sum(grepl("^failed", run_summary$status))
cat(
  "\n小亚群先验marker点图模块运行结束。\n",
  "成功绘图的大亚群：", completed_number, "\n",
  "失败的大亚群：", failed_number, "\n",
  "汇总目录：", global_output_dir, "\n",
  sep = ""
)
