# ==============================================================================
# 入口2：从大亚群注释对象开始，完成全部大亚群的小亚群分析
# ==============================================================================
# 使用方法：
#   1. 先完成大亚群注释，得到result/5.大亚群注释完成_sce.rds；
#   2. 修改config/04_小亚群分析参数.R；
#      如使用外部大亚群注释对象，再修改config/07_外部大亚群RDS输入参数.R；
#   3. 在RStudio中打开本文件并点击Source。
#
# 每个大亚群的PC自动使用co1/co2规则的推荐值；程序只会逐个询问resolution。
# 每个大亚群单独保存结果和检查点。某个大亚群报错后，修复问题再次Source本文件，
# 会从该大亚群继续，不重复已成功完成的大亚群。
# ==============================================================================


# 自动取得本入口脚本所在目录，不依赖可能被setwd()改变的工作目录。
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
      "请打开“运行_02_入口_从大亚群注释至小亚群完成.R”并点击Source，",
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
    "config/02_大亚群注释参数.R",
    "config/04_小亚群分析参数.R",
    "config/05_小亚群marker点图参数.R",
    "config/07_外部大亚群RDS输入参数.R",
    "R/00_载入依赖包.R",
    "R/09_大亚群注释.R",
    "R/10_小亚群自动分析函数.R",
    "R/11_小亚群marker点图函数.R"
  )
)
missing_pipeline_files <- required_pipeline_files[!file.exists(required_pipeline_files)]
if (length(missing_pipeline_files) > 0) {
  stop(
    "小亚群分析文件不完整。当前定位目录：", pipeline_root,
    "\n缺少文件：\n", paste(missing_pipeline_files, collapse = "\n")
  )
}

message("模块化流程目录：", pipeline_root)
source(file.path(pipeline_root, "config", "01_流程参数.R"), encoding = "UTF-8")
source(file.path(pipeline_root, "R", "00_载入依赖包.R"), encoding = "UTF-8")
source(file.path(pipeline_root, "config", "02_大亚群注释参数.R"), encoding = "UTF-8")
source(file.path(pipeline_root, "config", "04_小亚群分析参数.R"), encoding = "UTF-8")
source(file.path(pipeline_root, "config", "05_小亚群marker点图参数.R"), encoding = "UTF-8")
source(file.path(pipeline_root, "config", "07_外部大亚群RDS输入参数.R"), encoding = "UTF-8")

if (!requireNamespace("cluster", quietly = TRUE)) {
  stop("缺少R包cluster，无法计算轮廓系数。请安装后重新运行。")
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


# ----------------------------------------------------------------------------
# 读取大亚群注释完成对象并验证
# ----------------------------------------------------------------------------

uses_external_annotated_rds <- isTRUE(downstream_use_external_annotated_rds)
if (uses_external_annotated_rds) {
  external_file <- trimws(as.character(downstream_external_annotated_rds)[1])
  if (!nzchar(external_file)) {
    stop(
      "已经启用外部大亚群注释RDS，但尚未填写",
      "downstream_external_annotated_rds。"
    )
  }
  subcluster_input_file <- external_file
} else if (is.null(subcluster_input_file) ||
           !nzchar(as.character(subcluster_input_file)[1])) {
  subcluster_input_file <- file.path(result_dir, "5.大亚群注释完成_sce.rds")
}

# 如果使用默认大亚群对象，且注释配置比现有对象更新，则先从正式聚类对象快速重做
# 大亚群注释。配置没有变化时不会重复保存对象，避免小亚群检查点被无故判定失效。
default_annotation_file <- file.path(result_dir, "5.大亚群注释完成_sce.rds")
uses_default_annotation_file <- identical(
  normalizePath(subcluster_input_file, winslash = "/", mustWork = FALSE),
  normalizePath(default_annotation_file, winslash = "/", mustWork = FALSE)
)
annotation_config_file <- file.path(
  pipeline_root,
  "config",
  "02_大亚群注释参数.R"
)
annotation_refresh_needed <- !uses_external_annotated_rds &&
  uses_default_annotation_file && (
  !file.exists(default_annotation_file) ||
    file.mtime(annotation_config_file) > file.mtime(default_annotation_file)
)

if (annotation_refresh_needed) {
  clustered_sce_file <- file.path(result_dir, "4.正式聚类后_sce.rds")
  if (!file.exists(clustered_sce_file)) {
    stop(
      "需要从大亚群注释步骤开始，但没有找到正式聚类对象：",
      clustered_sce_file,
      "\n请先运行入口1。"
    )
  }
  message("检测到大亚群注释尚未完成或注释配置已更新，先重新执行大亚群注释。")
  sce <- readRDS(clustered_sce_file)

  final_pc_number <- NA_integer_
  final_resolution <- NA_real_
  clustering_parameter_file <- file.path(result_dir, "正式聚类参数.csv")
  if (file.exists(clustering_parameter_file)) {
    clustering_parameter <- read.csv(
      clustering_parameter_file,
      check.names = FALSE
    )
    if ("final_PC_number" %in% colnames(clustering_parameter)) {
      final_pc_number <- clustering_parameter$final_PC_number[1]
    }
    if ("final_resolution" %in% colnames(clustering_parameter)) {
      final_resolution <- clustering_parameter$final_resolution[1]
    }
  }

  source(
    file.path(pipeline_root, "R", "09_大亚群注释.R"),
    encoding = "UTF-8",
    local = .GlobalEnv
  )
  rm(sce)
  invisible(gc())
}

if (!file.exists(subcluster_input_file)) {
  stop(
    "没有找到大亚群注释完成对象：", subcluster_input_file,
    "\n请先运行“运行_01_入口_至大亚群注释.R”，",
    "或在config/04_小亚群分析参数.R中指定subcluster_input_file。"
  )
}

message("读取大亚群注释对象：", subcluster_input_file)
sce <- readRDS(subcluster_input_file)
if (!inherits(sce, "Seurat")) stop("subcluster_input_file读取后不是Seurat对象。")

configured_celltype_column <- if (uses_external_annotated_rds) {
  trimws(as.character(downstream_external_celltype_column)[1])
} else {
  "celltype"
}
if (!configured_celltype_column %in% colnames(sce@meta.data)) {
  stop(
    "输入Seurat对象的metadata中没有配置的大亚群注释列：",
    configured_celltype_column
  )
}
sce$celltype <- as.character(sce@meta.data[[configured_celltype_column]])
blank_celltype <- is.na(sce$celltype) | !nzchar(trimws(sce$celltype))
if (any(blank_celltype) && !isTRUE(downstream_allow_blank_celltype)) {
  stop(
    "大亚群注释列“", configured_celltype_column, "”中有",
    sum(blank_celltype), "个细胞的名称为空。请先补充注释，或在config/07中允许空值。"
  )
}
sce$celltype[blank_celltype] <- "unknown"

configured_rna_assay <- if (uses_external_annotated_rds) {
  trimws(as.character(downstream_external_rna_assay)[1])
} else {
  "RNA"
}
if (!"RNA" %in% Assays(sce)) {
  if (!configured_rna_assay %in% Assays(sce)) {
    stop(
      "输入Seurat对象中既没有RNA assay，也没有配置的assay：",
      configured_rna_assay
    )
  }
  message("把外部对象的assay“", configured_rna_assay, "”复制为RNA供下游使用。")
  sce[["RNA"]] <- sce[[configured_rna_assay]]
}

all_celltypes <- unique(as.character(sce$celltype))
all_celltypes <- all_celltypes[!is.na(all_celltypes) & nzchar(all_celltypes)]
if (identical(tolower(as.character(subcluster_target_celltypes)[1]), "all")) {
  target_celltypes <- setdiff(all_celltypes, subcluster_exclude_celltypes)
} else {
  target_celltypes <- resolve_requested_celltypes(
    subcluster_target_celltypes,
    all_celltypes
  )
}
if (length(target_celltypes) == 0) stop("没有需要进行小亚群分析的大亚群。")

subcluster_result_dir <- file.path(result_dir, "小亚群自动分析")
dir.create(subcluster_result_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(
  c(
    paste0("input_file=", normalizePath(
      subcluster_input_file,
      winslash = "/",
      mustWork = TRUE
    )),
    paste0("external_rds_enabled=", uses_external_annotated_rds),
    paste0("celltype_column_used=", configured_celltype_column),
    paste0("RNA_assay_available=", "RNA" %in% Assays(sce))
  ),
  file.path(subcluster_result_dir, "下游输入对象记录.txt")
)

# 每完成一个大亚群的小亚群分析，就把独立Seurat对象按大亚群名称复制到统一目录。
# 各大亚群原分析文件夹中的通用检查点文件仍然保留，确保断点续跑逻辑不受影响。
named_subcluster_object_dir <- file.path(
  subcluster_result_dir,
  "按大亚群命名的独立对象"
)
dir.create(named_subcluster_object_dir, recursive = TRUE, showWarnings = FALSE)

sync_named_subcluster_object <- function(current_celltype, current_folder) {
  source_file <- file.path(current_folder, "小亚群聚类与marker完成_sce.rds")
  if (!file.exists(source_file)) return(NA_character_)

  target_file <- file.path(
    named_subcluster_object_dir,
    paste0(basename(current_folder), "_小亚群聚类完成_sce.rds")
  )
  copy_needed <- !file.exists(target_file) ||
    file.mtime(source_file) > file.mtime(target_file)
  if (copy_needed && !file.copy(source_file, target_file, overwrite = TRUE)) {
    stop("无法保存按大亚群命名的独立对象：", target_file)
  }
  target_file
}

# 保存本次流程实际使用的完整先验marker目录、文献来源和排除项。
prior_marker_global_dir <- file.path(subcluster_result_dir, "先验marker点图汇总")
dir.create(prior_marker_global_dir, recursive = TRUE, showWarnings = FALSE)
all_prior_marker_catalogs <- lapply(
  names(subcluster_marker_panels),
  build_subcluster_marker_catalog
)
all_prior_marker_catalog <- do.call(rbind, all_prior_marker_catalogs)
rownames(all_prior_marker_catalog) <- NULL
write.csv(
  all_prior_marker_catalog,
  file.path(prior_marker_global_dir, "全部大亚群_小亚群先验marker目录.csv"),
  row.names = FALSE
)
write.csv(
  subcluster_marker_reference_sources,
  file.path(prior_marker_global_dir, "marker文献来源.csv"),
  row.names = FALSE
)
write.csv(
  subcluster_marker_excluded_entries,
  file.path(prior_marker_global_dir, "原始Excel中未进入正式点图的marker.csv"),
  row.names = FALSE
)

# 处理不同大亚群名称在替换特殊字符后发生重名的情况。
folder_names <- vapply(target_celltypes, safe_subcluster_name, character(1))
if (anyDuplicated(folder_names)) {
  duplicate_groups <- unique(folder_names[duplicated(folder_names)])
  for (duplicate_name in duplicate_groups) {
    indexes <- which(folder_names == duplicate_name)
    folder_names[indexes] <- paste0(folder_names[indexes], "__", seq_along(indexes))
  }
}
celltype_folder_map <- data.frame(
  celltype = target_celltypes,
  output_folder = folder_names,
  matched_marker_panel = vapply(
    target_celltypes,
    function(x) {
      answer <- tryCatch(resolve_marker_panel_key(x), error = function(e) NA_character_)
      if (length(answer) == 0) NA_character_ else answer
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)
write.csv(
  celltype_folder_map,
  file.path(subcluster_result_dir, "大亚群与输出目录对应关系.csv"),
  row.names = FALSE
)


# ----------------------------------------------------------------------------
# 小亚群断点状态
# ----------------------------------------------------------------------------

config_file <- file.path(pipeline_root, "config", "04_小亚群分析参数.R")
external_input_config_file <- file.path(
  pipeline_root,
  "config",
  "07_外部大亚群RDS输入参数.R"
)
analysis_function_file <- file.path(pipeline_root, "R", "10_小亚群自动分析函数.R")
input_info <- file.info(subcluster_input_file)
run_signature <- paste(
  unname(tools::md5sum(config_file)),
  unname(tools::md5sum(external_input_config_file)),
  unname(tools::md5sum(analysis_function_file)),
  normalizePath(subcluster_input_file, winslash = "/", mustWork = TRUE),
  input_info$size,
  as.numeric(input_info$mtime),
  sep = "|"
)

state_file <- file.path(subcluster_result_dir, "小亚群流程运行状态.rds")
completed_celltypes <- character()
previous_state <- NULL
subcluster_internal_resume_allowed <- FALSE
resume_mode <- tolower(as.character(subcluster_resume_mode)[1])
if (!resume_mode %in% c("auto", "restart")) {
  stop("subcluster_resume_mode只能填写'auto'或'restart'。")
}

if (resume_mode == "auto" && file.exists(state_file)) {
  previous_state <- tryCatch(readRDS(state_file), error = function(e) NULL)
  if (!is.null(previous_state) && identical(previous_state$run_signature, run_signature)) {
    subcluster_internal_resume_allowed <- TRUE
    completed_celltypes <- intersect(previous_state$completed_celltypes, target_celltypes)
    if (identical(previous_state$status, "error")) {
      message("检测到上次在大亚群“", previous_state$current_celltype, "”中断。")
      message("上次错误：", previous_state$error_message)
    }
  } else {
    message("检测到配置或输入对象已经变化，旧小亚群检查点自动失效，本次重新计算。")
  }
}

write_subcluster_state <- function(
    status,
    current_celltype = NA_character_,
    error_message = NA_character_) {
  state <- list(
    status = status,
    current_celltype = current_celltype,
    completed_celltypes = completed_celltypes,
    error_message = error_message,
    target_celltypes = target_celltypes,
    run_signature = run_signature,
    input_file = subcluster_input_file,
    updated_at = Sys.time()
  )
  saveRDS(state, state_file)
  invisible(state)
}

# 先验marker配置变化时，只重画DotPlot，不要求重新运行小亚群聚类。
prior_marker_code_files <- c(
  file.path(pipeline_root, "config", "05_小亚群marker点图参数.R"),
  file.path(pipeline_root, "R", "11_小亚群marker点图函数.R")
)
prior_marker_code_mtime <- max(file.info(prior_marker_code_files)$mtime)

run_embedded_prior_marker_plot <- function(
    current_celltype,
    current_folder,
    object = NULL) {
  panel_key <- tryCatch(
    resolve_marker_panel_key(current_celltype),
    error = function(e) e
  )
  if (inherits(panel_key, "error")) {
    warning(
      "大亚群“", current_celltype,
      "”的先验marker面板映射失败：", conditionMessage(panel_key)
    )
    return(data.frame(
      actual_celltype = current_celltype,
      marker_panel = NA_character_,
      status = "failed_name_mapping",
      message = conditionMessage(panel_key),
      stringsAsFactors = FALSE
    ))
  }
  if (is.na(panel_key)) {
    message("没有匹配到先验marker面板，跳过DotPlot：", current_celltype)
    return(data.frame(
      actual_celltype = current_celltype,
      marker_panel = NA_character_,
      status = "skipped_no_marker_panel",
      message = "可在config/05的subcluster_marker_manual_group_map中设置",
      stringsAsFactors = FALSE
    ))
  }

  object_file <- file.path(current_folder, "小亚群聚类与marker完成_sce.rds")
  output_dir <- file.path(current_folder, "先验marker点图")
  output_pdf <- file.path(
    output_dir,
    "先验marker_cluster水平_DotPlot_全部页面.pdf"
  )
  object_mtime <- if (file.exists(object_file)) {
    file.mtime(object_file)
  } else {
    as.POSIXct(0, origin = "1970-01-01")
  }
  newest_input_mtime <- max(prior_marker_code_mtime, object_mtime)

  if (file.exists(output_pdf) &&
      file.mtime(output_pdf) >= newest_input_mtime) {
    message("先验marker DotPlot已存在且为最新，跳过重画：", current_celltype)
    return(data.frame(
      actual_celltype = current_celltype,
      marker_panel = panel_key,
      status = "already_current",
      message = output_pdf,
      stringsAsFactors = FALSE
    ))
  }

  loaded_here <- is.null(object)
  if (loaded_here) {
    if (!file.exists(object_file)) {
      return(data.frame(
        actual_celltype = current_celltype,
        marker_panel = panel_key,
        status = "skipped_no_subcluster_object",
        message = paste0("未找到子对象：", object_file),
        stringsAsFactors = FALSE
      ))
    }
    object <- readRDS(object_file)
  }

  message("自动绘制先验marker DotPlot：", current_celltype, " -> ", panel_key)
  plot_result <- tryCatch(
    plot_prior_markers_for_celltype(
      object = object,
      actual_celltype = current_celltype,
      panel_key = panel_key,
      output_dir = output_dir
    ),
    error = function(e) e
  )
  if (loaded_here) {
    rm(object)
    invisible(gc())
  }

  if (inherits(plot_result, "error")) {
    warning(
      "大亚群“", current_celltype,
      "”的先验marker DotPlot绘制失败，但不影响聚类结果：",
      conditionMessage(plot_result)
    )
    return(data.frame(
      actual_celltype = current_celltype,
      marker_panel = panel_key,
      status = "failed_plotting",
      message = conditionMessage(plot_result),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    actual_celltype = current_celltype,
    marker_panel = panel_key,
    status = "completed",
    message = paste0(
      "cluster=", plot_result$cluster_number,
      "; marker覆盖=", plot_result$present_marker_number,
      "/", plot_result$configured_marker_number,
      "; pages=", plot_result$page_number
    ),
    stringsAsFactors = FALSE
  )
}


# ----------------------------------------------------------------------------
# 依次分析每个大亚群；每完成一个就保存独立结果
# ----------------------------------------------------------------------------

run_summary <- list()
prior_marker_run_results <- list()
for (celltype_index in seq_along(target_celltypes)) {
  current_celltype <- target_celltypes[celltype_index]
  current_folder <- file.path(subcluster_result_dir, folder_names[celltype_index])
  completed_object_file <- file.path(current_folder, "小亚群聚类与marker完成_sce.rds")
  mapping_file <- file.path(current_folder, "细胞与小亚群对应关系.csv")
  skipped_file <- file.path(current_folder, "细胞数不足_未重聚类.csv")

  checkpoint_is_complete <- current_celltype %in% completed_celltypes &&
    file.exists(mapping_file) &&
    (file.exists(completed_object_file) || file.exists(skipped_file))
  if (resume_mode == "auto" && checkpoint_is_complete) {
    message("跳过已完成的大亚群：", current_celltype)
    sync_named_subcluster_object(current_celltype, current_folder)
    prior_marker_run_results[[length(prior_marker_run_results) + 1L]] <-
      run_embedded_prior_marker_plot(
        current_celltype = current_celltype,
        current_folder = current_folder
      )
    next
  }

  message("\n", strrep("=", 72))
  message(
    "开始小亚群分析 [", celltype_index, "/", length(target_celltypes), "]：",
    current_celltype
  )
  message(strrep("=", 72))
  write_subcluster_state(status = "running", current_celltype = current_celltype)

  analysis_result <- tryCatch(
    analyze_one_parent_celltype(
      parent_object = sce,
      celltype = current_celltype,
      output_dir = current_folder
    ),
    error = function(e) e
  )

  if (inherits(analysis_result, "error")) {
    write_subcluster_state(
      status = "error",
      current_celltype = current_celltype,
      error_message = conditionMessage(analysis_result)
    )
    stop(
      "大亚群“", current_celltype, "”的小亚群分析失败。\n",
      "错误信息：", conditionMessage(analysis_result), "\n\n",
      "已完成的大亚群检查点均已保留。修复后重新Source入口2，",
      "会从该大亚群继续。",
      call. = FALSE
    )
  }

  if (identical(analysis_result$status, "skipped_too_few_cells")) {
    dir.create(current_folder, recursive = TRUE, showWarnings = FALSE)
    stale_files <- file.path(
      current_folder,
      c(
        "小亚群聚类与marker完成_sce.rds",
        "小亚群细胞数量汇总.csv",
        "正式小亚群聚类参数.csv",
        "01_PCA_Harmony完成_中间检查点.rds",
        "02_正式聚类完成_待marker_中间检查点.rds"
      )
    )
    unlink(stale_files[file.exists(stale_files)], force = TRUE)
    skipped_cells <- rownames(sce@meta.data)[as.character(sce$celltype) == current_celltype]
    mapping <- data.frame(
      cell = skipped_cells,
      celltype = current_celltype,
      subcluster = "未重聚类",
      subcluster_label = paste0(current_celltype, "_未重聚类"),
      stringsAsFactors = FALSE
    )
    write.csv(mapping, mapping_file, row.names = FALSE)
    write.csv(
      data.frame(
        celltype = current_celltype,
        cell_number = analysis_result$cell_number,
        minimum_required = subcluster_min_cells,
        status = "细胞数不足，未进行重聚类"
      ),
      skipped_file,
      row.names = FALSE
    )
    message("细胞数不足，已跳过重聚类：", analysis_result$cell_number, "个细胞。")
    prior_marker_run_results[[length(prior_marker_run_results) + 1L]] <- data.frame(
      actual_celltype = current_celltype,
      marker_panel = NA_character_,
      status = "skipped_too_few_cells",
      message = paste0("细胞数=", analysis_result$cell_number),
      stringsAsFactors = FALSE
    )
  } else {
    if (file.exists(skipped_file)) unlink(skipped_file, force = TRUE)
    sync_named_subcluster_object(current_celltype, current_folder)
    prior_marker_run_results[[length(prior_marker_run_results) + 1L]] <-
      run_embedded_prior_marker_plot(
        current_celltype = current_celltype,
        current_folder = current_folder,
        object = analysis_result$object
      )
    # 避免同时在内存中保留多个完整子对象；正式对象已经写入各自目录。
    analysis_result$object <- NULL
    run_summary[[current_celltype]] <- analysis_result
  }

  completed_celltypes <- unique(c(completed_celltypes, current_celltype))
  write_subcluster_state(status = "celltype_completed", current_celltype = current_celltype)
  message("大亚群小亚群分析完成并保存：", current_celltype)
  invisible(gc())
}

if (length(prior_marker_run_results) > 0) {
  prior_marker_run_summary <- do.call(rbind, prior_marker_run_results)
  rownames(prior_marker_run_summary) <- NULL
  write.csv(
    prior_marker_run_summary,
    file.path(prior_marker_global_dir, "先验marker点图自动运行汇总.csv"),
    row.names = FALSE
  )
}


# ----------------------------------------------------------------------------
# 将各大亚群的小cluster映射回原始大亚群对象
# ----------------------------------------------------------------------------

all_mapping <- list()
all_cluster_summary <- list()
all_parameter_summary <- list()

for (celltype_index in seq_along(target_celltypes)) {
  current_celltype <- target_celltypes[celltype_index]
  current_folder <- file.path(subcluster_result_dir, folder_names[celltype_index])
  mapping_file <- file.path(current_folder, "细胞与小亚群对应关系.csv")
  if (!file.exists(mapping_file)) {
    stop("汇总时缺少细胞映射表：", mapping_file)
  }
  current_mapping <- read.csv(mapping_file, check.names = FALSE, stringsAsFactors = FALSE)
  all_mapping[[length(all_mapping) + 1L]] <- current_mapping

  summary_file <- file.path(current_folder, "小亚群细胞数量汇总.csv")
  if (file.exists(summary_file)) {
    all_cluster_summary[[length(all_cluster_summary) + 1L]] <- read.csv(
      summary_file,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  } else {
    all_cluster_summary[[length(all_cluster_summary) + 1L]] <- data.frame(
      subcluster = "未重聚类",
      cell_number = nrow(current_mapping),
      celltype = current_celltype,
      cell_fraction = 1,
      stringsAsFactors = FALSE
    )
  }

  parameter_file <- file.path(current_folder, "正式小亚群聚类参数.csv")
  if (file.exists(parameter_file)) {
    all_parameter_summary[[length(all_parameter_summary) + 1L]] <- read.csv(
      parameter_file,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
}

mapping_table <- do.call(rbind, all_mapping)
if (anyDuplicated(mapping_table$cell)) stop("不同大亚群的细胞映射表中出现了重复细胞名。")
write.csv(
  mapping_table,
  file.path(subcluster_result_dir, "全部细胞与小亚群对应关系.csv"),
  row.names = FALSE
)

sce$subcluster <- NA_character_
sce$subcluster_label <- NA_character_
matched_cells <- intersect(rownames(sce@meta.data), mapping_table$cell)
match_index <- match(matched_cells, mapping_table$cell)
sce@meta.data[matched_cells, "subcluster"] <- mapping_table$subcluster[match_index]
sce@meta.data[matched_cells, "subcluster_label"] <- mapping_table$subcluster_label[match_index]

saveRDS(
  sce,
  file.path(subcluster_result_dir, "6.所有大亚群小亚群聚类完成_sce.rds")
)
if (isTRUE(subcluster_save_full_metadata_csv)) {
  write.csv(
    sce@meta.data,
    file.path(subcluster_result_dir, "全部小亚群注释后_metadata.csv")
  )
}

cluster_summary_table <- do.call(rbind, all_cluster_summary)
write.csv(
  cluster_summary_table,
  file.path(subcluster_result_dir, "全部小亚群细胞数量汇总.csv"),
  row.names = FALSE
)
if (length(all_parameter_summary) > 0) {
  parameter_summary_table <- do.call(rbind, all_parameter_summary)
  write.csv(
    parameter_summary_table,
    file.path(subcluster_result_dir, "全部大亚群正式聚类参数汇总.csv"),
    row.names = FALSE
  )
}

named_object_manifest <- data.frame(
  celltype = target_celltypes,
  object_file = vapply(
    seq_along(target_celltypes),
    function(celltype_index) {
      current_folder <- file.path(
        subcluster_result_dir,
        folder_names[celltype_index]
      )
      object_file <- sync_named_subcluster_object(
        target_celltypes[celltype_index],
        current_folder
      )
      if (is.na(object_file)) {
        NA_character_
      } else {
        normalizePath(object_file, winslash = "/", mustWork = TRUE)
      }
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)
write.csv(
  named_object_manifest,
  file.path(named_subcluster_object_dir, "大亚群独立对象文件清单.csv"),
  row.names = FALSE
)

write_subcluster_state(status = "completed", current_celltype = tail(target_celltypes, 1))
writeLines(capture.output(sessionInfo()), file.path(subcluster_result_dir, "sessionInfo.txt"))

cat(
  "\n全部目标大亚群的小亚群分析已完成。\n",
  "目标大亚群：", paste(target_celltypes, collapse = "、"), "\n",
  "合并对象：", file.path(subcluster_result_dir, "6.所有大亚群小亚群聚类完成_sce.rds"), "\n",
  "按大亚群命名的独立对象：", named_subcluster_object_dir, "\n",
  "先验marker点图汇总：", prior_marker_global_dir, "\n",
  "结果目录：", subcluster_result_dir, "\n",
  sep = ""
)
