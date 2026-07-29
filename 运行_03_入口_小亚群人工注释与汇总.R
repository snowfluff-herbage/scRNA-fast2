# ==============================================================================
# 入口3：汇总所有小亚群cluster，人工确认注释，并输出最终对象
# ==============================================================================
# 使用方法：
#   1. 先完成入口2；
#   2. Source本文件，程序会生成“小亚群人工注释总表.xlsx”；
#   3. 在Excel中修改黄色的“正式小亚群名称”列，保存并关闭；
#   4. 回到R控制台按回车，或再次Source本入口；
#   5. 所有注释完成文件统一保存到“小亚群注释完成文件”目录。
#
# 第一和第二候选名称根据先验marker DotPlot评分生成，只用于辅助判断。
# 正式名称允许任意修改，不会被候选名称库限制。
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
      "“运行_03_入口_小亚群人工注释与汇总.R”并点击Source。"
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
    "config/05_小亚群marker点图参数.R",
    "config/06_小亚群人工注释参数.R",
    "R/00_载入依赖包.R",
    "R/11_小亚群marker点图函数.R",
    "R/12_小亚群人工注释函数.R"
  )
)
missing_pipeline_files <- required_pipeline_files[!file.exists(required_pipeline_files)]
if (length(missing_pipeline_files) > 0) {
  stop(
    "入口3文件不完整。缺少：\n",
    paste(missing_pipeline_files, collapse = "\n")
  )
}

source(file.path(pipeline_root, "config", "01_流程参数.R"), encoding = "UTF-8")
source(file.path(pipeline_root, "R", "00_载入依赖包.R"), encoding = "UTF-8")
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop(
    "入口3需要R包openxlsx生成和读取注释总表。请先运行：\n",
    "install.packages(\"openxlsx\")"
  )
}
source(
  file.path(pipeline_root, "config", "05_小亚群marker点图参数.R"),
  encoding = "UTF-8"
)
source(
  file.path(pipeline_root, "config", "06_小亚群人工注释参数.R"),
  encoding = "UTF-8"
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


# ----------------------------------------------------------------------------
# 定位入口2结果并生成当前cluster结构
# ----------------------------------------------------------------------------

if (is.null(subcluster_annotation_input_root) ||
    !nzchar(as.character(subcluster_annotation_input_root)[1])) {
  subcluster_annotation_input_root <- file.path(
    result_dir,
    "小亚群自动分析"
  )
}
if (!dir.exists(subcluster_annotation_input_root)) {
  stop(
    "没有找到入口2的小亚群分析目录：",
    subcluster_annotation_input_root,
    "\n请先运行入口2。"
  )
}

folder_map_file <- file.path(
  subcluster_annotation_input_root,
  "大亚群与输出目录对应关系.csv"
)
if (!file.exists(folder_map_file)) {
  stop("缺少大亚群与输出目录对应关系：", folder_map_file)
}
folder_map <- read.csv(
  folder_map_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!all(c("celltype", "output_folder") %in% colnames(folder_map))) {
  stop("大亚群与输出目录对应关系.csv必须包含celltype和output_folder列。")
}

expected_annotation_table <- build_subcluster_annotation_table(
  subcluster_annotation_input_root,
  folder_map
)
if (nrow(expected_annotation_table) == 0) {
  stop("入口2结果中没有可供注释的小亚群cluster。")
}

if (is.null(subcluster_annotation_workbook) ||
    !nzchar(as.character(subcluster_annotation_workbook)[1])) {
  subcluster_annotation_workbook <- file.path(
    subcluster_annotation_input_root,
    "小亚群人工注释总表.xlsx"
  )
}

workbook_created <- FALSE
if (isTRUE(subcluster_annotation_rebuild_workbook) ||
    !file.exists(subcluster_annotation_workbook)) {
  create_subcluster_annotation_workbook(
    expected_annotation_table,
    subcluster_annotation_workbook
  )
  workbook_created <- TRUE
  message("已经生成小亚群人工注释总表：", subcluster_annotation_workbook)
} else {
  message("检测到已有小亚群人工注释总表，将保留现有填写内容：")
  message(subcluster_annotation_workbook)
}

if (isTRUE(subcluster_annotation_wait_for_edit) && interactive()) {
  readline(
    paste0(
      "请填写、保存并关闭“", subcluster_annotation_workbook,
      "”，完成后回到R控制台按回车继续："
    )
  )
} else if (workbook_created && !interactive()) {
  stop(
    "已生成小亚群人工注释总表：", subcluster_annotation_workbook, "\n",
    "请填写并保存后重新运行入口3。",
    call. = FALSE
  )
}

annotation_table <- read_and_validate_subcluster_annotations(
  subcluster_annotation_workbook,
  expected_annotation_table
)


# ----------------------------------------------------------------------------
# 将人工确认的注释写回合并对象和各大亚群独立对象
# ----------------------------------------------------------------------------

if (is.null(subcluster_annotation_output_dir) ||
    !nzchar(as.character(subcluster_annotation_output_dir)[1])) {
  subcluster_annotation_output_dir <- file.path(
    result_dir,
    "小亚群注释完成文件"
  )
}
dir.create(
  subcluster_annotation_output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

combined_input_file <- file.path(
  subcluster_annotation_input_root,
  "6.所有大亚群小亚群聚类完成_sce.rds"
)
if (!file.exists(combined_input_file)) {
  stop("缺少入口2生成的合并Seurat对象：", combined_input_file)
}

message("写入全部小亚群人工注释：", combined_input_file)
combined_object <- readRDS(combined_input_file)
if (!inherits(combined_object, "Seurat")) {
  stop("入口2合并对象读取后不是Seurat对象。")
}
combined_object <- apply_subcluster_annotations_to_object(
  combined_object,
  annotation_table,
  allow_unmapped = TRUE
)

combined_output_file <- file.path(
  subcluster_annotation_output_dir,
  "所有大亚群小亚群注释完成_sce.rds"
)
saveRDS(combined_object, combined_output_file)

annotation_output_table <- annotation_table
annotation_output_table$注释键 <- annotation_key(
  annotation_output_table$大亚群,
  annotation_output_table$cluster
)
write.csv(
  annotation_output_table,
  file.path(subcluster_annotation_output_dir, "小亚群人工注释对应关系.csv"),
  row.names = FALSE
)

if (isTRUE(subcluster_annotation_save_full_metadata_csv)) {
  write.csv(
    combined_object@meta.data,
    file.path(
      subcluster_annotation_output_dir,
      "所有细胞小亚群注释完成_metadata.csv"
    )
  )
}

named_object_input_dir <- file.path(
  subcluster_annotation_input_root,
  "按大亚群命名的独立对象"
)
manifest_rows <- list()

for (map_index in seq_len(nrow(folder_map))) {
  current_celltype <- as.character(folder_map$celltype[map_index])
  safe_name <- as.character(folder_map$output_folder[map_index])
  current_folder <- file.path(
    subcluster_annotation_input_root,
    as.character(folder_map$output_folder[map_index])
  )
  named_input_file <- file.path(
    named_object_input_dir,
    paste0(safe_name, "_小亚群聚类完成_sce.rds")
  )
  generic_input_file <- file.path(
    current_folder,
    "小亚群聚类与marker完成_sce.rds"
  )

  if (file.exists(named_input_file)) {
    current_object <- readRDS(named_input_file)
  } else if (file.exists(generic_input_file)) {
    current_object <- readRDS(generic_input_file)
  } else {
    selected_cells <- rownames(combined_object@meta.data)[
      as.character(combined_object$celltype) == current_celltype
    ]
    if (length(selected_cells) == 0) {
      warning("大亚群没有可保存的细胞，已跳过：", current_celltype)
      next
    }
    current_object <- subset(combined_object, cells = selected_cells)
  }

  current_annotation <- annotation_table[
    annotation_table$大亚群 == current_celltype,
    ,
    drop = FALSE
  ]
  current_object <- apply_subcluster_annotations_to_object(
    current_object,
    current_annotation,
    allow_unmapped = FALSE
  )
  Idents(current_object) <- "subcluster_annotation"

  current_output_file <- file.path(
    subcluster_annotation_output_dir,
    paste0(safe_name, "_小亚群注释完成_sce.rds")
  )
  saveRDS(current_object, current_output_file)

  umap_output_file <- NA_character_
  if ("umap" %in% Reductions(current_object)) {
    annotation_plot <- DimPlot(
      current_object,
      reduction = "umap",
      group.by = "subcluster_annotation",
      label = TRUE,
      repel = TRUE,
      shuffle = TRUE,
      raster = FALSE
    ) +
      ggtitle(paste0(current_celltype, " | 小亚群人工注释")) +
      NoLegend()
    umap_output_file <- file.path(
      subcluster_annotation_output_dir,
      paste0(safe_name, "_小亚群注释完成_UMAP.pdf")
    )
    ggsave(
      filename = umap_output_file,
      plot = annotation_plot,
      width = 10,
      height = 8
    )
  }

  manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
    大亚群 = current_celltype,
    细胞数 = ncol(current_object),
    小亚群名称数 = length(unique(as.character(current_object$subcluster_annotation))),
    注释完成对象 = normalizePath(
      current_output_file,
      winslash = "/",
      mustWork = TRUE
    ),
    注释UMAP = if (is.na(umap_output_file)) {
      NA_character_
    } else {
      normalizePath(umap_output_file, winslash = "/", mustWork = TRUE)
    },
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rm(current_object)
  invisible(gc())
}

if (length(manifest_rows) > 0) {
  output_manifest <- do.call(rbind, manifest_rows)
  write.csv(
    output_manifest,
    file.path(subcluster_annotation_output_dir, "小亚群注释完成文件清单.csv"),
    row.names = FALSE
  )
}

filled_workbook_copy <- file.path(
  subcluster_annotation_output_dir,
  "小亚群人工注释总表_已用于本次注释.xlsx"
)
if (!file.copy(
  subcluster_annotation_workbook,
  filled_workbook_copy,
  overwrite = TRUE
)) {
  warning("无法复制填写后的注释总表到最终输出目录。")
}

writeLines(
  capture.output(sessionInfo()),
  file.path(subcluster_annotation_output_dir, "sessionInfo.txt")
)

cat(
  "\n全部小亚群人工注释已经完成。\n",
  "注释总表：", subcluster_annotation_workbook, "\n",
  "合并对象：", combined_output_file, "\n",
  "全部注释完成文件：", subcluster_annotation_output_dir, "\n",
  sep = ""
)
