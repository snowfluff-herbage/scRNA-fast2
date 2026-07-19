# ==============================================================================
# 仅重新进行大亚群注释
# ==============================================================================
# 完整流程已经产生“4.正式聚类后_sce.rds”和marker结果后：
#   1. 修改config/02_大亚群注释参数.R；
#   2. Source本文件；
# 即可重新注释，不再重复QC、PCA、Harmony、PC梯度和resolution梯度。
# ==============================================================================

# 自动取得“本入口脚本”所在目录；不使用可能被setwd()改变的getwd()进行猜测。
command_file <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)

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
      "请打开“运行_02_仅重新进行大亚群注释.R”并点击Source，",
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
    "R/00_载入依赖包.R",
    "R/09_大亚群注释.R"
  )
)
missing_pipeline_files <- required_pipeline_files[!file.exists(required_pipeline_files)]
if (length(missing_pipeline_files) > 0) {
  stop(
    "已经定位到入口脚本，但模块化流程文件不完整。当前定位目录：",
    pipeline_root,
    "\n缺少文件：\n",
    paste(missing_pipeline_files, collapse = "\n"),
    "\n请保持config和R文件夹与入口脚本的相对位置不变。"
  )
}

message("模块化流程目录：", pipeline_root)

source(file.path(pipeline_root, "config", "01_流程参数.R"), encoding = "UTF-8")
source(file.path(pipeline_root, "R", "00_载入依赖包.R"), encoding = "UTF-8")
source(file.path(pipeline_root, "config", "02_大亚群注释参数.R"), encoding = "UTF-8")

clustered_sce_file <- file.path(result_dir, "4.正式聚类后_sce.rds")
if (!file.exists(clustered_sce_file)) {
  stop(
    "没有找到正式聚类对象：", clustered_sce_file,
    "。请先运行“运行_01_完整流程.R”。"
  )
}

sce <- readRDS(clustered_sce_file)

# 从完整流程保存的参数记录中恢复正式PC和resolution，仅用于最终运行信息输出。
clustering_parameter_file <- file.path(result_dir, "正式聚类参数.csv")
if (file.exists(clustering_parameter_file)) {
  clustering_parameter <- read.csv(clustering_parameter_file, check.names = FALSE)
  final_pc_number <- clustering_parameter$final_PC_number[1]
  final_resolution <- clustering_parameter$final_resolution[1]
}

source(
  file.path(pipeline_root, "R", "09_大亚群注释.R"),
  encoding = "UTF-8",
  local = .GlobalEnv
)
