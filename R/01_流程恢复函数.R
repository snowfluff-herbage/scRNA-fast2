# ==============================================================================
# 流程恢复函数：记录模块状态、保存必要对象、从失败模块继续
# ==============================================================================
# 本文件由主入口自动source，一般不需要用户修改。
# 恢复粒度是“功能块”：某个功能块报错后，下次运行会恢复上一功能块的对象，
# 并从报错功能块重新执行。这样可避免跳过报错前尚未完整写入的操作。

pipeline_checkpoint_dir <- file.path(result_dir, "流程检查点")
pipeline_state_file <- file.path(result_dir, "流程运行状态.rds")
dir.create(pipeline_checkpoint_dir, recursive = TRUE, showWarnings = FALSE)


# 根据模块名称返回检查点文件路径。
pipeline_checkpoint_path <- function(module_file) {
  module_id <- sub("_.*$", "", module_file)
  file.path(
    pipeline_checkpoint_dir,
    paste0(module_id, "_", tools::file_path_sans_ext(module_file), "_完成.rds")
  )
}


# 每个功能块只保存下一个功能块真正需要的对象，避免把全部工作空间写入磁盘。
save_pipeline_checkpoint <- function(module_file, envir = .GlobalEnv) {
  checkpoint_variables <- switch(
    module_file,
    "02_识别并读取输入数据.R" = c(
      "count_list", "sample_name", "input_data_type"
    ),
    "03_metadata与Seurat对象.R" = c(
      "scRNAlist", "sample_info", "sample_name"
    ),
    "04_逐样本质量控制.R" = c(
      "scRNAlist", "sample_info", "sample_name", "qc_summary", "mito_summary"
    ),
    "05_合并样本.R" = c(
      "sce", "sample_info", "sample_name"
    ),
    "06_标准化细胞周期与PCA.R" = c(
      "sce", "recommended_pc_number", "pc_number_candidates",
      "pc_selection_summary", "actual_pca_number"
    ),
    "07_Harmony与正式聚类.R" = c(
      "sce", "recommended_pc_number", "pc_number_candidates",
      "final_pc_number", "final_resolution", "pc.num"
    ),
    "08_寻找cluster_marker.R" = c(
      "sce", "markers", "top10_markers", "heatmap_genes",
      "final_pc_number", "final_resolution"
    ),
    "09_大亚群注释.R" = c(
      "sce", "final_pc_number", "final_resolution"
    ),
    character()
  )

  existing_variables <- checkpoint_variables[
    vapply(checkpoint_variables, exists, logical(1), envir = envir, inherits = FALSE)
  ]

  checkpoint <- list(
    module_file = module_file,
    completed_at = Sys.time(),
    objects = mget(existing_variables, envir = envir, inherits = FALSE)
  )

  saveRDS(checkpoint, pipeline_checkpoint_path(module_file))
  invisible(pipeline_checkpoint_path(module_file))
}


# 将某个已完成模块保存的对象恢复到全局环境，交给下一模块使用。
restore_pipeline_checkpoint <- function(module_file, envir = .GlobalEnv) {
  checkpoint_file <- pipeline_checkpoint_path(module_file)

  if (!file.exists(checkpoint_file)) {
    stop(
      "流程状态要求从检查点恢复，但没有找到文件：", checkpoint_file,
      "\n请把config/01_流程参数.R中的resume_mode改成'restart'，从头运行。"
    )
  }

  checkpoint <- readRDS(checkpoint_file)
  list2env(checkpoint$objects, envir = envir)
  message("已恢复检查点：", basename(checkpoint_file))
  invisible(checkpoint)
}


# 写入当前流程状态。错误信息只保存文字，不保存调用栈和巨大临时对象。
write_pipeline_state <- function(
    status,
    current_module = NA_character_,
    completed_modules = character(),
    error_message = NA_character_) {

  state <- list(
    status = status,
    current_module = current_module,
    completed_modules = completed_modules,
    error_message = error_message,
    updated_at = Sys.time(),
    raw_data_dir = raw_data_dir,
    metadata_file = metadata_file
  )

  saveRDS(state, pipeline_state_file)
  invisible(state)
}


# 报错发生在模块07时，PC或resolution可能已经由用户在控制台选择。
# 将这些小型参数保存到状态旁边，恢复时重新写回，不必再次询问。
save_error_parameters <- function(envir = .GlobalEnv) {
  parameter_names <- c("final_pc_number", "final_resolution")
  existing_parameters <- parameter_names[
    vapply(parameter_names, exists, logical(1), envir = envir, inherits = FALSE)
  ]
  parameter_values <- mget(existing_parameters, envir = envir, inherits = FALSE)
  saveRDS(
    parameter_values,
    file.path(pipeline_checkpoint_dir, "报错时已选择参数.rds")
  )
  invisible(parameter_values)
}


restore_error_parameters <- function(envir = .GlobalEnv) {
  parameter_file <- file.path(pipeline_checkpoint_dir, "报错时已选择参数.rds")
  if (file.exists(parameter_file)) {
    list2env(readRDS(parameter_file), envir = envir)
  }
  invisible(NULL)
}

