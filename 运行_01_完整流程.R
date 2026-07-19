# ==============================================================================
# Seurat V5模块化一站式流程——主入口
# ==============================================================================
# 使用方法：
#   1. 修改config/01_流程参数.R中的路径和参数；
#   2. 在RStudio中打开本文件并点击Source；
#   3. 各功能模块会按照正确顺序自动运行。
#
# 主脚本只负责加载配置和调度模块，具体分析代码保存在R文件夹中。
# ==============================================================================

# 自动取得“本入口脚本”所在目录。
# 不再使用getwd()作为备用目录，因为config文件中的setwd(project_dir)会改变工作目录；
# 如果重复运行旧写法，就可能把project_dir误认为模块化流程所在目录。
command_file <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)

if (length(command_file) > 0) {
  # 使用Rscript运行时，从--file=参数取得当前入口脚本的完整路径。
  entry_script <- sub("^--file=", "", command_file[1])
} else {
  # 在RStudio中点击Source或通过source()运行时，从调用栈中寻找ofile。
  # 使用最后一个有效ofile，使嵌套source时优先取得当前入口脚本，而不是外层脚本。
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
      "请在RStudio中打开“运行_01_完整流程.R”并点击Source，",
      "或者使用Rscript运行该文件。"
    )
  }
  entry_script <- source_files[length(source_files)]
}

entry_script <- normalizePath(entry_script, mustWork = TRUE)
pipeline_root <- dirname(entry_script)

# 定位完成后立即验证目录结构。验证失败时直接停止，不使用其他目录继续猜测。
required_pipeline_files <- file.path(
  pipeline_root,
  c(
    "config/01_流程参数.R",
    "config/02_大亚群注释参数.R",
    "R/00_载入依赖包.R",
    "R/01_流程恢复函数.R",
    "R/02_识别并读取输入数据.R"
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

# 先载入用户参数、依赖包、辅助函数和注释参数。
source(file.path(pipeline_root, "config", "01_流程参数.R"), encoding = "UTF-8")
source(file.path(pipeline_root, "R", "00_载入依赖包.R"), encoding = "UTF-8")
source(file.path(pipeline_root, "R", "01_辅助函数.R"), encoding = "UTF-8")
source(file.path(pipeline_root, "config", "02_大亚群注释参数.R"), encoding = "UTF-8")
source(file.path(pipeline_root, "R", "01_流程恢复函数.R"), encoding = "UTF-8")

# 按照分析顺序依次source。每个模块完成后，产生的对象会直接交给下一个模块。
pipeline_modules <- c(
  "02_识别并读取输入数据.R",
  "03_metadata与Seurat对象.R",
  "04_逐样本质量控制.R",
  "05_合并样本.R",
  "06_标准化细胞周期与PCA.R",
  "07_Harmony与正式聚类.R",
  "08_寻找cluster_marker.R",
  "09_大亚群注释.R"
)

resume_mode <- tolower(resume_mode)
if (!resume_mode %in% c("auto", "restart")) {
  stop("resume_mode只能填写'auto'或'restart'。")
}

start_index <- 1L
completed_modules <- character()
pipeline_already_completed <- FALSE


# ----------------------------------------------------------------------------
# 根据上次运行状态决定从哪个功能块开始
# ----------------------------------------------------------------------------
if (resume_mode == "restart") {
  message("resume_mode='restart'：忽略以前的运行状态，从第一个功能块开始。")
  write_pipeline_state(
    status = "initialized",
    current_module = pipeline_modules[1],
    completed_modules = character()
  )

} else if (file.exists(pipeline_state_file)) {
  previous_state <- tryCatch(
    readRDS(pipeline_state_file),
    error = function(e) NULL
  )

  if (is.null(previous_state)) {
    warning("流程运行状态.rds无法读取，本次从第一个功能块重新开始。")

  } else {
    # 防止更换raw_data_dir后误用旧项目的检查点。
    previous_raw_dir <- tolower(normalizePath(
      previous_state$raw_data_dir,
      winslash = "/",
      mustWork = FALSE
    ))
    current_raw_dir <- tolower(normalizePath(
      raw_data_dir,
      winslash = "/",
      mustWork = FALSE
    ))

    if (!identical(previous_raw_dir, current_raw_dir)) {
      stop(
        "检查点中的raw_data_dir与当前设置不同，不能自动恢复。\n",
        "检查点目录：", previous_state$raw_data_dir, "\n",
        "当前目录：", raw_data_dir, "\n",
        "请确认路径后，把resume_mode改成'restart'重新运行。"
      )
    }

    completed_modules <- intersect(
      previous_state$completed_modules,
      pipeline_modules
    )

    if (identical(previous_state$status, "completed")) {
      pipeline_already_completed <- TRUE
      start_index <- length(pipeline_modules) + 1L
      message(
        "检测到该项目的完整流程已经运行完成。",
        "如需从头重跑，请把resume_mode改成'restart'。"
      )

    } else if (identical(previous_state$status, "module_completed")) {
      completed_index <- match(previous_state$current_module, pipeline_modules)
      if (!is.na(completed_index)) {
        start_index <- completed_index + 1L
      }

    } else {
      # status为error或running都说明当前模块没有完整结束，从当前模块重跑。
      failed_index <- match(previous_state$current_module, pipeline_modules)
      if (!is.na(failed_index)) {
        start_index <- failed_index
        message(
          "检测到上次流程在以下功能块中断：", previous_state$current_module
        )
        if (!is.null(previous_state$error_message) &&
            !is.na(previous_state$error_message)) {
          message("上次错误：", previous_state$error_message)
        }
      }
    }
  }
}


# 极少数情况下，最后一个模块已经保存检查点，但R会话在写入completed状态前退出。
# 此时直接补写完成状态，不再错误尝试恢复一个不存在的“下一个模块”。
if (!pipeline_already_completed && start_index > length(pipeline_modules)) {
  pipeline_already_completed <- TRUE
  write_pipeline_state(
    status = "completed",
    current_module = tail(pipeline_modules, 1),
    completed_modules = pipeline_modules
  )
  message("所有功能块检查点均已存在，流程状态已补记为completed。")
}


# ----------------------------------------------------------------------------
# 恢复失败模块之前已经完整保存的对象
# ----------------------------------------------------------------------------
if (!pipeline_already_completed && start_index > 1L) {
  previous_module <- pipeline_modules[start_index - 1L]

  # 第二处参数覆盖区在原始矩阵读取完成后生效。
  source(
    file.path(pipeline_root, "config", "03_处理前参数覆盖.R"),
    encoding = "UTF-8",
    local = .GlobalEnv
  )

  restore_pipeline_checkpoint(previous_module, envir = .GlobalEnv)
  restore_error_parameters(envir = .GlobalEnv)

  message(
    "将从功能块“", pipeline_modules[start_index], "”继续运行，",
    "此前完成的功能块不会重复计算。"
  )
}


# ----------------------------------------------------------------------------
# 按顺序运行剩余功能块；每个功能块成功后立即保存检查点
# ----------------------------------------------------------------------------
if (!pipeline_already_completed && start_index <= length(pipeline_modules)) {
  modules_to_run <- pipeline_modules[
    seq.int(start_index, length(pipeline_modules))
  ]

  for (module_file in modules_to_run) {
    message("\n", strrep("=", 70))
    message("开始运行功能块：", module_file)
    message(strrep("=", 70))

    write_pipeline_state(
      status = "running",
      current_module = module_file,
      completed_modules = completed_modules
    )

    module_error <- tryCatch({
      source(
        file.path(pipeline_root, "R", module_file),
        encoding = "UTF-8",
        local = .GlobalEnv
      )

      # 原始矩阵读取完成后、建立Seurat对象和QC之前，载入第二处参数覆盖区。
      if (module_file == "02_识别并读取输入数据.R") {
        source(
          file.path(pipeline_root, "config", "03_处理前参数覆盖.R"),
          encoding = "UTF-8",
          local = .GlobalEnv
        )
      }

      save_pipeline_checkpoint(module_file, envir = .GlobalEnv)
      NULL

    }, error = function(e) e)

    if (inherits(module_error, "error")) {
      save_error_parameters(envir = .GlobalEnv)
      write_pipeline_state(
        status = "error",
        current_module = module_file,
        completed_modules = completed_modules,
        error_message = conditionMessage(module_error)
      )

      stop(
        "功能块运行失败：", module_file, "\n",
        "错误信息：", conditionMessage(module_error), "\n\n",
        "流程状态和上一功能块检查点已经保存。修正问题后重新Source本入口脚本，",
        "程序会自动从该功能块继续，不会重复此前已经完成的步骤。",
        call. = FALSE
      )
    }

    completed_modules <- unique(c(completed_modules, module_file))
    write_pipeline_state(
      status = "module_completed",
      current_module = module_file,
      completed_modules = completed_modules
    )
    message("功能块已完成并保存检查点：", module_file)
  }

  write_pipeline_state(
    status = "completed",
    current_module = tail(pipeline_modules, 1),
    completed_modules = pipeline_modules
  )
  message("全部模块已经运行完成。")
}
