# ==============================================================================
# 功能块03：生成/读取样本metadata，并建立逐样本Seurat对象
# 输入：count_list、sample_name
# 输出：sample_info、scRNAlist
# ==============================================================================

# ------------------------------------------------------------------------------
# 3.4.1 读取用户提交的样本级metadata
# ------------------------------------------------------------------------------

# metadata采用“一行代表一个样本”的格式，至少包含sample_id列。例如：
#
# | sample_id       | group  | batch  | patient  | tissue |
# | GSM123_sampleA  | Tumor  | Batch1 | Patient1 | Tumor  |
# | GSM456_sampleB  | Normal | Batch1 | Patient1 | Normal |
#
# sample_id必须与前面得到的sample_name完全一致，包括大小写、下划线和连接符。
# 后面的group、batch、patient、tissue只是示例；用户可以增加任意样本信息列。

# XLSX不依赖CSV的UTF-8/GBK编码，更适合在Windows Excel中填写中文分组信息。
# 本流程用openxlsx同时负责生成和读取样本metadata。
write_sample_metadata <- function(sample_info, metadata_file) {
  if (!grepl("\\.xlsx$", metadata_file, ignore.case = TRUE)) {
    stop(
      "metadata_file必须以.xlsx结尾：", metadata_file,
      "\n请在config/01_流程参数.R中修正metadata_file。"
    )
  }

  metadata_workbook <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(metadata_workbook, "sample_metadata")

  # 表头加深、冻结首行并添加筛选按钮，方便样本较多时在Excel中填写。
  metadata_header_style <- openxlsx::createStyle(
    fontColour = "#FFFFFF",
    fgFill = "#4472C4",
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    border = "Bottom"
  )
  openxlsx::writeData(
    metadata_workbook,
    sheet = "sample_metadata",
    x = sample_info,
    startRow = 1,
    startCol = 1,
    colNames = TRUE,
    rowNames = FALSE,
    withFilter = TRUE,
    keepNA = FALSE,
    headerStyle = metadata_header_style
  )
  openxlsx::freezePane(
    metadata_workbook,
    sheet = "sample_metadata",
    firstRow = TRUE
  )
  openxlsx::setColWidths(
    metadata_workbook,
    sheet = "sample_metadata",
    cols = seq_len(ncol(sample_info)),
    widths = "auto"
  )

  # overwrite=TRUE用于metadata缺少样本时，将补齐后的完整模板写回原文件。
  # 如果文件正被Excel占用，Windows可能拒绝覆盖，此时请关闭Excel后重试。
  tryCatch(
    openxlsx::saveWorkbook(
      metadata_workbook,
      file = metadata_file,
      overwrite = TRUE
    ),
    error = function(e) {
      stop(
        "无法写入metadata XLSX文件：", metadata_file,
        "\n请确认文件没有被Excel占用。",
        "\n原始错误：", conditionMessage(e)
      )
    }
  )

  invisible(metadata_file)
}

read_sample_metadata <- function(metadata_file) {
  if (!grepl("\\.xlsx$", metadata_file, ignore.case = TRUE)) {
    stop(
      "metadata_file必须是.xlsx文件：", metadata_file,
      "\n本版流程不再把CSV作为样本metadata输入。"
    )
  }

  metadata_sheets <- tryCatch(
    openxlsx::getSheetNames(metadata_file),
    error = function(e) {
      stop(
        "无法打开metadata XLSX文件：", metadata_file,
        "\n请确认文件已正常保存且没有损坏。",
        "\n原始错误：", conditionMessage(e)
      )
    }
  )

  if (length(metadata_sheets) == 0) {
    stop("metadata XLSX中没有可读取的工作表：", metadata_file)
  }

  # 自动生成的模板使用sample_metadata工作表。
  # 如果用户对工作表改了名，则读取第一个工作表。
  metadata_sheet <- if ("sample_metadata" %in% metadata_sheets) {
    "sample_metadata"
  } else {
    metadata_sheets[1]
  }

  sample_info <- tryCatch(
    openxlsx::read.xlsx(
      xlsxFile = metadata_file,
      sheet = metadata_sheet,
      colNames = TRUE,
      rowNames = FALSE,
      detectDates = FALSE,
      skipEmptyRows = TRUE,
      skipEmptyCols = FALSE,
      check.names = FALSE,
      na.strings = c("NA", "")
    ),
    error = function(e) {
      stop(
        "无法读取metadata XLSX工作表“", metadata_sheet, "”。",
        "\n请确认第一行是列名，且包含sample_id列。",
        "\n原始错误：", conditionMessage(e)
      )
    }
  )

  # 删除列名两端可能由Excel误输入的空格，但不改动用户的其他列名。
  colnames(sample_info) <- trimws(colnames(sample_info))
  if (!"sample_id" %in% colnames(sample_info)) {
    stop(
      "metadata XLSX工作表“", metadata_sheet,
      "”中没有识别到sample_id列：", metadata_file
    )
  }

  message("metadata读取成功，工作表：", metadata_sheet)
  sample_info
}

add_missing_sample_metadata <- function(sample_info, missing_sample) {
  if (nrow(sample_info) == 0) {
    sample_info <- data.frame(sample_id = character(), stringsAsFactors = FALSE)
  }

  if (!"sample_id" %in% colnames(sample_info)) {
    stop("用户metadata必须包含名为sample_id的列：", metadata_file)
  }

  standard_columns <- c("sample_id", "group", "batch", "patient", "tissue")
  for (current_column in standard_columns) {
    if (!current_column %in% colnames(sample_info)) {
      sample_info[[current_column]] <- ""
    }
  }

  missing_rows <- as.data.frame(
    lapply(sample_info, function(x) rep("", length(missing_sample))),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  missing_rows$sample_id <- missing_sample

  if ("group" %in% colnames(missing_rows)) {
    missing_rows$group <- "Group_to_fill"
  }
  if ("batch" %in% colnames(missing_rows)) {
    missing_rows$batch <- missing_sample
  }
  if ("patient" %in% colnames(missing_rows)) {
    missing_rows$patient <- missing_sample
  }
  if ("tissue" %in% colnames(missing_rows)) {
    missing_rows$tissue <- "Tissue_to_fill"
  }

  rbind(sample_info, missing_rows)
}

# 第一次运行且metadata文件不存在时，根据实际识别到的样本名生成模板。
if (!file.exists(metadata_file)) {
  sample_info_template <- data.frame(
    sample_id = sample_name,
    group = rep("Group_to_fill", length(sample_name)),
    batch = sample_name,
    patient = sample_name,
    tissue = rep("Tissue_to_fill", length(sample_name))
  )

  write_sample_metadata(sample_info_template, metadata_file)

  message("已经自动生成可填写的样本metadata模板：", metadata_file)

  if (interactive()) {
    # 在RStudio/R GUI中运行时，脚本暂停在这里。
    # 用户打开XLSX，填写并保存后回到R控制台按回车，脚本会在本次运行中继续读取。
    readline(paste0(
      "请填写并保存“", metadata_file,
      "”，完成后回到R控制台按回车继续："
    ))
  } else {
    # Rscript或服务器批处理模式无法等待用户现场填写，因此生成模板后安全停止。
    stop(
      "样本metadata模板已经生成：", metadata_file,
      "。当前为非交互式运行，请填写该XLSX文件后重新运行脚本。"
    )
  }
}


# 读取用户填写完成的metadata。
# 如果模板原本已经存在，也会直接从这里读取，不会再次覆盖用户填写内容。
repeat {
  sample_info <- read_sample_metadata(metadata_file)

  if (!"sample_id" %in% colnames(sample_info)) {
    stop("用户metadata必须包含名为sample_id的列：", metadata_file)
  }

  sample_info$sample_id <- trimws(as.character(sample_info$sample_id))
  sample_info <- sample_info[nzchar(sample_info$sample_id), , drop = FALSE]

  if (anyDuplicated(sample_info$sample_id)) {
    stop("用户metadata中的sample_id存在重复值，每个样本只能保留一行。")
  }

  missing_metadata <- setdiff(sample_name, sample_info$sample_id)
  if (length(missing_metadata) == 0) {
    break
  }

  sample_info <- add_missing_sample_metadata(sample_info, missing_metadata)
  write_sample_metadata(sample_info, metadata_file)

  message(
    "metadata中缺少以下样本，已自动补回模板：",
    paste(missing_metadata, collapse = "、")
  )

  if (interactive()) {
    readline(paste0(
      "请重新打开并填写“", metadata_file,
      "”中新补回的样本行，保存后回到R控制台按回车继续："
    ))
  } else {
    stop(
      "metadata中缺少样本，已自动补回模板：", metadata_file,
      "。当前为非交互式运行，请填写后重新运行脚本。"
    )
  }
}

# nCount_RNA和nFeature_RNA由CreateSeuratObject()根据表达矩阵自动计算。
# 如果用户metadata中也包含同名列，会覆盖真实QC信息，因此禁止使用这两个列名。
protected_metadata <- intersect(
  c("nCount_RNA", "nFeature_RNA"),
  colnames(sample_info)
)
if (length(protected_metadata) > 0) {
  stop(
    "用户metadata不能包含以下Seurat自动生成的列：",
    paste(protected_metadata, collapse = "、")
  )
}

# 按照count_list和sample_name的顺序重新排列metadata。
# metadata中存在但本次没有读取的额外样本不会加入Seurat对象。
unused_metadata <- setdiff(sample_info$sample_id, sample_name)
if (length(unused_metadata) > 0) {
  warning(
    "metadata中以下样本未在本次输入数据中找到，已忽略：",
    paste(unused_metadata, collapse = "、")
  )
}
sample_info <- sample_info[
  match(sample_name, sample_info$sample_id),
  ,
  drop = FALSE
]
rownames(sample_info) <- NULL

# 如果示例占位文字仍然存在，脚本继续运行但发出提醒。
# 这些占位值也会进入最终对象，建议正式分析前填写完整。
metadata_text <- unlist(sample_info, use.names = FALSE)
if (any(metadata_text %in% c("Group_to_fill", "Tissue_to_fill"))) {
  warning("样本metadata中仍存在Group_to_fill或Tissue_to_fill，请确认是否已经填写完整。")
}

# 保存本次真正读入并完成样本排序后的metadata，作为分析记录。
# 这个文件同样使用XLSX，便于直接在Excel中查看。
write_sample_metadata(
  sample_info,
  file.path(result_dir, "本次实际读入的样本metadata.xlsx")
)


# ------------------------------------------------------------------------------
# 3.4.2 建立Seurat对象并加入用户metadata的全部列
# ------------------------------------------------------------------------------

scRNAlist <- list()
for (i in seq_along(count_list)) {
  scRNAlist[[i]] <- CreateSeuratObject(
    counts = count_list[[i]],
    project = sample_name[i],
    min.cells = create_min_cells,
    min.features = create_min_features
  )

  # 添加样本前缀，使多个样本中相同的10X barcode在合并后仍保持唯一。
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = sample_name[i])

  # sample_info中一行代表一个样本，而Seurat metadata中一行代表一个细胞。
  # 因此把当前样本的一行信息复制到该样本的所有细胞，再通过AddMetaData加入。
  # 这种写法不需要提前指定列名，用户增加的全部metadata列都会自动保留。
  current_metadata <- sample_info[i, , drop = FALSE]
  cell_metadata <- current_metadata[
    rep(1, ncol(scRNAlist[[i]])),
    ,
    drop = FALSE
  ]
  rownames(cell_metadata) <- colnames(scRNAlist[[i]])

  scRNAlist[[i]] <- AddMetaData(
    object = scRNAlist[[i]],
    metadata = cell_metadata
  )

  rm(current_metadata, cell_metadata)
}

names(scRNAlist) <- sample_name
rm(count_list)
scRNAlist
