# ==============================================================================
# 入口3使用的小亚群人工注释函数
# 功能：汇总cluster、根据先验marker生成候选名称、创建XLSX、写回Seurat对象。
# ==============================================================================


annotation_key <- function(celltype, cluster) {
  paste(trimws(as.character(celltype)), trimws(as.character(cluster)), sep = "|||")
}


safe_annotation_filename <- function(x) {
  y <- trimws(as.character(x)[1])
  y <- gsub("[\\\\/:*?\"<>|]", "_", y)
  y <- gsub("[[:cntrl:]]", "_", y)
  y <- gsub("[. ]+$", "", y)
  if (!nzchar(y)) y <- "未命名大亚群"
  y
}


rank_prior_subtype_candidates <- function(plot_data_file) {
  empty_result <- data.frame(
    cluster = character(),
    prior_first = character(),
    prior_first_score = numeric(),
    prior_second = character(),
    prior_second_score = numeric(),
    stringsAsFactors = FALSE
  )
  if (!file.exists(plot_data_file)) return(empty_result)

  plot_data <- read.csv(
    plot_data_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  required_columns <- c(
    "cluster", "small_subtype", "preferred_gene", "avg.exp.scaled"
  )
  if (!all(required_columns %in% colnames(plot_data))) return(empty_result)

  plot_data$cluster <- as.character(plot_data$cluster)
  plot_data$small_subtype <- trimws(as.character(plot_data$small_subtype))
  plot_data$avg.exp.scaled <- suppressWarnings(
    as.numeric(plot_data$avg.exp.scaled)
  )
  plot_data <- plot_data[
    nzchar(plot_data$cluster) &
      nzchar(plot_data$small_subtype) &
      is.finite(plot_data$avg.exp.scaled),
    ,
    drop = FALSE
  ]
  if (nrow(plot_data) == 0) return(empty_result)

  score_rows <- lapply(
    split(
      seq_len(nrow(plot_data)),
      annotation_key(plot_data$cluster, plot_data$small_subtype)
    ),
    function(row_index) {
      current <- plot_data[row_index, , drop = FALSE]
      data.frame(
        cluster = current$cluster[1],
        small_subtype = current$small_subtype[1],
        prior_score = mean(current$avg.exp.scaled, na.rm = TRUE),
        marker_number = length(unique(current$preferred_gene)),
        mean_pct_exp = if ("pct.exp" %in% colnames(current)) {
          mean(suppressWarnings(as.numeric(current$pct.exp)), na.rm = TRUE)
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }
  )
  score_table <- do.call(rbind, score_rows)
  score_table$mean_pct_exp[!is.finite(score_table$mean_pct_exp)] <- NA_real_
  score_table <- score_table[
    order(
      score_table$cluster,
      -score_table$prior_score,
      -ifelse(is.na(score_table$mean_pct_exp), -Inf, score_table$mean_pct_exp),
      score_table$small_subtype
    ),
    ,
    drop = FALSE
  ]

  result_rows <- lapply(
    split(score_table, score_table$cluster),
    function(current) {
      current <- current[order(
        -current$prior_score,
        -ifelse(is.na(current$mean_pct_exp), -Inf, current$mean_pct_exp),
        current$small_subtype
      ), , drop = FALSE]
      data.frame(
        cluster = current$cluster[1],
        prior_first = current$small_subtype[1],
        prior_first_score = current$prior_score[1],
        prior_second = if (nrow(current) >= 2) current$small_subtype[2] else NA_character_,
        prior_second_score = if (nrow(current) >= 2) current$prior_score[2] else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  )
  result <- do.call(rbind, result_rows)
  rownames(result) <- NULL
  result
}


build_subcluster_annotation_table <- function(input_root, folder_map) {
  table_rows <- list()
  row_index <- 0L

  for (map_index in seq_len(nrow(folder_map))) {
    current_celltype <- as.character(folder_map$celltype[map_index])
    current_folder <- file.path(
      input_root,
      as.character(folder_map$output_folder[map_index])
    )
    summary_file <- file.path(current_folder, "小亚群细胞数量汇总.csv")
    mapping_file <- file.path(current_folder, "细胞与小亚群对应关系.csv")
    if (file.exists(summary_file)) {
      cluster_summary <- read.csv(
        summary_file,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    } else if (file.exists(mapping_file)) {
      mapping <- read.csv(
        mapping_file,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      cluster_summary <- as.data.frame(
        table(as.character(mapping$subcluster)),
        stringsAsFactors = FALSE
      )
      colnames(cluster_summary) <- c("subcluster", "cell_number")
      cluster_summary$cell_fraction <- cluster_summary$cell_number /
        sum(cluster_summary$cell_number)
    } else {
      stop("大亚群“", current_celltype, "”缺少小亚群汇总和细胞映射文件。")
    }

    required_summary_columns <- c("subcluster", "cell_number")
    if (!all(required_summary_columns %in% colnames(cluster_summary))) {
      stop("小亚群细胞数量汇总.csv缺少subcluster或cell_number列：", current_folder)
    }
    if (!"cell_fraction" %in% colnames(cluster_summary)) {
      cluster_summary$cell_fraction <- cluster_summary$cell_number /
        sum(cluster_summary$cell_number)
    }
    cluster_summary$subcluster <- as.character(cluster_summary$subcluster)

    panel_key <- tryCatch(
      resolve_marker_panel_key(current_celltype),
      error = function(e) NA_character_
    )
    plot_data_file <- file.path(
      current_folder,
      "先验marker点图",
      "先验marker_DotPlot底层数值.csv"
    )
    prior_rank <- rank_prior_subtype_candidates(plot_data_file)
    rank_match <- match(cluster_summary$subcluster, prior_rank$cluster)

    first_candidate <- prior_rank$prior_first[rank_match]
    first_score <- prior_rank$prior_first_score[rank_match]
    second_candidate <- prior_rank$prior_second[rank_match]
    second_score <- prior_rank$prior_second_score[rank_match]
    formal_name <- if (isTRUE(subcluster_annotation_prefill_prior)) {
      first_candidate
    } else {
      rep(NA_character_, nrow(cluster_summary))
    }

    row_index <- row_index + 1L
    table_rows[[row_index]] <- data.frame(
      大亚群 = current_celltype,
      cluster = cluster_summary$subcluster,
      细胞数 = as.integer(cluster_summary$cell_number),
      细胞比例 = as.numeric(cluster_summary$cell_fraction),
      先验marker面板 = if (is.na(panel_key)) "" else panel_key,
      先验推荐亚群名 = ifelse(is.na(first_candidate), "", first_candidate),
      先验推荐得分 = first_score,
      第二候选亚群名 = ifelse(is.na(second_candidate), "", second_candidate),
      第二候选得分 = second_score,
      正式小亚群名称 = ifelse(is.na(formal_name), "", formal_name),
      备注 = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  annotation_table <- do.call(rbind, table_rows)
  rownames(annotation_table) <- NULL
  annotation_table
}


build_prior_subtype_library <- function() {
  rows <- list()
  row_index <- 0L
  for (panel_key in names(subcluster_marker_panels)) {
    subtype_names <- names(subcluster_marker_panels[[panel_key]])
    for (subtype_name in subtype_names) {
      in_user <- !is.null(
        subcluster_marker_uploaded_panels[[panel_key]][[subtype_name]]
      )
      in_literature <- !is.null(
        subcluster_marker_literature_panels[[panel_key]][[subtype_name]]
      )
      source_label <- if (in_user && in_literature) {
        "用户Excel + 文献补充"
      } else if (in_user) {
        "用户Excel"
      } else {
        "文献补充"
      }
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        marker面板 = panel_key,
        候选小亚群名称 = subtype_name,
        来源 = source_label,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}


create_subcluster_annotation_workbook <- function(
    annotation_table,
    workbook_file) {
  workbook <- openxlsx::createWorkbook(creator = "SeuratV5模块化流程")
  openxlsx::addWorksheet(
    workbook,
    "注释填写表",
    gridLines = FALSE,
    tabColour = "#4472C4"
  )
  openxlsx::addWorksheet(
    workbook,
    "先验亚群名称库",
    gridLines = FALSE,
    tabColour = "#70AD47"
  )
  openxlsx::addWorksheet(
    workbook,
    "使用说明",
    gridLines = FALSE,
    tabColour = "#A5A5A5"
  )

  header_style <- openxlsx::createStyle(
    fontColour = "#FFFFFF",
    fgFill = "#4472C4",
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    border = "bottom",
    borderColour = "#D9E2F3"
  )
  editable_style <- openxlsx::createStyle(
    fgFill = "#FFF2CC",
    valign = "top",
    wrapText = TRUE
  )
  percent_style <- openxlsx::createStyle(numFmt = "0.0%")
  score_style <- openxlsx::createStyle(numFmt = "0.000")

  openxlsx::writeData(
    workbook,
    "注释填写表",
    annotation_table,
    withFilter = TRUE,
    headerStyle = header_style,
    keepNA = FALSE
  )
  openxlsx::freezePane(
    workbook,
    "注释填写表",
    firstRow = TRUE,
    firstCol = TRUE
  )
  openxlsx::setColWidths(
    workbook,
    "注释填写表",
    cols = seq_len(ncol(annotation_table)),
    widths = c(18, 10, 11, 11, 18, 28, 14, 28, 14, 30, 30)
  )
  if (nrow(annotation_table) > 0) {
    openxlsx::addStyle(
      workbook,
      "注释填写表",
      percent_style,
      rows = 2:(nrow(annotation_table) + 1),
      cols = which(colnames(annotation_table) == "细胞比例"),
      gridExpand = TRUE
    )
    openxlsx::addStyle(
      workbook,
      "注释填写表",
      score_style,
      rows = 2:(nrow(annotation_table) + 1),
      cols = which(colnames(annotation_table) %in% c("先验推荐得分", "第二候选得分")),
      gridExpand = TRUE
    )
    openxlsx::addStyle(
      workbook,
      "注释填写表",
      editable_style,
      rows = 2:(nrow(annotation_table) + 1),
      cols = which(colnames(annotation_table) %in% c("正式小亚群名称", "备注")),
      gridExpand = TRUE
    )
  }

  prior_library <- build_prior_subtype_library()
  openxlsx::writeData(
    workbook,
    "先验亚群名称库",
    prior_library,
    withFilter = TRUE,
    headerStyle = header_style
  )
  openxlsx::freezePane(
    workbook,
    "先验亚群名称库",
    firstRow = TRUE
  )
  openxlsx::setColWidths(
    workbook,
    "先验亚群名称库",
    cols = 1:3,
    widths = c(20, 34, 24)
  )

  instructions <- data.frame(
    项目 = c(
      "填写位置",
      "预填逻辑",
      "允许的修改",
      "不要修改",
      "空白名称",
      "保存后操作",
      "最终输出"
    ),
    说明 = c(
      "只需要修改“注释填写表”中黄色的“正式小亚群名称”和“备注”列。",
      "先验推荐来自对应大亚群的marker DotPlot平均标准化表达评分，仅用于辅助判读，不代表自动真值。",
      "正式名称可以改成任意中文或英文名称；“先验亚群名称库”中的候选名称也可以自行增加或删除。",
      "不要删除或改变“注释填写表”的大亚群、cluster行，也不要修改这两列的值。",
      if (isTRUE(subcluster_annotation_allow_blank)) {
        paste0("正式名称留空时会自动写为“", subcluster_annotation_unknown_label, "”。")
      } else {
        "正式名称不允许留空，入口3会要求补齐。"
      },
      "保存并关闭Excel后，重新运行入口3；交互运行时也可回到R控制台按回车继续。",
      "所有注释完成的独立对象、合并对象、UMAP和对应关系表会放入同一个输出文件夹。"
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  openxlsx::writeData(
    workbook,
    "使用说明",
    instructions,
    withFilter = FALSE,
    headerStyle = header_style
  )
  openxlsx::setColWidths(workbook, "使用说明", cols = 1, widths = 18)
  openxlsx::setColWidths(workbook, "使用说明", cols = 2, widths = 85)
  openxlsx::addStyle(
    workbook,
    "使用说明",
    openxlsx::createStyle(wrapText = TRUE, valign = "top"),
    rows = 2:(nrow(instructions) + 1),
    cols = 1:2,
    gridExpand = TRUE
  )

  dir.create(dirname(workbook_file), recursive = TRUE, showWarnings = FALSE)
  openxlsx::saveWorkbook(workbook, workbook_file, overwrite = TRUE)
  invisible(workbook_file)
}


read_and_validate_subcluster_annotations <- function(
    workbook_file,
    expected_table) {
  annotation_table <- tryCatch(
    openxlsx::read.xlsx(
      workbook_file,
      sheet = "注释填写表",
      check.names = FALSE,
      detectDates = FALSE
    ),
    error = function(e) {
      stop(
        "无法读取小亚群注释总表。请确认Excel已经保存并关闭：",
        conditionMessage(e)
      )
    }
  )
  required_columns <- c("大亚群", "cluster", "正式小亚群名称")
  missing_columns <- setdiff(required_columns, colnames(annotation_table))
  if (length(missing_columns) > 0) {
    stop("注释填写表缺少必要列：", paste(missing_columns, collapse = "、"))
  }

  annotation_table$大亚群 <- trimws(as.character(annotation_table$大亚群))
  annotation_table$cluster <- trimws(as.character(annotation_table$cluster))
  annotation_table$正式小亚群名称 <- trimws(
    as.character(annotation_table$正式小亚群名称)
  )
  annotation_table$正式小亚群名称[
    is.na(annotation_table$正式小亚群名称)
  ] <- ""
  current_keys <- annotation_key(annotation_table$大亚群, annotation_table$cluster)
  expected_keys <- annotation_key(expected_table$大亚群, expected_table$cluster)

  if (anyDuplicated(current_keys)) {
    duplicated_keys <- unique(current_keys[duplicated(current_keys)])
    stop("注释填写表中存在重复的大亚群/cluster：", paste(duplicated_keys, collapse = "、"))
  }
  missing_keys <- setdiff(expected_keys, current_keys)
  extra_keys <- setdiff(current_keys, expected_keys)
  if (length(missing_keys) > 0 || length(extra_keys) > 0) {
    stop(
      "注释填写表的cluster结构与当前入口2结果不一致。\n",
      if (length(missing_keys) > 0) {
        paste0("缺少：", paste(missing_keys, collapse = "、"), "\n")
      } else {
        ""
      },
      if (length(extra_keys) > 0) {
        paste0("多出：", paste(extra_keys, collapse = "、"), "\n")
      } else {
        ""
      },
      "如需按当前结果重建表格，请把",
      "subcluster_annotation_rebuild_workbook改为TRUE。"
    )
  }

  annotation_table <- annotation_table[
    match(expected_keys, current_keys),
    ,
    drop = FALSE
  ]
  blank_rows <- !nzchar(annotation_table$正式小亚群名称)
  if (any(blank_rows) && !isTRUE(subcluster_annotation_allow_blank)) {
    stop(
      "以下cluster尚未填写正式小亚群名称：",
      paste(
        annotation_key(
          annotation_table$大亚群[blank_rows],
          annotation_table$cluster[blank_rows]
        ),
        collapse = "、"
      )
    )
  }
  annotation_table$正式小亚群名称[blank_rows] <-
    subcluster_annotation_unknown_label
  annotation_table
}


apply_subcluster_annotations_to_object <- function(
    object,
    annotation_table,
    allow_unmapped = FALSE) {
  if (!all(c("celltype", "subcluster") %in% colnames(object@meta.data))) {
    stop("Seurat对象metadata中必须同时包含celltype和subcluster列。")
  }
  lookup_keys <- annotation_key(
    annotation_table$大亚群,
    annotation_table$cluster
  )
  lookup_values <- setNames(
    annotation_table$正式小亚群名称,
    lookup_keys
  )
  object_keys <- annotation_key(object$celltype, object$subcluster)
  final_names <- unname(lookup_values[object_keys])
  if (any(is.na(final_names)) && !isTRUE(allow_unmapped)) {
    missing_keys <- unique(object_keys[is.na(final_names)])
    stop(
      "以下对象中的大亚群/cluster没有注释对应关系：",
      paste(missing_keys, collapse = "、")
    )
  }
  object$subcluster_annotation <- final_names
  final_labels <- rep(NA_character_, length(final_names))
  mapped <- !is.na(final_names)
  final_labels[mapped] <- paste0(
    as.character(object$celltype)[mapped],
    "_",
    final_names[mapped]
  )
  object$subcluster_annotation_label <- final_labels
  object
}
