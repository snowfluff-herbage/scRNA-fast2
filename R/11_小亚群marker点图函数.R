# ==============================================================================
# 大亚群内小亚群先验marker的cluster水平DotPlot函数
# 本文件只定义函数，由入口2自动调用，也可由辅助脚本单独调用。
# ==============================================================================


safe_marker_filename <- function(x) {
  y <- trimws(as.character(x)[1])
  y <- gsub("[\\\\/:*?\"<>|]", "_", y)
  y <- gsub("[[:cntrl:]]", "_", y)
  y <- gsub("[. ]+$", "", y)
  if (!nzchar(y)) y <- "未命名"
  y
}


uploaded_marker_source_file <- function(panel_key) {
  source_map <- c(
    "T/NK" = "T_marker_31(1).xlsx",
    "B cell" = "B_markers(1).xlsx",
    "Plasma" = "PC_markers(1).xlsx",
    "Myeloid" = "Mye_markers(1).xlsx",
    "Fibroblast" = "Fib_markers(1).xlsx",
    "Endothelial" = "EC_markers.xlsx + Fib_markers(1).xlsx"
  )
  answer <- unname(source_map[panel_key])
  if (length(answer) == 0 || is.na(answer)) NA_character_ else answer
}


marker_exists_in_panel <- function(panel_set, panel_key, subtype, gene) {
  panel_group <- panel_set[[panel_key]]
  if (is.null(panel_group) || is.null(panel_group[[subtype]])) return(FALSE)
  gene %in% panel_group[[subtype]]
}


build_subcluster_marker_catalog <- function(panel_key) {
  panel <- subcluster_marker_panels[[panel_key]]
  if (is.null(panel) || length(panel) == 0) {
    stop("没有找到marker面板：", panel_key)
  }

  rows <- list()
  row_index <- 0L
  for (subtype_index in seq_along(panel)) {
    subtype <- names(panel)[subtype_index]
    genes <- unique(as.character(panel[[subtype]]))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    for (gene_index in seq_along(genes)) {
      gene <- genes[gene_index]
      in_uploaded <- marker_exists_in_panel(
        subcluster_marker_uploaded_panels,
        panel_key,
        subtype,
        gene
      )
      in_literature <- marker_exists_in_panel(
        subcluster_marker_literature_panels,
        panel_key,
        subtype,
        gene
      )
      origin <- if (in_uploaded && in_literature) {
        "用户Excel + 文献补充共同支持"
      } else if (in_uploaded) {
        "用户Excel"
      } else {
        "文献补充"
      }

      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        marker_panel = panel_key,
        small_subtype = subtype,
        gene = gene,
        marker_origin = origin,
        uploaded_source_file = if (in_uploaded) uploaded_marker_source_file(panel_key) else NA_character_,
        subtype_order = subtype_index,
        gene_order = gene_index,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}


# 使用名称中的关键词识别常见大亚群，不区分大小写。
# 长关键词优先，避免Fib中的字母B、Endo/Epi中的字母T被短关键词误识别。
canonical_large_group_keyword <- function(celltype_name) {
  original <- trimws(as.character(celltype_name)[1])
  if (is.na(original) || !nzchar(original)) return(NA_character_)
  compact <- gsub("[^a-z0-9]+", "", tolower(original))

  keyword_rules <- list(
    "Plasma" = "plasma",
    "Endothelial" = "endo",
    "Epithelial" = "epi",
    "Myeloid" = "mye",
    "Fibroblast" = "fib|stromal|caf",
    "B cell" = "b",
    "T/NK" = "t"
  )
  for (panel_key in names(keyword_rules)) {
    if (grepl(keyword_rules[[panel_key]], compact, perl = TRUE)) {
      return(panel_key)
    }
  }
  NA_character_
}


resolve_marker_panel_key <- function(actual_celltype) {
  actual_celltype <- trimws(as.character(actual_celltype)[1])

  manual_index <- which(
    tolower(names(subcluster_marker_manual_group_map)) == tolower(actual_celltype)
  )
  if (length(manual_index) == 1) {
    manual_key <- unname(subcluster_marker_manual_group_map[manual_index])
    if (!manual_key %in% names(subcluster_marker_panels)) {
      stop(
        "subcluster_marker_manual_group_map把“", actual_celltype,
        "”映射到了不存在的marker面板：", manual_key
      )
    }
    return(manual_key)
  }

  exact_index <- which(tolower(names(subcluster_marker_panels)) == tolower(actual_celltype))
  if (length(exact_index) == 1) return(names(subcluster_marker_panels)[exact_index])

  alias_matches <- names(subcluster_marker_big_group_aliases)[vapply(
    subcluster_marker_big_group_aliases,
    function(aliases) tolower(actual_celltype) %in% tolower(as.character(aliases)),
    logical(1)
  )]
  alias_matches <- intersect(alias_matches, names(subcluster_marker_panels))
  if (length(alias_matches) == 1) return(alias_matches)
  if (length(alias_matches) > 1) {
    stop("大亚群名称同时匹配了多个marker面板：", paste(alias_matches, collapse = "、"))
  }

  keyword_key <- canonical_large_group_keyword(actual_celltype)
  if (!is.na(keyword_key) && keyword_key %in% names(subcluster_marker_panels)) {
    return(keyword_key)
  }
  NA_character_
}


# 将用户填写的完整大亚群名称或关键词解析为对象中真实存在的celltype名称。
# 精确匹配优先；关键词命中多个真实名称时停止，避免静默选错。
resolve_requested_celltypes <- function(requested_celltypes, available_celltypes) {
  available_celltypes <- unique(trimws(as.character(available_celltypes)))
  available_celltypes <- available_celltypes[
    !is.na(available_celltypes) & nzchar(available_celltypes)
  ]
  resolved <- character()

  for (request in unique(trimws(as.character(requested_celltypes)))) {
    exact <- available_celltypes[
      tolower(available_celltypes) == tolower(request)
    ]
    if (length(exact) == 1) {
      resolved <- c(resolved, exact)
      next
    }

    request_key <- resolve_marker_panel_key(request)
    available_keys <- vapply(
      available_celltypes,
      function(x) {
        answer <- tryCatch(resolve_marker_panel_key(x), error = function(e) NA_character_)
        if (length(answer) == 0) NA_character_ else answer
      },
      character(1)
    )
    keyword_matches <- available_celltypes[
      !is.na(request_key) &
        !is.na(available_keys) &
        available_keys == request_key
    ]
    if (length(keyword_matches) == 1) {
      resolved <- c(resolved, keyword_matches)
      next
    }
    if (length(keyword_matches) > 1) {
      stop(
        "输入“", request, "”同时匹配多个大亚群：",
        paste(keyword_matches, collapse = "、"),
        "。请改用完整名称。"
      )
    }

    partial_matches <- available_celltypes[
      grepl(request, available_celltypes, ignore.case = TRUE, fixed = TRUE)
    ]
    if (length(partial_matches) == 1) {
      resolved <- c(resolved, partial_matches)
      next
    }
    stop(
      "无法把“", request, "”匹配到对象中的大亚群。可用名称：",
      paste(available_celltypes, collapse = "、")
    )
  }
  unique(resolved)
}


choose_marker_cluster_column <- function(object) {
  configured <- as.character(subcluster_marker_cluster_column)[1]
  if (!identical(tolower(configured), "auto")) {
    if (!configured %in% colnames(object@meta.data)) {
      stop("配置的cluster列不在子对象metadata中：", configured)
    }
    return(configured)
  }

  candidates <- c("subcluster", "seurat_clusters")
  available <- candidates[candidates %in% colnames(object@meta.data)]
  if (length(available) == 0) {
    stop("子对象metadata中没有subcluster或seurat_clusters列。")
  }
  available[1]
}


natural_cluster_levels <- function(values) {
  values <- unique(as.character(values))
  values <- values[!is.na(values) & nzchar(values)]
  numeric_values <- suppressWarnings(as.numeric(values))
  if (length(values) > 0 && all(!is.na(numeric_values))) {
    values[order(numeric_values)]
  } else {
    sort(values)
  }
}


match_catalog_genes_to_object <- function(catalog, object) {
  object_genes <- rownames(object)
  upper_object_genes <- toupper(object_genes)
  keep_lookup <- !duplicated(upper_object_genes)
  gene_lookup <- setNames(object_genes[keep_lookup], upper_object_genes[keep_lookup])

  catalog$matched_gene <- unname(gene_lookup[toupper(catalog$gene)])
  catalog$match_status <- ifelse(
    is.na(catalog$matched_gene),
    "missing",
    ifelse(catalog$matched_gene == catalog$gene, "exact", "case_insensitive")
  )
  catalog$present_in_object <- !is.na(catalog$matched_gene)
  catalog
}


expand_dotplot_data_with_catalog <- function(raw_dot_data, present_catalog) {
  rows <- list()
  row_index <- 0L
  raw_features <- as.character(raw_dot_data$features.plot)

  for (catalog_index in seq_len(nrow(present_catalog))) {
    current <- present_catalog[catalog_index, , drop = FALSE]
    matched_rows <- which(raw_features == current$matched_gene)
    if (length(matched_rows) == 0) next
    current_dot <- raw_dot_data[matched_rows, , drop = FALSE]
    current_dot$marker_panel <- current$marker_panel
    current_dot$small_subtype <- current$small_subtype
    current_dot$preferred_gene <- current$gene
    current_dot$matched_gene <- current$matched_gene
    current_dot$marker_origin <- current$marker_origin
    current_dot$subtype_order <- current$subtype_order
    current_dot$gene_order <- current$gene_order
    row_index <- row_index + 1L
    rows[[row_index]] <- current_dot
  }
  if (length(rows) == 0) return(data.frame())
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}


build_marker_dotplot_page <- function(
    plot_data,
    page_catalog,
    cluster_levels,
    actual_celltype,
    panel_key,
    page_number,
    total_pages) {

  subtype_levels <- unique(page_catalog$small_subtype[order(page_catalog$subtype_order)])
  gene_levels <- unique(page_catalog$gene[order(
    page_catalog$subtype_order,
    page_catalog$gene_order
  )])

  plot_data$cluster <- factor(as.character(plot_data$id), levels = cluster_levels)
  plot_data$small_subtype <- factor(plot_data$small_subtype, levels = subtype_levels)
  # 反转gene水平，使每个面板中配置最靠前的marker显示在最上方。
  plot_data$preferred_gene <- factor(plot_data$preferred_gene, levels = rev(gene_levels))

  ggplot(
    plot_data,
    aes(x = cluster, y = preferred_gene)
  ) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    facet_grid(
      rows = vars(small_subtype),
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    scale_size_continuous(
      name = "表达细胞比例 (%)",
      range = c(0, subcluster_marker_dot_scale),
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100)
    ) +
    scale_color_gradientn(
      name = "平均标准化表达",
      colors = subcluster_marker_colors
    ) +
    labs(
      title = paste0(actual_celltype, "：小亚群先验marker（cluster水平）"),
      subtitle = paste0(
        "marker面板：", panel_key,
        " | 第", page_number, "/", total_pages, "页"
      ),
      x = "小亚群cluster",
      y = "先验marker"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(face = "italic", size = 8),
      strip.placement = "outside",
      strip.background = element_rect(fill = "#EAF2F8", color = "#AAB7B8"),
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 8),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
}


plot_prior_markers_for_celltype <- function(
    object,
    actual_celltype,
    panel_key,
    output_dir) {

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!subcluster_marker_assay %in% Assays(object)) {
    stop("子对象中没有配置的assay：", subcluster_marker_assay)
  }
  DefaultAssay(object) <- subcluster_marker_assay

  # 自定义输入对象若仍有多个RNA data层，先合并后再画DotPlot。
  if (inherits(object[[subcluster_marker_assay]], "Assay5")) {
    data_layers <- Layers(object[[subcluster_marker_assay]], search = "^data")
    if (length(data_layers) > 1 ||
        (length(data_layers) == 1 && !"data" %in% data_layers)) {
      object[[subcluster_marker_assay]] <- JoinLayers(object[[subcluster_marker_assay]])
    }
  }

  cluster_column <- choose_marker_cluster_column(object)
  cluster_levels <- natural_cluster_levels(object@meta.data[[cluster_column]])
  if (length(cluster_levels) < 1) stop("cluster列中没有有效分组。")
  plot_width <- max(
    subcluster_marker_pdf_width,
    6 + 0.65 * length(cluster_levels)
  )
  object@meta.data[[cluster_column]] <- factor(
    as.character(object@meta.data[[cluster_column]]),
    levels = cluster_levels
  )

  catalog <- build_subcluster_marker_catalog(panel_key)
  catalog <- match_catalog_genes_to_object(catalog, object)
  write.csv(
    catalog,
    file.path(output_dir, "先验marker基因存在性检查.csv"),
    row.names = FALSE
  )

  subtype_total <- aggregate(
    gene ~ small_subtype,
    data = catalog,
    FUN = length
  )
  colnames(subtype_total)[2] <- "configured_marker_number"
  subtype_present <- aggregate(
    present_in_object ~ small_subtype,
    data = catalog,
    FUN = sum
  )
  colnames(subtype_present)[2] <- "present_marker_number"
  coverage <- merge(subtype_total, subtype_present, by = "small_subtype", all = TRUE)
  coverage$marker_coverage <- coverage$present_marker_number /
    coverage$configured_marker_number
  coverage <- coverage[match(unique(catalog$small_subtype), coverage$small_subtype), ]
  write.csv(
    coverage,
    file.path(output_dir, "各小亚群marker覆盖情况.csv"),
    row.names = FALSE
  )

  present_catalog <- catalog[catalog$present_in_object, , drop = FALSE]
  if (nrow(present_catalog) == 0) {
    stop(
      "marker面板中的基因在对象中全部缺失。",
      "如果rownames是Ensembl ID，请先转换为gene symbol。"
    )
  }

  raw_dot <- DotPlot(
    object = object,
    features = unique(present_catalog$matched_gene),
    assay = subcluster_marker_assay,
    group.by = cluster_column,
    scale = TRUE,
    dot.scale = subcluster_marker_dot_scale
  )
  raw_dot_data <- raw_dot$data
  plot_data <- expand_dotplot_data_with_catalog(raw_dot_data, present_catalog)
  if (nrow(plot_data) == 0) stop("DotPlot已经计算，但无法与marker目录匹配。")

  subtype_names <- unique(present_catalog$small_subtype[order(present_catalog$subtype_order)])
  page_number_by_subtype <- ceiling(
    seq_along(subtype_names) / max(1L, as.integer(subcluster_marker_max_panels_per_page))
  )
  subtype_pages <- split(subtype_names, page_number_by_subtype)

  page_plots <- list()
  page_heights <- numeric(length(subtype_pages))
  plot_data$page <- NA_integer_
  for (page_index in seq_along(subtype_pages)) {
    page_subtypes <- subtype_pages[[page_index]]
    page_catalog <- present_catalog[
      present_catalog$small_subtype %in% page_subtypes,
      ,
      drop = FALSE
    ]
    current_data <- plot_data[
      plot_data$small_subtype %in% page_subtypes,
      ,
      drop = FALSE
    ]
    plot_data$page[plot_data$small_subtype %in% page_subtypes] <- page_index

    page_plots[[page_index]] <- build_marker_dotplot_page(
      plot_data = current_data,
      page_catalog = page_catalog,
      cluster_levels = cluster_levels,
      actual_celltype = actual_celltype,
      panel_key = panel_key,
      page_number = page_index,
      total_pages = length(subtype_pages)
    )
    page_heights[page_index] <- max(
      6,
      subcluster_marker_base_height +
        nrow(page_catalog) * subcluster_marker_height_per_gene
    )

    page_stub <- paste0("先验marker_cluster水平_DotPlot_第", page_index, "页")
    ggsave(
      filename = file.path(output_dir, paste0(page_stub, ".pdf")),
      plot = page_plots[[page_index]],
      width = plot_width,
      height = page_heights[page_index],
      limitsize = FALSE
    )
    if (isTRUE(subcluster_marker_save_page_png)) {
      ggsave(
        filename = file.path(output_dir, paste0(page_stub, ".png")),
        plot = page_plots[[page_index]],
        width = plot_width,
        height = page_heights[page_index],
        dpi = subcluster_marker_png_dpi,
        limitsize = FALSE
      )
    }
  }

  # 同一大亚群的全部页面再保存到一个多页PDF中，便于连续查看。
  multi_page_pdf <- file.path(output_dir, "先验marker_cluster水平_DotPlot_全部页面.pdf")
  grDevices::pdf(
    file = multi_page_pdf,
    width = plot_width,
    height = max(page_heights),
    onefile = TRUE
  )
  tryCatch(
    {
      for (page_plot in page_plots) print(page_plot)
    },
    finally = grDevices::dev.off()
  )

  if (isTRUE(subcluster_marker_save_plot_data)) {
    selected_columns <- intersect(
      c(
        "marker_panel", "small_subtype", "preferred_gene", "matched_gene",
        "marker_origin", "id", "pct.exp", "avg.exp", "avg.exp.scaled", "page"
      ),
      colnames(plot_data)
    )
    output_data <- plot_data[, selected_columns, drop = FALSE]
    if ("id" %in% colnames(output_data)) colnames(output_data)[colnames(output_data) == "id"] <- "cluster"
    write.csv(
      output_data,
      file.path(output_dir, "先验marker_DotPlot底层数值.csv"),
      row.names = FALSE
    )
  }

  list(
    status = "completed",
    actual_celltype = actual_celltype,
    marker_panel = panel_key,
    cluster_column = cluster_column,
    cluster_number = length(cluster_levels),
    configured_marker_number = nrow(catalog),
    present_marker_number = nrow(present_catalog),
    marker_coverage = nrow(present_catalog) / nrow(catalog),
    page_number = length(page_plots),
    output_pdf = multi_page_pdf
  )
}
