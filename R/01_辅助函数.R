# ==============================================================================
# 辅助函数：由主流程自动载入，一般不需要用户修改
# ==============================================================================

# 自动判断features.tsv/genes.tsv中哪一列是真正的基因符号。
# 标准10X文件通常是“Ensembl ID、gene symbol、feature type”三列，gene symbol在第2列；
# 但部分GEO数据会删除Ensembl ID，变成“gene symbol、feature type”两列，
# 此时仍然固定gene.column=2就会把所有基因读成Gene Expression。
# 本函数优先寻找MT-/mt-、核糖体、血红蛋白和管家基因等典型基因符号；
# 如果没有命中，再根据每列唯一值数量排除Gene Expression这种类型列。
detect_10x_gene_column <- function(data_dir) {
  current_files <- list.files(data_dir, full.names = TRUE, recursive = FALSE)
  feature_files <- current_files[
    grepl(
      "(features?|genes?)\\.(tsv|txt|csv)(\\.gz)?$",
      basename(current_files), ignore.case = TRUE
    )
  ]

  # 如果目录中同时存在其他带feature/gene字样的文件，优先使用10X标准名称。
  standard_feature_file <- feature_files[
    tolower(basename(feature_files)) %in% c(
      "features.tsv.gz", "features.tsv",
      "genes.tsv.gz", "genes.tsv"
    )
  ]
  if (length(standard_feature_file) == 1) {
    feature_files <- standard_feature_file
  }

  # 文件尚未整理时可能找不到唯一features文件。
  # 暂时返回Seurat默认的第2列，后面Read10X失败后会先整理目录再重新判断。
  if (length(feature_files) != 1) {
    return(list(
      gene_column = 2L,
      feature_file = NA_character_,
      symbol_score = NA_real_,
      example = "未唯一找到features/genes文件，暂用第2列"
    ))
  }

  feature_table <- data.table::fread(
    feature_files[1],
    header = FALSE,
    data.table = FALSE
  )

  # 典型基因符号评分。正确的gene symbol列通常会命中大量这些基因；
  # Ensembl ID列和Gene Expression类型列一般不会命中。
  symbol_pattern <- paste0(
    "^(MT-|mt-|RPS|RPL|Rps|Rpl|HB[ABDEGMQZ]|Hba|Hbb)",
    "|^(MALAT1|Malat1|ACTB|Actb|GAPDH|Gapdh)$"
  )
  symbol_score <- vapply(
    feature_table,
    function(x) sum(grepl(symbol_pattern, as.character(x)), na.rm = TRUE),
    numeric(1)
  )
  unique_number <- vapply(
    feature_table,
    function(x) length(unique(as.character(x))),
    numeric(1)
  )

  if (max(symbol_score) > 0) {
    gene_column <- which.max(symbol_score)
  } else {
    # 没有典型基因命中时，先排除只有少量唯一值的feature type列。
    possible_gene_column <- which(unique_number > 100)
    if (2 %in% possible_gene_column) {
      gene_column <- 2L
    } else if (length(possible_gene_column) > 0) {
      gene_column <- possible_gene_column[1]
    } else {
      gene_column <- 1L
    }
  }

  return(list(
    gene_column = as.integer(gene_column),
    feature_file = basename(feature_files[1]),
    symbol_score = symbol_score[gene_column],
    example = paste(head(as.character(feature_table[[gene_column]]), 5), collapse = " | ")
  ))
}
# 该规则来自Harvard Chan Bioinformatics Core（HBC）的单细胞聚类培训流程：
#   co1：累计贡献比例超过90%，同时单个PC贡献比例小于5%的第一个位置；
#   co2：相邻PC贡献比例下降超过0.1%的最后一个位置再加1；
#   推荐PC数取co1和co2中的较小值。
#
# 为保持与用户提供代码一致，这里沿用stdev / sum(stdev)计算PC贡献比例。
# choose_pc是本脚本中保留的少量自定义函数之一，仅负责PC数量选择。
choose_pc <- function(seurat_obj) {
  # 读取PCA中每个主成分的标准差，并计算其相对贡献比例和累计比例。
  pct <- seurat_obj[["pca"]]@stdev /
    sum(seurat_obj[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)

  # 满足累计贡献比例>90且单个PC贡献比例<5的第一个PC编号。
  co1 <- which(cumu > 90 & pct < 5)[1]

  # 找出相邻两个PC贡献比例下降超过0.1%的位置。
  # 取最后一个符合位置并加1，作为第二个候选PC数量。
  co2_position <- which(
    (pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1
  )
  if (length(co2_position) > 0) {
    co2 <- sort(co2_position, decreasing = TRUE)[1] + 1
  } else {
    co2 <- NA_integer_
  }

  # 极少数数据可能无法满足co1或co2条件。
  # 如果其中一个为NA，使用另一个；如果都为NA，则退回manual_pc_number。
  available_pc <- c(co1, co2)
  available_pc <- available_pc[!is.na(available_pc)]
  if (length(available_pc) == 0) {
    recommended_pc <- manual_pc_number
  } else {
    recommended_pc <- min(available_pc)
  }

  # 推荐值不能超过实际已经计算的PCA维数，也不能小于1。
  recommended_pc <- max(
    1,
    min(recommended_pc, length(seurat_obj[["pca"]]@stdev))
  )

  cat("co1:", co1, "\n")
  cat("co2:", co2, "\n")
  cat("Recommended PC number:", recommended_pc, "\n")

  return(recommended_pc)
}

