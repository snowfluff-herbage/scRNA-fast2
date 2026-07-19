# ==============================================================================
# 功能块02：自动识别并读取H5、10X或CSV原始表达矩阵
# 依赖：config/01_流程参数.R、R/01_辅助函数.R
# 输出：count_list、sample_name
# ==============================================================================

####二、自动识别输入数据是H5、10X目录还是CSV表达矩阵####

# 自动识别只查看raw_data_dir本身以及它的一级子文件夹：
#   1. raw_data_dir下存在.h5或.hdf5文件，识别为h5；
#   2. raw_data_dir本身或一级样本文件夹内存在.mtx/.mtx.gz文件，识别为10x；
#   3. raw_data_dir下存在.csv或.csv.gz文件，识别为csv。
#
# 建议一个raw_data_dir中只保存一种输入格式。这样可避免把样本说明表等普通CSV
# 误认为表达矩阵，也能避免程序在多种格式同时存在时擅自选择错误数据。

top_level_files <- list.files(raw_data_dir, full.names = TRUE, recursive = FALSE)
h5_files <- top_level_files[
  grepl("\\.(h5|hdf5)$", basename(top_level_files), ignore.case = TRUE)
]
csv_files <- top_level_files[
  grepl("\\.csv(\\.gz)?$", basename(top_level_files), ignore.case = TRUE)
]
# 下面这个CSV是本脚本整理10X后自动生成的对应表，不是表达矩阵。
# 将它排除，避免脚本第二次运行时把“10x目录+对应表CSV”误判成两种输入格式。
csv_files <- csv_files[
  !basename(csv_files) %in% c(
    "10X文件与样本改名对应表.csv",
    "10X文件与样本前缀整理对应表.csv"
  )
]
# 10X辅助文件有时使用barcode.csv或features.csv；它们不是完整表达矩阵，
# 因此不参与CSV表达矩阵格式的自动判断。
csv_files <- csv_files[
  !grepl(
    "(barcodes?|features?|genes?)\\.csv(\\.gz)?$",
    basename(csv_files), ignore.case = TRUE
  )
]

# 寻找10X矩阵目录。既支持raw_data_dir本身就是一个10X样本，
# 也支持raw_data_dir/sample_A、raw_data_dir/sample_B这种多样本结构。
sample_dirs <- character()
root_files <- list.files(raw_data_dir, full.names = TRUE, recursive = FALSE)
if (any(grepl("\\.mtx(\\.gz)?$", basename(root_files), ignore.case = TRUE))) {
  sample_dirs <- c(sample_dirs, raw_data_dir)
}

first_level_dirs <- list.dirs(raw_data_dir, recursive = FALSE, full.names = TRUE)
first_level_dirs <- first_level_dirs[!startsWith(basename(first_level_dirs), ".")]
for (i in seq_along(first_level_dirs)) {
  current_files <- list.files(first_level_dirs[i], full.names = TRUE, recursive = FALSE)
  if (any(grepl("\\.mtx(\\.gz)?$", basename(current_files), ignore.case = TRUE))) {
    sample_dirs <- c(sample_dirs, first_level_dirs[i])
  }
}

input_data_type <- tolower(input_data_type)
if (!input_data_type %in% c("auto", "10x", "h5", "csv")) {
  stop("input_data_type只能填写auto、10x、h5或csv。")
}

if (input_data_type == "auto") {
  detected_type <- character()
  if (length(h5_files) > 0) detected_type <- c(detected_type, "h5")
  if (length(sample_dirs) > 0) detected_type <- c(detected_type, "10x")
  if (length(csv_files) > 0) detected_type <- c(detected_type, "csv")

  if (length(detected_type) == 0) {
    stop("没有识别到H5、10X矩阵目录或CSV表达矩阵，请检查raw_data_dir。")
  }
  if (length(detected_type) > 1) {
    stop(
      "raw_data_dir中同时识别到以下格式：", paste(detected_type, collapse = "、"),
      "。请把input_data_type手动改成需要读取的格式。"
    )
  }
  input_data_type <- detected_type
}

message("自动识别/指定的输入数据类型为：", input_data_type)



####三、按照识别结果分别读取H5、10X或CSV数据####

# 三种读取分支最终都会生成：
#   count_list  —— 每个元素是一个样本的基因×细胞原始计数矩阵；
#   sample_name —— 与count_list一一对应的样本名称。
# 后面再使用完全相同的CreateSeuratObject()代码建立Seurat对象。
count_list <- list()
sample_name <- character()


# ------------------------------------------------------------------------------
# 3.1 读取10x Genomics H5文件
# ------------------------------------------------------------------------------
if (input_data_type == "h5") {
  if (length(h5_files) == 0) {
    stop("已指定读取h5，但raw_data_dir下没有找到.h5或.hdf5文件。")
  }
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("读取H5需要hdf5r包，请先运行：install.packages('hdf5r')")
  }

  h5_files <- sort(h5_files)

  # 从文件名中删除常见GEO前缀和10X后缀，得到较简洁的样本名称。
  # make.unique()可以避免多个文件处理后得到完全相同的样本名。
  sample_name <- basename(h5_files)
  sample_name <- sub("\\.(h5|hdf5)$", "", sample_name, ignore.case = TRUE)
  sample_name <- sub("^GSM[0-9]+_", "", sample_name, ignore.case = TRUE)
  sample_name <- sub("_(filtered|raw)_(feature|gene)_bc_matrix$", "", sample_name, ignore.case = TRUE)
  sample_name <- make.unique(sample_name, sep = "_")

  for (i in seq_along(h5_files)) {
    message("正在读取H5样本：", sample_name[i])
    counts <- Read10X_h5(filename = h5_files[i], use.names = TRUE, unique.features = TRUE)

    # 含ADT或CRISPR等多模态数据时，Read10X_h5()会返回list。
    # 本流程分析RNA，因此优先提取名为Gene Expression的矩阵。
    if (is.list(counts)) {
      if ("Gene Expression" %in% names(counts)) {
        counts <- counts[["Gene Expression"]]
      } else {
        counts <- counts[[1]]
        warning("H5文件中未找到Gene Expression，已使用返回列表中的第一个矩阵。")
      }
    }
    count_list[[i]] <- counts
    rm(counts)
  }
}


# ------------------------------------------------------------------------------
# 3.2 读取10X文件夹；读取失败时可选择整理后重新读取
# ------------------------------------------------------------------------------

if (input_data_type == "10x") {
  if (length(sample_dirs) == 0) {
    stop("已指定读取10x，但没有找到含.mtx或.mtx.gz文件的目录。")
  }

  # 不建议同时把一个矩阵放在raw_data_dir根目录、另一些矩阵放在样本子目录。
  # 这种混合结构无法安全判断哪些文件属于同一个样本。
  if (raw_data_dir %in% sample_dirs && length(sample_dirs) > 1) {
    stop("检测到根目录和样本子目录中同时存在mtx文件，请统一为一种目录结构。")
  }

  # 按样本目录名称排序，使每次运行的读取顺序保持一致。
  sample_dirs <- sample_dirs[order(tolower(basename(sample_dirs)))]

  # 先对每个目录实际运行一次Read10X()。
  # 每个样本先自动判断gene.column，再实际读取并检验。
  read_test <- vector("list", length(sample_dirs))
  read_error <- character(length(sample_dirs))
  gene_column_10x <- integer(length(sample_dirs))
  gene_column_summary <- data.frame()

  for (i in seq_along(sample_dirs)) {
    gene_column_result <- detect_10x_gene_column(sample_dirs[i])
    gene_column_10x[i] <- gene_column_result$gene_column

    message(
      "正在检验10X目录：", sample_dirs[i],
      "；自动选择features第", gene_column_10x[i], "列作为基因名"
    )

    gene_column_summary <- rbind(
      gene_column_summary,
      data.frame(
        sample_id = basename(sample_dirs[i]),
        feature_file = gene_column_result$feature_file,
        selected_gene_column = gene_column_10x[i],
        typical_gene_symbol_score = gene_column_result$symbol_score,
        selected_column_example = gene_column_result$example
      )
    )

    read_test[[i]] <- try(
      Read10X(
        data.dir = sample_dirs[i],
        gene.column = gene_column_10x[i],
        unique.features = TRUE
      ),
      silent = TRUE
    )
    if (inherits(read_test[[i]], "try-error")) {
      read_error[i] <- as.character(read_test[[i]])
    }
  }

  failed_sample <- which(vapply(read_test, inherits, logical(1), what = "try-error"))

  if (length(failed_sample) > 0) {
    message("以下10X目录读取失败：")
    message(paste0("  ", sample_dirs[failed_sample], collapse = "\n"))

    organize_choice <- tolower(organize_10x_after_failure)
    if (organize_choice == "ask") {
      if (interactive()) {
        organize_choice <- tolower(trimws(readline(
          "是否开始整理10X目录和文件名？输入y开始整理，输入n停止："
        )))
      } else {
        stop(
          "10X读取失败；当前不是交互式R会话，无法等待键盘选择。",
          "如需自动整理，请把organize_10x_after_failure改成'yes'后重新运行。\n",
          paste(read_error[failed_sample], collapse = "\n")
        )
      }
    }

    if (organize_choice %in% c("yes", "y")) {
      message("开始整理10X样本文件夹和文件名……")

      # 数据库下载的10X文件经常全部平铺在一个目录中，例如：
      #   GSM123_sampleA_barcodes.tsv.gz
      #   GSM123_sampleA_features.tsv.gz
      #   GSM123_sampleA_matrix.mtx.gz
      #   GSM456_sampleB_barcodes.tsv.gz
      #   GSM456_sampleB_features.tsv.gz
      #   GSM456_sampleB_matrix.mtx.gz
      #
      # 下面分别删除文件名末尾的barcode、feature/gene和matrix部分，
      # 保留前面的GSM编号和样本名，并把这个共同前缀作为新样本文件夹名称。
      # 整理前先完成全部匹配；只有每个前缀恰好对应三个文件时才开始移动文件。
      file_manifest <- data.frame()

      for (i in seq_along(sample_dirs)) {
        current_dir <- sample_dirs[i]
        current_files <- list.files(current_dir, full.names = TRUE, recursive = FALSE)

        barcode_files <- current_files[
          grepl("barcodes?\\.(tsv|txt|csv)(\\.gz)?$", basename(current_files), ignore.case = TRUE)
        ]
        feature_files <- current_files[
          grepl("(features?|genes?)\\.(tsv|txt|csv)(\\.gz)?$", basename(current_files), ignore.case = TRUE)
        ]
        matrix_files <- current_files[
          grepl("\\.mtx(\\.gz)?$", basename(current_files), ignore.case = TRUE)
        ]

        # 分别从三类文件名中提取“类型关键词之前”的内容。
        # 末尾多余的下划线、横线、点号和空格会被删除，使文件夹名称更整洁。
        barcode_prefix <- sub(
          "[[:space:]_.-]*barcodes?\\.(tsv|txt|csv)(\\.gz)?$", "",
          basename(barcode_files), ignore.case = TRUE
        )
        feature_prefix <- sub(
          "[[:space:]_.-]*(features?|genes?)\\.(tsv|txt|csv)(\\.gz)?$", "",
          basename(feature_files), ignore.case = TRUE
        )
        matrix_prefix <- sub(
          "([[:space:]_.-]*matrix)?\\.mtx(\\.gz)?$", "",
          basename(matrix_files), ignore.case = TRUE
        )

        barcode_prefix <- sub("[[:space:]_.-]+$", "", barcode_prefix)
        feature_prefix <- sub("[[:space:]_.-]+$", "", feature_prefix)
        matrix_prefix <- sub("[[:space:]_.-]+$", "", matrix_prefix)

        # 如果三个文件本来已经是标准名称，文件名前面没有样本前缀，
        # 则只能在它们已经位于独立样本文件夹时沿用原文件夹名称。
        all_prefix <- c(barcode_prefix, feature_prefix, matrix_prefix)
        if (length(all_prefix) == 3 && all(all_prefix == "")) {
          if (normalizePath(current_dir) == normalizePath(raw_data_dir)) {
            stop(
              "raw_data_dir根目录中的10X文件没有样本前缀，无法据此建立样本文件夹。",
              "请先确认文件是否属于同一个样本。"
            )
          }
          barcode_prefix <- basename(current_dir)
          feature_prefix <- basename(current_dir)
          matrix_prefix <- basename(current_dir)
        }

        sample_prefix <- sort(unique(c(barcode_prefix, feature_prefix, matrix_prefix)))
        sample_prefix <- sample_prefix[nzchar(sample_prefix)]

        for (j in seq_along(sample_prefix)) {
          current_prefix <- sample_prefix[j]
          current_barcode <- barcode_files[barcode_prefix == current_prefix]
          current_feature <- feature_files[feature_prefix == current_prefix]
          current_matrix <- matrix_files[matrix_prefix == current_prefix]

          if (length(current_barcode) != 1 ||
              length(current_feature) != 1 ||
              length(current_matrix) != 1) {
            stop(
              "前缀“", current_prefix, "”没有恰好匹配一个barcode、",
              "一个feature/gene和一个matrix文件，请检查原始文件名。"
            )
          }

          file_manifest <- rbind(
            file_manifest,
            data.frame(
              old_sample_dir = current_dir,
              sample_prefix = current_prefix,
              old_barcode_file = current_barcode,
              old_feature_file = current_feature,
              old_matrix_file = current_matrix,
              new_sample_dir = file.path(raw_data_dir, current_prefix)
            )
          )
        }
      }

      if (nrow(file_manifest) == 0) {
        stop("没有从原始文件名中匹配到完整的10X样本。")
      }
      if (anyDuplicated(file_manifest$sample_prefix)) {
        stop("不同位置出现了相同的样本前缀，请先确认是否为重复文件。")
      }

      # 根据共同前缀建立样本文件夹，再把三个文件整理成Read10X标准名称。
      # 如果源文件已经gzip压缩，只移动并改名；如果未压缩，则真正进行gzip压缩。
      for (i in seq_len(nrow(file_manifest))) {
        target_dir <- file_manifest$new_sample_dir[i]
        dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)

        old_files <- c(
          file_manifest$old_barcode_file[i],
          file_manifest$old_feature_file[i],
          file_manifest$old_matrix_file[i]
        )
        new_files <- file.path(
          target_dir,
          c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
        )

        for (j in seq_along(old_files)) {
          # 文件已经位于目标目录且名称正确时，不重复处理。
          if (normalizePath(old_files[j], mustWork = TRUE) ==
              normalizePath(new_files[j], mustWork = FALSE)) next

          if (file.exists(new_files[j])) {
            stop("目标文件已经存在，为避免覆盖已停止：", new_files[j])
          }

          old_is_gz <- grepl("\\.gz$", old_files[j], ignore.case = TRUE)
          old_is_csv <- grepl("\\.csv(\\.gz)?$", old_files[j], ignore.case = TRUE)

          if (j %in% c(1, 2) && old_is_csv) {
            # barcode或feature为CSV时，转换为无表头的制表符文件并gzip压缩。
            temp_data <- data.table::fread(old_files[j], header = FALSE, data.table = FALSE)
            data.table::fwrite(
              temp_data,
              file = new_files[j],
              sep = "\t",
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE,
              compress = "gzip"
            )
            file.remove(old_files[j])
            rm(temp_data)
          } else if (old_is_gz) {
            if (!file.rename(old_files[j], new_files[j])) {
              stop("文件移动或改名失败：", old_files[j])
            }
          } else {
            if (!requireNamespace("R.utils", quietly = TRUE)) {
              stop("发现未压缩10X文件，请先安装R.utils：install.packages('R.utils')")
            }
            R.utils::gzip(
              old_files[j], destname = new_files[j],
              remove = TRUE, overwrite = FALSE
            )
          }
        }
      }

      sample_dirs <- unique(file_manifest$new_sample_dir)

      # 保存每个共同前缀、原始三个文件及新目录的对应关系，方便追溯。
      write.csv(
        file_manifest,
        file.path(raw_data_dir, "10X文件与样本前缀整理对应表.csv"),
        row.names = FALSE
      )

      # 整理完成后重新实际读取。若此处仍报错，说明问题不只是文件名，
      # 还可能是压缩文件损坏、矩阵维度不一致或features/barcodes内容不正确。
      read_test <- vector("list", length(sample_dirs))
      gene_column_10x <- integer(length(sample_dirs))
      gene_column_summary <- data.frame()

      for (i in seq_along(sample_dirs)) {
        gene_column_result <- detect_10x_gene_column(sample_dirs[i])
        gene_column_10x[i] <- gene_column_result$gene_column

        message(
          "整理后重新读取10X样本：", basename(sample_dirs[i]),
          "；自动选择features第", gene_column_10x[i], "列作为基因名"
        )

        gene_column_summary <- rbind(
          gene_column_summary,
          data.frame(
            sample_id = basename(sample_dirs[i]),
            feature_file = gene_column_result$feature_file,
            selected_gene_column = gene_column_10x[i],
            typical_gene_symbol_score = gene_column_result$symbol_score,
            selected_column_example = gene_column_result$example
          )
        )

        read_test[[i]] <- Read10X(
          data.dir = sample_dirs[i],
          gene.column = gene_column_10x[i],
          unique.features = TRUE
        )
      }
    } else {
      stop(
        "已取消整理，原文件没有被修改。第一次读取的错误为：\n",
        paste(read_error[failed_sample], collapse = "\n")
      )
    }
  }

  sample_name <- basename(sample_dirs)
  sample_name <- make.unique(sample_name, sep = "_")
  count_list <- read_test

  write.csv(
    gene_column_summary,
    file.path(result_dir, "10X基因名称列自动判断结果.csv"),
    row.names = FALSE
  )
}


# ------------------------------------------------------------------------------
# 3.3 读取CSV或CSV.GZ表达矩阵
# ------------------------------------------------------------------------------
if (input_data_type == "csv") {
  if (length(csv_files) == 0) {
    stop("已指定读取csv，但raw_data_dir下没有找到.csv或.csv.gz文件。")
  }

  csv_files <- sort(csv_files)
  sample_name <- basename(csv_files)
  sample_name <- sub("\\.csv(\\.gz)?$", "", sample_name, ignore.case = TRUE)
  sample_name <- sub("^GSM[0-9]+_", "", sample_name, ignore.case = TRUE)
  sample_name <- make.unique(sample_name, sep = "_")

  # CSV必须是“行=基因、列=细胞”的原始计数矩阵：
  #   第一列保存基因名，第一行保存细胞barcode，左上角可以为空或写gene。
  # fread()能够直接读取普通.csv和gzip压缩的.csv.gz，不需要先手动解压。
  for (i in seq_along(csv_files)) {
    message("正在读取CSV样本：", sample_name[i])
    count_table <- data.table::fread(
      csv_files[i], data.table = FALSE, check.names = FALSE
    )

    gene_name <- as.character(count_table[[1]])
    count_table <- count_table[, -1, drop = FALSE]
    count_matrix <- as.matrix(count_table)
    storage.mode(count_matrix) <- "numeric"

    if (anyNA(count_matrix)) {
      stop(
        "CSV转换为数值矩阵后出现NA：", basename(csv_files[i]),
        "。请确认第一列是基因名，其余各列均为原始计数。"
      )
    }

    rownames(count_matrix) <- make.unique(gene_name)
    count_list[[i]] <- Matrix::Matrix(count_matrix, sparse = TRUE)
    rm(count_table, count_matrix, gene_name)
  }
}


# ------------------------------------------------------------------------------
# 3.4 使用三种输入分支的读取结果建立Seurat对象
# ------------------------------------------------------------------------------

# 某些H5或10X目录同时含Gene Expression、抗体和CRISPR矩阵，读取结果是list。
# 本流程只分析转录组，所以统一优先选择Gene Expression。
for (i in seq_along(count_list)) {
  if (is.list(count_list[[i]])) {
    if ("Gene Expression" %in% names(count_list[[i]])) {
      count_list[[i]] <- count_list[[i]][["Gene Expression"]]
    } else {
      count_list[[i]] <- count_list[[i]][[1]]
      warning(sample_name[i], "未找到Gene Expression，已使用返回列表中的第一个矩阵。")
    }
  }
}
