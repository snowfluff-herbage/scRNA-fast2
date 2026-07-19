# ==============================================================================
# 单细胞标准聚类与大亚群注释流程（Seurat V5）
# ==============================================================================
# 本脚本按照“从头读取数据”的顺序编写，整体代码风格参考：
# xzx(从头）单细胞标准聚类注释流程SeuratV5全流程.R
#
# 脚本目前写到“大亚群注释”结束，不包含以下下游分析：
#   细胞比例比较、差异表达、富集分析、拟时序、细胞通讯、SCENIC、inferCNV、
#   CytoTRACE、细胞亚群再次聚类等。
#
# 运行顺序：
#   一、设置路径和参数
#   二、自动识别输入数据类型
#   三、分别读取H5、10X或CSV表达矩阵
#       如果10X读取失败，可按文件共同前缀建立样本目录并重新读取
#   四、逐样本质控
#   五、合并样本和补充metadata
#   六、标准化、高变基因、细胞周期和PCA
#   七、Harmony去批次
#   八、选择PC、聚类分辨率和UMAP/tSNE
#   九、寻找cluster marker
#   十、依据经典marker完成大亚群注释
#
# 注意：代码中不包含任何原队列的路径、样本名称、分组或cluster注释结果。
# ==============================================================================



####一、设置工作目录、载入R包和修改参数####

# Windows路径建议统一使用正斜杠“/”，避免反斜杠被R识别为转义符。
# 请把下面三个路径改成自己的真实路径。
project_dir <- "C:/singlecell/project"
raw_data_dir <- "C:/singlecell/project/rawdata"
result_dir <- "C:/singlecell/project/result"

dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
setwd(project_dir)
getwd()

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 30 * 1024^3)  #允许future使用最多30GB内存
set.seed(1234)

library(Seurat)
library(SeuratObject)
library(Matrix)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)
library(data.table)


# ----需要根据数据修改的参数----

# 输入数据类型。默认使用"auto"自动识别，也可以手动指定为：
#   "10x"：标准10X样本文件夹，或带共同样本前缀并平铺保存的三个矩阵文件；
#   "h5" ：一个或多个10x Genomics生成的.h5/.hdf5文件；
#   "csv"：一个或多个“基因×细胞”的.csv或.csv.gz表达矩阵。
# 如果raw_data_dir中同时放了多种格式，自动模式无法判断应读取哪一种，
# 此时请手动填写"10x"、"h5"或"csv"。
input_data_type <- "auto"

# 仅在10X目录第一次Read10X()读取失败时使用。
#   "ask"：在R控制台询问是否开始整理；
#   "yes"：不询问，直接整理；
#   "no" ：不修改任何文件，直接停止并显示读取错误。
# 整理时不再把样本命名为P1、P2、P3……，而是提取barcode、feature/gene、
# matrix前面的共同文件名前缀，以该前缀建立样本文件夹。例如：
#   GSM123_sampleA_barcodes.tsv.gz
#   GSM123_sampleA_features.tsv.gz
#   GSM123_sampleA_matrix.mtx.gz
# 会建立GSM123_sampleA文件夹，并把三个文件放入其中、改成10X标准文件名。
organize_10x_after_failure <- "ask"

# 用户填写的“样本级metadata”文件路径。
# 文件必须是CSV格式，并且至少包含一列名为sample_id。
# sample_id需要与本脚本读取后得到的sample_name完全一致。
# 其他列可以由用户自由增加，例如group、batch、patient、tissue、condition等；
# 这些列都会自动添加到最终sce@meta.data中。
metadata_file <- file.path(project_dir, "样本metadata.csv")

# 建立Seurat对象时的最低要求。
create_min_cells <- 3
create_min_features <- 100

# 正式QC阈值。请结合QC小提琴图修改，不要机械照搬。
min_features <- 200
max_features <- 8000
min_counts <- 500
max_counts <- 100000
max_percent_mito <- 20
max_percent_hb <- 5
min_percent_ribo <- 0

# 标准化、高变基因、PCA参数。
n_variable_features <- 3000
n_pcs_calculate <- 50

# 是否使用HBC培训流程中的co1/co2规则自动选择聚类所用PC数。
# TRUE：PCA完成后计算推荐PC数，作为PC梯度比较的参考；
# FALSE：不计算co1/co2，使用manual_pc_number作为参考值。
auto_choose_pc <- TRUE
manual_pc_number <- 20

# PC数量梯度。脚本会逐一使用这些PC数量完成邻居图、测试聚类和UMAP并绘图。
# 超过实际PCA/Harmony维度的数值会自动删除，co1/co2推荐值会自动加入梯度。
pc_number_list <- seq(10, 50, by = 5)

# 比较不同PC数量时暂时固定使用同一个resolution，避免同时改变两个参数。
pc_test_resolution <- 0.5

# 正式PC数量。默认NA表示绘图后由用户选择：
#   在RStudio交互运行时，控制台会要求输入候选PC数；
#   使用Rscript非交互运行时会在生成PC比较图和中间对象后停止，
#   用户查看图片并在这里填入候选值，重新运行后才会继续比较resolution。
final_pc_number <- NA_integer_

# Harmony使用哪一列作为批次信息。
# 默认每个10X捕获样本是一个批次，因此使用sample_id。
harmony_group <- "sample_id"

# 聚类时先计算多个分辨率，绘图后再选择其中一个作为正式cluster。
resolution_list <- seq(0.1, 1.0, by = 0.1)
# 默认NA表示等待用户查看resolution比较图后输入；也可以预先填写0.5等候选值。
final_resolution <- NA_real_

# 是否启用“正式处理前最终参数覆盖区”。
# 默认FALSE，使用上面的初始值；改成TRUE后，下方第二处参数值将覆盖这里的设置。
use_parameter_override <- FALSE



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


# ------------------------------------------------------------------------------
# 正式处理前最终参数覆盖区
# ------------------------------------------------------------------------------

# 这一处放在原始计数矩阵读取完成、CreateSeuratObject和QC正式开始之前。
# 如果前面只是使用通用默认值，可以把use_parameter_override改成TRUE，
# 然后只修改下面这一组参数；这里的值会覆盖脚本开头设置的同名初始值。
if (use_parameter_override) {
  # 建立Seurat对象时的最低要求。
  create_min_cells <- 3
  create_min_features <- 100

  # QC阈值。
  min_features <- 200
  max_features <- 8000
  min_counts <- 500
  max_counts <- 100000
  max_percent_mito <- 20
  max_percent_hb <- 5
  min_percent_ribo <- 0

  # 标准化、高变基因和PCA参数。
  n_variable_features <- 3000
  n_pcs_calculate <- 50
  auto_choose_pc <- TRUE
  manual_pc_number <- 20
  pc_number_list <- seq(10, 50, by = 5)
  pc_test_resolution <- 0.5
  final_pc_number <- NA_integer_

  # Harmony和聚类参数。
  harmony_group <- "sample_id"
  resolution_list <- seq(0.1, 1.0, by = 0.1)
  final_resolution <- NA_real_
}

# ------------------------------------------------------------------------------
# 3.4.1 读取用户提交的样本级metadata
# ------------------------------------------------------------------------------

# metadata采用“一行代表一个样本”的格式，至少包含sample_id列。例如：
#
# sample_id,group,batch,patient,tissue
# GSM123_sampleA,Tumor,Batch1,Patient1,Tumor
# GSM456_sampleB,Normal,Batch1,Patient1,Normal
#
# sample_id必须与前面得到的sample_name完全一致，包括大小写、下划线和连接符。
# 后面的group、batch、patient、tissue只是示例；用户可以增加任意样本信息列。

# 第一次运行且metadata文件不存在时，根据实际识别到的样本名生成模板。
if (!file.exists(metadata_file)) {
  sample_info_template <- data.frame(
    sample_id = sample_name,
    group = rep("Group_to_fill", length(sample_name)),
    batch = sample_name,
    patient = sample_name,
    tissue = rep("Tissue_to_fill", length(sample_name))
  )

  write.csv(
    sample_info_template,
    metadata_file,
    row.names = FALSE
  )

  message("已经自动生成可填写的样本metadata模板：", metadata_file)

  if (interactive()) {
    # 在RStudio/R GUI中运行时，脚本暂停在这里。
    # 用户打开CSV，填写并保存后回到R控制台按回车，脚本会在本次运行中继续读取。
    readline(paste0(
      "请填写并保存“", metadata_file,
      "”，完成后回到R控制台按回车继续："
    ))
  } else {
    # Rscript或服务器批处理模式无法等待用户现场填写，因此生成模板后安全停止。
    stop(
      "样本metadata模板已经生成：", metadata_file,
      "。当前为非交互式运行，请填写该CSV文件后重新运行脚本。"
    )
  }
}


# 读取用户填写完成的metadata。
# 如果模板原本已经存在，也会直接从这里读取，不会再次覆盖用户填写内容。
sample_info <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"sample_id" %in% colnames(sample_info)) {
  stop("用户metadata必须包含名为sample_id的列：", metadata_file)
}

sample_info$sample_id <- as.character(sample_info$sample_id)

if (anyDuplicated(sample_info$sample_id)) {
  stop("用户metadata中的sample_id存在重复值，每个样本只能保留一行。")
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

missing_metadata <- setdiff(sample_name, sample_info$sample_id)
if (length(missing_metadata) > 0) {
  stop(
    "以下已读取样本没有在用户metadata中找到：",
    paste(missing_metadata, collapse = "、"),
    "。请检查sample_id是否与样本名称完全一致。"
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
write.csv(
  sample_info,
  file.path(result_dir, "本次实际读入的样本metadata.csv"),
  row.names = FALSE
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



####四、逐样本计算QC指标和过滤低质量细胞####

# 推荐逐样本QC，而不是先合并再使用同一个阈值。
# 如果不同样本质量差异很大，可以把下面固定阈值改成每个样本独立阈值。

qc_summary <- data.frame()
mito_summary <- data.frame()

for (i in seq_along(scRNAlist)) {
  current_sample <- names(scRNAlist)[i]

  # 分别统计人类常见的“MT-”和小鼠常见的“mt-”线粒体基因命名。
  # 每个样本独立判断，哪一种匹配到的基因更多，就使用哪一种计算percent_mito。
  mito_genes_upper <- grep("^MT-", rownames(scRNAlist[[i]]), value = TRUE)
  mito_genes_lower <- grep("^mt-", rownames(scRNAlist[[i]]), value = TRUE)

  if (length(mito_genes_upper) >= length(mito_genes_lower) &&
      length(mito_genes_upper) > 0) {
    mito_genes <- mito_genes_upper
    mito_gene_style <- "MT-（人类常见命名）"
  } else if (length(mito_genes_lower) > 0) {
    mito_genes <- mito_genes_lower
    mito_gene_style <- "mt-（小鼠常见命名）"
  } else {
    mito_genes <- character()
    mito_gene_style <- "未识别到MT-/mt-线粒体基因"
  }

  if (length(mito_genes) > 0) {
    scRNAlist[[i]] <- PercentageFeatureSet(
      scRNAlist[[i]],
      features = mito_genes,
      col.name = "percent_mito"
    )
  } else {
    # 不再把percent_mito强行写成0。
    # 如果gene.column读错，或表达矩阵使用Ensembl ID，把它写成0会让线粒体QC失效，
    # 而且后续结果表面上仍能继续运行，因此这里直接停止并展示实际基因行名。
    stop(
      "样本 ", current_sample,
      " 没有识别到MT-/mt-线粒体基因，已停止QC，避免把percent_mito错误设为0。\n",
      "当前对象前20个基因名：",
      paste(head(rownames(scRNAlist[[i]]), 20), collapse = "、"),
      "\n请先查看结果目录中的“10X基因名称列自动判断结果.csv”。",
      "如果这些行名是Gene Expression，说明features文件的基因列选择仍不正确；",
      "如果是ENSG开头的Ensembl ID，则需要先把Ensembl ID转换为gene symbol。"
    )
  }

  # 核糖体和血红蛋白基因比例暂时沿用人类基因符号写法。
  scRNAlist[[i]] <- PercentageFeatureSet(
    scRNAlist[[i]], pattern = "^RP[SL]", col.name = "percent_ribo"
  )
  scRNAlist[[i]] <- PercentageFeatureSet(
    scRNAlist[[i]], pattern = "^HB[ABDEGMQZ]", col.name = "percent_hb"
  )

  # 在过滤细胞之前汇总每个样本的线粒体比例分布。
  # percent_cells_over_threshold表示有多少比例的细胞超过前面设置的max_percent_mito。
  percent_cells_over_threshold <- mean(
    scRNAlist[[i]]$percent_mito > max_percent_mito
  ) * 100

  if (median(scRNAlist[[i]]$percent_mito) > max_percent_mito) {
    mito_judgement <- "样本中位线粒体比例高于QC阈值，建议重点检查"
  } else if (percent_cells_over_threshold > 0) {
    mito_judgement <- "样本中位数未超阈值，但部分细胞超过QC阈值"
  } else {
    mito_judgement <- "当前样本没有细胞超过线粒体QC阈值"
  }

  mito_summary <- rbind(
    mito_summary,
    data.frame(
      sample_id = current_sample,
      mitochondrial_gene_style = mito_gene_style,
      mitochondrial_gene_number = length(mito_genes),
      mean_percent_mito = mean(scRNAlist[[i]]$percent_mito),
      median_percent_mito = median(scRNAlist[[i]]$percent_mito),
      Q1_percent_mito = as.numeric(quantile(scRNAlist[[i]]$percent_mito, 0.25)),
      Q3_percent_mito = as.numeric(quantile(scRNAlist[[i]]$percent_mito, 0.75)),
      max_percent_mito = max(scRNAlist[[i]]$percent_mito),
      QC_mito_threshold = max_percent_mito,
      percent_cells_over_threshold = percent_cells_over_threshold,
      judgement = mito_judgement
    )
  )

  # 先画过滤前QC图，再根据分布调整前面设置的阈值。
  p.qc.before <- VlnPlot(
    scRNAlist[[i]],
    features = c(
      "nFeature_RNA", "nCount_RNA",
      "percent_mito", "percent_ribo", "percent_hb"
    ),
    # 此时尚未运行NormalizeData()，RNA assay只有原始counts层。
    # 明确指定layer="counts"，可避免Seurat反复寻找不存在的data层并产生警告。
    assay = "RNA",
    layer = "counts",
    pt.size = 0,
    ncol = 5
  ) + plot_annotation(title = paste0(current_sample, "_QC过滤前"))

  ggsave(
    filename = file.path(result_dir, paste0(current_sample, "_QC过滤前.pdf")),
    plot = p.qc.before,
    width = 16,
    height = 4
  )

  cells_before <- ncol(scRNAlist[[i]])

  # subset中的条件必须同时满足。
  scRNAlist[[i]] <- subset(
    scRNAlist[[i]],
    subset = nFeature_RNA >= min_features &
      nFeature_RNA <= max_features &
      nCount_RNA >= min_counts &
      nCount_RNA <= max_counts &
      percent_mito <= max_percent_mito &
      percent_hb <= max_percent_hb &
      percent_ribo >= min_percent_ribo
  )

  cells_after <- ncol(scRNAlist[[i]])

  p.qc.after <- VlnPlot(
    scRNAlist[[i]],
    features = c(
      "nFeature_RNA", "nCount_RNA",
      "percent_mito", "percent_ribo", "percent_hb"
    ),
    assay = "RNA",
    layer = "counts",
    pt.size = 0,
    ncol = 5
  ) + plot_annotation(title = paste0(current_sample, "_QC过滤后"))

  ggsave(
    filename = file.path(result_dir, paste0(current_sample, "_QC过滤后.pdf")),
    plot = p.qc.after,
    width = 16,
    height = 4
  )

  qc_summary <- rbind(
    qc_summary,
    data.frame(
      sample_id = current_sample,
      cells_before_QC = cells_before,
      cells_after_QC = cells_after,
      median_nFeature = median(scRNAlist[[i]]$nFeature_RNA),
      median_nCount = median(scRNAlist[[i]]$nCount_RNA),
      median_percent_mito = median(scRNAlist[[i]]$percent_mito)
    )
  )
}

mito_summary
write.csv(
  mito_summary,
  file.path(result_dir, "各样本线粒体基因与比例判断_QC前.csv"),
  row.names = FALSE
)

qc_summary
write.csv(qc_summary, file.path(result_dir, "QC汇总表.csv"), row.names = FALSE)
saveRDS(scRNAlist, file.path(result_dir, "1.逐样本QC后对象列表.rds"))



####五、合并样本和检查metadata####

# 使用merge合并多个Seurat对象。
# Seurat V5会把不同样本保留为不同RNA layers，后面可以直接进行标准化和Harmony。
if (length(scRNAlist) == 1) {
  sce <- scRNAlist[[1]]
} else {
  sce <- merge(x = scRNAlist[[1]], y = scRNAlist[-1], project = "scRNA_project")
}

# 用户填写的每一列已经复制到每个细胞的sce@meta.data中。
# 同时在misc中保存原始的“一行一个样本”表，便于以后直接查看样本设计信息。
sce@misc$sample_metadata <- sample_info

rm(scRNAlist)

dim(sce)
table(sce$sample_id)

# 显示用户metadata中哪些列已经进入合并后的sce。
user_metadata_columns <- colnames(sample_info)
user_metadata_columns
head(sce@meta.data[, user_metadata_columns, drop = FALSE])

# group是常用但不是强制列；用户提供时再显示各组细胞数量。
if ("group" %in% colnames(sce@meta.data)) {
  table(sce$group)
}

# nCount与nFeature相关性图可辅助观察异常高UMI细胞。
p.scatter <- FeatureScatter(
  sce,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA",
  group.by = "sample_id"
)
ggsave(
  filename = file.path(result_dir, "QC_nCount与nFeature相关性.pdf"),
  plot = p.scatter,
  width = 8,
  height = 6
)

saveRDS(sce, file.path(result_dir, "2.合并与QC后_sce.rds"))



####六、标准化、高变基因、细胞周期评分和PCA####

DefaultAssay(sce) <- "RNA"

# merge后的Seurat V5对象通常为每个样本保留独立counts layer。
# 当前脚本使用RunHarmony直接校正PCA，不使用IntegrateLayers，因此这里先把RNA layers
# 合并为统一的counts层。这样CellCycleScoring、AddModuleScore和后续marker分析都能
# 使用明确的RNA data层，避免GetAssayData面对多个data layers时产生歧义。
sce[["RNA"]] <- JoinLayers(sce[["RNA"]])

# NormalizeData把原始UMI counts进行文库大小校正和log1p转换。
# JoinLayers后在统一的RNA counts层上进行标准化；批次信息仍保留在metadata中，
# 后面由Harmony根据harmony_group完成低维空间校正。
sce <- NormalizeData(
  sce,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


####细胞周期评分####

# 细胞周期评分不是必须删除细胞周期效应。
# 这里先计算S.Score、G2M.Score和Phase，供后续判断是否需要回归。
s.genes <- Seurat::cc.genes.updated.2019$s.genes
g2m.genes <- Seurat::cc.genes.updated.2019$g2m.genes
s.genes <- intersect(s.genes, rownames(sce))
g2m.genes <- intersect(g2m.genes, rownames(sce))

sce <- CellCycleScoring(
  object = sce,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE
)
sce$CC.Difference <- sce$S.Score - sce$G2M.Score

p.cycle <- VlnPlot(
  sce,
  features = c("S.Score", "G2M.Score"),
  group.by = "sample_id",
  ncol = 2,
  pt.size = 0
)
ggsave(
  filename = file.path(result_dir, "细胞周期评分.pdf"),
  plot = p.cycle,
  width = 10,
  height = 5
)


####高变基因、ScaleData和PCA####

sce <- FindVariableFeatures(
  sce,
  selection.method = "vst",
  nfeatures = n_variable_features
)

head(VariableFeatures(sce), 20)

p.variable <- VariableFeaturePlot(sce)
ggsave(
  filename = file.path(result_dir, "高变基因.pdf"),
  plot = p.variable,
  width = 8,
  height = 6
)

# 按照原脚本风格，对高变基因进行ScaleData。
# 这里不默认回归percent_mito或细胞周期，避免在未检查数据前过度校正生物学信号。
# 如果后续确认线粒体比例或细胞周期明显主导PCA，再自行添加vars.to.regress。
sce <- ScaleData(sce, features = VariableFeatures(sce))

sce <- RunPCA(
  sce,
  features = VariableFeatures(sce),
  npcs = n_pcs_calculate,
  verbose = FALSE
)

p.pca <- DimPlot(sce, reduction = "pca", group.by = "sample_id")
ggsave(
  filename = file.path(result_dir, "PCA_按样本.pdf"),
  plot = p.pca,
  width = 8,
  height = 6
)

p.elbow <- ElbowPlot(sce, ndims = n_pcs_calculate)
ggsave(
  filename = file.path(result_dir, "PCA_ElbowPlot.pdf"),
  plot = p.elbow,
  width = 7,
  height = 5
)


######按照co1/co2规则自动选择聚类使用的PC数######

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

if (auto_choose_pc) {
  recommended_pc_number <- choose_pc(sce)
} else {
  # 关闭co1/co2推荐时，使用manual_pc_number作为PC梯度的参考值。
  recommended_pc_number <- min(
    manual_pc_number,
    length(sce[["pca"]]@stdev)
  )
  message("已关闭co1/co2推荐，参考PC数为：", recommended_pc_number)
}

# 整理PC梯度：删除小于1或超过实际PCA维数的值，并自动加入co1/co2推荐值。
# 如果用户已经预先填写final_pc_number，也把它加入梯度并绘图，便于核对。
actual_pca_number <- length(sce[["pca"]]@stdev)
if (!is.na(final_pc_number) &&
    (final_pc_number < 1 || final_pc_number > actual_pca_number)) {
  stop("final_pc_number必须位于1到", actual_pca_number, "之间。")
}

pc_number_candidates <- c(pc_number_list, recommended_pc_number)
if (!is.na(final_pc_number)) {
  pc_number_candidates <- c(pc_number_candidates, final_pc_number)
}
pc_number_candidates <- sort(unique(as.integer(pc_number_candidates)))
pc_number_candidates <- pc_number_candidates[
  pc_number_candidates >= 1 & pc_number_candidates <= actual_pca_number
]

if (length(pc_number_candidates) == 0) {
  stop("PC梯度中没有可用数值，请检查pc_number_list和n_pcs_calculate。")
}

message(
  "本次将比较以下PC数量：",
  paste(pc_number_candidates, collapse = "、")
)

# 保存每个PC的贡献比例、累计比例、是否为推荐位置以及是否进入比较梯度。
pc_pct <- sce[["pca"]]@stdev / sum(sce[["pca"]]@stdev) * 100
pc_selection_summary <- data.frame(
  PC = seq_along(pc_pct),
  percent = pc_pct,
  cumulative_percent = cumsum(pc_pct),
  recommended_by_co1_co2 = seq_along(pc_pct) == recommended_pc_number,
  included_in_PC_gradient = seq_along(pc_pct) %in% pc_number_candidates
)
write.csv(
  pc_selection_summary,
  file.path(result_dir, "PC自动选择结果.csv"),
  row.names = FALSE
)



####七、Harmony去除样本或批次效应####

# Harmony根据PCA低维空间校正批次，不会直接改写RNA counts和data。
# harmony_group指定的metadata列必须存在，并且每个细胞都有对应值。
# 为了能够比较完整PC梯度，Harmony先使用候选梯度中的最大PC数；
# 后面的每次PC测试和正式聚类再分别使用1:当前候选PC数。
harmony_pc_num <- seq_len(max(pc_number_candidates))
sce <- RunHarmony(
  object = sce,
  group.by.vars = harmony_group,
  reduction.use = "pca",
  dims.use = harmony_pc_num,
  reduction.save = "harmony",
  verbose = TRUE
)

p.harmony <- DimPlot(sce, reduction = "harmony", group.by = harmony_group)
ggsave(
  filename = file.path(result_dir, "Harmony_按批次.pdf"),
  plot = p.harmony,
  width = 8,
  height = 6
)

saveRDS(sce, file.path(result_dir, "3.标准化_PCA_Harmony后_sce.rds"))



####八、选择PC、聚类分辨率和UMAP/tSNE####

######第一步：按照PC数量梯度分别聚类和绘图######

# PC比较阶段只改变PC数量，其他参数保持一致：
#   reduction固定为harmony；
#   resolution固定为pc_test_resolution；
#   UMAP的n.neighbors、min.dist和随机种子保持一致。
# 这样不同图片之间的主要区别来自纳入的PC数量，而不是resolution同时发生变化。

pc_gradient_plots <- list()
pc_gradient_summary <- data.frame()

for (i in seq_along(pc_number_candidates)) {
  current_pc_number <- pc_number_candidates[i]
  current_pc_dims <- seq_len(current_pc_number)

  message(
    "正在进行PC梯度比较：使用前", current_pc_number,
    "个Harmony维度，测试resolution = ", pc_test_resolution
  )

  # 每次使用临时对象完成邻居图、测试聚类和UMAP。
  # 这样不会在正式sce中同时保存大量测试图结构和UMAP结果，减少对象体积。
  sce_pc_test <- sce
  sce_pc_test <- FindNeighbors(
    sce_pc_test,
    reduction = "harmony",
    dims = current_pc_dims,
    verbose = FALSE
  )
  sce_pc_test <- FindClusters(
    sce_pc_test,
    resolution = pc_test_resolution,
    algorithm = 1,
    random.seed = 1234,
    verbose = FALSE
  )
  sce_pc_test <- RunUMAP(
    sce_pc_test,
    reduction = "harmony",
    dims = current_pc_dims,
    n.neighbors = 30,
    min.dist = 0.3,
    seed.use = 1234,
    verbose = FALSE
  )

  current_cluster_number <- length(unique(sce_pc_test$seurat_clusters))
  current_title <- paste0(
    "PC = 1:", current_pc_number,
    " | resolution = ", pc_test_resolution,
    " | clusters = ", current_cluster_number
  )
  if (current_pc_number == recommended_pc_number) {
    current_title <- paste0(current_title, " | co1/co2推荐")
  }

  pc_gradient_plots[[i]] <- DimPlot(
    sce_pc_test,
    reduction = "umap",
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE,
    shuffle = TRUE
  ) + ggtitle(current_title) + NoLegend()

  # 每个候选PC单独保存一张图，方便放大查看cluster边界和小群。
  ggsave(
    filename = file.path(
      result_dir,
      paste0("PC梯度_PC1-", current_pc_number, "_UMAP.pdf")
    ),
    plot = pc_gradient_plots[[i]],
    width = 8,
    height = 7
  )

  pc_gradient_summary <- rbind(
    pc_gradient_summary,
    data.frame(
      PC_number = current_pc_number,
      test_resolution = pc_test_resolution,
      cluster_number = current_cluster_number,
      recommended_by_co1_co2 = current_pc_number == recommended_pc_number
    )
  )

  rm(sce_pc_test)
  invisible(gc())
}

# 把所有PC梯度图合并到同一个PDF中，便于横向比较。
p.pc.gradient <- wrap_plots(pc_gradient_plots, ncol = 3)
ggsave(
  filename = file.path(result_dir, "PC梯度_全部UMAP比较.pdf"),
  plot = p.pc.gradient,
  width = 15,
  height = 5 * ceiling(length(pc_gradient_plots) / 3)
)

write.csv(
  pc_gradient_summary,
  file.path(result_dir, "PC梯度_聚类数量汇总.csv"),
  row.names = FALSE
)

# 保存绘图完成、尚未正式选择PC的中间对象。
# 该对象已经包含PCA和最大候选维度的Harmony结果，但不包含测试聚类产生的大量临时图。
saveRDS(
  sce,
  file.path(result_dir, "3.1_PC梯度比较完成_待选择PC_sce.rds")
)


######用户确定正式PC数量######

if (is.na(final_pc_number)) {
  if (interactive()) {
    pc_answer <- trimws(readline(paste0(
      "请查看result目录中的PC梯度图，并输入正式PC数量（候选值：",
      paste(pc_number_candidates, collapse = "、"),
      "；直接回车使用co1/co2推荐值", recommended_pc_number, "）："
    )))

    if (pc_answer == "") {
      final_pc_number <- recommended_pc_number
    } else {
      final_pc_number <- suppressWarnings(as.integer(pc_answer))
    }
  } else {
    stop(
      "PC梯度图已经全部生成。当前为非交互式运行，无法等待用户输入。",
      "请查看PC梯度_全部UMAP比较.pdf，将final_pc_number填写为以下候选值之一后重新运行：",
      paste(pc_number_candidates, collapse = "、")
    )
  }
}

if (is.na(final_pc_number) || !final_pc_number %in% pc_number_candidates) {
  stop(
    "final_pc_number必须是本次已经绘图的候选值之一：",
    paste(pc_number_candidates, collapse = "、")
  )
}

pc.num <- seq_len(final_pc_number)
message("正式聚类使用PC：1:", final_pc_number)


######第二步：固定正式PC数量后比较resolution######

# FindNeighbors只需要运行一次；不同resolution共用同一个邻居图。
# 一次输入多个resolution后，metadata会生成RNA_snn_res.0.1、0.2……等列。
sce <- FindNeighbors(
  sce,
  reduction = "harmony",
  dims = pc.num
)

sce <- FindClusters(
  sce,
  resolution = resolution_list,
  algorithm = 1,
  random.seed = 1234
)

# UMAP只依赖选择的reduction和PC，不依赖cluster分辨率，因此只需要计算一次。
sce <- RunUMAP(
  sce,
  reduction = "harmony",
  dims = pc.num,
  n.neighbors = 30,
  min.dist = 0.3,
  seed.use = 1234
)

sce <- RunTSNE(
  sce,
  reduction = "harmony",
  dims = pc.num,
  seed.use = 1234
)


####比较不同聚类分辨率####

# 如果安装了clustree，可以直接观察不同resolution之间cluster的分裂关系。
if (requireNamespace("clustree", quietly = TRUE)) {
  p.clustree <- clustree::clustree(sce@meta.data, prefix = "RNA_snn_res.")
  ggsave(
    filename = file.path(result_dir, "不同分辨率_clustree.pdf"),
    plot = p.clustree,
    width = 12,
    height = 9
  )
}

# 分别画出不同resolution的UMAP，便于人工选择正式分辨率。
resolution_plots <- list()
for (i in seq_along(resolution_list)) {
  resolution_column <- paste0("RNA_snn_res.", resolution_list[i])
  resolution_plots[[i]] <- DimPlot(
    sce,
    reduction = "umap",
    group.by = resolution_column,
    label = TRUE,
    repel = TRUE,
    shuffle = TRUE
  ) + ggtitle(paste0("resolution = ", resolution_list[i])) + NoLegend()
}

p.resolution <- wrap_plots(resolution_plots, ncol = 3)
ggsave(
  filename = file.path(result_dir, "不同聚类分辨率_UMAP.pdf"),
  plot = p.resolution,
  width = 15,
  height = 5 * ceiling(length(resolution_plots) / 3)
)

# 保存已经固定正式PC并完成resolution梯度计算的中间对象。
saveRDS(
  sce,
  file.path(result_dir, "3.2_resolution比较完成_待选择resolution_sce.rds")
)


#####正式选择聚类分辨率####

# PC数量确定后，才让用户根据UMAP和clustree选择正式resolution。
# 交互运行时可以直接在控制台输入；非交互运行时会生成图片和中间对象后停止。
if (is.na(final_resolution)) {
  if (interactive()) {
    resolution_answer <- trimws(readline(paste0(
      "请查看不同聚类分辨率_UMAP.pdf和clustree，并输入正式resolution（候选值：",
      paste(resolution_list, collapse = "、"), "）："
    )))
    final_resolution <- suppressWarnings(as.numeric(resolution_answer))
  } else {
    stop(
      "resolution梯度图已经生成。当前为非交互式运行，无法等待用户输入。",
      "请查看不同聚类分辨率_UMAP.pdf和clustree，将final_resolution填写为以下候选值之一后重新运行：",
      paste(resolution_list, collapse = "、")
    )
  }
}

# 使用数值差值匹配，避免小数在计算机内部表示造成0.3无法精确匹配的问题。
resolution_match <- which(abs(resolution_list - final_resolution) < 1e-10)
if (is.na(final_resolution) || length(resolution_match) != 1) {
  stop(
    "final_resolution必须是resolution_list中的候选值之一：",
    paste(resolution_list, collapse = "、")
  )
}
final_resolution <- resolution_list[resolution_match]

final_cluster_column <- paste0("RNA_snn_res.", final_resolution)
sce$seurat_clusters <- as.character(sce@meta.data[[final_cluster_column]])
Idents(sce) <- "seurat_clusters"

table(sce$seurat_clusters)

write.csv(
  data.frame(
    final_PC_number = final_pc_number,
    final_PC_range = paste0("1:", final_pc_number),
    final_resolution = final_resolution,
    clustering_algorithm = "Louvain_algorithm_1"
  ),
  file.path(result_dir, "正式聚类参数.csv"),
  row.names = FALSE
)

p.cluster <- DimPlot(
  sce,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  shuffle = TRUE,
  raster = FALSE
) + NoLegend()

p.sample <- DimPlot(
  sce,
  reduction = "umap",
  group.by = "sample_id",
  raster = FALSE
)

p.cluster.sample <- p.cluster + p.sample
ggsave(
  filename = file.path(result_dir, "正式聚类_UMAP.pdf"),
  plot = p.cluster.sample,
  width = 15,
  height = 7
)

p.tsne <- DimPlot(
  sce,
  reduction = "tsne",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) + NoLegend()
ggsave(
  filename = file.path(result_dir, "正式聚类_tSNE.pdf"),
  plot = p.tsne,
  width = 8,
  height = 7
)

saveRDS(sce, file.path(result_dir, "4.正式聚类后_sce.rds"))



####九、寻找每个cluster的marker基因####

# Seurat V5在FindAllMarkers前需要把各样本RNA layers重新合并。
# JoinLayers只合并表达矩阵层，不会删除Harmony降维结果和UMAP坐标。
sce[["RNA"]] <- JoinLayers(sce[["RNA"]])
DefaultAssay(sce) <- "RNA"
Idents(sce) <- "seurat_clusters"

markers <- FindAllMarkers(
  object = sce,
  assay = "RNA",
  test.use = "wilcox",
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)

write.csv(
  markers,
  file.path(result_dir, "所有cluster_marker.csv"),
  row.names = FALSE
)

# 每个cluster按照avg_log2FC选取前10个marker，用于快速查看cluster特征。
top10_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup()

write.csv(
  top10_markers,
  file.path(result_dir, "每个cluster_top10_marker.csv"),
  row.names = FALSE
)

# 热图使用RNA assay的scale.data，因此先对需要展示的marker进行ScaleData。
heatmap_genes <- unique(top10_markers$gene)
sce <- ScaleData(sce, assay = "RNA", features = heatmap_genes, verbose = FALSE)

p.heatmap <- DoHeatmap(
  sce,
  features = heatmap_genes,
  group.by = "seurat_clusters",
  raster = TRUE
) + NoLegend()

ggsave(
  filename = file.path(result_dir, "cluster_top10_marker热图.pdf"),
  plot = p.heatmap,
  width = 12,
  height = max(8, length(heatmap_genes) * 0.12)
)



####十、根据经典marker完成大亚群注释####

#####1.使用经典marker检查每个cluster####

# 以下marker只用于大亚群级别判断，不直接决定最终注释。
# 正式注释时需要同时结合FindAllMarkers结果、DotPlot、FeaturePlot和生物学背景。
genes_to_check <- c(
  # Epithelial
  "EPCAM", "KRT8", "KRT18", "KRT19",

  # Endothelial
  "PECAM1", "VWF", "CLDN5", "EMCN",

  # Fibroblast / stromal
  "COL1A1", "COL1A2", "DCN", "COL3A1",

  # T cells
  "CD3D", "CD3E", "TRAC", "CD4", "CD8A", "CD8B",

  # NK cells
  "NKG7", "GNLY", "KLRD1", "NCR1",

  # B cells
  "CD79A", "CD19", "MS4A1", "CD37",

  # Plasma cells
  "MZB1", "JCHAIN", "SDC1", "IGHG1",

  # Myeloid cells
  "LST1", "TYROBP", "FCER1G", "CTSS", "C1QA", "S100A8", "S100A9",

  # Dendritic cells
  "CD1C", "CLEC9A", "LILRA4", "LAMP3",

  # Mast cells
  "TPSAB1", "TPSB2", "CPA3", "MS4A2"
)

# 删除数据中不存在的基因，避免DotPlot因为个别基因缺失而警告。
genes_to_check <- intersect(genes_to_check, rownames(sce))

p.marker.dot <- DotPlot(
  sce,
  features = genes_to_check,
  assay = "RNA",
  group.by = "seurat_clusters"
) +
  coord_flip() +
  theme_bw() +
  ggtitle("大亚群经典marker")

ggsave(
  filename = file.path(result_dir, "大亚群经典marker_DotPlot.pdf"),
  plot = p.marker.dot,
  width = 12,
  height = max(8, length(genes_to_check) * 0.22)
)


#####2.按照原脚本方法直接进行大亚群注释####

# 先查看cluster编号和每个cluster的细胞数量。
# 然后结合上面的DotPlot、top10 marker热图和所有cluster_marker.csv，
# 把相应的cluster编号直接填入下面的c()中。
table(sce@meta.data$seurat_clusters)

# 例如某一类细胞对应cluster 0、3、7，可以写成：c(0, 3, 7)。
# 当前全部为c()，表示尚未指定cluster，不会把任何cluster错误分配给某类细胞。
cluster_T_NK <- c()
cluster_Plasma <- c()
cluster_B <- c()
cluster_Myeloid <- c()
cluster_Fibroblast <- c()
cluster_Endothelial <- c()
cluster_Mast <- c()
cluster_Epi <- c()

# 建立cluster与celltype的对应表；所有cluster初始均为unknown。
# 下面的写法与原脚本一致：根据ClusterID是否属于指定cluster向量，直接修改第2列。
celltype=data.frame(
  ClusterID=sort(unique(as.character(sce$seurat_clusters))),
  celltype='unknown'
)

celltype[celltype$ClusterID %in% cluster_T_NK,2]='T/NK'
celltype[celltype$ClusterID %in% cluster_Plasma,2]='Plasma'
celltype[celltype$ClusterID %in% cluster_B,2]='B'
celltype[celltype$ClusterID %in% cluster_Myeloid,2]='Myeloid'
celltype[celltype$ClusterID %in% cluster_Fibroblast,2]='Fibroblast'
celltype[celltype$ClusterID %in% cluster_Endothelial,2]='Endothelial'
celltype[celltype$ClusterID %in% cluster_Mast,2]='Mast'
celltype[celltype$ClusterID %in% cluster_Epi,2]='Epi'

celltype
table(celltype$celltype)


#####3.将大亚群注释写入Seurat对象####

# 先在metadata中增加celltype列，所有细胞初始标记为unknown。
sce@meta.data$celltype = "unknown"

# 按照原脚本方法逐行读取celltype对应表，根据seurat_clusters给每个细胞赋值。
for(i in 1:nrow(celltype)){
  sce@meta.data[
    which(sce@meta.data$seurat_clusters == celltype$ClusterID[i]),
    'celltype'
  ] <- celltype$celltype[i]
}

# 保存本次代码中实际使用的cluster与大亚群对应关系，方便复核；
# 该文件只是结果记录，不参与注释，也不需要人工填表后再读入。
write.csv(
  celltype,
  file.path(result_dir, "大亚群注释对应关系.csv"),
  row.names = FALSE
)

table(sce$celltype)
table(sce$celltype, sce$seurat_clusters)
Idents(sce) <- "celltype"


#####4.绘制大亚群注释结果####

p.celltype <- DimPlot(
  sce,
  reduction = "umap",
  group.by = "celltype",
  label = TRUE,
  repel = TRUE,
  shuffle = TRUE,
  raster = FALSE
) + NoLegend()

p.celltype.sample <- DimPlot(
  sce,
  reduction = "umap",
  group.by = "celltype",
  split.by = "sample_id",
  label = FALSE,
  raster = FALSE
)

ggsave(
  filename = file.path(result_dir, "大亚群注释_UMAP.pdf"),
  plot = p.celltype,
  width = 10,
  height = 8
)

ggsave(
  filename = file.path(result_dir, "各样本大亚群注释_UMAP.pdf"),
  plot = p.celltype.sample,
  width = max(12, 4 * length(unique(sce$sample_id))),
  height = 6
)

p.marker.celltype <- DotPlot(
  sce,
  features = genes_to_check,
  assay = "RNA",
  group.by = "celltype"
) +
  coord_flip() +
  theme_bw() +
  ggtitle("大亚群注释后marker复核")

ggsave(
  filename = file.path(result_dir, "大亚群注释后_marker复核.pdf"),
  plot = p.marker.celltype,
  width = 12,
  height = max(8, length(genes_to_check) * 0.22)
)


#####5.保存大亚群注释后的对象####

saveRDS(sce, file.path(result_dir, "5.大亚群注释完成_sce.rds"))
write.csv(sce@meta.data, file.path(result_dir, "大亚群注释后_metadata.csv"))
writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))

cat(
  "\n流程已运行到大亚群注释。\n",
  "Seurat对象：", file.path(result_dir, "5.大亚群注释完成_sce.rds"), "\n",
  "细胞数量：", ncol(sce), "\n",
  "基因数量：", nrow(sce), "\n",
  "正式PC范围：1:", final_pc_number, "\n",
  "正式resolution：", final_resolution, "\n",
  sep = ""
)


# ==============================================================================
# 脚本到此结束：暂不进行细胞亚群重聚类和其他下游分析
# ==============================================================================
