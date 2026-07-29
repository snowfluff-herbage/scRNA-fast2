# ==============================================================================
# 功能块00：载入整个流程需要的R包
# 由入口脚本自动source，一般不需要修改
# ==============================================================================

required_packages <- c(
  "Seurat", "SeuratObject", "Matrix", "dplyr",
  "ggplot2", "patchwork", "harmony", "data.table"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "缺少以下R包：", paste(missing_packages, collapse = "、"),
    "。请安装后重新运行主脚本。"
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(harmony)
  library(data.table)
})

