# ==============================================================================
# 小亚群自动重聚类参数
# 由“运行_02_入口_从大亚群注释至小亚群完成.R”读取；
# 修改本文件不会影响前面的大亚群流程。
# ==============================================================================

# 要分析的大亚群。
#   "all"：分析sce$celltype中除unknown以外的全部大亚群；
#   也可以填写名称向量，例如c("T/NK", "B cell", "Myeloid")。
#   名称也可只写兼容关键词B、Epi、Mye、T、Fib、Endo、Plasma，不区分大小写。
subcluster_target_celltypes <- "all"
subcluster_exclude_celltypes <- c("unknown", "Unknown", "未注释")

# 默认NULL时读取result_dir中的“5.大亚群注释完成_sce.rds”。
# 也可以指定另一个已经含有celltype列的Seurat对象RDS文件。
# 更推荐在config/07_外部大亚群RDS输入参数.R中集中设置外部对象、注释列和assay。
subcluster_input_file <- NULL

# 少于该细胞数的大亚群不再聚类，会在汇总表中记录为“细胞数不足”。
subcluster_min_cells <- 100

# 断点续跑方式：
#   "auto"    ：已成功完成的大亚群直接跳过；报错修复后从报错的大亚群继续；
#   "restart" ：全部大亚群重新计算。
# 当本配置文件或输入的大亚群对象发生变化时，auto会自动使旧检查点失效并重算。
subcluster_resume_mode <- "auto"


# ----------------------------------------------------------------------------
# 标准化、PCA和Harmony
# ----------------------------------------------------------------------------

subcluster_n_variable_features <- 3000
subcluster_n_pcs_calculate <- 50

# 这些基因保留在RNA表达矩阵中供后续查看和FindAllMarkers使用，
# 但从高变基因中排除，避免线粒体、核糖体或MALAT1主导小亚群聚类。
subcluster_exclude_variable_patterns <- c("^MT-", "^mt-", "^RP[SL]", "^Rp[sl]", "^MALAT1$")

# 是否在ScaleData时回归变量。默认不回归；需要时可填写：
# c("percent_mito", "nCount_RNA")
subcluster_vars_to_regress <- character()

# Harmony模式："always"、"auto"或"none"。
# always：只要指定批次列存在且包含至少两个非空水平就运行Harmony；
# auto  ：当前实现同样在存在至少两个批次时运行，并在结果中记录；
# none  ：不运行Harmony，后续直接使用PCA。
subcluster_harmony_mode <- "always"
subcluster_harmony_group <- "sample_id"

# 可按大亚群覆盖Harmony方式，适用于希望保留强样本异质性的亚群。
# 示例：subcluster_harmony_overrides <- list("Epithelial" = "none")
subcluster_harmony_overrides <- list()


# ----------------------------------------------------------------------------
# PC、resolution和聚类算法
# ----------------------------------------------------------------------------

# 每个大亚群分别计算co1和co2，正式PC自动取min(co1, co2)。
# 如果推荐PC低于subcluster_min_pc，则使用subcluster_min_pc；
# 如果超过实际PCA维数，则使用实际最大维数。运行时不再询问PC。
subcluster_min_pc <- 5

# 程序会为当前大亚群生成各候选resolution的UMAP和clustree，
# 然后在RStudio控制台逐个大亚群询问一次正式resolution。
subcluster_resolution_candidates <- c(0.2, 0.4, 0.6, 0.8, 1.0)

# 正式聚类算法固定，不再逐亚群提问：
# 1=Louvain；2=Louvain多层细化；3=SLM；4=Leiden。
subcluster_algorithm <- 1

# 默认留空，表示每个大亚群都现场询问resolution。
# 如需非交互批量运行，可按大亚群名称预先指定；未指定的亚群仍会询问。
# 示例：subcluster_resolution_overrides <- list("B cell" = 0.4, "Myeloid" = 0.6)
subcluster_resolution_overrides <- list()

# 自动推荐评分参数。仅用于比较计算可行性，不等于生物学真值；
# 最终仍应结合UMAP、marker和已知生物学复核。
subcluster_evaluation_max_cells <- 2000
subcluster_min_cluster_cells <- 20
subcluster_min_cluster_fraction <- 0.005
subcluster_score_weights <- c(
  silhouette = 0.50,
  stability = 0.35,
  cluster_size = 0.15
)


# ----------------------------------------------------------------------------
# 正式UMAP和marker参数
# ----------------------------------------------------------------------------

subcluster_umap_n_neighbors <- 30
subcluster_umap_min_dist <- 0.3

subcluster_marker_test <- "wilcox"
# 小亚群重聚类默认保留正、负marker，避免only.pos=TRUE遗漏有用信息。
subcluster_marker_only_pos <- FALSE
subcluster_marker_min_pct <- 0.10
subcluster_marker_logfc_threshold <- 0.25
subcluster_marker_padj_threshold <- 0.05
subcluster_top_marker_number <- 20
subcluster_heatmap_top_number <- 10

# 为避免超大对象写出巨大的metadata CSV，默认只保存小亚群汇总和映射表。
subcluster_save_full_metadata_csv <- FALSE

set.seed(1234)
