# ==============================================================================
# 入口4：单个大亚群重新聚类与注释参数
# ==============================================================================

# 默认NULL：在RStudio控制台显示全部大亚群并询问。
# 可填写完整名称或兼容关键词，例如"B"、"Epi"、"Mye"、"T"、"Fib"、
# "Endo"、"Plasma"，不区分大小写。
single_recluster_target_celltype <- NULL

# 交互运行时，每轮都会询问PC和resolution；以下值只作为回车时的默认值。
# PC为NA时使用当前大亚群co1/co2推荐值。
single_recluster_default_pc <- NA_integer_
single_recluster_default_resolution <- 0.4

# 非交互Rscript运行时必须给出有效PC和resolution，并自动把该组作为最终方案。
single_recluster_noninteractive_pc <- NA_integer_
single_recluster_noninteractive_resolution <- NA_real_

# 聚类算法：1=Louvain；2=Louvain多层细化；3=SLM；4=Leiden。
single_recluster_algorithm <- 1

# TRUE：保存每次试算的Seurat对象；FALSE：只保存UMAP、参数和cluster数量表。
# 试算很多时建议FALSE，以免占用大量磁盘。
single_recluster_save_trial_object <- FALSE

# TRUE：每次试算完成后尝试自动打开UMAP PDF；无法自动打开时仍会打印完整路径。
single_recluster_open_trial_pdf <- TRUE

# TRUE：如果上次已经确认最终PC/res并生成“待人工注释”对象，则直接恢复人工注释，
# 不重新做PCA/Harmony和试算。想重新开始试算时改为FALSE。
single_recluster_resume_pending_annotation <- TRUE

# 默认NULL时输出到result_dir/单大亚群重新聚类与注释。
single_recluster_output_root <- NULL

set.seed(1234)
