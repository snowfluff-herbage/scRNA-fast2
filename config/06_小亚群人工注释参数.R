# ==============================================================================
# 入口3：小亚群人工注释参数
# 由“运行_03_入口_小亚群人工注释与汇总.R”读取。
# ==============================================================================

# 默认NULL：读取result_dir/小亚群自动分析。
# 如需注释另一套入口2结果，可填写其“小亚群自动分析”目录。
subcluster_annotation_input_root <- NULL

# 默认NULL：注释总表保存在“小亚群自动分析”目录。
# 也可以指定其他.xlsx路径。
subcluster_annotation_workbook <- NULL

# TRUE：根据当前小亚群cluster重新生成注释总表，会覆盖旧表。
# FALSE：保留已经填写的注释表。若cluster结构变化，程序会提示后停止。
subcluster_annotation_rebuild_workbook <- FALSE

# TRUE：第一次生成总表时，将先验marker评分最高的亚群名称预填到
# “正式小亚群名称”列；这些只是建议，用户可以任意修改。
# FALSE：正式名称列保持空白，仅显示第一和第二候选名称。
subcluster_annotation_prefill_prior <- TRUE

# 正式名称留空时如何处理：
# TRUE：自动填写下面的“未注释”标签并继续；
# FALSE：只要存在空白正式名称就停止，要求补齐后再运行。
subcluster_annotation_allow_blank <- TRUE
subcluster_annotation_unknown_label <- "未注释"

# TRUE：在交互式RStudio中生成/发现注释表后暂停，等待用户填写、保存并关闭Excel。
# Rscript非交互运行时，新建总表后会停止；填写后再次运行入口3即可正式注释。
subcluster_annotation_wait_for_edit <- TRUE

# 最终文件统一保存到该目录。默认NULL时使用：
# result_dir/小亚群注释完成文件
subcluster_annotation_output_dir <- NULL

# 是否输出合并对象的完整metadata CSV。细胞很多时文件可能较大。
subcluster_annotation_save_full_metadata_csv <- FALSE

set.seed(1234)
