# ==============================================================================
# 下游小亚群分析：外部大亚群注释RDS输入参数
# 由入口2和入口4共同读取。
# ==============================================================================

# FALSE：使用流程默认生成的result/5.大亚群注释完成_sce.rds，
#        或config/04中的subcluster_input_file。
# TRUE ：直接使用下面指定的任意大亚群注释完成Seurat RDS文件。
downstream_use_external_annotated_rds <- TRUE

# 当downstream_use_external_annotated_rds为TRUE时必须填写。
# Windows路径请使用正斜杠，例如：
# downstream_external_annotated_rds <- "D:/singlecell/my_annotated_object.rds"
downstream_external_annotated_rds <- "C:/singlecell/project/result/5.大亚群注释完成_sce.rds"

# 外部Seurat对象中保存大亚群名称的metadata列。
# 如果不是celltype，例如为major_celltype，可直接改为相应列名。
downstream_external_celltype_column <- "celltype"

# 小亚群流程需要原始/标准化RNA assay。通常保持"RNA"。
# 如果外部对象的RNA assay使用其他名称，可在此指定；当对象中没有RNA assay时，
# 程序会把这里指定的assay复制为RNA供下游使用。
downstream_external_rna_assay <- "RNA"

# 是否允许大亚群名称为空。默认FALSE，发现NA或空名称时直接停止并提示。
downstream_allow_blank_celltype <- FALSE
