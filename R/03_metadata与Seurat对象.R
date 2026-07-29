# ==============================================================================
# 功能块03：生成/读取样本metadata，并建立逐样本Seurat对象
# 输入：count_list、sample_name
# 输出：sample_info、scRNAlist
# ==============================================================================

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




