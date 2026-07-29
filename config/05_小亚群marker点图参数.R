# ==============================================================================
# 大亚群内小亚群先验marker点图参数
# 由入口2自动读取；也可由“辅助_小亚群marker点图.R”单独读取并重画。
#
# marker分为两部分：
#   1. subcluster_marker_uploaded_panels：来自用户提供的6个Excel文件；
#   2. subcluster_marker_literature_panels：根据结直肠癌/肿瘤微环境单细胞文献补充。
# 最终DotPlot使用两部分的并集；不会根据先验marker自动强制给cluster命名。
# ==============================================================================


# ----要绘图的大亚群----
# "all"表示对入口2已经成功生成子对象、且存在marker面板的全部大亚群绘图。
# 也可以填写实际的大亚群名称，例如c("T/NK", "Myeloid", "B cell")。
subcluster_marker_target_celltypes <- "all"

# 默认NULL时读取result_dir/小亚群自动分析；也可指定另一个入口2结果目录。
subcluster_marker_input_root <- NULL

# 子对象中的cluster列。"auto"优先使用subcluster，其次使用seurat_clusters。
subcluster_marker_cluster_column <- "auto"

# RNA assay和DotPlot显示参数。
subcluster_marker_assay <- "RNA"
subcluster_marker_dot_scale <- 6
subcluster_marker_colors <- c("grey90", "#FEE08B", "#D73027")
subcluster_marker_max_panels_per_page <- 7
subcluster_marker_pdf_width <- 13
subcluster_marker_base_height <- 3
subcluster_marker_height_per_gene <- 0.23
subcluster_marker_save_page_png <- TRUE
subcluster_marker_png_dpi <- 220

# 是否同时输出每个大亚群的DotPlot底层数值表：
# cluster、小亚群先验名称、gene、表达细胞比例和平均标准化表达。
subcluster_marker_save_plot_data <- TRUE


# ----------------------------------------------------------------------------
# 大亚群名称匹配
# ----------------------------------------------------------------------------
# 左侧是marker面板中的标准名称，右侧是允许自动识别的实际sce$celltype名称。
# 如果你在config/02中使用了其他自由名称，可在相应向量中直接增加。
subcluster_marker_big_group_aliases <- list(
  "T/NK" = c("T/NK", "T_NK", "T cell and NK", "T细胞和NK细胞", "T/NK细胞"),
  "Plasma" = c("Plasma", "Plasma cell", "PC", "浆细胞"),
  "B cell" = c("B cell", "B", "B cells", "B细胞"),
  "Myeloid" = c("Myeloid", "Myeloid cell", "Mye", "髓系细胞"),
  "Fibroblast" = c("Fibroblast", "Fibroblasts", "CAF", "Stromal", "成纤维细胞", "基质细胞"),
  "Endothelial" = c("Endothelial", "Endothelial cell", "EC", "内皮细胞"),
  "Mast" = c("Mast", "Mast cell", "Mast cells", "肥大细胞"),
  "Epithelial" = c("Epithelial", "Epithelial cell", "Epi", "上皮细胞")
)

# 自由名称仍无法自动匹配时，在这里明确写“实际大亚群名称 = marker面板标准名称”。
# 示例：subcluster_marker_manual_group_map <- c("免疫T与NK" = "T/NK")
# 除别名外，名称中含B、Epi、Mye、T、Fib、Endo或Plasma时也会自动匹配，
# 且不区分大小写。较长关键词优先，避免Fib中的B等短关键词造成冲突。
subcluster_marker_manual_group_map <- character()


# ----------------------------------------------------------------------------
# 用户Excel中的原始先验marker
# ----------------------------------------------------------------------------

subcluster_marker_uploaded_panels <- list(
  "T/NK" = list(
    "CD4.c01.Tn(CCR7)" = c("CCR7", "SELL", "IL7R"),
    "CD4.c02.Tcm(ANXA1)" = c("ANXA1", "GPR183"),
    "CD4.c03.Th17(IL17A)" = c("IL17A", "FURIN", "RORA"),
    "CD4.c04.Tfh(CXCL13)" = c("CXCL13", "IL6ST", "MAF", "TOX"),
    "CD4.c05.Treg(FOXP3)" = c("IL2RA", "FOXP3", "BATF"),
    "CD4.c06.active(NR4A2)" = c("FOSB", "CD69", "JUN"),
    "CD4.c07.stress(HSPA1A)" = c("HSPA1A", "HSPA1B"),
    "CD4.c08.transiting" = c("RPL29", "RPL39"),
    "CD4.c09.Tex(LAG3)" = c("LAG3", "TNFRSF18"),
    "CD8.c10.Trm(XCL2)" = c("XCL1", "XCL2", "KLRD1", "ITGAE"),
    "CD8.c11.Tem(GZMK)" = c("GZMK", "CXCR4", "ENC1"),
    "CD8.c12.Teff(GZMB)" = c("GZMB", "GZMA", "GZMH", "NKG7"),
    "CD8.c13.Tex(HAVCR2)" = c("HAVCR2", "TIGIT", "LAYN"),
    "CD8.c14.MAIT(SLC4A10)" = c("KLRB1", "SLC4A10", "AQP3", "LTB"),
    "CD8.c15.stress(HSPA1A)" = c("HSPH1", "HSPD1"),
    "CD8.c16.active(CD69)" = c("JUND", "FOS", "MYADM"),
    "CD8.c17.transiting" = c("RPS4X", "RPS3A"),
    "Innate.Immune.Cell" = c("TRDC", "TRGC1", "NCAM1"),
    "T.c18.STMN1" = c("MKI67", "STMN1", "PCNA")
  ),

  "B cell" = list(
    "Bn.c01.TCL1A" = c("TCL1A", "IGHD", "IGHM", "FCER2"),
    "Bn.c02.NR4A2" = c("NR4A1", "NR4A2"),
    "Bm.c03.TXNIP" = c("CD27", "TNFRSF13B", "TXNIP"),
    "Bm.c04.GPR183" = c("GPR183", "CREM", "PTPRCAP"),
    "Bm.c05.HSPA1A" = c("HSPA1B", "HSPA1A", "DNAJB1"),
    "Bgc.c06.BCL6" = c("BCL6", "AICDA", "LMO2", "RGS13")
  ),

  "Plasma" = list(
    # 人类IgA恒定区只有IGHA1和IGHA2；原表IGHA3/IGHA4不进入正式marker面板。
    "PC.c01.IgA+" = c("IGHA1", "IGHA2"),
    "PC.c02.IgG+" = c("IGHG4", "IGHG3", "IGHG1", "IGHG2"),
    "PC.c03.early" = c("HLA-DQA1", "HLA-DPA1", "HLA-DRB1"),
    "PC.c04.cycling" = c("STMN1", "TYMS", "MKI67"),
    "PC.c05.IgM+" = c("IGHM"),
    # IGKV家族基因容易受BCR克隆型影响，解释时不能单独作为稳定细胞类型依据。
    "PC.c06.transisting" = c("IGKV1-9", "IGKV1-5", "IGKV1-6")
  ),

  "Myeloid" = list(
    "Macro.c01.FCN1" = c("FCN1", "VCAN", "S100A8"),
    "Macro.c02.C1QA" = c("C1QA", "C1QB", "APOE"),
    "Macro.c03.SPP1" = c("SPP1", "MARCO", "TREM2"),
    "Macro.c04.IL1A" = c("IL1A", "CCL20", "CXCL5"),
    "Macro.c05.HSPA1B" = c("HSPA1B", "HSPA1A", "JUN"),
    "Macro.c06.TYMS" = c("TK1", "TYMS", "MKI67"),
    "cDC.c07.CD1C" = c("CLEC10A", "CD1C", "HLA-DQA1"),
    "DC.c08.LAMP3" = c("LAMP3", "CCR7", "CCL22"),
    # IRF4保留在原始Excel中，但不作为正式pDC marker；见excluded_entries。
    "pDC.c09.LILRA4" = c("LILRA4", "JCHAIN")
  ),

  # Fib_markers.xlsx实际是广义stromal表，因此内皮相关cluster被拆到Endothelial面板。
  "Fibroblast" = list(
    "Stromal.c01.RGS5_pericyte" = c("RGS5", "NDUFA4L2", "MCAM", "NOTCH3"),
    "myoCAF.c03.SOX6" = c("SOX6", "CXCL14", "BMP4", "WNT5A"),
    "Stromal.c06.BMP5" = c("BMP5", "NSG1", "HSD17B2", "AGT"),
    "iCAF.c07.MMP1" = c("MMP3", "MMP1", "CXCL1"),
    "Stromal.c08.transiting" = c("RPL34", "RPS19", "RPS12"),
    "Stromal.c09.active" = c("JUN", "FOS", "FOSB"),
    "myoCAF.c10.CTHRC1" = c("CTHRC1", "COL8A1", "FBLN2"),
    "Stromal.c11.PI16" = c("PI16", "ADH1B", "SFRP1"),
    "Stromal.c12.ADAMDEC1" = c("ADAMDEC1", "CFD", "APOE"),
    "Stromal.c13.IL24" = c("IL24", "PI15", "IL1B"),
    "Stromal.c14.MYH11_SMC" = c("MYH11", "ACTA2", "ACTG2"),
    "SWC.c15.CDH19" = c("CDH19", "S100B", "CRYAB")
  ),

  # EC_markers.xlsx为内皮主面板，并合并Fib_markers.xlsx中的3组内皮marker。
  "Endothelial" = list(
    "EC.c01.GJA5_arterial" = c("GJA4", "GJA5", "SEMA3G", "FBLN5", "EFNB2"),
    "EC.c02.ACKR1_venous" = c(
      "ACKR1", "POSTN", "CPE", "IL1R1", "SELE", "SELP"
    ),
    "EC.c03.CAVIN1_venous" = c("CAVIN1", "CAVIN3", "IGFBP5", "PRCP"),
    "EC.c04.SELE_HEV" = c("SELE", "SERPINE1", "COTL1", "COL6A2"),
    "EC.c05.CD36_capillary" = c("CD36", "CD320", "BTNL9", "TMEM88"),
    "EC.c06.CA2_capillary" = c("CA2", "FOLH1", "KDR", "PIEZO2"),
    "EC.c07.LYVE1_lymphatic" = c("CCL21", "LYVE1", "RELN", "PDPN"),
    "EC.c08.ESM1_tip" = c(
      "ESM1", "PGF", "ANGPT2", "APLN", "PLVAP", "VWF", "PECAM1"
    ),
    "EC.c09.proliferating" = c("MKI67", "CCNB1"),
    "EC.c10.transiting" = c("RPS15", "RPL12", "TMSB10")
  )
)


# 原Excel中保留但不进入正式点图的条目，运行时会单独输出审计CSV。
subcluster_marker_excluded_entries <- data.frame(
  source_file = c("PC_markers(1).xlsx", "PC_markers(1).xlsx", "Mye_markers(1).xlsx"),
  original_cluster = c("PC.c01.IgA+", "PC.c01.IgA+", "pDC.c09.LILRA4"),
  gene = c("IGHA3", "IGHA4", "IRF4"),
  reason = c(
    "不是有效的人类IgA恒定区基因；人类正式基因符号为IGHA1和IGHA2",
    "不是有效的人类IgA恒定区基因；人类正式基因符号为IGHA1和IGHA2",
    "IRF4并非稳定特异的pDC标记；正式面板改用CLEC4C、TCF4、GZMB、IL3RA和IRF7"
  ),
  stringsAsFactors = FALSE
)


# ----------------------------------------------------------------------------
# 结直肠癌及肿瘤微环境单细胞文献补充marker
# ----------------------------------------------------------------------------

subcluster_marker_literature_panels <- list(
  "T/NK" = list(
    "CD4.c01.Tn(CCR7)" = c("TCF7", "LEF1", "MAL", "LTB"),
    "CD4.c02.Tcm(ANXA1)" = c("IL7R", "LTB", "MAL"),
    "CD4.c03.Th17(IL17A)" = c("KLRB1", "CCR6", "CCL20"),
    "CD4.c04.Tfh(CXCL13)" = c("PDCD1", "ICOS", "SH2D1A"),
    "CD4.c05.Treg(FOXP3)" = c("CTLA4", "TIGIT", "TNFRSF4"),
    "CD4.c06.active(NR4A2)" = c("NR4A1", "NR4A2", "DUSP1"),
    "CD4.c07.stress(HSPA1A)" = c("DNAJB1", "HSPA6"),
    "CD4.c09.Tex(LAG3)" = c("PDCD1", "TOX", "CTLA4", "CXCL13"),
    "CD8.c10.Trm(XCL2)" = c("CD69", "ZNF683", "CXCR6"),
    "CD8.c11.Tem(GZMK)" = c("CCL5", "EOMES", "CD44"),
    "CD8.c12.Teff(GZMB)" = c("PRF1", "GNLY", "CTSW", "FGFBP2"),
    "CD8.c13.Tex(HAVCR2)" = c("PDCD1", "TOX", "CXCL13", "CTLA4"),
    "CD8.c14.MAIT(SLC4A10)" = c("ZBTB16", "TRAV1-2", "NKG7"),
    "CD8.c15.stress(HSPA1A)" = c("HSPA1A", "HSPA1B", "DNAJB1"),
    "CD8.c16.active(CD69)" = c("CD69", "JUN", "FOSB"),
    "Innate.Immune.Cell" = c("TRGC2", "KLRD1", "NKG7"),
    "T.c18.STMN1" = c("TOP2A", "TYMS", "TUBA1B"),
    "NK.c19.FCGR3A" = c("KLRD1", "FCGR3A", "GNLY", "NKG7", "PRF1", "FGFBP2"),
    "gdT.c20.TRDC" = c("TRDC", "TRGC1", "TRGC2", "KLRD1", "NKG7")
  ),

  "B cell" = list(
    "Bn.c01.TCL1A" = c("IL4R", "HVCN1", "CD22", "MS4A1"),
    "Bn.c02.NR4A2" = c("CD83", "FOS", "JUN", "DUSP1"),
    "Bm.c03.TXNIP" = c("CD37", "CD44", "BANK1"),
    "Bm.c04.GPR183" = c("CD27", "CD22", "MS4A1"),
    "Bm.c05.HSPA1A" = c("HSPA6", "FOSB"),
    "Bgc.c06.BCL6" = c("MEF2B", "SERPINA9", "MKI67"),
    "B_atypical.c07.FCRL5" = c("FCRL4", "FCRL5", "ITGAX", "TBX21")
  ),

  "Plasma" = list(
    "PC.c01.IgA+" = c("JCHAIN", "MZB1", "DERL3", "SDC1"),
    "PC.c02.IgG+" = c("JCHAIN", "MZB1", "DERL3", "SDC1"),
    "PC.c03.early" = c("HLA-DRA", "CD74", "JCHAIN", "MZB1"),
    "PC.c04.cycling" = c("TOP2A", "PCNA", "TUBA1B"),
    "PC.c05.IgM+" = c("JCHAIN", "MZB1", "DERL3"),
    "PC.c06.transisting" = c("IGKC", "IGLC2", "CD79A", "MS4A1", "MZB1")
  ),

  "Myeloid" = list(
    "Macro.c01.FCN1" = c("S100A9", "CTSS", "LYZ", "THBS1"),
    "Macro.c02.C1QA" = c("C1QC", "APOC1", "MRC1", "CD68"),
    "Macro.c03.SPP1" = c("APOC1", "LGALS3", "CTSB", "CD68"),
    "Macro.c04.IL1A" = c("IL1B", "CXCL8", "NLRP3", "S100A9"),
    "Macro.c05.HSPA1B" = c("HSPA6", "DNAJB1", "FOSB"),
    "Macro.c06.TYMS" = c("TOP2A", "STMN1", "TUBA1B"),
    "cDC.c07.CD1C" = c("FCER1A", "CD1E", "HLA-DRA"),
    "DC.c08.LAMP3" = c("FSCN1", "CCL19", "CD83", "IDO1"),
    "pDC.c09.LILRA4" = c("GZMB", "TCF4", "CLEC4C", "IL3RA", "IRF7"),
    "cDC1.c10.CLEC9A" = c("CLEC9A", "XCR1", "BATF3", "CADM1", "DNASE1L3"),
    "Neutrophil.c11.FCGR3B" = c("FCGR3B", "CSF3R", "CXCR2", "S100A8", "S100A9")
  ),

  "Fibroblast" = list(
    "Stromal.c01.RGS5_pericyte" = c("CSPG4", "RBP1", "PDGFRB", "COL4A1"),
    "myoCAF.c03.SOX6" = c("F3", "COL1A1", "COL1A2", "DCN", "ACTA2", "TAGLN"),
    "Stromal.c06.BMP5" = c("COL15A1", "COL14A1", "DCN", "LUM"),
    "iCAF.c07.MMP1" = c("CXCL8", "IL6", "WNT5A", "MMP10"),
    "Stromal.c09.active" = c("JUNB", "DUSP1", "EGR1"),
    "myoCAF.c10.CTHRC1" = c(
      "POSTN", "MMP11", "SULF1", "COL1A1", "WNT5A",
      "COL11A1", "INHBA", "FAP"
    ),
    "Stromal.c11.PI16" = c("CFD", "COL14A1", "DPT", "OGN", "DCN"),
    "Stromal.c12.ADAMDEC1" = c("ADAM28", "CXCL14", "EDNRB", "PROCR"),
    "Stromal.c13.IL24" = c("CXCL2", "CCL2", "CXCL8"),
    "Stromal.c14.MYH11_SMC" = c("TAGLN", "CNN1", "DES", "RERGL"),
    "SWC.c15.CDH19" = c("SOX10", "PLP1", "S100A1", "MPZ"),
    "CAF.c16.FAP" = c("FAP", "PDPN", "COL11A1", "INHBA", "THY1")
  ),

  "Endothelial" = list(
    "EC.c01.GJA5_arterial" = c("SOX17", "KDR", "EFNB2"),
    "EC.c02.ACKR1_venous" = c("VWF", "PECAM1", "PLVAP", "NR2F2"),
    "EC.c03.CAVIN1_venous" = c("VWF", "PECAM1", "EMCN", "KDR", "RGCC"),
    "EC.c04.SELE_HEV" = c("ACKR1", "SELP", "ICAM1", "VCAM1", "PLVAP", "VWF"),
    "EC.c05.CD36_capillary" = c("CA4", "RGCC", "EMCN", "KDR", "PLVAP"),
    "EC.c06.CA2_capillary" = c("RGCC", "EMCN", "PECAM1", "VWF"),
    "EC.c07.LYVE1_lymphatic" = c("PROX1", "FLT4", "CCL21", "PDPN"),
    "EC.c08.ESM1_tip" = c("NID2", "KDR", "PLVAP", "VWF", "PECAM1"),
    "EC.c09.proliferating" = c("TOP2A", "TYMS", "STMN1", "PCNA")
  ),

  "Mast" = list(
    "Mast.c01.TPSAB1" = c("TPSAB1", "TPSB2", "CPA3", "MS4A2", "KIT", "HDC"),
    "Mast.c02.activated" = c("IL1RL1", "HPGDS", "LTC4S", "AREG", "FOS", "JUN")
  ),

  "Epithelial" = list(
    "Epi.c01.stem_TA" = c("LGR5", "OLFM4", "SMOC2", "SOX9", "MKI67"),
    "Epi.c02.enterocyte" = c("FABP1", "ALDOB", "APOA1", "KRT20", "CA1"),
    "Epi.c03.goblet" = c("MUC2", "SPINK4", "CLCA1", "FCGBP", "ZG16"),
    "Epi.c04.BEST4_OTOP2" = c("BEST4", "OTOP2", "CA7", "GUCA2A", "HES4"),
    "Epi.c05.tuft" = c("POU2F3", "TRPM5", "AVIL", "IL17RB", "SOX9"),
    "Epi.c06.enteroendocrine" = c("CHGA", "CHGB", "NEUROD1", "PAX6", "PCSK1N"),
    "Epi.c07.cycling" = c("MKI67", "TOP2A", "STMN1", "TYMS", "TUBA1B"),
    "Epi.c08.EMT_like" = c("VIM", "TGFBI", "LAMC2", "KRT17", "ITGA2", "S100A10"),
    "Epi.c09.secretory_tumor" = c("LYZ", "LCN2", "REG1A", "REG3A", "AGR2")
  )
)


# ----------------------------------------------------------------------------
# 合并原始marker和文献补充marker
# ----------------------------------------------------------------------------

merge_marker_panel_sets <- function(primary, supplemental) {
  big_groups <- unique(c(names(primary), names(supplemental)))
  result <- setNames(vector("list", length(big_groups)), big_groups)
  for (big_group in big_groups) {
    primary_group <- primary[[big_group]]
    supplemental_group <- supplemental[[big_group]]
    subtype_names <- unique(c(names(primary_group), names(supplemental_group)))
    result[[big_group]] <- setNames(vector("list", length(subtype_names)), subtype_names)
    for (subtype in subtype_names) {
      genes <- unique(c(primary_group[[subtype]], supplemental_group[[subtype]]))
      genes <- as.character(genes[!is.na(genes) & nzchar(genes)])
      result[[big_group]][[subtype]] <- genes
    }
  }
  result
}

subcluster_marker_panels <- merge_marker_panel_sets(
  subcluster_marker_uploaded_panels,
  subcluster_marker_literature_panels
)


# ----------------------------------------------------------------------------
# 文献来源记录：运行时会原样输出到结果目录
# ----------------------------------------------------------------------------

subcluster_marker_reference_sources <- data.frame(
  reference_id = c("CRC_immune_2021", "CRC_Tcell_2023", "CRC_stromal_2022",
                   "CRC_progression_2023", "CRC_BPC_2024", "CRC_EC_2024",
                   "CRC_epithelial_2022", "CRC_spatial_2021"),
  scope = c("T/NK、B、髓系", "T细胞", "CAF、巨噬细胞、DC", "髓系/免疫",
            "B/Plasma", "Endothelial", "Epithelial", "CRC免疫生态"),
  title = c(
    "Single-cell analyses reveal suppressive tumor microenvironment of human colorectal cancer",
    "Single-cell transcriptome analysis reveals T population heterogeneity and functions in colorectal cancer metastases",
    "Single-cell and spatial analysis reveal interaction of FAP+ fibroblasts and SPP1+ macrophages in colorectal cancer",
    "Dynamic heterogeneity of colorectal cancer during progression revealed clinical risk-associated cell types and regulations",
    "Single-cell transcriptome analysis reveals immunosuppressive landscape in overweight and obese colorectal cancer",
    "Endothelial cell heterogeneity in colorectal cancer: tip cells drive angiogenesis",
    "Single-cell Transcriptomics Reveals Early Molecular and Immune Alterations Underlying the Serrated Neoplasia Pathway",
    "Spatially organized multicellular immune hubs in human colorectal cancer"
  ),
  url = c(
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC8181206/",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC10394913/",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC8976074/",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC10290555/",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC10838453/",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC11342913/",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC9791140/",
    "https://pubmed.ncbi.nlm.nih.gov/34450029/"
  ),
  stringsAsFactors = FALSE
)

set.seed(1234)
