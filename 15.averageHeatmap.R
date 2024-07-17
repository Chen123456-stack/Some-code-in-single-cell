
#选取特点基因做热图
selected_genes_list <- list(
  'CD4 T cell' = c('CD3D','CD3E',"CD4", "IL7R", "CCR7"),
  'CD8 T cell' = c("CD8A", "CD8B", "IFNG",'GZMK','KLRG1'),
  'NK' = c("KLRF1", "NCR1", "GNLY",'NCAM1','XCL1'),
  'B cell' = c("MS4A1", "CD79A", "PAX5",'IGHD','CD22'),
  'Plasma cell' = c("IGHA1", "IGHG1", "IGHG4",'IGHG3','IGHG2'),
  'Myeloid cell' = c("CD14", "CD163", "C1QB",'C1QA','CLEC10A'),
  'pDC' = c("SCT","LILRA4","PACSIN1","CLEC4C","MAP1A"),
  'Neutrophil' = c("FCGR3B", "CXCR2", "CSF3R",'CXCR1','FFAR2'),
  'Basophil' = c("GATA2", "HDC", "GCSAML",'MS4A2','CPA3'),
  'Epithlial cell' = c("TMPRSS2", "EPCAM", "KRT18",'KRT8','LGALS4'),
  'Endothelial cell' = c("TM4SF1", "AQP1", "TSPAN7",'CLEC14A','TM4SF18'),
  'Hepatocyte' = c("ALB", "ALDOB", "FABP1",'ADH4','APOC3'),
  'Fibroblast' = c("COL1A2", "DCN", "PCDH18",'PLN','TCF21')
)

# 从每个 cluster 中选取特定基因
top5 <- cell.markers %>%
  filter((cluster == "CD4 T cell" & gene %in% selected_genes_list$`CD4 T cell`) |
           (cluster == "CD8 T cell" & gene %in% selected_genes_list$`CD8 T cell`) |
           (cluster == "NK" & gene %in% selected_genes_list$NK) |
           (cluster == "B cell" & gene %in% selected_genes_list$`B cell`) |
           (cluster == "Plasma cell" & gene %in% selected_genes_list$"Plasma cell") |
           (cluster == "Myeloid cell" & gene %in% selected_genes_list$"Myeloid cell") |
           (cluster == "pDC" & gene %in% selected_genes_list$pDC) |
           (cluster == "Neutrophil" & gene %in% selected_genes_list$Neutrophil) |
           (cluster == "Basophil" & gene %in% selected_genes_list$Basophil) |
           (cluster == "Epithlial cell" & gene %in% selected_genes_list$`Epithlial cell`) |
           (cluster == "Endothelial cell" & gene %in% selected_genes_list$"Endothelial cell") |
           (cluster == "Hepatocyte" & gene %in% selected_genes_list$Hepatocyte) |
           (cluster == "Fibroblast" & gene %in% selected_genes_list$Fibroblast))

table(sc_pbmc_int2$cell_type)
#选取特点基因做热图
selected_genes_list <- list(
  'CD4_C1_LEF1' = c("LEF1", "CCR7", "TRABD2A",'MAL','ACTN1'),
  'CD4_C2_FOXP3' = c("FOXP3", "IL2RA", "CTLA4",'TBC1D4','RTKN2'),
  'CD4_C3_IL7R' = c("CD4", "CD40LG", "GPR183",'IL6R','IL7R'),
  'CD8_C1_RGS1' = c("RGS1", "HLA-DQA1", "CTLA4",'PDCD1','DUSP4'),
  'CD8_C2_GPR183' = c("GPR183", "MYB", "SESN3",'FGFR1','PGM2L1'),
  'CD8_C3_CX3CR1' = c("FGFBP2","GZMH","CX3CR1","PRSS23","ADGRG1"),
  'CD8_C4_GZMK' = c("KLRB1", "KLRG1", "GZMK",'ZBTB16','SLC4A10'),
  'CD8_C5_LEF1' = c("LEF1", "TRGC2", "TRDC",'IKZF2','FOS'),
  'CD8_C6_MKI67' = c("MKI67", "TOP2A", "TYMS",'PCLAF','ASPM'),
  'gdT' = c("CD160", "TRGC2", "IKZF2",'CRTAM','TRDC'),
  'NK_C1_ITGAX' = c("ITGAX", "GNLY", "NCAM1",'XCL1','ITGAM'),
  'NK_C2_FCGR3A' = c("FGFBP2", "FCGR3A", "GNLY",'CX3CR1','GZMB'),
  'NK_C3_XCL1' = c("FCER1G", "KLRF1", "CEBPD",'CCL3','XCL1')
)

# 从每个 cluster 中选取特定基因
top5 <- cell.markers %>%
  filter((cluster == "CD4_C1_LEF1" & gene %in% selected_genes_list$CD4_C1_LEF1) |
           (cluster == "CD4_C2_FOXP3" & gene %in% selected_genes_list$CD4_C2_FOXP3) |
           (cluster == "CD4_C3_IL7R" & gene %in% selected_genes_list$CD4_C3_IL7R) |
           (cluster == "CD8_C1_RGS1" & gene %in% selected_genes_list$CD8_C1_RGS1) |
           (cluster == "CD8_C2_GPR183" & gene %in% selected_genes_list$CD8_C2_GPR183) |
           (cluster == "CD8_C3_CX3CR1" & gene %in% selected_genes_list$CD8_C3_CX3CR1) |
           (cluster == "CD8_C4_GZMK" & gene %in% selected_genes_list$CD8_C4_GZMK) |
           (cluster == "CD8_C5_LEF1" & gene %in% selected_genes_list$CD8_C5_LEF1) |
           (cluster == "CD8_C6_MKI67" & gene %in% selected_genes_list$CD8_C6_MKI67) |
           (cluster == "gdT" & gene %in% selected_genes_list$gdT) |
           (cluster == "NK_C1_ITGAX" & gene %in% selected_genes_list$NK_C1_ITGAX) |
           (cluster == "NK_C2_FCGR3A" & gene %in% selected_genes_list$NK_C2_FCGR3A)|
           (cluster == "NK_C3_XCL1" & gene %in% selected_genes_list$NK_C3_XCL1))


table(sc_pbmc_int2$cell_type)
#选取特点基因做热图
selected_genes_list <- list(
  'Mac_C1_CCL3' = c("CCL3", "IER3", "CCL3L3",'PRDM1','CXCL2'),
  'Mac_C2_CXCL12' = c("CD5L", "LILRB5", "CXCL12",'SELENOP','VCAM1'),
  'Mono_C1_CD14' = c("VCAN", "S100A12", "S100A8",'S100A9','CKAP4'),
  'Mono_C2_FCGR3A' = c("CDKN1C", "FCGR3A", "PELATON",'CX3CR1','SIGLEC10'),
  'DC_C1_CD1C' = c("CD1C", "FCER1A", "FCGR2B",'PLD4','SLC38A1'),
  'DC_C2_CLEC9A' = c("CLEC9A","CLNK","CADM1","XCR1","DNASE1L3"),
  'DC_C3_MKI67' = c("MKI67", "TYMS", "TOP2A",'PCLAF','TK1')
)

# 从每个 cluster 中选取特定基因
top5 <- cell.markers %>%
  filter((cluster == "Mac_C1_CCL3" & gene %in% selected_genes_list$Mac_C1_CCL3) |
           (cluster == "Mac_C2_CXCL12" & gene %in% selected_genes_list$Mac_C2_CXCL12) |
           (cluster == "Mono_C1_CD14" & gene %in% selected_genes_list$Mono_C1_CD14) |
           (cluster == "Mono_C2_FCGR3A" & gene %in% selected_genes_list$Mono_C2_FCGR3A) |
           (cluster == "DC_C1_CD1C" & gene %in% selected_genes_list$DC_C1_CD1C) |
           (cluster == "DC_C2_CLEC9A" & gene %in% selected_genes_list$DC_C2_CLEC9A) |
           (cluster == "DC_C3_MKI67" & gene %in% selected_genes_list$DC_C3_MKI67))

# 创建DataFrame
top5 <- data.frame(
  State = c(rep("Naive", 4), rep("Resident", 3), rep("Inhibitory", 5), rep("Cytokines", 11), 
            rep("Co-stimulatory", 3), rep("TF", 8), rep("Treg", 3), rep("Cell Type", 8)),
  gene = c("CCR7", "TCF7", "LEF1", "SELL", 
           "CD69", "RUNX3", "NR4A1", 
           "TIGIT", "CTLA4", "LAG3", "PDCD1", "HAVCR2", 
           "IL2", "IL17A", "IL17F", "LAMTOR3", "GNLY", "GZMB", "PRF1", "GZMK", "IFNG", "GZMA", "NKG7", 
           "TNFRSF14", "CD28", "ICOS", 
           "TBX21", "ZNF683", "ZEB2", "ID2", "EOMES", "HOPX", "HIF1A", "TOX", 
           "FOXP3", "IL2RA", "IKZF2", 
           "KLRF1", "KLRD1", "CD4", "CD3D", "CD3G", "CD3E", "CD8A", "CD8B")
)



NK <- c("IFNG", "CCL4", "CCL4L2", "CCL5", "IL32", "IL16", "CCL3", "XCL2", "XCL1", "FLT3LG", "S1PR1", "CXCR3", 
  "CXCR4", "S1PR5", "CX3CR1", "CXCR2", "S1PR4", "KLRK1", "KLRC2", "NCR1", "CD160", "NCR3", "SLAMF6", 
  "SLAMF7", "KLRC1", "CD300A", "TIGIT", "KIR3DL2", "KIR3DL1", "KIR2DL3", "KIR2DL1", "KLRB1", "SIGLEC7", 
  "HAVCR2", "TGFBR1", "IL12RB1", "IL2RG", "TGFBR3", "TGFBR2", "IL10RA", "IL18R1", "IL18RAP", "IL2RB", 
  "IL10RB", "FASLG", "GSDMD", "GZMB", "NKG7", "PRF1", "GZMA", "GZMH", "GZMK", "TNFSF10", "CD69", "TNFRSF18", 
  "ITGAM", "ITGAX", "NCAMI", "FCGR3A")

c1 <- c("ITGAM", "ITGAX", "NCAM1", "FCGR3A")

markers <- c('CXCR1','CXCR2','CXCR3','CXCR4','CXCR5','CXCR6','CCR1','CCR2','CCR3','CCR4',
             'CCR6','CCR7','CCR9','CCR10','XCR1','CXCR1')

markers <- c('HLA-A','HLA-B','HLA-C','HLA-E','HLA-F')
markers <- c('HLA-A','HLA-B','HLA-C','HLA-E','HLA-F',
             'HLA-DMA','HLA-DMB','HLA-DOA','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQA2','HLA-DQB1','HLA-DRA','HLA-DRB1','HLA-DRB5','HLA-DRB6')

t <- sc_pbmc_int2
sc_pbmc_int2 <- t

table(sc_pbmc_int2$cell_type)
sc_pbmc_int2 <- subset(sc_pbmc_int2,group%in%c('LI_AR','LI_NAR')&cell_type%in%c('14:Mac_C3_C1QA','15:Mac_C4_CXCL12',
                                                                                '16:Mac_C5_CD5L','17:Mac_C6_TIMD4'))



sc_pbmc_int2 <- subset(sc_pbmc_int2,group%in%c('LI_AR','LI_NAR')&cell_type%in%c("4:Mono_C1_VCAN",
                                                                                "5:Mono_C2_S100A12", "6:Mono_C3_FCGR3A","7:Mono_C4_ATF3","8:Mono_C5_FCN1",
                                                                                "9:Mono_C6_LRMDA","10:Mono_C7_RPS13","11:Mono_C8_RTF2","12:Mac_C1_HLA-DRA","13:Mac_C2_RGS1",
                                                                                '14:Mac_C3_C1QA','15:Mac_C4_CXCL12','16:Mac_C5_CD5L','17:Mac_C6_TIMD4'))


sc_pbmc_int2 <- subset(sc_pbmc_int2,group%in%c('LI_AR','LI_NAR')&cell_type%in%c('3:CD8_C1_GZMK',
                                                                                '4:CD8_C2_RGS1',
                                                                                '5:CD8_C3_GZMH',
                                                                                '6:CD8_C4_LEF1',
                                                                                '7:CD8_C5_KLRB1',
                                                                                '8:CD8_C6_TCF7',
                                                                                '9:CD8_C7_FOSB',
                                                                                '10:CD8_C8_MKI67',
                                                                                '11:CD8_C9_TIGIT'))




metadata <- sc_pbmc_int2@meta.data

# 创建一个新的列，将 cell_type 和 group 组合起来
metadata$cell_group <- paste(metadata$cell_type, metadata$group, sep = "_")

# 将新的元数据添加回 Seurat 对象
sc_pbmc_int2 <- AddMetaData(sc_pbmc_int2, metadata = metadata)

table(sc_pbmc_int2$cell_group)

# sc_pbmc_int2 <- CD8_Obj

mean_gene_exp <- as.matrix(
  data.frame(
    Seurat::AverageExpression(sc_pbmc_int2,
                              features = markers,
                              group.by = 'cell_group',
                              assays = 'RNA',
                              layer = 'data'
    )
  )
)

Idents(sc_pbmc_int2) <- sc_pbmc_int2$cell_group

#注意改列明的时候要保持一致
colnames(mean_gene_exp) <- levels(Seurat::Idents(sc_pbmc_int2))

#可以用这个方法
x <- table(sc_pbmc_int2$cell_group)
x
x <- as.data.frame(x)

colnames(mean_gene_exp) <- c(x$Var1)


# Z-score
htdf <- t(scale(t(mean_gene_exp), scale = TRUE, center = TRUE))
# htdf <- htdf[top5$gene, ]
htdf <- htdf[markers, ]




# color
htCol = c("#0099CC", "white", "#CC0033")
# htCol = c("blue", "white", "orange")
htRange = c(-2, 0, 2)
col_fun <- circlize::colorRamp2(htRange, htCol)



my.breaks <- c(seq(-1.5, -0.2, by=0.1), seq(0,1.5, by=0.1))
my.colors <- c(
  colorRampPalette(colors = c("#377EB8", "white"))(length(my.breaks)/2),
  colorRampPalette(colors = c("white", "#E41A1C"))(length(my.breaks)/2))
col_fun <- circlize::colorRamp2(my.breaks, my.colors)


annoColType = "light"
annoColTypeAlpha = 0

set.seed(666)
FJ <- c('#4C8562','#A24F1D','#65A8AE','#363C6C','#3D5C96','#5E5598','#67912F',
        '#806224','#D08685','#A7256F','#834626','#6376B9','#ABCD79')

anno_col <- FJ

print(c("Your cluster annotation color is:", anno_col))


names(anno_col) <- colnames(htdf)

clusterAnnoName = TRUE
# top annotation
column_ha <- ComplexHeatmap::HeatmapAnnotation(
  cluster = colnames(htdf),
  show_legend = FALSE,
  show_annotation_name = clusterAnnoName,
  col = list(cluster = anno_col)
)

# all genes
rowGene <- rownames(htdf)

# tartget gene
annoGene <- top5$gene

# get target gene index
index <- match(annoGene, rowGene)

# some genes annotation
geneMark <- ComplexHeatmap::rowAnnotation(
  gene = ComplexHeatmap::anno_mark(
    at = index,
    labels = annoGene,
    labels_gp = grid::gpar(
      fontface = "italic",
      fontsize = 10
    )
  )
)
right_annotation <- geneMark


column_split <- c(7,12)
class(column_split)
column_split <- as.integer(column_split)


library(ComplexHeatmap)
library(grid)

setwd('./1.nk and t/NK/')


my.breaks <- c(seq(-1.5, -0.2, by=0.1), seq(0,1.5, by=0.1))
my.colors <- c(
  colorRampPalette(colors = c("#377EB8", "white"))(length(my.breaks)/2),
  colorRampPalette(colors = c("white", "#E41A1C"))(length(my.breaks)/2))
col_fun <- circlize::colorRamp2(my.breaks, my.colors)

colnames(htdf)


###对坐标进行排序，按顺序进行排序
data <- colnames(htdf)

# 将数据转换为数据框，并拆分成数字和后缀部分
df <- data.frame(original = data) %>%
  mutate(number = as.numeric(sub(":.*", "", original)),
         suffix = ifelse(grepl("LI_AR", original), "LI_AR", "LI_NAR"))

# 按照数字和后缀排序
sorted_df <- df %>%
  arrange(number, suffix)

# 获取排序后的数据
sorted_data <- sorted_df$original

sorted_data


# 自定义排序向量
custom_order <- c(7, 12, 13, 4, 5, 6, 8, 9, 10, 11)

# 提取行名中排序顺序的部分
extract_order <- as.integer(sub(":.*", "", sorted_data))

# 根据自定义排序向量重新排序行名
ordered_row_names <- sorted_data[order(match(extract_order, custom_order))]

# 查看排序后的行名
ordered_row_names

# 重新排列data的列
htdf <- htdf[, sorted_data]



# 获取行数和列数,用于分割热图
n_rows <- nrow(t(htdf))
n_cols <- ncol(t(htdf))

# 创建行分割向量
split_vector <- rep(1, n_rows)
split_vector[7:n_rows] <- 2
split_vector[9:n_rows] <- 3
split_vector[17:n_rows] <- 4
# 创建列分割向量
split_vector1 <- rep(1, n_cols)
split_vector1[6:n_cols] <- 2
# split_vector1[8:n_cols] <- 3



dev.off()
pdf( file = "XCL.pdf", width = 20, height = 10 )
 ComplexHeatmap::Heatmap(
  t(htdf),
  name = "Z-score",
  cluster_columns = FALSE,
  cluster_rows = F,
  # row_title = "Cluster top Marker genes",
  right_annotation = NULL,
  show_row_names = TRUE,
  column_names_gp = grid::gpar(
    fontface = "italic",
    fontsize = 10
  ),
  row_names_gp = grid::gpar(
    # fontface = "italic",
    fontsize = 10
  ),
  row_names_side = "left",
  border = T,
  column_names_side = 'bottom',
  column_names_rot =90,
  top_annotation = NULL,  # 如果有列注释，可以在这里添加
  col = col_fun,
  # row_split = NULL,  # 如果有行分割信息，可以在这里添加
  width = unit(1, "cm") * ncol(t(htdf)), # 设置单元格宽度
  height = unit(1, "cm") * nrow(t(htdf)) # 设置单元格高度
  # row_split=split_vector
  # column_split = split_vector1 # 使用定义的 column_split
  # gap = unit(4, "mm")  # 设置列之间的间隙
)

dev.off()

###放置图例用
draw(t, 
     heatmap_legend_side = "top", 
     heatmap_legend_list = list(Legend(col_fun = col_fun, title = "Z-score", direction = "horizontal"))
)



dev.off()
pdf( file = "5.5-T_heatmap.pdf", width = 8, height = 15 )
  ComplexHeatmap::Heatmap(htdf,
    col_split = htdf[1, ],
    gap = 30,
  name = "Z-score",
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  row_title = "Cluster top Marker genes",
  # column_title = "Clusters",
  right_annotation = NULL,
  show_row_names = TRUE,
  row_names_gp = grid::gpar(
    fontface = "italic",
    fontsize = 10
  ),
  row_names_side = "right",
  border = FALSE,
  column_names_side = "top",
  column_names_rot = 45,
  top_annotation = column_ha,
  col = col_fun,
  row_split=top5$State
)
dev.off()





# 创建列注释数据框
signatureType_row <- paste0("Type", 1:13)
signatureType_col <- data.frame(Type = signatureType_row)
rownames(signatureType_col) <- colnames(htdf1)
signatureType_col$Type <- rownames(signatureType_col)
signatureType_col$Type <- factor(signatureType_col$Type)

class(htdf)
htdf1 <- as.data.frame(htdf)


my.breaks <- c(seq(-1, -0.5, by=0.01), seq(-0.4, 1, by=0.01))
my.breaks <- c(seq(-1.5, -0.2, by=0.1), seq(-0.1, 1.5, by=0.1))
my.colors <- c(
  colorRampPalette(colors = c("#377EB8", "white"))(length(my.breaks)/2),
  colorRampPalette(colors = c("white", "#E41A1C"))(length(my.breaks)/2))


my.colors <- c(
  colorRampPalette(colors = c("#0571B0", "white"))(length(my.breaks)/2),
  colorRampPalette(colors = c("white", "#CA0020"))(length(my.breaks)/2))

my.colors <- c(
  colorRampPalette(colors = c("#4DBBD5FF", "white"))(length(my.breaks)/2),
  colorRampPalette(colors = c("white", "#DC0000FF"))(length(my.breaks)/2))



dev.off()
pdf( file = "5.5-T_heatmap.pdf", width = 10, height = 20 )
pheatmap::pheatmap(t(htdf),
                   show_colnames = T,
                   show_rownames = T, annotation_names_col = F,
                   # annotation_col = signatureType_col,
                   #annotation_row = signatureType_row,
                   # gaps_row = c(3, 14, 17),
                   # gaps_col = c(1,2,3,4,5,6,7,8,9,10,11,12),
                   cluster_rows = F,
                   cluster_cols = F,
                   breaks = my.breaks,
                   color = my.colors,
                   cellwidth = 10, cellheight = 10,
                   angle_col=315,
                   border_color = "black",
                   fontsize = 8)

dev.off()

###########
t <- htdf
table(abs(t)>2)
t[t>=2]=2
t[t<=-2]=-2
table(is.na(t))

# 绘制热图
library(pheatmap)


# 构建行注释信息（行名与表达矩阵的行名row保持一致）
annotation_row = data.frame(
  Function = c(rep("Naive", 4), rep("Resident", 3), rep("Inhibitory", 5), rep("Cytokines", 11), 
           rep("Co-stimulatory", 3), rep("TF", 8), rep("Treg", 3), rep("Cell Type", 8)),
row.names = rownames(t))
head(annotation_row)

# 构建行注释信息（行名与表达矩阵的行名row保持一致）
annotation_col = data.frame(
  celltype = c(rep("CD4T", 3), rep("CD8T", 6), rep("gdT", 1), rep("NK", 3)),
  row.names = colnames(t))
head(annotation_col)

#修改注释标签的颜色
ann_colors = list(
  Function = c(Naive = "#E41A1C", Resident = "#377EB8", Inhibitory = "#4DAF4A", Cytokines = "#984EA3", 
            `Co-stimulatory` = "#FF7F00", TF = "#FFFF33", Treg = "#A65628", `Cell Type` = "#F781BF"),
  celltype = c(CD4T = "#999999", CD8T = "#66C2A5", gdT = "#FC8D62", NK = "#8DA0CB"))


dev.off()
pdf( file = "5.5-T_heatmap.pdf", width = 6, height = 8 )
pheatmap(t, 
         color = colorRampPalette(c("navy", "white", "#C9771D"))(100), 
         cluster_rows = FALSE, cluster_cols = FALSE,gaps_row = c(4,7,12,23,26,34,37,45),
         gaps_col = c(3,9,10),
         row_names_side = "left", 
         angle_col = c("315"),
         annotation_col = annotation_col,
         annotation_row = annotation_row)

dev.off()
