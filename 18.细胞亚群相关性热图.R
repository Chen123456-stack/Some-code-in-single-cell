library(Seurat)
library(tidyr)

# 假设sc_pbmc_int2已经是一个Seurat对象

# 查看meta.data的列名
colnames(sc_pbmc_int2@meta.data)

t <- sc_pbmc_int2

sc_pbmc_int2 <- t


sc_pbmc_int2 <- subset(sc_pbmc_int2,group%in%c('LI_NAR','PB_NAR'))
table(sc_pbmc_int2$cell_type,sc_pbmc_int2$group)

# 使用tidyr::unite将type和cell_type列合并为type_celltype列
sc_pbmc_int2@meta.data <- unite(sc_pbmc_int2@meta.data, 
                                "celltype_type", 
                                cell_type, group,
                                remove = FALSE)

# 查看type_celltype列的频率
table(sc_pbmc_int2@meta.data$celltype_type)

# 将identities设置为type_celltype
Idents(sc_pbmc_int2) <- sc_pbmc_int2@meta.data$celltype_type

# 计算平均表达量
exp <- AverageExpression(sc_pbmc_int2, group.by = 'celltype_type',layer = 'scale.data')$RNA

# 提取标准差最大的2000个基因
features <- names(tail(sort(apply(exp, 1, sd)), 5000))

# 选择这些基因的表达量数据
exp <- exp[which(row.names(exp) %in% features), ]

# 将exp转换为数据框
exp <- as.data.frame(exp)

# 计算Spearman相关系数
exp <- cor(exp, method = "spearman")

head(exp)
colnames(exp)


# 提取列名
colnames_exp <- colnames(exp)
head(colnames_exp)

# 提取LIVER和PBMC列名
liver_cols <- grep("PB-AR", colnames_exp, value = TRUE)
pbmc_cols <- grep("PB-NAR", colnames_exp, value = TRUE)

# 提取LIVER和PBMC后面的数字并排序
liver_order <- liver_cols[order(as.numeric(gsub(".*-PB-AR_", "", liver_cols)))]
pbmc_order <- pbmc_cols[order(as.numeric(gsub(".*-PB-NAR_", "", pbmc_cols)))]

# 合并排序后的列名
sorted_order <- c(liver_order, pbmc_order)

# 根据排序顺序重新排列数据框的列
exp_sorted <- exp[sorted_order, sorted_order]

# 检查重新排序后的列名
colnames(exp_sorted)



#行列注释
annotation_col = data.frame(
  celltype = c("1:CD4_C1_LEF1" ,  "2:CD4_C2_IL7R" ,  "3:CD8_C1_GZMK" ,  "4:CD8_C2_RGS1" ,  "5:CD8_C3_GZMH" ,  "6:CD8_C4_LEF1" , 
                "7:CD8_C5_KLRB1",  "8:CD8_C6_TCF7"  , "9:CD8_C7_FOSB" ,  "10:CD8_C8_MKI67", "11:CD8_C9_TIGIT", "12:NK_C1_FCGR3A",
                "13:NK_C2_XCL1" ,  "14:NK_C3_KLRF1" , "15:DN_LINC_PINT") ,
  group = c(rep("AR",15),rep("NAR",15))
)


row.names(annotation_col) <- colnames(exp_sorted)

annotation_row = data.frame(
  celltype = c("1:CD4_C1_LEF1" ,  "2:CD4_C2_IL7R" ,  "3:CD8_C1_GZMK" ,  "4:CD8_C2_RGS1" ,  "5:CD8_C3_GZMH" ,  "6:CD8_C4_LEF1" , 
               "7:CD8_C5_KLRB1",  "8:CD8_C6_TCF7"  , "9:CD8_C7_FOSB" ,  "10:CD8_C8_MKI67", "11:CD8_C9_TIGIT", "12:NK_C1_FCGR3A",
               "13:NK_C2_XCL1" ,  "14:NK_C3_KLRF1" , "15:DN_LINC_PINT") ,
  group = c(rep("AR",15),rep("NAR",15))
)
row.names(annotation_row) <- rownames(exp_sorted)

# 将 FJ 中的颜色按顺序填充到 annotation_colors 中的 celltype 和 type 类别




# 自定义 celltype 颜色顺序
annotation_colors <- list(
  celltype = c(
    "1:CD4_C1_LEF1" = FJ[1], "2:CD4_C2_IL7R" = FJ[2], "3:CD8_C1_GZMK" = FJ[3], "4:CD8_C2_RGS1" = FJ[4], "5:CD8_C3_GZMH" = FJ[5],
    "6:CD8_C4_LEF1" = FJ[6], "7:CD8_C5_KLRB1" = FJ[7], "8:CD8_C6_TCF7" = FJ[8], "9:CD8_C7_FOSB" = FJ[9], "10:CD8_C8_MKI67" = FJ[10],
    "11:CD8_C9_TIGIT" = FJ[11], "12:NK_C1_FCGR3A" = FJ[12], "13:NK_C2_XCL1" = FJ[13], "14:NK_C3_KLRF1" = FJ[14], "15:DN_LINC_PINT" = FJ[15]
  ),
  group = c("AR" = FJ[18], "NAR" = FJ[21])
)

# # 确认你的 annotation_col 包含 celltype 和 type 列
# annotation_col <- data.frame(
#   celltype = factor(colnames(exp_sorted), levels = c(
#     "1:CD4-C1-LEF1", "2:CD4-C2-IL7R", "3:CD8-C1-GZMK", "4:CD8-C2-RGS1",
#     "5:CD8-C3-GZMH", "6:CD8-C4-LEF1", "7:CD8-C5-KLRB1", "8:CD8-C6-TCF7",
#     "9:CD8-C7-FOSB", "10:CD8-C8-MKI67", "11:CD8-C9-TIGHT", "12:NK-C1-FCGR3A",
#     "13:NK-C2-XCL1", "14:NK-C3-KLRF1", "15:DN-LINC-PINT"
#   )),
#   type = rep(c("LIVER", "PBMC"), each = 15)
# )
# rownames(annotation_col) <- colnames(exp_sorted)


# 定义颜色断点
my.breaks <- c(seq(-1, -0.1, by = 0.01), seq(0.1, 1, by = 0.01))

# 定义颜色调色板
my.colors <- c(
  colorRampPalette(colors = c("#377EB8", "white"))(length(my.breaks)/2),
  colorRampPalette(colors = c("white", "#E41A1C"))(length(my.breaks)/2)
)


dev.off()
pdf( file = "T_PBMC_aso_heatmap.pdf", width = 12, height = 15 )
pheatmap::pheatmap(
  exp_sorted, 
  annotation_col = annotation_col,
  annotation_row = annotation_row, 
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = "#3B3838",
  gaps_col = c(15),
  gaps_row = c(15),
  color = my.colors,
  annotation_legend = TRUE,
  annotation_names_row = F,
  annotation_names_col = F,
  breaks = my.breaks,
  cellwidth = 12, cellheight = 12,
  annotation_colors = annotation_colors
)
dev.off()

colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)

#color = rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")),

# color = colorRampPalette(c("#0F0F74", "white","#FE3838"))(10),

# color = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")),




