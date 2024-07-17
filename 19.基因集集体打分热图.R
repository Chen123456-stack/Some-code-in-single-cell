setwd('./CD8')

CD8_data <- read_csv("CD81.csv")
marker.list = list()

for (col_name in colnames(CD8_data)) {
  genes <- CD8_data %>%
    pull(!!sym(col_name)) %>%
    na.omit() %>%
    as.vector()
  marker.list[[col_name]] <- genes
}

marker.list$Cytotoxicity

table(sc_pbmc_int2$cell_type)
# 子集化包含感兴趣的 CD8 细胞亚群
CD8_Obj <- subset(sc_pbmc_int2, cell_type %in% c('3:CD8_C1_GZMK', '4:CD8_C2_RGS1', '5:CD8_C3_GZMH', 
                                                 '6:CD8_C4_LEF1', '7:CD8_C5_KLRB1', '8:CD8_C6_TCF7', 
                                                 '9:CD8_C7_FOSB', '10:CD8_C8_MKI67', '11:CD8_C9_TIGIT'))

# 去除没有数据的细胞类型
CD8_Obj$cell_type <- droplevels(CD8_Obj$cell_type)

# 检查每个细胞类型的细胞数
table(CD8_Obj$cell_type)



CD8_Obj <-  AddModuleScore(CD8_Obj,
                           features = marker.list,
                           ctrl = 5,
                           name = "FunctionScore")
for(i in 1:length(marker.list)){
  colnames(CD8_Obj@meta.data)[colnames(CD8_Obj@meta.data) == paste0("FunctionScore", i)] <- names(marker.list)[i]
}
Idents(CD8_Obj) <- CD8_Obj$seurat_clusters
Differentiation <- c("Naive", "Activation/Effector function", "Exhaustion")
Function <- c("TCR Signaling", "Cytotoxicity", "Cytokine/Cytokine receptor",
              "Chemokine/Chemokine receptor", "Senescence", "Anergy",
              "NFKB Signaling", "Stress response", "MAPK Signaling", "Adhesion",
              "IFN Response")
Metabolism <- c("Oxidative phosphorylation", "Glycolysis", "Fatty acid metabolism")
Apoptosis <- c("Pro-apoptosis", "Anti-apoptosis")
MarkerNameVector <- c(Differentiation, Function, Metabolism, Apoptosis)
FunctionScoreMatrix <- matrix(0,
                              ncol = length(unique(CD8_Obj$seurat_clusters)),
                              nrow = length(marker.list))
colnames(FunctionScoreMatrix) <- paste0(rownames(table(CD8_Obj$cell_type)))
rownames(FunctionScoreMatrix) <- MarkerNameVector
for(ci in 1:ncol(FunctionScoreMatrix)){
  for(ri in 1:nrow(FunctionScoreMatrix)){
    FunctionVec <- as_tibble(CD8_Obj@meta.data) %>% pull(MarkerNameVector[ri])
    fv <- mean(FunctionVec[CD8_Obj$cell_type == levels(CD8_Obj$cell_type)[ci]])  #level比较重要，记得设置
    FunctionScoreMatrix[ri, ci] <- fv
  }
}

library(scales)
FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to=c(-1, 1)))
orderC = colnames(FunctionScoreMatrix)

FunctionScoreMatrix <- FunctionScoreMatrix[,orderC]
my.breaks <- c(seq(-1, 0, by=0.1), seq(0.1, 1, by=0.1))
my.colors <- c(
  colorRampPalette(colors = c("#377EB8", "white"))(length(my.breaks)/2),
  colorRampPalette(colors = c("white", "#E41A1C"))(length(my.breaks)/2))
## cellType_col <- data.frame(cell.type = CD8_Obj_CellType)
## rownames(cellType_col) <- colnames(FunctionScoreMatrix)
signatureType_row <- data.frame(Signature.type = c(
  rep("Differentiation", length(Differentiation)),
  rep("Function", length(Function)),
  rep("Metabolism", length(Metabolism)),
  rep("Apoptosis", length(Apoptosis))))
rownames(signatureType_row) <- MarkerNameVector


#6DCCFD #FD9AA0

library(pheatmap)

dev.off()
pdf( file = "T_heatmap.pdf", width = 7, height = 8)
pheatmap::pheatmap(FunctionScoreMatrix,
         show_colnames = T,
         show_rownames = T,
         ## annotation_col = cellType_col,
         annotation_row = signatureType_row,
         gaps_row = c(3, 14, 17),
         cluster_rows = F,
         cluster_cols = F,
         breaks = my.breaks,
         color = my.colors,
         cellwidth = 12, cellheight = 12,
         legend_position = "left",
         angle_col=90,
         border_color = "black",
         fontsize = 8,
         width = 5,
         silent = TRUE,
         height = 3.8)
dev.off()


colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)



##CD4
CD8_data <- read_csv("CD4.csv")
marker.list = list()
for (col_name in colnames(CD8_data)) {
  genes <- CD8_data %>% pull(!!sym(col_name)) %>% na.omit()
  marker.list[[col_name]] <- genes
}

marker.list <- marker.list[1:17]

table(sc_pbmc_int2$cell_type)
# 子集化包含感兴趣的 CD8 细胞亚群
CD8_Obj <- subset(sc_pbmc_int2, cell_type %in% c('1:CD4_C1_LEF1', '2:CD4_C2_IL7R'))

# 去除没有数据的细胞类型
CD8_Obj$cell_type <- droplevels(CD8_Obj$cell_type)

# 检查每个细胞类型的细胞数
table(CD8_Obj$cell_type)



CD8_Obj <-  AddModuleScore(CD8_Obj,
                           features = marker.list,
                           ctrl = 5,
                           name = "FunctionScore")
for(i in 1:length(marker.list)){
  colnames(CD8_Obj@meta.data)[colnames(CD8_Obj@meta.data) == paste0("FunctionScore", i)] <- names(marker.list)[i]
}
Idents(CD8_Obj) <- CD8_Obj$seurat_clusters
colnames(CD8_data)

Differentiation <- c("Naive", "Activation/Effector function", "Exhaustion")
Function <- c( "TCR signaling", "Cytotoxicity" ,"Cytokine/Cytokine receptor","Chemokine/Chemokine receptor", "Stress response",             
               "Adhesion", "IFN response" ,"Treg signature" ,"Costimulatory molecules" )
Metabolism <- c("Oxidative phosphorylation", "Glycolysis", "Lipid metabolism")
Apoptosis <- c("Pro-apoptosis", "Anti-apoptosis")
MarkerNameVector <- c(Differentiation, Function, Metabolism, Apoptosis)
FunctionScoreMatrix <- matrix(0,
                              ncol = length(unique(CD8_Obj$seurat_clusters)),
                              nrow = length(marker.list))
colnames(FunctionScoreMatrix) <- paste0(rownames(table(CD8_Obj$cell_type)))
rownames(FunctionScoreMatrix) <- MarkerNameVector
for(ci in 1:ncol(FunctionScoreMatrix)){
  for(ri in 1:nrow(FunctionScoreMatrix)){
    FunctionVec <- as_tibble(CD8_Obj@meta.data) %>% pull(MarkerNameVector[ri])
    fv <- mean(FunctionVec[CD8_Obj$cell_type == levels(CD8_Obj$cell_type)[ci]])
    FunctionScoreMatrix[ri, ci] <- fv
  }
}

library(scales)
FunctionScoreMatrix <- t(apply(FunctionScoreMatrix, 1, rescale, to=c(-1, 1)))
orderC = colnames(FunctionScoreMatrix)

FunctionScoreMatrix <- FunctionScoreMatrix[,orderC]
my.breaks <- c(seq(-0.5, 0, by=0.1), seq(0.1, 0.5, by=0.1))
my.colors <- c(
  colorRampPalette(colors = c("#377EB8", "white"))(length(my.breaks)/2),
  colorRampPalette(colors = c("white", "#E41A1C"))(length(my.breaks)/2))
## cellType_col <- data.frame(cell.type = CD8_Obj_CellType)
## rownames(cellType_col) <- colnames(FunctionScoreMatrix)
signatureType_row <- data.frame(Signature.type = c(
  rep("Differentiation", length(Differentiation)),
  rep("Function", length(Function)),
  rep("Metabolism", length(Metabolism)),
  rep("Apoptosis", length(Apoptosis))))
rownames(signatureType_row) <- MarkerNameVector


#6DCCFD #FD9AA0

library(pheatmap)

dev.off()
pdf( file = "T_heatmap.pdf", width = 7, height = 8)
pheatmap::pheatmap(FunctionScoreMatrix,
                   show_colnames = T,
                   show_rownames = T,
                   ## annotation_col = cellType_col,
                   annotation_row = signatureType_row,
                   gaps_row = c(3, 14, 17),
                   cluster_rows = F,
                   cluster_cols = F,
                   breaks = my.breaks,
                   color = my.colors,
                   cellwidth = 12, cellheight = 12,
                   angle_col=315,
                   border_color = "black",
                   fontsize = 8,
                   width = 5,
                   height = 3.8)
dev.off()











###基因集打分小提琴图
CD8_data <- read_csv("CD4.csv")
marker.list = list()
for (col_name in colnames(CD8_data)) {
  genes <- CD8_data %>% pull(!!sym(col_name)) %>% na.omit()
  marker.list[[col_name]] <- genes
}

# marker.list <- marker.list[1:17]

table(sc_pbmc_int2$cell_type)
# 子集化包含感兴趣的 CD8 细胞亚群
CD8_Obj <- subset(sc_pbmc_int2, cell_type %in% c('1:CD4_C1_LEF1'))

# 去除没有数据的细胞类型
CD8_Obj$cell_type <- droplevels(CD8_Obj$cell_type)

# 检查每个细胞类型的细胞数
table(CD8_Obj$cell_type)



CD8_Obj <-  AddModuleScore(CD8_Obj,
                           features = marker.list,
                           ctrl = 5,
                           name = "FunctionScore")
for(i in 1:length(marker.list)){
  colnames(CD8_Obj@meta.data)[colnames(CD8_Obj@meta.data) == paste0("FunctionScore", i)] <- names(marker.list)[i]
}
colnames(CD8_Obj@meta.data)

CD8_Obj <- subset(CD8_Obj, group %in% c('LI_NAR','PB_NAR'))

# 最关键的作图
Idents(CD8_Obj) <- CD8_Obj$group
my_comparisons <- list( c("LI_NAR", "PB_NAR"))

VlnPlot(CD8_Obj, features= "Chemokine/Chemokine receptor" , group.by = "group", cols = c('#20466D', '#4C8562'), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.size = 0)+
  scale_y_continuous(limits = c(-0.3,1))+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = " Chemokine/Chemokine receptor", x="") + theme(legend.position="right") +  
  # stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  scale_fill_manual(values = c('#20466D', '#4C8562'), labels = c("LIVER", "PBMC"))+
  stat_compare_means(comparisons = my_comparisons,label="p.signif",size=5,method = 'wilcox.test',show.legend = T,label.x=1.5)
  
  ggsave(filename = 'Activation_Effector function.pdf',width = 8,height = 6,dpi = 300)
  
  
  
  
  
  
  
  # 给定的列表，在基因列表后面加上后缀
  marker_list <- c('3:CD8_C1_GZMK','4:CD8_C2_RGS1','5:CD8_C3_GZMH','6:CD8_C4_LEF1',
                   '7:CD8_C5_KLRB1','8:CD8_C6_TCF7','9:CD8_C7_FOSB','10:CD8_C8_MKI67',
                   '11:CD8_C9_TIGIT')
  
  # 生成含有LI_AR和LI_NAR变体的引用列表
  formatted_markers <- c()
  for (marker in marker_list) {
    formatted_marker_AR <- paste0(marker, '_LI_AR')
    formatted_marker_NAR <- paste0(marker, '_LI_NAR')
    formatted_markers <- c(formatted_markers, formatted_marker_AR, formatted_marker_NAR)
  }
  
  # 输出结果
  formatted_markers
