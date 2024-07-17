# devtools::install_github("YosefLab/VISION@v2.1.0") 
# devtools::install_github("wu-yc/scMetabolism")


library(scMetabolism)
library(tidyverse)
library(rsvd)
library(Seurat)
library(pheatmap)
library(ComplexHeatmap)
library(ggsci)

table(sc_pbmc_int2$cell_type, sc_pbmc_int2$group)
sc_NK = subset(sc_pbmc_int2, cell_type %in% c('1:DC_C1_CD1C',"15:Mac_C3_CXCL12",'17:Mac_C5_TIMD4') )
table(sc_NK$cell_type, sc_NK$group)

sc_NK$cell_type <- droplevels(sc_NK$cell_type)


Idents(sc_NK) <- "cell_type"
countexp.Seurat <- sc.metabolism.Seurat(obj = sc_NK,  #Seuratde单细胞object
                                        method = "AUCell", 
                                        imputation = F, 
                                        ncores = 10, 
                                        metabolism.type = "KEGG")

saveRDS(countexp.Seurat,file = 'NKscMetabolish_reactome.rds')



##将metabolism赋值到metadata里面
#提取score结果
score <- countexp.Seurat@assays$METABOLISM$score
score[1:4,1:4]
#将score中barcode的点转为下划线
score_change <- score %>% 
  select_all(~str_replace_all(., "\\.", "-"))  #基因ID不规范会报错,下划线替换-
#确定细胞barcode椅子
identical(colnames(score_change) , rownames(countexp.Seurat@meta.data))
#[1] TRUE
countexp.Seurat@meta.data <- cbind(countexp.Seurat@meta.data,t(score_change) )



#可以直接使用Seurat的相关函数
p1 <- FeaturePlot(countexp.Seurat,features = "Bile acid and bile salt metabolism",reduction = 'tsne',split.by = 'group')
p1
p2 <- VlnPlot(countexp.Seurat,features = "Primary bile acid biosynthesis")
p1 + p2


VlnPlot( countexp.Seurat, features = 'Primary bile acid biosynthesis', 
         cols = c("olivedrab", "salmon"),
         pt.size = 0,
         group.by = "cell_type",
         split.plot = TRUE,
         split.by = "group" 
)

rownames(countexp.Seurat@assays$METABOLISM$score)


#######常规作图

DimPlot.metabolism(obj = countexp.Seurat, 
                   pathway = "Bile acid and bile salt metabolism", 
                   dimention.reduction.type = "tsne", 
                   dimention.reduction.run = F, size = 0.5)



table(countexp.Seurat@assays$METABOLISM)

input.pathway<-c("Primary bile acid biosynthesis","Oxidative phosphorylation", 
                 "Citrate cycle (TCA cycle)")

input.pathway <- rownames(countexp.Seurat@assays$METABOLISM$score)[1:25]

DotPlot.metabolism(obj = countexp.Seurat, 
                   pathway = input.pathway, 
                   phenotype = c("group"),  #更改phenotype 参数
                   norm = "y")

BoxPlot.metabolism(obj = countexp.Seurat, 
                   pathway = "Primary bile acid biosynthesis",
                     # input.pathway[1:4], 
                   phenotype = "group", 
                   ncol = 2) +
  scale_fill_nejm()
ggsave(filename = 'NK_boxplot_kegg.pdf',width = 15,height = 10)
ggsave(filename = 'NK_boxplot_reactome.pdf',width = 15,height = 10)

######组间比较
features = colnames(countexp.Seurat@meta.data[38:48])
col <- pal_npg("nrc")(4)
DotPlot( object = countexp.Seurat, 
         group.by = "cell_type",
         split.by = 'group',
         cols = col,
         features = features )+
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))





######热图
#可以计算celltype均值，然后绘制
colnames(countexp.Seurat@meta.data)
df <- countexp.Seurat@meta.data

df$new_column <- paste(df$cell_type, df$group, sep="")


#19列开始是代谢通路的得分，按照celltype计算均值
avg_df = aggregate(df[,36:ncol(df)],
                   list(df$new_column),
                   mean)
library(dplyr)

avg_df <- as.data.frame(avg_df)
class(df)
#热图需要转为矩阵
avg_df <- avg_df %>% 
  select(1:87) %>% #展示前20个
  column_to_rownames("Group.1") 
avg_df[1:4,1:4]

rownames <-  rownames(avg_df)
# 假设你的行名存储在一个字符向量中
rownames <- c("13:Mac_C1_CCL3LI_AR", "13:Mac_C1_CCL3LI_NAR", "13:Mac_C1_CCL3PB_AR", "13:Mac_C1_CCL3PB_NAR",
              "5:Mono_C1_VCANLI_AR", "5:Mono_C1_VCANLI_NAR", "5:Mono_C1_VCANPB_AR", "5:Mono_C1_VCANPB_NAR",
              "6:Mono_C2_S100A12LI_AR", "6:Mono_C2_S100A12LI_NAR", "6:Mono_C2_S100A12PB_AR", "6:Mono_C2_S100A12PB_NAR",
              "8:Mono_C4_ATF3LI_AR", "8:Mono_C4_ATF3LI_NAR", "8:Mono_C4_ATF3PB_AR", "8:Mono_C4_ATF3PB_NAR")



# # 定义自定义排序函数
# custom_order <- function(x) {
#   # 提取后缀并按照要求排序
#   suffix <- gsub(".*_(\\w+)$", "\\1", x)
#   order_index <- case_when(
#     grepl("LI_AR$", x) ~ 1,
#     grepl("LI_NAR$", x) ~ 2,
#     grepl("PB_AR$", x) ~ 3,
#     grepl("PB_NAR$", x) ~ 4,
#     TRUE ~ 5  # 处理未匹配到的情况
#   )
# 
#   # 返回排序后的索引
#   return(order_index)
# }
# 
# # 根据自定义排序函数对行名进行排序
# sorted_rownames <- rownames[order(custom_order(rownames))]
# 
# # 输出排序后的行名
# sorted_rownames



# 定义自定义排序函数
custom_order <- function(x) {
  # 定义排序顺序
  order <- c("LI_AR", "LI_NAR", "PB_AR", "PB_NAR")

  # 提取后缀并映射到排序顺序
  suffix <- gsub(".*_(\\w+)$", "\\1", x)
  order_index <- match(suffix, order)

  # 返回排序后的索引
  return(order_index)
}

# 根据自定义排序函数对行名进行排序
sorted_rownames <- rownames[order(custom_order(rownames))]

# 输出排序后的行名
sorted_rownames

avg_df_sorted <- avg_df[match(sorted_rownames, rownames(avg_df)),]

avg_df_sorted$new_column <- NULL

# 提取行名包含 LI_AR 和 LI_NAR 的行
li_ar_nar_indices <- grep("LI_AR|LI_NAR", rownames(avg_df_sorted))
li_ar_nar_df <- avg_df_sorted[li_ar_nar_indices, ]

# 提取行名包含 PB_AR 和 PB_NAR 的行
pb_ar_nar_indices <- grep("PB_AR|PB_NAR", rownames(avg_df_sorted))
pb_ar_nar_df <- avg_df_sorted[pb_ar_nar_indices, ]


library(SCP)


dev.off()
pdf( file = "DM_meta_PBMC_heatmap_kegg2.pdf", width = 10, height = 30 )
# pdf( file = "NK_Metabolish_reactome.pdf", width = 8, height = 10 )

pheatmap(t(pb_ar_nar_df), 
         show_colnames = T,
         scale='row', 
         cluster_rows = T,
         angle_col='45',
         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
         cluster_cols = F)
dev.off()





######其他热图
dev.off()
pdf( file = "NK_heatmap_kegg.pdf", width = 12, height = 14 )
exp <- apply(avg_df, 2, scale)
rownames(exp) <- rownames(avg_df)
# 组件
h_state <- Heatmap(t(exp),
                   column_title = "state_gsva",
                   col = colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
                   name= "gsva ",
                   show_row_names = TRUE,
                   column_names_rot=30,
                   show_column_names = TRUE)

h_state
dev.off()


########自定义作图



#提取坐标（umap图）
umap.loc<-countexp.Seurat@reductions$umap@cell.embeddings

row.names(umap.loc)<-colnames(countexp.Seurat)
signature_exp<-countexp.Seurat@assays$METABOLISM$score

input.pathway <- rownames(countexp.Seurat@assays$METABOLISM$score)[1:5]
input.pathway <- input.pathway

signature_ggplot<-data.frame(umap.loc, t(signature_exp[input.pathway,]),
                             countexp.Seurat@active.ident,countexp.Seurat@meta.data$group)

colnames(signature_ggplot)
library(wesanderson)
pal <- wes_palette("Zissou1", 100, type = "continuous")

####R包自带作图
library(ggplot2)
plot <- ggplot(data=signature_ggplot, aes(x=umap_1, y=umap_2, color = signature_ggplot[,1])) +  #this plot is great
  geom_point(size = 1) +
  scale_fill_gradientn(colours = pal) +
  scale_color_gradientn(colours = pal) +
  labs(color = input.pathway) +
  #xlim(0, 2)+ ylim(0, 2)+
  xlab("UMAP 1") +ylab("UMAP 2") +
  theme(aspect.ratio=1)+
  #theme_bw()
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_wrap(~ countexp.Seurat.meta.data.group)
plot

#####nature medicine 作图

#添加标签
class_avg <- signature_ggplot %>%
  group_by(countexp.Seurat.active.ident) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2)
  )
colnames(input.pathway)
colnames(signature_ggplot)
ggplot(data=signature_ggplot, aes(x=umap_1, y=umap_2))  +
  labs(color = "Fructose.and.mannose.metabolism" ) +
  # ggrepel::geom_label_repel(aes(label = countexp.Seurat.active.ident),
  #                           data = class_avg,
  #                           label.size = 1,
  #                           segment.color = NA)+
  geom_point(aes(colour  = signature_ggplot[,3])) + viridis::scale_color_viridis(option="A") +
  theme(legend.position = "none") + theme_bw()+ theme (panel.grid=element_blank ())
  # + facet_wrap(~ countexp.Seurat.meta.data.group)

#######分组dotplot图
DotPlot.metabolism_chen <- function(obj, input.pathway, phenotype, group,norm = "y"){
   input.norm = norm
  input.pathway <- input.pathway
  input.parameter<-phenotype
  input.group<-group  ###添加分组信息
  
  
  metadata<-countexp.Seurat@meta.data
  metabolism.matrix <- countexp.Seurat@assays$METABOLISM$score
  
  cat("\nPlease Cite: \nYingcheng Wu, Qiang Gao, et al. Cancer Discovery. 2021. \nhttps://pubmed.ncbi.nlm.nih.gov/34417225/   \n\n")
  
  
  metadata[,input.parameter]<-as.character(metadata[,input.parameter])
  metadata[,input.group]<-as.character(metadata[,input.group])
  metabolism.matrix_sub<-t(metabolism.matrix[input.pathway,])
  
  #arrange large table
  gg_table<-c()
  for (i in 1:length(input.pathway)){
    gg_table<-rbind(gg_table, cbind(metadata[,input.parameter], input.pathway[i], metabolism.matrix_sub[,i], metadata[,input.group]))
  }
  gg_table<-data.frame(gg_table)
  
  
  ####添加分组信息
  gg_table[5] = paste(gg_table[,1], gg_table[,4], sep = '.')
  temp <- gg_table[1]
  gg_table[1] <- gg_table[5]
  gg_table[5] <- temp
  
  
  #get median value
  gg_table_median<-c()
  input.group.x<-unique(as.character(gg_table[,1]))
  input.group.y<-unique(as.character(gg_table[,2]))
  
  
  for (x in 1:length(input.group.x)){
    for (y in 1:length(input.group.y)){
      gg_table_sub<-subset(gg_table, gg_table[,1] == input.group.x[x] & gg_table[,2] == input.group.y[y])
      gg_table_median<-rbind(gg_table_median, cbind(input.group.x[x], input.group.y[y], median(as.numeric(as.character(gg_table_sub[,3])))))
      
    }
  }
  gg_table_median<-data.frame(gg_table_median)
  gg_table_median[,3]<-as.numeric(as.character(gg_table_median[,3]))
  
  
  #normalize
  gg_table_median_norm<-c()
  input.group.x<-unique(as.character(gg_table[,1]))
  input.group.y<-unique(as.character(gg_table[,2]))
  
  
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
  if (input.norm == "y")
    for (y in 1:length(input.group.y)){
      gg_table_median_sub<-subset(gg_table_median, gg_table_median[,2] == input.group.y[y])
      norm_value<- range01(as.numeric(as.character(gg_table_median_sub[,3])))
      gg_table_median_sub[,3]<-norm_value
      gg_table_median_norm<-rbind(gg_table_median_norm, gg_table_median_sub)
    }
  
  if (input.norm == "x")
    for (x in 1:length(input.group.x)){
      gg_table_median_sub<-subset(gg_table_median, gg_table_median[,1] == input.group.x[x])
      norm_value<- range01(as.numeric(as.character(gg_table_median_sub[,3])))
      gg_table_median_sub[,3]<-norm_value
      gg_table_median_norm<-rbind(gg_table_median_norm, gg_table_median_sub)
    }
  
  if (input.norm == "na") gg_table_median_norm<-gg_table_median
  
  
  gg_table_median_norm<-data.frame(gg_table_median_norm)
  gg_table_median_norm[,3]<-as.numeric(as.character(gg_table_median_norm[,3]))
  
  
  
  
  library(wesanderson)
  pal <- wes_palette("Zissou1", 100, type = "continuous")
  
  ggplot(data=gg_table_median_norm, aes(x=gg_table_median_norm[,1], y=gg_table_median_norm[,2], color = gg_table_median_norm[,3])) +
    geom_point(data=gg_table_median_norm, aes(size = gg_table_median_norm[,3])) + #geom_line() +
    #theme_bw()+theme(aspect.ratio=0.5, axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("Metabolic Pathway")+ xlab(input.parameter)+
    theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1), #aspect.ratio=1,
                     panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    scale_color_gradientn(colours = pal) +
    labs(color = "Value", size = "Value") +
    #facet_wrap(~tissueunique, ncol = 1) +
    #theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    NULL
}


input.pathway <- rownames(countexp.Seurat@assays$METABOLISM$score)[1:15]
DotPlot.metabolism_chen(obj = countexp.Seurat, 
                   input.pathway = input.pathway, 
                   phenotype = "cell_type",  #更改phenotype 参数
                   group='group',
                   norm = "y")
ggsave(filename = 'NK_Dotplot_kegg.pdf',width = 8,height = 10)
ggsave(filename = 'NK_Dotplot_reactome.pdf',width = 8,height = 10)
