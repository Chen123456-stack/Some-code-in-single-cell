library(SeuratData)
library(Seurat)
library(ggplot2)
library(ggSCvis)
library(dplyr)
library(scRNAtoolVis)


#####dotplot#####
#查看细胞数量
Idents(sc_pbmc_int2)
table(sc_pbmc_int2@meta.data$cell_type)
unique(sc_pbmc_int2@meta.data$cell_type)
#细胞排序
sc_pbmc_int2@meta.data$cell_type <- factor(sc_pbmc_int2@meta.data$cell_type,
                                           levels = c( "T",
                                                       'NK',
                                                       "B",
                                                       "Plasma cell", 
                                                       "Monocyte",
                                                       'Macrophage',
                                                       "cDC2",
                                                       "cDC1",
                                                       "pDC",
                                                       "Neutrophil",
                                                       'Basophil',
                                                       "Magakaryocyte",
                                                       "Endothelial cell",
                                                       "Hepatocyte",
                                                       'Fibroblast'
                                           ))
# sc_pbmc_int2@meta.data$cell_type <- factor(sc_pbmc_int2@meta.data$cell_type,
#                                            levels = c(
#                                              'CD4Tn',
#                                              'CD4Tm',
#                                              'CD4Treg',
#                                              'CD8Tn',
#                                              'CD8Te',
#                                              'CD8Tex',
#                                              'CD8Tm',
#                                              'CD8pro',
#                                              'CD8MAIT',
#                                              'DN',
#                                              'NK_cycle',
#                                              'NK_rest'))

#获取topmarkers
top10 = cell.markers %>% # 挑选每个群差异靠前的基因，非必要，可以用别的方法筛选或指定
  filter( avg_log2FC > 0 ) %>% # 是否只看上调的，都可以，看自己的需求
  group_by(cluster) %>% # 按cluster分组
  filter( p_val_adj < 0.05 ) %>% # 按校正P值筛选，看自己的需求
  top_n( 10, avg_log2FC )%>% # 提取各cluster种logFC最大的3个，看前多少就设置多少
  group_by() # 取消分组

top5 = cell.markers %>% # 挑选每个群差异靠前的基因，非必要，可以用别的方法筛选或指定
  filter( avg_log2FC > 0 ) %>% # 是否只看上调的，都可以，看自己的需求
  group_by(cluster) %>% # 按cluster分组
  filter( p_val_adj < 0.05 ) %>% # 按校正P值筛选，看自己的需求
  top_n( 5, avg_log2FC )%>% # 提取各cluster种logFC最大的3个，看前多少就设置多少
  group_by()  # 取消分组

top3 = cell.markers %>% # 挑选每个群差异靠前的基因，非必要，可以用别的方法筛选或指定
  filter( avg_log2FC > 0 ) %>% # 是否只看上调的，都可以，看自己的需求
  group_by(cluster) %>% # 按cluster分组
  filter( p_val_adj < 0.05 ) %>% # 按校正P值筛选，看自己的需求
  top_n( 3, avg_log2FC ) %>% # 提取各cluster种logFC最大的3个，看前多少就设置多少
  group_by() # 取消分组

#jjDotplot作图
jjDotPlot(object = sc_pbmc_int2,
          gene = top5$gene,
          id = 'cell_type',
          ytree = F)

ggsave(filename = 'Dotplot_pbmc.pdf',width = 12,height = 8)

#批量火山图

#提取设置好细胞类型的颜色
colour=c("#F08080","#1E90FF","#7CFC00", "#808000","#FF00FF","#FA8072"
         ,"#7B68EE","#800080","#A0522D","#D2B48C","#D2691E")

#可自己设置基因
mygene <- c('CCL3','CCL5','CXCL9')
#可top基因
jjVolcano(diffData = final_result,#批量差异基因
          log2FC.cutoff = 0.25,
          col.type = "updown",
          topGeneN = 5,
          #myMarkers = mygene,
          legend.position=c(0.9,0.9),
          # tile.col=colour)
          tile.col = corrplot::COL2('RdBu', 20)[5:16])

ggsave(filename = 'Volcano_pbmc.pdf',width = 12,height = 8)



FeaturePlot(sc_pbmc_int2, features = c( "CD3E", "KLRF1", "CD79A",'LYZ','FCGR3B','ALB'),ncol = 3,
            reduction = 'umap',cols = c("grey", "#4C8562"))#min.cutoff阈值，只有表达高于此细胞群体10%的表达值的细胞才会显示

# FeaturePlot作图，cluster marker基因
FeaturePlot(sc_pbmc_int2, features = c( "CD3E",'CD8A', "KLRF1", "CD79A",'LYZ','FCGR3B','CD1C','CLEC9A','ALB'),
            reduction = 'umap',cols = c("grey", "#4C8562"))#min.cutoff阈值，只有表达高于此细胞群体10%的表达值的细胞才会显示
ggsave(filename = 'featureplot.pdf',width = 14,height = 8)

FeaturePlot(sc_pbmc_int2, features = c( "CD4", "CD8A", "NCAM1",'FCGR3A'),
            reduction = 'umap',raster=FALSE,cols = c("grey", "red"),min.cutoff = "q2")#min.cutoff阈值，只有表达高于此细胞群体10%的表达值的细胞才会显示
ggsave(filename = 'featureplot.pdf',width = 14,height = 12)