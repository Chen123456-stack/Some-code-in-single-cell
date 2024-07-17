
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(Seurat)
library(scater)
library(RColorBrewer)
library("pheatmap")


#差异基因热图展示

sc_pbmc_int2 <- subset(sc_pbmc_int2, cell_type == "NK_rest")

sce <- as.SingleCellExperiment(sc_pbmc_int2)

# 标记基因
marker_sign <- c('XCL1','XCL2','CCL3','CCL4','CCL4L2','APOA2','APOC3','APOC1','APOH','APOE','APOB','IFNG','CCL18','TNF','CCL3L3','CXCR6',
                 'IL1B','CCR5','CCR1','GZMK','GZMB','CX3CR1','S100A8','CXCR2','S100A9','CXCR1'
)

celltype_mean <- aggregateAcrossCells(as(sce, "SingleCellExperiment"),  
                                      ids = sce$type, 
                                      statistics = "mean",
                                      use.assay.type = "counts", 
                                      subset.row = marker_sign
)
result <- celltype_mean@assays@data@listData$counts
class(result)

col<- colorRampPalette(c("blue", "white", "red"))(256)

#对数据进行改变，使热图更好看
TF_matrix=t(scale(t(result)))#标准化，可以考虑是否使用
TF_matrix <- result
TF_matrix[TF_matrix>3]=3
TF_matrix[TF_matrix< -3] = -3

#热图
p <- pheatmap(TF_matrix, scale='none',cluster_col = FALSE,
              color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
              cluster_row = T,
              cutree_rows =3, cellwidth = 34, cellheight = 14,angle_col=0)
ggsave(p,filename = 'pheatmap.pdf',width = 10,height = 10)

#热图
sscVis::plotMatrix.simple(result,
                          out.prefix=sprintf("%s.group.reference","./Fig1"),
                          show.number=F,
                          waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                          exp.name=expression(italic(AverageExpression)),
                          z.hi=2,
                          palatte=viridis::viridis(7),
                          pdf.width = 6, pdf.height = 10)











