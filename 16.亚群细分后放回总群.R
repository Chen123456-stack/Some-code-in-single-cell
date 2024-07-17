
library(dplyr)


sc_pbmc_int2@meta.data$cellid <- rownames(sc_pbmc_int2@meta.data)
M@meta.data$cellid <- rownames(M@meta.data)
M$new_celltype <- M$cell_type



colnames(M@meta.data)
colnames(sc_pbmc_int2@meta.data)

metadata<- left_join(x=sc_pbmc_int2@meta.data[,c(36,37)],y= M@meta.data[,c(24,36,37)], by= "cellid")

metadata$subcell_type2 <- ifelse(is.na(metadata$cell_type), 
                                as.character(metadata$subcell_type),
                                as.character(metadata$new_celltype))


sc_pbmc_int2@meta.data$subcell_type2 <- metadata$subcell_type2


Idents(sc_pbmc_int2) <- "subcell_type2"

p4 = UMAPPlot( sc_pbmc_int2, alpha=0.8,
               label =F,cols = zzm60colors,repel=T,
               group.by = 'subcell_type2')+
  ggtitle("UMAP")+
  # theme_bw( )+
  theme(
    plot.title = element_text( face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank()   # 去掉次要网格线
  ) ; p4

dir.create('./cellchat')
setwd('./cellchat/')

saveRDS(sc_pbmc_int2,'cellchat.rds')

table(sc_pbmc_int2$group,sc_pbmc_int2$subcell_type2)


table(sc_pbmc_int2$subcell_type2)
sc_pbmc_int2 <- sc_pbmc_int2[, !(sc_pbmc_int2$subcell_type2 %in% c('T','NK'))]

colnames(sc_pbmc_int2@meta.data)
colnames(kup@meta.data)

table(sc_pbmc_int2$cell_type)

##后续进一步调整分群

sc_pbmc_int2 <- FindNeighbors(sc_pbmc_int2, reduction = "harmony", dims = 1:30) 
sc_pbmc_int2 <- FindClusters(sc_pbmc_int2,resolution = 0.1)
sc_pbmc_int2 <- RunUMAP(sc_pbmc_int2, dims = 1:30,reduction = "harmony", n.neighbors =30L,
                        umap.method = "umap-learn",metric = "correlation",min.dist = 1) 

sc_pbmc_int2 <- RunTSNE(sc_pbmc_int2, dims = 1:25, method = "FIt-SNE", reduction = "harmony", 
                        perplexity = 30, fitsnePath = "/home/data/t010560/FItSNE_local/FIt-SNE")


p4 = TSNEPlot( sc_pbmc_int2, alpha=0.8,
               label =T,cols = mycolor,repel=T,raster.dpi=c(512,512),order=T,shuffle = FALSE,
               group.by = 'cell_type')+
  ggtitle("TSNE")+
  # theme_bw( )+
  theme(
    plot.title = element_text( face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank()   # 去掉次要网格线
  ) ; p4

p4 = UMAPPlot( sc_pbmc_int2, alpha=0.8,
               label =T,cols = mycolor,repel=T,raster.dpi=c(1024,1024),order=F,shuffle = FALSE,
               group.by = 'subcell_type')+
  ggtitle("UMAP")+
  # theme_bw( )+
  theme(
    plot.title = element_text( face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank()   # 去掉次要网格线
  ) ; p4

