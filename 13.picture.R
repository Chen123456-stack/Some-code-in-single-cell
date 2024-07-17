#一些作图的标准及参数设置
setwd('~/picture')
#umap图
p3 = UMAPPlot( sc_pbmc_int2,
               label =T,pt.size=0.8,
               group.by = 'seurat_clusters')+
  ggtitle("UMAP")+
  theme_bw( )+
  theme(
    plot.title = element_text( face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank()   # 去掉次要网格线
  ) +
  scale_color_manual(values = colorRampPalette(pal_igv("default")(50))(46)); p3

# table(sc_pbmc_int2$group,sc_pbmc_int2$seurat_clusters)
# table(sc_pbmc_int2$group,sc_pbmc_int2$cell_type)

p4 = UMAPPlot( sc_pbmc_int2,
               label =T,pt.size=0.8,
               group.by = 'cell_type')+
  ggtitle("UMAP")+
  theme_bw( )+
  theme(
    plot.title = element_text( face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank()   # 去掉次要网格线
  ) + scale_color_manual(values = colorRampPalette(pal_igv("default")(30))(25)); p4

p3 + p4

ggsave( filename = "5.6-细胞分群注释umap.pdf", width = 18, height = 8 )

#type
p = UMAPPlot( sc_pbmc_int2,
               label =F,
               order=F,shuffle = T,
               group.by = 'type')+
  ggtitle("UMAP")+
  theme_bw( )+
  theme(
    plot.title = element_text( face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank()   # 去掉次要网格线
  ) + scale_color_manual(values = colorRampPalette(pal_igv("default")(30))(25)); p

ggsave( filename = "t.pdf", width = 9, height = 8 ,dpi = 300)
































