

CD8_Obj <- sc_pbmc_int2
CD8_Obj$cell_type <- droplevels(CD8_Obj$cell_type)
table(CD8_Obj$group,CD8_Obj$cell_type)



Idents(CD8_Obj) <- CD8_Obj$cell_type
cell.markers <- FindAllMarkers(object = CD8_Obj, 
                               only.pos = T, # 是否只保留表达相对上调的基因，设置FALSE则会保留下调的
                               test.use = "wilcox", # 默认使用 wilcox 非参数检验，其它选项可以查看说明
                               slot = "data", # 需要注意的是，默认使用 data，而不是 counts
                               min.pct = 0.25, # 设置表达比例的阈值，没有统一标准，差异基因很多的情况下可以把阈值调高，差异基因有1000个够用
                               logfc.threshold = 0.25 # 设置 log2FC 即差异倍率的阈值，没有统一标准，差异基因很多的情况下可以把阈值调高
)


# 假设你感兴趣的亚群是 Cluster 1
cluster1_markers <- cell.markers %>% filter(cluster == "13:Mac_C2_RGS1")

# 找到在 Cluster 1 中高表达，而在其他亚群中表达很低的基因
# 这里我们假设其他亚群中的表达值低于 min.pct（最小表达比例阈值）
exclusive_cluster1_markers <- cluster1_markers %>%
  rowwise() %>%
  filter(all(avg_log2FC > 0.25 & pct.1 > 0.25 & pct.2 < 0.1))

# 查看结果
print(exclusive_cluster1_markers)


saveRDS( cell.markers, file = "M_cellMarkers.rds" )

library( xlsx )
write.xlsx(cell.markers, file = "M_CellMarkers.xlsx" )


colorsForDataType <- c("#6DCCDD", "#EDCAE0", "#F494BE", "#F9B26C", "#A6ADCC", "#C4DA5D")

markers <- c("GZMK", "GZMB",'GZMA', "PRF1",  
             'RGS1',"FAS", "FASLG", "PDCD1","LAG3", "CTLA4",
             "FGFBP2", "GZMH", "GNLY",'CX3CR1','FCGR3A',
             "CCR7", "SELL", "TCF7",
             'KLRB1','IL7R','KLRG1','CD160',
             'FOSB','HSPA1A','HSPA1B','IFNG',
             'TOP2A','MKI67','TYMS','MCM2',
             'GPR183','MYB','CD28','CD27',
             "CD44", "CD69",'EOMES','ITGA1','IL7R','CXCR4')

markers <- c('CCL1','CCL2',"CCL3", "CCL3L3",'CCL4', "CCL5",  'CCL7','CCL8','CCL18',
             'CXCL1','CXCL2','CXCL3','CXCL6','CXCL8','CXCL9','CXCL10')


markers <- c('CCR1','CCR2',"CCR3", "CCR4",'CCR5', "CCR6",  'CXCR1','CXCR2','CXCR3',
             'CXCR4','CXCR5','CXCR6','CXCR8','CXCR9','CXCR10')


markers <- c('CCR7','SELL', 'LEF1', 'TCF7',
            'CD69','CD38','NKG7','KLRB1','FCGR3A','CX3CR1','FGFBP2','JUN','FOS','FOSB','NR4A1','EOMES',
            'CD3D','CD3E','CD3G','CD247','NFATC2','NFATC1','NFATC3','NFATC3','IRF4',
             'GZMA','GZMB','GZMH','GZMK','GNLY','PRF1','TNF','TNFSF9','IFNG','IFI44L',
             'CCR1','CCR4', 'CCR5', 'CCR7', 'CXCR3', 'CXCR4','CXCR5','CXCR6','RGS1',
             'CCL3','CCL3L3','CCL4','CCL4L1','CCL4L2','CCL5', 'CXCL13',
            'MKI67','TOP2A',
            'PDCD1','HAVCR2','CTLA4','TIGIT')


markers <- c("C1QA","C1QB","C1QC",'CD68','TIMD4','CD163','CD5L','MARCO','VSIG4','MRC1',
             "VCAN","FCN1",'S100A8','S100A9','S100A10',"S100A12",
             'CD14','FCGR3A','RGS1','CD69','CXCR4',
             'HLA-A','HLA-B','HLA-C','HLA-E','HLA-F',
             'HLA-DMA','HLA-DMB','HLA-DOA','HLA-DPA1','HLA-DPB1','HLA-DQA1','HLA-DQA2','HLA-DQB1','HLA-DRA','HLA-DRB1','HLA-DRB5','HLA-DRB6'
             )


markers <- c('CCL3','CCL3L3','CCL4','CCL4L2','CCL5', 'XCL1','XCL2','IFNG','IL16','IL32','IL18','IL18R1','KLRC2',
             'S1PR1', 'CXCR3','CXCR4','CXCR6','CX3CR1','GZMA','GZMB','GZMK','PRF1','NKG7','CD160','CD69','ITGAE','ITGA1','RGS1','NCAM1',
             'FCGR3A','B3GAT1','CREM','NR4A1','NR4A2','NR4A3','HSPA1A','DNAJB1','NFKBIA')

markers <- c('CCL3','CCL3L3','CCL4','CCL4L2','CCL5', 'XCL1','XCL2','IFNG','HSPA1A','CD160','CD69','RGS1')

#TIGIT,TIM3
#'CCL3','CCL3L3','CCL4','CCL4L2','CCL2','CXCL2','CXCL10','CXCL12','CXCL16',
#'CCR1','CCR2','CCR5',,'CXCR6'
#'IFI27','IFIT1','IFITM2','IFITM3','IFI6','IFIT3','IFIT2','IFI44'

markers <- c('CXCR1','CXCR2','CXCR3','CXCR4','CXCR5','CXCR6','CCR1','CCR2','CCR3','CCR4',
             'CCR6','CCR7','CCR8','CCR9','CCR10','XCR1','XCR2')


markers <- c(
'GZMA','GZMB','GZMM','GZMK','GZMH','PRF1','GNLY','TNFSF10',#cytotoxic
  'KLRK1','KLRC2','NCR1','CD160','NCR3','SLAMF6','SLAMF7', #activite
'KLRC1','CD300A','TIGIT','KIR3DL2','KIR3DL1','KIR2DL3','KIR2DL1','KLRB1','SIGLEC7','HAVCR2'#inhibit
)

markers <- c('FCGR3A','NCAM1')

### 'HBA1','CD3G','CD3D','DSTN','ZBTB38','CD2','PPDPF','PRDM1','LGALS1','S100A4','ITGB1','IL32','VIM','LINCO1871','S100A6','PTMS','GZMH','CD3E','CCL5','KLRC2'

markers <- c('CCL3')
table(t$cell_type)


CD8_Obj <- sc_pbmc_int2

table(sc_pbmc_int2$cell_type,sc_pbmc_int2$group)


t <- CD8_Obj

CD8_Obj <- subset(CD8_Obj,group%in%c('LI_AR','LI_NAR'))

CD8_Obj <- subset(CD8_Obj,group%in%c('LI_AR','LI_NAR')&cell_type%in%c( '3:CD8_C1_GZMK',
                                                                       '4:CD8_C2_RGS1',
                                                                       '5:CD8_C3_GZMH',
                                                                       '6:CD8_C4_LEF1',
                                                                       '7:CD8_C5_KLRB1',
                                                                       '8:CD8_C6_TCF7',
                                                                       '9:CD8_C7_FOSB',
                                                                       '10:CD8_C8_MKI67',
                                                                       '11:CD8_C9_TIGIT'))

table(CD8_Obj$cell_type,CD8_Obj$group)

CD8_Obj$cell_type <- droplevels(CD8_Obj$cell_type)

metadata <- CD8_Obj@meta.data

# 创建一个新的列，将 cell_type 和 group 组合起来
metadata$cell_group <- paste(metadata$subcell_type2, metadata$group, sep = "_")

# 将新的元数据添加回 Seurat 对象
CD8_Obj <- AddMetaData(CD8_Obj, metadata = metadata)

table(CD8_Obj$cell_group)
# 
# CD8_Obj <- scNK
Idents(CD8_Obj) <- CD8_Obj$cell_group
gene <- intersect(markers, rownames(CD8_Obj))
p <- DotPlot(CD8_Obj, features =gene)+ coord_flip()+
  scale_fill_gradient2(low =  "#377EB8", high ="#E41A1C")
p
data <- p$data[,c('id','features.plot','pct.exp','avg.exp.scaled')]
table(data$id)


data$id <- factor(data$id,levels = rev(c('1:CD4_C1_LEF1_LI_AR','1:CD4_C1_LEF1_LI_NAR','2:CD4_C2_IL7R_LI_AR','2:CD4_C2_IL7R_LI_NAR',
                                         "3:CD8_C1_GZMK_LI_AR"  ,  "3:CD8_C1_GZMK_LI_NAR" ,  "4:CD8_C2_RGS1_LI_AR"  ,  "4:CD8_C2_RGS1_LI_NAR" ,  "5:CD8_C3_GZMH_LI_AR",
                                         "5:CD8_C3_GZMH_LI_NAR" ,  "6:CD8_C4_LEF1_LI_AR" ,   "6:CD8_C4_LEF1_LI_NAR" ,  "7:CD8_C5_KLRB1_LI_AR" ,  "7:CD8_C5_KLRB1_LI_NAR" ,
                                        "8:CD8_C6_TCF7_LI_AR" ,   "8:CD8_C6_TCF7_LI_NAR" ,  "9:CD8_C7_FOSB_LI_AR" ,   "9:CD8_C7_FOSB_LI_NAR"  , "10:CD8_C8_MKI67_LI_AR" ,
                                          "10:CD8_C8_MKI67_LI_NAR" ,"11:CD8_C9_TIGIT_LI_AR",  "11:CD8_C9_TIGIT_LI_NAR",
                                          '12:NK_C1_FCGR3A_LI_AR','12:NK_C1_FCGR3A_LI_NAR','13:NK_C2_XCL1_LI_NAR','13:NK_C2_XCL1_LI_AR',
                                        '14:NK_C3_KLRF1_LI_NAR','14:NK_C3_KLRF1_LI_AR','15:DN_LINC-PINT_LI_AR','15:DN_LINC-PINT_LI_NAR')))


###对坐标进行排序，按顺序进行排序
x <- data$id

# 将数据转换为数据框，并拆分成数字和后缀部分
df <- data.frame(original = x) %>%
  mutate(number = as.numeric(sub(":.*", "", original)),
         suffix = ifelse(grepl("LI_AR", original), "LI_AR", "LI_NAR"))

# 按照数字和后缀排序
sorted_df <- df %>%
  arrange(number, suffix)

# 获取排序后的数据
sorted_data <- sorted_df$original

sorted_data

data$id <- factor(data$id, levels=rev(c(sorted_data)))

# 
# data$id <- factor(data$id, levels = c('6:Mono_C2_S100A12_PB_AR','7:Mono_C3_FCGR3A_PB_AR','8:Mono_C4_ATF3_PB_AR','9:Mono_C5_FCN1_PB_AR',
#                                       '10:Mono_C6_LRMDA_PB_AR','11:Mono_C7_RPS13_PB_AR','12:Mono_C8_RTF2_PB_AR','13:Mac_C1_CCL3_PB_AR',
#                                       '6:Mono_C2_S100A12_PB_NAR','7:Mono_C3_FCGR3A_PB_NAR','8:Mono_C4_ATF3_PB_NAR','9:Mono_C5_FCN1_PB_NAR',
#                                       '10:Mono_C6_LRMDA_PB_NAR','11:Mono_C7_RPS13_PB_NAR','12:Mono_C8_RTF2_PB_NAR','13:Mac_C1_CCL3_PB_NAR'))
# 
# labels <- c('6:Mono_C2_S100A12_PB_AR','7:Mono_C3_FCGR3A_PB_AR','8:Mono_C4_ATF3_PB_AR','9:Mono_C5_FCN1_PB_AR',
#             '10:Mono_C6_LRMDA_PB_AR','11:Mono_C7_RPS13_PB_AR','12:Mono_C8_RTF2_PB_AR','13:Mac_C1_CCL3_PB_AR',
#             '6:Mono_C2_S100A12_PB_NAR','7:Mono_C3_FCGR3A_PB_NAR','8:Mono_C4_ATF3_PB_NAR','9:Mono_C5_FCN1_PB_NAR',
#             '10:Mono_C6_LRMDA_PB_NAR','11:Mono_C7_RPS13_PB_NAR','12:Mono_C8_RTF2_PB_NAR','13:Mac_C1_CCL3_PB_NAR')
# 
# # 替换标签中的 "PB" 为 "LI"
# labels_modified <- gsub("PB", "LI", labels)

labels <- c('14:NK_C3_KLRF1_LI_NAR','14:NK_C3_KLRF1_LI_AR','13:NK_C2_XCL1_LI_NAR','13:NK_C2_XCL1_LI_AR','12:NK_C1_FCGR3A_LI_NAR','12:NK_C1_FCGR3A_LI_AR')


labels <- c("3:CD8_C1_GZMK_LI_AR"  ,  "3:CD8_C1_GZMK_LI_NAR" ,  "4:CD8_C2_RGS1_LI_AR"  ,  "4:CD8_C2_RGS1_LI_NAR" ,  "5:CD8_C3_GZMH_LI_AR",
                                                           "5:CD8_C3_GZMH_LI_NAR" ,  "6:CD8_C4_LEF1_LI_AR" ,   "6:CD8_C4_LEF1_LI_NAR" ,  "7:CD8_C5_KLRB1_LI_AR" ,  "7:CD8_C5_KLRB1_LI_NAR" ,
                                                           "8:CD8_C6_TCF7_LI_AR" ,   "8:CD8_C6_TCF7_LI_NAR" ,  "9:CD8_C7_FOSB_LI_AR" ,   "9:CD8_C7_FOSB_LI_NAR"  , "10:CD8_C8_MKI67_LI_AR" ,
                                                           "10:CD8_C8_MKI67_LI_NAR" ,"11:CD8_C9_TIGIT_LI_AR",  "11:CD8_C9_TIGIT_LI_NAR")
#替换标签中的 "PB" 为 "LI"
labels_modified <- gsub("LI", "PB", labels)
# 
data$id <- factor(data$id, levels = rev(labels_modified))




# #反向排序id
data <- data %>%
  mutate(id = factor(id, levels = rev(unique(id))))

plotx <- ggplot(data, aes(x = id, y = features.plot)) +        ## global aes
  geom_point(aes(fill = avg.exp.scaled, size = pct.exp),
             color = 'black',
             shape = 21,
             stroke = 0.01)  +    ## geom_point for circle illusion
  #scale_fill_gradientn(colours=rev(color),limits=c(0,max(data$avg.exp)))+       ## color of the corresponding aes
  # scale_x_discrete(breaks=0:13, labels=paste0("CD8_c", 0:13)) +
  xlab("") + ylab("") +
  scale_fill_gradient2(low =  "#377EB8", high ="#E41A1C")+
  #scale_fill_gradient2(high = colorsForDataType[3], mid = "white", low = "#6DCCFF")+
  # scale_fill_gradientn(
  #   ## limits = c(-1.5, 1.5),
  #   colors = c("#5DBCFF", "#6DCCFF", "white", colorsForDataType[3], "#F484AE"))+
  scale_size(range = c(0, 6), limits = c(0, 100), breaks = c(0,20,40,60,80,100))+   
  theme_bw()+
  # coord_flip() +
  theme(
    text = element_text(size = 8),
    panel.grid.major = element_line(colour = "grey90", size=0.2),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text.x=element_text(angle = 315, vjust = 2, hjust = 0,size=10),
    axis.text.y=element_text(angle =0, vjust = 1, hjust = 0,size=10),
    legend.position="right",
    legend.title=element_text(size=10)) +
  guides(size = guide_legend(title.position="top",
                             title.hjust = 0.5,
                             ncol = 1,
                             byrow = T,
                             override.aes = list(stroke = 0.4)),
         fill = guide_colourbar(title.position = "top", title.hjust = 0.5))

plotx


VlnPlot(CD8_Obj,features = c('CCL3','CCL5'),pt.size = 0,split.by = 'group')

ggsave( "PBMC_bubbleplot.pdf", plotx, width =1000/3.4,dpi = 600,
        height =160/1.7, units = "mm")

#常用参数
ggsave( "CCL3_bubbleplot.pdf", plotx, width =1200/3.4,dpi = 600,
        height = 400/1.7, units = "mm")

ggsave( "PB.pdf", plotx, width =1200/3.4,dpi = 600,
        height = 400/1.7, units = "mm")

ggsave( "CCL3_bubbleplot2.pdf", plotx, width =550/3.4,dpi = 600,
        height = 400/1.7, units = "mm")

ggsave( "CCL3_bubbleplot2.pdf", plotx, width =800,dpi = 600,
       height = 300, units = "mm")


ggsave( "CCL3_M.pdf", plotx, width =2000/3.4,dpi = 600,
        height = 100/1.7, units = "mm")
