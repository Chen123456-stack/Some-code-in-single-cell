# 上皮细胞拟时序分析
.libPaths(c('~/R/x86_64-pc-linux-gnu-library/4.3',
            '/refdir/Rlib',
            '/usr/local/lib/R/library'))

setwd('~/T集中抽样2w 可用/change/')
library(Seurat)
#没有安装包的需要先安装这些包
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(reshape2)
library(Biobase)
library(ggsci)
library(ggpubr)
library(data.table)
library(viridis)
library(igraph)
library(Matrix)
library(VGAM)
#BiocManager::install("monocle")
#install.packages("./monocle_2.24.0.tar.gz", repos = NULL, type = "source")
# 加载Biobase包
library(monocle)#没有安装包的需要先安装这些包

setwd('./T集中抽样2w 可用/change/')


dir.create('./2.monocle')
setwd('2.monocle')

set.seed(12345)


####创建monocle对象####
# 01选择需要构建细胞分化轨迹的细胞类型（subset提取感兴趣的名字，上面已经准备好名字了）
table(sc_pbmc_int2$cell_type,sc_pbmc_int2$group)
seurat=subset(sc_pbmc_int2,cell_type%in%c("CD8_C1_RGS1",'CD8_C2_GPR183','CD8_C3_CX3CR1',
                                          'CD8_C4_GZMK','CD8_C5_LEF1','CD8_C6_MKI67'))
# seurat <- sc_pbmc_int2
unique(seurat$cell_type)
colnames(seurat@meta.data)
table(seurat$cell_type,seurat$group)
# 02构建monocle独有的celldata对象
#count矩阵
data=LayerData(seurat, assay = 'RNA', layer = 'counts')#使用counts表达值(第一个准备文件：基因count矩阵)
#metadata
pd <- new("AnnotatedDataFrame", data = seurat@meta.data[,c(1,8,9,24)])
#gene信息
fData=data.frame(gene_short_name=rownames(data),row.names = row.names(data))#构建一个含有基因名字的数据框
fd <- new("AnnotatedDataFrame", data = fData)
#构建对象
cds <- newCellDataSet(#as.matrix(data),
                      as(as.matrix(data),'sparseMatrix'),#若数据太大转换为稀疏矩阵
                      phenoData = pd, 
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily=negbinomial.size())#创建一个monocle的对象
cds #cellData对象；monocle独有




# 03相当于归一化Size factors帮助标准化基因在细胞间的差异, "dispersion"值帮助后面的差异表达分析执行。
cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds,cores=10)



# 04数据降维，特征基因的选择
# 过滤低表达基因；差异分析获得top1000排序基因（可视化用来排序的基因）
###注意--differentialGeneTest这一步有些久
cds <- detectGenes(cds, min_expr = 0.1)#统计，过滤低表达基因（表达量过滤低表达基因）
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >= 10))#（数量过滤低表达基因）
print(head(pData(cds)))

######monocle差异分析选择拟时序相关基因####
#====================method1====================================================
#monocle计算高变基因
cds_monocle <- cds
disp_table <- dispersionTable(cds_monocle)
unsup_clustering_genes <- subset(disp_table, 
                                 mean_expression >= 0.01 
                                 & dispersion_empirical >= 1 * dispersion_fit)[1:2000,]#meanexpression筛选标准按自己情况
head(unsup_clustering_genes)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
plot_ordering_genes(cds)
write.csv(unsup_clustering_genes,'degForCellOrdering.csv')
ggsave(filename = 'monocle2_ordering_gene.pdf')
save(cds, file = 'cds_monocle.RData')
#====================method2====================================================
#使用seurat计算
cds_seurat <- cds
seurat_variable_genes <- VariableFeatures(FindVariableFeatures(sc_pbmc_int2, nfeatures = 4000,
                                                               assay = "RNA"), assay = "RNA")
cds <- setOrderingFilter(cds, seurat_variable_genes)
plot_ordering_genes(cds)
save(cds, file = 'cds_seurat.RData')


#====================method3====================================================
#differentialGeneTest
cds_DGT <- cds
diff_test_res <- differentialGeneTest(cds_DGT,fullModelFormulaStr = "~cell_type",cores=10)
ordering_genes <- row.names (subset(diff_test_res, qval < 0.1))
#选择或者不选择
# ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:3000]#或者选择前多少的基因
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)
save(cds, file = "cds_DGT.RData")



#################################################################################


# 降维；绘制细胞轨迹（绘制state与细胞类型的轨迹图）
cds <- reduceDimension(cds, reduce_method = 'DDRTree',
                       # max_components = 2,
                       # residualModelFormulaStr = "~orig.ident",
                       verbose = T)#降维

cds <- orderCells(cds)#排序

# 后面根据推测选择起点
# cds <- orderCells(cds,root_state =3)


####运行完orderCells，保存数据####
save.image(file = "orderCells_cds_done.RData")  # 保存所有对象到my_workspace.RData


save.image(file = "after_beam.RData")  # 保存所有对象到my_workspace.RData

# load("orderCells_cds_done.RData")

p2=plot_cell_trajectory(cds, color_by = "group") +#可视化细胞类型
  theme(text = element_text(size = 18))  # 设置字体大小为 18
p2


p1=plot_cell_trajectory(cds, color_by = "State") +#可视化鉴定的State
  theme(text = element_text(size = 18))  # 设置字体大小为 18
p1
ggsave(p1,filename = 'monocle2_state_trajectory.pdf', width = 10, height = 6.5)

p2=plot_cell_trajectory(cds, color_by = "cell_type", use_color_gradient=F) +#可视化细胞类型
  theme(text = element_text(size = 18))  # 设置字体大小为 18
p2
ggsave(p2,filename = 'monocle2_celltype_trajectory.pdf', width = 10, height = 6.5)

#4.按发育时间分类画轨迹图，可以与state的图对应上
p3=plot_cell_trajectory(cds, color_by = "Pseudotime")+
  theme(text = element_text(size = 18))  # 设置字体大小为 18
p3
ggsave(filename = 'monocle2_Pseudotime.pdf',  width = 10, height = 6.5) 

p4=plot_cell_trajectory(cds, color_by = "orig.ident")+
  theme(text = element_text(size = 18))  # 设置字体大小为 18
p4
ggsave(filename = 'monocle2_orig.ident.pdf', width = 10, height = 6.5) 

p5=plot_complex_cell_trajectory(cds,x=1,y=2,color_by = 'cell_type')
p5

p6=ggplot(pData(cds),aes(Pseudotime,colour = cds$cell_type,fill=cds$cell_type))+
            geom_density(bw=0.5,size =1,alpha = 0.5)+theme_classic()
p6
ggsave(filename = 'monocle2_density.pdf', width = 10, height = 6.5) 

#标记想观察基因
genes= c(ordering_genes)[1:4]
p1 = plot_genes_in_pseudotime(cds[genes],color_by="State")
p2 = plot_genes_in_pseudotime(cds[genes],color_by="cell_type")
p3 = plot_genes_in_pseudotime(cds[genes],color_by="Pseudotime")
ggsave("7Trajectory_Pseudotime.pdf",plot=p1|p2|p3, width=10, height=6.5)

p1 = plot_genes_jitter(cds[genes],grouping="State",color_by="State")
p2 = plot_genes_violin(cds[genes],grouping="State",color_by="State")
p3 = plot_genes_in_pseudotime(cds[genes],color_by="State")
ggsave("8Trajectory_jitter.pdf",plot=p1|p2|p3, width = 10, height = 6.5)

#单独基因的轨迹表达情况
pData(cds)$GZMK =log2(exprs(cds)["GZMK",] +1)
p1 = plot_cell_trajectory(cds,color_by='GZMK') +
  scale_color_continuous(type ="viridis")
p1

pData(cds)$CXCL9 =log2(exprs(cds)["CXCL9",] +1)
p2 = plot_cell_trajectory(cds,color_by='CXCL9') +
   scale_color_continuous(type ="viridis")
p2

pData(cds)$MARCO =log2(exprs(cds)["MARCO",] +1)
p3 = plot_cell_trajectory(cds,color_by='MARCO') +
  scale_color_continuous(type ="viridis")
p3

pData(cds)$MARCO =log2(exprs(cds)["MARCO",] +1)
p4 = plot_cell_trajectory(cds,color_by='MARCO') +
  scale_color_continuous(type ="viridis")
p4
ggsave("9Trajectory_Expression.pdf",plot = p1|p2|p3|p4, width=14,height=6.5)


#指定细胞分化轨迹的起始点，即设置root cell,这个是指定细胞类型。
#自定义了一个名为GM_state的函数，函数主要计算包含特定类型细胞最多的state，并返回state编号
#即特定细胞类型最多的state为指定起点
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$cell_type)[,"Mono_C1_CD14"]#细胞类型填写这里即可
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  } else {return (1)}
}
#有时间点采样，就将imm_cell_type改成hours；0时刻为起始状态
cds <- orderCells(cds, root_state = GM_state(cds)) #指定分化起始点，排序细胞


#其他绘图。。。等等要明白理解的
#把每个分支单独画出来
p1 <- plot_cell_trajectory(cds, color_by = "cell_type")+
  facet_wrap(~group, nrow = 2)#facet_wrap函数可以把每个state单独highlight
p1
ggsave(filename = 'monocle2_facet_celltype.pdf',width = 10, height = 6.5) 

# 按celltype单独画出来
p2=plot_cell_trajectory(cds, color_by = "cell_type")+
  facet_wrap(~cell_type, nrow = 5) +
  theme(text = element_text(size = 25))  # 设置字体大小为 18
p2
ggsave(file='monocle2_celltype分开.pdf', width = 12, height = 18)

p2=plot_cell_trajectory(cds, color_by = "group")+
  facet_wrap(~group, nrow = 5) +
  theme(text = element_text(size = 25))  # 设置字体大小为 18
p2

## 指定特定基因在各个state的表达(可辅助判断细胞起点)
#注意：可以可视化一些增殖基因在哪些state高表达（一般干性细胞高表达或细胞周期marker）
#然后可以假设其为分化起点。（当然也可以看自身感兴趣的）
blast_genes <- row.names(subset(fData(cds),gene_short_name %in% c("CCNB2", "CCNE1", "CCND1")))#获取特定基因向量
plot_genes_jitter(cds[blast_genes,],grouping = "State",min_expr = 0.1)
ggsave(filename = 'monocle2_cycling_gene_state.pdf')
blast_genes <- row.names(subset(fData(cds),gene_short_name %in% c("CCL5", "TCEA1", "RB1CC1")))#获取特定基因向量
plot_genes_jitter(cds[blast_genes,],grouping = "State",min_expr = 0.1)#可视化基因表达，散点图
ggsave(filename = 'monocle2_specific_gene_state.pdf') 
plot_genes_jitter(cds[blast_genes,],grouping = "cell_type",min_expr = 0.1)#可视化基因表达，散点图
ggsave(filename = 'monocle2_specific_gene_celltype.pdf') 



# 3.鉴定伪时间相关的基因，即分化过程相关 ----------------------------------------------------
#鉴定伪时间相关的基因，即分化过程（state）相关的差异基因（可以指定计算个数，20可为100）
##根据Pseudotime找差异基因
diff_test_res=differentialGeneTest(cds[ordering_genes,],
                                   fullModelFormulaStr = "~sm.ns(Pseudotime)",cores=10)#鉴定伪时间相关的基因，即分化过程相关
write.csv(diff_test_res,file='diff_test_res.csv') 





##根据State找差异基因
diff_test_res=differentialGeneTest(cds[ordering_genes,],
                                   fullModelFormulaStr = "~sm.ns(State)",cores=10)#鉴定伪时间相关的基因，即分化过程相关
write.csv(diff_test_res,file='diff_test_res.csv')

save(diff_test_res,BEAM_res,file = 'res.RData')

# 
# 
# ## 注意可以将top3g改成任何自身感兴趣的基因，进行可视化
# top5g=rownames(diff_test_res[order(diff_test_res$qval),])[1:5] # top 5个显著基因
# plot_genes_in_pseudotime(cds[top5g,], color_by="cell_type", ncol = 1) +# top 5个显著基因，可视化，pseudotime作为横坐标
#   theme(text = element_text(size = 15))  # 设置字体大小为 18
# ggsave(filename = 'monocle2_significant_gene.pdf') 

####伪时间相关基因热图####
#伪时间基因排序
diff_res = diff_test_res[order(diff_test_res$qval),]
diff_res = diff_res[1:10,]

class(rownames(diff_res))

gc()

dev.off()
pdf( file = "T_heatmap.pdf", width = 7, height = 12)
p <- plot_pseudotime_heatmap(cds[rownames(diff_res),], #挑选的20个计算了分化差异的基因
                        num_clusters =4,
                        cores = 10, #指定线程数
                        show_rownames = T,
                        cluster_rows = T,
                        return_heatmap = T,
                        use_gene_short_name = T)#对伪时间相关的基因绘制热图，查看具备相似表达模式的基因簇
table(p$tree_row$order,p$tree_row$labels)

p
dev.off()


gene <- c('CCL3',	'CCL3L3',	'CCL4',	'CCL4L2',	'CCL5',	'CD69',	'CCR7',	'CD2',	'CD38',	'CD74',	'CD8A',	'CREM',	'CTLA4',	'CXCR4',
          'CXCR6',	'FOS'	,'FOSB'	,'GNLY',	'GZMA',	'HLA-DQA1',	'HLA-DRA',	'HMGB2'	,'IFI16',	'IFNG',	'IL10',	'IL2RB'	,	'JUN',
          'LEF1',	'MKI67',	'NR4A1',	'NR4A3',	'PDCD1',	'RGS1',	'RGS2',	'SELL',	'TCF7',	'TIGIT',	'TOP2A',	'TOX'	,'TOX2',	'GZMK',	'HLA-DRB1',
          'CCR5',	'CXCL13',	'PIM2',	'CD27',	'GZMB',	'CXCR2',	'HAVCR2',	'BTLA',	'LAG3',	'EOMES'	,'HOPX',
          'ARNT',	'ETV1',	'IRF8',	'KLRB1',	'FOXP1')

gene <- c('S1PR1',  'CD74', 'HLA-DQA1', 'CXCR4', 'FOS', 'GZMB', 'IFNG', 'PDCD1', 'GZMK','CCL3','CCL3L3','CD69',
          'FCMR','TOX', 'MKI67', 'HLA-DRA', 'CD2', 'CD38', 'CD8A', 'CREM', 'CTLA4', 'CXCR6', 'GZMA', 'CCR7',
          'PCED1B', 'HLA-DQA1', 'HLA-DRB1', 'IFI16', 'IL10', 'IL2RB', 'IL7R', 'JUN', 'LEF1', 'MKI67', 'NR4A1', 'FOSB'	,
          'SNX9','NR4A3', 'RGS1', 'RGS2', 'SELL', 'TCF7', 'TIGIT', 'TOP2A', 'TOX2', 'GZMK', 'CCR5', 'PRF1',
          'EEF1G','CXCL13', 'S100A9', 'PIM2', 'CD27', 'CXCR2', 'HAVCR2', 'BTLA', 'CCR6','KLF10','DNAJC2','BCL11B','ELK3','NFKB1',
          'PLAC8','LAG3', 'EOMES', 'HOPX', 'ARNT', 'ETV1', 'IRF8', 'FOXP1','LEF1',)


gene <- c("IL7R" , "CCR7" , "SELL" , "FOXO1", "KLF2" , "KLF3" , "LEF1" , "TCF7" , "ACTN1", "FOXP1",##naive
          "FAS"  ,  "FASLG" , "CD44" ,  "CD69" ,  "CD38"  , "NKG7"  , "KLRB1" , "KLRD1" , "KLRF1" , "KLRG1" ,#effect
          "KLRK1" , "FCGR3A", "CX3CR1", "CD300A" ,"FGFBP2", "ID2"  ,  "ID3"   , "PRDM1" , "RUNX3" , "TBX21" ,
          "ZEB2" ,  "BATF"  , "IRF4"  , "NR4A1" , "NR4A2" , "NR4A3" , "PBX3" ,  "ZNF683", "HOPX"  , "FOS"  , 
          "FOSB"  , "JUN"   , "JUNB"  , "JUND"  , "STAT1",  "STAT2" , "STAT5A", "STAT6" , "STAT4",  "EOMES" ,
          "GZMA" , "GZMB"  ,   "GZMH"  ,   "GZMK"  ,   "GNLY"   ,  "PRF1"  ,   "IFNG"  ,   "TNF"   ,#Cytotoxicity
          "SERPINB1", "SERPINB6", "SERPINB9", "CTSA"  ,   "CTSB"  ,   "CTSC"  ,   "CTSD"  ,   "CTSW"  ,   "CST3"  ,
          "CST7"  ,   "CSTB"  ,   "LAMP1" ,   "LAMP3"   , "CAPN2", 
          "CCR4"  , "CCR5"  ,  "CXCR3" , "CXCR4",  "CXCR5" , "CXCR6"  ,"CCL3" ,  "CCL4"  , "CCL4L1","CCL4L2","CCL5" ,#Chemokine
          "CXCL13" ,"CXCL8" , "XCL1" ,"XCL2",
          "NT5E"  ,  "IZUMO1R", "NRP1"  ,  "DGKA"  ,  "CBLB"  ,  "RNF128" , "ITCH"  ,  "NFATC2" , "EGR2"  ,#Anergy
          "EGR3" ,     "TOB1",
          "PDCD1" , "LAYN" ,  "HAVCR2", "LAG3"  , "CD244",  "CTLA4"  ,"LILRB1", "TIGIT" , "TOX"  ,  "VSIR"  , #exhust
          "BTLA" ,  "ENTPD1", "CD160" , "LAIR1" )

duplicated_genes <- gene[duplicated(gene)]
duplicated_genes



source('./add.flag.R')
add.flag(p,kept.labels = gene,repel.degree = 0.2)


####获取gene，重新画图
clusters <- cutree(p$tree_row, k = 4)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)

clustering$gene <- rownames(clustering)

custom_order <- c(2, 1, 3, 4)

# 2. 使用 match 函数生成排序索引
order_genecluster <- match(clustering$Gene_Clusters, custom_order)

# 3. 根据排序索引重新排列数据框
clustering_sorted <- clustering[order(order_genecluster), ]

# 打印排序后的数据框
print(clustering_sorted$gene)



gene <- c("CCR7", "SELL", "LEF1", "TCF7", "ACTN1", "FOXP1", "NT5E", "DGKA", "RNF128", "IL7R", "KLF2", 
           "KLF3", "CD44", "KLRB1", "FCGR3A", "CX3CR1", "CD300A", "FGFBP2", "ZNF683", "STAT6", "LILRB1", "VSIR", 
           "LAIR1", "GZMH", "GNLY", "CTSA", "CTSW", "CST3", "LAMP1", "CAPN2", "CXCR5", "FOXO1", 
           "FASLG", "CD38", "KLRK1", "PRDM1", "IRF4", "NR4A2", "NR4A3", "PBX3", "JUND", "STAT1", "STAT2", 
          "PDCD1", "LAYN", "HAVCR2", "LAG3", "CTLA4", "TIGIT", "TOX", "BTLA", "IFNG", 
           "SERPINB9", "CTSD", "CXCR6", "CCL3", "CCL4L2", "CXCL13", "CXCL8", "CBLB", "EGR2", "TOB1", "FAS", 
           "CD69", "NKG7", "KLRD1", "KLRF1", "KLRG1", "ID2", "ID3", "RUNX3", "TBX21", "ZEB2", "BATF", "NR4A1", 
           "HOPX", "FOS", "FOSB", "JUN", "JUNB", "STAT4", "EOMES", "CD244", "CD160", "GZMA", "GZMB", "GZMK", 
           "PRF1", "TNF", "SERPINB1", "CTSB", "CTSC", "CST7", "CSTB", "CCR4", "CCR5", "CXCR3", "CXCR4", "CCL4", 
           "CCL5", "XCL1", "XCL2", "NRP1", "NFATC2", "EGR3",'MKI67','TOP2A','RGS1')








#更改颜色
plot_pseudotime_heatmap(cds[rownames(diff_res),],
                        num_clusters = 4,
                        cores = 10,
                        show_rownames = T,
                        return_heatmap = T,
                        hmcols = colorRampPalette(viridis(4))(1000))#对伪时间相关的基因绘制热图，查看具备相似表达模式的基因簇


hp.genes <-p$tree_row$labels[p$tree_row$order]
Time_diff_sig <- diff_test_res[hp.genes, c("gene _short_name", "pval", "gval")]
write.csv(Time_diff_sig,"13Time_diff sig.csv",row.names = F)
                             

# 4.鉴定分支依赖基因--------------------------------------------------------------
### 分支比较后，分支依赖基因的鉴定（热图重现）
# 基因簇按照自身需求，设置；图形上侧显示两个分支（朝不同方向分化）

BEAM_res <- BEAM(cds[ordering_genes,], branch_point = 1, cores = 10)#分支依赖基因鉴定，这里只选择了20个基因进行演示

BEAM_res = BEAM_res[order(BEAM_res$qval),]

BEAM_res2 <- subset(BEAM_res, gene_short_name %in% gene)

BEAM_res2 = BEAM_res[1:20,]


#BEAM #鉴定分支依赖的基因
plot_genes_branched_heatmap(cds[rownames(BEAM_res2),],
                            branch_point = 1,
                            num_clusters = 3,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T
)#绘制分支热图
ggsave(filename = 'monocle2_beam_heatmap.pdf') 

### 选取BEAM的三个基因做例证(感兴趣基因可以替换rownames(BEAM_res)[]可视化)
plot_genes_branched_pseudotime(cds[rownames(BEAM_res),][7:10],
                               branch_point = 1,
                               color_by = "cell_type",ncol = 1)#绘制分支pseudotime散点图
ggsave(filename = 'monocle2_beam_gene.pdf') 
plot_genes_branched_pseudotime(cds[rownames(BEAM_res)[1:5],],
                               branch_point = 1,
                               color_by = "State",ncol = 1)#绘制分支pseudotime散点图
gene <- c('GZMB','PRF1','GZMA','NKG7','GNLY','KLRG1')
plot_genes_branched_pseudotime(cds[gene,],
                               branch_point = 1,
                               color_by = "group",ncol = 2)+ scale_fill_manual(values = FJ)+
                              scale_color_manual(values = FJ)+
guides(color = guide_legend(override.aes = list(shape = 15, size = 4)))+
  labs(color = "Celltype")+
  theme_bw()

FJ <- c('#4C8562','#A24F1D','#65A8AE','#363C6C','#3D5C96','#5E5598','#67912F',
        '#D08685','#806224','#A7256F','#834626','#6376B9','#ABCD79','#DBB263',
        '#B5A0C6','#4A277B','#BF9B25','#A72023','#A3BCD3','#C66E26','#528A31',
        '#535353')




#######提取坐标自己画图
plot(cds@reducedDimS)

t <- pData(cds)

plotdf2=as.data.frame(t(cds@reducedDimS))

colnames(plotdf2)=c("component1","component2")

x=combine(t,plotdf2)
colnames(x)

identical(rownames(x), rownames(data1))
t =combine(x,data1)

p <- ggplot(data = x, aes(x = component1, y = component2, color = Naive)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = FJ) +
  scale_fill_manual(values = FJ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    legend.position = 'top',
    legend.text = element_text(size = 8),  # 设置图例文本字体大小
    legend.title = element_text(size = 10)  # 设置图例标题字体大小
  ) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 4)))+
  labs(color = "Celltype")

p
ggsave(filename ='celltype.pdf',width = 10,height = 6 )


setwd('./CD8T umap聚类取细胞2w/change/')

# 使用 viridis 配色方案
#  '#4C8562','#A24F1D'

p <- ggplot(data = t, aes(x = component1, y = component2, color = Exhausted3)) +
  geom_point(size = 1) +
  scale_color_gradientn(colors = colorRampPalette(c( "#F6B04A","#7B2063"))(100)) +
  theme_minimal()+
  theme(
    plot.title = element_text( face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank(),   # 去掉次要网格线
    legend.position = 'top',
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    legend.background = element_rect(fill = "transparent", colour = NA),
    axis.text.x = element_blank(), # 去掉X轴标签
    axis.text.y =element_blank()
  )

p
ggsave(filename ='Pseudotime.pdf',width = 10,height = 6 )





#旋转180°
x_rotated <- x
x_rotated$component1 <- -x$component1
x_rotated$component2 <- -x$component2

p <- ggplot(data = x_rotated, aes(x = component1, y = component2, color = Pseudotime)) +
  geom_point(size = 0.5) +
  scale_color_gradientn(colors = colorRampPalette(c("#F6B04A", "#7B2063"))(100)) +
  theme_minimal() +
  theme(
    plot.title = element_text( face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank(),   # 去掉次要网格线
    legend.position = 'top',
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    legend.background = element_rect(fill = "transparent", colour = NA),
    axis.text.x = element_blank(), # 去掉X轴标签
    axis.text.y =element_blank()
  )

p

ggsave(filename ='Pseudotime.pdf',width = 8,height = 6 )


p <- ggplot(data = x_rotated, aes(x = component1, y = component2, color = cell_type)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = FJ) +
  scale_fill_manual(values = FJ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    legend.position = 'top',
    legend.text = element_text(size = 8),  # 设置图例文本字体大小
    legend.title = element_text(size = 10)  # 设置图例标题字体大小
  ) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 4)))+
  labs(color = "Celltype")

p
ggsave(filename ='celltype.pdf',width = 8,height = 6 )

##提取非排斥患者
x_filtered <- subset(x_rotated, group %in% c("LI_NAR", "PB_NAR"))



p <- ggplot(data = x_rotated, aes(x = component1, y = component2, color = Pseudotime)) +
  geom_point(size = 0.5) +
  scale_color_gradientn(colors = colorRampPalette(c("#F6B04A", "#7B2063"))(100)) +
  theme_minimal() +
  theme(
    plot.title = element_text( face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank(),   # 去掉次要网格线
    legend.position = 'top',
    legend.key = element_rect(fill = "transparent", colour = "transparent"),
    legend.background = element_rect(fill = "transparent", colour = NA),
    axis.text.x = element_blank(), # 去掉X轴标签
    axis.text.y =element_blank()
  )+
  facet_wrap(~group)

p

ggsave(filename ='Pseudotime.pdf',width = 8,height = 6 )


###分组画图
p <- ggplot(data = x_filtered, aes(x = component1, y = component2, color = cell_type)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = FJ) +
  scale_fill_manual(values = FJ) +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 4)))+
  labs(color = "Celltype")+
  facet_wrap(~type)+
  theme_bw()+
  theme(
    axis.text.x = element_blank(), 
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y =element_blank()
  )

p
ggsave(filename ='split_celltype.pdf',width = 12,height = 6 )






