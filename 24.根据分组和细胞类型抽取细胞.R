#加载R包
.libPaths(c('/mnt/hpc/home/lvjunqiang/miniconda3/envs/R/lib/R/library',
            "/mnt/hpc/home/lvjunqiang/R/x86_64-pc-linux-gnu-library/4.3",
            "/opt/R/4.3.1/lib/R/library" ))

setwd("/mnt/hpc/home/lvjunqiang/M")
library(Seurat)
####增加线程
library( future )
availableCores()
nbrOfWorkers()

plan("multisession", workers = 20 ) ###增加线程，冲冲冲！！！
plan()
#没有安装包的需要先安装这些包
library(VGAM)
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

#BiocManager::install("monocle")
#install.packages("./monocle_2.24.0.tar.gz", repos = NULL, type = "source")
# 加载Biobase包
library(monocle)#没有安装包的需要先安装这些包
sc_pbmc_int2 <- readRDS('M.rds')

# setwd('3.monocleT')

set.seed(12345)

table(sc_pbmc_int2$cell_type,sc_pbmc_int2$group)
seurat=subset(sc_pbmc_int2,cell_type%in%c("5:Mono_C1_VCAN", "6:Mono_C2_S100A12", "7:Mono_C3_FCGR3A",
                                          "8:Mono_C4_ATF3", "9:Mono_C5_FCN1", "10:Mono_C6_LRMDA", "11:Mono_C7_RPS13", 
                                          "12:Mono_C8_RTF2","13:Mac_C1_CCL3",'14:Mac_C2_C1QA','15:Mac_C3_CXCL12',
                                          '16:Mac_C4_CD5L','17:Mac_C5_TIMD4'))

table(seurat$cell_type)
seurat$cell_type <- droplevels(seurat$cell_type)
table(seurat$cell_type)

# 根据细胞比例和分组信息选择50,000个细胞的名称


Idents(seurat) <- seurat$cell_type

# 获取细胞类型和分组信息
cell_types <- seurat@meta.data$cell_type
groups <- seurat@meta.data$group

# 创建一个包含细胞类型和分组的联合标识符
combined_labels <- interaction(cell_types, groups, sep = "_")

# 获取每种组合的细胞数量
combined_counts <- table(combined_labels)


# 计算每种组合需要抽取的数量
total_cells <- 20000
combined_proportions <- combined_counts / sum(combined_counts)
cells_to_sample <- round(combined_proportions * total_cells)

# 打印每种组合需要抽取的数量
print(cells_to_sample)

# 计算每种组合需要抽取的数量
total_cells <- 20000
combined_proportions <- combined_counts / sum(combined_counts)

cells_to_sample <- round(combined_proportions * total_cells)
table(cells_to_sample)
# 获取UMAP坐标
umap_coords <- Embeddings(seurat, reduction = "umap")

# 从每种组合中根据UMAP坐标距离中心点选取细胞
selected_cells <- unlist(lapply(names(cells_to_sample), function(combined_label) {
  split_label <- strsplit(combined_label, "_")
  cell_type <- paste(split_label[[1]][1:(length(split_label[[1]])-2)], collapse="_")
  group <- paste(split_label[[1]][(length(split_label[[1]])-1):length(split_label[[1]])], collapse="_")

  
  print(paste("split_label:", split_label))
  print(paste("cell_type:", cell_type))
  print(paste("group:", group))
  
  # 使用WhichCells函数来筛选细胞
  cells_in_combination <- WhichCells(seurat, idents = cell_type)
  cells_in_group <- cells_in_combination[seurat@meta.data[cells_in_combination, "group"] == group]
  
  if (length(cells_in_group) < 2) {
    warning(paste("Not enough cells for cell type", cell_type, "and group", group))
    return(NULL)
  }
  
  # 计算每种组合细胞的 UMAP 坐标中心点
  center <- colMeans(umap_coords[cells_in_group, , drop = FALSE])

  
  # 计算每个细胞到中心点的距离
  distances <- apply(umap_coords[cells_in_group, ], 1, function(coord) {
    sum((coord - center)^2)
  })
  
  # 根据距离对细胞排序
  ordered_cells <- cells_in_group[order(distances)]
  
  # 选择前n个细胞，如果组合中的细胞数量不足，则选择所有细胞
  cells_to_select <- min(length(ordered_cells), cells_to_sample[combined_label])
  ordered_cells[1:cells_to_select]
}))

# 检查抽取的细胞总数
length(selected_cells)

# 创建一个包含这些细胞的子集Seurat对象
seurat_sub <- subset(seurat, cells = selected_cells)

table(seurat_sub$cell_type,seurat_sub$group)

dim(seurat_sub)

table(seurat$cell_type,seurat$group)

p1 <- DimPlot(seurat,reduction = 'umap',group.by = 'cell_type')
p2 <- DimPlot(seurat_sub,reduction = 'umap',group.by = 'cell_type')

p1+p2
seurat <- seurat_sub

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


#====================method3====================================================
#differentialGeneTest
cds_DGT <- cds
diff_test_res <- differentialGeneTest(cds_DGT,fullModelFormulaStr = "~cell_type",cores=10)
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
#选择或者不选择
# ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:3000]#或者选择前多少的基因
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)
save(cds, file = "cds_DGT1.RData")

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



