
#===============================================================================
#                             一、软件安装
#===============================================================================
setwd('D:\\KS项目\\公众号文章\\cellchat V2更新教程')
setwd('./cellchat')
#===============================================================================
#安装包并加载单细胞转录组数据
#如果安装出现如下的报错
# ERROR: dependency 'presto' is not available for package 'CellChat'
devtools::install_github("immunogenomics/presto")
#可能会解决你的问题


#安装cellchat V2
devtools::install_github("jinworks/CellChat")
#依赖包要求：
# Install NMF (>= 0.23.0) using install.packages('NMF'). Please check here for other solutions if you encounter any issue. You might can set Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=TRUE) if it throws R version error.
# Install circlize (>= 0.4.12) using devtools::install_github("jokergoo/circlize") if you encounter any issue.
# Install ComplexHeatmap using devtools::install_github("jokergoo/ComplexHeatmap") if you encounter any issue.


library(CellChat)
library(Seurat)

sceV5 <- load("D:/KS项目/公众号文章/cellchat V2更新教程/sceV5.RData")

#===============================================================================
#                             二、基本分析
#===============================================================================
#===============================================================================
#input data- normalized gene exp matix, cell label
#如果你的数据是counts，那么请使用cellchat中的normalizeData函数进行处理
#input数据的构建和之前的cellchat版本一致

#我们的数据分为两组，对照组和实验组,分别提取数据
table(sc_pbmc_int2$group)

# 增加GROUP，将其分为排斥和非排斥组
sc_pbmc_int2$GROUP <- NA

# 将每个样本分配到新的组别
sc_pbmc_int2$GROUP[sc_pbmc_int2$group %in% c("LI_AR", "PB_AR")] <- "AR"
sc_pbmc_int2$GROUP[sc_pbmc_int2$group %in% c("LI_NAR", "PB_NAR")] <- "NAR"

# 检查新的组别
table(sc_pbmc_int2$GROUP)


sceV5 <- sc_pbmc_int2
unique(sceV5$group)
# [1] MDA5 HD 
HD <- subset(sceV5,group%in%c('PB_NAR','LI_NAR'))
MDA <- subset(sceV5,group%in%c('PB_AR','LI_AR'))
# rm(sceV5)

#我们首先分析单个样本、分析HD组的通讯情况
# DefaultAssay(HD) <- 'SCT'
HD_input <- GetAssayData(HD, layer = 'data',assay='RNA')
HD_meta <- HD@meta.data[,c("group","cell_type")]
colnames(HD_meta) <-  c("group","labels")
#保证meta的行名和input matrix的列名是一致的
identical(colnames(HD_input),rownames(HD_meta)) 
# [1] TRUE

#creat cellchat obj
HD.cellchat <- createCellChat(object = HD_input, meta = HD_meta, group.by = "labels")
# cellchat = createCellChat(HD, assay = 'SCT',group.by = 'cell_type')


#Add cell information
# A = createCellChat(object = HD_input)
# cellchat <- addMeta(A, meta = HD_meta)
# cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(HD.cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(HD.cellchat@idents)) # number of cells in each cell group

#===============================================================================
#设置 ligand-receptor interaction 数据库
#在使用CellChat推断细胞之间的通讯之前，需要设置配体-受体相互作用数据库，并识别过表达的配体或受体。
#CellChatDB数据库整理有人和鼠的数据
#CellChatDB v2包含约3300个经过验证的分子相互作用，包括约40%的分泌自分泌/旁分泌信号相互作用，约17%的细胞外基质(ECM)受体相互作用，约13%的细胞-细胞接触相互作用和约30%的非蛋白信号。
#相比于·之前的版本，CellChatDB v2包含有更多的受配体数据
#cellchat的CellChatDB数据库也支持用户自行修改

#我们示例数据是人的，这里我们选择人的数据库
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
# showDatabaseCategory(CellChatDB)

# use a subset of CellChatDB for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
# showDatabaseCategory(CellChatDB.use)

# cellchat obj中设置使用的数据库
HD.cellchat@DB <- CellChatDB

#===============================================================================
#预处理细胞-细胞通讯分析的表达数据
#为了推断细胞状态特异性的通信，CellChat在一个细胞组中识别过表达的配体或受体，然后在配体或受体过表达时识别过表达的配体-受体相互作用。
#subsetData可以选择expression data of signaling genes in CellChatDB.use，
#比如选择"Cell-Cell Contact"等。但是需要注意的是，即使你使用全部的数据库，这一步也是需要跑的
HD.cellchat <- subsetData(HD.cellchat) 
future::plan("multisession", workers = 10) # 并行计算
HD.cellchat <- identifyOverExpressedGenes(HD.cellchat)
HD.cellchat <- identifyOverExpressedInteractions(HD.cellchat)
#默认情况下，cellchat使用object@data.signaling进行网络推断
#但是他也提供了projectData函数，将基因投射到PPI，两者各有优势，原作者说PPI并不会带来假通讯
#或者说很低。具体选取哪种方式，可以做个对比，我们这里使用默认的
# cellchat <- projectData(cellchat, PPI.human)

#===============================================================================
#推断细胞通讯网络
#推断的配体-受体对的数量显然取决于计算每个细胞组平均基因表达的方法。默认情况下，CellChat使用一种称为“trimean”的统计稳健均值方法，该方法产生的交互比其他方法少。然而，我们发现CellChat在预测更强的相互作用方面表现良好，这对于缩小相互作用范围以进行进一步的实验验证非常有帮助。
#上面翻译的这段话什么意思呢，就是说要推断网络，首先要计算celltype的平均基因能表达，而且computeCommunProb推断的时候会设置平均值cut，那么每种方法的cut是不一样的
#cellchat提供了多个方法，默认是一种叫做triMean的方法。这种方法呢产生的interaction会少一点，但是都是interaction较强的
#那么，假设你的推断中没有出现你生物学意义观测的interaction，可以选择换一种方法。例如truncatedMean，选择较低的cut值
HD.cellchat <- computeCommunProb(HD.cellchat,type = "triMean",raw.use = T)
# triMean is used for calculating the average gene expression per cell group. 
# [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2023-12-29 16:33:46]"
# |==========================================================================| 100%
# [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2023-12-29 16:36:20]"


#如果在某些细胞组中只有少数细胞，用户可以过滤掉细胞间通信。默认情况下，每个细胞组中进行细胞间通信所需的最小细胞数为10
#这个过滤是有意义的，毕竟如果只有少出几颗细胞，那么结果太具有偶然性了
HD.cellchat <- filterCommunication(HD.cellchat, min.cells = 100)
# table(HD$cell_type)

#在信号通路水平推断细胞间的通讯
#CellChat通过汇总与每个信号通路相关的所有配体-受体相互作用的通信概率来计算信号通路水平上的通信概率。
#NB: The inferred intercellular communication network of each ligand-receptor pair and each signaling pathway is stored in the slot ‘net’ and ‘netP’, respectively.
HD.cellchat <- computeCommunProbPathway(HD.cellchat)

#互作网络整合,可以设置soure和target，我们这里默认了，就是全部的
HD.cellchat <- aggregateNet(HD.cellchat)

#至此，整个推断就结束了，过程很简单，比比较快

#===============================================================================
#数据提取，subsetCommunication函数
df.net <- subsetCommunication(HD.cellchat) #这里返回的是所有细胞类型受配体数据，基本上所有的信息
#那么也可以设置soure/target来获取特定的互作，但是没啥必要，只要有全部的，subset数据就可以了


#===============================================================================
#                             三、基本可视化
#===============================================================================
#===============================================================================
#可视化1---所有celltype互作数目-------------------------------------------------
#figure1
#circle展示互作数目
groupSize <- as.numeric(table(HD.cellchat@idents)) 
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(HD.cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(HD.cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#figure2
#heatmap展示互作数据
pheatmap::pheatmap(HD.cellchat@net$count, border_color = "black", 
                   cluster_cols = F, fontsize = 10, cluster_rows = F,
                   display_numbers = T,number_color="black",number_format = "%.0f")

#figure3
#贝克图，展示每一个celltype作为source与其他celltype的互作
#指定顺序和指定颜色
# celltype_order <- c( B, Basophil, cDC1, cDC2, Endothelial cell, Fibroblast,  Macrophage,Magakaryocyte,Monocyte, Neutrophil
#                      NK,pDC,Plasma cell,T )
# color.use <- c("#843C39", "#8C6D31", "#E6550D", "#3182BD", "#54990F","#E7BA52", "#31A354", "#E41A1C", "#6BAED6")
# mat <- as.data.frame(HD.cellchat@net$weight)
# mat <- mat[celltype_order,]#行排序
# mat <- mat[,celltype_order] %>% as.matrix()

mat <- HD.cellchat@net$weight
par(mfrow = c(4,4), xpd=TRUE,mar = c(1, 1, 1, 1))
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, arrow.size=0.05,
                   arrow.width=1, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i]
                   # ,color.use = color.use
                   )
}

#cellchat提供的可视化真的是眼花缭乱，有很多，在实际应用中不必每一种都放在文章中
#这里简单的过一遍即可，至于需要呈现哪一种，可以自己抉择，主要是能够展示自己的结果

#按照通路来展示互作
#你可以查看下在你的富集结果中，都有哪些通路.可视化自己感兴趣、与研究相关的即可
HD.cellchat@netP$pathways
#展示的方式也很多

#figure4
#cicrle图

pathways.show <- c("MHC-II") 
netVisual_aggregate(HD.cellchat, signaling = pathways.show)


#figure5
#弦图
netVisual_aggregate(HD.cellchat, signaling = pathways.show, layout = "chord")

#figure6
#热图
#展示特定通路
p1 = netVisual_heatmap(HD.cellchat, signaling = pathways.show,color.heatmap = "Reds")
#展示所有
p2 = netVisual_heatmap(HD.cellchat,color.heatmap = "Reds")
p1+p2


#figure7
#分组展示弦图，也就是说我们的celltype有很多，假设你有很多亚群一起做的分析
#那么就可以设置分组，让你的互作更加清晰，例如这里的例子我们将celltype分为免疫细胞和其他细胞两组
unique(HD.cellchat@idents)
group.cellType <- c(rep("immune", 5), rep("other", 4))
names(group.cellType) <- c("Kers","Mast","Men","Mon","Tcell","ECs","Fibs","lang","SMCs")
group.cellType
netVisual_chord_cell(HD.cellchat, signaling = pathways.show, 
                     group = group.cellType, title.name = paste0(pathways.show, " signaling network"))


#figure8
#计算每个配体-受体对整体信号通路的贡献，并可视化由单个配体-受体对介导的细胞-细胞通信
netAnalysis_contribution(HD.cellchat, signaling = pathways.show)

#figure9
#可视化单个配体-受体对介导的细胞-细胞通讯
#提取特定通路中的受配体对
pairLR <- extractEnrichedLR(HD.cellchat, signaling = pathways.show, geneLR.return = FALSE)
# pairLR <- extractEnrichedLR(HD.cellchat, signaling = HD.cellchat@netP$pathways, geneLR.return = FALSE)
pairLR[1,]
# [1] "HLA-DPA1_CD4"
netVisual_individual(HD.cellchat, signaling = pathways.show, pairLR.use = "HLA-DPA1_CD4", layout = "circle")



#可视化2---可视化细胞通讯受配体对-------------------------------------------------
#上面那么多的可视化，其实都是一类，就是有没有互作，或者互作多少
#那么第二个重要的可视化就是受配体对了，毕竟我们分析了通讯最终的目的还是想只是什么介导的

#figure10
#经典的气泡图（这个数据可自己提取，然后用ggplot作图，这里我们使用cellchat自带的）
#这里是用netVisual_bubble函数，可以指定sourece和target细胞，如果不指定，那就是所有细胞
netVisual_bubble(HD.cellchat, sources.use = "Tcell",  remove.isolate = FALSE)
#当然了，通路和受配体都是可以设定的，查看帮助函数

#figure11
#和弦图展示受配体互作
netVisual_chord_gene(HD.cellchat, sources.use = "Tcell", lab.cex = 0.5,legend.pos.y =50)



#===============================================================================
#                             四、系统分析及可视化
#===============================================================================
#===============================================================================
#cell-cell通信网络系统分析，以便于对结果更好的解释

#识别细胞群的信号作用(例如，主要的发送者，接受者)以及主要的信号作用
#识别细胞间通信网络中的主要的信号发送者、接收者、中介者和影响者。

#计算并可视化网络中心性得分
#figure12
HD.cellchat <- netAnalysis_computeCentrality(HD.cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(HD.cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

#使用散点图在二维空间中可视化占优势的发送者(sourece)和接收器(target)。
#figure13
gg1 <- netAnalysis_signalingRole_scatter(HD.cellchat)#所有通路
gg2 <- netAnalysis_signalingRole_scatter(HD.cellchat, signaling = c("CXCL", "CCL"))#感兴趣通路
gg1 + gg2

#识别对某些细胞群的输出或输入信号贡献最大的信号
#这就可以回答哪些信号对某些细胞群的输出或输入信号贡献最大的问题。
#figure14
ht1 <- netAnalysis_signalingRole_heatmap(HD.cellchat, pattern = "outgoing",width=8,height = 12)
ht2 <- netAnalysis_signalingRole_heatmap(HD.cellchat, pattern = "incoming",width=8,height = 12)
ht1 + ht2


#===============================================================================
#探索多种细胞类型和信号通路如何协调在一起
#除了探索单个通路的详细通信外，一个重要的问题是多个细胞群和信号通路如何协调发挥作用。
#CellChat采用模式识别方法来识别全局通信模式。
#首先确定pattern数量，根据NMF R包中实现的两个指标来推断模式的数量，包括Cophenetic和Silhouette。
#对于一定数量的模式，一个合适的模式数量是Cophenetic值和Silhouette值开始突然下降。
library(NMF)
library(ggalluvial)
#figure15
selectK(HD.cellchat, pattern = "outgoing")


#识别特定细胞群和一般通信模式的主要信号
#figure16
nPatterns = 2
HD.cellchat <- identifyCommunicationPatterns(HD.cellchat, pattern = "outgoing", k = nPatterns)

#figure17,冲击图
netAnalysis_river(HD.cellchat, pattern = "outgoing")

#figure18,气泡图
netAnalysis_dot(HD.cellchat, pattern = "outgoing")

#incoming的鉴定方式和outgoing一样，只不过将参数pattern修改为incoming

#===============================================================================
#多组分析时候应用可能更好
#CellChat能够量化所有重要信号通路之间的相似性，然后根据它们的细胞通信网络相似性对它们进行分组。
#可以根据功能或结构相似性进行分组。
#基于功能相似性来识别信号组
HD.cellchat <- computeNetSimilarity(HD.cellchat, type = "functional")
HD.cellchat <- netEmbedding(HD.cellchat, type = "functional")
HD.cellchat <- netClustering(HD.cellchat, type = "functional")
# figure19
netVisual_embedding(HD.cellchat, type = "functional", label.size = 3.5)


#基于结构相似性识别信号组
HD.cellchat <- computeNetSimilarity(HD.cellchat, type = "structural")
HD.cellchat <- netEmbedding(HD.cellchat, type = "structural")
HD.cellchat <- netClustering(HD.cellchat, type = "structural")
#figure20
netVisual_embedding(HD.cellchat, type = "structural", label.size = 3.5)

#保存文件
saveRDS(HD.cellchat, file = "PB_NAR.cellchat.rds")


#===============================================================================
#                             五、一体化函数分析另一组cellchat
#===============================================================================
#一体化函数
KS_cellchat <- function(input_obj,
                        assay= NULL,
                        group.by = NULL,
                        workers,
                        species=c('human','mouse'),
                        CellChatDB.use=NULL,#a character vector, which is a subset of c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact","Non-protein Signaling")
                        PPIuse=F,
                        type="triMean",# c("triMean", "truncatedMean", "thresholdedMean", "median")
                        min.cells = 10
){
  
  
  
  cellchat.obj = createCellChat(input_obj, assay = assay, group.by = group.by)
  
  if(species=='human'){
    
    CellChatDB <- CellChatDB.human
    ppi = PPI.human
  }
  
  if(species =="mouse"){
    
    CellChatDB <- CellChatDB.mouse
    ppi = PPI.mouse
  }
  
  
  if(is.null(CellChatDB.use)){
    
    cellchat.obj@DB <- CellChatDB
    
  }else{
    
    CellChatDB <- subsetDB(CellChatDB, search = CellChatDB.use, key = "annotation")
    cellchat.obj@DB <- CellChatDB
  }
  
  cellchat.obj <- subsetData(cellchat.obj) 
  future::plan("multisession", workers = workers) # 并行计算
  cellchat.obj <- identifyOverExpressedGenes(cellchat.obj)
  cellchat.obj <- identifyOverExpressedInteractions(cellchat.obj)
  
  if(PPIuse==F){
    
    cellchat.obj <- computeCommunProb(cellchat.obj, type = type)
    
  }else{
    
    cellchat.obj <- projectData(cellchat.obj, ppi)
    cellchat.obj <- computeCommunProb(cellchat.obj, raw.use=F, type = type)
  }
  
  
  cellchat.obj <- filterCommunication(cellchat.obj, min.cells = min.cells)
  cellchat.obj <- computeCommunProbPathway(cellchat.obj)
  
  #互作网络整合,可以设置soure和target，我们这里默认了，就是全部的
  cellchat.obj <- aggregateNet(cellchat.obj)
  
  return(cellchat.obj)
  
}


#===============================================================================
#分析MDA组的通讯
table(sc_pbmc_int2$subcell_type2,sc_pbmc_int2$group)


sceV5 <- subset(sc_pbmc_int2,type%in%c('LIVER')&subcell_type2%in%c('1:CD4_C1_LEF1','2:CD4_C2_IL7R', '3:CD8_C1_GZMK', '4:CD8_C2_RGS1','5:CD8_C3_GZMH',
                                                               '6:CD8_C4_LEF1', '7:CD8_C5_KLRB1', '8:CD8_C6_TCF7', '9:CD8_C7_FOSB','10:CD8_C8_MKI67',
                                                               '11:CD8_C9_TIGIT', '12:NK_C1_FCGR3A','13:NK_C2_XCL1', '14:NK_C3_KLRF1', 
                                                               '8:Mono_C4_ATF3',
                                                               '13:Mac_C1_CCL3'))


# sceV5 <- sc_pbmc_int2
HD <- subset(sceV5,GROUP=='NAR')
MDA <- subset(sceV5,GROUP=='AR')
table(HD$cell_type)
table(MDA$cell_type)

HD <- subset(sc_pbmc_int2,group=='LI_NAR'&subcell_type2%in%c('4:CD8_C2_RGS1','10:CD8_C8_MKI67',
                                                      '13:NK_C2_XCL1',
                                                      '8:Mono_C4_ATF3',
                                                      '13:Mac_C1_CCL3'))
MDA <- subset(sc_pbmc_int2,group=='LI_AR'&subcell_type2%in%c('4:CD8_C2_RGS1','10:CD8_C8_MKI67',
                                                             '13:NK_C2_XCL1',
                                                             '8:Mono_C4_ATF3',
                                                             '13:Mac_C1_CCL3'))
table(HD$group,HD$subcell_type2)
table(MDA$group,MDA$subcell_type2)

gc()
options (future.globals.maxSize = 2200 * 1024^4) # 设置环境大小限制，稍大于提示信息要求的就行
plan("multisession", workers = 10 ) 
plan()

HD.cellchat <- KS_cellchat(HD, assay = 'RNA', group.by = "subcell_type2",
                            workers=10, species='human',PPIuse=TRUE)

gc()
options (future.globals.maxSize = 2200 * 1024^4) # 设置环境大小限制，稍大于提示信息要求的就行
plan("multisession", workers = 10 ) 
plan()
MDA.cellchat <- KS_cellchat(MDA, assay = 'RNA', group.by = "subcell_type2",
                            workers=10, species='human',PPIuse=TRUE)

#运行日志
# [1] "Create a CellChat object from a Seurat object"
# The `meta.data` slot in the Seurat object is used as cell meta information 
# Set cell identities for the new CellChat object 
# The cell groups used for CellChat analysis are  ECs, Fibs, Kers, lang, Mast, Men, Mon, SMCs, Tcell 
# triMean is used for calculating the average gene expression per cell group. 
# [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2024-01-02 10:42:49]"
# |===========================================================================| 100%
# [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2024-01-02 10:50:02]"

setwd('./new')
saveRDS(HD.cellchat, file = "LI_NAR.cellchatPPI.rds")
saveRDS(MDA.cellchat, file = "LI_AR.cellchatPPI.rds")
#===============================================================================
#                             六、cellchat多组比较分析
#===============================================================================
#merge cellchat obj of different group
object.list <- list(NAR = HD.cellchat, AR = MDA.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
# The cell barcodes in merged 'meta' is   
# Error in `.rowNamesDF<-`(x, value = value) : 'row.names'的长度不对

# MDA.cellchat@meta <- MDA.cellchat@meta[,c(9,38)]
# colnames(MDA.cellchat@meta) <- c("group","labels")
# 
# HD.cellchat@meta <- HD.cellchat@meta[,c(9,38)]
# colnames(HD.cellchat@meta) <- c("group","labels")
#===============================================================================
#基本比较及可视化，例如互作数目的多少等等
#多组比较主要回答一下的生物学问题，在实际中，我们设置对照组，实验组，肯定是想知道处理之后有什么变化
# 细胞间通讯是否增强
# 哪些细胞类型之间的相互作用发生了显著变化
# 主要source和target是如何从一种情况变化到另一种情况

#figure21
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
ggsave(filename = 'figure21.pdf',width = 12,height = 8)

#figure22
#不同细胞群之间相互作用数量或相互作用强度的差异.红色表示与第一个数据集相比，互作增强，蓝色对应的是减弱
dev.off()
pdf( file = "figure22.pdf", width = 12, height = 8 )
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE)
netVisual_diffInteraction(cellchat, weight.scale = TRUE, measure = "weight")
dev.off()

#还可以使用热图更详细地显示相互作用的不同数量或相互作用强度。顶部的彩色条形图表示热图(incoming signaling)中显示的列值之和。右边的彩色条形图表示值行之和(outcoming signaling。
#figure23
dev.off()
pdf( file = "figure23.pdf", width = 12, height = 8 )
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2
dev.off()




#很显然，以上的比较只适用于配对的2组数据
#如果有更多的数据集进行比较，可以直接显示每个数据集中任意两个细胞群之间的相互作用数量或相互作用强度。
#figure24，其实就是每组的互作数据展示在一起而已
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = TRUE, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


#figure25,修改下netVisual_circle的参数，显示互作数量，让图更加直观。
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = TRUE, label.edge= T, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}


#===============================================================================
#比较主要的信号
#比较二维空间中out和in的相互作用强度，可以随时识别在不同数据集之间发送或接收信号发生重大变化的细胞群。
#figure26
dev.off()
pdf( file = "figure26MHCI.pdf", width = 16, height = 8 )
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i])+
  scale_y_continuous(limits = c(0,15)) + scale_x_continuous(limits = c(0,15))
  }
patchwork::wrap_plots(plots = gg)
dev.off()

+scale_y_continuous(limits = c(0,4)) + scale_x_continuous(limits = c(0,5))

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Error in netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i],  :
#                                              Please run `netAnalysis_computeCentrality` to compute the network centrality scores!

##计算网络
object.list[[1]] <- netAnalysis_computeCentrality(object.list[[1]], slot.name = "netP")
selectK(object.list[[1]], pattern = "outgoing")
nPatterns = 2
object.list[[1]] <- identifyCommunicationPatterns(object.list[[1]], pattern = "outgoing", k = nPatterns)


object.list[[2]] <- netAnalysis_computeCentrality(object.list[[2]], slot.name = "netP")
selectK(object.list[[2]], pattern = "outgoing")
nPatterns = 2
object.list[[2]] <- identifyCommunicationPatterns(object.list[[2]], pattern = "outgoing", k = nPatterns)


#在识别了主要的信号celltype变化后，我们还可以比较两组之间celltype特定信号的变化
#figure27
netAnalysis_signalingChanges_scatter(cellchat, idents.use = '12:NK_C1_FCGR3A',
                                     # signaling.exclude = "MHC-II",
                                     comparison = c(2, 1))


#===============================================================================
#CellChat可以识别具有较大(或较小)差异的信号网络，信号组，
#以及基于它们在多种生物条件下的细胞-细胞通信网络的保守和特异性信号通路。
#按照功能识别
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, umap.method = 'uwot',type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
#figure28
netVisual_embeddingPairwise(cellchat, type = "functional", label.size = 3.5)



#通过比较每个信号通路的信息流/相互作用强度，我们可以识别信号通路，(i)关闭，(ii)减少，(iii)打开或(iv)增加，通过改变它们在一种条件下的信息流与另一种条件相比。
#比较每个通路整体信息流
#这个柱形图对重要的信号通路进行了排序。顶部红色的信号通路在HD中富集，这些绿色的信号通路在MDA中富集。
#figure29
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T,
               cutoff.pvalue = 0.01,thresh = 0.01,do.stat = TRUE,do.flip=TRUE,axis.gap = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
ggsave(filename = 'figure29.pdf',width = 12,height = 12)


#比较与每个细胞群相关的传出(或传入)信号
#figure30
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
pathway.union <- object.list[[i]]@netP$pathways
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 12)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 6, height = 12)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

pathway.union <- object.list[[i+1]]@netP$pathways
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 12)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 6, height = 12)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))


#===============================================================================
#比较不同组细胞互作配体-受体对
##figure31
netVisual_bubble(cellchat,sources.use = c('12:NK_C1_FCGR3A','13:NK_C2_XCL1', '14:NK_C3_KLRF1'),
                 targets.use=c('8:Mono_C4_ATF3','13:Mac_C1_CCL3'), comparison = c(1, 2), angle.x = 45)
#===============================================================================
#差异比较
#通过差异表达分析识别功能障碍信号.基于差异表达分析(DEA)来识别上调和下调的信号配体-受体对。
#具体而言，就是先做差异分析，然后根据发送细胞中配体和接收细胞中受体的折叠变化获得上调和下调的信号。
pos.dataset = "AR"#设置数据集。这里表示logFC相比于其他数据集上调的数据集
features.name = paste0(pos.dataset, ".merged")
cellchat <- identifyOverExpressedGenes(cellchat, 
                                       group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name, 
                                       only.pos = FALSE, 
                                       thresh.pc = 0.1, thresh.fc = 0.05,
                                       thresh.p = 0.05, group.DE.combined =TRUE) 

#将差异表达分析的结果映射到推断的细胞-细胞通信上，以便轻松地管理感兴趣的配体-受体对
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
#提取ligand-receptor对MDA中上调的配体
net.up <- subsetCommunication(cellchat, 
                              net = net, 
                              datasets = "AR",
                              ligand.logFC = 0.05, 
                              receptor.logFC = NULL)
#提取ligand-receptor对HD中上调的配体和上调的受体
net.down <- subsetCommunication(cellchat,
                                net = net, 
                                datasets = "NAR",
                                ligand.logFC = -0.05, 
                                receptor.logFC = NULL)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

pairLR.use.up = net.up[, "interaction_name", drop = F]

my_palette <- colorRampPalette(c("darkblue","yellow","red"))(n=100)
#原始
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, 
                        sources.use = c('12:NK_C1_FCGR3A','13:NK_C2_XCL1', '14:NK_C3_KLRF1'),
                        targets.use=c('13:Mac_C1_CCL3'),
                        comparison = c(1, 2),  
                        angle.x = 45, 
                        remove.isolate = TRUE,
                        title.name = paste0("Up-regulated signaling in ", 
                                            names(object.list)[2]),
                        line.on = TRUE,
                        min.quantile=0.1,
                        max.quantile=1,
                        dot.size.min=2,
                        dot.size.max=4,
                        color.text=c('black','black'),
                        font.size = 14) +    
  scale_color_gradientn('Mean expression', colors=my_palette) 

  
  gg1

 table( cellchat@idents$AR)
#修改
 HLA <- pairLR.use.up[grep(c("^CCL|CXCL"), pairLR.use.up$interaction_name),] 
 HLA <- as.data.frame(HLA)
 colnames(HLA) <- c('interaction_name')

 
 gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, 
                        sources.use = c('12:NK_C1_FCGR3A','13:NK_C2_XCL1', '14:NK_C3_KLRF1'),
                        targets.use=c('13:Mac_C1_CCL3','8:Mono_C4_ATF3'),
                        comparison = c(1, 2),  
                        angle.x =315, 
                        vjust.x = 1,
                        hjust.x =0,
                        remove.isolate = TRUE,
                        line.on = TRUE,
                        max.quantile=1,
                        dot.size.min=2,
                        dot.size.max=6,
                        color.text=c('black','black'),
                        font.size = 10)+
  theme(legend.position = 'right')+
   coord_flip() + 
  scale_color_gradientn('Mean expression', colors=my_palette) 
   # scale_y_discrete(position = "top")

gg1



# 
# library(jjAnno)
# p2 <- annoSegment(object = gg1,
#                   annoPos = 'top', 
#                   aesGroup = T,
#                   aesGroName = 'source',
#                   yPosition = c(48),
#                   lwd = 18,
#                   segWidth = 0.6,
#                   pCol=c('#4C8562','#A24F1D','#65A8AE'),
#                   addText=T,
#                   textSize =8,
#                   textHVjust = 0.05,
#                   textCol = c("black","black",'black'))

# '1:CD4_C1_LEF1','2:CD4_C2_IL7R', 
#'3:CD8_C1_GZMK', '4:CD8_C2_RGS1','5:CD8_C3_GZMH','6:CD8_C4_LEF1', '7:CD8_C5_KLRB1', '8:CD8_C6_TCF7', '9:CD8_C7_FOSB','10:CD8_C8_MKI67', '11:CD8_C9_TIGIT'
#'12:NK_C1_FCGR3A','13:NK_C2_XCL1', '14:NK_C3_KLRF1'
ggsave(filename = 'T.pdf',width = 12,height = 6,dpi=300)
ggsave(filename = 'NK.pdf',width = 12,height = 6,dpi=300)




pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down,  
                        sources.use = c('8:Mono_C4_ATF3','13:Mac_C1_CCL3'),
                        targets.use=c('4:CD8_C2_RGS1','5:CD8_C3_GZMH','3:CD8_C1_GZMK'),
                        comparison = c(1, 2),  
                        angle.x = 90, remove.isolate = TRUE,
                        title.name = paste0("Down-regulated signaling in ", 
                                            names(object.list)[2]))
gg2
gg1 + gg2 
##figure32
ggsave(filename = 'figure32MHCII.pdf',width = 12,height = 8)

#===============================================================================
#不同通路的互作比较
#figure33
pathways.show='MHC-I'
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}


#figure34 弦图
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}


save(cellchat, file = 'cellchat_merge.RData')

