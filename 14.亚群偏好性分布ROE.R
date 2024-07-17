#+++, Ro/e > 1; 
#++, 0.8 < Ro/e ≤ 1; 
#+, 0.2 ≤ Ro/e ≤ 0.8; 
#+/−, 0 < Ro/e < 0.2; 
#−, Ro/e = 0
#remotes::install_github("Japrin/sscVis", quiet = TRUE)
#remotes::install_github("Japrin/sscClust", quiet = TRUE)
#remotes::install_github("Japrin/Startrac", quiet = TRUE)

sc_pbmc_int2 <- scpbmc_int2 

rm(list=ls())
gc()
library(sscVis)
library(Seurat)
library(SeuratObject)
library(data.table)
library(readr)
library(dplyr)
library(sscClust)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(Startrac)


setwd("/mnt/data/home/tycloud/LNM_doublet")
sce <- readRDS("/mnt/data/home/tycloud/LNM_doublet/sce_thy_malig_UCell_k_NMF.rds")
sce <- sc_pbmc_int2
cols <- c( '#7cbb5f', '#368650', '#E6550D', '#5e4d9a', '#78c2ed', '#866017',
           '#E41A1C', '#e0dfed', '#279ad7', '#ef7b77')
DimPlot(sce,reduction = "umap",group.by = "cell_type",
        split.by = "orig.ident",label =F,cols=cols)+DimPlot(sce,reduction = "umap",group.by = "orig.ident",
                                                            label =F)
DotPlot(sce, group.by = "mincelltype"  ,features = markers, col.min = 0)+coord_flip()

#数据来源文章：Lineage tracking reveals dynamic relationships of T cells in colorectal cancer
#都可以从对应的NCBI网页下载：
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108989


print(getwd())
#输出路径及文件前缀
out.prefix = "~/home/data/t010560/God/mt25/1.nk and t"
#自定义主题
if(T){
  text.size = 8
  text.angle = 45
  text.hjust = 1
  legend.position = "right"
  mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                   # axis.ticks = element_line(color = "black"),
                   axis.title = element_text(size = text.size,color ="black"), 
                   axis.text = element_text(size=text.size,color = "black"),
                   axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
                   #axis.line = element_line(color = "black"),
                   #axis.ticks = element_line(color = "black"),
                   #panel.grid.minor.y = element_blank(),
                   #panel.grid.minor.x = element_blank(),
                   panel.grid=element_blank(), # 去网格线
                   legend.position = legend.position,
                   legend.text = element_text(size= text.size),
                   legend.title= element_text(size= text.size)
                   # panel.border = element_rect(size = 0.7, lineorig.ident = "solid", colour = "black")
                   # strip.background = element_rect(color="black",size= 1, lineorig.ident="solid") # fill="#FC4E07", 
  )
}

#### 1.读入数据
meta.tb <- sc_pbmc_int2@meta.data
head(meta.tb)
saveRDS(meta.tb,"sce_meta.data.rds")


#### 2.数据预处理、计算Ro/e
### 2.1 数据预处理
cellInfo.tb = meta.tb %>% data.table()
cellInfo.tb$cell_type = as.character(cellInfo.tb$cell_type)###行
cellInfo.tb$group = cellInfo.tb$group %>% as.factor()
loc.avai.vec <- levels(cellInfo.tb[["group"]])

### 2.2 计算Ro/e值
startrac.dist <- unclass(Startrac::calTissueDist(cellInfo.tb,byPatient = F,
                                                 colname.cluster="cell_type",
                                                 colname.patient = "patient",
                                                 colname.tissue = "group"))
startrac.dist <- startrac.dist[,loc.avai.vec]



### 2.3 做分类变量切分
cuts <- c(0.0000000, 0.0000001, 0.8, 1.5, 2, Inf)

# cuts <- c(0.0000000, 0.0000001, 0.5, 1.0, 1.5, Inf)
startrac.dist.bin.values <- factor(c("-", "+/-", "+", "++", "+++"),levels=c("-", "+/-", "+", "++", "+++"))
startrac.dist.bin <- matrix(startrac.dist.bin.values[findInterval(startrac.dist, cuts)],
                            ncol=ncol(startrac.dist))

colnames(startrac.dist.bin) <- colnames(startrac.dist)
rownames(startrac.dist.bin) <- rownames(startrac.dist)


#输出路径及文件前缀
out.prefix = "~/t"
#### 3. 可视化
#### 3. 可视化
### 3.1 plotMatrix.simple可视化
sscVis::plotMatrix.simple(startrac.dist,
                          out.prefix=sprintf("%s.startrac.dist",out.prefix),
                          show.number=F,
                          clust.row=T,
                          #waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                          exp.name=expression(italic(R)[o/e]),
                          z.hi=2,
                          #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
                          palatte=viridis::viridis(7),
                          pdf.width = 4.5, pdf.height = 4.5)

### 3.2 pheatmap可视化
#限定范围
startrac.dist1 = ifelse(startrac.dist>2,2,startrac.dist)

#不限定
startrac.dist1 <- startrac.dist

plot_2 = pheatmap::pheatmap(startrac.dist,
                            color = viridis::viridis(7),
                            cluster_rows = T,
                            cluster_cols = F,
                            fontsize = 10,
                            border_color = "white",
                            treeheight_col = 0,
                            treeheight_row = 0
)


print(plot_2)



### 3.3 ggplot2可视化
plot_data = reshape2::melt(startrac.dist1)
head(plot_data)
colnames(plot_data) = c("mincelltype", "group", "Roe")
startrac.dist.bin2 = reshape2::melt(startrac.dist.bin)
plot_data$bin = startrac.dist.bin2$value %>% factor(levels=c("-", "+/-", "+", "++", "+++"))

plot_data$mincelltype = factor(plot_data$mincelltype,
                               levels = rev(plot_2$tree_row$labels[plot_2$tree_row$order]))
head(plot_data)
#plot_data<-plot_data[22:98,]
##############===========================================
p3 <- ggplot(data = plot_data, aes(x = group, y = mincelltype, fill = Roe)) +
  geom_tile() +
  theme_minimal() +
  scale_x_discrete(expand = c(0.0, 0.0)) +
  scale_y_discrete(expand = c(0.0, 0.0)) +
  scale_fill_gradientn(name = "Ro/e",
                       colours = viridis::viridis(7)) +
  labs(x = NULL, y = NULL, title = "Ro/e") +
  mytheme
p3

table(sc_pbmc_int2$cell_type)

### 3.4 等级可视化
plot_data$mincelltype=factor(plot_data$mincelltype,
                             levels = rev(c(  "1:CD4_C1_LEF1" ,  "2:CD4_C2_IL7R" ,  "3:CD8_C1_GZMK" ,  "4:CD8_C2_RGS1" ,  "5:CD8_C3_GZMH" ,  "6:CD8_C4_LEF1" , 
                                          "7:CD8_C5_KLRB1",  "8:CD8_C6_TCF7"  , "9:CD8_C7_FOSB" ,  "10:CD8_C8_MKI67", "11:CD8_C9_TIGIT", "12:NK_C1_FCGR3A",
                                          "13:NK_C2_XCL1" ,  "14:NK_C3_KLRF1" , "15:DN_LINC-PINT")))

plot_data$group=factor(plot_data$group,
                             levels = (c('LI_AR','PB_AR','LI_NAR','PB_NAR')))


###绘制带数字的热图
plot_data$Roe <- round(plot_data$Roe, 2)
ggplot(plot_data,aes(group,mincelltype))+
  geom_point(aes(fill=Roe),shape=22, color='white',size=14)+
  scale_fill_gradientn(colours =colorRampPalette(c("#3162A0","white","#C7252C"))(7))+
  labs(x = NULL,y = NULL,fill="Enrichment(Ro/e)")+
  guides(fill=guide_colorbar(barheight = 8))+
  theme_bw()+  theme(panel.grid.major= element_blank(),
                     panel.grid.minor= element_blank())+
  theme(axis.text=element_text(colour='black',size=10))+
  theme(axis.text.x =element_text(angle =0,hjust =0.9,vjust = 0.9,size=10))+
  geom_text(aes(label = Roe), color = "black", size = 3, vjust = 0.3, hjust = 0.5)
ggsave(filename = 'Roe.pdf',width = 5.5,height = 8)





p4 <- ggplot(data = plot_data, aes(x = group, y = mincelltype, fill = bin)) +
  geom_tile() +
  scale_fill_manual(
    name = NULL,
    values = c(
      "-" = "white",
      "+/-" = "#fde6ce",
      "+" = "#fcc08b",
      "++" = "#f5904a",
      "+++" = "#e6540d"
    ),
    labels = c(0,0.8, 1.5, 2, ">2")
  ) +
  geom_text(
    aes(label = bin),
    color = "black",
    size = 3,
    show.legend = T
  ) +
  guides(fill = guide_legend(title = "Ro/e",
                             override.aes = list(label = c( "+/-", "+", "++","+++"
                             )))) + theme_minimal() +
  scale_x_discrete(expand = c(0.0, 0.0)) +
  scale_y_discrete(expand = c(0.0, 0.0)) +
  labs(x = NULL, y = NULL, title = "Ro/e") +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, "cm")) + mytheme;p4

ggsave("ROe_c1_c7_origident.pdf",plot=p4,width = 5,height = 8)

######type

#### 2.数据预处理、计算Ro/e
### 2.1 数据预处理
cellInfo.tb = meta.tb %>% data.table()
cellInfo.tb$mincelltype = as.character(cellInfo.tb$mincelltype)###行
cellInfo.tb$type = cellInfo.tb$type %>% as.factor()
loc.avai.vec <- levels(cellInfo.tb[["type"]])

### 2.2 计算Ro/e值
startrac.dist <- unclass(Startrac::calTissueDist(cellInfo.tb,byPatient = F,
                                                 colname.cluster="mincelltype",
                                                 colname.patient = "patient",
                                                 colname.tissue = "type"))
startrac.dist <- startrac.dist[,loc.avai.vec]



### 2.3 做分类变量切分
cuts <- c(0.0000000, 0.0000001, 0.8, 1.5, 2, Inf)
startrac.dist.bin.values <- factor(c("-", "+/-", "+", "++", "+++"),levels=c("-", "+/-", "+", "++", "+++"))
startrac.dist.bin <- matrix(startrac.dist.bin.values[findInterval(startrac.dist, cuts)],
                            ncol=ncol(startrac.dist))

colnames(startrac.dist.bin) <- colnames(startrac.dist)
rownames(startrac.dist.bin) <- rownames(startrac.dist)



#### 3. 可视化
#### 3. 可视化
### 3.1 plotMatrix.simple可视化
sscVis::plotMatrix.simple(startrac.dist,
                          out.prefix=sprintf("%s.startrac.dist",out.prefix),
                          show.number=F,
                          clust.row=T,
                          #waterfall.row=T,par.warterfall = list(score.alpha = 2,do.norm=T),
                          exp.name=expression(italic(R)[o/e]),
                          z.hi=2,
                          #palatte=rev(brewer.pal(n = 7,name = "RdYlBu")),
                          palatte=viridis::viridis(7),
                          pdf.width = 4.5, pdf.height = 4.5)

### 3.2 pheatmap可视化
startrac.dist1 = ifelse(startrac.dist>2,2,startrac.dist)

plot_2 = pheatmap::pheatmap(startrac.dist1,
                            color = viridis::viridis(7),
                            cluster_rows = T,
                            cluster_cols = F,
                            fontsize = 10,
                            border_color = "white",
                            treeheight_col = 0,
                            treeheight_row = 0
)


print(plot_2)



### 3.3 ggplot2可视化
plot_data = reshape2::melt(startrac.dist1)
head(plot_data)
colnames(plot_data) = c("mincelltype", "type", "Roe")
startrac.dist.bin2 = reshape2::melt(startrac.dist.bin)
plot_data$bin = startrac.dist.bin2$value %>% factor(levels=c("-", "+/-", "+", "++", "+++"))

plot_data$mincelltype = factor(plot_data$mincelltype,
                               levels = rev(plot_2$tree_row$labels[plot_2$tree_row$order]))
head(plot_data)
#plot_data<-plot_data[8:21,]
##############===========================================
p3 <- ggplot(data = plot_data, aes(x = type, y = mincelltype, fill = Roe)) +
  geom_tile() +
  theme_minimal() +
  scale_x_discrete(expand = c(0.0, 0.0)) +
  scale_y_discrete(expand = c(0.0, 0.0)) +
  scale_fill_gradientn(name = "Ro/e",
                       colours = viridis::viridis(7)) +
  labs(x = NULL, y = NULL, title = "Ro/e") +
  mytheme
p3



### 3.4 等级可视化
plot_data$mincelltype=factor(plot_data$mincelltype,
                             levels = c( "Normal","c1", "c2" ,"c3" ,
                                         "c4", "c5" , "c6", "c7" ))
p4 <- ggplot(data = plot_data, aes(x = type, y = mincelltype, fill = bin)) +
  geom_tile() +
  scale_fill_manual(
    name = "Ro/e",
    values = c(
      "-" = "white",
      "+/-" = "#fde6ce",
      "+" = "#fcc08b",
      "++" = "#f5904a",
      "+++" = "#e6540d"
    ),
    labels = c(0,0.8, 1.5, 2, ">2")
  ) +
  geom_text(
    aes(label = bin),
    color = "black",
    size = 3,
    show.legend = T
  ) +
  guides(fill = guide_legend(title = "Ro/e",
                             override.aes = list(label = c( "+/-", "+", "++", "+++")))) + theme_minimal() +
  scale_x_discrete(expand = c(0.0, 0.0)) +
  scale_y_discrete(expand = c(0.0, 0.0)) +
  labs(x = NULL, y = NULL, title = "Ro/e") +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.length.y = unit(0, "cm")) + mytheme;p4

ggsave("ROe_normal_c1_c7_type.pdf",plot=p4,width = 3,height = 5)
