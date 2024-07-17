#=================================================================================
#                   ggplot修饰monocle2热图
#=================================================================================

#一、颜色修改------------------------------------------------------------------

setwd('D:/KS项目/公众号文章/ggplot修饰monocle热图')
library(monocle)
library(monocle)
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)
library(dplyr)
library(tidytree)
library(viridis)
library(scales)
load("D:/KS项目/公众号文章/ggplot修饰monocle热图/cds.RData")#拟时cds
#查看随着pseudotime变化的基因
# cds_DGT_pseudotimegenes <- differentialGeneTest(mouse_monocle,fullModelFormulaStr = "~sm.ns(Pseudotime)")
#这里完全是为了展示较少的基因所以控制了阈值，实际需要展示什么基因，自己选择
cds_DGT_pseudotimegenes_sig <- subset(ds_DGT_pseudotimegenes, qval < 0.01)

mouse_data <- readRDS("D:/KS项目/公众号文章/mouse_data.rds")
# marker <- FindAllMarkers(mouse_data, only.pos = T,logfc.threshold = 0.5)
# marker <- marker[which(marker$p_val_adj<0.05),]
top15 <- marker %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
top15_ordergene <- cds_DGT_pseudotimegenes_sig[top15$gene, ]
Time_genes <- top15_ordergene %>% pull(gene_short_name) %>% as.character()
Time_genes <- unique(Time_genes)

p <- plot_pseudotime_heatmap(mouse_monocle[Time_genes,], 
                             num_cluster = 4, 
                             show_rownames = T, 
                             return_heatmap = T,
                             hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))

p <- plot_pseudotime_heatmap(mouse_monocle[Time_genes,], 
                             num_cluster = 4, 
                             show_rownames = T, 
                             return_heatmap = T,
                             hmcols =colorRampPalette(rev(brewer.pal(9, "PRGn")))(100))

p <- plot_pseudotime_heatmap(mouse_monocle[Time_genes,], 
                             num_cluster = 4, 
                             show_rownames = T, 
                             return_heatmap = F,
                             hmcols = viridis(256))


#二、数据提取，pheatmap重现-----------------------------------------------------------------
#提取数据
cds_subset <- cds
newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime), 
                                       max(pData(cds_subset)$Pseudotime),
                                       length.out = 100)) 
m <- genSmoothCurves(cds[gene,], 
                     trend_formula = '~sm.ns(Pseudotime, df=3)',  
                     relative_expr = T, new_data = newdata)
m=m[!apply(m,1,sum)==0,]
m <- log10(m+1) #log处理数据
#数据缩放
m=m[!apply(m,1,sd)==0,]
m=Matrix::t(scale(Matrix::t(m),center=TRUE))
m[m > 3] = 3#热图最大值
m[m <- 3] = -3#热图最小值，可自行调整

##排序设置
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
library(pheatmap)
#先做一个普通热图

p1 <- pheatmap(m, 
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               clustering_method = "ward.D2",
               cutree_rows=3,
               filename=NA,
               border_color = NA,
               fontsize_row = 8,
               color=colorRampPalette(c("navy","white","firebrick3"))(100),
               clustering_callback = callback)


#列注释
annotation_col = data.frame(
  pseudotime = rescale(newdata$Pseudotime, to = c(-1, 1)))
row.names(annotation_col) <- colnames(m)


#行注释
annotation_row <- data.frame(Cluster=factor(cutree(p1$tree_row, 4)))
row.names(annotation_row) <- rownames(m)

rowcolor <- c("#85B22E","#E29827","#922927",'#57C3F3') 
names(rowcolor) <- c("1","2","3","4") #类型颜色

#注释颜色修改
ann_colors <- list(pseudotime=viridis(100),
                   Cluster=rowcolor) #颜色设置

library(colorRamps)
bks <- seq(-3.1, 3.1, by = 0.1)
hmcols <- blue2green2red(length(bks) - 1) 
rm(hmcols)



p1 <- pheatmap(m,
               useRaster = TRUE,
               cluster_cols = FALSE,
               cluster_rows = TRUE,
               show_rownames = TRUE,
               show_colnames = FALSE,
               clustering_method = "ward.D2",
               cutree_rows = 4,
               filename = NA,
               border_color = NA,
               fontsize_row = 8,
               color = colorRampPalette(c("navy","white","firebrick3"))(100),
               clustering_callback = callback)

dev.off()
pdf( file = "monocle_heatmap.pdf", width = 10, height =10)
p2 <- pheatmap(m, 
               useRaster = T,
               cluster_cols=FALSE, 
               cluster_rows=T, 
               show_rownames=T, 
               show_colnames=F, 
               clustering_method = "ward.D2",
               cutree_rows=4,
               filename=NA,
               border_color = NA,
               fontsize_row = 6,
               color=hmcols,
               # annotation_col = annotation_col,
               annotation_colors=ann_colors,
               annotation_row = annotation_row,
               clustering_callback = callback,
               annotation_legend = F,
               annotation_names_col = F,
               annotation_names_row = F,
               main="Pseudotime")

dev.off()

plotdf=pData(cds)
ggplot(plotdf, aes(x=Pseudotime,y=cell_type,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(),linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )
ggsave("tmp2.pdf",width = 16,height = 7,units = "cm",dpi = 300)




library(grid)

标签基因


genes <- c("CCR7", "SELL", "LEF1", "TCF7")
#三、数据提取，ggplot重现------------------------------------------------------
#数据就是前面提取的数据
#但是需要将数据转化为ggplot的长格式
heat_gg <- m
heat_gg <- as.data.frame(heat_gg)
heat_gg <- heat_gg%>% mutate(gene=row.names(.)) %>% melt()#转化为ggplot画图需要的长列表

p3 <- ggplot(heat_gg,aes(x=variable,y=gene,fill=value))+
  geom_raster()+
  scale_fill_gradient2(low="#003366", high="#990033", mid="white",
                       name='Pseudotime')+
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = 'black',
                                   size = 8))+
  scale_y_discrete(position = "right")

#层次聚类
library(ggtree)
library(tidytree)
d <- m%>% as.data.frame()
hc <- hclust(dist(d), method = "ward.D2")
clus <-cutree(hc, 4)
d1 = data.frame(label=names(clus),member=factor(clus))
ph <- ggtree(as.phylo(hc))

#行注释
cluster <- d1 %>%
  ggplot(aes(x=1, y=label,fill=member))+
  geom_tile() + 
  scale_y_discrete(position="right") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank())+
  labs(fill = "cluster")+
  scale_fill_manual(values = c("#708090",'#68A180','#F3B1A0', '#D6E7A3'))


#列注释
group <- colnames(m) %>% as.data.frame() %>% 
  mutate(group=newdata$Pseudotime) %>%
  mutate(p="group") %>%
  ggplot(aes(.,y=1,fill=group))+
  geom_tile() + 
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text = element_blank(),
        panel.grid = element_blank())+
  scale_fill_gradientn(colours = c("#85B22E","#E29827",'#57C3F3',"#922927"))+
  scale_x_discrete(expand = c(0,0))+
  geom_segment(aes(x = 5, y = 1, xend = 95, yend = 1),
               arrow = arrow(length = unit(0.1, "inches"),
                             type = 'closed'))+
  theme(plot.margin = margin(0,0,0,0))



library(aplot)
p4 <- p3 %>%insert_left(cluster, width = 0.04)%>%
  insert_top(group, height =0.02)%>%
  insert_left(ph,width=0.1)

p4


#山脊图-------------------------------------------------------------------------
#这里拼图各显神通吧，AI修饰最好
library(ggridges)
plotdf=pData(mouse_monocle)
plotdf$celltype <- factor(plotdf$celltype, 
                          levels = c("PMN(4)","PMN(1)","PMN(3)","PMN(5)",
                                     "PMN(2)","PMN(6)","PMN(0)","PMN(7)"))
#细胞拟时顺序是按照下面的图来设定的
# plot_cell_trajectory(mouse_monocle, cell_size = 2.2, color_by = "celltype") +
#   facet_wrap(~celltype, nrow = 2)

p5 <- ggplot(plotdf, aes(x=Pseudotime,y=celltype,fill=celltype))+
  geom_density_ridges(scale=1) +
  scale_y_discrete(position = 'right')+
  theme_minimal()+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(colour = 'black', size=8))+
  scale_x_continuous(position = 'top')


p6 <- ggplot(heat_gg,aes(x=variable,y=gene,fill=value))+
  geom_raster()+
  scale_fill_gradient2(low="#003366", high="#990033", mid="white",
                       name='Pseudotime')+
  theme(axis.title = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = 'black',
                                   size = 8))+
  scale_y_discrete(position = "right")

library(patchwork)
p5+p6+plot_layout(ncol = 1,heights=c(1,5),guides = 'collect')


#四、ggplot重做基因拟时细胞散点图------------------------------------------------------

genes <- c("Anxa1", "Ncf1","Ltf","Camp","Ngp","Chil3")
genes_exp <- list()

for(i in 1:length(genes)){
  A <- log2(exprs(mouse_monocle)[genes[i],]+1)
  A <- as.data.frame(A)
  genes_exp[[i]] <- A
}

gene_exp <- do.call(cbind, genes_exp)
colnames(gene_exp) <- genes
#将上述几个基因的拟时表达添加到monocle
pData(mouse_monocle) = cbind(pData(mouse_monocle), gene_exp)

#提取作图数据，只需要游基因表达和拟时即可
data <- pData(mouse_monocle)
colnames(data)
#选择需要的列即可，我这里的origin.ident就是分组
data<-data %>% select("celltype","Pseudotime","Anxa1", "Ncf1",
                      "Ltf","Camp","Ngp","Chil3","S100a8","S100a9")

features <- c("Anxa1", "Ncf1","Ltf","Camp","Ngp","Chil3","S100a8","S100a9")
plist <- list()
#ggplot作图
for (i in 1:length(features)){
  df <- data[, colnames(data)%in% c("celltype","Pseudotime",features[i])]
  colnames(df) <- c("celltype","Pseudotime",'gene')
  p <- ggplot(df, aes(x=Pseudotime, 
                      y=gene)) + 
    geom_point(aes(color=celltype), size=0.8)+
    geom_smooth(method = "loess",level = 0.95,
                formula = y~x, color='black',se=F)+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black",size = 1),
          axis.title = element_blank(),
          axis.text  = element_text(size = 15, colour = 'black'),
          legend.position = 'none')+
    ggtitle(features[i])
  plist[[i]] <- p
}

library(Seurat)
CombinePlots(plist, ncol = 2)



#===============================================================================

#山脊图
library(monocle)
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)

test=readRDS("test_monocle.rds")
plotdf=pData(test)

ggplot(plotdf, aes(x=Pseudotime,y=celltype,fill=celltype))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )
ggsave("tmp1.pdf",width = 13,height = 7,units = "cm")

ggplot(plotdf, aes(x=Pseudotime,y=celltype,fill = stat(x))) +
  geom_density_ridges_gradient(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_fill_gradientn(name="Pseudotime",colors = colorRampPalette(rev(brewer.pal(10, "Spectral")))(99))+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )
ggsave("tmp2.pdf",width = 13,height = 7,units = "cm")









