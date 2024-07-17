
library(dplyr)
library(monocle)
library(ggplot2)
library(Seurat)
setwd("D:/KS项目/公众号文章/拟时基因趋势_通路趋势")

#===============================================================================
#                           基因随拟时表达变化
#===============================================================================
colnames(pData(cds))#mouse_monocle是已经做好拟时的monocle2对象
#将需要展示的基因表达量添加到cds， 这里我使用的是log2处理，也可以使用log10进行标准化
genes <- c("TIGHT",'CTLA4','HAVCR2')
genes_exp <- list()


# "GZMA",'GZMB','PRF1','NKG7','KLRG1','GNLY','IFNG'

for(i in 1:length(genes)){
  A <- log2(exprs(cds)[genes[i],]+1)
  A <- as.data.frame(A)
  genes_exp[[i]] <- A
}

gene_exp <- do.call(cbind, genes_exp)
colnames(gene_exp) <- genes
#将上述几个基因的拟时表达添加到monocle
pData(cds) = cbind(pData(cds), gene_exp)


#提取作图数据，只需要游基因表达和拟时即可
data <- pData(cds)
colnames(data)
#选择需要的列即可，我这里的origin.ident就是分组
data<-data %>% select('cell_type',"group",'type',"Pseudotime","TIGIT",'CTLA4','HAVCR2')


##提取非排斥患者
data <- subset(data, group %in% c("LI_AR", "LI_NAR"))
table(data$cell_type)

data <- subset(data, cell_type!='10:CD8_C8_MKI67')

data$cell_type <- droplevels(data$cell_type)

table(data$cell_type)
#ggplot作图

library(data.table)
# 将 data 转换为 data.table
data <- as.data.table(data)

# 使用 data.table 的 melt 方法
data_long_m <- melt(data, id.vars = c("group", "Pseudotime"),
                    measure.vars = 5:7,
                    variable.name = "gene",
                    value.name = "value")

# #使用分屏，应该就是文献种的办法
# #首先将data宽数据转化为长数据
# data_long_m<-melt(data, id.vars = c("type", "Pseudotime"), #需保留的不参与聚合的变量列名
#                   measure.vars = 4:6,#选择需要转化的列
#                   variable.name = c('gene'),#聚合变量的新列名
#                   value.name = 'value')#聚合值的新列名
colnames(data_long_m)



ggplot(data_long_m, aes(x=Pseudotime, y=value, color=group))+
  geom_smooth(aes(fill=group))+ #平滑的填充
  xlab('pseudotime') + 
  ylab('Relative expression') +
  facet_wrap(~gene, scales = "free_y")+ #分面，y轴用各自数据
  theme(axis.text = element_text(color = 'black',size = 12),
        axis.title = element_text(color = 'black',size = 14),
        strip.text = element_text(color = 'black',size = 14))+ #分面标题
  scale_color_manual(name=NULL, values = c("#089E86","#3D5387"))+#修改颜色
  scale_fill_manual(name=NULL, values = c("#089E86","#3D5387"))#修改颜色

ggsave("Exhuasted.pdf",  dpi = 600, width = 12, height =4,)

#===============================================================================
#                           通路随拟时表达变化
#===============================================================================
library(dplyr)
library(monocle)
library(ggplot2)
library(Seurat)
library(fgsea)
library(nichenetr)
table(seurat$group,seurat$cell_type)
dim(seurat)
mouse_data <- load(./mouse_data.rds)
#将拟时添加到seurat对象
mouse_sc <- seurat[, colnames(cds@reducedDimS)]
mouse_sc$pseudotime <- cds@phenoData@data$Pseudotime

#文章里面使用的是代谢通路，当让实际情况中并不一定是代谢通路
#所有通路基因的gmt文件可以去GSEA官网下载，这里我们以代谢通路为例子

pathway_file <- "~/Migsdb/CELL_PROLIFERATION_GO_0008283.v2023.2.Hs.gmt"
pathways <- gmtPathways(pathway_file)
pathway_names <- names(pathways)
#选择需要的通路查看基因（自己关注的或者感兴趣的通路）
#这里我随便找两个举例子
Glutathione <- pathways[["Glutathione metabolism"]]
Caffeine <- pathways[["Caffeine metabolism"]]


# #我用的是鼠的示例数据，所以需要将人的代谢基因转化为鼠的。如果是人的
# #就不必操作了
# Glutathione <- as.data.frame(Glutathione)
# Glutathione$gene = Glutathione$Glutathione %>% convert_human_to_mouse_symbols()
# #去除NA
# Glutathione <- na.omit(Glutathione)
# Glutathione_gene <- Glutathione$gene
# 
# 
# Caffeine <- as.data.frame(Caffeine)
# Caffeine$gene = Caffeine$Caffeine %>% convert_human_to_mouse_symbols()
# #去除NA
# Caffeine <- na.omit(Caffeine)
# Caffeine_gene <- Caffeine$gene

#获取评分基因表
CD8_data <- read_csv("CD81.csv")
marker.list = list()
for (col_name in colnames(CD8_data)) {
  genes <- CD8_data %>% pull(!!sym(col_name)) %>% na.omit()
  marker.list[[col_name]] <- genes
}

marker.list <- marker.list[1:19]

#计算代谢通路评分，添加到seurat对象
DefaultAssay(mouse_sc) <- "RNA"
mouse_sc <- AddModuleScore(mouse_sc,
                           features=pathways,
                           pool = rownames(mouse_sc), k=F, nbin=24,
                           name = c('Proliferation'))


mouse_sc <-  AddModuleScore(mouse_sc,
                           features = marker.list,
                           ctrl = 5,
                           name = "FunctionScore")
for(i in 1:length(marker.list)){
  colnames(mouse_sc@meta.data)[colnames(mouse_sc@meta.data) == paste0("FunctionScore", i)] <- names(marker.list)[i]
}

#提取作图数据
data1 <- mouse_sc@meta.data
colnames(data1)
data1 <- data1[, c("group",'type', "pseudotime",'Proliferation1', "Naive" ,"Activation/Effector function", "Exhaustion","TCR Signaling" ,
                   "Cytotoxicity" ,"Cytokine/Cytokine receptor","Chemokine/Chemokine receptor", "Senescence"      )]



data1 <- subset(data1, group %in% c("LI_NAR", "PB_NAR"))

data1_long   <-melt(data1, id.vars = c("type", "pseudotime"), #需保留的不参与聚合的变量列名
                    measure.vars = 4:6,#选择需要转化的列
                    variable.name = c('pathway'),#聚合变量的新列名
                    value.name = 'value')#聚合值的新列名


data1_new <- data1_long
levels(data1_new$pathway) <- c("Naive1","Effector2",'Exhausted3')


#作图和上述一样
ggplot(data1_new, aes(x=pseudotime, y=value, color=type))+
  geom_smooth(aes(fill=type))+ #平滑的填充
  xlab('pseudotime') + 
  ylab('Relative expression') +
  facet_wrap(~levels(pathway), scales = "free_y")+ #分面，y轴用各自数据
  theme(axis.text = element_text(color = 'black',size = 12),
        axis.title = element_text(color = 'black',size = 12),
        strip.text = element_text(color = 'black',size = 14))+ #分面标题
  scale_color_manual(name=NULL, values = c("#089E86","#3D5387"))+#修改颜色
  scale_fill_manual(name=NULL, values = c("#089E86","#3D5387"))#修改颜色


####将基因集打分得情况应用到拟时序分析得图中

#提取拟时序数据
t <- pData(cds)
#提取细胞坐标位置
plotdf2=as.data.frame(t(cds@reducedDimS))

colnames(plotdf2)=c("component1","component2")
#确定两个dataframe是否相同得rownames
identical(rownames(t), rownames(plotdf2))

#合并
x=combine(t,plotdf2)
colnames(x)

#data1包括打分信息，确定两个dataframe得rownames是否相同
identical(rownames(x), rownames(data1))

t =combine(x,data1)
colnames(t)

t <- subset(t, group %in% c("LI_NAR", "PB_NAR"))

x_rotated <- t
x_rotated$component1 <- -t$component1
x_rotated$component2 <- -t$component2

#合并完之后进行画图
colnames(x_rotated)[colnames(x_rotated) == "Proliferation1"] <- "Proliferation"

my.breaks <- c(seq(0, 0.15, by=0.01), seq(0.5, 0.5, by=0.0001))
my.colors <- c(
  colorRampPalette(colors = c("#E0DDD7", "#A24F1D"))(length(my.breaks)/2),
  colorRampPalette(colors = c("#A24F1D", "#A24F1D"))(length(my.breaks)/2))

ggplot(data = x_rotated, aes(x = component1, y = component2, color = `Proliferation`)) +
  geom_point(size = 0.5) +
  scale_color_gradientn(colors = colorRampPalette(c("#E0DDD7", "#A24F1D"))(3)) +
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
ggsave(filename ='Proliferation.pdf',width = 8,height = 6,dpi = 300)






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

p <- ggplot(data = t, aes(x = component1, y = component2, color = Proliferation1)) +
  geom_point(size = 2) +
  scale_color_gradientn(colors = colorRampPalette(c( "#7B2063","green"))(100)) +
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






