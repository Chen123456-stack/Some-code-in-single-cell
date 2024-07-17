#细胞比例统计分析

#单个柱状图
setwd('./CD8/')
table(sc_pbmc_int2$cell_type,sc_pbmc_int2$group)#查看各组细胞数

sc_pbmc_int2 <- subset(sc_pbmc_int2,group%in%c('LI_NAR','PB_NAR')&cell_type%in%c('3:CD8_C1_GZMK', '4:CD8_C2_RGS1','5:CD8_C3_GZMH','6:CD8_C4_LEF1',
                                                                                 '7:CD8_C5_KLRB1','8:CD8_C6_TCF7','9:CD8_C7_FOSB','10:CD8_C8_MKI67',
                                                                                 '11:CD8_C9_TIGIT'))
sc_pbmc_int2$cell_type <- droplevels(sc_pbmc_int2$cell_type)

Idents(sc_pbmc_int2) <- sc_pbmc_int2$cell_type
table(sc_pbmc_int2$orig.ident)#查看各组细胞数
table(sc_pbmc_int2$cell_type,sc_pbmc_int2$group)#查看各组细胞数
prop.table(table(Idents(sc_pbmc_int2)))
table(Idents(sc_pbmc_int2), sc_pbmc_int2$orig.ident)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(sc_pbmc_int2), sc_pbmc_int2$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- data.frame(Cellratio)
library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]
cellper$orig.ident <- rownames(cellper)

###添加分组信息
meta <- combine(meta1,meta2)
cellper$type <- meta[rownames(cellper),'type']#R添加列
cellper$group <- meta[rownames(cellper),'group']#R添加列
## 将数据框转换为长格式
df_long <- melt(cellper, id.vars = c("group",'type','orig.ident'), 
                variable.name = "cell_type", value.name = "Proportion")
# 将 Proportion 列乘以 100
df_long$Proportion <- df_long$Proportion * 100

list = list(c("LIVER","PBMC"))
cell_type_data <- subset(df_long, cell_type == "11:CD8_C9_TIGIT")
cell_type_data$orig <- str_extract(cell_type_data$orig.ident, "^[^_]+")
cell_type_data$type <- factor(cell_type_data$type, levels = c("LIVER", "PBMC"))

p <- ggplot(cell_type_data,aes(type,Proportion,color=type, fill=type))+
  geom_boxplot(aes(fill=type),
               alpha=0.1, width=0.2)+
  geom_jitter(position = position_jitterdodge(jitter.height=0.01, # 散点抖动高度
                                              jitter.width = 0.001, # 散点抖动宽度
                                              dodge.width = 0.75),shape=1,size=3)+ # x轴方向上的闪避量
  scale_color_manual(values = c('#20466D', '#4C8562','#745EA9','#FABB83'))+
  scale_fill_manual(values = c('#20466D', '#4C8562','#745EA9','#FABB83'))+
  theme_bw()+
  labs(x = "", y = "Percentage in CD8 T cell") +
  scale_y_continuous() +  # 去掉轴的扩展，使数据紧贴轴线
  theme(panel.grid = element_blank(),  # 去掉网格线
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(),
        legend.position = "none",
        axis.line.x = element_line(),  # 保留X轴线
        axis.line.y = element_line(),  # 保留Y轴线
        axis.title = element_text( size = 10),
        panel.border = element_blank(),  # 去掉面板边框
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10)) +  # 去掉面板背景
  #geom_line(aes(group = orig), position = position_dodge(0.00001),color="grey80", linewidth = 0.9) + # 添加连接线，表示同一时间点的配对数据
  stat_compare_means(comparisons = list,label="p.signif",size=5,method = 'wilcox.test',show.legend = T,label.x=1.5)+
  ggtitle("11:CD8_C9_TIGIT") 
p


ggsave(filename = '细胞比例_11:CD8_C9_TIGIT_NAR.pdf',height = 4.5,width = 4)


