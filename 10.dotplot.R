
#将cell_type和group赋值为新的一列，好进一步计算cell_type和group组合的结果
sc_pbmc_int2 <- AddMetaData(sc_pbmc_int2,metadata = paste(sc_pbmc_int2$cell_type, 
                           sc_pbmc_int2$individual, sep = "_"), col.name = "GROUP")

table(sc_pbmc_int2$GROUP)

#提取想要观察的差异基因
expr <- AverageExpression(sc_pbmc_int2, assays = "RNA",features = 'CCL3L3',group.by = 'GROUP', slot = "data")[[1]]
expr <- as.matrix(expr)
expr<- t(expr)


# 提取 cell_type 和 group
cell_type_group <- rownames(expr)
split_cell_type_group <- strsplit(cell_type_group, "-")
cell_type <- sapply(split_cell_type_group, "[", 1)
group <- sapply(split_cell_type_group, function(x) paste(x[c(2, 3)], collapse = "_"))

# 将提取出来的列添加到矩阵中
expr <- cbind(cell_type = cell_type, group = group, expr)

expr <- as.data.frame(expr)
exprs <- expr
#celltype转换为因子，固定顺序
exprs$cell_type <- factor(exprs$cell_type,
                          levels = c('T', 'NK','Neutrophil','B','Plasma cell','Basophil',
                                     'Monocyte','Macrophage','cDC1','cDC2','pDC','Hepatocyte',
                                     'Endothelial cell','Fibroblast', 'Magakaryocyte' ))

exprs$CCL3L3 <- as.numeric(exprs$CCL3L3)
p1 <- ggplot(exprs, aes(x = cell_type, y = group)) + # 建立映射
  geom_point(aes(size =CCL3L3), color = '#00bbb1') # 绘制散点，固定颜色为蓝色
p1


# 图表个性化调整:
p2 <- p1 +
  scale_size_continuous(range = c(0,8)) + # 气泡大小范围调整
  theme_classic() +
  theme(legend.text = element_text(size = 14), # 图例文本字号
        legend.title = element_text(size = 16), # 图例标题字号
        axis.text = element_text(size = 14), # 坐标轴标签字号
        axis.title = element_text(size = 16), # 坐标轴标题字号
        axis.text.x = element_text(angle = 60, hjust = 1)) + # x轴标签旋转60°
  labs(x = 'cell_type', y = 'Group', # xy轴标题修改
       size = 'mean_CCL3') + # 图例标题修改
  coord_flip() # 坐标轴翻转
p2

ggsave(filename = 'CCR1.pdf',width = 12,height = 8)



ggplot(exprs, aes(cell_type, group, size=CCL3L3, fill=CCL3L3)) + 
  geom_point(shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill=NA))) + 
  scale_size_continuous(range = c(0,8))+
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major.x = element_line(color = "grey80"),
    panel.grid.major.y = element_line(color = "grey80"),
    axis.title = element_blank(),
    axis.text.y = element_text(color='black',size=12),
    axis.text.x = element_text(color='black',size=12, angle = 90, hjust = 1, vjust = 0.5))+
  scale_fill_gradientn(colours = c('#5749a0', '#0f7ab0', '#00bbb1',
                                   '#bef0b0', '#fdf4af', '#f9b64b',
                                   '#ec840e', '#ca443d', '#a51a49')) + # 图例标题修改
  coord_flip() # 坐标轴翻转
