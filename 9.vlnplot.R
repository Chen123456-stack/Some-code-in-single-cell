
library(Seurat)
library(tidyverse)
library(gghalves)
library(gridExtra)

table(sc_pbmc_int2$cell_type,sc_pbmc_int2$group)
sc_NK = subset(sc_pbmc_int2, cell_type %in% c('NK_rest')&type%in%c('LIVER') )
gene <- c('CCR1','CCR2','CX3CR1','XCL1','XCL2')

#考虑是否用data做
exprs <- data.frame(FetchData(object = sc_NK, vars = c("cell_type",gene,"group"),layer = 'data'))##vars对象可更改，如果想在不同的cluster上进行比对，
                                                                                   # 可以加入srueat_clusters(具体变量最好是head（combined），看一下自己的数据rowname
exprs$Proj <- "Seurat"

df <- gather(exprs, gene, expression, c(-cell_type,-group,-Proj))


##多个基因合在一张图上
p <- ggplot() +
  geom_half_violin(data = df[df$group == 'LI_AR',] %>% filter(gene=="CCR1"),
                   aes(x = gene, y = expression, fill = group),
                   color = 'black',
                   scale = 'width') +
  geom_half_violin(data = df[df$group == 'LI_NAR',]%>% filter(gene=="CCR1"),
                   aes(x = gene, y = expression, fill = group),
                   color = 'black',
                   scale = 'width',
                   side = 'r')+
  # geom_half_violin(data = df[df$group == 'PB_AR',] %>% filter(gene=="CCR2"),
  #                  aes(x = gene, y = expression, fill = group),
  #                  color = 'black',
  #                  scale = 'width') +
  # geom_half_violin(data = df[df$group == 'PB_NAR',]%>% filter(gene=="CCR2"),
  #                  aes(x = gene, y = expression, fill = group),
  #                  color = 'black',
  #                  scale = 'width',
  #                  side = 'r')+
  # geom_half_violin(data = df[df$group == 'PB_AR',] %>% filter(gene=="CX3CR1"),
  #                  aes(x = gene, y = expression, fill = group),
  #                  color = 'black',
  #                  scale = 'width') +
  # geom_half_violin(data = df[df$group == 'PB_NAR',]%>% filter(gene=="CX3CR1"),
  #                  aes(x = gene, y = expression, fill = group),
  #                  color = 'black',
  #                  scale = 'width',
  #                  side = 'r')+
  geom_half_violin(data = df[df$group == 'LI_AR',] %>% filter(gene=="XCL1"),
                   aes(x = gene, y = expression, fill = group),
                   color = 'black',
                   scale = 'width') +
  geom_half_violin(data = df[df$group == 'LI_NAR',]%>% filter(gene=="XCL1"),
                   aes(x = gene, y = expression, fill = group),
                   color = 'black',
                   scale = 'width',
                   side = 'r')+
  geom_half_violin(data = df[df$group == 'LI_AR',] %>% filter(gene=="XCL2"),
                   aes(x = gene, y = expression, fill = group),
                   color = 'black',
                   scale = 'width') +
  geom_half_violin(data = df[df$group == 'LI_NAR',]%>% filter(gene=="XCL2"),
                   aes(x = gene, y = expression, fill = group),
                   color = 'black',
                   scale = 'width',
                   side = 'r')
p
p1 <- p +
  theme_bw() +
  scale_fill_manual(values = c("#3274A1","#E1812C")) +
  labs(x = 'Gene', y = 'Expression Level')
p1

ggsave(filename = 'NKrest_PB.pdf',width = 8,height = 8)


#单个For循环绘制
# 创建空的图表列表
plot_list <- list()

# 循环替换基因并创建半小提琴图层
for (gene in c('XCL1','XCL2')) {
  # 创建半小提琴图层
  violin_layer1 <- ggplot() + geom_half_violin(data = exprs[exprs$group == "LI_AR", ],
                                               aes(x = Proj, y = !!sym(gene), fill = group),
                                               color = 'black',
                                               scale = 'width')
  violin_layer2 <- geom_half_violin(data = exprs[exprs$group == 'LI_NAR',],
                                    aes(x = Proj, y = !!sym(gene), fill = group),
                                    color = 'black',
                                    scale = 'width',
                                    side = 'r')
  # 添加图层到图表列表中
  plot_list[[gene]] <- violin_layer1 + violin_layer2 +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank()) +
    scale_fill_manual(values = c("#3274A1","#E1812C")) +
    labs(x = gene ,y = 'Expression Level')
}
p3 <- do.call(grid.arrange, c(plot_list, nrow =1, ncol = 2))

ggplot2::ggsave(filename = paste0("p3.pdf"),
                p3,
                height=15,
                width=30,
                dpi=300,
                units="cm",
                bg = "#FFFFFF")


##for循环绘制普通提琴图

plot_list <- list()
my_comparisons <- list( c( "LI_NAR","LI_AR"))
for (gene in c('XCL1','XCL2')) {
  gene_sym <- rlang::sym(gene)
  # 创建半小提琴图层
  violin_layer1 <- VlnPlot(sc_NK, features = gene, split.by = 'group', pt.size = 0)+
    stat_compare_means(comparisons = my_comparisons)+  #加统计分析
    ylim(-0.5, 6)
  # 添加图层到图表列表中
  plot_list[[gene]] <- violin_layer1 +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank()) +
    scale_fill_manual(values = c("#3274A1","#E1812C")) +
    labs(x = gene, y = 'Expression Level')
}
p3 <- do.call(grid.arrange, c(plot_list, nrow =1, ncol = 2))
ggplot2::ggsave(filename = paste0("p3.pdf"),
                p3,
                height=15,
                width=30,
                dpi=300,
                units="cm",
                bg = "#FFFFFF")

#对分组进行排序
Idents(sc_NK) <- sc_NK$group
sc_NK@meta.data$group <- factor(sc_NK@meta.data$group,
                                           levels = c('LI_AR','LI_NAR'))

VlnPlot(sc_NK, features = 'XCL1', split.by = 'group', pt.size = 0)+
  stat_compare_means(comparisons = my_comparisons)+  #加统计分析
  ylim(-1, 6)

