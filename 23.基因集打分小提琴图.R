setwd('./CD8')

CD8_data <- read_csv("CD4.csv")
marker.list = list()
for (col_name in colnames(CD8_data)) {
  genes <- CD8_data %>% pull(!!sym(col_name)) %>% na.omit()
  marker.list[[col_name]] <- genes
}

# marker.list <- marker.list[1:18]

table(sc_pbmc_int2$cell_type)
# 子集化包含感兴趣的 CD8 细胞亚群
CD8_Obj <- subset(sc_pbmc_int2, cell_type %in% c('1:CD4_C1_LEF1', '2:CD4_C2_IL7R'))

# 去除没有数据的细胞类型
CD8_Obj$cell_type <- droplevels(CD8_Obj$cell_type)

# 检查每个细胞类型的细胞数
table(CD8_Obj$cell_type)



CD8_Obj <-  AddModuleScore(CD8_Obj,
                           features = marker.list,
                           ctrl = 5,
                           name = "FunctionScore")
for(i in 1:length(marker.list)){
  colnames(CD8_Obj@meta.data)[colnames(CD8_Obj@meta.data) == paste0("FunctionScore", i)] <- names(marker.list)[i]
}

colnames(CD8_Obj@meta.data)

CD8_Obj <- subset(CD8_Obj, group %in% c('LI_NAR','PB_NAR'))

subCD8_Obj <- subset(CD8_Obj, cell_type %in% c('1:CD4_C1_LEF1'))

# 最关键的作图
Idents(subCD8_Obj) <- subCD8_Obj$group
my_comparisons <- list( c("LI_NAR", "PB_NAR"))

VlnPlot(subCD8_Obj, features= "Chemokine/Chemokine receptor" , group.by = "group", cols = c('#20466D', '#4C8562'), pt.size = 0) +  theme_classic() + 
  theme(text = element_text(size=20, colour = "black")) + RotatedAxis() + 
  geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.size = 0)+
 scale_y_continuous(limits = c(-0.25,1))+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank())+ 
  labs(title = "", y = "Chemokine/Chemokine receptor", x="") + theme(legend.position="right") +  
  # stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1),  geom = "pointrange", color = "black")+
  scale_fill_manual(values = c('#20466D', '#4C8562'), labels = c("LIVER", "PBMC"))+
  stat_compare_means(comparisons = my_comparisons,label="p.signif",size=5,method = 'wilcox.test',show.legend = T,label.x=1.5)

ggsave(filename = '1:CD4_C1_LEF1_Chemokine_Chemokine receptor.pdf',width = 6,height = 6,dpi = 300)







library(Seurat)
library(ggplot2)
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )
  return(p)
}

## main function
StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

modify_vlnplot(sc_pbmc_int2, c('CD163'), pt.size=0, cols=my36colors)

#配色方案
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')











