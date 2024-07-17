# 检查每个细胞类型的细胞数
table(sc_pbmc_int2$cell_type)

# 去除没有数据的细胞类型
sc_pbmc_int2$cell_type <- droplevels(sc_pbmc_int2$cell_type)

# 检查每个细胞类型的细胞数
table(sc_pbmc_int2$type)



library(ggrepel)
library(ggsci)
library(tidyverse)
library(ggplot2)

#分组饼图

  # 1. Get plot data

  plot_data <- sc_pbmc_int2@meta.data %>%
    dplyr::select(type, cell_type)
  
  plot_data <- table(plot_data[['cell_type']], plot_data[['type']]) %>%
    as.data.frame() %>%
    group_by(Var2) %>%
    mutate(Total = sum(Freq),
           Proportion = round((Freq / Total)*100,2),
           labels = scales::percent(Freq / Total)) %>%
    mutate(text_y = Total - (Freq/2)) %>%
    dplyr::rename('type' = "Var2") %>%
    dplyr::rename('cell_type' = "Var1")
  
  # 3.plot
  p <- ggplot(plot_data, aes(x = 2, y = Proportion, fill = cell_type)) +
    geom_bar(width = 1, stat = "identity", color = "black", size = 0.5) +  # 添加黑色边框，设置边框宽度
    coord_polar(theta = "y", start = 0) +
    theme_void() +
    facet_wrap(~ type, strip.position = "top") +
    #labs(fill = 'cell_type') +
    xlab("") +
    scale_color_manual(values = FJ) +
    scale_fill_manual(values = FJ) +
    xlim(0.5, 2.5)  + 
    geom_label_repel(aes(label = labels,y=Proportion/2), 
                              size = 2, 
                              nudge_x = 0.6,
                              nudge_y = 0.4,
                              show.legend = FALSE, 
                              segment.color = "grey50",
                              max.overlaps = 50)

  
  p
  
  ggsave(filename = 'pie_type.pdf',width = 15,height = 10  )
  
  #不分组饼图
  plot_data <- sc_pbmc_int2@meta.data %>%
    dplyr::select( cell_type)
  
  plot_data <- table(plot_data[['cell_type']]) %>%
    as.data.frame() %>%
    mutate(Total = sum(Freq),
           Proportion = round((Freq / Total)*100,2),
           labels = scales::percent(Freq / Total)) %>%
    mutate(text_y = Total - (Freq/2)) %>%
    dplyr::rename('cell_type' = "Var1")
  

  # 3.plot
  p <- ggplot(plot_data, aes(x = 2, y = Proportion, fill = cell_type)) +
    geom_bar(width = 1, stat = "identity", color = "black", size = 0.5) +  # 添加黑色边框，设置边框宽度
    coord_polar(theta = "y", start = 0) +
    theme_void() +
    # facet_wrap(~ type, strip.position = "top") +
    labs(fill = 'cell_type') +
    xlab("") +
    scale_color_manual(values = FJ) +
    scale_fill_manual(values = FJ) 
  
  
    # + 
    # geom_text_repel(aes(label = labels,y=Proportion/2), 
    #                  size = 2, 
    #                  # nudge_x = 0.6,
    #                  # nudge_y = 0.4,
    #                  show.legend = FALSE, 
    #                  segment.color = "grey50",
    #                  max.overlaps = 50)
  
  
  p

  ggsave(filename = 'pie.pdf',width = 15,height = 10  )
  
  ######
  
  names <- table(sc_pbmc_int2$cell_type) %>% names()
  ratio <- table(sc_pbmc_int2$cell_type) %>% as.numeric()
  pielabel <- paste0(names," (", round(ratio/sum(ratio)*100,2), "%)")
  
  selected_labels <- c("CD4 T cell", "CD8 T cell", "NK",'B cell','Plasma cell','Myeloid cell','pDC','Neutrophil')
  # 逆时针旋转 45 度
  p5 <- pie(ratio, labels = selected_labels, radius = 1.0, clockwise = F, 
            main = "cell_type", col = FJ, start = 90)
  p5
  
  p5 <- pie(ratio, labels = selected_labels, radius = 1.0, clockwise = FALSE, 
            main = "cell_type", col = FJ, start = 90)  # 90度为从正上方开始

  
  #细胞数量
  
  # 将样本、细胞类型、细胞个数制成表格
  sample_table <- as.data.frame(table(sc_pbmc_int2@meta.data$individual,sc_pbmc_int2@meta.data$cell_type))
  # 定义列名
  names(sample_table ) <- c("Samples","celltype","CellNumber")
  # 规定样本顺序
  #sample_table$Samples<-factor(sample_table$Samples,levels = c("1NT_P21","2NT_P24","3NT_P25","1T_C21","2T_C24","3T_C25"))
  # 指定颜色
  #colors <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#1F78B4", "#33A02C")
  # 横轴为样本,纵轴为细胞比例
  p1 <-ggplot(sample_table,aes(x=Samples,weight=CellNumber,fill=celltype))+
    geom_bar(position="fill",width = 0.7,size = 0.5,colour = '#222222')+
    scale_fill_manual(values = FJ) + 
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          axis.line.x = element_line(colour = "black") ,
          axis.line.y = element_line(colour = "black") ,
          plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
    )+labs(y="Percentage")+RotatedAxis()
  p1
  
  
  p2 <- ggplot(sample_table,aes(x=Samples,weight=CellNumber,fill=celltype))+
    geom_bar(position="stack",width = 0.7,size = 0.5,colour = '#222222')+
    scale_fill_manual(values = FJ) + 
    theme_bw()+
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          axis.title.x = element_text(),
          axis.line.x = element_line(colour = "black") ,
          axis.line.y = element_line(colour = "black") ,
          plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
    )+labs(y="Cell number")+RotatedAxis()
  p2   

  ggsave(filename = 'cellnumber.pdf',width = 12,height = 6)  
  
  
  
  
  
  
  
  