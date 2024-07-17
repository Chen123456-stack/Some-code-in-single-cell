#monocle中的可视化函数
plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "cell_type")

#提取数据=======================================================================
data_df <- t(reducedDimS(cds)) %>% as.data.frame() %>% #提取坐标
  select_(Component_1 = 1, Component_2 = 2) %>% #重命名
  rownames_to_column("cells") %>% #rownames命名
  mutate(pData(cds)$State) %>% #添加State
  mutate(pData(cds)$Pseudotime, 
         pData(cds)$orig.ident, 
         pData(cds)$cell_type,
         pData(cds)$subcell_type,
         pData(cds)$group)#将这些需要作图的有用信息都添加上

colnames(data_df) <- c("cells","Component_1","Component_2","State",
                       "Pseudotime","orig.ident","celltype",'subcelltype','group')
#==============================================================================
#轨迹数据提取---完全摘录于monocle包原函数
dp_mst <- minSpanningTree(cds)
reduced_dim_coords <- reducedDimK(cds)
ica_space_df <- Matrix::t(reduced_dim_coords) %>% as.data.frame() %>% 
  select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>% 
  mutate(sample_name = rownames(.), sample_state = rownames(.))

#构建一个做轨迹线图的数据
edge_df <- dp_mst %>% igraph::as_data_frame() %>% 
  select_(source = "from", target = "to") %>% 
  left_join(ica_space_df %>% select_(source = "sample_name", 
                                     source_prin_graph_dim_1 = "prin_graph_dim_1", 
                                     source_prin_graph_dim_2 = "prin_graph_dim_2"), by = "source") %>% 
  left_join(ica_space_df %>% select_(target = "sample_name", 
                                     target_prin_graph_dim_1 = "prin_graph_dim_1", 
                                     target_prin_graph_dim_2 = "prin_graph_dim_2"), by = "target")

#==============================================================================
#计算细胞比例
Cellratio <- prop.table(table(data_df$State, data_df$group), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c('State',"group","Freq")
#==============================================================================

#ggplot作图
library(ggplot2)
library(tidydr)
library(ggforce)
library(ggrastr)

g <- ggplot() + 
  geom_point_rast(data = data_df, aes(x = Component_1, 
                                      y = Component_2,
                                      color =Pseudotime)) + #散点图
  scale_color_viridis()+#密度色
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", 
                          xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), 
               linewidth = 1, 
               linetype = "solid", na.rm = TRUE, data = edge_df)+#添加轨迹线
  theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+#坐标轴主题修改
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_arc(arrow = arrow(length = unit(0.15, "inches"), #曲线箭头
                         type = "closed",angle=315),
           aes(x0=0,y0=-3,r=5, start=-0.4*pi, end=0.4*pi),lwd=1)+
  geom_arc_bar(data=subset(Cellratio,State=='1'),stat = "pie",#添加饼图
               aes(x0=-15,y0=0,r0=0,r=2.5,amount=Freq,fill=group))+
  geom_arc_bar(data=subset(Cellratio,State=='2'),stat = "pie",
               aes(x0=2,y0=9,r0=0,r=2.5,amount=Freq,fill=group))+
  geom_arc_bar(data=subset(Cellratio,State=='3'),stat = "pie",
               aes(x0=15,y0=-3,r0=0,r=2.5,amount=Freq,fill=group))+
  scale_fill_manual(values = mycolor)
g

c('#67A9CC','#DA8A87')


