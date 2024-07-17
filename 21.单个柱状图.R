library(tidyverse)
library(gapminder)
library(ggpubr)
library(rstatix)
library(ggprism)


#算法有错误，需要在调整
table(sc_pbmc_int2$cell_type)
t <- table(sc_pbmc_int2$orig.ident, sc_pbmc_int2$cell_type, sc_pbmc_int2$group)
t <- as.data.frame(t)
colnames(t) <- c('orig.ident','cell_type', 'patient','count')
t_filtered <- subset(t, count != 0)
df <- t_filtered

# 添加细胞大类列
df <- df %>%
  mutate(cell_category = case_when(
    grepl("CD4", cell_type) ~ "CD4",
    grepl("CD8", cell_type) ~ "CD8",
    grepl("NK", cell_type) ~ "NK",
    TRUE ~ "Other"
  ))

# 只保留 CD4, CD8 和 NK 数据
df <- df %>%
  filter(cell_category %in% c("CD4", "CD8", "NK"))

# 计算每个大类中的总数
total_counts <- df %>%
  group_by(patient, cell_category) %>%
  summarise(total_count = sum(count))

# 计算每个亚群占对应大类总数的比例
df_proportions <- df %>%
  left_join(total_counts, by = c("patient", "cell_category")) %>%
  mutate(proportion = count / total_count)

df_proportions$proportion <- df_proportions$proportion*100



# 查看结果
df_proportions

table(df_proportions$cell_type)
list = list(c("LI_AR","PB_AR"),c('LI_NAR','PB_NAR'))

cell_type_data <- subset(df_proportions, cell_type == "3:CD8_C1_GZMK")
cell_type_data <- subset(cell_type_data, !grepl("^LJH_L", orig.ident))


cell_type_data$orig <- str_extract(cell_type_data$orig.ident, "^[^_]+")

cell_type_data$patient <- factor(cell_type_data$patient, levels = c("LI_AR", "PB_AR", "LI_NAR", "PB_NAR"))

p <- ggplot(cell_type_data,aes(patient,proportion,color=patient, fill=patient))+
  geom_boxplot(aes(fill=patient),
               alpha=0.1)+
  geom_jitter(position = position_jitterdodge(jitter.height=0.01, # 散点抖动高度
                                              jitter.width = 0.001, # 散点抖动宽度
                                              dodge.width = 0.75),shape=1,size=3)+ # x轴方向上的闪避量
  scale_color_manual(values = c('#20466D', '#4C8562','#745EA9','#FABB83'))+
  scale_fill_manual(values = c('#20466D', '#4C8562','#745EA9','#FABB83'))+
  theme_bw()+
  labs(x = "", y = "Percentage in CD8 T cell") +
  scale_y_continuous() +  # 去掉轴的扩展，使数据紧贴轴线
  theme(panel.grid = element_blank(),  # 去掉网格线
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(),
        legend.position = "none",
        axis.line.x = element_line(),  # 保留X轴线
        axis.line.y = element_line(),  # 保留Y轴线
        axis.title = element_text( size = 16),
        panel.border = element_blank(),  # 去掉面板边框
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 16)) +  # 去掉面板背景
  geom_line(aes(group = orig), position = position_dodge(0.00001),color="grey80", linewidth = 0.9) + # 添加连接线，表示同一时间点的配对数据
  stat_compare_means(comparisons = list,label="p.signif",size=5,method = 'wilcox.test',show.legend = T,label.x=1.5)+
  ggtitle("4:CD8_C2_RGS1") 
p
ggsave(filename = '4:CD8_C2_RGS1_group.pdf',width = 4,height = 5)

# aes(group = patient)

# 
# 
# mytheme <- theme_prism() + #这里用一个Graphpad Prism 风格主题
#   theme(strip.text = element_text(size = 18),
#         axis.line = element_line(color = "black",size = 0.4),
#         axis.text.y = element_text(color = "black",size = 18),
#         axis.text.x = element_text(color = "black",size = 16),
#         axis.title = element_text(color = "black",size = 20),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(size = 0.2,color = "#e5e5e5"),
#         legend.position = "none")
# df_proportions$proportion <- df_proportions$proportion * 100
# cell_type_data <- subset(df_proportions, cell_type == "1:CD4_C1_LEF1")



