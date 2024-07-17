
Idents(sc_pbmc_int2) <- sc_pbmc_int2$cell_type
t <- table(Idents(sc_pbmc_int2), sc_pbmc_int2$group)
t <- as.data.frame(t)
colnames(t) <- c('cell_type', 'patient','count')
t <- t %>% filter(count != 0)
df <- t

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

# 查看结果
df_proportions



# 重新整理数据以适应绘图需求
df_plot <- df_proportions %>%
  mutate(patient = factor(patient, levels = unique(patient))) %>%
  mutate(cell_category = factor(cell_category, levels = c("CD4", "CD8", "NK")))



head(df_plot)


ggplot(df_plot, aes(x = interaction(patient, cell_category), y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_discrete(labels = function(x) sub("\\..*", "", x)) +
  facet_wrap(~ patient, scales = "free_x", nrow = 1) +
  labs(x = "Patient", y = "Proportion", title = "") +
  scale_fill_manual(values = c(FJ)) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(), # 去掉X轴标签
    axis.text.y = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggplot(df_plot, aes(x = interaction(patient, cell_category), y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", color = "black",size=0.9) +  # 增加黑框
  scale_x_discrete(labels = function(x) sub("\\..*", "", x)) +
  facet_wrap(~ patient, scales = "free_x", nrow = 1) +
  labs(x = "Patient", y = "Proportion", title = "") +
  scale_fill_manual(values = c(FJ)) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(), # 去掉X轴标签
    axis.text.y = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave( filename = '细胞比例柱状图.pdf', width = 10, height = 8 )


FJ <- c('#4C8562','#A24F1D','#65A8AE','#363C6C','#3D5C96','#5E5598','#67912F',
        '#806224','#D08685','#A7256F','#834626','#6376B9',
        '#B5A0C6','#4A277B','#BF9B25','#A72023','#A3BCD3','#C66E26','#528A31',
        '#535353')

mycolor <- c('#2BA077','#745EA9','#3084BB','#AA6168','#9DD097',
             '#845550','#CEB8B5','#B4C380','#CE7BB3','#C1B2D4',
             '#DA3235','#B0C9E6','#FCC4C4','#A7DEDB','#D6D995',
             '#F78838','#F4B2CC','#FABB83','#23C0CE','#B5A0C6','#4A277B')



t <- c("#F08080", "#228B22","#C71585", "#FF7F0E", "#FFBB78", "#2CA02C", "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
  "#8C564B","#C49C94", "#E377C2", "#F7B6D2")  # Example colors for CD4 cell types
colors_CD8 <- c()  # Example colors for CD8 cell types
colors_NK <- c()  

##用于生成图例

# Install and load required packages
library(ggplot2)
library(gridExtra)

# Create separate data frames for each cell category
df_plot_CD4 <- df_plot %>%
  filter(cell_category == "CD4")
df_plot_CD8 <- df_plot %>%
  filter(cell_category == "CD8")
df_plot_NK <- df_plot %>%
  filter(cell_category == "NK")

# Define color palettes for each cell category
colors_CD4 <- c('#4C8562','#A24F1D')  # Example colors for CD4 cell types
colors_CD8 <- c('#65A8AE','#363C6C','#3D5C96','#5E5598','#67912F','#806224','#D08685','#A7256F','#834626')  # Example colors for CD8 cell types
colors_NK <- c('#6376B9','#B5A0C6','#4A277B')  # Example colors for NK cell types

# Create plots for each cell category
plot_CD4 <- ggplot(df_plot_CD4, aes(x = interaction(patient, cell_type), y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", color = "black",size=0.9) +
  scale_x_discrete(labels = function(x) sub("\\..*", "", x)) +
  facet_wrap(~ patient, scales = "free_x", nrow = 1) +
  labs(x = "Patient", y = "Proportion", title = "CD4 Cell Type Proportions") +
  scale_fill_manual(values = colors_CD4) +  # Colors for CD4 cell types
  guides(fill = guide_legend(title = "CD4 T cell", ncol = 2)) +  # Adjust legend position and title
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_hline(yintercept = 1, linetype = "dashed")

plot_CD8 <- ggplot(df_plot_CD8, aes(x = interaction(patient, cell_type), y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", color = "black",size=0.9) +
  scale_x_discrete(labels = function(x) sub("\\..*", "", x)) +
  facet_wrap(~ patient, scales = "free_x", nrow = 1) +
  labs(x = "Patient", y = "Proportion", title = "CD8 Cell Type Proportions") +
  scale_fill_manual(values = colors_CD8) +  # Colors for CD8 cell types
  guides(fill = guide_legend(title = "CD8 T cell", ncol = 2,direction = "vertical" )) +  # Adjust legend position and title
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_hline(yintercept = 1, linetype = "dashed")



plot_NK <- ggplot(df_plot_NK, aes(x = interaction(patient, cell_type), y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack", color = "black",size=0.9) +
  scale_x_discrete(labels = function(x) sub("\\..*", "", x)) +
  facet_wrap(~ patient, scales = "free_x", nrow = 1) +
  labs(x = "Patient", y = "Proportion", title = "CD8 Cell Type Proportions") +
  scale_fill_manual(values = colors_NK) +  # Colors for CD8 cell types
  guides(fill = guide_legend(title = "NK", ncol = 2)) +  # Adjust legend position and title
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  geom_hline(yintercept = 1, linetype = "dashed")

                                
p = plot_CD4+plot_CD8+plot_NK + plot_layout( guides='collect',ncol = 1 ) # 拼图 patchwork
p
ggsave( p, filename = '图例.pdf', width = 10, height = 20 )

# "1:CD4_C1_LEF1" = "#4C8562", "2:CD4_C2_IL7R" = "#A24F1D", 
# "3:CD8_C1_GZMK" = "#A72023", "4:CD8_C2_RGS1" = "#3D5C96", "5:CD8_C3_GZMH" = "#5E5598",
# "6:CD8_C4_LEF1" = "#C66E26", "7:CD8_C5_KLRB1" = "#535353", "8:CD8_C6_TCF7" = "#806224",
# "9:CD8_C7_FOSB" = "#D08685", "10:CD8_C8_MKI67" = "#A3BCD3", "11:CD8_C9_TIGHT" = "#528A31",
# "12:NK_C1_FCGR3A" = "#927c9a", "13:NK_C2_XCL1" = "#3674a2", "14:NK_C3_KLRF1" = "#A3BCD3"

