results_list_li <- list()#创建空列表
table(sc_pbmc_int2$cell_type,sc_pbmc_int2$group)
#通过for循环实现
for (i in c('3:CD8_C1_GZMK', '4:CD8_C2_RGS1','5:CD8_C3_GZMH','6:CD8_C4_LEF1',
            '7:CD8_C5_KLRB1','8:CD8_C6_TCF7','9:CD8_C7_FOSB','10:CD8_C8_MKI67',
            '11:CD8_C9_TIGIT'
)) {
  # 寻找差异基因
  group_degs <- FindMarkers(
    subset(sc_pbmc_int2, cell_type == i), 
    logfc.threshold = 0.25,
    min.pct = 0.1,
    ident.1 = "LI_NAR",
    ident.2 = "PB_NAR",
    group.by = "group"
  ) %>% 
    mutate(gene = rownames(.), cluster = i)
  # 将结果添加到列表中
  results_list_li[[i]] <- group_degs
}
# 将列表中的数据框合并为一个数据框
final_result <- bind_rows(results_list_li)


cell.markers <- final_result
cell.markers$change <- ifelse(cell.markers$avg_log2FC > 0, "UP", "DOWN")

# 筛选差异基因
lsxx = cell.markers %>% 
  #filter( cluster == "7:CD8_C5_KLRB1" ) %>% #筛选某一亚群差异基因
  filter( p_val_adj < 0.05 ) %>% #筛选p
  filter( avg_log2FC >  0.585 | avg_log2FC < -0.585 )%>%   #0.585
  arrange( desc( avg_log2FC ) ) #对avg_log2FC进行降序排序

#分别筛选出上调和下调的gene
gene_up = subset(lsxx, change == 'UP')
gene_down = subset(lsxx, change == 'DOWN')

head(gene_up)


##升高基因做富集分析
gene_clusters <- split(gene_up$gene, gene_up$cluster)

# 定义一个函数进行GO富集分析
perform_enrichGO <- function(gene_list) {
  enrichGO(gene         = gene_list,
           OrgDb        = "org.Hs.eg.db",
           keyType      = "SYMBOL",  # 使用基因符号进行富集分析
           ont          = "BP",      # BP: Biological Process, CC: Cellular Component, MF: Molecular Function
           pAdjustMethod = "BH",     # 使用 Benjamini-Hochberg 方法进行p值校正
           pvalueCutoff  = 0.05,     # p值阈值
           qvalueCutoff  = 0.05)     # q值阈值
}

# 对每个基因集进行富集分析
go_results_up <- lapply(gene_clusters, perform_enrichGO)



##降低基因做富集分析
gene_clusters <- split(gene_down$gene, gene_down$cluster)

go_results_down <- lapply(gene_clusters, perform_enrichGO)


##过滤
filter_and_label <- function(enrich_result, label) {
  result_df <- enrich_result@result %>%
    filter(p.adjust < 0.05) %>%
    mutate(change = label)
  return(result_df)
}

# 过滤并添加标签
go_results_up_filtered <- lapply(go_results_up, filter_and_label, label = "up")
go_results_down_filtered <- lapply(go_results_down, filter_and_label, label = "down")

##合并富集分析结果

merge_results <- function(cluster_name) {
  up_result <- go_results_up_filtered[[cluster_name]]
  down_result <- go_results_down_filtered[[cluster_name]]
  combined_result <- rbind(up_result, down_result)
  return(combined_result)
}

# 获取所有cluster的名字
clusters <- intersect(names(go_results_up_filtered), names(go_results_down_filtered))

# 合并结果
combined_go_results <- lapply(clusters, merge_results)
names(combined_go_results) <- clusters

# 将所有cluster的结果合并到一个数据框中
all_combined_results <- do.call(rbind, lapply(names(combined_go_results), function(cluster_name) {
  result_df <- combined_go_results[[cluster_name]]
  result_df$cluster <- cluster_name
  return(result_df)
}))



# 选择每个cluster上调和下调的前10个通路
select_top_shared_terms <- function(df, top_n = 10) {
  df %>%
    group_by(cluster, change) %>%
    arrange(cluster, change, p.adjust) %>%
    slice(1:top_n) %>%
    ungroup()
}

top_shared_terms_per_cluster <- select_top_shared_terms(all_combined_results, top_n = 10)

table(top_shared_terms_per_cluster$cluster)


select_shared_terms_across_clusters <- function(df) {
  up_terms <- df %>% filter(change == "up") %>% select(cluster, Description)
  down_terms <- df %>% filter(change == "down") %>% select(cluster, Description)
  
  # 获取上调通路的交集
  shared_up_terms <- up_terms %>%
    group_by(Description) %>%
    filter(n_distinct(cluster) > 1) %>%
    ungroup()
  
  # 获取下调通路的交集
  shared_down_terms <- down_terms %>%
    group_by(Description) %>%
    filter(n_distinct(cluster) > 1) %>%
    ungroup()
  
  # 合并上调和下调的共享通路
  shared_terms <- bind_rows(shared_up_terms, shared_down_terms)
  
  return(shared_terms)
}

# 获取不同cluster之间的相同通路
shared_terms_across_clusters <- select_shared_terms_across_clusters(top_shared_terms_per_cluster)

# 将共享的通路与原始数据合并，保留p值和其他信息
shared_terms_detailed <- inner_join(shared_terms_across_clusters, all_combined_results, 
                                    by = c("cluster", "Description"))

shared_terms_detailed$'-log10pvalue' <- ifelse(shared_terms_detailed$change == 'up', -log10(shared_terms_detailed$pvalue), -(-log10(shared_terms_detailed$pvalue)))







# 选择每个cluster上调和下调的前10个通路
select_top_terms <- function(df, top_n = 10) {
  df_up <- df %>% filter(change == "up") %>% head(top_n)
  df_down <- df %>% filter(change == "down") %>% head(top_n)
  return(rbind(df_up, df_down))
}

top_terms_per_cluster <- lapply(names(combined_go_results), function(cluster_name) {
  top_terms <- select_top_terms(combined_go_results[[cluster_name]])
  top_terms$cluster <- cluster_name
  return(top_terms)
})

top_terms_df <- do.call(rbind, top_terms_per_cluster)

top_terms_df$'-log10pvalue' <- ifelse(top_terms_df$change == 'up', -log10(top_terms_df$pvalue), -(-log10(top_terms_df$pvalue)))



plot_data$Roe <- round(plot_data$Roe, 2)
ggplot(shared_terms_detailed,aes(cluster,Description))+
  geom_point(aes(fill=`-log10pvalue`),shape=22, color='white',size=14)+
  scale_fill_gradientn(colours =colorRampPalette(c("#3162A0","white","#C7252C"))(40))+
  labs(x = NULL,y = NULL,fill="Enrichment(Ro/e)")+
  guides(fill=guide_colorbar(barheight = 8))+
  theme_bw()+  theme(panel.grid.major= element_blank(),
                     panel.grid.minor= element_blank())+
  theme(axis.text=element_text(colour='black',size=10))+
  theme(axis.text.x =element_text(angle =0,hjust =0.9,vjust = 0.9,size=10))
# geom_text(aes(label = `-log10pvalue`), color = "black", size = 3, vjust = 0.3, hjust = 0.5)
ggsave(filename = 'Roe.pdf',width = 5.5,height = 8)



# 绘制dotplot
dotplot <- ggplot(shared_terms_detailed, aes(x = cluster, y = Description)) +
  geom_point(aes(size = Count, color = `-log10pvalue`)) +
  # facet_wrap(~ change) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "-log10(p-value)") +
  theme_minimal() +
  labs(title = "GO Enrichment Dotplot",
       x = "Cluster",
       y = "GO Term",
       color = "Adjusted p-value",
       size = "Gene Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 显示dotplot
print(dotplot)

# 保存dotplot
ggsave("GO_enrichment_dotplot.png", plot = dotplot, width = 12, height = 8)








