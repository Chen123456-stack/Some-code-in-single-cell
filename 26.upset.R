head(final_result)

final_result$change <- ifelse(final_result$avg_log2FC > 0, "UP", "DOWN")

#分别筛选出上调和下调的gene
gene_up = subset(final_result, change == 'UP')
gene_down = subset(final_result, change == 'DOWN')

head(gene_up)

library(dplyr)
library(ggvenn)

# 创建一个列表，键为cluster，值为上调基因
upregulated_genes <- gene_up %>%
  filter(change == "UP") %>%
  group_by(cluster) %>%
  summarize(genes = list(gene))

downregulated_genes <- gene_down %>%
  filter(change == "DOWN") %>%
  group_by(cluster) %>%
  summarize(genes = list(gene))

# 转换为一个命名列表
upregulated_genes_list <- setNames(upregulated_genes$genes, upregulated_genes$cluster)

upregulated_genes_list <- upregulated_genes_list[c(3,4,5,7,9)]

downregulated_genes_list <- setNames(downregulated_genes$genes, downregulated_genes$cluster)

downregulated_genes_list <- downregulated_genes_list[c(3,4,5,7,9)]


example <- fromList(upregulated_genes_list)#Upset 自带函数转化数据结构  
example <- fromList(downregulated_genes_list)#Upset 自带函数转化数据结构  

setsBarColors <-  c('#4C8562','#A24F1D','#65A8AE')

dev.off()
pdf("naive_down_veen.pdf", width = 10, height = 8)
 upset(example,
      nsets=length(example),
      sets = c('3:CD8_C1_GZMK', '4:CD8_C2_RGS1','5:CD8_C3_GZMH',
               '7:CD8_C5_KLRB1','9:CD8_C7_FOSB'),
      nintersects=30,
      keep.order = TRUE,
      point.size = 3,
      line.size = 1,
      number.angles = 0,
      text.scale = c(1.5, 1.2, 1.2, 1.5, 1.5, 1.25), 
      order.by="freq",
      matrix.color="#4285F4",
      main.bar.color = 'black',
      sets.bar.color=setsBarColors,
      queries = list(list(query=intersects,params=list('3:CD8_C1_GZMK', '4:CD8_C2_RGS1','5:CD8_C3_GZMH',
                                                       '7:CD8_C5_KLRB1','9:CD8_C7_FOSB'),color='#A72023',active=T)))
dev.off()



# 获取上调基因的交集列表
upregulated_genes_intersection <- Reduce(intersect, upregulated_genes_list)


go_up <- enrichGO(OrgDb="org.Hs.eg.db", # 种属基因数据
                  gene = upregulated_genes_intersection, 
                  keyType='SYMBOL',
                  ont ="BP", # CC BP MF
                  pvalueCutoff = 1, # 设置校准P值筛选阈值
                  qvalueCutoff = 1, # 设置Q值筛选阈值
                  readable = TRUE # 自动将 ID 转换为 SYMBOL
) %>% 
  # clusterProfiler::simplify( . ) %>% # 删减重复条目，可选，相似结果比较多时好用
  enrichplot::pairwise_termsim() # emap 网络图需要做此转换




# 获取下调基因的交集列表
downregulated_genes_intersection <- Reduce(intersect, downregulated_genes_list)

#下调基因富集
go_down <- enrichGO(OrgDb="org.Hs.eg.db", # 种属基因数据
                    gene = downregulated_genes_intersection, 
                    keyType='SYMBOL',
                    ont ="BP", # CC BP MF
                    pvalueCutoff = 1, 
                    qvalueCutoff = 1, 
                    readable = TRUE # 自动将 ID 转换为 SYMBOL
) %>% 
  # clusterProfiler::simplify( . )%>% # 删减重复条目，可选，相似结果比较多时好用
  enrichplot::pairwise_termsim() # emap 网络图需要做此转换


#筛选p.adjust<0.05通路
go_up <- go_up@result %>%
  filter(p.adjust<0.05) %>%
  mutate(change='up')

go_down <- go_down@result %>%
  filter(p.adjust<0.5) %>%
  mutate(change='down')

# go_up_sig <- go_up[c(1:20),]
# go_down_sig <- go_down[c(1:20),]



go_up_sig$'-log10pvalue' <-  -log10(go_up_sig$pvalue)
df$Description <- str_wrap(df$Description, width =20)


#以富集结果表Top20为例：
go_up <- go_down[c(1,3:10,12,13,14,16,17,18,19,21,22,24,26),]
go_up <- go_up %>%
  arrange(desc(Count))


go_up$Description <- str_wrap(go_up$Description, width =40)
#指定绘图顺序（转换为因子）：
go_up$Description <- factor(go_up$Description,levels=rev(go_up$Description))
print(levels(go_up$Description))  # 检查因子顺序



#Top20富集数目条形图
mytheme <- theme(axis.title = element_text(size = 13),
                 axis.text = element_text(size = 11),
                 plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11)) #自定义主题


p <- ggplot(data = go_up,
              aes(x = Count, y = Description, fill = -log10(pvalue)))+
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_distiller(palette = 'Blues',direction = 3) +
  labs(x = "Number of Gene",
       y = "pathway",
       title = "GO enrichment") +
  theme_bw() +
  mytheme
p

ggsave(filename = "naive_down_enrich.pdf", width = 10, height = 10)


###naive上调top20画热图
subset_go_up = subset(gene_down, cluster%in% c('6:CD8_C4_LEF1','8:CD8_C6_TCF7', '11:CD8_C9_TIGIT'))

subset_go_up <- subset_go_up %>%
  filter(gene %in% downregulated_genes_intersection)




# 使用tidyr::unite将type和cell_type列合并为type_celltype列
sc_pbmc_int<- subset(sc_pbmc_int2,group%in%c('LI_NAR','PB_NAR')&cell_type%in%c('3:CD8_C1_GZMK', '4:CD8_C2_RGS1','5:CD8_C3_GZMH',
                                                                               '7:CD8_C5_KLRB1','9:CD8_C7_FOSB'))

sc_pbmc_int <- subset(sc_pbmc_int,orig.ident%in%c('HZY_L', 'LXY_L', 'LXY_P', 'LYT_L', 'LYT_P',
                                                  'RSQ_L', 'RSQ_P', 'SME_P', 'WMX_L', 'WMX_P', 'XXY_L', 'XXY_P',
                                                  'ZSQ_L', 'ZSQ_P'))


table(sc_pbmc_int$cell_type,sc_pbmc_int$individual)

sc_pbmc_int$cell_type <- droplevels(sc_pbmc_int$cell_type)


# 将identities设置为type_celltype
Idents(sc_pbmc_int) <- sc_pbmc_int@meta.data$individual
# 计算平均表达量
exp <- AverageExpression(sc_pbmc_int, group.by = 'individual',layer = 'data')$RNA


head(upregulated_genes_intersection)
head(downregulated_genes_intersection)

cc <- c(upregulated_genes_intersection,downregulated_genes_intersection)
# 提取这些基因的表达量
selected_genes_expression <- exp[rownames(exp) %in% cc, ]






library(pheatmap)
library(RColorBrewer)

# 绘制热图
p <-pheatmap(selected_genes_expression, 
         cluster_rows = T, 
         cluster_cols = F, 
         show_rownames = T, 
         show_colnames = TRUE, 
         scale = "row", 
         color = colorRampPalette(c("navy","white","firebrick3"))(100),
         main = "Gene Expression")

p


###热图添加标签函数
add.flag <- function(pheatmap,
                     kept.labels,
                     repel.degree) {
  
  # repel.degree = number within [0, 1], which controls how much 
  #                space to allocate for repelling labels.
  ## repel.degree = 0: spread out labels over existing range of kept labels
  ## repel.degree = 1: spread out labels over the full y-axis
  
  heatmap <- pheatmap$gtable
  
  new.label <- heatmap$grobs[[which(heatmap$layout$name == "row_names")]] 
  
  # keep only labels in kept.labels, replace the rest with ""
  new.label$label <- ifelse(new.label$label %in% kept.labels, 
                            new.label$label, "")
  
  # calculate evenly spaced out y-axis positions
  repelled.y <- function(d, d.select, k = repel.degree){
    # d = vector of distances for labels
    # d.select = vector of T/F for which labels are significant
    
    # recursive function to get current label positions
    # (note the unit is "npc" for all components of each distance)
    strip.npc <- function(dd){
      if(!"unit.arithmetic" %in% class(dd)) {
        return(as.numeric(dd))
      }
      
      d1 <- strip.npc(dd$arg1)
      d2 <- strip.npc(dd$arg2)
      fn <- dd$fname
      return(lazyeval::lazy_eval(paste(d1, fn, d2)))
    }
    
    full.range <- sapply(seq_along(d), function(i) strip.npc(d[i]))
    selected.range <- sapply(seq_along(d[d.select]), function(i) strip.npc(d[d.select][i]))
    
    return(unit(seq(from = max(selected.range) + k*(max(full.range) - max(selected.range)),
                    to = min(selected.range) - k*(min(selected.range) - min(full.range)), 
                    length.out = sum(d.select)), 
                "npc"))
  }
  new.y.positions <- repelled.y(new.label$y,
                                d.select = new.label$label != "")
  new.flag <- segmentsGrob(x0 = new.label$x,
                           x1 = new.label$x + unit(0.15, "npc"),
                           y0 = new.label$y[new.label$label != ""],
                           y1 = new.y.positions)
  
  # shift position for selected labels
  new.label$x <- new.label$x + unit(0.2, "npc")
  new.label$y[new.label$label != ""] <- new.y.positions
  
  # add flag to heatmap
  heatmap <- gtable::gtable_add_grob(x = heatmap,
                                     grobs = new.flag,
                                     t = 4, 
                                     l = 4
  )
  
  # replace label positions in heatmap
  heatmap$grobs[[which(heatmap$layout$name == "row_names")]] <- new.label
  
  # plot result
  grid.newpage()
  grid.draw(heatmap)
  
  # return a copy of the heatmap invisibly
  invisible(heatmap)
}


# 示例特定基因列表
#naive上下调基因
genes <- c( 'S100A8','S100A9','CSF3R','RPS10','RPS13','RPS23','RPL11','RPL26','NDUFC1','NDUFA4',
            'ISCU','UQCRFS1','NDUFB9','COX7A2','ATP5PD','NDUFB11','COX7A2','G0S2','APOO',
            'CD8A','RBP4','APOC3','APOC1','HSPA1B','TF','HLA-DRB1','TGFB1',
            'IL2RB','CCL4','CXCR3','FOSB','CCL4','EGR1','CREM','CCL5','RGS1','CD63','EOMES','PIM1','CD27',
            'CD83','TOX','JAK3','IKZF3','CD2','GPR18','CD74','NR4A3','KLRD1','NKG7','IL21R','CD226')




#effector上下调基因
genes <- c( 'HSPA1B','HSPA1A','EGR1','IFNG','APOA2','CCL3','NR4A3','CRTAM',
            'IL16','IFI16','ATF4','LY9','RORA','CD83','IL27RA','IL6ST','EOMES','RGS1','PIM1','GPR18',
            'LY9','CD27','BCL3','CCL4','ATF2','JUNB','FOS','NR4A1','FASLG',
            'RPL11','RPL26','RPL30','COX7C','UQCRFS1','COX7A2','NDUFB9','NDUFA1','NDUFA4','ATP5PD','ATP5ME',
            'COX7A2','COX6C','COX5B','COX6B1','FOXP1','DNAJC10','S100A8','S100A10','CTSS','LYN','PDCD2')



dev.off()
pdf("effector_UP_and_down_heatmap.pdf", width = 10, height = 10)
add.flag(p,
         kept.labels = genes,
         repel.degree = 0.2)
dev.off()



