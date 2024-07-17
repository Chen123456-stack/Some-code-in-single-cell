######不同亚群之间富集分析差异 同时画GO和KEGG

library(viridis)



mycol3 <- c('#6BA5CE', '#F5AA5F')
cmap <- c("viridis", "magma", "inferno", "plasma", "cividis", "rocket", "mako", "turbo")


table(CD8_Obj$cell_type)

##只保留上调即可

Idents(CD8_Obj) <- CD8_Obj$cell_type
cell.markers <- FindAllMarkers(object = CD8_Obj, 
                               only.pos = T, # 是否只保留表达相对上调的基因，设置FALSE则会保留下调的
                               test.use = "wilcox", # 默认使用 wilcox 非参数检验，其它选项可以查看说明
                               slot = "data", # 需要注意的是，默认使用 data，而不是 counts
                               min.pct = 0.25, # 设置表达比例的阈值，没有统一标准，差异基因很多的情况下可以把阈值调高，差异基因有1000个够用
                               logfc.threshold = 0.25 # 设置 log2FC 即差异倍率的阈值，没有统一标准，差异基因很多的情况下可以把阈值调高
)


# 筛选差异基因
lsxx = cell.markers %>% 
  filter( cluster == "14:NK_C3_KLRF1" ) %>% #筛选某一亚群差异基因
  filter( p_val_adj < 0.05 ) %>% #筛选p
  filter( avg_log2FC >  0.585 | avg_log2FC < -0.585 )%>%   #0.585
  arrange( desc( avg_log2FC ) ) #对avg_log2FC进行降序排序


### 1.3 富集分析----
# 基于 clusterProfiler 包

##BP
#go富集分析
#上调基因富集
go_up <- enrichGO(OrgDb="org.Hs.eg.db", # 种属基因数据
                  gene = lsxx$gene, 
                  keyType='SYMBOL',
                  ont ="BP", # CC BP MF
                  pvalueCutoff = 1, # 设置校准P值筛选阈值
                  qvalueCutoff = 1, # 设置Q值筛选阈值
                  readable = TRUE # 自动将 ID 转换为 SYMBOL
) 
# %>% 
  # clusterProfiler::simplify( . ) %>% # 删减重复条目，可选，相似结果比较多时好用
  # enrichplot::pairwise_termsim() # emap 网络图需要做此转换



allgg <- as.data.frame(go_up)

go_up <- allgg %>%
  filter(p.adjust<0.05)

go_up_sig <- go_up[c(1,3,6),c(2,5,8,9)]




annot = bitr( lsxx$gene, fromType = "SYMBOL",  toType = c("ENTREZID"),
              OrgDb = "org.Hs.eg.db", drop = F ) %>% 
  drop_na() %>% 
  distinct( SYMBOL, .keep_all = T ) 
#将ENTREZID转到基因列表中
lsxx = lsxx %>% 
  mutate( ENTREZID = annot$ENTREZID[ match( gene, annot$SYMBOL ) ] ) %>% 
  na.omit( ENTREZID ) 




#kegg富集分析
kegg_up <- enrichKEGG(gene = lsxx$ENTREZID, 
                      keyType = "kegg",
                      organism   = 'hsa', # 种属设置
                      pvalueCutoff = 1, # 设置1意味着不筛选，保留所有结果，可以在后续处理中筛选
                      qvalueCutoff = 1,
                      use_internal_data = F ) %>% # 是否使用本地数据，默认 F 
  setReadable( ., OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')  %>% # 将 ENTRZID 转换成 SYMBOL
  pairwise_termsim() # emap 网络图需要做此转换

kegg_up <- as.data.frame(kegg_up)

kegg_up <- kegg_up %>%
  filter(p.adjust<0.1)

kegg_up <- kegg_up[c(1,2,3),c(4,7,10,11)]

kegg_up$term <- 'KEGG'
go_up_sig$term <- 'GO'



g <- combine(kegg_up,go_up_sig)

g$term <- factor(g$term,levels = c('KEGG','GO'))

g$Description <- factor(g$Description)


# 先对数据进行排序，确保按顺序排列
g <- g[order(g$term, g$Count),]

# 将 Description 重新因子化，按排序后的顺序
g$Description <- factor(g$Description, levels = unique(g$Description))

g$geneID <- sapply(g$geneID, function(x) {
  genes <- str_split(x, "/")[[1]]
  paste(genes[1:min(5, length(genes))], collapse = "/")
})


# 画图
p <- ggplot(data = g, aes(x = Count, y = Description, fill=term)) +
  geom_bar(width = 0.5, stat = 'identity') +
  theme_bw() + 
  scale_x_continuous(expand = c(0, 0.5)) +
  scale_fill_manual(values = alpha(mycol3, 0.66))

p <- p + theme(axis.text.y = element_blank()) + 
  geom_text(data = g,
            aes(x = 0.1, y = Description, label = Description),
            size = 4.8,
            hjust = 0) 

p <- p+  geom_text(data = g,
                aes(x = 0.1, y = Description, label = geneID , color=-log10(pvalue)),
                size = 4,
                fontface = 'italic',
                hjust = 0,
                vjust = 2.7) +
  scale_colour_viridis(option=cmap[5], direction=-1) 

 
p


ggsave(filename = '13.pdf',width = 6,height = 4,dpi = 300)

# Basic plot
p <- ggplot(data = g, aes(x = Count, y = Description, fill=term)) +
  geom_bar(width = 0.5, stat = 'identity') +
  theme_bw() + 
  scale_x_continuous(expand = c(0, 0.5)) +
  scale_fill_manual(values = alpha(mycol3, 0.66))

# Adding text labels
p <- p + 
  theme(axis.text.y = element_blank()) + 
  geom_text(data = g,
            aes(x = 0.1, y = Description, label = Description, color = -log10(pvalue)),
            size = 4,
            fontface = 'italic',
            hjust = 0,
            vjust = 2.7) +
  scale_colour_viridis_c(option = "E", direction = -1)  # Use appropriate viridis option

# Print the plot
print(p)

ggsave(filename = '12.pdf',width = 6,height = 4,dpi = 300)
