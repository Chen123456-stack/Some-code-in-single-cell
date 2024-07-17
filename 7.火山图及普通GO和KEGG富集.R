#火山图及普通富集分析

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse )
library(clusterProfiler )
library(enrichplot)
library(cols4all)
library(xlsx)
#准备基因列表
#基因列表准备均按logfc.threshold = 0.01,min.cells.group = 0.01准备

dir.create('./CD4_C1_LEF1')
setwd('./CD4')

##批量获取差异基因
###寻找AR/NAR各个亚群差异基因
results_list_li <- list()#创建空列表
table(sc_pbmc_int2$cell_type,sc_pbmc_int2$group)
#通过for循环实现
for (i in c('3:CD8_C1_GZMK'
  # 'NK_rest', 'NK_cycle'
            # , 'CD4Tn','CD4Tm', 'CD4Treg', 'CD8Tn', 'CD8Te', 'CD8Tm', 'CD8pro', 'CD8MAIT', 'DN'
)) {
  # 寻找差异基因
  group_degs <- FindMarkers(
    subset(sc_pbmc_int2, cell_type == i), 
    logfc.threshold = 0.25,
    min.pct = 0.1,
    ident.1 = "LI_AR",
    ident.2 = "LI_NAR",
    group.by = "group"
  ) %>% 
    mutate(gene = rownames(.), cluster = i)
  # 将结果添加到列表中
  results_list_li[[i]] <- group_degs
}
# 将列表中的数据框合并为一个数据框
final_result <- bind_rows(results_list_li)
#性别相关基因
# gender_related_genes <- c("DDX3Y", "TXLNGY", "RPS4Y1",'USP9Y','EIF1AY','KDM5D',
#                           'UTY','XIST','ZFY','LINC00278','ENSG00000247970',
#                           'PRKY','KDM5D',	'NLGN4Y','TTTY14','TBL1Y','TMSB4Y',
#                           'ZNF839P1',	'TTTY15',	'AC009977.1',	'AC011297.1','HLA-DRB6')  # 用实际的性别相关基因名称替换这些示例基因

# # 过滤性别相关基因
# final_result <- final_result %>% filter(!gene %in% gender_related_genes)
#导出基因列表
# write.xlsx(final_result, file = "NK_liver.xlsx" )
# 
# #分别获取单个亚群的基因列表
# NK_cycle_degs <- final_result %>%
#   filter(cluster == "NK_cycle")
# 
# NK_rest_degs <- final_result %>%
#   filter(cluster == "NK_rest")


##标准
NK_degs <- group_degs
NK_degs$difference <- NK_degs$pct.1 - NK_degs$pct.2#获取Δpct

NK_degs <- NK_degs[NK_degs$gene != 'ENSG00000276241', ]


log2FC = 1
padj = 0.05 
NK_degs$threshold="ns";
NK_degs[which(NK_degs$avg_log2FC  > log2FC & NK_degs$p_val_adj <padj),]$threshold="up";
NK_degs[which(NK_degs$avg_log2FC  < (-log2FC) & NK_degs$p_val_adj < padj),]$threshold="down";
NK_degs$threshold=factor(NK_degs$threshold, levels=c('down','ns','up'))

ggplot(NK_degs, aes(x=difference, y=avg_log2FC,color = threshold)) + 
  geom_point(size=2) + 
  scale_color_manual(values=c( "#E64B35" , "gray", "#00468B"),
                     limits = c("up", "ns", "down")) + 
  geom_label_repel(data=subset(NK_degs, avg_log2FC >= 1 & difference >= 0.1 & p_val_adj <= 0.05), 
                   aes(label=gene),  #添加label
                   color="black", #设置label中标签的颜色
                   segment.colour = "black",#设置label框的颜色
                   label.padding = 0.1, 
                   max.overlaps = 10,
                   segment.size = 0.3,  #框的大小
                   size=4)+
  geom_label_repel(data=subset(NK_degs, avg_log2FC <= -1 & difference <= -0.1 & p_val_adj <= 0.05), 
                   aes(label=gene), label.padding = 0.1, 
                   color="black",
                   segment.colour = "black",
                   segment.size = 0.3, size=4,
                   max.overlaps = 10)+
  geom_label_repel( data = subset(NK_degs, gene %in% c('CCR9','IL18','CD27','CCR5','MKI67','IL1B','ITGA1',
                                                       'CCR5','RGS1','LEF1','CCL3','CCL3L3','CCL4','CCR1',
                                                      'IL2RA','GZMB','CD38')),
                   aes( label = gene ),size = 4,label.padding = 0.1,  color = "black" ,max.overlaps = 5)+
   # geom_label_repel( data = subset(NK_degs, gene %in% c('CCR9','IL18','CD27','CCR5','MKI67','IRF7','IL1B','ITGA1',
   #                                                      'CCR5','RGS1','LEF1','TCF7',
   #                                                     'CXCR4','IL2RA','GZMA','CXCR6','GZMK','GZMB','LDLR','CD38')),
   #                  aes( label = gene ),size = 4,label.padding = 0.1,  color = "black" ,max.overlaps = 5)+
  geom_vline(xintercept = 0.0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  theme_bw()+
  theme(legend.title =element_blank(),
        legend.justification.right = 3,
    legend.position = 'right')

ggsave(filename = '13:NK_C2_XCL1_12:NK_C1_FCGR3A.pdf',width = 10,height = 9)


#暂时舍弃，也可以画

# 
# # 制作差异火山图
# #对角线火山图，更适合单细胞数据，只考虑
# NK_degs <- final_result
# NK_degs$difference <- NK_degs$pct.1 - NK_degs$pct.2#获取Δpct
# #选择展示标签的基因
# NK_degs_sig <- NK_degs[which(NK_degs$p_val_adj<0.05 & abs(NK_degs$avg_log2FC) >1),]
# NK_degs_sig$gene <- rownames(NK_degs_sig)
# 
# #作图
# ggplot(NK_degs, aes(x=difference, y=avg_log2FC)) + 
#   geom_point(size=2, color="grey60") + 
#   geom_point(data=NK_degs[which(NK_degs$p_val_adj<0.05 & NK_degs$avg_log2FC>0.1),],
#              aes(x=difference, y=avg_log2FC),
#              size=2, color="red")+
#   geom_point(data=NK_degs[which(NK_degs$p_val_adj<0.05 & NK_degs$avg_log2FC< -0.1),],
#              aes(x=difference, y=avg_log2FC),
#              size=2, color="blue")+
#   #仅展示手动选择的部分基因
#   # geom_label_repel( data = subset(NK_degs, gene %in% c('CCR1','CCR2','CCR9','IL18','CD27','CCR5','MKI67','IRF7','IL1B',
#   #                                                      'PIM1','CXCR3','IL2RA','S100A9','GZMK','GZMB','CCL3','LDLR','CD38')),
#   #                   aes( label = gene ),size = 4.5, color = "black" )+
#   #展示筛选出的基因集
#   geom_label_repel( data = NK_degs_sig,aes( label = gene ),size = 4.5, color = "black" )+
#   theme_classic()+
#   theme(axis.text.x = element_text(colour = 'black',size = 12),
#         axis.text.y = element_text(colour = 'black',size = 12),
#         axis.title = element_text(colour = 'black',size = 15),
#         axis.line = element_line(color = 'black', size = 1))+
#   geom_hline(yintercept = 0,lty=2,lwd = 1)+
#   geom_vline(xintercept = 0,lty=2,lwd = 1)+
#   ylab("Log-fold Change")+
#   xlab("Delta Percent")
# ggsave(filename = 'vol.pdf',width = 6,height = 6)



##传统火山图
cell.markers <- NK_rest_degs
colnames( cell.markers )
table( abs(cell.markers$avg_log2FC) > 2 )
#去除一些高表达没有意义的基因，影响火山图美观
rows_to_keep <- !(rownames(cell.markers) %in% c("IFI44L", "MX1", "ENSG00000262202", "XAF1",
                                                'ENSG00000289826','ENSG00000278847'))
cell.markers <- cell.markers[rows_to_keep, ]
#筛选出部分基因用于集体展示标签
plotdt = cell.markers %>% mutate( gene = ifelse( abs( avg_log2FC ) >= 0 & p_val_adj < 0.05, gene, NA ) )
dim(plotdt)
#作图
  ggplot( cell.markers, aes( x = avg_log2FC, y = -log10( p_val_adj ), 
                                size = pct.1,
                                color = avg_log2FC ) ) +
  geom_point( ) +
  ggtitle( label = "NK_rest", subtitle = "AR vs NAR" ) +
  #展示筛选出的批量基因标签
  geom_text_repel(plotdt, aes( label = gene ), size = 4.5, color = "black" ) +
  #展示自己选择标签并进行位置调节
  geom_label( data = subset(cell.markers, gene %in% c('CCR2','CCR9','IL18','CD27','IRF7','IL1B',
                                                      'PIM1','CXCR3','GZMK','GZMB','CCL3','LDLR','CD38')),
              aes( label = gene ), size = 3, color = "black" )+
  geom_label( data = subset(cell.markers, gene %in% c('S100A9')),
              aes( label = gene ), size = 3, color = "black",vjust=-1.5,hjust=0.1)+
  geom_label( data = subset(cell.markers, gene %in% c('CCR1')),
              aes( label = gene ), size = 3, color = "black",vjust=1.5,hjust=0.1)+
  geom_label( data = subset(cell.markers, gene %in% c('IL2RA','MKI67')),
              aes( label = gene ), size = 3, color = "black",vjust=-0.5,hjust=0)+
  geom_label( data = subset(cell.markers, gene %in% c('CCR5')),
              aes( label = gene ), size = 3, color = "black",vjust=0,hjust=1.5)+
  theme_bw( )+
  theme(
    plot.title = element_text( face = "bold", hjust = 0.5),
    plot.background = element_rect(fill = "transparent",colour = NA)
  ) +
  scale_color_gradient2( low = "blue", high = "red", 
                         mid = "grey", midpoint = 0 ) +
  scale_size( range = c(1,3) )
ggsave(filename = 'NKdim cluster huoshan.pdf',width =8, height = 6 )

###GO/KEGG富集分析#########
#根据差异基因avg_log2FC正负添加上下调情况
cell.markers <- group_degs
cell.markers$change <- ifelse(cell.markers$avg_log2FC > 0, "UP", "DOWN")

# 筛选差异基因
lsxx = cell.markers %>% 
  #filter( cluster == "CD4+ T cells" ) %>% #筛选某一亚群差异基因
  filter( p_val_adj < 0.05 ) %>% #筛选p
  filter( avg_log2FC >  1 | avg_log2FC < -1 )%>%   #0.585
  arrange( desc( avg_log2FC ) ) #对avg_log2FC进行降序排序

#分别筛选出上调和下调的gene
gene_up = subset(lsxx, change == 'UP')
gene_down = subset(lsxx, change == 'DOWN')


### 1.3 富集分析----
# 基于 clusterProfiler 包

##BP
#go富集分析
#上调基因富集
go_up <- enrichGO(OrgDb="org.Hs.eg.db", # 种属基因数据
                       gene = gene_up$gene, 
                       keyType='SYMBOL',
                       ont ="BP", # CC BP MF
                       pvalueCutoff = 1, # 设置校准P值筛选阈值
                       qvalueCutoff = 1, # 设置Q值筛选阈值
                       readable = TRUE # 自动将 ID 转换为 SYMBOL
) %>% 
  # clusterProfiler::simplify( . ) %>% # 删减重复条目，可选，相似结果比较多时好用
  enrichplot::pairwise_termsim() # emap 网络图需要做此转换

#下调基因富集
go_down <- enrichGO(OrgDb="org.Hs.eg.db", # 种属基因数据
                         gene = gene_down$gene, 
                         keyType='SYMBOL',
                         ont ="BP", # CC BP MF
                         pvalueCutoff = 1, 
                         qvalueCutoff = 1, 
                         readable = TRUE # 自动将 ID 转换为 SYMBOL
) %>% 
  # clusterProfiler::simplify( . )%>% # 删减重复条目，可选，相似结果比较多时好用
  enrichplot::pairwise_termsim() # emap 网络图需要做此转换

##CC
go_up <- enrichGO(OrgDb="org.Hs.eg.db", # 种属基因数据
                  gene = gene_up$gene, 
                  keyType='SYMBOL',
                  ont ="CC", # CC BP MF
                  pvalueCutoff = 1, # 设置校准P值筛选阈值
                  qvalueCutoff = 1, # 设置Q值筛选阈值
                  readable = TRUE # 自动将 ID 转换为 SYMBOL
) %>% 
  # clusterProfiler::simplify( . ) %>% # 删减重复条目，可选，相似结果比较多时好用
  enrichplot::pairwise_termsim() # emap 网络图需要做此转换

#下调基因富集
go_down <- enrichGO(OrgDb="org.Hs.eg.db", # 种属基因数据
                    gene = gene_down$gene, 
                    keyType='SYMBOL',
                    ont ="CC", # CC BP MF
                    pvalueCutoff = 1, 
                    qvalueCutoff = 1, 
                    readable = TRUE # 自动将 ID 转换为 SYMBOL
) %>% 
  # clusterProfiler::simplify( . )%>% # 删减重复条目，可选，相似结果比较多时好用
  enrichplot::pairwise_termsim() # emap 网络图需要做此转换


##MF
go_up <- enrichGO(OrgDb="org.Hs.eg.db", # 种属基因数据
                  gene = gene_up$gene, 
                  keyType='SYMBOL',
                  ont ="MF", # CC BP MF
                  pvalueCutoff = 1, # 设置校准P值筛选阈值
                  qvalueCutoff = 1, # 设置Q值筛选阈值
                  readable = TRUE # 自动将 ID 转换为 SYMBOL
) %>% 
  clusterProfiler::simplify( . ) %>% # 删减重复条目，可选，相似结果比较多时好用
  enrichplot::pairwise_termsim() # emap 网络图需要做此转换

#下调基因富集
go_down <- enrichGO(OrgDb="org.Hs.eg.db", # 种属基因数据
                    gene = gene_down$gene, 
                    keyType='SYMBOL',
                    ont ="MF", # CC BP MF
                    pvalueCutoff = 1, 
                    qvalueCutoff = 1, 
                    readable = TRUE # 自动将 ID 转换为 SYMBOL
) %>% 
  # clusterProfiler::simplify( . )%>% # 删减重复条目，可选，相似结果比较多时好用
  enrichplot::pairwise_termsim() # emap 网络图需要做此转换


dotplot(go_up)
dotplot(go_down)

#筛选p.adjust<0.05通路
go_up <- go_up@result %>%
  filter(p.adjust<0.05) %>%
  mutate(change='up')

go_down <- go_down@result %>%
  filter(p.adjust<0.5) %>%
  mutate(change='down')

go_up_sig <- go_up[c(1:10),]
go_down_sig <- go_down[c(1:10),]
#合并上调与下调
df <- rbind(go_up_sig, go_down_sig)

 # #将下调对应的pvalue指定负数，以便绘图：
df$'-log10pvalue' <- ifelse(df$change == 'up', -log10(df$pvalue), -(-log10(df$pvalue)))



# head(df)
# 
# #指定因子；调整绘图顺序按照上下调分别排序：
# level_up <- df$Description[df$change=='up']
# level_down <- df$Description[df$change=='down']
# level <- c(level_up,level_down)
# unique_levels <- unique(level) 
# 
# df$Description <- factor(df$Description,
#                          levels = rev( unique(level)))
# 
# ###选取top5可视化
# top10_df <- df %>% group_by(change) %>%arrange('-log10pvalue') %>%slice_head(n = 10)
# head(df)

###较为好看的富集展示图

df$Description <- str_wrap(df$Description, width =20)
p <- ggplot(df, aes(reorder(Description, `-log10pvalue`), `-log10pvalue`,fill=change))+ 
  geom_bar(stat = 'identity',alpha = 0.7) + 
  scale_fill_manual(breaks=c("up","no","down"),values = c("red","gray","#08519C"))+
  labs(x = "Pathways")+
  coord_flip()+
  theme_bw()+
  theme(
    axis.text.y = element_blank()
  )+
  geom_text(data = df[1:10,], aes(x = Description, y = -0.1, label = Description),size= 4, hjust = 1, color = 'black')+
  geom_text(data = df[11:20,], aes(x = Description, y = 0.1, label = Description), size= 4,hjust = 0, color = 'black')

p

ggsave(filename = '11:CD8_C9_TIGIT_Go.pdf',width = 10,height = 10)

#重新绘制富集条形图（代码同上）：
ggplot(df, aes(reorder(Description, `-log10pvalue`), `-log10pvalue`,fill=change)) + 
  geom_bar(stat = 'identity',alpha = 1) + 
  scale_fill_manual(breaks=c("up","down"),values = c("red","#08519C"))+
  labs(x = "Description")+
  coord_flip()+
  theme_bw()+
  theme(axis.text.y = element_text(face = "bold", size = 12))+ 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40) )
ggsave(filename = 'Macro_GO.pdf',width = 10,height = 6)

p1 <- ggplot(top10_df,aes(x = `-log10pvalue`,y = Description,fill = 'change')) +
  geom_col() +
  theme_classic()
p1
#添加不同颜色
#部分色板，感兴趣可自行替换尝试效果；
#army_rose，sunset，temps，blue_red2，div_bu_wh_rd，div_rainbow，div_tq_wh_pk
p <- p1 + scale_fill_continuous_c4a_div('sunset', mid = 0)+
  scale_y_discrete(labels = function(x) str_wrap(x, width = 50))+ 
  theme(axis.text.x = element_text(size = 10))+
  theme(axis.text.y = element_text(size = 12))#大于0上调；小于0下调
p
ggsave(p,filename = 'NKdim cluster go.pdf',width = 10, height = 6 )




# symbol 转 ENTREZID
annot = bitr( lsxx$gene, fromType = "SYMBOL",  toType = c("ENTREZID"),
              OrgDb = "org.Hs.eg.db", drop = F ) %>% 
  drop_na() %>% 
  distinct( SYMBOL, .keep_all = T ) 
#将ENTREZID转到基因列表中
lsxx = lsxx %>% 
  mutate( ENTREZID = annot$ENTREZID[ match( gene, annot$SYMBOL ) ] ) %>% 
  na.omit( ENTREZID ) 

#分别筛选出上调和下调的gene
gene_up = subset(lsxx, change == 'UP')
gene_down = subset(lsxx, change == 'DOWN')


#kegg富集分析
kegg_up <- enrichKEGG(gene = gene_up$ENTREZID, 
                           keyType = "kegg",
                           organism   = 'hsa', # 种属设置
                           pvalueCutoff = 1, # 设置1意味着不筛选，保留所有结果，可以在后续处理中筛选
                           qvalueCutoff = 1,
                           use_internal_data = F ) %>% # 是否使用本地数据，默认 F 
  setReadable( ., OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')  %>% # 将 ENTRZID 转换成 SYMBOL
  pairwise_termsim() # emap 网络图需要做此转换

kegg_down <- enrichKEGG(gene = gene_down$ENTREZID, 
                             keyType = "kegg",
                             organism   = 'hsa', # 种属设置
                             pvalueCutoff = 1, # 设置1意味着不筛选，保留所有结果，可以在后续处理中筛选
                             qvalueCutoff = 1,
                             use_internal_data = F ) %>% # 是否使用本地数据，默认 F 
  setReadable( ., OrgDb = 'org.Hs.eg.db', keyType = 'ENTREZID')  %>% # 将 ENTRZID 转换成 SYMBOL
  pairwise_termsim() # emap 网络图需要做此转换

dotplot(kegg_up)
dotplot(kegg_down)

kegg_up <- kegg_up@result %>%
  filter(p.adjust<0.1) %>%
  mutate(change='up')
kegg_down <- kegg_down@result %>%
  filter(p.adjust<0.1) %>%
  mutate(change='down')
#合并上调与下调
df <- rbind(up_kegg,down_kegg)

#将下调对应的pvalue指定负数，以便绘图：
df$'-log10pvalue' <- ifelse(df$change == 'up', -log10(df$pvalue), -(-log10(df$pvalue)))
head(df)

#指定因子；调整绘图顺序按照上下调分别排序：
level_up <- df$Description[df$change=='up']
level_down <- df$Description[df$change=='down']
level <- c(level_up,level_down)
unique_levels <- unique(level) 
df$Description <- factor(df$Description,
                         levels = rev( unique(level)))
###选取top5可视化
top10_df <- df %>%group_by(change) %>%arrange('-log10pvalue') %>% slice_head(n = 10)
head(df)

#重新绘制富集条形图（代码同上）：
p1 <- ggplot(top10_df,aes(x = `-log10pvalue`,y = Description,fill = `-log10pvalue`)) +
  geom_col() +
  theme_classic()
p1

p1 + scale_fill_continuous_c4a_div('sunset', mid = 0) #大于0上调；小于0下调


#
save( go_data, kegg_data, file = '1.3-常规富集分析.Rdata' )




## === 其他数据库通路
###其他数据库的分析，如HALLMARK，可以先在MsigDB中将该基因集下载下来
geneset <- read.gmt("data/MsigDB/v7.4/h.all.v7.4.symbols.gmt")
##geneset一共两列，一列是通路，一列是基因名；所有通路前都有HALLMARK，后续可以去掉，比较直观
geneset$term <- gsub(pattern = "HALLMARK","", geneset$term)
geneset$term <- str_to_title(geneset$term)  #大小写转换

my_path <- enricher(gene=DEG, pvalueCutoff = 1, qvalueCutoff = 1, TERM2GENE=geneset)
#使用enricher这个函数，gene是我们的差异表达基因，TERM2GENE是指定手动输入的集合
p1 <- barplot(my_path, showCategory=15,color = "pvalue")
#根据p值选择top15
p1
ggsave("result/enrich_HALLMARK.png", plot = p1, width = 12, height = 10)

my_path <- data.frame(my_path)
write.csv(my_path,"result/enrich_HALLMARK.csv")






















































