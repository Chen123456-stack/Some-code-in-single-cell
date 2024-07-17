####富集分析及基因集评分

####批量GSEA分析#######
library(Seurat)
library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(ggplot2)
library(RColorBrewer)
library(Scillus)
library(GseaVis)
library(ggpubr)
setwd('../')
setwd('./1.NKT')
dir.create('1.GSEA')
setwd('./1.GSEA')


##########GSEA通路批量富集分析（fgsea富集方法）#######
#####msigdbr准备相关基因集#####
#查看基因集的物种
msigdbr_show_species()
#查看人类基因集
m_df = msigdbr(species = "Homo sapiens")
head(m_df)

#查看基因集类别
print(m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat),n=23)

#提取数据集
gene_sets <- msigdbr(species = "Homo sapiens", category = "H",subcategory = 'GO:BP') 
             # %>% dplyr::filter(gs_name=="REACTOME_TCR_SIGNALING")#特定提取某个数据集，基因集打分用
             # gene_sets <-list(TCR_signaling=unique(gene_sets$gene_symbol)) #基因集打分用

#使用 fgsea富集分析需要准备的基因集数据框
gene_sets <- split(x = gene_sets$gene_symbol, f = gene_sets$gs_name)


#使用 clusterProfiler富集分析需要准备的基因集数据框
gene_sets <- gene_sets %>%
  dplyr::distinct(gs_name, gene_symbol) %>%
  as.data.frame()



###########准备富集数据集，将其准备为列表的形式
#####手动下载基因集#######
###fgsea基因集准备
###去掉基因集通路中的下划线和空格
h.sets <- as.list(read_lines('/home/data/t010560/Migsdb/c2.cp.kegg_medicus.v2023.2.Hs.symbols.gmt'))
for(i in 1:length(h.sets)){
  oneline =str_split(h.sets[[i]],'\t')[[1]]
  names(h.sets)[i]= str_replace_all(str_replace(oneline[1],"HALLMARK_",""),'_','')
  h.sets[[i]]= oneline[3:length(oneline)]
}

###不去基因集通路中的下划线和空格
h.sets <- as.list(read_lines('/home/data/t010560/Migsdb/h.all.v2023.2.Hs.symbols.gmt'))
for(i in 1:length(h.sets)){
  oneline =str_split(h.sets[[i]],'\t')[[1]]
  names(h.sets)[i] = oneline[1]
  h.sets[[i]]= oneline[3:length(oneline)]
}

####clusterprofile基因集准备
h.sets <-read.gmt("/home/data/t010560/Migsdb/h.all.v2023.2.Hs.symbols.gmt")
head(h.sets)


########准备seurat数据
immune <- sc_pbmc_int2
Idents(immune) <- "cell_type"  #设置默认ident，必备
table(immune$cell_type,immune$group)

#批量/单独获取差异基因
DEG <- find_diff_genes(dataset = immune, 
                       clusters = c('7:Mono_C4_ATF3'
                                    # ,'NK_cycle', 'CD4Tn','CD4Tm', 'CD4Treg', 'CD8Tn',
                                    # 'CD8Te', 'CD8Tm', 'CD8pro', 'CD8MAIT', 'DN'
                                    ),
                       comparison = c("group", "LI_NAR", "LI_AR"),
                       logfc.threshold = 0.01,   # threshold of 0 is used for GSEA
                       min.cells.group = 0.01)   # To include clusters with only 1 cell

##细胞亚群间的差异
DEG = FindMarkers( sc_pbmc_int2, 
                          logfc.threshold = 0.01,
                          min.pct = 0.01, # 表达比例的阈值设置
                          only.pos = FALSE,
                          ident.1 = "Kup_C1_CXCL2", ident.2 = "Kup_C2_LRMDA", # ident 1 vs 2
                          group.by = "cell_type" ) %>% # 指定用来分组的信息
  # filter( p_val_adj < 0.05 ) %>% # 需要按校准P筛选可以加上这条
  mutate( gene = rownames(.) ) # 将行名新增为列，方便检索基因和后续作图
####差异基因clusterprofile gsea富集##########
#准备基因列表（只需要基因名和avg_log2FC即可）
DEG <- final_result
genelist <- DEG$avg_log2FC
names(genelist) <- DEG$feature
genelist <- sort(genelist,decreasing = T)#按照avg_log2FC进行排序
class(genelist)
head(genelist)
#gsea富集分析
GSEA_enrichment <- GSEA(genelist,#待富集的基因列表
                        TERM2GENE = h.sets,  #自己选择基因集
                        pvalueCutoff = 1,  #指定 p 值阈值（可指定 1 以输出全部）
                        pAdjustMethod = 'BH',
                        minGSSize=5,
                        maxGSSize = 800,
                        eps = 0,
                        seed = 123,
                        nPermSimple = 10000)
#提取富集结果，做柱状图
df <- GSEA_enrichment@result
# %>% filter(p.adjust<0.1)
write.csv(df,'GSEA_7:Mono_C4_ATF3.csv')


#筛选并按NES排序
df_sig <- df %>% filter(qvalue <0.8) %>%arrange(desc(NES))
#选择部分进行作图
df_sig <- df_sig[c(1,2,3,4,5,6,8,9,10,30,31,32),]
#根据NES定义上调下调
df_sig$label <- ifelse(df_sig$NES > 0, "up", "down")
#作柱状图
ggplot(df_sig, aes(reorder(Description, NES), NES,fill=label)) + 
  geom_bar(stat = 'identity',alpha = 1) + 
  scale_fill_manual(breaks=c("up","down"),values = c("red","#08519C"))+
  labs(x = "Description")+
  coord_flip()+
  theme_bw()
ggsave(filename = '1.pdf',width = 8,height = 6)

#棒棒糖图
ggdotchart(df_sig, x = "Description", y = "NES", 
           color = "qvalue",
           sorting = "descending",
           add = "segments",
           rotate = TRUE,
           dot.size = 8,
           gradient.cols = c("darkgreen", "grey"),
           label=df_sig$size,
           font.label = list(color = "white", size = 8, vjust = 0.5),              
           ggtheme = theme_pubr(),
           ylab = F)+
  geom_hline(yintercept = 0, linetype=2, linewidth=0.5)+
  theme(legend.position = "bottom", legend.title = element_blank(), legend.text = element_text(size = 5))
ggsave(filename = '1.pdf',width = 8,height = 6)
#gsea通路图GseaVis包
#选择展示通路
setid <- c('HALLMARK_INTERFERON_ALPHA_RESPONSE',
           'HALLMARK_INTERFERON_GAMMA_RESPONSE',
           'HALLMARK_IL6_JAK_STAT3_SIGNALING',
           'HALLMARK_MTORC1_SIGNALING',
           'HALLMARK_ALLOGRAFT_REJECTION'
)
#一个/多个通路同时展示
gseaNb(object  = GSEA_enrichment,
       geneSetID = setid,
       curveCol = ggsci::pal_npg()(9),
       subPlot = 2,
       # legend.position='right',
       addPval = T,
       pvalX = 1,
       pvalY = 1.2)

#多个通路分开展示
lapply(setid, function(x){
  gseaNb(object = GSEA_enrichment,
         geneSetID = x,
         addPval = T,
         pvalX = 0.8,pvalY = 0.8,
         pCol = 'black',
         pHjust = 0)
}) -> gseaList
# 合并图片
cowplot::plot_grid(plotlist = gseaList,ncol = 2,align = 'hv')
ggsave(filename = '2.pdf',width = 20,height = 22)



###批量gsea富集（fgsea富集）####
gsea_res <- test_GSEA(DEG,pathway = h.sets)###可选择自己想要富集的基因集

#dotplot作图（多个细胞的多个富集结果）
plot_GSEA(gsea_res, p_cutoff = 0.05, colors = c("#0570b0", "grey", "#d7301f"))

ggsave(filename = 'gsea_liver.pdf',width = 12,height = 10)


####调整后的plot，设置NES为颜色，-log10adjustP为大小

plot_GSEA <- function(gsea_res, 
                      p_cutoff = 0.05,
                      colors = c('#0570b0','grey','#d7301f')) {
  
  gsea_res %>%
    filter(.data$padj <= p_cutoff) %>%
    # mutate(color = -log10(.data$padj) * sign(.data$NES)) %>%
    ggplot() + 
    geom_point(aes(x = factor(.data$pathway), 
                   y = factor(.data$cluster),
                   size = -log10(.data$padj), 
                   color = sign(.data$NES))) +
    scale_size(name = bquote(-log[10]~"Adj. p-value")) +
    scale_color_gradient2(name = 'NES',
                          low = colors[1],
                          mid = colors[2],
                          high = colors[3],
                          midpoint = 0) +
    coord_flip() +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
}

#####保留HALLMARK
plot_GSEA <- function(gsea_res, 
                      p_cutoff = 0.05,
                      colors = c('#0570b0','grey','#d7301f')) {
  
  t <- gsea_res %>%
    filter(.data$padj <= p_cutoff) %>%
    mutate(color = -log10(.data$padj) * sign(.data$NES)) %>%
    ggplot() + 
    geom_point(aes(x = factor(.data$pathway), 
                   y = factor(.data$cluster),
                   size = abs(.data$NES), 
                   color = .data$color)) +
    scale_size(name = "Normalized\nEnrichment\nScore Size") +
    scale_color_gradient2(name = bquote(-log[10]~"Adj. p-value"), 
                          low = colors[1], 
                          mid = colors[2], 
                          high = colors[3], 
                          midpoint = 0) +
    coord_flip() +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
}





####KS科研 clusterProfiler gsea go/kegg富集分析#####


#加载函数

#参数解释---GSEA分析是基于已经完成差异分析结果，且纳入所有基因
#KS_GSEA作用是进行转录组数据GSEA分析，提供clusterProfiler和fgsea两种R包分析方式
#KS_GSEA适用于human和mouse两个物种，支持KEGG、GO（BP）的GSEA分析
#parameter
# @gene：做好的差异分析文件的gene symbol。
# @LogFC：做好的差异分析文件的变化倍数 logFC。
# @analysis：选择参数('GO',"KEGG")
# @package：选择需要使用的包('clusterProfiler','fgsea')
# @OrgDb：物种选择，人："org.Hs.eg.db"， 小鼠："org.Mm.eg.db"

KS_GSEA <- function(gene,
                    LogFC,
                    analysis=c('GO',"KEGG"),
                    package=c('clusterProfiler','fgsea'),
                    OrgDb=c("org.Hs.eg.db", "org.Mm.eg.db")){
  if(OrgDb=="org.Hs.eg.db"){
    organism = "hsa"
    species = "Homo sapiens"
  }
  
  if(OrgDb=="org.Mm.eg.db"){
    organism = "mmu"
    species = "Mus musculus"
  }
  
  if(package=="clusterProfiler"){
    entrezID <- bitr(gene, fromType = "SYMBOL",toType = "ENTREZID", OrgDb = OrgDb)#genesymbol转化为ID
    
    genelist <- LogFC
    names(genelist) <- gene
    genelist <- genelist[names(genelist) %in% entrezID[,1]]
    names(genelist) <- entrezID[match(names(genelist),entrezID[,1]),2]
    genelist <- sort(genelist,decreasing = T)
    
    #install.packages('R.utils')
    R.utils::setOption( "clusterProfiler.download.method",'auto')
    
    if(analysis=="KEGG"){
      KEGG_gesa <- gseKEGG(geneList = genelist,
                           organism = organism,#不同物种查询：https://www.genome.jp/kegg/catalog/org_list.html
                           minGSSize = 10,
                           maxGSSize = 500,
                           pvalueCutoff = 0.05,
                           pAdjustMethod = "BH",
                           verbose = FALSE,
                           eps = 0)
      
      KEGG_gesa_result = setReadable(KEGG_gesa, OrgDb = OrgDb, keyType = "ENTREZID")
      
      KEGG_gesa_table <- KEGG_gesa_result@result
      write.csv(KEGG_gesa_table, file = './KEGG_gesa_table.csv')
      
      
      return(KEGG_gesa_result)
      
    }
    else{
      GO_gesa <- gseGO(geneList = genelist,
                       ont = "BP",
                       OrgDb=OrgDb,
                       minGSSize = 10,
                       maxGSSize = 500,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                       verbose = FALSE)
      
      GO_gesa_result = setReadable(GO_gesa,OrgDb = OrgDb, keyType = "ENTREZID")
      
      GO_gesa_table <- GO_gesa_result@result 
      write.csv(GO_gesa_table, file = './GO_gesa_table.csv')
      
      return(GO_gesa_result)
      
    }
    
    
  }
  
  if(package=="fgsea"){
    if(analysis=="KEGG"){
      
      geneset_KEGG = msigdbr(species = species,#mouse:Mus musculus
                             category = "C2", 
                             subcategory = "CP:KEGG") %>% dplyr::select(gs_name,gene_symbol)
      geneset_KEGG$gs_name <- gsub('KEGG_','',geneset_KEGG$gs_name)#去除前缀KEGG_
      geneset_KEGG$gs_name <- tolower(geneset_KEGG$gs_name)#将大写换为小写
      geneset_KEGG$gs_name <- gsub('_',' ',geneset_KEGG$gs_name)#将_转化为空格
      library(Hmisc)
      geneset_KEGG$gs_name <- capitalize(geneset_KEGG$gs_name)#首字母大写
      GSEA_geneset <- geneset_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)
      
    }else{
      
      geneset_GO = msigdbr(species = species,#mouse:Mus musculus
                           category = "C5", 
                           subcategory = "GO:BP") %>% dplyr::select(gs_name,gene_symbol)
      
      geneset_GO$gs_name <- gsub('GOBP_','',geneset_GO$gs_name)#去除前缀KEGG_
      geneset_GO$gs_name <- tolower(geneset_GO$gs_name)#将大写换为小写
      geneset_GO$gs_name <- gsub('_',' ',geneset_GO$gs_name)#将_转化为空格
      library(Hmisc)
      geneset_GO$gs_name <- capitalize(geneset_GO$gs_name)#首字母大写
      GSEA_geneset <- geneset_GO %>% split(x = .$gene_symbol, f = .$gs_name)
    }
    
    df <- data.frame(logFC = LogFC, gene = gene)
    df <- df[order(df$logFC,decreasing = T),]
    ranks <- df$logFC
    names(ranks) <- df$gene
    
    ## GSEA分析
    GSEA_df <- fgsea(pathways = GSEA_geneset, 
                     stats = ranks,
                     minSize=10,
                     maxSize=500,
                     eps=0.0)
    library(data.table)
    fwrite(GSEA_df, file="./GSEA_df.txt", sep="\t", sep2=c("", " ", ""))
    
    return(GSEA_df)
    
  }
  
}

#参数解释---KS_GSEA_plot函数GSEA分析结果可视化是基于完成GSEA分析得到的文件
#KS_GSEA_plot作用是进行GSEA结果的可视化，可选择需要的通路进行展示

#parameter
# @inputType：输入文件是基于哪种包分析的结果('clusterProfiler','fgsea')
# @analysis：执行的分析是哪种？('GO',"KEGG")
# @data：GSEA分析得到的结果文件，如果是clusterProfiler包，输入的是一个 large gseaResult对象，如果是fgsea包， 输入是一个"data.table" "data.frame"文件
# @term: 需要展示的term，一张图只能展示一个
# @gene：做好的差异分析文件的gene symbol。
# @LogFC：做好的差异分析文件的变化倍数 logFC。
# @OrgDb：物种选择，人："org.Hs.eg.db"， 小鼠："org.Mm.eg.db"

KS_GSEA_plot <- function(inputType=c('clusterProfiler','fgsea'),
                         analysis=c('GO',"KEGG"),
                         data,
                         term,
                         gene,
                         LogFC,
                         OrgDb
){
  
  if(OrgDb=="org.Hs.eg.db"){
    species = "Homo sapiens"
  }
  
  if(OrgDb=="org.Mm.eg.db"){
    species = "Mus musculus"
  }
  
  
  
  if(inputType=='clusterProfiler'){
    #clusterprofile
    gseaScores <- getFromNamespace("gseaScores", "DOSE")
    
    # define function
    gsInfo <- function(object, terms) {
      geneList <- object@geneList
      
      gsea_result_table <- object@result
      site=which(gsea_result_table$Description==terms, arr.ind = TRUE)
      genesetid=gsea_result_table$ID[site]
      
      geneSetID <- object@result[genesetid, "ID"]
      
      geneSet <- object@geneSets[[geneSetID]]
      exponent <- object@params[["exponent"]]
      df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
      df$ymin <- 0
      df$ymax <- 0
      pos <- df$position == 1
      h <- diff(range(df$runningScore))/20
      df$ymin[pos] <- -h
      df$ymax[pos] <- h
      df$geneList <- geneList
      
      df$Description <- object@result[geneSetID, "Description"]
      return(df)
    }
    
    gsea_plotData = gsInfo(data, terms = term)
    
    gsea_plotData <- gsea_plotData %>%
      mutate("gene_name" = data@gene2Symbol) %>%
      filter(position == 1) 
    
    colnames(gsea_plotData)[2] <- 'y'
    
    test = gsea_plotData
    
    #NES和adjust Pvalue
    gesa_table <- data@result
    terms = gesa_table[gesa_table$Description==term,]
    labels_NES = round(terms$NES, 4)
    labels_FDR = format(terms$qvalue, scientific = T,digits = 2)
    
    
    
  }
  
  if(inputType=='fgsea'){
    if(analysis=="KEGG"){
      geneset_KEGG = msigdbr(species = species,#mouse:Mus musculus
                             category = "C2", 
                             subcategory = "CP:KEGG") %>% dplyr::select(gs_name,gene_symbol)
      geneset_KEGG$gs_name <- gsub('KEGG_','',geneset_KEGG$gs_name)#去除前缀KEGG_
      geneset_KEGG$gs_name <- tolower(geneset_KEGG$gs_name)#将大写换为小写
      geneset_KEGG$gs_name <- gsub('_',' ',geneset_KEGG$gs_name)#将_转化为空格
      library(Hmisc)
      geneset_KEGG$gs_name <- capitalize(geneset_KEGG$gs_name)#首字母大写
      GSEA_geneset <- geneset_KEGG %>% split(x = .$gene_symbol, f = .$gs_name)
      
    }else{
      geneset_GO = msigdbr(species = species,#mouse:Mus musculus
                           category = "C5", 
                           subcategory = "GO:BP") %>% dplyr::select(gs_name,gene_symbol)
      
      geneset_GO$gs_name <- gsub('GOBP_','',geneset_GO$gs_name)#去除前缀KEGG_
      geneset_GO$gs_name <- tolower(geneset_GO$gs_name)#将大写换为小写
      geneset_GO$gs_name <- gsub('_',' ',geneset_GO$gs_name)#将_转化为空格
      library(Hmisc)
      geneset_GO$gs_name <- capitalize(geneset_GO$gs_name)#首字母大写
      GSEA_geneset <- geneset_GO %>% split(x = .$gene_symbol, f = .$gs_name)
      
      
    }
    
    
    df <- data.frame(logFC = LogFC, gene = gene)
    df <- df[order(df$logFC,decreasing = T),]
    ranks <- df$logFC
    names(ranks) <- df$gene
    
    
    
    ToPlot  <- sapply(data$pathway, function(Pathway){
      pathway <- GSEA_geneset[[Pathway]]
      stats <- ranks
      rnk <- rank(-stats)
      ord <- order(rnk)
      statsAdj <- stats[ord]
      statsAdj <- sign(statsAdj)*(abs(statsAdj)^1)
      statsAdj <- statsAdj/max(abs(statsAdj))
      pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
      pathway <- sort(pathway)
      gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, returnAllExtremes = T)
      bottoms <- gseaRes$bottoms
      tops <- gseaRes$tops
      n <- length(statsAdj)
      xs <- as.vector(rbind(pathway - 1, pathway))
      ys <- as.vector(rbind(bottoms, tops))
      toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0), pathway = Pathway)
      return(list(toPlot))
    })
    
    #点图数据
    gsea_plotData <- ToPlot[[term]]
    
    test <-  gsea_plotData[1:floor(nrow(gsea_plotData)/length(1)),]
    #NES和adjust Pvalue
    terms = data[data$pathway==term,]
    labels_NES = round(terms$NES, 4)
    labels_FDR = format(terms$padj, scientific = T,digits = 2)
    
  }
  
  
  
  test$xend = test$x+1
  test$xend2 = ''
  test$xend3 = ''
  
  for (i in 1:nrow(test)) {
    test$xend2[i] = test$xend[i+1] - test$xend[i]
  }
  
  for (i in 1:ncol(test)){
    test[,i][is.na(test[,i])] <- 5
  } 
  
  test$xend2 <- as.numeric(test$xend2)
  
  for (i in 1:nrow(test)) {
    test$xend3[i] = (test$x[i] + test$xend2[i])
  }
  
  test$xend3 <- as.numeric(test$xend3)
  
  
  if(analysis=="KEGG"){
    title_terms = paste0("KEGG:",term)
  }
  
  if(analysis=="GO"){
    title_terms = paste0("GO:",term)
  }
  
  
  if(labels_NES>0){
    fill_color = "#90191B"
  }
  
  if(labels_NES<0){
    fill_color = "#3F90C9"
  }
  
  
  p=ggplot(test) + 
    geom_rect(aes(xmin = x-1,xmax = xend3, ymin = -0.04 , ymax = 0.04, fill=x), lwd=4)+
    scale_fill_gradientn(colours = colorRampPalette(c("#90191B","white","#3F90C9"))(100))+
    geom_rect(aes(xmin = x,xmax = xend, ymin = -0.04, ymax = 0.04, fill=x), color="black", size=0.5)+
    geom_point(data=gsea_plotData, aes(x = x, y = y),fill=fill_color, shape=21, size=4) + 
    geom_hline(yintercept = 0, linetype=3, lwd = 1) +
    scale_x_continuous(expand = c(0.01,0.01))+
    ylab("Enrichment Score") +
    xlab('')+
    labs(title=title_terms)+
    theme_bw() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(colour = "black",size = 1),
          axis.text.x=element_blank(),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.title.y = element_text(size = 15,colour = 'black'),
          legend.position="none",
          plot.title = element_text(size=15,hjust =0.5, face = 'bold'),
          axis.ticks.x = element_blank())
  
  if(labels_NES>0){
    lable_y = max(gsea_plotData$y)-0.2
    
    a=nrow(gsea_plotData)
    
    lable_x = gsea_plotData$x[a]-gsea_plotData$x[ceiling(a/1.8)]
    
    ylim = expand_limits(y=c(-0.2, NA))
  }
  
  if(labels_NES<0){
    
    lable_y = min(gsea_plotData$y)+0.2
    a=nrow(gsea_plotData)
    
    lable_x = gsea_plotData$x[a]-gsea_plotData$x[ceiling(a/2)]
    ylim = expand_limits(y=c(NA, 0.2))
    
  }
  
  p+ylim+annotate(geom = 'text', label=paste("NES =",labels_NES, "\n", "FDR =", labels_FDR), 
                  x=lable_x, y=lable_y, size=4)
  
  
}


#差异基因-----------------------------------------------------------------------
load("D:/KS项目/公众号文章/GSEA分析及可视化函数/sce_mar.RData")
df <- FindMarkers(Macrophage, min.pct = 0, logfc.threshold = 0,
                  group.by = "group",ident.1 ="BM",ident.2="GM")

df <- NK_cycle_degs
df$gene = rownames(df)


#GSEA分析-----------------------------------------------------------------------
GSEA_CP <- KS_GSEA(gene = df$gene,
                   LogFC = df$avg_log2FC,
                   analysis = "KEGG",
                   package = 'clusterProfiler',
                   OrgDb = 'org.Hs.eg.db')
class(GSEA_CP)
# [1] "gseaResult"
# attr(,"package")
# [1] "DOSE"



#GSEA可视化-----------------------------------------------------------------------
#NES>0
p1=KS_GSEA_plot(inputType = "clusterProfiler",
                analysis = "KEGG",
                data = GSEA_CP,
                term = 'Cytokine-cytokine receptor interaction',
                gene = df$gene,
                LogFC  = df$avg_log2FC,
                OrgDb = 'org.Hs.eg.db')

#NES<0
p2=KS_GSEA_plot(inputType = "clusterProfiler",
                analysis = "KEGG",
                data = GSEA_CP,
                term = 'Ribosome',
                gene = df$gene,
                LogFC  = df$avg_log2FC,
                OrgDb = 'org.Hs.eg.db')

#批量循环
#NES>0
pathway <- c('Oxidative phosphorylation',
             "Parkinson disease",
             "Biosynthesis of amino acids",
             "Cardiac muscle contraction")

pathway_list <- list()
for (i in 1:length(pathway)) {
  p = KS_GSEA_plot(inputType = "clusterProfiler",
                   analysis = "KEGG",
                   data = GSEA_CP,
                   term = pathway[i],
                   gene = df$gene,
                   LogFC  = df$avg_log2FC,
                   OrgDb = 'org.Hs.eg.db')
  pathway_list[[i]] <- p
  
}
#组合图
CombinePlots(pathway_list, ncol = 2)


#棒棒糖图，这个图很好看
library(ggpubr)
df <- read.csv('KEGG_gesa_table.csv', header = T)
df_sig <- df %>%
  filter(qvalue <0.01) %>%
  arrange(desc(qvalue))

ggdotchart(df_sig[0:15,], x = "Description", y = "NES", 
           color = "qvalue",
           sorting = "descending",
           add = "segments",
           rotate = TRUE,
           dot.size = 8,
           gradient.cols = c("darkgreen", "grey"),
           font.label = list(color = "white", size = 8, vjust = 0.5),              
           ggtheme = theme_pubr(),
           ylab = F)+
  geom_hline(yintercept = 0, linetype=2, linewidth=0.5)


#####基因集打分，以M1,M2巨噬细胞打分为例#######
#手动下载基因集
###导入M1
m1m2_pws <- read_lines("./GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP.v2023.2.Hs.gmt") %>%
  lapply(str_split, "\\t") %>% 
  unlist(recursive = F) %>% 
  lapply(function(x) setNames(list(x[-c(1:2)]), x[1])) %>% 
  unlist(recursive = F)

### 导入M2
m1m2_pws <- append(m1m2_pws, read_lines("./GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN.v2023.2.Hs.gmt") %>%
                     lapply(str_split, "\\t") %>% 
                     unlist(recursive = F) %>% 
                     lapply(function(x) setNames(list(x[-c(1:2)]), x[1])) %>% 
                     unlist(recursive = F))

### 运用addmodulescore函数，加入富集分数到metadata中
imm_anno <- AddModuleScore(object = sc_pbmc_int2, features = m1m2_pws, name = c("M1up", "M1dn"), nbin = 12)

#提琴图
VlnPlot(imm_anno, features = c("M1up1", "M1dn2"), 
        group.by = "cell_type", pt.size = 0)

#箱线图
data<- FetchData(imm_anno,vars = c("cell_type","M1_UP1"))#取得基因集打分后得值
p<-ggplot(data,aes(cell_type,M1_UP1))+
   geom_boxplot()+theme_bw()+RotatedAxis()

data<- FetchData(imm_anno,vars = c("cell_type","M2_UP2"))
p<-ggplot(data,aes(cell_type,M2_UP2))+
   geom_boxplot()+theme_bw()+RotatedAxis()

#提取单细胞数据中的某个基因比如CCL3/CCL3L3绘制箱线图以及进行差异分析（相对于小提琴图更加直观）
#不分组差异分析
#设置要对比的对象
my_comparisons <- list( c("CXCL2_MACRO", "DSCAM_MACRO"), 
                        c("CD5L_MACRO", "CXCL2_MACRO"),c("CCL5_MACRO", "CXCL2_MACRO"))
data<- FetchData(imm_anno,vars = c("cell_type","CCL3",'group'),layer='data')#取得单个基因
p<-ggplot(data,aes(cell_type,CCL3,color=cell_type))+
  geom_boxplot()+theme_bw()+RotatedAxis()+
  stat_compare_means(comparisons = my_comparisons)
p
ggsave(filename = 'CCL3.pdf',width =8,height = 6 )

#分组差异分析
data<- FetchData(imm_anno,vars = c("cell_type","CCL3L3",'group'))
p<-ggplot(data,aes(cell_type,CCL3L3,color=group))+
  geom_boxplot(aes(fill=group),
               alpha=0.1)+
  # geom_jitter(position = position_jitterdodge(jitter.height=0.75, # 散点抖动高度
  #                                             jitter.width = 0.5, # 散点抖动宽度
  #                                             dodge.width = 0.75))+ # x轴方向上的闪避量
  scale_color_manual(values = pal_npg('nrc')(9)[c(1,3)])+
  scale_fill_manual(values = pal_npg('nrc')(9)[c(1,3)])+
  theme_bw()+
  theme(panel.grid = element_blank())+
  stat_compare_means(aes(group = group),label="p.signif",
                     show.legend = T)
p
ggsave(filename = 'CCL3_group.pdf',width =8,height = 6 )

#featureplot作图
FeaturePlot(imm_anno,'M1_UP1',cols=rev(brewer.pal(10, name = "RdBu")))
FeaturePlot(imm_anno,'M2_UP2',cols=rev(brewer.pal(10, name = "RdBu")))

#stat_compare_means函数添加P值
my_comparisons <- list( c("CXCL2_MACRO", "DSCAM_MACRO"), c("CD5L_MACRO", "DSCAM_MACRO"),c("CCL5_MACRO", "DSCAM_MACRO"))

Idents(imm_anno) <- 'cell_type'

VlnPlot(imm_anno, features =c("M2_UP2") , pt.size = 0) +
  NoLegend() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_compare_means(comparisons = my_comparisons)+
  ylim(-0.15, 0.15)

####其他基因集打分####

imm_flam <- sc_pbmc_int2
table(sc_pbmc_int2$cell_type,sc_pbmc_int2$group)
sc_pbmc_int2 <- subset(sc_pbmc_int2,type%in%c('LIVER')&cell_type%in%c('4:Mono_C1_VCAN','5:Mono_C2_S100A12','6:Mono_C3_FCGR3A','7:Mono_C4_ATF3',
                                                                      '8:Mono_C5_FCN1','9:Mono_C6_LRMDA','10:Mono_C7_RPS13','11:Mono_C8_RTF2',
                                                                      '12:Mac_C1_HLA-DRA','13:Mac_C2_RGS1','14:Mac_C3_C1QA','15:Mac_C4_CXCL12',
                                                                      '16:Mac_C5_CD5L','17:Mac_C6_TIMD4'))

#对细胞类型进行排序
sc_pbmc_int2$cell_type <- factor(sc_pbmc_int2$cell_type,
                          levels = c('cDC1', 'cDC2',  'Mono_FCGR3A' ,'Mono_CD14',  'Mono_Inter', 'Mono_macro','Macro'))
#HALLMARK_INFLAMMATORY_RESPONSE
sc_pbmc_int2 <- t

  inflam <- read_lines("~/Migsdb/phagocytosis.v2023.2.Hs.gmt") %>%
  lapply(str_split, "\\t") %>% 
  unlist(recursive = F) %>% 
  lapply(function(x) setNames(list(x[-c(1:2)]), x[1])) %>% 
  unlist(recursive = F)
inflam[[2]] <- NULL
imm_flam <- AddModuleScore(object = sc_pbmc_int2, features = inflam,
                           name = c("ANTIGEN_PROCESSING_AND_PRESENTATION"),  
                           seed = 1, search = FALSE, assay = NULL,
                           nbin = 21, ctrl = 5, pool = NULL)
#普通提琴图
# c("#089E86","#3D5387")

#提琴图加箱线图
VlnPlot(imm_flam, features = c("ANTIGEN_PROCESSING_AND_PRESENTATION1"), split.by  = 'group',
        group.by = "cell_type", pt.size = 0, cols = c("#089E86","#3D5387"),fill=NULL)+
  # geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.size = 0)+
  theme(legend.position = "none") +
  labs(x = NULL,y = "MHC-II mediated antigen processing presentation",title = NULL) + 
  theme(axis.text.x = element_text(angle=45,size=12))+
  theme(axis.text.y = element_text(size=12),
          axis.title.y = element_text(size = 12))
ggsave(filename = 'MHC-II mediated antigen processing presentation.pdf',width =10,height = 6) 



data<- FetchData(imm_flam,vars = c("cell_type","ANTIGEN_PROCESSING_AND_PRESENTATION1",
                                   'orig.ident','group'))
##ggplot 小提琴

#不同cell_type小提琴图
ggplot(data, aes(x = cell_type, y = ANTIGEN_PROCESSING_AND_PRESENTATION1)) + 
  geom_violin(fill = NA,color = "#4A856D",size=0.8, scale = "width") +  # 小提琴图不填充颜色，保留边框颜色
  geom_boxplot(width = 0.15, fill = NA, color = "#4A856D", outlier.size = 0) + 
  theme_bw() + 
  theme(legend.position = "none", 
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),   # 去掉次要网格线
        axis.text.x = element_text(angle = 30, vjust = 0.5)) +
  labs(x = NULL, y = "Kupffer score", title = NULL) +
  scale_fill_manual(values = mycolor)

# 绘制分组的小提琴图和箱线图
ggplot(data, aes(x = cell_type, y = ANTIGEN_PROCESSING_AND_PRESENTATION1, fill = group)) + 
  geom_violin(aes(color = group), fill = NA, size = 0.8, scale = "width", position = position_dodge(0.8)) +  # 小提琴图不填充颜色，保留边框颜色
  # geom_boxplot(width = 0.15, aes(color = group), fill = NA, outlier.size = 0, position = position_dodge(0.8)) + 
  theme_bw() + 
  theme(legend.position = "top", 
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),   # 去掉次要网格线
        axis.text.x = element_text(angle = 30, vjust = 0.5)) +
  labs(x = NULL, y = "Kupffer score", title = NULL) +
  scale_fill_manual(values = mycolor) +
  scale_color_manual(values = mycolor)


# 
# VlnPlot(imm_flam, features = c('ANTIGEN_PROCESSING_AND_PRESENTATION1'),split.by = 'group', 
#         group.by = "cell_type", pt.size = 0)+
#   labs(x = NULL,y = "MHC-II_ANTIGEN_PROCESSING",title = NULL) 
ggsave(filename = 'Kupffer.pdf',width = 12,height = 4) 
ggsave(filename = 'Class I MHC Mediated Antigen Processing Presentation.pdf',width = 8,height = 6) 

#显著性提琴图
VlnPlot(imm_flam, features =c("INFLAMMATORY_RESPONSE1") , pt.size = 0) +
  NoLegend() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_compare_means(comparisons = my_comparisons)+
  ylim(-0.3, 0.5)
ggsave(filename = 'T Cell Proliferation.pdf',width = 10,height = 8) 


#箱线图
data<- FetchData(imm_flam,vars = c("cell_type","ANTIGEN_PROCESSING_AND_PRESENTATION1",
                                   'orig.ident','group'))



data$cell_type <- factor(data$cell_type,levels = rev(c("1:DC_C1_CD1C","2:DC_C2_CLEC9A","3:DC_C3_MKI67","4:Mono_C1_VCAN",
                                                                   "5:Mono_C2_S100A12", "6:Mono_C3_FCGR3A","7:Mono_C4_ATF3","8:Mono_C5_FCN1",
                                                                   "9:Mono_C6_LRMDA","10:Mono_C7_RPS13","11:Mono_C8_RTF2","12:Mac_C1_HLA-DRA",
                                                                   "13:Mac_C2_RGS1","14:Mac_C3_C1QA","15:Mac_C4_CXCL12","16:Mac_C5_CD5L","17:Mac_C6_TIMD4")))


#fill=group,

#填充箱线图
ggplot(data,aes(cell_type,ANTIGEN_PROCESSING_AND_PRESENTATION1,pt.size=0)) + 
  geom_boxplot(outlier.shape =1,color =c("#3D5387"),width=0.6,outlier.size = 0.8,size=0.8) + 
  theme_bw() + 
  # geom_jitter(position = position_jitter(width = 0.1), color = "black",size=0) +
  labs(x = "", y = "Phagocytosis score") +
  theme(legend.position = "top",
        # panel.grid.major = element_blank(),  # 去掉主网格线
        # panel.grid.minor = element_blank()
        ) + # 去掉主网格线
  theme(axis.text.x = element_text(angle=0,vjust = 0.5))+
  scale_fill_manual(values = mycolor)+
  coord_flip()

ggsave(filename = 'Kupffer.pdf',width =400/3.4,
       height = 400/1.7, units = "mm",dpi = 600) 
ggsave(filename = 'Phagocytosis.pdf',width =400/3.4,
       height = 400/1.7, units = "mm",dpi = 600) 

#非填充箱线图
ggplot(data, aes(cell_type, ANTIGEN_PROCESSING_AND_PRESENTATION1, fill = group, color = group)) + 
  geom_boxplot(outlier.shape = 1, width = 0.6, outlier.size = 0.8, size = 0.8, fill = NA) + 
  theme_bw() + 
  labs(x = "", y = "Kupffer score") +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 315, vjust = 0.5)) +
  scale_fill_manual(values = c("#089E86","#3D5387")) +
  scale_color_manual(values = c("#089E86","#3D5387"))

ggsave(filename = 'Kupffer.pdf',width = 4,height = 10,dpi = 600) 




ggsave(filename = 'Class I MHC Mediated Antigen Processing Presentation group.pdf',width = 10,height = 6)

#求打分平均值
result <- data %>%
  group_by(orig.ident, cell_type,group) %>%
  summarize(ANTIGEN_PROCESSING_AND_PRESENTATION1 = mean(ANTIGEN_PROCESSING_AND_PRESENTATION1))


library(ggpubr)
ggplot(result,aes(cell_type,ANTIGEN_PROCESSING_AND_PRESENTATION1)) + 
  geom_boxplot(outlier.shape = 20,width=0.5,color = "#4A856D") + 
  theme_bw() +
  geom_jitter(position = position_jitter(width = 0.001), color = "#4A856D",size=1) +
  labs(x = "", y = "Kupffer score") +
  theme(legend.position = "top",
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),      # 去掉面板边框
        axis.line = element_line(color = "black")) + 
  theme(axis.text.x = element_text(angle=315,vjust = 0.5,size=12))+
  scale_fill_manual(values = mycolor)




#HALLMARK_ALLOGRAFT_REJECTION
allo <- read_lines("../../.././HALLMARK_ALLOGRAFT_REJECTION.v2023.2.Hs.gmt") %>%
  lapply(str_split, "\\t") %>% 
  unlist(recursive = F) %>% 
  lapply(function(x) setNames(list(x[-c(1:2)]), x[1])) %>% 
  unlist(recursive = F)
imm_allo <- AddModuleScore(object = sc_pbmc_int2, features = allo,
                           name = c("ALLOGRAFT_REJECTION"),  
                           seed = 1, search = FALSE, assay = NULL,
                           nbin = 21, ctrl = 5, pool = NULL)

#普通提琴图
VlnPlot(imm_allo, features = c("ALLOGRAFT_REJECTION1"), 
        group.by = "cell_type", pt.size = 0)

#显著性提琴图
VlnPlot(imm_allo, features =c("ALLOGRAFT_REJECTION1") , pt.size = 0) +
  NoLegend() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_compare_means(comparisons = my_comparisons)+
  ylim(-0.1, 0.5)
ggsave(filename = 'ALLOGRAFT_REJECTION.pdf',width = 8,height = 6) 



###批量获取差异基因(自己写的)#####
library(dplyr)
results_list <- list()

sc_pbmc_int2 <- sc
table(sc_pbmc_int2$cell_type,sc_pbmc_int2$group)

for (i in c('NK_rest', 'NK_cycle'
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
  results_list[[i]] <- group_degs
}
# 将列表中的数据框合并为一个数据框
final_result <- bind_rows(results_list)
#过滤性别相关基因
gender_related_genes <- c("DDX3Y", "TXLNGY", "RPS4Y1",'USP9Y','EIF1AY','KDM5D',
                          'UTY','XIST','ZFY','LINC00278','ENSG00000247970',
                          'PRKY','KDM5D',	'NLGN4Y','TTTY14','TBL1Y','TMSB4Y',
                          'ZNF839P1',	'TTTY15',	'AC009977.1',	'AC011297.1','HLA-DRB6')  # 用实际的性别相关基因名称替换这些示例基因

final_result <- final_result %>%
  filter(!gene %in% gender_related_genes)

#导出为excel
write.xlsx(final_result, file = "NK_liver.xlsx" )

#取相应子集
NK_cycle_degs <- final_result %>%
  filter(cluster == "NK_cycle")

NK_rest_degs <- final_result %>%
  filter(cluster == "NK_rest")

























