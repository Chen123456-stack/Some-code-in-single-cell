


t <- sc_pbmc_int2
sc_pbmc_int2 <- t

sc_pbmc_int2 <- subset(sc_pbmc_int2,cell_type%in%c("4:Mono_C1_VCAN",
                                                   "5:Mono_C2_S100A12", "6:Mono_C3_FCGR3A","7:Mono_C4_ATF3","8:Mono_C5_FCN1",
                                                   "9:Mono_C6_LRMDA","10:Mono_C7_RPS13","11:Mono_C8_RTF2","12:Mac_C1_HLA-DRA",
                                                   "13:Mac_C2_RGS1","14:Mac_C3_C1QA","15:Mac_C4_CXCL12","16:Mac_C5_CD5L","17:Mac_C6_TIMD4"))

Idents(sc_pbmc_int2)
Idents(sc_pbmc_int2) = "cell_type" # 把celltype设置为主要的分类标签

#
cell.markers <- FindAllMarkers(object = sc_pbmc_int2, 
                               only.pos = FALSE, # 是否只保留表达相对上调的基因，设置FALSE则会保留下调的
                               test.use = "wilcox", # 默认使用 wilcox 非参数检验，其它选项可以查看说明
                               slot = "data", # 需要注意的是，默认使用 data，而不是 counts
                               min.pct = 0.25, # 设置表达比例的阈值，没有统一标准，差异基因很多的情况下可以把阈值调高，差异基因有1000个够用
                               logfc.threshold = 0.25 # 设置 log2FC 即差异倍率的阈值，没有统一标准，差异基因很多的情况下可以把阈值调高
)

library( xlsx )
write.xlsx(cell.markers, file = "5.8_CellMarkers.xlsx" )

degs <- cell.markers

head(degs)
names(degs)
dim(degs)

markers=degs
{
  
  all_degs=markers
  df=all_degs
  ##筛选阈值确定：p＜0.05，|log2FC|＞1
  p_val_adj = 0.05
  avg_log2FC = 1
  head(all_degs)
  
  
  #根据阈值添加上下调分组标签：
  df$direction <- case_when(
    df$avg_log2FC > avg_log2FC & df$p_val_adj < p_val_adj ~ "up",
    df$avg_log2FC < -avg_log2FC & df$p_val_adj < p_val_adj ~ "down",
    TRUE ~ 'none'
  )
  head(df)
  print(getwd())
  df=df[df$direction!="none",]
  head(df)
  dim(df)
  df$mygroup=paste(df$group,df$direction,sep = "_")
  head(df)
  table(df$mygroup)
  dim(df)
  
  print(getwd())
  
  dir.create("./degs_enrichments0.3fc")
  setwd("./degs_enrichments0.3fc")
  ##########################----------------------enrichment analysis==================================================
  #https://mp.weixin.qq.com/s/WyT-7yKB9YKkZjjyraZdPg
  
  
  {
    library(clusterProfiler)
    library(org.Hs.eg.db) #人
    library(org.Mm.eg.db) #鼠
    library(ggplot2)
    # degs_for_nlung_vs_tlung$gene=rownames(degs_for_nlung_vs_tlung)
    head(df)
    #df$cluster=df$mygroup
    
    df=df %>%dplyr::group_by(cluster)%>%
      filter(p_val_adj <0.05,
             avg_log2FC >0)  #对所有上调基因进行富集
    
    
    sce.markers=df
    head(sce.markers)
    print(getwd())
    ids <- suppressWarnings(bitr(sce.markers$gene, 'SYMBOL', 'ENTREZID', 'org.Hs.eg.db'))
    head(ids)
    
    head(sce.markers)
    tail(sce.markers)
    dim(sce.markers)
    sce.markers=merge(sce.markers,ids,by.x='gene',by.y='SYMBOL')
    head(sce.markers)
    dim(sce.markers)
    sce.markers$group=sce.markers$cluster
    sce.markers=sce.markers[sce.markers$group!="none",]
    dim(sce.markers)
    head(sce.markers)
    #sce.markers=openxlsx::read.xlsx("~/silicosis/spatial_transcriptomicsharmony_cluster_0.5res_gsea/sce.markers_for_each_clusterfor_enrichment.xlsx")
    #sce.markers$cluster=sce.markers$mygroup
    dim(sce.markers)
    head(sce.markers)
    gcSample=split(sce.markers$ENTREZID, sce.markers$cluster)
    
    library(clusterProfiler)
    gcSample # entrez id , compareCluster 
    names(gcSample)
    print("===========开始go============")
    xx <- compareCluster(gcSample, fun="enrichGO",OrgDb="org.Hs.eg.db" ,   #'org.Hs.eg.db',
                         pvalueCutoff=0.05) #organism="hsa",
    
    xx.BP <- compareCluster(gcSample, fun="enrichGO",OrgDb="org.Hs.eg.db" ,   #'org.Hs.eg.db',
                            pvalueCutoff=0.05,readable=TRUE,
                            ont="BP") #organism="hsa",
    p=clusterProfiler::dotplot(object = xx.BP,showCategory = 20,
                               label_format =60) 
    p=p+ theme(axis.text.x = element_text(angle = 90, 
                                          vjust = 0.5, hjust=0.5))
    p
    ggsave('degs_compareCluster-BP_enrichment--3.pdf',plot = p,width = 13,height = 40,limitsize = F)
    ggsave('degs_compareCluster-BP_enrichment--3.png',plot = p,width = 13,height = 40,limitsize = F)
    
    xx.CC <- compareCluster(gcSample, fun="enrichGO",OrgDb="org.Hs.eg.db" ,   #'org.Hs.eg.db',
                            pvalueCutoff=0.05,readable=TRUE,
                            ont="CC") #organism="hsa",
    p=clusterProfiler::dotplot(object = xx.CC,showCategory = 20,
                               label_format =60) 
    p=p+ theme(axis.text.x = element_text(angle = 90, 
                                          vjust = 0.5, hjust=0.5))
    p
    ggsave('degs_compareCluster-CC_enrichment--3.pdf',plot = p,width = 13,height = 40,limitsize = F)
    ggsave('degs_compareCluster-CC_enrichment--3.png',plot = p,width = 13,height = 40,limitsize = F)
    
    xx.MF <- compareCluster(gcSample, fun="enrichGO",OrgDb="org.Hs.eg.db" ,   #'org.Hs.eg.db',
                            pvalueCutoff=0.05,readable=TRUE,
                            ont="MF") #organism="hsa",
    p=clusterProfiler::dotplot(object = xx.MF,showCategory = 20,
                               label_format =60) 
    p=p+ theme(axis.text.x = element_text(angle = 90, 
                                          vjust = 0.5, hjust=0.5))
    p
    ggsave('degs_compareCluster-MF_enrichment--3.pdf',plot = p,width = 13,height = 40,limitsize = F)
    ggsave('degs_compareCluster-MF_enrichment--3.png',plot = p,width = 13,height = 40,limitsize = F)
    
    
    
    print(getwd())
    
    .libPaths()
    
    
    print("===========开始 kegg============")
    gg<-clusterProfiler::compareCluster(gcSample,fun = "enrichKEGG",  #readable=TRUE,
                                        keyType = 'kegg',  #KEGG 富集
                                        organism='hsa',#"rno",
                                        pvalueCutoff = 0.05 #指定 p 值阈值（可指定 1 以输出全部
    )
    
    p=clusterProfiler::dotplot(object = xx,showCategory = 20,
                               label_format =100) 
    p=p+ theme(axis.text.x = element_text(angle = 90, 
                                          vjust = 0.5, hjust=0.5))
    p
    ggsave('degs_compareCluster-GO_enrichment--3.pdf',plot = p,width = 13,height = 40,limitsize = F)
    ggsave('degs_compareCluster-GO_enrichment--3.png',plot = p,width = 13,height = 40,limitsize = F)
    
    xx
    write.csv(xx,file = "compareCluster-GO_enrichment.csv")
    
    p=clusterProfiler::dotplot(gg,showCategory = 20,
                               label_format = 40) 
    p4=p+ theme(axis.text.x = element_text(angle = 90, 
                                           vjust = 0.5, hjust=0.5))
    p4
    print(paste("保存位置",getwd(),sep = "  :   "))
    ggsave('degs_compareCluster-KEGG_enrichment-2.pdf',plot = p4,width = 13,height = 25,limitsize = F)
    ggsave('degs_compareCluster-KEGG_enrichment-2.png',plot = p4,width = 13,height = 25,limitsize = F)
    
    gg
    openxlsx::write.xlsx(gg,file = "compareCluster-KEGG_enrichment.xlsx")
    
    getwd()
    openxlsx::write.xlsx(sce.markers,file = "sce.markers_for_each_clusterfor_enrichment.xlsx")
    
    
    #   save(xx.BP,xx.CC,xx.MF, gg, file = "~/silicosis/spatial_transcriptomicsharmony_cluster_0.5res_gsea/xx.Rdata")
    
    xx.BP@compareClusterResult$oncology="BP"
    xx.CC@compareClusterResult$oncology="CC"
    xx.MF@compareClusterResult$oncology="MF"
    
    
    xx.all=do.call(rbind,list(xx.BP@compareClusterResult,
                              xx.CC@compareClusterResult,
                              xx.MF@compareClusterResult))
    head(xx.all)    
    openxlsx::write.xlsx(xx.all,file = "enrichments_all.xlsx")    
    print(getwd())
    save(xx.BP,xx.CC,xx.MF,xx.all,sce.markers, gg, file = "./xx.Rdata")
    
    
    if (F) { 
      #大鼠
      print("===========开始go============")
      xx <-clusterProfiler::compareCluster(gcSample, fun="enrichGO",OrgDb="org.Rn.eg.db",
                                           readable=TRUE,
                                           ont = 'ALL',  #GO Ontology，可选 BP、MF、CC，也可以指定 ALL 同时计算 3 者
                                           pvalueCutoff=0.05) #organism="hsa", #'org.Hs.eg.db',
      
      print("===========开始 kegg============")
      gg<-clusterProfiler::compareCluster(gcSample,fun = "enrichKEGG",
                                          keyType = 'kegg',  #KEGG 富集
                                          organism="rno",
                                          pvalueCutoff = 0.05 #指定 p 值阈值（可指定 1 以输出全部
      )
      
      p=dotplot(xx) 
      p2=p+ theme(axis.text.x = element_text(angle = 90, 
                                             vjust = 0.5, hjust=0.5))
      p2
      ggsave('degs_compareCluster-GO_enrichment-2.pdf',plot = p2,width = 6,height = 20,limitsize = F)
      xx
      openxlsx::write.xlsx(xx,file = "compareCluster-GO_enrichment.xlsx")
      # -----  
      p=dotplot(gg) 
      p4=p+ theme(axis.text.x = element_text(angle = 90, 
                                             vjust = 0.5, hjust=0.5))
      p4
      print(paste("保存位置",getwd(),sep = "  :   "))
      ggsave('degs_compareCluster-KEGG_enrichment-2.pdf',plot = p4,width = 6,height = 12,limitsize = F)
      gg
      openxlsx::write.xlsx(gg,file = "compareCluster-KEGG_enrichment.xlsx")   
    }  
  }
  
  
}
