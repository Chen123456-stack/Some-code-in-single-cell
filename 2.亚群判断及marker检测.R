
#亚群判定及检测
########淋巴系粗检测
gene = list(
  CD4=c("CD4", "CD3D","CD3E","IL7R","CD3G",'CCNB1','CCR7'),
  CD8=c("CD8A")  ,
  NK=c('NKG7',"GNLY", "NCAM1","FCGR3A",'KLRD1','GZMB','KLRC1','KLRF1'),
  B = c('MS4A1','CD19','CD79A', 'CD79B', 'CD38'))

#####淋巴细胞亚群细分#####
gene= list(
  B=c('MS4A1','CD19','CD79A'),
  Tcell=c('CD3D','CD3E','CD3G','CD4','CD8A'),
  Native=c('CCR7','LEF1','SELL','TCF7'),
  Cytotoxic_effector=c('GNLY','IFNG','NKG7','PRF1',"GZMH","GZMB",'GZMA','GZMK'),
  Immune_checkpoint=c('HAVCR2','LAG3','PDCD1','CTLA4','TIGIT','BTLA','KLRC1'),
  Tem=c('ANXA1','ANKRD28','IL7R','CD69','CD40LG','ZNF683'),
  CD4Treg=c('FOXP3','IL2RA','IKZF2'),
  NK=c('NCR1','NCAM1','TYROBP','FGFBP2','KLRD1','KLRF1','KLRB1','CX3CR1','FCGR3A',
       'XCL1','XCL2'),
  Tγδ=c('TRGC1','TRDC'),
  Proliferating=c('MKI67','TOP2A'))


#########髓系细胞marker#####
gene = list(
  Mac=c("C1QA","C1QB","C1QC",'CD68','CD163','VCAM1','CD5L','MARCO',"SELENOP","RNASE1","DAB2","LGMN","PLTP","MAF","SLCO2B1"),
  mono=c("VCAN","FCN1","CD300E","S100A12","EREG","APOBEC3A","STXBP2","ASGR1","CCR2","NRG1"),
  CD14 = c('CD14'),#Mono_CD14
  CD16 = c('FCGR3A'),#Mono_FCGR3A
  neutrophils =  c("FCGR3B",'CSF3R',"CXCR2","SLC25A37","G0S2","CXCR1","ADGRG3","PROK2","STEAP4","CMTM2" ),
  pDC = c("GZMB","SCT","CLIC3","LRRC26","LILRA4","PACSIN1","CLEC4C","MAP1A","PTCRA","C12orf75"),
  DC1 = c("CLEC9A","XCR1","CLNK","CADM1","ENPP1","SNX22","NCALD","DBN1","HLA-DOB","PPY"),
  DC2=c( "CD1C","CLEC10A","FCER1A","CD1E","HLA-DQA2","HLA-DQB1","HLA-DPB1","CD2","GPAT3","CCND2","ENHO","PKIB","CD1B"),
  DC3 =  c("HMSD","ANKRD33B","LAD1","CCR7","LAMP3","CCL19","CCL22","INSM1","TNNT2","TUBB2B")
)

gene = list(pDC = c("CLEC4C","IRF7","TCF4","GZMB"),
            cDC1 = c("XCR1","CLNK","CLEC9A"),
            cDC2 = c("FCER1A","HLA-DPB1","HLA-DQB1","CD1E","CD1C","CLEC10A","HLA-DQA2"),
            DC3 = c("CCL19","LAMP3","IDO1","IDO2","LAD1","FSCN1","CCR7","LY75","CCL22","CD40","BIRC3","NFKB2"),
            Macrophages = c("APOC1","HLA-DRB5","C1QA","C1QB"),
            RTMs = c("THBS1"),#Resident tissue macrophages
            Lam = c("APOE"),#Lipid associated macrophages
            Monocytes = c("LYZ","HLA-DRB1","TIMP1","S100A11","CXCL8","IL1B","PTGS2","S100A9","S100A8","MMP19"),
            Mono_C = c('CD14'),#Mono_CD14
            Mono_F = c('FCGR3A'),#Mono_FCGR3A
            Mast = c('TPSAB1' , 'TPSB2'))


genes_to_check =gene

genes_to_check = lapply(genes_to_check, str_to_upper)
p_all_markers=DotPlot(sc_pbmc_int2, 
                      features = genes_to_check,
                      group.by = c('seurat_clusters'),
                      scale = T,assay='RNA' )+
  theme(axis.text.x=element_text(angle=45,hjust = 1))
p_all_markers

ggsave(filename = 'lymphocyte_check1.pdf',width = 30,height = 15)
ggsave(filename = 'lymphocyte_check2.pdf',width = 30,height = 15)
ggsave(filename = 'myeloids_check.pdf',width = 30,height = 15)
ggsave(filename = 'myeloids_check2.pdf',width = 30,height = 15)


#####张泽民T细胞分群######

#########CD8
gene= list(
  CD8TN = c('CCR7','LEF1', 'SELL', 'TCF7', 'CD27', 'S1PR1'),
  CD8TCM=c( 'IL7R', 'CD28', 'PRF1', 'GZMA', 'CCL5', 'GPR183'),
  CD8TEMRATEFF = c('KLRG1', 'CX3CR1', 'FCGR3A', 'FGFBP2', 'GZMH', 'TBX21', 'EOMES','S1PR5'),
  CD8TEM= c('GZMK', 'CXCR4', 'CXCR3', 'CD44'),
  CD8TRM=c('CD6', 'XCL1', 'XCL2', 'MYADM', 'CAPG', 'NR4A1','NR4A2','NR4A3', 'CD69', 'ITGAE'),
  CD8IEL = c('CD160', 'KIR2DL4','TMIGD2', 'KLRC1','KLRC2','KLRC3', 'IKZF2', 'ENTPD1'),
  CD8TEX =c('HAVCR2', 'CXCL13', 'PDCD1', 'LAYN', 'TOX', 'IFNG', 'GZMB', 'MIR155HG', 'TNFRSF9'),
  MAIT = c('SLC4A10', 'KLRB1', 'ZBTB16', 'NCR3', 'RORC', 'RORA'))


####CD4
gene= list(
  CD4+TN <- c('CCR7','LEF1', 'SELL', 'TCF7', 'CD27','CD28', 'S1PR1'),
  CD4+Blood-TCM <- c('CCR7', 'SELL', 'PTGER2', 'ICAM2', 'ANXA1', 'ANXA2', 'S1PR1'),
  CD4+TEMRA/TEFF <- c('KLRG1','CX3CR1', 'NKG7', 'PRF1', 'GNLY', 'GZMH', 'TBX21', 'CTSW', 'S1PR1', 'S1PR5'),
  CD4+Normal-TCM <- c('CCR7','TCF7','RGS1', 'CD69'),
  CD4+TRM <- c('CD69', 'KLRB1', 'PTGER4', 'IL7R', 'CXCR6', 'NR4A1','NR4A2','NR4A3', 'MYADM'),
  TFH <- c('CXCR5', 'BCL6', 'ICA1', 'TOX', 'TOX2', 'IL6ST', 'MAGEH1', 'BTLA', 'ICOS', 'PDCD1', 'CD200'),
  CD4+TEM/TH1-ike<- c('GZMK', 'GZMA', 'CCL5', 'IFNG', 'RUNX3', 'EOMES', 'CXCR3', 'CXCR4', 'CD44'),
  TH17R <- c('IL23R', 'RORC', 'IL17A', 'FURIN', 'CTSH', 'CCR6', 'KLRB1.', 'CAPG', 'ITGAE'),
  TH1-likeTcells <- c('CXCL13', 'IFNG', 'CXCR3', 'BHLHE40', 'GZMB', 'PDCD1', 'HAVCR2', 'ICOS', 'IGFLR1', 'ITGAE'),
  Blood-Treg<- c('FOXP3', 'IL2RA', 'IL10RA', 'IKZF2', 'RTKN2', 'CDC25B', 'S1PR4'),
  TFR <- c('FOXP3','IL2RA','CXCR5','PDCD1','IL10','CCR4','CD69'),
  Tumour-Treg<- c('FOXP3', 'CCR8', 'TNFRSF18', 'LAYN', 'TNFRSF9', 'IKZF2', 'RTKN2', 'CTLA4', 'BATF', 'IL21R'))



