setwd('~/NSCLC/Test/')
library(clusterProfiler)
library(org.Hs.eg.db)
library(readr)

Fibro_marker <-readRDS('./res/data/OV_SC_Fibro_marker.rds')


############################################################
#
#--------------------富集分析————————————————————————————
#
############################################################

Overlap_genes_ENTREZID = bitr(Fibro_marker$gene, #数据集
                              fromType="SYMBOL", #输入为SYMBOL格式
                              toType="ENTREZID",  # 转为ENTERZID格式
                              OrgDb="org.Hs.eg.db") #人类 数据库

ego_ALL <- enrichGO(gene = Overlap_genes_ENTREZID$ENTREZID,
                    OrgDb = org.Hs.eg.db, #没有organism="human"，改为OrgDb=org.Hs.eg.db
                    #keytype = 'ENSEMBL',
                    ont = "ALL", #也可以是 CC  BP  MF中的一种
                    pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                    pvalueCutoff = 1, #P值会过滤掉很多，可以全部输出
                    qvalueCutoff = 1,
                    readable = TRUE) 

pdf('./res/fig/02_Fig2B_GOEnrichment_Plot.pdf',height = 9,width = 6)
dotplot(ego_ALL, showCategory = 10, split="ONTOLOGY") + 
  facet_grid(ONTOLOGY~., scale="free")
dev.off()

kk <-  enrichKEGG(
  gene          = Overlap_genes_ENTREZID$ENTREZID,
  keyType     = "kegg",
  organism   = 'hsa',
  pvalueCutoff      = 1,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 1
)

str(kk)
head(kk)
kk<-as.data.frame(kk)
gene_symbols <- bitr(kk$geneID, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)

pdf('./res/fig/02_Fig2B_KEGGEnrichment_Plot.pdf',height = 6,width = 6)
dotplot(kk, showCategory = 15)
dev.off()

write_csv(as.data.frame(ego_ALL),'./res/data/02_ego_ALL.csv')
write_csv(as.data.frame(kk),'./res/data/02_kk.csv')

