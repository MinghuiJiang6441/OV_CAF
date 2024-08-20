setwd('~/NSCLC/Test/')
library(Seurat)
library(SingleR)
library(celldex)
library(ggplot2)
library(reshape2)
suppressMessages(library(tidyverse))
suppressMessages(library(viridis))
suppressMessages(library(ggsci))
library(reshape2)


source('~/NSCLC/Test/SciBet.R')
source('~/NSCLC/Test/RcppExports.R')
OV_seurat <-readRDS('./data/OV_seurat.rds')

saveRDS(OV_seurat,'./data/OV_seurat.rds')

############################################################
#
#--------------细胞注释————————————————————————————————
#
############################################################

ref <- celldex::HumanPrimaryCellAtlasData()

OV_seurat <- JoinLayers(OV_seurat)

DefaultAssay(OV_seurat) <- "RNA"
seurat_matrix <- GetAssayData(OV_seurat, slot = "data")

# 运行SingleR进行注释
singleR_results <- SingleR(test = seurat_matrix, ref = ref, 
                           labels = ref$label.main)

# 将SingleR结果添加到Seurat对象中
OV_seurat$SingleR.labels <- singleR_results$labels

rm(seurat_matrix)

# 查看注释结果
head(OV_seurat$SingleR.labels)
DimPlot(OV_seurat, reduction = "umap", group.by = "SingleR.labels")

OV_seurat <- AddMetaData(OV_seurat, metadata = singleR_results$labels, 
                            col.name = "SingleR_Labels")

############################################################
#
#--------------------marker基因可视化————————————————————————————
#
############################################################


# 筛选每个细胞簇的标记基因
Idents(OV_seurat) <-"SingleR_Labels" 
markers <- FindAllMarkers(OV_seurat, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 1, 
                          test.use = "wilcox")
# 过滤标记基因
markers <- markers[markers$p_val_adj < 0.05, ]

head(markers)
# 选取每组前5个基因
top_markers <- markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)

pdf('./res/fig/01_Fig1F_DotPlot.pdf',width = 18,height = 6)
DotPlot(OV_seurat, features = unique(top_markers$gene) )+ RotatedAxis()+ 
  scale_color_gradient(low = "white", high = "red")
dev.off()

# 过滤出Fibroblasts群集的基因
fibroblasts_markers <- markers[markers$cluster == "Fibroblasts", ]
# 按照avg_log2FC排序，并选择前10个基因
top10_fibroblasts_genes <- fibroblasts_markers[order(-fibroblasts_markers$avg_log2FC), ][1:10, ]
# 查看前10个基因
top10_fibroblasts_genes

pdf('./res/fig/01_Fig1G_FeaturePlot_fibroblasts_genes.pdf',width = 11,height = 11)
FeaturePlot(OV_seurat, features = top10_fibroblasts_genes$gene[1:9], cols = c("lightgrey", "purple"), reduction = "umap")
dev.off()

# 筛选满足条件的基因
filtered_genes <- fibroblasts_markers[
  abs(fibroblasts_markers$avg_log2FC) > 1 & 
    fibroblasts_markers$pct.1 > 0.25 & 
    fibroblasts_markers$p_val_adj < 0.05, 
]

# 查看筛选后的基因

dim(filtered_genes)
#[1] 588   7

OV_SC_Fibro_marker <-filtered_genes 
saveRDS(OV_SC_Fibro_marker,'./res/data/OV_SC_Fibro_marker.rds')

############################################################
#
#--------------Fibro亚型鉴定——————————————————————————————
#
############################################################

Idents(OV_seurat) <- 'SingleR.labels'
fibroblasts_subset <- subset(OV_seurat, idents = "Fibroblasts")

# 创建样本信息的映射表
sample_info <- data.frame(
  sample = c("GSM5599220", "GSM5599221", "GSM5599222", "GSM5599223", "GSM5599224", 
             "GSM5599225", "GSM5599226", "GSM5599227", "GSM5599228", "GSM5599229", 
             "GSM5599230", "GSM5599231"),
  type = c("normal", "normal", "normal", "normal", "normal", 
           "tumor", "tumor", "tumor", "tumor", "tumor", 
           "tumor", "tumor")
)

# 将样本信息添加到fibroblasts_subset的元数据中
fibroblasts_subset$sample_type <- sample_info$type[match(fibroblasts_subset$sample, sample_info$sample)]

# 查看添加后的元数据
table(fibroblasts_subset@meta.data$sample_type)

fibroblasts_subset <- SCTransform(fibroblasts_subset, verbose = TRUE,assay = 'RNA',vars.to.regress = c('nFeature_RNA','nCount_RNA')) %>% RunPCA()

pc.num = 1:25
fibroblasts_subset<- RunTSNE(fibroblasts_subset, reduction="pca", dims=pc.num,seed.use = 123) %>% 
  RunUMAP(reduction="pca", dims=pc.num,seed.use = 123) %>%
  FindNeighbors(reduction="pca", dims=pc.num) 

  fibroblasts_subset<-  FindClusters( fibroblasts_subset, resolution=1,random.seed = 12345)

marker<- c("DCN","COL1A1",
           "SLPI","PI16",'CLU', #NMF
           "IGF1",'C7','APOD',#NAF
           'CXCL1','CXCL6','CCL2','CCL11','CXCL14',#CAF
           'CXCL5','CXCL8','CSF3','MMP3','MMP1',#CAF
           'TAGLN','ACTA2',"ACTG2", #VSMC
           "ACTA2","MYH11","MYL9",#Myofibroblasts
           "POSTN","CXCL14", #eCAF
           "CFD","CD34",#iCAF
           "RGS5","MCAM","COL4A1","COL4A2", #Pericyte
           'PECAM1', 'VWF'
           
)
DotPlot(fibroblasts_subset,group.by = 'seurat_clusters',features = unique(marker), cluster.idents = T,col.min = 0) + coord_flip()


DimPlot(object = fibroblasts_subset, reduction = "tsne", group.by = c('seurat_clusters'),label = T ,seed=1,label.size = 3)

Pericyte <-1
Myofibroblasts <- c(0,6)
NAF <- 2
CAF <- c(3,4,5)

current.cluster.ids <- c(Pericyte,Myofibroblasts,NAF,CAF)

new.cluster.ids <- c(rep("Pericyte",length(Pericyte)),
                     rep("Myofibroblasts",length(Myofibroblasts)),
                     rep("NAF",length(NAF)),
                     rep("CAF",length(CAF)))

fibroblasts_subset@meta.data$Subcluster <- plyr::mapvalues(x = fibroblasts_subset@active.ident, from = current.cluster.ids, to = new.cluster.ids)
head(fibroblasts_subset@meta.data)

pdf('./res/fig/01Fig1G_CAF_Dimplot.pdf',width = 5,height = 5)
DimPlot(object = fibroblasts_subset, reduction = "tsne",cols = rev(colors), group.by = c('Subcluster'),label = F ,seed=1,label.size = 3)
dev.off()

saveRDS(fibroblasts_subset,'./res/data/fibroblasts_subset.rds')

FEcell_counts <- FetchData(fibroblasts_subset, vars = c("sample_type", "Subcluster",'sample')) %>%
  mutate(tissue_type = factor(sample_type, levels = c("tumor", "normal")))

colors <- c(
  "#C13F99", "#6A2473", "#331F40", "#FFDCE4", "#F2B28D",
  "#F24B59", "#A63F48", "#9373D9", "#68A6A6", "#83A660"
)

pdf('./res/fig/01_Fig1H_CAF_sampletype.pdf',width = 5,height = 1)
FEcell_counts %>%
  #filter(tissue_type == "Tumor") %>%
  ggplot(aes(x = sample_type, fill = Subcluster)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = rev(colors)) +
  coord_flip() +
  scale_y_reverse()
dev.off()

Idents(fibroblasts_subset) <-'Subcluster'
marker <- FindAllMarkers(object = fibroblasts_subset,assay = "RNA", only.pos = T, min.pct = 0.25, min.diff.pct = 0.25)



CAF_marker <- subset(marker,cluster=='CAF')

saveRDS(CAF_marker,'./res/data/CAF_marker.rds')
############################################################
#
#--------------marker基因可视化————————————————————————————————
#
############################################################

vln.df=as.data.frame(fibroblasts_subset[["RNA"]]$data[CAF_marker$gene[1:6],])
vln.df$gene=rownames(vln.df)
vln.df=melt(vln.df,id="gene")
colnames(vln.df)[c(2,3)]=c("CB","exp")

fibroblasts_subset@meta.data$CB <- rownames(fibroblasts_subset@meta.data)
anno=fibroblasts_subset@meta.data[,c("CB","Subcluster")]
vln.df=inner_join(vln.df,anno,by="CB")
vln.df$gene=factor(vln.df$gene,levels = CAF_marker$gene[1:6]) #为了控制画图的基因顺序

pdf('./res/fig/01Fig1I_CAF_marler_vlinplot.pdf',width = 3,height = 6)
vln.df%>%ggplot(aes(Subcluster,exp))+geom_violin(aes(fill=gene),scale = "width")+
  facet_grid(vln.df$gene~.,scales = "free_y")+
  scale_fill_brewer(palette = "Set3",direction = 1)+
  scale_x_discrete("")+scale_y_continuous("")+
  theme_bw()+
  theme(
    axis.text.x.bottom = element_text(angle = 45,hjust = 1,vjust = 1),
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
    legend.position = "none"
  )
dev.off()