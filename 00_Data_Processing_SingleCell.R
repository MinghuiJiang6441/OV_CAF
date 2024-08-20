setwd('~/NSCLC/Test/')
library(Seurat)
library(dplyr)
library(gridExtra)
library(paletteer)
library(harmony)

unzip('~/Test/GSE184880_RAW.zip', exdir = '~/Test')

files <- list.files('~/Test/GSE184880_RAW', full.names = TRUE)
files

# 定义文件路径
file_paths <- list.files('/data/home/jiangminghui/Test/GSE184880_RAW', full.names = TRUE)


# 创建一个空列表来存储Seurat对象
seurat_list <- list()

# 循环读取每个样本的数据
for (i in seq(1, length(file_paths), by = 3)) {
  # 读取barcodes, genes, 和 matrix文件
  barcodes <- read.table(file_paths[i], header = FALSE)
  genes <- read.table(file_paths[i + 1], header = FALSE)
  matrix <- Matrix::readMM(file_paths[i + 2])
  
  # 设置matrix的行名和列名
  rownames(matrix) <- make.unique(genes$V2)
  colnames(matrix) <- barcodes$V1
  
  # 创建Seurat对象
  seurat_obj <- CreateSeuratObject(counts = matrix, project = "GSE184880", min.cells = 3, min.features = 200)
  
  # 添加样本名称
  sample_name <- strsplit(basename(file_paths[i]), "_")[[1]][1]
  seurat_obj$sample <- sample_name
  
  # 将Seurat对象添加到列表中
  seurat_list[[sample_name]] <- seurat_obj
}

# 合并所有Seurat对象
OV_seurat <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = names(seurat_list), project = "GSE184880")

# 计算线粒体基因比例
OV_seurat[["percent.mt"]] <- PercentageFeatureSet(OV_seurat, pattern = "^MT-")

# 绘制质量控制图
VlnPlot(OV_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# 设置过滤标准
OV_seurat <- subset(OV_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)


# 归一化数据
OV_seurat <- NormalizeData(OV_seurat) %>%FindVariableFeatures(selection.method = "vst",verbose = TRUE) 


top10 <- head(VariableFeatures(OV_seurat), 10)# 筛选前10个高变基因
top10
plot1 <- VariableFeaturePlot(OV_seurat) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE,xnudge = 0,ynudge = 0) # 额外标记前10个HVG。

pdf('./res/fig/00_Fig1A_VariableFeatures.pdf',width = 6,height = 7)
grid.arrange(plot1, plot2, ncol = 1)
dev.off()

OV_seurat<- ScaleData(OV_seurat) %>%RunPCA()

OV_seurat <-RunHarmony(OV_seurat,reduction = 'pca',lambda = 1,group.by.vars = 'sample')


VizDimLoadings(OV_seurat, dims = 1:10, reduction = "pca")

DimPlot(OV_seurat, reduction = "pca")

DimHeatmap(OV_seurat, dims = 1:10, cells = 500, balanced = TRUE)
OV_seurat <- JackStraw(OV_seurat, num.replicate = 100)
OV_seurat <- ScoreJackStraw(OV_seurat, dims = 1:5)
JackStrawPlot(OV_seurat,xmax = 0.01, ymax = 0.1,dims = 1:5) # 将重要的PC的P值进行可视化，通常情况下此步骤可以跳跃，选择不同数量的PC结果会有很大的不同，且建议使用更多的PC进行下游分析，如果选择过少，则会丢失较多的特征值，对分析产生不利影响。通常默认选择20个，可选择10~50个，结果通常产生不了太大变化。
ElbowPlot(OV_seurat)
pdf('./res/fig/00_Fig1C_JackStrawPlot.pdf',width = 5,height = 5)
JackStrawPlot(OV_seurat)
dev.off()

pdf('./res/fig/00_Fig1D_ElbowPlot.pdf',width = 5,height = 5)
ElbowPlot(OV_seurat)
dev.off()

# 聚类分析
OV_seurat <- FindNeighbors(OV_seurat, dims = 1:10,reduction = 'pca')
OV_seurat <- FindClusters(OV_seurat, resolution = 0.5)

# UMAP降维
OV_seurat <- RunUMAP(OV_seurat, dims = 1:10)

d_palettes<- palettes_d_names

mycol<-paletteer_d( "ggsci::default_igv",n=51)
set.seed(123) # 为了可重复性设置随机种子
random_colors <- sample(mycol, 12)

pdf('./res/fig/00_Fig1E_DimPlot.pdf',width = 8,height = 4)
# 可视化结果
DimPlot(OV_seurat, reduction = "umap", raster = F,group.by = c('seurat_clusters',"sample"),cols =mycol )
dev.off()

saveRDS(OV_seurat,'./data/OV_seurat.rds')

