library(easyTCGA)
library(GSVA)
library(GSEABase)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)

setwd('~/NSCLC/Test/')

getmrnaexpr("TCGA-OV")
OV_SC_Fibro_marker <-readRDS('./res/data/OV_SC_Fibro_marker.rds')

load(file = "~/NSCLC/Test/output_mRNA_lncRNA_expr/TCGA-OV_mrna_expr_count.rdata")
#load('~/NSCLC/Test/output_mRNA_lncRNA_expr/TCGA-OV_clinicalSE.rdata')
#load('~/NSCLC/Test/output_mRNA_lncRNA_expr/TCGA-OV_gene_info.rdata')
OV_clin <- read.delim('./data/OV_clin.txt', header = TRUE, row.names = 1, sep = '\t')
OV_count <-mrna_expr_count
# 创建一个映射表，将hgnc_id映射到gene_name
############################################################
#
#--------------生存分析————————————————————————————————
#
############################################################
colnames(OV_count) %in%rownames(OV_clin) %>% any()
rownames(OV_clin)[1:5]
colnames(OV_count)[1:5]
colnames(OV_count)<-substr(colnames(OV_count),start = 1,stop = 12)

bothsamples <- intersect(colnames(OV_count),rownames(OV_clin))
length(bothsamples)
OV_count <-OV_count[,bothsamples]
OV_clin<-OV_clin[bothsamples,]
colnames(OV_count) %in%rownames(OV_clin) %>% any()


OV_OS_clin<- data.frame(OS.time = OV_clin$A1_OS, OS = OV_clin$A2_Event)
head(OV_OS_clin)



OV_OS_clin$OS.time <- as.numeric(OV_OS_clin$OS.time)
OV_OS_clin$OS<- ifelse(OV_OS_clin$OS== 'Alive',0,1)



gene_set <- GeneSet(CAF_marker$gene, setName = "CAF_sig")

# 将GeneSet对象转换为GeneSetCollection
gene_set_collection <- GeneSetCollection(list(gene_set))

ssgsea_scores <- gsva(as.matrix(OV_count), gene_set_collection, method = "ssgsea")

group <- ifelse(ssgsea_scores > median(ssgsea_scores), "High", "Low") %>%as.character()

# 将分组信息添加到临床数据中
OV_OS_clin$CAF_group <- group

OV_OS_clin <- OV_OS_clin[!is.na(OV_OS_clin$OS.time), ]

table(OV_OS_clin$OS)
surv_time <- OV_OS_clin$OS.time
surv_status <- OV_OS_clin$OS
# 创建生存对象
surv_object <- Surv(time = surv_time, event = surv_status)
# 进行生存分析
fit <- survfit(surv_object ~ OV_OS_clin$CAF_group, data = OV_OS_clin)

ggsurvplot(fit, data = OV_OS_clin, pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, legend.labs = c("High", "Low"),
           palette = c("#E7B800", "#2E9FDF"))

#p值不显著
save(OV_clin,file = './data/OV_clin.rdata')
save(OV_count,file = './data/OV_count.rdata')
save(OV_OS_clin,file = './data/OV_OS_clin.rdata')

############################################################
#
#--------------差异表达————————————————————————————————
#
############################################################
DEG<- read.table('./data/OV_DEG.txt',header = T,sep = '\t')

head(DEG)
colnames(DEG) <-c('gene','ESNG',"T",'N','Log2.Fold.Change','adjp')

range(DEG$Log2.Fold.Change.)


cut_off_padj =0.01 #设置padj的阈值

cut_off_log2FC =1 #设置log2FC的阈值

DEG$Sig = ifelse(DEG$adjp < cut_off_padj &    #根据阈值筛选差异显著的上下调基因，与差异不显著的基因
                      abs(DEG$Log2.Fold.Change) >= cut_off_log2FC,  #abs绝对值
                    ifelse(DEG$Log2.Fold.Change > cut_off_log2FC ,'Up','Down'),'no')
head(DEG)

# 创建火山图
volcano_plot <- ggplot(DEG, aes(x = Log2.Fold.Change, y = -log10(adjp))) +
  geom_point(aes(color = Sig), alpha = 0.8, size = 2, shape = 17) +  # shape = 17 表示三角形
  scale_color_manual(values =c("Up" = "#FF5A33", "Down" = "#20854e")) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-log10 Adjusted p-value") +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(),  # 去掉主网格线
        panel.grid.minor = element_blank())  # 去掉次网格线

# 显示火山图
print(volcano_plot)
ggsave('./res/fig/02_Fig2A_OV_TCGA_DEG.pdf',volcano_plot,width = 5,height = 5)















