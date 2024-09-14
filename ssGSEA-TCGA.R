######ssGSEA
rm(list=ls())
library(tidyverse)
library(GSVA)
library(dplyr)
setwd("E:/lab/CA3/ssGSEA")


DEG<-read.csv("lihc_tcga_exp.csv",comment.char="!",
              stringsAsFactors=F,
              sep=",",
              header=T,row.names = 1)
DEG$ENSEMBL=rownames(DEG)
df <- bitr(unique(DEG$ENSEMBL), fromType = "ENSEMBL",
           toType = c( "SYMBOL"),
           OrgDb = org.Hs.eg.db)
deg=DEG
deg=merge(deg,df,by='ENSEMBL')
deg<-deg[,-1]
# 找到重复的SYMBOL值
duplicate_symbols <- deg$SYMBOL[duplicated(deg$SYMBOL)]

# 删除重复的行
deg <- deg[!deg$SYMBOL %in% duplicate_symbols, ]
row.names(deg) <- deg$SYMBOL
deg <- deg[, -which(names(deg) == "SYMBOL")]
write.csv(deg, "lihc_tcga_exp-419.csv")



deg<-read.csv("lihc_tcga_exp-419.csv",comment.char="!",
              stringsAsFactors=F,
              sep=",",
              header=T,row.names = 1)

# 然后重新执行选择列的代码
selected_columns <- grep("01", colnames(deg), value = TRUE)
tumor <- deg[, selected_columns]
# 或者
#tumor <- deg %>% select(all_of(selected_columns))
#tumor=deg %>% select(matches("01"))

write.csv(tumor, "lihc_tcga_exp-369.csv")



## 去除低表达的探针
dat<- tumor[!apply(tumor,1,sum)==0,]
## 过滤重复的probe
boxplot(dat[,1:369],las=2)
gs = clusterProfiler::read.gmt("KEGG_NITROGEN_METABOLISM.v2023.2.Hs.gmt")
a = as.matrix(dat)
gs = split(gs$gene, gs$term)
## build GSVA parameter object
gsvapar <- gsvaParam(a, gs, maxDiff=TRUE)

## estimate GSVA enrichment scores for the three sets
gsva_es <- gsva(gsvapar)


write.csv(gsva_es, "TCGA-LIHC-ssGSEA-CA3.csv")



#基因集需要是list为对象。默认情况下，kcdf="Gaussian"，适用于输入表达式值连续的情况，如对数尺度的微阵列荧光单元、RNA-seq log-CPMs、log-RPKMs或log-TPMs。当输入表达式值是整数计数时，比如那些从RNA-seq实验中得到的值，那么这个参数应该设置为kcdf="Poisson"
ssgsea<-read.csv("TCGA-LIHC-ssGSEA-CA3.csv",comment.char="!",
                 stringsAsFactors=F,
                 sep=",",
                 header=T,row.names = 1)
# Min-Max标准化是指对原始数据进行线性变换，将值映射到[0，1]之间
# 这里是将每个样本中不同的免疫细胞比例标准化到0-1之间
ssgsea.1 <- as.data.frame(t(ssgsea))
# 最小-最大归一化函数
min_max_normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

# 对 Score 列进行归一化
ssgsea.1$KEGG_NITROGEN_METABOLISM <- min_max_normalize(ssgsea.1$KEGG_NITROGEN_METABOLISM)

write.csv(ssgsea.1, "TCGA-LIHC-ssGSEA-CA3-2.csv")




