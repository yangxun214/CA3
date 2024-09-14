######ssGSEA
rm(list=ls())
library(tidyverse)
library(GSVA)
library(dplyr)
library(org.Hs.eg.db)

setwd("E:/lab/CA3/ssGSEA/pr")

DEG<-read.csv("LIHC-CPTAC_Log Ratio150.csv",comment.char="!",
              stringsAsFactors=F,
              sep=",",
              header=T,row.names = 1)

deg<-read.csv("LIHC-CPTAC_Log Ratio150.csv",comment.char="!",
              stringsAsFactors=F,
              sep=",",
              header=T,row.names = 1)

tumor <- deg


## 去除低表达的探针
dat<- tumor[!apply(tumor,1,sum)==0,]
## 过滤重复的probe
boxplot(dat[,1:150],las=2)
gs = clusterProfiler::read.gmt("HP_INCREASED_SERUM_LACTATE.v2023.2.Hs.gmt")
a = as.matrix(dat)
gs = split(gs$gene, gs$term)
## build GSVA parameter object
gsvapar <- gsvaParam(a, gs, maxDiff=TRUE)

## estimate GSVA enrichment scores for the three sets
gsva_es <- gsva(gsvapar)


write.csv(gsva_es, "pr-LIHC-ssGSEA-lactate.csv")



#基因集需要是list为对象。默认情况下，kcdf="Gaussian"，适用于输入表达式值连续的情况，如对数尺度的微阵列荧光单元、RNA-seq log-CPMs、log-RPKMs或log-TPMs。当输入表达式值是整数计数时，比如那些从RNA-seq实验中得到的值，那么这个参数应该设置为kcdf="Poisson"
ssgsea<-read.csv("pr-LIHC-ssGSEA-NITROGEN.csv",comment.char="!",
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
ssgsea.1$HP_INCREASED_SERUM_LACTATE <- min_max_normalize(ssgsea.1$HP_INCREASED_SERUM_LACTATE)

write.csv(ssgsea.1, "TCGA-LIHC-ssGSEA-CA3-2.csv")




