
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})

library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
library(org.Hs.eg.db)
library(patchwork)
setwd("E:/lab/CA3/CA3-GSEA/fudan")
#输入文件
expFile = read.csv("exp-Tumor.csv",header=TRUE,sep=",")
       
sgene="CA3"       #输入进行单基因GSEA的基因名称

rt=expFile
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#按基因表达分组
group <- ifelse(data[c(sgene),]>median(data[c(sgene),]), "High", "Low")   
group <- factor(group,levels = c("High","Low"))

#差异分析
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(High-Low,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
deg=topTable(fit2,adjust='fdr',number=nrow(data))
Diff=deg
#保存单基因分组的所有基因差异结果
DIFFOUT=rbind(id=colnames(Diff),Diff)
write.table(DIFFOUT,file=paste0("1.","DIFF_all-CA3.xls"),sep="\t",quote=F,col.names=F)




logFC_t=0
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)

deg$symbol=rownames(deg)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
DEG<- deg[deg$P.Value < 0.05, ]
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
data_all_sort <- DEG %>% 
  arrange(desc(logFC))

geneList = data_all_sort$logFC #把foldchange按照从大到小提取出来
names(geneList) <- data_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID
#names(geneList) <- data_all_sort$symbol
head(geneList)


#GSEA分析——GO
Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont="all", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)

Go_gse <- as.data.frame(Go_gseresult)

#基因名称转换，例如将基因 ENTREZ ID 转换为 SYMBOL ID
GO <- Go_gse
rownames(data_all_sort) <- as.character(data_all_sort$ENTREZ)

for (i in 1:nrow(GO)) {
  new <- c()
  old <- unlist(strsplit(GO[i,'core_enrichment'], '/'))
  for (id in old) new <- c(new, data_all_sort[id,'symbol'])
  GO[i,'core_enrichment'] <- paste(new, collapse = '/')
}


write.table (GO, file ="fudan_Go_gseresult.csv", sep =",", row.names =TRUE)


# 假设 GO 是你的数据框
matching_rows <- grep("oxidoreductase activity", GO$Description, ignore.case = TRUE)

print(matching_rows)
integer(0)
# 假设 GO 是你的数据框
matching_rows <- grep("cell redox homeostasis", GO$Description, ignore.case = TRUE)

print(matching_rows)
integer(0)

dev.off()
gseaplot2(Go_gseresult,1009,color="red",pvalue_table = T)

# 获取 NES 值
nes <- Go_gseresult$NES[3737]


#==================KEGG================================
#开始富集分析
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               nPerm        = 10000,
               minGSSize    = 10,
               maxGSSize    = 200,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none" )
class(kk2)
colnames(kk2@result)
kegg_result <- as.data.frame(kk2)
rownames(kk2@result)[head(order(kk2@result$enrichmentScore))]
af=as.data.frame(kk2@result)
rownames(df) <- as.character(df$ENTREZ)

for (i in 1:nrow(af)) {
  new <- c()
  old <- unlist(strsplit(af[i,'core_enrichment'], '/'))
  for (id in old) new <- c(new, df[id,'SYMBOL'])
  af[i,'core_enrichment'] <- paste(new, collapse = '/')
}


write.table(af,file=paste0("fudan_gseKEGG-CA3.xls"),sep="\t",quote=F,col.names=T)

#排序后分别取GSEA结果的前5个和后5个
num=5
pdf(paste0("2.","down_GSEA-CA3.pdf"),width = 8,height = 8)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[head(order(kk2@result$enrichmentScore),num)])
dev.off()
pdf(paste0("2.","up_GSEA-CA3.pdf"),width = 8,height = 8)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[tail(order(kk2@result$enrichmentScore),num)])
dev.off()
#排序后取前5个和后5个一起展示
num=5
pdf(paste0("2.","CA3-all_GSEA.pdf"),width = 10,height =10)
p<-gseaplot2(kk2, geneSetID = rownames(kk2@result)[c(head(order(kk2@result$enrichmentScore),num),tail(order(kk2@result$enrichmentScore),num))])
p

dev.new()

#单独展示,自行保存
gseaplot2(kk2,
          title = "Th1 and Th2 cell differentiation",  #设置标题
          "hsa04658", #绘制hsa04658通路的结果，通路名称与编号对应
          color="red", #线条颜色
          base_size = 20, #基础字体的大小
          subplots = 1:3, 
          pvalue_table = T) # 显示p值

