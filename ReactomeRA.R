library(ReactomePA)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
BiocManager::install("org.Hs.eg.db")
BiocManager::install("org.Mm.eg.db")
genelist_input <- fread(file="C:/Users/ShiqiHu/Desktop/cluster1.txt", header = T, sep='\t', data.table = F)
genename <- as.character(genelist_input[,1]) 
gene_map <- select(org.Mm.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))
###or
gene_map <- AnnotationDbi::select(org.Mm.eg.db, keys=genename, keytype="SYMBOL", columns=c("SYMBOL","ENTREZID"))
##
colnames(gene_map)[1]<-"Gene"
write.csv(as.data.frame(gene_map),"C:/Users/ShiqiHu/Desktop/cluster1entrezid.csv",row.names =F)


aaa<-inner_join(gene_map,genelist_input,by = "Gene")
aaa<-aaa[,-1]
aaa<-na.omit(aaa)
aaa$logFC<-sort(aaa$logFC,decreasing = T)


geneList = aaa[,2]
names(geneList) = as.character(aaa[,1])
geneList= sort(geneList, decreasing = TRUE)

geneList

#GSEA——GO
Go_gseresult <- gseGO(geneList, 'org.Mm.eg.db', keyType = "ENTREZID", ont="all", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
#GSEA——KEGG
KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
#GSEA——Reactome
Go_Reactomeresult <- gsePathway(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
#
write.table (Go_gseresult, file ="C:/Users/ShiqiHu/Desktop/Go_gseresult.csv", sep =",", row.names =TRUE)
write.table (KEGG_gseresult, file ="C:/Users/ShiqiHu/Desktop/KEGG_gseresult.csv", sep =",", row.names =TRUE)
write.table (Go_Reactomeresult, file ="C:/Users/ShiqiHu/Desktop/Go_Reactomeresult.csv", sep =",", row.names =TRUE)

ridgeplot(Go_gseresult,10)
gseaplot(Go_gseresult,1,pvalue_table = TRUE)
gseaplot2(Go_Reactomeresult,212,pvalue_table = TRUE)
gseaplot2(Go_gseresult, 1:4, pvalue_table = TRUE)





devtools::install_github("nicolash2/gggsea")
BiocManager::install ("GSVA")
BiocManager :: install ("msigdb")
#加载包
library(clusterProfiler)
library(gggsea)
library(ggplot2)
library(GSVA)
library(enrichplot)
library('GSEABase')
library(fgsea)
library(org.Hs.eg.db)
library(tidyverse)
library(dplyr)
library(msigdb) 
library(GSEABase) 
aaa<-inner_join(gene_map,genelist_input,by = "Gene")
aaa<-aaa[,-1]
aaa<-na.omit(aaa)
aaa$logFC<-sort(aaa$logFC,decreasing = T)


geneList = aaa[,2]
names(geneList) = as.character(aaa[,1])
geneList= sort(geneList, decreasing = TRUE)


#读取基因集（本示例使用的 MSigDB 数据库中的 hallmark）
#使用clusterProfiler中的read.gmt函数读取下载的gmt基因集文件
hallmark <- read.gmt("h.all.v2022.1.Hs.symbols.gmt")



msigdb.mm  <-  getMsigdb(org = 'mm',id = c("SYM", "EZID"))
msigdb.mm  <-  appendKEGG(msigdb.mm)
listCollections(msigdb.mm)
listSubCollections(msigdb.mm)
hallmarks  <-  subsetCollection(msigdb.mm, 'h')
msigdb_ids <- geneIds(hallmarks)

gsea.re1<- GSEA(genename,  
                TERM2GENE = hallmarks,  
                pvalueCutoff = 1) 
