BiocManager::install("hgu133plus2.db")
library("annotate")
library("hgu133plus2.db")
BiocManager::install ("GO.db")
library(GO.db)
install.packages("dplyr")
library(dplyr)
library(AnnotationDbi)
BiocManager::install("clusterProfiler")
library(clusterProfiler)
BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(ggplot2)

data <- read.csv("C:/Users/Khalila/Downloads/GO2.csv")
View(data)

cluster_0 <- data$X[data$Cluster_HDBSCAN == 0]
mapped0 <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = cluster_0,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)
gene_list0 <- unique(mapped0$ENTREZID)
ego0 <- enrichGO(
  gene          = gene_list0,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)
p0 <- barplot(ego0, showCategory = 10, title = "GO Enrichment - Cellular Component (Cluster 0)")
p0 + theme(axis.text.y = element_text(size = 9))

cluster_1 <- data$X[data$Cluster_HDBSCAN == 1]
mapped1 <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = cluster_1,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)
gene_list1 <- unique(mapped1$ENTREZID)
ego1 <- enrichGO(
  gene          = gene_list1,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "cc",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)
p1<- barplot(ego1, showCategory = 10, title = "GO Enrichment - Celluler Component (Cluster 1)")
p1 + theme(axis.text.y = element_text(size = 9))

##

cluster_0_bi <- data$X[data$Bicluster == 0]
mapped0_bi <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = cluster_0_bi,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)
gene_list0_bi <- unique(mapped0_bi$ENTREZID)
ego0_bi <- enrichGO(
  gene          = gene_list0_bi,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)
p0_bi <- barplot(ego0_bi, showCategory = 10, title = "GO Enrichment - Cellular Component (Bicluster 0)")
p0_bi + theme(axis.text.y = element_text(size = 9))

cluster_1_bi <- data$X[data$Bicluster == 1]
mapped1_bi <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = cluster_1_bi,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)
gene_list1_bi <- unique(mapped1_bi$ENTREZID)
ego1_bi <- enrichGO(
  gene          = gene_list1_bi,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)
p1_bi <- barplot(ego1_bi, showCategory = 10, title = "GO Enrichment - Cellular Component (Bicluster 1)")
p1_bi + theme(axis.text.y = element_text(size = 9))

cluster_2_bi <- data$X[data$Bicluster == 2]
mapped2_bi <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = cluster_2_bi,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)
gene_list2_bi <- unique(mapped2_bi$ENTREZID)
ego2_bi <- enrichGO(
  gene          = gene_list2_bi,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)
p2_bi <- barplot(ego2_bi, showCategory = 10, title = "GO Enrichment - Cellular Component (Bicluster 2)")
p2_bi + theme(axis.text.y = element_text(size = 9))

cluster_3_bi <- data$X[data$Bicluster == 3]
mapped3_bi <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = cluster_3_bi,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)
gene_list3_bi <- unique(mapped3_bi$ENTREZID)
ego3_bi <- enrichGO(
  gene          = gene_list3_bi,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)
p3_bi <- barplot(ego3_bi, showCategory = 10, title = "GO Enrichment - Cellular Component (Bicluster 3)")
p3_bi + theme(axis.text.y = element_text(size = 9))

cluster_4_bi <- data$X[data$Bicluster == 4]
mapped4_bi <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = cluster_4_bi,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)
gene_list4_bi <- unique(mapped4_bi$ENTREZID)
ego4_bi <- enrichGO(
  gene          = gene_list4_bi,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)
p4_bi <- barplot(ego4_bi, showCategory = 10, title = "GO Enrichment - Cellular Component (Bicluster 4)")
p4_bi + theme(axis.text.y = element_text(size = 9))

cluster_5_bi <- data$X[data$Bicluster == 5]
mapped5_bi <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = cluster_5_bi,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)
gene_list5_bi <- unique(mapped5_bi$ENTREZID)
ego5_bi <- enrichGO(
  gene          = gene_list5_bi,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)
p5_bi <- barplot(ego5_bi, showCategory = 10, title = "GO Enrichment - Biological Process (Bicluster 5)")
p5_bi + theme(axis.text.y = element_text(size = 9))

cluster_6_bi <- data$X[data$Bicluster == 6]
mapped6_bi <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = cluster_6_bi,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)
gene_list6_bi <- unique(mapped6_bi$ENTREZID)
ego6_bi <- enrichGO(
  gene          = gene_list6_bi,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)
p6_bi <- barplot(ego6_bi, showCategory = 10, title = "GO Enrichment - Biological Process (Bicluster 6)")
p6_bi + theme(axis.text.y = element_text(size = 9))

cluster_7_bi <- data$X[data$Bicluster == 7]
mapped7_bi <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = cluster_7_bi,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
)
gene_list7_bi <- unique(mapped7_bi$ENTREZID)
ego7_bi <- enrichGO(
  gene          = gene_list7_bi,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)
p7_bi <- barplot(ego7_bi, showCategory = 10, title = "GO Enrichment - Biological Process (Bicluster 7)")
p7_bi + theme(axis.text.y = element_text(size = 9))

