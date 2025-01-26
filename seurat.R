######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

#install.packages("clustree")
#install.packages("Seurat")
#install.packages("harmony")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("GSVA")
#BiocManager::install("GSEABase")
#BiocManager::install("limma")
#BiocManager::install("SingleR")
#BiocManager::install("celldex")
#BiocManager::install("monocle")

#install.packages("devtools")
#devtools::install_github('immunogenomics/presto')




#######################01.????ǰ?ڴ????ͽ???#######################
#???ð?
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)
library(clustree)
library(harmony)

logFCfilter=1          #logFC?Ĺ???????
adjPvalFilter=0.05     #????????pvalue?Ĺ???????

#???ù???Ŀ¼
workDir="D:/216ST/05.Seurat"
setwd(workDir)

#??ȡ????
dirs=list.dirs(workDir)
dirs_sample=dirs[-1]
names(dirs_sample)=gsub(".+\\/(.+)", "\\1", dirs_sample)
counts <- Read10X(data.dir = dirs_sample)
pbmc = CreateSeuratObject(counts, min.cells=3, min.features=100)

#ʹ??PercentageFeatureSet??????????��???????İٷֱ?
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
#???ƻ?????????С????ͼ
pdf(file="01.featureViolin.pdf", width=10, height=6.5)
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
dev.off()
pbmc=subset(x = pbmc, subset = nFeature_RNA > 50 & percent.mt < 15)    #?????ݽ??й???

#???Ʋ????????????Ե?ͼ??
pdf(file="01.featureCor.pdf", width=13, height=7)
plot1 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size=1.5)
plot2 <- FeatureScatter(object = pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt",,pt.size=1.5)
CombinePlots(plots = list(plot1, plot2))
dev.off()

#?????ݽ??б?׼??
pbmc <- NormalizeData(object=pbmc, normalization.method="LogNormalize", scale.factor=10000)
#??ȡϸ????????ϵ???ϴ??Ļ???
pbmc <- FindVariableFeatures(object=pbmc, selection.method="vst", nfeatures=1500)
#????????????ͼ
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file="01.featureVar.pdf", width=10, height=6)
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
dev.off()



#######################02.PCA???ɷַ???#######################
##PCA????
scale.genes <-  rownames(pbmc)
pbmc=ScaleData(pbmc, features = scale.genes)        #PCA??ά֮ǰ?ı?׼Ԥ???���??
pbmc=RunPCA(object= pbmc, npcs=20, pc.genes=VariableFeatures(object=pbmc))     #PCA????
pbmc=RunHarmony(pbmc, "orig.ident")

#????ÿ??PCA?ɷֵ?????????
pdf(file="02.pcaGene.pdf", width=10, height=8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca", nfeatures=20)
dev.off()

#????PCAͼ??
pdf(file="02.PCA.pdf", width=7.5, height=5)
DimPlot(object=pbmc, reduction="pca")
dev.off()

#PCA????ͼ
pdf(file="02.pcaHeatmap.pdf", width=10, height=8)
DimHeatmap(object=pbmc, dims=1:4, cells=500, balanced=TRUE, nfeatures=30, ncol=2)
dev.off()

#?õ?ÿ??PC??pֵ?ֲ?
pbmc <- JackStraw(object=pbmc, num.replicate=100)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:20)
pdf(file="02.pcaJackStraw.pdf",width=8,height=6)
JackStrawPlot(object=pbmc, dims=1:20)
dev.off()



#######################03.ϸ????????????marker????#######################
##????????
pcSelect=20
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)     #?????ڽӾ???

#??ϸ??????,??ϸ????׼ģ?黯
pbmc <- FindClusters(pbmc, resolution=seq(0.5, 1.2, by=0.1))
pbmc <- FindClusters(object = pbmc, resolution=0.6)
#????????ͼ??
pdf(file="03.cluster.pdf", width=7, height=6)
pbmc <-RunUMAP(object = pbmc, dims = 1:pcSelect)        #UMAP????
DimPlot(pbmc, reduction = "umap", pt.size = 2, label = TRUE)      #UMAP???ӻ?
#pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect)             #TSNE????
#TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)     #TSNE???ӻ?
dev.off()
write.table(pbmc$seurat_clusters,file="03.Cluster.txt",quote=F,sep="\t",col.names=F)

##????ÿ???????Ĳ???????
pbmc.markers <- FindAllMarkers(object = pbmc,
                               only.pos = FALSE,
                               min.pct = 0.25,
                               logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.markers,file="03.clusterMarkers.txt",sep="\t",row.names=F,quote=F)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10, file = "FindAllMarkers_top10.csv")
#????marker??ÿ??????????ͼ
pdf(file="03.clusterHeatmap.pdf",width=15, height=15)
DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()
dev.off()



#######################04.SingleR R??ע??ϸ??????#######################
pbmc_for_SingleR <- GetAssayData(pbmc, layer="data")
clusters<-pbmc@meta.data$seurat_clusters
ref1=get(load("ref_Human_all.RData"))
ref2=get(load("ref_Hematopoietic.RData"))
ref3=get(load("DatabaseImmuneCellExpressionData.Rdata"))
ref4=get(load("BlueprintEncode_bpe.se_human.RData"))
ref5=get(load("HumanPrimaryCellAtlas_hpca.se_human.RData"))
ref6=get(load("MonacoImmuneData.Rdata"))
ref7=get(load("NovershternHematopoieticData.Rdata"))
singler=SingleR(test=pbmc_for_SingleR, ref =list(ref1, ref2, ref3, ref4, ref5, ref6, ref7),
                labels=list(ref1$label.main,ref2$label.main,ref3$label.main,ref4$label.main,ref5$label.main,ref6$label.main,ref7$label.main), clusters = clusters)

singler$labels=gsub("_", " ", singler$labels)
singler$labels=gsub("T cells, CD4\\+", "CD4+ T-cells", singler$labels)
singler$labels=gsub("T cells, CD8\\+", "CD8+ T-cells", singler$labels)
clusterAnn=as.data.frame(singler)
clusterAnn=cbind(id=row.names(clusterAnn), clusterAnn)
clusterAnn=clusterAnn[,c("id", "labels")]
write.table(clusterAnn,file="04.clusterAnn.txt",quote=F,sep="\t", row.names=F)
clusterAnn$labels <- c("T Cells",
                       "B Cells", 
                       "CD8+ T Cells", 
                       "B Cells", 
                       'B Cells', 
                       'Enterocytes', 
                       'Fibroblasts', 
                       'Epithelial Cells', 
                       'Dendritic Cells', 
                       'Proliferating Cells', 
                       'Macrophages', 
                       'Endothelial Cells', 
                       'Mast Cells', 
                       'Goblet Cells', 
                       'Epithelial Cells', 
                       'Neutrophils', 
                       'Enterocytes', 
                       'Schwann Cells', 
                       'Inflammatory Epithelial Cells', 
                       'Pericytes or Smooth Muscle Cells', 
                       'Enteroendocrine Cells')
write.table(clusterAnn,file="04.clusterAnn_GPT.txt",quote=F,sep="\t", row.names=F)

cellAnn=c()
for(i in 1:length(pbmc$seurat_clusters)){
  index=pbmc$seurat_clusters[i]
  cellAnn=c(cellAnn, clusterAnn[index,2])
}
cellAnnOut=cbind(names(pbmc$seurat_clusters), cellAnn)
colnames(cellAnnOut)=c("id", "labels")
write.table(cellAnnOut, file="04.cellAnn.txt", quote=F, sep="\t", row.names=F)

#ϸ??ע?ͺ??Ŀ??ӻ?
newLabels=clusterAnn$labels
names(newLabels)=levels(pbmc)
pbmc=RenameIdents(pbmc, groups)
pdf(file="04.cellAnn.pdf", width=15, height=8)
DimPlot(pbmc, reduction = "umap", pt.size = 2, label = TRUE)     #UMAP???ӻ?
#TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)           #TSNE???ӻ?
dev.off()
#???????ӻ?
Type=gsub("(.*?)\\..*", "\\1", colnames(pbmc))
names(Type)=colnames(pbmc)
pbmc=AddMetaData(object=pbmc, metadata=Type, col.name="Type")
pdf(file="04.group.cellAnn.pdf", width=30, height=8)
DimPlot(pbmc, reduction = "umap", pt.size = 2, label = TRUE, split.by="Type")       #UMAP???ӻ?
#TSNEPlot(object = pbmc, pt.size = 1, label = TRUE, split.by="Type")     #TSNE???ӻ?
dev.off()

celltype = data.frame(ClusterID=rownames(clusterAnn), celltype=clusterAnn$labels, stringsAsFactors = F)
pbmc@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  pbmc@meta.data[which(pbmc@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
Idents(pbmc)=pbmc$Type

#细胞比例柱状图

df <- table(pbmc@active.ident,pbmc@meta.data$orig.ident) %>% melt()
colnames(df) <- c("Cluster","Sample","Number")
df$Cluster <- factor(df$Cluster)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

sample_color <- col_vector[1:9]
ggplot(data = df, aes(x =Number, y = Sample, fill =  Cluster)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=col_vector[1:20]) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="Ratio",y="")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  )


#ϸ?????͵Ĳ???????
pbmc.markers=FindAllMarkers(object = pbmc,
                            only.pos = FALSE,
                            min.pct = 0.25,
                            logfc.threshold = 1)
sig.cellMarkers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.cellMarkers,file="04.cellMarkers.txt",sep="\t",row.names=F,quote=F)
write.csv(pbmc.markers, 'ALL.csv')
#ϸ??????????????
groups=gsub("(.*?)\\..*", "\\1", colnames(pbmc))
groups=paste0(groups, "_", cellAnn)
names(groups)=colnames(pbmc)
pbmc=AddMetaData(object=pbmc, metadata=groups, col.name="Type")
for(cellName in c("Goblet Cells",                 
                  "Inflammatory Epithelial Cells",
                  "Proliferating Cells",
                  "Neutrophils",
                  "Dendritic Cells",
                  "Enteroendocrine Cells")        ){
  conName=paste0("Control_", cellName)
  treatName=paste0("UC_", cellName)
  
  if( (length(groups[groups==conName])>10) & (length(groups[groups==treatName])>10) ){
    pbmc.markers=FindMarkers(pbmc, ident.1=treatName, ident.2=conName, logfc.threshold=0.1)
    sig.markersGroup=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
    sig.markersGroup=cbind(Gene=row.names(sig.markersGroup), sig.markersGroup)
    write.csv(sig.markersGroup,file=paste0("05.", cellName, ".diffGene.csv"))
  }
}


#???浥ϸ?????ݵĶ???
save(pbmc, cellAnn, file="Seurat.Rdata")

#load("Seurat.Rdata")

###################################05.monocle R??ϸ???켣????###################################
#׼??ϸ???켣??????Ҫ???ļ?
monocle.matrix=as.matrix(pbmc@assays$RNA$data)
monocle.sample=pbmc@meta.data
monocle.geneAnn=data.frame(gene_short_name=row.names(monocle.matrix), row.names=row.names(monocle.matrix))
monocle.clusterAnn=clusterAnn
monocle.markers=sig.markers

#??Seurat????ת??Ϊmonocle??Ҫ??ϸ????????ϸ??ע?ͱ????ͻ???ע?ͱ???
data <- as(as.matrix(monocle.matrix), 'sparseMatrix')
pd<-new("AnnotatedDataFrame", data = monocle.sample)
fd<-new("AnnotatedDataFrame", data = monocle.geneAnn)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd)
names(pData(cds))[names(pData(cds))=="seurat_clusters"]="Cluster"
pData(cds)[,"Cluster"]=paste0("cluster",pData(cds)[,"Cluster"])

#????ϸ??????????
clusterAnn=as.character(monocle.clusterAnn[,2])
names(clusterAnn)=paste0("cluster",monocle.clusterAnn[,1])
pData(cds)$cell_type2 <- plyr::revalue(as.character(pData(cds)$Cluster),clusterAnn)

#ϸ???켣????????
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, as.vector(sig.markers$gene))
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree')
cds <- orderCells(cds)
#????ϸ???켣״̬??ͼ??
pdf(file="06.trajectory.State.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "State")
dev.off()
#?????ֻ?ʱ????ϸ???켣ͼ
pdf(file="06.trajectory.Pseudotime.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "Pseudotime")
dev.off()
#????ϸ?????Ƶ?ϸ???켣ͼ
pdf(file="06.trajectory.cellType.pdf",width=6.5,height=6)
plot_cell_trajectory(cds,color_by = "cell_type2")
dev.off()
#??????????ϸ???켣ͼ
pdf(file="06.trajectory.cluster.pdf",width=6.5,height=6)
plot_cell_trajectory(cds, color_by = "Cluster")
dev.off()


######Video source: https://ke.biowolf.cn
######??????ѧ??: https://www.biowolf.cn/
######΢?Ź??ںţ?biowolf_cn
######???????䣺biowolf@foxmail.com
######????΢??: 18520221056

pbmc_epi <- subset(pbmc, Type=="UC_Epithelial Cells")
pbmc_t_Bcell <- subset(pbmc, subset=group=="Treat_B-cells")

pbmc_epi <- RunPCA(pbmc_epi, features = VariableFeatures(pbmc_epi)) 
plot1 <- DimPlot(pbmc_epi, reduction = "pca", group.by="orig.ident") 
plot2 <- ElbowPlot(pbmc_epi, ndims=40, reduction="pca") 
plotc <- plot1+plot2
plotc
ggsave("cluster/pca.pdf", plot = plotc, width = 8, height = 4) 

# 选取平缓的elbow，不用更改
pc.num=1:10

#reduction = "harmony" 
pbmc_epi <- FindNeighbors(pbmc_epi, dims = pc.num) 

# 聚类
pbmc_epi <- FindClusters(pbmc_epi)
table(pbmc_epi@meta.data$seurat_clusters)
metadata <- pbmc_epi@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'cluster/cell_cluster.csv',row.names = F)

#tSNE
pbmc_epi <- RunTSNE(pbmc_epi,dims = pc.num,dim.embed=2)
embed_tsne <- Embeddings(pbmc_epi, 'tsne')
write.csv(embed_tsne,'cluster/embed_tsne.csv') 
plot1 = DimPlot(pbmc_epi, reduction = "tsne",label = T) 
plot1
ggsave("cluster/tsne.pdf", plot = plot1, width = 8, height = 7)


#UMAP可视化
pbmc_epi <- RunUMAP(pbmc_epi, dims = pc.num)
embed_umap <- Embeddings(pbmc_epi, 'umap')
write.csv(embed_umap,'cluster/embed_umap.csv') 
plot2 = DimPlot(pbmc_epi, reduction = "umap",label = T) 
plot2
ggsave("pbmc_epi/UMAP.pdf", plot = plot2, width = 8, height = 7)

#Find markers for each cluster
pbmc.markers=FindAllMarkers(object = pbmc_epi,
                            only.pos = FALSE,
                            min.pct = 0.25,
                            logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,file="pbmc_epi/03.top10Markers.csv")

#????marker??ÿ??????????ͼ
pdf(file="03.clusterHeatmap.pdf",width=15, height=15)
DoHeatmap(object = pbmc_epi, features = top10$gene) + NoLegend()
dev.off()


pbmc_for_SingleR <- GetAssayData(pbmc_epi, layer="data")
clusters<-pbmc_epi@meta.data$seurat_clusters
ref1=get(load("ref_Human_all.RData"))
ref2=get(load("ref_Hematopoietic.RData"))
ref3=get(load("DatabaseImmuneCellExpressionData.Rdata"))
ref4=get(load("BlueprintEncode_bpe.se_human.RData"))
ref5=get(load("HumanPrimaryCellAtlas_hpca.se_human.RData"))
ref6=get(load("MonacoImmuneData.Rdata"))
ref7=get(load("NovershternHematopoieticData.Rdata"))
singler=SingleR(test=pbmc_for_SingleR, ref =list(ref1, ref2, ref3, ref4, ref5, ref6, ref7),
                labels=list(ref1$label.main,ref2$label.main,ref3$label.main,ref4$label.main,ref5$label.main,ref6$label.main,ref7$label.main), clusters = clusters)

singler$labels=gsub("_", " ", singler$labels)
singler$labels=gsub("T cells, CD4\\+", "CD4+ T-cells", singler$labels)
singler$labels=gsub("T cells, CD8\\+", "CD8+ T-cells", singler$labels)
clusterAnn=as.data.frame(singler)
clusterAnn=cbind(id=row.names(clusterAnn), clusterAnn)
clusterAnn=clusterAnn[,c("id", "labels")]
write.table(clusterAnn,file="04.clusterAnn.txt",quote=F,sep="\t", row.names=F)
clusterAnn$labels <- c('Enterocytes', 'Goblet Cells', 'Progenitor Cells', 'Tuft Cells')
write.table(clusterAnn,file="04.clusterAnn_GPT.txt",quote=F,sep="\t", row.names=F)

cellAnn=c()
for(i in 1:length(pbmc_epi$seurat_clusters)){
  index=pbmc_epi$seurat_clusters[i]
  cellAnn=c(cellAnn, clusterAnn[index,2])
}
cellAnnOut=cbind(names(pbmc_epi$seurat_clusters), cellAnn)
colnames(cellAnnOut)=c("id", "labels")
write.table(cellAnnOut, file="04.cellAnn.txt", quote=F, sep="\t", row.names=F)

#ϸ??ע?ͺ??Ŀ??ӻ?
newLabels=clusterAnn$labels
names(newLabels)=levels(pbmc_epi)
pbmc_epi=RenameIdents(pbmc_epi, newLabels)
pdf(file="pbmc_epi/04.cellAnn.pdf", width=10, height=6)
DimPlot(pbmc_epi, reduction = "umap", pt.size = 2, label = TRUE)     #UMAP???ӻ?
#TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)           #TSNE???ӻ?
dev.off()





pbmc_t_Bcell <- RunPCA(pbmc_t_Bcell, features = VariableFeatures(pbmc_t_Bcell)) 
plot1 <- DimPlot(pbmc_t_Bcell, reduction = "pca", group.by="orig.ident") 
plot2 <- ElbowPlot(pbmc_t_Bcell, ndims=40, reduction="pca") 
plotc <- plot1+plot2
plotc
ggsave("cluster/pca.pdf", plot = plotc, width = 8, height = 4) 

# 选取平缓的elbow，不用更改
pc.num=1:10

#reduction = "harmony" 
pbmc_t_Bcell <- FindNeighbors(pbmc_t_Bcell, dims = pc.num) 

# 聚类
pbmc_t_Bcell <- FindClusters(pbmc_t_Bcell)
table(pbmc_t_Bcell@meta.data$seurat_clusters)
metadata <- pbmc_t_Bcell@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)
write.csv(cell_cluster,'cluster/cell_cluster.csv',row.names = F)

#tSNE
pbmc_t_Bcell <- RunTSNE(pbmc_t_Bcell,dims = pc.num,dim.embed=2)
embed_tsne <- Embeddings(pbmc_t_Bcell, 'tsne')
write.csv(embed_tsne,'cluster/embed_tsne.csv') 
plot1 = DimPlot(pbmc_t_Bcell, reduction = "tsne",label = T) 
plot1
ggsave("cluster/tsne.pdf", plot = plot1, width = 8, height = 7)


#UMAP可视化
pbmc_t_Bcell <- RunUMAP(pbmc_t_Bcell, dims = pc.num)
embed_umap <- Embeddings(pbmc_t_Bcell, 'umap')
write.csv(embed_umap,'cluster/embed_umap.csv') 
plot2 = DimPlot(pbmc_t_Bcell, reduction = "umap",label = T) 
plot2
ggsave("pbmc_t_Bcell/UMAP.pdf", plot = plot2, width = 8, height = 7)

#Find markers for each cluster
pbmc.markers=FindAllMarkers(object = pbmc_t_Bcell,
                            only.pos = FALSE,
                            min.pct = 0.25,
                            logfc.threshold = logFCfilter)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10,file="pbmc_t_Bcell/03.top10Markers.csv")

#????marker??ÿ??????????ͼ
pdf(file="03.clusterHeatmap.pdf",width=15, height=15)
DoHeatmap(object = pbmc_t_Bcell, features = top10$gene) + NoLegend()
dev.off()


pbmc_for_SingleR <- GetAssayData(pbmc_t_Bcell, layer="data")
clusters<-pbmc_t_Bcell@meta.data$seurat_clusters
ref1=get(load("ref_Human_all.RData"))
ref2=get(load("ref_Hematopoietic.RData"))
ref3=get(load("DatabaseImmuneCellExpressionData.Rdata"))
ref4=get(load("BlueprintEncode_bpe.se_human.RData"))
ref5=get(load("HumanPrimaryCellAtlas_hpca.se_human.RData"))
ref6=get(load("MonacoImmuneData.Rdata"))
ref7=get(load("NovershternHematopoieticData.Rdata"))
singler=SingleR(test=pbmc_for_SingleR, ref =list(ref1, ref2, ref3, ref4, ref5, ref6, ref7),
                labels=list(ref1$label.main,ref2$label.main,ref3$label.main,ref4$label.main,ref5$label.main,ref6$label.main,ref7$label.main), clusters = clusters)

singler$labels=gsub("_", " ", singler$labels)
singler$labels=gsub("T cells, CD4\\+", "CD4+ T-cells", singler$labels)
singler$labels=gsub("T cells, CD8\\+", "CD8+ T-cells", singler$labels)
clusterAnn=as.data.frame(singler)
clusterAnn=cbind(id=row.names(clusterAnn), clusterAnn)
clusterAnn=clusterAnn[,c("id", "labels")]
write.table(clusterAnn,file="04.clusterAnn.txt",quote=F,sep="\t", row.names=F)
clusterAnn$labels <- c("Plasma Cells",
                       "Memory B Cells", 
                       "Plasma Cells", 
                       "Germinal Center B Cells", 
                       'Plasma Cells', 
                       'Proliferating B Cells', 
                       'Plasma Cells', 
                       'Plasma Cells', 
                       'Plasma Cells', 
                       'Plasma Cells', 
                       'Plasma Cells', 
                       'Plasma Cells', 
                       'Memory B Cells', 
                       'Plasma Cells', 
                       'Activated B Cells',
                       'Plasma Cells')
write.table(clusterAnn,file="04.clusterAnn_GPT.txt",quote=F,sep="\t", row.names=F)

cellAnn=c()
for(i in 1:length(pbmc_t_Bcell$seurat_clusters)){
  index=pbmc_t_Bcell$seurat_clusters[i]
  cellAnn=c(cellAnn, clusterAnn[index,2])
}
cellAnnOut=cbind(names(pbmc_t_Bcell$seurat_clusters), cellAnn)
colnames(cellAnnOut)=c("id", "labels")
write.table(cellAnnOut, file="04.cellAnn.txt", quote=F, sep="\t", row.names=F)

#ϸ??ע?ͺ??Ŀ??ӻ?
newLabels=clusterAnn$labels
names(newLabels)=levels(pbmc_t_Bcell)
pbmc_t_Bcell=RenameIdents(pbmc_t_Bcell, newLabels)
pdf(file="04.cellAnn.pdf", width=10, height=6)
DimPlot(pbmc_t_Bcell, reduction = "umap", pt.size = 2, label = TRUE)     #UMAP???ӻ?
#TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)           #TSNE???ӻ?
dev.off()
