## install.packages("GeneNMF") ## 需要4.3.0以上R版本哦！
library(GeneNMF)
library(Seurat)
library(ggplot2)
library(patchwork)
library(Matrix)
library(RcppML)
library(viridis)
library(msigdbr)
library(fgsea)
library(readr)


###数据准备###

###数据准备###
#scRNA <- readRDS("scRNA.RDS") #单细胞数据
#提取肿瘤上皮
#pbmc@meta.data$Type <- gsub("Control", "CT", pbmc@meta.data$Type)
setwd("D:/216ST/05.Seurat/GeneNMF/UC")
scRNA_epi=subset(pbmc,celltype=='Schwann Cells')
scRNA_epi=subset(scRNA_epi,Type=="UC")
#合并前各个样本是分开的
scRNA_epi <- JoinLayers(scRNA_epi)
Layers(scRNA_epi, assay = "RNA") #合并后便解决了后面assay、slot的读取问题

seu.list <- SplitObject(scRNA_epi,split.by = "orig.ident") #先按照样本分开，我们有6个肿瘤组织
do.call(rbind,lapply(seu.list, dim)) #查看每个样本的基因数和细胞数
geneNMF.programs <- multiNMF(seu.list, assay="RNA", k=4:9, min.exp = 0.05,seed=2025)
geneNMF.metaprograms<-getMetaPrograms(geneNMF.programs,
                                        metric="cosine",
                                        weight.explained=0.5,
                                        nMP=10)

str(geneNMF.metaprograms$metaprograms.genes)

ph<-plotMetaPrograms(geneNMF.metaprograms,
                       similarity.cutoff=c(0.1,1))
ph

#看看每个MP的稳定性
geneNMF.metaprograms$metaprograms.metrics
##每个基因对MP的贡献度
geneNMF.metaprograms$metaprograms.genes.weights$MP6
##基因列表


top_p<-lapply(geneNMF.metaprograms$metaprograms.genes,
               function(program){
                 runGSEA(genes=program,universe=rownames(scRNA_epi),category="C5")
               })
#runGSEA()函数可用于链接msigDB数据库并评估检测到的MPs参与啥通路

head(top_p$MP1)

#nature绘图
J<-geneNMF.metaprograms[["programs.similarity"]]
tree<-geneNMF.metaprograms[["programs.tree"]]
cl_members<-geneNMF.metaprograms[["programs.clusters"]]
labs.order<-labels(as.dendrogram(tree))

## 获取样本顺序，照着跑就行
downsample=1000
#downsample, to avoid overloading the graphics
if(length(cl_members)>downsample){
  sid.keep<-downsampleMin(cl_members,size=downsample)
  labs.order<-labs.order[labs.order%in%sid.keep]
  
  cl_members<-cl_members[labs.order]
  J<-J[labs.order,labs.order]
  
  #disable 'tree' for downsampled heatmap
  showtree<-FALSE
}


cl_names<-names(cl_members)
cl_members[!is.na(cl_members)]<-paste0("MP",cl_members[!is.na(cl_members)])
names(cl_members)<-cl_names

#Recover order of MP clusters
cluster.order<-unique(cl_members[labs.order])
nMP<-length(cluster.order)

#Gaps for heatmap
diffs<-diff(as.integer(as.factor(cl_members)))
gaps<-which(diffs!=0)

## 制作注释文件，照着跑就行
#Annotation column
annotation_col<-as.data.frame(cl_members)
colnames(annotation_col)<-"Metaprogram"
annotation_col[["Metaprogram"]]<-factor(cl_members,levels=cluster.order)

### 这里0，1也是可以的，无所谓
similarity.cutoff=c(0.1,1)

#把最大最小值过滤一下，照着跑就行
J[J<similarity.cutoff[1]]<-similarity.cutoff[1]
J[J>similarity.cutoff[2]]<-similarity.cutoff[2]

showtree=TRUE

if(!showtree){
  tree<-FALSE
}

library(ComplexHeatmap)
pdf("heatmap.pdf", width = 13, height = 10)
ComplexHeatmap::pheatmap(J,scale='none',
                         annotation_names_col=FALSE,
                         annotation_names_row=FALSE,
                         color=c('white',viridis(100,option="A",direction=-1)),
                         show_rownames=F,
                         show_colnames=F,
                         cluster_cols=tree,
                         cluster_rows=tree,
                         annotation_col=annotation_col,
                         annotation_row=annotation_col)
dev.off()

J2=as.data.frame(J)
library(tidyr)
library(tidyverse)
# 将矩阵转换为长格式
  long_matrix <- J2 %>%
  rownames_to_column("Var1") %>%
  pivot_longer(cols = -Var1, 
               names_to = "Var2", 
               values_to = "value")

# 获取行和列的聚类顺序*
row_order <- tree$order
col_order <-tree$order

# 重新排序长矩阵*
long_matrix$Var1 <- factor(long_matrix$Var1, levels = rownames(J)[row_order])
long_matrix$Var2 <- factor(long_matrix$Var2, levels = colnames(J)[col_order])


library(ggplot2)
library(viridis)

ggplot(long_matrix, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() + # 使用热图图形*
  scale_fill_gradient(low = "white", high = viridis(100, option = "A", direction = -1)) + # 手动设置渐变色*
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(), axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))+
  xlab("Program") + ylab("Program") + # 设置轴标签*
  scale_y_discrete(limits = rev(levels(long_matrix$Var1)))


top_p<- lapply(geneNMF.metaprograms$metaprograms.genes,
               function(program){
                 runGSEA(genes=program,universe=rownames(scRNA_epi),category="C5",subcategory = 'BP') 
               }) #runGSEA()函数可用于链接msigDB数据库并评估检测到的MPs与数据库中signature的重叠情况*
  
  
  # 提取每个MP的前5个term并生成文本（仅显示MP标题和term）*
  top5_terms_text <- lapply(top_p, function(program_result) {
  # 按照 pval 排序并选择前5个路径*
      top5 <- program_result[order(program_result$pval), ][1:5, ]
      
  # 格式化文本，直接列出term*
        term_text <- top5$pathway
        
        return(term_text)
  })

# 使用gridExtra包中的textGrob来绘制文本*
library(gridExtra)
library(grid)

# 创建文本图形*
  plots <- lapply(names(top5_terms_text), function(mp_name) {
    # 每个MP的名称 + 对应的前5个term*
      term_text <- paste0(mp_name,': ',top5_terms_text[[mp_name]], collapse = "\n")
      # 将文本作为textGrob格式绘制*
        text_grob <- textGrob(term_text, gp = gpar(fontsize = 10))
        return(text_grob)
  })

dev.off()
# 绘制每个MP的前5个term文本*
grid.arrange(grobs = plots, ncol = 1)

# 提取每个MP的前5个term并生成文本（仅显示term，不显示MP名称）*
  top5_terms_text <- lapply(top_p, function(program_result) {
    # 按照 pval 排序并选择前5个路径*
      top5 <- program_result[order(program_result$pval), ][1:5, ]
      
      # 格式化文本，直接列出term（只包含pathway）*
        term_text <- top5$pathway
        
        return(term_text)
  })

# 使用gridExtra包中的textGrob来绘制文本*
library(gridExtra)
library(grid)

# 颜色设置：为每个MP分配不同的颜色*
  color_palette <- RColorBrewer::brewer.pal(length(top5_terms_text), "Set3") # 使用Set3颜色*
  
  # 创建文本图形*
  plots <- lapply(seq_along(top5_terms_text), function(i) {
    # 获取当前MP的路径*
      term_text <- paste(top5_terms_text[[i]], collapse = "\n")
      
      # 为每个MP分配颜色*
        term_color <- color_palette[i]
        
        # 使用 textGrob 绘制路径名，使用不同颜色*
          text_grob <- textGrob(term_text, gp = gpar(fontsize = 10, col = term_color))
          return(text_grob)
  })

# 绘制每个MP的前5个term文本，每个MP的颜色不同*
  grid.arrange(grobs = plots, ncol = 1)



library(dplyr)

# 将每个MP的富集结果合并为一个大的数据框*
  combined_results <- lapply(names(top_p), function(mp_name) {
    # 提取每个MP的富集结果*
      program_result <- top_p[[mp_name]]
      
      # 添加MP的名称列*
        program_result$MP <- mp_name
        
        # 返回带有MP名称的数据框*
          return(program_result)
  }) %>% bind_rows() # 将所有MP的结果合并*
  
  # 查看合并后的数据框*
  head(combined_results)



# 提取每个MP的前5个pathway*
  top5_combined_results <- lapply(names(top_p), function(mp_name) {
    # 提取每个MP的富集结果*
      program_result <- top_p[[mp_name]]
      
      # 按pval排序并选择前5个pathway*
        top5 <- program_result[order(program_result$pval), ][1:5, ]
        
        # 添加MP名称列*
          top5$MP <- mp_name
          
          return(top5)
  }) %>% bind_rows() # 合并所有MP的前5个pathway*
  
  # 查看前5个pathway的数据框*
  head(top5_combined_results)
  top5_combined_results$overlapGenes <- NULL
  top5_combined_results <- as.data.frame(top5_combined_results)
  write.csv(top5_combined_results, 'top5_MP_GSEA.csv')


library(ggplot2)

# 绘制条形图展示不同MP的前5个富集pathway*
 p1<- ggplot(top5_combined_results, aes(x = reorder(pathway, pval), y = -log10(pval), fill = MP)) +
  geom_bar(stat = "identity", show.legend = TRUE) + # 画条形图*
  coord_flip() + # 翻转坐标轴，横向展示*
  theme_minimal() + # 简约主题*
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # x轴标签旋转*
  labs(x = "Pathway", y = "-log10(p-value)", title = "Top 5 Enriched Pathways for Each MP") + # 添加标题和标签*
  scale_fill_brewer(palette = "Set3") + # 设置颜色*
  facet_wrap(~ MP, scales = "free_y") # 分面展示每个MP，y轴独立*
  
  #补充画个评分图
  
  mp.genes <- geneNMF.metaprograms$metaprograms.genes
seu2 <- AddModuleScore(pbmc,
                       features = mp.genes,
                       name = names(mp.genes))

mp.genes <- do.call("cbind", mp.genes)
mp.genes <- as.data.frame(mp.genes)
write.csv(mp.genes, "mp.genes.csv")


plist <- c()
for (i in 16:length(seu2@meta.data)){
plist[[i-15]] <- FeaturePlot(seu2, features =colnames(seu2@meta.data)[i]) +
   scale_color_gradientn(colors = RColorBrewer::brewer.pal(11, "Spectral"))
}

pdf("featureplot.pdf", width = 10, height = 10)
wrap_plots(plist)
dev.off()

pbmc=RenameIdents(pbmc, Type)
pbmc.markers=FindMarkers(pbmc, ident.1='UC', ident.2='Control', logfc.threshold=0.1)
sig.markers=pbmc.markers[(abs(as.numeric(as.vector(pbmc.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(pbmc.markers$p_val_adj))<adjPvalFilter),]

setwd('D:/216ST/GeneNMF/Enterocytes/UC')
mp.genes<-read.csv("mp.genes.csv")
mp.genes
sig.genes<-read.csv("05.Enterocytes.diffGene.csv")
sig.genes
