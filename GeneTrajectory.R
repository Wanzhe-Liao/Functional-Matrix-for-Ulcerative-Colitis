library(GeneTrajectory)
require(Seurat)
require(scales)
require(ggplot2)
require(viridis)
require(dplyr)
require(GeneTrajectory)
require(Matrix)
require(plot3D)


assay <- "RNA"
DefaultAssay(pbmc) <- assay
pbmc <- FindVariableFeatures(pbmc, nfeatures = 500)

all_genes <- na.omit(pbmc@assays$RNA@meta.data$var.features)[1:10]
expr_percent <- apply(as.matrix(pbmc@assays$RNA$data[all_genes, ]) > 0, 1, sum)/ncol(pbmc)
genes <- all_genes[which(expr_percent > 0.01 & expr_percent < 0.5)]
length(genes)

# Compute the Diffusion Map cell embedding
pbmc <- GeneTrajectory::RunDM(object=pbmc,
  reduction = "pca",
  dims = 1:20,
  K = 10,
  sigma = NULL,
  n.components = 30,
  t = 1,
  dist.mat = NULL,
  reduction.key = "DM_"
)
# Calculate cell-cell graph distances over a cell-cell kNN graph
cell.graph.dist <- GetGraphDistance(pbmc, K = 10)
# Coarse-grain the cell graph by grouping cells into `N`=500 "meta-cells"
cg_output <- CoarseGrain(pbmc, cell.graph.dist, genes, N = 500)

# Create a virtualenv using reticulate
if(!reticulate::virtualenv_exists('gene_trajectory')){
  reticulate::virtualenv_create('gene_trajectory', packages=c('gene_trajectory'))
}
reticulate::use_virtualenv('gene_trajectory')
# Import the function to compute gene-gene distances
cal_ot_mat_from_numpy <- reticulate::import('gene_trajectory.compute_gene_distance_cmd')$cal_ot_mat_from_numpy
# Compute gene-gene distances 
gene.dist.mat <- cal_ot_mat_from_numpy(ot_cost = cg_output[["graph.dist"]], gene_expr = cg_output[["gene.expression"]], num_iter_max = 50000, show_progress_bar = TRUE)
rownames(gene.dist.mat) <- cg_output[["features"]]
colnames(gene.dist.mat) <- cg_output[["features"]]

gene_embedding <- GetGeneEmbedding(gene.dist.mat, K = 5)$diffu.emb

# Extract 3 gene trajectories
gene_trajectory <- ExtractGeneTrajectory(gene_embedding, gene.dist.mat, N = 3, t.list = c(4,7,7), K = 5)
table(gene_trajectory$selected)
# Visualize gene trajectories
par(mar = c(1.5,1.5,1.5,1.5))
scatter3D(gene_embedding[,1],
          gene_embedding[,2],
          gene_embedding[,3],
          bty = "b2", colvar = as.integer(as.factor(gene_trajectory$selected))-1,
          main = "trajectory", pch = 19, cex = 1, theta = 45, phi = 0,
          col = ramp.col(c(hue_pal()(3))))

# Extract the ordered list of genes along each gene trajectory
gene_list <- list()
for (i in 1:3){
  gene_trajectory_sub <- gene_trajectory[which(gene_trajectory$selected == paste0("Trajectory-", i)),]
  genes <- rownames(gene_trajectory_sub)[order(gene_trajectory_sub[, paste0("Pseudoorder", i)])]
  gene_list[[i]] <- genes
}

library(SeuratWrappers)
pbmc <- RunALRA(pbmc)

dim(gene.dist.mat)

pbmc <- AddGeneBinScore(pbmc, gene_trajectory, N.bin = 5, trajectories = 1:3, assay = "alra", reverse = c(F, F, T))

# Visualize gene bin plots for each gene trajectory
FeaturePlot(pbmc, pt.size = 0.05, features = paste0("Trajectory",1,"_genes", 1:5), ncol = 5, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes() & theme(title = element_text(size = 10))
pbmc <- AddGeneBinScore(pbmc, gene_trajectory, N.bin = 5, trajectories = 1:3, assay = "alra", reverse = c(F, F, T))
FeaturePlot(pbmc, pt.size = 0.05, features = paste0("Trajectory",2,"_genes", 1:5), ncol = 5, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes() & theme(title = element_text(size = 10))
FeaturePlot(pbmc, pt.size = 0.05, features = paste0("Trajectory",3,"_genes", 1:5), ncol = 5, order = T) &
  scale_color_gradientn(colors = rev(brewer_pal(palette = "RdYlBu")(10))) & NoLegend() & NoAxes() & theme(title = element_text(size = 10))