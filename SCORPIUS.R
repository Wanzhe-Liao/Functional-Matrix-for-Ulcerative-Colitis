#install.packages("SCORPIUS")
library(SCORPIUS)
expression_data <- t(as.matrix(pbmc_epi@assays$RNA@layers$data))
a <- pbmc_epi@assays$RNA@features@.Data
b <- pbmc_epi@assays$RNA@cells@.Data
colnames(expression_data) <- rownames(a)
rownames(expression_data) <- rownames(b)
expression_scale <- t(as.matrix(pbmc_epi@assays$RNA@layers$scale.data))
dim(expression_scale)
c <- rownames(a[which(a[,3] == "TRUE"),])
colnames(expression_scale) <- c
rownames(expression_scale) <- rownames(b)
group_name =  as.factor(as.character(pbmc_epi@active.ident))
table(group_name)
space <- reduce_dimensionality(expression_scale, "spearman")
draw_trajectory_plot(space, group_name, contour = TRUE)
traj <- infer_trajectory(space)
draw_trajectory_plot(space, group_name, traj$path, contour = TRUE)

gimp <- gene_importances(expression_scale, 
                         traj$time,   num_permutations = 24, #置换次数  
                         num_threads = 128, #线程数  
                         ntree = 10000,
                         ntree_perm = 1000) #每次随机重排后的树的数量 


gimp$qvalue <- p.adjust(gimp$pvalue, "BH", length(gimp$pvalue))
gene_sel <- gimp$gene[gimp$qvalue < .05]
expr_sel <- scale_quantile(expression_data[,gene_sel])
dim(expr_sel)
time <- traj$time
draw_trajectory_heatmap(expr_sel, time)
draw_trajectory_heatmap(expr_sel, time,                         
                        progression_group=group_name)









pbmc_epi <- subset(pbmc, subset=Type=="Treat")
