setwd('D:/216ST/GeneNMF/Schwann Cells/UC')

a <- read.csv("05.Goblet Cells.diffGene.csv")
a
mp.genes <- read.csv("mp.genes.csv")
head(mp.genes)
mp.genes$X <- NULL

# 保留原始结构，交集后用 NA 填充
filtered_mp.genes <- as.data.frame(lapply(mp.genes, function(column) {
  result <- intersect(column, a$Gene)
  # 用 NA 填充到与原列长度一致
  c(result, rep(NA, length(column) - length(result)))
}))
filtered_mp.genes <- filtered_mp.genes[apply(filtered_mp.genes, 1, function(y) any(!is.na(y))),]
# 查看结果
print(filtered_mp.genes)
write.csv(filtered_mp.genes,"filtered_mp_genes.csv")

