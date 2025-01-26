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


#正常组样本分布
x <- table(pbmc@active.ident, pbmc@meta.data$orig.ident)
colnames(x) <- c('CT','CT','CT','CT','CT','CT','CT','CT','UC','UC','UC','UC')
x3= t(t(x)/rowSums(t(x)))
x4 = as.data.frame(as.table(t(x3)))
colnames(x4) = c("sample","celltype","Freq")
library(stringr)
x4$group = str_remove(x4$sample,pattern = '[0-9]')


top<-function(x){
  return(mean(x)+sd(x)/sqrt(length(x)))
}
bottom<-function(x){
  return(mean(x)-sd(x)/sqrt(length(x)))
}
table(x4$group)
# Normal  Tumor 
# 21     21 

#### 根据个人更改！！！
dose_CT <-x4[which(x4$group=="CT"),]
dose_UC<-x4[which(x4$group=="UC"),]
############################################
ggplot(data=dose_CT,aes(x=celltype,y=Freq,fill=celltype))+
  stat_summary(geom = "bar",fun = "mean",
               position = position_dodge(0.9))+
  stat_summary(geom = "errorbar",
               fun.min = bottom,
               fun.max = top,
               position = position_dodge(0.9),
               width=0.2)+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="Celltype",y="Proportion")+
  geom_point(data=dose_sham,aes(celltype,Freq),size=3,pch=19)+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) +ggtitle("CT")

#肿瘤组样本分布
ggplot(data=dose_UC,aes(x=celltype,y=Freq,fill=celltype))+
  stat_summary(geom = "bar",fun = "mean",
               position = position_dodge(0.9))+
  stat_summary(geom = "errorbar",
               fun.min = bottom,
               fun.max = top,
               position = position_dodge(0.9),
               width=0.2)+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="Celltype",y="Proportion")+
  geom_point(data=dose_MCAO,aes(celltype,Freq),size=3,pch=19)+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) +ggtitle("UC")


分组堆叠柱状图
table(pbmc$orig.ident)#查看各组细胞数
# Normal1 Normal2 Normal3  Tumor1  Tumor2  Tumor3 
# 525      51     706    3747      21    2113

prop.table(table(Idents(pbmc)))

table(Idents(pbmc), pbmc$orig.ident)#各组不同细胞群细胞数


Cellratio <- prop.table(table(Idents(pbmc), pbmc$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio

Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
ggplot(data = Cellratio, aes(x =Var2, y = Freq, fill =  Var1)) +
  geom_bar(stat = "identity", width=0.8,position="fill")+
  scale_fill_manual(values=col_vector[1:20]) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Ratio")+
  ####用来将y轴移动位置
  theme(axis.text.y = element_text(size=12, colour = "black"))+
  theme(axis.text.x = element_text(size=12, colour = "black"))+
  theme(
    axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
  ) 


#差异统计图
table(pbmc$orig.ident)#查看各组细胞数
prop.table(table(Idents(pbmc)))
table(Idents(pbmc), pbmc$orig.ident)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(pbmc), pbmc$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- data.frame(Cellratio)
library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]


###添加分组信息
sample <- c('CT1','CT2','CT3','CT4','CT5','CT6','CT7','CT8','UC1','UC2','UC3','UC4')
group <- c('CT','CT','CT','CT','CT','CT','CT','CT','UC','UC','UC','UC')
samples <- data.frame(sample, group)#创建数据框
rownames(cellper)<-c('CT1','CT2','CT3','CT4','CT5','CT6','CT7','CT8','UC1','UC2','UC3','UC4')
rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']#R添加列
cellper$group <- samples[rownames(cellper),'group']#R添加列


#组间分析
pplist = list()
sce_groups = levels(pbmc@active.ident)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(cowplot)
for(group_ in sce_groups){
  cellper_  = cellper %>% select(one_of(c('sample','group',group_)))#选择一组数据
  colnames(cellper_) = c('sample','group','percent')#对选择数据列命名
  cellper_$percent = as.numeric(cellper_$percent)#数值型数据
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75), 
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))#上下分位数
  print(group_)
  print(cellper_$median)
  
  pp1 = ggplot(cellper_,aes(x=group,y=percent)) + #ggplot作图
    geom_jitter(shape = 21,aes(fill=group),width = 0.25) + 
    stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') + 
    labs(title = group_,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)
  
  ###组间t检验分析
  labely = max(cellper_$percent)
  compare_means(percent ~ group,  data = cellper_)
  my_comparisons <- list( c("CT", "UC") )
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}

library(cowplot)
#patch <- c()
#for (i in levels(df[,"Cluster"])){
#patch <- patch + plot_grid(pplist[[i]])
#}
wrap_plots(pplist)
