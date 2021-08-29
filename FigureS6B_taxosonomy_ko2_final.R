# Ying Zhu, so_zy2003@126.com
# 2021-07



rm(list=ls()) # clean enviroment object
library("reshape2", quietly=T, warn.conflicts=F)
library(ggalluvial)

# Set ggplot2 drawing parameter
main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=7),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=7),
                   text=element_text(family="sans", size=7))

# Design of experiment
design = read.table("metadata.txt", header=T, row.names= 1, sep="\t") 

#  raw reads count of each ko in each sample

ko2 <-read.delim("ko2.txt",  header=T, sep="\t", comment.char = "#")

library(tidyr)
ko2 <-  separate(ko2,KEGG_Pathways,into=c("L1", "L2"), sep=";")
rownames(ko2) <-ko2$L2 


kegg <- ko2[,-1]

#numeric matrix
kegg <- kegg[,-(152:153)]
tax_count_sum <- kegg

# normalization
per = t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 # normalization to total 100


# plot for each age group
mean_sort = per[(order(-rowSums(per))), ] # decrease sort
colSums(mean_sort)#验证百分比总和为1

# the top 10 ko
mean_sort=as.data.frame(mean_sort)
other = colSums(mean_sort[42:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(42-1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[42] = c("Low Abundance")
topN=rownames(mean_sort)


mean_sort <- mean_sort[,rownames(design)]


mat=mean_sort

mat_t = t(mat)

mat_t2 = merge(design[c("Day")], mat_t, by="row.names")

mat_t2 = mat_t2[,-1]

mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]

geno = mat_mean$Day
colnames(mat_mean_final) = geno

mat_mean_final = as.data.frame(mat_mean_final)

mat_mean_final$L2 = rownames(mat_mean_final)


data_all = as.data.frame(melt(mat_mean_final, id.vars=c("L2")))

data_all <- na.omit(data_all)
p = ggplot(data_all, aes(x=variable, y = value, fill = L2 )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_y_continuous(labels = scales::percent) + 
  labs(x="Day", y="Relative Abundance (%)")+main_theme
tax.stack=p

# 绘制冲击图alluvium
p = ggplot(data = data_all, aes(x = variable, y = value, 
                                alluvium = L2, stratum=L2)) +
  geom_alluvium(aes(fill = L2), alpha = 0.75) +
  geom_stratum(aes(fill=L2))+
  labs(x="Day", y="Relative Abundance (%)")+
  main_theme 

tax.alluvium=p

tax.alluvium

# ggsave("tax_alluvium_ko2_all41_group.pdf", tax.alluvium, width = 12, height = 5)

