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

#  Design of experiment
design = read.table("metadata.txt", header=T, row.names= 1, sep="\t") 

# raw reads count of each ko in each sample

ko3 <-read.delim("ko3.txt",  header=T, sep="\t", comment.char = "#")

library(tidyr)
ko3 <-  separate(ko3,KEGG_Pathways,into=c("L1", "L2", "L3"), sep=";")
rownames(ko3) <-ko3$L3 


kegg <- ko3[,-1]


kegg <- kegg[,-(152:154)]
tax_count_sum <- kegg


per = t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 # normalization to total 100


# plot for each age group
mean_sort = per[(order(-rowSums(per))), ] # decrease sort
colSums(mean_sort)

# the top 10
mean_sort=as.data.frame(mean_sort)
other = colSums(mean_sort[40:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(40-1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[40] = c("Low Abundance")
# write.table(mean_sort, file="Top10ko3.txt", append = F, sep="\t", quote=F, row.names=T, col.names=T)
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
mat_mean_final$L3 = rownames(mat_mean_final)


data_all = as.data.frame(melt(mat_mean_final, id.vars=c("L3")))


p = ggplot(data_all, aes(x=variable, y = value, fill = L3 )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_y_continuous(labels = scales::percent) + 
  labs(x="Day", y="Relative Abundance (%)")+main_theme
tax.stack=p

# alluvium
p = ggplot(data = data_all, aes(x = variable, y = value, alluvium = L3, stratum=L3)) +
  geom_alluvium(aes(fill = L3), alpha = 0.75) +
  geom_stratum(aes(fill=L3))+
  labs(x="Day", y="Relative Abundance (%)")+
  main_theme 

tax.alluvium=p

tax.alluvium

# ggsave("tax_alluvium_ko3_top40_group.pdf", tax.alluvium, width = 12, height = 5)

