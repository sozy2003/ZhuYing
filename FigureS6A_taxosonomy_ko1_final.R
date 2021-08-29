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

ko1 <-read.delim("ko1.txt",  header=T, sep="\t", comment.char = "#")

rownames(ko1) <-ko1$KEGG_Pathways


kegg <- ko1[,-1]

# generate pure numeric matrix
kegg <- kegg[,-(152:153)]
tax_count_sum <- kegg


per = t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 # normalization to total 100



mean_sort = per[(order(-rowSums(per))), ] # decrease sort
colSums(mean_sort)#verify

#filter the top ko
mean_sort=as.data.frame(mean_sort)
topN=rownames(mean_sort)


mean_sort <- mean_sort[,rownames(design)]

# plot for each age group

mat=mean_sort

mat_t = t(mat)

mat_t2 = merge(design[c("Day")], mat_t, by="row.names")

mat_t2 = mat_t2[,-1]

mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]

geno = mat_mean$Day
colnames(mat_mean_final) = geno

mat_mean_final = as.data.frame(mat_mean_final)

mat_mean_final$L1 = rownames(mat_mean_final)


data_all = as.data.frame(melt(mat_mean_final, id.vars=c("L1")))

data_all <- na.omit(data_all)
p = ggplot(data_all, aes(x=variable, y = value, fill = L1 )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_y_continuous(labels = scales::percent) + 
  labs(x="Day", y="Relative Abundance (%)")+main_theme
tax.stack=p

# alluvium
p = ggplot(data = data_all, aes(x = variable, y = value, 
                                alluvium = L1, stratum=L1)) +
  geom_alluvium(aes(fill = L1), alpha = 0.75) +
  geom_stratum(aes(fill=L1))+
  labs(x="Day", y="Relative Abundance (%)")+
  main_theme 

tax.alluvium=p

tax.alluvium

# ggsave("tax_alluvium_ko1_all8_group.pdf", tax.alluvium, width = 12, height = 5)



