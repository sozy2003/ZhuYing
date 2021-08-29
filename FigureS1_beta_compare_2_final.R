# Ying Zhu, so_zy2003@126.com
# 2021-07

rm(list=ls()) # clean enviroment object

# Install related packages
if (FALSE){
  source("https://bioconductor.org/biocLite.R")
  biocLite(c("ggplot2","grid","scales","vegan","agricolae","ggrepel","dplyr"))
  install.packages("devtools", repo="http://cran.us.r-project.org")
  library(devtools)
  install_github("vqv/ggbiplot")
}



library("ggplot2") # load related packages
library("grid")
library("scales")
library("vegan")
library("agricolae")
library("ggbiplot")
library("dplyr")
library("ggrepel")

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


otu_table = read.delim("otu_table.txt", row.names= 1,  header=T, sep="\t")

calc_dis = function(subset){
  
 
  idx = grepl(subset, design$Site)
  sub_design = design[idx,]
idx=rownames(sub_design) %in% colnames(otu_table)
  sub_design=sub_design[idx,]
  count = otu_table[, rownames(sub_design)]
  norm = t(t(count)/colSums(count,na=T)) * 100 # normalization to total 100
  norm=as.data.frame(norm)
  
 # day 1 as an example
  # change the number for the interesting time point
  final=rownames(sub_design[sub_design$Day %in% 1,])
  ck=norm[,final] 
  ck_mean= rowMeans(ck)
  norm$ck_mean=ck_mean 
  
 
  bray_curtis = vegdist(t(norm), method = "bray")
  bray_curtis= as.matrix(bray_curtis)
  

  dat = t(as.data.frame(c("sampleA","sampleB","0.15","group","genosite")))
  colnames(dat) = c("sampleA","sampleB","distance","group","type")
  rownames(dat) = c("test")
  
  
  for (i in sort(unique(sub_design$Day))){
    group = rownames(sub_design[sub_design$Day %in% i,])
    for (m in 1:length(group)) {
      x = c(group[m],"ck_mean",bray_curtis[group[m],"ck_mean"],i,subset)
      dat=rbind(dat,x)
    }
  }
  dat = as.data.frame(dat[-1,], stringsAsFactors=F) # 删除首行框架数据
  # dat = dat[dat$distance != 0,] # 
 
  dat$distance=round(as.numeric(as.character(dat$distance)), digits=3)

  dat$group = factor(dat$group, levels=unique(dat$group))
  return(dat)
}


dat = calc_dis("Emei")

all = dat



#day 1 
all_day1 <- all
all_day1$control <- 1


# we cant get data for other time points as day 1 

#combine the data
all_beta_compare <- rbind(all_day1,all_day3, all_day6,all_day9,
                          all_day12,all_day15,all_day18,all_day21,
                          all_day24,all_day27,all_day30,all_day37,
                          all_day44)

all_beta_compare$within <- ifelse(all_beta_compare$group==all_beta_compare$control, "within","between")
all_beta <- all_beta_compare

all_beta_compare$within <- factor(all_beta_compare$within, level = c("within","between"))
all_beta_compare$control <- factor(all_beta_compare$control, 
                                       level=unique(sort(all_beta_compare$group)))
all_beta_compare$group <- factor(all_beta_compare$group, 
                                     level=unique(sort(all_beta_compare$group)))

# boxplot
p = ggplot(all_beta_compare, aes(x=group, y=distance, 
                                     group=group, col=within)) +
  geom_boxplot(alpha=1, outlier.size=0, size=0.5, width=0.5) +
  labs(x="Age", y="Bray-Curtis distance") +
  scale_color_manual( values=c("red", "black"))+
  theme_classic()+
  geom_jitter(position=position_jitter(0.17), size=0.3, alpha=0.5)+
  facet_wrap(~control)
p
all_beta_compare_p <- p

# ggsave(paste("boxplot_compare_all_day",".pdf", sep=""), 
#        all_beta_compare_p, width = 10, height =8 )







