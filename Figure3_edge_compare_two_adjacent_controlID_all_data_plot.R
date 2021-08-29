# Ying Zhu, so_zy2003@126.com
# 2021-07

# clean environment variables

rm(list=ls())

# taxonomy for each ASV
taxonomy = read.delim("taxonomy.txt", row.names= 1,header=F, sep="\t")
colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species")
taxonomy <- taxonomy[-1,]

#read the files which is produced in script edge_compare_two_adjacent_controlID_all_data.R
d3_d1<- read.table("d3-d1-210706.txt", header = T, sep = "\t")
#filtering taxonomy
taxonomy_d3_d1 <- taxonomy[rownames(taxonomy)%in%rownames(d3_d1),]
d3_d1_tax = merge(d3_d1, taxonomy_d3_d1, by="row.names")
d3_d1_tax$compares <- "d3_d1"

d6_d3 <- read.table("d6-d3-210706.txt", header = T, sep = "\t")
taxonomy_d6_d3 <- taxonomy[rownames(taxonomy)%in%rownames(d6_d3),]
d6_d3_tax = merge(d6_d3, taxonomy_d6_d3, by="row.names")
d6_d3_tax$compares <- "d6_d3"


d24_d21 <- read.table("d24-d21-210706.txt", header = T, sep = "\t")
taxonomy_d24_d21 <- taxonomy[rownames(taxonomy)%in%rownames(d24_d21),]
d24_d21_tax = merge(d24_d21, taxonomy_d24_d21, by="row.names")
d24_d21_tax$compares <- "d24_d21"


d27_d24 <- read.table("d27-d24-210706.txt", header = T, sep = "\t")
taxonomy_d27_d24 <- taxonomy[rownames(taxonomy)%in%rownames(d27_d24),]
d27_d24_tax = merge(d27_d24, taxonomy_d27_d24, by="row.names")
d27_d24_tax$compares <- "d27_d24"


d30_d27 <- read.table("d30-d27-210706.txt", header = T, sep = "\t")
taxonomy_d30_d27 <- taxonomy[rownames(taxonomy)%in%rownames(d30_d27),]
d30_d27_tax = merge(d30_d27, taxonomy_d30_d27, by="row.names")
d30_d27_tax$compares <- "d30_d27"

 
d44_d37 <- read.table("d44-d37-210706.txt", header = T, sep = "\t")
taxonomy_d44_d37 <- taxonomy[rownames(taxonomy)%in%rownames(d44_d37),]
d44_d37_tax = merge(d44_d37, taxonomy_d44_d37, by="row.names")
d44_d37_tax$compares <- "d44_d37"


#combine all the comparison
compares_2points <- rbind(d3_d1_tax,  d6_d3_tax, 
                          d24_d21_tax, d27_d24_tax,d30_d27_tax,d44_d37_tax)

# add a column for comparison level
compares_2points $compares <- factor(compares_2points $compares,
                                     levels=c("d3_d1","d6_d3", 
                                              "d24_d21", "d27_d24","d30_d27","d44_d37"
                                     ))
# add a column for ASV ID 
compares_2points$ID <- compares_2points$Row.names

library(tidyr)
temp <-  separate(compares_2points,ID,into=c("ASV", "number"), sep="_")

num <- sort(unique(as.numeric(temp$number)))

compares_2points$ID <- factor(compares_2points$ID, level=paste("ASV_", num, sep=""))

colnames(compares_2points)[colnames(compares_2points)=="regulated"] <- "level"

compares_2points$level[compares_2points$level=="up"]<- "Enriched"
compares_2points$level[compares_2points$level=="down"]<- "Depleted"

# plotting
library(ggplot2)
p <- ggplot(compares_2points , aes(x=compares, y=ID, shape=level,
                                   size=logCPM, col=class))+ 
  geom_point()+theme_bw()+xlab("")+ylab("")+
  scale_shape_manual(values=c(25, 17))+
  theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))

edge_compare_p <- p
edge_compare_p

# ggsave(edge_compare_p, file="compare_two_adjacent_controlID.pdf",
#        height =18, width =15, units="cm")


