# Ying Zhu, so_zy2003@126.com
# 2021-07



# clean enviroment object
rm(list=ls()) 

# load related packages
library("ggplot2") 
library("vegan")

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


design$phase <- design$Day

design$phase[design$Day<12] <- "p1"
design$phase[design$Day%in%c(12,15,18,21,24)] <- "p2"
design$phase[design$Day>24] <- "p3"

design$Day2 <- factor(design$Day2, 
                      levels = c("d1", "d3","d6", "d9", "d12", "d15", 
                                 "d18", "d21", "d24","d27", "d30", "d37", "d44"))

#  PCoA bray_curtis
distance = read.table("unifrac.txt", sep="\t", header=T, check.names=F)


rownames(distance)=distance[,1]
distance=distance[,-1]
distance = distance[rownames(design), rownames(design)] 

# subset and reorder distance matrix

# Classical multidimensional scaling (MDS) of a data matrix. Also known as principal coordinates analysis
pcoa = cmdscale(distance, k=4, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points) # get coordinate string, format to dataframme
colnames(points) = c("x", "y", "z","a") 
eig = pcoa$eig


points = cbind(points, design[match(rownames(points), rownames(design)), ])


# age effect
# plot PCo 1 and 2
p = ggplot(points, aes(x=x, y=y, color=Day))
p = p + geom_point(alpha=.7, size=2) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title=" Weighted Unifrac PCoA") + main_theme
q= p + scale_color_gradientn(colours=rainbow(7))
q
beta_pcoa_day_Weighted_unifrac1v2 <- q
# ggsave("beta_pcoa_day_Weighted_unifrac1v2_20210607.pdf", beta_pcoa_day_Weighted_unifrac1v2, width = 4, height = 2.5)


# plot PCo 1 and 3

p = ggplot(points, aes(x=x, y=z, color=Day,shape=Stage))
p = p + geom_point(alpha=.7, size=2) +
  scale_color_gradientn(colours=rainbow(7)) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 3 (", format(100 * eig[3] / sum(eig), digits=4), "%)", sep=""),
       title=" PCoA") + main_theme
p


#three stages
p = ggplot(points, aes(x=x, y=y, color=phase))
p = p + geom_point(alpha=.7, size=2) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title="Weighted Unifrac PCoA") + main_theme+stat_ellipse(level=0.95)
p
beta_pcoa_phase_Weighted_unifrac <- p
# ggsave("beta_pcoa_phase_Weighted_unifrac1v2_20210607.pdf", beta_pcoa_phase_Weighted_unifrac, width = 4, height = 2.5)



# sex effect
p = ggplot(points, aes(x=x, y=y, color=Sex))
p = p + geom_point(alpha=.7, size=2) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title="Bray_Curtis PCoA") + main_theme+stat_ellipse(level=0.95)
p

# family group
p = ggplot(points, aes(x=x, y=y, color=Group))
p = p + geom_point(alpha=.7, size=2) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title="Bray_Curtis PCoA") + main_theme+stat_ellipse(level=0.95)
p

# individual ID
p = ggplot(points, aes(x=x, y=y, color=Individual))
p = p + geom_point(alpha=.7, size=2) +
  labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title="Bray_Curtis PCoA") + main_theme+stat_ellipse(level=0.95)
p



#Permutational Multivariate Analysis of Variance Using Distance Matrices (PERMANOVA)

#permanova  置换多元方差分析
# Permutational Multivariate Analysis of Variance Using Distance Matrices (PERMANOVA)
adonis(distance~design$Day2, permutations = 999)#significant
#R2=0.272, P=0.001

adonis(distance~design$Group, permutations = 999)

adonis(distance~design$phase, permutations = 999)#significant

adonis(distance~design$Individual, permutations = 999)
adonis(distance~design$Sex, permutations = 999)

# the final model
adonis_day2 <- adonis(distance~design$phase+design$Day2,
                      permutations = 999)

adonis_day2$aov.tab  


