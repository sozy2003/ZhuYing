# Ying Zhu, so_zy2003@126.com
# 2021-07

# clean environment variables
rm(list=ls()) 

# load related packages
if(!require(ggplot2))install.packages("ggplot2")
library("ggplot2") 

# Set ggplot2 drawing parameter
main_theme <-  theme(panel.background=element_blank(),
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

#   Design of experiment
design <-  read.table("metadata.txt", header=T, row.names= 1, sep="\t") 

#alpha diversity
alphadiv <-  read.table("vegan.txt", header=T, row.names= 1, sep="\t")


#reorder the rows
alphadiv <-  alphadiv[rownames(design), ]

#conbine the two data frames
alpha = cbind(alphadiv, design)

# plot the alpha diversity for total data 
#boxplot to see the pattern, shannon index
p=ggplot(alpha,aes(x=as.factor(Day),y=shannon))+
           geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent")+
           labs(x="Day", y="Shannon") + theme_classic()
p=p+geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7)
              
p=p+main_theme
p.shannon.box=p
# ggsave(paste("alpha_box_shannon",".pdf", sep=""), p.shannon.box, width = 4, height =2.5 )

#scatter point for shannon index
p=ggplot(alpha)+
  geom_point(aes(x=Day,y=shannon),alpha=1)+
  geom_smooth(aes(x=Day,y=shannon),method = "lm", formula = y ~ poly(x,3))+
  labs(x="Day", y="Shannon Index") 

p=p+main_theme
p.shannon.point=p
p.shannon.point
# ggsave(paste("alpha_point_shannon",".pdf", sep=""), p.shannon.point, width = 4, height =2.5 )

#--------------以下是尝试--------------------
p=ggplot(alpha, aes(x=Day,y=shannon, color=Individual))+
  geom_point(alpha=1)+
  geom_smooth(method = "lm", formula = y ~ poly(x,3), se=F)+
  labs(x="Day", y="Shannon Index") 
p

# individual
p=ggplot(alpha, aes(x=Day,y=shannon, color=Individual))+
  geom_point(alpha=1)+
  geom_smooth(method = "lm", formula = y ~ x, se=F)+
  labs(x="Day", y="Shannon Index") 
p
# ggsave(paste("shannon_all_individuals_linear",".pdf", sep=""),
#        p, width = 5, height =4 )

# individual 
p=ggplot(alpha, aes(x=Day,y=shannon, color=Individual))+
  geom_point(alpha=1)+
  geom_smooth(method = "lm", formula = y ~ poly(x,3), se=F)+
  labs(x="Day", y="Shannon Index") 
p
# ggsave(paste("shannon_all_individuals_polynomial_linear",".pdf",
#              sep=""), p, width = 5, height =4 )

# individual facet
p=ggplot(alpha, aes(x=Day,y=shannon))+
  geom_point(alpha=1)+
geom_smooth(method = "lm", formula = y ~ poly(x,3), se=F)+
  labs(x="Day", y="Shannon Index")+
  facet_wrap(~Individual, scales="free", nrow=4)
p
# ggsave(paste("shannon_all_individuals_facet",".pdf", sep=""),
#        p, width = 8, height =6 )



#scatter point, richness
p=ggplot(alpha)+
  geom_point(aes(x=Day,y=richness),alpha=1)+
  geom_smooth(aes(x=Day,y=richness),method = "lm", formula = y ~ poly(x,3))+
  labs(x="Day", y="Richness Index") 
# p=p+theme(axis.text=element_text(angle=45,vjust=1, hjust=1))
p=p+main_theme
p.richness.point=p
# ggsave(paste("alpha_point_richenss",".pdf", sep=""), p.richness.point, width = 4, height =2.5 )


# fit the model
if(!require("mgcv"))install.packages("mgcv")
library(mgcv)
model.gam <- gam(shannon~s(Day),data=alpha)
summary(model.gam) 


# use the mean value for each day
library(dplyr)
alpha.mean <- alpha[,c(1:6,10)]%>%
  group_by(Day)%>%
  summarise_all(mean)
alpha.mean <- as.data.frame(alpha.mean)

#plot based on mean values 
p=ggplot(alpha.mean,aes(x=Day,y=shannon))+
  labs(x="Day", y="mean value of shannon index") + theme_classic()+
  geom_point(size=1,alpha=1)+main_theme
p1=p+geom_smooth(method = "lm", formula = y ~ x)# simple linear model
p2=p+geom_smooth(method = "lm", formula = y ~ poly(x,3))# polynomial model
p3=p+geom_smooth(method=gam, formula = y~s(x))#generalized additive models
p1
p2
p3
# ggsave(paste("alpha_mean_point_shannon_linear",".pdf", sep=""), 
#        p1, width = 4, height =2.5 )
# 
# ggsave(paste("alpha_mean_point_shannon",".pdf", sep=""), 
#        p2, width = 4, height =2.5 )
# 

#plot richness
p=ggplot(alpha.mean,aes(x=Day,y=richness))+
  labs(x="Day", y="richness") + theme_classic()+
  geom_point(size=1,alpha=1)+main_theme
p1=p+geom_smooth()
p2=p+geom_smooth(method = "lm", formula = y ~ poly(x,3))
p3=p+geom_smooth(method=gam, formula = y~s(x))
p.mean.richness.point=p2

p1
p2
p3
# ggsave(paste("alpha_mean_point_richness",".pdf", sep=""), p.mean.richness.point, width = 4, height =2.5 )
