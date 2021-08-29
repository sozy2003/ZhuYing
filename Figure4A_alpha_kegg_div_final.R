# Ying Zhu, so_zy2003@126.com
# 2021-07


# clean environment variables

rm(list=ls()) 

# load related packages
library("ggplot2") 

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
design = read.table("metadata_weight.txt", header=T, row.names= 1, sep="\t") 

#aggregate a column for three stages
design$phase <- design$Day

design$phase[design$Day<12] <- "p1"  #stage one in Figure 1
design$phase[design$Day%in%c(12,15,18,21,24)] <- "p2"  #stage two in Figure 1
design$phase[design$Day>24] <- "p3" #stage three in Figure 1

#alpha diversity
alphadiv= read.table("vegan_picrust.txt", header=T, row.names= 1, sep="\t")


alphadiv = alphadiv[rownames(design), ]

alpha = cbind(alphadiv, design)

# plot the alpha diversity for total data 
#boxplot
p=ggplot(alpha,aes(x=as.factor(Day),y=shannon))+
  geom_boxplot(alpha=1, outlier.size=0, size=0.7, width=0.5, fill="transparent")+
  labs(x="Day", y="Shannon") + theme_classic()
p=p+geom_jitter(position=position_jitter(0.17), size=1, alpha=0.7)

p=p+theme(axis.text=element_text(angle=45,vjust=1, hjust=1))
p.shannon.box=p
# ggsave(paste("picrust_alpha_box_shannon",".pdf", sep=""), p.shannon.box, width = 4, height =2.5 )

#scatter point
p=ggplot(alpha)+
  geom_point(aes(x=Day,y=shannon),alpha=1)+
  geom_smooth(aes(x=Day,y=shannon),method = "lm", formula = y ~ poly(x,2))+
labs(x="Day", y="Shannon Index") + theme_classic()
# p=p+theme(axis.text=element_text(angle=45,vjust=1, hjust=1))
p=p+main_theme
p.shannon.point=p
 # ggsave(paste("picrust_alpha_point_shannon",".pdf", sep=""), p.shannon.point, width = 4, height =2.5 )

#scatter point ,richness
p=ggplot(alpha)+
  geom_point(aes(x=Day,y=richness),alpha=1)+
  geom_smooth(aes(x=Day,y=richness),method = "lm", formula = y ~ poly(x,3))
labs(x="Day", y="Shannon") + theme_classic()
p=p+theme(axis.text=element_text(angle=45,vjust=1, hjust=1))
p=p+main_theme
p.shannon.point=p
# ggsave(paste("picrust_alpha_point_richness",".pdf", sep=""), p.shannon.point, width = 4, height =2.5 )



#fit the model

#load the packages
if(!require("mgcv"))install.packages("mgcv")
library(mgcv)


#group by time points
library(dplyr)
alpha.mean <- alpha[,c(1:6,10)]%>%
  group_by(Day)%>%
  summarise_all(mean)
alpha.mean <- as.data.frame(alpha.mean)

#shannon
p=ggplot(alpha.mean,aes(x=Day,y=shannon))+
  labs(x="Day", y="shannon_index") + theme_classic()+
  geom_point(size=1,alpha=1)
p1=p+geom_smooth()
p2=p+geom_smooth(method = "lm", formula = y ~ poly(x,3))
p3=p+geom_smooth(method = "lm", formula = y ~ poly(x,2))
p3=p+geom_smooth(method=gam, formula = y~s(x))
p.mean.shannon.point=p2
# ggsave(paste("picrust_alpha_mean_point_shannon",".pdf", sep=""), p.mean.shannon.point, width = 4, height =2.5 )



#richness
p=ggplot(alpha.mean,aes(x=Day,y=richness))+
  labs(x="Day", y="richness") + theme_classic()+
  geom_point(size=1,alpha=1)
p1=p+geom_smooth()
p2=p+geom_smooth(method = "lm", formula = y ~ poly(x,3))
p3=p+geom_smooth(method=gam, formula = y~s(x))
p.mean.richness.point=p2
# ggsave(paste("picrust_alpha_mean_point_richness",".pdf", sep=""), p.mean.richness.point, width = 4, height =2.5 )


# check the relationship between signal independent variable 
#and dependent variable to choose possible variables in the final models
dat <- alpha
lm_day <- lm(shannon~Day, dat)#
summary(lm_day)# age significant

lm_stage <- lm(shannon~Stage, dat)
summary(lm_stage)#diet significant

lm_phase <- lm(shannon~phase, dat)
summary(lm_phase)#stage significant

lm_sex <- lm(shannon~Sex, dat)
summary(lm_sex)

lm_group <- lm(shannon~Group, dat)
summary(lm_group)#group non significant


lm_food <- lm(shannon~food, dat)
summary(lm_group)#food, non significant


lm_weight <- lm(shannon~weight, dat)
summary(lm_weight)#?

#age和stage, phase都显著的预测因子

# linear mixed-effects models

library(nlme)

# two stages
lme_1 <- lme(shannon~Day+Stage, 
             random=~1|Individual,method ="ML",data=dat)
summary(lme_1)

lme_2 <- lme(shannon~poly(Day,2)+Stage, 
             random=~1|Individual,method ="ML",data=dat)
summary(lme_2)


lme_3 <- lme(shannon~poly(Day,3)+Stage, 
             random=~1|Individual,method ="ML",data=dat)
summary(lme_3)

anova(lme_1,lme_2, lme_3)
# choose lme_2


# three stages
lme_1 <- lme(shannon~Day+phase, 
             random=~1|Individual,method ="ML",data=dat)
summary(lme_1)

lme_2 <- lme(shannon~poly(Day,2)+phase, 
             random=~1|Individual,method ="ML",data=dat)
summary(lme_2)


lme_3 <- lme(shannon~poly(Day,3)+phase, 
             random=~1|Individual,method ="ML",data=dat)
summary(lme_3)

anova(lme_1,lme_2, lme_3)
# choose lme_2 which has the smallest AIC as the final model 



# three stages+diet change
lme_1 <- lme(shannon~Day+phase+Stage, 
             random=~1|Individual,method ="ML",data=dat)
summary(lme_1)

lme_2 <- lme(shannon~poly(Day,2)+phase+Stage, 
             random=~1|Individual,method ="ML",data=dat)
summary(lme_2)


lme_3 <- lme(shannon~poly(Day,3)+phase+Stage, 
             random=~1|Individual,method ="ML",data=dat)
summary(lme_3)

anova(lme_1,lme_2, lme_3)
# choose lme_2 which has the smallest AIC as the final model 
