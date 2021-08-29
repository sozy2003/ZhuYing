# Ying Zhu, so_zy2003@126.com
# 2021-07


# clean environment variables
rm(list=ls()) 

# load related packages
library("ggplot2") 

# Set ggplot2 drawing parameter, such as axis line and text size, lengend and title size, and so on.
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


betadiv= read.table("picrustbray_curtis.txt", header=T, row.names= 1, sep="\t")


betadiv = betadiv[rownames(design), rownames(design)]

betadis=as.matrix(betadiv)

betadis[upper.tri(betadis,diag = T)] <- NA

day <- unique(sort(design$Day))


p=list() 
for (i in day){
  p[[i]] <- betadis[design$Day==i, design$Day==i]
  nrow <- dim(p[[i]])[1]
  ncol <- dim(p[[i]])[2]
  row <- rep(rownames(p[[i]]), ncol)
  col <- rep(colnames(p[[i]]),each=nrow)
  p[[i]] <- data.frame(row, col, value=as.numeric(p[[i]]))
  p[[i]]$day <- rep(i,nrow(p[[i]]))
  
}
p

#combine all the data

frame <- rbind(p[[1]],p[[3]],p[[6]],p[[9]],p[[12]],p[[15]],
               p[[18]],p[[21]], p[[24]], p[[27]],p[[30]],p[[37]],
               p[[44]])


frame <- na.omit(frame)

#boxplot
p=ggplot(frame,aes(x=as.factor(day),y=value))+
  geom_boxplot()+
  labs(x="Day", y="bray_curtis") + theme_classic()
p=p+geom_jitter(position=position_jitter(0.17), size=0.001, alpha=0.3,color='blue')+
  geom_smooth(method = "lm", formula = y ~ poly(x,3))
p=p+theme(axis.text=element_text(angle=45,vjust=1, hjust=1))
p.boxplot=p
# ggsave(paste("kegg_beta_box_bray_curtis",".pdf", sep=""), p.boxplot, width = 4, height =2.5 )



#scatter point
p=ggplot(frame,aes(x=day,y=value))+
  labs(x="Day", y="Bray-Curtis Distance") + theme_classic()+
  geom_point(size=1,alpha=1)
p=p+geom_smooth(method = "lm", formula = y ~ poly(x,3))

# p=p+theme(axis.text=element_text(angle=45,vjust=1, hjust=1))
p=p+main_theme
p.point=p
 # ggsave(paste("kegg_beta_scatter_bray_curtis",".pdf", sep=""), p.point, width = 4, height =2.5 )



#GAM fitting
if(!require("mgcv"))install.packages("mgcv")
library(mgcv)

#mean value
library(dplyr)
temp <- frame%>%
  group_by(day)%>%
  summarize(mean(value))

temp <- as.data.frame(temp)
colnames(temp) <- c("day", "mean")

#scatter point
p=ggplot(temp,aes(x=day,y=mean))+
  labs(x="Day", y="bray_curtis_mean") + theme_classic()+
  geom_point(size=1,alpha=1)
p3=p+geom_smooth(method = "lm", formula = y ~ poly(x,3))+main_theme
p2=p+geom_smooth(method = "lm", formula = y ~ poly(x,2))+main_theme
p1=p+geom_smooth(method = "lm", formula = y ~ x)+main_theme
p1
p2
p3
p.mean.point=p3
# ggsave(paste("kegg_beta_mean_scatter_bray_curtis",".pdf", sep=""), p.mean.point, width = 4, height =2.5 )

#----------------fit the model
# linear model
frame$stage <- frame$day
frame$stage[frame$stage<=24] <- "one"
frame$stage[frame$stage!="one"] <- "two"

lm1 <- lm(value~day+stage, frame)
summary(lm1) #R2=0.01

lm2 <- lm(value~poly(day,2)+stage, frame)
summary(lm2) #R2=0.01033

lm3 <- lm(value~poly(day,3)+stage, frame)
summary(lm3) #R2=0.068 stagetwo

anova(lm1,lm2,lm3)
#选择lm3
# coef_beta_kegg_age<- summary(lm3)$coefficient
# write.table(coef_beta_kegg_age, file="coef_beta_kegg_age.txt", sep = '\t')

#--------------final model----linear mixed-effects models
if(!require(nlme))install.packages("nlme")
library(nlme)

head(frame)
# add a column for individual ID

library(tidyr)

frame$temp <- frame$col
frame1 <- separate(frame,temp, sep="D", into=c("Individual", "day2"))
frame1 <- frame1[,colnames(frame1)!="day2"]

colnames(frame1) <- c("row", "col", "beta_distance", "Day", "Stage","Individual" )


# two stages
lme_1 <- lme(beta_distance~Day+Stage, 
             random=~1|Individual,method ="ML",data=frame1)
summary(lme_1)#stage, day


lme_2 <- lme(beta_distance~poly(Day,2)+Stage, 
             random=~1|Individual,method ="ML",data=frame1)
summary(lme_2)#stage


lme_3 <- lme(beta_distance~poly(Day,3)+Stage, 
             random=~1|Individual,method ="ML",data=frame1)
summary(lme_3)#stage

anova(lme_1,lme_2,lme_3)
# choose lme_3



# three stages
frame1$phase[frame1$Day<12] <- "p1"  #stage one in Figure 1
frame1$phase[frame1$Day%in%c(12,15,18,21,24)] <- "p2"  #stage two in Figure 1
frame1$phase[frame1$Day>24] <- "p3" #stage three in Figure 1

lme_1 <- lme(beta_distance~Day+phase, 
             random=~1|Individual,method ="ML",data=frame1)
summary(lme_1)


lme_2 <- lme(beta_distance~poly(Day,2)+phase, 
             random=~1|Individual,method ="ML",data=frame1)
summary(lme_2)# stage 2


lme_3 <- lme(beta_distance~poly(Day,3)+phase, 
             random=~1|Individual,method ="ML",data=frame1)
summary(lme_3)#stage 3 

anova(lme_1,lme_2,lme_3)
# choose lme_3



#three stages+diet changes

if(!require(lme4))install.packages(lme4)
if(!require(lmerTest))install.packages("lmerTest")
library(lme4)
library(lmerTest)
library(car)

lme_1 <- lmer(beta_distance~Day+phase+Stage+(1|Individual),data=frame1)
summary(lme_1)
vif(lme_1)

lme_2 <- lmer(beta_distance~poly(Day,2)+phase+Stage+(1|Individual), 
              data=frame1)
summary(lme_2)


lme_3 <- lmer(beta_distance~poly(Day,3)+phase+Stage+(1|Individual), 
              data=frame1)
summary(lme_3)#stage 3 

vif(lme_3)

anova(lme_1,lme_2,lme_3)
# choose lme_3
