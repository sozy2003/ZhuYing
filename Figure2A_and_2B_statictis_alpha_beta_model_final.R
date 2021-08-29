# Ying Zhu, so_zy2003@126.com
# 2021-07

# clean environment variables
rm(list=ls()) # clean enviroment object


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

# Design of experiment
design = read.table("metadata_weight.txt", header=T, row.names= 1, sep="\t") 

#aggregate a column for three stages
design$phase <- design$Day

design$phase[design$Day<12] <- "p1"  #stage one in Figure 1
design$phase[design$Day%in%c(12,15,18,21,24)] <- "p2"  #stage two in Figure 1
design$phase[design$Day>24] <- "p3" #stage three in Figure 1

# alpha diversity
vegan = read.table("vegan.txt", header=T, row.names= 1, sep="\t") 
vegan <- vegan[rownames(design),]

rownames(design)==rownames(vegan)

# generate data frame for fitting
dat <- data.frame(vegan, group=design$Group, day=design$Day,day2=design$Day2, 
                  sex=design$Sex, stage=design$Stage, individual=design$Individual,
                  weight=design$weight, food=design$food, phase=design$phase )


#----------linear model-------------
# age
lm_day <- lm(shannon~day, dat)#
summary(lm_day)# significant

#two stages
lm_stage <- lm(shannon~stage, dat)
summary(lm_stage)#stage two # significant

#three stages
lm_phase <- lm(shannon~phase, dat)
summary(lm_phase)#stage three  # significant

# sex
lm_sex <- lm(shannon~sex, dat)
summary(lm_sex)

# family group
lm_group <- lm(shannon~group, dat)
summary(lm_group)

# food
lm_food <- lm(shannon~food, dat)
summary(lm_food)

#weight
lm_weight <- lm(shannon~weight, dat)
summary(lm_weight)# significant

# age and stage/ phase are possible predictors. 


#--------linear mixed effects models
library(nlme)

# two stages
lme_1 <- lme(shannon~day+stage, 
             random=~1|individual,method ="ML",data=dat)
summary(lme_1) 

lme_2 <- lme(shannon~day+I(day^2)+stage, 
             random=~1|individual,method ="ML",data=dat)
summary(lme_2) 

lme_3 <- lme(shannon~day+I(day^2)+I(day^3)+stage, 
             random=~1|individual,method ="ML",data=dat)
summary(lme_3)


anova(lme_1,lme_2,lme_3) # choose lme_3 as the final model



# three stages
lme_1 <- lme(shannon~day+phase, 
             random=~1|individual,method ="ML",data=dat)
summary(lme_1)

lme_2 <- lme(shannon~day+I(day^2)+phase, 
             random=~1|individual,method ="ML",data=dat)
summary(lme_2) 

lme_3 <- lme(shannon~day+I(day^2)+I(day^3)+phase, 
             random=~1|individual,method ="ML",data=dat)
summary(lme_3) #day, day^2, day^3


lme_4 <- lme(shannon~day+I(day^2)+I(day^3)+stage, 
             random=~1|individual,method ="ML",data=dat)
summary(lme_4) #day, day^2, day^3

anova(lme_1,lme_2,lme_3) # choose lme_3 as the final model



# three stages +diet change
lme_1 <- lme(shannon~day+phase+stage, 
             random=~1|individual,method ="ML",data=dat)
summary(lme_1)

lme_2 <- lme(shannon~day+I(day^2)+phase+stage, 
             random=~1|individual,method ="ML",data=dat)
summary(lme_2) 

lme_3 <- lme(shannon~day+I(day^2)+I(day^3)+phase+stage, 
             random=~1|individual,method ="ML",data=dat)
summary(lme_3) #day, day^2, day^3

vif(lme_3)


anova(lme_1,lme_2,lme_3) # choose lme_3 as the final model


#-----------------full model--------
lme_full <- lme(shannon~day+phase+stage+sex+group, 
             random=~1|individual,method ="ML",data=dat)
summary(lme_full) #day, day^2, day^3

lme_del1 <- lme(shannon~day+phase+stage+sex, 
                random=~1|individual,method ="ML",data=dat)
summary(lme_del1) #day, day^2, day^3

lme_del2 <- lme(shannon~day+phase+stage, 
                random=~1|individual,method ="ML",data=dat)
summary(lme_del2) #day, day^2, day^3

lme_del3 <- lme(shannon~day+stage, 
                random=~1|individual,method ="ML",data=dat)
summary(lme_del3) #day, day^2, day^3

lme_del4 <- lme(shannon~day+stage, 
                random=~1|individual,method ="ML",data=dat)
summary(lme_del4) #day, day^2, day^3

anova(lme_full, lme_del1,lme_del2, lme_del3, lme_del4)


#-------------------beta diversity models----------
#beta  diversity
betadiv= read.table("bray_curtis.txt", header=T, row.names= 1, sep="\t")

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



frame <- rbind(p[[1]],p[[3]],p[[6]],p[[9]],p[[12]],p[[15]],
               p[[18]],p[[21]], p[[24]], p[[27]],p[[30]],p[[37]],
               p[[44]])

#去除所有NA的值
frame <- na.omit(frame)

#建模
lm_day <- lm(value~day+I(day^2)+I(day^3), frame)
summary(lm_day)


#----------------fit the model


#--------------final model----linear mixed-effects models
if(!require(nlme))install.packages("nlme")
library(nlme)

# add a column for two stages
frame$stage <- frame$day
frame$stage[frame$stage<=24] <- "one"
frame$stage[frame$stage!="one"] <- "two"


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
summary(lme_1)#stage 


lme_2 <- lme(beta_distance~poly(Day,2)+Stage, 
             random=~1|Individual,method ="ML",data=frame1)
summary(lme_2)#day^2


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



# three stages and diet changes
lme_1 <- lme(beta_distance~Day+phase+Stage,
             random=~1|Individual,method ="ML",data=frame1)
summary(lme_1)
# nlme packages could not deal with this model. 
# 
# lme_2 <- lme(beta_distance~poly(Day,2)+phase+Stage, 
#              random=~1|Individual,method ="ML",data=frame1)
# summary(lme_2)# stage 2
# 
# 
# lme_3 <- lme(beta_distance~poly(Day,3)+phase+Stage, 
#              random=~1|Individual,method ="ML",data=frame1)
# summary(lme_3)#stage 3 
# 
# anova(lme_1,lme_2,lme_3)
# # choose lme_3


if(!require(lme4))install.packages(lme4)
if(!require(lmerTest))install.packages("lmerTest")
library(lme4)
library(lmerTest)
library(car)


lme_1 <- lmer(beta_distance~Day+(1|phase)+Stage+(1|Individual),data=frame1)
summary(lme_1)
vif(lme_1)

lme_2 <- lmer(beta_distance~poly(Day,2)+(1|phase)+Stage+(1|Individual), 
             data=frame1)
summary(lme_2)# stage 2


lme_3 <- lmer(beta_distance~poly(Day,3)+(1|phase)+Stage+(1|Individual), 
             data=frame1)
summary(lme_3)#stage 3 

lme_3 <- lmer(beta_distance~poly(Day,3)+(Stage|phase)+Stage+(1|Individual), 
              data=frame1)
summary(lme_3)#stage 3 



vif(lme_3)

anova(lme_1,lme_2,lme_3)
# choose lme_3


