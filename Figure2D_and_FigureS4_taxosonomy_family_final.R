# Ying Zhu, so_zy2003@126.com
# 2021-07

# clean environment variables
rm(list=ls()) 

#install and load the packages
if(!require(reshape2))install.packages("reshape2")
if(!require(ggalluvial))install.packages("ggalluvial")

library("reshape2", quietly=T, warn.conflicts=F)
library("ggalluvial")

# Set ggplot2 ploting parameter
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

#  raw reads count of each OTU in each sample
otu_table = read.delim("otu_table.txt", row.names= 1,  header=T, sep="\t")

# taxonomy for each OTU
taxonomy = read.delim("taxonomy.txt", row.names= 1,header=F, sep="\t")
colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species")

taxonomy <- taxonomy[-1,]

tax_count = merge(taxonomy, otu_table, by="row.names")

# group by family
tax_count_sum = aggregate(tax_count[,-(1:8)], by=tax_count[6], FUN=sum) # mean

rownames(tax_count_sum) = tax_count_sum$family
# generate numeric matrix
tax_count_sum = tax_count_sum[,-1]
# normalization to total 100
per = t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 


# descending sort based on abundance 
mean_sort = per[(order(-rowSums(per))), ] # decrease sort
colSums(mean_sort)

# top10,including Low abundance
mean_sort=as.data.frame(mean_sort)
other = colSums(mean_sort[10:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(10-1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[10] = c("Low Abundance")

# write.table(mean_sort, file="Top10phylum_ProClass.txt", append = F, sep="\t", quote=F, row.names=T, col.names=T)

topN=rownames(mean_sort)


mean_sort <- mean_sort[,rownames(design)]


# plot for each time points
mat=mean_sort

mat_t = t(mat)

mat_t2 = merge(design[c("Day")], mat_t, by="row.names")

mat_t2 = mat_t2[,-1]

mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final = do.call(rbind, mat_mean)[-1,]

geno = mat_mean$Day
colnames(mat_mean_final) = geno

mat_mean_final = as.data.frame(mat_mean_final)

mat_mean_final$family = rownames(mat_mean_final)

data_all = as.data.frame(melt(mat_mean_final, id.vars=c("family")))


#plotting
p = ggplot(data_all, aes(x=variable, y = value, fill = family )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+ylab("Percentage (%)")+main_theme
tax.stack=p
tax.stack
# ggsave("tax_stack_family_top9_group.pdf", tax.stack, width = 8, height = 5)


# alluvium
p = ggplot(data = data_all, aes(x = variable, y = value, alluvium = family, stratum=family)) +
  geom_alluvium(aes(fill = family), alpha = 0.75) +
  geom_stratum(aes(fill=family))+
  labs(x="Day", y="Relative Abundance (%)")+
main_theme 
tax.alluvium=p
tax.alluvium
# ggsave("tax_alluvium_family_top9.pdf", tax.alluvium, width = 8, height = 5)


#boxplot for the top 10 families
mat_t2_melt <- melt(mat_t2, id.vars = "Day")
colnames(mat_t2_melt)[2] <- "Family"
top10family_relative <- ggplot(mat_t2_melt,aes(x=as.factor(Day),y=value, fill=Family))+
  geom_boxplot()+
 geom_jitter(size=0.5, alpha=0.5)+
  facet_wrap(~Family)+theme_classic()+
  labs(x="Day", y="Relative Abundance (%)")
top10family_relative
# ggsave("top10family_relative_boxplot.pdf",top10family_relative, width = 12, height = 8 )


# fit the model
mat_t3 <-  merge(design[c("Day", "Stage", "Individual", "phase", "food", "Group", "Sex")], mat_t, by="row.names")
rownames(mat_t3) <- mat_t3$Row.names
mat_t3 <- mat_t3[,-1]


# what is influence of the factors on abundance of dominant families 
#age, linear model
lm1 <- lm(Enterobacteriaceae~Day, mat_t3)
summary(lm1)#系数负值R2=0.1391 significant

lm2<- lm(Clostridiaceae_1~Day, mat_t3)
summary(lm2)#系数负值，R2=0.02302 significant 


lm3 <- lm(Carnobacteriaceae~Day, mat_t3)
summary(lm3)#系数正值R2=0.1231

lm4<- lm(Peptostreptococcaceae~Day,  mat_t3)
summary(lm4)# #系数正值R2= 0.2293 

#three stages
table(mat_t3$phase)
lm1 <- lm(Enterobacteriaceae~phase, mat_t3)
summary(lm1)#significant

lm2<- lm(Clostridiaceae_1~phase, mat_t3)
summary(lm2)#non significant 

lm3 <- lm(Carnobacteriaceae~phase, mat_t3)
summary(lm3)#significant

lm4<- lm(Peptostreptococcaceae~Day,  mat_t3)
summary(lm4)#significant 


# sex
lm1 <- lm(Enterobacteriaceae~Sex, mat_t3)
summary(lm1)#

lm2<- lm(Clostridiaceae_1~Sex, mat_t3)
summary(lm2)#non significant 

lm3 <- lm(Carnobacteriaceae~Sex, mat_t3)
summary(lm3)#non significant

lm4<- lm(Peptostreptococcaceae~Sex,  mat_t3)
summary(lm4)#non significant 

# family group
lm1 <- lm(Enterobacteriaceae~Group, mat_t3)
summary(lm1)##non significant

lm2<- lm(Clostridiaceae_1~Group, mat_t3)
summary(lm2)#non significant 

lm3 <- lm(Carnobacteriaceae~Group, mat_t3)
summary(lm3)#non significant

lm4<- lm(Peptostreptococcaceae~Group,  mat_t3)
summary(lm4)#non significant 


# food
lm1 <- lm(Enterobacteriaceae~food, mat_t3)
summary(lm1)##significant

lm2<- lm(Clostridiaceae_1~food, mat_t3)
summary(lm2)#non significant 

lm3 <- lm(Carnobacteriaceae~food, mat_t3)
summary(lm3)#significant

lm4<- lm(Peptostreptococcaceae~food,  mat_t3)
summary(lm4)#significant 

#three possible predictors: phase, age and food
# only incluede age in the model due to the strong relationship between age and food

# linear mixed-effects models
if(!require(nlme))install.packages("nlme")
library(nlme)
# two stages
lme_1 <- lme(Enterobacteriaceae~Day+Stage, 
                 random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_1)#stage


lme_2 <- lme(Clostridiaceae_1~Day+Stage, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_2)

lme_3 <- lme(Carnobacteriaceae~Day+Stage, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_3)#stage


lme_4 <- lme(Peptostreptococcaceae~Day+Stage, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_4)#Day


#three stages
lme_1 <- lme(Enterobacteriaceae~Day+phase, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_1)


lme_2 <- lme(Clostridiaceae_1~Day+phase, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_2)#Day, negative

lme_3 <- lme(Carnobacteriaceae~Day+phase, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_3)#stage three


lme_4 <- lme(Peptostreptococcaceae~Day+phase, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_4)#Day, positive



#three stages+diet change
lme_1 <- lme(Enterobacteriaceae~Day+phase+Stage, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_1)


lme_2 <- lme(Clostridiaceae_1~Day+phase+Stage, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_2)

lme_3 <- lme(Carnobacteriaceae~Day+phase+Stage, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_3)#stage three


lme_4 <- lme(Peptostreptococcaceae~Day+phase+Stage, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_4)#Day, positive, phase three?



