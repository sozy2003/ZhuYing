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

# Design of experiment
design = read.table("metadata_weight.txt", header=T, row.names= 1, sep="\t") 

#aggregate a column for three stages
design$phase <- design$Day

design$phase[design$Day<12] <- "p1"  #stage one in Figure 1
design$phase[design$Day%in%c(12,15,18,21,24)] <- "p2"  #stage two in Figure 1
design$phase[design$Day>24] <- "p3" #stage three in Figure 1

# raw reads count of each ASV in each sample
otu_table <-  read.delim("otu_table.txt", row.names= 1,  header=T, sep="\t")

#  taxonomy for each OTU
taxonomy <-  read.delim("taxonomy.txt", row.names= 1,header=F, sep="\t")
colnames(taxonomy) <-  c("kingdom","phylum","class","order","family","genus","species")

taxonomy <- taxonomy[-1,]

# summary for each interesting tax

# select p__Proteobacteria line
idx <-  taxonomy$phylum == "Proteobacteria"

taxonomy$full <- as.character(taxonomy$phylum) 

taxonomy[idx,]$full <- as.character(taxonomy[idx,]$class)
# add annotation for otu table
tax_count <-  merge(taxonomy, otu_table, by="row.names")

# group by column "full"
tax_count_sum <-  aggregate(tax_count[,-(1:9)], by=tax_count[9], FUN=sum) # mean
# rownames
rownames(tax_count_sum) <-  tax_count_sum$full
# generate numeric matrix 
tax_count_sum <-  tax_count_sum[,-1]
# normalization to total 100 
per <-  t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 


# descending sort based on abundance 
mean_sort <- per[(order(-rowSums(per))), ] # decrease sort
colSums(mean_sort)

# top10,including Low abundance
mean_sort<-as.data.frame(mean_sort)
other <- colSums(mean_sort[10:dim(mean_sort)[1], ])
mean_sort <- mean_sort[1:(10-1), ]
mean_sort <- rbind(mean_sort,other)
rownames(mean_sort)[10] <- c("Low Abundance")

# write.table(mean_sort, file="Top10phylum_ProClass.txt", append = F, sep="\t", quote=F, row.names=T, col.names=T)
# top 10
topN<-rownames(mean_sort)

mean_sort <- mean_sort[,rownames(design)]

# plot for samples
mean_sort$phylumpro <- rownames(mean_sort)

data_all <- as.data.frame(melt(mean_sort, id.vars=c("phylumpro")))

data_all <- merge(data_all, design[c("Day", "Individual")], by.x="variable", by.y = "row.names")

p <- ggplot(data_all, aes(x=variable, y = value, fill = phylumpro )) + 
  geom_bar(stat = "identity",position="fill", width=1)+ 
  scale_y_continuous(labels = scales::percent) + 
 facet_grid( ~ Day, scales = "free_x", switch = "x") +  main_theme +
 theme(axis.ticks.x = element_blank(), legend.position="top", axis.text.x = element_blank(), strip.background = element_blank())+
  xlab("Groups")+ylab("Percentage (%)")           

p
# ggsave("tax_stack_phylumpro_sample.pdf", p, width = 10, height = 6)


# plot for each time points
mat <-  mean_sort[,1:(dim(mean_sort)[2]-1)]

mat_t <-  t(mat)

mat_t2 <-  merge(design[c("Day")], mat_t, by="row.names")

rownames(mat_t2) <- mat_t2$Row.names
mat_t2 <-  mat_t2[,-1]

mat_mean <-  aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean
mat_mean_final <-  do.call(rbind, mat_mean)[-1,]

geno <-  mat_mean$Day
colnames(mat_mean_final) <-  geno

mat_mean_final <-  as.data.frame(mat_mean_final)

mat_mean_final$phylumpro <-  rownames(mat_mean_final)


data_all <-  as.data.frame(melt(mat_mean_final, id.vars=c("phylumpro")))


p <-  ggplot(data_all, aes(x=variable, y = value, fill = phylumpro )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+ylab("Percentage (%)")+main_theme+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
tax.stack<- p
tax.stack
# ggsave("tax_stack_phylumpro_top9_group.pdf", tax.stack, width = 8, height = 5)


# alluvium
p = ggplot(data = data_all, aes(x = variable, y = value, alluvium = phylumpro, stratum=phylumpro)) +
  geom_alluvium(aes(fill = phylumpro), alpha = 0.75) +
  geom_stratum(aes(fill=phylumpro))+
  labs(x="Day", y="Relative Abundance (%)")+
  main_theme 
tax.alluvium=p

tax.alluvium
# ggsave("tax_alluvium_phylumpro_top9.pdf", tax.alluvium, width = 8, height = 5)


#scatter point plot based on mean value
p=ggplot(data_all ,aes(x=as.numeric(variable),y=value,color=phylumpro))+
  labs(x="Day", y="Relative Abundance") + main_theme+
  geom_point(size=2,alpha=1)+
  scale_x_continuous(breaks = 1:nlevels(data_all$variable), 
                     labels = levels(data_all$variable))
# smooth fit
p1=p+geom_smooth(method = "lm", formula = y ~ poly(x,1))
p2=p+geom_smooth(method = "lm", formula = y ~ poly(x,2))
p3=p+geom_smooth(method = "lm", formula = y ~ poly(x,3))
p1
p2
p3

# ggsave("tax_mean_phy_top10_linear.pdf", p1, width = 8, height = 5)
# ggsave("tax_mean_phy_top10_polylinear.pdf", p3, width = 8, height = 5)


# boxplot
mat_t2_melt <- melt(mat_t2, id.vars = "Day")
colnames(mat_t2_melt)[2] <- "Phylum"
top10phylum_relative <- ggplot(mat_t2_melt,aes(x=as.factor(Day),y=value, fill=Phylum))+
  geom_boxplot()+
  geom_jitter(size=0.5, alpha=0.5)+
  facet_wrap(~Phylum)+theme_classic()+
  labs(x="Day", y="Relative Abundance (%)")
top10phylum_relative
# ggsave("top10phylum_relative_boxplot.pdf",top10phylum_relative, width = 12, height = 8 )


#check the Fusobacteria 
mat_t2_melt_fuso <- mat_t2_melt[mat_t2_melt$Phylum=="Fusobacteria",]
ggplot(mat_t2_melt_fuso,aes(x=as.factor(Day),y=value))+
  geom_boxplot()+
  geom_jitter(size=0.5, alpha=0.5)+
  # facet_wrap(~Phylum)+theme_classic()+
  labs(x="Day", y="Relative Abundance (%)")


# fit the model
mat_t3 <-  merge(design[c("Day", "Stage", "Individual", "phase", "food", "Group", "Sex")], mat_t, by="row.names")
rownames(mat_t3) <- mat_t3$Row.names
mat_t3 <- mat_t3[,-1]

# what is influence of the factors on abundance of dominant phyla 
#age, linear model
lm1 <- lm(Gammaproteobacteria~Day, mat_t3)
summary(lm1)# significant

lm2<- lm(Firmicutes~Day, mat_t3)
summary(lm2)# significant

# polynomial model
lm3 <- lm(Gammaproteobacteria~poly(Day,2), mat_t3)
summary(lm3)

lm4<- lm(Firmicutes~poly(Day,2),  mat_t2)
summary(lm4)

anova(lm1, lm3) # polynomial model did not improve the fitting
anova(lm2, lm4) # polynomial model did not improve the fitting

lm5 <- lm(Gammaproteobacteria~poly(Day,3), mat_t3)
summary(lm5)

lm6<- lm(Firmicutes~poly(Day,3),  mat_t2)
summary(lm6)

anova(lm1, lm5) # polynomial model did not improve the fitting
anova(lm2, lm6) # polynomial model did not improve the fitting


#three stages
table(mat_t3$phase)
lm1 <- lm(Gammaproteobacteria~phase, mat_t3)
summary(lm1)# significant

lm2<- lm(Firmicutes~phase, mat_t3)
summary(lm2)# significant

#food
lm1 <- lm(Gammaproteobacteria~food, mat_t3)
summary(lm1)# significant

lm2<- lm(Firmicutes~food, mat_t3)
summary(lm2)# significant


#sex
lm1 <- lm(Gammaproteobacteria~Sex, mat_t3)
summary(lm1)# nonsignificant

lm2<- lm(Firmicutes~Sex, mat_t3)
summary(lm2)# significant


#family group
lm1 <- lm(Gammaproteobacteria~Group, mat_t3)
summary(lm1)# nonsignificant

lm2<- lm(Firmicutes~Group, mat_t3)
summary(lm2)# significant


#relation among three possible predictors: phase, age and food

# load the packages
if(!require(Hmisc))install.packages("Hmisc")
library("Hmisc")

#age~food, 相关strong positive relationship between age and food
rcorr(mat_t3$Day, mat_t3$food, type="spearman")#正相关
ggplot(mat_t3, aes(x=Day, y=food))+
  geom_point()

# only incluede age in the model due to the strong relationship between age and food

#------------------final models
# linear mixed-effects models to control individual effects
library(nlme)

# two stages
#Fusobacteria
lme_fuso <- lme(Fusobacteria~Day+Stage, 
                random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_fuso) # 

colnames(mat_t3)

#Gammaproteobacteria
lme_gamma <- lme(Gammaproteobacteria~Day+Stage, 
                random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_gamma) #day P=0.05


# three stages
#Fusobacteriaceae
lme_fuso <- lme(Fusobacteria~Day+phase, 
                random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_fuso)

colnames(mat_t3)

#Gammaproteobacteria
lme_gamma <- lme(Gammaproteobacteria~Day+phase, 
                 random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_gamma)

# chosse lme_fuso and lme_gamma (three stages)




# three stages and diet changes
#Fusobacteriaceae
lme_fuso <- lme(Fusobacteria~Day+phase+Stage, 
                random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_fuso)

colnames(mat_t3)

#Gammaproteobacteria
lme_gamma <- lme(Gammaproteobacteria~Day+phase+Stage, 
                 random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_gamma)

# chosse lme_fuso and lme_gamma (three stages+diet change)??
