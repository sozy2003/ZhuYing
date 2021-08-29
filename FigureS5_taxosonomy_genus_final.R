# Ying Zhu, so_zy2003@126.com
# 2021-07

rm(list=ls()) # clean enviroment object
library("reshape2", quietly=T, warn.conflicts=F)
library(ggalluvial)

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
design = read.table("metadata_weight.txt", header=T, row.names= 1, sep="\t") 

design$phase <- design$Day

design$phase[design$Day<12] <- "p1"  #stage one in Figure 1
design$phase[design$Day%in%c(12,15,18,21,24)] <- "p2"  #stage two in Figure 1
design$phase[design$Day>24] <- "p3" #stage three in Figure 1

#raw reads count of each OTU in each sample
otu_table = read.delim("otu_table.txt", row.names= 1,  header=T, sep="\t")

# taxonomy for each OTU
taxonomy = read.delim("taxonomy.txt", row.names= 1,header=F, sep="\t")
colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species")

taxonomy <- taxonomy[-1,]


tax_count = merge(taxonomy, otu_table, by="row.names")

# genus summary
tax_count_sum = aggregate(tax_count[,-(1:8)], by=tax_count[7], FUN=sum) # mean

rownames(tax_count_sum) = tax_count_sum$genus

tax_count_sum = tax_count_sum[,-1]

per = t(t(tax_count_sum)/colSums(tax_count_sum,na=T)) * 100 # normalization to total 100



mean_sort = per[(order(-rowSums(per))), ] # decrease sort
colSums(mean_sort)#验证百分比总和为1

# filtering 
mean_sort=as.data.frame(mean_sort)
other = colSums(mean_sort[10:dim(mean_sort)[1], ])
mean_sort = mean_sort[1:(10-1), ]
mean_sort = rbind(mean_sort,other)
rownames(mean_sort)[10] = c("Low Abundance")

# write.table(mean_sort, file="Top10phylum_ProClass.txt", append = F, sep="\t", quote=F, row.names=T, col.names=T)

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

mat_mean_final$genus = rownames(mat_mean_final)


data_all = as.data.frame(melt(mat_mean_final, id.vars=c("genus")))


p = ggplot(data_all, aes(x=variable, y = value, fill = genus )) + 
  geom_bar(stat = "identity",position="fill", width=0.7)+ 
  scale_y_continuous(labels = scales::percent) + 
  xlab("Groups")+ylab("Percentage (%)")+main_theme+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
tax.stack=p
# ggsave("tax_stack_genus_top9_group.pdf", tax.stack, width = 8, height = 5)


# alluvium
p = ggplot(data = data_all, aes(x = variable, y = value, alluvium = genus, stratum=genus)) +
  geom_alluvium(aes(fill = genus), alpha = 0.75) +
  geom_stratum(aes(fill=genus))+
  labs(x="Day", y="Relative Abundance (%)")+
  main_theme 
tax.alluvium=p
# ggsave("tax_alluvium_genus_top9.pdf", tax.alluvium, width = 8, height = 5)


# facet
mat_t2_melt <- melt(mat_t2, id.vars = "Day")
top10genus_relative <- ggplot(mat_t2_melt,aes(x=as.factor(Day),y=value, fill=variable))+
  geom_boxplot()+
  geom_jitter(size=0.5, alpha=0.5)+
  facet_wrap(~variable)+theme_classic()+
  labs(x="Day", y="Relative Abundance (%)")
# ggsave("top10genus_relative_boxplot.pdf",top10genus_relative, width = 12, height = 8 )


#---------------linear model------------
#age effect
colnames(mat_t2)[2] <- "Escherichia"
lm1 <- lm(Escherichia~Day, mat_t2)
summary(lm1)# significant R2=0.1529 

lm2<- lm(Clostridium_sensu_stricto~Day, mat_t2)
summary(lm2)#significant，R2=0.02302


lm3 <- lm(Plesiomonas~Day, mat_t2)
summary(lm3)#significant R2=0.02302

lm4<- lm(Catellicoccus~Day,  mat_t2)
summary(lm4)# #significant R2=0.1233


# two stages
mat_t3 = merge(mat_t,design[c("Day", "Stage", "Individual", "phase")], by="row.names")
rownames(mat_t3) <- mat_t3$Row.names
mat_t3 <- mat_t3[,-1]
colnames(mat_t3)[1] <- "Escherichia"

lm1 <- lm(Escherichia~Day+Stage, mat_t3)
summary(lm1)#stage significant

lm2<- lm(Clostridium_sensu_stricto~Day+Stage, mat_t3)
summary(lm2)


lm3 <- lm(Plesiomonas~Day+Stage, mat_t3)
summary(lm3)

lm4<- lm(Catellicoccus~Day+Stage,  mat_t3)
summary(lm4)#stage significant

#---------------linear mixed models

library(nlme)

# two stages
lme_1 <- lme(Escherichia~Day+Stage, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_1)#stage significant


lme_2 <- lme(Clostridium_sensu_stricto~Day+Stage, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_2)

lme_3 <- lme(Plesiomonas~Day+Stage, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_3)


lme_4 <- lme(Catellicoccus~Day+Stage, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_4) #stage significant

# three stages
lme_1 <- lme(Escherichia~Day+phase, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_1)


lme_2 <- lme(Clostridium_sensu_stricto~Day+phase, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_2)#Day negative
 
lme_3 <- lme(Plesiomonas~Day+phase, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_3) #Day, negative


lme_4 <- lme(Catellicoccus~Day+phase, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_4) #stage three 





# three stages+diet change
lme_1 <- lme(Escherichia~Day+phase+Stage, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_1)


lme_2 <- lme(Clostridium_sensu_stricto~Day+phase+Stage, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_2)

lme_3 <- lme(Plesiomonas~Day+phase+Stage, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_3)


lme_4 <- lme(Catellicoccus~Day+phase+Stage, 
             random=~1|Individual,method ="ML",data=mat_t3)
summary(lme_4) #stage three
