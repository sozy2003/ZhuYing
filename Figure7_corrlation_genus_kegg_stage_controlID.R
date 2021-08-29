# Ying Zhu, so_zy2003@126.com
# 2021-07



# clean environment object
rm(list=ls()) 

#load variables
load(file = "kegg_stage_dif_210706.Rdata")
load(file = "genus_stage_dif_210706.Rdata")

#load packages
if(!require(Hmisc))install.packages("Hmisc")
library("Hmisc")
# differential genus
sub_g.t <- as.data.frame(t(genus_stage_dif))
#differential kegg pathways
count.t <- as.data.frame(t(kegg_stage_dif))


#将所有成对的比较
filter <- list()
#设置菌属数目
m <- ncol(sub_g.t)
for (j in 1:m){
  p <- apply(count.t,2, function(i){
    rcorr(i, sub_g.t[,j],type="spearman" )#"spearman"
  })
  q <- unlist(p)
  # generate matrix
  q <- matrix(q, ncol=12,byrow=T)
  genus_kegg_cor <- as.data.frame(q[,c(2,6,10)])
  genus_kegg_cor$kegg <- rownames(kegg_stage_dif) #
  colnames(genus_kegg_cor) <- c("R2", "n", "P", "KEGG")
  dim(genus_kegg_cor)
  #delete missing data
  genus_kegg_cor <- na.omit(genus_kegg_cor)
  dim(genus_kegg_cor)
  # filter[[j]] <- genus_kegg_cor[genus_kegg_cor$P<0.05,]
  filter[[j]] <- genus_kegg_cor 
}


#correlation between 21 genera and 10 kegg pathways
relation <- list()

for (i in 1:m){
  relation[[i]]<-as.data.frame(filter[i])
  relation[[i]]$genus<-rownames(genus_stage_dif)[i]#
}

relation.all <- do.call(rbind, relation)
dim(relation.all)

# generate matrix
genus_kegg_relation_all <- matrix(relation.all[,1], nrow=ncol(count.t), byrow=F)
# change to data frame
genus_kegg_relation_all <- as.data.frame(genus_kegg_relation_all)
colnames(genus_kegg_relation_all) <- rownames(genus_group)
rownames(genus_kegg_relation_all) <- rownames(kegg_group)


#dif_genus_stage$genus=="Unassigned"
whi <- which(dif_genus_stage$genus=="Unassigned")
dif_genus_stage$genus[whi] <- paste(dif_genus_stage$family[whi],"_Unassigned", sep = "")
dif_genus_stage <- dif_genus_stage[order(dif_genus_stage$regulated),]

#reorder the column
genus_kegg_relation_all <- genus_kegg_relation_all[ ,dif_genus_stage$genus]
# reorder the row
genus_group$genus <- rownames(genus_group)
genus_group <- genus_group[match(dif_genus_stage$genus,genus_group$genus),]


# row annotation
annotation_row <- kegg_group

rownames(annotation_row)==rownames(genus_kegg_relation_all)

#column annotation
annotation_col <- cbind(phylum=genus_group$phylum, 
                        
                                regulate=dif_genus_stage$regulated)
annotation_col <- annotation_col [,-2] 
annotation_col <- as.data.frame(annotation_col)

rownames(annotation_col) <- colnames(genus_kegg_relation_all)
colnames(annotation_col) <- "Phylum"



# load package
library(pheatmap)


p <- pheatmap(genus_kegg_relation_all, cluster_cols = F, 
              cluster_rows = T, 
              # scale = "row",
              annotation_row = annotation_row,
              annotation_col = annotation_col,
              show_rownames = T,
              filename = "heatmap_genus_kegg_stage_relation_controlID.pdf", 
              width = 20, height = 6)
