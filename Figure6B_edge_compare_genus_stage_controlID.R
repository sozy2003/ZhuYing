# Ying Zhu, so_zy2003@126.com
# 2021-07

# clear variables
rm(list = ls()) 
options(stringsAsFactors = F)

# #load data
filter_count <- read.table("count_6Genus.txt", row.names= 1,  header=T, sep="\t")#get the count matrix
colnames(filter_count)


# Design of experiment
metadata = read.table("metadata.txt", header=T, row.names= 1, sep="\t") 

#rearrange the column of filter_count according to the row names of metadata
filter_count = filter_count[, rownames(metadata)]

exprSet <- filter_count
exprSet <- exprSet[,!colnames(exprSet)%in%c("A4D24","B2D24","C1D24")]
dim(exprSet)
exprSet[1:6,1:6]

#cross filtering
metadata <- metadata[rownames(metadata)%in%colnames(exprSet),]
dim(metadata)

#define variable for inviduals and stages
#variable including two stages
group_list <- metadata$Stage #stage one: loach paste; stage two: loach paste+fresh loach
group_list <- as.factor(group_list)
table(group_list)

#variable including all the individals
individual_list <- as.factor(metadata$Individual)

#check the data
data.frame(Sampe=rownames(metadata), individual_list,group_list)


# install and load edgeR package
if(!require(edgeR))install.packages("edgeR")
library(edgeR)

# the design matrix is formed from an  additive model for paired designs.
# we are interested in the average effect of the diet over a population of chicks. 
design <- model.matrix(~0+individual_list+group_list)
rownames(design) <- colnames(exprSet)
colnames(design) #stage one as a control
# This creats a design matrix with 15 columns: 14 for the individuals and 1 more for the differences between stages.

# construct DGEList object
DEG <- DGEList(counts=exprSet)

# check the library size 
DEG$samples
DEG$samples$lib.size <- colSums(DEG$counts)
DEG$samples

#calculate the overall count for each sample
library(dplyr)
apply(exprSet,2,sum)
# This is just for checking data.

# normalization
DEG <- calcNormFactors(DEG)

# caculate model parameters
DEG <- estimateGLMCommonDisp(DEG,design)
DEG <- estimateGLMTrendedDisp(DEG, design)
DEG <- estimateGLMTagwiseDisp(DEG, design)

# fit glms model
fit <- glmFit(DEG, design)

  # differential analysis
  lrt <- glmLRT(fit)
  
# filtering the results
  DEG_edgeR <- as.data.frame(topTags(lrt, n=nrow(DEG)))
  
  fc_cutoff <- 2
  fdr <- 0.05
  
  DEG_edgeR$regulated <- "normal"
  
  loc_up <- intersect(which(DEG_edgeR$logFC>log2(fc_cutoff)),which(DEG_edgeR$FDR<fdr))
  loc_down <- intersect(which(DEG_edgeR$logFC< (-log2(fc_cutoff))),which(DEG_edgeR$FDR<fdr))
  
  DEG_edgeR$regulated[loc_up] <- "up"
  DEG_edgeR$regulated[loc_down] <- "down"
  

#select the significantly differential ASVs
deg <- DEG_edgeR[DEG_edgeR$regulated!="normal",]
dim(deg)
deg
table(deg$regulated)

deg$whole <- rownames(deg)

# separate the column "whole" to kingdom, phylum,class, order, family and genus
library(tidyr)
dif_g <-  separate(deg,whole,
                   into=c("kindom", "phylum", "class", "order", "family", "genus"),
                   sep="\\|")

colnames(dif_g)

#filtering raw data for plotting based on results of above differential analysis 
sub_g <- exprSet[rownames(exprSet)%in%rownames(dif_g),]
dim(sub_g)

# rename the row names of sub_g
genus_name <- as.data.frame(rownames(sub_g))
dim(genus_name)

colnames(genus_name)[1] <- "whole"
genus_name<-  separate(genus_name,whole,
                       into=c("kindom", "phylum", "class", "order", "family", "genus"),
                       sep="\\|")
colnames(genus_name)


whi <- which(genus_name$genus=="Unassigned")
genus_name$genus[whi] <- paste(genus_name$family[whi],"_Unassigned", sep = "")

rownames(sub_g) <- genus_name$genus

sub_g <- sub_g[,rownames(metadata)]

# calculate the mean count for each day

sampFile = as.data.frame(metadata$Day,row.names = row.names(metadata))
colnames(sampFile)[1] = "group"
mat_t = t(sub_g)
mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,-1]
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean

otu_norm_group = do.call(rbind, mat_mean)[-1,]
colnames(otu_norm_group) = mat_mean$group

library(dplyr)
otu_norm_group <- apply(otu_norm_group, 2, as.numeric)
rownames(otu_norm_group) <- colnames(mat_t)

# column annotation for heatmap 
stage <- c("one", "one","one","one","one","one","one","one",
           "one","two","two","two","two")
stage <- as.data.frame(stage)
rownames(stage) <- colnames(otu_norm_group)

rownames(otu_norm_group)==rownames(sub_g)

#row annotation for heatmap
rowgroup <- genus_name$phylum
rowgroup <- as.data.frame(rowgroup)
colnames(rowgroup) <- "phylum"

rownames(rowgroup) <- rownames(otu_norm_group)



# load package pheatmap
if(!require(pheatmap))install.packages("pheatmap")
library(pheatmap)

# pheatmap(otu_norm_group,scale="row",cluster_cols = F, cluster_rows = T,
#          annotation_col=stage, annotation_row=rowgroup,
#          filename = "heatmap_df_stage_controlID.pdf", width = 10, height = 5)


# save variables for the relation analysis 
genus_stage_dif <- as.data.frame(otu_norm_group)

genus_group <- rowgroup
genus_group$phylum <- genus_name$phylum
dif_genus_stage <- dif_g
# save(genus_stage_dif, genus_group, dif_genus_stage,
#      file="genus_stage_dif_210706.Rdata")

