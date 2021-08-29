# Ying Zhu, so_zy2003@126.com
# 2021-07



# clear variables
rm(list = ls()) 
options(stringsAsFactors = F)

#load otu table, keep the ASVs whose mean abundance is larger than 0.01%
filter_count <- read.table("count_8OTU.txt", row.names= 1,  header=T, sep="\t")
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

#define variable for individuals and stages
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



# calculate the mean count for stage group

sub_g = t(t(exprSet)/colSums(exprSet,na=T)) * 100 # normalization to total 100

sampFile = as.data.frame(metadata$Stage,row.names = row.names(metadata))
colnames(sampFile)[1] = "group"
mat_t = t(sub_g)
mat_t2 = merge(sampFile, mat_t, by="row.names")
mat_t2 = mat_t2[,-1]
mat_mean = aggregate(mat_t2[,-1], by=mat_t2[1], FUN=mean) # mean

otu_norm_group = do.call(rbind, mat_mean)[-1,]
colnames(otu_norm_group) = mat_mean$group

#
temp <- merge(otu_norm_group,sub_g,by="row.names")
rownames(temp) <-temp$Row.names 
temp <- temp[,-1]

two_one_compare <-merge(DEG_edgeR,temp, by="row.names") 
two_one_compare <- two_one_compare[(order(two_one_compare$PValue)),]

colnames(two_one_compare)[1:9] <- c("ID", "log2FC", "log2CPM","LR", "PValue",
                                    "FDR", "level", "MeanA", "MeanB")

two_one_compare$level[two_one_compare$level=="normal"] <- "NotSig"
two_one_compare$level[two_one_compare$level=="up"] <- "Enriched"
two_one_compare$level[two_one_compare$level=="down"] <- "Depleted"

two_one_compare <- two_one_compare[,colnames(two_one_compare)!="LR"]

# save results
# write.table(two_one_compare,"two-one_controlID.txt",row.names= F, sep="\t")
