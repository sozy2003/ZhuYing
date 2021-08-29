# Ying Zhu, so_zy2003@126.com
# 2021-07


rm(list=ls()) # clean environment object

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

# Design of experiment
design = read.table("metadata.txt", header=T, row.names= 1, sep="\t") 

#  raw reads count of each OTU in each sample
otu_table = read.delim("otu_table.txt", row.names= 1,  header=T, sep="\t")

# taxonomy for each OTU, tab seperated
taxonomy = read.delim("taxonomy.txt", row.names= 1,header=F, sep="\t")
colnames(taxonomy) = c("kingdom","phylum","class","order","family","genus","species")

taxonomy <- taxonomy[-1,]


tax_count = merge(taxonomy, otu_table, by="row.names")



# summary phyla
tax_count_sum = aggregate(tax_count[,-(1:8)], by=tax_count[3], FUN=sum) # mean

rownames(tax_count_sum) = tax_count_sum$phylum

tax_count_sum = tax_count_sum[,-1]


phylum<- as.matrix(rowSums(tax_count_sum), ncol=1)

phylum <- as.data.frame(phylum)
colnames(phylum)[1] <- "count"
phylum$relative <- phylum$count/sum(phylum)
sum(phylum$relative)


# calculte the sum of all the sample for each OTU

otu_sum <- as.matrix(rowSums(otu_table)/sum(otu_table),ncol=1)
otu_sum <- as.data.frame(otu_sum)
colnames(otu_sum)[1] <- "relative"
#top 20
sum(otu_sum$relative[1:20])

#annotation
otu_sum_tax <- merge(taxonomy, otu_sum, by="row.names")
sum(otu_sum_tax$relative)
#select phyla Firmicutes and Proteobacteria 
firm <- otu_sum_tax[otu_sum_tax$phylum=="Firmicutes",]
firm <- firm[order(-firm$relative),]
rownames(firm) <- firm$Row.names

proto <- otu_sum_tax[otu_sum_tax$phylum=="Proteobacteria",]
proto  <- proto [order(-proto $relative),]
rownames(proto) <- proto$Row.names
#Ratios of the relative 
#abundances of the 20 most abundant OTUs classified as Proteobacteria  and Firmicutes 
#to the total relative abundances of Proteobacteria and Firmicutes, respectively. 

firm_ratio <- firm$relative/phylum$relative[rownames(phylum)=="Firmicutes"]
firm_ratio <- as.matrix(firm_ratio,ncol=1)
rownames(firm_ratio) <-rownames(firm)
colnames(firm_ratio) <- "relative"
firm_ratio <- as.data.frame(firm_ratio)
firm_ratio$asv <- firm$Row.names


proto_ratio <- proto$relative/phylum$relative[rownames(phylum)=="Proteobacteria"]
proto_ratio <- as.matrix(proto_ratio,ncol=1)
rownames(proto_ratio) <-rownames(proto)
colnames(proto_ratio) <- "relative"
proto_ratio  <- as.data.frame(proto_ratio )
proto_ratio$asv <- proto$Row.names

top_firm <- firm_ratio[1:20,]
top_proto <-proto_ratio[1:20,]

# write.table(top_firm, file="top_firm.txt", 
            sep="\t")

# write.table(top_proto, file="top_proto.txt", 
            sep="\t")


sig_otu <- c("ASV_10", "ASV_11","ASV_12",
             "ASV_15", "ASV_16", "ASV_17", "ASV_18", 
             "ASV_19")
top_firm <- read.table("top_firm.txt", sep="\t", header=T,row.names = 1)
top_firm$sig <- top_firm$asv
top_firm$sig[top_firm$sig%in%sig_otu] <- "sig"
top_firm$sig[top_firm$sig!="sig"] <- "nonsig"
ratio_firm_p <- ggplot(data=top_firm)+
  geom_bar(aes(x=reorder(asv, -relative), y=relative, fill=sig), stat = "identity")+
  main_theme+
  labs(x="", y="Mean firm. OTU abundance/ Total Fir. abundance")+
  theme(axis.text=element_text(angle=45,vjust=1, hjust=1))
ratio_firm_p
# ggsave("ratio_firm_p.pdf", ratio_firm_p, width = 20, height = 10, units = "cm")

#proteo
top_proto <- read.table("top_proto.txt", sep="\t", header=T,row.names = 1)
colnames(top_proto) <- c("relative","asv")
top_proto$sig <- top_proto$asv
top_proto$sig[top_proto$sig%in%sig_otu] <- "sig"
top_proto$sig[top_proto$sig!="sig"] <- "nonsig"
ratio_proto_p <- ggplot(data=top_proto)+
  geom_bar(aes(x=reorder(asv, -relative), y=relative,fill=sig), stat = "identity")+
  main_theme+labs(x="", y="Mean pro. OTU abundance/ Total Pro. abundance")+
  theme(axis.text=element_text(angle=45,vjust=1, hjust=1))
ratio_proto_p
# ggsave("ratio_proto_p.pdf", ratio_proto_p, width = 20, height = 10, units = "cm")




