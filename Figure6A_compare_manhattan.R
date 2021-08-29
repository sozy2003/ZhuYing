# Ying Zhu, so_zy2003@126.com
# 2021-07

# clear variables
rm(list = ls()) 
options(stringsAsFactors = F)

#install and load packages

site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
package_list = c("ggplot2","dplyr")

for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos=site)
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}


# theme for plotting
main_theme = theme(panel.background=element_blank(), panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"), axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"), axis.text=element_text(color="black", size=7),
                   legend.position="right", legend.background=element_blank(), legend.key=element_blank(), legend.text= element_text(size=7),
                   text=element_text(family="sans", size=7))

# differential results
x = read.table("two-one_controlID.txt", header=T, row.names= 1, sep="\t", stringsAsFactors = F)
# subset
x = x[,1:7]
x = na.omit(x)

# add taxonomy annotation 
taxonomy = read.table("taxonomy.txt", sep = "\t", row.names=1, header=T, stringsAsFactors = F)
taxonomy = taxonomy[rownames(x),]
x = cbind(x, taxonomy)

# add  new columns
x$neglogp = -log10(x$PValue)

x$otu=rownames(x)
x = arrange(x, Kingdom, Phylum, Class, Order, Family, Genus, Species)
x$otu = factor(x$otu, levels=x$otu)   # set x order
x$num = 1:dim(x)[1]

# top 10 phylum 
per= read.delim("sum_p.txt", sep = "\t", row.names=1, header=T)
mean = rowMeans(per)
per = as.data.frame(mean[order(mean, decreasing = T)])
top_tax=head(rownames(per), n=10)

# Low Abundance
x$tax = x$Phylum

if (length(unique(x$tax)) > length(top_tax)){
  x[!(x$tax %in% top_tax),]$tax = "Low Abundance" # no level can get value
}

# label order
label = unique(x$tax)
label = label[!(label %in% "Low Abundance")] # delete low abundance
x$tax = factor(x$tax, levels = c(label, "Low Abundance"))
#
temp = x[x$tax %in% label, c("tax","num")]
mat_mean = aggregate(temp[,-1], by=temp[1], FUN=mean) # mean


# adjust value for -log10(x$PValue)
if (max(x$neglogp)>20){
  x[x$neglogp>20,]$neglogp = 20
}

# Manhattan plot
FDR = min(x$neglogp[x$level!="NotSig"])
p = ggplot(x, aes(x=num, y=neglogp, color=tax, size=log2CPM, shape=level)) +
  geom_point(alpha=.7) +
  geom_hline(yintercept=FDR, linetype=2, color="lightgrey") +
  scale_shape_manual(values=c(25, 17, 20))+
  scale_size(breaks=c(5, 10, 15)) +
  labs(x="Phylum of Features", y="-log10(P)") +
  labs(title=paste(gsub(".txt", "", basename("stage three vs stage one and two .txt")), sep=" ")) +
  main_theme +
  # theme(legend.position="top") +
  scale_x_continuous(breaks=mat_mean$x, labels=mat_mean$tax) +
  theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1)) + ylim(0,20)
p
# ggsave(file=paste("two-one.manhattan.p_controlID.pdf", sep=""), 
#        p, width = 183, height = 100, useDingbats=F, units = "mm")

