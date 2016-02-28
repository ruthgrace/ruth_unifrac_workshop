#!/usr/bin/env Rscript
options(error=recover)

#this script prints out PDF pcoa plots and distance matrices, given an OTU table, phylogenetic tree, and metadata

source("UniFrac.r")
library(ape)
library(phangorn)
library(vegan)

otu.tab <- read.table("data/tongue_tongue_data/tongue_tongue_data.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

otu.tab <- t(as.matrix(otu.tab))

#sort taxa from most to least abundant
taxaOrder <- rev(order(apply(otu.tab,2,sum)))
otu.tab <- otu.tab[,taxaOrder]

# read and root tree (rooted tree is required)
tree <- read.tree("data/tongue_tongue_data/tongue_tongue_subtree.tre")
tree <- midpoint(tree)

#filter out all OTUs which are < 1% abundance in every sample
samples.one.percent <- 0.01 * apply(otu.tab,1,sum)

otus.greater.than.one.percent <- apply(otu.tab,2,function(x) { return(length(which(x > samples.one.percent))) })

otu.tab <- otu.tab[,which(otus.greater.than.one.percent > 0)]

tree$tip.label <- gsub("'","",tree$tip.label)
absent <- tree$tip.label[!(tree$tip.label %in% colnames(otu.tab))]
if (length(absent) != 0) {
		tree <- drop.tip(tree, absent)
}


groups <- as.factor(substr(rownames(otu.tab),0,4))

#rarefy data for unweighted unifrac
otu.tab.rarefy <- rrarefy(otu.tab, min(apply(otu.tab,1,sum)))

#get rid of zero sum columns
otu.tab.rarefy.sum <- apply(otu.tab.rarefy,2,sum)
otu.tab.rarefy.notzero <- otu.tab.rarefy[, otu.tab.rarefy.sum > 0]
otu.tab.rarefy <- otu.tab.rarefy.notzero

#get rid of zero sum OTUs in tree
tree.rarefy <- tree

tree.rarefy$tip.label <- gsub("'","",tree.rarefy$tip.label)
absent <- tree.rarefy$tip.label[!(tree.rarefy$tip.label %in% colnames(otu.tab.rarefy))]
if (length(absent) != 0) {
		tree.rarefy <- drop.tip(tree.rarefy, absent)
}


#calculate distance matrix
unweighted <- getDistanceMatrix(otu.tab.rarefy,tree.rarefy,method="unweighted",verbose=TRUE)
write.table(unweighted,file="output/tongue_tongue_unweighted_distance_matrix.txt",sep="\t",quote=FALSE)
all_distance_matrices <- getDistanceMatrix(otu.tab,tree,method="all",verbose=TRUE)
weighted <- all_distance_matrices[["weighted"]]
write.table(weighted,file="output/tongue_tongue_weighted_distance_matrix.txt",sep="\t",quote=FALSE)
information <- all_distance_matrices[["information"]]
write.table(information,file="output/tongue_tongue_information_distance_matrix.txt",sep="\t",quote=FALSE)
ratio_no_log <- all_distance_matrices[["ratio_no_log"]]
write.table(ratio_no_log,file="output/tongue_tongue_ratio_no_log_distance_matrix.txt",sep="\t",quote=FALSE)

unweighted.pcoa <- pcoa(unweighted)
weighted.pcoa <- pcoa(weighted)
information.pcoa <- pcoa(information)
ratio_no_log.pcoa <- pcoa(ratio_no_log)


#function to get variance explained for the PCOA component labels
getVarExplained <- function(vector) {
	rawVarEx <- apply(vector,2,function(x) sd(x)*sd(x))
	totalVarExplained <- sum(rawVarEx)
	varEx <- rawVarEx/totalVarExplained
	return(varEx)
}


unweighted.varEx <- getVarExplained(unweighted.pcoa$vectors)
weighted.varEx <- getVarExplained(weighted.pcoa$vectors)
information.varEx <- getVarExplained(information.pcoa$vectors)
ratio_no_log.varEx <- getVarExplained(ratio_no_log.pcoa$vectors)

#get vector version of distance matrices for correlation plots below
unweighted.vector <- unlist(unweighted[lower.tri(unweighted,diag=TRUE)])
weighted.vector <- unlist(weighted[lower.tri(weighted,diag=TRUE)])
information.vector <- unlist(information[lower.tri(information,diag=TRUE)])
ratio_no_log.vector <- unlist(ratio_no_log[lower.tri(ratio_no_log,diag=TRUE)])

pdf("output/tongue_tongue_pcoa_plots.pdf")

#plot pcoa plots
plot(unweighted.pcoa$vectors[,1],unweighted.pcoa$vectors[,2], col=groups,main="Unweighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(unweighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(unweighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
plot(weighted.pcoa$vectors[,1],weighted.pcoa$vectors[,2], col=groups,main="Weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(weighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(weighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend(0.2,0.32,levels(groups),col=palette(),pch=19)
plot(information.pcoa$vectors[,1],information.pcoa$vectors[,2], col=groups,main="Information UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(information.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(information.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
plot(ratio_no_log.pcoa$vectors[,1],ratio_no_log.pcoa$vectors[,2], col=groups,main="Ratio UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(ratio_no_log.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(ratio_no_log.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)

#plot correlation between different UniFrac modes
plot(unweighted.vector,information.vector,main="unweighted vs. information UniFrac")
plot(weighted.vector,information.vector,main="weighted vs. information UniFrac")
plot(unweighted.vector,weighted.vector,main="unweighted vs. weighted UniFrac")
plot(ratio_no_log.vector,information.vector,main="ratio vs. information UniFrac")
plot(ratio_no_log.vector,weighted.vector,main="ratio vs. weighted UniFrac")

dev.off()
