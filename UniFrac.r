# ruth's unifrac method


library(phangorn)
library(ape)
library(zCompositions)

gm_mean = function(x, na.rm=TRUE){
  exp(mean(log(x), na.rm=na.rm) )
}

# gm_vector = function(node_count, other_counts, geometric_mean) {
#   gm_vec <- array(0,length(node_count))
#   # a count of -1 means that the OTU has been amalgamated to other OTUs in the node_count
#   zero_means <- list()
#   counter = 1
#   for (i in 1:length(node_count)) {
#   # CHECK IF THIS IS CORRECT - indexing for other_counts
#     gm_vec[i] <- gm_mean(c(node_count[i],other_counts[i,][which(other_counts[i,]!=-1)]))
#     if (gm_vec[i] == 0) {
#       gm_vec[i] = geometric_mean[i]
#     }
#   }
#   return(gm_vec)
# }

gm_from_props = function(rownum, colnums) {
  otuPropsPerNode.adjustedZeros <- get("otuPropsPerNode.adjustedZeros",envir = .GlobalEnv)
  return(log2(otuPropsPerNode.adjustedZeros[rownum,colnums[1]]) - mean(log2(otuPropsPerNode.adjustedZeros[rownum,colnums])))
}

make_otu_props_positive = function() {
  otuPropsPerNode <- get("otuPropsPerNode",envir = .GlobalEnv)
  otuPropsPerNode.adjustedZeros <- get("otuPropsPerNode.adjustedZeros",envir = .GlobalEnv)
  otuPropsPerNode <- abs(otuPropsPerNode)
  otuPropsPerNode.adjustedZeros <- abs(otuPropsPerNode.adjustedZeros)
  assign("otuPropsPerNode", otuPropsPerNode, envir = .GlobalEnv)
  assign("otuPropsPerNode.adjustedZeros", otuPropsPerNode.adjustedZeros, envir = .GlobalEnv)
}

get_children = function(root) {
  unifrac.tree <- get("unifrac.tree",envir=.GlobalEnv)
  children <- unifrac.tree$edge[which(unifrac.tree$edge[,1]==root),2]
  return(children)
}

get_subtree = function() {
  otuPropsPerNode.adjustedZeros <- get("otuPropsPerNode.adjustedZeros",envir = .GlobalEnv)
  return(colnames(otuPropsPerNode.adjustedZeros)[which(otuPropsPerNode.adjustedZeros[1,] < 0)])
}

build_weights = function(root) {
  # make all weight values positive
  make_otu_props_positive()
  # if root is leaf, then weight has already been calculated and we can skip this
  if (root > length(unifrac.tree$tip.label)) {
    
    children <- get_children(root)
    
    # construct values for left child
    build_weights(children[1])
    leftChildSubtree <- get_subtree()
    
    # construct values for right child
    build_weights(children[2])
    rightChildSubtree <- get_subtree()
    
    otuPropsPerNode <- get("otuPropsPerNode",envir = .GlobalEnv)
    otuPropsPerNode.adjustedZeros <- get("otuPropsPerNode.adjustedZeros",envir = .GlobalEnv)
    weightsPerNode <- get("weightsPerNode",envir = .GlobalEnv)
    unifrac.tree <- get("unifrac.tree",envir=.GlobalEnv)
    unifrac.treeLeaves <- get("unifrac.treeLeaves",envir=.GlobalEnv)
    
    # the negative values in leftChildProportions are for all the nodes in the subtree of the left child
    # the negative values in rightChildProportions are for all the nodes in the subtree of the right child
    
    leafColumns <- colnames(otuPropsPerNode)[which(colnames(otuPropsPerNode) %in% unifrac.treeLeaves)]
    
    # make both left and right child subtrees negative
    otuPropsPerNode[,leftChildSubtree] <- (-1)*abs(otuPropsPerNode[,leftChildSubtree])
    otuPropsPerNode.adjustedZeros[,leftChildSubtree] <- (-1)*abs(otuPropsPerNode.adjustedZeros[,leftChildSubtree])
    
    # get columns for use in weight calculation (leaves not in right or left subtree plus root)
    columnsForWeightCalculation <- leafColumns
    columnsForWeightCalculation <- columnsForWeightCalculation[which(!(columnsForWeightCalculation %in% leftChildSubtree))]
    columnsForWeightCalculation <- columnsForWeightCalculation[which(!(columnsForWeightCalculation %in% rightChildSubtree))]
    columnsForWeightCalculation <- c(root, columnsForWeightCalculation)
    
    # calculate proportion of current node
    otuPropsPerNode[,root] <- abs(otuPropsPerNode[,children[1]]) + abs(otuPropsPerNode[,children[2]])
    otuPropsPerNode.adjustedZeros[,root] <- abs(otuPropsPerNode.adjustedZeros[,children[1]]) + abs(otuPropsPerNode.adjustedZeros[,children[2]])
    assign("otuPropsPerNode", otuPropsPerNode, envir = .GlobalEnv)
    assign("otuPropsPerNode.adjustedZeros", otuPropsPerNode.adjustedZeros, envir = .GlobalEnv)
    
    # calculate weight of current node
    weightsPerNode[,root] <- sapply(c(1:nrow(otuPropsPerNode.adjustedZeros)),function(x) { gm_from_props(x, columnsForWeightCalculation)})
  }
  # make root node values negative
  otuPropsPerNode[,root] <- (-1)*abs(otuPropsPerNode[,root])
  otuPropsPerNode.adjustedZeros[,root] <- (-1)*abs(otuPropsPerNode.adjustedZeros[,root])
  
  assign("otuPropsPerNode", otuPropsPerNode, envir = .GlobalEnv)
  assign("otuPropsPerNode.adjustedZeros", otuPropsPerNode.adjustedZeros, envir = .GlobalEnv)
  assign("weightsPerNode", weightsPerNode, envir = .GlobalEnv)
}

calculateDistanceMatrix <- function(weights, method, otuTable, verbose, pruneTree, normalize) {
  unifrac.tree <- get("unifrac.tree",envir=.GlobalEnv)
  otuPropsPerNode <- get("otuPropsPerNode",envir=.GlobalEnv)
  
  nSamples <- nrow(otuTable)
  distanceMatrix <- matrix(NA,nrow=nSamples,ncol=nSamples)
  rownames(distanceMatrix) <- rownames(otuTable)
  colnames(distanceMatrix) <- rownames(otuTable)
  
  branchLengths <- unifrac.tree$edge.length
  
  weightColnames <- as.numeric(colnames(weights))
  
  weights <- weights[,match(unifrac.tree$edge[,2], weightColnames)]
  
	for (i in 1:nSamples) {
		for (j in i:nSamples) {

				if (method == "weighted" || method == "information" || method == "ratio" || method == "ratio_no_log") {
					# the formula is sum of (proportional branch lengths * | proportional abundance for sample A - proportional abundance for sample B| )
					if (pruneTree==TRUE){
						includeBranchLengths <- which( (otuPropsPerNode[i,] > 0) | (otuPropsPerNode[j,] > 0) )
						if (normalize==TRUE && (method != "ratio" || method != "ratio_no_log")) {
							distance <- sum( branchLengths[includeBranchLengths] * abs(weights[i,includeBranchLengths] - weights[j,includeBranchLengths]) )/sum( branchLengths[includeBranchLengths]* (weights[i,includeBranchLengths] + weights[j,includeBranchLengths]) )
						}
						else {
							distance <- sum( branchLengths[includeBranchLengths] * abs(weights[i,includeBranchLengths] - weights[j,includeBranchLengths]) )/sum( branchLengths[includeBranchLengths])
						}
					}
					else {
						distance <- sum( branchLengths * abs(weights[i,] - weights[j,]) )/sum(branchLengths)
						if (normalize==TRUE) {
							distance <- sum( branchLengths * abs(weights[i,] - weights[j,]) )/sum(branchLengths * (weights[i,] + weights[j,]))
						}
					}

				}
			else {
				if (method!="unweighted") {
					warning(paste("Invalid method",method,", using unweighted Unifrac instead"))
				}
				# the formula is sum of (branch lengths * (1 if one sample has counts and not the other, 0 otherwise) )
				#	i call the (1 if one sample has counts and not the other, 0 otherwise) xorBranchLength
				xorBranchLength <- as.numeric(xor( weights[i,] > 0, weights[j,] > 0))
				if (pruneTree==TRUE) {
					includeBranchLengths <- which( (weights[i,] > 0) | (weights[j,] > 0) )
					distance <- sum( branchLengths[includeBranchLengths] *  xorBranchLength[includeBranchLengths])/sum(branchLengths[includeBranchLengths])
				}
				else {
					distance <- sum( branchLengths *  xorBranchLength)/sum(branchLengths)
				}

			}
			distanceMatrix[i,j] <- distance
			distanceMatrix[j,i] <- distance

		}
	}

	if(verbose) {	print("done")	}

	return(distanceMatrix)

}

#valid methods are unweighted, weighted, information, and ratio. Any other method will result in a warning and the unweighted analysis
#pruneTree option prunes the tree for each comparison to exclude branch lenxgths not present in both samples
#normalize divides the value at each node by sum of weights to guarantee output between 0 and 1 (breaks the triangle inequality)

#otuTable must have samples as rows, OTUs as columns
#tree must be phylo tree object from ape package (can use read.tree method to read in a newick tree to this format)

getDistanceMatrix <- function(otuTable,tree,method="weighted",verbose=FALSE,pruneTree=FALSE,normalize=TRUE)  {

	if (length(which(is.na(otuTable))) > 0) {
		stop("OTU count table has NA")
	}

	if (!is.rooted(tree)) {
		tree <- midpoint(tree)
		if(verbose) { print("Rooting tree by midpoint") }
	}

	if (attributes(tree)$order!="postorder") {
		tree <- reorder(tree,order="postorder")
		if (verbose) { print("Reordering tree as postorder for distance calculation algorithm") }
	}
  
  #make globally available copy of tree
  assign("unifrac.tree", tree, envir = .GlobalEnv)
  unifrac.tree <- get("unifrac.tree", envir = .GlobalEnv)
  
  #make globally available copy of leaves
  assign("unifrac.treeLeaves", c(1:length(tree$tip.label)), envir = .GlobalEnv)
  unifrac.treeLeaves <- get("unifrac.treeLeaves", envir = .GlobalEnv)
  
	# get proportions
	readsPerSample <- apply(otuTable,1,sum)
	otu.prop <- otuTable/readsPerSample
	otu.prop <- as.matrix(otu.prop)
	rownames(otu.prop) <- rownames(otuTable)
	colnames(otu.prop) <- colnames(otuTable)
	if(verbose) {	print("calculated proportional abundance")	}

	# add priors to zeros based on bayesian approach
	otuTable.adjustedZeros <- cmultRepl(otuTable, method="CZM", output="counts")
  # make any negative numbers as close to zero as possible - this is probably due to a precision error.
  otuTable.adjustedZeros <- apply(otuTable.adjustedZeros,2,function(x) { x[which(x < 0)] <- .Machine$double.eps; return(x) })
  readsPerSample <- apply(otuTable.adjustedZeros,1,sum)
  otu.prop.adjustedZeros <- otuTable.adjustedZeros/readsPerSample
  otu.prop.adjustedZeros <- as.matrix(otu.prop.adjustedZeros)
	rownames(otu.prop.adjustedZeros) <- rownames(otuTable)
	colnames(otu.prop.adjustedZeros) <- colnames(otuTable)
  
	##get cumulative proportional abundance for the nodes (nodes are ordered same as in the phylo tree representation)

  otuPropsPerNode <- matrix(NA, ncol=length(tree$edge.length)+1, nrow=length(rownames(otuTable)))
  otuPropsPerNode.adjustedZeros <- matrix(NA, ncol=length(tree$edge.length)+1, nrow=length(rownames(otuTable)))
  weightsPerNode <- matrix(NA, ncol=length(tree$edge.length)+1, nrow=length(rownames(otuTable)))
  
	#each row is a sample
  rownames(otuPropsPerNode) <- rownames(otuTable)
  rownames(otuPropsPerNode.adjustedZeros) <- rownames(otuTable)
  rownames(weightsPerNode) <- rownames(otuTable)
  #each column is a node in the tree
  colnames(otuPropsPerNode) <- c(1:(length(tree$edge.length)+1))
  colnames(otuPropsPerNode.adjustedZeros) <- c(1:(length(tree$edge.length)+1))
  colnames(weightsPerNode) <- c(1:(length(tree$edge.length)+1))
  
  leafNodes <- c(1:length(tree$tip.label))
  leafOrder <- match(tree$tip.label,colnames(otu.prop))
  
  otuPropsPerNode[,leafNodes] <- otu.prop[,leafOrder]
  otuPropsPerNode.adjustedZeros[,leafNodes] <- otu.prop.adjustedZeros[,leafOrder]
  weightsPerNode[,leafNodes] <- log2(otu.prop.adjustedZeros[,leafOrder])
  weightsPerNode[,leafNodes] <- t(apply(weightsPerNode[,leafNodes],1,function(x) { return(x - mean(x))}))
  
  # the tree is in postorder, so the last edges belong to the root
  root <- tree$edge[nrow(tree$edge),1]

  if(verbose) {	print("calculating weights...")	}
  
  assign("otuPropsPerNode", otuPropsPerNode, envir = .GlobalEnv)
  assign("otuPropsPerNode.adjustedZeros", otuPropsPerNode.adjustedZeros, envir = .GlobalEnv)
  assign("weightsPerNode", weightsPerNode, envir = .GlobalEnv)
  
  build_weights(root)
  
  otuPropsPerNode <- get("otuPropsPerNode",envir = .GlobalEnv)
  otuPropsPerNode.adjustedZeros <- get("otuPropsPerNode.adjustedZeros",envir = .GlobalEnv)
  weightsPerNode <- get("weightsPerNode",envir = .GlobalEnv)
    
  # build_weights makes everything negative, make everything positive again
  otuPropsPerNode[,c(1:ncol(otuPropsPerNode))] <- abs(otuPropsPerNode[,c(1:ncol(otuPropsPerNode))])
  otuPropsPerNode.adjustedZeros[,c(1:ncol(otuPropsPerNode.adjustedZeros))] <- abs(otuPropsPerNode.adjustedZeros[,c(1:ncol(otuPropsPerNode.adjustedZeros))])
  
  assign("otuPropsPerNode", otuPropsPerNode, envir = .GlobalEnv)
  assign("otuPropsPerNode.adjustedZeros", otuPropsPerNode.adjustedZeros, envir = .GlobalEnv)
  assign("weightsPerNode", weightsPerNode, envir = .GlobalEnv)
  
  if(verbose) {	print("calculating pairwise distances...")	}
  
  #convert table according to weight
  if (method == "all") {
    returnList <- list()
    # unweighted
    print("calculating unweighted distance matrix")
    weights <- otuPropsPerNode
    returnList[["unweighted"]] <- calculateDistanceMatrix(weights, "unweighted", otuTable, verbose, pruneTree, normalize)
    # weighted
    print("calculating weighted distance matrix")
    weights <- otuPropsPerNode
    returnList[["weighted"]] <- calculateDistanceMatrix(weights, "weighted", otuTable, verbose, pruneTree, normalize)
    # information
    print("calculating information distance matrix")
    weights <- otuPropsPerNode.adjustedZeros*log2(otuPropsPerNode.adjustedZeros)
    returnList[["information"]] <- calculateDistanceMatrix(weights, "information", otuTable, verbose, pruneTree, normalize)
    # ratio
    print("calculating ratio distance matrix")
    weights <- weightsPerNode
    returnList[["ratio"]] <- calculateDistanceMatrix(weights, "ratio", otuTable, verbose, pruneTree, normalize)
    # ratio no log
    print("calculating no log ratio distance matrix")
    weights <- weightsPerNode
    weights <- 2^weights
    returnList[["ratio_no_log"]] <- calculateDistanceMatrix(weights, "ratio_no_log", otuTable, verbose, pruneTree, normalize)
    return(returnList)
  } else if (method=="information") {
    weights <- otuPropsPerNode.adjustedZeros*log2(otuPropsPerNode.adjustedZeros)
    return(calculateDistanceMatrix(weights, method, otuTable, verbose, pruneTree, normalize))
  } else if (method == "ratio_no_log") {
    weights <- weightsPerNode
    weights <- 2^weights
    return(calculateDistanceMatrix(weights, method, otuTable, verbose, pruneTree, normalize))
  } else if (method == "ratio") {
    weights <- weightsPerNode
    return(calculateDistanceMatrix(weights, method, otuTable, verbose, pruneTree, normalize))
  }
  else {
    weights <- otuPropsPerNode
    return(calculateDistanceMatrix(weights, method, otuTable, verbose, pruneTree, normalize))
  }
}
