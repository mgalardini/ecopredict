library(ape)
library(geiger)
library(nlme)
library(phytools)

args = commandArgs(trailingOnly=TRUE)

# tree
# phenotypic data
# dists
# pic

pdata <- read.csv(file.path(args[2]), row.names = 1, sep = '\t')
tree <- read.tree(file.path(args[1]))
tree <- midpoint.root(tree)
tree <- multi2di(tree)
ldist <- cophenetic.phylo(tree)

ndist <- dist.nodes(tree)
tree <- compute.brlen(tree, tree$edge.length)
npdist <- dist.nodes(tree)

vpic <- list()
for (i in 1:length(pdata)) {
  cond <- colnames(pdata)[i]
  vec <- pdata[, cond]
  names(vec) <- rownames(pdata)
  vpic[[i]] <- pic(vec, tree, scaled=T)
}
vpic <- data.frame(vpic)
colnames(vpic) <- colnames(pdata)
npdist <- npdist[rownames(vpic), rownames(vpic)]

write.csv(npdist, file.path(args[3]))
write.csv(vpic, file.path(args[4]))
