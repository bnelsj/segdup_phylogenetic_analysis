library(ape)
library(phytools)

args = commandArgs(trailingOnly=T)
input = args[1]
output = args[2]
root = args[3]

tree = read.tree(input,)
rooted = root(tree, root, resolve.root=T)
rooted_mid = midpoint.root(rooted)

pdf(output, width=8, height=6)
plot(rooted_mid, cex=0.5)
dev.off()
