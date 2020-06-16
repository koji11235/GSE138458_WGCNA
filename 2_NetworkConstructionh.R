library(tidyverse)
library(WGCNA)

setwd("~/Script")

lnames = load("../Output/SLE_GSE138458-01-dataInput.RData")

powers = c(c(1:10), seq(from = 12, to=20, by=2))
#Dockerのメモリを2GB->20GBに変更
sft = pickSoftThreshold(datExpr$E, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 14


###########################################################################
# 2.c.2 Block-wise network construction and module detection
for (deepSplit in c(0:4)){
  bwnet = blockwiseModules(datExpr$E, maxBlockSize = 16000,
                           power = softPower, TOMType = "unsigned", minModuleSize = 50,
                           reassignThreshold = 0, mergeCutHeight = 0.25,
                           numericLabels = TRUE,
                           saveTOMs = TRUE,
                           deepSplit = deepSplit,
                           saveTOMFileBase = paste0("../Output/network_power14/SLE_GSE138458_TOM-blockwise","_softpower_",softPower,"_deepsplit_",deepSplit),
                           verbose = 3)
  print(table(bwnet$colors))
  
  sizeGrWindow(12, 9)
  # Convert labels to colors for plotting
  mergedColors = labels2colors(bwnet$colors)
  # Plot the dendrogram and the module colors underneath
  pdf(paste0("../Output/network_power14/DendroAndColors","_deepsplit_",deepSplit,".pdf"))
  plotDendroAndColors(bwnet$dendrograms[[1]], mergedColors[bwnet$blockGenes[[1]]],
                      "Module colors",
                      main="Cluster Dendrogram: Block 1",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
  plotDendroAndColors(bwnet$dendrograms[[2]], mergedColors[bwnet$blockGenes[[2]]],
                      "Module colors",
                      main="Cluster Dendrogram: Block 2",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  plotDendroAndColors(bwnet$dendrograms[[3]], mergedColors[bwnet$blockGenes[[3]]],
                      "Module colors",
                      main="Cluster Dendrogram: Block 3",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
  moduleLabels = bwnet$colors
  moduleColors = labels2colors(bwnet$colors)
  MEs = bwnet$MEs;
  geneTrees = bwnet$dendrograms;
  save(MEs, moduleLabels, moduleColors, geneTrees,
       file = paste0("../Output/network_power14/SLE_GSE138458-02-networkConstruction_deepsplit_",deepSplit,".RData"))
}
###########################################################################
##  adjacencyの計算がメモリ的に無理
adjacency = adjacency(datExpr, power = softPower)


## 2.b.3 Topological Overlap Matrix (TOM)
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


## 2.b.4 Clustering using TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


lnames = load("../Output/network_power14/SLE_GSE138458-02-networkConstruction_deepsplit_0.RData")
geneTree=geneTrees[[3]]
lname =load("../Output/network_power14/SLE_GSE138458_TOM-blockwise_softpower_14_deepsplit_0-block.3.RData")
dissTOM = 1-TOM
a
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 4, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


