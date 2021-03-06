library(tidyverse)
library(WGCNA)

setwd("~/Script")

lnames = load("../Data/SLE_GSE138458-01-dataInput.RData")

powers = c(c(1:10), seq(from = 12, to=20, by=2))
#Dockerのメモリを2GB->20GBに変更
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
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

lnames = load("../Data/SLE_GSE138458-01-dataInput.RData")
#[1] "datExpr"   "datTraits"
softPower = 14
adjacency = adjacency(datExpr, power = softPower)


TOM = TOMsimilarity(adjacency,TOMType = "signed");
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average");


for (deepSplit in c(0:4)){
  minModuleSize = 30;
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = deepSplit, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize);
  print(table(dynamicMods))
  
  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  # Plot the dendrogram and colors underneath
  sizeGrWindow(8,6)
  pdf(paste0("../Output/001_NetworkConstruction/NetworkConstruction","_deepsplit",deepSplit, "_", "signed",".pdf"))
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  #dev.off()
  
  # Calculate eigengenes
  MEList = moduleEigengenes(datExpr, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs);
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");
  # Plot the result
  #sizeGrWindow(7, 6)
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  moduleColors = dynamicColors
  colorOrder = c("grey", standardColors(50))
  moduleLabels = match(moduleColors, colorOrder)-1
  
  # Save module colors and labels for use in subsequent parts
  save(MEs, moduleLabels, moduleColors, geneTree, file = paste0("../Data/ConstructedNetwork/SLE_GSE138458-01-networkConstruction", "_", "signed", "_",deepSplit,".RData"))
  
  
  # Merge Modules 
  
  MEDissThres = 0.15
  # Plot the cut line into the dendrogram
  abline(h=MEDissThres, col = "red")
  # Call an automatic merging function
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors;
  print(table(mergedColors))
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs;
  
  #sizeGrWindow(12, 9)
  #pdf(file = paste0("../Output/geneDendro-", deepSplit,".pdf"), wi = 9, he = 6)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
  # Rename to moduleColors
  moduleColors = mergedColors
  # Construct numerical labels corresponding to the colors
  colorOrder = c("grey", standardColors(50));
  moduleLabels = match(moduleColors, colorOrder)-1;
  MEs = mergedMEs;
  # Save module colors and labels for use in subsequent parts
  save(MEs, moduleLabels, moduleColors, geneTree, file = paste0("../Data/ConstructedNetwork/SLE_GSE138458-01-networkConstruction", "_", "signed", "_",deepSplit, "_", "merged",".RData"))
}


###########################################################################

lnames = load("../Data/SLE_GSE138458-01-dataInput.RData")
#[1] "datExpr"   "datTraits"
softPower = 14
adjacency = adjacency(datExpr, power = softPower)


TOM = TOMsimilarity(adjacency,TOMType = "unsigned");
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average");

for (deepSplit in c(0:4)){
  minModuleSize = 30;
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = deepSplit, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize);
  print(table(dynamicMods))
  
  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  # Plot the dendrogram and colors underneath
  sizeGrWindow(8,6)
  pdf(paste0("../Output/NetworkConstruction","_deepsplit",deepSplit, "_", "unsigned",".pdf"))
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  #dev.off()
  
  # Calculate eigengenes
  MEList = moduleEigengenes(datExpr, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1-cor(MEs);
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");
  # Plot the result
  #sizeGrWindow(7, 6)
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  # Save module colors and labels for use in subsequent parts
  moduleColors = dynamicColors
  colorOrder = c("grey", standardColors(50))
  moduleLabels = match(moduleColors, colorOrder)-1
  
  save(MEs, moduleLabels, moduleColors, geneTree, file = paste0("../Data/ConstructedNetwork/SLE_GSE138458-01-networkConstruction", "_", "unsigned", "_",deepSplit,".RData"))
  
  
  # Merge Modules 
  
  MEDissThres = 0.15
  # Plot the cut line into the dendrogram
  abline(h=MEDissThres, col = "red")
  # Call an automatic merging function
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors;
  print(table(mergedColors))
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs;
  
  #sizeGrWindow(12, 9)
  #pdf(file = paste0("../Output/geneDendro-", deepSplit,".pdf"), wi = 9, he = 6)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
  # Rename to moduleColors
  moduleColors = mergedColors
  # Construct numerical labels corresponding to the colors
  colorOrder = c("grey", standardColors(50));
  moduleLabels = match(moduleColors, colorOrder)-1;
  MEs = mergedMEs;
  # Save module colors and labels for use in subsequent parts
  save(MEs, moduleLabels, moduleColors, geneTree, file = paste0("../Data/ConstructedNetwork/SLE_GSE138458-01-networkConstruction", "_", "unsigned", "_",deepSplit, "_", "merged",".RData"))
}






