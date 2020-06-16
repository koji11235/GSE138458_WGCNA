setwd("~/Project/20190709_SystemicSclerosis/3.Script")

library(limma)
data <- read.ilmn(""../1.Data/GSE58095_non-normalized.txt"", sep = ""\t"", probeid = ""ID_REF"", expr = ""SAMPLE"")
normdata <- backgroundCorrect(data, method = ""normexp"", offset=20) #offset??????"
#normalizeBetweenArrays?̏o?͂?log2?????ďo?Ă????炵??: https://support.bioconductor.org/p/62765/
"Expdata_quontile <-  normalizeBetweenArrays(normdata, method =""quantile"" )#""scale"")"
"#plotDensities(Expdata_quontile, col=""black"")"
"save(Expdata_quontile,file=""Expdata_quontile.Rdata"")"

boxplot(Expdata_quontile$E)
#Expdata$E
"#PCA <- prcomp(t(Expdata$E), scale=T)"
"#SSC_col=rgb(0.5, 0.5, 0.5, alpha=0.5) #SSC?Q?̐F?̐ݒ?"
"#HC_col=rgb(1, 0, 0, alpha=0.5) #HC?Q?̐F?̐ݒ?"
"#Color = c(rep(SSC_col, 59), rep(HC_col, 43)) #?e?Qn=3?̂Ƃ??̃T???v???̐F?w??"
"#plot(PCA$x[,1], PCA$x[,2], col=Color, pch=16, cex=2.5, xlab=""PCA1"",ylab=""PCA2"")"
#plot(hclust(dist(t(Expdata$E))))







#series matrix???p????WGCNA
#?l?b?g???[?N?̍쐬#
"setwd(""~/Project/20190709_SystemicSclerosis/3.Script"")"
library(WGCNA)
options(stringsAsFactors = FALSE)

"load(""../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/Expdata_quontile.Rdata"")"

#WGCNAdata <- log2(Expdata$E)
datExpr <- Expdata_quontile$E
keep<-rowSums(datExpr>log2(50))>=102
"datExpr<-datExpr[keep,]%>% t()"

"#lnames <- load(""dataInput.Rdata"")"

#soft-thresholding Power ???̐ݒ?
"powers <- c(c(1:10), seq(from = 12, to = 26, by = 2))"

#?l?b?g???[?N?g?|???W?[????
#???ꂼ???̃??ɑ΂???Scale-free fit index??Mean connectivity???Z?o????
#unsigned?i?????ցE?t???ցj?Csigned?i?????ւ̂݁j???ꂼ???Z?o
"sft_unsigned <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = ""unsigned"")"
"sft_signed <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = ""signed"")"

#Scale-free fit index??Mean connectivity???}??
"sizeGrWindow(9, 5)"
"pdf(""Scale_free_fit_index_and_Mean_connectivity_unsigned.pdf"", width = 9, height = 5)"
"par(mfrow = c(1, 2))"
cex1 = 0.9
"plot(sft_unsigned$fitIndices[,1], -sign(sft_unsigned$fitIndices[,3])*sft_unsigned$fitIndices[,2],"
"     xlab = ""Soft Threshold (power)"", ylab = ""Scale Free Topology Model Fit, signed R^2"", type = ""n"","
"     main = paste(""Scale independence (unsigned)""))"
"text(sft_unsigned$fitIndices[,1], -sign(sft_unsigned$fitIndices[,3])*sft_unsigned$fitIndices[,2],"
"     labels = powers, cex = cex1, col = ""red"")"
"abline(h = 0.90, col = ""red"") #R^2?J?b?g?I?t?l0.90?ɐԐ?"

"plot(sft_unsigned$fitIndices[,1], sft_unsigned$fitIndices[,5],"
"     xlab = ""Soft Threshold (power)"", ylab = ""Mean Connectivity"", type = ""n"","
"     main = paste(""Mean connectivity (unsigned)""))"
"text(sft_unsigned$fitIndices[,1], sft_unsigned$fitIndices[,5], labels = powers, cex = cex1, col = ""red"")"
dev.off()

"sizeGrWindow(9, 5)"
"pdf(""Scale_free_fit_index_and_Mean_connectivity_signed.pdf"", width = 9, height = 5)"
"par(mfrow = c(1, 2))"
cex1 = 0.9
"plot(sft_signed$fitIndices[,1], -sign(sft_signed$fitIndices[,3])*sft_signed$fitIndices[,2],"
"     xlab = ""Soft Threshold (power)"", ylab = ""Scale Free Topology Model Fit, signed R^2"", type = ""n"","
"     main = paste(""Scale independence (signed)""))"
"text(sft_signed$fitIndices[,1], -sign(sft_signed$fitIndices[,3])*sft_signed$fitIndices[,2],"
"     labels = powers, cex = cex1, col = ""red"")"
"abline(h = 0.90, col = ""red"") #R^2?J?b?g?I?t?l0.90?ɐԐ?"

"plot(sft_signed$fitIndices[,1], sft_signed$fitIndices[,5],"
"     xlab = ""Soft Threshold (power)"", ylab = ""Mean Connectivity"", type = ""n"","
"     main = paste(""Mean connectivity (signed)""))"
"text(sft_signed$fitIndices[,1], sft_signed$fitIndices[,5], labels = powers, cex = cex1, col = ""red"")"
dev.off()

#Scalefree plot??R^2?????Ȃ????œK?ȃ????T??
#unsigned?̏ꍇ
"pdf(""Scalefree_plot_unsigned.pdf"", width = 12, height = 9)"
"par(mfrow=c(2,2))"
"for(i in c(c(6:10), seq(from = 12, to = 22, by = 2))){"
softPower_unsigned <- i
"adjacency_unsigned <- adjacency(datExpr, power = softPower_unsigned, type = ""unsigned"")"
diag(adjacency_unsigned) <- 0
"connectivity_unsigned <- apply(adjacency_unsigned, 1, sum)"
"scaleFreePlot(connectivity_unsigned, truncated = T)"
"mtext(paste(""beta ="", i), side = 3, line = 0, adj = 0)"
}
dev.off()

#???̑I?ѕ???scale free fit index??0.90?ɒB?????ŏ??l?Cmean connectivity???Ⴗ???Ă͂????Ȃ??Cscale free R^2??1?ɋ߂?
#?? = 14??scale free fit index??0.90?ɒB?????ŏ??l?ł????Cscale free R^2??0.95
#unsigned?̓? = 14???p???ĉ??͂???
softPower_unsigned <- 12
"adjacency_unsigned <- adjacency(datExpr, power = softPower_unsigned, type = ""unsigned"")"

#signed?̏ꍇ
"pdf(""Scalefree_plot_signed.pdf"", width = 12, height = 9)"
"par(mfrow=c(2,2))"
"for(i in c(seq(from = 12, to = 26, by = 2))){"
  softPower_signed <- i
"  adjacency_signed <- adjacency(datExpr, power = softPower_signed, type = ""signed"")"
  diag(adjacency_signed) <- 0
"  connectivity_signed <- apply(adjacency_signed, 1, sum)"
"  scaleFreePlot(connectivity_signed, truncated = T)"
"  mtext(paste(""beta ="", i), side = 3, line = 0, adj = 0)"
}
dev.off()

#?? = 20??scale free fit index??0.90?ɒB?????ŏ??l?ł????Cscale free R^2??0.89
#signed?̓? = 20???p???ĉ??͂???
softPower_signed <- 22
"adjacency_signed <- adjacency(datExpr, power = softPower_signed, type = ""signed"")"

#adjacency?i?אڍs???l?j??topological overlap?iTOM?j?ɕϊ?
#unsigned
"TOM_unsigned <- TOMsimilarity(adjacency_unsigned, TOMType = ""unsigned"")"
dissTOM_unsigned <- 1-TOM_unsigned

#signed
"TOM_signed <- TOMsimilarity(adjacency_signed, TOMType = ""signed"")"
dissTOM_signed <- 1-TOM_signed

#1-TOM???K?w?I?N???X?^?????O
"geneTree_unsigned <- hclust(as.dist(dissTOM_unsigned), method = ""average"")"
"geneTree_signed <- hclust(as.dist(dissTOM_signed), method = ""average"")"

#?f???h???O???????}??
"sizeGrWindow(12, 9)"
"pdf(""Gene_clustering_dissTOM_unsigned.pdf"", width = 12, height = 9)"
"plot(geneTree_unsigned, xlab = """", sub = """", main = ""Gene clustering on TOM-based dissimilarity (unsigned)"","
"     labels = FALSE, hang = 0.04)"
dev.off()
"pdf(""Gene_clustering_dissTOM_signed.pdf"", width = 12, height = 9)"
"plot(geneTree_signed, xlab = """", sub = """", main = ""Gene clustering on TOM-based dissimilarity (signed)"","
"     labels = FALSE, hang = 0.04)"
dev.off()


#???W???[?????o??eigengene?̎Z?o#
#???W???[???T?C?Y?̍ŏ??l???ݒ?
minModuleSize <- 30

#deepSplit??0?`4?????ꂼ?ꎩ?????͂??C???L?̃R?[?h???J???Ԃ????s

#########????????#########
for(i in 0:4){
  
deepSplit <- i #???̐??l??0?`4?ɕς???

##unsigned_merge?Ȃ?
#?t?H???_?쐬?ƈړ?
"ifelse(!dir.exists(paste(""~/Project/20190709_SystemicSclerosis/3.Script/unsigned_"", deepSplit, sep = """")), "
"       dir.create(paste(""~/Project/20190709_SystemicSclerosis/3.Script/unsigned_"", deepSplit, sep = """")), FALSE)"
"setwd(paste(""~/Project/20190709_SystemicSclerosis/3.Script/unsigned_"", deepSplit, sep = """"))"

#Dynamic tree cut???p???ă??W???[?????o
"dynamicMods_unsigned <- cutreeDynamic(dendro = geneTree_unsigned, distM = dissTOM_unsigned,"
"                                      deepSplit = deepSplit, pamRespectsDendro = FALSE, pamStage = FALSE,"
                                      minClusterSize = minModuleSize)
#???W???[??No.???F?ɕϊ????C?e???W???[???̃v???[?u?????ۑ?
dynamicColors_unsigned <- labels2colors(dynamicMods_unsigned)
module_unsigned <- as.data.frame(table(dynamicColors_unsigned))
"write.csv(module_unsigned, paste(""module_unsigned_"",deepSplit,"".csv"", sep = """"))"

#module eigengene?iME?j?̎Z?o
"MEList_unsigned <- moduleEigengenes(datExpr, colors = dynamicColors_unsigned)"
MEs_unsigned <- MEList_unsigned$eigengenes

#ME???p???????W???[???̃N???X?^?????O
MEDiss_unsigned <- 1-cor(MEs_unsigned)
"METree_unsigned <- hclust(as.dist(MEDiss_unsigned), method = ""average"")"
"pdf(paste(""Clustering_ME_unsigned_"", deepSplit, "".pdf"",sep = """"), width = 7, height = 6)"


"plot(METree_unsigned, main = paste(""Clustering of module eigengenes (unsigned, deepSplit "", deepSplit, "")"", sep = """"),"
"     xlab = """", sub = """")"




#ME???߂??i?????????Ă????j???W???[?????T??
#ME?f???h???O?????̃J?b?g?I?t?l???ݒ肷??
MEDissThres <- 0.25 #0.15?`0.25?͈̔͂Ńf???h???O???????݂ēK?X????
#?J?b?g???C????????
"abline(h = MEDissThres, col = ""red"")"
dev.off()

#?f?[?^??rename???ĕۑ?
moduleColors_unsigned <- dynamicColors_unsigned
"colorOrder <- c(""grey"", standardColors(100))"
"moduleLabels_unsigned <- match(moduleColors_unsigned, colorOrder)-1"
"save(MEs_unsigned, moduleLabels_unsigned, moduleColors_unsigned, geneTree_unsigned, "
"     file = paste(""networkConstruction_StepByStep_unsigned_"", deepSplit, "".Rdata"", sep = """"))"

#ME?̕ۑ?
"MEs_unsigned2 <- cbind(rownames(datExpr), MEs_unsigned)"
"write.csv(MEs_unsigned2, paste(""MEs_unsigned_"",deepSplit,"".csv"", sep = """"), row.names = F)"

#intramodular connectivity (kWithin)?̎Z?o
"kwithin_unsigned <- intramodularConnectivity(adjacency_unsigned, moduleColors_unsigned)"

#?v???[?u?Ckwithin?Ccolor?̕\???쐬???C?ۑ?
"binded_unsigned <- cbind(rownames(kwithin_unsigned), kwithin_unsigned, moduleColors_unsigned)"
"colnames(binded_unsigned)[1] <- ""ProbeID"""
"write.csv(binded_unsigned, paste(""Probe_kWitnin_color_unsigned_"", deepSplit, "".csv"", sep=""""), row.names=F)"


##unsigned_merge????
#?J?b?g?I?t?l?ȉ??̃N???X?^?[??merge
"merge_unsigned <- mergeCloseModules(datExpr, dynamicColors_unsigned, cutHeight = MEDissThres, verbose = 3)"
mergedColors_unsigned <- merge_unsigned$colors
#merge???̊e???W???[???̃v???[?u?????ۑ?
merged_module_unsigned <- as.data.frame(table(mergedColors_unsigned))
"write.csv(merged_module_unsigned, paste(""module_unsigned_"",deepSplit,""_merged.csv"", sep = """"))"

#merge???????W???[????ME
mergedMEs_unsigned <- merge_unsigned$newMEs
#?f???h???O??????merge?O???̐F?t?????????W???[?????}??
"pdf(paste(""Gene_dendrogram_and_module_colors_unsigned_"", deepSplit, "".pdf"",sep = """"), width = 9, height = 6)"
"plotDendroAndColors(geneTree_unsigned, cbind(dynamicColors_unsigned, mergedColors_unsigned),"
"                    c(""Dynamic Tree Cut"", ""Merged dynamic""),"
"                    dendroLabels = FALSE, hang = 0.03,"
"                    addGuide = TRUE, guideHang = 0.05,"
"                    main = paste(""Gene dendrogram and module colors (unsigned, deepSplit "", deepSplit, "")"", sep = """"))"
dev.off()

#merge????ME???p???????W???[???̃N???X?^?????O
MEDiss_merged_unsigned <- 1-cor(mergedMEs_unsigned)
"METree_merged_unsigned <- hclust(as.dist(MEDiss_merged_unsigned), method = ""average"")"
"pdf(paste(""Clustering_ME_unsigned_"", deepSplit, ""_merged.pdf"",sep = """"), width = 7, height = 6)"
"plot(METree_merged_unsigned, main = paste(""Clustering of module eigengenes (unsigned, deepSplit "", deepSplit, "", merged)"", sep = """"),"
"     xlab = """", sub = """")"
dev.off()

#?f?[?^??rename???ĕۑ?
moduleColors_unsigned <- mergedColors_unsigned
"colorOrder <- c(""grey"", standardColors(100))"
"moduleLabels_unsigned <- match(moduleColors_unsigned, colorOrder)-1"
"save(mergedMEs_unsigned, moduleLabels_unsigned, moduleColors_unsigned, geneTree_unsigned, "
"     file = paste(""networkConstruction_StepByStep_unsigned_"", deepSplit, ""_merged.Rdata"", sep = """"))"

#ME?̕ۑ?
"mergedMEs_unsigned2 <- cbind(rownames(datExpr), mergedMEs_unsigned)"
"write.csv(mergedMEs_unsigned2, paste(""MEs_unsigned_"",deepSplit,""_merged.csv"", sep = """"), row.names = F)"

#intramodular connectivity (kWithin)?̎Z?o
"kwithin_unsigned <- intramodularConnectivity(adjacency_unsigned, moduleColors_unsigned)"

#?v???[?u?Ckwithin?Ccolor?̕\???쐬???C?ۑ?
"binded_unsigned <- cbind(rownames(kwithin_unsigned), kwithin_unsigned, moduleColors_unsigned)"
"colnames(binded_unsigned)[1] <- ""ProbeID"""
"write.csv(binded_unsigned, paste(""Probe_kWitnin_color_unsigned_"", deepSplit, ""_merged.csv"", sep=""""), row.names=F)"


##signed_merge?Ȃ?
#?t?H???_?̍쐬?ƈړ?
"ifelse(!dir.exists(paste(""~/Project/20190709_SystemicSclerosis/3.Script/signed_"", deepSplit, sep = """")), "
"       dir.create(paste(""~/Project/20190709_SystemicSclerosis/3.Script/signed_"", deepSplit, sep = """")), FALSE)"
"setwd(paste(""~/Project/20190709_SystemicSclerosis/3.Script/signed_"", deepSplit, sep = """"))"

#Dynamic tree cut???p???ă??W???[?????o
"dynamicMods_signed <- cutreeDynamic(dendro = geneTree_signed, distM = dissTOM_signed,"
"                                      deepSplit = deepSplit, pamRespectsDendro = FALSE, pamStage = FALSE,"
                                      minClusterSize = minModuleSize)
#???W???[??No.???F?ɕϊ????C?e???W???[???̃v???[?u?????ۑ?
dynamicColors_signed <- labels2colors(dynamicMods_signed)
module_signed <- as.data.frame(table(dynamicColors_signed))
"write.csv(module_signed, paste(""module_signed_"",deepSplit,"".csv"", sep = """"))"

#module eigengene?iME?j?̎Z?o
"MEList_signed <- moduleEigengenes(datExpr, colors = dynamicColors_signed)"
MEs_signed <- MEList_signed$eigengenes

#ME???p???????W???[???̃N???X?^?????O
MEDiss_signed <- 1-cor(MEs_signed)
"METree_signed <- hclust(as.dist(MEDiss_signed), method = ""average"")"
"pdf(paste(""Clustering_ME_signed_"", deepSplit, "".pdf"",sep = """"), width = 7, height = 6)"
"plot(METree_signed, main = paste(""Clustering of module eigengenes (signed, deepSplit "", deepSplit, "")"", sep = """"),"
"     xlab = """", sub = """")"
#ME???߂??i?????????Ă????j???W???[?????T??
#ME?f???h???O?????̃J?b?g?I?t?l???ݒ肷??
MEDissThres <- 0.25 #0.15?`0.25?͈̔͂Ńf???h???O???????݂ēK?X????
#?J?b?g???C????????
"#abline(h = MEDissThres, col = ""red"")"
dev.off()

#?f?[?^??rename???ĕۑ?
moduleColors_signed <- dynamicColors_signed
"colorOrder <- c(""grey"", standardColors(100))"
"moduleLabels_signed <- match(moduleColors_signed, colorOrder)-1"
"save(MEs_signed, moduleLabels_signed, moduleColors_signed, geneTree_signed, "
"     file = paste(""networkConstruction_StepByStep_signed_"", deepSplit, "".Rdata"", sep = """"))"

#ME?̕ۑ?
"MEs_signed2 <- cbind(rownames(datExpr), MEs_signed)"
"write.csv(MEs_signed2, paste(""MEs_signed_"",deepSplit,"".csv"", sep = """"), row.names = F)"

#intramodular connectivity (kWithin)?̎Z?o
"kwithin_signed <- intramodularConnectivity(adjacency_signed, moduleColors_signed)"

#?v???[?u?Ckwithin?Ccolor?̕\???쐬???C?ۑ?
"binded_signed <- cbind(rownames(kwithin_signed), kwithin_signed, moduleColors_signed)"
"colnames(binded_signed)[1] <- ""ProbeID"""
"write.csv(binded_signed, paste(""Probe_kWitnin_color_signed_"", deepSplit, "".csv"", sep=""""), row.names=F)"


##signed_merge????
#?J?b?g?I?t?l?ȉ??̃N???X?^?[??merge
"merge_signed <- mergeCloseModules(datExpr, dynamicColors_signed, cutHeight = MEDissThres, verbose = 3)"
mergedColors_signed <- merge_signed$colors
#merge???̊e???W???[???̃v???[?u?????ۑ?
merged_module_signed <- as.data.frame(table(mergedColors_signed))
"write.csv(merged_module_signed, paste(""module_signed_"",deepSplit,""_merged.csv"", sep = """"))"

#merge???????W???[????ME
mergedMEs_signed <- merge_signed$newMEs
#?f???h???O??????merge?O???̐F?t?????????W???[?????}??
"pdf(paste(""Gene_dendrogram_and_module_colors_signed_"", deepSplit, "".pdf"",sep = """"), width = 9, height = 6)"
"plotDendroAndColors(geneTree_signed, cbind(dynamicColors_signed, mergedColors_signed),"
"                    c(""Dynamic Tree Cut"", ""Merged dynamic""),"
"                    dendroLabels = FALSE, hang = 0.03,"
"                    addGuide = TRUE, guideHang = 0.05,"
"                    main = paste(""Gene dendrogram and module colors (signed, deepSplit "", deepSplit, "")"", sep = """"))"
dev.off()

#merge????ME???p???????W???[???̃N???X?^?????O
MEDiss_merged_signed <- 1-cor(mergedMEs_signed)
"METree_merged_signed <- hclust(as.dist(MEDiss_merged_signed), method = ""average"")"
"pdf(paste(""Clustering_ME_signed_"", deepSplit, ""_merged.pdf"",sep = """"), width = 7, height = 6)"
"plot(METree_merged_signed, main = paste(""Clustering of module eigengenes (signed, deepSplit "", deepSplit, "", merged)"", sep = """"),"
"     xlab = """", sub = """")"
dev.off()

#?f?[?^??rename???ĕۑ?
moduleColors_signed <- mergedColors_signed
"colorOrder <- c(""grey"", standardColors(100))"
"moduleLabels_signed <- match(moduleColors_signed, colorOrder)-1"
"save(mergedMEs_signed, moduleLabels_signed, moduleColors_signed, geneTree_signed, "
"     file = paste(""networkConstruction_StepByStep_signed_"", deepSplit, ""_merged.Rdata"", sep = """"))"

#ME?̕ۑ?
"mergedMEs_signed2 <- cbind(rownames(datExpr), mergedMEs_signed)"
"write.csv(mergedMEs_signed2, paste(""MEs_signed_"",deepSplit,""_merged.csv"", sep = """"), row.names = F)"

#intramodular connectivity (kWithin)?̎Z?o
"kwithin_signed <- intramodularConnectivity(adjacency_signed, moduleColors_signed)"

#?v???[?u?Ckwithin?Ccolor?̕\???쐬???C?ۑ?
"binded_signed <- cbind(rownames(kwithin_signed), kwithin_signed, moduleColors_signed)"
"colnames(binded_signed)[1] <- ""ProbeID"""
"write.csv(binded_signed, paste(""Probe_kWitnin_color_signed_"", deepSplit, ""_merged.csv"", sep=""""), row.names=F)"
}
#########?????܂?#########