#series matrix???p????WGCNA
#???W???[????trait?̊֘A??#
library(WGCNA)
library(tidyverse)

"setwd(""~/Project/20190709_SystemicSclerosis/3.Script"")"

options(stringsAsFactors = FALSE)
"lnames <- load(""../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/Expdata_quontile.Rdata"")"
datExpr<-Expdata_quontile$E
keep<-rowSums(datExpr>log2(50))>=102
"datExpr<-datExpr[keep,]%>%t()"
"annot <- read.csv(""../1.Data/Annotation/GeneAnnotation_GPL10558.csv"", header = T, stringsAsFactors = F)"

#datTraits?̐??̂??킩?????ȉ??A?T???v???̃A?m?e?[?V?????????p
"datTraits<-read.table(""../1.Data/GSE58095_all_sample_annotatiions.txt"",header = TRUE, skip=4,sep=""\t"",nrows = 102)"
"choosen<-c(#7,12,"
"           13)#,15,16)"
"datTraits<-datTraits[,13] %>% as.data.frame()"
"names(datTraits)<-""Skin_score"""
#deepSplit??0?`4?????ꂼ?????͂??C???L?̃R?[?h???J???Ԃ????s

#########????????#########
for(i in 0:4){
  deepSplit <- i #???̐??l??0?`4?ɕς???

##unsigned_merge?Ȃ?
#?t?H???_?ړ?
"#setwd(paste(""//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_"", deepSplit, sep = """"))"
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
"lnames <- load(paste(""../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/unsigned_"", deepSplit,""/networkConstruction_StepByStep_unsigned_"", deepSplit, "".Rdata"", sep = """"))"



#???`?q?i?v???[?u?j???ƃT???v?????????`
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

#ME???Čv?Z
"MEs0 <- moduleEigengenes(datExpr, moduleColors_unsigned)$eigengenes"
MEs <- orderMEs(MEs0)
"moduleTraitCor <- cor(MEs, datTraits, use = ""p"",method = ""spearman"")"
"moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)"

#???W???[????trait?̊֘A?????q?[?g?}?b?v?Ŏ???
"textMatrix <- paste(signif(moduleTraitCor, 2), ""\n("","
"                    signif(moduleTraitPvalue, 1), "")"", sep = """")"
dim(textMatrix) <- dim(moduleTraitCor)
#dev.new()
"#pdf(paste(""../2.Output/Signed_Unsigned_Deepsplit/unsigned_"", deepSplit,""/Module_trait_relationships_unsigned_"", deepSplit, "".pdf"",sep = """"), width = 10, height = 18)"
"par(mar = c(6, 8.5, 3, 3))"
"labeledHeatmap(Matrix = moduleTraitCor,"
"               xLabels = names(datTraits),"
"               yLabels = names(MEs),"
"               ySymbols = names(MEs),"
"               colorLabels = FALSE,"
"               colors = blueWhiteRed(50),"
"               textMatrix = textMatrix,"
"               setStdMargins = FALSE,"
"               cex.text = 0.5,"
"               zlim = c(-1, 1),"
"               main = paste(""Module-trait relationships (unsigned, deepSplit "", deepSplit, "")"", sep = """"))"
#dev.off()

#???W???[????trait?̑??֌W???y??p-value???\?ɂ??ĕۑ?
colnames(textMatrix) <- names(datTraits)
rownames(textMatrix) <- names(MEs)
"write.csv(textMatrix, paste(""../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/unsigned_"", deepSplit,""/Module_trait_correlation_unsigned_"", deepSplit, "".csv"", sep = """"))"

#???ڂ???trait?i???̏ꍇ??Fibrosis?j?Ɋւ???Gene significance (GS)??Module membership (MM)???Z?o????
skinscore <- as.data.frame(datTraits$Skin_score)
"names(skinscore) <-""skinscore"""
"modNames <- substring(names(MEs), 3)"

"geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = ""p""))"
"MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))"
"names(geneModuleMembership) <- paste(""MM"", modNames, sep = """")"
"names(MMPvalue) <- paste(""p.MM"", modNames, sep = """")"

"geneTraitSignificance <- as.data.frame(cor(datExpr, skinscore, use = ""p""))"
"GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))"
"names(geneTraitSignificance) <- paste(""GS."", names(skinscore), sep = """")"
"names(GSPvalue) <- paste(""p.GS."", names(skinscore), sep = """")"

#Probe ID??gene name?ɕϊ?
probes <- colnames(datExpr)
"probes2annot <- match(probes, annot$ID)"
#?A?m?e?[?V?????????Ă??Ȃ??v???[?u????0?ł??邩?m?F
sum(is.na(probes2annot)) 

#?A?m?e?[?V?????????y??skinscore?Ɋւ???GS?̕\???쐬
"geneInfo0 <- data.frame(Probe_ID = probes,"
"                        geneSymbol = annot$Symbol[probes2annot],"
"                        Entrez_gene_ID = annot$Entrez_Gene_ID[probes2annot],"
"                        moduleColor = moduleColors_unsigned,"
"                        geneTraitSignificance,"
                        GSPvalue)
#skinscore??p?l?ŕ??בւ?
"modOrder <- order(-abs(cor(MEs, skinscore, use = ""p"")))"
#MM???????ǉ?
for(mod in 1:ncol(geneModuleMembership))
{
  oldNames <- names(geneInfo0)
"  geneInfo0 <- data.frame(geneInfo0, geneModuleMembership[,modOrder[mod]],"
"                          MMPvalue[, modOrder[mod]])"
"  names(geneInfo0) <- c(oldNames, paste(""MM."", modNames[modOrder[mod]], sep = """"),"
"                        paste(""p.MM."", modNames[modOrder[mod]], sep = """"))"
}

#IPA annotation???????ǉ?
"#setwd(""E:/IPA_gene_annotation"")"
"#IPAannot <- read.csv(""../1.Data/Annotation/GeneAnnotation_GPL10558.csv"", header = T, stringsAsFactors = F)"
"#IPAannot <- IPAannot[,-c(2:4)]"
"#geneInfo0 <- merge(geneInfo0, IPAannot, by.x = 1, by.y = 1, all.x = T)"

#module color?y??GS?ŕ??בւ?
"geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.skinscore))"
"geneInfo <- geneInfo0[geneOrder, ]"

#?A?m?e?[?V?????????y??skinscore?Ɋւ???GS/MM?̕\???????o??
"#setwd(paste(""//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_"", deepSplit, sep = """"))"
"write.csv(geneInfo, paste(""../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/unsigned_"",deepSplit,""/geneInfo_unsigned_"", deepSplit, "".csv"",sep = """" ), row.names = FALSE)"


##unsigned_merge????
#?t?H???_?ړ?
"#setwd(paste(""//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_"", deepSplit, sep = """"))"
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
"lnames <- load(paste(""../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/unsigned_"",deepSplit,""/networkConstruction_StepByStep_unsigned_"", deepSplit, ""_merged.Rdata"", sep = """"))"

#???`?q?i?v???[?u?j???ƃT???v?????????`
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

#ME???Čv?Z
"MEs0 <- moduleEigengenes(datExpr, moduleColors_unsigned)$eigengenes"
MEs <- orderMEs(MEs0)
"moduleTraitCor <- cor(MEs, datTraits, use = ""p"")"
"moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)"

#???W???[????trait?̊֘A?????q?[?g?}?b?v?Ŏ???
"textMatrix <- paste(signif(moduleTraitCor, 2), ""\n("","
"                    signif(moduleTraitPvalue, 1), "")"", sep = """")"
dim(textMatrix) <- dim(moduleTraitCor)
"#pdf(paste(""../2.Output/Signed_Unsigned_Deepsplit/unsigned_"",deepSplit,""/Module_trait_relationships_unsigned_"", deepSplit, ""_merged.pdf"",sep = """"), width = 10, height = 18)"
"par(mar = c(6, 8.5, 3, 3))"
"labeledHeatmap(Matrix = moduleTraitCor,"
"               xLabels = names(datTraits),"
"               yLabels = names(MEs),"
"               ySymbols = names(MEs),"
"               colorLabels = FALSE,"
"               colors = blueWhiteRed(50),"
"               textMatrix = textMatrix,"
"               setStdMargins = FALSE,"
"               cex.text = 0.5,"
"               zlim = c(-1, 1),"
"               main = paste(""Module-trait relationships (unsigned, deepSplit "", deepSplit, "", merged)"", sep = """"))"
#dev.off()

#???W???[????trait?̑??֌W???y??p-value???\?ɂ??ĕۑ?
colnames(textMatrix) <- names(datTraits)
rownames(textMatrix) <- names(MEs)
"write.csv(textMatrix, paste(""../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/unsigned_"",deepSplit,""/Module_trait_correlation_unsigned_"", deepSplit, ""_merged.csv"", sep = """"))"

#???ڂ???trait?i???̏ꍇ??skinscore?j?Ɋւ???Gene significance (GS)??Module membership (MM)???Z?o????
skinscore <- as.data.frame(datTraits$Skin_score)
"names(skinscore) <-""skinscore"""
"modNames <- substring(names(MEs), 3)"

"geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = ""p""))"
"MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))"
"names(geneModuleMembership) <- paste(""MM"", modNames, sep = """")"
"names(MMPvalue) <- paste(""p.MM"", modNames, sep = """")"

"geneTraitSignificance <- as.data.frame(cor(datExpr, skinscore, use = ""p""))"
"GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))"
"names(geneTraitSignificance) <- paste(""GS."", names(skinscore), sep = """")"
"names(GSPvalue) <- paste(""p.GS."", names(skinscore), sep = """")"

#Probe ID??gene name?ɕϊ?
probes <- colnames(datExpr)
"probes2annot <- match(probes, annot$ID)"
#?A?m?e?[?V?????????Ă??Ȃ??v???[?u????0?ł??邩?m?F
sum(is.na(probes2annot)) 

#?A?m?e?[?V?????????y??skinscore?Ɋւ???GS?̕\???쐬
"geneInfo0 <- data.frame(Probe_ID = probes,"
"                        geneSymbol = annot$Symbol[probes2annot],"
"                        Entrez_gene_ID = annot$Entrez_Gene_ID[probes2annot],"
"                        moduleColor = moduleColors_unsigned,"
"                        geneTraitSignificance,"
                        GSPvalue)
#skinscore??p?l?ŕ??בւ?
"modOrder <- order(-abs(cor(MEs, skinscore, use = ""p"")))"
#MM???????ǉ?
for(mod in 1:ncol(geneModuleMembership))
{
  oldNames <- names(geneInfo0)
"  geneInfo0 <- data.frame(geneInfo0, geneModuleMembership[,modOrder[mod]],"
"                          MMPvalue[, modOrder[mod]])"
"  names(geneInfo0) <- c(oldNames, paste(""MM."", modNames[modOrder[mod]], sep = """"),"
"                        paste(""p.MM."", modNames[modOrder[mod]], sep = """"))"
}

#IPA annotation???????ǉ?
"#setwd(""E:/IPA_gene_annotation"")"
"#IPAannot <- read.csv(""GPL14951_IPAannot.csv"", header = T, stringsAsFactors = F)"
"#IPAannot <- IPAannot[,-c(2:4)]"
"#geneInfo0 <- merge(geneInfo0, IPAannot, by.x = 1, by.y = 1, all.x = T)"

#module color?y??GS?ŕ??בւ?
"geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.skinscore))"
"geneInfo <- geneInfo0[geneOrder, ]"
#?A?m?e?[?V?????????y??skinscore?Ɋւ???GS/MM?̕\???????o??
"#setwd(paste(""../2.Output/Signed_Unsigned_Deepsplit/unsigned_"", deepSplit, sep = """"))"
"write.csv(geneInfo, paste(""../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/unsigned_"", deepSplit,""/geneInfo_unsigned_"", deepSplit, ""_merged.csv"",sep = """" ), row.names = FALSE)"

##signed_merge?Ȃ?
#?t?H???_?ړ?
"#setwd(paste(""//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_"", deepSplit, sep = """"))"
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
"lnames <- load(paste(""../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/signed_"",deepSplit,""/networkConstruction_StepByStep_signed_"", deepSplit, "".Rdata"", sep = """"))"

#???`?q?i?v???[?u?j???ƃT???v?????????`
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

#ME???Čv?Z
"MEs0 <- moduleEigengenes(datExpr, moduleColors_signed)$eigengenes"
MEs <- orderMEs(MEs0)
"moduleTraitCor <- cor(MEs, datTraits, use = ""p"")"
"moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)"

#???W???[????trait?̊֘A?????q?[?g?}?b?v?Ŏ???
"textMatrix <- paste(signif(moduleTraitCor, 2), ""\n("","
"                    signif(moduleTraitPvalue, 1), "")"", sep = """")"
dim(textMatrix) <- dim(moduleTraitCor)
"#pdf(paste(""../2.Output/Signed_Unsigned_Deepsplit/signed_"",deepSplit,""/Module_trait_relationships_signed_"", deepSplit, "".pdf"",sep = """"), width = 10, height = 18)"
"par(mar = c(6, 8.5, 3, 3))"
"labeledHeatmap(Matrix = moduleTraitCor,"
"               xLabels = names(datTraits),"
"               yLabels = names(MEs),"
"               ySymbols = names(MEs),"
"               colorLabels = FALSE,"
"               colors = blueWhiteRed(50),"
"               textMatrix = textMatrix,"
"               setStdMargins = FALSE,"
"               cex.text = 0.5,"
"               zlim = c(-1, 1),"
"               main = paste(""Module-trait relationships (signed, deepSplit "", deepSplit, "")"", sep = """"))"
#dev.off()

#???W???[????trait?̑??֌W???y??p-value???\?ɂ??ĕۑ?
colnames(textMatrix) <- names(datTraits)
rownames(textMatrix) <- names(MEs)
"write.csv(textMatrix, paste(""../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/signed_"",deepSplit,""/Module_trait_correlation_signed_"", deepSplit, "".csv"", sep = """"))"

#???ڂ???trait?i???̏ꍇ??skinscore?j?Ɋւ???Gene significance (GS)??Module membership (MM)???Z?o????
skinscore <- as.data.frame(datTraits$Skin_score)
"names(skinscore) <-""skinscore"""
"modNames <- substring(names(MEs), 3)"

"geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = ""p""))"
"MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))"
"names(geneModuleMembership) <- paste(""MM"", modNames, sep = """")"
"names(MMPvalue) <- paste(""p.MM"", modNames, sep = """")"

"geneTraitSignificance <- as.data.frame(cor(datExpr, skinscore, use = ""p""))"
"GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))"
"names(geneTraitSignificance) <- paste(""GS."", names(skinscore), sep = """")"
"names(GSPvalue) <- paste(""p.GS."", names(skinscore), sep = """")"

#Probe ID??gene name?ɕϊ?
probes <- colnames(datExpr)
"probes2annot <- match(probes, annot$ID)"
#?A?m?e?[?V?????????Ă??Ȃ??v???[?u????0?ł??邩?m?F
sum(is.na(probes2annot)) 

#?A?m?e?[?V?????????y??skinscore?Ɋւ???GS?̕\???쐬
"geneInfo0 <- data.frame(Probe_ID = probes,"
"                        geneSymbol = annot$Symbol[probes2annot],"
"                        Entrez_gene_ID = annot$Entrez_Gene_ID[probes2annot],"
"                        moduleColor = moduleColors_signed,"
"                        geneTraitSignificance,"
                        GSPvalue)
#skinscore??p?l?ŕ??בւ?
"modOrder <- order(-abs(cor(MEs, skinscore, use = ""p"")))"
#MM???????ǉ?
for(mod in 1:ncol(geneModuleMembership))
{
  oldNames <- names(geneInfo0)
"  geneInfo0 <- data.frame(geneInfo0, geneModuleMembership[,modOrder[mod]],"
"                          MMPvalue[, modOrder[mod]])"
"  names(geneInfo0) <- c(oldNames, paste(""MM."", modNames[modOrder[mod]], sep = """"),"
"                        paste(""p.MM."", modNames[modOrder[mod]], sep = """"))"
}

#IPA annotation???????ǉ?
"#setwd(""E:/IPA_gene_annotation"")"
"#IPAannot <- read.csv(""GPL14951_IPAannot.csv"", header = T, stringsAsFactors = F)"
"#IPAannot <- IPAannot[,-c(2:4)]"
"#geneInfo0 <- merge(geneInfo0, IPAannot, by.x = 1, by.y = 1, all.x = T)"

#module color?y??GS?ŕ??בւ?
"geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.skinscore))"
"geneInfo <- geneInfo0[geneOrder, ]"
#?A?m?e?[?V?????????y??skinscore?Ɋւ???GS/MM?̕\???????o??
"#setwd(paste(""//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_"", deepSplit, sep = """"))"
"write.csv(geneInfo, paste(""../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/signed_"",deepSplit,""/geneInfo_signed_"", deepSplit, "".csv"",sep = """" ), row.names = FALSE)"


##signed_merge????
#?t?H???_?ړ?
"#setwd(paste(""//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_"", deepSplit, sep = """"))"
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
"lnames <- load(paste(""../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/signed_"",deepSplit,""/networkConstruction_StepByStep_signed_"", deepSplit, ""_merged.Rdata"", sep = """"))"

#???`?q?i?v???[?u?j???ƃT???v?????????`
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

#ME???Čv?Z
"MEs0 <- moduleEigengenes(datExpr, moduleColors_signed)$eigengenes"
MEs <- orderMEs(MEs0)
"moduleTraitCor <- cor(MEs, datTraits, use = ""p"")"
"moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)"

#???W???[????trait?̊֘A?????q?[?g?}?b?v?Ŏ???
"textMatrix <- paste(signif(moduleTraitCor, 2), ""\n("","
"                    signif(moduleTraitPvalue, 1), "")"", sep = """")"
dim(textMatrix) <- dim(moduleTraitCor)
"#pdf(paste(""../2.Output/Signed_Unsigned_Deepsplit/signed_"",deepSplit,""/Module_trait_relationships_signed_"", deepSplit, ""_merged.pdf"",sep = """"), width = 10, height = 18)"
"par(mar = c(6, 8.5, 3, 3))"
"labeledHeatmap(Matrix = moduleTraitCor,"
"               xLabels = names(datTraits),"
"               yLabels = names(MEs),"
"               ySymbols = names(MEs),"
"               colorLabels = FALSE,"
"               colors = blueWhiteRed(50),"
"               textMatrix = textMatrix,"
"               setStdMargins = FALSE,"
"               cex.text = 0.5,"
"               zlim = c(-1, 1),"
"               main = paste(""Module-trait relationships (signed, deepSplit "", deepSplit, "", merged)"", sep = """"))"
#dev.off()

#???W???[????trait?̑??֌W???y??p-value???\?ɂ??ĕۑ?
colnames(textMatrix) <- names(datTraits)
rownames(textMatrix) <- names(MEs)
"write.csv(textMatrix, paste(""../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/signed_"",deepSplit,""/Module_trait_correlation_signed_"", deepSplit, ""_merged.csv"", sep = """"))"

#???ڂ???trait?i???̏ꍇ??skinscore?j?Ɋւ???Gene significance (GS)??Module membership (MM)???Z?o????
skinscore <- as.data.frame(datTraits$Skin_score)
"names(skinscore) <-""skinscore"""
"modNames <- substring(names(MEs), 3)"

"geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = ""p""))"
"MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))"
"names(geneModuleMembership) <- paste(""MM"", modNames, sep = """")"
"names(MMPvalue) <- paste(""p.MM"", modNames, sep = """")"

"geneTraitSignificance <- as.data.frame(cor(datExpr, skinscore, use = ""p""))"
"GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))"
"names(geneTraitSignificance) <- paste(""GS."", names(skinscore), sep = """")"
"names(GSPvalue) <- paste(""p.GS."", names(skinscore), sep = """")"

#Probe ID??gene name?ɕϊ?
probes <- colnames(datExpr)
"probes2annot <- match(probes, annot$ID)"
#?A?m?e?[?V?????????Ă??Ȃ??v???[?u????0?ł??邩?m?F
sum(is.na(probes2annot)) 

#?A?m?e?[?V?????????y??skinscore?Ɋւ???GS?̕\???쐬
"geneInfo0 <- data.frame(Probe_ID = probes,"
"                        geneSymbol = annot$Symbol[probes2annot],"
"                        Entrez_gene_ID = annot$Entrez_Gene_ID[probes2annot],"
"                        moduleColor = moduleColors_signed,"
"                        geneTraitSignificance,"
                        GSPvalue)
#skinscore??p?l?ŕ??בւ?
"modOrder <- order(-abs(cor(MEs, skinscore, use = ""p"")))"
#MM???????ǉ?
for(mod in 1:ncol(geneModuleMembership))
{
  oldNames <- names(geneInfo0)
"  geneInfo0 <- data.frame(geneInfo0, geneModuleMembership[,modOrder[mod]],"
"                          MMPvalue[, modOrder[mod]])"
"  names(geneInfo0) <- c(oldNames, paste(""MM."", modNames[modOrder[mod]], sep = """"),"
"                        paste(""p.MM."", modNames[modOrder[mod]], sep = """"))"
}

#IPA annotation???????ǉ?
"#setwd(""E:/IPA_gene_annotation"")"
"#IPAannot <- read.csv(""GPL14951_IPAannot.csv"", header = T, stringsAsFactors = F)"
"#IPAannot <- IPAannot[,-c(2:4)]"
"#geneInfo0 <- merge(geneInfo0, IPAannot, by.x = 1, by.y = 1, all.x = T)"

#module color?y??GS?ŕ??בւ?
"geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.skinscore))"
"geneInfo <- geneInfo0[geneOrder, ]"
#?A?m?e?[?V?????????y??skinscore?Ɋւ???GS/MM?̕\???????o??
"#setwd(paste(""//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_"", deepSplit, sep = """"))"
"write.csv(geneInfo, paste(""../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/signed_"",deepSplit,""/geneInfo_signed_"", deepSplit, ""_merged.csv"",sep = """" ), row.names = FALSE)"
}
#########?????܂?#########


#???ڂ???trait?ƍł????ւ??郂?W???[???i???̏ꍇcyan?j?ɂ????āC??GS??MM?̈??`?q(probe)?????肷??
"trait <- ""skinscore"""
"module <- ""tan"""
"column <- match(module, modNames)"
moduleGenes <- moduleColors_unsigned == module

#???ڃ??W???[???ɂ?????GS-MM?v???b?g
"#pdf(paste(""GS_MM_unsigned_"", deepSplit, ""_"", trait, ""_"", module, "".pdf"",sep = """"), width = 7, height = 7)"
"par(mfrow = c(1, 1))"
"verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),"
"                   abs(geneTraitSignificance[moduleGenes, 1]),"
"                       xlab = paste(""Module Membership in"", module, ""module""),"
"                       ylab = paste(""Gene significance for"", trait),"
"                       main = paste(""Module membership vs. gene significance"", ""(unsigned_"", deepSplit, "")\n"", sep = """"),"
"                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)"
#dev.off()