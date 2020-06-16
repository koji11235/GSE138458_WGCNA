library(tidyverse)
library(WGCNA)
library(clusterProfiler)
library(reactomePA)
library(org.Hs.eg.db)

setwd("~/Script")

lnames = load("../Data/SLE_GSE138458-01-dataInput.RData")
print(lnames)
###########################################################################

lnames = load(file = "../Data/SLE_GSE138458-01-networkConstruction_signed_4_merged.RData")
print(lnames)

datTraits <- datTraits %>% 
  dplyr::select(-c(case_control, sledai_group, array_ID))

## 3.a Quantifying module–trait associations
# Define numbers of genes and samples
nGenes = ncol(datExpr$E);
nSamples = nrow(datExpr$E);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr$E, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
pdf(file = paste0("Module-trait_relationships.pdf"))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

pvalue_threshold = 0.05

sum(moduleTraitPvalue[,"case_control_levels"] < pvalue_threshold)
significant_module <- moduleTraitPvalue[moduleTraitPvalue[,"case_control_levels"] < pvalue_threshold,] %>% 
  rownames() %>% 
  str_sub(start = 3)



gene_annotation = read_tsv("../Data/GPL10558-50081.txt", comment = "#", col_names = TRUE)

gene_info <- tibble(
  probeID = colnames(datExpr$E),
  module = moduleColors
) %>% 
  left_join(gene_annotation %>% dplyr::select(c(ID, Entrez_Gene_ID, Symbol, Synonyms)), by = c("probeID" = "ID"))


library(clusterProfiler)
library(reactomePA)
library(org.Hs.eg.db)

groupby_enrichment <- gene_info %>% 
  group_by(module) %>% 
  nest() %>% 
  filter(module != "grey", module %in% significant_module) %>% 
  #head(2) %>% 
  mutate(GOenrich = map(data, .f = function(data){enrichGO(data$Entrez_Gene_ID, OrgDb = org.Hs.eg.db, ont = "ALL", readable = TRUE)})#,
         #GOenrich_readable = map(data, .f = function(data){enrichGO(data$Entrez_Gene_ID, OrgDb = org.Hs.eg.db, ont = "ALL", readable = TRUE)})
         )

map2(groupby_enrichment$GOenrich , groupby_enrichment$module,
    .f = function(GOenrich, module){
      bar = barplot(GOenrich)
      ggsave(filename = paste0("../Output/GOenrichment/",module,"_barplot",".pdf"), plot = bar,
             width = 300,height = 200,units = "mm")
      GOenrich@result %>% 
        write_csv(paste0("../Output/GOenrichment/",module,"_result",".csv"))
      }
    )


###########################################################################

setwd("~/Script")

lnames = load("../Data/SLE_GSE138458-01-dataInput.RData")
print(lnames)
datTraits <- datTraits %>% 
  dplyr::select(-c(case_control, sledai_group, array_ID))

gene_annotation = read_tsv("../Data/GPL10558-50081.txt", comment = "#", col_names = TRUE)


constructedNetwork_paths <- list.files(path = "../Data/", pattern = "SLE_GSE138458-01-networkConstruction_.*.RData", full.names = TRUE)
constructedNetwork_paths[[1]]
for (constructedNetwork in constructedNetwork_paths){
  condition = basename(constructedNetwork) %>% str_extract(pattern = "(un)?signed_[0-9](_merged)?")
  print(condition)
  outputdir_base = paste0("../Output/001_NetworkConstruction/", condition, "/")
  if (!dir.exists(outputdir_base)){
    dir.create(outputdir_base)}
  lnames = load(file = constructedNetwork)
  print(lnames)
  
  
  
  ## 3.a Quantifying module–trait associations
  # Define numbers of genes and samples
  nGenes = ncol(datExpr$E);
  nSamples = nrow(datExpr$E);
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(datExpr$E, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, datTraits, use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
  
  sizeGrWindow(10,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  png(file = paste0(outputdir_base,"Module-trait_relationships.png"), width = 720, height = 480)
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
  dev.off()
  
  pvalue_threshold = 0.05
  
  sum(moduleTraitPvalue[,"case_control_levels"] < pvalue_threshold)
  significant_module <- moduleTraitPvalue[moduleTraitPvalue[,"case_control_levels"] < pvalue_threshold,] %>% 
    rownames() %>% 
    str_sub(start = 3)
  
  gene_info <- tibble(
    probeID = colnames(datExpr$E),
    module = moduleColors
  ) %>% 
    left_join(gene_annotation %>% dplyr::select(c(ID, Entrez_Gene_ID, Symbol, Synonyms)), by = c("probeID" = "ID"))
  
  
  # Enrichment 
  
  groupby_enrichment <- gene_info %>% 
    group_by(module) %>% 
    nest() %>% 
    filter(module != "grey", module %in% significant_module) %>% 
    #head(2) %>% 
    mutate(GOenrich = map(data, .f = function(data){enrichGO(data$Entrez_Gene_ID, OrgDb = org.Hs.eg.db, ont = "ALL", readable = TRUE)})#,
           #GOenrich_readable = map(data, .f = function(data){enrichGO(data$Entrez_Gene_ID, OrgDb = org.Hs.eg.db, ont = "ALL", readable = TRUE)})
    )
  
  map2(groupby_enrichment$GOenrich , groupby_enrichment$module,
       .f = function(GOenrich, module){
         bar = barplot(GOenrich)
         ggsave(filename = paste0(outputdir_base, module,"_barplot",".pdf"), plot = bar,
                width = 300,height = 200,units = "mm")
         GOenrich@result %>% 
           write_csv(paste0(outputdir_base, module,"_result",".csv"))
       }
  )
}










