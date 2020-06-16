library(tidyverse)
library(WGCNA)

setwd("~/20200421_SLE_WGCNA/Script")

datExpr0 <- read_tsv("../Data/GSE138458_series_matrix.txt",comment = "!")

datExpr0 <- datExpr0 %>% 
  as.data.frame() %>% 
  column_to_rownames("ID_REF") %>% t()

gsg =   goodSamplesGenes(datExpr0,verbose = 3)
gsg$allOK
#[1] FALSE

#datExpr0 %>% 
#  dplyr::select(colnames(datExpr0)[1:ncol(datExpr0)][gsg$goodGenes]) %>% 
#  filter(gsg$goodSamples)

datExpr0 <- datExpr0[gsg$goodSamples,gsg$goodGenes]
#dim(datExpr0)
#[1]   330 47323

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

abline(h = 100, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 100, minSize = 1)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


## Create Trait data #######
traitData <-read_tsv("../Data/GSE138458_series_matrix.txt",skip = 38,n_max = 27,col_names = FALSE) #%>% t()
colnames(traitData)<- traitData[27,] %>% as.character()
traitData <- traitData[1:3,] %>%
  rename(Characteristic="ID_REF") %>%
  mutate(Characteristic=c("case_control", "study_group","sledai_group"))%>%
  tidyr::pivot_longer(
    cols = -Characteristic,
    names_to = "ID_REF",
    values_to = "value",
    names_prefix = ""
  ) %>%
  pivot_wider(
    names_from=Characteristic,
    values_from=value
  ) %>% 
  dplyr::select(-study_group) %>% 
  mutate(sledai_group_levels=sledai_group_map[sledai_group],
         case_control_levels=case_control_map[case_control])

traitData %>%
  write_csv("../Data/Traits.csv")
#######


## 1.c Loading clinical trait data
traitData<- read_csv("../Data/Traits.csv",)
dim(traitData)
names(traitData)

datTraits <- traitData %>% 
  filter(ID_REF %in% rownames(datExpr)) %>% 
  as.data.frame() %>% 
  column_to_rownames("ID_REF")

collectGarbage();


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
#traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
# plotDendroAndColors(sampleTree2, traitColors,
#                     groupLabels = names(datTraits),
#                     main = "Sample dendrogram and trait heatmap")
plot(sampleTree2, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

save(datExpr, datTraits, file = "../Output/SLE_GSE138458-01-dataInput.RData")





