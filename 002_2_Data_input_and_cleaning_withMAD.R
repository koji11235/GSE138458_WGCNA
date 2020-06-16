library(WGCNA)
library(tidyverse)
library(stats)

#############################################################
## Gene Filtering ##
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/faq.html
#Should I filter probesets or genes?
# Probesets or genes may be filtered by mean expression or variance 
# or their robust analogs such as median and median absolute deviation, MAD) 
# since low-expressed or non-varying genes usually represent noise. 
# Whether it is better to filter by mean expression or variance is a matter of debate; 
# both have advantages and disadvantages, but more importantly, 
# they tend to filter out similar sets of genes since mean and variance are usually related.

lname <- load("../Data/datExpr_filtered.Rdata")

#############################################################
# datExpr$E
# mean(datExpr$E)
# rowMeans(datExpr$E)
# #keep<-rowSums(datExpr>log2(50))>=102
# keep<-rowMeans(datExpr$E) >=mean(datExpr$E)
# sum(keep)
# datExpr$E <- datExpr$E[keep,] %>% t()
# dim(datExpr)
#############################################################
# filtering by mad 

dim(datExpr$E)
gene_mad = apply(datExpr$E, MARGIN = 1, FUN = mad)
hist(gene_mad)
gene_mad_rank = rank(-gene_mad, ties.method = "first")
keep = gene_mad_rank < 6000
hist(gene_mad[keep])

datExpr$E <- datExpr$E[keep,] %>% t()
dim(datExpr)



gsg =   goodSamplesGenes(datExpr$E,verbose = 3)
gsg$allOK
#[1] TRUE
datExpr0 <- datExpr
save(datExpr0, file="../Data/datExpr0.Rdata")


#############################################################
lname <- load("../Data/datExpr0.Rdata")

#datExpr0 %>% 
#  dplyr::select(colnames(datExpr0)[1:ncol(datExpr0)][gsg$goodGenes]) %>% 
#  filter(gsg$goodSamples)

#datExpr0 <- datExpr0[gsg$goodSamples,gsg$goodGenes]
#dim(datExpr0)
#[1]   330 47323

sampleTree = hclust(dist(datExpr0$E), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
cutHeight = 83
abline(h = cutHeight, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = cutHeight, minSize = 1)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0$E[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


## Create Trait data #######
traitData <-read_tsv("../Data/GSE138458_series_matrix.txt",skip = 38,n_max = 28,col_names = FALSE) #%>% t()
colnames(traitData)<- traitData[27,] %>% as.character()
sledai_group_map <- c(0,1,2)
names(sledai_group_map) <- c("sledai group: no val", "sledai group: High", "sledai group: Low")
case_control_map <- c(0, 1)
names(case_control_map) <- c("case/control: Control", "case/control: SLE Case")

traitData <- traitData[c(1:3, 11),] %>%
  dplyr::rename(Characteristic="ID_REF") %>%
  mutate(Characteristic=c("case_control", "study_group","sledai_group", "array_ID"))%>%
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
  filter(array_ID %in% rownames(datExpr)) %>% 
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

save(datExpr, datTraits, file = "../Data/SLE_GSE138458-01-dataInput.RData")

#############################################################




