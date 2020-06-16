setwd("~/Script")

library(tidyverse)
library(limma)

read_tsv("../Data/GSE138458_non-normalized.txt")%>% 
  rename_at(vars(starts_with("200")), function(x){paste0("SAMPLE",x)}) %>% 
  write_tsv("../Data/GSE138458_non-normalized_renamed.txt")
  
data <- read.ilmn("../Data/GSE138458_non-normalized_renamed.txt", sep = "\t", probeid ="ID_REF", expr = "SAMPLE")
boxplot(log2(data$E[,1:20]))
length(colnames(datExpr$E)) %/% 50
pdf(file = "../Output/boxplot_rawdata.pdf")
for (i in 1:(length(colnames(data$E)) %/% 50 + 1)){
  #for (i in 1:2){
  print(i)
  start = 50 * (i - 1) + 1 
  end = 50 * i
  if (end > length(colnames(data$E))){
    end = length(colnames(data$E))
  }
  boxplot(log2(data$E[,start:end]))
}
dev.off()

normdata <- backgroundCorrect(data, method = "normexp", offset=20) 
normdata$E
boxplot(log2(normdata$E[,1:10]))
pdf(file = "../Output/boxplot_normdata.pdf")
for (i in 1:(length(colnames(normdata$E)) %/% 50 + 1)){
  #for (i in 1:2){
  print(i)
  start = 50 * (i - 1) + 1 
  end = 50 * i
  if (end > length(colnames(normdata$E))){
    end = length(colnames(normdata$E))
  }
  boxplot(log2(normdata$E[,start:end]))
}
dev.off()

matcor<-cor(log2(normdata$E))
range(as.vector (matcor))

heatmap.2(matcor,trace="none",col=heat.colors(40),Rowv=FALSE,
          #ColSideColors=geno, RowSideColors=treatment,
          cexCol=1,cexRow=1,labCol="")
matcor
rm_matcor <- rowMeans(matcor)
rm_matcor[rm_matcor < 0.94]
hist(rm_matcor, breaks=seq(0.70,1,0.01))
quantile(rm_matcor, 0.1)
keep <- names(rm_matcor[rm_matcor > quantile(rm_matcor, 0.1)])
Not_kept <- names(rm_matcor[rm_matcor <= quantile(rm_matcor, 0.1)])
heatmap.2(matcor[keep, keep],trace="none",col=heat.colors(40),
          #ColSideColors=geno, RowSideColors=treatment,
          cexCol=1,cexRow=1,labCol="")
range(as.vector (matcor[keep, keep]))
datExpr <-  normalizeBetweenArrays(normdata, method ="quantile")
save(datExpr, file="../Data/datExpr.Rdata")
datExpr <- normalizeBetweenArrays(normdata[,keep], method ="quantile")
save(datExpr, file="../Data/datExpr_filtered.Rdata")
#6 samples were determined to be outliers and were dropped prior to normalization:
#200319680117_E 200308380073_J 200308380073_K 200308380073_L 200308380104_E 200667730029_E
#Nor kept
names(rm_matcor[rm_matcor <= quantile(rm_matcor, 0.1)])
# [1] "200319680117_E" "200319680117_G" "200319680117_H" "200319680117_K" "200793290138_C"
# [6] "200793290138_E" "200793290138_K" "200532570064_B" "200319680110_A" "200319680135_D"
# [11] "200319680135_H" "200319680135_J" "200363680045_G" "200363680045_K" "200316700105_L"
# [16] "200319680125_C" "200319680125_G" "200319680125_K" "200308380073_J" "200308380073_K"
# [21] "200308380073_L" "200308380104_E" "200308380104_I" "200319680102_J" "200319680092_K"
# [26] "200667730029_E" "200667730029_F" "200363680050_B" "200793810013_A" "200793810013_G"
# [31] "200319680127_A" "200319680127_C" "200363680043_J" "200319680077_L"

traitData<- read_csv("../Data/Traits.csv")
traitData %>% 
  filter(array_ID  %in% keep)  %>% 
  pull(sledai_group) %>% 
  table()


matcor <- cor(datExpr$E)
range(as.vector (matcor))

heatmap.2(matcor,trace="none",col=heat.colors(40),
          #ColSideColors=geno, RowSideColors=treatment,
          cexCol=1,cexRow=1,labCol="")


lname <-load("../Data/datExpr_filtered.Rdata")
length(colnames(datExpr$E)) %/% 50
pdf(file = "../Output/boxplot_quantile.pdf")
for (i in 1:(length(colnames(datExpr$E)) %/% 50 + 1)){
#for (i in 1:2){
  print(i)
  start = 50 * (i - 1) + 1 
  end = 50 * i
  if (end > length(colnames(datExpr$E))){
    end = length(colnames(datExpr$E))
  }
  boxplot(datExpr$E[,start:end])
}
dev.off()

#boxplot(datExpr$E[,1:50])
