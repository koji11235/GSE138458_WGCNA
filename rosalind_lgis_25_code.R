001_choose_deepsplit_signed.R





library(tidyverse)

setwd("~/Project/20190709_SystemicSclerosis/3.Script")

output_file_base_path<-"../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/"





####################################################
#gene number of each grey module

col_types=list("_","c","d")
result_df <- data.frame(matrix(rep(NA, 4), nrow=1))[numeric(0), ]
colnames(result_df) <- c("deepsplit","sign_unsigned", "merged_unmerged","freq_grey_module")

condition_deepsplit<-0:4
condition_signed_unsigned<-c("signed","unsigned")
condition_merged<-c("","_merged")
condition<-list(rep(condition_deepsplit,each=4),rep(condition_signed_unsigned,each=2,5),rep(condition_merged,10))
for(i in 1:length(condition[[1]])){
  deepSplit <- condition[[1]][[i]]
  signed_unsigned<-condition[[2]][[i]]
  merged_unmerged<-condition[[3]][[i]]
  
  dir_path<-paste(output_file_base_path, signed_unsigned, "_", deepSplit,"/", sep = "")
  
  module_data<-read_csv(paste0(dir_path,"module_",signed_unsigned,"_",deepSplit,merged_unmerged,".csv"),col_types = col_types)
  
  varname <-colnames(module_data)[1]
  freq_grey_module<-module_data %>% filter_(paste(varname,"==","'grey'")) %>% .$Freq #%>% .[[1]]
  result_df<-rbind(result_df,data.frame(deepSplit=deepSplit,signed_unsigned=signed_unsigned,merged_unmerged=str_sub(merged_unmerged,2),freq_grey_module=freq_grey_module))
}
result_df
result_df %>% 
  ggplot(aes(x=deepSplit,y=freq_grey_module))+
  geom_line(aes(color=interaction(signed_unsigned,merged_unmerged,sep = "_")))+
  facet_grid(signed_unsigned~merged_unmerged)+
  guides(color=FALSE)
ggsave(file = "freq_grey_by_deepsplit.png",dpi=300)

#########################################
#number of modules of each condition
col_types=list("_","c","d")
result_df <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
colnames(result_df) <- c("deepSplit","signed_unsigned", "merged_unmerged","color","freq")

condition_deepsplit<-0:4
condition_signed_unsigned<-c("signed","unsigned")
condition_merged<-c("","_merged")
condition<-list(rep(condition_deepsplit,each=4),rep(condition_signed_unsigned,each=2,5),rep(condition_merged,10))

for(i in 1:length(condition[[1]])){
  deepSplit <- condition[[1]][[i]]
  signed_unsigned<-condition[[2]][[i]]
  merged_unmerged<-condition[[3]][[i]]
  
  dir_path<-paste(output_file_base_path, signed_unsigned, "_", deepSplit,"/", sep = "")
  
  module_data<-read_csv(paste0(dir_path,"module_",signed_unsigned,"_",deepSplit,merged_unmerged,".csv"),col_types = col_types)
  colnames(module_data)<-c("color","freq")
  result_df<-rbind(result_df,
                   data.frame(deepSplit=deepSplit,
                              signed_unsigned=signed_unsigned,
                              merged_unmerged=str_sub(merged_unmerged,2),
                              color=module_data[,1], 
                              freq=module_data[,2]))
}
result_df
result_df %>% 
  group_by(signed_unsigned,deepSplit,merged_unmerged) %>% 
  count() %>% 
  ggplot(aes(x=deepSplit,y=n))+
  geom_line(aes(color=interaction(signed_unsigned,merged_unmerged,sep = "_")))+
  facet_grid(signed_unsigned~merged_unmerged)+
  guides(color=FALSE)+
  ggtitle("Number of modules")
ggsave(file = "freq_grey_by_deepsplit.png",dpi=300)




###################################################
#Num of gene in each module
result_df %>% 
  group_by(signed_unsigned,deepSplit,merged_unmerged) %>% 
  summarize(genenum_mean=mean(freq,na.rm = T),
            genenum_median=median(freq,na.rm = T),
            sd=sd(freq)) %>% 
  filter(signed_unsigned=="signed",merged_unmerged=="") %>% 
  ggplot(aes(x=deepSplit))+
    geom_line(aes(y=genenum_median))


result_df %>% 
  filter(signed_unsigned=="signed",merged_unmerged=="") %>% 
  ggplot(aes(x=as.factor(deepSplit)))+
  geom_violin(trim=T,aes(y=freq))
    


result_df %>% 
  filter(signed_unsigned=="signed",merged_unmerged=="",color!="grey") %>% 
  ggplot(aes(x=freq,fill=color))+
  geom_histogram()+
  facet_grid(deepSplit~.)+
  guides(fill=F)+xlim(0,1000)


result_df %>% 
  filter(signed_unsigned=="signed",merged_unmerged=="") %>% 
  ggplot(aes(x=reorder(color,freq),y=log10(freq),fill=color))+
  geom_bar(stat = "identity")+
  facet_grid(.~deepSplit,scales='free',  space = "free")+
  guides(fill=F)+
  coord_flip()






001_DiomensionReductionPlot.R





library(tidyverse)
library(stringr)
library(Rtsne)

setwd( "C:/Users/ywt2100/Documents/Project/20190709_SystemicSclerosis/3.Script")
lname<-load("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/datExpr_filtered.Rdata")
datExpr
#datExpr<-Expdata_quontile$E
#keep<-rowSums(datExpr>log2(50))>=102
#datExpr<-datExpr[keep,]%>%t()
#save(datExpr,file = "../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/datExpr_filtered.Rdata")

sample_annotation<-read_tsv("../1.Data/GSE58095_all_sample_annotatiions.txt",skip = 4,n_max = 102)
sample_annotation$`Characteristics: Gender (1=male)`

group<-sample_annotation$`Characteristics: Group`
group[is.na(group)]<-"SSc"

#PCA

result <- prcomp(datExpr , scale=T)
summary(result)
#biplot(result)
result$rotation
result$x %>% 
  as.data.frame() %>%
  mutate(group=group) %>% 
  ggplot(aes(x=PC1,y=PC2))+geom_point(aes(color=group))+
  ggtitle("PCA")

result$x %>% 
  as.data.frame() %>%
  mutate(group=group) %>% 
  ggplot(aes(x=PC1))+geom_histogram(aes(fill=group))

importance_df<-summary(result)$importance[3,] %>% as.data.frame()%>% rownames_to_column(var = "PC")%>% rownames_to_column()
colnames(importance_df)<-c("index","PC","Cumulative_Proportion")
importance_df %>% 
  ggplot(aes(x=as.integer(index),y=Cumulative_Proportion))+geom_bar(stat = "identity")
importance_df %>% 
  ggplot(aes(x=as.integer(index),y=Cumulative_Proportion))+geom_line(stat = "identity")




result <- prcomp(datExpr[group=="SSc",] , scale=T)
summary(result)
#biplot(result)
result$rotation
result$x %>% 
  as.data.frame() %>%
  #mutate(group=c(rep("SSc",59),rep("Control",43))) %>% 
  ggplot(aes(x=PC1,y=PC2))+geom_point()


######################################################
#MDS

data.dist<-dist(datExpr, method = "euclidean", diag = TRUE, upper = TRUE)
data.cmd<-cmdscale(data.dist, eig=TRUE)

data.cmd$points %>% 
  as.data.frame() %>% 
  mutate(group=group) %>% 
  ggplot(aes(x=V1,y=V2))+
  geom_point(aes(color=group))+
  ggtitle("metric MDS")




data.dist<-dist(datExpr[group=="SSc",], method = "euclidean", diag = TRUE, upper = TRUE)
data.cmd<-cmdscale(data.dist, eig=TRUE)

data.cmd$points %>% 
  as.data.frame() %>% 
  ggplot(aes(x=V1,y=V2))+
  geom_point()+
  ggtitle("metric MDS")

######################################################
#tSNE0

set.seed(1) # ?Č????̊m??
# verbose=TRUE ?œr???o?߂??o?͂????܂?
data.tsne <- Rtsne(datExpr, check_duplicates = FALSE, verbose=TRUE, perplexity = 20)
#colors = rainbow(length(unique(iris$Species)))
#plot(data.tsne$Y, t='n', main="Rtsne")
#text(data.tsne$Y)

data.tsne$Y %>% as.data.frame() %>% 
  mutate(group=group) %>% 
  ggplot(aes(x=V1,y=V2))+
  geom_point(aes(color=group))+
  ggtitle("tSNE")



d<-datExpr[group=="SSc",]
#d<-datExpr[group=="Control",]


set.seed(1) # ?Č????̊m??
# verbose=TRUE ?œr???o?߂??o?͂????܂?
data.tsne <- Rtsne(as.matrix(d), check_duplicates = FALSE, verbose=TRUE, perplexity = 10)
#colors = rainbow(length(unique(iris$Species)))
#plot(data.tsne$Y, t='n', main="Rtsne")
#text(data.tsne$Y)

data.tsne$Y %>% as.data.frame() %>% 
  ggplot(aes(x=V1,y=V2))+
  geom_point()+
  ggtitle("tSNE")

######################################################



set.seed(1) # ?Č????̊m??
# verbose=TRUE ?œr???o?߂??o?͂????܂?
data.tsne <- Rtsne(datExpr, check_duplicates = FALSE, verbose=TRUE, perplexity = 20)
#colors = rainbow(length(unique(iris$Species)))
#plot(data.tsne$Y, t='n', main="Rtsne")
#text(data.tsne$Y)

data.tsne$Y %>% as.data.frame() %>% 
  mutate(group=group,gender=sample_annotation$`Characteristics: Gender (1=male)`) %>% 
  ggplot(aes(x=V1,y=V2))+
  geom_point(aes(color=as.factor(gender)))+
  facet_wrap(~group)+
  ggtitle("tSNE")












#########################################
#reference

library(readr)
library(Rtsne)
# The competition datafiles are in the directory ../input
# Read competition data files:
train <- read_csv("../input/train.csv")
test <- read_csv("../input/test.csv")
train$label <- as.factor(train$label)

# shrinking the size for the time limit
numTrain <- 10000
set.seed(1)
rows <- sample(1:nrow(train), numTrain)
train <- train[rows,]
# using tsne
set.seed(1) # for reproducibility
tsne <- Rtsne(train[,-1], dims = 2, perplexity=30, verbose=TRUE, max_iter = 500)
# visualizing
colors = rainbow(length(unique(train$label)))
names(colors) = unique(train$label)
plot(tsne$Y, t='n', main="tsne")
text(tsne$Y, labels=train$label, col=colors[train$label])

# compare with pca
pca = princomp(train[,-1])$scores[,1:2]
plot(pca, t='n', main="pca")
text(pca, labels=train$label,col=colors[train$label])





001_GOEnrichment_.R





library(tidyverse)
library(stringr)
setwd("~/Project/20190709_SystemicSclerosis/3.Script")
rm(list = ls(all.names = TRUE))


output_file_base_path<-"../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/"
if (!dir.exists(paste0(output_file_base_path,"01_Compare_deepsplit_signed","/GOEnrichment"))){
  dir.create(paste0(output_file_base_path,"01_Compare_deepsplit_signed","/GOEnrichment"))
}
condition_deepsplit<-0:4
condition_signed_unsigned<-c("signed","unsigned")
condition_merged<-c("","_merged")

condition<-list(rep(condition_deepsplit,each=4),rep(condition_signed_unsigned,each=2,5),rep(condition_merged,10))
for(i in 1:length(condition[[1]])){
  deepSplit <- condition[[1]][[i]]
  signed_unsigned<-condition[[2]][[i]]
  merged_unmerged<-condition[[3]][[i]]
  
  dir_path<-paste(output_file_base_path, signed_unsigned, "_", deepSplit,"/", sep = "")
  
  go_file_path<-paste(dir_path,"GOEnrichmentTable_",signed_unsigned,"_",deepSplit,merged_unmerged,".csv",sep="")
  print(go_file_path)
  go_data<-read_csv(go_file_path)
  #print(head(go_data))
  go_data$termName<-str_wrap(go_data$termName, width = 40)
  
  go_data %>% mutate()
  tmp<-go_data%>% distinct(module) 
  tmp<-tmp%>% mutate(index=row.names(tmp))%>% 
    mutate(col= if_else(as.double(index) < (nrow(tmp)/2), true = 1, false = 2))
  
  ylim<-max(-log10(go_data$BonferoniP)) +1
  go_data %>% 
    left_join(tmp,by="module") %>% 
    #filter(col==half) %>% 
    filter(BonferoniP<=0.05) %>% 
    ggplot(aes(x=reorder(termName,-BonferoniP),y=-log10(BonferoniP)))+
    geom_bar(aes(fill=module),stat = "identity",width = 0.6)+#c("#0086ab","#f79646"))+
    facet_grid(module~.,  scales='free_y',  space = "free")+
    labs(y="-log(p-value)",x="")+
    theme_classic()+
    theme(
      axis.text=element_text(size=8),
      panel.background = element_rect(fill = "transparent",color = NA),
      panel.grid.minor = element_line(color = NA), 
      panel.grid.major = element_line(color = NA),
      plot.background = element_rect(fill = "transparent",color = NA) )+
    scale_y_continuous(expand = c(0, 0),limits =c(0,ylim) )+
    coord_flip()+
    guides(fill=FALSE)
  savefile<-paste(output_file_base_path,"01_Compare_deepsplit_signed","/GOEnrichment","/GOEnrichment_",signed_unsigned,"_",deepSplit,merged_unmerged,".png",sep="")
  ggsave(file = savefile, dpi = 320, width = 150, height = 300, units = "mm")
}

#???????ŕ`?悷???o?[?W????
if(FALSE){
for(i in 1:length(condition[[1]])){
  deepSplit <- condition[[1]][[i]]
  signed_unsigned<-condition[[2]][[i]]
  merged_unmerged<-condition[[3]][[i]]
  go_file_path<-paste("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/",signed_unsigned,"_",deepSplit,"/GOEnrichmentTable_",signed_unsigned,"_",deepSplit,merged_unmerged,".csv",sep="")
  print(go_file_path)
  go_data<-read_csv(go_file_path)
  #print(head(go_data))
  go_data$termName<-str_wrap(go_data$termName, width = 40)
  
  go_data %>% mutate()
  tmp<-go_data%>% distinct(module) 
  tmp<-tmp%>% mutate(index=row.names(tmp))%>% 
    mutate(col= if_else(as.double(index) < (nrow(tmp)/2), true = 1, false = 2))
  
  ylim<-max(-log10(go_data$BonferoniP)) +1
  for (half in 1:2){
    go_data %>% 
      left_join(tmp,by="module") %>% 
      filter(col==half) %>% 
      filter(BonferoniP<=0.05) %>% 
      ggplot(aes(x=reorder(termName,desc(BonferoniP)),y=-log10(BonferoniP)))+
      geom_bar(aes(fill=module),stat = "identity",width = 0.6)+#c("#0086ab","#f79646"))+
      facet_grid(module~.,  scales='free_y',  space = "free")+
      labs(y="-log(p-value)",x="")+
      theme_classic()+
      theme(
        axis.text=element_text(size=8),
        panel.background = element_rect(fill = "transparent",color = NA),
        panel.grid.minor = element_line(color = NA), 
        panel.grid.major = element_line(color = NA),
        plot.background = element_rect(fill = "transparent",color = NA) )+
      scale_y_continuous(expand = c(0, 0),limits =c(0,ylim) )+
      coord_flip()+
      guides(fill=FALSE)
    savefile<-paste(dir_path,"GOEnrichment/GOEnrichment_",signed_unsigned,"_",deepSplit,merged_unmerged,"_",half,".png",sep="")
    ggsave(file = savefile, dpi = 320, width = 150, height = 300, units = "mm")
  }
}
}







001_GOEnrichment_withBP.R





#series matrix???p????WGCNA
#???W???[????GO enrichment????#
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
library(WGCNA)
library(anRichment)
library(tidyverse)
options(stringsAsFactors = FALSE)

output_file_base_path<-"../2.Output/02_signed_2_unmerged/"

lnames <- load("../2.Output/01_Compare_deepsplit_signed/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/Expdata_quontile.Rdata")
annot <- read.csv("../1.Data/Annotation/GeneAnnotation_GPL10558.csv", header = T, stringsAsFactors = F)
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor
#?????????܂?????
deepSplit <- 2
signed_unsigned<-"signed"
merged_unmerged<-""

dir_path<-paste(output_file_base_path, signed_unsigned, "_", deepSplit,"/", sep = "")

lnames <- load(paste(dir_path,"networkConstruction_StepByStep_",signed_unsigned, "_", deepSplit, merged_unmerged, ".Rdata", sep = ""))
keep<-rowSums(Expdata_quontile$E>log2(50))>=102
Expdata<-Expdata_quontile$E[keep,]%>% t()
#Probe ID??Gene ID?ɕϊ?
probes <- colnames(Expdata)
probes2annot <- match(probes, annot$ID)
#Entrez Gene ID?iLocusLinkID?j???ǂݍ???
allLLIDs <- annot$Entrez_Gene_ID[probes2annot]

#GO enrichment????
#########################
#?e???W???[???ɂ???top 10 term???Ԃ?
GOenr <- GOenrichmentAnalysis(get(paste("moduleColors_",signed_unsigned,sep = "")), allLLIDs, organism = "human", nBestP = 10,
                              ontologies = "BP")

tab <- GOenr$bestPTerms$BP$enrichment
######################################
#???L??depreciated?炵???̂?anRichment package???g???ď???????

GOcollection = buildGOcollection(organism = "human")
GOenrichment = enrichmentAnalysis(
  classLabels = moduleColors_signed, identifiers = allLLIDs,
  refCollection = GOcollection,
  useBackground = "given",
  threshold = 1e-4,
  thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "grey")

GO.BPcollection = subsetCollection(GOcollection, tags = "GO.BP")
GO.BPenrichment = enrichmentAnalysis(
  classLabels = moduleColors_signed, identifiers = allLLIDs,
  refCollection = GO.BPcollection,
  useBackground = "given",
  threshold = 1e-4,
  thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "grey")
tag<-GO.BPenrichment$enrichmentTable
#GO enrichment???͌??ʂ??????o??
write.csv(tab, paste(dir_path,"GOEnrichmentTable_", signed_unsigned,"_", deepSplit, merged_unmerged,"withBP", ".csv",sep = "" ), quote = TRUE, row.names = FALSE)


go_file_path<-paste(dir_path,"GOEnrichmentTable_",signed_unsigned,"_",deepSplit,merged_unmerged,"withBP",".csv",sep="")
go_data<-read_csv(go_file_path)
#go_data$termName<-str_wrap(go_data$termName, width = 40)

go_data %>% mutate()
tmp<-go_data%>% distinct(module) 
tmp<-tmp%>% mutate(index=row.names(tmp))%>% 
  mutate(col= if_else(as.double(index) < (nrow(tmp)/2), true = 1, false = 2))

ylim<-max(-log10(go_data$BonferoniP)) +1
go_data %>% 
  mutate(termName=str_wrap(termName, width = 50)) %>% 
  left_join(tmp,by="module") %>% 
  filter(module %in% significant_module) %>% 
  filter(BonferoniP<=0.05) %>% 
  ggplot(aes(x=reorder(termName,-BonferoniP),y=-log10(BonferoniP)))+
  geom_bar(aes(fill=module),stat = "identity",width = 0.6)+
  facet_grid(module~.,  scales='free_y')+
  labs(y="-log(Bonferoni p-value)",x="")+
  theme_classic()+
  theme(
    axis.text=element_text(size=8),
    panel.background = element_rect(fill = "transparent",color = NA),
    panel.grid.minor = element_line(color = NA), 
    panel.grid.major = element_line(color = NA),
    plot.background = element_rect(fill = "transparent",color = NA) )+
  scale_y_continuous(expand = c(0, 0),limits =c(0,ylim) )+
  scale_fill_identity()+
  coord_flip()+
  guides(fill=FALSE)
savefile<-paste0(output_file_base_path,"GOEnrichment_withBP_ver3",".png")
ggsave(file = savefile, dpi = 320, width = 150, height = 300, units = "mm")
text_size=12
go_data %>% 
  mutate(termName=str_wrap(termName, width = 50)) %>% 
  left_join(tmp,by="module") %>% 
  filter(module %in% significant_module) %>% 
  filter(BonferoniP<=0.05) %>% 
  mutate(color=module %>% 
           str_replace_all("lightyellow","yellow") %>% 
           str_replace_all("lightcyan","cyan")) %>% 
  ggplot(aes(x=reorder(termName,BonferoniP),y=-log10(BonferoniP)))+
  geom_bar(aes(fill=color),stat = "identity",width = 0.6)+
  facet_grid(.~module,  scales='free_x')+
  labs(y="-log(Bonferoni p-value)",x="")+
  theme_classic()+
  theme(
    axis.text=element_text(),#size=text_size),
    axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5),
    panel.background = element_rect(fill = "transparent",color = NA),
    panel.grid.minor = element_line(color = NA), 
    panel.grid.major = element_line(color = NA),
    plot.background = element_rect(fill = "transparent",color = NA) )+
  scale_y_continuous(expand = c(0, 0),limits =c(0,27) )+
  scale_fill_identity()+
  #coord_flip()+
  guides(fill=FALSE)
savefile<-paste0(output_file_base_path,"GOEnrichment_withBP_ver6",".png")
ggsave(file = savefile, dpi = 320, width = 300, height = 140, units = "mm")

go_data %>% 
  left_join(tmp,by="module") %>% 
  filter(module %in% significant_module) %>% 
  filter(BonferoniP<=0.05) %>% 
  mutate(termName=termName %>% str_replace_all("endoplasmic reticulum","ER"),
         termName=str_wrap(termName, width =45 ),
         color=module %>% 
           str_replace_all("lightyellow","yellow") %>% 
           str_replace_all("lightcyan","cyan")) %>% 
  group_by(module )%>%
  top_n(n = 5, wt = -BonferoniP ) %>% 
  ggplot(aes(x=reorder(termName,BonferoniP),y=-log10(BonferoniP)))+
  geom_bar(aes(fill=color),stat = "identity",width = 0.6)+
  facet_grid(.~module,  scales='free_x')+
  labs(y="-log(Bonferoni p-value)",x="")+
  theme_classic()+
  theme(
    axis.text=element_text(size=8),#size=text_size),
    axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")  )+
  scale_y_continuous(expand = c(0, 0),limits =c(0,27) )+
  scale_fill_identity()+
  #coord_flip()+
  guides(fill=FALSE)
savefile<-paste0(output_file_base_path,"GOEnrichment_withBP_ver7",".png")
ggsave(file = savefile, dpi = 320, width = 330, height = 140, units = "mm",  bg = "transparent")






001_MakeWGCNAList.R





library(WGCNA)
library(anRichment)
library(tidyverse)
library(stringr)
options(stringsAsFactors = FALSE)
setwd("~/Project/20190709_SystemicSclerosis/3.Script/")
rm(list = ls(all.names = TRUE))

output_file_base_path<-"../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/"

annot<-read_csv("../1.Data/Annotation/GeneAnnotation_GPL10558.csv")
GeneID2Symbol<-annot[,c("Entrez_Gene_ID","Symbol")] %>% 
  filter(!is.na(Entrez_Gene_ID)) %>% 
  distinct(Entrez_Gene_ID,.keep_all =TRUE)

load("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/Expdata_quontile.Rdata")
datExpr<-Expdata_quontile$E
keep<-rowSums(datExpr>log2(50))>=102
datExpr<-datExpr[keep,]%>% t()

probes <- colnames(datExpr)
probes2annot <- match(probes, annot$ID)
#Entrez Gene ID?iLocusLinkID?j???ǂݍ???
allLLIDs <- annot$Entrez_Gene_ID[probes2annot]


if (!dir.exists(paste0(output_file_base_path,"WGCNAgenelist"))){
  dir.create(paste0(output_file_base_path,"WGCNAgenelist"))
}

condition_deepsplit<-0:4
condition_signed_unsigned<-c("signed","unsigned")
condition_merged<-c("","_merged")
condition<-list(rep(condition_deepsplit,each=4),rep(condition_signed_unsigned,each=2,5),rep(condition_merged,10))
for(i in 1:length(condition[[1]])){
  deepSplit <- condition[[1]][[i]]
  signed_unsigned<-condition[[2]][[i]]
  merged_unmerged<-condition[[3]][[i]]
  
  dir_path<-paste(output_file_base_path, signed_unsigned, "_", deepSplit,"/", sep = "")
  
  lnames <- load(paste(dir_path,"networkConstruction_StepByStep_",signed_unsigned,"_", deepSplit, merged_unmerged,".Rdata", sep = ""))

  WGCNAgenelist <- data.frame(Entrez_Gene_ID=allLLIDs, Module=get(paste("moduleColors_",signed_unsigned,sep="")))%>%
    distinct(.keep_all = TRUE) %>%
    left_join(GeneID2Symbol,by = "Entrez_Gene_ID")%>%
    .[,c(3,2,1)] %>% 
    write_csv(paste(output_file_base_path,"WGCNAgenelist/WGCNAgenelist_",signed_unsigned,"_", deepSplit, merged_unmerged,".csv", sep = ""))
  
}






001_modified_series_matrix_module_GOEnrichment.R





#series matrix???p????WGCNA
#???W???[????GO enrichment????#
setwd("~/Project/20190709_SystemicSclerosis/3.Script")
library(WGCNA)
options(stringsAsFactors = FALSE)

output_file_base_path<-"../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/"

lnames <- load("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/Expdata_quontile.Rdata")
annot <- read.csv("../1.Data/Annotation/GeneAnnotation_GPL10558.csv", header = T, stringsAsFactors = F)


condition_deepsplit<-0:4
condition_signed_unsigned<-c("signed","unsigned")
condition_merged<-c("","_merged")

condition<-list(rep(condition_deepsplit,each=4),rep(condition_signed_unsigned,each=2,5),rep(condition_merged,10))
for(i in 1:length(condition[[1]])){
  deepSplit <- condition[[1]][[i]]
  signed_unsigned<-condition[[2]][[i]]
  merged_unmerged<-condition[[3]][[i]]
  
  dir_path<-paste(output_file_base_path, signed_unsigned, "_", deepSplit,"/", sep = "")
  
  lnames <- load(paste(dir_path,"networkConstruction_StepByStep_",signed_unsigned, "_", deepSplit, merged_unmerged, ".Rdata", sep = ""))
  #Probe ID??Gene ID?ɕϊ?
  probes <- colnames(datExpr)
  probes2annot <- match(probes, annot$ID)
  #Entrez Gene ID?iLocusLinkID?j???ǂݍ???
  allLLIDs <- annot$Entrez_Gene_ID[probes2annot]
  
  #GO enrichment????
  #?e???W???[???ɂ???top 10 term???Ԃ?
  GOenr <- GOenrichmentAnalysis(get(paste("moduleColors_",signed_unsigned,sep = "")), allLLIDs, organism = "human", nBestP = 10)
  tab <- GOenr$bestPTerms[[4]]$enrichment
  
  #GO enrichment???͌??ʂ??????o??
  write.csv(tab, paste(dir_path,"GOEnrichmentTable_", signed_unsigned,"_", deepSplit, merged_unmerged, ".csv",sep = "" ), quote = TRUE, row.names = FALSE)
  


}



if(FALSE){
#deepSplit??0?`4?????ꂼ?????͂??Cfor()?ɂ??艺?L?̃R?[?h???J???Ԃ????s
for(i in 0:4){
  deepSplit <- i #???̐??l??0?`4?ɕς???
  
  ##unsigned_merge?Ȃ?
  #?t?H???_?ړ?
  #setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, sep = ""))
  #?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
  lnames <- load(paste("unsigned_",deepSplit,"/networkConstruction_StepByStep_unsigned_", deepSplit, ".Rdata", sep = ""))
  #Probe ID??Gene ID?ɕϊ?
  probes <- colnames(datExpr)
  probes2annot <- match(probes, annot$ID)
  #Entrez Gene ID?iLocusLinkID?j???ǂݍ???
  allLLIDs <- annot$Entrez_Gene_ID[probes2annot]
  
  #GO enrichment????
  #?e???W???[???ɂ???top 10 term???Ԃ?
  GOenr <- GOenrichmentAnalysis(moduleColors_unsigned, allLLIDs, organism = "human", nBestP = 10)
  tab <- GOenr$bestPTerms[[4]]$enrichment
  
  #GO enrichment???͌??ʂ??????o??
  write.csv(tab, paste("unsigned_",deepSplit,"/GOEnrichmentTable_unsigned_", deepSplit, ".csv",sep = "" ), quote = TRUE, row.names = FALSE)
  
  
  ##unsigned_merge????
  #?t?H???_?ړ?
  #setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, sep = ""))
  #?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
  lnames <- load(paste("unsigned_",deepSplit,"/networkConstruction_StepByStep_unsigned_", deepSplit, "_merged.Rdata", sep = ""))
  #Probe ID??Gene ID?ɕϊ?
  probes <- colnames(datExpr)
  probes2annot <- match(probes, annot$ID)
  #Entrez Gene ID?iLocusLinkID?j???ǂݍ???
  allLLIDs <- annot$Entrez_Gene_ID[probes2annot]
  
  #GO enrichment????
  #?e???W???[???ɂ???top 10 term???Ԃ?
  GOenr <- GOenrichmentAnalysis(moduleColors_unsigned, allLLIDs, organism = "human", nBestP = 10)
  tab <- GOenr$bestPTerms[[4]]$enrichment
  
  #GO enrichment???͌??ʂ??????o??
  write.csv(tab, paste("unsigned_",deepSplit,"/GOEnrichmentTable_unsigned_", deepSplit, "_merged.csv",sep = "" ), quote = TRUE, row.names = FALSE)
  
  
  ##signed_merge?Ȃ?
  #?t?H???_?ړ?
  #setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
  #?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
  lnames <- load(paste("signed_",deepSplit,"/networkConstruction_StepByStep_signed_", deepSplit, ".Rdata", sep = ""))
  #Probe ID??Gene ID?ɕϊ?
  probes <- colnames(datExpr)
  probes2annot <- match(probes, annot$ID)
  #Entrez Gene ID?iLocusLinkID?j???ǂݍ???
  allLLIDs <- annot$Entrez_Gene_ID[probes2annot]
  
  #GO enrichment????
  #?e???W???[???ɂ???top 10 term???Ԃ?
  GOenr <- GOenrichmentAnalysis(moduleColors_signed, allLLIDs, organism = "human", nBestP = 10)
  tab <- GOenr$bestPTerms[[4]]$enrichment
  
  #GO enrichment???͌??ʂ??????o??
  write.csv(tab, paste("signed_",deepSplit,"/GOEnrichmentTable_signed_", deepSplit, ".csv",sep = "" ), quote = TRUE, row.names = FALSE)
  
  
  ##signed_merge????
  #?t?H???_?ړ?
  #setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
  #?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
  lnames <- load(paste("signed_",deepSplit,"/networkConstruction_StepByStep_signed_", deepSplit, "_merged.Rdata", sep = ""))
  #Probe ID??Gene ID?ɕϊ?
  probes <- colnames(datExpr)
  probes2annot <- match(probes, annot$ID)
  #Entrez Gene ID?iLocusLinkID?j???ǂݍ???
  allLLIDs <- annot$Entrez_Gene_ID[probes2annot]
  
  #GO enrichment????
  #?e???W???[???ɂ???top 10 term???Ԃ?
  GOenr <- GOenrichmentAnalysis(moduleColors_signed, allLLIDs, organism = "human", nBestP = 10)
  tab <- GOenr$bestPTerms[[4]]$enrichment
  
  #GO enrichment???͌??ʂ??????o??
  write.csv(tab, paste("signed_",deepSplit,"/GOEnrichmentTable_signed_", deepSplit, "_merged.csv",sep = "" ), quote = TRUE, row.names = FALSE)
}
}





001_Normarize_Expdata.R





setwd("~/Project/20190709_SystemicSclerosis/3.Script")

library(limma)
data <- read.ilmn("../1.Data/GSE58095_non-normalized.txt", sep = "\t", probeid = "ID_REF", expr = "SAMPLE")
normdata <- backgroundCorrect(data, method = "normexp", offset=20) #offset??????
#normalizeBetweenArrays?̏o?͂?log2?????ďo?Ă????炵??: https://support.bioconductor.org/p/62765/
Expdata_quontile <-  normalizeBetweenArrays(normdata, method ="quantile" )#"scale")
#plotDensities(Expdata_quontile, col="black")
save(Expdata_quontile,file="Expdata_quontile.Rdata")

boxplot(Expdata_quontile$E)
#Expdata$E
#PCA <- prcomp(t(Expdata$E), scale=T)
#SSC_col=rgb(0.5, 0.5, 0.5, alpha=0.5) #SSC?Q?̐F?̐ݒ?
#HC_col=rgb(1, 0, 0, alpha=0.5) #HC?Q?̐F?̐ݒ?
#Color = c(rep(SSC_col, 59), rep(HC_col, 43)) #?e?Qn=3?̂Ƃ??̃T???v???̐F?w??
#plot(PCA$x[,1], PCA$x[,2], col=Color, pch=16, cex=2.5, xlab="PCA1",ylab="PCA2")
#plot(hclust(dist(t(Expdata$E))))





001_PathwayAnalysis_with_ReactomePA.R





#Ref: https://bioconductor.org/packages/release/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html

#BiocManager::install("ReactomePA")
library(ReactomePA)
library(tidyverse)
setwd("~/Project/20190709_SystemicSclerosis/3.Script")


base_path<-"../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/02_signed_2_unmerged/Modules/Significant/all/"
colors<-list.files(base_path) %>% str_sub(end=-5)
colors<-colors[colors!="significant"]
#color<-colors[[1]]
module_data<-read_csv(paste0(base_path,color,".csv"))



for (color in colors){
  print(color)
  assign(color,read_csv(paste0(base_path,color,".csv")))
  if (nrow(get(color))==0){
    next
  }
  x <- enrichPathway(organism ="human",gene=get(color) %>% .$Gene, pvalueCutoff=0.1, readable=T)
  if (nrow(as.data.frame(x))==0){
    next
  }
  pdf(file=paste0("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/02_signed_2_unmerged/PA/","PA_",color,".pdf"), width=11.6,height=8.2)
  print(barplot(x, showCategory=8))
  print(dotplot(x, showCategory=15))
  print(emapplot(x))
  print(cnetplot(x, categorySize="pvalue"))
  dev.off()
}






001_RetrieveSignificantModules.R





library(tidyverse)
library(dplyr)
setwd("\\\\cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

#MEs_df<-read_csv("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/signed_3/MEs_signed_3.csv") %>% as.data.frame()
#group<-c(rep("SSc",59),rep("Control",43))
#result<-matrix(0,1,ncol(MEs_df)-1)
#colnames(result) <- colnames(MEs_df)[2:ncol(MEs_df)]

#???????ݒ?
signed_unsigned<-"signed"
deepSplit<-2
merged_unmerged<-""
file_path<-paste0("../2.Output/02_signed_2_unmerged/signed_2/","/MEs_",signed_unsigned,"_",deepSplit,".csv")
MEs_df<-read_csv(file_path) %>% as.data.frame()

sample_annotation<-read_tsv("../1.Data/GSE58095_all_sample_annotatiions.txt",skip = 4,n_max = 102)
group<-sample_annotation$`Characteristics: Group`
group[is.na(group)]<-"SSc"

r<-rep(1,ncol(MEs_df)-1)
###################################
for (i in 2:(ncol(MEs_df))){
  r[[i-1]]<- t.test(MEs_df[group=="SSc",i], MEs_df[group=="Control",i], var.equal = F, alternative = "two.sided", paired = F)$p.value
}

color_pvalue_df<-tibble(ModuleColor=colnames(MEs_df[2:ncol(MEs_df)]),pvalue=r) %>% 
  mutate(qvalue=pvalue %>% p.adjust( method = "BH"))
color_pvalue_df %>% 
  mutate(ModuleColor=str_sub(ModuleColor,3,-1)) %>% 
  #filter(pvalue<0.01) %>% 
  write_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv")
####################################
color_pvalue_df<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv")

n_significant_module<-color_pvalue_df %>% filter(pvalue<0.01) %>% nrow()
color_code_df<-tibble(module=significant_module,
                      color_code=col2hex(significant_module %>% str_replace_all("white","antiquewhite")%>% str_replace_all("lightyellow","lightyellow2")))


color_pvalue_df<-color_pvalue_df %>% 
  filter(pvalue<0.01) %>% 
  arrange(pvalue) %>% 
  left_join(color_code_df,by=c("ModuleColor"="module")) 
color_pvalue_df%>% 
  ggplot(aes(x=reorder(ModuleColor,-pvalue),y=-log10(pvalue)))+
  #geom_bar(aes(fill=ModuleColor),stat = "identity")+
  geom_bar(aes(),stat = "identity",fill=color_pvalue_df %>% arrange(-pvalue) %>% .$color_code)+
  labs(y="-log(p-value)",x="")+
  theme_classic()+
  theme(
    axis.text=element_text(size=8),
    panel.background = element_rect(fill = "transparent",color = NA),
    panel.grid.minor = element_line(color = NA), 
    panel.grid.major = element_line(color = NA),
    plot.background = element_rect(fill = "transparent",color = NA) )+
  scale_y_continuous(expand = c(0, 0))+
  coord_flip()+
  #ggtitle(paste0(signed_unsigned," ",deepSplit," ",merged_unmerged,"\n","significant / all = ",n_significant_module," / ",nrow(color_pvalue_df))) +
  guides(fill=FALSE)
ggsave(paste0("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/01_Compare_deepsplit_signed/",
              "significant_modules_",n_significant_module,"_of_",nrow(color_pvalue_df),".png"),dpi=300)






001_Retrieve_Genes_of_each_Module.R





library(tidyverse)
library(dplyr)

setwd("~/Project/20190709_SystemicSclerosis/3.Script")
signed_unsigned<-"signed"
deepSplit<-2
merged_unmerged<-""
dir_path<-paste0("../2.Output/01_Compare_deepsplit_signed/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/",signed_unsigned,"_", deepSplit,"/")

#module color?????[?h
lnames <- load(paste0(dir_path,"networkConstruction_StepByStep_signed_", deepSplit, merged_unmerged,".Rdata"))
moduleColors_signed
#datExpr?????[?h
load("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/datExpr_filtered.Rdata")

gene_annotation <- read_csv("../1.Data/Annotation/GeneAnnotation_GPL10558.csv")
cols<-c("ID","Symbol","Entrez_Gene_ID")


sample_annotation<-read_tsv("../1.Data/GSE58095_all_sample_annotatiions.txt",skip = 4,n_max = 102)
group<-sample_annotation$`Characteristics: Group`
group[is.na(group)]<-"SSc"


#Symbol, ID, module color???܂ރf?[?^?t???[?????쐬
symbol_color_df<-tibble(ID=datExpr %>% t() %>% rownames(),
       module_color=moduleColors_signed) %>% 
  left_join(gene_annotation[,cols],by = "ID")

#ILMNID, logFC, p.value, TestExpression, ControlExpression???܂ރf?[?^?t???[?????쐬
#?e?X?̌v?Z?p?֐??????`
calculate_logFC = function(x) {
  #print(x)
  #print(class(x))
  x<-x[2:103]
  x<-as.double(x)
  #print(x)
  logFC<-mean(x[group=="SSc"],na.rm=T)-mean(x[group=="Control"],na.rm=T)
  return(logFC)
}
calculate_pvalue = function(x) {
  x<-x[2:103]
  x<-as.double(x)
  return(t.test(x[group=="SSc"],x[group=="Control"])$p.value)
}

logFC_df<-datExpr %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var="ILMNID") 
logFC_df<-logFC_df %>% 
  mutate(logFC=apply(logFC_df,1,calculate_logFC),
         p.value=apply(logFC_df,1,calculate_pvalue),
         TestExpression=apply(logFC_df,1,function(x){mean(as.double(x[2:103])[group=="SSc"],na.rm=T)}),
         ControlExpression=apply(logFC_df,1,function(x){mean(as.double(x[2:103])[group=="Control"],na.rm=T)})) %>%
  dplyr::select(ILMNID,logFC,p.value,TestExpression,ControlExpression) 

symbol_color_df<-left_join(logFC_df,symbol_color_df,by = c("ILMNID"="ID"))
symbol_color_df %>% 
  write_csv("../2.Output/02_signed_2_unmerged/logFC_ModuleColor.csv")



#?L?ӂȃ??W???[???J???[???擾
file_path<-paste0("../2.Output/02_signed_2_unmerged/",signed_unsigned,"_",deepSplit,"/MEs_",signed_unsigned,"_",deepSplit,".csv")
MEs_df<-read_csv(file_path) %>% as.data.frame()
r<-rep(1,ncol(MEs_df)-1)

for (i in 2:(ncol(MEs_df))){
  r[[i-1]]<- t.test(MEs_df[group=="SSc",i], MEs_df[group=="Control",i], var.equal = F, alternative = "two.sided", paired = F)$p.value
}
color_pvalue_df<-tibble(ModuleColor=colnames(MEs_df[2:ncol(MEs_df)]),pvalue=r)
write_csv(color_pvalue_df,"../2.Output/02_signed_2_unmerged/color_pvalue.csv")
significant_module<-color_pvalue_df %>% filter(pvalue<0.01) %>% .$ModuleColor %>% str_sub(3)

colors<-symbol_color_df$module_color %>% unique()


#?L?ӂ??L?ӂłȂ????Ńf?B???N?g??????????
#?e???W???[????GeneID?ƃV???{?????ۑ?
base_path<-"../2.Output/02_signed_2_unmerged/Modules/"
sig_dir_path<-paste0(base_path,"Significant/")
if (!dir.exists(sig_dir_path)){
  dir.create(sig_dir_path)
}
nsig_dir_path<-paste0(base_path,"NoSignificant/")
if (!dir.exists(nsig_dir_path)){
  dir.create(nsig_dir_path)
}

for (color in colors){
  if(color %in% significant_module){
    result_dir_path<-paste0(base_path,"Significant/all/")
  }else{
    result_dir_path<-paste0(base_path,"NoSignificant/all/")
  }
  print(result_dir_path)
  if (!dir.exists(result_dir_path)){
    dir.create(result_dir_path)
  }
  result<-symbol_color_df %>% 
    filter(module_color==color) %>% #,!is.na(Entrez_Gene_ID)) %>% 
    select(-module_color) 
  #colnames(result)<-c("Gene","logFC","p.value","Test Expression","Control Expression","Symbol")
  write_csv(result,paste0(result_dir_path,color,".csv"))
}

#######################
#logFC?̌v?Z

calculate_logFC = function(x) {
  #print(x)
  #print(class(x))
  x<-x[2:103]
  x<-as.double(x)
  #print(x)
  logFC<-mean(x[group=="SSc"],na.rm=T)-mean(x[group=="Control"],na.rm=T)
  return(logFC)
}
calculate_pvalue = function(x) {
  x<-x[2:103]
  x<-as.double(x)
  return(t.test(x[group=="SSc"],x[group=="Control"])$p.value)
}
test<-datExpr %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var="ILMNID") 
test<-test %>% 
  mutate(logFC=apply(test,1,calculate_logFC),
         p.value=apply(test,1,calculate_pvalue),
         TestExpression=apply(test,1,function(x){mean(as.double(x[2:103])[group=="SSc"],na.rm=T)}),
         ControlExpression=apply(test,1,function(x){mean(as.double(x[2:103])[group=="Control"],na.rm=T)})) %>%
  dplyr::select(ILMNID,logFC,p.value,TestExpression,ControlExpression)



test %>% ggplot(aes(logFC,-log10(p.value)))+
  geom_point(aes(color=(abs(logFC) > 0.5 & p.value < 0.01)))+
  theme_classic()+
  theme(
    legend.position="none",
    panel.background = element_rect(fill = "transparent",color = NA),
    panel.grid.minor = element_line(color = NA), 
    panel.grid.major = element_line(color = NA),
    plot.background = element_rect(fill = "transparent",color = NA) )+
  scale_color_manual(values = c("#bfbfbf","#0086ab","#f79646","#9bbb59","#da6272","#777777"))+
  ggtitle("SSc vs Control\nthreshold: abs(logFC) > 0.5 & p.value < 0.01")
ggsave("../2.Output/02_signed_2_unmerged/VolcanoPlot_SSc_Vs_Control.png",dpi = 300)




datExpr %>% 
  as.data.frame() %>% 
  #rownames_to_column(var = "sample") %>% 
  mutate(group=group) %>% 
  gather(key = ID, value = Expression, -group) %>% 
  #arrange(Expression) %>% 
  #spread(key = group,value = Expression)
  group_by(group,ID) %>% 
  summarise(Exp=mean(Expression)) %>% 
  spread(key = group,value = Exp) %>% 
  mutate(logFC=log2(SSc/Control))
  





datExpr %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "sample") %>%  
  gather(key=ID,value =Expression,  -sample)%>% 
  ggplot(aes(sample,Expression))+geom_boxplot()








001_series_matrix_modified.R





#series matrix???p????WGCNA
#?l?b?g???[?N?̍쐬#
setwd("~/Project/20190709_SystemicSclerosis/3.Script")
library(WGCNA)
library(tidyverse)
options(stringsAsFactors = FALSE)


output_file_base_path<-"../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/"
load("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/Expdata_quontile.Rdata")
png_resolution<-list(height = 1440, width = 1440, res = 216)


#WGCNAdata <- log2(Expdata$E)
datExpr <- Expdata_quontile$E
keep<-rowSums(datExpr>log2(50))>=102
datExpr<-datExpr[keep,]%>% t()

#lnames <- load("dataInput.Rdata")

#soft-thresholding Power ???̐ݒ?
powers <- c(c(1:10), seq(from = 12, to = 26, by = 2))

#?l?b?g???[?N?g?|???W?[????
#???ꂼ???̃??ɑ΂???Scale-free fit index??Mean connectivity???Z?o????
#unsigned?i?????ցE?t???ցj?Csigned?i?????ւ̂݁j???ꂼ???Z?o
sft_unsigned <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "unsigned")
sft_signed <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")

#Scale-free fit index??Mean connectivity???}??
sizeGrWindow(9, 5)
png(paste(output_file_base_path,"Scale_free_fit_index_and_Mean_connectivity_unsigned.png"),height = png_resolution$height,width = png_resolution$width,res = png_resolution$res)#, width = 9, height = 5)
par(mfrow = c(1, 2))
cex1 = 0.9
plot(sft_unsigned$fitIndices[,1], -sign(sft_unsigned$fitIndices[,3])*sft_unsigned$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence (unsigned)"))
text(sft_unsigned$fitIndices[,1], -sign(sft_unsigned$fitIndices[,3])*sft_unsigned$fitIndices[,2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red") #R^2?J?b?g?I?t?l0.90?ɐԐ?

plot(sft_unsigned$fitIndices[,1], sft_unsigned$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity (unsigned)"))
text(sft_unsigned$fitIndices[,1], sft_unsigned$fitIndices[,5], labels = powers, cex = cex1, col = "red")
dev.off()

sizeGrWindow(9, 5)
png(paste(output_file_base_path,"Scale_free_fit_index_and_Mean_connectivity_signed.png"),height = png_resolution$height,width = png_resolution$width,res = png_resolution$res)#, width = 9, height = 5)
par(mfrow = c(1, 2))
cex1 = 0.9
plot(sft_signed$fitIndices[,1], -sign(sft_signed$fitIndices[,3])*sft_signed$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence (signed)"))
text(sft_signed$fitIndices[,1], -sign(sft_signed$fitIndices[,3])*sft_signed$fitIndices[,2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red") #R^2?J?b?g?I?t?l0.90?ɐԐ?

plot(sft_signed$fitIndices[,1], sft_signed$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity (signed)"))
text(sft_signed$fitIndices[,1], sft_signed$fitIndices[,5], labels = powers, cex = cex1, col = "red")
dev.off()

#Scalefree plot??R^2?????Ȃ????œK?ȃ????T??
#unsigned?̏ꍇ
png(paste(output_file_base_path,"Scalefree_plot_unsigned.png"),height = png_resolution$height,width = png_resolution$width,res = png_resolution$res)#, width = 12, height = 9)
par(mfrow=c(2,2))
for(i in c(c(6:10), seq(from = 12, to = 22, by = 2))){
  softPower_unsigned <- i
  adjacency_unsigned <- adjacency(datExpr, power = softPower_unsigned, type = "unsigned")
  diag(adjacency_unsigned) <- 0
  connectivity_unsigned <- apply(adjacency_unsigned, 1, sum)
  scaleFreePlot(connectivity_unsigned, truncated = T)
  mtext(paste("beta =", i), side = 3, line = 0, adj = 0)
}
dev.off()

#???̑I?ѕ???scale free fit index??0.90?ɒB?????ŏ??l?Cmean connectivity???Ⴗ???Ă͂????Ȃ??Cscale free R^2??1?ɋ߂?
#?? = 14??scale free fit index??0.90?ɒB?????ŏ??l?ł????Cscale free R^2??0.95
#unsigned?̓? = 14???p???ĉ??͂???
softPower_unsigned <- 12
adjacency_unsigned <- adjacency(datExpr, power = softPower_unsigned, type = "unsigned")

#signed?̏ꍇ
png(paste(output_file_base_path,"Scalefree_plot_signed.png"),height = png_resolution$height,width = png_resolution$width,res = png_resolution$res)#, width = 12, height = 9)
par(mfrow=c(2,2))
for(i in c(seq(from = 12, to = 26, by = 2))){
  softPower_signed <- i
  adjacency_signed <- adjacency(datExpr, power = softPower_signed, type = "signed")
  diag(adjacency_signed) <- 0
  connectivity_signed <- apply(adjacency_signed, 1, sum)
  scaleFreePlot(connectivity_signed, truncated = T)
  mtext(paste("beta =", i), side = 3, line = 0, adj = 0)
}
dev.off()

#?? = 20??scale free fit index??0.90?ɒB?????ŏ??l?ł????Cscale free R^2??0.89
#signed?̓? = 20???p???ĉ??͂???
softPower_signed <- 22
adjacency_signed <- adjacency(datExpr, power = softPower_signed, type = "signed")

#adjacency?i?אڍs???l?j??topological overlap?iTOM?j?ɕϊ?
#unsigned
TOM_unsigned_file_path<-paste(output_file_base_path,"TOM_unsigned.Rdata",sep="")
dissTOM_unsigned_file_path<-paste(output_file_base_path,"dissTOM_unsigned.Rdata",sep="")
if ( !(file.exists(TOM_unsigned_file_path) & file.exists(dissTOM_unsigned_file_path))){
  TOM_unsigned <- TOMsimilarity(adjacency_unsigned, TOMType = "unsigned")
  dissTOM_unsigned <- 1-TOM_unsigned
  save(TOM_unsigned,file = TOM_unsigned_file_path)
  save(dissTOM_unsigned,file = dissTOM_unsigned_file_path)
}else{
  print("loading files...")
  load(TOM_unsigned_file_path)
  load(TOM_unsigned_file_path)
}
#signed
TOM_signed_file_path<-paste(output_file_base_path,"TOM_signed.Rdata",sep="")
dissTOM_signed_file_path<-paste(output_file_base_path,"dissTOM_signed.Rdata",sep="")
if ( !(file.exists(TOM_signed_file_path) & file.exists(dissTOM_signed_file_path))){
  TOM_signed <- TOMsimilarity(adjacency_signed, TOMType = "signed")
  dissTOM_signed <- 1-TOM_signed
  save(TOM_signed,file = paste(output_file_base_path,"TOM_signed.Rdata",sep=""))
  save(dissTOM_signed,file = paste(output_file_base_path,"dissTOM_signed.Rdata",sep=""))
}else{
  print("loading files...")
  load(TOM_signed_file_path)
  load(TOM_signed_file_path)
}
#1-TOM???K?w?I?N???X?^?????O
geneTree_unsigned <- hclust(as.dist(dissTOM_unsigned), method = "average")
geneTree_signed <- hclust(as.dist(dissTOM_signed), method = "average")

#?f???h???O???????}??
sizeGrWindow(12, 9)
png(paste(output_file_base_path,"Gene_clustering_dissTOM_unsigned.png"),height = png_resolution$height,width = png_resolution$width,res = png_resolution$res)#, width = 12, height = 9)
plot(geneTree_unsigned, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity (unsigned)",
     labels = FALSE, hang = 0.04)
dev.off()
png(paste(output_file_base_path,"Gene_clustering_dissTOM_signed.png"),height = png_resolution$height,width = png_resolution$width,res = png_resolution$res)#, width = 12, height = 9)
plot(geneTree_signed, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity (signed)",
     labels = FALSE, hang = 0.04)
dev.off()


#???W???[?????o??eigengene?̎Z?o#
#???W???[???T?C?Y?̍ŏ??l???ݒ?
minModuleSize <- 30


#deepSplit??0?`4??signed unsigned, merged unmerged?̂??ׂĂ̏????????? 
condition_deepsplit<-0:4
condition_signed_unsigned<-c("signed","unsigned")
condition_merged<-c("","_merged")
condition<-list(rep(condition_deepsplit,each=4),rep(condition_signed_unsigned,each=2,5),rep(condition_merged,10))


#########????????#########
for(i in 1:length(condition[[1]])){
  deepSplit <- condition[[1]][[i]]
  signed_unsigned<-condition[[2]][[i]]
  #merged_unmerged<-condition[[3]][[i]]
  

  
  ##signed_merge?Ȃ?
  #?t?H???_?̍쐬?ƈړ?
  dir_path<-paste(output_file_base_path, signed_unsigned, "_", deepSplit,"/", sep = "")
  ifelse(!dir.exists(dir_path), 
         dir.create(dir_path), FALSE)
  #setwd(dir_name)
  
  #Dynamic tree cut???p???ă??W???[?????o
  #dynamicMods_signed <- cutreeDynamic(dendro = geneTree_signed, distM = dissTOM_signed,
                                      #deepSplit = deepSplit, pamRespectsDendro = FALSE, pamStage = FALSE,
                                      #minClusterSize = minModuleSize)
  geneTree<-get(paste("geneTree_",signed_unsigned,sep=""))
  dynamicMods<-paste("dynamicMods_",signed_unsigned,sep="") %>% 
    assign(cutreeDynamic(dendro = geneTree, distM = get(paste("dissTOM_",signed_unsigned,sep="")),
                       deepSplit = deepSplit, pamRespectsDendro = FALSE, pamStage = FALSE,
                       minClusterSize = minModuleSize))
  
  #???W???[??No.???F?ɕϊ????C?e???W???[???̃v???[?u?????ۑ?
  #dynamicColors_signed <- labels2colors(dynamicMods_signed)
  dynamicColors<-paste("dynamicColors_",signed_unsigned,sep="") %>% 
    assign(labels2colors(dynamicMods))
  
  #module_signed <- as.data.frame(table(dynamicColors_signed))
  module<-paste("module_",signed_unsigned,sep="") %>% 
    assign(as.data.frame(table(dynamicColors)))
  write.csv(module, paste(dir_path,"module_",signed_unsigned,"_",deepSplit,".csv", sep = ""))
  
  #module eigengene?iME?j?̎Z?o
  #MEList_signed <- moduleEigengenes(datExpr, colors = dynamicColors_signed)
  MEList<-paste("MEList_",signed_unsigned,sep="") %>% 
    assign(moduleEigengenes(datExpr, colors = dynamicColors))
  #MEs_signed <- MEList_signed$eigengenes
  MEs<-paste("MEs_",signed_unsigned,sep="") %>% 
    assign(MEList %>% .$eigengenes)
  
  #ME???p???????W???[???̃N???X?^?????O
  #MEDiss_signed <- 1-cor(MEs_signed)
  MEDiss <- paste("MEDiss_",signed_unsigned,sep="") %>% 
    assign(1-cor(MEs))
  
  #METree_signed <- hclust(as.dist(MEDiss_signed), method = "average")
  METree <- paste("METree",signed_unsigned,sep="") %>% 
    assign(hclust(as.dist(MEDiss), method = "average"))
  
  png(paste(dir_path,"Clustering_ME_",signed_unsigned,"_", deepSplit, ".png",sep = ""),height = png_resolution$height,width = png_resolution$width,res = png_resolution$res)#, width = 7, height = 6)
  plot(METree, main = paste("Clustering of module eigengenes (",signed_unsigned," deepSplit, ", deepSplit, ")", sep = ""),
       xlab = "", sub = "")
  #ME???߂??i?????????Ă????j???W???[?????T??
  #ME?f???h???O?????̃J?b?g?I?t?l???ݒ肷??
  MEDissThres <- 0.25 #0.15?`0.25?͈̔͂Ńf???h???O???????݂ēK?X????
  #?J?b?g???C????????
  #abline(h = MEDissThres, col = "red")
  dev.off()
  
  #?f?[?^??rename???ĕۑ?
  #moduleColors_signed <- dynamicColors_signed
  moduleColors <- paste("modulecolors_",signed_unsigned,sep = "") %>% 
    assign(dynamicColors)
  colorOrder <- c("grey", standardColors(100))
  #moduleLabels_signed <- match(moduleColors_signed, colorOrder)-1
  moduleLabels <- paste("moduleLabels_",signed_unsigned,sep="") %>% 
    assign(match(moduleColors, colorOrder)-1)
  #save(MEs_signed, moduleLabels_signed, moduleColors_signed, geneTree_signed, 
       #file = paste("networkConstruction_StepByStep_signed_", deepSplit, ".Rdata", sep = ""))
  save(MEs, moduleLabels, moduleColors, geneTree, 
       file = paste(dir_path,"networkConstruction_StepByStep_",signed_unsigned,"_", deepSplit, ".Rdata", sep = ""))
  
  #ME?̕ۑ?
  #MEs_signed2 <- cbind(rownames(datExpr), MEs_signed)
  MEs2 <- paste("MEs_",signed_unsigned,"2",sep = "") %>% 
    assign(cbind(rownames(datExpr), MEs))
  write.csv(MEs2, paste("MEs_",signed_unsigned,"_",deepSplit,".csv", sep = ""), row.names = F)
  
  #intramodular connectivity (kWithin)?̎Z?o
  #kwithin_signed <- intramodularConnectivity(adjacency_signed, moduleColors_signed)
  kwithin <- intramodularConnectivity(get(paste("adjacency_",signed_unsigned,sep = "")), moduleColors)
  
  #?v???[?u?Ckwithin?Ccolor?̕\???쐬???C?ۑ?
  #binded_signed <- cbind(rownames(kwithin_signed), kwithin_signed, moduleColors_signed)
  binded <- cbind(rownames(kwithin), kwithin, moduleColors)
  colnames(binded_signed)[1] <- "ProbeID"
  write.csv(binded, paste(dir_path,"Probe_kWitnin_color_",signed_unsigned, "_", deepSplit, ".csv", sep=""), row.names=F)

  
  
  
###############################################################################
  ##signed_merge????
  #?J?b?g?I?t?l?ȉ??̃N???X?^?[??merge
  #merge_signed <- mergeCloseModules(datExpr, dynamicColors_signed, cutHeight = MEDissThres, verbose = 3)
  merged <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  mergedColors <- merged$colors
  #merge???̊e???W???[???̃v???[?u?????ۑ?
  merged_module <- as.data.frame(table(mergedColors))
  write.csv(merged_module, paste(dir_path,"module_",signed_unsigned,"_",deepSplit,"_merged.csv", sep = ""))
  
  #merge???????W???[????ME
  mergedMEs <- merge$newMEs
  #?f???h???O??????merge?O???̐F?t?????????W???[?????}??
  png(paste("Gene_dendrogram_and_module_colors_signed_", deepSplit, ".png",sep = ""),height = png_resolution$height,width = png_resolution$width,res = png_resolution$res)#, width = 9, height = 6)
  plotDendroAndColors(geneTree_signed, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = paste("Gene dendrogram and module colors (",signed_unsigned,", deepSplit ", deepSplit, ")", sep = ""))
  dev.off()
  
  #merge????ME???p???????W???[???̃N???X?^?????O
  MEDiss_merged <- 1-cor(mergedMEs)
  METree_merged <- hclust(as.dist(MEDiss_merged), method = "average")
  png(paste("Clustering_ME_",signed_unsigned, "_", deepSplit, "_merged.png",sep = ""),height = png_resolution$height,width = png_resolution$width,res = png_resolution$res)#, width = 7, height = 6)
  plot(METree_merged, main = paste("Clustering of module eigengenes (",signed_unsigned, ", deepSplit ", deepSplit, ", merged)", sep = ""),
       xlab = "", sub = "")
  dev.off()
  
  #?f?[?^??rename???ĕۑ?
  moduleColors <- mergedColors
  colorOrder <- c("grey", standardColors(100))
  moduleLabels <- match(moduleColors, colorOrder)-1
  save(mergedMEs, moduleLabels, moduleColors, geneTree, 
       file = paste("networkConstruction_StepByStep_",signed_unsigned, "_", deepSplit, "_merged.Rdata", sep = ""))
  
  #ME?̕ۑ?
  mergedMEs2 <- cbind(rownames(datExpr), mergedMEs)
  write.csv(mergedMEs2, paste("MEs_",signed_unsigned, "_",deepSplit,"_merged.csv", sep = ""), row.names = F)
  
  #intramodular connectivity (kWithin)?̎Z?o
  kwithin <- intramodularConnectivity(adjacency, moduleColors)
  
  #?v???[?u?Ckwithin?Ccolor?̕\???쐬???C?ۑ?
  binded <- cbind(rownames(kwithin), kwithin, moduleColors)
  colnames(binded)[1] <- "ProbeID"
  write.csv(binded, paste("Probe_kWitnin_color_",signed_unsigned, "_", deepSplit, "_merged.csv", sep=""), row.names=F)
}
#########?????܂?#########







001_series_matrix_module_CellTypeEnrichment_mod.R





#series matrix???p????WGCNA
#???W???[????Cell-Type enrichment????#
library(WGCNA)
library(anRichment)
#library(anRichmentMethods)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(stringr)
options(stringsAsFactors = FALSE)
#lnames <- load("Expdata_quontile.Rdata")
#setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
annot <- read_csv("../1.Data/Annotation/GeneAnnotation_GPL10558.csv")


#?q?g?̍זE?ʈ??`?q???X?g?̍쐬
#???肽??Collection?????w?肷??
collection <- newCollection()
#Organism???w?肷??
organism <- "Human"

#?t?H???_????Genelist???ꊇ?œǂݍ??܂??邽?߂̐ݒ?
#setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data/Human_Liver_CellType")
#?t?@?C???i?[?p?X???w??
files <- list.files(path = "../1.Data/Cell_Enrivhment/Skin+Blood/Jensen_TISSUES/",full.names = T)
#files <- list.files(path = "C:/Users/ywt2100/Desktop/human",full.names = T)
files_names <- list.files("../1.Data/Cell_Enrivhment/Skin+Blood/Jensen_TISSUES/" , full.names = F, pattern="csv") 
#files_names <- list.files("C:/Users/ywt2100/Desktop/human" , full.names = F, pattern="csv") 
files_names <- gsub(".csv", "", files_names) 

#Genelist??Geneset?ɂ??邽?߂̊֐????ǂݍ??܂???
Togeneset <- function (X, Y){
  colnames(X) <- "Gene" #???????????I??Gene?ɕύX
  Symbol = X$Gene  #?֋X???̃??l?[??
  Entrez.0 = convert2entrez(organism = organism, symbol = Symbol)
  Entrez = unique(Entrez.0)#???`?q??(Symbol)??entrez?`???ɕϊ?????
  print(table(is.finite(Entrez)))  #entrez?`???ɕϊ??ł??????`?q???̊m?F
  return(newGeneSet(
    geneEntrez = Entrez,
    geneEvidence = "IEP",
    geneSource = "",
    ID = Y,
    name = Y,
    description = Y,
    source = "",
    organism = organism,
    internalClassification = Y,
    groups = Y,
    lastModified = "")
  )
}
#?convert2entrez

#?t?H???_????Genelist???ꊇ??geneset?ɕϊ????Ă??????L??Collection?ɓ?????
for(i in 1:length(files)){
  genes <- read_csv(files[i])#?????????t?@?C?????ǂݍ???
  genes2<- genes[1]#$Symbol
  assign(files_names[i], genes2)
  genes <- unique(genes)#?d???L???ꍇ????????
  genes_set <- Togeneset(genes, files_names[i])#genelist??geneset??
  assign(paste(files_names[i]), genes_set)
  collection <- addToCollection(collection, genes_set)#Geneset?????ꂽ??collection?ɓ?????
  rm(genes, genes_set)
}

#???`?q?Z?b?g?ۑ?
#setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
#save(collection, file="Human_skin_blood_CellType_collection.Rdata")
save(collection, file="test.Rdata")

lname<-load("../1.Data/Cell_Enrichment/Human_skin_blood_CellType_collection.Rdata")


#CelltypeEnrichment?֐????ǂݍ??܂???
###################?֐?????????###################
CelltypeEnrichment <- function (X){
  names(X) <- c("Gene", "Group") #???????????I??Gene??Group?ɕύX
  entrez = X$Gene  #?֋X???̃??l?[??
  Group = X$Group?@#?֋X???̃??l?[??
  print(table(Group))  #???`?q???X?g?̐??̊m?F?p?B
  
  #????????anRichment?p?b?P?[?W?̒???enrichmentAnalysis?֐????g?p
  analysisresult = enrichmentAnalysis(?@?@
    classLabels = Group, identifiers = entrez,?@#classLabels???????̂??̂??????̉??͑Ώۈ??`?q?Q?Ƃ??ĔF??
    refCollection = my_collection,  #reference?Ƃ??Đݒ肷???f?[?^?Z?b?g?B
    #useGroups=files_names %>% as.vector(),
    useBackground = "given",?@?@#???̓o?b?N?O???E???h?̐ݒ??i?d?v?j?B?????̓C???v?b?g?????S???`?q?̂????ǂݍ??߂????̂??g?p
    threshold = 1,?@?@#?G?????b?`?????g???͌??ʂ??Ƃ??ďo?͂?????csv?t?@?C???ł́Cp?l?̂??????l
    thresholdType = "Bonferroni",?@#???Lp?l?ɕt?????āCBonferroni?␳????p?l??臒l?Ƃ???
    getOverlapEntrez = FALSE,?@?@#?o??csv?t?@?C???ɂāC?I?[?o?[???b?v???????`?q????entrez?`???ŏo?͂????C
    getOverlapSymbols = TRUE,    #?o??csv?t?@?C???ɂāC?I?[?o?[???b?v???????`?q????symbol?`???ŏo?͂????C
    maxReportedOverlapGenes = 10000,?@?@#???L?I?[?o?[???b?v???????`?q???ǂ̂??炢?\???????邩???ݒ??B?S?Ăق????̂?10,000??
    removeDuplicatesInDifferentClasses =FALSE,?@#???????W???[?????ɓ??????`?q?????????ꍇ???C???̂܂܏???????
    entrySeparator = ",",?@?@#?I?[?o?[???b?v???????`?q?Q?ɂ??āC","?ŕ????ďo?͂??Ă????B?ǂ??ł??悵
    ignoreLabels = "grey", #grey???W???[???̈??`?q?Q?ɂ??Ă͖????i???͂??Ȃ??j
    combineEnrichmentTables = FALSE) 
  
  
  countsInDataSet <- analysisresult$effectiveBackgroundSize  #???̓o?b?N?O???E???h?̈??`?q?????m??
  print(table(countsInDataSet))?@?@?@?@?@?@?@?@?@?@?@?@?@?@?@?@#???????o??
  
  
  Resulttable <- analysisresult$enrichmentTable
  Resulttable <- separate(Resulttable, overlapGenes, into=as.character(c(1:1000)), sep=",")
  Resulttable2 <- subset(analysisresult$enrichmentTable, analysisresult$enrichmentTable$pValue < 1) #?o?͗p?ɐ??`
  options(warn=-1) #???̃R?[?h??warning???o???B???????????????R?[?h
  #Overlap???`?q?Q???C1???`?q???ƂɃZ???ɕ????ĕ\???????邽?߂̃R?[?h(1000?͔C?ӁB)???̍??Ƃ?warning???o?邪???Ə??͖????Ȃ?
  write.table(Resulttable2, file = "WGCNA_enrichment_result.csv",row.names = FALSE,sep=",") #???ʂ?csv?t?@?C???ɏo??
  list <- by(Resulttable, Resulttable$class, data.frame) #?O???[?v???ƂɃt?@?C???`???????X?g??
  sapply(1:dim(list), function(x){write.csv(list[x], file=paste0("Result_", dimnames(list)[[1]][x], ".csv"), row.names=FALSE)})?@#?O???[?v???ƂɌ??ʂ??o??
  
  #?????????C???͌??ʃf?[?^(Resulttable)?̂????O???t?ɕK?v?Ȃ??̂??Ԃ????ʂ?
  Rank_all <- sapply(1:dim(list), function(x){list(list[[x]][1:5,c(4,6,9)])})#?O???[?v???Ƃ?TOP10?̕K?v?ȗ????????o??
  names(Rank_all) <- sapply(1:dim(list), function(x){names(Rank_all) <- names(list[x])})#???L?ŃO???[?v???????????̂ł??????????l?[??
  for (i in 1:dim(list)) Rank_all[[i]]$nCommonGenes <- paste("(",Rank_all[[i]]$nCommonGenes,")") #?O???t?̃??x???p?̖??O???`
  for (i in 1:dim(list)) Rank_all[[i]]$pValue <- -(log10(as.numeric(Rank_all[[i]]$pValue)))?@#?O???t?̃??x???p?̖??O???`
  Rank_all <- na.omit(Rank_all) #???ʂ?TOP10?ɖ????Ȃ??ꍇ??NA????
  for (i in 1:dim(list)) Rank_all[[i]] <- transform(Rank_all[[i]], "Rename"=(paste(Rank_all[[i]]$dataSetName,Rank_all[[i]]$nCommonGenes,sep="?@")))#?O???t?̃??x???p?̖??O???`
  
  #list Rename
  names(Rank_all) <- sapply(1:dim(list), function(x){names(Rank_all) <- names(list[x])})
  
  
  ##????????ggplot2???g?????}???̎w??
  for (i in 1:dim(list)) ggsave(file=paste0("Result_", names(Rank_all)[i], ".pdf"),plot = ggplot(Rank_all[[i]], aes(x=reorder(Rank_all[[i]]$Rename, Rank_all[[i]]$pValue), y=Rank_all[[i]]$pValue)) +  
                                  geom_bar(stat="identity", width=.5,fill="black") +   
                                  coord_flip() +                                     
                                  xlab("Cell-type\n(#Overlap genes)") + 
                                  ylab("-log(P-value)"))
  
  #Cell-type enrichment p-value?̃q?[?g?}?b?v???}??
  CellTyperesult <- read.csv("WGCNA_enrichment_result.csv")
  logP <- -log10(CellTyperesult$pValue)
  CellTypeP <- data.frame(Module = CellTyperesult$class, CellType = CellTyperesult$dataSetID, logP = logP)
  ghm <- ggplot(CellTypeP, aes(x = CellType, y = Module, fill = logP))
  ghm <- ghm + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  ghm <- ghm + scale_x_discrete(limits = c("Hep3_cluster5", "Hep5_cluster14", "Hep4_cluster6", "Hep6_cluster15", "Hep1_cluster1", "Hep2_cluster3", 
                                           "PeriportalLSEC_cluster11", "CentralVenousLSEC_cluster12", "PortalEndothelialCell_cluster13",
                                           "Cholangiocyte_cluster17", "HepaticStellateCell_cluster20",
                                           "InflammatoryMacs_cluster4", "NonInflammatoryMacs_cluster10",
                                           "CD3abTcell_cluster2", "gdTcell1_cluster9", "gdTcell2_cluster18", "NKlikeCell_cluster8",
                                           "MatureBcell_cluster16", "PlasmaCell_cluster7", "ErythroidCell_cluster19"))
  ghm <- ghm + geom_tile(aes(fill = logP))
  ghm <- ghm + xlab("Cell Type") + ylab("Module")
  ghm <- ghm + geom_text(aes(label = round(CellTypeP$logP, 1)), size = 2)
  ghm <- ghm + scale_fill_gradient(low = "white", high = "red", limits = c(0, 15))
  pdf("CellType_pvalue_heatmap.pdf", width = 5, height = 7)
  plot(ghm)
  dev.off()
  
}
###################?֐??????܂?###################
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
load("../2.Output/01_Compare_deepsplit_signed/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/Expdata_quontile.Rdata")
datExpr <- Expdata_quontile$E
keep<-rowSums(datExpr>log2(50))>=102
datExpr<-datExpr[keep,]%>% t()

#collection$groups<-files_names %>% as.list()

#???ۂɉ???

#deepSplit??0?`4?????ꂼ?????͂??Cfor()?ɂ??艺?L?̃R?[?h???J???Ԃ????s
for(i in 0:4){
  deepSplit <- i #???̐??l??0?`4?ɕς???
  
  ##unsigned_merge?Ȃ?
  #?t?H???_?ړ?
  #setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, sep = ""))
  #?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
  lnames <- load(paste("unsigned_",deepSplit,"/networkConstruction_StepByStep_unsigned_", deepSplit, ".Rdata", sep = ""))
  #Probe ID??Entrez Gene ID?ɕϊ?????
  probes <- colnames(datExpr)
  probes2annot <- match(probes, annot$ID)
  #Entrez Gene ID?iLocusLinkID?j???ǂݍ???
  allLLIDs <- annot$Entrez_Gene_ID[probes2annot]
  #WGCNA?ɂ??????`?q?ƃ??W???[???̕\???쐬
  WGCNAgenelist <- data.frame(allLLIDs, moduleColors_unsigned)
  #?t?H???_?쐬?ƈړ?
  ifelse(!dir.exists(paste("unsigned_", deepSplit, "/CellType", sep = "")), 
         dir.create(paste("unsigned_", deepSplit, "/CellType", sep = "")), FALSE)
  #setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, "/CellType", sep = ""))
  #???W???[????Cell-Type enrichment????
  CelltypeEnrichment(WGCNAgenelist)


  ##unsigned_merge????
  #?t?H???_?ړ?
  #setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, sep = ""))
  #?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
  lnames <- load(paste("unsigned_",deepSplit,"/networkConstruction_StepByStep_unsigned_", deepSplit, "_merged.Rdata", sep = ""))
  #Probe ID??Entrez Gene ID?ɕϊ?????
  probes <- colnames(datExpr)
  probes2annot <- match(probes, annot$ID)
  #Entrez Gene ID?iLocusLinkID?j???ǂݍ???
  allLLIDs <- annot$Entrez_Gene_ID[probes2annot]
  #WGCNA?ɂ??????`?q?ƃ??W???[???̕\???쐬
  WGCNAgenelist <- data.frame(allLLIDs, moduleColors_unsigned)
  #?t?H???_?쐬?ƈړ?
  ifelse(!dir.exists(paste("unsigned_", deepSplit, "/CellType_merged", sep = "")), 
         dir.create(paste("unsigned_", deepSplit, "/CellType_merged", sep = "")), FALSE)
  #setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, "/CellType_merged", sep = ""))
  #???W???[????Cell-Type enrichment????
  CelltypeEnrichment(WGCNAgenelist)


  ##signed_merge?Ȃ?
  #?t?H???_?ړ?
  #setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
  #?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
  lnames <- load(paste("signed_",deepSplit,"/networkConstruction_StepByStep_signed_", deepSplit, ".Rdata", sep = ""))
  #Probe ID??Entrez Gene ID?ɕϊ?????
  probes <- colnames(datExpr)
  probes2annot <- match(probes, annot$ID)
  #Entrez Gene ID?iLocusLinkID?j???ǂݍ???
  allLLIDs <- annot$Entrez_Gene_ID[probes2annot]
  #WGCNA?ɂ??????`?q?ƃ??W???[???̕\???쐬
  WGCNAgenelist <- data.frame(allLLIDs, moduleColors_signed) 
  #?t?H???_?쐬?ƈړ?
  ifelse(!dir.exists(paste("signed_", deepSplit, "/CellType", sep = "")), 
         dir.create(paste("signed_", deepSplit, "/CellType", sep = "")), FALSE)
  #setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, "/CellType", sep = ""))
  #???W???[????Cell-Type enrichment????
  CelltypeEnrichment(WGCNAgenelist)
  
  
  ##signed_merge????
  #?t?H???_?ړ?
  #setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
  #?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
  lnames <- load(paste("signed_",deepSplit,"/networkConstruction_StepByStep_signed_", deepSplit, "_merged.Rdata", sep = ""))
  #Probe ID??Entrez Gene ID?ɕϊ?????
  probes <- colnames(datExpr)
  probes2annot <- match(probes, annot$ID)
  #Entrez Gene ID?iLocusLinkID?j???ǂݍ???
  allLLIDs <- annot$Entrez_Gene_ID[probes2annot]
  #WGCNA?ɂ??????`?q?ƃ??W???[???̕\???쐬
  WGCNAgenelist <- data.frame(allLLIDs, moduleColors_signed)
  #?t?H???_?쐬?ƈړ?
  ifelse(!dir.exists(paste("signed_", deepSplit, "/CellType_merged", sep = "")), 
         dir.create(paste("signed_", deepSplit, "/CellType_merged", sep = "")), FALSE)
  #setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, "/CellType_merged", sep = ""))
  #???W???[????Cell-Type enrichment????
  CelltypeEnrichment(WGCNAgenelist)
}
###########################################################
deepSplit <- 2 #???̐??l??0?`4?ɕς???

##signed_merge?Ȃ?
#?t?H???_?ړ?
#setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("../2.Output/02_signed_2_unmerged/","signed_",deepSplit,"/networkConstruction_StepByStep_signed_", deepSplit, ".Rdata", sep = ""))
#Probe ID??Entrez Gene ID?ɕϊ?????
probes <- colnames(datExpr )
probes2annot <- match(probes, annot$ID)
#Entrez Gene ID?iLocusLinkID?j???ǂݍ???
allLLIDs <- annot$Entrez_Gene_ID[probes2annot]
#WGCNA?ɂ??????`?q?ƃ??W???[???̕\???쐬
WGCNAgenelist <- data.frame(allLLIDs, moduleColors_signed) 
#?t?H???_?쐬?ƈړ?
ifelse(!dir.exists(paste("../2.Output/02_signed_2_unmerged/","signed_", deepSplit, "/CellType", sep = "")), 
       dir.create(paste("../2.Output/02_signed_2_unmerged/","signed_", deepSplit, "/CellType", sep = "")), FALSE)
#setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, "/CellType", sep = ""))
#???W???[????Cell-Type enrichment????
CelltypeEnrichment(WGCNAgenelist)







001_series_matrix_module_GOEnrichment.R





#series matrix???p????WGCNA
#???W???[????GO enrichment????#
setwd("~/Project/20190709_SystemicSclerosis/3.Script")
library(WGCNA)
options(stringsAsFactors = FALSE)
lnames <- load("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/Expdata_quontile.Rdata")
#setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
annot <- read.csv("../1.Data/Annotation/GeneAnnotation_GPL10558.csv", header = T, stringsAsFactors = F)

#deepSplit??0?`4?????ꂼ?????͂??Cfor()?ɂ??艺?L?̃R?[?h???J???Ԃ????s
for(i in 0:4){
deepSplit <- i #???̐??l??0?`4?ɕς???

##unsigned_merge?Ȃ?
#?t?H???_?ړ?
#setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("unsigned_",deepSplit,"/networkConstruction_StepByStep_unsigned_", deepSplit, ".Rdata", sep = ""))
#Probe ID??Gene ID?ɕϊ?
probes <- colnames(datExpr)
probes2annot <- match(probes, annot$ID)
#Entrez Gene ID?iLocusLinkID?j???ǂݍ???
allLLIDs <- annot$Entrez_Gene_ID[probes2annot]

#GO enrichment????
#?e???W???[???ɂ???top 10 term???Ԃ?
GOenr <- GOenrichmentAnalysis(moduleColors_unsigned, allLLIDs, organism = "human", nBestP = 10)
tab <- GOenr$bestPTerms[[4]]$enrichment

#GO enrichment???͌??ʂ??????o??
write.csv(tab, paste("unsigned_",deepSplit,"/GOEnrichmentTable_unsigned_", deepSplit, ".csv",sep = "" ), quote = TRUE, row.names = FALSE)


##unsigned_merge????
#?t?H???_?ړ?
#setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("unsigned_",deepSplit,"/networkConstruction_StepByStep_unsigned_", deepSplit, "_merged.Rdata", sep = ""))
#Probe ID??Gene ID?ɕϊ?
probes <- colnames(datExpr)
probes2annot <- match(probes, annot$ID)
#Entrez Gene ID?iLocusLinkID?j???ǂݍ???
allLLIDs <- annot$Entrez_Gene_ID[probes2annot]

#GO enrichment????
#?e???W???[???ɂ???top 10 term???Ԃ?
GOenr <- GOenrichmentAnalysis(moduleColors_unsigned, allLLIDs, organism = "human", nBestP = 10)
tab <- GOenr$bestPTerms[[4]]$enrichment

#GO enrichment???͌??ʂ??????o??
write.csv(tab, paste("unsigned_",deepSplit,"/GOEnrichmentTable_unsigned_", deepSplit, "_merged.csv",sep = "" ), quote = TRUE, row.names = FALSE)


##signed_merge?Ȃ?
#?t?H???_?ړ?
#setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("signed_",deepSplit,"/networkConstruction_StepByStep_signed_", deepSplit, ".Rdata", sep = ""))
#Probe ID??Gene ID?ɕϊ?
probes <- colnames(datExpr)
probes2annot <- match(probes, annot$ID)
#Entrez Gene ID?iLocusLinkID?j???ǂݍ???
allLLIDs <- annot$Entrez_Gene_ID[probes2annot]

#GO enrichment????
#?e???W???[???ɂ???top 10 term???Ԃ?
GOenr <- GOenrichmentAnalysis(moduleColors_signed, allLLIDs, organism = "human", nBestP = 10)
tab <- GOenr$bestPTerms[[4]]$enrichment

#GO enrichment???͌??ʂ??????o??
write.csv(tab, paste("signed_",deepSplit,"/GOEnrichmentTable_signed_", deepSplit, ".csv",sep = "" ), quote = TRUE, row.names = FALSE)


##signed_merge????
#?t?H???_?ړ?
#setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("signed_",deepSplit,"/networkConstruction_StepByStep_signed_", deepSplit, "_merged.Rdata", sep = ""))
#Probe ID??Gene ID?ɕϊ?
probes <- colnames(datExpr)
probes2annot <- match(probes, annot$ID)
#Entrez Gene ID?iLocusLinkID?j???ǂݍ???
allLLIDs <- annot$Entrez_Gene_ID[probes2annot]

#GO enrichment????
#?e???W???[???ɂ???top 10 term???Ԃ?
GOenr <- GOenrichmentAnalysis(moduleColors_signed, allLLIDs, organism = "human", nBestP = 10)
tab <- GOenr$bestPTerms[[4]]$enrichment

#GO enrichment???͌??ʂ??????o??
write.csv(tab, paste("signed_",deepSplit,"/GOEnrichmentTable_signed_", deepSplit, "_merged.csv",sep = "" ), quote = TRUE, row.names = FALSE)
}





001_series_matrix_module_trait_relationship.R





#series matrix???p????WGCNA
#???W???[????trait?̊֘A??#
library(WGCNA)
library(tidyverse)

setwd("~/Project/20190709_SystemicSclerosis/3.Script")

options(stringsAsFactors = FALSE)
lnames <- load("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/Expdata_quontile.Rdata")
datExpr<-Expdata_quontile$E
keep<-rowSums(datExpr>log2(50))>=102
datExpr<-datExpr[keep,]%>%t()
annot <- read.csv("../1.Data/Annotation/GeneAnnotation_GPL10558.csv", header = T, stringsAsFactors = F)

#datTraits?̐??̂??킩?????ȉ??A?T???v???̃A?m?e?[?V?????????p
datTraits<-read.table("../1.Data/GSE58095_all_sample_annotatiions.txt",header = TRUE, skip=4,sep="\t",nrows = 102)
choosen<-c(#7,12,
           13)#,15,16)
datTraits<-datTraits[,13] %>% as.data.frame()
names(datTraits)<-"Skin_score"
#deepSplit??0?`4?????ꂼ?????͂??C???L?̃R?[?h???J???Ԃ????s

#########????????#########
for(i in 0:4){
  deepSplit <- i #???̐??l??0?`4?ɕς???

##unsigned_merge?Ȃ?
#?t?H???_?ړ?
#setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/unsigned_", deepSplit,"/networkConstruction_StepByStep_unsigned_", deepSplit, ".Rdata", sep = ""))



#???`?q?i?v???[?u?j???ƃT???v?????????`
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

#ME???Čv?Z
MEs0 <- moduleEigengenes(datExpr, moduleColors_unsigned)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, datTraits, use = "p",method = "spearman")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

#???W???[????trait?̊֘A?????q?[?g?}?b?v?Ŏ???
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
#dev.new()
#pdf(paste("../2.Output/Signed_Unsigned_Deepsplit/unsigned_", deepSplit,"/Module_trait_relationships_unsigned_", deepSplit, ".pdf",sep = ""), width = 10, height = 18)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships (unsigned, deepSplit ", deepSplit, ")", sep = ""))
#dev.off()

#???W???[????trait?̑??֌W???y??p-value???\?ɂ??ĕۑ?
colnames(textMatrix) <- names(datTraits)
rownames(textMatrix) <- names(MEs)
write.csv(textMatrix, paste("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/unsigned_", deepSplit,"/Module_trait_correlation_unsigned_", deepSplit, ".csv", sep = ""))

#???ڂ???trait?i???̏ꍇ??Fibrosis?j?Ɋւ???Gene significance (GS)??Module membership (MM)???Z?o????
skinscore <- as.data.frame(datTraits$Skin_score)
names(skinscore) <-"skinscore"
modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

geneTraitSignificance <- as.data.frame(cor(datExpr, skinscore, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(skinscore), sep = "")
names(GSPvalue) <- paste("p.GS.", names(skinscore), sep = "")

#Probe ID??gene name?ɕϊ?
probes <- colnames(datExpr)
probes2annot <- match(probes, annot$ID)
#?A?m?e?[?V?????????Ă??Ȃ??v???[?u????0?ł??邩?m?F
sum(is.na(probes2annot)) 

#?A?m?e?[?V?????????y??skinscore?Ɋւ???GS?̕\???쐬
geneInfo0 <- data.frame(Probe_ID = probes,
                        geneSymbol = annot$Symbol[probes2annot],
                        Entrez_gene_ID = annot$Entrez_Gene_ID[probes2annot],
                        moduleColor = moduleColors_unsigned,
                        geneTraitSignificance,
                        GSPvalue)
#skinscore??p?l?ŕ??בւ?
modOrder <- order(-abs(cor(MEs, skinscore, use = "p")))
#MM???????ǉ?
for(mod in 1:ncol(geneModuleMembership))
{
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(geneInfo0, geneModuleMembership[,modOrder[mod]],
                          MMPvalue[, modOrder[mod]])
  names(geneInfo0) <- c(oldNames, paste("MM.", modNames[modOrder[mod]], sep = ""),
                        paste("p.MM.", modNames[modOrder[mod]], sep = ""))
}

#IPA annotation???????ǉ?
#setwd("E:/IPA_gene_annotation")
#IPAannot <- read.csv("../1.Data/Annotation/GeneAnnotation_GPL10558.csv", header = T, stringsAsFactors = F)
#IPAannot <- IPAannot[,-c(2:4)]
#geneInfo0 <- merge(geneInfo0, IPAannot, by.x = 1, by.y = 1, all.x = T)

#module color?y??GS?ŕ??בւ?
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.skinscore))
geneInfo <- geneInfo0[geneOrder, ]

#?A?m?e?[?V?????????y??skinscore?Ɋւ???GS/MM?̕\???????o??
#setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, sep = ""))
write.csv(geneInfo, paste("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/unsigned_",deepSplit,"/geneInfo_unsigned_", deepSplit, ".csv",sep = "" ), row.names = FALSE)


##unsigned_merge????
#?t?H???_?ړ?
#setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/unsigned_",deepSplit,"/networkConstruction_StepByStep_unsigned_", deepSplit, "_merged.Rdata", sep = ""))

#???`?q?i?v???[?u?j???ƃT???v?????????`
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

#ME???Čv?Z
MEs0 <- moduleEigengenes(datExpr, moduleColors_unsigned)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

#???W???[????trait?̊֘A?????q?[?g?}?b?v?Ŏ???
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
#pdf(paste("../2.Output/Signed_Unsigned_Deepsplit/unsigned_",deepSplit,"/Module_trait_relationships_unsigned_", deepSplit, "_merged.pdf",sep = ""), width = 10, height = 18)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships (unsigned, deepSplit ", deepSplit, ", merged)", sep = ""))
#dev.off()

#???W???[????trait?̑??֌W???y??p-value???\?ɂ??ĕۑ?
colnames(textMatrix) <- names(datTraits)
rownames(textMatrix) <- names(MEs)
write.csv(textMatrix, paste("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/unsigned_",deepSplit,"/Module_trait_correlation_unsigned_", deepSplit, "_merged.csv", sep = ""))

#???ڂ???trait?i???̏ꍇ??skinscore?j?Ɋւ???Gene significance (GS)??Module membership (MM)???Z?o????
skinscore <- as.data.frame(datTraits$Skin_score)
names(skinscore) <-"skinscore"
modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

geneTraitSignificance <- as.data.frame(cor(datExpr, skinscore, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(skinscore), sep = "")
names(GSPvalue) <- paste("p.GS.", names(skinscore), sep = "")

#Probe ID??gene name?ɕϊ?
probes <- colnames(datExpr)
probes2annot <- match(probes, annot$ID)
#?A?m?e?[?V?????????Ă??Ȃ??v???[?u????0?ł??邩?m?F
sum(is.na(probes2annot)) 

#?A?m?e?[?V?????????y??skinscore?Ɋւ???GS?̕\???쐬
geneInfo0 <- data.frame(Probe_ID = probes,
                        geneSymbol = annot$Symbol[probes2annot],
                        Entrez_gene_ID = annot$Entrez_Gene_ID[probes2annot],
                        moduleColor = moduleColors_unsigned,
                        geneTraitSignificance,
                        GSPvalue)
#skinscore??p?l?ŕ??בւ?
modOrder <- order(-abs(cor(MEs, skinscore, use = "p")))
#MM???????ǉ?
for(mod in 1:ncol(geneModuleMembership))
{
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(geneInfo0, geneModuleMembership[,modOrder[mod]],
                          MMPvalue[, modOrder[mod]])
  names(geneInfo0) <- c(oldNames, paste("MM.", modNames[modOrder[mod]], sep = ""),
                        paste("p.MM.", modNames[modOrder[mod]], sep = ""))
}

#IPA annotation???????ǉ?
#setwd("E:/IPA_gene_annotation")
#IPAannot <- read.csv("GPL14951_IPAannot.csv", header = T, stringsAsFactors = F)
#IPAannot <- IPAannot[,-c(2:4)]
#geneInfo0 <- merge(geneInfo0, IPAannot, by.x = 1, by.y = 1, all.x = T)

#module color?y??GS?ŕ??בւ?
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.skinscore))
geneInfo <- geneInfo0[geneOrder, ]
#?A?m?e?[?V?????????y??skinscore?Ɋւ???GS/MM?̕\???????o??
#setwd(paste("../2.Output/Signed_Unsigned_Deepsplit/unsigned_", deepSplit, sep = ""))
write.csv(geneInfo, paste("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/unsigned_", deepSplit,"/geneInfo_unsigned_", deepSplit, "_merged.csv",sep = "" ), row.names = FALSE)

##signed_merge?Ȃ?
#?t?H???_?ړ?
#setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/signed_",deepSplit,"/networkConstruction_StepByStep_signed_", deepSplit, ".Rdata", sep = ""))

#???`?q?i?v???[?u?j???ƃT???v?????????`
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

#ME???Čv?Z
MEs0 <- moduleEigengenes(datExpr, moduleColors_signed)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

#???W???[????trait?̊֘A?????q?[?g?}?b?v?Ŏ???
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
#pdf(paste("../2.Output/Signed_Unsigned_Deepsplit/signed_",deepSplit,"/Module_trait_relationships_signed_", deepSplit, ".pdf",sep = ""), width = 10, height = 18)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships (signed, deepSplit ", deepSplit, ")", sep = ""))
#dev.off()

#???W???[????trait?̑??֌W???y??p-value???\?ɂ??ĕۑ?
colnames(textMatrix) <- names(datTraits)
rownames(textMatrix) <- names(MEs)
write.csv(textMatrix, paste("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/signed_",deepSplit,"/Module_trait_correlation_signed_", deepSplit, ".csv", sep = ""))

#???ڂ???trait?i???̏ꍇ??skinscore?j?Ɋւ???Gene significance (GS)??Module membership (MM)???Z?o????
skinscore <- as.data.frame(datTraits$Skin_score)
names(skinscore) <-"skinscore"
modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

geneTraitSignificance <- as.data.frame(cor(datExpr, skinscore, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(skinscore), sep = "")
names(GSPvalue) <- paste("p.GS.", names(skinscore), sep = "")

#Probe ID??gene name?ɕϊ?
probes <- colnames(datExpr)
probes2annot <- match(probes, annot$ID)
#?A?m?e?[?V?????????Ă??Ȃ??v???[?u????0?ł??邩?m?F
sum(is.na(probes2annot)) 

#?A?m?e?[?V?????????y??skinscore?Ɋւ???GS?̕\???쐬
geneInfo0 <- data.frame(Probe_ID = probes,
                        geneSymbol = annot$Symbol[probes2annot],
                        Entrez_gene_ID = annot$Entrez_Gene_ID[probes2annot],
                        moduleColor = moduleColors_signed,
                        geneTraitSignificance,
                        GSPvalue)
#skinscore??p?l?ŕ??בւ?
modOrder <- order(-abs(cor(MEs, skinscore, use = "p")))
#MM???????ǉ?
for(mod in 1:ncol(geneModuleMembership))
{
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(geneInfo0, geneModuleMembership[,modOrder[mod]],
                          MMPvalue[, modOrder[mod]])
  names(geneInfo0) <- c(oldNames, paste("MM.", modNames[modOrder[mod]], sep = ""),
                        paste("p.MM.", modNames[modOrder[mod]], sep = ""))
}

#IPA annotation???????ǉ?
#setwd("E:/IPA_gene_annotation")
#IPAannot <- read.csv("GPL14951_IPAannot.csv", header = T, stringsAsFactors = F)
#IPAannot <- IPAannot[,-c(2:4)]
#geneInfo0 <- merge(geneInfo0, IPAannot, by.x = 1, by.y = 1, all.x = T)

#module color?y??GS?ŕ??בւ?
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.skinscore))
geneInfo <- geneInfo0[geneOrder, ]
#?A?m?e?[?V?????????y??skinscore?Ɋւ???GS/MM?̕\???????o??
#setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
write.csv(geneInfo, paste("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/signed_",deepSplit,"/geneInfo_signed_", deepSplit, ".csv",sep = "" ), row.names = FALSE)


##signed_merge????
#?t?H???_?ړ?
#setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/signed_",deepSplit,"/networkConstruction_StepByStep_signed_", deepSplit, "_merged.Rdata", sep = ""))

#???`?q?i?v???[?u?j???ƃT???v?????????`
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

#ME???Čv?Z
MEs0 <- moduleEigengenes(datExpr, moduleColors_signed)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

#???W???[????trait?̊֘A?????q?[?g?}?b?v?Ŏ???
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
#pdf(paste("../2.Output/Signed_Unsigned_Deepsplit/signed_",deepSplit,"/Module_trait_relationships_signed_", deepSplit, "_merged.pdf",sep = ""), width = 10, height = 18)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships (signed, deepSplit ", deepSplit, ", merged)", sep = ""))
#dev.off()

#???W???[????trait?̑??֌W???y??p-value???\?ɂ??ĕۑ?
colnames(textMatrix) <- names(datTraits)
rownames(textMatrix) <- names(MEs)
write.csv(textMatrix, paste("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/signed_",deepSplit,"/Module_trait_correlation_signed_", deepSplit, "_merged.csv", sep = ""))

#???ڂ???trait?i???̏ꍇ??skinscore?j?Ɋւ???Gene significance (GS)??Module membership (MM)???Z?o????
skinscore <- as.data.frame(datTraits$Skin_score)
names(skinscore) <-"skinscore"
modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

geneTraitSignificance <- as.data.frame(cor(datExpr, skinscore, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(skinscore), sep = "")
names(GSPvalue) <- paste("p.GS.", names(skinscore), sep = "")

#Probe ID??gene name?ɕϊ?
probes <- colnames(datExpr)
probes2annot <- match(probes, annot$ID)
#?A?m?e?[?V?????????Ă??Ȃ??v???[?u????0?ł??邩?m?F
sum(is.na(probes2annot)) 

#?A?m?e?[?V?????????y??skinscore?Ɋւ???GS?̕\???쐬
geneInfo0 <- data.frame(Probe_ID = probes,
                        geneSymbol = annot$Symbol[probes2annot],
                        Entrez_gene_ID = annot$Entrez_Gene_ID[probes2annot],
                        moduleColor = moduleColors_signed,
                        geneTraitSignificance,
                        GSPvalue)
#skinscore??p?l?ŕ??בւ?
modOrder <- order(-abs(cor(MEs, skinscore, use = "p")))
#MM???????ǉ?
for(mod in 1:ncol(geneModuleMembership))
{
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(geneInfo0, geneModuleMembership[,modOrder[mod]],
                          MMPvalue[, modOrder[mod]])
  names(geneInfo0) <- c(oldNames, paste("MM.", modNames[modOrder[mod]], sep = ""),
                        paste("p.MM.", modNames[modOrder[mod]], sep = ""))
}

#IPA annotation???????ǉ?
#setwd("E:/IPA_gene_annotation")
#IPAannot <- read.csv("GPL14951_IPAannot.csv", header = T, stringsAsFactors = F)
#IPAannot <- IPAannot[,-c(2:4)]
#geneInfo0 <- merge(geneInfo0, IPAannot, by.x = 1, by.y = 1, all.x = T)

#module color?y??GS?ŕ??בւ?
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.skinscore))
geneInfo <- geneInfo0[geneOrder, ]
#?A?m?e?[?V?????????y??skinscore?Ɋւ???GS/MM?̕\???????o??
#setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
write.csv(geneInfo, paste("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/signed_",deepSplit,"/geneInfo_signed_", deepSplit, "_merged.csv",sep = "" ), row.names = FALSE)
}
#########?????܂?#########


#???ڂ???trait?ƍł????ւ??郂?W???[???i???̏ꍇcyan?j?ɂ????āC??GS??MM?̈??`?q(probe)?????肷??
trait <- "skinscore"
module <- "tan"
column <- match(module, modNames)
moduleGenes <- moduleColors_unsigned == module

#???ڃ??W???[???ɂ?????GS-MM?v???b?g
#pdf(paste("GS_MM_unsigned_", deepSplit, "_", trait, "_", module, ".pdf",sep = ""), width = 7, height = 7)
par(mfrow = c(1, 1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Gene significance for", trait),
                       main = paste("Module membership vs. gene significance", "(unsigned_", deepSplit, ")\n", sep = ""),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
#dev.off()





001_series_matrix_module_trait_relationship_modified.R





rm(list = ls(all.names = TRUE))

#series matrix???p????WGCNA
#???W???[????trait?̊֘A??#
library(WGCNA)
library(tidyverse)

setwd("~/Project/20190709_SystemicSclerosis/3.Script")

Expdata_path<-"../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/Expdata_quontile.Rdata"
gene_annotation_file_path<-"../1.Data/Annotation/GeneAnnotation_GPL10558.csv"
sample_annotation_file_path<-"../1.Data/GSE58095_all_sample_annotatiions.txt"

#?f?t?H???g height=480, width=480, res=72 ?̔{?????p?????Ɨǂ???
#png_resolution<-list(height = 4800, width = 4800, res = 300)
png_resolution<-list(height = 1440, width = 1440, res = 216)
#png_resolution<-c(1440, 1440, 216)

options(stringsAsFactors = FALSE)
lnames <- load(Expdata_path)
datExpr<-Expdata_quontile$E
keep<-rowSums(datExpr>log2(50))>=102
datExpr<-datExpr[keep,]%>%t()
annot <- read.csv(gene_annotation_file_path, header = T, stringsAsFactors = F)

#datTraits?̐??̂??킩?????ȉ??A?T???v???̃A?m?e?[?V?????????p
datTraits<-read.table(sample_annotation_file_path,header = TRUE, skip=4,sep="\t",nrows = 102)
choosen<-c(7,#12,
  13)#,15,16)
datTraits<-datTraits[,choosen] %>% as.data.frame()
names(datTraits)<-c("Sex","Skin_score")


#deepSplit??0?`4??signed unsigned, merged unmerged?̂??ׂĂ̏????????? 
condition_deepsplit<-0:4
condition_signed_unsigned<-c("signed","unsigned")
condition_merged<-c("","_merged")
condition<-list(rep(condition_deepsplit,each=4),rep(condition_signed_unsigned,each=2,5),rep(condition_merged,10))


for(i in 1:length(condition[[1]])){
  deepSplit <- condition[[1]][[i]]
  signed_unsigned<-condition[[2]][[i]]
  merged_unmerged<-condition[[3]][[i]]
  
  base_file_path<-paste("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/",signed_unsigned,"_", deepSplit,"/", sep = "")
  
  networkconstruction_file_path<-paste(base_file_path,"networkConstruction_StepByStep_",signed_unsigned, "_", deepSplit, merged_unmerged,".Rdata", sep = "")
  png_file_path<-paste(base_file_path,"Module_trait_relationships_",signed_unsigned,"_",  deepSplit, merged_unmerged, ".png",sep = "")
  csv_file_path<-paste(base_file_path,"Module_trait_correlation_",signed_unsigned,"_", deepSplit, merged_unmerged, ".csv", sep = "")
  geneinfo_file_path<-paste(base_file_path,"geneInfo_",signed_unsigned, "_", deepSplit, merged_unmerged, ".csv",sep = "" )
  
  lnames <- load(networkconstruction_file_path)
  
  
  #???`?q?i?v???[?u?j???ƃT???v?????????`
  nGenes <- ncol(datExpr)
  nSamples <- nrow(datExpr)
  
  moduleColor<-get(paste("moduleColors_",signed_unsigned,sep = ""))
  
  #ME???Čv?Z
  MEs0 <- moduleEigengenes(datExpr, moduleColor)$eigengenes
  MEs <- orderMEs(MEs0)
  moduleTraitCor <- cor(MEs, datTraits, use = "p",method = "spearman")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  
  #???W???[????trait?̊֘A?????q?[?g?}?b?v?Ŏ???
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) <- dim(moduleTraitCor)
  png(png_file_path,height = png_resolution$height,width = png_resolution$width,res = png_resolution$res)#, width = 10, height = 18)
  #do.call("png",png_resolution)
  par(mar = c(6, 8.5, 3, 3))
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(datTraits),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1, 1),
                 main = paste("Module-trait relationships (",signed_unsigned,", deepSplit ", deepSplit," ",str_sub(merged_unmerged,2), ")", sep = ""))
  dev.off()
  
  ###########################
  moduleTraitCor %>% 
    as.data.frame() %>%
    rownames_to_column("Color") %>% 
    gather(key=cor_to,value=correlation,Sex,Skin_score) %>% 
    arrange(desc(correlation)) %>% 
    #mutate(cor_to = "Skin skore") %>% 
    ggplot(aes(x=cor_to,y=reorder(Color,correlation)))+
    geom_tile(aes(fill=correlation))+
    geom_text(aes(label = round(correlation, 2)))+
    scale_fill_gradient2(low = "green", high = "red") +
    labs(x="",y="")+
    ggtitle(paste("Module-trait relationships (",signed_unsigned,", deepSplit ", deepSplit," ",str_sub(merged_unmerged,2), ")", sep = ""))
  ggsave(paste("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/01_Compare_deepsplit_signed/TraitCorrelation/","Module_trait_relationships_byggplot_",signed_unsigned,"_",  deepSplit, merged_unmerged, ".png",sep = ""),dpi=300)
  
  ###########################
  
  #???W???[????trait?̑??֌W???y??p-value???\?ɂ??ĕۑ?
  colnames(textMatrix) <- names(datTraits)
  rownames(textMatrix) <- names(MEs)
  write.csv(textMatrix, csv_file_path)
  

  #???ڂ???trait?i???̏ꍇ??Fibrosis?j?Ɋւ???Gene significance (GS)??Module membership (MM)???Z?o????
  skinscore <- as.data.frame(datTraits$Skin_score)
  names(skinscore) <-"skinscore"
  modNames <- substring(names(MEs), 3)
  
  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) <- paste("MM", modNames, sep = "")
  names(MMPvalue) <- paste("p.MM", modNames, sep = "")
  
  geneTraitSignificance <- as.data.frame(cor(datExpr, skinscore, use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) <- paste("GS.", names(skinscore), sep = "")
  names(GSPvalue) <- paste("p.GS.", names(skinscore), sep = "")
  
  #Probe ID??gene name?ɕϊ?
  probes <- colnames(datExpr)
  probes2annot <- match(probes, annot$ID)
  #?A?m?e?[?V?????????Ă??Ȃ??v???[?u????0?ł??邩?m?F
  sum(is.na(probes2annot)) 
  
  #?A?m?e?[?V?????????y??skinscore?Ɋւ???GS?̕\???쐬
  geneInfo0 <- data.frame(Probe_ID = probes,
                          geneSymbol = annot$Symbol[probes2annot],
                          Entrez_gene_ID = annot$Entrez_Gene_ID[probes2annot],
                          moduleColor = moduleColor,
                          geneTraitSignificance,
                          GSPvalue)
  #skinscore??p?l?ŕ??בւ?
  modOrder <- order(-abs(cor(MEs, skinscore, use = "p")))
  #MM???????ǉ?
  for(mod in 1:ncol(geneModuleMembership))
  {
    oldNames <- names(geneInfo0)
    geneInfo0 <- data.frame(geneInfo0, geneModuleMembership[,modOrder[mod]],
                            MMPvalue[, modOrder[mod]])
    names(geneInfo0) <- c(oldNames, paste("MM.", modNames[modOrder[mod]], sep = ""),
                          paste("p.MM.", modNames[modOrder[mod]], sep = ""))
  }
  
  #IPA annotation???????ǉ?
  #setwd("E:/IPA_gene_annotation")
  #IPAannot <- read.csv("../1.Data/Annotation/GeneAnnotation_GPL10558.csv", header = T, stringsAsFactors = F)
  #IPAannot <- IPAannot[,-c(2:4)]
  #geneInfo0 <- merge(geneInfo0, IPAannot, by.x = 1, by.y = 1, all.x = T)
  
  #module color?y??GS?ŕ??בւ?
  geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.skinscore))
  geneInfo <- geneInfo0[geneOrder, ]
  
  #?A?m?e?[?V?????????y??skinscore?Ɋւ???GS/MM?̕\???????o??
  #setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, sep = ""))
  write.csv(geneInfo, geneinfo_file_path, row.names = FALSE)
  
}


deepSplit <- 2
signed_unsigned<-"signed"
merged_unmerged<-""

base_file_path<-paste("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/",signed_unsigned,"_", deepSplit,"/", sep = "")

networkconstruction_file_path<-paste(base_file_path,"networkConstruction_StepByStep_",signed_unsigned, "_", deepSplit, merged_unmerged,".Rdata", sep = "")
lname<-load(networkconstruction_file_path)

trait <- "skinscore"
module <- "plum1"

modNames<-moduleColors_signed %>% unique()
column <- match(module, modNames)
moduleGenes <- moduleColors_signed == module

#???ڃ??W???[???ɂ?????GS-MM?v???b?g
png(paste("GS_MM_",signed_unsigned,"_", deepSplit, "_", trait, "_", module, ".png",sep = ""),height = png_resolution$height,width = png_resolution$width,res = png_resolution$res)
par(mfrow = c(1, 1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", trait),
                   main = paste("Module membership vs. gene significance", "(signed_", deepSplit, ")\n", sep = ""),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
############################################################################################################






001_series_matrix_WGCNA.R





#series matrix???p????WGCNA
#?l?b?g???[?N?̍쐬#
setwd("~/Project/20190709_SystemicSclerosis/3.Script")
library(WGCNA)
options(stringsAsFactors = FALSE)

load("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/Expdata_quontile.Rdata")

#WGCNAdata <- log2(Expdata$E)
datExpr <- Expdata_quontile$E
keep<-rowSums(datExpr>log2(50))>=102
datExpr<-datExpr[keep,]%>% t()

#lnames <- load("dataInput.Rdata")

#soft-thresholding Power ???̐ݒ?
powers <- c(c(1:10), seq(from = 12, to = 26, by = 2))

#?l?b?g???[?N?g?|???W?[????
#???ꂼ???̃??ɑ΂???Scale-free fit index??Mean connectivity???Z?o????
#unsigned?i?????ցE?t???ցj?Csigned?i?????ւ̂݁j???ꂼ???Z?o
sft_unsigned <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "unsigned")
sft_signed <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")

#Scale-free fit index??Mean connectivity???}??
sizeGrWindow(9, 5)
pdf("Scale_free_fit_index_and_Mean_connectivity_unsigned.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
cex1 = 0.9
plot(sft_unsigned$fitIndices[,1], -sign(sft_unsigned$fitIndices[,3])*sft_unsigned$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence (unsigned)"))
text(sft_unsigned$fitIndices[,1], -sign(sft_unsigned$fitIndices[,3])*sft_unsigned$fitIndices[,2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red") #R^2?J?b?g?I?t?l0.90?ɐԐ?

plot(sft_unsigned$fitIndices[,1], sft_unsigned$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity (unsigned)"))
text(sft_unsigned$fitIndices[,1], sft_unsigned$fitIndices[,5], labels = powers, cex = cex1, col = "red")
dev.off()

sizeGrWindow(9, 5)
pdf("Scale_free_fit_index_and_Mean_connectivity_signed.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))
cex1 = 0.9
plot(sft_signed$fitIndices[,1], -sign(sft_signed$fitIndices[,3])*sft_signed$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2", type = "n",
     main = paste("Scale independence (signed)"))
text(sft_signed$fitIndices[,1], -sign(sft_signed$fitIndices[,3])*sft_signed$fitIndices[,2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red") #R^2?J?b?g?I?t?l0.90?ɐԐ?

plot(sft_signed$fitIndices[,1], sft_signed$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity (signed)"))
text(sft_signed$fitIndices[,1], sft_signed$fitIndices[,5], labels = powers, cex = cex1, col = "red")
dev.off()

#Scalefree plot??R^2?????Ȃ????œK?ȃ????T??
#unsigned?̏ꍇ
pdf("Scalefree_plot_unsigned.pdf", width = 12, height = 9)
par(mfrow=c(2,2))
for(i in c(c(6:10), seq(from = 12, to = 22, by = 2))){
softPower_unsigned <- i
adjacency_unsigned <- adjacency(datExpr, power = softPower_unsigned, type = "unsigned")
diag(adjacency_unsigned) <- 0
connectivity_unsigned <- apply(adjacency_unsigned, 1, sum)
scaleFreePlot(connectivity_unsigned, truncated = T)
mtext(paste("beta =", i), side = 3, line = 0, adj = 0)
}
dev.off()

#???̑I?ѕ???scale free fit index??0.90?ɒB?????ŏ??l?Cmean connectivity???Ⴗ???Ă͂????Ȃ??Cscale free R^2??1?ɋ߂?
#?? = 14??scale free fit index??0.90?ɒB?????ŏ??l?ł????Cscale free R^2??0.95
#unsigned?̓? = 14???p???ĉ??͂???
softPower_unsigned <- 12
adjacency_unsigned <- adjacency(datExpr, power = softPower_unsigned, type = "unsigned")

#signed?̏ꍇ
pdf("Scalefree_plot_signed.pdf", width = 12, height = 9)
par(mfrow=c(2,2))
for(i in c(seq(from = 12, to = 26, by = 2))){
  softPower_signed <- i
  adjacency_signed <- adjacency(datExpr, power = softPower_signed, type = "signed")
  diag(adjacency_signed) <- 0
  connectivity_signed <- apply(adjacency_signed, 1, sum)
  scaleFreePlot(connectivity_signed, truncated = T)
  mtext(paste("beta =", i), side = 3, line = 0, adj = 0)
}
dev.off()

#?? = 20??scale free fit index??0.90?ɒB?????ŏ??l?ł????Cscale free R^2??0.89
#signed?̓? = 20???p???ĉ??͂???
softPower_signed <- 22
adjacency_signed <- adjacency(datExpr, power = softPower_signed, type = "signed")

#adjacency?i?אڍs???l?j??topological overlap?iTOM?j?ɕϊ?
#unsigned
TOM_unsigned <- TOMsimilarity(adjacency_unsigned, TOMType = "unsigned")
dissTOM_unsigned <- 1-TOM_unsigned

#signed
TOM_signed <- TOMsimilarity(adjacency_signed, TOMType = "signed")
dissTOM_signed <- 1-TOM_signed

#1-TOM???K?w?I?N???X?^?????O
geneTree_unsigned <- hclust(as.dist(dissTOM_unsigned), method = "average")
geneTree_signed <- hclust(as.dist(dissTOM_signed), method = "average")

#?f???h???O???????}??
sizeGrWindow(12, 9)
pdf("Gene_clustering_dissTOM_unsigned.pdf", width = 12, height = 9)
plot(geneTree_unsigned, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity (unsigned)",
     labels = FALSE, hang = 0.04)
dev.off()
pdf("Gene_clustering_dissTOM_signed.pdf", width = 12, height = 9)
plot(geneTree_signed, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity (signed)",
     labels = FALSE, hang = 0.04)
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
ifelse(!dir.exists(paste("~/Project/20190709_SystemicSclerosis/3.Script/unsigned_", deepSplit, sep = "")), 
       dir.create(paste("~/Project/20190709_SystemicSclerosis/3.Script/unsigned_", deepSplit, sep = "")), FALSE)
setwd(paste("~/Project/20190709_SystemicSclerosis/3.Script/unsigned_", deepSplit, sep = ""))

#Dynamic tree cut???p???ă??W???[?????o
dynamicMods_unsigned <- cutreeDynamic(dendro = geneTree_unsigned, distM = dissTOM_unsigned,
                                      deepSplit = deepSplit, pamRespectsDendro = FALSE, pamStage = FALSE,
                                      minClusterSize = minModuleSize)
#???W???[??No.???F?ɕϊ????C?e???W???[???̃v???[?u?????ۑ?
dynamicColors_unsigned <- labels2colors(dynamicMods_unsigned)
module_unsigned <- as.data.frame(table(dynamicColors_unsigned))
write.csv(module_unsigned, paste("module_unsigned_",deepSplit,".csv", sep = ""))

#module eigengene?iME?j?̎Z?o
MEList_unsigned <- moduleEigengenes(datExpr, colors = dynamicColors_unsigned)
MEs_unsigned <- MEList_unsigned$eigengenes

#ME???p???????W???[???̃N???X?^?????O
MEDiss_unsigned <- 1-cor(MEs_unsigned)
METree_unsigned <- hclust(as.dist(MEDiss_unsigned), method = "average")
pdf(paste("Clustering_ME_unsigned_", deepSplit, ".pdf",sep = ""), width = 7, height = 6)


plot(METree_unsigned, main = paste("Clustering of module eigengenes (unsigned, deepSplit ", deepSplit, ")", sep = ""),
     xlab = "", sub = "")




#ME???߂??i?????????Ă????j???W???[?????T??
#ME?f???h???O?????̃J?b?g?I?t?l???ݒ肷??
MEDissThres <- 0.25 #0.15?`0.25?͈̔͂Ńf???h???O???????݂ēK?X????
#?J?b?g???C????????
abline(h = MEDissThres, col = "red")
dev.off()

#?f?[?^??rename???ĕۑ?
moduleColors_unsigned <- dynamicColors_unsigned
colorOrder <- c("grey", standardColors(100))
moduleLabels_unsigned <- match(moduleColors_unsigned, colorOrder)-1
save(MEs_unsigned, moduleLabels_unsigned, moduleColors_unsigned, geneTree_unsigned, 
     file = paste("networkConstruction_StepByStep_unsigned_", deepSplit, ".Rdata", sep = ""))

#ME?̕ۑ?
MEs_unsigned2 <- cbind(rownames(datExpr), MEs_unsigned)
write.csv(MEs_unsigned2, paste("MEs_unsigned_",deepSplit,".csv", sep = ""), row.names = F)

#intramodular connectivity (kWithin)?̎Z?o
kwithin_unsigned <- intramodularConnectivity(adjacency_unsigned, moduleColors_unsigned)

#?v???[?u?Ckwithin?Ccolor?̕\???쐬???C?ۑ?
binded_unsigned <- cbind(rownames(kwithin_unsigned), kwithin_unsigned, moduleColors_unsigned)
colnames(binded_unsigned)[1] <- "ProbeID"
write.csv(binded_unsigned, paste("Probe_kWitnin_color_unsigned_", deepSplit, ".csv", sep=""), row.names=F)


##unsigned_merge????
#?J?b?g?I?t?l?ȉ??̃N???X?^?[??merge
merge_unsigned <- mergeCloseModules(datExpr, dynamicColors_unsigned, cutHeight = MEDissThres, verbose = 3)
mergedColors_unsigned <- merge_unsigned$colors
#merge???̊e???W???[???̃v???[?u?????ۑ?
merged_module_unsigned <- as.data.frame(table(mergedColors_unsigned))
write.csv(merged_module_unsigned, paste("module_unsigned_",deepSplit,"_merged.csv", sep = ""))

#merge???????W???[????ME
mergedMEs_unsigned <- merge_unsigned$newMEs
#?f???h???O??????merge?O???̐F?t?????????W???[?????}??
pdf(paste("Gene_dendrogram_and_module_colors_unsigned_", deepSplit, ".pdf",sep = ""), width = 9, height = 6)
plotDendroAndColors(geneTree_unsigned, cbind(dynamicColors_unsigned, mergedColors_unsigned),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = paste("Gene dendrogram and module colors (unsigned, deepSplit ", deepSplit, ")", sep = ""))
dev.off()

#merge????ME???p???????W???[???̃N???X?^?????O
MEDiss_merged_unsigned <- 1-cor(mergedMEs_unsigned)
METree_merged_unsigned <- hclust(as.dist(MEDiss_merged_unsigned), method = "average")
pdf(paste("Clustering_ME_unsigned_", deepSplit, "_merged.pdf",sep = ""), width = 7, height = 6)
plot(METree_merged_unsigned, main = paste("Clustering of module eigengenes (unsigned, deepSplit ", deepSplit, ", merged)", sep = ""),
     xlab = "", sub = "")
dev.off()

#?f?[?^??rename???ĕۑ?
moduleColors_unsigned <- mergedColors_unsigned
colorOrder <- c("grey", standardColors(100))
moduleLabels_unsigned <- match(moduleColors_unsigned, colorOrder)-1
save(mergedMEs_unsigned, moduleLabels_unsigned, moduleColors_unsigned, geneTree_unsigned, 
     file = paste("networkConstruction_StepByStep_unsigned_", deepSplit, "_merged.Rdata", sep = ""))

#ME?̕ۑ?
mergedMEs_unsigned2 <- cbind(rownames(datExpr), mergedMEs_unsigned)
write.csv(mergedMEs_unsigned2, paste("MEs_unsigned_",deepSplit,"_merged.csv", sep = ""), row.names = F)

#intramodular connectivity (kWithin)?̎Z?o
kwithin_unsigned <- intramodularConnectivity(adjacency_unsigned, moduleColors_unsigned)

#?v???[?u?Ckwithin?Ccolor?̕\???쐬???C?ۑ?
binded_unsigned <- cbind(rownames(kwithin_unsigned), kwithin_unsigned, moduleColors_unsigned)
colnames(binded_unsigned)[1] <- "ProbeID"
write.csv(binded_unsigned, paste("Probe_kWitnin_color_unsigned_", deepSplit, "_merged.csv", sep=""), row.names=F)


##signed_merge?Ȃ?
#?t?H???_?̍쐬?ƈړ?
ifelse(!dir.exists(paste("~/Project/20190709_SystemicSclerosis/3.Script/signed_", deepSplit, sep = "")), 
       dir.create(paste("~/Project/20190709_SystemicSclerosis/3.Script/signed_", deepSplit, sep = "")), FALSE)
setwd(paste("~/Project/20190709_SystemicSclerosis/3.Script/signed_", deepSplit, sep = ""))

#Dynamic tree cut???p???ă??W???[?????o
dynamicMods_signed <- cutreeDynamic(dendro = geneTree_signed, distM = dissTOM_signed,
                                      deepSplit = deepSplit, pamRespectsDendro = FALSE, pamStage = FALSE,
                                      minClusterSize = minModuleSize)
#???W???[??No.???F?ɕϊ????C?e???W???[???̃v???[?u?????ۑ?
dynamicColors_signed <- labels2colors(dynamicMods_signed)
module_signed <- as.data.frame(table(dynamicColors_signed))
write.csv(module_signed, paste("module_signed_",deepSplit,".csv", sep = ""))

#module eigengene?iME?j?̎Z?o
MEList_signed <- moduleEigengenes(datExpr, colors = dynamicColors_signed)
MEs_signed <- MEList_signed$eigengenes

#ME???p???????W???[???̃N???X?^?????O
MEDiss_signed <- 1-cor(MEs_signed)
METree_signed <- hclust(as.dist(MEDiss_signed), method = "average")
pdf(paste("Clustering_ME_signed_", deepSplit, ".pdf",sep = ""), width = 7, height = 6)
plot(METree_signed, main = paste("Clustering of module eigengenes (signed, deepSplit ", deepSplit, ")", sep = ""),
     xlab = "", sub = "")
#ME???߂??i?????????Ă????j???W???[?????T??
#ME?f???h???O?????̃J?b?g?I?t?l???ݒ肷??
MEDissThres <- 0.25 #0.15?`0.25?͈̔͂Ńf???h???O???????݂ēK?X????
#?J?b?g???C????????
#abline(h = MEDissThres, col = "red")
dev.off()

#?f?[?^??rename???ĕۑ?
moduleColors_signed <- dynamicColors_signed
colorOrder <- c("grey", standardColors(100))
moduleLabels_signed <- match(moduleColors_signed, colorOrder)-1
save(MEs_signed, moduleLabels_signed, moduleColors_signed, geneTree_signed, 
     file = paste("networkConstruction_StepByStep_signed_", deepSplit, ".Rdata", sep = ""))

#ME?̕ۑ?
MEs_signed2 <- cbind(rownames(datExpr), MEs_signed)
write.csv(MEs_signed2, paste("MEs_signed_",deepSplit,".csv", sep = ""), row.names = F)

#intramodular connectivity (kWithin)?̎Z?o
kwithin_signed <- intramodularConnectivity(adjacency_signed, moduleColors_signed)

#?v???[?u?Ckwithin?Ccolor?̕\???쐬???C?ۑ?
binded_signed <- cbind(rownames(kwithin_signed), kwithin_signed, moduleColors_signed)
colnames(binded_signed)[1] <- "ProbeID"
write.csv(binded_signed, paste("Probe_kWitnin_color_signed_", deepSplit, ".csv", sep=""), row.names=F)


##signed_merge????
#?J?b?g?I?t?l?ȉ??̃N???X?^?[??merge
merge_signed <- mergeCloseModules(datExpr, dynamicColors_signed, cutHeight = MEDissThres, verbose = 3)
mergedColors_signed <- merge_signed$colors
#merge???̊e???W???[???̃v???[?u?????ۑ?
merged_module_signed <- as.data.frame(table(mergedColors_signed))
write.csv(merged_module_signed, paste("module_signed_",deepSplit,"_merged.csv", sep = ""))

#merge???????W???[????ME
mergedMEs_signed <- merge_signed$newMEs
#?f???h???O??????merge?O???̐F?t?????????W???[?????}??
pdf(paste("Gene_dendrogram_and_module_colors_signed_", deepSplit, ".pdf",sep = ""), width = 9, height = 6)
plotDendroAndColors(geneTree_signed, cbind(dynamicColors_signed, mergedColors_signed),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = paste("Gene dendrogram and module colors (signed, deepSplit ", deepSplit, ")", sep = ""))
dev.off()

#merge????ME???p???????W???[???̃N???X?^?????O
MEDiss_merged_signed <- 1-cor(mergedMEs_signed)
METree_merged_signed <- hclust(as.dist(MEDiss_merged_signed), method = "average")
pdf(paste("Clustering_ME_signed_", deepSplit, "_merged.pdf",sep = ""), width = 7, height = 6)
plot(METree_merged_signed, main = paste("Clustering of module eigengenes (signed, deepSplit ", deepSplit, ", merged)", sep = ""),
     xlab = "", sub = "")
dev.off()

#?f?[?^??rename???ĕۑ?
moduleColors_signed <- mergedColors_signed
colorOrder <- c("grey", standardColors(100))
moduleLabels_signed <- match(moduleColors_signed, colorOrder)-1
save(mergedMEs_signed, moduleLabels_signed, moduleColors_signed, geneTree_signed, 
     file = paste("networkConstruction_StepByStep_signed_", deepSplit, "_merged.Rdata", sep = ""))

#ME?̕ۑ?
mergedMEs_signed2 <- cbind(rownames(datExpr), mergedMEs_signed)
write.csv(mergedMEs_signed2, paste("MEs_signed_",deepSplit,"_merged.csv", sep = ""), row.names = F)

#intramodular connectivity (kWithin)?̎Z?o
kwithin_signed <- intramodularConnectivity(adjacency_signed, moduleColors_signed)

#?v???[?u?Ckwithin?Ccolor?̕\???쐬???C?ۑ?
binded_signed <- cbind(rownames(kwithin_signed), kwithin_signed, moduleColors_signed)
colnames(binded_signed)[1] <- "ProbeID"
write.csv(binded_signed, paste("Probe_kWitnin_color_signed_", deepSplit, "_merged.csv", sep=""), row.names=F)
}
#########?????܂?#########







002_Compare_expression_dif_vs_lim.R





library(tidyverse)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")


# check SUZ12 expression control vs lSSc vs dSSc

Expdat <- read_csv("../1.Data/GSE58095_main/Expdata_GSE58095.csv")
sample_annotation <- read_tsv("../1.Data/GSE58095_main/GSE58095_all_sample_annotations.txt", skip = 4, n_max = 102)
gene_annotation <- read_csv("../1.Data/GSE58095_main/Annotation/GeneAnnotation_GPL10558.csv")

target_gene = "SUZ12"

target_probeid <- gene_annotation %>% 
  filter(Symbol == target_gene) %>% 
  pull(ID)

Expdat %>% 
  filter(ProbeID == target_probeid)


sample_annotation_Early <- sample_annotation %>% 
  #drop_na(`Characteristics: Group`, `Characteristics: Diffuse vs limited`)
  filter(`Characteristics: Time point` =="Early") %>% 
  rename(Patient=`Characteristics: Patient`,
         Group=`Characteristics: Group`,
         Dif_vs_lim=`Characteristics: Diffuse vs limited`)


group_list <- sample_annotation_Early$Group
names(group_list) <- sample_annotation_Early$Patient
Dif_vs_lim_list <- sample_annotation_Early$Dif_vs_lim
names(Dif_vs_lim_list) <- sample_annotation_Early$Patient

my_group<-sample_annotation %>% 
  mutate(`Characteristics: Group`=if_else(!is.na(`Characteristics: Time point`) & `Characteristics: Time point`=="Late",
                                          as.character(group_list[`Characteristics: Patient`]),
                                          `Characteristics: Group`),
         `Characteristics: Diffuse vs limited`=if_else(!is.na(`Characteristics: Time point`) & `Characteristics: Time point`=="Late",
                                                       as.character(Dif_vs_lim_list[`Characteristics: Patient`]),
                                                       `Characteristics: Diffuse vs limited`)) %>% 
  mutate(my_group=if_else(`Characteristics: Group`=="SSc",
                        str_c(`Characteristics: Group`, "_", `Characteristics: Diffuse vs limited`),
                        `Characteristics: Group`) )%>% 
  pull(my_group)
  
  

Expdat %>% 
  filter(ProbeID == target_probeid) %>% 
  gather(key="Patient",
         value="expression",-ProbeID) %>% 
  mutate(Group=my_group) %>% 
  group_by(Group) %>% 
  summarise(mean = mean(expression))
  






002_MakeModuleSummary.R





library(tidyverse)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")


# read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv") %>% 
#   dplyr::select(ProbeID ,kTotal,kWithin,kOut,kDiff) %>% 
#   left_join(read_csv("../2.Output/02_signed_2_unmerged/logFC_ModuleColor.csv") ,by=c("ProbeID"="ILMNID")) %>% 
#   write_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv")



read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv") %>% 
  filter(module_color=="lightyellow") %>% 
  arrange(desc(kWithin)) %>% 
  mutate(dense_rank=dense_rank(kWithin))->a
  filter(dense_rank(kWithin)>0.2)
  
kWithinRank_Lightyellow<-read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv") %>% 
  filter(module_color=="lightyellow") %>% 
  arrange(desc(kWithin)) %>% 
  mutate(
    percent_rank=percent_rank(desc(kWithin)),
    dense_rank=dense_rank(desc(kWithin)))
  
SSc_variant<-read_csv("../2.Output/12_Comparison_withGWAS/1.Data/SSc_variant.csv") %>% 
  distinct(entrez_id) %>% .$entrez_id

kWithinRank_Lightyellow %>% 
  mutate(is_GWASgene= Entrez_Gene_ID %in% SSc_variant) %>% 
  write_csv("../2.Output/24_Lightyellow_BN/Lightyellow_kWithinGWAS.csv")


significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor



geneid2KwithinModulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv")

geneid2KwithinModulecolor %>% 
  filter(module_color=="yellow") %>% 
  arrange(desc(kWithin))%>% 
  mutate(percent_rank=percent_rank(desc(kWithin)),
    dense_rank=dense_rank(desc(kWithin)),
    is_GWASgene= Entrez_Gene_ID %in% SSc_variant) %>% 
  write_csv("../2.Output/02_signed_2_unmerged/Yellow_summary.csv")
  
color="lightcyan"


geneid2KwithinModulecolor %>% 
  filter(module_color==color) %>% 
  arrange(desc(kWithin))%>% 
  mutate(percent_rank=percent_rank(desc(kWithin)),
         dense_rank=dense_rank(desc(kWithin)),
         is_GWASgene= Entrez_Gene_ID %in% SSc_variant) %>% 
  write_csv(paste0("../2.Output/02_signed_2_unmerged/ModuleSummary/",color,"_summary.csv"))

all_module <- read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  .$ModuleColor



for( color in all_module){
  
  geneid2KwithinModulecolor %>% 
    filter(module_color==color) %>% 
    arrange(desc(kWithin))%>% 
    mutate(percent_rank=percent_rank(desc(kWithin)),
           dense_rank=dense_rank(desc(kWithin)),
           is_GWASgene= Entrez_Gene_ID %in% SSc_variant) %>% 
    write_csv(paste0("../2.Output/02_signed_2_unmerged/ModuleSummary/",color,"_summary.csv"))
}







004_EnrichmentAnalysis_by_ImmuneCellCollection.R





#m(list = ls(all.names = TRUE))
library(WGCNA)
library(anRichment)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
options(stringsAsFactors = FALSE)

##############???͂Ɏg?p???????`?q?R???N?V???????쐬????
#?܂??֋X???̋?Collection??????
collection <- newCollection()

#Collection?̂??ƂɂȂ??e???`?q???X?g?i?????C?w?b?_?[?????ň??`?q??symbol?@entrez?̏ꍇ?͕ϊ??s?v?Ȃ̂ŉ??L?X?N???v?g?ύX)
#?t?H???_????Genelist???ꊇ?œǂݍ??܂??邽?߂̐ݒ?
#genelist?̃t?@?C?????????̂܂?enrichment???͎??̈??`?q?O???[?v???ɔ??f???????̂ŁC?K?v?ł????΂??ꂼ?ꃊ?l?[??
setwd("E:/01_Analysis/01_INF/05_Cell Specific gene/human immune cell/pSI/P??0.1/") #?t?@?C???i?[?p?X???w??
files <- list.files(path = "E:/01_Analysis/01_INF/05_Cell Specific gene/human immune cell/pSI/P??0.1/", full.names = T)
files_names <- list.files("E:/01_Analysis/01_INF/05_Cell Specific gene/human immune cell/pSI/P??0.1/", full.names = F, pattern="csv") 
files_names <-  gsub(".csv", "", files_names) 
files_names
cell_names <- str_replace(files_names, pattern="_expression_gene", replacement="")

#genelist?̗R???ƂȂ?Organism???w?肷??
organism <- "Human"

#Genelist??Geneset?ɂ??邽?߂̊֐????ǂݍ??܂???
Togeneset <- function (X, Y){
  colnames(X) <- "Gene" #???????????I??Gene?ɕύX
  Symbol = X$Gene  #?֋X???̃??l?[??
  #????entrez?`???̈??`?q???????????ꍇ?́C???L??Entrez = X$Gene?ɂ??āC??2?s?̃X?N???v?g?폜
  Entrez.0 = convert2entrez(organism = organism, symbol = Symbol) 
  Entrez = unique(Entrez.0)#???`?q??(Symbol)??entrez?`???ɕϊ?????
  print(table(is.finite(Entrez)))  #entrez?`???ɕϊ??ł??????`?q???̊m?F
  newGeneSet(
    geneEntrez = Entrez,
    geneEvidence = "IEP",
    geneSource = "",
    ID = Y,
    name = Y,
    description = Y,
    source = "",
    organism = organism,
    internalClassification = Y,
    groups = Y,
    lastModified = "")
  
}

#?t?H???_????Genelist???S??Togeneset?֐???geneset?ɕϊ????Ă??????L??Collection?ɓ?????
for(i in 1:length(files)){
  genes <- read.csv(files[i], fileEncoding="UTF-8")#?????????t?@?C?????ǂݍ???
  genes2<- genes$Symbol
  assign(cell_names[i], genes2)
  genes <- unique(genes)#?d???L???ꍇ????????
  genes_set <- Togeneset(genes, cell_names[i])#genelist??geneset??
  assign(paste(cell_names[i]), genes_set)
  collection <- addToCollection(collection, genes_set)#Geneset?????ꂽ??collection?ɓ?????
  #?????L??"Human_CNS_5cells_mouse_collection"?͎????̂??肽??collection?????ݒ?
  rm(genes, genes_set)
}

#doublecollecion <- newCollection()
#doublecollecion <- addToCollection(doublecollecion, Human_CNS_5cells_collection,Human_INF_36cells_collection_N)

##?K?v?ł????Ώ??L?ō쐬?????R???N?V???????ۑ??B
#setwd("C:/Users/mjd9761/Documents/Enrichment/celltype_enrichment/Genecollectiondata")
#save(collection, file="Human_CNS_5cells_collection.Rdata")
#???x?R???N?V???????????΁C?????????̓??[?h???Ă??ł??g????
#load("Human_CNS_5cells_collection.Rdata")

######?I?v?V????
#???????R???N?V?????????????????????^?R???N?V???????쐬???āC?ꊇ?ŃG?????b?`?????g???͂??\?B
#CNS?̈??`?q???X?g?Ɖ??ǍזE?̈??`?q???X?g???ʁX?ō쐬???????ǁC?????ɃG?????b?`?????g???͂??????Ƃ??ȂǁB
#?ȉ??????̗?
#CNS_INF_collecion <- newCollection()
#CNS_INF_collecion <- addToCollection(doublecollecion, Human_CNS_5cells_collection,Human_INF_36cells_collection_N)

#enrichment???͂Ɏg?p????collection?Ɋ܂܂??????`?q???X?g?????o?^????
#???Lenrichment?֐??ɂ????āCTOP?????܂ł?p-value???ʂ??Z?o???邩?ɉe???B?C?ӂŕύX??
nlist_collection <- length(collection)

#?זEenrichment???͂Ńq?[?g?}?b?v???쐬?????ۂɍזE?̕??я????ݒ?
cellorder <- c("Colony.Forming.Unit.Granulocyte", "Colony.Forming.Unit.Megakaryocytic", "Colony.Forming.Unit.Monocyte", "Common.myeloid.progenitor",
               "Basophil", "Eosinophil", "Monocyte", "Neutrophil", "Neutrophilic.Metamyelocyte", "mDC", "pDC", "Granulocyte.monocyte.progenitor",
               "Mature.NKcell.CD56negaCD16negaCD3nega", "Mature.NKcell.CD56negaCD16posiCD3nega", "Mature.NKcell.CD56posiCD16posiCD3nega",
               "Pro.Bcell", "Early.Bcell", "Naive.Bcell", "Mature.Bcell", "Mature.Bcell.class.able.to.switch", "Mature.Bcell.class.switched",
               "Naive.CD4posi.Tcell", "Naive.CD8posi.Tcell", "CD4posi.Tcm", "CD4posi.Tem", "CD8posi.Tcm", "CD8posi.Tem", "CD8posi.Tem_RA", "NKT",
               "Erythroid.CD34negaCD71lowGlyAposi", "Erythroid.CD34negaCD71negaGlyAposi", "Erythroid.CD34negaCD71posiGlyAnega", 
               "Erythroid.CD34negaCD71posiGlyAposi", "Erythroid.CD34posiCD71posiGlyAnega", "Megakaryocyte.erythroid.progenitor",
               "Megakaryocyte", "HSC.CD133posiCD34dim", "HSC.CD38negaCD34posi")

#CelltypeEnrichment?֐????ǂݍ??܂???
CelltypeEnrichment <- function (X){
  names(X) <- c("Gene", "Group") #???????????I??Gene??Group?ɕύX
  symbol = X$Gene  #?֋X???̃??l?[??
  Group = X$Group?@#?֋X???̃??l?[??
  print(table(Group))  #???`?q???X?g?̐??̊m?F?p?B
  entrez = convert2entrez(organism = "human", symbol = symbol)?@#???`?q??(Symbol)??entrez?`???ɕϊ?????
  na.omit(entrez)
  print(table(is.finite(entrez)))  #entrez?`???ɕϊ??ł??????`?q???̊m?F
  #????????anRichment?p?b?P?[?W?̒???enrichmentAnalysis?֐????g?p
  analysisresult = enrichmentAnalysis(?@?@
    classLabels = Group, identifiers = entrez,?@#classLabels???????̂??̂??????̉??͑Ώۈ??`?q?Q?Ƃ??ĔF???B
    refCollection = collection,  #reference?Ƃ??Đݒ肷???f?[?^?Z?b?g?B
    useBackground = "given",?@?@#???̓o?b?N?O???E???h?̐ݒ??i?d?v?j?B?????̓C???v?b?g?????S???`?q?̂????ǂݍ??߂????̂??g?p?B
    threshold = 1,?@?@#?G?????b?`?????g???͌??ʂ??Ƃ??ďo?͂?????csv?t?@?C???ł́Cp?l?̂??????l?B?????͍L???Ƃ??Ă????B
    thresholdType = "Bonferroni",?@#???Lp?l?ɕt?????āCBonferroni?␳????p?l??臒l?Ƃ???
    getOverlapEntrez = FALSE,?@?@#?o??csv?t?@?C???ɂāC?I?[?o?[???b?v???????`?q????entrez?`???ŏo?͂????C
    getOverlapSymbols = TRUE,    #?o??csv?t?@?C???ɂāC?I?[?o?[???b?v???????`?q????symbol?`???ŏo?͂????C
    maxReportedOverlapGenes = 10000,?@?@#???L?I?[?o?[???b?v???????`?q???ǂ̂??炢?\???????邩???ݒ??B?S?Ăق????̂?10,000??
    removeDuplicatesInDifferentClasses =FALSE,?@#???????W???[?????ɓ??????`?q?????????ꍇ???C???̂܂܏???????
    entrySeparator = ",",?@?@#?I?[?o?[???b?v???????`?q?Q?ɂ??āC","?ŕ????ďo?͂??Ă????B?ǂ??ł??悵
    ignoreLabels = "grey", #grey???W???[???̈??`?q?Q?ɂ??Ă͖????i???͂??Ȃ??j
    combineEnrichmentTables = FALSE) 
  
  
  countsInDataSet<- analysisresult$effectiveBackgroundSize  #???̓o?b?N?O???E???h?̈??`?q?????m??
  print(table(countsInDataSet))?@?@?@?@?@?@?@?@?@?@?@?@?@?@?@?@#???????o??
  
  
  Resulttable <- analysisresult$enrichmentTable
  Resulttable <- separate(Resulttable,overlapGenes, into=as.character(c(1:1000)), sep=",")
  Resulttable2 <- subset(analysisresult$enrichmentTable, analysisresult$enrichmentTable$pValue < 0.05) #?o?͗p?ɐ??`
  options(warn=-1) #???̃R?[?h??warning???o???B???????????????R?[?h
  #Overlap???`?q?Q???C1???`?q???ƂɃZ???ɕ????ĕ\???????邽?߂̃R?[?h(1000?͔C?ӁB)???̍??Ƃ?warning???o?邪???Ə??͖????Ȃ?
  write.table(Resulttable2, file = "WGCNA_enrichment_result.csv",row.names = FALSE,sep=",") #???ʂ?csv?t?@?C???ɏo??
  list <- by(Resulttable, Resulttable$class, data.frame) #?O???[?v???ƂɃt?@?C???`???????X?g??
  sapply(1:dim(list), function(x){write.csv(list[x], file=paste0("Result_", dimnames(list)[[1]][x], ".csv"), row.names=FALSE)})?@#?O???[?v???ƂɌ??ʂ??o??
  
  #?????????C???͌??ʃf?[?^(Resulttable)?̂????O???t?ɕK?v?Ȃ??̂??Ԃ????ʂ?
  Rank_all <- sapply(1:dim(list), function(x){list(list[[x]][1:nlist_collection,c(4,6,9)])})#?O???[?v???ƂɕK?v?ȗ????????o??
  names(Rank_all) <- sapply(1:dim(list), function(x){names(Rank_all) <- names(list[x])})#???L?ŃO???[?v???????????̂ł??????????l?[??
  for (i in 1:dim(list)) Rank_all[[i]]$nCommonGenes <- paste("(",Rank_all[[i]]$nCommonGenes,")") #?O???t?̃??x???p?̖??O???`
  for (i in 1:dim(list)) Rank_all[[i]]$pValue <- -(log10(as.numeric(Rank_all[[i]]$pValue)))?@#?O???t?̃??x???p?̖??O???`
  Rank_all <- na.omit(Rank_all) #???ʂ?TOP10?ɖ????Ȃ??ꍇ??NA????
  for (i in 1:dim(list)) Rank_all[[i]] <- transform(Rank_all[[i]], "Rename"=(paste(Rank_all[[i]]$dataSetName,Rank_all[[i]]$nCommonGenes,sep="?@")))#?O???t?̃??x???p?̖??O???`
  
  #list Rename
  names(Rank_all) <- sapply(1:dim(list), function(x){names(Rank_all) <- names(list[x])})
  
  
  ##????????ggplot2???g?????}???̎w??
  for (i in 1:dim(list)) ggsave(file=paste0("Result_", names(Rank_all)[i], ".pdf"),plot = ggplot(Rank_all[[i]], aes(x=reorder(Rank_all[[i]]$Rename, Rank_all[[i]]$pValue), y=Rank_all[[i]]$pValue)) +  
                                  geom_bar(stat="identity", width=.5,fill="black") +   
                                  coord_flip() +                                     
                                  xlab("Cell-type\n(#Overlap genes)") + 
                                  ylab("-log(P-value)"))
  
  #Cell-type enrichment p-value?̃q?[?g?}?b?v???}??
  CellTyperesult <- read.csv("WGCNA_enrichment_result.csv")
  logP <- -log10(CellTyperesult$pValue)
  CellTypeP <- data.frame(Module = CellTyperesult$class, CellType = CellTyperesult$dataSetID, logP = logP)
  ghm <- ggplot(CellTypeP, aes(x = CellType, y = Module, fill = logP))
  ghm <- ghm + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  ghm <- ghm + scale_x_discrete(limits = cellorder)
  ghm <- ghm + geom_tile(aes(fill = logP))
  ghm <- ghm + xlab("Cell Type") + ylab("Module")
  ghm <- ghm + geom_text(aes(label = round(CellTypeP$logP, 0)), size = 1)
  ghm <- ghm + scale_fill_gradient(low = "white", high = "red", limits = c(0, 120))#
  pdf("CellType_pvalue_heatmap.pdf", width = 7, height = 7)
  plot(ghm)
  dev.off()
  
}

###################?֐??????܂?###################

#???ۂɉ???
setwd("C:/Users/ana1484/Desktop/shiota/")
files <- list.files(path = "//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/tmp/MakeWGCNAlist/WGCNAList", full.names = T)
files_names <- list.files(path = "//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/tmp/MakeWGCNAlist/WGCNAList", full.names = F, pattern="csv") 
files_names <-  gsub(".csv", "", files_names) 
#files_names

for (i in 1:length(files)){
  #print(file)
  setwd("C:/Users/ana1484/Desktop/shiota/")
  filepath<-files[[i]]
  dirname<-files_names[[i]]
  if (!file.exists(paste("./",dirname,sep = ""))){
    dir.create(paste("./",dirname, sep = ""))
    setwd(paste("C:/Users/ana1484/Desktop/shiota/",dirname, sep = ""))
    WGCNAgenelist <- read.csv(file = filepath, header = TRUE)
    CelltypeEnrichment(WGCNAgenelist)
  }
  #print(head(WGCNAgenelist))
  #WGCNAgenelist2 <- read.table("WGCNA_gene_module.txt", header = TRUE)
}




#save(collection, file = "Human_Immune_Cell_Collection.RData")






004_EnrichmentAnalysis_by_SkinAndBloodCollection.R





rm(list = ls(all.names = TRUE))
library(WGCNA)
library(anRichment)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
options(stringsAsFactors = FALSE)

load("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/tmp/MakeWGCNAlist/Human_skin_blood_CellType_collection.Rdata")

nlist_collection <- length(collection)

#CelltypeEnrichment?֐????ǂݍ??܂???
CelltypeEnrichment <- function (X){
  names(X) <- c("Gene", "Group") #???????????I??Gene??Group?ɕύX
  symbol = X$Gene  #?֋X???̃??l?[??
  Group = X$Group?@#?֋X???̃??l?[??
  print(table(Group))  #???`?q???X?g?̐??̊m?F?p?B
  entrez = convert2entrez(organism = "human", symbol = symbol)?@#???`?q??(Symbol)??entrez?`???ɕϊ?????
  na.omit(entrez)
  print(table(is.finite(entrez)))  #entrez?`???ɕϊ??ł??????`?q???̊m?F
  #????????anRichment?p?b?P?[?W?̒???enrichmentAnalysis?֐????g?p
  analysisresult = enrichmentAnalysis(?@?@
    classLabels = Group, identifiers = entrez,?@#classLabels???????̂??̂??????̉??͑Ώۈ??`?q?Q?Ƃ??ĔF???B
    refCollection = collection,  #reference?Ƃ??Đݒ肷???f?[?^?Z?b?g?B
    useBackground = "given",?@?@#???̓o?b?N?O???E???h?̐ݒ??i?d?v?j?B?????̓C???v?b?g?????S???`?q?̂????ǂݍ??߂????̂??g?p?B
    threshold = 1,?@?@#?G?????b?`?????g???͌??ʂ??Ƃ??ďo?͂?????csv?t?@?C???ł́Cp?l?̂??????l?B?????͍L???Ƃ??Ă????B
    thresholdType = "Bonferroni",?@#???Lp?l?ɕt?????āCBonferroni?␳????p?l??臒l?Ƃ???
    getOverlapEntrez = FALSE,?@?@#?o??csv?t?@?C???ɂāC?I?[?o?[???b?v???????`?q????entrez?`???ŏo?͂????C
    getOverlapSymbols = TRUE,    #?o??csv?t?@?C???ɂāC?I?[?o?[???b?v???????`?q????symbol?`???ŏo?͂????C
    maxReportedOverlapGenes = 10000,?@?@#???L?I?[?o?[???b?v???????`?q???ǂ̂??炢?\???????邩???ݒ??B?S?Ăق????̂?10,000??
    removeDuplicatesInDifferentClasses =FALSE,?@#???????W???[?????ɓ??????`?q?????????ꍇ???C???̂܂܏???????
    entrySeparator = ",",?@?@#?I?[?o?[???b?v???????`?q?Q?ɂ??āC","?ŕ????ďo?͂??Ă????B?ǂ??ł??悵
    ignoreLabels = "grey", #grey???W???[???̈??`?q?Q?ɂ??Ă͖????i???͂??Ȃ??j
    combineEnrichmentTables = FALSE) 
  
  
  countsInDataSet<- analysisresult$effectiveBackgroundSize  #???̓o?b?N?O???E???h?̈??`?q?????m??
  print(table(countsInDataSet))?@?@?@?@?@?@?@?@?@?@?@?@?@?@?@?@#???????o??
  
  
  Resulttable <- analysisresult$enrichmentTable
  Resulttable <- separate(Resulttable,overlapGenes, into=as.character(c(1:1000)), sep=",")
  Resulttable2 <- subset(analysisresult$enrichmentTable, analysisresult$enrichmentTable$pValue < 0.05) #?o?͗p?ɐ??`
  options(warn=-1) #???̃R?[?h??warning???o???B???????????????R?[?h
  #Overlap???`?q?Q???C1???`?q???ƂɃZ???ɕ????ĕ\???????邽?߂̃R?[?h(1000?͔C?ӁB)???̍??Ƃ?warning???o?邪???Ə??͖????Ȃ?
  write.table(Resulttable2, file = "WGCNA_enrichment_result.csv",row.names = FALSE,sep=",") #???ʂ?csv?t?@?C???ɏo??
  list <- by(Resulttable, Resulttable$class, data.frame) #?O???[?v???ƂɃt?@?C???`???????X?g??
  sapply(1:dim(list), function(x){write.csv(list[x], file=paste0("Result_", dimnames(list)[[1]][x], ".csv"), row.names=FALSE)})?@#?O???[?v???ƂɌ??ʂ??o??
  
  #?????????C???͌??ʃf?[?^(Resulttable)?̂????O???t?ɕK?v?Ȃ??̂??Ԃ????ʂ?
  Rank_all <- sapply(1:dim(list), function(x){list(list[[x]][1:nlist_collection,c(4,6,9)])})#?O???[?v???ƂɕK?v?ȗ????????o??
  names(Rank_all) <- sapply(1:dim(list), function(x){names(Rank_all) <- names(list[x])})#???L?ŃO???[?v???????????̂ł??????????l?[??
  for (i in 1:dim(list)) Rank_all[[i]]$nCommonGenes <- paste("(",Rank_all[[i]]$nCommonGenes,")") #?O???t?̃??x???p?̖??O???`
  for (i in 1:dim(list)) Rank_all[[i]]$pValue <- -(log10(as.numeric(Rank_all[[i]]$pValue)))?@#?O???t?̃??x???p?̖??O???`
  Rank_all <- na.omit(Rank_all) #???ʂ?TOP10?ɖ????Ȃ??ꍇ??NA????
  for (i in 1:dim(list)) Rank_all[[i]] <- transform(Rank_all[[i]], "Rename"=(paste(Rank_all[[i]]$dataSetName,Rank_all[[i]]$nCommonGenes,sep="?@")))#?O???t?̃??x???p?̖??O???`
  
  #list Rename
  names(Rank_all) <- sapply(1:dim(list), function(x){names(Rank_all) <- names(list[x])})
  
  
  ##????????ggplot2???g?????}???̎w??
  for (i in 1:dim(list)) ggsave(file=paste0("Result_", names(Rank_all)[i], ".pdf"),plot = ggplot(Rank_all[[i]], aes(x=reorder(Rank_all[[i]]$Rename, Rank_all[[i]]$pValue), y=Rank_all[[i]]$pValue)) +  
                                  geom_bar(stat="identity", width=.5,fill="black") +   
                                  coord_flip() +                                     
                                  xlab("Cell-type\n(#Overlap genes)") + 
                                  ylab("-log(P-value)"))
  
  #Cell-type enrichment p-value?̃q?[?g?}?b?v???}??
  CellTyperesult <- read.csv("WGCNA_enrichment_result.csv")
  logP <- -log10(CellTyperesult$pValue)
  CellTypeP <- data.frame(Module = CellTyperesult$class, CellType = CellTyperesult$dataSetID, logP = logP)
  ghm <- ggplot(CellTypeP, aes(x = CellType, y = Module, fill = logP))
  ghm <- ghm + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  #ghm <- ghm + scale_x_discrete(limits = cellorder)
  ghm <- ghm + geom_tile(aes(fill = logP))
  ghm <- ghm + xlab("Cell Type") + ylab("Module")
  ghm <- ghm + geom_text(aes(label = round(CellTypeP$logP, 0)), size = 1)
  ghm <- ghm + scale_fill_gradient(low = "white", high = "red", limits = c(0, 120))
  pdf("CellType_pvalue_heatmap.pdf", width = 7, height = 7)
  plot(ghm)
  dev.off()
  
}

###################?֐??????܂?###################

#???ۂɉ???
setwd("C:/Users/ana1484/Desktop/shiota/")
files <- list.files(path = "//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/tmp/MakeWGCNAlist/WGCNAList", full.names = T)
files_names <- list.files(path = "//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/tmp/MakeWGCNAlist/WGCNAList", full.names = F, pattern="csv") 
files_names <-  gsub(".csv", "", files_names) 
#files_names

for (i in 1:length(files)){
  #print(file)
  setwd("C:/Users/ana1484/Desktop/shiota/")
  filepath<-files[[i]]
  dirname<-paste("C:/Users/ana1484/Desktop/shiota/",files_names[[i]],"_Skin",sep = "")
  if (!file.exists(dirname)){
    dir.create(dirname)
    setwd(dirname)
    WGCNAgenelist <- read.csv(file = filepath, header = TRUE)
    CelltypeEnrichment(WGCNAgenelist)
  }
  #print(head(WGCNAgenelist))
  #WGCNAgenelist2 <- read.table("WGCNA_gene_module.txt", header = TRUE)
}







004_Make_Human_Skin_Blood_Cell_collection.R





setwd("~/Project/20190709_SystemicSclerosis/MakeWGCNAlist")
library(WGCNA)
library(anRichment)
library(tidyverse)
library(stringr)

collection <- newCollection()

#?t?H???_????Genelist???ꊇ?œǂݍ??܂??邽?߂̐ݒ?
#?t?@?C???i?[?p?X???w??
EnrichmentCollectionPath<-"./Cell_Enrivhment/Skin+Blood/Jensen_TISSUES/"
files <- list.files(path = EnrichmentCollectionPath,full.names = T)
files_names <- list.files(EnrichmentCollectionPath , full.names = F, pattern="csv")
files_names <- gsub(".csv", "", files_names) 

#Genelist??Geneset?ɂ??邽?߂̊֐????ǂݍ??܂???
Togeneset <- function (X, Y){
  #Organism???w?肷??
  organism <- "Human"
  colnames(X) <- "Gene" #???????????I??Gene?ɕύX
  Symbol = X$Gene  #?֋X???̃??l?[??
  Entrez.0 = convert2entrez(organism = organism, symbol = Symbol)
  Entrez = unique(Entrez.0)#???`?q??(Symbol)??entrez?`???ɕϊ?????
  print(table(is.finite(Entrez)))  #entrez?`???ɕϊ??ł??????`?q???̊m?F
  newGeneSet(
    geneEntrez = Entrez,
    geneEvidence = "IEP",
    geneSource = "",
    ID = Y,
    name = Y,
    description = Y,
    source = "",
    organism = organism,
    internalClassification = Y,
    groups = Y,
    lastModified = "")
}
#?convert2entrez

#?t?H???_????Genelist???ꊇ??geneset?ɕϊ????Ă??????L??Collection?ɓ?????
for(i in 1:length(files)){
  genes_set <- read_csv(files[i]) %>% 
    .[1] %>% 
    unique() %>% 
    Togeneset(files_names[i])#genelist??geneset??
  collection <- addToCollection(collection, genes_set)#Geneset?????ꂽ??collection?ɓ?????
  rm(genes_set)
}

#???`?q?Z?b?g?ۑ?
save(collection, file="Human_skin_blood_CellType_collection.Rdata")






005_BN_BMA.R





library(tidyverse)
library(tidygraph)
library(ggraph)
library(igraph)
library(bnlearn)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

#BMA
#c2PB
BMA_pre=read.table("../2.Output/05_BayesianNetwork/BN_SSc/result/BMA_pre_re.csv",sep=",",header = F,row.names = NULL,stringsAsFactors = F)
data=read.table("../2.Output/05_BayesianNetwork/ME_SkinScore_forBN.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
rownames(BMA_pre)=rownames(data)
colnames(BMA_pre)=rownames(data)
for (i in 1:nrow(BMA_pre)){
  for (j in 1:ncol(BMA_pre)){
    if (BMA_pre[i,j]>30 && BMA_pre[j,i]>30){
      if (BMA_pre[i,j]>BMA_pre[j,i]){
        BMA_pre[j,i]=0
      }else if (BMA_pre[i,j]<BMA_pre[j,i]){
        BMA_pre[i,j]=0
      }else{
        message(print("node",i,"and node",j,"is undirected"))
      }
    }
  }
} #?Ō???eles?͂Ȃ?

#???Ɗ??S?ɓƗ??????m?[?h???폜
d_node=c()
for (i in 1:ncol(BMA_pre)){
  if (sum(BMA_pre[i,])==0 && sum(BMA_pre[,i])==0){
    d_node=c(d_node,i)
  }
} #d_node is empty 

BMA=BMA_pre
for (i in 1:nrow(BMA)){
  for (j in 1:ncol(BMA)){
    if (BMA[i,j]>1){
      BMA[i,j]=1
    }
  }
}
#DAG?̊m?F
dag=empty.graph(rownames(BMA))
amat(dag)=as.matrix(BMA) #the specified network contains cycles.

for (i in 3:ncol(BMA)){
  print(i)
  dag=empty.graph(row.names(BMA)[1:i])
  amat(dag)=as.matrix(BMA[1:i,1:i])
} #42 stop

Node <- c()
BMA_1 <- BMA
for (i in 3:ncol(BMA)){
  tryCatch({
    dag=empty.graph(rownames(BMA_1)[1:i])
    amat(dag)=as.matrix(BMA_1[1:i,1:i])
  },
  error=function(e){
    node <<- c(node,i)
    BMA_1 <<- BMA_1[-i,-i]
  })
  if(i==ncol(BMA_1)) break
}
node <- sapply(1:length(node), function(x){node[x]+x-1})


dag=empty.graph(rownames(BMA_final))
amat(dag)=as.matrix(BMA_final)
#ok
pdf("../2.Output/05_BayesianNetwork/BMA.pdf",paper="a4",width=9.5,height = 7,pointsize = 10)
g=Rgraphviz::layoutGraph(bnlearn::as.graphNEL(dag))
graph::parRenderInfo(g)=list(graph=list(main="c2PB BMA"),nodes=list(fontsize=30))
Rgraphviz::renderGraph(g)
dev.off()

#???????킩?????̂ŉ??߂ł??郂?W???[?????I??
#module-module????
rm(list=ls())
data=read.table("Superior temporal gyrus/exp_header_rownames.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
ME=as.data.frame(t(data[-c(1:3),]))

library(WGCNA)
module_cor = cor(ME, use = "p",method ="pearson")
#pvalue?Z?o
nSamples=dim(ME)[1]
module_score_Pvalue = corPvalueStudent(module_cor, nSamples)

module_cor_1=module_cor
for (i in 1:dim(module_score_Pvalue)[1]){
  for (j in 1:dim(module_score_Pvalue)[2]){
    if (module_score_Pvalue[i,j]>=0.05){
      module_cor_1[i,j]=0
    }
  }
}

for (i in 1:dim(module_cor_1)[1]){
  module_cor_1[i,i]=0
}

install.packages("sna")
library(sna)
gplot(module_cor_1,usearrows = F,displaylabels = T)




######################################################################################
#BMA
#c2PB
BMA_pre=read.table("../2.Output/05_BayesianNetwork/BN_SSc/result/BMA_pre_re.csv",sep=",",header = F,row.names = NULL,stringsAsFactors = F)
data=read.table("../2.Output/05_BayesianNetwork/ME_SkinScore_forBN.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
rownames(BMA_pre)=rownames(data)
colnames(BMA_pre)=rownames(data)
for (i in 1:nrow(BMA_pre)){
  for (j in 1:ncol(BMA_pre)){
    if (BMA_pre[i,j]>30 && BMA_pre[j,i]>30){
      if (BMA_pre[i,j]>BMA_pre[j,i]){
        BMA_pre[j,i]=0
      }else if (BMA_pre[i,j]<BMA_pre[j,i]){
        BMA_pre[i,j]=0
      }else{
        message(print("node",i,"and node",j,"is undirected"))
      }
    }
  }
} #?Ō???eles?͂Ȃ?

#???Ɗ??S?ɓƗ??????m?[?h???폜
d_node=c()
for (i in 1:ncol(BMA_pre)){
  if (sum(BMA_pre[i,])==0 && sum(BMA_pre[,i])==0){
    d_node=c(d_node,i)
  }
} #d_node is empty 

BMA=BMA_pre
for (i in 1:nrow(BMA)){
  for (j in 1:ncol(BMA)){
    if (BMA[i,j]>1){
      BMA[i,j]=1
    }
  }
}
#DAG?̊m?F
dag=empty.graph(rownames(BMA))
amat(dag)=as.matrix(BMA) #the specified network contains cycles.

for (i in 3:ncol(BMA)){
  print(i)
  dag=empty.graph(row.names(BMA)[1:i])
  amat(dag)=as.matrix(BMA[1:i,1:i])
} #42 stop

node=c()
for (i in 1:ncol(BMA)){
  tryCatch({
    BMA_1=BMA[-i,-i]
    dag=empty.graph(rownames(BMA_1))
    amat(dag)=as.matrix(BMA_1)
    node=c(node,i)
    rm(BMA_1)
  },error=function(e){print(i)})
}
node #22 36 42

in_36=which(BMA[,36]==1) # 12 35 39
out_42=which(BMA[42,]==1) #5  7 10 12 20 28 34 35 41

for (j in in_36){
  BMA_2=BMA[-j,-j]
  node=c()
  for (i in 1:dim(BMA_2)[1]){
    tryCatch({
      BMA_1=BMA_2[-i,-i]
      dag=empty.graph(rownames(BMA_1))
      amat(dag)=as.matrix(BMA_1)
      node=c(node,i)
      rm(BMA_1)
    },error=function(e){print(i)})
  }
  assign(paste0("node_",j),node)
  rm(BMA_2)
  rm(node)
}

BMA_pre[c(22,35,36,42),c(22,35,36,42)]
BMA_pre[c(12,22,36,42),c(12,22,36,42)]

#36->22????
BMA_final=BMA
BMA_final[36,22]=0
dag=empty.graph(rownames(BMA_final))
amat(dag)=as.matrix(BMA_final)
#ok
pdf("../2.Output/05_BayesianNetwork/BMA.pdf",paper="a4",width=9.5,height = 7,pointsize = 10)
g=Rgraphviz::layoutGraph(bnlearn::as.graphNEL(dag))
graph::parRenderInfo(g)=list(graph=list(main="c3PB BMA"),nodes=list(fontsize=30))
Rgraphviz::renderGraph(g)
dev.off()

#???????킩?????̂ŉ??߂ł??郂?W???[?????I??
#module-module????
rm(list=ls())
data=read.table("Superior temporal gyrus/exp_header_rownames.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
ME=as.data.frame(t(data[-c(1:3),]))

library(WGCNA)
module_cor = cor(ME, use = "p",method ="pearson")
#pvalue?Z?o
nSamples=dim(ME)[1]
module_score_Pvalue = corPvalueStudent(module_cor, nSamples)

module_cor_1=module_cor
for (i in 1:dim(module_score_Pvalue)[1]){
  for (j in 1:dim(module_score_Pvalue)[2]){
    if (module_score_Pvalue[i,j]>=0.05){
      module_cor_1[i,j]=0
    }
  }
}

for (i in 1:dim(module_cor_1)[1]){
  module_cor_1[i,i]=0
}

install.packages("sna")
library(sna)
gplot(module_cor_1,usearrows = F,displaylabels = T)







005_BN_SelectSampler.R





####?㑤?????̃T???v?????p????BN####
rm(list=ls())
library(tidyverse)
#setwd("\\\\cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/02_CNS/001_AD/2019XX_AD_patient_GSE84422_reanalysis/analysis/module analysis")
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
if(0){
  ME=read.table("mergedMEs_signed_ver.2.csv",sep=",",header = T,row.names = 1,stringsAsFactors = F)
  ME_1=cbind(data.frame("region"=rownames(ME)),ME)
  ME_1=subset(ME_1,grepl("X05",ME_1$region))
  
  Trait=read.table("Traits.csv",sep=",",header = T,row.names = 1)
  Trait_1=cbind(data.frame("region"=rownames(Trait)),Trait)
  Trait_1=subset(Trait_1,grepl("X05",Trait_1$region))
  
  MergedData=merge(Trait_1,ME_1,by.x = "region",by.y = "region")
  
  rownames(MergedData)=MergedData$region
  MergedData=MergedData[,-1]
  colnames(MergedData)[2]="Cerad"
  colnames(MergedData)[3]="Braak"
  
  MergedData_1=as.data.frame(t(MergedData))
  MergedData_1=MergedData_1[-which(rownames(MergedData_1)=="grey"),]
  
  dir.create("Superior temporal gyrus")
  write.table(MergedData_1,"Superior temporal gyrus/exp_header_rownames.txt",sep="\t",row.names = T,col.names = T)
  write.table(MergedData_1,"Superior temporal gyrus/exp.txt",sep="\t",row.names = F,col.names = F)
}

ME_skinscore_file_path<-"../../tmp/ME_SkinScore_forBN_nafilled.txt"
result_dir_path<-"../../tmp/BN_SSc/result/"
BN_result_dir_path<-"../../tmp/BN_SSc/BN_results/"
MergedData_1<-read_tsv(ME_skinscore_file_path) %>% 
  as.data.frame() %>% 
  mutate(variable=str_replace_all(variable,"ME","")) %>% 
  column_to_rownames("variable")

#LL??Mcmc?ɑ΂??郊?X?g???쐬
library(bnlearn)
library(Rgraphviz)
sampler=c("1PB","2PB","3PB","4PB","c2PB","c3PB","c4PB","M.1PB","M.2PB","M.3PB","M.4PB","M.c2PB","M.c3PB","M.c4PB","M.REV50","M.STR","REV50","STR")
#
tmp<-tibble(
  module=rownames(MergedData_1),
  module_function=c("Myofibroblast","ECM","Protein lipidation","Mitochondria","Immune Response","ECM","Apoptosis",
                    "ER, RNA","Blood Vessel Development","Blood Vessel Development","Keratinization","Metabolism, Cell Cycle",
                    "ER, RNA","Nucleic Acid Metabolism","Histone, Chromatin","","")) 
tmp<-tibble(
  module=rownames(MergedData_1),
  module_function=c("MFB","FB,ECM","MFB,Blood,Protein lipidation","DC,Mitochon","Immune ,B?","FB,ECM","MP, Apoptosis",
                    "DC,ER, RNA","FB,Blood Vessel Development","FB,MP,Blood Vessel Development","Keratinization","Epithelium, Cell Cycle",
                    "DC,ER, RNA","DC, Nucleic Acid Metabolism","Epithelium, Cell Cycle","","")) 
module_function<-tmp %>% 
  mutate(node_label=str_c(module,module_function,sep="\n"))
#node=rownames(MergedData_1)
node=module_function$node_label

LL_list=list()
Mcmc_list=list()
for(i in c(1:length(sampler))){
  #filename?̒??`
  file_in_LL=paste(result_dir_path,"LL_",sampler[i],".csv",sep="")
  file_in_Mcmc=paste(result_dir_path,"Mcmc_",sampler[i],".csv",sep="")
  
  #LL
  LL=read.table(file_in_LL,sep=",",header=F,row.names = NULL,stringsAsFactors = F)
  LL=data.frame("iteration"=c(1:1000),LL)
  LL_list=c(LL_list,list(LL))
  
  #Mcmc
  Mcmc=read.table(file_in_Mcmc,sep=",",header=F,row.names = NULL,stringsAsFactors = F)
  rownames(Mcmc)=c(node)
  colnames(Mcmc)=c(node)
  Mcmc_mat=as.matrix(Mcmc)
  dag=empty.graph(c(node))
  amat(dag)=Mcmc_mat
  Mcmc_list=c(Mcmc_list,list(dag))
}
#rename
names(LL_list)=sampler
names(Mcmc_list)=sampler
#save
save(LL_list,file=paste0(result_dir_path,"LL_list.rda"))
save(Mcmc_list,file=paste0(result_dir_path,"Mcmc_list.rda"))

#LL ?o??
dir.create(BN_result_dir_path)
pdf(paste0(BN_result_dir_path,"LL.pdf"),paper="a4",width=9.5,height = 7,pointsize = 10)
par(mfrow=c(1,2))
for(i in c(1:length(sampler))){
  main=sampler[i]
  plot(LL_list[[i]],xlab="iteration",ylab="Log Likelihood",main=main)
}
dev.off()

#max_LL???o
max_LL=sapply(1:length(sampler),function(x){max(LL_list[[x]]$V1)})
max_LL=data.frame("sampler"=c(sampler),"max_LL"=max_LL)
write.table(max_LL,paste0(BN_result_dir_path,"max_LL.csv"),sep=",",row.names = F,col.names = T)

#DAG ?o??
pdf(paste0(BN_result_dir_path,"DAG_ver3.pdf"),paper="a4",width=9.5,height = 7,pointsize = 10)
par(mfrow=c(1,2))
for(i in c(1:length(sampler))){
  main=sampler[i]
  g=Rgraphviz::layoutGraph(bnlearn::as.graphNEL(Mcmc_list[[i]]))
  graph::parRenderInfo(g)=list(graph=list(main=main),nodes=list(fontsize=30))
  Rgraphviz::renderGraph(g)
}
dev.off()

###########BMA#################
#c3PB
BMA_pre=read.table("Superior temporal gyrus/GSE84422_STG_BN/result/BMA_pre_re.csv",sep=",",header = F,row.names = NULL,stringsAsFactors = F)
data=read.table("Superior temporal gyrus/exp_header_rownames.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
rownames(BMA_pre)=rownames(data)
colnames(BMA_pre)=rownames(data)
for (i in 1:nrow(BMA_pre)){
  for (j in 1:ncol(BMA_pre)){
    if (BMA_pre[i,j]>30 && BMA_pre[j,i]>30){
      if (BMA_pre[i,j]>BMA_pre[j,i]){
        BMA_pre[j,i]=0
      }else if (BMA_pre[i,j]<BMA_pre[j,i]){
        BMA_pre[i,j]=0
      }else{
        message(print("node",i,"and node",j,"is undirected"))
      }
    }
  }
} #?Ō???eles?͂Ȃ?

#???Ɗ??S?ɓƗ??????m?[?h???폜
d_node=c()
for (i in 1:ncol(BMA_pre)){
  if (sum(BMA_pre[i,])==0 && sum(BMA_pre[,i])==0){
    d_node=c(d_node,i)
  }
} #d_node is empty 

BMA=BMA_pre
for (i in 1:nrow(BMA)){
  for (j in 1:ncol(BMA)){
    if (BMA[i,j]>1){
      BMA[i,j]=1
    }
  }
}
#DAG?̊m?F
dag=empty.graph(rownames(BMA))
amat(dag)=as.matrix(BMA) #the specified network contains cycles.

for (i in 3:ncol(BMA)){
  print(i)
  dag=empty.graph(row.names(BMA)[1:i])
  amat(dag)=as.matrix(BMA[1:i,1:i])
} #42 stop

node=c()
for (i in 1:ncol(BMA)){
  tryCatch({
    BMA_1=BMA[-i,-i]
    dag=empty.graph(rownames(BMA_1))
    amat(dag)=as.matrix(BMA_1)
    node=c(node,i)
    rm(BMA_1)
  },error=function(e){print(i)})
}
node #22 36 42

in_36=which(BMA[,36]==1) # 12 35 39
out_42=which(BMA[42,]==1) #5  7 10 12 20 28 34 35 41

for (j in in_36){
  BMA_2=BMA[-j,-j]
  node=c()
  for (i in 1:dim(BMA_2)[1]){
    tryCatch({
      BMA_1=BMA_2[-i,-i]
      dag=empty.graph(rownames(BMA_1))
      amat(dag)=as.matrix(BMA_1)
      node=c(node,i)
      rm(BMA_1)
    },error=function(e){print(i)})
  }
  assign(paste0("node_",j),node)
  rm(BMA_2)
  rm(node)
}

BMA_pre[c(22,35,36,42),c(22,35,36,42)]
BMA_pre[c(12,22,36,42),c(12,22,36,42)]

#36->22????
BMA_final=BMA
BMA_final[36,22]=0
dag=empty.graph(rownames(BMA_final))
amat(dag)=as.matrix(BMA_final)
#ok
pdf("Superior temporal gyrus/BN_results/BMA.pdf",paper="a4",width=9.5,height = 7,pointsize = 10)
g=Rgraphviz::layoutGraph(bnlearn::as.graphNEL(dag))
graph::parRenderInfo(g)=list(graph=list(main="c3PB BMA"),nodes=list(fontsize=30))
Rgraphviz::renderGraph(g)
dev.off()

#???????킩?????̂ŉ??߂ł??郂?W???[?????I??
#module-module????
rm(list=ls())
data=read.table("Superior temporal gyrus/exp_header_rownames.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
ME=as.data.frame(t(data[-c(1:3),]))

library(WGCNA)
module_cor = cor(ME, use = "p",method ="pearson")
#pvalue?Z?o
nSamples=dim(ME)[1]
module_score_Pvalue = corPvalueStudent(module_cor, nSamples)

module_cor_1=module_cor
for (i in 1:dim(module_score_Pvalue)[1]){
  for (j in 1:dim(module_score_Pvalue)[2]){
    if (module_score_Pvalue[i,j]>=0.05){
      module_cor_1[i,j]=0
    }
  }
}

for (i in 1:dim(module_cor_1)[1]){
  module_cor_1[i,i]=0
}

install.packages("sna")
library(sna)
gplot(module_cor_1,usearrows = F,displaylabels = T)






005_BN_SelectSampler_2.R





library(tidyverse)
library(tidygraph)
library(ggraph)
library(igraph)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
ME_skinscore_file_path<-"../../tmp/ME_SkinScore_forBN_nafilled.txt"
#result_dir_path<-"../../tmp/BN_SSc/result/"
#BN_result_dir_path<-"../../tmp/BN_SSc/BN_results/"
BN_result_files<-list.files("../../tmp/BN_SSc/result/",pattern ="Mcmc_*" ,full.names = T)
BN_result_files<-BN_result_files[-which(BN_result_files %in% "../../tmp/BN_SSc/result/Mcmc_list.rda" )]

MergedData_1<-read_tsv(ME_skinscore_file_path) %>% 
  as.data.frame() %>% 
  mutate(variable=str_replace_all(variable,"ME","")) %>% 
  column_to_rownames("variable")
if(0){
tmp<-tibble(
  module=rownames(MergedData_1),
  module_function=c("MFB,","FB,ECM","MFB,Protein lipidation,Blood","DC,Mitochon","Immune ,B?","FB,ECM","MP,Apoptosis",
                    "DC,ER, RNA","FB,Blood Vessel Development","FB,MP,Blood Vessel Development","Keratinization","Epithelium, Cell Cycle",
                    "DC,ER, RNA","DC, Nucleic Acid Metabolism","Epithelium, Cell Cycle","","")) 
}

module2function<-tibble(
  module=rownames(MergedData_1),
  celltype=c("MFB","FB", "MFB" ,"DC","Immune","FB" ,"MP","DC","FB","FB_MP","Keratinocyte","Epithelium","DC","DC","Epithelium","",""),
  module_function=c("","ECM","Protein lipidation","Mitochondria","Immune response","ECM","Apoptosis","ER_RNA","Blood Vessel Development","Blood Vessel Development",
                    "Keratinization","Cell Cycle","ER_RNA","Nucleic Acid Metabolism","Cell Cycle","",""))

immune_cells<-c("DC","Immune","MP")
fibroblast<-c("MFB","FB","FB_MP")
traits<-c("Skin_score_at_biopsy_site","total_skin_score")

BN_result_file<-BN_result_files[2]
for (BN_result_file in BN_result_files){
  print(BN_result_file)
  g<-read_csv(BN_result_file,col_names = rownames(MergedData_1))%>% 
    mutate(from=rownames(MergedData_1)) %>% 
    gather(key=to,value=edge,-from) %>% 
    filter(edge==1) %>% 
    select(-edge) %>% 
    as_tbl_graph() %N>% 
    left_join(module2function,by=c("name"="module")) %>% 
    mutate(module_color=if_else(name %in% traits,"red",name),    #name %>% str_replace("Skin_score_at_biopsy_site","red") %>% str_replace("total_skin_score" ,"red"),
           celltype_large=if_else(celltype %in% immune_cells,"immune cell",if_else(celltype %in% fibroblast,"fibroblast","")),
           celltype_color=if_else(celltype_large=="immune cell","#f79646",if_else(celltype_large=="fibroblast","#9bbb59","grey")) )%>% 
    ggraph(layout = "sugiyama")+
    geom_edge_diagonal2(arrow = arrow(length = unit(4, 'mm'),type = "closed"), 
                   end_cap = circle(10, 'mm'))+
    geom_node_point(aes(color=celltype_color),size=20)+
    geom_node_point(aes(color=module_color),size=10)+
    geom_node_text(aes(label=str_c(name,"\n",module_function)))+
    scale_color_identity()+
    theme_graph()+
    theme(plot.margin= unit(c(0,0,0,0), "lines"))
  savefile<-paste0("../2.Output/05_BayesianNetwork/Image/",str_sub(basename(BN_result_file),6,-5),".png")
  ggsave(plot = g,filename = savefile,dpi=300)
}
#"#0086ab","#f79646","#9bbb59","#da6272","#777777","#bfbfbf"

BN_result_file<-"../../tmp/BN_SSc/result/Mcmc_c2PB.csv"

graph<-read_csv(BN_result_file,col_names = rownames(MergedData_1))%>% 
  mutate(from=rownames(MergedData_1)) %>% 
  gather(key=to,value=edge,-from) %>% 
  filter(edge==1) %>% 
  select(-edge) %>% 
  as_tbl_graph() %N>% 
  left_join(module2function,by=c("name"="module")) %>% 
  mutate(module_color=if_else(name %in% traits,"red",name),    #name %>% str_replace("Skin_score_at_biopsy_site","red") %>% str_replace("total_skin_score" ,"red"),
         celltype_large=if_else(celltype %in% immune_cells,"immune cell",if_else(celltype %in% fibroblast,"fibroblast","")),
         celltype_color=if_else(celltype_large=="immune cell","#f79646",if_else(celltype_large=="fibroblast","#9bbb59","grey")) )

layout<-create_layout(graph,layout = "sugiyama")
#####################
layout$x<--layout$x
layout$y<--layout$y

#tmp<-layout[layout$name=="white",c("x","y")]
#layout[layout$name=="white",c("x","y")]<-layout[layout$name=="yellow",c("x","y")]
#layout[layout$name=="yellow",c("x","y")]<-tmp

layout[layout$name=="yellow",c("x","y")]<-c(0.5,2.3)
layout[layout$name=="tan",c("x","y")]<-c(-0.75,2)
layout[layout$name=="white",c("x","y")]<-c(1,2)
layout[layout$name=="turquoise",c("x","y")]<-c(0,1.7)

layout[layout$name=="darkmagenta",c("x","y")]<-c(-1,1.1)
layout[layout$name=="greenyellow",c("x","y")]<-c(0,1.1)

layout[layout$name=="royalblue",c("x","y")]<-c(-0.5,0.8)
layout[layout$name=="orange",c("x","y")]<-c(0.5,0.8)

layout[layout$name=="lightyellow",c("x","y")]<-c(0,0.5)

layout[layout$name=="midnightblue",c("x","y")]<-c(0.5,0.2)
layout[layout$name=="violet",c("x","y")]<-c(1.25,0.2)

layout[layout$name=="darkred",c("x","y")]<-c(0.5,-0.4)

layout[layout$name=="plum1",c("x","y")]<-c(-0,-0.7)

layout[layout$name=="darkgrey",c("x","y")]<-c(-1,-1)
layout[layout$name=="lightcyan",c("x","y")]<-c(1,-1)

layout[layout$name=="total_skin_score",c("x","y")]<-c(-0.5,-1.3)
layout[layout$name=="Skin_score_at_biopsy_site",c("x","y")]<-c(0.5,-1.6)
####################################################################################

ggraph(layout)+
  geom_edge_diagonal(arrow = arrow(length = unit(4, 'mm'),type = "closed"), 
                      end_cap = circle(10, 'mm'),color="#0086ab",
                      width=0.9)+
  #geom_node_point(aes(color=celltype_color),size=20)+
  geom_node_point(aes(color=module_color%>% str_replace_all("white","antiquewhite")),size=20)+
  geom_node_text(aes(label=name %>% str_replace_all("_"," ")),color="#404040")+
  scale_color_identity()+
  theme_graph()+
  theme(plot.margin= unit(c(0,0,0,0), "lines"),
        panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))#+
  xlim(c(-1.25,1.5))
savefile<-paste0("../2.Output/05_BayesianNetwork/Image/",str_sub(basename(BN_result_file),6,-5),"_modified4",".png")
ggsave(filename = savefile,width = 190, height = 190, units = "mm",bg="transparent")
####################################################################################

graph<-read_csv(BN_result_file,col_names = rownames(MergedData_1))%>% 
  mutate(from=rownames(MergedData_1)) %>% 
  gather(key=to,value=edge,-from) %>% 
  filter(edge==1) %>% 
  select(-edge) %>% 
  as_tbl_graph() %N>% 
  left_join(module2function,by=c("name"="module")) %>% 
  mutate(module_color=if_else(name %in% traits,"red",name),    #name %>% str_replace("Skin_score_at_biopsy_site","red") %>% str_replace("total_skin_score" ,"red"),
         celltype_large=if_else(celltype %in% immune_cells,"immune cell",if_else(celltype %in% fibroblast,"fibroblast","")),
         celltype_color=if_else(celltype_large=="immune cell","#f79646",if_else(celltype_large=="fibroblast","#9bbb59","grey")) )
ggraph(graph ,layout = 'sugiyama')+
  geom_edge_diagonal(arrow = arrow(length = unit(4, 'mm'),type = "closed"), 
                     end_cap = circle(10, 'mm'),color="#0086ab",
                     width=0.9)+
  #geom_node_point(aes(color=celltype_color),size=20)+
  geom_node_point(aes(color=module_color%>% str_replace_all("white","antiquewhite")),size=20)+
  geom_node_text(aes(label=name %>% str_replace_all("_"," ")),color="#404040")+
  scale_color_identity()+
  theme_graph()+
  theme(plot.margin= unit(c(0,0,0,0), "lines"),
        panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))#+
  xlim(c(-1.25,1.5))
  ####################################################################################







005_Correlation_for_BaysianNetwork.R





library(WGCNA)
library(caret)
library(GGally)
library(ggraph)
library(tidygraph)
library(tidyverse)

correlation_method<-"spearman"
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
lnames<-load("../2.Output/02_signed_2_unmerged/signed_2/networkConstruction_StepByStep_signed_2.Rdata")
sample_annotation_file_path<-"../1.Data/GSE58095/GSE58095_all_sample_annotations.txt"
#MEs_signed
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor
sample_annotation <- read_tsv(sample_annotation_file_path,skip = 4,n_max = 102) %>% 
  rename_all(function(x) x %>% str_replace_all(" ","_")%>% 
               str_replace_all(":","_") %>% 
               str_replace_all("\\(","") %>% 
               str_replace_all("\\)","") %>% 
               str_replace_all("=","is"))%>% 
  dplyr::select(-c(title,Characteristics__Patient,Characteristics__GSE47162_Sample_name,Characteristics__Time_point,
            molecule,label,description,platform)) %>% 
  mutate(Characteristics__Group=Characteristics__Group %>% replace_na("SSc"))
  

tmp <- dummyVars(~., data=sample_annotation %>% dplyr::select(-"Sample_name"))
sample_annotation.dummy <- as.data.frame(predict(tmp, sample_annotation)) %>% 
  mutate(Sample_name=sample_annotation$Sample_name) %>%
  as_tibble()%>% 
  rename_at(vars(starts_with("Characteristics__")),function(x) str_sub(x,18,-1)) %>% 
  dplyr::select(Sample_name,everything())


MEs_withTrait<-MEs_signed %>% as.data.frame() %>% 
  rownames_to_column("Sample_name") %>% 
  as_tibble() %>% 
  mutate(Sample_name=str_c("SAMPLE",Sample_name)) %>% 
  left_join(sample_annotation.dummy,by="Sample_name")

ME_Trait_cor<-MEs_withTrait%>% 
  dplyr::select(-c("Sample_name","GroupSSc","GroupControl") )%>% stats::cor( use = "complete.obs",method = correlation_method) %>% 
  as.data.frame() %>% 
  rownames_to_column("vs1") %>% 
  as_tibble() %>% 
  gather(key="vs2",value="correlation",-vs1) #

tmp<-MEs_withTrait%>% 
  dplyr::select(-c("Sample_name","GroupSSc","GroupControl") )%>% stats::cor( use = "complete.obs",method = correlation_method) 
tmp[upper.tri(tmp)] <-NA
tmp%>% 
  as.data.frame() %>% 
  rownames_to_column("vs1") %>% 
  as_tibble() %>% 
  gather(key="vs2",value="correlation",-vs1) %>% 
  filter(vs1!=vs2) %>% 
  filter(!near(correlation,-1.0)) %>% 
  arrange(-abs(correlation)) %>% 
  #write_csv("../2.Output/05_BayesianNetwork/ME_Trait_correlation.csv")
  ggplot(aes())+
  geom_histogram(aes(x=correlation))
  

ME_Trait_cor %>% 
  ggplot(aes(x=vs1,y=vs2,fill=correlation))+
  geom_tile()+
  scale_fill_gradient2() + 
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  xlab("")+
  ylab("")

#MEs_withTrait %>% select(-"Sample_name") %>% ggpairs(aes_string(alpha=0.5))
target<-c(str_c("ME",significant_module),c("total_skin_score","Skin_score_at_biopsy_site"))
ME_Trait_graph<-ME_Trait_cor %>% 
  filter(vs1 %in% target & vs2 %in% target) %>% 
  filter(vs1!=vs2) %>% 
  rename(vs1="from",vs2="to") %>% 
  as_tbl_graph()

ME_Trait_graph %>% 
  activate("edges") %>% 
  filter(abs(correlation)>0.6) %>% 
  activate("nodes") %>%
  mutate(color=recode(name, `total_skin_score` = "TRred",`Skin_score_at_biopsy_site` = "TRred", .default = name),
         color=str_sub(color,3,-1)) %>% 
  mutate(label=recode(name, `total_skin_score` = "TRtotal_skin_score",`Skin_score_at_biopsy_site` = "TRSkin_score_at_biopsy_site", .default = name),
         label=str_sub(label,3,-1)) %>% 
  ggraph(layout = "kk")+
  geom_edge_link0(aes(width=abs(correlation),color=correlation>0),alpha=0.8)+
  geom_node_point(aes(color="#404040"),size=21)+
  geom_node_point(aes(color=color),size=20)+
  geom_node_text(aes(label=label),color="#404040")+
  scale_color_identity()+
  theme_graph()+
  theme(legend.position ="none")+
  labs(title="total skinscore correlation P<0.01 | abs(correlation)>0.6")
ggsave(filename = "../2.Output/05_BayesianNetwork/correlation_network.png")

#baysian network ?`???p?e?L?X?g?t?@?C???쐬
#?s???ϐ??A?????T???v????tab???؂??t?@?C???? cserv48?̂ǂ????ɒu??
MEs_withTrait %>% 
  dplyr::select(Sample_name,tidyselect::one_of(target)) %>% 
  replace_na(list(total_skin_score=0,Skin_score_at_biopsy_site=0)) %>% 
  gather(key=variable,value=value,-Sample_name) %>% 
  spread(key=Sample_name,value=value) %>% 
  write_tsv("../2.Output/05_BayesianNetwork/ME_SkinScore_forBN_nafilled.txt")
  

###########################################################################
#?t?я????̉???
sample_annotation.dummy%>% select(-c("Sample_name","GroupSSc","GroupControl") )%>% 
  stats::cor( use = "complete.obs",method = correlation_method) %>% 
  as.data.frame() %>% 
  ggpairs(aes_string(alpha=0.5))

tmp<-sample_annotation.dummy%>% select(-c("Sample_name","GroupSSc","GroupControl") )%>% 
  stats::cor( use = "complete.obs",method = correlation_method) %>% 
  as.data.frame()
tmp[upper.tri(tmp)] <-NA

tmp%>% 
  as.data.frame() %>% 
  rownames_to_column("vs1") %>% 
  as_tibble() %>% 
  gather(key="vs2",value="correlation",-vs1) %>% 
  filter(vs1!=vs2) %>% 
  filter(!near(correlation,-1.0)) %>% 
  arrange(-abs(correlation)) %>% 
  write_csv("../2.Output/05_BayesianNetwork/Trait_correlation.csv")
  ggplot(aes())+
  geom_histogram(aes(x=correlation))






005_DrawBaysianNetworkPlot.R





library(tidyverse)
library(tidygraph)
library(ggraph)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

BN_result_file<-"../../tmp/BN_SSc/result/Mcmc_c2PB.csv"

graph<-read_csv(BN_result_file,col_names = rownames(MergedData_1))%>% 
  mutate(from=rownames(MergedData_1)) %>% 
  gather(key=to,value=edge,-from) %>% 
  filter(edge==1) %>% 
  select(-edge) %>% 
  as_tbl_graph() %N>% 
  left_join(module2function,by=c("name"="module")) %>% 
  mutate(module_color=if_else(name %in% traits,"red",name),    #name %>% str_replace("Skin_score_at_biopsy_site","red") %>% str_replace("total_skin_score" ,"red"),
         celltype_large=if_else(celltype %in% immune_cells,"immune cell",if_else(celltype %in% fibroblast,"fibroblast","")),
         celltype_color=if_else(celltype_large=="immune cell","#f79646",if_else(celltype_large=="fibroblast","#9bbb59","grey")) )

layout<-create_layout(graph,layout = "kk")
#####################
layout$x<--layout$x
layout$y<--layout$y

#tmp<-layout[layout$name=="white",c("x","y")]
#layout[layout$name=="white",c("x","y")]<-layout[layout$name=="yellow",c("x","y")]
#layout[layout$name=="yellow",c("x","y")]<-tmp

layout[layout$name=="yellow",c("x","y")]<-c(0.5,2.3)
layout[layout$name=="tan",c("x","y")]<-c(-0.75,2)
layout[layout$name=="white",c("x","y")]<-c(1,2)
layout[layout$name=="turquoise",c("x","y")]<-c(0,1.7)

layout[layout$name=="darkmagenta",c("x","y")]<-c(-1,1.1)
layout[layout$name=="greenyellow",c("x","y")]<-c(0,1.1)

layout[layout$name=="royalblue",c("x","y")]<-c(-0.5,0.8)
layout[layout$name=="orange",c("x","y")]<-c(0.5,0.8)

layout[layout$name=="lightyellow",c("x","y")]<-c(0,0.5)

layout[layout$name=="midnightblue",c("x","y")]<-c(0.5,0.2)
layout[layout$name=="violet",c("x","y")]<-c(1.25,0.2)

layout[layout$name=="darkred",c("x","y")]<-c(0.5,-0.4)

layout[layout$name=="plum1",c("x","y")]<-c(-0,-0.7)

layout[layout$name=="darkgrey",c("x","y")]<-c(-1,-1)
layout[layout$name=="lightcyan",c("x","y")]<-c(1,-1)

layout[layout$name=="total_skin_score",c("x","y")]<-c(-0.5,-1.3)
layout[layout$name=="Skin_score_at_biopsy_site",c("x","y")]<-c(0.5,-1.6)
#####################

ggraph(layout)+
  geom_edge_diagonal(arrow = arrow(length = unit(4, 'mm'),type = "closed"), 
                     end_cap = circle(10, 'mm'),color="#0086ab",
                     width=0.9)+
  #geom_node_point(aes(color=celltype_color),size=20)+
  geom_node_point(aes(color=module_color%>% str_replace_all("white","antiquewhite")),size=20)+
  geom_node_text(aes(label=name %>% str_replace_all("_"," ")),color="#404040")+
  scale_color_identity()+
  theme_graph()+
  theme(plot.margin= unit(c(0,0,0,0), "lines"),
        panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))+
  xlim(c(-1.25,1.5))
savefile<-paste0("../2.Output/05_BayesianNetwork/Image/",str_sub(basename(BN_result_file),6,-5),"_modified4",".png")
ggsave(filename = savefile,width = 190, height = 190, units = "mm",bg="transparent")






005_DrawBNPlot_afterBMA.R





library(tidyverse)
library(tidygraph)
library(ggraph)
library(igraph)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

ME_skinscore_file_path<-"../2.Output/05_BayesianNetwork/ME_SkinScore_forBN_nafilled.txt"
BMA_pre_file_path<-"../2.Output/05_BayesianNetwork/BN_SSc/result/BMA_pre_re.csv"
traits<-c("Skin_score_at_biopsy_site","total_skin_score")

ME_skinscore_df<-read_tsv(ME_skinscore_file_path) %>% 
  as.data.frame() %>% 
  mutate(variable=str_replace_all(variable,"ME","")) %>% 
  column_to_rownames("variable")
BMA_pre=read.table(BMA_pre_file_path,sep=",",header = F,row.names = NULL,stringsAsFactors = F)

my_scale<-function(variable,sc_min,sc_max){
  tmp<-(variable-min(variable))*((sc_max-sc_min)/(max(variable)-min(variable)))
  return(tmp+sc_min)
}

#BMA_graph %N>% as_tibble() %>% .$centrality %>% my_scale(5,30)
#############################################################################
#PageRank ?l?m?F?p
BMA.mat<-read_csv(BMA_pre_file_path,col_names = rownames(ME_skinscore_df))%>% 
  mutate(from=rownames(ME_skinscore_df)) %>% 
  as.data.frame() %>% 
  column_to_rownames("from") %>% 
  as.matrix()

BMA.mat.t<-BMA.mat %>% t()
diag(BMA.mat.t)[which(rowSums(BMA.mat.t)==0)] <- 100 #rowsum??0?ł????ꍇ?ɑΊp??????1

g <- graph_from_adjacency_matrix(BMA.mat.t/100, mode="directed",weighted = T)

page_rank(g)$vector

#BMA.mat.t<-BMA.mat %>% t()
diag(BMA.mat)[which(rowSums(BMA.mat)==0)] <- 100 #rowsum??0?ł????ꍇ?ɑΊp??????1

g <- graph_from_adjacency_matrix(BMA.mat/100, mode="directed",weighted = T)

page_rank(g)$vector

((BMA.mat/100)/rowSums(BMA.mat/100)) %>% t()->m
mat<-BMA.mat/100

m<-t(mat/rowSums(mat))
n<-nrow(m)
c=0.85
M<-(c*m)+((1-c)*matrix(1/n,nrow=n,ncol=n))
ev=eigen(M)$vectors[,1]
pr<-ev/sum(ev)
names(pr)<-rownames(mat)

#############################################################################

#centrality_pagerank?o?[?W????
BMA_graph.withPageRank<-read_csv(BMA_pre_file_path,col_names = rownames(ME_skinscore_df))%>% 
  mutate(from=rownames(ME_skinscore_df)) %>% 
  gather(key=to,value=edge,-from) %>% 
  filter(edge>=30) %>% 
  as_tbl_graph() %N>%
  
  mutate(centrality = centrality_pagerank(weights = edge,directed = T), #centrality_eigen
         size=my_scale(centrality,5,30)) %>%
  mutate(color=if_else(name %in% traits,"red",name))

#?ŗL?x?N?g?????S???o?[?W????
BMA_graph<-read_csv(BMA_pre_file_path,col_names = rownames(ME_skinscore_df))%>% 
  mutate(from=rownames(ME_skinscore_df)) %>% 
  gather(key=to,value=edge,-from) %>% 
  filter(edge>=30) %>% 
  as_tbl_graph() %N>%
  #centrality_pagerank
  mutate(centrality = centrality_eigen(weights = edge,directed = F), #centrality_eigen
         size=my_scale(centrality,5,30)) %>%
  mutate(color=if_else(name %in% traits,"red",name))

#############################################################################

ggraph(BMA_graph,layout = "sugiyama")+
  geom_edge_diagonal(aes(edge_width=edge,edge_alpha=edge,edge_colour=edge),
                     arrow = arrow(length = unit(4, 'mm'),type = "closed"),end_cap = circle(10, 'mm')#,color="#0086ab"
                     )+
  #geom_node_point(aes(color=celltype_color),size=20)+
  geom_node_point(aes(color=color%>% str_replace_all("white","antiquewhite")),size=20)+
  geom_node_text(aes(label=name %>% str_replace_all("_"," ")),color="#404040")+
  scale_color_identity()+
  scale_edge_width_continuous(range=c(0.5,2))+
  scale_edge_color_gradient(low = "red",high = "#0086ab")+
  theme_graph()+
  theme(plot.margin= unit(c(0,0,0,0), "lines"),
        panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))#+
  xlim(c(-1.25,1.5))
savefile<-paste0("../2.Output/05_BayesianNetwork/Image/","BNplot_c2PB_afterBMA",".png")
ggsave(filename = savefile,width = 190, height = 190, units = "mm",bg="transparent")


BMA_graph %E>% as_tibble() %>% 
  ggplot(aes(x=edge))+
  geom_histogram()



#####################
layout<-create_layout(BMA_graph,layout = "sugiyama")

layout[layout$name=="yellow",c("x","y")]<-c(0.5,2.3)
layout[layout$name=="tan",c("x","y")]<-c(-0.75,2)
layout[layout$name=="white",c("x","y")]<-c(1,2)
layout[layout$name=="turquoise",c("x","y")]<-c(0,1.7)

layout[layout$name=="darkmagenta",c("x","y")]<-c(-1,1.1)
layout[layout$name=="greenyellow",c("x","y")]<-c(0,1.1)

layout[layout$name=="royalblue",c("x","y")]<-c(-0.5,0.8)
layout[layout$name=="orange",c("x","y")]<-c(0.5,0.8)

layout[layout$name=="lightyellow",c("x","y")]<-c(0,0.5)

layout[layout$name=="midnightblue",c("x","y")]<-c(0.5,0.2)
layout[layout$name=="violet",c("x","y")]<-c(1.25,0.2)

layout[layout$name=="darkred",c("x","y")]<-c(0.5,-0.4)

layout[layout$name=="plum1",c("x","y")]<-c(-0,-0.7)

layout[layout$name=="darkgrey",c("x","y")]<-c(-1,-1)
layout[layout$name=="lightcyan",c("x","y")]<-c(1,-1)

layout[layout$name=="total_skin_score",c("x","y")]<-c(-0.5,-1.3)
layout[layout$name=="Skin_score_at_biopsy_site",c("x","y")]<-c(0.5,-1.6)
####################################################################################

ggraph(layout)+
  geom_edge_diagonal(aes(edge_width=edge,
                         edge_alpha=edge,
                         edge_colour=edge,
                         start_cap = circle(20+1, unit = "dida"), 
                         end_cap = circle(20+1, unit = "dida")),
                     arrow = arrow(length = unit(4, "mm"),type = "closed")
  )+
  #geom_node_point(aes(color=celltype_color),size=20)+
  geom_node_point(aes(color=color%>% str_replace_all("white","antiquewhite")),size=20)+
  geom_node_text(aes(label=name %>% str_replace_all("_"," ")),color="#404040")+
  scale_color_identity()+
  scale_edge_width_continuous(range=c(0.1,2),guide=FALSE)+
  scale_edge_color_gradient(low = "white",high = "#0086ab",guide=FALSE)+
  scale_edge_alpha(range = c(0.1, 1),guide=F)+
  #scale_size(range = c(5,30))+
  scale_size_identity()+
  #guides(alpha=FALSE)+
  theme_graph()+
  theme(plot.margin= unit(c(0,0,0,0), "lines"),
        panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))+
  xlim(c(-1.25,1.5))
savefile<-paste0("../2.Output/05_BayesianNetwork/BMA/","BNplot_c2PB_afterBMA",".png")
ggsave(filename = savefile,width = 190, height = 190, units = "mm",bg="transparent")

####################################################################################
#eigen vector centrality added
ggraph(layout)+
  geom_edge_diagonal(aes(edge_width=edge,
                         edge_alpha=edge,
                         edge_colour=edge,
                         start_cap = circle(node1.size+1, unit = "dida"), 
                         end_cap = circle(node2.size+1, unit = "dida")),
                     arrow = arrow(length = unit(4, "mm"),type = "closed")
  )+
  #geom_node_point(aes(color=celltype_color),size=20)+
  geom_node_point(aes(size=size,color=color%>% str_replace_all("white","antiquewhite")))+
  geom_node_text(aes(label=name %>% str_replace_all("_"," "),size=4.5),color="#404040")+
  scale_color_identity()+
  scale_edge_width_continuous(range=c(0.1,2),guide=FALSE)+
  scale_edge_color_gradient(low = "white",high = "#0086ab",guide=FALSE)+
  scale_edge_alpha(range = c(0.1, 1),guide=F)+
  #scale_size(range = c(5,30))+
  scale_size_identity()+
  #guides(alpha=FALSE)+
  theme_graph()+
  theme(plot.margin= unit(c(0,0,0,0), "lines"),
        panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))+
  xlim(c(-1.25,1.5))+
  ylim(c(-1.5,2.5))

savefile<-paste0("../2.Output/05_BayesianNetwork/BMA/","BNplot_c2PB_afterBMA_withEigenCent","_ver3",".png")
ggsave(filename = savefile,width = 120, height = 120, units = "mm",bg="transparent")

ggraph(layout)+
  geom_edge_diagonal(aes(#edge_width=edge,edge_alpha=edge,edge_colour=edge,
                         #start_cap = circle(node1.size, unit = "mm"), 
                         end_cap = circle(node2.size, unit = "points")),
                     arrow = arrow(length = unit(4, "mm"),type = "closed")
  )



#####################
layout.withPageRank<-create_layout(BMA_graph.withPageRank,layout = "sugiyama")

layout.withPageRank[layout.withPageRank$name=="yellow",c("x","y")]<-c(0.5,2.3)
layout.withPageRank[layout.withPageRank$name=="tan",c("x","y")]<-c(-0.75,2)
layout.withPageRank[layout.withPageRank$name=="white",c("x","y")]<-c(1,2)
layout.withPageRank[layout.withPageRank$name=="turquoise",c("x","y")]<-c(0,1.7)

layout.withPageRank[layout.withPageRank$name=="darkmagenta",c("x","y")]<-c(-1,1.1)
layout.withPageRank[layout.withPageRank$name=="greenyellow",c("x","y")]<-c(0,1.1)

layout.withPageRank[layout.withPageRank$name=="royalblue",c("x","y")]<-c(-0.5,0.8)
layout.withPageRank[layout.withPageRank$name=="orange",c("x","y")]<-c(0.5,0.8)

layout.withPageRank[layout.withPageRank$name=="lightyellow",c("x","y")]<-c(0,0.5)

layout.withPageRank[layout.withPageRank$name=="midnightblue",c("x","y")]<-c(0.5,0.2)
layout.withPageRank[layout.withPageRank$name=="violet",c("x","y")]<-c(1.25,0.2)

layout.withPageRank[layout.withPageRank$name=="darkred",c("x","y")]<-c(0.5,-0.4)

layout.withPageRank[layout.withPageRank$name=="plum1",c("x","y")]<-c(-0,-0.7)

layout.withPageRank[layout.withPageRank$name=="darkgrey",c("x","y")]<-c(-1,-1)
layout.withPageRank[layout.withPageRank$name=="lightcyan",c("x","y")]<-c(1,-1)

layout.withPageRank[layout.withPageRank$name=="total_skin_score",c("x","y")]<-c(-0.5,-1.3)
layout.withPageRank[layout.withPageRank$name=="Skin_score_at_biopsy_site",c("x","y")]<-c(0.5,-1.6)
####################################################################################
ggraph(layout.withPageRank)+
  geom_edge_diagonal(aes(edge_width=edge,
                         edge_alpha=edge,
                         edge_colour=edge,
                         start_cap = circle(node1.size+1, unit = "dida"), 
                         end_cap = circle(node2.size+1, unit = "dida")),
                     arrow = arrow(length = unit(4, "mm"),type = "closed")
  )+
  #geom_node_point(aes(color=celltype_color),size=20)+
  geom_node_point(aes(size=size,color=color%>% str_replace_all("white","antiquewhite")))+
  geom_node_text(aes(label=name %>% str_replace_all("_"," ")),color="#404040")+
  scale_color_identity()+
  scale_edge_width_continuous(range=c(0.1,2),guide=FALSE)+
  scale_edge_color_gradient(low = "white",high = "#0086ab",guide=FALSE)+
  scale_edge_alpha(range = c(0.1, 1),guide=F)+
  #scale_size(range = c(5,30))+
  scale_size_identity()+
  #guides(alpha=FALSE)+
  theme_graph()+
  theme(plot.margin= unit(c(0,0,0,0), "lines"),
        panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))+
  xlim(c(-1.25,1.5))

savefile<-paste0("../2.Output/05_BayesianNetwork/BMA/","BNplot_c2PB_afterBMA_withPageRank",".png")
ggsave(filename = savefile,width = 190, height = 190, units = "mm",bg="transparent")
####################################################################################
#eigen vector centrality added final version

ggraph(layout)+
  geom_edge_diagonal(aes(edge_width=edge,
                         edge_alpha=edge,
                         edge_colour=edge,
                         start_cap = circle(node1.size+1, unit = "dida"), 
                         end_cap = circle(node2.size+1, unit = "dida")),
                     arrow = arrow(length = unit(3, "mm"),type = "closed")
  )+
  #geom_node_point(aes(color=celltype_color),size=20)+
  geom_node_point(aes(size=size,color=color%>% str_replace_all("white","antiquewhite")))+
  geom_node_text(aes(label=name %>% str_replace_all("_"," "),size=4.5),color="#404040")+
  scale_color_identity()+
  scale_edge_width_continuous(range=c(0.05,1),guide=FALSE)+
  scale_edge_color_gradient(low = "white",high = "#0086ab",guide=FALSE)+
  scale_edge_alpha(range = c(0.01, 1),guide=F)+
  #scale_size(range = c(5,30))+
  scale_size_identity()+
  #guides(alpha=FALSE)+
  theme_graph()+
  theme(plot.margin= unit(c(0,0,0,0), "lines"),
        panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))+
  xlim(c(-1.5,1.5))+
  ylim(c(-1.5,2.7))

savefile<-paste0("../2.Output/05_BayesianNetwork/BMA/","BNplot_c2PB_afterBMA_withEigenCent","_ver4",".png")
ggsave(filename = savefile,width = 130, height = 130, units = "mm",bg="transparent")


#BMA final version

ggraph(layout)+
  geom_edge_diagonal(aes(edge_width=edge,
                         edge_alpha=edge,
                         edge_colour=edge,
                         start_cap = circle(16, unit = "dida"), 
                         end_cap = circle(16, unit = "dida")),
                     arrow = arrow(length = unit(3, "mm"),type = "closed")
  )+
  #geom_node_point(aes(color=celltype_color),size=20)+
  geom_node_point(aes(color=color%>% str_replace_all("white","antiquewhite")),size=15)+
  geom_node_text(aes(label=name %>% str_replace_all("_"," "),size=4.5),color="#404040")+
  scale_color_identity()+
  scale_edge_width_continuous(range=c(0.01,1),guide=FALSE)+
  scale_edge_color_gradient(low = "white",high = "#0086ab",guide=FALSE)+
  scale_edge_alpha(range = c(0.01, 1),guide=F)+
  #scale_size(range = c(5,30))+
  scale_size_identity()+
  #guides(alpha=FALSE)+
  theme_graph()+
  theme(plot.margin= unit(c(0,0,0,0), "lines"),
        panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))+
  xlim(c(-1.5,1.5))+
  ylim(c(-1.5,2.7))

savefile<-paste0("../2.Output/05_BayesianNetwork/BMA/","BNplot_c2PB_afterBMA","_ver2",".png")
ggsave(filename = savefile,width = 130, height = 130, units = "mm",bg="transparent")






005_ggdag_test.R





install.packages("ggdag")
library(tidyverse)
library(ggdag)
smoking_ca_dag <- dagify(cardiacarrest ~ cholesterol,
                         cholesterol ~ smoking + weight,
                         smoking ~ unhealthy,
                         weight ~ unhealthy,
                         labels = c("cardiacarrest" = "Cardiac\n Arrest", 
                                    "smoking" = "Smoking",
                                    "cholesterol" = "Cholesterol",
                                    "unhealthy" = "Unhealthy\n Lifestyle",
                                    "weight" = "Weight"),
                         latent = "unhealthy",
                         exposure = "smoking",
                         outcome = "cardiacarrest") %>% 
  tidy_dagitty()

smoking_ca_dag

ggdag_paths(smoking_ca_dag, text = FALSE, use_labels = "label", shadow = TRUE) +
  theme_dag(base_size = 14) +
  theme(legend.position = "none", strip.text = element_blank()) + 
  # set node aesthetics
  scale_color_manual(values = "#0072B2", na.value = "grey80") + 
  # set label aesthetics
  scale_fill_manual(values = "#0072B2", na.value = "grey80") + 
  # set arrow aesthetics
  ggraph::scale_edge_color_manual(values = "#0072B2", na.value = "grey80") +
  scale_y_reverse()+
  ggtitle("Open paths from smoking to cardiac arrest")

ggdag_paths(graph)
#ox-LDL is not gene?
yoshida_target<-c(7124,3596,7040,652,2309, 51341,  960 )
yoshida_target<-c("TNF??", "IL-13", "TFG-??1", "BMP4", "ox-LDL", "FOXO3a", "ZBTB7A", "CD44")

read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  filter(Entrez_Gene_ID %in% yoshida_target ) 


read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  mutate(Synonyms=Synonyms %>% str_split("; ")) %>% 
  unnest(Synonyms) %>% 
  filter(Symbol %in% yoshida_target | Synonyms %in% yoshida_target) ->a
  







005_reference_GSE84422_STG_modules_BN_.R





####?㑤?????̃T???v?????p????BN####
rm(list=ls())
setwd("\\\\cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/02_CNS/001_AD/2019XX_AD_patient_GSE84422_reanalysis/analysis/module analysis")

ME=read.table("mergedMEs_signed_ver.2.csv",sep=",",header = T,row.names = 1,stringsAsFactors = F)
ME_1=cbind(data.frame("region"=rownames(ME)),ME)
ME_1=subset(ME_1,grepl("X05",ME_1$region))

Trait=read.table("Traits.csv",sep=",",header = T,row.names = 1)
Trait_1=cbind(data.frame("region"=rownames(Trait)),Trait)
Trait_1=subset(Trait_1,grepl("X05",Trait_1$region))

MergedData=merge(Trait_1,ME_1,by.x = "region",by.y = "region")

rownames(MergedData)=MergedData$region
MergedData=MergedData[,-1]
colnames(MergedData)[2]="Cerad"
colnames(MergedData)[3]="Braak"

MergedData_1=as.data.frame(t(MergedData))
MergedData_1=MergedData_1[-which(rownames(MergedData_1)=="grey"),]

dir.create("Superior temporal gyrus")
write.table(MergedData_1,"Superior temporal gyrus/exp_header_rownames.txt",sep="\t",row.names = T,col.names = T)
write.table(MergedData_1,"Superior temporal gyrus/exp.txt",sep="\t",row.names = F,col.names = F)

#LL??Mcmc?ɑ΂??郊?X?g???쐬
library(bnlearn)
sampler=c("1PB","2PB","3PB","4PB","c2PB","c3PB","c4PB","M.1PB","M.2PB","M.3PB","M.4PB","M.c2PB","M.c3PB","M.c4PB","M.REV50","M.STR","REV50","STR")
node=rownames(MergedData_1)
LL_list=list()
Mcmc_list=list()
for(i in c(1:length(sampler))){
  #filename?̒??`
  file_in_LL=paste("Superior temporal gyrus/GSE84422_STG_BN/result/LL_",sampler[i],".csv",sep="")
  file_in_Mcmc=paste("Superior temporal gyrus/GSE84422_STG_BN/result/Mcmc_",sampler[i],".csv",sep="")
  
  #LL
  LL=read.table(file_in_LL,sep=",",header=F,row.names = NULL,stringsAsFactors = F)
  LL=data.frame("iteration"=c(1:1000),LL)
  LL_list=c(LL_list,list(LL))
  
  #Mcmc
  Mcmc=read.table(file_in_Mcmc,sep=",",header=F,row.names = NULL,stringsAsFactors = F)
  rownames(Mcmc)=c(node)
  colnames(Mcmc)=c(node)
  Mcmc_mat=as.matrix(Mcmc)
  dag=empty.graph(c(node))
  amat(dag)=Mcmc_mat
  Mcmc_list=c(Mcmc_list,list(dag))
}
#rename
names(LL_list)=sampler
names(Mcmc_list)=sampler
#save
save(LL_list,file="Superior temporal gyrus/GSE84422_STG_BN/result/LL_list.rda")
save(Mcmc_list,file="Superior temporal gyrus/GSE84422_STG_BN/result/Mcmc_list.rda")

#LL ?o??
dir.create("Superior temporal gyrus/BN_results")
pdf("Superior temporal gyrus/BN_results/LL.pdf",paper="a4",width=9.5,height = 7,pointsize = 10)
par(mfrow=c(1,2))
for(i in c(1:length(sampler))){
  main=sampler[i]
  plot(LL_list[[i]],xlab="iteration",ylab="Log Likelihood",main=main)
}
dev.off()

#max_LL???o
max_LL=sapply(1:length(sampler),function(x){max(LL_list[[x]]$V1)})
max_LL=data.frame("sampler"=c(sampler),"max_LL"=max_LL)
write.table(max_LL,"Superior temporal gyrus/BN_results/max_LL.csv",sep=",",row.names = F,col.names = T)

#DAG ?o??
pdf(paste0("Superior temporal gyrus/BN_results/DAG.pdf"),paper="a4",width=9.5,height = 7,pointsize = 10)
par(mfrow=c(1,2))
for(i in c(1:length(sampler))){
  main=sampler[i]
  g=Rgraphviz::layoutGraph(bnlearn::as.graphNEL(Mcmc_list[[i]]))
  graph::parRenderInfo(g)=list(graph=list(main=main),nodes=list(fontsize=30))
  Rgraphviz::renderGraph(g)
}
dev.off()

###########BMA#################
#c3PB
BMA_pre=read.table("Superior temporal gyrus/GSE84422_STG_BN/result/BMA_pre_re.csv",sep=",",header = F,row.names = NULL,stringsAsFactors = F)
data=read.table("Superior temporal gyrus/exp_header_rownames.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
rownames(BMA_pre)=rownames(data)
colnames(BMA_pre)=rownames(data)
for (i in 1:nrow(BMA_pre)){
  for (j in 1:ncol(BMA_pre)){
    if (BMA_pre[i,j]>30 && BMA_pre[j,i]>30){
      if (BMA_pre[i,j]>BMA_pre[j,i]){
        BMA_pre[j,i]=0
      }else if (BMA_pre[i,j]<BMA_pre[j,i]){
        BMA_pre[i,j]=0
      }else{
        message(print("node",i,"and node",j,"is undirected"))
      }
    }
  }
} #?Ō???eles?͂Ȃ?

#???Ɗ??S?ɓƗ??????m?[?h???폜
d_node=c()
for (i in 1:ncol(BMA_pre)){
  if (sum(BMA_pre[i,])==0 && sum(BMA_pre[,i])==0){
    d_node=c(d_node,i)
  }
} #d_node is empty 

BMA=BMA_pre
for (i in 1:nrow(BMA)){
  for (j in 1:ncol(BMA)){
    if (BMA[i,j]>1){
      BMA[i,j]=1
    }
  }
}
#DAG?̊m?F
dag=empty.graph(rownames(BMA))
amat(dag)=as.matrix(BMA) #the specified network contains cycles.

for (i in 3:ncol(BMA)){
  print(i)
  dag=empty.graph(row.names(BMA)[1:i])
  amat(dag)=as.matrix(BMA[1:i,1:i])
} #42 stop

node=c()
for (i in 1:ncol(BMA)){
  tryCatch({
    BMA_1=BMA[-i,-i]
    dag=empty.graph(rownames(BMA_1))
    amat(dag)=as.matrix(BMA_1)
    node=c(node,i)
    rm(BMA_1)
  },error=function(e){print(i)})
}
node #22 36 42

in_36=which(BMA[,36]==1) # 12 35 39
out_42=which(BMA[42,]==1) #5  7 10 12 20 28 34 35 41

for (j in in_36){
  BMA_2=BMA[-j,-j]
  node=c()
  for (i in 1:dim(BMA_2)[1]){
    tryCatch({
      BMA_1=BMA_2[-i,-i]
      dag=empty.graph(rownames(BMA_1))
      amat(dag)=as.matrix(BMA_1)
      node=c(node,i)
      rm(BMA_1)
    },error=function(e){print(i)})
  }
  assign(paste0("node_",j),node)
  rm(BMA_2)
  rm(node)
}

BMA_pre[c(22,35,36,42),c(22,35,36,42)]
BMA_pre[c(12,22,36,42),c(12,22,36,42)]

#36->22????
BMA_final=BMA
BMA_final[36,22]=0
dag=empty.graph(rownames(BMA_final))
amat(dag)=as.matrix(BMA_final)
#ok
pdf("Superior temporal gyrus/BN_results/BMA.pdf",paper="a4",width=9.5,height = 7,pointsize = 10)
g=Rgraphviz::layoutGraph(bnlearn::as.graphNEL(dag))
graph::parRenderInfo(g)=list(graph=list(main="c3PB BMA"),nodes=list(fontsize=30))
Rgraphviz::renderGraph(g)
dev.off()

#???????킩?????̂ŉ??߂ł??郂?W???[?????I??
#module-module????
rm(list=ls())
data=read.table("Superior temporal gyrus/exp_header_rownames.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
ME=as.data.frame(t(data[-c(1:3),]))

library(WGCNA)
module_cor = cor(ME, use = "p",method ="pearson")
#pvalue?Z?o
nSamples=dim(ME)[1]
module_score_Pvalue = corPvalueStudent(module_cor, nSamples)

module_cor_1=module_cor
for (i in 1:dim(module_score_Pvalue)[1]){
  for (j in 1:dim(module_score_Pvalue)[2]){
    if (module_score_Pvalue[i,j]>=0.05){
      module_cor_1[i,j]=0
    }
  }
}

for (i in 1:dim(module_cor_1)[1]){
  module_cor_1[i,i]=0
}

install.packages("sna")
library(sna)
gplot(module_cor_1,usearrows = F,displaylabels = T)






007.R











007_Barplot_RankModulePreservation.R





library(tidyverse)
library(WGCNA)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")


lname<- load("../2.Output/07_ModulePreservation/Testset_GSE45485/modulePreservationGSE58095_GSE45485.Rdata")
mp$preservation
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor

mp_quality<-read_csv("../2.Output/07_ModulePreservation/Testset_GSE45485/Table_GSE58095_GSE45485_preservation_quality.csv")

#############################################
#test
mp_quality %>% 
  ggplot(aes(x=module %>% reorder(Zsummary.pres), y= Zsummary.pres))+
    geom_bar(aes(fill=module),stat = "identity")+
    theme_minimal()+theme(panel.grid=element_blank())+
    theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
          panel.border = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
          #legend.box.background = element_rect(fill = "transparent"),
          #legend.position = c(1,0.1),
          #legend.justification = c(0,0),
          #axis.text=element_text(size=16),
          axis.text.y = element_text(size=12,lineheight = 1.5),
          axis.title.x = element_text(size=12, color= "#5f5f5f"),
          plot.margin=unit(c(1,1,1,1),"cm"))+
    scale_fill_identity()+   #,"#9bbb59"
    scale_y_continuous(expand = c(0, 0))+
    coord_flip()+
    xlab("")
##############################################

my_barplot_fromMP<-function(mp_quality){
  mp_quality %>% 
    filter(module %in% significant_module) %>% 
    arrange(desc(Zsummary.pres)) %>% 
    head(5) %>% 
    ggplot(aes(x=module %>% reorder(-Zsummary.pres), y= Zsummary.pres))+
    geom_bar(aes(fill=module),stat = "identity")+
    geom_hline(yintercept = 10,color= "#da6272",alpha=0.7,size=1,linetype="dashed")+
    theme_minimal()+theme(panel.grid=element_blank())+
    theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
          panel.border = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
          #legend.box.background = element_rect(fill = "transparent"),
          #legend.position = c(1,0.1),
          #legend.justification = c(0,0),
          #axis.text=element_text(size=16),
          axis.text.x = element_text(size=16,lineheight = 1.5,angle = 30,hjust = 1),
          axis.title.y = element_text(size=16, color= "#5f5f5f"),
          plot.margin=unit(c(1,1,1,1),"cm"))+
    scale_fill_identity()+   #,"#9bbb59"
    scale_y_continuous(expand = c(0, 0))+
    #coord_flip()+
    xlab("")+
    ylab("Preservation Zsummary")
}


mp_quality %>% 
  filter(module %in% significant_module) %>% 
  arrange(desc(Zsummary.pres)) %>% 
  head(5) %>% 
  ggplot(aes(x=module %>% reorder(-Zsummary.pres), y= Zsummary.pres))+
  geom_bar(aes(fill=module),stat = "identity")+
  geom_hline(yintercept = 10,color= "#da6272",alpha=0.7,size=1,linetype="dashed")+
  theme_minimal()+theme(panel.grid=element_blank())+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent"),
        #legend.position = c(1,0.1),
        #legend.justification = c(0,0),
        #axis.text=element_text(size=16),
        axis.text.x = element_text(size=16,lineheight = 1.5,angle = 30,hjust = 1),
        axis.title.y = element_text(size=16, color= "#5f5f5f"),
        plot.margin=unit(c(1,1,1,1),"cm"))+
  scale_fill_identity()+   #,"#9bbb59"
  scale_y_continuous(expand = c(0, 0))+
  #coord_flip()+
  xlab("")+
  ylab("Preservation Zsummary")
ggsave("../2.Output/07_ModulePreservation/Zsummary_Rank_bar_ver2.png",
       dpi = 320, width = 130, height = 130, units = "mm",  bg = "transparent")


mp_quality<-read_csv("../2.Output/18_Comparison_withMP/Testset_GSE81071_CLE/Table_GSE58095_GSE81071_preservation_quality.csv")
my_barplot_fromMP(mp_quality)
ggsave("../2.Output/18_Comparison_withMP/Zsummary_Rank_bar_vsCLE_ver1.png",
       dpi = 320, width = 130, height = 130, units = "mm",  bg = "transparent")

#########################################################################################
#pres ?? Qual?ǂ????g?????̊m?F
ref = 1
test = 2

modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
# leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
# Text labels for points
text = modColors[plotMods];
# Auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
mp$preservation$observed[[ref]][[test]][, 0:2]

mains = c("Preservation Median rank", "Preservation Zsummary");
# Start the plot
sizeGrWindow(10, 5);
#pdf(fi="Plots/BxHLiverFemaleOnly-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
for (p in 1:2)
{
  min = min(plotData[, p], na.rm = TRUE);
  max = max(plotData[, p], na.rm = TRUE);
  # Adjust ploting ranges appropriately
  if (p==2)
  {
    if (min > -max/10) min = -max/10
    ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
  } else
    ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
  plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
       main = mains[p],
       cex = 2.4,
       ylab = mains[p], xlab = "Module size", log = "x",
       ylim = ylim,
       xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
  labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
  # For Zsummary, add threshold lines
  if (p==2)
  {
    abline(h=0)
    abline(h=2, col = "blue", lty = 2)
    abline(h=10, col = "darkgreen", lty = 2)
  }
}
# If plotting into a file, close it
dev.off();






007_CalculateModulePreservation.R





#http://www.iu.a.u-tokyo.ac.jp/~kadota/r.html
#AnnotationDBi?̃_?E?????[?h?ߒ???vtcr,backports??install::packages??install?ł?????????????
#https://cran.r-project.org/web/packages/backports/index.html
#????Cran????tar?t?@?C???_?E?????[?h???Ă??Ė??????????ꂽ
#patchwork?̃_?E?????[?h?͈ȉ?????master???_?E?????[?h????devtools::install_local("path/to/masterfile")
#https://github.com/thomasp85/patchwork #proxy???ʂ???devtools::install_github("thomasp85/patchwork")??OK?Ȃ͂?
#patchwork reference: https://gotellilab.github.io/GotelliLabMeetingHacks/NickGotelli/ggplotPatchwork.html
library("hgug4112a.db", character.only=T)
library(tidyverse)
library(WGCNA)
library(limma)
library(patchwork)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

refset <- "GSE58095"
testset<-"GSE22356"

output_dir_path<-paste0("../2.Output/07_ModulePreservation/Testset_",testset,"_PBMC/")

datExprRef_color_Rdata_path <- paste0("../2.Output/07_ModulePreservation/Referenceset_GSE58095/datExprRef_color_",refset,".Rdata")
datExprTest_color_path<-paste0(output_dir_path,"datExprTest_color_",testset,".Rdata")
modulepreservation_result_path<-paste0(output_dir_path,"modulePreservation",refset,"_", testset,".Rdata")
preservation_quality_table_path<-paste(output_dir_path,"Table_", refset, "_", testset, "_preservation_quality.csv", sep = "")

ZS_result_image_path<-paste0(output_dir_path,"modulePreservation_",refset,"_", testset, "_Zsummary.png")
MR_result_image_path<-paste0(output_dir_path,"modulePreservation_",refset,"_", testset, "_medianRank.png")
ZS_MR_result_image_path <- paste0(output_dir_path, "modulePreservation_",refset,"_", testset, "_Zsummary-medianRank.png")


significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor

color_value<-c("#da6272","#0086ab","#777777","#f79646","#0086ab","#9bbb59","#bfbfbf")
#######################################################################################
#module preservation ?v?Z
lname_Ref<-load(datExprRef_color_Rdata_path) #"datExprRef" "colorsRef"
lname_Test<-load(datExprTest_color_Rdata_path) #"datExprTest" "colorsTest" 

setLabels <- c(refset, testset)
multiExpr <- list(refset = list(data = datExprRef), testset = list(data = datExprTest))
multiColor <- list(refset = colorsRef)

#module preservation???Z?o
#????networkType??signed???w??
system.time({
  mp <- modulePreservation(multiExpr, multiColor,
                           networkType = "signed",
                           referenceNetworks = 1,
                           nPermutations = 200,
                           randomSeed = 1,
                           quickCor = 0,
                           verbose = 3)
})
save(mp,file = modulepreservation_result_path)
#######################################################################################
#visualize the result

load(modulepreservation_result_path)
#mp
ref <- 1
test <- 2
statsObs <- cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ <- cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])
#preservation medianRank??Zsummary???\??
#preservation??quality?????r
preservation_quality <- cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
                              signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))%>% 
  rownames_to_column("module") %>% 
  left_join(mp$preservation$Z[[ref]][[test]] %>% dplyr::select(moduleSize) %>% rownames_to_column("module"),by="module")
write_csv(preservation_quality,preservation_quality_table_path )

#Zsummary.pres
g_zs<-read_csv(preservation_quality_table_path) %>% 
  filter(!(module %in% c("gold","grey"))) %>% 
  ggplot(aes(x=moduleSize,y=Zsummary.pres))+
  #geom_rect(xmin=-Inf, xmax=Inf, ymin=10, ymax=Inf,alpha=0.8,fill="#da6272")+
  #geom_rect(xmin=-Inf, xmax=Inf, ymin=2, ymax=10,alpha=0.8,fill="#E9A1AA")+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = 2)+
  geom_hline(yintercept = 10)+
  geom_point(aes(color=module),size=5)+
  geom_text(aes(label=module),nudge_y = -0.4)+
  scale_color_identity()+
  scale_x_log10(limits=c(4,450))+
  labs(x="Module Size",y="Preservation Zsummary",title = "Preservation Zsummary")+
  theme_classic()+
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"))
g_zs
ggsave(g_zs,filename = ZS_result_image_path, 
       dpi=2000, width = 140, height = 140, units = "mm",  bg = "transparent")


#medianRank.pres
g_mr<-read_csv(preservation_quality_table_path) %>% 
  filter(!(module %in% c("gold","grey"))) %>% 
  ggplot(aes(x=moduleSize,y=medianRank.pres))+
  #geom_hline(yintercept = 2)+
  #geom_hline(yintercept = 10)+
  geom_point(aes(color=module),size=5)+
  geom_text(aes(label=module),nudge_y = -1)+
  scale_color_identity()+
  scale_x_log10(limits=c(4,450))+
  scale_y_reverse(limits=c(40,-0))+
  #ylim(-5,45)+
  labs(x="Module Size",y="Preservation Median Rank",title = "Preservation Median Rank")+
  theme_classic()+
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"))
g_mr
ggsave(g_mr,filename = MR_result_image_path, 
       dpi=2000, width = 140, height = 140, units = "mm",  bg = "transparent")

g_zs + g_mr
ggsave(g_zs + g_mr, filename =  ZS_MR_result_image_path,
       dpi=2000, width = 280, height = 140, units = "mm",  bg = "transparent")

#########################################################################################
#significant_module
#Zsummary.pres
g_sm_zs<-read_csv(preservation_quality_table_path) %>% 
  filter(!(module %in% c("gold","grey"))) %>% 
  filter(module %in% significant_module) %>% 
  ggplot(aes(x=moduleSize,y=Zsummary.pres))+
  #geom_rect(xmin=-Inf, xmax=Inf, ymin=10, ymax=Inf,alpha=0.8,fill="#da6272")+
  #geom_rect(xmin=-Inf, xmax=Inf, ymin=2, ymax=10,alpha=0.8,fill="#E9A1AA")+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = 2)+
  geom_hline(yintercept = 10)+
  geom_point(aes(color=module),size=5)+
  geom_text(aes(label=module),nudge_y = -0.4)+
  scale_color_identity()+
  scale_x_log10(limits=c(4,450))+
  labs(x="Module Size",y="Preservation Zsummary",title = "Preservation Zsummary")+
  theme_classic()+
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"))
g_sm_zs
ggsave(g_sm_zs,filename = paste0("../2.Output/07_ModulePreservation/modulePreservation_",refset,"_", testset,"_SignificantModule", "_Zsummary.png"), 
       dpi=2000, width = 140, height = 140, units = "mm",  bg = "transparent")


#medianRank.pres
g_sm_mr<-read_csv(preservation_quality_table_path) %>% 
  filter(!(module %in% c("gold","grey"))) %>% 
  filter(module %in% significant_module) %>% 
  ggplot(aes(x=moduleSize,y=medianRank.pres))+
  #geom_hline(yintercept = 2)+
  #geom_hline(yintercept = 10)+
  geom_point(aes(color=module),size=5)+
  geom_text(aes(label=module),nudge_y = -1)+
  scale_color_identity()+
  scale_x_log10(limits=c(4,450))+
  scale_y_reverse(limits=c(40,-0))+
  #ylim(-5,45)+
  labs(x="Module Size",y="Preservation Median Rank",title = "Preservation Median Rank")+
  theme_classic()+
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"))
g_sm_mr
ggsave(g_sm_mr,filename = paste0(output_dir_path,"modulePreservation_",refset,"_", testset,"_SignificantModule", "_medianRank.png"), 
       dpi=2000, width = 140, height = 140, units = "mm",  bg = "transparent")

g_sm_zs + g_sm_mr
#output_image_path<-paste0("../2.Output/07_ModulePreservation/modulePreservation_",refset,"_", testset, "_Zsummary-medianRank.png")
ggsave(g_sm_zs + g_sm_mr,filename = paste0(output_dir_path,"modulePreservation_",refset,"_", testset, "_SignificantModule","_Zsummary-medianRank.png"),
       dpi=2000, width = 280, height = 140, units = "mm",  bg = "transparent")







007_CalculateModulePreservation_test.R





#http://www.iu.a.u-tokyo.ac.jp/~kadota/r.html
#AnnotationDBi?̃_?E?????[?h?ߒ???vtcr,backports??install::packages??install?ł?????????????
#https://cran.r-project.org/web/packages/backports/index.html
#????Cran????tar?t?@?C???_?E?????[?h???Ă??Ė??????????ꂽ
#patchwork?̃_?E?????[?h?͈ȉ?????master???_?E?????[?h????devtools::install_local("path/to/masterfile")
#https://github.com/thomasp85/patchwork #proxy???ʂ???devtools::install_github("thomasp85/patchwork")??OK?Ȃ͂?
#patchwork reference: https://gotellilab.github.io/GotelliLabMeetingHacks/NickGotelli/ggplotPatchwork.html
library("hgug4112a.db", character.only=T)
library(tidyverse)
library(WGCNA)
library(limma)
library(patchwork)
#library(cowplot)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

refset <- "GSE58095"
testset<-"GSE9285"

GSE9285_annotation_path<-"../1.Data/GSE9285/Annotation_GSE9285.txt"
ProcessedData_path<-"../1.Data/GSE9285/ProcessedData/"
testExpdata_path<-paste0("../1.Data/GSE9285/","Expdata.csv")
refExpdata_path<-"../2.Output/01_Compare_deepsplit_signed/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/Expdata_quontile.Rdata"

geneid2modulecolor_path<-"../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv"

geneid2modulecolor<-read_csv(geneid2modulecolor_path) %>% 
  select(Entrez_Gene_ID,module) %>% 
  distinct(Entrez_Gene_ID,module)
color_value<-c("#da6272","#0086ab","#777777","#f79646","#0086ab","#9bbb59","#bfbfbf")

#######################################################################################
#Refset?̉??H

GSE58095_annotation_path<-"../1.Data/GSE58095/Annotation/GeneAnnotation_GPL10558.csv"
GSE58095_annotation<-read_csv(GSE58095_annotation_path) %>% 
  select(ID,Entrez_Gene_ID) %>% 
  left_join(geneid2modulecolor,by="Entrez_Gene_ID") %>% 
  filter(!is.na(Entrez_Gene_ID),!is.na(module))
lname<-load(refExpdata_path)
#Expdata_quontile
datExprRef_pre <- Expdata_quontile$E
keep<-rowSums(datExprRef_pre>log2(50))>=102
datExprRef_pre<-datExprRef_pre[keep,] %>% 
  as.data.frame() %>% 
  rownames_to_column("ILMNID") %>% 
  left_join(GSE58095_annotation,by=c("ILMNID"="ID")) %>% 
  mutate(exprSignal_median=select(., matches("(?:\\d){3}")) %>% apply( 1, median))%>% 
  arrange(Entrez_Gene_ID %>% as.numeric(),desc(exprSignal_median))%>% 
  distinct(Entrez_Gene_ID,.keep_all = T) %>% 
  filter(!is.na(Entrez_Gene_ID)) %>% 
  column_to_rownames("Entrez_Gene_ID") 
datExprRef<-datExprRef_pre%>% 
  select(-c(ILMNID,exprSignal_median,module))%>% t()
dim(datExprRef)
colorsRef <- datExprRef_pre$module 
datExprRef %>% as.tibble()
#colnames(datExprRef) %>% length()
#datExprRef_pre %>% rownames() %>% length()
#(colnames(datExprRef) ==datExprRef_pre %>% rownames() ) %>% all()
#######################################################################################
#Testset?̉??H
datExprTest_pre<-read_csv(testExpdata_path) %>% 
  left_join(geneid2modulecolor,by=c("ENTREZID"="Entrez_Gene_ID")) %>% 
  distinct(ENTREZID,.keep_all = T)
table(datExprTest_pre$module)
table(geneid2modulecolor$module)
#loss????geneid?̊???
((table(geneid2modulecolor$module)-table(datExprTest_pre$module))/table(geneid2modulecolor$module)) %>% round(2)

#datExprTest<-datExprTest_pre[,2:length(datExprTest_pre)-1] %>% t()
datExprTest<-datExprTest_pre%>% select(starts_with("GSM"))%>% t()
dim(datExprTest)
colnames(datExprTest) <- datExprTest_pre$ENTREZID
datExprTest[1:10, 1:10]
colorsTest <- datExprTest_pre$module 


#######################################################################################
#module preservation ?v?Z
datExprRef_color_path <- "../2.Output/07_ModulePreservation/datExprRef_color_GSE58095.Rdata"
datExprTest_color_path<-"..//2.Output/07_ModulePreservation/datExprTest_color_GSE9285.Rdata"

lname<-load(datExprRef_color_path)
load(datExprTest_color_path)

setLabels <- c(refset, testset)
multiExpr <- list(refset = list(data = datExprRef), testset = list(data = datExprTest))
multiColor <- list(refset = colorsRef)

#module preservation???Z?o
#????networkType??signed???w??
system.time({7
  mp <- modulePreservation(multiExpr, multiColor,
                           networkType = "signed",
                           referenceNetworks = 1,
                           nPermutations = 200,
                           randomSeed = 1,
                           quickCor = 0,
                           verbose = 3)
})
modulepreservation_result_path<-paste("../2.Output/07_ModulePreservation/modulePreservation",refset,"_", testset,"2",".Rdata")
save(mp,file = modulepreservation_result_path)
#######################################################################################
#visualize the result

#reference: tutorial 
#######################################################################################
#https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/ModulePreservation/Tutorials/MiniTutorial-MouseLiver.pdf
if(0){
  load(modulepreservation_result_path)
  mp
  ref <- 1
  test <- 2
  statsObs <- cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
  statsZ <- cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])
  #preservation medianRank??Zsummary???\??
  #preservation??quality?????r
  preservation_quality <- cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
                                signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))
  write.csv(preservation_quality, paste("Table_", refset, "_", testset, "_preservation_quality.csv", sep = ""))
  
  # Module labels and module sizes are also contained in the results
  modColors = rownames(mp$preservation$observed[[ref]][[test]])
  moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
  # leave grey and gold modules out
  plotMods = !(modColors %in% c("grey"));
  # Text labels for points
  text = modColors[plotMods];
  # Auxiliary convenience variable
  plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
  # Main titles for the plot
  mains = c("Preservation Median rank", "Preservation Zsummary");
  # Start the plot
  sizeGrWindow(10, 5);
  output_image_path<-paste0("../2.Output/07_ModulePreservation/modulePreservation_",refset,"_", testset, "_Zsummary-medianRank.png")
  png(filename = output_image_path, width=960, height=480)
  par(mfrow = c(1,2))
  par(mar = c(4.5,4.5,2.5,1))
  for (p in 1:2)
  {
    min = min(plotData[, p], na.rm = TRUE);
    max = max(plotData[, p], na.rm = TRUE);
    # Adjust ploting ranges appropriately
    if (p==2)
    {
      if (min > -max/10) min = -max/10
      ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
    } else
      ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
    plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
         main = mains[p],
         cex = 2.4,
         ylab = mains[p], xlab = "Module size", log = "x",
         ylim = ylim,
         xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
    labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
    # For Zsummary, add threshold lines
    if (p==2)
    {
      abline(h=0)
      abline(h=2, col = "blue", lty = 2)
      abline(h=10, col = "darkgreen", lty = 2)
    }
  }
  # If plotting into a file, close it
  dev.off()
}
#######################################################################################
preservation_quality_table_path<-paste("../2.Output/07_ModulePreservation/","Table_", refset, "_", testset, "_preservation_quality2.csv", sep = "")

load(modulepreservation_result_path)
#mp
ref <- 1
test <- 2
statsObs <- cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ <- cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])
#preservation medianRank??Zsummary???\??
#preservation??quality?????r
preservation_quality <- cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
                              signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))%>% 
  rownames_to_column("module") %>% 
  left_join(mp$preservation$Z[[ref]][[test]] %>% select(moduleSize) %>% rownames_to_column("module"),by="module")
write_csv(preservation_quality,preservation_quality_table_path )

#Zsummary.pres
g_zs<-read_csv(preservation_quality_table_path) %>% 
  filter(!(module %in% c("gold","grey"))) %>% 
  ggplot(aes(x=moduleSize,y=Zsummary.pres))+
  #geom_rect(xmin=-Inf, xmax=Inf, ymin=10, ymax=Inf,alpha=0.8,fill="#da6272")+
  #geom_rect(xmin=-Inf, xmax=Inf, ymin=2, ymax=10,alpha=0.8,fill="#E9A1AA")+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = 2)+
  geom_hline(yintercept = 10)+
  geom_point(aes(color=module),size=5)+
  geom_text(aes(label=module),nudge_y = -0.4)+
  scale_color_identity()+
  scale_x_log10(limits=c(4,450))+
  labs(x="Module Size",y="Preservation Zsummary",title = "Preservation Zsummary")+
  theme_classic()+
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"))
g_zs
#output_image_path<-paste0("../2.Output/07_ModulePreservation/modulePreservation_",refset,"_", testset, "_Zsummary.png")

ggsave(g_zs,filename = paste0("../2.Output/07_ModulePreservation/modulePreservation_",refset,"_", testset, "_Zsummary.png"), 
       dpi=2000, width = 140, height = 140, units = "mm",  bg = "transparent")


#medianRank.pres
g_mr<-read_csv(preservation_quality_table_path) %>% 
  filter(!(module %in% c("gold","grey"))) %>% 
  ggplot(aes(x=moduleSize,y=medianRank.pres))+
  #geom_hline(yintercept = 2)+
  #geom_hline(yintercept = 10)+
  geom_point(aes(color=module),size=5)+
  geom_text(aes(label=module),nudge_y = -1)+
  scale_color_identity()+
  scale_x_log10(limits=c(4,450))+
  scale_y_reverse(limits=c(40,-0))+
  #ylim(-5,45)+
  labs(x="Module Size",y="Preservation Median Rank",title = "Preservation Median Rank")+
  theme_classic()+
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"))
#cowplot::plot_grid(g_zs,g_mr)
g_mr
#output_image_path<-paste0("../2.Output/07_ModulePreservation/modulePreservation_",refset,"_", testset, "_medianRank.png")
ggsave(g_mr,filename = paste0("../2.Output/07_ModulePreservation/modulePreservation_",refset,"_", testset, "_medianRank.png"), 
       dpi=2000, width = 140, height = 140, units = "mm",  bg = "transparent")

g_zs + g_mr
#output_image_path<-paste0("../2.Output/07_ModulePreservation/modulePreservation_",refset,"_", testset, "_Zsummary-medianRank.png")
ggsave(g_zs + g_mr,filename = paste0("../2.Output/07_ModulePreservation/modulePreservation_",refset,"_", testset, "_Zsummary-medianRank.png"),
       dpi=2000, width = 280, height = 140, units = "mm",  bg = "transparent")










007_ModulePreservation_test.R





#http://www.iu.a.u-tokyo.ac.jp/~kadota/r.html
#AnnotationDBi?̃_?E?????[?h?ߒ???vtcr,backports??install::packages??install?ł?????????????
#https://cran.r-project.org/web/packages/backports/index.html
#????Cran????tar?t?@?C???_?E?????[?h???Ă??Ė??????????ꂽ
library("hgug4112a.db", character.only=T)
library(tidyverse)
library(WGCNA)
library(limma)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

GSE9285_annotation_path<-"../1.Data/GSE9285/Annotation_GSE9285.txt"
ProcessedData_path<-"../1.Data/GSE9285/ProcessedData/"
testExpdata_path<-paste0("../1.Data/GSE9285/","Expdata.csv")
refExpdata_path<-"../2.Output/01_Compare_deepsplit_signed/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/Expdata_quontile.Rdata"

#out_f <- "../1.Data/GSE9285/Annotation_GSE9285.txt"                   #?o?̓t?@?C???????w?肵??out_f?Ɋi?[
#param1 <- "hgug4112a.db"               #?p?b?P?[?W?????w??
param <- c("REFSEQ", "SYMBOL","ENTREZID","ENSEMBL")        #?~?????A?m?e?[?V???????????w???i?L?q?\?ȃ??X?g?͈ȉ???keytypes(dbname)?????擾?\?j


####################################################################################################
#Annotation ?t?@?C???̍쐬
#?K?v?ȃp?b?P?[?W?̃C???X?g?[??(???p???????p?b?P?[?W?????ɑ??݂??Ă?????2???ڈȍ~?͕K?v?Ȃ?)
#source("http://bioconductor.org/biocLite.R")#?w?肵???p?b?P?[?W?̃C???X?g?[??
#biocLite(param1, suppressUpdates=TRUE) #?w?肵???p?b?P?[?W?̃C???X?g?[??
#BiocManager::install(param1)
#?K?v?ȃp?b?P?[?W?????[?h
#library(param1, character.only=T)      #param1?Ŏw?肵???p?b?P?[?W?̓ǂݍ???

#?O????(???o?\?ȏ????̃??X?g?A?b?v?????уL?[???\??)
#hoge <- eval(parse(text=param1))       #?I?u?W?F?N?g???̕ύX?ihoge?Ƃ??Ď??舵???j
#cols(hoge)                             #hoge???Ɋ܂܂??Ă??????????\??
keytypes(hgug4112a.db)                         #hoge???????ۂɎg?p?????L?[???\??
param_key <- keys(hgug4112a.db, keytype="PROBEID")#param1?Ŏw?肵?????̂??L?[?ɂ??Ă???

#?{??
out <- AnnotationDbi::select(hgug4112a.db, keys=param_key, keytype="PROBEID", columns=param)#?A?m?e?[?V???????????o???ʂ?out?Ɋi?[
write.table(out, GSE9285_annotation_path, sep="\t", append=F, quote=F, row.names=F)#out?̒??g???w?肵???t?@?C?????ŕۑ?

####################################################################################################
#Expdata?̍쐬

ProcessedData_files<-list.files(ProcessedData_path,full.names = T)
annotation<-read_tsv(GSE9285_annotation_path) %>% 
  select(-REFSEQ) %>% 
  distinct()

result<-tibble()
for (file in ProcessedData_files){
  result<-result %>% 
    bind_rows(read_delim(file,delim = "\t",skip = ) %>% 
                mutate(sample=basename(file) %>% str_split("_") %>% .[[1]] %>% .[1]))
}

RG<- result %>% 
  spread(key = sample,value = VALUE) %>% 
  inner_join(annotation,by=c(`Reporter Identifier`="PROBEID")) %>% 
  select(`Reporter Identifier`,"SYMBOL","ENTREZID","ENSEMBL",everything())%>% 
  mutate(exprSignal_median=dplyr::select(.,starts_with("GSM")) %>%  apply( 1, median,na.rm=T)) %>% #select(exprSignal_median,everything())->a
  arrange(ENTREZID,desc(exprSignal_median)) %>% 
  distinct(ENTREZID,.keep_all = T) %>% #select(ENTREZID,exprSignal_median)
  filter(!is.na(ENTREZID)) %>% 
  column_to_rownames(var = "ENTREZID") %>% 
  select(starts_with("GSM")) %>% 
  as.matrix() %>% 
  backgroundCorrect.matrix(method="normexp")

Expdata<-normalizeBetweenArrays(RG, method="quantile")
write_csv(Expdata %>% as.data.frame() %>% rownames_to_column(var="ENTREZID"),testExpdata_path)
####################################################################################################
#test by saitoh san's data
if(0){
  setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190403_NAFLD_NASH_patient_WGCNA_GSE49541/analysis/WGCNA")
  dat0_path<-"../signed_4/geneID_color_kWithin_normalizedData_logFC_qvalue.csv"
  dat0 <- read_csv(dat0_path) %>% 
    arrange(Entrez_gene_ID %>% as.numeric(),desc(kWithin))%>% 
    distinct(Entrez_gene_ID,.keep_all = T)
  table(duplicated(dat0$Entrez_gene_ID))
  dat0[, 6:75]
  datExprRef <- dat0[, 6:75]
  datExprRef <- t(datExprRef)
  dim(datExprRef)
  colnames(datExprRef) <- dat0$Entrez_gene_ID
  datExprRef[1:10, 1:10]
  colorsRef <- dat0$moduleColor
  
  
  setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190606_NAFLD_NASH_patient_E_MEXP_3291/analysis/module_preservation")
  saitohsandata_path<-paste0("Sample_normalized_data_logFC_qvalue.csv")
  
  data<-read_csv(saitohsandata_path)%>% 
    mutate(exprSignal_median=select(., c(5:49)) %>% apply( 1, median)) %>% 
    arrange(Entrez_gene_ID %>% as.numeric(),desc(exprSignal_median))%>% 
    #select(exprSignal_median,everything())
    distinct(Entrez_gene_ID,.keep_all = T)
  table(duplicated(data$Entrez_gene_ID))
  data[130:140, 1:5]
  datExprTest <- data[, 5:49]
  datExprTest <- t(datExprTest)
  colnames(datExprTest) <- data$Entrez_gene_ID
  datExprTest[1:5, 130:140]
  dim(datExprTest)
  
  ref2test <- match(colnames(datExprRef), colnames(datExprTest))
  table(is.finite(ref2test))
  refset <- "GSE49541"
  testset <- "E_MEXP_3291"
  setLabels <- c(refset, testset)
  multiExpr <- list(refset = list(data = datExprRef), testset = list(data = datExprTest))
  multiColor <- list(refset = colorsRef)
  
  system.time({
    mp <- modulePreservation(multiExpr, multiColor,
                             networkType = "signed",
                             referenceNetworks = 1,
                             nPermutations = 200,
                             randomSeed = 1,
                             quickCor = 0,
                             verbose = 3)
  })
}
####################################################################################################
geneid2modulecolor_path<-"../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv"

geneid2modulecolor<-read_csv(geneid2modulecolor_path) %>% 
  select(Entrez_Gene_ID,module) %>% 
  #filter(module!="grey") %>% 
  distinct(Entrez_Gene_ID,module)


#datExprTest_pre<-read_csv(testExpdata_path) %>% 
#  mutate(exprSignal_median=select(., c(5:79)) %>% apply( 1, median)) %>% #select(exprSignal_median,everything())
#  arrange(ENTREZID  %>% as.numeric(),desc(exprSignal_median))%>% 
#select(exprSignal_median,everything())
#  distinct(ENTREZID ,.keep_all = T)%>% 
#  left_join(geneid2modulecolor,by=c("ENTREZID"="Entrez_Gene_ID")) %>% 
#  select(module, exprSignal_median,everything())

datExprTest_pre<-read_csv(testExpdata_path) %>% 
  left_join(geneid2modulecolor,by=c("ENTREZID"="Entrez_Gene_ID")) %>% 
  distinct(ENTREZID,.keep_all = T)
table(datExprTest_pre$module)
table(geneid2modulecolor$module)
#loss????geneid?̊???
((table(geneid2modulecolor$module)-table(datExprTest_pre$module))/table(geneid2modulecolor$module)) %>% round(2)

#datExprTest<-datExprTest_pre[,2:length(datExprTest_pre)-1] %>% t()
datExprTest<-datExprTest_pre%>% select(starts_with("GSM"))%>% t()
dim(datExprTest)
colnames(datExprTest) <- datExprTest_pre$ENTREZID
datExprTest[1:10, 1:10]
colorsTest <- datExprTest_pre$module 

#############################

GSE58095_annotation_path<-"../1.Data/GSE58095/Annotation/GeneAnnotation_GPL10558.csv"
GSE58095_annotation<-read_csv(GSE58095_annotation_path) %>% 
  select(ID,Entrez_Gene_ID) %>% 
  left_join(geneid2modulecolor,by="Entrez_Gene_ID") %>% 
  filter(!is.na(Entrez_Gene_ID),!is.na(module))
lname<-load(refExpdata_path)
#Expdata_quontile
datExprRef_pre <- Expdata_quontile$E
keep<-rowSums(datExprRef_pre>log2(50))>=102
datExprRef_pre<-datExprRef_pre[keep,] %>% 
  as.data.frame() %>% 
  rownames_to_column("ILMNID") %>% 
  left_join(GSE58095_annotation,by=c("ILMNID"="ID")) %>% 
  mutate(exprSignal_median=select(., matches("(?:\\d){3}")) %>% apply( 1, median))%>% 
  arrange(Entrez_Gene_ID %>% as.numeric(),desc(exprSignal_median))%>% 
  distinct(Entrez_Gene_ID,.keep_all = T) %>% 
  filter(!is.na(Entrez_Gene_ID)) %>% 
  column_to_rownames("Entrez_Gene_ID") 
datExprRef<-datExprRef_pre%>% 
  select(-c(ILMNID,exprSignal_median,module))%>% t()
dim(datExprRef)
colorsRef <- datExprRef_pre$module 
datExprRef %>% as.tibble()
#colnames(datExprRef) %>% length()
#datExprRef_pre %>% rownames() %>% length()
#(colnames(datExprRef) ==datExprRef_pre %>% rownames() ) %>% all()

refset <- "GSE58095"
testset<-"GSE9285"
setLabels <- c(refset, testset)
multiExpr <- list(refset = list(data = datExprRef), testset = list(data = datExprTest))
multiColor <- list(refset = colorsRef)

#module preservation???Z?o
#????networkType??signed???w??
system.time({
  mp <- modulePreservation(multiExpr, multiColor,
                           networkType = "signed",
                           referenceNetworks = 1,
                           nPermutations = 200,
                           randomSeed = 1,
                           quickCor = 0,
                           verbose = 3)
})






007_Preprocessing.R





#http://www.iu.a.u-tokyo.ac.jp/~kadota/r.html
#AnnotationDBi?̃_?E?????[?h?ߒ???vtcr,backports??install::packages??install?ł?????????????
#https://cran.r-project.org/web/packages/backports/index.html
#????Cran????tar?t?@?C???_?E?????[?h???Ă??Ė??????????ꂽ
library("hgug4112a.db", character.only=T)
library(tidyverse)
library(WGCNA)
library(limma)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

#ReferenceSet related
refExpdata_path<-"../2.Output/01_Compare_deepsplit_signed/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/Expdata_quontile.Rdata"
GSE58095_annotation_path<-"../1.Data/GSE58095/Annotation/GeneAnnotation_GPL10558.csv"
probeid_kwithin_path<-"../2.Output/02_signed_2_unmerged/signed_2/Probe_kWitnin_color_signed_2.csv"
datExprRef_color_path <- "../2.Output/07_ModulePreservation/datExprRef_color_GSE58095.Rdata"

#testset related
GSE9285_annotation_path<-"../1.Data/GSE9285/Annotation_GSE9285.txt"
GSE9285_ProcessedData_path<-"../1.Data/GSE9285/ProcessedData/"
testExpdata_path<-paste0("../1.Data/GSE9285/","Expdata.csv")
datExprTest_color_path<-"..//2.Output/07_ModulePreservation/datExprTest_color_GSE9285_2.Rdata"

geneid2modulecolor_path<-"../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv"

geneid2modulecolor<-read_csv(geneid2modulecolor_path) %>% 
  select(Entrez_Gene_ID,module) %>% 
  distinct(Entrez_Gene_ID,module)
#######################################################################################
#Refset?̉??H


GSE58095_annotation<-read_csv(GSE58095_annotation_path) %>% 
  select(ID,Entrez_Gene_ID) %>% 
  left_join(geneid2modulecolor,by="Entrez_Gene_ID") %>% 
  filter(!is.na(Entrez_Gene_ID),!is.na(module))
GSE58095_probeid2kwithin<-read_csv(probeid_kwithin_path) %>% 
  select(ProbeID,kWithin)

lname<-load(refExpdata_path)
#Expdata_quontile
datExprRef_pre <- Expdata_quontile$E
keep<-rowSums(datExprRef_pre>log2(50))>=102
datExprRef_pre <- datExprRef_pre[keep,] %>% 
  as.data.frame() %>% 
  rownames_to_column("ILMNID") %>% 
  left_join(GSE58095_probeid2kwithin %>% left_join(GSE58095_annotation,by=c("ProbeID"="ID")),by=c("ILMNID"="ProbeID")) %>% 
  filter(!is.na(Entrez_Gene_ID)) %>% 
  #mutate(exprSignal_median=select(., matches("(?:\\d){3}")) %>% apply( 1, median))%>% #exprSignal_median?????GeneID?̏d???????????Ƃ?
  arrange(Entrez_Gene_ID %>% as.numeric(),desc(kWithin))%>% 
  distinct(Entrez_Gene_ID,.keep_all = T) %>% 
  column_to_rownames("Entrez_Gene_ID") 
datExprRef<-datExprRef_pre%>% 
  select(-c(ILMNID,kWithin,module))%>% t()
dim(datExprRef)
colorsRef <- datExprRef_pre$module 
#datExprRef %>% as.tibble()
#colnames(datExprRef) %>% length()
#datExprRef_pre %>% rownames() %>% length()
#(colnames(datExprRef) ==datExprRef_pre %>% rownames() ) %>% all()
save(list = c("datExprRef","colorsRef"),file = datExprRef_color_path)
####################################################################################################
#testset?̉??H

param <- c("REFSEQ", "SYMBOL","ENTREZID","ENSEMBL")        #?~?????A?m?e?[?V???????????w???i?L?q?\?ȃ??X?g?͈ȉ???keytypes(dbname)?????擾?\?j
keytypes(hgug4112a.db)                         #hoge???????ۂɎg?p?????L?[???\??
param_key <- keys(hgug4112a.db, keytype="PROBEID")#param1?Ŏw?肵?????̂??L?[?ɂ??Ă???

#?{??
out <- AnnotationDbi::select(hgug4112a.db, keys=param_key, keytype="PROBEID", columns=param)#?A?m?e?[?V???????????o???ʂ?out?Ɋi?[
write.table(out, GSE9285_annotation_path, sep="\t", append=F, quote=F, row.names=F)#out?̒??g???w?肵???t?@?C?????ŕۑ?

ProcessedData_files<-list.files(GSE9285_ProcessedData_path,full.names = T)
GSE9285_annotation<-read_tsv(GSE9285_annotation_path) %>% 
  select(-REFSEQ) %>% 
  distinct()

ProcessedData<-tibble()
for (file in ProcessedData_files){
  ProcessedData<-ProcessedData %>% 
    bind_rows(read_delim(file,delim = "\t",skip = ) %>% 
                mutate(sample=basename(file) %>% str_split("_") %>% .[[1]] %>% .[1]))
}


RG<- ProcessedData %>% 
  spread(key = sample,value = VALUE) %>% 
  column_to_rownames(var = "Reporter Identifier") %>% 
  as.matrix() %>% 
  backgroundCorrect.matrix(method="normexp")

RGnormalized<-normalizeBetweenArrays(RG, method="quantile") 
datExprTest_pre<-RGnormalized%>% as.data.frame() %>% 
  rownames_to_column("ReporterIdentifier") %>% 
  left_join(GSE9285_annotation %>% 
              left_join(geneid2modulecolor,by=c("ENTREZID"="Entrez_Gene_ID")),by=c("ReporterIdentifier"="PROBEID")) %>% 
  #select(`Reporter Identifier`,"SYMBOL","ENTREZID","ENSEMBL",everything())%>% 
  filter(!is.na(ENTREZID))  %>% 
  filter(!is.na(module)) %>% 
  mutate(exprSignal_median=dplyr::select(.,starts_with("GSM")) %>%  apply( 1, median,na.rm=T)) %>% #select(exprSignal_median,everything())->a
  arrange(ENTREZID,desc(exprSignal_median)) %>% 
  distinct(ENTREZID,.keep_all = T) 


table(datExprTest_pre$module)
table(geneid2modulecolor$module)
#loss????geneid?̊???
((table(geneid2modulecolor$module)-table(datExprTest_pre$module))/table(geneid2modulecolor$module)) %>% round(2)

#datExprTest<-datExprTest_pre[,2:length(datExprTest_pre)-1] %>% t()
datExprTest<-datExprTest_pre%>% select(starts_with("GSM"))%>% t()
dim(datExprTest)
colnames(datExprTest) <- datExprTest_pre$ENTREZID
datExprTest[1:10, 1:10]
colorsTest <- datExprTest_pre$module 

save(list = c("datExprTest","colorsTest"),file = datExprTest_color_path)
####################################################################################################

#testset?̉??H backgroundCorrect.matrix?Ȃ?ver

param <- c("REFSEQ", "SYMBOL","ENTREZID","ENSEMBL")        #?~?????A?m?e?[?V???????????w???i?L?q?\?ȃ??X?g?͈ȉ???keytypes(dbname)?????擾?\?j
keytypes(hgug4112a.db)                         #hoge???????ۂɎg?p?????L?[???\??
param_key <- keys(hgug4112a.db, keytype="PROBEID")#param1?Ŏw?肵?????̂??L?[?ɂ??Ă???

#?{??
out <- AnnotationDbi::select(hgug4112a.db, keys=param_key, keytype="PROBEID", columns=param)#?A?m?e?[?V???????????o???ʂ?out?Ɋi?[
write.table(out, GSE9285_annotation_path, sep="\t", append=F, quote=F, row.names=F)#out?̒??g???w?肵???t?@?C?????ŕۑ?

ProcessedData_files<-list.files(GSE9285_ProcessedData_path,full.names = T)
GSE9285_annotation<-read_tsv(GSE9285_annotation_path) %>% 
  select(-REFSEQ) %>% 
  distinct()

ProcessedData<-tibble()
for (file in ProcessedData_files){
  ProcessedData<-ProcessedData %>% 
    bind_rows(read_delim(file,delim = "\t",skip = ) %>% 
                mutate(sample=basename(file) %>% str_split("_") %>% .[[1]] %>% .[1]))
}

ProcessedData %>% 
  ggplot(aes(x=sample,y=VALUE))+
  geom_boxplot()

#RG<- ProcessedData %>% 
#  spread(key = sample,value = VALUE) %>% 
#  column_to_rownames(var = "Reporter Identifier") %>% 
#  as.matrix() %>% 
#  backgroundCorrect.matrix(method="normexp")
RGnormalized<-ProcessedData %>% 
  spread(key = sample,value = VALUE) %>% 
  column_to_rownames(var = "Reporter Identifier") %>% 
  as.matrix() %>% 
  normalizeBetweenArrays(method="quantile")

RGnormalized %>% as.data.frame() %>% rownames_to_column("ID") %>% 
  gather(key="sample",value="VALUE",-ID)%>% 
  ggplot(aes(x=sample,y=VALUE))+
  geom_boxplot()


datExprTest_pre<-RGnormalized%>% as.data.frame() %>% 
  rownames_to_column("ReporterIdentifier") %>% 
  left_join(GSE9285_annotation %>% 
              left_join(geneid2modulecolor,by=c("ENTREZID"="Entrez_Gene_ID")),by=c("ReporterIdentifier"="PROBEID")) %>% 
  #select(`Reporter Identifier`,"SYMBOL","ENTREZID","ENSEMBL",everything())%>% 
  filter(!is.na(ENTREZID))  %>% 
  filter(!is.na(module)) %>% 
  mutate(exprSignal_median=dplyr::select(.,starts_with("GSM")) %>%  apply( 1, median,na.rm=T)) %>% #select(exprSignal_median,everything())->a
  arrange(ENTREZID,desc(exprSignal_median)) %>% 
  distinct(ENTREZID,.keep_all = T) 


table(datExprTest_pre$module)
table(geneid2modulecolor$module)
#loss????geneid?̊???
((table(geneid2modulecolor$module)-table(datExprTest_pre$module))/table(geneid2modulecolor$module)) %>% round(2)

#datExprTest<-datExprTest_pre[,2:length(datExprTest_pre)-1] %>% t()
datExprTest<-datExprTest_pre%>% select(starts_with("GSM"))%>% t()
dim(datExprTest)
colnames(datExprTest) <- datExprTest_pre$ENTREZID
datExprTest[1:10, 1:10]
colorsTest <- datExprTest_pre$module 

save(list = c("datExprTest","colorsTest"),file = datExprTest_color_path)






007_Preprocessing_test.R





#http://www.iu.a.u-tokyo.ac.jp/~kadota/r.html
#AnnotationDBi?̃_?E?????[?h?ߒ???vtcr,backports??install::packages??install?ł?????????????
#https://cran.r-project.org/web/packages/backports/index.html
#????Cran????tar?t?@?C???_?E?????[?h???Ă??Ė??????????ꂽ
library("hgug4112a.db", character.only=T)
library(tidyverse)
library(WGCNA)
library(limma)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

#ReferenceSet related
refExpdata_path<-"../2.Output/01_Compare_deepsplit_signed/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/Expdata_quontile.Rdata"
GSE58095_annotation_path<-"../1.Data/GSE58095/Annotation/GeneAnnotation_GPL10558.csv"
probeid_kwithin_path<-"../2.Output/02_signed_2_unmerged/signed_2/Probe_kWitnin_color_signed_2.csv"
datExprRef_color_path <- "../2.Output/07_ModulePreservation/datExprRef_color_GSE58095.Rdata"

#testset related
GSE9285_annotation_path<-"../1.Data/GSE9285/Annotation_GSE9285.txt"
GSE9285_ProcessedData_path<-"../1.Data/GSE9285/ProcessedData/"
testExpdata_path<-paste0("../1.Data/GSE9285/","Expdata.csv")
datExprTest_color_path<-"..//2.Output/07_ModulePreservation/datExprTest_color_GSE9285.Rdata"

geneid2modulecolor_path<-"../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv"


#######################################################################################
#Refset?̉??H


GSE58095_annotation<-read_csv(GSE58095_annotation_path) %>% 
  select(ID,Entrez_Gene_ID) %>% 
  left_join(geneid2modulecolor,by="Entrez_Gene_ID") %>% 
  filter(!is.na(Entrez_Gene_ID),!is.na(module))
GSE58095_probeid2kwithin<-read_csv(probeid_kwithin_path) %>% 
  select(ProbeID,kWithin)

lname<-load(refExpdata_path)
#Expdata_quontile
datExprRef_pre <- Expdata_quontile$E
keep<-rowSums(datExprRef_pre>log2(50))>=102
datExprRef_pre <- datExprRef_pre[keep,] %>% 
  as.data.frame() %>% 
  rownames_to_column("ILMNID") %>% 
  left_join(GSE58095_probeid2kwithin %>% left_join(GSE58095_annotation,by=c("ProbeID"="ID")),by=c("ILMNID"="ProbeID")) %>% 
  filter(!is.na(Entrez_Gene_ID)) %>% 
  #mutate(exprSignal_median=select(., matches("(?:\\d){3}")) %>% apply( 1, median))%>% #exprSignal_median?????GeneID?̏d???????????Ƃ?
  arrange(Entrez_Gene_ID %>% as.numeric(),desc(kWithin))%>% 
  distinct(Entrez_Gene_ID,.keep_all = T) %>% 
  column_to_rownames("Entrez_Gene_ID") 
datExprRef<-datExprRef_pre%>% 
  select(-c(ILMNID,kWithin,module))%>% t()
dim(datExprRef)
colorsRef <- datExprRef_pre$module 
#datExprRef %>% as.tibble()
#colnames(datExprRef) %>% length()
#datExprRef_pre %>% rownames() %>% length()
#(colnames(datExprRef) ==datExprRef_pre %>% rownames() ) %>% all()
save(list = c("datExprRef","colorsRef"),file = datExprRef_color_path)



####################################################################################################
#testset?̉??H

param <- c("REFSEQ", "SYMBOL","ENTREZID","ENSEMBL")        #?~?????A?m?e?[?V???????????w???i?L?q?\?ȃ??X?g?͈ȉ???keytypes(dbname)?????擾?\?j
keytypes(hgug4112a.db)                         #hoge???????ۂɎg?p?????L?[???\??
param_key <- keys(hgug4112a.db, keytype="PROBEID")#param1?Ŏw?肵?????̂??L?[?ɂ??Ă???

#?{??
out <- AnnotationDbi::select(hgug4112a.db, keys=param_key, keytype="PROBEID", columns=param)#?A?m?e?[?V???????????o???ʂ?out?Ɋi?[
write.table(out, GSE9285_annotation_path, sep="\t", append=F, quote=F, row.names=F)#out?̒??g???w?肵???t?@?C?????ŕۑ?

ProcessedData_files<-list.files(GSE9285_ProcessedData_path,full.names = T)
GSE9285_annotation<-read_tsv(GSE9285_annotation_path) %>% 
  select(-REFSEQ) %>% 
  distinct()

ProcessedData<-tibble()
for (file in ProcessedData_files){
  ProcessedData<-ProcessedData %>% 
    bind_rows(read_delim(file,delim = "\t",skip = ) %>% 
                mutate(sample=basename(file) %>% str_split("_") %>% .[[1]] %>% .[1]))
}

#RG<- ProcessedData %>% 
#  spread(key = sample,value = VALUE) %>% 
#  inner_join(annotation,by=c(`Reporter Identifier`="PROBEID")) %>% 
#  select(`Reporter Identifier`,"SYMBOL","ENTREZID","ENSEMBL",everything())%>% 
#  mutate(exprSignal_median=dplyr::select(.,starts_with("GSM")) %>%  apply( 1, median,na.rm=T)) %>% #select(exprSignal_median,everything())->a
#  arrange(ENTREZID,desc(exprSignal_median)) %>% 
#  distinct(ENTREZID,.keep_all = T) %>% #select(ENTREZID,exprSignal_median)
#  filter(!is.na(ENTREZID)) %>% 
#  column_to_rownames(var = "ENTREZID") %>% 
#  select(starts_with("GSM")) %>% 
#  as.matrix() %>% 
#  backgroundCorrect.matrix(method="normexp")

RG<- ProcessedData %>% 
  spread(key = sample,value = VALUE) %>% 
  column_to_rownames(var = "Reporter Identifier") %>% 
  as.matrix() %>% 
  backgroundCorrect.matrix(method="normexp")

RGnormalized<-normalizeBetweenArrays(RG, method="quantile") 
datExprTest_pre<-RGnormalized%>% as.data.frame() %>% 
  rownames_to_column("ReporterIdentifier") %>% 
  left_join(GSE9285_annotation %>% 
              left_join(geneid2modulecolor,by=c("ENTREZID"="Entrez_Gene_ID")),by=c("ReporterIdentifier"="PROBEID")) %>% 
  #select(`Reporter Identifier`,"SYMBOL","ENTREZID","ENSEMBL",everything())%>% 
  filter(!is.na(ENTREZID))  %>% 
  filter(!is.na(module)) %>% 
  mutate(exprSignal_median=dplyr::select(.,starts_with("GSM")) %>%  apply( 1, median,na.rm=T)) %>% #select(exprSignal_median,everything())->a
  arrange(ENTREZID,desc(exprSignal_median)) %>% 
  distinct(ENTREZID,.keep_all = T) 

  #select(starts_with("GSM")) %>% 
#write_csv(Expdata %>% as.data.frame() %>% rownames_to_column(var="ENTREZID"),testExpdata_path)

#datExprTest_pre<-RGfiltered%>% 
#  left_join(geneid2modulecolor,by=c("ENTREZID"="Entrez_Gene_ID")) %>% 
#  distinct(ENTREZID,.keep_all = T)
nrow(RGfiltered)
testExpdata

table(datExprTest_pre$module)
table(geneid2modulecolor$module)
#loss????geneid?̊???
((table(geneid2modulecolor$module)-table(datExprTest_pre$module))/table(geneid2modulecolor$module)) %>% round(2)

#datExprTest<-datExprTest_pre[,2:length(datExprTest_pre)-1] %>% t()
datExprTest<-datExprTest_pre%>% select(starts_with("GSM"))%>% t()
dim(datExprTest)
colnames(datExprTest) <- datExprTest_pre$ENTREZID
datExprTest[1:10, 1:10]
colorsTest <- datExprTest_pre$module 

save(list = c("datExprTest","colorsTest"),file = datExprTest_color_path)
####################################################################################################





007_testdataPreprocessing_fromDL_test.R





#library("hgug4112a.db", character.only=T)
library(tidyverse)
#library(WGCNA)
library(limma)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

target_file_path<-"../1.Data/GSE9285/test_Targets.txt"
target_processed_file_path<-"../1.Data/GSE9285/test_Targets_processed.txt"
raw_data_dir_path<-"../1.Data/GSE9285/RAWfiles/DecompressedData"
raw_data_fileformat<-".gpr"

targets<-read_tsv(target_file_path) %>% 
  mutate(FileName=str_c(Sample,raw_data_fileformat),
         )%>% # .$FileName
  rename(Name="Group")
write_tsv(targets,target_processed_file_path)

targets <- readTargets(target_processed_file_path,row.names="Name")
RG <- read.maimages(targets, path = raw_data_dir_path, source="genepix", green.only=FALSE)
#RG_test1 <- read.maimages(targets, source="agilent",green.only=TRUE)
#boxplot(RG$E,ylim=c(0,300))
#?plotMD(RG)
#Reference: limma user guide 6-3
RG.b <- backgroundCorrect(RG,method="minimum")
plotDensities(RG.b)
MA.p <-normalizeWithinArrays(RG.b)
plotDensities(MA.p)
MA.q <- normalizeBetweenArrays(RG.b, method="quantile")
plotDensities(MA.q, col="black")

boxplot(as.data.frame(MA.q) %>% select(starts_with("GSM")),main="Quantile normalization")#,
        #col= batch,names=1:nsampl)

library(GEOquery)

gse <- getGEO(filename="../1.Data/GSE9285/GSE9285_family.soft.gz",GSEMatrix = T)
show(gse)
gse@gsms[1]$GSM236379@dataTable@table->a

#gse_raw<-getGEOSuppFiles(filename="../1.Data/GSE9285/RAWfiles/GSE9285_RAW.tar")
gse_raw<-getGEO(filename="../1.Data/GSE9285/RAWfiles/GSE9285_RAW.tar")



########################################################################################
#GSE45485
series_matrix_path<-"../1.Data/GSE45485/GSE45485_series_matrix.txt.gz"
datExprTest_color_path<-"../2.Output/07_ModulePreservation/Testset_GSE45485/datExprTest_color_GSE45485.Rdata"
geneid2modulecolor_path<-"../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv"
geneid2modulecolor<-read_csv(geneid2modulecolor_path) %>% 
  select(Entrez_Gene_ID,module) %>% 
  distinct(Entrez_Gene_ID,module)

gse45485 <- getGEO(filename="../1.Data/GSE45485/GSE45485_family.soft.gz",GSEMatrix = T)
show(gse45485)
gse45485@gsms$GSM1104220@dataTable

testset_annotation<-gse45485@gpls$GPL6480@dataTable@table

testset_annotation<-testset_annotation %>% #colnames()
  select("ID","GENE") %>% 
  left_join(geneid2modulecolor,by=c("GENE"="Entrez_Gene_ID")) %>% 
  filter(!is.na(GENE) & !is.na(module))
  

#!Sample_data_processing
#Data for both channels were Lowess-normalized and then the log(2) ratio was taken
#read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45485/matrix/", comment="!")
SM<-read_tsv(series_matrix_path, comment="!")
SM%>% 
  gather(key="sample",value="value",-ID_REF) %>% 
  ggplot(aes(x=sample,y=value))+
  geom_boxplot()
ggsave(filename = "../2.Output/07_ModulePreservation/Testset_GSE45485/Boxplot_GSE45485_SM.png")

SM.q<-SM%>% 
  column_to_rownames(var = "ID_REF") %>% 
  as.matrix() %>% 
  normalizeBetweenArrays(method="quantile")

SM.q%>% as.data.frame() %>% 
  rownames_to_column("ID_REF") %>% 
  gather(key="sample",value="value",-ID_REF) %>% 
  ggplot(aes(x=sample,y=value))+
  geom_boxplot()
ggsave(filename = "../2.Output/07_ModulePreservation/Testset_GSE45485/Boxplot_GSE45485_SM_quantile.png")

datExprTest_pre<-SM.q %>% as.data.frame() %>%
  rownames_to_column(var="ID") %>% 
  left_join(testset_annotation,by="ID") %>% 
  filter(!is.na(GENE))  %>% 
  filter(!is.na(module)) %>% 
  mutate(exprSignal_median=dplyr::select(.,starts_with("GSM")) %>%  apply( 1, median,na.rm=T)) %>% #select(exprSignal_median,everything())->a
  arrange(GENE,desc(exprSignal_median)) %>% 
  distinct(GENE,.keep_all = T) 

datExprTest<-datExprTest_pre%>% 
  column_to_rownames(var="GENE") %>% 
  select(starts_with("GSM"))%>% t()
dim(datExprTest)
datExprTest[1:10, 1:10]
colorsTest <- datExprTest_pre$module 

save(list = c("datExprTest","colorsTest"),file = datExprTest_color_path)


########################################################################################
#Functionize
testset<-"GSE22356"

Entrez_Gene_ID_exp<-"ENTREZ_GENE_ID"
Probe_ID_exp<-"ID"

input_dir_path<-"../1.Data/GSE22356_PBMC/"
output_dir_path<-paste0("../2.Output/07_ModulePreservation/Testset_",testset,"_PBMC/")

series_matrix_path<-paste0(input_dir_path,testset,"_series_matrix.txt.gz")
gse_family_soft_path<-paste0(input_dir_path,testset,"_family.soft.gz")

Boxplot_SM_image_path <- paste0(output_dir_path,"Boxplot_",testset,"_SM.png")
Boxplot_SM.q_image_path <- paste0(output_dir_path,"Boxplot_",testset,"_SM_quantile.png")
datExprTest_color_path<-paste0(output_dir_path,"datExprTest_color_",testset,".Rdata")

geneid2modulecolor_path<-"../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv"

if (!file.exists(output_dir_path)){
  dir.create(output_dir_path)
}

geneid2modulecolor<-read_csv(geneid2modulecolor_path) %>% 
  select(Entrez_Gene_ID,module) %>% 
  distinct(Entrez_Gene_ID,module)

gse <- getGEO(filename=gse_family_soft_path)
gplid<-Meta(gse)$platform_id
gpl_data<-gse@gpls[gplid][[1]]@dataTable@table %>% as.tibble() %>% 
  select(eval(enquo(Probe_ID_exp)),eval(enquo(ENTREZ_GENE_ID_exp)))%>% 
  rename(Entrez_Gene_ID=eval(enquo(Entrez_Gene_ID_exp)),
         Probe_ID=eval(enquo(Probe_ID_exp))) %>% 
  unnest(Entrez_Gene_ID = strsplit(Entrez_Gene_ID, " /// ")) %>% #Probe?ɑ΂???GeneID???????o?^?????Ă????̂ł??̑Ή?" /// "?ŋ??؂????Ă?
  mutate(Entrez_Gene_ID=as.numeric(Entrez_Gene_ID))

testset_annotation<-gpl_data %>% 
  left_join(geneid2modulecolor,by="Entrez_Gene_ID") %>% 
  filter(!is.na(Entrez_Gene_ID) & !is.na(module))


#!Sample_data_processing
#Raw data was normalized using RMA and presented as non-log transformed signal intensity
SM<-read_tsv(series_matrix_path, comment="!")
SM%>% 
  gather(key="sample",value="value",-ID_REF) %>% 
  ggplot(aes(x=sample,y=value))+
  geom_boxplot()
ggsave(filename = Boxplot_SM_image_path)

SM.q<-SM%>% 
  column_to_rownames(var = "Probe_ID") %>% 
  as.matrix() %>% 
  normalizeBetweenArrays(method="quantile")

SM.q%>% as.data.frame() %>% 
  rownames_to_column("Probe_ID") %>% 
  gather(key="sample",value="value",-Probe_ID) %>% 
  ggplot(aes(x=sample,y=value))+
  geom_boxplot()
ggsave(filename =Boxplot_SM.q_image_path )

datExprTest_pre<-SM.q %>% as.data.frame() %>%
  rownames_to_column(var="Probe_ID") %>% 
  left_join(testset_annotation,by="Probe_ID") %>% 
  filter(!is.na(Entrez_Gene_ID))  %>% 
  filter(!is.na(module)) %>% 
  mutate(exprSignal_median=dplyr::select(.,starts_with("GSM")) %>%  apply( 1, median,na.rm=T)) %>% #select(exprSignal_median,everything())->a
  arrange(Entrez_Gene_ID,desc(exprSignal_median)) %>% 
  distinct(Entrez_Gene_ID,.keep_all = T) 

datExprTest<-datExprTest_pre%>% 
  column_to_rownames(var="Entrez_Gene_ID") %>% 
  select(starts_with("GSM"))%>% t()
dim(datExprTest)
datExprTest[1:10, 1:10]
colorsTest <- datExprTest_pre$module 

save(list = c("datExprTest","colorsTest"),file = datExprTest_color_path)






007_test_by_saitosans_data.R





library(tidyverse)
library(WGCNA)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190403_NAFLD_NASH_patient_WGCNA_GSE49541/analysis/WGCNA")

dat0_path<-"../signed_4/geneID_color_kWithin_normalizedData_logFC_qvalue.csv"
dat0 <- read_csv(dat0_path) %>% 
  arrange(Entrez_gene_ID %>% as.numeric(),desc(kWithin))%>% 
  distinct(Entrez_gene_ID,.keep_all = T)
table(duplicated(dat0$Entrez_gene_ID))
dat0[, 6:75]
datExprRef <- dat0[, 6:75]
datExprRef <- t(datExprRef)
dim(datExprRef)
colnames(datExprRef) <- dat0$Entrez_gene_ID
datExprRef[1:10, 1:10]
colorsRef <- dat0$moduleColor
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190606_NAFLD_NASH_patient_E_MEXP_3291/analysis/module_preservation")
saitohsandata_path<-paste0("Sample_normalized_data_logFC_qvalue.csv")

data<-read_csv(saitohsandata_path)%>% 
  mutate(exprSignal_median=select(., c(5:49)) %>% apply( 1, median)) %>% 
  arrange(Entrez_gene_ID %>% as.numeric(),desc(exprSignal_median))%>% 
  #select(exprSignal_median,everything())
  distinct(Entrez_gene_ID,.keep_all = T)
table(duplicated(data$Entrez_gene_ID))
data[130:140, 1:5]
datExprTest <- data[, 5:49]
datExprTest <- t(datExprTest)
colnames(datExprTest) <- data$Entrez_gene_ID
datExprTest[1:5, 130:140]
dim(datExprTest)

ref2test <- match(colnames(datExprRef), colnames(datExprTest))
table(is.finite(ref2test))
refset <- "GSE49541"
testset <- "E_MEXP_3291"
setLabels <- c(refset, testset)
multiExpr <- list(refset = list(data = datExprRef), testset = list(data = datExprTest))
multiColor <- list(refset = colorsRef)

system.time({
  mp <- modulePreservation(multiExpr, multiColor,
                           networkType = "signed",
                           referenceNetworks = 1,
                           nPermutations = 200,
                           randomSeed = 1,
                           quickCor = 0,
                           verbose = 3)
})





008_pathway_analysis_KEGG.R





library(tidyverse)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
c("#0086ab","#f79646","#9bbb59","#da6272","#777777","#bfbfbf")

####################################################################
#KEGGgraph test
library(KEGGgraph)
library(KEGG.db)
library("XML")
mapkKGML <- system.file("extdata/hsa04010.xml", package="KEGGgraph")
result <- xmlParse(file =mapkKGML)

####################################################################
library(clusterProfiler)
library(pathview)
library(tidyverse)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
output_dir<-"//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/2.Output/08_PathwayAnalysis/"

search_kegg_organism('hsa', by='kegg_code')#organismcode Ref: https://www.kegg.jp/kegg-bin/get_htext#D1
module<-read_csv(file = "../2.Output/02_signed_2_unmerged/logFC_ModuleColor.csv") 
write_keggEnrichResult<-function(KEresult,module_color){
  if (!dir.exists(paste0(output_dir,"/",module_color))){
    dir.create(paste0(output_dir,"/",module_color))
  }
  write_csv(KEresult,path = paste0(output_dir,module_color,"/","KEGG_PathwayEnrichment_",module_color,".csv"))
}
get_pathview<-function(logFC.namedlist, pathway.id.list,module_color, output_dir){
  current_dir<-getwd()#save first directory path
  if (!dir.exists(paste0(output_dir,"/",module_color))){
    dir.create(paste0(output_dir,"/",module_color))
  }
  if(length(pathway.id.list)==0){return(vector())}
  setwd(paste0(output_dir,"/",module_color))
  pv <- vector("list", length(pathway.id.list))
  names(pv)<-pathway.id.list
  for (i in seq_len(length(pathway.id.list))) {
    pv[[i]] <- pathview(gene.data = logFC.namedlist, pathway.id = pathway.id.list[i], 
                        species = "hsa", gene.idtype = "entrez",
                        limit = list(gene=2, cpd=1),
                        bins = list(gene=10,cpd=10), 
                        low =list(gene = "#0086ab", cpd = "blue"),#?????̐?#0086ab  #sns blue #34419A
                        mid = list(gene = "#FBFDC7", cpd ="gray"), #?x?[?W?? DDD9C3 sns yellow #FBFDC7
                        high = list(gene = "#E46C0A", cpd = "yellow"), #?????̃I?????W #E46C0A #sns red #AE0926
                        out.suffix = paste0(module_color))
  }
  setwd(current_dir)#return to first directory
  return(pv)
}

PA_summary.byModule<-module %>% 
  group_by(module_color) %>% 
  nest()%>% 
  mutate(kk=map(data,~{enrichKEGG(gene =.$Entrez_Gene_ID,organism = 'hsa')@result %>%  as.tibble()}))%>% 
  filter(module_color!="grey") %>% 
  mutate(logFC.namedlist=map(data,
                             function(df){
                               logFC<-df$logFC
                               names(logFC) <- df$Entrez_Gene_ID
                               return(logFC)
                            }),
         pathway.id.list=map(kk,~{filter(.,qvalue<=0.10) %>% .$ID}))

map2(.x = PA_summary.byModule$kk,.y = PA_summary.byModule$module_color,
     .f = write_keggEnrichResult)

PA_summary.byModule %>% 
  mutate(pathview=dplyr::select(.,logFC.namedlist,pathway.id.list,module_color) %>% mutate(output_dir=output_dir)%>% pmap(#.l = list(logFC.namedlist=logFC.namedlist,pathway.id.list=pathway.id.list, module_color=module_color,output_dir=rep(output_dir,nrow(tmp3))),
    .f = get_pathview))#->tmp4






008_pathway_analysis_Reactome.R





####################################################################
library(tidyverse)
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db) 
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
output_dir<-"//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/2.Output/08_PathwayAnalysis/"
c("#0086ab","#f79646","#9bbb59","#da6272","#777777","#bfbfbf")

data<-read_csv("../2.Output/02_signed_2_unmerged/logFC_ModuleColor.csv")
#clusterProfiler::bitr(data$Entrez_Gene_ID,fromType = "ENTREZID" ,toType = "ENSEMBL", OrgDb="org.Hs.eg.db")

tmp<-data %>% 
  nest(-module_color) %>% 
  #.[1:5,] %>% 
  mutate(enrchRsl=map(data,.f=~{enrichPathway(.$Entrez_Gene_ID,pvalueCutoff=0.05, readable=T)}),
         enrchRsl.filtered=map(enrchRsl,
                               .f = ~{as_tibble(.) %>% 
                                   filter(qvalue<=0.10) %>% 
                                   mutate(rank = dense_rank(qvalue))}))

map(tmp$reactomeEnrichment,.f=~{filter(as_tibble(.), qvalue<=0.10)})

tmp2<-tmp %>% 
  mutate(enrchRsl.filtered=map(enrchRsl,
                             .f = ~{as_tibble(.) %>% 
                                 filter(qvalue<=0.10) %>% 
                                 mutate(rank = dense_rank(qvalue))}))

#https://reactome.org/ContentService/exporter/diagram/R-HSA-8876198.png?flg=DENND2C,DENND4C,DENND6A,GDI2,RAB3GAP2,RAB3IP,RAB8B
download_ReactomeDiagram<-function(pathwayid,flag.list,output_dir,rank){
  query_url<-paste0("https://reactome.org/ContentService/exporter/diagram/",pathwayid,".png?quality=10&flg=",str_c(str_split(flag.list,"/")[[1]],collapse=","))
  filename<-paste0(output_dir,"Reactome_",formatC(rank,width=2,flag="0"),"_",pathwayid,".png")
  try(download.file(query_url,filename, mode = 'wb'))
}



map2(.x=tmp$module_color,.y=tmp$enrchRsl.filtered,
     .f = function(module_color,enrchRsl.filtered){
       output_dir<-paste0("../2.Output/08_PathwayAnalysis_/",module_color,"/")
       if (!dir.exists(output_dir)){
         dir.create(output_dir)
       }
       write_csv(enrchRsl.filtered,paste0(output_dir,"Reactome_PathwayEnrichment_",module_color,".csv"))
       enrchRsl_tmp<-enrchRsl.filtered %>% 
         filter(rank<=10) %>% 
         mutate(output_dir=output_dir) %>% 
         dplyr::select(ID,geneID,output_dir,rank) %>% 
         dplyr::rename(pathwayid="ID",
                flag.list="geneID")
       pmap(enrchRsl_tmp,download_ReactomeDiagram)
       
     })

####################################################################


a %>% as.tibble() %>% 
  dplyr::select(ID,geneID)

data %>% 
  filter(module_color=="turquoise") %>% 
  .$Entrez_Gene_ID %>% 
  enrichPathway(pvalueCutoff=0.05, readable=T) %>% 
  .@result #%>% #as_tibble() %>% 
  
barplot(a,showCategory=8)
dotplot(a,showCategory=8)
emapplot(a)
cnetplot(a, categorySize="pvalue")
a@result %>% as_tibble() %>% .$p.adjust

data %>% 
  filter(!is.na(Entrez_Gene_ID) & module_color!="grey") %>%
  dplyr::select(Entrez_Gene_ID,module_color) %>% 
  nest(-module_color) %>% 
  mutate(geneid.list=map(data,.f = ~{.$Entrez_Gene_ID %>% as.character()})) %>% 
  .$geneid.list->test
names(test)<-data %>% 
  dplyr::select(Entrez_Gene_ID,module_color) %>% 
  nest(-module_color) %>% 
  mutate(geneid.list=map(data,.f = ~{.$Entrez_Gene_ID %>% as.character()})) %>% 
  filter(module_color!="grey") %>% 
  .$module_color
xx.formula <- compareCluster(Entrez_Gene_ID~module_color, data=data, fun = 'groupGO', OrgDb='org.Hs.eg.db')



significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor
compareModule_res <- compareCluster(test[names(test) %in% significant_module], fun = 'enrichPathway')
dotplot(compareModule_res)

test[names(test) %in% significant_module]

####################################################################

library(DOSE)
data("geneList")
yellow<-read_csv("../2.Output/08_PathwayAnalysis_/yellow/Reactome_PathwayEnrichment_yellow.csv")

viewPathway(yellow[2,]$Description )+
  geom_edge_link()
vp<-viewPathway("E2F mediated regulation of DNA replication", readable=TRUE, foldChange=geneList)
vp %>% class()
vp$data


vp+
  geom_edge_link()+
  geom_node_point(aes(color=color %>% as.numeric()),size=5)+
  geom_node_text(aes(label=name))
a<-get_edges(vp)
E(vp)

as_long_data_frame(vp)
as_tbl_graph(vp$data) %>% 
  ggraph()+
  geom_edge_link()
reactome.db
library(GOSemSim)
library(DOSE)
library(igraph)
library(tidygraph)
library(ggraph)

viewPathway2 <- function(pathName,
                        organism="human",
                        readable=TRUE,
                        foldChange=NULL,
                        keyType = "ENTREZID",
                        layout = "kk", ...){
  
  ## call pathways via imported from graphite has the following issue:
  ##
  ## Error: processing vignette 'ReactomePA.Rnw' failed with diagnostics:
  ## no item called "package:graphite" on the search list
  ## Execution halted
  ##
  
  pkg <- "graphite"
  require(pkg, character.only=TRUE)
  
  # convertion to the names that graphite::pathways understands
  org2org <- list(arabidopsis="athaliana",
                  bovine="btaurus",
                  canine="cfamiliaris",
                  chicken="ggallus",
                  ecolik12="ecoli",
                  fly="dmelanogaster",
                  human="hsapiens",
                  mouse="mmusculus",
                  pig="sscrofa",
                  rat="rnorvegicus",
                  celegans="celegans",
                  xenopus="xlaevis",
                  yeast="scerevisiae",
                  zebrafish="drerio")
  
  if(!(organism %in% names(org2org))){
    cat(paste(c("the list of supported organisms:",names(org2org)), collapse='\n'))
    stop(sprintf("organism %s is not supported", organism))
  }
  pathways <- eval(parse(text="pathways"))
  p <- pathways(org2org[[organism]], 'reactome')[[pathName]]
  
  if (readable) {
    p <- convertIdentifiers(p, "symbol")
    if (!is.null(foldChange)) {
      stopifnot(!any(duplicated(names(foldChange)))) # can't have two value for one gene
      OrgDb <- getDb(organism)
      names(foldChange) <- DOSE::EXTID2NAME(OrgDb, names(foldChange), keyType)
      
    }
  } else {
    if (!is.null(foldChange)) {
      p <- convertIdentifiers(p, "entrez")
    }
  }
  
  g <- pathwayGraph(p)
  gg <- igraph.from.graphNEL(g)
  gg <- as.undirected(gg)
  gg <- setting.graph.attributes(gg)
  V(gg)$name <- sub("[^:]+:", "", V(gg)$name)
  
  if (!is.null(foldChange)) {
    ## gg <- scaleNodeColor(gg, foldChange)
    fc <- foldChange[V(gg)$name]
    V(gg)$color <- fc
    ## palette <- enrichplot:::fc_palette(fc)
    
  }
  ## netplot(gg, foldChange=foldChange, ...)
  ggraph(as_tbl_graph(gg), layout=layout) +
    geom_edge_link(alpha=.8, colour='darkgrey') +
    geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~size)) +
    scale_color_continuous(low="red", high="blue", name = "fold change", na.value = "#E5C494") +
    geom_node_text(aes_(label=~name), repel=TRUE) +
    ## scale_color_gradientn(name = "fold change", colors=palette, na.value = "#E5C494") +
    scale_size(guide = FALSE) + theme_void()
}
getDb <- function(db.name.or.index = NULL, debug = FALSE) {
  if (debug) {
    print(paste("IN:", match.call()[[1]]))
  }
  
  .separator <- .Platform$file.sep # Platform dependent path separator.
  
  # LOAD DATABASE INFO  #######################################################
  
  # Get package path.
  if (require(strvalidator)) {
    packagePath <- path.package("strvalidator", quiet = FALSE)
  } else {
    warning("Package path for strvalidator not found!")
    return(NULL)
  }
  subFolder <- "extdata"
  fileName <- "database.txt"
  
  # Create complete file path.
  filePath <- paste(packagePath, subFolder, fileName, sep = .separator)
  
  .db <- read.delim(
    file = filePath, header = TRUE, sep = "\t", quote = "\"",
    dec = ".", fill = TRUE, stringsAsFactors = FALSE
  )
  
  # Available databases. Must match else if construct.
  databases <- unique(.db$Database)
  
  # Check if NULL
  if (is.null(db.name.or.index)) {
    db <- databases
    
    # String provided.
  } else {
    
    # Check if number or string.
    if (is.numeric(db.name.or.index)) {
      
      # Set index to number.
      index <- db.name.or.index
    } else {
      
      # Find matching database index (case insensitive)
      index <- match(toupper(db.name.or.index), toupper(databases))
    }
    
    # No matching database.
    if (is.na(index)) {
      db <- NA
      
      # Assign matching database information.
    } else {
      db <- .db[.db$Database == databases[index], ]
    }
  }
  
  return(db)
}
vp2<-viewPathway2("E2F mediated regulation of DNA replication",organism="human", readable=TRUE, foldChange=geneList)

ggraph(as_tbl_graph(gg), layout="kk") +
  geom_edge_link(alpha=.8, colour='darkgrey') +
  geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~size)) +
  scale_color_continuous(low="red", high="blue", name = "fold change", na.value = "#E5C494") +
  geom_node_text(aes_(label=~name), repel=TRUE) +
  ## scale_color_gradientn(name = "fold change", colors=palette, na.value = "#E5C494") +
  scale_size(guide = FALSE) + theme_graph()



##' @importFrom igraph V
##' @importFrom igraph V<-
##' @importFrom igraph E
##' @importFrom igraph E<-
setting.graph.attributes <- function(g, node.size=8,
                                     node.color="#B3B3B3",
                                     edege.width=2,
                                     edege.color="#8DA0CB") {
  V(g)$size <- node.size
  V(g)$color <- node.color
  V(g)$label <- V(g)$name
  
  E(g)$width <- edege.width
  E(g)$color <- edege.color
  
  return(g)
}






008_pathway_analysis_test.R





library(tidyverse)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
c("#0086ab","#f79646","#9bbb59","#da6272","#777777","#bfbfbf")

####################################################################
#KEGGgraph test
library(KEGGgraph)
library(KEGG.db)
library("XML")
mapkKGML <- system.file("extdata/hsa04010.xml", package="KEGGgraph")
result <- xmlParse(file =mapkKGML)


####################################################################
#BiocManager::install(c("clusterProfiler","pathview"))
library(clusterProfiler)
library(pathview)
search_kegg_organism('hsa', by='kegg_code')#organismcode Ref: https://www.kegg.jp/kegg-bin/get_htext#D1


module<-read_csv(file = "../2.Output/02_signed_2_unmerged/logFC_ModuleColor.csv") %>% 
  filter(module_color=="yellow" & !is.na(Entrez_Gene_ID))
module.gene<-module%>% 
  .$Entrez_Gene_ID
module.logFC<-module$logFC
as.vector(module$logFC,)

kk <- enrichKEGG(gene = module.gene, organism = 'hsa')
head(kk, n=10)


module.gene<-module$Entrez_Gene_ID
names(module.logFC) <- module$Entrez_Gene_ID
pv<-pathview(gene.data = module.logFC, 
         pathway.id = "hsa03015", 
         species = "hsa", 
         limit = list(gene=2, cpd=1),
         bins = list(gene=10,cpd=10), 
         low =list(gene = "#0086ab", cpd = "blue"),#?????̐?#0086ab  #sns blue #34419A
         mid = list(gene = "#FBFDC7", cpd ="gray"), #?x?[?W?? DDD9C3 sns yellow #FBFDC7
         high = list(gene = "#E46C0A", cpd = "yellow"),#?????̃I?????W #E46C0A #sns red #AE0926
         out.suffix = "test_suffix")
gene_in_pathway<-kk@result %>% filter(ID=="hsa03015") %>% .$geneID %>% 
  str_split("/") %>% .[[1]]
module %>% 
  filter(Entrez_Gene_ID %in% gene_in_pathway) %>% 
  dplyr::select(Entrez_Gene_ID, Symbol, logFC)

kk@result 

####################################################################
library(clusterProfiler)
library(pathview)
search_kegg_organism('hsa', by='kegg_code')#organismcode Ref: https://www.kegg.jp/kegg-bin/get_htext#D1


module<-read_csv(file = "../2.Output/02_signed_2_unmerged/logFC_ModuleColor.csv") 

tmp<-module %>% 
  group_by(module_color) %>% 
  nest()

tmp2<-tmp %>% 
  mutate(kk=map(data,~{enrichKEGG(gene =.$Entrez_Gene_ID,organism = 'hsa')@result %>%  as.tibble()}))
 
map2(.x = tmp2$kk,.y = tmp2$module_color,
     .f = ~{write_csv(.x,path = paste0("../2.Output/08_PathwayAnalysis/","KEGG_PathwayEnrichment_",.y,".csv"))})

tmp3<-tmp2 %>% 
  filter(module_color!="grey") %>% 
  mutate(logFC.namedlist=map(data,
                             function(df){
                               logFC<-df$logFC
                               names(logFC) <- df$Entrez_Gene_ID
                               return(logFC)}),
         pathway.id.list=map(kk,~{filter(.,qvalue<=0.10) %>% .$ID}))



get_pathview<-function(logFC.namedlist, pathway.id.list,module_color, output_dir){
  current_dir<-getwd()#save first directory path
  if (!dir.exists(paste0(output_dir,"/",module_color))){
    dir.create(paste0(output_dir,"/",module_color))
  }
  if(length(pathway.id.list)==0){return(vector())}
  setwd(paste0(output_dir,"/",module_color))
  pv <- vector("list", length(pathway.id.list))
  names(pv)<-pathway.id.list
  for (i in seq_len(length(pathway.id.list))) {
    pv[[i]] <- pathview(gene.data = logFC.namedlist, pathway.id = pathway.id.list[i], 
                        species = "hsa", gene.idtype = "entrez",
                        limit = list(gene=2, cpd=1),
                        bins = list(gene=10,cpd=10), 
                        low =list(gene = "#0086ab", cpd = "blue"),#?????̐?#0086ab  #sns blue #34419A
                        mid = list(gene = "#FBFDC7", cpd ="gray"), #?x?[?W?? DDD9C3 sns yellow #FBFDC7
                        high = list(gene = "#E46C0A", cpd = "yellow"), #?????̃I?????W #E46C0A #sns red #AE0926
                        out.suffix = paste0(module_color))
  }
  setwd(current_dir)#return to first directory
  return(pv)
}
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
output_dir<-"//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/2.Output/08_PathwayAnalysis/"
tmp3[39:nrow(tmp3),] %>% #filter(module_color!="sienna3") %>% 
  mutate(pathview=dplyr::select(.,logFC.namedlist,pathway.id.list,module_color) %>% mutate(output_dir=output_dir)%>% pmap(#.l = list(logFC.namedlist=logFC.namedlist,pathway.id.list=pathway.id.list, module_color=module_color,output_dir=rep(output_dir,nrow(tmp3))),
                       .f = get_pathview))->tmp4

tmp3$module_color
tmp3[39:nrow(tmp3),]

tmp3 %>% filter(module_color=="midnightblue") %>% .$pathway.id.list

tmp3 %>% 
  mutate(pathview=map(.x = tmp3,.f = function(df){}))

setwd("../2.Output/08_PathwayAnalysis/")
pmap(.x = tmp3$logFC.namedlist,.y = tmp3$pathway.id.list,
     .f = ~{})


a$pathway.id.list[[1]]

a$kk %>% .[[1]] %>% as.tibble()

tmp %>% 
  mutate(kk=map(data,myfunc))
tmp %>% 
  mutate(geneid=map(data,.f = ~{.$Entrez_Gene_ID }))

tmp$data[[1]]$Entrez_Gene_ID


myfunc<-function(df){enrichKEGG(gene =df$Entrez_Gene_ID,organism = 'hsa')}
myfunc(tmp$data[[1]])






008_reactomefi_test.R





library(reactomefi)
library(tidyverse)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
genes<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  filter(module=="yellow") %>% 
  .$Entrez_Gene_ID %>% as.character()
genes
network <- reactomefi::ReactomeFINetwork("2012", genes)
annotate(network, "Pathway")
plot(network)





008_retrieve_ReactomePathwayNetwork.R





library(tidyverse)
library(GOSemSim)
library(DOSE)
library(graphite)
library(org.Mm.eg.db)
library(igraph)
library(tidygraph)
library(ggraph)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

#yellow<-read_csv("../2.Output/08_PathwayAnalysis_/yellow/Reactome_PathwayEnrichment_yellow.csv")
#pathName<-yellow[2,]$Description
#foldChange<-NULL
data("geneList")

foldChange<-geneList
pathName<-"E2F mediated regulation of DNA replication"
organism<-"human"
org2org <- list(arabidopsis="athaliana",
                bovine="btaurus",
                canine="cfamiliaris",
                chicken="ggallus",
                ecolik12="ecoli",
                fly="dmelanogaster",
                human="hsapiens",
                mouse="mmusculus",
                pig="sscrofa",
                rat="rnorvegicus",
                celegans="celegans",
                xenopus="xlaevis",
                yeast="scerevisiae",
                zebrafish="drerio")

p<-graphite::pathways(org2org[[organism]], 'reactome')[[pathName]]
p <- convertIdentifiers(p, "entrez")

#undirected?ɂ????????Ă邯??directed?̂܂܂̕????????ʑ???
#directed??undirected?????????O???t?̕`?悪?ł????΂悢?
gg <- pathwayGraph(p) %>% 
  igraph.from.graphNEL() %>% 
  as.undirected() %>% 
  setting.graph.attributes() %>% 
  as_tbl_graph() %N>% 
  mutate(name=str_sub(name,start = 10,end = -1),
         label=mapIds(org.Hs.eg.db,name, 'SYMBOL',"ENTREZID") %>% as.vector())

if (!is.null(foldChange)) {
  #names(foldChange) <- mapIds(org.Hs.eg.db, names(foldChange), 'SYMBOL',"ENTREZID") #%>% as.vector()
  ## gg <- scaleNodeColor(gg, foldChange)
  fc <- foldChange[V(gg)$name]
  V(gg)$color <- fc
  ## palette <- enrichplot:::fc_palette(fc)
  
}

ggraph(gg, layout="kk") +
  geom_edge_link(alpha=.8, colour='darkgrey') +
  geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~size)) +
  scale_color_continuous(low="red", high="blue", name = "fold change", na.value = "#E5C494") +
  geom_node_text(aes_(label=~label), repel=TRUE) +
  ## scale_color_gradientn(name = "fold change", colors=palette, na.value = "#E5C494") +
  scale_size(guide = FALSE) + 
  theme_graph()



#######################################################################
#reference::https://github.com/YuLab-SMU/ReactomePA/blob/master/R/viewPathway.R
viewPathway2 <- function(pathName,
                         organism="human",
                         readable=TRUE,
                         foldChange=NULL,
                         keyType = "ENTREZID",
                         layout = "kk", ...){
  
  ## call pathways via imported from graphite has the following issue:
  ##
  ## Error: processing vignette 'ReactomePA.Rnw' failed with diagnostics:
  ## no item called "package:graphite" on the search list
  ## Execution halted
  ##
  
  pkg <- "graphite"
  require(pkg, character.only=TRUE)
  
  # convertion to the names that graphite::pathways understands
  org2org <- list(arabidopsis="athaliana",
                  bovine="btaurus",
                  canine="cfamiliaris",
                  chicken="ggallus",
                  ecolik12="ecoli",
                  fly="dmelanogaster",
                  human="hsapiens",
                  mouse="mmusculus",
                  pig="sscrofa",
                  rat="rnorvegicus",
                  celegans="celegans",
                  xenopus="xlaevis",
                  yeast="scerevisiae",
                  zebrafish="drerio")
  
  if(!(organism %in% names(org2org))){
    cat(paste(c("the list of supported organisms:",names(org2org)), collapse='\n'))
    stop(sprintf("organism %s is not supported", organism))
  }
  pathways <- eval(parse(text="pathways"))
  p <- pathways(org2org[[organism]], 'reactome')[[pathName]]
  
  if (readable) {
    p <- convertIdentifiers(p, "symbol")
    if (!is.null(foldChange)) {
      stopifnot(!any(duplicated(names(foldChange)))) # can't have two value for one gene
      OrgDb <- getDb(organism)
      names(foldChange) <- DOSE::EXTID2NAME(OrgDb, names(foldChange), keyType)
      
    }
  } else {
    if (!is.null(foldChange)) {
      p <- convertIdentifiers(p, "entrez")
    }
  }
  
  g <- pathwayGraph(p)
  gg <- igraph.from.graphNEL(g)
  gg <- as.undirected(gg)
  gg <- setting.graph.attributes(gg)
  V(gg)$name <- sub("[^:]+:", "", V(gg)$name)
  
  if (!is.null(foldChange)) {
    ## gg <- scaleNodeColor(gg, foldChange)
    fc <- foldChange[V(gg)$name]
    V(gg)$color <- fc
    ## palette <- enrichplot:::fc_palette(fc)
    
  }
  ## netplot(gg, foldChange=foldChange, ...)
  ggraph(as_tbl_graph(gg), layout=layout) +
    geom_edge_link(alpha=.8, colour='darkgrey') +
    geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~size)) +
    scale_color_continuous(low="red", high="blue", name = "fold change", na.value = "#E5C494") +
    geom_node_text(aes_(label=~name), repel=TRUE) +
    ## scale_color_gradientn(name = "fold change", colors=palette, na.value = "#E5C494") +
    scale_size(guide = FALSE) + theme_void()
}
getDb <- function(db.name.or.index = NULL, debug = FALSE) {
  if (debug) {
    print(paste("IN:", match.call()[[1]]))
  }
  
  .separator <- .Platform$file.sep # Platform dependent path separator.
  
  # LOAD DATABASE INFO  #######################################################
  
  # Get package path.
  if (require(strvalidator)) {
    packagePath <- path.package("strvalidator", quiet = FALSE)
  } else {
    warning("Package path for strvalidator not found!")
    return(NULL)
  }
  subFolder <- "extdata"
  fileName <- "database.txt"
  
  # Create complete file path.
  filePath <- paste(packagePath, subFolder, fileName, sep = .separator)
  
  .db <- read.delim(
    file = filePath, header = TRUE, sep = "\t", quote = "\"",
    dec = ".", fill = TRUE, stringsAsFactors = FALSE
  )
  
  # Available databases. Must match else if construct.
  databases <- unique(.db$Database)
  
  # Check if NULL
  if (is.null(db.name.or.index)) {
    db <- databases
    
    # String provided.
  } else {
    
    # Check if number or string.
    if (is.numeric(db.name.or.index)) {
      
      # Set index to number.
      index <- db.name.or.index
    } else {
      
      # Find matching database index (case insensitive)
      index <- match(toupper(db.name.or.index), toupper(databases))
    }
    
    # No matching database.
    if (is.na(index)) {
      db <- NA
      
      # Assign matching database information.
    } else {
      db <- .db[.db$Database == databases[index], ]
    }
  }
  
  return(db)
}
vp2<-viewPathway2("E2F mediated regulation of DNA replication",organism="human", readable=TRUE, foldChange=geneList)

ggraph(as_tbl_graph(gg), layout="kk") +
  geom_edge_link(alpha=.8, colour='darkgrey') +
  geom_node_point(aes_(color=~as.numeric(as.character(color)), size=~size)) +
  scale_color_continuous(low="red", high="blue", name = "fold change", na.value = "#E5C494") +
  geom_node_text(aes_(label=~name), repel=TRUE) +
  ## scale_color_gradientn(name = "fold change", colors=palette, na.value = "#E5C494") +
  scale_size(guide = FALSE) + theme_graph()







009_CalculateCorrelationILD.R





library(WGCNA)
library(caret)
library(GGally)
library(ggraph)
library(tidygraph)
library(tidyverse)
library(polycor)

correlation_method<-"spearman"
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
lnames<-load("../2.Output/02_signed_2_unmerged/signed_2/networkConstruction_StepByStep_signed_2.Rdata")
sample_annotation_file_path<-"../1.Data/GSE58095/GSE58095_all_sample_annotations.txt"
output_file_base_path<-"../2.Output/09_ILDCorrelation/"

sample_annotation <- read_tsv(sample_annotation_file_path,skip = 4,n_max = 102) %>% 
  rename_all(function(x) x %>% str_replace_all(" ","_")%>% 
               str_replace_all(":","_") %>% 
               str_replace_all("\\(","") %>% 
               str_replace_all("\\)","") %>% 
               str_replace_all("=","is"))%>% 
  dplyr::select(-c(title,Characteristics__Patient,Characteristics__GSE47162_Sample_name,Characteristics__Time_point,
                   molecule,label,description,platform)) %>% 
  mutate(Characteristics__Group=Characteristics__Group %>% replace_na("SSc"),
         characteristics__FVC_less_than_70=characteristics__FVC_less_than_70 %>% replace_na("No ILD"),
         Characteristics__On_immunosurressive_agents_1isyes = Characteristics__On_immunosurressive_agents_1isyes %>% replace_na(0))


tmp <- dummyVars(~., data=sample_annotation %>% dplyr::select(-"Sample_name"))
sample_annotation.dummy <- as.data.frame(predict(tmp, sample_annotation)) %>% 
  mutate(Sample_name=sample_annotation$Sample_name) %>%
  as_tibble()%>% 
  rename_at(vars(starts_with("Characteristics__")),function(x) str_sub(x,18,-1)) %>% 
  dplyr::select(Sample_name,everything())


MEs_withTrait<-MEs_signed %>% as.data.frame() %>% 
  rownames_to_column("Sample_name") %>% 
  as_tibble() %>% 
  mutate(Sample_name=str_c("SAMPLE",Sample_name)) %>% 
  left_join(sample_annotation.dummy,by="Sample_name")

ME_Trait_cor<-MEs_withTrait%>% 
  dplyr::select(-c("Sample_name","GroupSSc","GroupControl") )%>% stats::cor( use = "complete.obs",method = correlation_method) %>% 
  as.data.frame() %>% 
  rownames_to_column("vs1") %>% 
  as_tibble() %>% 
  gather(key="vs2",value="correlation",-vs1) #

MEs_withTrait%>% 
  dplyr::select(-c("Sample_name","GroupSSc","GroupControl") )%>% stats::cor( use = "complete.obs",method = correlation_method) %>% 
  as.data.frame() %>% 
  rownames_to_column("vs1") %>% 
  as_tibble() %>% 
  select(vs1,FVC_less_than_70ILD,On_immunosurressive_agents_1isyes) %>% 
  filter(grepl("^ME", vs1)) %>% 
  gather(key="feature",value="correlation",-vs1)%>% 
  ggplot(aes(x=feature,y=vs1,fill=correlation))+
  geom_tile()+
  scale_fill_gradient2() + 
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  xlab("")+
  ylab("")


MEs_withTrait %>% 
  group_by(FVC_less_than_70ILD) %>% 
  nest()

MEs_withTrait %>% 
  gather(key=feature,value=value,-c("FVC_less_than_70ILD","Sample_name")) %>% 
  group_by(feature) %>% 
  filter(grepl("^ME", feature)) %>% 
  nest() %>% 
  mutate(ILD_cor=map(data,
                     function(df){
                       tmp<-cor.test(df$FVC_less_than_70ILD,df$value)
                       return(list(tmp$estimate,tmp$p.value))}),
    cor=ILD_cor[[1]][[1]],
    p.value=ILD_cor[[1]][[2]]) %>% 
  filter(p.value<0.05)

#a %>% 
#  unnest(ILD_cor) %>% 
#  group_by(feature) %>% 
#  mutate(col=seq_along(feature)) %>% 
#  spread(key=col,value=ILD_cor)

#######################################################################
#ttest ILD vs Not ILD
MEs_withTrait %>% 
  select(-Sample_name) %>% 
  gather(key=feature,value=value,-FVC_less_than_70ILD) %>% 
  filter(grepl("^ME", feature)) %>% 
  group_by(feature) %>% 
  nest() %>% 
  mutate(t.test=map(data,
                    function(df){
                      ILD<-df %>% filter(FVC_less_than_70ILD == 1) %>% .$value
                      NotILD<-df %>% filter(FVC_less_than_70ILD == 0) %>% .$value
                      t.result<-t.test(ILD, NotILD, var.equal=T)
                      return(c(t.result$p.value,t.result$estimate["mean of x"],t.result$estimate["mean of y"]))
                    }),
         t.test.pvalue=t.test[[1]][[1]],
         t.test.est.ild=t.test[[1]][[2]],
         t.test.est.notild=t.test[[1]][[3]]) %>% 
  filter(t.test.pvalue<=0.05) %>% 
  arrange(t.test.pvalue) %>% {
  ggplot(.,aes(x=feature %>% str_sub(3,-1) %>% reorder(desc(t.test.pvalue)),y=-log10(t.test.pvalue)))+
  geom_bar(aes(fill=feature %>% str_sub(3,-1)),stat = "identity")+
  scale_fill_identity()+
  scale_y_continuous(expand = c(0, 0),limits =c(0,max(-log10(.$t.test.pvalue))*1.1))+
  coord_flip()+
  theme_classic()+
  theme(plot.margin= unit(c(0,0,0,0), "lines"),
        panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.title.y=element_blank())+
      labs(y="-log(P-value)")
}
savefile<-paste0(output_file_base_path,"ILDvsNotILD_Ttest_pvalue_barplot",".png")
ggsave(file = savefile, dpi = 320, width = 100, height = 150, units = "mm",  bg = "transparent")


##############################################################################################
#ttest on Immunosuppressive vs not
MEs_withTrait %>% 
  select(-Sample_name) %>% 
  gather(key=feature,value=value,-On_immunosurressive_agents_1isyes) %>% 
  filter(grepl("^ME", feature)) %>% 
  group_by(feature) %>% 
  nest() %>% 
  mutate(t.test=map(data,
                    function(df){
                      ILD<-df %>% filter(On_immunosurressive_agents_1isyes == 1) %>% .$value
                      NotILD<-df %>% filter(On_immunosurressive_agents_1isyes == 0) %>% .$value
                      t.result<-t.test(ILD, NotILD, var.equal=T)
                      return(c(t.result$p.value,t.result$estimate["mean of x"],t.result$estimate["mean of y"]))
                    }),
         t.test.pvalue=t.test[[1]][[1]],
         t.test.est.ild=t.test[[1]][[2]],
         t.test.est.notild=t.test[[1]][[3]]) %>% 
  filter(t.test.pvalue<=0.05) %>% 
  arrange(t.test.pvalue) %>% {
    ggplot(.,aes(x=feature %>% str_sub(3,-1) %>% reorder(desc(t.test.pvalue)),y=-log10(t.test.pvalue)))+
      geom_bar(aes(fill=feature %>% str_sub(3,-1)),stat = "identity")+
      scale_fill_identity()+
      scale_y_continuous(expand = c(0, 0),limits =c(0,max(-log10(.$t.test.pvalue))*1.1))+
      coord_flip()+
      theme_classic()+
      theme(plot.margin= unit(c(0,0,0,0), "lines"),
            panel.background = element_blank(),#fill = "transparent"), # bg of the panel
            panel.border = element_blank(),
            plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
            panel.grid.major = element_blank(), # get rid of major grid
            panel.grid.minor = element_blank(), # get rid of minor grid
            legend.background = element_rect(fill = "transparent"), # get rid of legend bg
            legend.box.background = element_rect(fill = "transparent"),
            axis.title.y=element_blank())+
      labs(y="-log(P-value)")
  }
savefile<-paste0(output_file_base_path,"ImmunosuppressivevsNot_Ttest_pvalue_barplot",".png")
ggsave(file = savefile, dpi = 320, width = 100, height = 150, units = "mm",  bg = "transparent")

##############################################################################################
#correlation
MEs_withTrait %>% 
  select(-Sample_name) %>% 
  gather(key=feature,value=value,-c("FVC_less_than_70ILD","On_immunosurressive_agents_1isyes")) %>% 
  filter(grepl("^ME", feature)) %>% 
  group_by(feature) %>% 
  nest() %>% 
  mutate(ILD.cor.result=map(data,
                    function(df){
                      polychor.result<-polychor(df[["value"]], df[["FVC_less_than_70ILD"]],std.err = T)
                      return(polychor.result)
                    }),
         ILD.cor=ILD.cor.result[[1]]$rho,
         ILD.pvalue=pchisq(ILD.cor.result[[1]]$chisq, ILD.cor.result[[1]]$df, lower.tail = FALSE),
         ImmunSup.cor.result=map(data,
                            function(df){
                              polychor.result<-polychor(df[["value"]], df[["On_immunosurressive_agents_1isyes"]],std.err = T)
                              return(polychor.result)
                            }),
         ImmunSup.cor=ImmunSup.cor.result[[1]]$rho,
         ImmunSup.pvalue=pchisq(ImmunSup.cor.result[[1]]$chisq, ImmunSup.cor.result[[1]]$df, lower.tail = FALSE)
         ) %>% 
  filter(ILD.pvalue<=0.05 | ImmunSup.pvalue<=0.05) %>% 
  arrange(ILD.pvalue) %>% {
    ggplot(.,aes(x=feature %>% str_sub(3,-1) %>% reorder(desc(t.test.pvalue)),y=-log10(t.test.pvalue)))+
      geom_point(aes(fill=feature %>% str_sub(3,-1)),stat = "identity")+
      scale_fill_identity()+
      scale_y_continuous(expand = c(0, 0),limits =c(0,max(-log10(.$t.test.pvalue))*1.1))+
      coord_flip()+
      theme_classic()+
      theme(plot.margin= unit(c(0,0,0,0), "lines"),
            panel.background = element_blank(),#fill = "transparent"), # bg of the panel
            panel.border = element_blank(),
            plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
            panel.grid.major = element_blank(), # get rid of major grid
            panel.grid.minor = element_blank(), # get rid of minor grid
            legend.background = element_rect(fill = "transparent"), # get rid of legend bg
            legend.box.background = element_rect(fill = "transparent"),
            axis.title.y=element_blank())+
      labs(y="-log(P-value)")
  }

MEs_withTrait %>% 
  select(-Sample_name) %>% 
  gather(key=feature,value=value,-c("FVC_less_than_70ILD","On_immunosurressive_agents_1isyes")) %>% 
  filter(grepl("^ME", feature)) %>% 
  group_by(feature) %>% 
  nest() %>% 
  mutate(ILD.cor.result=map(data,
                            function(df){
                              cor.result<-cor.test(df[["value"]], df[["FVC_less_than_70ILD"]],method = "spearman")
                              return(cor.result)
                            }),
         ILD.cor=ILD.cor.result[[1]]$estimate,
         ILD.pvalue=ILD.cor.result[[1]]$p.value,
         ImmunSup.cor.result=map(data,
                                 function(df){
                                   cor.result<-cor.test(df[["value"]], df[["On_immunosurressive_agents_1isyes"]],method = "spearman")
                                   return(cor.result)
                                 }),
         ImmunSup.cor=ImmunSup.cor.result[[1]]$estimate,
         ImmunSup.pvalue=ImmunSup.cor.result[[1]]$p.value	
         ) %>% 
  filter(ILD.pvalue<=0.05 | ImmunSup.pvalue<=0.05) %>% 
  arrange(ILD.pvalue) %>%
  select(-c("data","ILD.cor.result","ImmunSup.cor.result")) %>% 
  gather(key=cor,value=value,-feature) %>% 
  separate(col = cor,into = c("ILD_Immune","cor_pvalue")) %>% 
  spread(key=cor_pvalue,value=value) %>% 
  ggplot(aes(x=feature %>% str_sub(3,-1) %>% reorder(cor)))+
  geom_point(aes(fill=feature %>% str_sub(3,-1),size=-log10(pvalue),color=cor,y=ILD_Immune))+
  scale_fill_identity()+
  scale_color_gradient2()+
  scale_size_continuous(range=c(5,20))+
  coord_flip()+
  theme_classic()+
  theme(plot.margin= unit(c(0,0,0,0), "lines"),
        panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.title.y=element_blank())+
  labs(x="",y="")

savefile<-paste0(output_file_base_path,"Correlation_dotplot",".png")
ggsave(file = savefile, dpi = 320, width = 100, height = 150, units = "mm",  bg = "transparent")

###############################################################################################
tmp<-MEs_withTrait%>% 
  dplyr::select(-c("Sample_name","GroupSSc","GroupControl") )%>% stats::cor( use = "complete.obs",method = correlation_method) 
tmp[upper.tri(tmp)] <-NA
tmp%>% 
  as.data.frame() %>% 
  rownames_to_column("vs1") %>% 
  as_tibble() %>% 
  gather(key="vs2",value="correlation",-vs1) %>% 
  filter(vs1!=vs2) %>% 
  filter(!near(correlation,-1.0)) %>% 
  arrange(-abs(correlation)) %>% 
  #write_csv("../2.Output/05_BayesianNetwork/ME_Trait_correlation.csv")
  ggplot(aes())+
  geom_histogram(aes(x=correlation))


ME_Trait_cor %>% 
  ggplot(aes(x=vs1,y=vs2,fill=correlation))+
  geom_tile()+
  scale_fill_gradient2() + 
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  xlab("")+
  ylab("")

#MEs_withTrait %>% select(-"Sample_name") %>% ggpairs(aes_string(alpha=0.5))
target<-c(str_c("ME",significant_module),c("total_skin_score","Skin_score_at_biopsy_site"))
ME_Trait_graph<-ME_Trait_cor %>% 
  filter(vs1 %in% target & vs2 %in% target) %>% 
  filter(vs1!=vs2) %>% 
  rename(vs1="from",vs2="to") %>% 
  as_tbl_graph()

ME_Trait_graph %>% 
  activate("edges") %>% 
  filter(abs(correlation)>0.6) %>% 
  activate("nodes") %>%
  mutate(color=recode(name, `total_skin_score` = "TRred",`Skin_score_at_biopsy_site` = "TRred", .default = name),
         color=str_sub(color,3,-1)) %>% 
  mutate(label=recode(name, `total_skin_score` = "TRtotal_skin_score",`Skin_score_at_biopsy_site` = "TRSkin_score_at_biopsy_site", .default = name),
         label=str_sub(label,3,-1)) %>% 
  ggraph(layout = "kk")+
  geom_edge_link0(aes(width=abs(correlation),color=correlation>0),alpha=0.8)+
  geom_node_point(aes(color="#404040"),size=21)+
  geom_node_point(aes(color=color),size=20)+
  geom_node_text(aes(label=label),color="#404040")+
  scale_color_identity()+
  theme_graph()+
  theme(legend.position ="none")+
  labs(title="total skinscore correlation P<0.01 | abs(correlation)>0.6")
ggsave(filename = "../2.Output/05_BayesianNetwork/correlation_network.png")





010_Atg5Pathway_Visualization.R





library(tidygraph)
library(ReactomePA)

reactome_result<-read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv") %>% 
  filter(moduleColors_signed==module) %>% .$Entrez_Gene_ID %>% 
  enrichPathway(pvalueCutoff=0.05, readable=T)

bp<-barplot(reactome_result,color = "pvalue",showCategory = 5)
bp+
  ylab("# of Gene")+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))
ggsave("../2.Output/10_YellowModuleAnalysis/Reactome_barplot.png",width=250,height=150,dpi="retina",units = "mm",bg="transparent")
target<-reactome_result
target@result<-reactome_result@result %>% 
  filter(ID=="R-HSA-168928")

cnetplot(reactome_result,showCategory = 2)+ 
  theme(legend.position = 'none',
        panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))
ggsave("../2.Output/10_YellowModuleAnalysis/aatg5_cnetplot_ver1.png",height=150,width=150,dpi="retina",units="mm",bg="transparent")







010_ExosomeMarker_test.R






setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")


exosome_marker<-read_tsv("../1.Data/ExoCarta/ExoCarta_top100_protein_details_5.txt") %>% .$`Gene Symbol`

exosome_marker.entrez<-read_tsv("../1.Data/ExoCarta/ExoCarta_top100_protein_details_5.txt") %>%
  mutate(Entrez_Gene_ID=mapIds(org.Hs.eg.db, 
                               read_tsv("../1.Data/ExoCarta/ExoCarta_top100_protein_details_5.txt") %>% .$`Gene Symbol`,
                               'ENTREZID',"SYMBOL") %>% as.vector()) %>% 
  .$Entrez_Gene_ID
genesymbol.yellow<-read_csv("../2.Output/02_signed_2_unmerged/Modules/Significant/all/yellow.csv") %>% 
  filter(!is.na(Symbol)) %>% 
  distinct(Symbol) %>% .$Symbol   
read_csv("../2.Output/02_signed_2_unmerged/Modules/Significant/all/yellow.csv") %>% 
  mutate(is_exomarker=Symbol %in% exosome_marker) %>% 
  filter(is_exomarker)

read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  filter(Symbol %in% exosome_marker) %>% 
  .$module %>% 
  table()

read_csv("../2.Output/02_signed_2_unmerged/Modules/Significant/all/yellow.csv") %>% 
  mutate(is_exomarker=Entrez_Gene_ID %in% exosome_marker.entrez) %>% 
  filter(is_exomarker)






010_GSEA_YellowModuleConcentration.R





library(org.Hs.eg.db)
library(tidyverse)
library(clusterProfiler)
library(ReactomePA)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

read_csv("../2.Output/02_signed_2_unmerged/Modules/Significant/all/yellow.csv") %>% 
  .$Entrez_Gene_ID %>% 
  as.character()->yellow.entrez

genelist<-read_csv("../2.Output/02_signed_2_unmerged/logFC_ModuleColor.csv") %>% .$logFC
names(genelist)<-read_csv("../2.Output/02_signed_2_unmerged/logFC_ModuleColor.csv") %>% .$Entrez_Gene_ID
gsego.result<-gseGO(sort(genelist,decreasing = TRUE), 
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "BP",
                    pvalueCutoff = 0.05, 
                    pAdjustMethod = "BH"
)
gsego.result.simple<-simplify(gsego.result)


ridgeplot(gsego.result)
dotplot(gsego.result)
as.data.frame(gsego.result) %>% as_tibble() %>% 
  mutate(core_enrichment.list=core_enrichment %>% str_split("/")) %>% 
  mutate(concentration=map(core_enrichment.list,
                           function(elist){
                             yellow.num<-(elist %in% yellow.entrez) %>% sum()
                             total<-length(elist)
                             #return(yellow.num/total)
                             return(c(yellow.num,total))
                           })
  ) %>% 
  mutate(yellow.num=map(concentration,~{paste(as.character(.[[1]]),"/",as.character(.[[2]]))}) %>% as.character(),
         concentration=map(concentration,~{.[[1]]/.[[2]]})%>% as.numeric()) %>% 
  filter(concentration>0.05 & qvalues<0.1) %>% 
  arrange(qvalues) %>% 
  dplyr::select(-core_enrichment.list) %>% 
  write_csv("../2.Output/10_YellowModuleAnalysis/GSEA_GO_highYellowModuleConc.csv")

ridgeplot(gsego.result.simple)
dotplot(gsego.result.simple)
as.data.frame(gsego.result.simple) %>% as_tibble() %>% 
  mutate(core_enrichment.list=core_enrichment %>% str_split("/")) %>% 
  mutate(concentration=map(core_enrichment.list,
                           function(elist){
                             yellow.num<-(elist %in% yellow.entrez) %>% sum()
                             total<-length(elist)
                             #return(yellow.num/total)
                             return(c(yellow.num,total))
                           })
  ) %>% 
  mutate(yellow.num=map(concentration,~{paste(as.character(.[[1]]),"/",as.character(.[[2]]))})%>% as.character(),
         concentration=map(concentration,~{.[[1]]/.[[2]]}) %>% as.numeric()) %>% 
  filter(concentration>0.05 & qvalues<0.1) %>% 
  arrange(qvalues) %>% 
  dplyr::select(-core_enrichment.list) %>% 
  write_csv("../2.Output/10_YellowModuleAnalysis/GSEA_GO_simplified_highYellowModuleConc.csv")

goid.simple.yellowConcentrated<-read_csv("../2.Output/10_YellowModuleAnalysis/GSEA_GO_simplified_highYellowModuleConc.csv") %>% .$ID
gsego.result.simple@result<-gsego.result.simple@result %>% 
  filter(ID %in% goid.simple.yellowConcentrated) %>% 
  mutate(rowname=ID) %>% 
  mutate(core_enrichment.list=core_enrichment %>% str_split("/")) %>% 
  mutate(concentration=map(core_enrichment.list,
                           function(elist){
                             yellow.num<-(elist %in% yellow.entrez) %>% sum()
                             total<-length(elist)
                             return(yellow.num/total)
                             #return(c(yellow.num,total))
                           }) %>% as.numeric(),
         concentration_color=-concentration
  ) %>% 
  arrange(-concentration) %>% 
  as.data.frame() %>% 
  column_to_rownames("rowname")

rp<-ridgeplot(gsego.result.simple,fill = "concentration_color")
rp$data <-rp$data %>% 
  arrange(concentration_color) %>% 
  mutate(category=ordered(category,levels=gsego.result.simple@result$Description%>% rev()) ) 
rp

ggsave("../2.Output/10_YellowModuleAnalysis/GSEA_GO_simplified_highYellowModuleConc.png",width=400,height=300,units = "mm")

ridgeplot(gsego.result)

goid.yellowConcentrated<-read_csv("../2.Output/10_YellowModuleAnalysis/GSEA_highYellowModuleConc.csv") %>% .$ID
gsego.result@result<-gsego.result@result %>% 
  filter(ID %in% goid.yellowConcentrated) %>% 
  mutate(rowname=ID) %>% 
  mutate(core_enrichment.list=core_enrichment %>% str_split("/")) %>% 
  mutate(concentration=map(core_enrichment.list,
                           function(elist){
                             yellow.num<-(elist %in% yellow.entrez) %>% sum()
                             total<-length(elist)
                             return(yellow.num/total)
                             #return(c(yellow.num,total))
                           }) %>% as.numeric(),
         concentration_color=-concentration
  ) %>% 
  arrange(-concentration) %>% 
  filter(concentration>0.12) %>% 
  as.data.frame() %>% 
  column_to_rownames("rowname")
gsego.result@result %>% as_tibble() %>% 
  dplyr::select(-core_enrichment.list) %>% 
  write_csv("../2.Output/10_YellowModuleAnalysis/GSEA_highYellowModuleConc_.csv")

rp<-ridgeplot(gsego.result,fill = "concentration_color")
rp$data <-rp$data %>% 
  arrange(concentration_color) %>% 
  mutate(category=ordered(category,levels=gsego.result@result$Description%>% rev()) ) 
rp

ggsave("../2.Output/10_YellowModuleAnalysis/Image/GSEA_GO_highYellowModuleConc.png",width=250,height=200,units = "mm")

ridgeplot(gsego.result)

#############################################################

gsereactome.result<-gsePathway(sort(genelist,decreasing = TRUE), 
                               organism = "human",
                               pvalueCutoff = 0.05, 
                               pAdjustMethod = "BH"
                               )


ridgeplot(gsereactome.result)
dotplot(gsereactome.result)
as.data.frame(gsereactome.result) %>% as_tibble() %>% 
  mutate(core_enrichment.list=core_enrichment %>% str_split("/")) %>% 
  mutate(concentration=map(core_enrichment.list,
                           function(elist){
                             yellow.num<-(elist %in% yellow.entrez) %>% sum()
                             total<-length(elist)
                             #return(yellow.num/total)
                             return(c(yellow.num,total))
                           })
  ) %>% 
  mutate(yellow.num=map(concentration,~{paste(as.character(.[[1]]),"/",as.character(.[[2]]))}) %>% as.character(),
         concentration=map(concentration,~{.[[1]]/.[[2]]})%>% as.numeric()) %>% 
  filter(concentration>0.05 & qvalues<0.1) %>% 
  arrange(qvalues) %>% 
  dplyr::select(-core_enrichment.list) %>% 
  write_csv("../2.Output/10_YellowModuleAnalysis/GSEA_Reactome_highYellowModuleConc.csv")

id.yellowConcentrated<-read_csv("../2.Output/10_YellowModuleAnalysis/GSEA_Reactome_highYellowModuleConc.csv") %>% .$ID
gsereactome.result@result<-gsereactome.result@result %>% 
  filter(ID %in% id.yellowConcentrated) %>% 
  mutate(rowname=ID) %>% 
  mutate(core_enrichment.list=core_enrichment %>% str_split("/")) %>% 
  mutate(concentration=map(core_enrichment.list,
                           function(elist){
                             yellow.num<-(elist %in% yellow.entrez) %>% sum()
                             total<-length(elist)
                             return(yellow.num/total)
                             #return(c(yellow.num,total))
                           }) %>% as.numeric(),
         concentration_color=-concentration
  ) %>% 
  arrange(-concentration) %>% 
  as.data.frame() %>% 
  column_to_rownames("rowname")

rp.reactome<-ridgeplot(gsereactome.result,fill = "concentration_color")
rp.reactome$data <-rp.reactome$data %>% 
  arrange(concentration_color) %>% 
  mutate(category=ordered(category,levels=gsereactome.result@result$Description%>% rev()) ) 
rp.reactome

ggsave("../2.Output/10_YellowModuleAnalysis/GSEA_Reactome_highYellowModuleCOnc.png",width=400,height=300,units = "mm")


#############################################################

gseKEGG.result<-gseKEGG(sort(genelist,decreasing = TRUE), 
                               organism = "human",
                               pvalueCutoff = 0.05, 
                               pAdjustMethod = "BH"
)


ridgeplot(gseKEGG.result)
dotplot(gseKEGG.result)
as.data.frame(gseKEGG.result) %>% as_tibble() %>% 
  mutate(core_enrichment.list=core_enrichment %>% str_split("/")) %>% 
  mutate(concentration=map(core_enrichment.list,
                           function(elist){
                             yellow.num<-(elist %in% yellow.entrez) %>% sum()
                             total<-length(elist)
                             #return(yellow.num/total)
                             return(c(yellow.num,total))
                           })
  ) %>% 
  mutate(yellow.num=map(concentration,~{paste(as.character(.[[1]]),"/",as.character(.[[2]]))}) %>% as.character(),
         concentration=map(concentration,~{.[[1]]/.[[2]]})%>% as.numeric()) %>% 
  filter(concentration>0.05 & qvalues<0.1) %>% 
  arrange(qvalues) %>% 
  dplyr::select(-core_enrichment.list) %>% 
  write_csv("../2.Output/10_YellowModuleAnalysis/GSEA_KEGG_highYellowModuleConc.csv")

id.yellowConcentrated<-read_csv("../2.Output/10_YellowModuleAnalysis/GSEA_KEGG_highYellowModuleConc.csv") %>% .$ID
gseKEGG.result@result<-gseKEGG.result@result %>% 
  filter(ID %in% id.yellowConcentrated) %>% 
  mutate(rowname=ID) %>% 
  mutate(core_enrichment.list=core_enrichment %>% str_split("/")) %>% 
  mutate(concentration=map(core_enrichment.list,
                           function(elist){
                             yellow.num<-(elist %in% yellow.entrez) %>% sum()
                             total<-length(elist)
                             return(yellow.num/total)
                             #return(c(yellow.num,total))
                           }) %>% as.numeric(),
         concentration_color=-concentration
  ) %>% 
  arrange(-concentration) %>% 
  as.data.frame() %>% 
  column_to_rownames("rowname")

rp.KEGG<-ridgeplot(gseKEGG.result,fill = "concentration_color")
rp.KEGG$data <-rp.KEGG$data %>% 
  arrange(concentration_color) %>% 
  mutate(category=ordered(category,levels=gseKEGG.result@result$Description%>% rev()) ) 
rp.KEGG

ggsave("../2.Output/10_YellowModuleAnalysis/GSEA_KEGG_highYellowModuleCOnc.png",width=400,height=300,units = "mm")


library(GO.db)
select(GO.db,keys = "GO:0030162",columns = )






010_NetworkPlot_Plum1.R





library(WGCNA)
library(tidyverse)
library(ggraph)
library(tidygraph)
library(gplots)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
lname_TOM<-load(file = "../2.Output/02_signed_2_unmerged/TOM.RData")
lname_signed2<-load("../2.Output/02_signed_2_unmerged/signed_2/networkConstruction_StepByStep_signed_2.Rdata")
lname_probes<-load("../2.Output/03_Network_Plot/probes.Rdata")
gene_annotation <- read_csv("../1.Data/GSE58095/Annotation/GeneAnnotation_GPL10558.csv")

kwithin_df<-read_csv("../2.Output/02_signed_2_unmerged/signed_2/Probe_kWitnin_color_signed_2.csv")

modules<-"plum1"
kWithin_top20_probe<-read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv") %>% 
  filter(moduleColors_signed==modules) %>% 
  mutate(kWithin_percent_rank=kWithin %>% percent_rank(),
         percent_rank_top20p=kWithin_percent_rank>0.8,
         kWithin_dense_rank=dense_rank(dplyr::desc(kWithin)),
         dense_rank_top20p=kWithin_dense_rank<=(max(kWithin_dense_rank)*0.2)) %>% 
  filter(dense_rank_top20p) %>% 
  dplyr::select(-percent_rank_top20p) %>% .$ProbeID



#TO_threshold<-0.05
TO_percent_threshold<-0.9
modTOM<-TOM[moduleColors_signed %in% modules,moduleColors_signed %in% modules]
dimnames(modTOM)<-list(probes[moduleColors_signed %in% modules],probes[moduleColors_signed %in% modules])
modTOM[upper.tri(modTOM,diag = TRUE)]<-NA
t_graph<-modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>% 
  mutate(percent=percent_rank(TO),
         TO_top=percent>TO_percent_threshold) %>% 
  #filter(TO_top & TO>TO_threshold & (from %in% kWithin_top20_probe) & (to %in% kWithin_top20_probe)) %>% 
  #filter(percent >= TO_percent_threshold) %>% 
  dplyr::select(-TO) %>% 
  as_tbl_graph(directed=FALSE) %>% 
  activate("nodes") %>% 
  left_join(gene_annotation %>% dplyr::select(ID,Symbol,Entrez_Gene_ID),by=c("name"="ID")) %>% 
  left_join(kwithin_df %>% dplyr::select(ProbeID,moduleColors_signed),by=c("name"="ProbeID")) %E>%
  filter(percent >= TO_percent_threshold) 

t_graph%>% 
  ggraph(layout = "kk")+
  geom_edge_link(aes(),colour=modules,width=1)+
  geom_node_point(aes(color=moduleColors_signed),size=10)+
  geom_node_text(aes(label=Symbol),color="black",size=2)+
  theme_graph()+
  scale_color_identity()+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))
ggsave("../2.Output/19_ImmuneCell-Fibroblast_Signaling/Plum1-Yellow_LeukocyteTransendothelialMigration/TOM_Network_Plum1.png",
       height = 200,width = 200,dpi = 450,units = "mm",bg="transparent")


geneid_hsa04670<-read_csv("../2.Output/11_Plum1ModuleAnalysis/Kegg/KEGG_PathwayEnrichment_plum1.csv") %>% 
  filter(ID=="hsa04670") %>% .$geneID %>% 
  str_split("/") %>%  .[[1]]

kWithin_top50p_probe<-read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv") %>% 
  filter(moduleColors_signed==modules) %>% 
  mutate(kWithin_percent_rank=kWithin %>% percent_rank(),
         percent_rank_top20p=kWithin_percent_rank>0.5,
         kWithin_dense_rank=dense_rank(dplyr::desc(kWithin)),
         dense_rank_top20p=kWithin_dense_rank<=(max(kWithin_dense_rank)*0.5)) %>% 
  filter(dense_rank_top20p) %>% 
  dplyr::select(-percent_rank_top20p) %>% .$ProbeID


t_graph %N>% 
  mutate(isin_hsa04670 = Entrez_Gene_ID %in% geneid_hsa04670,
         color= if_else(isin_hsa04670, "#da5019", moduleColors_signed))%N>% 
  left_join(kwithin_df, by=c("name" = "ProbeID")) %N>% 
  mutate(alpha = if_else(name %in% kWithin_top50p_probe,1 ,0.5 )) %>% 
  ggraph(layout = "kk")+
  geom_edge_link0(aes(),colour=modules,width=1)+
  geom_node_point(aes(color=color, size = kWithin))+
  geom_node_text(aes(label=Symbol),color="black", size= 4)+
  theme_graph()+
  scale_alpha_identity()+
  scale_color_identity()+
  scale_size_continuous(range = c(3,15))+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        legend.position = 'none')

ggsave("../2.Output/19_ImmuneCell-Fibroblast_Signaling/Plum1-Yellow_LeukocyteTransendothelialMigration/TOM_Network_Plum1_withLTM_tver2.png",
       height = 200,width = 200,dpi = 600,units = "mm",bg="transparent")











010_NetworkPlot_Yellow.R





library(WGCNA)
library(tidyverse)
library(ggraph)
library(tidygraph)
library(gplots)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
load("../2.Output/01_Compare_deepsplit_signed/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/Expdata_quontile.Rdata")#Expdata
#Expdata_quontile$E
keep<-rowSums(Expdata_quontile$E>log2(50))>=102
Expdata<-Expdata_quontile$E[keep,]%>% t()
#probes<-colnames(Expdata)
#save(probes,file = "../2.Output/03_Network_Plot/probes.Rdata")
load("../2.Output/03_Network_Plot/probes.Rdata")
lname<-load("../2.Output/02_signed_2_unmerged/signed_2/networkConstruction_StepByStep_signed_2.Rdata")
kwithin_df<-read_csv("../2.Output/02_signed_2_unmerged/signed_2/Probe_kWitnin_color_signed_2.csv")
#TOM <- TOMsimilarityFromExpr(Expdata, networkType = "signed", power = 20, TOMType = "signed")
#save(TOM, file = "../2.Output/02_signed_2_unmerged/TOM.RData")
# ?v?Z?ς݂̏ꍇ??TOM?f?[?^???ǂݍ???
l<-load(file = "../2.Output/02_signed_2_unmerged/TOM.RData")

gene_annotation <- read_csv("../1.Data/GSE58095/Annotation/GeneAnnotation_GPL10558.csv")
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor

# ME???Čv?Z
MEs <- moduleEigengenes(Expdata, moduleColors_signed)$eigengenes
row.names(MEs) <- row.names(Expdata)
module_color <- str_sub(colnames(MEs), start = 3) #?擪2????ME???폜
len <- length(module_color)

modules<-significant_module %>% head(3)
modules<-c("royalblue","lightcyan")
modules<-"yellow"

kWithin_top20_probe<-read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv") %>% 
  filter(moduleColors_signed=="yellow") %>% 
  mutate(kWithin_percent_rank=kWithin %>% percent_rank(),
         percent_rank_top20p=kWithin_percent_rank>0.8,
         kWithin_dense_rank=dense_rank(dplyr::desc(kWithin)),
         dense_rank_top20p=kWithin_dense_rank<=(max(kWithin_dense_rank)*0.2)) %>% 
  filter(dense_rank_top20p) %>% 
  dplyr::select(-percent_rank_top20p) %>% .$ProbeID


goid_list=c("GO:0070062", "GO:0071985", "GO:1904896", "GO:1904903"	)

exosome_related_gene<-select(org.Hs.eg.db,keys=goid_list,keytype="GO",columns = c("SYMBOL","ENTREZID")) %>% distinct(ENTREZID) %>% .$ENTREZID


threshold<-0.1
modTOM<-TOM[moduleColors_signed %in% modules,moduleColors_signed %in% modules]
dimnames(modTOM)<-list(colnames(Expdata)[moduleColors_signed %in% modules],colnames(Expdata)[moduleColors_signed %in% modules])
modTOM[upper.tri(modTOM,diag = TRUE)]<-NA
t_graph<-modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>% 
  mutate(percent=percent_rank(TO),
         TO_top=percent>0.99) %>% 
  filter(TO_top & TO>threshold & (from %in% kWithin_top20_probe) & (to %in% kWithin_top20_probe)) %>% 
  dplyr::select(-TO) %>% 
  as_tbl_graph(directed=FALSE) %>% 
  activate("nodes") %>% 
  left_join(gene_annotation %>% dplyr::select(ID,Symbol,Entrez_Gene_ID),by=c("name"="ID")) %>% 
  left_join(kwithin_df %>% dplyr::select(ProbeID,moduleColors_signed),by=c("name"="ProbeID")) %>% 
  mutate(exosome_related=Entrez_Gene_ID %in% exosome_related_gene)

t_graph%>% 
  ggraph(layout = "kk")+
  geom_edge_link(aes(),alpha=0.8,colour="lightgray")+
  geom_node_point(aes(color=moduleColors_signed),size=10)+
  geom_node_text(aes(label=Symbol),color="black",size=2)+
  theme_graph()+
  scale_color_identity()

#Exosome related
t_graph%>% 
  ggraph(layout = "kk")+
  geom_edge_link(aes(),alpha=0.8,colour="lightgray")+
  geom_node_point(aes(color=exosome_related),size=10)+
  geom_node_text(aes(label=Symbol),color="black",size=2)+
  theme_graph()#+
scale_color_identity()

#?}???쐬?p
color_code_df %>% 
  mutate(x=1:nrow(color_code_df),y=1) %>% 
  ggplot(aes(x,y,color=module))+
  geom_point(size=10)+
  scale_color_manual(values=color_code_df$color_code)+
  theme_classic()


ggsave(filename = "../2.Output/10_YellowModuleAnalysis/Networkplot_ver1.png",width = 1000, height = 1000, units = "mm")

##########################################################################
read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv") %>% 
  mutate(kWithin_percent_rank=kWithin %>% percent_rank(),
         top20p=kWithin_percent_rank>0.8) %>% 
  filter(top20p & moduleColors_signed=="yellow") %>% 
  dplyr::select(-top20p)

modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>% 
  mutate(percent=percent_rank(TO),
         top20=percent>0.99) %>% 
  ggplot(aes(x=TO,fill=top20))+
  geom_histogram(binwidth = 0.0005)

t_graph<-modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>% 
  mutate(percent=percent_rank(TO),
         TO_top=percent>0.99) %>% 
  filter(TO_top  & (from %in% kWithin_top20_probe) & (to %in% kWithin_top20_probe)) %>% 
  dplyr::select(-TO) %>% 
  as_tbl_graph(directed=FALSE) %>% 
  activate("nodes") %>% 
  left_join(gene_annotation %>% select(ID,Symbol),by=c("name"="ID")) %>% 
  left_join(kwithin_df %>% select(ProbeID,moduleColors_signed),by=c("name"="ProbeID"))

t_graph%>% 
  ggraph(layout = "kk")+
  geom_edge_link(aes(),alpha=0.8,colour="lightgray")+
  geom_node_point(aes(color=moduleColors_signed),size=10)+
  geom_node_text(aes(label=Symbol),color="black",size=2)+
  theme_graph()+
  scale_color_identity()

############################################################################

############################################
#Turquoise + yellow??EXOSOME?Ɋ֘A????GENEID?̃l?b?g???[?N

#"GO:0070062" Extracellular Exosome
#"GO:0071985" multivesicular body sorting pathway
#"GO:1904896"	ESCRT complex disassembly
#"GO:1904903"	ESCRT III complex disassembly
goid_list=c("GO:0070062", "GO:0071985", "GO:1904896", "GO:1904903"	)

exosome_related_gene<-select(org.Hs.eg.db,keys=goid_list,keytype="GO",columns = c("SYMBOL","ENTREZID")) %>% distinct(ENTREZID) %>% .$ENTREZID
  #filter(ENTREZID %in% yellow.entrez) %>% 
  #distinct(ENTREZID)
modules<-c("yellow","turquoise")

kWithin_top20_probe.yellow<-read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv") %>% 
  filter(moduleColors_signed=="yellow") %>% 
  mutate(kWithin_percent_rank=kWithin %>% percent_rank(),
         percent_rank_top20p=kWithin_percent_rank>0.8,
         kWithin_dense_rank=dense_rank(dplyr::desc(kWithin)),
         dense_rank_top20p=kWithin_dense_rank<=(max(kWithin_dense_rank)*0.2)) %>% 
  filter(dense_rank_top20p) %>% 
  dplyr::select(-percent_rank_top20p) %>% .$ProbeID

kWithin_top20_probe.turquoise<-read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv") %>% 
  filter(moduleColors_signed=="turquoise") %>% 
  mutate(kWithin_percent_rank=kWithin %>% percent_rank(),
         percent_rank_top20p=kWithin_percent_rank>0.8,
         kWithin_dense_rank=dense_rank(dplyr::desc(kWithin)),
         dense_rank_top20p=kWithin_dense_rank<=(max(kWithin_dense_rank)*0.2)) %>% 
  filter(dense_rank_top20p) %>% 
  dplyr::select(-percent_rank_top20p) %>% .$ProbeID

kWithin_top20_probe<-c(kWithin_top20_probe.yellow,kWithin_top20_probe.turquoise)

threshold<-0.1
modTOM<-TOM[moduleColors_signed %in% modules,moduleColors_signed %in% modules]
dimnames(modTOM)<-list(colnames(Expdata)[moduleColors_signed %in% modules],colnames(Expdata)[moduleColors_signed %in% modules])
modTOM[upper.tri(modTOM,diag = TRUE)]<-NA
t_graph<-modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>% 
  mutate(percent=percent_rank(TO),
         TO_top=percent>0.995) %>% 
  filter(TO_top & TO>threshold & (from %in% kWithin_top20_probe) & (to %in% kWithin_top20_probe)) %>% 
  dplyr::select(-TO) %>% 
  as_tbl_graph(directed=FALSE) %>% 
  activate("nodes") %>% 
  left_join(gene_annotation %>% dplyr::select(ID,Symbol,Entrez_Gene_ID),by=c("name"="ID")) %>% 
  left_join(kwithin_df %>% dplyr::select(ProbeID,moduleColors_signed),by=c("name"="ProbeID")) %>% 
  mutate(exosome_related=Entrez_Gene_ID %in% exosome_related_gene)

t_graph%>% 
  ggraph(layout = "kk")+
  geom_edge_link(aes(),alpha=0.8,colour="lightgray")+
  geom_node_point(aes(color=exosome_related),size=13)+
  geom_node_point(aes(color=moduleColors_signed),size=10)+
  geom_node_text(aes(label=Symbol),color="black",size=2)+
  theme_graph()#+
scale_color_identity()





010_NetworkPlot_Yellow_withATG5pathway.R





library(WGCNA)
library(tidyverse)
library(ggraph)
library(tidygraph)
library(gplots)
library(clusterProfiler)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
lname_TOM<-load(file = "../2.Output/02_signed_2_unmerged/TOM.RData")
lname_signed2<-load("../2.Output/02_signed_2_unmerged/signed_2/networkConstruction_StepByStep_signed_2.Rdata")
lname_probes<-load("../2.Output/03_Network_Plot/probes.Rdata")
gene_annotation <- read_csv("../1.Data/GSE58095/Annotation/GeneAnnotation_GPL10558.csv")

kwithin_df<-read_csv("../2.Output/02_signed_2_unmerged/signed_2/Probe_kWitnin_color_signed_2.csv")

module<-"yellow"
kWithin_top20_probe<-read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv") %>% 
  filter(moduleColors_signed==module) %>% 
  mutate(kWithin_percent_rank=kWithin %>% percent_rank(),
         percent_rank_top20p=kWithin_percent_rank>0.8,
         kWithin_dense_rank=dense_rank(dplyr::desc(kWithin)),
         dense_rank_top20p=kWithin_dense_rank<=(max(kWithin_dense_rank)*0.2)) %>% 
  filter(dense_rank_top20p) %>% 
  dplyr::select(-percent_rank_top20p) %>% .$ProbeID



TO_threshold<-0.1
TO_percent_threshold<-0.95
modTOM<-TOM[moduleColors_signed %in% module,moduleColors_signed %in% module]
dimnames(modTOM)<-list(probes[moduleColors_signed %in% module],probes[moduleColors_signed %in% module])
modTOM[upper.tri(modTOM,diag = TRUE)]<-NA
t_graph<-modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>% 
  as_tbl_graph(directed=FALSE) %N>% 
  left_join(gene_annotation %>% dplyr::select(ID,Symbol,Entrez_Gene_ID),by=c("name"="ID")) %>% 
  left_join(kwithin_df %>% dplyr::select(ProbeID,moduleColors_signed,kWithin),by=c("name"="ProbeID")) #%>% 
  #filter((name %in% kWithin_top20_probe))%E>%
  #mutate(percent=percent_rank(TO),
  #       TO_top=percent>TO_percent_threshold) #%>% 
  #filter(TO_top)

t_graph%>% 
  ggraph(layout = "fr")+
  geom_edge_link(aes(),colour=module,width=1)+
  geom_node_point(aes(color=moduleColors_signed),size=10)+
  geom_node_text(aes(label=Symbol),color="black",size=2)+
  theme_graph()+
  scale_color_identity()
ggsave("../2.Output/11_Plum1ModuleAnalysis/TOM_Network_Plum1.png")


library(ReactomePA)
target_pathwayID<-"R-HSA-168928"
reactome_result<-read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv") %>% 
  filter(moduleColors_signed==module) %>% .$Entrez_Gene_ID %>% 
  enrichPathway(pvalueCutoff=0.05, readable=F)

gene_in_target_pathway<-reactome_result@result %>% as_tibble() %>% 
  filter(ID==target_pathwayID) %>% .$geneID %>% 
  str_split("/") %>% .[[1]]

#f6ca06 yellow
t_graph %N>%
  mutate(is_inPathway=Entrez_Gene_ID %>% as.character() %in% gene_in_target_pathway,
         node_color=if_else(is_inPathway,"#da5019","#f6ca06"))%>% 
  filter(kWithin>=10,
         !is.na(Symbol)) %E>% 
  #arrange(desc(TO)) %>% 
  #top_n(120) %>% 
  filter(percent>0.8) %>% 
  ggraph(layout = "fr")+
  geom_edge_link(aes(),colour="#f6ca06",width=1)+
  geom_node_point(aes(color=node_color),size=10)+
  geom_node_text(aes(label=Symbol),color="black",size=2)+
  theme_graph()+
  scale_color_identity()

kwithin_df %>% 
  filter(moduleColors_signed==module) %>% 
  ggplot(aes(x=kWithin))+
  geom_histogram()
  
t_graph%N>%
  mutate(is_inPathway=Entrez_Gene_ID %>% as.character() %in% gene_in_target_pathway,
         node_color=if_else(is_inPathway,"#da5019","#f6ca06"))%>% 
  filter(#percent_rank(kWithin)>0.8,
         !is.na(Symbol),
         is_inPathway
         ) %E>% 
  filter(percent_rank(TO)>0.9) %>% 
  ggraph(layout = "nicely")+
  #geom_edge_link(aes(),colour="#f6ca06",width=1)+
  geom_node_point(aes(color=node_color),size=10)+
  geom_node_text(aes(label=Symbol),color="black",size=2)+
  theme_graph()+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
              panel.border = element_blank(),
              plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
              panel.grid.major = element_blank(), # get rid of major grid
              panel.grid.minor = element_blank(), # get rid of minor grid
              legend.background = element_rect(fill = "transparent"), # get rid of legend bg
              legend.box.background = element_rect(fill = "transparent"))+
  scale_color_identity()
ggsave()

modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>% 
  #filter(percent_rank(TO)>0.90) %>% 
  as_tbl_graph(directed=FALSE) %N>% 
  left_join(gene_annotation %>% dplyr::select(ID,Symbol,Entrez_Gene_ID),by=c("name"="ID")) %>% 
  left_join(kwithin_df %>% dplyr::select(ProbeID,moduleColors_signed,kWithin),by=c("name"="ProbeID")) %>% 
  mutate(is_inPathway=Entrez_Gene_ID %>% as.character() %in% gene_in_target_pathway,
         node_color=if_else(is_inPathway,"#da5019","#f6ca06"))%>% 
  filter(percent_rank(kWithin)>0.8,
         !is.na(Symbol),
         #is_inPathway
  ) %E>% 
  filter(percent_rank(TO)>0.90) %>% 
  ggraph(layout = "nicely")+
  geom_edge_link(aes(),colour="#f6ca06",width=1)+
  geom_node_point(aes(color=node_color),size=20)+
  geom_node_text(aes(label=Symbol),color="black",size=4)+
  theme_graph()+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
              panel.border = element_blank(),
              plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
              panel.grid.major = element_blank(), # get rid of major grid
              panel.grid.minor = element_blank(), # get rid of minor grid
              legend.background = element_rect(fill = "transparent"), # get rid of legend bg
              legend.box.background = element_rect(fill = "transparent"))+
  scale_color_identity()
ggsave("../2.Output/10_YellowModuleAnalysis/Networkplot_ver2.png",height = 500,width = 500,dpi = 450,units = "mm",bg="transparent")

#####################################################
library(tweenr)
library(igraph)
igraph_layouts <- c('star', 'circle', 'gem', 'dh', 'graphopt', 'grid', 'mds', 
                    'randomly', 'fr', 'kk', 'drl', 'lgl')
igraph_layouts <- sample(igraph_layouts)
graph <- t_graph %N>%
  mutate(is_inPathway=Entrez_Gene_ID %>% as.character() %in% gene_in_target_pathway,
         node_color=if_else(is_inPathway,"#da5019","#f6ca06"))%>% 
  filter(kWithin>=10,
         !is.na(Symbol),
         Symbol!="CEBPG") %E>% 
  #arrange(desc(TO)) %>% 
  #top_n(120) %>% 
  filter(percent>0.9) 
V(graph)$degree <- degree(graph)
layouts <- lapply(igraph_layouts, create_layout, graph = graph)
layouts_tween <- tween_states(c(layouts, layouts[1]), tweenlength = 1, 
                              statelength = 1, ease = 'cubic-in-out', 
                              nframes = length(igraph_layouts) * 16 + 8)
title_transp <- tween_t(c(0, 1, 0, 0, 0), 16, 'cubic-in-out')[[1]]
for (i in seq_len(length(igraph_layouts))) {
  #tmp_layout <- layouts_tween[layouts_tween$.frame == i, ]
  
  layout <- igraph_layouts[ceiling(i / 16)]
  tmp_layout<-create_layout(graph,layout=layout)
  title_alpha <- title_transp[i %% 16]
  p <- ggraph(tmp_layout, 'manual') +
    geom_edge_link(aes(),colour="#f6ca06",width=1)+
    geom_node_point(aes(color=node_color),size=10)+
    geom_node_text(aes(label=Symbol),color="black",size=2)+
    theme_graph()+
    scale_color_identity()+
    theme(legend.position = 'none', 
          plot.title = element_text(colour = alpha('black', title_alpha)))
  plot(p)
}
#####################################################

#exosome

read_csv("../2.Output/02_signed_2_unmerged/Modules/Significant/all/yellow.csv") %>% 
  .$Entrez_Gene_ID %>% 
  as.character()->yellow.entrez

genelist<-read_csv("../2.Output/02_signed_2_unmerged/logFC_ModuleColor.csv") %>% .$logFC
names(genelist)<-read_csv("../2.Output/02_signed_2_unmerged/logFC_ModuleColor.csv") %>% .$Entrez_Gene_ID
gsego.result<-gseGO(sort(genelist,decreasing = TRUE), 
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "BP",
                    pvalueCutoff = 0.05, 
                    pAdjustMethod = "BH"
)

target_goid<-c("GO:1904896", "GO:1904903", "GO:0071985")
GSEA_highYellowModuleConc<-read_csv("../2.Output/10_YellowModuleAnalysis/GSEA/GSEA_highYellowModuleConc.csv")
target_go_geneid<-GSEA_highYellowModuleConc %>% 
  filter(ID %in% target_goid) %>% 
  transmute(core_enrichment=core_enrichment %>% str_split("\\/")) %>% 
  unnest() %>% .$core_enrichment


modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>% 
  #filter(percent_rank(TO)>0.90) %>% 
  as_tbl_graph(directed=FALSE) %N>% 
  left_join(gene_annotation %>% dplyr::select(ID,Symbol,Entrez_Gene_ID),by=c("name"="ID")) %>% 
  left_join(kwithin_df %>% dplyr::select(ProbeID,moduleColors_signed,kWithin),by=c("name"="ProbeID")) %>% 
  mutate(is_inPathway=Entrez_Gene_ID %>% as.character() %in% target_go_geneid,
         node_color=if_else(is_inPathway,"#da5019","#f6ca06"))%>% 
  filter(percent_rank(kWithin)>0.8,
         !is.na(Symbol),
         #is_inPathway
  ) %E>% 
  filter(percent_rank(TO)>0.90) %>% 
  ggraph(layout = "nicely")+
  geom_edge_link(aes(),colour="#f6ca06",width=1)+
  geom_node_point(aes(color=node_color),size=20)+
  geom_node_text(aes(label=Symbol),color="black",size=4)+
  theme_graph()+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))+
  scale_color_identity()
ggsave("../2.Output/10_YellowModuleAnalysis/Networkplot_exosome_ver1.png",height = 500,width = 500,dpi = 450,units = "mm",bg="transparent")











010_PathwayAnalysis_yellow.R





library(tidyverse)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

c("#777777","#f79646","#0086ab","#da6272","#9bbb59","#bfbfbf")


yellow_reactome<-read_csv("../2.Output/10_YellowModuleAnalysis/Reactome/Reactome_PathwayEnrichment_yellow.csv")
yellow_kegg <- read_csv("../2.Output/10_YellowModuleAnalysis/KEGG/KEGG_PathwayEnrichment_yellow.csv")

yellow_reactome %>% 
  mutate(Description=Description %>% str_wrap(width=30)) %>% 
  separate(GeneRatio, into = c("numerator","denominator"),sep="/") %>% 
  mutate(numerator= numerator %>% as.numeric(),
         denominator = denominator %>% as.numeric(),
         GeneRatio= numerator / denominator )%>% 
  head(4) %>% 
  ggplot(aes(x=Description %>% reorder(-pvalue), y= -log10(qvalue)))+
  geom_bar(aes(fill=GeneRatio),stat = "identity")+
  theme_minimal()+theme(panel.grid=element_blank())+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent"),
        #legend.position = c(1,0.1),
        #legend.justification = c(0,0),
        axis.text=element_text(size=12))+
  scale_fill_gradient2(low = "#f79646",high = "#0086ab")+   #,"#9bbb59"
  scale_y_continuous(expand = c(0, 0))+
  coord_flip()+
  xlab("")

ggsave("../2.Output/10_YellowModuleAnalysis/Yellow_reactomePA.png",
       dpi = 320, width = 140, height = 110, units = "mm",  bg = "transparent")


yellow_reactome %>% 
  mutate(Description=Description %>% str_wrap(width=30)) %>% 
  separate(GeneRatio, into = c("numerator","denominator"),sep="/") %>% 
  mutate(numerator= numerator %>% as.numeric(),
         denominator = denominator %>% as.numeric(),
         GeneRatio= numerator / denominator )%>% 
  head(4) %>% 
  mutate(color=c("#7fc2d5","#0086ab", "#7fc2d5", "#7fc2d5")) %>% 
  ggplot(aes(x=Description %>% reorder(-pvalue), y= -log10(qvalue)))+
  geom_bar(aes(fill=color),stat = "identity")+
  theme_minimal()+theme(panel.grid=element_blank())+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent"),
        #legend.position = c(1,0.1),
        #legend.justification = c(0,0),
        axis.text=element_text(size=16),
        plot.margin=unit(c(1,1,1,1),"cm"))+
  scale_fill_identity()+   #,"#9bbb59"
  scale_y_continuous(expand = c(0, 0))+
  coord_flip()+
  xlab("")

ggsave("../2.Output/10_YellowModuleAnalysis/Yellow_reactomePA_ver2.png",
       dpi = 320, width = 160, height = 150, units = "mm",  bg = "transparent")


######################################################
colors=c(rep("#5f5f5f",2),rep("#0086ab",3))#[order(c(rep("#0086ab",3),rep("#7fc2d5",3)))]
character_size= 12

yellow_reactome %>% 
  mutate(Description=Description %>% str_wrap(width=30)) %>% 
  separate(GeneRatio, into = c("numerator","denominator"),sep="/") %>% 
  mutate(numerator= numerator %>% as.numeric(),
         denominator = denominator %>% as.numeric(),
         GeneRatio= numerator / denominator )%>% 
  head(5) %>% 
  mutate(color=c(rep("#0086ab",3),rep("#7fc2d5",2))) %>% 
  ggplot(aes(x=Description %>% reorder(-pvalue), y= -log10(qvalue)))+
  geom_bar(aes(fill=color),stat = "identity")+
  theme_minimal()+theme(panel.grid=element_blank())+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent"),
        #legend.position = c(1,0.1),
        #legend.justification = c(0,0),
        #axis.text=element_text(size=16),
        axis.text.y = element_text(size=character_size,color = colors,lineheight = 1.5),
        axis.title.x = element_text(size=character_size, color= "#5f5f5f"),
        plot.margin=unit(c(1,1,1,1),"cm"))+
  scale_fill_identity()+   #,"#9bbb59"
  scale_y_continuous(expand = c(0, 0))+
  coord_flip()+
  xlab("")

ggsave("../2.Output/10_YellowModuleAnalysis/Yellow_reactomePA_ver7.png",
       dpi = 320, width = 150, height = 100, units = "mm",  bg = "transparent")

yellow_kegg %>% 
  mutate(Description=Description %>% str_wrap(width=30)) %>% 
  separate(GeneRatio, into = c("numerator","denominator"),sep="/") %>% 
  mutate(numerator= numerator %>% as.numeric(),
         denominator = denominator %>% as.numeric(),
         GeneRatio= numerator / denominator )%>% 
  head(5) %>% 
  mutate(color=c(rep("#0086ab",3),rep("#7fc2d5",2))) %>% 
  ggplot(aes(x=Description %>% reorder(-pvalue), y= -log10(qvalue)))+
  geom_bar(aes(fill=color),stat = "identity")+
  theme_minimal()+theme(panel.grid=element_blank())+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent"),
        #legend.position = c(1,0.1),
        #legend.justification = c(0,0),
        #axis.text=element_text(size=16),
        axis.text.y = element_text(size=character_size,color = colors,lineheight = 1.5),
        axis.title.x = element_text(size=character_size, color= "#5f5f5f"),
        plot.margin=unit(c(1,1,1,1),"cm"))+
  scale_fill_identity()+   #,"#9bbb59"
  scale_y_continuous(expand = c(0, 0))+
  coord_flip()+
  xlab("")

ggsave("../2.Output/10_YellowModuleAnalysis/Yellow_keggPA_ver3.png",
       dpi = 320, width = 150, height = 90, units = "mm",  bg = "transparent")






010_test.R





library(org.Hs.eg.db)
library(tidyverse)
library(clusterProfiler)
library(ReactomePA)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
geneid_list<-read_csv("../2.Output/10_YellowModuleAnalysis/GSEA_highYellowModuleConc.csv") %>% 
  .$core_enrichment %>% .[[1]] %>% 
  str_split("/") %>% .[[1]]
read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  filter(Entrez_Gene_ID %in% geneid_list)

GSEA_result<-read_csv("../2.Output/10_YellowModuleAnalysis/GSEA_highYellowModuleConc.csv")
geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv")
for (i in 1:nrow(GSEA_result)){
  print(GSEA_result[i,2])
  geneid_list<-GSEA_result %>% 
    .$core_enrichment %>% .[[i]] %>% 
    str_split("/") %>% .[[1]]
  geneid2modulecolor %>% 
    filter(Entrez_Gene_ID %in% geneid_list) %>% 
    .$module %>% 
    table() %>% 
    print()
}

read_csv("../2.Output/02_signed_2_unmerged/signed_2/Probe_kWitnin_color_signed_2.csv") %>% 
  left_join(read_csv("../1.Data/GSE58095/Annotation/GeneAnnotation_GPL10558.csv") %>% dplyr::select(ID,Entrez_Gene_ID),by=c("ProbeID"="ID")) %>% 
  write_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv")

read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv") %>% 
  mutate(kWithin_percent_rank=kWithin %>% percent_rank(),
         top20p=kWithin_percent_rank>0.8) %>% 
  arrange(desc(kWithin_percent_rank)) %>% 
  ggplot(aes(x=kWithin,fill=top20p))+
  geom_histogram()

geneid2KwithinModulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv") %>% 
  mutate(Entrez_Gene_ID=Entrez_Gene_ID %>% as.character())

for (i in 1:nrow(GSEA_result)){
  print(GSEA_result[i,2])
  geneid_list<-GSEA_result %>% 
    .$core_enrichment %>% .[[i]] %>% 
    str_split("/") %>% .[[1]]
  geneid2KwithinModulecolor %>% 
    filter(Entrez_Gene_ID %in% geneid_list) %>% 
    arrange(desc(kWithin))%>% 
    print()
}
geneid2KwithinModulecolor%>% 
  filter(Entrez_Gene_ID %in% geneid_list)



library(org.Hs.eg.db)
read_csv("../2.Output/02_signed_2_unmerged/Modules/Significant/all/yellow.csv") %>% 
  .$Entrez_Gene_ID %>% 
  as.character()->yellow.entrez
columns(org.Hs.eg.db)
select(org.Hs.eg.db,keys=GSEA_result$ID[3],keytype="GO",columns = c("SYMBOL","ENTREZID")) %>% 
  filter(ENTREZID %in% yellow.entrez)
keys(org.Hs.eg.db)
org.Hs.eg.db

(yellow.entrez %in%geneid_list) %>% 
  sum()
################################################################################

"GO:0070062"

select(org.Hs.eg.db,keys="GO:0070062",keytype="GO",columns = c("SYMBOL","ENTREZID")) %>% 
  filter(ENTREZID %in% yellow.entrez) %>% 
  distinct(ENTREZID)

kWithin_top20_probe.yellow<-read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv") %>% 
  filter(moduleColors_signed=="yellow") %>% 
  mutate(kWithin_percent_rank=kWithin %>% percent_rank(),
         percent_rank_top20p=kWithin_percent_rank>0.8,
         kWithin_dense_rank=dense_rank(dplyr::desc(kWithin)),
         dense_rank_top20p=kWithin_dense_rank<=(max(kWithin_dense_rank)*0.2)) %>% 
  filter(dense_rank_top20p) %>% 
  dplyr::select(-percent_rank_top20p) %>% .$ProbeID

kWithin_top20_probe.turquoise<-read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv") %>% 
  filter(moduleColors_signed=="turquoise") %>% 
  mutate(kWithin_percent_rank=kWithin %>% percent_rank(),
         percent_rank_top20p=kWithin_percent_rank>0.8,
         kWithin_dense_rank=dense_rank(dplyr::desc(kWithin)),
         dense_rank_top20p=kWithin_dense_rank<=(max(kWithin_dense_rank)*0.2)) %>% 
  filter(dense_rank_top20p) %>% 
  dplyr::select(-percent_rank_top20p) %>% .$ProbeID

geneid2KwithinModulecolor %>% 
  filter(Entrez_Gene_ID %in% (select(org.Hs.eg.db,keys="GO:0070062",keytype="GO",columns = c("SYMBOL","ENTREZID")) %>% .$ENTREZID)) %>% 
  arrange(desc(kWithin)) %>% #.$moduleColors_signed %>% table()
  filter(moduleColors_signed%in% c("turquoise","yellow")) %>% 
  left_join(read_csv("../1.Data/GSE58095/Annotation/GeneAnnotation_GPL10558.csv") %>% dplyr::select(ID,Symbol),by=c("ProbeID"="ID")) %>% 
  mutate(top20=ProbeID %in% c(kWithin_top20_probe.yellow,kWithin_top20_probe.turquoise))->a
  filter(ProbeID %in% c(kWithin_top20_probe.yellow,kWithin_top20_probe.turquoise))
  filter(moduleColors_signed=="yellow")
  


genelist<-read_csv("../2.Output/02_signed_2_unmerged/logFC_ModuleColor.csv") %>% .$logFC
names(genelist)<-read_csv("../2.Output/02_signed_2_unmerged/logFC_ModuleColor.csv") %>% .$Entrez_Gene_ID
gsego.result<-gseGO(sort(genelist,decreasing = TRUE), 
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "CC",
                    pvalueCutoff = 0.05, 
                    pAdjustMethod = "BH"
)

gsego.cc.result<-gsego.result

as.data.frame(gsego.cc.result) %>% as_tibble() %>% 
  mutate(core_enrichment.list=core_enrichment %>% str_split("/")) %>% 
  mutate(concentration=map(core_enrichment.list,
                           function(elist){
                             yellow.num<-(elist %in% yellow.entrez) %>% sum()
                             total<-length(elist)
                             #return(yellow.num/total)
                             return(c(yellow.num,total))
                           })
  ) %>% 
  mutate(yellow.num=map(concentration,~{paste(as.character(.[[1]]),"/",as.character(.[[2]]))}) %>% as.character(),
         concentration=map(concentration,~{.[[1]]/.[[2]]})%>% as.numeric()) %>% 
  filter(concentration>0.05 & qvalues<0.1) %>% 
  arrange(qvalues) %>% 
  dplyr::select(-core_enrichment.list) %>% 
  write_csv("../2.Output/10_YellowModuleAnalysis/GSEA_GO_CC_highYellowModuleConc.csv")

goid.yellowConcentrated<-read_csv("../2.Output/10_YellowModuleAnalysis/GSEA_GO_CC_highYellowModuleConc.csv") %>% .$ID
gsego.cc.result@result<-gsego.cc.result@result %>% 
  filter(ID %in% goid.yellowConcentrated) %>% 
  mutate(rowname=ID) %>% 
  mutate(core_enrichment.list=core_enrichment %>% str_split("/")) %>% 
  mutate(concentration=map(core_enrichment.list,
                           function(elist){
                             yellow.num<-(elist %in% yellow.entrez) %>% sum()
                             total<-length(elist)
                             return(yellow.num/total)
                             #return(c(yellow.num,total))
                           }) %>% as.numeric(),
         concentration_color=-concentration
  ) %>% 
  arrange(-concentration) %>% 
  filter(concentration>0.12) %>% 
  as.data.frame() %>% 
  column_to_rownames("rowname")
gsego.cc.result@result %>% as_tibble() %>% 
  dplyr::select(-core_enrichment.list) %>% 
  write_csv("../2.Output/10_YellowModuleAnalysis/GSEA_CC_highYellowModuleConc_.csv")

rp<-ridgeplot(gsego.cc.result,fill = "concentration_color")
rp$data <-rp$data %>% 
  arrange(concentration_color) %>% 
  mutate(category=ordered(category,levels=gsego.cc.result@result$Description%>% rev()) ) 
rp

ggsave("../2.Output/10_YellowModuleAnalysis/Image/GSEA_GO_highYellowModuleConc.png",width=250,height=200,units = "mm")

read_csv("../2.Output/02_signed_2_unmerged/Modules/Significant/all/turquoise.csv") %>% 
  .$Entrez_Gene_ID %>% 
  as.character()->turquoise.entrez

ego<-enrichGO(gene = c(yellow.entrez,turquoise.entrez),
              OrgDb = org.Hs.eg.db,
              keyType = "ENTREZID",
              ont = "CC",
              pvalueCutoff = 0.1, 
              pAdjustMethod = "BH",
              readable = TRUE
)

############################################
#"GO:0070062" Extracellular Exosome
#"GO:0071985" multivesicular body sorting pathway
#"GO:1904896"	ESCRT complex disassembly
#"GO:1904903"	ESCRT III complex disassembly
goid_list=c("GO:0070062", "GO:0071985", "GO:1904896", "GO:1904903"	)

select(org.Hs.eg.db,keys=goid_list,keytype="GO",columns = c("SYMBOL","ENTREZID")) %>% 
  filter(ENTREZID %in% yellow.entrez) %>% 
  distinct(ENTREZID)


############################################
#MSigDB canpnical pathway
library(DOSE)
data(geneList)
gene <- names(geneList)[abs(geneList) > 2]
gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c5 <- read.gmt(gmtfile)
egmt.test <- enricher(gene, TERM2GENE=c5)
head(as.data.frame(egmt.test))


gene.yellow<-read_csv("../2.Output/02_signed_2_unmerged/Modules/Significant/all/yellow.csv") %>% 
  distinct(Entrez_Gene_ID) %>% 
  filter(!is.na(Entrez_Gene_ID)) %>% 
  .$Entrez_Gene_ID %>% as.character() 

gene.yellow<-read_csv("../2.Output/02_signed_2_unmerged/Modules/Significant/all/plum1.csv") %>% 
  distinct(Entrez_Gene_ID) %>% 
  filter(!is.na(Entrez_Gene_ID)) %>% 
  .$Entrez_Gene_ID %>% as.character() 

gene.yellow<-read_csv("../2.Output/02_signed_2_unmerged/Modules/Significant/all/yellow.csv") %>% 
  distinct(Entrez_Gene_ID) %>% 
  filter(!is.na(Entrez_Gene_ID)) %>% 
  .$Entrez_Gene_ID %>% as.character() 
#gmtfile <- system.file("extdata", "c5.cc.v5.0.entrez.gmt", package="clusterProfiler")
c2 <- read.gmt("../1.Data/MSigDB_CP/c2.cp.v7.0.entrez.gmt")
egmt <- enricher(unique(gene.yellow), TERM2GENE=c2,
                 pvalueCutoff = 0.5,qvalueCutoff = 0.5)
head(as.data.frame(egmt))

egmt %>% as_tibble()->a






010_YellowKwithinTop.R





library(org.Hs.eg.db)
library(tidyverse)
library(clusterProfiler)
library(ReactomePA)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

lname<-load("../2.Output/02_signed_2_unmerged/signed_2/networkConstruction_StepByStep_signed_2.Rdata")

#read_csv("../2.Output/02_signed_2_unmerged/Modules/Significant/all/yellow.csv")
read_csv("../2.Output/02_signed_2_unmerged/signed_2/Probe_kWitnin_color_signed_2.csv") %>% 
  filter(moduleColors_signed=="yellow") %>% 
  arrange(desc(kWithin)) %>% 
  filter(percent_rank(kWithin) >= 0.8) %>% .$ProbeID->yellow_top20p
  #filter(top_n(kWithin)/n(kWithin)>0.2)

read_csv("../2.Output/02_signed_2_unmerged/signed_2/Probe_kWitnin_color_signed_2.csv") %>% 
  filter(moduleColors_signed=="yellow") %>% 
  arrange(desc(kWithin)) %>% 
  mutate(percent=percent_rank(kWithin)) %>% 
  filter(percent> 0.8)->a2

read_csv("../2.Output/02_signed_2_unmerged/Modules/Significant/all/yellow.csv") %>% 
  filter(ILMNID %in% yellow_top20p) %>% 
  .$Entrez_Gene_ID->yellow_top20p.entrez

read_csv("../2.Output/02_signed_2_unmerged/Modules/Significant/all/yellow.csv") %>% 
  filter(ILMNID %in% yellow_top20p) %>% 
  .$Symbol->yellow_top20p.symbol

ego<-enrichGO(gene = yellow_top20p.entrez,
         OrgDb = org.Hs.eg.db,
         keyType = "ENTREZID",
         ont = "ALL",
         pvalueCutoff = 0.05, 
         pAdjustMethod = "BH",
         readable = TRUE
         )


ep<-enrichPathway(gene = yellow_top20p.entrez,
                  organism = "human",
                  pvalueCutoff = 0.05, 
                  pAdjustMethod = "BH",
                  readable = TRUE
                  )
ekegg<-enrichKEGG(gene = yellow_top20p.entrez,
                  organism = "hsa",
                  pvalueCutoff = 0.05, 
                  pAdjustMethod = "BH"
                  )

clusterProfiler::dotplot(ego)
emapplot(ego)
cnetplot(ego)

clusterProfiler::dotplot(ep)
emapplot(ep)
cnetplot(ep)

clusterProfiler::dotplot(ekegg)
emapplot(ekegg)
cnetplot(ekegg)

#library(RDAVIDWebService)
#enrichDAVID(yellow_top20p.entrez)
genelist<-read_csv("../2.Output/02_signed_2_unmerged/logFC_ModuleColor.csv") %>% .$logFC
names(genelist)<-read_csv("../2.Output/02_signed_2_unmerged/logFC_ModuleColor.csv") %>% .$Entrez_Gene_ID
gsego.result<-gseGO(sort(genelist,decreasing = TRUE), 
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "BP",
                    pvalueCutoff = 0.05, 
                    pAdjustMethod = "BH"
                    )

ridgeplot(gsego.result)
dotplot(gsego.result)
as.data.frame(gsego.result) %>% as_tibble() %>% 
  mutate(core_enrichment.list=core_enrichment %>% str_split("/")) %>% 
  mutate(percent=map(core_enrichment.list,function(elist){
    return(length(elist))
    return((elist %in% yellow.entrez) %>% sum())
    yellow.num<-(elist %in% yellow.entrez) %>% sum()
    total<-length(elist)
    return(list(yellow.num,total))
  })) %>% #.$percent
  separate(col = percent,into = c("a","b"))->s

  #dplyr::select(core_enrichment.list)


read_csv("../2.Output/02_signed_2_unmerged/Modules/Significant/all/yellow.csv") %>% 
  .$Entrez_Gene_ID %>% 
  as.character()->yellow.entrez

for (i in 1:nrow(a)){
  num<-(a$core_enrichment.list[[i]] %in% yellow.entrez) %>% sum()
  (num/length(a$core_enrichment.list[[i]])) %>% 
  print()
}



















012_Barplot_enrichresult.R





library(tidyverse)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

c("#777777","#f79646","#0086ab","#da6272","#9bbb59","#bfbfbf")

#kegg_pathway<-read_csv("../2.Output/11_Plum1ModuleAnalysis/Kegg/KEGG_PathwayEnrichment_plum1.csv")
lname<-load("../2.Output/12_Comparison_withGWAS/1.Data/enrichGO_SScDM_ALL.Rdata")
lname<-load("../2.Output/12_Comparison_withGWAS/1.Data/enrichGO_SScDM_BP.Rdata")
lname<-load("../2.Output/12_Comparison_withGWAS/1.Data/enrichGO_SScSLE_BP.Rdata")


SScDM_SSc<-part2GOenrichment_SScDM$GO[[1]]@result
colors=c(rep("#0086ab",5))#c(rep("#0086ab",3),rep("#7fc2d5",2))
character_size= 12
##################################################################
SScDM_SScDM %>% 
  mutate(Description=Description %>% str_wrap(width=30)) %>% 
  separate(GeneRatio, into = c("numerator","denominator"),sep="/") %>% 
  mutate(numerator= numerator %>% as.numeric(),
         denominator = denominator %>% as.numeric(),
         GeneRatio= numerator / denominator )%>% 
  head(4) %>% 
  mutate(color=c(rep("#0086ab",4))) %>% 
  ggplot(aes(x=Description %>% reorder(-pvalue), y= -log10(qvalue)))+
  geom_bar(aes(fill=color),stat = "identity")+
  theme_minimal()+theme(panel.grid=element_blank())+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent"),
        #legend.position = c(1,0.1),
        #legend.justification = c(0,0),
        #axis.text=element_text(size=16),
        axis.text.y = element_text(size=character_size,color = "#5f5f5f",lineheight = 1.5),
        axis.title.x = element_text(size=character_size, color= "#5f5f5f"),
        plot.margin=unit(c(1,1,1,1),"cm"))+
  scale_fill_identity()+   #,"#9bbb59"
  scale_y_continuous(expand = c(0, 0))+
  coord_flip()+
  xlab("")
ggsave("../2.Output/12_Comparison_withGWAS/2.Image/Barplot_SScDM_SSc.png",
       dpi = 320, width = 150, height = 100, units = "mm",  bg = "transparent")

#################################################################
plotbar_by_enrichresult<-function(enrichresult,savefile,barcolor="#0086ab"){
  topn=4
  enrichresult %>% 
    mutate(Description=Description %>% str_wrap(width=42)) %>% 
    separate(GeneRatio, into = c("numerator","denominator"),sep="/") %>% 
    mutate(numerator= numerator %>% as.numeric(),
           denominator = denominator %>% as.numeric(),
           GeneRatio= numerator / denominator )%>% 
    head(topn) %>% 
    mutate(color=barcolor) %>% 
    ggplot(aes(x=Description %>% reorder(-pvalue), y= -log10(qvalue)))+
    geom_bar(aes(fill=color),stat = "identity")+
    theme_minimal()+theme(panel.grid=element_blank())+
    theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
          panel.border = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
          #legend.box.background = element_rect(fill = "transparent"),
          #legend.position = c(1,0.1),
          #legend.justification = c(0,0),
          #axis.text=element_text(size=16),
          axis.text.y = element_text(size=character_size,color = "#5f5f5f",lineheight = 1.5),
          axis.title.x = element_text(size=character_size, color= "#5f5f5f"),
          plot.margin=unit(c(1,1,1,1),"cm"))+
    scale_fill_identity()+   #,"#9bbb59"
    scale_y_continuous(expand = c(0, 0))+
    coord_flip()+
    xlab("")
  ggsave(savefile,dpi = 320, width = 150, height = 100, units = "mm",  bg = "transparent")
}

SScDM_SSc<-simplify(part2GOenrichment_SScDM$GO[[1]])@result
savefile<-"../2.Output/12_Comparison_withGWAS/2.Image/Barplot_SScDM_SSc.png"
plotbar_by_enrichresult(SScDM_SSc,savefile,barcolor = alpha('#0086ab',0.3) )

SScDM_SScDM<-simplify(part2GOenrichment_SScDM$GO[[2]])@result
savefile<-"../2.Output/12_Comparison_withGWAS/2.Image/Barplot_SScDM_SScDM_ver2.png"
plotbar_by_enrichresult(SScDM_SScDM,savefile, barcolor = "#d3cdeb" )

SScDM_DM<-simplify(part2GOenrichment_SScDM$GO[[3]])@result
savefile<-"../2.Output/12_Comparison_withGWAS/2.Image/Barplot_SScDM_DM_ver2.png"
plotbar_by_enrichresult(SScDM_DM,savefile,barcolor= alpha("#da6272",0.3))
#alpha('#0086ab',0.3),  alpha("#da6272",0.3))


SScSLE_SSc<-part2GOenrichment_SScSLE$GO[[1]]@result
savefile<-"../2.Output/12_Comparison_withGWAS/2.Image/Barplot_SScSLE_SSc.png"
plotbar_by_enrichresult(SScSLE_SSc,savefile, barcolor = alpha('#0086ab',0.3) )

SScSLE_SScSLE<-part2GOenrichment_SScSLE$GO[[2]]@result
savefile<-"../2.Output/12_Comparison_withGWAS/2.Image/Barplot_SScSLE_SScSLE_ver2.png"
plotbar_by_enrichresult(SScSLE_SScSLE,savefile, barcolor = "#d0e6db" )

SScSLE_SLE<-part2GOenrichment_SScSLE$GO[[4]]@result
savefile<-"../2.Output/12_Comparison_withGWAS/2.Image/Barplot_SScSLE_SLE_ver2.png"
plotbar_by_enrichresult(SScSLE_SLE,savefile,barcolor= alpha('#9bbb59',0.3))






012_GWASgene_in_GSEA.R





library(tidyverse)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

gwas_gene<-read_csv("../2.Output/12_GWAS_compare2SLE_?畆?؉?/GWASgene_in_Module_gwascatalog.csv") %>% .$Entrez_Gene_ID
gse_GO_all<-read_csv("../2.Output/12_GWAS_compare2SLE_?畆?؉?/gse_GO_all.csv") %>% 
  mutate(core_enrichment_list=str_split(core_enrichment,"\\/")) 
gse_GO_all<-gse_GO_all %>% 
  mutate(gwas_gene=map(.x = gse_GO_all$core_enrichment_list,.f = ~intersect(x = .x,y = gwas_gene)),)
gse_GO_all %>% 
  mutate(is_empty=map_lgl(gse_GO_all$gwas_gene,.f=~(length(.x)==0))) %>% 
  filter(!is_empty) %>% 
  mutate(gwas_gene=map_chr(gwas_gene,.f=~str_c(.x,collapse = "/"))) %>% 
  dplyr::select(-core_enrichment_list) %>% 
  write_csv("../2.Output/12_GWAS_compare2SLE_?畆?؉?/GSEA_withGWASgene.csv")
 






012_GWASgene_in_Module.R





library(tidyverse)
library(org.Hs.eg.db)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
geneid2Kwithin<-read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv")
geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") 
SSc_GWAS_catalog<-read_csv("../2.Output/12_GWAS_compare2SLE_?畆?؉?/efotraits_EFO_0000717_SystemicScleoderma.csv")
NatCom_GWAS<-read_csv("../2.Output/12_GWAS_compare2SLE_?畆?؉?/EL_et_al_2019_NatCom/table2.csv")
NatCom_GWAS_locus<-NatCom_GWAS %>% 
  dplyr::select(2) %>% 
  mutate(locus=`Credible set locus` %>% str_split("-") ) %>% 
  unnest() %>% 
  mutate(locus =locus %>% str_remove_all("\\?")) %>% 
  .$locus
NatCom_GWAS_locus<-NatCom_GWAS_locus[!is.na(NatCom_GWAS_locus)]

SSc_GWAS_gene_list<-SSc_GWAS_catalog %>% 
  dplyr::select(1:2,`Mapped gene`) %>% 
  mutate(mapped_gene_mod=`Mapped gene` %>% 
           str_remove(",") %>% str_trim(side="both")) %>%
  .$mapped_gene_mod %>% 
  unique() 

SSc_GWAS_EntrezGeneID<-mapIds(org.Hs.eg.db, SSc_GWAS_gene_list, 'ENTREZID', 'SYMBOL')
SSc_GWAS_EntrezGeneID<-SSc_GWAS_EntrezGeneID[!is.na(SSc_GWAS_EntrezGeneID)]
tibble(Symbol=names(SSc_GWAS_EntrezGeneID),
       GeneID=SSc_GWAS_EntrezGeneID) %>% 
  write_csv("../2.Output/12_GWAS_compare2SLE_?畆?؉?/SSc_GWAS_Genes.csv")

geneid2modulecolor%>% 
  dplyr::select(1:4) %>% 
  left_join(geneid2Kwithin %>% dplyr::select(1:4), by= c("probe"="ProbeID")) %>% 
  mutate(Entrez_Gene_ID =Entrez_Gene_ID %>% as.character()) %>% 
  filter(Symbol %in% SSc_GWAS_gene_list)

geneid2modulecolor%>% 
  dplyr::select(1:4) %>% 
  left_join(geneid2Kwithin %>% dplyr::select(1:4), by= c("probe"="ProbeID")) %>% 
  mutate(Entrez_Gene_ID =Entrez_Gene_ID %>% as.character()) %>% 
  filter(Entrez_Gene_ID  %in% SSc_GWAS_EntrezGeneID) %>% 
  write_csv("../2.Output/12_GWAS_compare2SLE_?畆?؉?/GWASgene_in_Module_gwascatalog.csv")


geneid2modulecolor%>% 
  dplyr::select(1:4) %>% 
  left_join(geneid2Kwithin %>% dplyr::select(1:4), by= c("probe"="ProbeID")) %>% 
  mutate(Entrez_Gene_ID =Entrez_Gene_ID %>% as.character()) %>% 
  filter(Symbol %in% NatCom_GWAS_locus)%>% 
  write_csv("../2.Output/12_GWAS_compare2SLE_?畆?؉?/GWASgene_in_Module_NatCom2019.csv")


Gwas_gene<-read_csv("../2.Output/SSc_GWAS_Genes.csv") %>% 
  .$Symbol       

read_tsv("../2.Output/Participating Molecules [R-HSA-2028269].tsv",skip = 2) %>% 
  separate(col = `ADP(3-) [ChEBI:456216]`, into = c("a","Symbol"),sep = " ") %>% 
  filter(Symbol %in% Gwas_gene)

Gwas_gene_ID<-read_csv("../2.Output/SSc_GWAS_Genes.csv") %>% 
  .$GeneID %>% 
  as.character()

read_csv("../2.Output/Hippo.csv",col_names =F) %>% 
  mutate(X1 = X1 %>% str_replace_all("\\?\\?","")) %>% 
  filter(X1 %in% Gwas_gene_ID)

read_tsv("../2.Output/Participating Molecules [R-HSA-112310].tsv") %>% 
  separate(col = MoleculeName, into = c("a","Symbol"),sep = " ") %>% 
  filter(Symbol %in% Gwas_gene)






012_Make_GWASgene_Module_Graph.R





library(tidyverse)
library(tidygraph)
library(ggraph)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

entrez2part_2dis<-read_csv("../2.Output/12_Comparison_withGWAS/1.Data/Entrezid2partOfVenn_2diseases.csv")
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor
geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  filter(module!="grey") %>% 
  filter(module %in% significant_module) %>% 
  distinct(Entrez_Gene_ID,.keep_all = T)

graph<-entrez2part_2dis %>% 
  dplyr::rename(from="module",to="Entrez_Gene_ID") %>%
  as_tbl_graph(directed=T)

nodeID2module<-graph %N>%
  as.data.frame() %>% 
  filter(name %in% significant_module) %>% 
  rownames_to_column("nodeID") %>% as_tibble() %>% 
  mutate(nodeID=nodeID %>% as.integer())

graph<-graph %N>% 
  left_join(geneid2modulecolor %>% mutate(Entrez_Gene_ID=Entrez_Gene_ID %>% as.character()),by=c("name"="Entrez_Gene_ID")) %>% 
  mutate(color=if_else(is.na(Symbol), name, module),
         label=if_else(is.na(Symbol), name, Symbol),
         size=if_else(is.na(Symbol), 20, 10)) %>% 
  activate("edges") %>% 
  left_join(nodeID2module %>% dplyr::rename(module="name"),by=c("from"="nodeID"))

graph %>% 
  ggraph(layout='circlepack')+
  geom_edge_link0(aes(color=module))+
  geom_node_point(aes(color=color,size=size))+
  geom_node_text(aes(label=label, size=size/4))+
  scale_edge_color_identity()+
  scale_color_identity()+
  scale_size_identity()+
  theme_graph()+
  theme(plot.margin= unit(c(0,0,0,0), "lines"),
        panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))
ggsave("../2.Output/12_Comparison_withGWAS/3.Graph/GwasGene_Module.png",
       width = 190, height = 190, units = "mm",bg="transparent")
###########################################
#processing for cytoscape
graph %N>% as_tibble() %>% 
  left_join(entrez2part_2dis %>% 
              mutate(Entrez_Gene_ID=Entrez_Gene_ID %>% as.character()) %>% 
              dplyr::select(Entrez_Gene_ID, SSc,SLE, DM, key_SSc_SLE,key_SSc_DM),by=c("name"="Entrez_Gene_ID")) %>% 
  write_csv("../2.Output/12_Comparison_withGWAS/3.Graph/1.NodeData4cytoscape.csv")

nodeID2name<-graph %N>% 
  as_tibble() %>% 
  dplyr::select(name) %>% 
  as.data.frame() %>% 
  rownames_to_column("nodeID") %>% 
  as_tibble() %>% 
  mutate(nodeID=nodeID %>% as.integer())

graph %E>%
  as_tibble() %>% 
  left_join(nodeID2name,by=c("from"="nodeID")) %>% 
  dplyr::select(-from) %>% rename(from=name)%>% 
  left_join(nodeID2name,by=c("to"="nodeID")) %>% 
  dplyr::select(-to) %>% rename(to=name) %>% 
  dplyr::select(from, to, everything())%>% 
  distinct(from,to,.keep_all = T) %>% 
  write_csv("../2.Output/12_Comparison_withGWAS/3.Graph/2.EdgeData4cytoscape.csv")

#SSc vs DM
graph %N>% 
  left_join(entrez2part_2dis %>% 
              mutate(Entrez_Gene_ID=Entrez_Gene_ID %>% as.character()) %>% 
              dplyr::select(Entrez_Gene_ID, SSc,SLE, DM, key_SSc_SLE,key_SSc_DM) %>% 
              distinct(Entrez_Gene_ID,.keep_all=T),
            by=c("name"="Entrez_Gene_ID")) %>% 
  filter((name %in% significant_module) | !is.na(key_SSc_DM)) %>% 
  as_tibble()  %>% 
  write_csv("../2.Output/12_Comparison_withGWAS/3.Graph/3.NodeData4cytoscape_SScDM.csv")

graph %N>% 
  left_join(entrez2part_2dis %>% 
              mutate(Entrez_Gene_ID=Entrez_Gene_ID %>% as.character()) %>% 
              dplyr::select(Entrez_Gene_ID, SSc,SLE, DM, key_SSc_SLE,key_SSc_DM) %>% 
              distinct(Entrez_Gene_ID,.keep_all=T),
            by=c("name"="Entrez_Gene_ID")) %>% 
  filter((name %in% significant_module) | !is.na(key_SSc_DM)) %E>%
  as_tibble() %>% 
  left_join(nodeID2name,by=c("from"="nodeID")) %>% 
  dplyr::select(-from) %>% rename(from=name)%>% 
  left_join(nodeID2name,by=c("to"="nodeID")) %>% 
  dplyr::select(-to) %>% rename(to=name) %>% 
  dplyr::select(from, to, everything())%>% 
  distinct(from,to,.keep_all = T) %>% 
  write_csv("../2.Output/12_Comparison_withGWAS/3.Graph/4.EdgeData4cytoscape_SScDM.csv")

#SSc vs SLE
graph %N>% 
  left_join(entrez2part_2dis %>% 
              mutate(Entrez_Gene_ID=Entrez_Gene_ID %>% as.character()) %>% 
              dplyr::select(Entrez_Gene_ID, SSc,SLE, DM, key_SSc_SLE,key_SSc_DM) %>% 
              distinct(Entrez_Gene_ID,.keep_all=T),
            by=c("name"="Entrez_Gene_ID")) %>% 
  filter((name %in% significant_module) | !is.na(key_SSc_SLE)) %>% 
  as_tibble()  %>% 
  write_csv("../2.Output/12_Comparison_withGWAS/3.Graph/5.NodeData4cytoscape_SScSLE.csv")

graph %N>% 
  left_join(entrez2part_2dis %>% 
              mutate(Entrez_Gene_ID=Entrez_Gene_ID %>% as.character()) %>% 
              dplyr::select(Entrez_Gene_ID, SSc,SLE, DM, key_SSc_SLE,key_SSc_DM) %>% 
              distinct(Entrez_Gene_ID,.keep_all=T),
            by=c("name"="Entrez_Gene_ID")) %>% 
  filter((name %in% significant_module) | !is.na(key_SSc_SLE)) %E>%
  as_tibble() %>% 
  left_join(nodeID2name,by=c("from"="nodeID")) %>% 
  dplyr::select(-from) %>% rename(from=name)%>% 
  left_join(nodeID2name,by=c("to"="nodeID")) %>% 
  dplyr::select(-to) %>% rename(to=name) %>% 
  dplyr::select(from, to, everything())%>% 
  distinct(from,to,.keep_all = T) %>% 
  write_csv("../2.Output/12_Comparison_withGWAS/3.Graph/6.EdgeData4cytoscape_SScSLE.csv")

##################







012_Retrieve_GWASgene.R





#install_github("ramiromagno/gwasrapidd")
library(gwasrapidd)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

###Example from GitHub######
studies <- get_studies(efo_trait = 'triple-negative breast cancer')
studies@studies[1:4]
############################

ssc_variant <- get_variants(efo_id = "EFO_0000717") #efo_trait= "systemic scleroderma"
ssc_variant@entrez_ids %>% 
  left_join(ssc_variant@variants ,by = "variant_id") %>% 
  write_csv("../2.Output/12_Comparison_withGWAS/1.Data/SSc_variant.csv")

sle_study <- get_studies(efo_id = "EFO_0002690")#efo_trait= "systemic lupus erythematosus"
sle_study@studies[1:4]

sle_variant <- get_variants(efo_id = "EFO_0002690")
sle_variant@entrez_ids %>% 
  left_join(sle_variant@variants ,by = "variant_id") %>% 
  write_csv("../2.Output/12_Comparison_withGWAS/1.Data/SLE_variant.csv")
sle_variant@variants[1:4]
sle_variant@genomic_contexts

dm_variant <- get_variants(efo_id = "EFO_0000398") #efo_trait= "dermatomyositis"
dm_variant@entrez_ids %>% 
  left_join(dm_variant@variants ,by = "variant_id") %>% 
  write_csv("../2.Output/12_Comparison_withGWAS/1.Data/DM_variant.csv")






012_SScvariant2Module.R





library(tidyverse)

SSc_variant<-read_csv("../2.Output/12_Comparison_withGWAS/1.Data/SSc_variant.csv")%>% 
  mutate(entrez_id=entrez_id %>% as.character())
geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  mutate(Entrez_Gene_ID=Entrez_Gene_ID %>% as.character())

SSc_variant %>% 
  left_join(geneid2modulecolor,by=c("entrez_id"="Entrez_Gene_ID")) %>% 
  write_csv("../2.Output/12_Comparison_withGWAS/1.Data/SSc_variant2module.csv")






012_VennDiagrams_part_components.R





library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library()

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

#GO enrichment result was saved at 
#"../2.Output/12_Comparison_withGWAS/1.Data/enrichGO_ALL.Rdata"
#"../2.Output/12_Comparison_withGWAS/1.Data/enrichGO_BP.Rdata"

#create dataframe for GO enrichment
SSc_variant<-read_csv("../2.Output/12_Comparison_withGWAS/1.Data/SSc_variant.csv")%>% 
  distinct(entrez_id,.keep_all = T)
DM_variant<-read_csv("../2.Output/12_Comparison_withGWAS/1.Data/DM_variant.csv")%>% 
  distinct(entrez_id,.keep_all = T)
SLE_variant<-read_csv("../2.Output/12_Comparison_withGWAS/1.Data/SLE_variant.csv")%>% 
  distinct(entrez_id,.keep_all = T)

entrez2part<-tibble(Entrez_Gene_ID= c(SSc_variant$entrez_id,SLE_variant$entrez_id,DM_variant$entrez_id) %>% unique()) %>% 
  mutate(SSc=if_else(Entrez_Gene_ID %in% SSc_variant$entrez_id, "SSc",""),
         SLE=if_else(Entrez_Gene_ID %in% SLE_variant$entrez_id, "SLE",""),
         DM=if_else(Entrez_Gene_ID %in% DM_variant$entrez_id, "DM",""),
         key=str_c(SSc,SLE,DM,sep = "-") %>% 
           str_replace_all(pattern = "^-*|-*$" ,replacement = "") %>% 
           str_replace_all(pattern = "--" ,replacement = "-")
           )
entrez2part %>% 
  group_by(key) %>% 
  summarise(count=n())


#do GO enrichment test
#ont = "All"
part2GOenrichment<-entrez2part %>% 
  group_by(key) %>%
  nest() %>% 
  mutate(GO=map(data,.f = function(data){
    enrichGO(data$Entrez_Gene_ID,org.Hs.eg.db,keyType = "ENTREZID",ont = "ALL")
    }))
save(part2GOenrichment, file = "../2.Output/12_Comparison_withGWAS/1.Data/enrichGO_ALL.Rdata")
map2(.x=part2GOenrichment$GO,.y=part2GOenrichment$key,function(GO, key){
  barplot(GO)+ ggtitle(key)
  ggsave(paste0("../2.Output/12_Comparison_withGWAS/2.Image/partwise_GO/All/",key,".png"))})

#ont = "BP"
part2GOenrichment_BP<-entrez2part %>% 
  group_by(key) %>%
  nest() %>% 
  mutate(GO=map(data,.f = function(data){
    enrichGO(data$Entrez_Gene_ID,org.Hs.eg.db,keyType = "ENTREZID",ont = "BP")
  }))
save(part2GOenrichment_BP, file = "../2.Output/12_Comparison_withGWAS/1.Data/enrichGO_BP.Rdata")
map2(.x=part2GOenrichment_BP$GO,.y=part2GOenrichment_BP$key,function(GO, key){
  barplot(GO)+ ggtitle(key)
  ggsave(paste0("../2.Output/12_Comparison_withGWAS/2.Image/partwise_GO/BP/",key,".png"))})

########################################################
#Only two diseases like SSc vs SLE, SSc vs DM
entrez2part_2dis<-tibble(Entrez_Gene_ID= c(SSc_variant$entrez_id,SLE_variant$entrez_id,DM_variant$entrez_id) %>% unique()) %>% 
  mutate(SSc=if_else(Entrez_Gene_ID %in% SSc_variant$entrez_id, "SSc",""),
         SLE=if_else(Entrez_Gene_ID %in% SLE_variant$entrez_id, "SLE",""),
         DM=if_else(Entrez_Gene_ID %in% DM_variant$entrez_id, "DM",""),
         key_SSc_SLE=str_c(SSc,SLE,sep = "-") %>% 
           str_replace_all(pattern = "^-*|-*$" ,replacement = ""),
         key_SSc_DM=str_c(SSc,DM,sep = "-") %>% 
           str_replace_all(pattern = "^-*|-*$" ,replacement = "")
  )

#SSc vs SLE All
part2GOenrichment_SScSLE<-entrez2part_2dis %>% 
  group_by(key_SSc_SLE) %>%
  nest() %>% 
  mutate(GO=map(data,.f = function(data){
    enrichGO(data$Entrez_Gene_ID,org.Hs.eg.db,keyType = "ENTREZID",ont = "ALL")
  }))
save(part2GOenrichment_SScSLE, file = "../2.Output/12_Comparison_withGWAS/1.Data/enrichGO_SScSLE_ALL.Rdata")
map2(.x=part2GOenrichment_SScSLE$GO,.y=part2GOenrichment_SScSLE$key_SSc_SLE,function(GO, key){
  barplot(GO)+ ggtitle(key)
  ggsave(paste0("../2.Output/12_Comparison_withGWAS/2.Image/partwise_GO/All_SScSLE/",key,".png"))})

#SSc vs SLE BP
part2GOenrichment_SScSLE<-entrez2part_2dis %>% 
  group_by(key_SSc_SLE) %>%
  nest() %>% 
  mutate(GO=map(data,.f = function(data){
    enrichGO(data$Entrez_Gene_ID,org.Hs.eg.db,keyType = "ENTREZID",ont = "BP")
  }))
save(part2GOenrichment_SScSLE, file = "../2.Output/12_Comparison_withGWAS/1.Data/enrichGO_SScSLE_BP.Rdata")
map2(.x=part2GOenrichment_SScSLE$GO,.y=part2GOenrichment_SScSLE$key_SSc_SLE,function(GO, key){
  barplot(GO)+ ggtitle(key)
  ggsave(paste0("../2.Output/12_Comparison_withGWAS/2.Image/partwise_GO/BP_SScSLE/",key,".png"))})


#SSc vs DM All
part2GOenrichment_SScDM<-entrez2part_2dis %>% 
  group_by(key_SSc_DM) %>%
  nest() %>% 
  mutate(GO=map(data,.f = function(data){
    enrichGO(data$Entrez_Gene_ID,org.Hs.eg.db,keyType = "ENTREZID",ont = "ALL")
  }))
save(part2GOenrichment_SScDM, file = "../2.Output/12_Comparison_withGWAS/1.Data/enrichGO_SScDM_ALL.Rdata")
map2(.x=part2GOenrichment_SScDM$GO,.y=part2GOenrichment_SScDM$key_SSc_DM,function(GO, key){
  barplot(GO)+ ggtitle(key)
  ggsave(paste0("../2.Output/12_Comparison_withGWAS/2.Image/partwise_GO/All_SScDM/",key,".png"))})

#SSc vs DM BP
part2GOenrichment_SScDM<-entrez2part_2dis %>% 
  group_by(key_SSc_DM) %>%
  nest() %>% 
  mutate(GO=map(data,.f = function(data){
    enrichGO(data$Entrez_Gene_ID,org.Hs.eg.db,keyType = "ENTREZID",ont = "BP")
  }))
save(part2GOenrichment_SScDM, file = "../2.Output/12_Comparison_withGWAS/1.Data/enrichGO_SScDM_BP.Rdata")
map2(.x=part2GOenrichment_SScDM$GO,.y=part2GOenrichment_SScDM$key_SSc_DM,function(GO, key){
  barplot(GO)+ ggtitle(key)
  ggsave(paste0("../2.Output/12_Comparison_withGWAS/2.Image/partwise_GO/BP_SScDM/",key,".png"))})

#########################################################
geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv")
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor


entrez2part %>% 
  left_join(geneid2modulecolor,by="Entrez_Gene_ID") %>% 
  filter(!is.na(module)) %>% 
  filter(module %in% significant_module) %>% 
  write_csv("../2.Output/12_Comparison_withGWAS/1.Data/Entrezid2partOfVenn.csv")


entrez2part_2dis%>% 
  left_join(geneid2modulecolor,by="Entrez_Gene_ID") %>% 
  filter(!is.na(module)) %>% 
  filter(module %in% significant_module) %>% 
  #dplyr::rename(from="module",to="Symbol") %>%
  write_csv("../2.Output/12_Comparison_withGWAS/1.Data/Entrezid2partOfVenn_2diseases.csv")






012_Vennduagram.R





library(VennDiagram)

c("#da6272","#0086ab","#777777","#f79646","#0086ab","#9bbb59","#bfbfbf")

Data_SSc<-read_csv("../2.Output/12_Comparison_withGWAS/1.Data/SSc_variant.csv") %>% 
  distinct(entrez_id) %>% unlist()
Data_SLE<-read_csv("../2.Output/12_Comparison_withGWAS/1.Data/SLE_variant.csv") %>% 
  distinct(entrez_id) %>% unlist()
Data_DM<-read_csv("../2.Output/12_Comparison_withGWAS/1.Data/DM_variant.csv") %>% 
  distinct(entrez_id) %>% unlist()

venn.diagram(
    x = list(
      Data_SSc ,
      Data_SLE,
      Data_DM
    ),
    category.names = c("SSc" , "SLE" , "DM"),
    filename = '../2.Output/12_Comparison_withGWAS/2.Image/GWASgene_VennPlot_ver2.png',
    output = TRUE ,
    imagetype="png" ,
    height = 1440 ,
    width = 1440 ,
    resolution = 900,
    compression = "lzw",
    lwd = 1,
    col=c("#0086ab","#9bbb59","#da6272"),
    fill = c(alpha('#0086ab',0.3),  alpha('#9bbb59',0.3), alpha("#da6272",0.3)),
    cex = 0.5,
    fontfamily = "Helvetica",
    cat.cex = 0.5,
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "Helvetica",
    cat.col = c( '#0086ab', '#9bbb59', "#da6272"),
    rotation = 1
  )
#################################################################
#SSc vs SLE

venn.diagram(
  x = list(
    Data_SSc ,
    Data_SLE
  ),
  category.names = c("SSc" , "SLE" ),
  filename = '../2.Output/12_Comparison_withGWAS/2.Image/GWASgene_VennPlot_SScSLE.png',
  output = TRUE ,
  imagetype="png" ,
  height = 2880 ,
  width = 2880 ,
  resolution = 1800,
  compression = "lzw",
  lwd = 1,
  col=c("#0086ab","#9bbb59"),
  fill = c(alpha('#0086ab',0.3),  alpha('#9bbb59',0.3)),
  cex = 0.5,
  fontfamily = "Helvetica",
  cat.cex = 0.5,
  cat.default.pos = "outer",
  cat.pos = c(27, -27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "Helvetica",
  cat.col = c( '#0086ab','#9bbb59'),
  #rotation = 1
  inverted=T
)



#################################################################
#SSc vs DM

venn.diagram(
  x = list(
    Data_SSc ,
    Data_DM
  ),
  category.names = c("SSc" , "DM" ),
  filename = '../2.Output/12_Comparison_withGWAS/2.Image/GWASgene_VennPlot_SScDM.png',
  output = TRUE ,
  imagetype="png" ,
  height = 1440 ,
  width = 1440 ,
  resolution = 900,
  compression = "lzw",
  lwd = 1,
  col=c("#0086ab","#da6272"),
  fill = c(alpha('#0086ab',0.3),  alpha("#da6272",0.3)),
  cex = 0.5,
  fontfamily = "Helvetica",
  cat.cex = 0.5,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "Helvetica",
  cat.col = c( '#0086ab',"#da6272"),
  #rotation = 1
  inverted=F
)






###sample from R Graph gallery##########
# Libraries
# library(tidyverse)
# library(hrbrthemes)
# library(tm)
# library(proustr)
# 
# # Load dataset from github
# data <- read.table("https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/14_SeveralIndepLists.csv", header=TRUE) 
# to_remove <- c("_|[0-9]|\\.|function|^id|script|var|div|null|typeof|opts|if|^r$|undefined|false|loaded|true|settimeout|eval|else|artist")
# data <- data %>% filter(!grepl(to_remove, word)) %>% filter(!word %in% stopwords('fr')) %>% filter(!word %in% proust_stopwords()$word)
# 
# # library
# library(VennDiagram)
# 
# #Make the plot
# venn.diagram(
#   x = list(
#     data %>% filter(artist=="booba") %>% select(word) %>% unlist() , 
#     data %>% filter(artist=="nekfeu") %>% select(word) %>% unlist() , 
#     data %>% filter(artist=="georges-brassens") %>% select(word) %>% unlist()
#   ),
#   category.names = c("Booba (1995)" , "Nekfeu (663)" , "Brassens (471)"),
#   filename = 'IMG/venn.png',
#   output = TRUE ,
#   imagetype="png" ,
#   height = 480 , 
#   width = 480 , 
#   resolution = 300,
#   compression = "lzw",
#   lwd = 1,
#   col=c("#440154ff", '#21908dff', '#fde725ff'),
#   fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),
#   cex = 0.5,
#   fontfamily = "sans",
#   cat.cex = 0.3,
#   cat.default.pos = "outer",
#   cat.pos = c(-27, 27, 135),
#   cat.dist = c(0.055, 0.055, 0.085),
#   cat.fontfamily = "sans",
#   cat.col = c("#440154ff", '#21908dff', '#fde725ff'),
#   rotation = 1
# )

#####################################






014_IPA_Barplot.R





library(tidyverse)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

plotbar_by_enrichresult<-function(enrichresult,savefile,barcolor="#0086ab"){
  topn=4
  enrichresult %>% 
    mutate(Description=Description %>% str_wrap(width=42)) %>% 
    separate(GeneRatio, into = c("numerator","denominator"),sep="/") %>% 
    mutate(numerator= numerator %>% as.numeric(),
           denominator = denominator %>% as.numeric(),
           GeneRatio= numerator / denominator )%>% 
    head(topn) %>% 
    mutate(color=barcolor) %>% 
    ggplot(aes(x=Description %>% reorder(-pvalue), y= -log10(qvalue)))+
    geom_bar(aes(fill=color),stat = "identity")+
    theme_minimal()+theme(panel.grid=element_blank())+
    theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
          panel.border = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
          #legend.box.background = element_rect(fill = "transparent"),
          #legend.position = c(1,0.1),
          #legend.justification = c(0,0),
          #axis.text=element_text(size=16),
          axis.text.y = element_text(size=character_size,color = "#5f5f5f",lineheight = 1.5),
          axis.title.x = element_text(size=character_size, color= "#5f5f5f"),
          plot.margin=unit(c(1,1,1,1),"cm"))+
    scale_fill_identity()+   #,"#9bbb59"
    scale_y_continuous(expand = c(0, 0))+
    coord_flip()+
    xlab("")
  ggsave(savefile,dpi = 320, width = 150, height = 100, units = "mm",  bg = "transparent")
}


IPR_result<-read_tsv("../2.Output/14_IPA_Result/2.Result/lightcyan.txt",skip = 2)

lightcyan_show<-c(1,2,7, 27)
barcolor<-"#0086ab"
barcolor=c(rep("#7fc2d5",3),rep("#0086ab",1))
IPR_result[lightcyan_show,] %>% 
  mutate(color=barcolor) %>% 
  ggplot(aes(x=`Ingenuity Canonical Pathways` %>% reorder(`-log(p-value)`),y=`-log(p-value)`))+
  geom_bar(aes(fill=color),stat = "identity")+
  theme_minimal()+theme(panel.grid=element_blank())+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent"),
        #legend.position = c(1,0.1),
        #legend.justification = c(0,0),
        #axis.text=element_text(size=16),
        axis.text.y = element_text(size=character_size,color = "#5f5f5f",lineheight = 1.5),
        axis.title.x = element_text(size=character_size, color= "#5f5f5f"),
        plot.margin=unit(c(1,1,1,1),"cm"))+
  scale_fill_identity()+   #,"#9bbb59"
  scale_y_continuous(expand = c(0, 0))+
  coord_flip()+
  xlab("")

ggsave("../2.Output/14_IPA_Result/Lightcyan/IPA_result_lightsyan_ver1.png",
       dpi = 320, width = 150, height = 150, units = "mm",  bg = "transparent")

IPR_result_turquoise<-read_tsv("../2.Output/14_IPA_Result/2.Result/turquoise_CP.txt",skip = 2)

show<-c(1,2,6,18)
barcolor<-"#0086ab"
barcolor=c(rep("#7fc2d5",2),rep("#0086ab",2))
IPR_result_turquoise[show,] %>% 
  mutate(color=barcolor) %>% 
  ggplot(aes(x=`Ingenuity Canonical Pathways` %>% reorder(`-log(p-value)`),y=`-log(p-value)`))+
  geom_bar(aes(fill=color),stat = "identity")+
  theme_minimal()+theme(panel.grid=element_blank())+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent"),
        #legend.position = c(1,0.1),
        #legend.justification = c(0,0),
        #axis.text=element_text(size=16),
        axis.text.y = element_text(size=character_size,color = "#5f5f5f",lineheight = 1.5),
        axis.title.x = element_text(size=character_size, color= "#5f5f5f"),
        plot.margin=unit(c(1,1,1,1),"cm"))+
  scale_fill_identity()+   #,"#9bbb59"
  scale_y_continuous(expand = c(0, 0))+
  coord_flip()+
  xlab("")

ggsave("../2.Output/14_IPA_Result/Turquoise/IPA_result_Turquoise_ver1.png",
       dpi = 320, width = 150, height = 150, units = "mm",  bg = "transparent")

IPR_result_lightyellow<-read_tsv("../2.Output/14_IPA_Result/2.Result/Lightyellow_CP.txt",skip = 2)

show<-c(1,2,3,4)
focus= c("Leukocyte Extravasation Signaling")
character_size=12
IPR_result_turquoise[show,] %>% 
  mutate(color=if_else(`Ingenuity Canonical Pathways` %in% focus,"#0086ab", "#7fc2d5" ),
         `Ingenuity Canonical Pathways`= `Ingenuity Canonical Pathways` %>% str_wrap(35)) %>% 
  ggplot(aes(x=`Ingenuity Canonical Pathways` %>% reorder(`-log(p-value)`),y=`-log(p-value)`))+
  geom_bar(aes(fill=color),stat = "identity")+
  theme_minimal()+theme(panel.grid=element_blank())+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent"),
        #legend.position = c(1,0.1),
        #legend.justification = c(0,0),
        #axis.text=element_text(size=16),
        axis.text.y = element_text(size=character_size,color = "#5f5f5f",lineheight = 1.5),
        axis.title.x = element_text(size=character_size, color= "#5f5f5f"),
        plot.margin=unit(c(1,1,1,1),"cm"))+
  scale_fill_identity()+   #,"#9bbb59"
  scale_y_continuous(expand = c(0, 0))+
  coord_flip()+
  xlab("")

ggsave("../2.Output/14_IPA_Result/Lightyellow/IPA_result_ver2.png",
       dpi = 320, width = 150, height = 60, units = "mm",  bg = "transparent")








014_IPA_Barplot_functionized.R





library(tidyverse)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")


plot_IPAresult<-function(IPA_result,show= c(1,2,3,4), focus="",character_size=12){0
  
  IPA_result[show,] %>% 
    mutate(color=if_else(`Ingenuity Canonical Pathways` %in% focus,"#0086ab", "#7fc2d5" ),
           `Ingenuity Canonical Pathways`= `Ingenuity Canonical Pathways` %>% str_wrap(35)) %>% 
    ggplot(aes(x=`Ingenuity Canonical Pathways` %>% reorder(`-log(p-value)`),y=`-log(p-value)`))+
    geom_bar(aes(fill=color),stat = "identity")+
    theme_minimal()+theme(panel.grid=element_blank())+
    theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
          panel.border = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
          #legend.box.background = element_rect(fill = "transparent"),
          #legend.position = c(1,0.1),
          #legend.justification = c(0,0),
          #axis.text=element_text(size=16),
          axis.text.y = element_text(size=character_size,color = "#5f5f5f",lineheight = 1.5),
          axis.title.x = element_text(size=character_size, color= "#5f5f5f"),
          plot.margin=unit(c(1,1,1,1),"cm"))+
    scale_fill_identity()+   #,"#9bbb59"
    scale_y_continuous(expand = c(0, 0))+
    coord_flip()+
    xlab("")
}



IPA_result<-read_tsv("../2.Output/14_IPA_Result/2.Result/Plum1_CP.txt",skip = 2)
show<-c(1,2,3,8)
focus= c("Leukocyte Extravasation Signaling")
plot_IPAresult(IPA_result,show=show, focus=focus)
ggsave("../2.Output/14_IPA_Result/Plum1/IPA_result_Plum1_ver1.png",
       dpi = 320, width = 150, height = 60, units = "mm",  bg = "transparent")

IPA_result<-read_tsv("../2.Output/14_IPA_Result/2.Result/Yellow_CP.txt",skip = 2)
show<-c(1,2,3,68)
focus= c("CXCR4 Signaling")
plot_IPAresult(IPA_result,show=show, focus=focus)
ggsave("../2.Output/14_IPA_Result/Plum1/IPA_result_yellow_ver1.png",
       dpi = 320, width = 150, height = 60, units = "mm",  bg = "transparent")



IPA_result<-read_tsv("../2.Output/14_IPA_Result/2.Result/lightcyan.txt",skip = 2)
show<-c(1,2,7, 27)
focus= c("Leukocyte Extravasation Signaling")
plot_IPAresult(IPA_result,show=show, focus=focus,character_size=16)
ggsave("../2.Output/14_IPA_Result/Lightcyan/IPA_result_lightcyan_ver2.png",
       dpi = 320, width = 180, height = 120, units = "mm",  bg = "transparent")






015_bnlearn_test.R





#https://fisproject.jp/2018/01/bayesian-network-in-r-introduction/

library(bnlearn)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

data(learning.test)
df <- learning.test
tibble::glimpse(df)

res <- gs(df)
plot(res, main = "the Grow-Shrink (GS)")

res2 <- iamb(df)
plot(res2, main = "the Incremental Association")

res3 <- hc(df)
plot(res3, main = "the hill-climbing (HC)")

fitted <- bn.fit(res3, data = df, method = "bayes")
#plot(fitted)
particles <- cpdist(fitted,
                    nodes = "A",
                    evidence = (B == "C"))


prop.table(table(particles))

cpquery(fitted,
        event = (A == "b"),
        evidence = (D == "c"))



library(bnlearn)

data(marks)

ug=empty.graph(names(marks))
plot(ug)
arcs(ug)=matrix(c("MECH","VECT"),ncol=2, byrow=T)
ug
dag=empty.graph(names(marks))
mat=matrix(c(0,1,1,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0),nrow=5,dimnames = list(nodes(dag),nodes(dag)))
amat(dag)=mat
plot(dag)
nbr(dag,"ANL")
chld=children(dag,"VECT")
par=parents(dag,"VECT")
o.par=base::sapply(chld,parents,x=dag)
unique(c(chld,par,o.par[o.par != "VECT"]))

mb(dag,"VECT")

score(dag,data=marks, type="loglik-g")
dag.eq=reverse.arc(dag, "STAT", "ANL")
plot(dag.eq)
score(dag.eq,data=marks, type="loglik-g")
all.equal(cpdag(dag),cpdag(dag.eq))

dag2=drop.arc(dag, from="STAT",to="ANL")
dag3=drop.arc(dag, from="ALG",to="VECT")
plot(dag2)
plot(dag3)
vstructs(dag2)
vstructs(dag3)

all.equal(cpdag(dag2),cpdag(dag3))

hl2=list(arcs=vstructs(dag2,arcs = T),lwd=4, col="black")
hl3=list(arcs=vstructs(dag3,arcs = T),lwd=4, col="black")
graphviz.plot(dag2, highlight = hl2, layout="fdp",main="dag2")
graphviz.plot(dag3, highlight = hl3, layout="fdp",main="dag2")
graphviz.plot(cpdag(dag2), highlight = hl2, layout="fdp",main="cpdag(dag2)")

bn.gs=gs(marks)
bn.hc=hc(marks)
plot(bn.hc)

bn.gs$arcs


spec = paste0("[PKC][PKA|PKC][praf|PKC:PKA][pmek|PKC:PKA:praf]",
             "[p44.42|pmek:PKA][pakts473|p44.42:PKA][P38|PKC:PKA]",
             "[pjnk|PKC:PKA][plcg][PIP3|plcg][PIP2|plcg:PIP3]")

net = bnlearn::model2network(spec)
class(net)

graphviz.plot(net, shape = "ellipse")

protein_dataset = read.table("../2.Output/15_BaysianNetwork_inModule/_bnlearn_test/sachs.data.txt",header = T)
print(hc(protein_dataset))
test=hc(protein_dataset)
plot(test)
spec="[PIP2][p44.42][PKC][PIP3|PIP2][pakts473|p44.42][pjnk|PKC][plcg|PIP3][PKA|p44.42:pakts473][P38|PKC:pjnk][pmek|plcg:P38][praf|pmek]"
net = bnlearn::model2network(spec)
class(net)
graphviz.plot(net, shape = "ellipse")
test$optimized















015_BN_inYellowModule.R





library(tidyverse)
library(tidygraph)
library(ggraph)
library(ReactomePA)
library(clusterProfiler)
library(bnlearn)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

lname<-load("../2.Output/01_Compare_deepsplit_signed/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/datExpr_filtered.Rdata")
#datExpr
module_<-read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv") %>% 
  filter(moduleColors_signed=="yellow") %>% 
  mutate(Entrez_Gene_ID = Entrez_Gene_ID %>% as.character())

reactome_result<-read_csv("../2.Output/02_signed_2_unmerged/logFC_ModuleColor.csv") %>% 
  filter(module_color=="yellow") %>% 
  .$Entrez_Gene_ID %>% 
  enrichPathway(pvalueCutoff=0.05, readable=F)
reactome_result@result %>% 
  write_csv("../2.Output/10_YellowModuleAnalysis/Reactome/Reactome_PathwayEnrichment_yellow_readableF.csv")


kegg_result<-read_csv("../2.Output/10_YellowModuleAnalysis/KEGG/KEGG_PathwayEnrichment_yellow.csv") %>% 
  head(3)
reactome_result<-read_csv("../2.Output/10_YellowModuleAnalysis/Reactome/Reactome_PathwayEnrichment_yellow_readableF.csv") %>% 
  head(3)


kegg<-kegg_result$geneID
names(kegg)<- kegg_result$ID
kegg<- map(kegg,.f=~str_split(.x, "\\/")[[1]])0

reactome<-reactome_result$geneID
names(reactome)<- reactome_result$Description 
reactome<- map(reactome,.f=~str_split(.x, "\\/")[[1]])
#######################
i=1
kegg[[i]]
names(kegg)[[i]]
for (i in seq_along(kegg)){
module_<-module_ %>% 
  mutate(!!names(kegg)[[i]]:= Entrez_Gene_ID %in% !!kegg[[i]])
}


module %>% 
  mutate(!!names(kegg)[[i]]:= Entrez_Gene_ID %in% !!kegg[[i]]) %>% 
  filter(!!as.name(names(kegg[[i]])))
  

datExpr %>% as_tibble() %>% 
  dplyr::select(one_of(module %>% filter(Entrez_Gene_ID %in% !!kegg[[i]]) %>% .$ProbeID)) %>% 
  #t() %>% 
  prcomp() %>% .$x %>% .[,1]
a$x[,1]

module %>% 
  filter(Entrez_Gene_ID %in% !!kegg[[i]]) %>% .$ProbeID

tibble(index=1:102) %>% 
  mutate(!!names(kegg)[[i]]:=datExpr %>% as_tibble() %>% 
           dplyr::select(one_of(module %>% filter(Entrez_Gene_ID %in% !!kegg[[i]]) %>% .$ProbeID)) %>% 
           #t() %>% 
           prcomp() %>% .$x %>% .[,1])
############################################

add_pc1<-function(df,groupname,geneid_list){
  df %>% 
    mutate(!!groupname:= datExpr %>% 
             as_tibble() %>% 
             dplyr::select(one_of(module %>% filter(Entrez_Gene_ID %in% !!geneid_list) %>% .$ProbeID)) %>% 
             prcomp() %>% .$x %>% .[,1]) %>% 
    return()
}
add_pc1(tibble(index=1:102), names(kegg)[[i]],  kegg[[i]])

basedf=tibble(index=1:102)
for (i in seq_along(kegg)){
  basedf<-basedf %>% 
    add_pc1(names(kegg)[[i]],  kegg[[i]])
}

for (i in seq_along(reactome)){
  basedf<-basedf %>% 
    add_pc1(names(reactome)[[i]],  reactome[[i]])
}
df<-basedf %>% 
  dplyr::select(-index)

res <- gs(df)
plot(res, main = "the Grow-Shrink (GS)")

cs.res<-hc(df)
plot(cs.res)


spec = paste0("[Mitophagy - animal][Translation|Mitophagy - animal]",
   "[DDX58/IFIH1-mediated induction of interferon-alpha/beta|Mitophagy - animal]",
   "[mRNA surveillance pathway|Mitophagy - animal:DDX58/IFIH1-mediated induction of interferon-alpha/beta]",
   "[Interleukin-1 family signaling|Translation:DDX58/IFIH1-mediated induction of interferon-alpha/beta]",
   "[RNA transport|mRNA surveillance pathway:Translation:Interleukin-1 family signaling]")

net = bnlearn::model2network(spec)
class(net)

graphviz.plot(net, shape = "ellipse")



graph<-tibble(from=bnlearn::arcs(net)[,"from"],
                 to=bnlearn::arcs(net)[,"to"]) %>% 
  as_tbl_graph(directed = TRUE)

graph %>% 
  ggraph(layout = "sugiyama")+
  geom_edge_link(arrow = arrow(length = unit(4, 'mm'),type = "open"), 
                                                end_cap = circle(10, 'mm'),color="#7fc2d5",
                                                width=0.9)+
  geom_node_point(size=20,color="#0086ab")+
  geom_node_text(aes(label=name))+
  theme_graph()+
  theme(plot.margin= unit(c(0,0,0,0), "lines"),
        panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))+
  xlim(c(-1.25,2))

layout<-create_layout(graph,layout = "sugiyama")
layout[layout$name=="Mitophagy - animal",c("x","y")]<-c(0.75,4)
layout[layout$name=="RNA transport",c("x","y")]<-c(0.75,1)
  
  
ggraph(layout)+
  geom_edge_link(arrow = arrow(length = unit(4, 'mm'),type = "open"),
                     aes(start_cap = label_rect(node1.witdh), 
                         end_cap = label_rect(node2.name)),
                 color="#7fc2d5",
                 width=0.9)+
  # geom_edge_diagonal(aes(start_cap = label_rect(node1.name), 
  #                        end_cap = label_rect(node2.name)),
  #                    arrow = arrow(length = unit(4, 'mm'),type = "open"),color="#7fc2d5"
  # )+
  #geom_node_point(size=20,color="#0086ab")+
  geom_node_label(aes(label=name %>% str_wrap(30)),
                  color="white",
                  size=4,
                  label.padding = unit(0.75, "lines"),
                  label.r = unit(0.5, "lines"),
                  #fill="#7fc2d5"
                  fill="#4bacc6")+
  theme_graph()+
  theme(plot.margin= unit(c(0,0,0,0), "lines"),
        panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))+
  xlim(c(-0.5,2))+
  ylim(c(0.5,4.5))
ggsave(filename = "../2.Output/15_BaysianNetwork_inModule/Yellow_BayessianNetPlot_ver2.png",
       width = 150, height = 100, units = "mm",bg="transparent")


graph %E>% 
  as_tibble() %>% 
  write_csv("../2.Output/15_BaysianNetwork_inModule/EdgeData.csv")

graph %N>% 
  as.data.frame() %>% 
  rename(label="name") %>% 
  rownames_to_column("nodeID") %>% 
  write_csv("../2.Output/15_BaysianNetwork_inModule/NodeData.csv")






016_Heatmap_and_Yellow_Barplot.R





library(tidyverse)
library(reshape2)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
c("#0086ab","#f79646","#9bbb59","#da6272","#777777","#bfbfbf")
#Create summarized data file 
# results<-list.files("../2.Output/16_Celltype_by_BaseSpace/Result/",full.names = T)
# 
# #Check whether Cell type data
# for (result in results){
#   if (read_csv(result,col_names = F)[[2,1]]!="Cell type by body systems"){
#     print(paste("warning!",result,"Not Cell type Data"))
#   }else{
#     print("OK")
#   }
# }
# 
# summary_df<-tibble()
# for (r in results){
#   print(r)
#   summary_df<-dplyr::bind_rows(summary_df,
#                                read_csv(r,skip = 3) %>% 
#                                  mutate(`P-Value`=`P-Value` %>% as.double(),
#                                         module=r %>% basename() %>% str_extract("_[[:alnum:]]{1,}\\.") %>% str_sub(2,-2))
#   )
# }
# 
# summary_df %>% 
#   write_csv("../2.Output/16_Celltype_by_BaseSpace/basespace_celltype_result_summary.csv")

summary_df<-read_csv("../2.Output/16_Celltype_by_BaseSpace/basespace_celltype_result_summary.csv")

summary_df %>% 
  #filter(`P-Value`<0.001) %>% 
  group_by(module) %>% 
  #top_n(n = -10,wt=`P-Value`) %>% 
  arrange(`P-Value`) %>% 
  group_map(~head(.x,20L),keep = T)

summary_df %>% 
  filter(`P-Value`<0.001) %>% 
  group_by(module) %>% 
  top_n(n = -10,wt=`P-Value`) %>% 
  arrange(`P-Value`) %>% 
  #arrange(Name)
  ggplot(aes(x=module,y=Name))+
  geom_tile(aes(fill=`P-Value`))


top_10_each_module<-summary_df %>% 
  filter(`P-Value`<0.001) %>% 
  group_by(module) %>% 
  top_n(n = -10,wt=`P-Value`) %>% 
  arrange(`P-Value`)
data4heatmap<-top_10_each_module %>% 
  dplyr::select(c(`Body System`,"Name",`P-Value`, "module")) %>% 
  unite(key,`Body System`,Name,sep=".") %>% 
  spread(key=module, value = `P-Value`,fill = 1) %>% 
  as.data.frame() %>% 
  column_to_rownames("key") %>% 
  as.matrix() 

df<-melt(data4heatmap)
colnames(df) <- c("celltype", "module", "P-Value")
head(df)

clr<--log10(data4heatmap)%>% 
  heatmap(scale = "none")

celltype.idx  <- rownames(data4heatmap)[clr$rowInd]
module.idx <- colnames(data4heatmap)[clr$colInd]

df$celltype  <- factor(df$celltype, levels = w)
df$module <- factor(df$module, levels = module.idx)

ggplot(df, aes(x=module, y=celltype, fill=-log10(`P-Value`)))+
  geom_tile()+
  theme_classic()+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 90))+
  scale_fill_gradient(low = "white", high = "red")
ggsave("../2.Output/16_Celltype_by_BaseSpace/celltype_heatmap_allmodule.png",
       height = 400,width = 400,units = "mm")

df %>% 
  mutate(key=celltype) %>% 
  separate(celltype,into = c("Body System", "Name"),sep = "\\.",remove = F) %>% .$Name

label<-celltype.idx %>% str_extract("\\..{1,}") %>% str_sub(2,-1)

df %>% 
  mutate(key=celltype) %>% 
  separate(celltype,into = c("Body System", "Name"),sep = "\\.",remove = F) %>% 
  ggplot(aes(x=module, y=celltype))+
  geom_point(aes(size=-log10(`P-Value`),color=`Body System`))+
  theme_classic()+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 90))#+ 
scale_y_discrete(labels= label )
#scale_color_gradient(low = "white", high = "red")
ggsave("../2.Output/16_Celltype_by_BaseSpace/celltype_heatmap_allmodule_ver2_exact.png",
       height = 400,width = 400,units = "mm")


df %>% 
  mutate(key=celltype,
         `P-Value`=replace(`P-Value`, `P-Value`==1,NA)) %>% 
  separate(celltype,into = c("Body System", "Name"),sep = "\\.",remove = F) %>% 
  ggplot(aes(x=module, y=celltype))+
  geom_point(aes(size=-log10(`P-Value`),color=`Body System`))+
  theme_classic()+
  #theme_linedraw()
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_line(color="grey"), # get rid of major grid
        #panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 90))+ 
  scale_y_discrete(labels= label )
ggsave("../2.Output/16_Celltype_by_BaseSpace/celltype_heatmap_allmodule_ver5.png",
       height = 400,width = 400,units = "mm")

############################################################
#significant modules

significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor
data4heatmap_s<-data4heatmap[,colnames(data4heatmap) %in% significant_module]
data4heatmap_s<-data4heatmap_s[rowSums(data4heatmap_s)!=ncol(data4heatmap_s),]


df_s<-melt(data4heatmap_s)
colnames(df_s) <- c("celltype", "module", "P-Value")
#df_s<-df_s %>% 
#  filter(`P-Value`<= 1e-10)
head(df_s)

clr<--log10(data4heatmap_s)%>% 
  heatmap(scale = "none")

celltype.idx  <- rownames(data4heatmap_s)[clr$rowInd]
module.idx <- colnames(data4heatmap_s)[clr$colInd]

df_s$celltype  <- factor(df_s$celltype, levels = celltype.idx)
df_s$module <- factor(df_s$module, levels = module.idx)


label<-celltype.idx %>% str_extract("\\..{1,}") %>% str_sub(2,-1)

df_s %>% 
  mutate(key=celltype,
         `P-Value`=replace(`P-Value`, `P-Value`==1,NA)) %>% 
  separate(celltype,into = c("Body System", "Name"),sep = "\\.",remove = F) %>% 
  ggplot(aes(x=module, y=celltype))+
  geom_point(aes(size=-log10(`P-Value`),color=`Body System`))+
  theme_classic()+
  #theme_linedraw()
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_line(color="grey"), # get rid of major grid
        #panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 90))+ 
  scale_y_discrete(labels= label )
ggsave("../2.Output/16_Celltype_by_BaseSpace/celltype_heatmap_significantmodule_ver3.png",
       height = 400,width = 400,units = "mm")


df_s %>% 
  mutate(key=celltype,
         `P-Value`=replace(`P-Value`, `P-Value`==1,NA)) %>% 
  separate(celltype,into = c("Body System", "Name"),sep = "\\.",remove = F) %>% 
  ggplot(aes(x=module, y=Name))+
  geom_point(aes(size=-log10(`P-Value`),color=`Body System`))+
  theme_classic()+
  #theme_linedraw()
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_line(color="grey"), # get rid of major grid
        #panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 90))+ 
  facet_wrap(facets = .~`Body System`, scales = "free",drop = T)
ggsave("../2.Output/16_Celltype_by_BaseSpace/celltype_heatmap_significantmodule_facet_ver4.png",
       height = 200,width = 500,units = "mm")

df_s %>% 
  mutate(`P-Value`=replace(`P-Value`, `P-Value`==1,NA),
         # module=module %>% as.character()
  ) %>% 
  separate(celltype,into = c("Body System", "Name"),sep = "\\.",remove = F) %>% 
  filter(!( `Body System` %in% c("Visual System"))) %>% 
  filter(!is.na(`P-Value`)) %>% 
  ggplot(aes(x=module, y=Name))+
  geom_point(aes(size=-log10(`P-Value`),color=`Body System`))+
  theme_classic()+
  #theme_linedraw()
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_line(color="grey"), # get rid of major grid
        #panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 90))+ 
  #scale_y_discrete(labels= label )+
  facet_wrap(facets = .~ `Body System`, scales = "free",drop = T,)
ggsave("../2.Output/16_Celltype_by_BaseSpace/celltype_heatmap_significantmodule_facet_ver5.png",
       height = 500,width = 500,units = "mm",limitsize = FALSE)

df_s %>% 
  mutate(`P-Value`=replace(`P-Value`, `P-Value`==1,NA),
         # module=module %>% as.character()
  ) %>% 
  separate(celltype,into = c("Body System", "Name"),sep = "\\.",remove = F) %>% 
  filter(( `Body System` %in% c("Immune System","Cardiovascular System","Integumentary System"))) %>% 
  filter(!is.na(`P-Value`)) %>% 
  ggplot(aes(x=module, y=Name %>% str_wrap(width=45)))+
  geom_point(aes(size=-log10(`P-Value`),color=`Body System`))+
  theme_classic()+
  #theme_linedraw()
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_line(color="grey"), # get rid of major grid
        #panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 90),
        legend.position  = "None")+ 
  ylab(label = "")+
  #scale_y_discrete(labels= label )+
  facet_wrap(facets = .~ `Body System`, scales = "free",drop = T,)
ggsave("../2.Output/16_Celltype_by_BaseSpace/celltype_heatmap_significantmodule_facet_ver6.png",
       height = 200,width = 400,units = "mm",limitsize = FALSE)

df_s %>% 
  mutate(`P-Value`=replace(`P-Value`, `P-Value`==1,NA),
         # module=module %>% as.character()
  ) %>% 
  separate(celltype,into = c("Body System", "Name"),sep = "\\.",remove = F) %>% 
  filter(( `Body System` %in% c("Immune System","Cardiovascular System","Integumentary System"))) %>% 
  #filter(`P-Value`<=1e-10) %>% 
  filter(!is.na(`P-Value`)) %>% 
  ggplot(aes(x=module, y=Name %>% str_wrap(width=45)))+
  geom_point(aes(size=-log10(`P-Value`),color=`Body System`))+
  theme_classic()+
  #theme_linedraw()
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_line(color="grey"), # get rid of major grid
        #panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 90),
        legend.position  = "None")+ 
  ylab(label = "")+
  #scale_y_discrete(labels= label )+
  facet_wrap(facets =  `Body System`~., scales = "free",drop = T,ncol=1)
ggsave("../2.Output/16_Celltype_by_BaseSpace/celltype_heatmap_significantmodule_facet_ver7.png",
       height = 300,width = 200,units = "mm",limitsize = FALSE)


df_s %>% 
  mutate(`P-Value`=replace(`P-Value`, `P-Value`==1,NA),
         # module=module %>% as.character()
  ) %>% 
  separate(celltype,into = c("Body System", "Name"),sep = "\\.",remove = F) %>% 
  filter(( `Body System` %in% c("Immune System","Cardiovascular System","Integumentary System"))) %>% 
  filter(!is.na(`P-Value`)) %>% 
  ggplot(aes(x=module, y=Name %>% str_wrap(width=45)))+
  geom_tile(aes(fill=-log10(`P-Value`)))+
  theme_classic()+
  #theme_linedraw()
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        #panel.grid.major = element_blank(),
        panel.grid.major =element_line(color="grey"), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 90),
        legend.position  = "None")+ 
  ylab(label = "")+
  #scale_y_discrete(labels= label )+
  scale_fill_gradient(low = "white",high="red")+
  facet_wrap(facets = .~ `Body System`, scales = "free",drop = T,)
ggsave("../2.Output/16_Celltype_by_BaseSpace/celltype_heatmap_significantmodule_facet_ver8.png",
       height = 200,width = 400,units = "mm",limitsize = FALSE)


df_s %>% 
  as_tibble() %>% 
  filter(module=="yellow"&`P-Value`<0.01) %>% 
  arrange(`P-Value`)%>% 
  separate(celltype,into = c("Body System", "Name"),sep = "\\.",remove = F) %>% 
  filter(( `Body System` %in% c("Immune System","Cardiovascular System","Integumentary System"))) %>% 
  mutate(colorcode=as.character(map(`Body System`,.f = function(x){
    system2color<-c("Immune System"="#f79646","Cardiovascular System"="#da6272","Integumentary System"="#9bbb59")
    return(system2color[[x]])
  }))) %>% 
  head(8) %>% 
  ggplot(aes(y=-log10(`P-Value`), x=reorder(Name %>% str_wrap(width=45),`P-Value`)))+
  geom_bar(aes(fill=colorcode),stat = "identity")+
  theme_minimal()+theme(panel.grid=element_blank())+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent"),
        #legend.position = c(1,0.1),
        #legend.justification = c(0,0),
        #axis.text=element_text(size=16),
        axis.text.x = element_text(size=12,color = "#5f5f5f",lineheight = 1.5,angle=30,hjust = 1),
        axis.title.x = element_text(size=12, color= "#5f5f5f"),
        plot.margin=unit(c(1,1,1,1),"cm"))+
  scale_fill_identity()+   #,"#9bbb59"
  scale_y_continuous(expand = c(0, 0))+
  #coord_flip()+
  xlab("")

df_s %>% 
  as_tibble() %>% 
  filter(module=="yellow"&`P-Value`<0.01) %>% 
  arrange(`P-Value`)%>% 
  separate(celltype,into = c("Body System", "Name"),sep = "\\.",remove = F) %>% 
  filter(( `Body System` %in% c("Immune System","Cardiovascular System","Integumentary System"))) %>% 
  mutate(colorcode=as.character(map(`Body System`,.f = function(x){
    system2color<-c("Immune System"="#f79646","Cardiovascular System"="#da6272","Integumentary System"="#9bbb59")
    return(system2color[[x]])
  }))) %>% 
  head(5) %>% 
  ggplot(aes(y=-log10(`P-Value`), x=reorder(Name %>% str_replace_all("of","of\n"),-`P-Value`)))+
  geom_bar(aes(fill=colorcode),stat = "identity")+
  theme_minimal()+theme(panel.grid=element_blank())+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent"),
        #legend.position = c(1,0.1),
        #legend.justification = c(0,0),
        #axis.text=element_text(size=16),
        axis.text.y = element_text(size=12,color = "#5f5f5f",lineheight =1.5),
        axis.text.x = element_text(size=12,color = "#5f5f5f",lineheight = 1.5),
        axis.title.x = element_text(size=12, color= "#5f5f5f"),
        plot.margin=unit(c(1,1,1,1),"cm"))+
  scale_fill_identity()+   #,"#9bbb59"
  scale_y_continuous(expand = c(0, 0))+
  xlab("")+
  ylab("-log10(P-Value)")+
  coord_flip()
ggsave("../2.Output/16_Celltype_by_BaseSpace/celltype_yellow_barplot_ver1.png",
       dpi = 320, width = 150, height = 120, units = "mm",  bg = "transparent")






016_Heatmap_test.R





library(tidyverse)
library(reshape2)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

results<-list.files("../2.Output/16_Celltype_by_BaseSpace/Result/",full.names = T)
read_csv(results[1],col_names = F)[[2,1]]=="Cell type by body systems"

#Check whether Cell type data
for (result in results){
  if (read_csv(result,col_names = F)[[2,1]]!="Cell type by body systems"){
    print(paste("warning!",result,"Not Cell type Data"))
  }else{
    print("OK")
  }
}

read_csv(results[1],skip = 3) %>% 
  mutate(module=results[1] %>% basename() %>% str_extract("_[[:alnum:]]{1,}\\.") %>% str_sub(2,-2))

summary_df<-tibble()
for (r in results){
  print(r)
  summary_df<-dplyr::bind_rows(summary_df,
    read_csv(r,skip = 3) %>% 
      mutate(`P-Value`=`P-Value` %>% as.double(),
             module=r %>% basename() %>% str_extract("_[[:alnum:]]{1,}\\.") %>% str_sub(2,-2))
  )
}

summary_df %>% 
  write_csv("../2.Output/16_Celltype_by_BaseSpace/basespace_celltype_result_summary.csv")


summary_df %>% 
  #filter(`P-Value`<0.001) %>% 
  group_by(module) %>% 
  #top_n(n = -10,wt=`P-Value`) %>% 
  arrange(`P-Value`) %>% 
  group_map(~head(.x,20L),keep = T)

summary_df %>% 
  filter(`P-Value`<0.001) %>% 
  group_by(module) %>% 
  top_n(n = -10,wt=`P-Value`) %>% 
  arrange(`P-Value`) %>% 
  #arrange(Name)
  ggplot(aes(x=module,y=Name))+
  geom_tile(aes(fill=`P-Value`))


top_10_each_module<-summary_df %>% 
  filter(`P-Value`<0.001) %>% 
  group_by(module) %>% 
  top_n(n = -10,wt=`P-Value`) %>% 
  arrange(`P-Value`)
data4heatmap<-top_10_each_module %>% 
  dplyr::select(c(`Body System`,"Name",`P-Value`, "module")) %>% 
  unite(key,`Body System`,Name,sep=".") %>% 
  spread(key=module, value = `P-Value`,fill = 1) %>% 
  as.data.frame() %>% 
  column_to_rownames("key") %>% 
  as.matrix() 

df<-melt(data4heatmap)
colnames(df) <- c("celltype", "module", "P-Value")
head(df)

clr<--log10(data4heatmap)%>% 
  heatmap(scale = "none")
  
celltype.idx  <- rownames(data4heatmap)[clr$rowInd]
module.idx <- colnames(data4heatmap)[clr$colInd]

df$celltype  <- factor(df$celltype, levels = w)
df$module <- factor(df$module, levels = module.idx)

ggplot(df, aes(x=module, y=celltype, fill=-log10(`P-Value`)))+
  geom_tile()+
  theme_classic()+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 90))+
  scale_fill_gradient(low = "white", high = "red")
ggsave("../2.Output/16_Celltype_by_BaseSpace/celltype_heatmap_allmodule.png",
       height = 400,width = 400,units = "mm")

df %>% 
  mutate(key=celltype) %>% 
  separate(celltype,into = c("Body System", "Name"),sep = "\\.",remove = F) %>% .$Name

label<-celltype.idx %>% str_extract("\\..{1,}") %>% str_sub(2,-1)

df %>% 
  mutate(key=celltype) %>% 
  separate(celltype,into = c("Body System", "Name"),sep = "\\.",remove = F) %>% 
  ggplot(aes(x=module, y=celltype))+
  geom_point(aes(size=-log10(`P-Value`),color=`Body System`))+
  theme_classic()+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 90))#+ 
  scale_y_discrete(labels= label )
  #scale_color_gradient(low = "white", high = "red")
ggsave("../2.Output/16_Celltype_by_BaseSpace/celltype_heatmap_allmodule_ver2_exact.png",
       height = 400,width = 400,units = "mm")


df %>% 
  mutate(key=celltype,
         `P-Value`=replace(`P-Value`, `P-Value`==1,NA)) %>% 
  separate(celltype,into = c("Body System", "Name"),sep = "\\.",remove = F) %>% 
  ggplot(aes(x=module, y=celltype))+
  geom_point(aes(size=-log10(`P-Value`),color=`Body System`))+
  theme_classic()+
  #theme_linedraw()
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_line(color="grey"), # get rid of major grid
        #panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 90))+ 
  scale_y_discrete(labels= label )
ggsave("../2.Output/16_Celltype_by_BaseSpace/celltype_heatmap_allmodule_ver5.png",
       height = 400,width = 400,units = "mm")

############################################################
#significant modules

significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor
data4heatmap_s<-data4heatmap[,colnames(data4heatmap) %in% significant_module]
data4heatmap_s<-data4heatmap_s[rowSums(data4heatmap_s)!=ncol(data4heatmap_s),]


df_s<-melt(data4heatmap_s)
colnames(df_s) <- c("celltype", "module", "P-Value")
#df_s<-df_s %>% 
#  filter(`P-Value`<= 1e-10)
head(df_s)

clr<--log10(data4heatmap_s)%>% 
  heatmap(scale = "none")

celltype.idx  <- rownames(data4heatmap_s)[clr$rowInd]
module.idx <- colnames(data4heatmap_s)[clr$colInd]

df_s$celltype  <- factor(df_s$celltype, levels = celltype.idx)
df_s$module <- factor(df_s$module, levels = module.idx)


label<-celltype.idx %>% str_extract("\\..{1,}") %>% str_sub(2,-1)

df_s %>% 
  mutate(key=celltype,
         `P-Value`=replace(`P-Value`, `P-Value`==1,NA)) %>% 
  separate(celltype,into = c("Body System", "Name"),sep = "\\.",remove = F) %>% 
  ggplot(aes(x=module, y=celltype))+
  geom_point(aes(size=-log10(`P-Value`),color=`Body System`))+
  theme_classic()+
  #theme_linedraw()
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_line(color="grey"), # get rid of major grid
        #panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 90))+ 
  scale_y_discrete(labels= label )
ggsave("../2.Output/16_Celltype_by_BaseSpace/celltype_heatmap_significantmodule_ver3.png",
       height = 400,width = 400,units = "mm")


df_s %>% 
  mutate(key=celltype,
         `P-Value`=replace(`P-Value`, `P-Value`==1,NA)) %>% 
  separate(celltype,into = c("Body System", "Name"),sep = "\\.",remove = F) %>% 
  ggplot(aes(x=module, y=Name))+
  geom_point(aes(size=-log10(`P-Value`),color=`Body System`))+
  theme_classic()+
  #theme_linedraw()
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_line(color="grey"), # get rid of major grid
        #panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 90))+ 
  facet_wrap(facets = .~`Body System`, scales = "free",drop = T)
ggsave("../2.Output/16_Celltype_by_BaseSpace/celltype_heatmap_significantmodule_facet_ver4.png",
       height = 200,width = 500,units = "mm")

df_s %>% 
  mutate(`P-Value`=replace(`P-Value`, `P-Value`==1,NA),
         # module=module %>% as.character()
         ) %>% 
  separate(celltype,into = c("Body System", "Name"),sep = "\\.",remove = F) %>% 
  filter(!( `Body System` %in% c("Visual System"))) %>% 
  filter(!is.na(`P-Value`)) %>% 
  ggplot(aes(x=module, y=Name))+
  geom_point(aes(size=-log10(`P-Value`),color=`Body System`))+
  theme_classic()+
  #theme_linedraw()
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_line(color="grey"), # get rid of major grid
        #panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 90))+ 
  #scale_y_discrete(labels= label )+
  facet_wrap(facets = .~ `Body System`, scales = "free",drop = T,)
ggsave("../2.Output/16_Celltype_by_BaseSpace/celltype_heatmap_significantmodule_facet_ver5.png",
       height = 500,width = 500,units = "mm",limitsize = FALSE)

df_s %>% 
  mutate(`P-Value`=replace(`P-Value`, `P-Value`==1,NA),
         # module=module %>% as.character()
  ) %>% 
  separate(celltype,into = c("Body System", "Name"),sep = "\\.",remove = F) %>% 
  filter(( `Body System` %in% c("Immune System","Cardiovascular System","Integumentary System"))) %>% 
  filter(!is.na(`P-Value`)) %>% 
  ggplot(aes(x=module, y=Name %>% str_wrap(width=45)))+
  geom_point(aes(size=-log10(`P-Value`),color=`Body System`))+
  theme_classic()+
  #theme_linedraw()
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_line(color="grey"), # get rid of major grid
        #panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 90),
        legend.position  = "None")+ 
  ylab(label = "")+
  #scale_y_discrete(labels= label )+
  facet_wrap(facets = .~ `Body System`, scales = "free",drop = T,)
ggsave("../2.Output/16_Celltype_by_BaseSpace/celltype_heatmap_significantmodule_facet_ver6.png",
       height = 200,width = 400,units = "mm",limitsize = FALSE)

df_s %>% 
  mutate(`P-Value`=replace(`P-Value`, `P-Value`==1,NA),
         # module=module %>% as.character()
  ) %>% 
  separate(celltype,into = c("Body System", "Name"),sep = "\\.",remove = F) %>% 
  filter(( `Body System` %in% c("Immune System","Cardiovascular System","Integumentary System"))) %>% 
  #filter(`P-Value`<=1e-10) %>% 
  filter(!is.na(`P-Value`)) %>% 
  ggplot(aes(x=module, y=Name %>% str_wrap(width=45)))+
  geom_point(aes(size=-log10(`P-Value`),color=`Body System`))+
  theme_classic()+
  #theme_linedraw()
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_line(color="grey"), # get rid of major grid
        #panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 90),
        legend.position  = "None")+ 
  ylab(label = "")+
  #scale_y_discrete(labels= label )+
  facet_wrap(facets =  `Body System`~., scales = "free",drop = T,ncol=1)
ggsave("../2.Output/16_Celltype_by_BaseSpace/celltype_heatmap_significantmodule_facet_ver7.png",
       height = 300,width = 200,units = "mm",limitsize = FALSE)


df_s %>% 
  mutate(`P-Value`=replace(`P-Value`, `P-Value`==1,NA),
         # module=module %>% as.character()
  ) %>% 
  separate(celltype,into = c("Body System", "Name"),sep = "\\.",remove = F) %>% 
  filter(( `Body System` %in% c("Immune System","Cardiovascular System","Integumentary System"))) %>% 
  filter(!is.na(`P-Value`)) %>% 
  ggplot(aes(x=module, y=Name %>% str_wrap(width=45)))+
  geom_tile(aes(fill=-log10(`P-Value`)))+
  theme_classic()+
  #theme_linedraw()
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        #panel.grid.major = element_blank(),
        panel.grid.major =element_line(color="grey"), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 90),
        legend.position  = "None")+ 
  ylab(label = "")+
  #scale_y_discrete(labels= label )+
  scale_fill_gradient(low = "white",high="red")+
  facet_wrap(facets = .~ `Body System`, scales = "free",drop = T,)
ggsave("../2.Output/16_Celltype_by_BaseSpace/celltype_heatmap_significantmodule_facet_ver8.png",
       height = 200,width = 400,units = "mm",limitsize = FALSE)






018_CalculateModulePreservation.R





library(tidyverse)
library(WGCNA)
library(limma)
library(patchwork)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

refset <- "GSE58095"
testset<-"GSE128314"

output_dir_path<-paste0("../2.Output/18_Comparison_withOtherDisease/Testset_",testset,"/")

datExprRef_color_Rdata_path <- paste0("../2.Output/07_ModulePreservation/Referenceset_GSE58095/datExprRef_color_",refset,".Rdata")
datExprTest_color_Rdata_path<-paste0(output_dir_path,"datExprTest_color_",testset,".Rdata")
modulepreservation_result_path<-paste0(output_dir_path,"modulePreservation",refset,"_", testset,".Rdata")
preservation_quality_table_path<-paste(output_dir_path,"Table_", refset, "_", testset, "_preservation_quality.csv", sep = "")

ZS_result_image_path<-paste0(output_dir_path,"modulePreservation_",refset,"_", testset, "_Zsummary.png")
MR_result_image_path<-paste0(output_dir_path,"modulePreservation_",refset,"_", testset, "_medianRank.png")
ZS_MR_result_image_path <- paste0(output_dir_path, "modulePreservation_",refset,"_", testset, "_Zsummary-medianRank.png")


significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor

color_value<-c("#da6272","#0086ab","#777777","#f79646","#0086ab","#9bbb59","#bfbfbf")
#######################################################################################
#module preservation ?v?Z
lname_Ref<-load(datExprRef_color_Rdata_path) #"datExprRef" "colorsRef"
lname_Test<-load(datExprTest_color_Rdata_path) #"datExprTest" "colorsTest" 

setLabels <- c(refset, testset)
multiExpr <- list(refset = list(data = datExprRef), testset = list(data = datExprTest))
multiColor <- list(refset = colorsRef)

#module preservation???Z?o
#????networkType??signed???w??
system.time({
  mp <- modulePreservation(multiExpr, multiColor,
                           networkType = "signed",
                           referenceNetworks = 1,
                           nPermutations = 200,
                           randomSeed = 1,
                           quickCor = 0,
                           verbose = 3)
})
save(mp,file = modulepreservation_result_path)
#######################################################################################
#visualize the result

load(modulepreservation_result_path)
#mp
ref <- 1
test <- 2
statsObs <- cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ <- cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])
#preservation medianRank??Zsummary???\??
#preservation??quality?????r
preservation_quality <- cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
                              signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))%>% 
  rownames_to_column("module") %>% 
  left_join(mp$preservation$Z[[ref]][[test]] %>% dplyr::select(moduleSize) %>% rownames_to_column("module"),by="module")
write_csv(preservation_quality,preservation_quality_table_path )

#Zsummary.pres
g_zs<-read_csv(preservation_quality_table_path) %>% 
  filter(!(module %in% c("gold","grey"))) %>% 
  ggplot(aes(x=moduleSize,y=Zsummary.pres))+
  #geom_rect(xmin=-Inf, xmax=Inf, ymin=10, ymax=Inf,alpha=0.8,fill="#da6272")+
  #geom_rect(xmin=-Inf, xmax=Inf, ymin=2, ymax=10,alpha=0.8,fill="#E9A1AA")+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = 2)+
  geom_hline(yintercept = 10)+
  geom_point(aes(color=module),size=5)+
  geom_text(aes(label=module),nudge_y = -0.4)+
  scale_color_identity()+
  scale_x_log10(limits=c(4,450))+
  labs(x="Module Size",y="Preservation Zsummary",title = "Preservation Zsummary")+
  theme_classic()+
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"))
g_zs
ggsave(g_zs,filename = ZS_result_image_path, 
       dpi=2000, width = 140, height = 140, units = "mm",  bg = "transparent")


#medianRank.pres
g_mr<-read_csv(preservation_quality_table_path) %>% 
  filter(!(module %in% c("gold","grey"))) %>% 
  ggplot(aes(x=moduleSize,y=medianRank.pres))+
  #geom_hline(yintercept = 2)+
  #geom_hline(yintercept = 10)+
  geom_point(aes(color=module),size=5)+
  geom_text(aes(label=module),nudge_y = -1)+
  scale_color_identity()+
  scale_x_log10(limits=c(4,450))+
  scale_y_reverse(limits=c(40,-0))+
  #ylim(-5,45)+
  labs(x="Module Size",y="Preservation Median Rank",title = "Preservation Median Rank")+
  theme_classic()+
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"))
g_mr
ggsave(g_mr,filename = MR_result_image_path, 
       dpi=2000, width = 140, height = 140, units = "mm",  bg = "transparent")

g_zs + g_mr
ggsave(g_zs + g_mr, filename =  ZS_MR_result_image_path,
       dpi=2000, width = 280, height = 140, units = "mm",  bg = "transparent")

#########################################################################################
#significant_module
#Zsummary.pres
g_sm_zs<-read_csv(preservation_quality_table_path) %>% 
  filter(!(module %in% c("gold","grey"))) %>% 
  filter(module %in% significant_module) %>% 
  ggplot(aes(x=moduleSize,y=Zsummary.pres))+
  #geom_rect(xmin=-Inf, xmax=Inf, ymin=10, ymax=Inf,alpha=0.8,fill="#da6272")+
  #geom_rect(xmin=-Inf, xmax=Inf, ymin=2, ymax=10,alpha=0.8,fill="#E9A1AA")+
  geom_hline(yintercept = 0)+
  geom_hline(yintercept = 2)+
  geom_hline(yintercept = 10)+
  geom_point(aes(color=module),size=5)+
  geom_text(aes(label=module),nudge_y = -0.4)+
  scale_color_identity()+
  scale_x_log10(limits=c(4,450))+
  labs(x="Module Size",y="Preservation Zsummary",title = "Preservation Zsummary")+
  theme_classic()+
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"))
g_sm_zs
ggsave(g_sm_zs,filename = paste0("../2.Output/07_ModulePreservation/modulePreservation_",refset,"_", testset,"_SignificantModule", "_Zsummary.png"), 
       dpi=2000, width = 140, height = 140, units = "mm",  bg = "transparent")


#medianRank.pres
g_sm_mr<-read_csv(preservation_quality_table_path) %>% 
  filter(!(module %in% c("gold","grey"))) %>% 
  filter(module %in% significant_module) %>% 
  ggplot(aes(x=moduleSize,y=medianRank.pres))+
  #geom_hline(yintercept = 2)+
  #geom_hline(yintercept = 10)+
  geom_point(aes(color=module),size=5)+
  geom_text(aes(label=module),nudge_y = -1)+
  scale_color_identity()+
  scale_x_log10(limits=c(4,450))+
  scale_y_reverse(limits=c(40,-0))+
  #ylim(-5,45)+
  labs(x="Module Size",y="Preservation Median Rank",title = "Preservation Median Rank")+
  theme_classic()+
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"))
g_sm_mr
ggsave(g_sm_mr,filename = paste0(output_dir_path,"modulePreservation_",refset,"_", testset,"_SignificantModule", "_medianRank.png"), 
       dpi=2000, width = 140, height = 140, units = "mm",  bg = "transparent")

g_sm_zs + g_sm_mr
#output_image_path<-paste0("../2.Output/07_ModulePreservation/modulePreservation_",refset,"_", testset, "_Zsummary-medianRank.png")
ggsave(g_sm_zs + g_sm_mr,filename = paste0(output_dir_path,"modulePreservation_",refset,"_", testset, "_SignificantModule","_Zsummary-medianRank.png"),
       dpi=2000, width = 280, height = 140, units = "mm",  bg = "transparent")







018_CalculateModulePreservation_functionized.R





library(tidyverse)
library(WGCNA)
library(limma)
library(patchwork)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

refset <- "GSE58095"
testset<-"GSE81071"

output_dir_path<-paste0("../2.Output/18_Comparison_withOtherDisease/Testset_",testset,"/")

datExprRef_color_Rdata_path <- paste0("../2.Output/07_ModulePreservation/Referenceset_GSE58095/datExprRef_color_",refset,".Rdata")
datExprTest_color_Rdata_path<-paste0(output_dir_path,"datExprTest_color_",testset,".Rdata")

calculate_modulepresearvation<-function(refset,testset,output_dir_path,datExprRef_color_Rdata_path,datExprTest_color_Rdata_path){
  modulepreservation_result_path<-paste0(output_dir_path,"modulePreservation",refset,"_", testset,".Rdata")
  preservation_quality_table_path<-paste(output_dir_path,"Table_", refset, "_", testset, "_preservation_quality.csv", sep = "")
  
  ZS_result_image_path<-paste0(output_dir_path,"modulePreservation_",refset,"_", testset, "_Zsummary.png")
  MR_result_image_path<-paste0(output_dir_path,"modulePreservation_",refset,"_", testset, "_medianRank.png")
  ZS_MR_result_image_path <- paste0(output_dir_path, "modulePreservation_",refset,"_", testset, "_Zsummary-medianRank.png")
  
  
  significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
    filter(pvalue<=0.01) %>% 
    .$ModuleColor
  
  color_value<-c("#da6272","#0086ab","#777777","#f79646","#0086ab","#9bbb59","#bfbfbf")
  #######################################################################################
  #module preservation ?v?Z
  lname_Ref<-load(datExprRef_color_Rdata_path) #"datExprRef" "colorsRef"
  lname_Test<-load(datExprTest_color_Rdata_path) #"datExprTest" "colorsTest" 
  
  setLabels <- c(refset, testset)
  multiExpr <- list(refset = list(data = datExprRef), testset = list(data = datExprTest))
  multiColor <- list(refset = colorsRef)
  
  #module preservation???Z?o
  #????networkType??signed???w??
  system.time({
    mp <- modulePreservation(multiExpr, multiColor,
                             networkType = "signed",
                             referenceNetworks = 1,
                             nPermutations = 200,
                             randomSeed = 1,
                             quickCor = 0,
                             verbose = 3)
  })
  save(mp,file = modulepreservation_result_path)
  #######################################################################################
  #visualize the result
  
  load(modulepreservation_result_path)
  #mp
  ref <- 1
  test <- 2
  statsObs <- cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
  statsZ <- cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])
  #preservation medianRank??Zsummary???\??
  #preservation??quality?????r
  preservation_quality <- cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
                                signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))%>% 
    rownames_to_column("module") %>% 
    left_join(mp$preservation$Z[[ref]][[test]] %>% dplyr::select(moduleSize) %>% rownames_to_column("module"),by="module")
  write_csv(preservation_quality,preservation_quality_table_path )
  
  #Zsummary.pres
  g_zs<-read_csv(preservation_quality_table_path) %>% 
    filter(!(module %in% c("gold","grey"))) %>% 
    ggplot(aes(x=moduleSize,y=Zsummary.pres))+
    #geom_rect(xmin=-Inf, xmax=Inf, ymin=10, ymax=Inf,alpha=0.8,fill="#da6272")+
    #geom_rect(xmin=-Inf, xmax=Inf, ymin=2, ymax=10,alpha=0.8,fill="#E9A1AA")+
    geom_hline(yintercept = 0)+
    geom_hline(yintercept = 2)+
    geom_hline(yintercept = 10)+
    geom_point(aes(color=module),size=5)+
    geom_text(aes(label=module),nudge_y = -0.4)+
    scale_color_identity()+
    scale_x_log10(limits=c(4,450))+
    labs(x="Module Size",y="Preservation Zsummary",title = "Preservation Zsummary")+
    theme_classic()+
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "transparent"))
  g_zs
  ggsave(g_zs,filename = ZS_result_image_path, 
         dpi=2000, width = 140, height = 140, units = "mm",  bg = "transparent")
  
  
  #medianRank.pres
  g_mr<-read_csv(preservation_quality_table_path) %>% 
    filter(!(module %in% c("gold","grey"))) %>% 
    ggplot(aes(x=moduleSize,y=medianRank.pres))+
    #geom_hline(yintercept = 2)+
    #geom_hline(yintercept = 10)+
    geom_point(aes(color=module),size=5)+
    geom_text(aes(label=module),nudge_y = -1)+
    scale_color_identity()+
    scale_x_log10(limits=c(4,450))+
    scale_y_reverse(limits=c(40,-0))+
    #ylim(-5,45)+
    labs(x="Module Size",y="Preservation Median Rank",title = "Preservation Median Rank")+
    theme_classic()+
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "transparent"))
  g_mr
  ggsave(g_mr,filename = MR_result_image_path, 
         dpi=2000, width = 140, height = 140, units = "mm",  bg = "transparent")
  
  g_zs + g_mr
  ggsave(g_zs + g_mr, filename =  ZS_MR_result_image_path,
         dpi=2000, width = 280, height = 140, units = "mm",  bg = "transparent")
  
  #########################################################################################
  #significant_module
  #Zsummary.pres
  g_sm_zs<-read_csv(preservation_quality_table_path) %>% 
    filter(!(module %in% c("gold","grey"))) %>% 
    filter(module %in% significant_module) %>% 
    ggplot(aes(x=moduleSize,y=Zsummary.pres))+
    #geom_rect(xmin=-Inf, xmax=Inf, ymin=10, ymax=Inf,alpha=0.8,fill="#da6272")+
    #geom_rect(xmin=-Inf, xmax=Inf, ymin=2, ymax=10,alpha=0.8,fill="#E9A1AA")+
    geom_hline(yintercept = 0)+
    geom_hline(yintercept = 2)+
    geom_hline(yintercept = 10)+
    geom_point(aes(color=module),size=5)+
    geom_text(aes(label=module),nudge_y = -0.4)+
    scale_color_identity()+
    scale_x_log10(limits=c(4,450))+
    labs(x="Module Size",y="Preservation Zsummary",title = "Preservation Zsummary")+
    theme_classic()+
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "transparent"))
  g_sm_zs
  ggsave(g_sm_zs,filename = paste0("../2.Output/07_ModulePreservation/modulePreservation_",refset,"_", testset,"_SignificantModule", "_Zsummary.png"), 
         dpi=2000, width = 140, height = 140, units = "mm",  bg = "transparent")
  
  
  #medianRank.pres
  g_sm_mr<-read_csv(preservation_quality_table_path) %>% 
    filter(!(module %in% c("gold","grey"))) %>% 
    filter(module %in% significant_module) %>% 
    ggplot(aes(x=moduleSize,y=medianRank.pres))+
    #geom_hline(yintercept = 2)+
    #geom_hline(yintercept = 10)+
    geom_point(aes(color=module),size=5)+
    geom_text(aes(label=module),nudge_y = -1)+
    scale_color_identity()+
    scale_x_log10(limits=c(4,450))+
    scale_y_reverse(limits=c(40,-0))+
    #ylim(-5,45)+
    labs(x="Module Size",y="Preservation Median Rank",title = "Preservation Median Rank")+
    theme_classic()+
    theme(panel.background = element_blank(),
          panel.border = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "transparent"))
  g_sm_mr
  ggsave(g_sm_mr,filename = paste0(output_dir_path,"modulePreservation_",refset,"_", testset,"_SignificantModule", "_medianRank.png"), 
         dpi=2000, width = 140, height = 140, units = "mm",  bg = "transparent")
  
  g_sm_zs + g_sm_mr
  #output_image_path<-paste0("../2.Output/07_ModulePreservation/modulePreservation_",refset,"_", testset, "_Zsummary-medianRank.png")
  ggsave(g_sm_zs + g_sm_mr,filename = paste0(output_dir_path,"modulePreservation_",refset,"_", testset, "_SignificantModule","_Zsummary-medianRank.png"),
         dpi=2000, width = 280, height = 140, units = "mm",  bg = "transparent")
  
}



refset <- "GSE58095"
testset<-"GSE128314"
calculate_modulepresearvation(refset,testset,output_dir_path,datExprRef_color_Rdata_path,datExprTest_color_Rdata_path)

#########################################################

refset <- "GSE58095"
testset<-"GSE32245"
suffix<- "_Dermatomyositis"

output_dir_path<-paste0("../2.Output/18_Comparison_withOtherDisease/Testset_",testset,suffix,"/")

datExprRef_color_Rdata_path <- paste0("../2.Output/07_ModulePreservation/Referenceset_GSE58095/datExprRef_color_",refset,".Rdata")
datExprTest_color_Rdata_path<-paste0(output_dir_path,"datExprTest_color_",testset,".Rdata")

# Error in .checkExpr(multiData, verbose, indent) : 
#   The submitted 'multiExpr' data contain genes or samples
# with zero variance or excessive counts of missing entries.
# Please use the function goodSamplesGenes on each set to identify the problematic
# genes and samples, and remove them before running modulePreservation. 
#NA??0?Œu???????H?????̂?????????Ǝv?????ǂƂ肠???????点???
lname_Test<-load(datExprTest_color_Rdata_path) #"datExprTest" "colorsTest" 
datExprTest[is.na(datExprTest)]<-0
datExprTest_color_Rdata_path<-paste0(output_dir_path,"datExprTest_color_",testset,"_NArepalaced0",".Rdata")

save(list = c("datExprTest","colorsTest"),file = datExprTest_color_Rdata_path)

calculate_modulepresearvation(refset,testset,output_dir_path,datExprRef_color_Rdata_path,datExprTest_color_Rdata_path)


#NA??1%?ȏ??܂ރT???v???????????đS17?T???v???ɂ??????̂Ń??g???C
refset <- "GSE58095"
testset<-"GSE32245"
suffix<- "_Dermatomyositis"

output_dir_path<-paste0("../2.Output/18_Comparison_withOtherDisease/Testset_",testset,suffix,"/")

datExprRef_color_Rdata_path <- paste0("../2.Output/07_ModulePreservation/Referenceset_GSE58095/datExprRef_color_",refset,".Rdata")

datExprTest_color_NAdeleted_file_path<-paste0( testset_dir_path,"datExprTest_color_",gseid,"_NAdeleted",".Rdata")

calculate_modulepresearvation(refset,testset,output_dir_path,datExprRef_color_Rdata_path,datExprTest_color_NAdeleted_file_path)


#############################################################################
refset <- "GSE58095"
testset<-"GSE81071"
suffix<- "_CLE"

output_dir_path<-paste0("../2.Output/18_Comparison_withMP/Testset_",testset,suffix,"/")


datExprRef_color_Rdata_path <- paste0("../2.Output/07_ModulePreservation/Referenceset_GSE58095/datExprRef_color_",refset,".Rdata")
datExprTest_color_Rdata_path<-paste0(output_dir_path,"datExprTest_color_",testset,".Rdata")

calculate_modulepresearvation(refset,testset,output_dir_path,datExprRef_color_Rdata_path,datExprTest_color_Rdata_path)


















018_testdataPreprocessing_GSE128314.R





library(tidyverse)
#library(WGCNA)
library(limma)
library(GEOquery)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

########################################################################################
#GSE128314
gseid<-"GSE128314"
series_matrix_path<-paste0("../1.Data/ModulePresearvation/GSE128314_Dermatomyositis/",gseid,"_series_matrix.txt.gz")
soft_file_path<-paste0("../1.Data/ModulePresearvation/GSE128314_Dermatomyositis/",gseid,"_family.soft.gz")
testset_dir_path<-paste0("../2.Output/18_Comparison_withOtherDisease/Testset_",gseid,"/")
datExprTest_color_file_path<-paste0( testset_dir_path,"datExprTest_color_",gseid,".Rdata")
if (!dir.exists(testset_dir_path)){
  dir.create(testset_dir_path)
}

geneid2modulecolor_path<-"../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv"
geneid2modulecolor<-read_csv(geneid2modulecolor_path) %>% 
  select(Entrez_Gene_ID,module) %>% 
  distinct(Entrez_Gene_ID,module)

geo_data <- getGEO(filename=soft_file_path,GSEMatrix = T)
show(geo_data)
#geo_data@gsms$GSM1104220@dataTable
geo_data@gpls$GPL17586

testset_annotation<-geo_data@gpls$GPL17586@dataTable@table

#testset_annotation<-testset_annotation %>% #colnames()
#  dplyr::select("ID","GENE") %>% 
#  left_join(geneid2modulecolor,by=c("GENE"="Entrez_Gene_ID")) %>% 
#  filter(!is.na(GENE) & !is.na(module))
#a<-testset_annotation %>% 
#  dplyr::select("gene_assignment") %>% 
#  separate(col = "gene_assignment",into =c("1","2","3","4","geneid","extra") ,
#           sep=" // ",extra = "warn",fill="warn") %>% 
#  separate(col = "geneid",into = c("Entrez_Gene_ID","other_id"),
#           sep=" /// ",extra = "warn",fill="warn")
testset_annotation<-testset_annotation %>%
  dplyr::select("ID","gene_assignment") %>% 
    separate(col = "gene_assignment",into =c("1","2","3","4","geneid","extra") ,
             sep=" // ",extra = "warn",fill="warn") %>% 
    separate(col = "geneid",into = c("Entrez_Gene_ID","other_id"),
             sep=" /// ",extra = "warn",fill="warn") %>% 
  dplyr::select("ID","Entrez_Gene_ID") %>% 
  mutate(Entrez_Gene_ID=Entrez_Gene_ID %>% as.numeric()) %>% 
  filter(Entrez_Gene_ID !="---") %>% 
  left_join(geneid2modulecolor,by="Entrez_Gene_ID") %>% 
  filter(!is.na(module))
  


#!Sample_data_processing
#Data for both channels were Lowess-normalized and then the log(2) ratio was taken
#read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45485/matrix/", comment="!")

##in series matrix annotation# 
#Data summarized and normalized using the RMA routine as implemented in the oligo R library
##in oligo userguide# https://www.bioconductor.org/packages/release/bioc/vignettes/oligo/inst/doc/oug.pdf
#Note that rma's output is in the log2 scale
SM<-read_tsv(series_matrix_path, comment="!")
SM%>% 
  gather(key="sample",value="value",-ID_REF) %>% 
  ggplot(aes(x=sample,y=value))+
  geom_boxplot()
ggsave(filename = paste0(testset_dir_path,"/Boxplot_",gseid,"_SM.png"))

consition=c(replicate("Dermatomyositis",n = 8),replicate("Healthy",n = 5))

# data4MAplot<-SM %>%
#  as.data.frame() %>%
#  column_to_rownames("ID_REF") %>%
#  as.matrix() %>%
#  t() %>%
#  as_tibble() %>%
#  mutate(condition=consition) %>%
#  group_by(condition) %>%
#  dplyr::summarise_each(dplyr::funs(mean)) %>%
#  as.data.frame() %>%
#  column_to_rownames("condition") %>%
#  as.matrix() %>%
#  t() %>%
#  as_tibble() %>%
#  rownames_to_column(var = "Entrez_Gene_ID") %>%
#  mutate(M=Dermatomyositis - Healthy,
#        A=(Dermatomyositis + Healthy)/2)


Dermatomyositis=c("GSM3671280", "GSM3671281", "GSM3671282", "GSM3671283", "GSM3671284", "GSM3671285", "GSM3671286", "GSM3671287")
Healthy=c("GSM3671288", "GSM3671289", "GSM3671290", "GSM3671291", "GSM3671292")

data4MAplot<-SM %>%
  mutate(Dermatomyositis = dplyr::select(.,one_of(Dermatomyositis)) %>% 
           rowMeans(na.rm = T),
         Healthy = dplyr::select(.,one_of(Healthy)) %>% 
           rowMeans(na.rm = T),
         M=Dermatomyositis - Healthy,
         A=(Dermatomyositis + Healthy)/2)


data4MAplot %>% 
  ggplot(aes(x=A,y=M))+
  geom_point()
ggsave(filename = paste0(testset_dir_path,"/MAplot_",gseid,"_SM.png"))



SM.q<-SM%>% 
  column_to_rownames(var = "ID_REF") %>% 
  as.matrix() %>% 
  normalizeBetweenArrays(method="quantile")

SM.q%>% as.data.frame() %>% 
  rownames_to_column("ID_REF") %>% 
  gather(key="sample",value="value",-ID_REF) %>% 
  ggplot(aes(x=sample,y=value))+
  geom_boxplot()
ggsave(filename = paste0(testset_dir_path,"/Boxplot_",gseid,"_SM_quantile.png"))

datExprTest_pre<-SM.q %>% as.data.frame() %>%
  rownames_to_column(var="ID") %>% 
  left_join(testset_annotation,by="ID") %>% 
  filter(!is.na(Entrez_Gene_ID))  %>% 
  filter(!is.na(module)) %>% 
  mutate(exprSignal_median=dplyr::select(.,starts_with("GSM")) %>%  apply( 1, median,na.rm=T)) %>% #select(exprSignal_median,everything())->a
  arrange(Entrez_Gene_ID,desc(exprSignal_median)) %>% 
  distinct(Entrez_Gene_ID,.keep_all = T) 

datExprTest<-datExprTest_pre%>% 
  column_to_rownames(var="Entrez_Gene_ID") %>% 
  select(starts_with("GSM"))%>% t()
dim(datExprTest)
datExprTest[1:10, 1:10]
colorsTest <- datExprTest_pre$module 

save(list = c("datExprTest","colorsTest"),file = datExprTest_color_file_path)

#MA plot?`??!!
datExprTest %>% rownames()
consition=c(replicate("Dermatomyositis",n = 8),replicate("Healthy",n = 5))
data4MAplot<-datExprTest %>% 
  as_tibble() %>% 
  mutate(condition=consition) %>% 
  group_by(condition) %>% 
  dplyr::summarise_each(dplyr::funs(mean)) %>% 
  as.data.frame() %>% 
  column_to_rownames("condition") %>% 
  as.matrix() %>% 
  t() %>% 
  as_tibble() %>% 
  rownames_to_column(var = "Entrez_Gene_ID") %>% 
  mutate(M=Dermatomyositis - Healthy,
         A=(Dermatomyositis + Healthy)/2)
data4MAplot %>% 
  ggplot(aes(x=A,y=M))+
  geom_point()
ggsave(filename = paste0(testset_dir_path,"/MAplot_",gseid,"_SM_quantile.png"))


########################################################################################





018_testdataPreprocessing_GSE32245_Dermatomyositis.R





rm(list=ls())

library(tidyverse)
#library(WGCNA)
library(limma)
library(GEOquery)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

########################################################################################
#GSE32245
gseid<-"GSE32245"
suffix<-"_Dermatomyositis"
#series_matrix_path<-paste0("../1.Data/ModulePresearvation/",gseid, suffix,"/",gseid,"_series_matrix.txt.gz")
series_matrix_path<-paste0("../1.Data/ModulePresearvation/",gseid, suffix,"/",gseid,"-GPL9258","_series_matrix.txt.gz")



soft_file_path<-paste0("../1.Data/ModulePresearvation/",gseid, suffix,"/",gseid,"_family.soft.gz")
testset_dir_path<-paste0("../2.Output/18_Comparison_withOtherDisease/Testset_",gseid, suffix,"/")
datExprTest_color_file_path<-paste0( testset_dir_path,"datExprTest_color_",gseid,".Rdata")
if (!dir.exists(testset_dir_path)){
  dir.create(testset_dir_path)
}

geneid2modulecolor_path<-"../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv"
geneid2modulecolor<-read_csv(geneid2modulecolor_path) %>% 
  select(Entrez_Gene_ID,module) %>% 
  distinct(Entrez_Gene_ID,module)

#Make ProbeID to Gene ID to Module color dataframe (testset_annotation)
#####################################################################

geo_data <- getGEO(filename=soft_file_path,GSEMatrix = T)
show(geo_data)
#geo_data@gsms$GSM1104220@dataTable
geo_data@gpls$GPL9258

testset_annotation<-geo_data@gpls$GPL9258@dataTable@table
Entrez_Gene_ID_exp<-"CompositeSequence BioSequence [GeneId]"

#testset_annotation<-testset_annotation %>% #colnames()
#  dplyr::select("ID","GENE") %>% 
#  left_join(geneid2modulecolor,by=c("GENE"="Entrez_Gene_ID")) %>% 
#  filter(!is.na(GENE) & !is.na(module))

# testset_annotation<-testset_annotation %>%
#   dplyr::select("ID","gene_assignment") %>% 
#   separate(col = "gene_assignment",into =c("1","2","3","4","geneid","extra") ,
#            sep=" // ",extra = "warn",fill="warn") %>% 
#   separate(col = "geneid",into = c("Entrez_Gene_ID","other_id"),
#            sep=" /// ",extra = "warn",fill="warn") %>% 
#   dplyr::select("ID","Entrez_Gene_ID") %>% 
#   mutate(Entrez_Gene_ID=Entrez_Gene_ID %>% as.numeric()) %>% 
#   filter(Entrez_Gene_ID !="---") %>% 
#   left_join(geneid2modulecolor,by="Entrez_Gene_ID") %>% 
#   filter(!is.na(module))


# testset_annotation<-testset_annotation %>% #colnames()
#    dplyr::select("ID",Entrez_Gene_ID_exp) %>%
#    left_join(geneid2modulecolor,by=c(Entrez_Gene_ID_exp="Entrez_Gene_ID")) %>%
#    filter(!is.na(GENE) & !is.na(module))

testset_annotation<-eval(parse(text = paste0('testset_annotation %>% 
  dplyr::select(\'ID\',Entrez_Gene_ID_exp) %>%
  left_join(geneid2modulecolor,by=c(\'',Entrez_Gene_ID_exp,'\'=\'Entrez_Gene_ID\')) %>%
  filter(!is.na(\`',Entrez_Gene_ID_exp,'\`) & !is.na(module))'))) %>% 
  mutate(ID=ID %>% as.character())%>% 
  dplyr::rename(Entrez_Gene_ID=Entrez_Gene_ID_exp)

# tmp<-testset_annotation %>% #colnames()
#   dplyr::select("ID",Entrez_Gene_ID_exp)
# 
# eval(parse(text = paste0('left_join(tmp,geneid2modulecolor,by=c(\'',Entrez_Gene_ID_exp,'\'=\'Entrez_Gene_ID\'))'))) %>% 
#   filter(!is.na(GENE) & !is.na(module))

#####################################################################

f<-file(series_matrix_path,"r")

sample_annotation<-list()
characteristic_id=1
for(i in 1:100){
  line<-readLines(con=f,1)
  elements<-str_split(line, "\t")[[1]]
  rowname<-elements[1]
  #print(rowname)
  if (rowname == "!Sample_characteristics_ch2"){
    #print("!!!!!")
    #print(line)
    sample_annotation[[paste0(rowname,"_",characteristic_id)]]= elements[2:length(elements)] %>% str_replace_all("\"","")
    characteristic_id<- characteristic_id + 1
  }
  if (rowname=="\"ID_REF\""){
    sample_annotation[["sample_ID"]]=elements[2:length(elements)] %>% str_replace_all("\"","")
    break()
  }
}
sample_annotation<-sample_annotation %>% 
  as_tibble()
Healthy <- sample_annotation %>% 
  filter(`!Sample_characteristics_ch2_4` == "disease state: Healthy") %>% .$sample_ID
Dermatomyositis<- sample_annotation %>% 
  filter(`!Sample_characteristics_ch2_4` != "disease state: Healthy") %>% .$sample_ID

#!Sample_data_processing
#Data for both channels were Lowess-normalized and then the log(2) ratio was taken
#read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45485/matrix/", comment="!")

##in series matrix annotation# 
#Data summarized and normalized using the RMA routine as implemented in the oligo R library
##in oligo userguide# https://www.bioconductor.org/packages/release/bioc/vignettes/oligo/inst/doc/oug.pdf
#Note that rma's output is in the log2 scale
SM<-read_tsv(series_matrix_path, comment="!") %>% 
  mutate_at(vars(starts_with("GSM")),as.double)
data4boxplot<-SM%>% 
  gather(key="sample",value="value",-ID_REF) 
data4boxplot%>% 
  ggplot(aes(x=sample,y=value))+
  geom_boxplot()
ggsave(filename = paste0(testset_dir_path,"/Boxplot_",gseid,"_SM.png"))


# data4MAplot<-SM %>%
#  as.data.frame() %>%
#  column_to_rownames("ID_REF") %>%
#  as.matrix() %>%
#  t() %>%
#  as_tibble() %>%
#  mutate(condition=consition) %>%
#  group_by(condition) %>%
#  dplyr::summarise_each(dplyr::funs(mean)) %>%
#  as.data.frame() %>%
#  column_to_rownames("condition") %>%
#  as.matrix() %>%
#  t() %>%
#  as_tibble() %>%
#  rownames_to_column(var = "Entrez_Gene_ID") %>%
#  mutate(M=Dermatomyositis - Healthy,
#        A=(Dermatomyositis + Healthy)/2)

data4MAplot<-SM %>%
  mutate(Dermatomyositis = dplyr::select(.,one_of(Dermatomyositis)) %>% 
           rowMeans(na.rm = T),
         Healthy = dplyr::select(.,one_of(Healthy)) %>% 
           rowMeans(na.rm = T),
         M=Dermatomyositis - Healthy,
         A=(Dermatomyositis + Healthy)/2)


data4MAplot %>% 
  ggplot(aes(x=A,y=M))+
  geom_point()
ggsave(filename = paste0(testset_dir_path,"/MAplot_",gseid,"_SM.png"))



SM.q<-SM%>% 
  column_to_rownames(var = "ID_REF") %>% 
  as.matrix() %>% 
  normalizeBetweenArrays(method="quantile")

SM.q%>% as.data.frame() %>% 
  rownames_to_column("ID_REF") %>% 
  gather(key="sample",value="value",-ID_REF) %>% 
  ggplot(aes(x=sample,y=value))+
  geom_boxplot()
ggsave(filename = paste0(testset_dir_path,"/Boxplot_",gseid,"_SM_quantile.png"))

datExprTest_pre<-SM.q %>% as.data.frame() %>%
  rownames_to_column(var="ID") %>% 
  left_join(testset_annotation,by="ID") %>% 
  filter(!is.na(Entrez_Gene_ID))  %>% 
  filter(!is.na(module)) %>% 
  mutate(exprSignal_median=dplyr::select(.,starts_with("GSM")) %>%  apply( 1, median,na.rm=T)) %>% #select(exprSignal_median,everything())->a
  arrange(Entrez_Gene_ID,desc(exprSignal_median)) %>% 
  distinct(Entrez_Gene_ID,.keep_all = T) 

datExprTest<-datExprTest_pre%>% 
  column_to_rownames(var="Entrez_Gene_ID") %>% 
  select(starts_with("GSM"))%>% t()
dim(datExprTest)
datExprTest[1:10, 1:10]
colorsTest <- datExprTest_pre$module 

save(list = c("datExprTest","colorsTest"),file = datExprTest_color_file_path)
########################################################################################
# Error in .checkExpr(multiData, verbose, indent) : 
#   The submitted 'multiExpr' data contain genes or samples
# with zero variance or excessive counts of missing entries.
# Please use the function goodSamplesGenes on each set to identify the problematic
# genes and samples, and remove them before running modulePreservation. 
#NA???????????炵???̂ł??̑Ή?
#?????1/3?ȏ?NA??Probe?͏??? ????NA??0?u???Ƃ??

#????????Null?????܂??ɂ??????T???v?????܂܂??Ă邩?炻?????????O?????ׂ????Ǝv??
#PCA???`???A?O???l???ۂ??T???v???????o???ď??O????
#GSM798993, GSM799022, 

lname<-load(datExprTest_color_file_path)
pca_ExprTest<-prcomp(x=datExprTest,scale. = F)
pca_ExprTest$x %>% 
  as.data.frame() %>% 
  rownames_to_column("sample_ID") %>% 
  ggplot(aes(x=PC1,y=PC2))+
  geom_point()+
  geom_text(aes(label=sample_ID))


datExprTest %>% 
  apply(MARGIN = 1,function(x){return(sum(is.na(x)))})
datExprTest %>% 
  apply(MARGIN = 1,function(x){return((sum(is.na(x))/length(x))*100)})
# GSM798989 GSM798992 GSM798993 GSM798994 GSM798995 GSM798998 GSM799000 
# 101        15      4019        19       546        54        38 
# GSM799002 GSM799003 GSM799006 GSM799011 GSM799012 GSM799017 GSM799019 
# 7        11        26        38        32       410        23 
# GSM799022 GSM799024 GSM799025 GSM799028 GSM799029 GSM799030 GSM799033 
# 1128        28       105       464        70        20        20 
# GSM799037 GSM799038 GSM799040 GSM799042 
# 17       408        15        26 
#1% ?ȏ?Null?̃T???v???????????ĉ??͂ɉ????Ȃ?
keep<-datExprTest %>% 
  apply(MARGIN = 1,function(x){return((sum(is.na(x))/length(x))*100)})<=1
kept<-datExprTest[keep,]
dim(kept)
#[1]   17 9231
kept %>% 
  apply(MARGIN = 2,function(x){any(is.na(x))}) %>% 
  sum()
#190
sample_annotation %>% 
  filter(!(sample_ID  %in% rownames( kept))) %>% 
  .$`!Sample_characteristics_ch2_4`
#removed sample disease state
# [1] "disease state: active"   "disease state: inactive"
# [3] "disease state: active"   "disease state: active"  
# [5] "disease state: inactive" "disease state: active"  
# [7] "disease state: active"   "disease state: Healthy" 

datExprTest<-datExprTest[keep,]
datExprTest[is.na(datExprTest)]<-0
datExprTest_color_NAdeleted_file_path<-paste0( testset_dir_path,"datExprTest_color_",gseid,"_NAdeleted",".Rdata")
save(list = c("datExprTest","colorsTest"),file = datExprTest_color_NAdeleted_file_path)

  
########################################################################################
# #MA plot?`??!!
# datExprTest %>% rownames()
# consition=c(replicate("Dermatomyositis",n = 8),replicate("Healthy",n = 5))
# data4MAplot<-datExprTest %>% 
#   as_tibble() %>% 
#   #mutate(condition=consition) %>% 
#   #left_join(sample_annotation)
#   group_by(condition) %>% 
#   dplyr::summarise_each(dplyr::funs(mean)) %>% 
#   as.data.frame() %>% 
#   column_to_rownames("condition") %>% 
#   as.matrix() %>% 
#   t() %>% 
#   as_tibble() %>% 
#   rownames_to_column(var = "Entrez_Gene_ID") %>% 
#   mutate(M=Dermatomyositis - Healthy,
#          A=(Dermatomyositis + Healthy)/2)
# data4MAplot %>% 
#   ggplot(aes(x=A,y=M))+
#   geom_point()
# ggsave(filename = paste0(testset_dir_path,"/MAplot_",gseid,"_SM_quantile.png"))


########################################################################################





018_testdataPreprocessing_GSE81071_CLE.R





rm(list=ls())

library(tidyverse)
#library(WGCNA)
library(limma)
library(GEOquery)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

########################################################################################
#GSE32245
gseid<-"GSE81071"
suffix<-"_CLE"
#series_matrix_path<-paste0("../1.Data/ModulePresearvation/",gseid, suffix,"/",gseid,"_series_matrix.txt.gz")
series_matrix_path<-paste0("../1.Data/ModulePresearvation/",gseid, suffix,"/",gseid,"_series_matrix.txt.gz")



soft_file_path<-paste0("../1.Data/ModulePresearvation/",gseid, suffix,"/",gseid,"_family.soft.gz")
testset_dir_path<-paste0("../2.Output/18_Comparison_withMP/Testset_",gseid, suffix,"/")
datExprTest_color_file_path<-paste0( testset_dir_path,"datExprTest_color_",gseid,".Rdata")
if (!dir.exists(testset_dir_path)){
  dir.create(testset_dir_path)
}

geneid2modulecolor_path<-"../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv"
geneid2modulecolor<-read_csv(geneid2modulecolor_path) %>% 
  select(Entrez_Gene_ID,module) %>% 
  distinct(Entrez_Gene_ID,module)

#Make ProbeID to Gene ID to Module color dataframe (testset_annotation)
#####################################################################

geo_data <- getGEO(filename=soft_file_path,GSEMatrix = T)
show(geo_data)
#geo_data@gsms$GSM1104220@dataTable
geo_data@gpls$GPL19983

testset_annotation<-geo_data@gpls$GPL19983@dataTable@table
Entrez_Gene_ID_exp<-"ENTREZ_GENE_ID"

#testset_annotation<-testset_annotation %>% #colnames()
#  dplyr::select("ID","GENE") %>% 
#  left_join(geneid2modulecolor,by=c("GENE"="Entrez_Gene_ID")) %>% 
#  filter(!is.na(GENE) & !is.na(module))

# testset_annotation<-testset_annotation %>%
#   dplyr::select("ID","gene_assignment") %>% 
#   separate(col = "gene_assignment",into =c("1","2","3","4","geneid","extra") ,
#            sep=" // ",extra = "warn",fill="warn") %>% 
#   separate(col = "geneid",into = c("Entrez_Gene_ID","other_id"),
#            sep=" /// ",extra = "warn",fill="warn") %>% 
#   dplyr::select("ID","Entrez_Gene_ID") %>% 
#   mutate(Entrez_Gene_ID=Entrez_Gene_ID %>% as.numeric()) %>% 
#   filter(Entrez_Gene_ID !="---") %>% 
#   left_join(geneid2modulecolor,by="Entrez_Gene_ID") %>% 
#   filter(!is.na(module))


# testset_annotation<-testset_annotation %>% #colnames()
#    dplyr::select("ID",Entrez_Gene_ID_exp) %>%
#    left_join(geneid2modulecolor,by=c(Entrez_Gene_ID_exp="Entrez_Gene_ID")) %>%
#    filter(!is.na(GENE) & !is.na(module))

testset_annotation<-eval(parse(text = paste0('testset_annotation %>% 
  dplyr::select(\'ID\',Entrez_Gene_ID_exp) %>%
  left_join(geneid2modulecolor,by=c(\'',Entrez_Gene_ID_exp,'\'=\'Entrez_Gene_ID\')) %>%
  filter(!is.na(\`',Entrez_Gene_ID_exp,'\`) & !is.na(module))'))) %>% 
  mutate(ID=ID %>% as.character())%>% 
  dplyr::rename(Entrez_Gene_ID=Entrez_Gene_ID_exp)

# tmp<-testset_annotation %>% #colnames()
#   dplyr::select("ID",Entrez_Gene_ID_exp)
# 
# eval(parse(text = paste0('left_join(tmp,geneid2modulecolor,by=c(\'',Entrez_Gene_ID_exp,'\'=\'Entrez_Gene_ID\'))'))) %>% 
#   filter(!is.na(GENE) & !is.na(module))

#####################################################################

f<-file(series_matrix_path,"r")

sample_annotation<-list()
characteristic_id=1
for(i in 1:100){
  line<-readLines(con=f,1)
  elements<-str_split(line, "\t")[[1]]
  rowname<-elements[1]
  #print(rowname)
  if (rowname == "!Sample_characteristics_ch1"){
    #print("!!!!!")
    #print(line)
    sample_annotation[[paste0(rowname,"_",characteristic_id)]]= elements[2:length(elements)] %>% str_replace_all("\"","")
    characteristic_id<- characteristic_id + 1
  }
  if (rowname=="\"ID_REF\""){
    sample_annotation[["sample_ID"]]=elements[2:length(elements)] %>% str_replace_all("\"","")
    break()
  }
}
sample_annotation<-sample_annotation %>% 
  as_tibble()
Healthy <- sample_annotation %>% 
  #filter((`!Sample_characteristics_ch1_2` == "disease state: healthy")| (`!Sample_characteristics_ch1_1` == "disease state: Normal")) %>% .$sample_ID
  filter((`!Sample_characteristics_ch1_2` == "disease state: healthy")) %>% .$sample_ID
  
CLE<- sample_annotation %>% 
  #filter((`!Sample_characteristics_ch1_2` == "disease state: sCLE")| (`!Sample_characteristics_ch1_1` == "disease state: SCLE")) %>% .$sample_ID
  filter((`!Sample_characteristics_ch1_2` == "disease state: sCLE")) %>% .$sample_ID

#!Sample_data_processing
#Data for both channels were Lowess-normalized and then the log(2) ratio was taken
#read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE45nnn/GSE45485/matrix/", comment="!")

##in series matrix annotation# 
#Data summarized and normalized using the RMA routine as implemented in the oligo R library
##in oligo userguide# https://www.bioconductor.org/packages/release/bioc/vignettes/oligo/inst/doc/oug.pdf
#Note that rma's output is in the log2 scale
SM<-read_tsv(series_matrix_path, comment="!") %>% 
  mutate_at(vars(starts_with("GSM")),as.double) %>% 
  dplyr::select(one_of(c("ID_REF",Healthy,CLE)))
data4boxplot<-SM%>% 
  gather(key="sample",value="value",-ID_REF) 
data4boxplot%>% 
  ggplot(aes(x=sample,y=value))+
  geom_boxplot()
ggsave(filename = paste0(testset_dir_path,"/Boxplot_",gseid,"_SM.png"))


# data4MAplot<-SM %>%
#  as.data.frame() %>%
#  column_to_rownames("ID_REF") %>%
#  as.matrix() %>%
#  t() %>%
#  as_tibble() %>%
#  mutate(condition=consition) %>%
#  group_by(condition) %>%
#  dplyr::summarise_each(dplyr::funs(mean)) %>%
#  as.data.frame() %>%
#  column_to_rownames("condition") %>%
#  as.matrix() %>%
#  t() %>%
#  as_tibble() %>%
#  rownames_to_column(var = "Entrez_Gene_ID") %>%
#  mutate(M=Dermatomyositis - Healthy,
#        A=(Dermatomyositis + Healthy)/2)

data4MAplot<-SM %>%
  mutate(CLE = dplyr::select(.,one_of(CLE)) %>% 
           rowMeans(na.rm = T),
         Healthy = dplyr::select(.,one_of(Healthy)) %>% 
           rowMeans(na.rm = T),
         M=CLE - Healthy,
         A=(CLE + Healthy)/2)


data4MAplot %>% 
  ggplot(aes(x=A,y=M))+
  geom_point()
ggsave(filename = paste0(testset_dir_path,"/MAplot_",gseid,"_SM.png"))



SM.q<-SM%>% 
  column_to_rownames(var = "ID_REF") %>% 
  as.matrix() %>% 
  normalizeBetweenArrays(method="quantile")

SM.q%>% as.data.frame() %>% 
  rownames_to_column("ID_REF") %>% 
  gather(key="sample",value="value",-ID_REF) %>% 
  ggplot(aes(x=sample,y=value))+
  geom_boxplot()
ggsave(filename = paste0(testset_dir_path,"/Boxplot_",gseid,"_SM_quantile.png"))

datExprTest_pre<-SM.q %>% as.data.frame() %>%
  rownames_to_column(var="ID") %>% 
  left_join(testset_annotation,by="ID") %>% 
  filter(!is.na(Entrez_Gene_ID))  %>% 
  filter(!is.na(module)) %>% 
  mutate(exprSignal_median=dplyr::select(.,starts_with("GSM")) %>%  apply( 1, median,na.rm=T)) %>% #select(exprSignal_median,everything())->a
  arrange(Entrez_Gene_ID,desc(exprSignal_median)) %>% 
  distinct(Entrez_Gene_ID,.keep_all = T) 

datExprTest<-datExprTest_pre%>% 
  column_to_rownames(var="Entrez_Gene_ID") %>% 
  select(starts_with("GSM"))%>% t()
dim(datExprTest)
datExprTest[1:10, 1:10]
colorsTest <- datExprTest_pre$module 

save(list = c("datExprTest","colorsTest"),file = datExprTest_color_file_path)
########################################################################################






019_gene_in_SignalingReceptorBinding.R





library(tidyverse)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

signaling_receptor_binding<-read_tsv("../2.Output/19_ImmuneCell-Fibroblast_Signaling/GeneOntology/GO0005102_signaling_receptor_binding.txt",col_names = F) %>% 
  dplyr::rename(UniProtKB=X1,Symbol=X2) %>% 
  .$Symbol
read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  mutate(in_signaling_receptor_binding= Symbol  %in% signaling_receptor_binding ) %>% 
  filter(in_signaling_receptor_binding) %>%
  filter(module != "grey") %>% 
  arrange(module) %>% 
  write_csv("../2.Output/19_ImmuneCell-Fibroblast_Signaling/GeneOntology/Gene_in_SignalingReceptorBinding.csv")
  
kWithin_top20_ofGroup<-read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv") %>% 
  group_by(moduleColors_signed) %>% 
  arrange(desc(kWithin)) %>% 
  top_n(n() * 0.2,wt = kWithin ) %>% .$Entrez_Gene_ID

read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  mutate(in_signaling_receptor_binding= Symbol  %in% signaling_receptor_binding ,
         kWithin_top20= as.character(Entrez_Gene_ID) %in% as.character(kWithin_top20_ofGroup)) %>% 
  filter(in_signaling_receptor_binding) %>%
  filter(module != "grey") %>% 
  arrange(module) %>% 
  write_csv("../2.Output/19_ImmuneCell-Fibroblast_Signaling/GeneOntology/Gene_in_SignalingReceptorBinding.csv")



#######################################
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor

read_csv("../2.Output/19_ImmuneCell-Fibroblast_Signaling/GeneOntology/Gene_in_SignalingReceptorBinding.csv") %>% 
  filter(module %in% significant_module) %>% 
  write_csv("../2.Output/19_ImmuneCell-Fibroblast_Signaling/GeneOntology/Gene_in_SignalingReceptorBinding_withUniprot_sigMod.csv")





#######################################







019_KeggPathwayMapping.R





library(tidyverse)
library(gplots)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
target_module<- c("yellow", "turquoise", "plum1")

geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv")

gene_in_pathway<-read_csv("../2.Output/19_ImmuneCell-Fibroblast_Signaling/LeukocyteTransendothelialMigration/Kegg_Pathway_Mapping/hsa04670.csv",col_names = F) %>% 
  mutate(Entrez_Gene_ID = X1 %>% str_remove("\\?\\?")) %>% 
  .$Entrez_Gene_ID
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor

gene_in_module<-geneid2modulecolor %>% 
  filter(module %in% target_module) %>% 
  dplyr::select(Entrez_Gene_ID, module) %>% 
  .$Entrez_Gene_ID 

gene_in_sigmodule<-geneid2modulecolor %>% 
  filter(module %in% significant_module) %>% 
  dplyr::select(Entrez_Gene_ID, module) %>% 
  .$Entrez_Gene_ID 
SSc_GWAS_gene<-read_csv("../2.Output/12_GWAS_compare2SLE_?畆?؉?/SSc_GWAS_Genes.csv")
  
tibble(Entrez_Gene_ID=gene_in_pathway,
       map_color="white") %>% 
  filter(!(Entrez_Gene_ID %in% gene_in_module)) %>% 
  bind_rows(geneid2modulecolor %>% 
              filter(module %in% target_module) %>% 
              dplyr::select(Entrez_Gene_ID, module) %>% 
              dplyr::rename(map_color=module) %>% 
              mutate(Entrez_Gene_ID = Entrez_Gene_ID%>% as.character())) %>% 
  write_csv("../2.Output/19_ImmuneCell-Fibroblast_Signaling/Plum1-Yellow_LeukocyteTransendothelialMigration/Kegg_Pathway_Mapping/kegg_pathway_mapper.csv")



geneid2modulecolor %>% 
  filter(Entrez_Gene_ID %in% gene_in_pathway,
         (module %in% significant_module) )%>% 
  dplyr::select(Entrez_Gene_ID,  module) %>% 
  write_csv("../2.Output/19_ImmuneCell-Fibroblast_Signaling/Plum1-Yellow_LeukocyteTransendothelialMigration/Kegg_Pathway_Mapping/kegg_pathway_mapper_sigmodule.csv")
#lightyellow ?Ƃ?Lightcyan?Ƃ??͐F???F???ł????炵??
geneid2modulecolor %>% 
  filter(Entrez_Gene_ID %in% gene_in_pathway,
         (module %in% significant_module) )%>% 
  mutate(color_code = col2hex(module)) %>% 
  dplyr::select(Entrez_Gene_ID, color_code, module, Symbol) %>% 
  write_csv("../2.Output/19_ImmuneCell-Fibroblast_Signaling/Plum1-Yellow_LeukocyteTransendothelialMigration/Kegg_Pathway_Mapping/kegg_pathway_mapper_sigmodule_namemod.csv")

tibble(Entrez_Gene_ID=gene_in_pathway,
       map_color="white") %>% 
  filter(!(Entrez_Gene_ID %in% gene_in_sigmodule)) %>% 
  bind_rows(geneid2modulecolor %>% 
              filter(module %in% significant_module) %>% 
              dplyr::select(Entrez_Gene_ID, module) %>% 
              dplyr::rename(map_color=module) %>% 
              mutate(Entrez_Gene_ID = Entrez_Gene_ID%>% as.character())) %>% 
  mutate(fill_color_code = col2hex(map_color),
         isGWASgene = Entrez_Gene_ID %in% SSc_GWAS_gene$GeneID,
         edge_color_code= if_else(isGWASgene, "red", "black"),
         fill_color_code_edge_color_code=str_c(fill_color_code,edge_color_code, sep=",")) %>% 
  dplyr::select(Entrez_Gene_ID,fill_color_code_edge_color_code) %>% 
  write_csv("../2.Output/19_ImmuneCell-Fibroblast_Signaling/kegg_pathway_mapper_sigmod_ver7.csv")








019_PathwayBar_Plum1.R





library(tidyverse)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

c("#777777","#f79646","#0086ab","#da6272","#9bbb59","#bfbfbf")

kegg_pathway<-read_csv("../2.Output/11_Plum1ModuleAnalysis/Kegg/KEGG_PathwayEnrichment_plum1.csv")
kegg_pathway %>% 
  head(5) %>% 
  ggplot(aes(x=Description %>% reorder(-qvalue), y= -log10(qvalue)))+
  geom_bar(stat = "identity",fill="plum1")+
  theme_minimal()+theme(panel.grid=element_blank())+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"),
        legend.margin=element_blank())+
  scale_y_continuous(expand = c(0, 0))+
  coord_flip()+
  xlab("")

kegg_pathway %>% 
  head(5) %>% 
  ggplot(aes(x=Description %>% reorder(-qvalue), y= -log10(qvalue)))+
  geom_bar(stat = "identity",fill="plum1")+
  theme_minimal()+theme(panel.grid=element_blank())+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))+
  scale_y_continuous(expand = c(0, 0))+
  coord_flip()+
  xlab("")


kegg_pathway %>% 
  separate(GeneRatio, into = c("numerator","denominator"),sep="/") %>% 
  mutate(numerator= numerator %>% as.numeric(),
         denominator = denominator %>% as.numeric(),
         GeneRatio= numerator / denominator )%>% 
  head(5) %>% 
  ggplot(aes(x=Description %>% reorder(-qvalue), y= -log10(qvalue)))+
  geom_bar(aes(fill=GeneRatio),stat = "identity")+
  theme_minimal()+theme(panel.grid=element_blank())+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent"),
        legend.position = c(1,0.1),
        legend.justification = c(1,0),
        axis.text=element_text(size=12))+
  scale_fill_gradient2(low = "#f79646",high = "#0086ab")+   #,"#9bbb59"
  scale_y_continuous(expand = c(0, 0))+
  coord_flip()+
  xlab("")

ggsave("../2.Output/19_ImmuneCell-Fibroblast_Signaling/Plum1-Yellow_LeukocyteTransendothelialMigration/plum1_keggPA.png",
       dpi = 320, width = 140, height = 110, units = "mm",  bg = "transparent")









019_PPI_Networl_yellow_centered.R





library(tidyverse)
library(tidygraph)
library(ggraph)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  dplyr::select(c("Entrez_Gene_ID", "Symbol",  "module")) %>% 
  filter(module %in% significant_module)
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor

interaction<-read_tsv("../2.Output/19_ImmuneCell-Fibroblast_Signaling/2.PPI_string/inSRB_sigMod/string_interactions.tsv") %>% 
  dplyr::select(c(`#node1`, "node2","node1_external_id", "node2_external_id")) %>% 
  dplyr::rename(node1=`#node1`)
mapping<-read_tsv("../2.Output/19_ImmuneCell-Fibroblast_Signaling/2.PPI_string/inSRB_sigMod/string_mapping.tsv")
stringid2modulecolor<-mapping %>% 
  left_join(geneid2modulecolor,by=c("queryItem"="Entrez_Gene_ID")) %>% 
  dplyr::select(c("queryItem", "stringId" , "Symbol", "module","preferredName")) %>% 
  dplyr::rename(Entrez_Gene_ID=queryItem) %>% 
  distinct(stringId, .keep_all=T )

centered_module="yellow"


PPI_graph<-interaction %>% 
  as_tbl_graph() %E>% 
  left_join(stringid2modulecolor %>% dplyr::rename_all(function(x){paste0(x,"_node1")}), by=c("node1_external_id"="stringId_node1")) %>% 
  left_join(stringid2modulecolor %>% dplyr::rename_all(function(x){paste0(x,"_node2")}), by=c("node1_external_id"="stringId_node2")) %N>%
  left_join(stringid2modulecolor, by=c("name"="preferredName")) 
  

tmp<-interaction %>% 
  left_join(stringid2modulecolor %>% dplyr::rename_all(function(x){paste0(x,"_node1")}), by=c("node1_external_id"="stringId_node1")) %>% 
  left_join(stringid2modulecolor %>% dplyr::rename_all(function(x){paste0(x,"_node2")}), by=c("node2_external_id"="stringId_node2")) %>% 
  filter((module_node1==centered_module)|(module_node2==centered_module))
yellow_contact_gene<-unique(c(tmp$node1,tmp$node2))



PPI_graph %>% 
  ggraph(layout = "graphopt")+
  geom_edge_link()+
  geom_node_point(aes(color=module),size=10)+
  geom_node_text(aes(label=name ))+
  scale_color_identity()
  




PPI_graph %>% 
  #activate("edges") %>% filter(module_node1 ==centered_module | module_node2 ==centered_module)%N>% 
  activate("nodes") %>% filter(name %in% yellow_contact_gene) %>% 
  ggraph(layout = "graphopt")+
  geom_edge_link()+
  geom_node_point(aes(color=module),size=10)+
  geom_node_text(aes(label=name ))+
  scale_color_identity()



centered_module_id<-PPI_graph %>% 
  #activate("edges") %>% filter(module_node1 ==centered_module | module_node2 ==centered_module)%N>% 
  activate("nodes") %>% filter(name %in% yellow_contact_gene) %N>% 
  as_data_frame() %>% 
  rownames_to_column("ID") %>% 
  filter(module==centered_module) %>% 
  .$name


full_graph<-create_complete(length(centered_module_id)) %N>% 
  mutate(name=centered_module_id) %E>% 
  mutate(visible=F)


set.seed(28)
PPI_graph %>% 
  #activate("edges") %>% filter(module_node1 ==centered_module | module_node2 ==centered_module)%N>% 
  activate("nodes") %>% filter(name %in% yellow_contact_gene) %>% 
  activate("edges") %>% mutate(visible=T) %>% 
  graph_join(full_graph,by = "name") %>% 
  activate("edges") %>% mutate(edge_color=if_else(visible, "grey66", "transparent")) %>% 
  ggraph(layout = "graphopt")+
  geom_edge_link0(aes(color=edge_color))+
  geom_node_point(aes(color=module),size=10)+
  geom_node_text(aes(label=name ))+
  scale_color_identity()+
  scale_edge_color_identity()+
  theme_graph()+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))


ggsave("../2.Output/19_ImmuneCell-Fibroblast_Signaling/Plum1-Yellow_LeukocyteTransendothelialMigration/PPI_Network_yellow_centered.png",
       height = 150,width = 200,dpi = 450,units = "mm",bg="transparent")

for (i in 1:100){
set.seed(i)
g<-PPI_graph %>% 
  #activate("edges") %>% filter(module_node1 ==centered_module | module_node2 ==centered_module)%N>% 
  activate("nodes") %>% filter(name %in% yellow_contact_gene) %>% 
  activate("edges") %>% mutate(visible=T) %>% 
  graph_join(full_graph,by = "name") %>% 
  activate("edges") %>% mutate(edge_color=if_else(visible, "grey", "transparent")) %>% 
  ggraph(layout = "graphopt")+
  geom_edge_link(aes(color=edge_color))+
  geom_node_point(aes(color=module),size=10)+
  geom_node_text(aes(label=name ))+
  scale_color_identity()+
  scale_edge_color_identity()+
  theme_graph()+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))+
  labs(title = i)
print(g)

}













019_process_string_interactionfile.R





library(tidyverse)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

mapping<-read_tsv("../2.Output/19_ImmuneCell-Fibroblast_Signaling/PPI_string/inSRB_sigMod/string_mapping.tsv")

interaction<-read_tsv("../2.Output/19_ImmuneCell-Fibroblast_Signaling/PPI_string/inSRB_sigMod/string_interactions_more.tsv")
nodes<-read_csv("../2.Output/19_ImmuneCell-Fibroblast_Signaling/PPI_string/inSRB_sigMod/string_interaction_nodes.csv")
geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv")

significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor

nodes %>% 
  left_join(mapping,by = c("name"="preferredName")) %>% 
  left_join(geneid2modulecolor %>% filter(module %in% significant_module), by = c("queryItem"="Entrez_Gene_ID")) %>% 
  dplyr::rename(Entrez_Gene_ID=queryItem) %>% 
  distinct(Entrez_Gene_ID,.keep_all = T) %>% 
  dplyr::select(-c("selected", "name", "shared name")) %>% 
  write_csv("../2.Output/19_ImmuneCell-Fibroblast_Signaling/PPI_string/inSRB_sigMod/string_interaction_nodes_withModCol.csv")






020_CalculatePhyper_PTPN11targetgeneInModule.R





library(tidyverse)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor
PTPN11_targets<-readxl::read_excel("../2.Output/20_PTPN11/PTPN11_targets/1-s2.0-S2211124718302195-mmc3.xlsx",sheet = 2) %>% 
  dplyr::select(Protein, `Student's t-test Significant PDGF_PDGF-Shp2i`)%>% 
  filter(!is.na(`Student's t-test Significant PDGF_PDGF-Shp2i`)) %>% 
  distinct(Protein) %>% 
  write_csv("../2.Output/20_PTPN11/PTPN11_targets/1.PTPN11_targets_filtered_uniprotID.csv",col_names = FALSE)

#-> https://www.uniprot.org/uploadlists/

#read_tsv("../2.Output/20_PTPN11/PTPN11_targets/2.UniprotID2EntrezGeneID.tab")

#-> http://refdic.rcai.riken.jp/tools/matchom.cgi

PTPN11_targets_humanEntrezID<-read_tsv("../2.Output/20_PTPN11/PTPN11_targets/3.MouseEntrezID2HumanEntrezID.txt", skip = 3)$`Homo sapiens`


#allgene_in_Analysis<-
#read_csv("../1.Data/GSE58095_main/Annotation/GeneAnnotation_GPL10558.csv")
geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  dplyr::select(Entrez_Gene_ID, Symbol, module)

geneid2modulecolor %>% 
  mutate(is_PTPN11target=Entrez_Gene_ID %in% PTPN11_targets_humanEntrezID) %>% 
  filter(is_PTPN11target) %>% .$module %>% 
  table()

nrow(geneid2modulecolor)# 12875
length(PTPN11_targets_humanEntrezID)# 1362

# hitInSample = num_PTPN11_target  
# hitInPop = length(PTPN11_targets_humanEntrezID) 
# failInPop = nrow(geneid2modulecolor) - hitInPop 
# sampleSize = genenum_module

phyper_result<-geneid2modulecolor%>% 
  mutate(is_PTPN11target=Entrez_Gene_ID %in% PTPN11_targets_humanEntrezID) %>% 
  group_by(module) %>% 
  nest() %>% 
  mutate(genenum_module=nrow(data[[1]]),
         num_PTPN11_target=nrow(data[[1]] %>% filter(is_PTPN11target)),
         p.value=phyper(num_PTPN11_target-1, 
                        length(PTPN11_targets_humanEntrezID),
                        nrow(geneid2modulecolor)-length(PTPN11_targets_humanEntrezID), 
                        genenum_module, lower.tail= FALSE)) %>% 
  arrange(p.value)
#  filter(p.value<=0.1)

phyper_result %>% 
  dplyr::select(-data) %>% 
  write_csv("../2.Output/20_PTPN11/PTPN11_targets/phyper_results.csv")


significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor
character_size=12
read_csv("../2.Output/20_PTPN11/3.PTPN11_targets/phyper_results.csv") %>% 
  filter(module %in% significant_module) %>% 
  head(5) %>% 
  ggplot(aes(x=module %>% reorder(-p.value), y=-log10(p.value)))+
  geom_bar(aes(fill=module),stat="identity")+
  theme_minimal()+theme(panel.grid=element_blank())+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent"),
        #legend.position = c(1,0.1),
        #legend.justification = c(0,0),
        #axis.text=element_text(size=16),
        axis.text.y = element_text(size=character_size,color = "#5f5f5f",lineheight = 1.5),
        axis.title.x = element_text(size=character_size, color= "#5f5f5f"),
        plot.margin=unit(c(1,1,1,1),"cm"))+
  scale_fill_identity()+   #,"#9bbb59"
  scale_y_continuous(expand = c(0, 0))+
  coord_flip()+
  xlab("")

ggsave("../2.Output/20_PTPN11/PTPN11_targets/phyper_barplot_ver1.png",
       dpi = 320, width = 150, height = 150, units = "mm",  bg = "transparent")

#vertical
character_size=12
read_csv("../2.Output/20_PTPN11/3.PTPN11_targets//phyper_results.csv") %>% 
  filter(module %in% significant_module) %>% 
  head(4) %>% 
  ggplot(aes(x=module %>% reorder(p.value), y=-log10(p.value)))+
  geom_bar(aes(fill=module),stat="identity")+
  geom_hline(yintercept = 2,color= "#da6272",alpha=0.7,size=1,linetype="dashed")+
  
  theme_minimal()+theme(panel.grid=element_blank())+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent"),
        #legend.position = c(1,0.1),
        #legend.justification = c(0,0),
        #axis.text=element_text(size=character_size,angle = 30,hjust = 1),
        axis.text.x = element_text(size=character_size,color = "#5f5f5f",lineheight = 1.5,angle = 30,hjust = 1),
        axis.title.y = element_text(size=character_size, color= "#5f5f5f"),
        plot.margin=unit(c(1,1,1,1),"cm"))+
  scale_fill_identity()+   #,"#9bbb59"
  scale_y_continuous(expand = c(0, 0))+
  xlab("")

ggsave("../2.Output/20_PTPN11/3.PTPN11_targets/phyper_barplot_ver3.png",
       dpi = 320, width = 80, height = 100, units = "mm",  bg = "transparent")

read_csv("../2.Output/20_PTPN11/PTPN11_targets/phyper_results.csv") %>% 
  filter(module %in% significant_module) %>% 
  head(5) %>% 
  ggplot(aes(x=module %>% reorder(p.value), y=-log10(p.value)))+
  geom_bar(aes(fill=module),stat="identity")+
  geom_hline(yintercept = 2,color="#F79646",size=1,alpha=0.5)+
  theme_minimal()+theme(panel.grid=element_blank())+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent"),
        #legend.position = c(1,0.1),
        #legend.justification = c(0,0),
        #axis.text=element_text(size=character_size,angle = 30,hjust = 1),
        axis.text.x = element_text(size=character_size,color = "#5f5f5f",lineheight = 1.5,angle = 30,hjust = 1),
        axis.title.y = element_text(size=character_size, color= "#5f5f5f"),
        plot.margin=unit(c(1,1,1,1),"cm"))+
  scale_fill_identity()+   #,"#9bbb59"
  scale_y_continuous(expand = c(0, 0))+
  xlab("")

ggsave("../2.Output/20_PTPN11/PTPN11_targets/phyper_barplot_ver3.png",
       dpi = 320, width = 150, height = 100, units = "mm",  bg = "transparent")


############################################
#http://mengnote.blogspot.com/2012/12/calculate-correct-hypergeometric-p.html
x=1; # x could be 0~5 
hitInSample = 5  # could be 0~5
hitInPop = 13 
failInPop = 52 - hitInPop 
sampleSize = 5
phyper(hitInSample-1, hitInPop, failInPop, sampleSize, lower.tail= FALSE)
phyper(hitInSample, hitInPop, failInPop, sampleSize, lower.tail= TRUE)






020_PTPN11_IPA_Metacore_REsult_bar.R





library(tidyverse)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")


IPR_result_turquoise<-read_tsv("../2.Output/20_PTPN11/PTPN11_targets/IPA_result/PTPN11target_CP.txt",skip = 2)
character_size=12
show<-c(1,2,3,5)
barcolor<-"#0086ab"
barcolor=c(rep(c("#7fc2d5","#0086ab"),2))#,rep("#0086ab",2))
IPR_result_turquoise[show,] %>% 
  mutate(color=barcolor,
         `Ingenuity Canonical Pathways` =`Ingenuity Canonical Pathways` %>% str_wrap(30)) %>% 
  ggplot(aes(x=`Ingenuity Canonical Pathways` %>% reorder(`-log(p-value)`),y=`-log(p-value)`))+
  geom_bar(aes(fill=color),stat = "identity")+
  theme_minimal()+theme(panel.grid=element_blank())+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent"),
        #legend.position = c(1,0.1),
        #legend.justification = c(0,0),
        #axis.text=element_text(size=16),
        axis.text.y = element_text(size=character_size,color = "#5f5f5f",lineheight = 1.5),
        axis.title.x = element_text(size=character_size, color= "#5f5f5f"),
        plot.margin=unit(c(1,1,1,1),"cm"))+
  scale_fill_identity()+   #,"#9bbb59"
  scale_y_continuous(expand = c(0, 0))+
  coord_flip()+
  xlab("")

ggsave("../2.Output/20_PTPN11/PTPN11_targets/IPA_result/IPA_result_ver1.png",
       dpi = 320, width = 150, height = 80, units = "mm",  bg = "transparent")


metacore_result_turquoise<-read_csv("../2.Output/20_PTPN11/PTPN11_targets/Metacore_result/Turquoise_PTPN11_target.csv",skip = 2)
character_size=8
show<-c(1,2,4,5)
barcolor<-"#0086ab"
#barcolor=c(rep(c("#7fc2d5","#0086ab"),2))#,rep("#0086ab",2))

focus=c("Oxidative stress_ROS-mediated MAPK activation via canonical pathways", "Development_Angiopoietin - Tie2 signaling")
metacore_result_turquoise[show,] %>% 
  mutate(color=if_else(Maps %in% focus,"#0086ab", "#7fc2d5" ),
         Maps =Maps %>% str_wrap(35)) %>% 
  ggplot(aes(x=Maps %>% reorder(-pValue),y=-log10(pValue)))+
  geom_bar(aes(fill=color),stat = "identity")+
  theme_minimal()+theme(panel.grid=element_blank())+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent"),
        #legend.position = c(1,0.1),
        #legend.justification = c(0,0),
        #axis.text=element_text(size=16),
        axis.text.y = element_text(size=character_size,color = "#5f5f5f",lineheight = 1.5),
        axis.title.x = element_text(size=character_size, color= "#5f5f5f"),
        plot.margin=unit(c(1,1,1,1),"cm"))+
  scale_fill_identity()+   #,"#9bbb59"
  scale_y_continuous(expand = c(0, 0))+
  coord_flip()+
  xlab("")

ggsave("../2.Output/20_PTPN11/PTPN11_targets/Metacore_result/Metacore_result_ver1.png",
       dpi = 320, width = 150, height = 80, units = "mm",  bg = "transparent")

metacore_result_turquoise<-read_csv("../2.Output/20_PTPN11/PTPN11_targets/Metacore_result/Turquoise_PTPN11_target.csv",skip = 2)
character_size=12
show<-c(1,2,4,5)
barcolor<-"#0086ab"
#barcolor=c(rep(c("#7fc2d5","#0086ab"),2))#,rep("#0086ab",2))

focus=c("ROS-mediated MAPK activation via canonical pathways", "Angiopoietin - Tie2 signaling")
metacore_result_turquoise[show,] %>% 
  separate(col=Maps,into = c("Category", "Maps"),sep="_") %>% 
  mutate(color=if_else(Maps %in% focus,"#0086ab", "#7fc2d5" ),
         Maps =Maps %>% str_wrap(35)) %>% 
  ggplot(aes(x=Maps %>% reorder(-pValue),y=-log10(pValue)))+
  geom_bar(aes(fill=color),stat = "identity")+
  theme_minimal()+theme(panel.grid=element_blank())+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent",colour ="transparent",size=0), # get rid of legend bg
        #legend.box.background = element_rect(fill = "transparent"),
        #legend.position = c(1,0.1),
        #legend.justification = c(0,0),
        #axis.text=element_text(size=16),
        axis.text.y = element_text(size=character_size,color = "#5f5f5f",lineheight = 1.5),
        axis.title.x = element_text(size=character_size, color= "#5f5f5f"),
        plot.margin=unit(c(1,1,1,1),"cm"))+
  scale_fill_identity()+   #,"#9bbb59"
  scale_y_continuous(expand = c(0, 0))+
  coord_flip()+
  xlab("")

ggsave("../2.Output/20_PTPN11/PTPN11_targets/Metacore_result/Metacore_result_ver2.png",
       dpi = 320, width = 150, height = 80, units = "mm",  bg = "transparent")






020_TOMnetwork_withPTPN11target.R





library(tidyverse)
library(WGCNA)
library(ggraph)
library(tidygraph)
library(gplots)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
lname_TOM<-load(file = "../2.Output/02_signed_2_unmerged/TOM.RData")
lname_signed2<-load("../2.Output/02_signed_2_unmerged/signed_2/networkConstruction_StepByStep_signed_2.Rdata")
lname_probes<-load("../2.Output/03_Network_Plot/probes.Rdata")
gene_annotation <- read_csv("../1.Data/GSE58095_main/Annotation/GeneAnnotation_GPL10558.csv")

kwithin_df<-read_csv("../2.Output/02_signed_2_unmerged/signed_2/Probe_kWitnin_color_signed_2.csv")
geneid2KwithinModulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv")
geneid2modulecolor<- read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv")
PTPN11_targets_humanEntrezID<-read_tsv("../2.Output/20_PTPN11/PTPN11_targets/3.MouseEntrezID2HumanEntrezID.txt", skip = 3)$`Homo sapiens`



modules<-"turquoise"
kWithin_top20_probe<-read_csv("../2.Output/02_signed_2_unmerged/geneid2KwithinModulecolor.csv") %>% 
  filter(moduleColors_signed==modules) %>% 
  mutate(kWithin_percent_rank=kWithin %>% percent_rank(),
         percent_rank_top20p=kWithin_percent_rank>0.8,
         kWithin_dense_rank=dense_rank(dplyr::desc(kWithin)),
         dense_rank_top20p=kWithin_dense_rank<=(max(kWithin_dense_rank)*0.2)) %>% 
  filter(dense_rank_top20p) %>% 
  dplyr::select(-percent_rank_top20p) %>% .$ProbeID



#TO_threshold<-0.05
TO_denserank_threshold<-0.95
modTOM<-TOM[moduleColors_signed %in% modules,moduleColors_signed %in% modules]
dimnames(modTOM)<-list(probes[moduleColors_signed %in% modules],probes[moduleColors_signed %in% modules])
modTOM[upper.tri(modTOM,diag = TRUE)]<-NA
t_graph<-modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>% 
  filter(!is.na(TO)) %>% 
  mutate(dense_rank=dense_rank(dplyr::desc(TO)))%>% 
  as_tbl_graph(directed=FALSE) %>% 
  activate("nodes") %>% 
  left_join(gene_annotation %>% dplyr::select(ID,Symbol,Entrez_Gene_ID),by=c("name"="ID")) %>% 
  left_join(kwithin_df %>% dplyr::select(ProbeID,moduleColors_signed,kWithin),by=c("name"="ProbeID")) %>% 
  mutate(dense_rank=dense_rank(dplyr::desc(kWithin)) ) %>% 
  filter(dense_rank <= max(dense_rank)*0.2)%>% 
  activate("edges") %>%
  filter(dense_rank <= max(dense_rank)*0.01) 

t_graph%>% 
  ggraph(layout = "kk")+
  geom_edge_link(aes(),colour=modules,width=1)+
  geom_node_point(aes(color=moduleColors_signed),size=10)+
  geom_node_text(aes(label=Symbol),color="black",size=2)+
  theme_graph()+
  scale_color_identity()+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))

ggsave("../2.Output/20_PTPN11/TOM_Network_Turquoise.png",
       height = 200,width = 200,dpi = 450,units = "mm",bg="transparent")


t_graph %>% 
  activate("nodes") %>% 
  mutate(is_PTPN11target= Entrez_Gene_ID %in% PTPN11_targets_humanEntrezID,
         color= if_else(is_PTPN11target, "#da5019" , moduleColors_signed))%>% 
  ggraph(layout = "kk")+
  geom_edge_link(aes(),colour=modules,width=1)+
  geom_node_point(aes(color=color),size=10)+
  geom_node_text(aes(label=Symbol),color="black",size=2)+
  theme_graph()+
  scale_color_identity()+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))

ggsave("../2.Output/20_PTPN11/TOM_Network_Turquoise_withPTPN11target.png",
       height = 200,width = 200,dpi = 450,units = "mm",bg="transparent")

t_graph %>% 
  activate("nodes") %>% 
  mutate(is_PTPN11target= Entrez_Gene_ID %in% PTPN11_targets_humanEntrezID,
         color= if_else(is_PTPN11target, "#da5019" , moduleColors_signed))%>%
  activate("nodes") %>% 
  as.data.frame() %>% 
  rownames_to_column("nodeID") %>% 
  rename(probeID=name) %>% 
  write_csv("../2.Output/20_PTPN11/TOMturquoise_NodeData4sytescape.csv")
t_graph %>% 
  activate("nodes") %>% 
  mutate(is_PTPN11target= Entrez_Gene_ID %in% PTPN11_targets_humanEntrezID,
         color= if_else(is_PTPN11target, "#da5019" , moduleColors_signed))%>%
  activate("edges") %>% 
  as_tibble() %>% 
  mutate(color="Turquoise") %>% 
  write_csv("../2.Output/20_PTPN11/TOMturquoise_EdgeData4sytescape.csv")

###########################################################################

kwithin2isPTPN11target<-geneid2KwithinModulecolor %>% 
  left_join(geneid2modulecolor %>% dplyr::select(-Entrez_Gene_ID), by= c("ProbeID"="probe")) %>% 
  mutate(is_PTPN11target=Entrez_Gene_ID %in% PTPN11_targets_humanEntrezID) %>% 
  filter(moduleColors_signed=="turquoise") %>% 
  arrange(desc(kWithin))

kWithin_PTPN11target<-kwithin2isPTPN11target %>% 
  filter(is_PTPN11target) %>% 
  .$kWithin
kWithin_NotPTPN11target<-kwithin2isPTPN11target %>% 
  filter(!is_PTPN11target) %>% 
  .$kWithin


t.test(kWithin_PTPN11target, kWithin_NotPTPN11target, var.equal=F)

kwithin2isPTPN11target %>% 
  filter(is_PTPN11target) %>% 
  write_tsv("../2.Output/20_PTPN11/Turquoise_PTPN11targetgenes.tsv")






020_Turquoise_Centered_PPI_Network.R





library(tidyverse)
library(tidygraph)
library(ggraph)
library(gplots)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  dplyr::select(c("Entrez_Gene_ID", "Symbol",  "module")) %>% 
  filter(module %in% significant_module)
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor

interaction<-read_tsv("../2.Output/19_ImmuneCell-Fibroblast_Signaling/2.PPI_string/inSRB_sigMod/string_interactions.tsv") %>% 
  dplyr::select(c(`#node1`, "node2","node1_external_id", "node2_external_id")) %>% 
  dplyr::rename(node1=`#node1`)
mapping<-read_tsv("../2.Output/19_ImmuneCell-Fibroblast_Signaling/2.PPI_string/inSRB_sigMod/string_mapping.tsv")
stringid2modulecolor<-mapping %>% 
  left_join(geneid2modulecolor,by=c("queryItem"="Entrez_Gene_ID")) %>% 
  dplyr::select(c("queryItem", "stringId" , "module","preferredName")) %>% 
  dplyr::rename(Entrez_Gene_ID=queryItem) %>% 
  distinct(stringId, .keep_all=T )

centered_module="turquoise"


PPI_graph<-interaction %>% 
  as_tbl_graph() %E>% 
  left_join(stringid2modulecolor %>% dplyr::rename_all(function(x){paste0(x,"_node1")}), by=c("node1_external_id"="stringId_node1")) %>% 
  left_join(stringid2modulecolor %>% dplyr::rename_all(function(x){paste0(x,"_node2")}), by=c("node1_external_id"="stringId_node2")) %N>%
  left_join(stringid2modulecolor, by=c("name"="preferredName")) 


tmp<-interaction %>% 
  left_join(stringid2modulecolor %>% dplyr::rename_all(function(x){paste0(x,"_node1")}), by=c("node1_external_id"="stringId_node1")) %>% 
  left_join(stringid2modulecolor %>% dplyr::rename_all(function(x){paste0(x,"_node2")}), by=c("node2_external_id"="stringId_node2")) %>% 
  filter((module_node1==centered_module)|(module_node2==centered_module))
contact_gene<-unique(c(tmp$node1,tmp$node2))



PPI_graph %>% 
  ggraph(layout = "graphopt")+
  geom_edge_link()+
  geom_node_point(aes(color=module),size=10)+
  geom_node_text(aes(label=name ))+
  scale_color_identity()





PPI_graph %>% 
  #activate("edges") %>% filter(module_node1 ==centered_module | module_node2 ==centered_module)%N>% 
  activate("nodes") %>% filter(name %in% contact_gene) %>% 
  ggraph(layout = "graphopt")+
  geom_edge_link()+
  geom_node_point(aes(color=module),size=10)+
  geom_node_text(aes(label=name ))+
  scale_color_identity()



centered_module_id<-PPI_graph %>% 
  #activate("edges") %>% filter(module_node1 ==centered_module | module_node2 ==centered_module)%N>% 
  activate("nodes") %>% filter(name %in% contact_gene) %N>% 
  as.data.frame() %>% 
  rownames_to_column("ID") %>% 
  filter(module==centered_module) %>% 
  .$name


full_graph<-create_complete(length(centered_module_id)) %N>% 
  mutate(name=centered_module_id) %E>% 
  mutate(visible=F)


set.seed(28)
PPI_graph %>% 
  #activate("edges") %>% filter(module_node1 ==centered_module | module_node2 ==centered_module)%N>% 
  activate("nodes") %>% filter(name %in% contact_gene) %>% 
  activate("edges") %>% mutate(visible=T) %>% 
  graph_join(full_graph,by = "name") %>% 
  activate("edges") %>% mutate(edge_color=if_else(visible, "grey66", "transparent")) %>% 
  ggraph(layout = "graphopt")+
  geom_edge_link0(aes(color=edge_color))+
  geom_node_point(aes(color=module),size=10)+
  geom_node_text(aes(label=name ))+
  scale_color_identity()+
  scale_edge_color_identity()+
  theme_graph()+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))


ggsave("../2.Output/20_PTPN11/PPI_Network_turquoise_centered.png",
       height = 150,width = 200,dpi = 450,units = "mm",bg="transparent")

PPI_graph %>% 
  #activate("edges") %>% filter(module_node1 ==centered_module | module_node2 ==centered_module)%N>% 
  activate("nodes") %>% filter(name %in% contact_gene) %>% 
  activate("edges") %>% mutate(visible=T) %>% 
  graph_join(full_graph,by = "name") %>% 
  activate("edges") %>% mutate(edge_color=if_else(visible, "grey66", "transparent")) %N>%
  as.data.frame() %>%
  rename(label=name) %>% 
  rownames_to_column("nodeID") %>% 
  mutate(colorcode=col2hex(module)) %>% 
  write_csv("../2.Output/20_PTPN11/Turquoise_PPI_NodeData4cytoscape.csv")
PPI_graph %>% 
  #activate("edges") %>% filter(module_node1 ==centered_module | module_node2 ==centered_module)%N>% 
  activate("nodes") %>% filter(name %in% contact_gene) %>% 
  activate("edges") %>% mutate(visible=T) %>% 
  graph_join(full_graph,by = "name") %>% 
  activate("edges") %>% mutate(edge_color=if_else(visible, "grey66", "transparent")) %E>% 
  as_tibble() %>% 
  write_csv("../2.Output/20_PTPN11/Turquoise_PPI_EdgeData4cytoscape.csv")


centered_gene="PTPN11"

tmp<-interaction %>% 
  filter((node1==centered_gene)|( node2 ==centered_gene))
contact_gene<-unique(c(tmp$node1,tmp$node2))

PPI_graph %>% 
  activate("nodes") %>%
  filter(name %in% contact_gene) %>% 
  activate("nodes") %>%
  as.data.frame() %>%
  rename(label=name) %>% 
  rownames_to_column("nodeID") %>% 
  mutate(colorcode=col2hex(module)) %>% 
  write_csv("../2.Output/20_PTPN11/PTPN_PPI_NodeData4cytoscape.csv")

PPI_graph %>% 
  activate("nodes") %>%
  filter(name %in% contact_gene) %>% 
  activate("edges") %>%
  as_tibble() %>% 
  write_csv("../2.Output/20_PTPN11/PTPN_PPI_EdgeData4cytoscape.csv")







021_ITGAM_PPI_Network.R





library(tidyverse)
library(org.Hs.eg.db)
library(tidygraph)
library(ggraph)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  dplyr::select(c("Entrez_Gene_ID", "Symbol",  "module")) %>% 
  mutate(Entrez_Gene_ID=Entrez_Gene_ID %>% as.character()) %>% 
  filter(module %in% significant_module) %>% 
  distinct(Entrez_Gene_ID,.keep_all = T)

ppi<-read_tsv("../2.Output/21_ITGAM/string_interactions_ITGAM.tsv")

ppi_graph<-ppi %>% 
  as_tbl_graph() %>% 
  activate("nodes") %>% 
  mutate(Entrez_Gene_ID=mapIds(org.Hs.eg.db, name, 'ENTREZID',"SYMBOL") %>% as.vector()) %>%
  left_join(geneid2modulecolor,by=c("Entrez_Gene_ID"="Entrez_Gene_ID")) %>% 
  mutate(color=if_else(is.na(module), "grey", module))

ppi_graph %N>%
  as.data.frame() %>%
  rename(label=name) %>% 
  rownames_to_column("nodeID") %>% 
  write_csv("../2.Output/21_ITGAM/ITGAM_PPI_NodeData4cytoscape.csv")
ppi_graph %E>% 
  as_tibble() %>% 
  write_csv("../2.Output/21_ITGAM/ITGAM_PPI_EdgeData4cytoscape.csv")



ppi_graph %>% 
  ggraph(layout = "kk")+
  geom_edge_link0()+
  geom_node_point(aes(color=color),size=20)+
  geom_node_text(aes(label=name))+
  scale_color_identity()





















021_Lightcyan_Centered_PPI_Network.R





library(tidyverse)
library(tidygraph)
library(ggraph)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  dplyr::select(c("Entrez_Gene_ID", "Symbol",  "module")) %>% 
  filter(module %in% significant_module)
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor

interaction<-read_tsv("../2.Output/19_ImmuneCell-Fibroblast_Signaling/2.PPI_string/inSRB_sigMod/string_interactions.tsv") %>% 
  dplyr::select(c(`#node1`, "node2","node1_external_id", "node2_external_id")) %>% 
  dplyr::rename(node1=`#node1`)
mapping<-read_tsv("../2.Output/19_ImmuneCell-Fibroblast_Signaling/2.PPI_string/inSRB_sigMod/string_mapping.tsv")
stringid2modulecolor<-mapping %>% 
  left_join(geneid2modulecolor,by=c("queryItem"="Entrez_Gene_ID")) %>% 
  dplyr::select(c("queryItem", "stringId" , "module","preferredName")) %>% 
  dplyr::rename(Entrez_Gene_ID=queryItem) %>% 
  distinct(stringId, .keep_all=T )

centered_module="lightcyan"


PPI_graph<-interaction %>% 
  as_tbl_graph() %E>% 
  left_join(stringid2modulecolor %>% dplyr::rename_all(function(x){paste0(x,"_node1")}), by=c("node1_external_id"="stringId_node1")) %>% 
  left_join(stringid2modulecolor %>% dplyr::rename_all(function(x){paste0(x,"_node2")}), by=c("node1_external_id"="stringId_node2")) %N>%
  left_join(stringid2modulecolor, by=c("name"="preferredName")) 


tmp<-interaction %>% 
  left_join(stringid2modulecolor %>% dplyr::rename_all(function(x){paste0(x,"_node1")}), by=c("node1_external_id"="stringId_node1")) %>% 
  left_join(stringid2modulecolor %>% dplyr::rename_all(function(x){paste0(x,"_node2")}), by=c("node2_external_id"="stringId_node2")) %>% 
  filter((module_node1==centered_module)|(module_node2==centered_module))
contact_gene<-unique(c(tmp$node1,tmp$node2))



PPI_graph %>% 
  ggraph(layout = "graphopt")+
  geom_edge_link()+
  geom_node_point(aes(color=module),size=10)+
  geom_node_text(aes(label=name ))+
  scale_color_identity()





PPI_graph %>% 
  #activate("edges") %>% filter(module_node1 ==centered_module | module_node2 ==centered_module)%N>% 
  activate("nodes") %>% filter(name %in% contact_gene) %>% 
  ggraph(layout = "graphopt")+
  geom_edge_link()+
  geom_node_point(aes(color=module),size=10)+
  geom_node_text(aes(label=name ))+
  scale_color_identity()



centered_module_id<-PPI_graph %>% 
  #activate("edges") %>% filter(module_node1 ==centered_module | module_node2 ==centered_module)%N>% 
  activate("nodes") %>% filter(name %in% contact_gene) %N>% 
  as.data.frame() %>% 
  rownames_to_column("ID") %>% 
  filter(module==centered_module) %>% 
  .$name


full_graph<-create_complete(length(centered_module_id)) %N>% 
  mutate(name=centered_module_id) %E>% 
  mutate(visible=F)


set.seed(28)
PPI_graph %>% 
  #activate("edges") %>% filter(module_node1 ==centered_module | module_node2 ==centered_module)%N>% 
  activate("nodes") %>% filter(name %in% contact_gene) %>% 
  activate("edges") %>% mutate(visible=T) %>% 
  graph_join(full_graph,by = "name") %>% 
  activate("edges") %>% mutate(edge_color=if_else(visible, "grey66", "transparent")) %>% 
  ggraph(layout = "graphopt")+
  geom_edge_link0(aes(color=edge_color))+
  geom_node_point(aes(color=module),size=10)+
  geom_node_text(aes(label=name ))+
  scale_color_identity()+
  scale_edge_color_identity()+
  theme_graph()+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))


ggsave("../2.Output/19_ImmuneCell-Fibroblast_Signaling/PPI_Network_lightcyan_centered.png",
       height = 150,width = 200,dpi = 450,units = "mm",bg="transparent")

for (i in 1:100){
  set.seed(i)
  g<-PPI_graph %>% 
    #activate("edges") %>% filter(module_node1 ==centered_module | module_node2 ==centered_module)%N>% 
    activate("nodes") %>% filter(name %in% yellow_contact_gene) %>% 
    activate("edges") %>% mutate(visible=T) %>% 
    graph_join(full_graph,by = "name") %>% 
    activate("edges") %>% mutate(edge_color=if_else(visible, "grey", "transparent")) %>% 
    ggraph(layout = "graphopt")+
    geom_edge_link(aes(color=edge_color))+
    geom_node_point(aes(color=module),size=10)+
    geom_node_text(aes(label=name ))+
    scale_color_identity()+
    scale_edge_color_identity()+
    theme_graph()+
    theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
          panel.border = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          panel.grid.major = element_blank(), # get rid of major grid
          panel.grid.minor = element_blank(), # get rid of minor grid
          legend.background = element_rect(fill = "transparent"), # get rid of legend bg
          legend.box.background = element_rect(fill = "transparent"))+
    labs(title = i)
  print(g)
  
}













022_RetrieveFibrosisCausalGene.R





library(tidyverse)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

gene_attribute_matrix <- read_tsv("../2.Output/22_OMIM_CausalGene/gene_attribute_matrix.txt")
geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  mutate(Entrez_Gene_ID=Entrez_Gene_ID %>% as.character())

gamat_colname<-gene_attribute_matrix %>% colnames()
gamat_1strow<-gene_attribute_matrix[1,] %>% as.character()
colnames(gene_attribute_matrix)<-map2_chr(gamat_colname,gamat_1strow,~str_c(.x,.y,sep="_"))

gene_attribute_matrix<-gene_attribute_matrix[3:nrow(gene_attribute_matrix),] %>% 
  dplyr::rename(GeneSymbol=`#_#`,OMIMID=`#_1_#`,GeneID=Disease_OMIMID)


fibrosis_related<-gene_attribute_matrix %>% 
  dplyr::select(GeneSymbol,OMIMID,GeneID,contains("fibrosis")) %>% 
  filter_at(.vars = vars(contains("OMIM:")), any_vars(.>=1))

gene_attribute_matrix %>% 
  dplyr::select(GeneSymbol,OMIMID,GeneID,contains("collagen")) %>% 
  filter_at(.vars = vars(contains("OMIM:")), any_vars(.>=1))

fibrosis_related %>% 
  left_join(geneid2modulecolor,by = c("GeneID"="Entrez_Gene_ID")) %>% 
  dplyr::select(GeneSymbol, OMIMID, GeneID, module, everything())  %>% 
  write_csv("../2.Output/22_OMIM_CausalGene/fibrosis_related.csv")






025_TF_target_enricjment.R





library(tidyverse)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")


TF_targets<-read_tsv("../2.Output/25_TranscriptionFactor_ofModule/ENCODE_Transcription_Factor_Targets/gene_attribute_matrix.txt",comment = "#")
geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  dplyr::select(Entrez_Gene_ID, Symbol, module)


nested<-TF_targets %>% 
  select(-X2) %>% 
  filter(`GeneID/GeneID` %in% geneid2modulecolor$Entrez_Gene_ID) %>% 
  gather(key="TF_GeneID",value="is_target",-c(GeneSym,`GeneID/GeneID`)) %>% 
  filter(is_target==1) %>% 
  select(-is_target) %>% 
  group_by(TF_GeneID) %>% 
  nest()%>% 
  mutate(targets_EntrezID=data[[1]] %>% .$`GeneID/GeneID` %>% list())
  
phyper_result<-map(nested$targets_EntrezID,
    .f = function(TFtargets){
      phyper_result=geneid2modulecolor%>% 
        filter(module!="grey") %>% 
        mutate(is_TFtarget=Entrez_Gene_ID %in% TFtargets) %>% 
        group_by(module) %>% 
        nest() %>% 
        mutate(genenum_module=nrow(data[[1]]),
               num_TF_target_inModule=nrow(data[[1]] %>% filter(is_TFtarget)),
               total_TF_target=length(TFtargets)) %>% 
        filter(num_TF_target_inModule!=0) %>% 
        mutate(target_gene_symbol=data[[1]] %>% filter(is_TFtarget) %>% head(100) %>% pull(Symbol) %>% str_c(collapse =","),#str_c( pull( filter(data[[1]],is_TFtarget) ,Symbol) ,collapse =","),
               p.value=phyper(num_TF_target_inModule-1, 
                              length(TFtargets),
                              nrow(geneid2modulecolor)-length(TFtargets), 
                              genenum_module, lower.tail= FALSE))%>% 
        arrange(p.value)
      return(phyper_result)
    }
  )
phyper_result

result=tibble()

for (i in seq_along(phyper_result)){
  add<-phyper_result[[i]] %>% 
    as_tibble() %>% 
    dplyr::select(-data) %>% 
    mutate(q.value=p.adjust(phyper_result[[i]]$p.value, method = "BH"),
           TF_GeneID=nested$TF_GeneID[[i]])
  result<-result %>% 
    bind_rows(add)
}

geneid2TFsymbol<-read_tsv("../2.Output/25_TranscriptionFactor_ofModule/ENCODE_Transcription_Factor_Targets/gene_attribute_matrix.txt")
geneid2TFsymbol=tibble(
  Symbol=geneid2TFsymbol %>% colnames() %>% .[4:ncol(geneid2TFsymbol)],
  Entrez_Gene_ID=geneid2TFsymbol[2,] %>% as.character()%>% .[4:ncol(geneid2TFsymbol)])

geneid2TFsymbol %>% tail()


result %>% 
  filter(q.value<0.01) %>% 
  left_join(geneid2TFsymbol, by=c("TF_GeneID"="Entrez_Gene_ID")) %>% 
  rename(TF_Symbol="Symbol") %>% 
  left_join(geneid2modulecolor %>% 
              mutate(Entrez_Gene_ID=Entrez_Gene_ID %>% as.character()) %>% 
              rename(TF_module="module"),by=c("TF_GeneID"="Entrez_Gene_ID")) %>% 
  arrange(module) %>% 
  dplyr::select(module,TF_Symbol,TF_GeneID , TF_module, p.value,q.value, genenum_module, num_TF_target_inModule, total_TF_target, target_gene_symbol) %>% 
  write_csv("../2.Output/25_TranscriptionFactor_ofModule/encode_TFtargets_summary_ver3.csv")
  #filter(module=="lightyellow")






025_TF_target_enricjment_ChIP-Atlas.R





library(tidyverse)
library(org.Hs.eg.db)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

YAP_targets<-read_tsv("../2.Output/25_TranscriptionFactor_ofModule/ChIP-Atlas/YAP1.tsv")
YAP_targets<-tibble(Symbol=YAP_targets$Target_genes,
       Entrez_Gene_ID=mapIds(org.Hs.eg.db, Symbol, 'ENTREZID', 'SYMBOL'))
YAP_targets %>% 
  filter(is.na(Entrez_Gene_ID)) %>% 
  write_csv("../2.Output/25_TranscriptionFactor_ofModule/NA_genes.csv")

###############################


geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  dplyr::select(Entrez_Gene_ID, Symbol, module)
YAP_targets_entrezid<-YAP_targets %>% 
  filter(!is.na(Entrez_Gene_ID) & (Entrez_Gene_ID %in% geneid2modulecolor$Entrez_Gene_ID)) %>% 
  pull(Entrez_Gene_ID)



geneid2modulecolor %>% 
  mutate(is_target=Entrez_Gene_ID %in% YAP_targets_entrezid) %>% 
  filter(is_target) %>% .$module %>% 
  table()

nrow(geneid2modulecolor)# 12875
length(YAP_targets_entrezid)# 1362

# hitInSample = num_PTPN11_target  
# hitInPop = length(PTPN11_targets_humanEntrezID) 
# failInPop = nrow(geneid2modulecolor) - hitInPop 
# sampleSize = genenum_module
phyper_result<-geneid2modulecolor%>% 
  mutate(is_target=Entrez_Gene_ID %in% YAP_targets_entrezid) %>% 
  group_by(module) %>% 
  filter(module!="grey") %>% 
  nest() %>% 
  mutate(genenum_module=nrow(data[[1]]),
         num_targets=nrow(data[[1]] %>% filter(is_target))) %>% 
  filter(num_targets!=0) %>% 
  mutate(target_gene_symbol=data[[1]] %>% filter(is_target) %>% head(100) %>% pull(Symbol) %>% str_c(collapse =","),
         p.value=phyper(num_targets-1, 
                        length(YAP_targets_entrezid),
                        nrow(geneid2modulecolor)-length(YAP_targets_entrezid), 
                        genenum_module, lower.tail= FALSE) ) %>% 
  arrange(p.value)
phyper_result<-phyper_result %>% 
  dplyr::select(-data) %>% 
  as_tibble() %>% 
  mutate(q.value=p.adjust(phyper_result$p.value, method = "BH"))  %>% 
  dplyr::select(module, genenum_module, num_targets, p.value, q.value, target_gene_symbol) %>% 
  write_csv("../2.Output/25_TranscriptionFactor_ofModule/ChIP-Atlas_YAP1_summary.csv")


########################################
TAZ_targets<-read_tsv("../2.Output/25_TranscriptionFactor_ofModule/ChIP-Atlas/TAZ.tsv")
TAZ_targets<-tibble(Symbol=TAZ_targets$Target_genes,
                    Entrez_Gene_ID=mapIds(org.Hs.eg.db, Symbol, 'ENTREZID', 'SYMBOL'))
TAZ_targets %>% 
  filter(is.na(Entrez_Gene_ID)) %>% 
  write_csv("../2.Output/25_TranscriptionFactor_ofModule/NA_genes.csv")

geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  dplyr::select(Entrez_Gene_ID, Symbol, module)
TAZ_targets_entrezid<-TAZ_targets %>% 
  filter(!is.na(Entrez_Gene_ID) & (Entrez_Gene_ID %in% geneid2modulecolor$Entrez_Gene_ID)) %>% 
  pull(Entrez_Gene_ID)



geneid2modulecolor %>% 
  mutate(is_target=Entrez_Gene_ID %in% TAZ_targets_entrezid) %>% 
  filter(is_target) %>% .$module %>% 
  table()

nrow(geneid2modulecolor)# 12875
length(TAZ_targets_entrezid)# 1362

# hitInSample = num_PTPN11_target  
# hitInPop = length(PTPN11_targets_humanEntrezID) 
# failInPop = nrow(geneid2modulecolor) - hitInPop 
# sampleSize = genenum_module

phyper_result<-geneid2modulecolor%>% 
  mutate(is_target=Entrez_Gene_ID %in% TAZ_targets_entrezid) %>% 
  group_by(module) %>% 
  filter(module!="grey") %>% 
  nest() %>% 
  mutate(genenum_module=nrow(data[[1]]),
         num_targets=nrow(data[[1]] %>% filter(is_target))) %>% 
  filter(num_targets!=0) %>% 
  mutate(target_gene_symbol=data[[1]] %>% filter(is_target) %>% head(100) %>% pull(Symbol) %>% str_c(collapse =","),
         p.value=phyper(num_targets-1, 
                        length(TAZ_targets_entrezid),
                        nrow(geneid2modulecolor)-length(TAZ_targets_entrezid), 
                        genenum_module, lower.tail= FALSE) ) %>% 
  arrange(p.value)
phyper_result<-phyper_result %>% 
  dplyr::select(-data) %>% 
  as_tibble() %>% 
  mutate(q.value=p.adjust(phyper_result$p.value, method = "BH"))  %>% 
  dplyr::select(module, genenum_module, num_targets, p.value, q.value, target_gene_symbol) %>% 
  write_csv("../2.Output/25_TranscriptionFactor_ofModule/ChIP-Atlas_TAZ1_summary.csv")



##########################################
#SRX757326
YAP_targets<-read_tsv("../2.Output/25_TranscriptionFactor_ofModule/ChIP-Atlas/YAP1.tsv")%>% 
  filter(`SRX757326|H2052`>0)
  
YAP_targets<-tibble(Symbol=YAP_targets$Target_genes,
                    Entrez_Gene_ID=mapIds(org.Hs.eg.db, Symbol, 'ENTREZID', 'SYMBOL'))


geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  dplyr::select(Entrez_Gene_ID, Symbol, module)
YAP_targets_entrezid<-YAP_targets %>% 
  filter(!is.na(Entrez_Gene_ID) & (Entrez_Gene_ID %in% geneid2modulecolor$Entrez_Gene_ID)) %>% 
  pull(Entrez_Gene_ID)



geneid2modulecolor %>% 
  mutate(is_target=Entrez_Gene_ID %in% YAP_targets_entrezid) %>% 
  filter(is_target) %>% .$module %>% 
  table()

nrow(geneid2modulecolor)# 12875
length(YAP_targets_entrezid)# 1362

# hitInSample = num_PTPN11_target  
# hitInPop = length(PTPN11_targets_humanEntrezID) 
# failInPop = nrow(geneid2modulecolor) - hitInPop 
# sampleSize = genenum_module
phyper_result<-geneid2modulecolor%>% 
  mutate(is_target=Entrez_Gene_ID %in% YAP_targets_entrezid) %>% 
  group_by(module) %>% 
  filter(module!="grey") %>% 
  nest() %>% 
  mutate(genenum_module=nrow(data[[1]]),
         num_targets=nrow(data[[1]] %>% filter(is_target))) %>% 
  filter(num_targets!=0) %>% 
  mutate(target_gene_symbol=data[[1]] %>% filter(is_target) %>% head(100) %>% pull(Symbol) %>% str_c(collapse =","),
         p.value=phyper(num_targets-1, 
                        length(YAP_targets_entrezid),
                        nrow(geneid2modulecolor)-length(YAP_targets_entrezid), 
                        genenum_module, lower.tail= FALSE) ) %>% 
  arrange(p.value)
phyper_result<-phyper_result %>% 
  dplyr::select(-data) %>% 
  as_tibble() %>% 
  mutate(q.value=p.adjust(phyper_result$p.value, method = "BH"))  %>% 
  dplyr::select(module, genenum_module, num_targets, p.value, q.value, target_gene_symbol) %>% 
  write_csv("../2.Output/25_TranscriptionFactor_ofModule/ChIP-Atlas_YAP1_summary_SRX757326.csv")

##########################################
#SRX1038944|IMR-90
YAP_targets<-read_tsv("../2.Output/25_TranscriptionFactor_ofModule/ChIP-Atlas/YAP1.tsv")%>% 
  filter(`SRX1038944|IMR-90`>0 | `SRX1038945|IMR-90`>0)

YAP_targets<-tibble(Symbol=YAP_targets$Target_genes,
                    Entrez_Gene_ID=mapIds(org.Hs.eg.db, Symbol, 'ENTREZID', 'SYMBOL'))


geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  dplyr::select(Entrez_Gene_ID, Symbol, module)
YAP_targets_entrezid<-YAP_targets %>% 
  filter(!is.na(Entrez_Gene_ID) & (Entrez_Gene_ID %in% geneid2modulecolor$Entrez_Gene_ID)) %>% 
  pull(Entrez_Gene_ID)



geneid2modulecolor %>% 
  mutate(is_target=Entrez_Gene_ID %in% YAP_targets_entrezid) %>% 
  filter(is_target) %>% .$module %>% 
  table()

nrow(geneid2modulecolor)# 12875
length(YAP_targets_entrezid)# 1362

# hitInSample = num_PTPN11_target  
# hitInPop = length(PTPN11_targets_humanEntrezID) 
# failInPop = nrow(geneid2modulecolor) - hitInPop 
# sampleSize = genenum_module
phyper_result<-geneid2modulecolor%>% 
  mutate(is_target=Entrez_Gene_ID %in% YAP_targets_entrezid) %>% 
  group_by(module) %>% 
  filter(module!="grey") %>% 
  nest() %>% 
  mutate(genenum_module=nrow(data[[1]]),
         num_targets=nrow(data[[1]] %>% filter(is_target))) %>% 
  filter(num_targets!=0) %>% 
  mutate(target_gene_symbol=data[[1]] %>% filter(is_target) %>% head(100) %>% pull(Symbol) %>% str_c(collapse =","),
         p.value=phyper(num_targets-1, 
                        length(YAP_targets_entrezid),
                        nrow(geneid2modulecolor)-length(YAP_targets_entrezid), 
                        genenum_module, lower.tail= FALSE) ) %>% 
  arrange(p.value)
phyper_result<-phyper_result %>% 
  dplyr::select(-data) %>% 
  as_tibble() %>% 
  mutate(q.value=p.adjust(phyper_result$p.value, method = "BH"))  %>% 
  dplyr::select(module, genenum_module, num_targets, p.value, q.value, target_gene_symbol) %>% 
  write_csv("../2.Output/25_TranscriptionFactor_ofModule/ChIP-Atlas_YAP1_summary_SRX1038944.csv")






025_TF_target_enricjment_test.R





library(tidyverse)

setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")

TF_targets<-read_tsv("../2.Output/25_TranscriptionFactor_ofModule/ENCODE_Transcription_Factor_Targets/gene_attribute_matrix.txt",comment = "#")
geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  dplyr::select(Entrez_Gene_ID, Symbol, module)


targets_EntrezID<-TF_targets %>% 
  dplyr::select( GeneSym ,`GeneID/GeneID` ,`1820`) %>% 
  filter(`GeneID/GeneID` %in% geneid2modulecolor$Entrez_Gene_ID) %>% 
  filter(`1820`==1) %>% 
  .$`GeneID/GeneID`


#allgene_in_Analysis<-
#read_csv("../1.Data/GSE58095_main/Annotation/GeneAnnotation_GPL10558.csv")


geneid2modulecolor %>% 
  mutate(is_TFtarget=Entrez_Gene_ID %in% targets_EntrezID) %>% 
  filter(is_TFtarget) %>% .$module %>% 
  table()

nrow(geneid2modulecolor)# 12875
length(targets_EntrezID)#  6801

# hitInSample = num_TF_target  
# hitInPop = length(PTPN11_targets_humanEntrezID) 
# failInPop = nrow(geneid2modulecolor) - hitInPop 
# sampleSize = genenum_module

phyper_result<-geneid2modulecolor%>% 
  mutate(is_TFtarget=Entrez_Gene_ID %in% targets_EntrezID) %>% 
  group_by(module) %>% 
  nest() %>% 
  mutate(genenum_module=nrow(data[[1]]),
         num_TF_target=nrow(data[[1]] %>% filter(is_TFtarget)),
         p.value=phyper(num_TF_target-1, 
                        length(targets_EntrezID),
                        nrow(geneid2modulecolor)-length(targets_EntrezID), 
                        genenum_module, lower.tail= FALSE)) %>% 
  arrange(p.value)


######################################################
TF_targets<-read_tsv("../2.Output/25_TranscriptionFactor_ofModule/ENCODE_Transcription_Factor_Targets/gene_attribute_matrix.txt",comment = "#")
geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  dplyr::select(Entrez_Gene_ID, Symbol, module)

targets_EntrezID<-TF_targets %>% 
  dplyr::select( GeneSym ,`GeneID/GeneID` ,`1820`) %>% 
  filter(`GeneID/GeneID` %in% geneid2modulecolor$Entrez_Gene_ID) %>% 
  filter(`1820`==1) %>% 
  .$`GeneID/GeneID`

nrow(geneid2modulecolor)# 12875
length(targets_EntrezID)# 1362

# hitInSample = num_TF_target  
# hitInPop = length(PTPN11_targets_humanEntrezID) 
# failInPop = nrow(geneid2modulecolor) - hitInPop 
# sampleSize = genenum_module

phyper_result<-geneid2modulecolor%>% 
  mutate(is_TFtarget=Entrez_Gene_ID %in% targets_EntrezID) %>% 
  group_by(module) %>% 
  nest() %>% 
  mutate(genenum_module=nrow(data[[1]]),
         num_TF_target=nrow(data[[1]] %>% filter(is_TFtarget)),
         p.value=phyper(num_TF_target-1, 
                        length(targets_EntrezID),
                        nrow(geneid2modulecolor)-length(targets_EntrezID), 
                        genenum_module, lower.tail= FALSE)) %>% 
  arrange(p.value)

######################################################
TF_targets<-read_tsv("../2.Output/25_TranscriptionFactor_ofModule/ENCODE_Transcription_Factor_Targets/gene_attribute_matrix.txt",comment = "#")
geneid2modulecolor<-read_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv") %>% 
  dplyr::select(Entrez_Gene_ID, Symbol, module)


nested<-TF_targets %>% 
  select(-X2) %>% 
  filter(`GeneID/GeneID` %in% geneid2modulecolor$Entrez_Gene_ID) %>% 
  gather(key="TF_GeneID",value="is_target",-c(GeneSym,`GeneID/GeneID`)) %>% 
  filter(is_target==1) %>% 
  select(-is_target) %>% 
  group_by(TF_GeneID) %>% 
  nest()%>% 
  mutate(targets_EntrezID=data[[1]] %>% .$`GeneID/GeneID` %>% list())

phyper_result<-map(nested$targets_EntrezID,
                   .f = function(TFtargets){
                     phyper_result=geneid2modulecolor%>% 
                       mutate(is_TFtarget=Entrez_Gene_ID %in% TFtargets) %>% 
                       group_by(module) %>% 
                       nest() %>% 
                       mutate(genenum_module=nrow(data[[1]]),
                              num_TF_target_inModule=nrow(data[[1]] %>% filter(is_TFtarget)),
                              total_TF_target=length(TFtargets),
                              target_gene_symbol=data[[1]] %>% filter(is_TFtarget) %>% .$Symbol %>% str_c(collapse =","),
                              p.value=phyper(num_TF_target_inModule-1, 
                                             length(TFtargets),
                                             nrow(geneid2modulecolor)-length(TFtargets), 
                                             genenum_module, lower.tail= FALSE))%>% 
                       arrange(p.value)
                     return(phyper_result)
                   }
)
phyper_result
TFtargets=nested$targets_EntrezID[[1]]
geneid2modulecolor%>% 
  filter(module!="grey") %>% 
  mutate(is_TFtarget=Entrez_Gene_ID %in% TFtargets) %>% 
  group_by(module) %>% 
  nest() %>% 
  mutate(genenum_module=nrow(data[[1]]),
         num_TF_target_inModule=nrow(data[[1]] %>% filter(is_TFtarget)),
         total_TF_target=length(TFtargets),
         target_gene_symbol=data[[1]] %>% filter(is_TFtarget) %>% .$Symbol %>% str_c(collapse =","))



phyper_result[[1]] %>% 
  as_tibble() %>% 
  mutate(q.value=p.adjust(phyper_result[[1]]$p.value, method = "BH"),
         TF_GeneID=nested$TF_GeneID[[1]])

result=tibble()

for (i in seq_along(phyper_result)){
  add<-phyper_result[[i]] %>% 
    as_tibble() %>% 
    dplyr::select(-data) %>% 
    mutate(q.value=p.adjust(phyper_result[[i]]$p.value, method = "BH"),
           TF_GeneID=nested$TF_GeneID[[i]])
  result<-result %>% 
    bind_rows(add)
}

result %>% 
  filter(q.value<0.01) %>% 
  group_by(module) %>% 
  nest()->a

geneid2TFsymbol<-read_tsv("../2.Output/25_TranscriptionFactor_ofModule/ENCODE_Transcription_Factor_Targets/gene_attribute_matrix.txt")
geneid2TFsymbol=tibble(
  Symbol=geneid2TFsymbol %>% colnames() %>% .[4:ncol(geneid2TFsymbol)],
  Entrez_Gene_ID=geneid2TFsymbol[2,] %>% as.character()%>% .[4:ncol(geneid2TFsymbol)])

geneid2TFsymbol %>% tail()


result %>% 
  filter(q.value<0.01) %>% 
  left_join(geneid2TFsymbol, by=c("TF_GeneID"="Entrez_Gene_ID")) %>% 
  arrange(module) %>% 
  write_csv("../2.Output/25_TranscriptionFactor_ofModule/encode_TFtargets_summary.csv")
#filter(module=="lightyellow")






24_CellTypeEnrichment_biperticegraph.R





library(tidyverse)
library(tidygraph)
library(igraph)
library(ggraph)
library(gplots)

significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor

celltype_enrichment_result<-read_csv("../2.Output/04_CelltypeEnrichment/ARCHS4_Tissues_ver2/WGCNA_enrichment_result_ver2.csv") 

celltype_enrichment_result %>% 
  filter(!(dataSetID %in% c("BLASTOCYST"))) %>% 
  filter(FDR<0.10) %>% 
  filter(class %in% significant_module) %>% 
  group_by(class) %>% 
  group_modify(~ head(.x, 5)) %>% 
  dplyr::select(c(1,3),FDR) %>% 
  rename(from="class",to="dataSetID") %>% 
  as_tbl_graph() %>% #directed = TRUE) %>% 
  activate("nodes") %>% 
  mutate(type=if_else(name %in% celltype_enrichment_result$class,1,0),
         label_module=if_else(name %in% celltype_enrichment_result$class,name,"red"),
         label_celltype=if_else(name %in% celltype_enrichment_result$class,"",name),
         point_alpha=if_else(name %in% celltype_enrichment_result$class,1,0)) %>% 
  create_layout(layout = "bipartite")->layout
layout$x<--layout$x
layout %>% 
  ggraph(layout = "bipartite") +
  geom_edge_link0(width=1,color="#0086ab",alpha=0.8)+
  geom_node_label(aes(label=label_celltype %>% str_replace_all("_"," ")%>% str_to_title()),size=5,color="grey20")+
  geom_node_point(aes(alpha=point_alpha,color=label_module),size=20)+
  geom_node_text(aes(label=label_module,alpha=point_alpha),color="grey20")+
  coord_flip()+
  ylim(c(-0.5,1.5))+
  scale_alpha_identity()+
  scale_color_identity()+
  theme_graph()+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))
ggsave(filename = "../2.Output/04_CelltypeEnrichment/result_ver6.png",height = 250,width = 170,units = "mm",bg="transparent")

################################################
#?q?[?g?}?b?v???`??
color_key<-colorRampPalette(c("white",  "red"))

tmp<-celltype_enrichment_result%>% 
  filter(class %in% significant_module) %>% 
  group_by(class) %>% 
  dplyr::select(c(1,3),FDR) %>% 
  
  rename(from="class",to="dataSetID") %>% 
  as_tibble() %>% 
  mutate(to=str_replace_all(to,"_"," ") %>% str_to_title()) %>% 
  tidyr::spread(key="to",value="FDR") %>% 
  as.data.frame() %>% 
  column_to_rownames("from") %>% 
  as.matrix() %>% t()
tmp[is.na(tmp)]<-1
png(filename = "../2.Output/04_CelltypeEnrichment/Heatmap_ver4.png",height=1440, width=1600, res=216)

heatmap.2(-log10(tmp),scale = "none", col =color_key,trace = "none", density.info = "none",margins = c(10, 15),
          breaks=seq(0,10,1)) 
dev.off()

#ARCHS?S???o?[?W????
tmp2<-read_csv("../2.Output/04_CelltypeEnrichment/ARCHS4_Tissues/WGCNA_enrichment_result.csv") %>% 
  filter(class %in% significant_module) %>% 
  group_by(class) %>% 
  dplyr::select(c(1,3),FDR) %>% 
  rename(class="from",dataSetID="to") %>% 
  as_tibble() %>% 
  mutate(to=str_replace_all(to,"_"," ") %>% str_to_title()) %>% 
  tidyr::spread(key="to",value="FDR") %>% 
  as.data.frame() %>% 
  column_to_rownames("from") %>% 
  as.matrix() %>% t()
tmp2[is.na(tmp2)]<-1

-log10(tmp2) %>% 
  heatmap.2(scale = "column", col = bluered(100), 
            trace = "none", density)






BN_SelectSampler_2.R





library(tidyverse)
library(tidygraph)
library(ggraph)
library(igraph)
library(bnlearn)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
ME_skinscore_file_path<-"../../tmp/ME_SkinScore_forBN_nafilled.txt"
#result_dir_path<-"../../tmp/BN_SSc/result/"
#BN_result_dir_path<-"../../tmp/BN_SSc/BN_results/"
BN_result_files<-list.files("../../tmp/BN_SSc/result/",pattern ="Mcmc_*" ,full.names = T)
BN_result_files<-BN_result_files[-which(BN_result_files %in% "../../tmp/BN_SSc/result/Mcmc_list.rda" )]

MergedData_1<-read_tsv(ME_skinscore_file_path) %>% 
  as.data.frame() %>% 
  mutate(variable=str_replace_all(variable,"ME","")) %>% 
  column_to_rownames("variable")
if(0){
tmp<-tibble(
  module=rownames(MergedData_1),
  module_function=c("MFB,","FB,ECM","MFB,Protein lipidation,Blood","DC,Mitochon","Immune ,B?","FB,ECM","MP,Apoptosis",
                    "DC,ER, RNA","FB,Blood Vessel Development","FB,MP,Blood Vessel Development","Keratinization","Epithelium, Cell Cycle",
                    "DC,ER, RNA","DC, Nucleic Acid Metabolism","Epithelium, Cell Cycle","","")) 
}

module2function<-tibble(
  module=rownames(MergedData_1),
  celltype=c("MFB","FB", "MFB" ,"DC","Immune","FB" ,"MP","DC","FB","FB_MP","Keratinocyte","Epithelium","DC","DC","Epithelium","",""),
  module_function=c("","ECM","Protein lipidation","Mitochondria","Immune response","ECM","Apoptosis","ER_RNA","Blood Vessel Development","Blood Vessel Development",
                    "Keratinization","Cell Cycle","ER_RNA","Nucleic Acid Metabolism","Cell Cycle","",""))

immune_cells<-c("DC","Immune","MP")
fibroblast<-c("MFB","FB","FB_MP")
traits<-c("Skin_score_at_biopsy_site","total_skin_score")

BN_result_file<-BN_result_files[2]
for (BN_result_file in BN_result_files){
  print(BN_result_file)
  g<-read_csv(BN_result_file,col_names = rownames(MergedData_1))%>% 
    mutate(from=rownames(MergedData_1)) %>% 
    gather(key=to,value=edge,-from) %>% 
    filter(edge==1) %>% 
    select(-edge) %>% 
    as_tbl_graph() %N>% 
    left_join(module2function,by=c("name"="module")) %>% 
    mutate(module_color=if_else(name %in% traits,"red",name),    #name %>% str_replace("Skin_score_at_biopsy_site","red") %>% str_replace("total_skin_score" ,"red"),
           celltype_large=if_else(celltype %in% immune_cells,"immune cell",if_else(celltype %in% fibroblast,"fibroblast","")),
           celltype_color=if_else(celltype_large=="immune cell","#f79646",if_else(celltype_large=="fibroblast","#9bbb59","grey")) )%>% 
    ggraph(layout = "sugiyama")+
    geom_edge_diagonal2(arrow = arrow(length = unit(4, 'mm'),type = "closed"), 
                   end_cap = circle(10, 'mm'))+
    geom_node_point(aes(color=celltype_color),size=20)+
    geom_node_point(aes(color=module_color),size=10)+
    geom_node_text(aes(label=str_c(name,"\n",module_function)))+
    scale_color_identity()+
    theme_graph()+
    theme(plot.margin= unit(c(0,0,0,0), "lines"))
  savefile<-paste0("../2.Output/05_BayesianNetwork/Image/",str_sub(basename(BN_result_file),6,-5),".png")
  ggsave(plot = g,filename = savefile,dpi=300)
}
#"#0086ab","#f79646","#9bbb59","#da6272","#777777","#bfbfbf"

BN_result_file<-"../../tmp/BN_SSc/result/Mcmc_c2PB.csv"

graph<-read_csv(BN_result_file,col_names = rownames(MergedData_1))%>% 
  mutate(from=rownames(MergedData_1)) %>% 
  gather(key=to,value=edge,-from) %>% 
  filter(edge==1) %>% 
  select(-edge) %>% 
  as_tbl_graph() %N>% 
  left_join(module2function,by=c("name"="module")) %>% 
  mutate(module_color=if_else(name %in% traits,"red",name),    #name %>% str_replace("Skin_score_at_biopsy_site","red") %>% str_replace("total_skin_score" ,"red"),
         celltype_large=if_else(celltype %in% immune_cells,"immune cell",if_else(celltype %in% fibroblast,"fibroblast","")),
         celltype_color=if_else(celltype_large=="immune cell","#f79646",if_else(celltype_large=="fibroblast","#9bbb59","grey")) )

layout<-create_layout(graph,layout = "sugiyama")
#####################
layout$x<--layout$x
layout$y<--layout$y

#tmp<-layout[layout$name=="white",c("x","y")]
#layout[layout$name=="white",c("x","y")]<-layout[layout$name=="yellow",c("x","y")]
#layout[layout$name=="yellow",c("x","y")]<-tmp

layout[layout$name=="yellow",c("x","y")]<-c(0.5,2.3)
layout[layout$name=="tan",c("x","y")]<-c(-0.75,2)
layout[layout$name=="white",c("x","y")]<-c(1,2)
layout[layout$name=="turquoise",c("x","y")]<-c(0,1.7)

layout[layout$name=="darkmagenta",c("x","y")]<-c(-1,1.1)
layout[layout$name=="greenyellow",c("x","y")]<-c(0,1.1)

layout[layout$name=="royalblue",c("x","y")]<-c(-0.5,0.8)
layout[layout$name=="orange",c("x","y")]<-c(0.5,0.8)

layout[layout$name=="lightyellow",c("x","y")]<-c(0,0.5)

layout[layout$name=="midnightblue",c("x","y")]<-c(0.5,0.2)
layout[layout$name=="violet",c("x","y")]<-c(1.25,0.2)

layout[layout$name=="darkred",c("x","y")]<-c(0.5,-0.4)

layout[layout$name=="plum1",c("x","y")]<-c(-0,-0.7)

layout[layout$name=="darkgrey",c("x","y")]<-c(-1,-1)
layout[layout$name=="lightcyan",c("x","y")]<-c(1,-1)

layout[layout$name=="total_skin_score",c("x","y")]<-c(-0.5,-1.3)
layout[layout$name=="Skin_score_at_biopsy_site",c("x","y")]<-c(0.5,-1.6)
####################################################################################

ggraph(layout)+
  geom_edge_diagonal(arrow = arrow(length = unit(4, 'mm'),type = "closed"), 
                      end_cap = circle(10, 'mm'),color="#0086ab",
                      width=0.9)+
  #geom_node_point(aes(color=celltype_color),size=20)+
  geom_node_point(aes(color=module_color%>% str_replace_all("white","antiquewhite")),size=20)+
  geom_node_text(aes(label=name %>% str_replace_all("_"," ")),color="#404040")+
  scale_color_identity()+
  theme_graph()+
  theme(plot.margin= unit(c(0,0,0,0), "lines"),
        panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))#+
  xlim(c(-1.25,1.5))
savefile<-paste0("../2.Output/05_BayesianNetwork/Image/",str_sub(basename(BN_result_file),6,-5),"_modified4",".png")
ggsave(filename = savefile,width = 190, height = 190, units = "mm",bg="transparent")
####################################################################################

graph<-read_csv(BN_result_file,col_names = rownames(MergedData_1))%>% 
  mutate(from=rownames(MergedData_1)) %>% 
  gather(key=to,value=edge,-from) %>% 
  filter(edge==1) %>% 
  select(-edge) %>% 
  as_tbl_graph() %N>% 
  left_join(module2function,by=c("name"="module")) %>% 
  mutate(module_color=if_else(name %in% traits,"red",name),    #name %>% str_replace("Skin_score_at_biopsy_site","red") %>% str_replace("total_skin_score" ,"red"),
         celltype_large=if_else(celltype %in% immune_cells,"immune cell",if_else(celltype %in% fibroblast,"fibroblast","")),
         celltype_color=if_else(celltype_large=="immune cell","#f79646",if_else(celltype_large=="fibroblast","#9bbb59","grey")) )
ggraph(graph ,layout = 'sugiyama')+
  geom_edge_diagonal(arrow = arrow(length = unit(4, 'mm'),type = "closed"), 
                     end_cap = circle(10, 'mm'),color="#0086ab",
                     width=0.9)+
  #geom_node_point(aes(color=celltype_color),size=20)+
  geom_node_point(aes(color=module_color%>% str_replace_all("white","antiquewhite")),size=20)+
  geom_node_text(aes(label=name %>% str_replace_all("_"," ")),color="#404040")+
  scale_color_identity()+
  theme_graph()+
  theme(plot.margin= unit(c(0,0,0,0), "lines"),
        panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))#+
  xlim(c(-1.25,1.5))
  ####################################################################################
  #BMA
  #c2PB
  BMA_pre=read.table("../2.Output/05_BayesianNetwork/BN_SSc/result/BMA_pre_re.csv",sep=",",header = F,row.names = NULL,stringsAsFactors = F)
  data=read.table("../2.Output/05_BayesianNetwork/ME_SkinScore_forBN.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
  rownames(BMA_pre)=rownames(data)
  colnames(BMA_pre)=rownames(data)
  for (i in 1:nrow(BMA_pre)){
    for (j in 1:ncol(BMA_pre)){
      if (BMA_pre[i,j]>30 && BMA_pre[j,i]>30){
        if (BMA_pre[i,j]>BMA_pre[j,i]){
          BMA_pre[j,i]=0
        }else if (BMA_pre[i,j]<BMA_pre[j,i]){
          BMA_pre[i,j]=0
        }else{
          message(print("node",i,"and node",j,"is undirected"))
        }
      }
    }
  } #?Ō???eles?͂Ȃ?
  
  #???Ɗ??S?ɓƗ??????m?[?h???폜
  d_node=c()
  for (i in 1:ncol(BMA_pre)){
    if (sum(BMA_pre[i,])==0 && sum(BMA_pre[,i])==0){
      d_node=c(d_node,i)
    }
  } #d_node is empty 
  
  BMA=BMA_pre
  for (i in 1:nrow(BMA)){
    for (j in 1:ncol(BMA)){
      if (BMA[i,j]>1){
        BMA[i,j]=1
      }
    }
  }
  #DAG?̊m?F
  dag=empty.graph(rownames(BMA))
  amat(dag)=as.matrix(BMA) #the specified network contains cycles.
  
  for (i in 3:ncol(BMA)){
    print(i)
    dag=empty.graph(row.names(BMA)[1:i])
    amat(dag)=as.matrix(BMA[1:i,1:i])
  } #42 stop
  
  node=c()
  for (i in 1:ncol(BMA)){
    tryCatch({
      BMA_1=BMA[-i,-i]
      dag=empty.graph(rownames(BMA_1))
      amat(dag)=as.matrix(BMA_1)
      node=c(node,i)
      rm(BMA_1)
    },error=function(e){print(i)})
  }
  node #22 36 42
  
  in_36=which(BMA[,36]==1) # 12 35 39
  out_42=which(BMA[42,]==1) #5  7 10 12 20 28 34 35 41
  
  for (j in in_36){
    BMA_2=BMA[-j,-j]
    node=c()
    for (i in 1:dim(BMA_2)[1]){
      tryCatch({
        BMA_1=BMA_2[-i,-i]
        dag=empty.graph(rownames(BMA_1))
        amat(dag)=as.matrix(BMA_1)
        node=c(node,i)
        rm(BMA_1)
      },error=function(e){print(i)})
    }
    assign(paste0("node_",j),node)
    rm(BMA_2)
    rm(node)
  }
  
  BMA_pre[c(22,35,36,42),c(22,35,36,42)]
  BMA_pre[c(12,22,36,42),c(12,22,36,42)]
  
  #36->22????
  BMA_final=BMA
  BMA_final[36,22]=0
  dag=empty.graph(rownames(BMA_final))
  amat(dag)=as.matrix(BMA_final)
  #ok
  pdf("../2.Output/05_BayesianNetwork/BMA.pdf",paper="a4",width=9.5,height = 7,pointsize = 10)
  g=Rgraphviz::layoutGraph(bnlearn::as.graphNEL(dag))
  graph::parRenderInfo(g)=list(graph=list(main="c3PB BMA"),nodes=list(fontsize=30))
  Rgraphviz::renderGraph(g)
  dev.off()
  
  #???????킩?????̂ŉ??߂ł??郂?W???[?????I??
  #module-module????
  rm(list=ls())
  data=read.table("Superior temporal gyrus/exp_header_rownames.txt",sep="\t",header = T,row.names = 1,stringsAsFactors = F)
  ME=as.data.frame(t(data[-c(1:3),]))
  
  library(WGCNA)
  module_cor = cor(ME, use = "p",method ="pearson")
  #pvalue?Z?o
  nSamples=dim(ME)[1]
  module_score_Pvalue = corPvalueStudent(module_cor, nSamples)
  
  module_cor_1=module_cor
  for (i in 1:dim(module_score_Pvalue)[1]){
    for (j in 1:dim(module_score_Pvalue)[2]){
      if (module_score_Pvalue[i,j]>=0.05){
        module_cor_1[i,j]=0
      }
    }
  }
  
  for (i in 1:dim(module_cor_1)[1]){
    module_cor_1[i,i]=0
  }
  
  install.packages("sna")
  library(sna)
  gplot(module_cor_1,usearrows = F,displaylabels = T)
  






celltype_enrichment_WGCNAgenelist.R





rm(list = ls(all.names = TRUE))
library(WGCNA)
#setwd("~/01. R????/03.Analysis/01. SLE/human_SLE/WGCNA/Enrichment") #?t?@?C???i?[?p?X???w??
setwd("~/Project/20190709_SystemicSclerosis/3.Script")

#source("https://bioconductor.org/biocLite.R");
#biocLite(c("AnnotationDBI", "GO.db", "org.Hs.eg.db", "org.Mm.eg.db",
#           "XML", "WGCNA", "TxDb.Hsapiens.UCSC.hg19.knownGene", 
#           "TxDb.Mmusculus.UCSC.mm10.knownGene")); 


#source("https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/GeneAnnotation/installAnRichment.R");
#install.packages("anRichment_0.97-1.tar.gz", repos = NULL, type = "source")
#install.packages("anRichmentMethods_0.87-1.tar.gz", repos = NULL, type = "source")

#installAnRichment(); 

library(ggplot2)
library(WGCNA)


library(anRichment)
library(tidyr)
library(ggplot2)
options(stringsAsFactors = FALSE)
##############???͂Ɏg?p???????`?q?R???N?V???????쐬????
#?܂??֋X???̋?Collection??????
collection <- newCollection()

#Collection?̂??ƂɂȂ??e???`?q???X?g?i?????C?w?b?_?[?????ň??`?q??symbol?@entrez?̏ꍇ?͕ϊ??s?v?Ȃ̂ŉ??L?X?N???v?g?ύX)
#?t?H???_????Genelist???ꊇ?œǂݍ??܂??邽?߂̐ݒ?
#genelist?̃t?@?C?????????̂܂?enrichment???͎??̈??`?q?O???[?v???ɔ??f???????̂ŁC?K?v?ł????΂??ꂼ?ꃊ?l?[??
#setwd("~/01. R????/03.Analysis/01. SLE/human_SLE/WGCNA/Enrichment") #?t?@?C???i?[?p?X???w??
#files <- list.files(path = "../1.Data/Cell_Enrivhment/Skin+Blood/Jensen_TISSUES/",full.names = T)
#files_names <- list.files("../1.Data/Cell_Enrivhment/Skin+Blood/Jensen_TISSUES/" , full.names = F,pattern="csv") 
#files_names <-  gsub(".csv", "", files_names) 

files <- list.files(path = "C:/Users/ywt2100/Desktop/human",full.names = T)
files_names <- list.files("C:/Users/ywt2100/Desktop/human" , full.names = F, pattern="csv") 
files_names <- gsub(".csv", "", files_names) 



#genelist?̗R???ƂȂ?Organism???w?肷??
organism <- "Human"

#Genelist??Geneset?ɂ??邽?߂̊֐????ǂݍ??܂???
Togeneset <- function (X, Y){
  colnames(X) <- "Gene" #???????????I??Gene?ɕύX
  Symbol = X$Gene  #?֋X???̃??l?[??
  #????entrez?`???̈??`?q???????????ꍇ?́C???L??Entrez = X$Gene?ɂ??āC??2?s?̃X?N???v?g?폜
  Entrez.0 = convert2entrez(organism = organism, symbol = Symbol) 
  Entrez = unique(Entrez.0)#???`?q??(Symbol)??entrez?`???ɕϊ?????
  print(table(is.finite(Entrez)))  #entrez?`???ɕϊ??ł??????`?q???̊m?F
  newGeneSet(
    geneEntrez = Entrez,
    geneEvidence = "IEP",
    geneSource = "",
    ID = Y,
    name = Y,
    description = Y,
    source = "",
    organism = organism,
    internalClassification = Y,
    groups = Y,
    lastModified = "")
}

#?t?H???_????Genelist???S??Togeneset?֐???geneset?ɕϊ????Ă??????L??Collection?ɓ?????
for(i in 1:length(files)){
  genes <- read.csv(files[i], fileEncoding="UTF-8")#?????????t?@?C?????ǂݍ???
  genes2<- genes$Symbol
  assign(files_names[i], genes2)
  genes <- unique(genes)#?d???L???ꍇ????????
  genes_set <- Togeneset(genes, files_names[i])#genelist??geneset??
  assign(paste(files_names[i]), genes_set)
  collection <- addToCollection(collection, genes_set)#Geneset?????ꂽ??collection?ɓ?????
  #?????L??"Human_CNS_5cells_mouse_collection"?͎????̂??肽??collection?????ݒ?
  rm(genes, genes_set)
}

#doublecollecion <- newCollection()
#doublecollecion <- addToCollection(doublecollecion, Human_CNS_5cells_collection,Human_INF_36cells_collection_N)

##?K?v?ł????Ώ??L?ō쐬?????R???N?V???????ۑ??B
#setwd("C:/Users/mjd9761/Documents/Enrichment/celltype_enrichment/Genecollectiondata")
save(collection, file="Human_immune_test.Rdata")
#???x?R???N?V???????????΁C?????????̓??[?h???Ă??ł??g????
#load("Human_CNS_5cells_collection.Rdata")

######?I?v?V????
#???????R???N?V?????????????????????^?R???N?V???????쐬???āC?ꊇ?ŃG?????b?`?????g???͂??\?B
#CNS?̈??`?q???X?g?Ɖ??ǍזE?̈??`?q???X?g???ʁX?ō쐬???????ǁC?????ɃG?????b?`?????g???͂??????Ƃ??ȂǁB
#?ȉ??????̗?
#CNS_INF_collecion <- newCollection()
#CNS_INF_collecion <- addToCollection(doublecollecion, Human_CNS_5cells_collection,Human_INF_36cells_collection_N)


#enrichment???͂Ɏg?p????collection?Ɋ܂܂??????`?q???X?g?????o?^????
#???Lenrichment?֐??ɂ????āCTOP?????܂ł?p-value???ʂ??Z?o???邩?ɉe???B?C?ӂŕύX??
nlist_collection <- length(collection)

#CelltypeEnrichment?֐????ǂݍ??܂???
CelltypeEnrichment <- function (X){
  names(X) <- c("Gene", "Group") #???????????I??Gene??Group?ɕύX
  #symbol = X$Gene  #?֋X???̃??l?[??
  X<-X %>% filter(!is.na(Gene))
  entrez=X$Gene
  Group = X$Group?@#?֋X???̃??l?[??
  print(table(Group))  #???`?q???X?g?̐??̊m?F?p?B
  #entrez = convert2entrez(organism = "Human", symbol = symbol)?@#???`?q??(Symbol)??entrez?`???ɕϊ?????
  print(table(is.finite(entrez)))  #entrez?`???ɕϊ??ł??????`?q???̊m?F
  #????????anRichment?p?b?P?[?W?̒???enrichmentAnalysis?֐????g?p
  analysisresult = enrichmentAnalysis(?@?@
    classLabels = Group, identifiers = entrez,?@#classLabels???????̂??̂??????̉??͑Ώۈ??`?q?Q?Ƃ??ĔF???B
    refCollection = collection,  #reference?Ƃ??Đݒ肷???f?[?^?Z?b?g?B
    useBackground = "given",?@?@#???̓o?b?N?O???E???h?̐ݒ??i?d?v?j?B?????̓C???v?b?g?????S???`?q?̂????ǂݍ??߂????̂??g?p?B
    threshold = 1,?@?@#?G?????b?`?????g???͌??ʂ??Ƃ??ďo?͂?????csv?t?@?C???ł́Cp?l?̂??????l?B?????͍L???Ƃ??Ă????B
    thresholdType = "Bonferroni",?@#???Lp?l?ɕt?????āCBonferroni?␳????p?l??臒l?Ƃ???
    getOverlapEntrez = FALSE,?@?@#?o??csv?t?@?C???ɂāC?I?[?o?[???b?v???????`?q????entrez?`???ŏo?͂????C
    getOverlapSymbols = TRUE,    #?o??csv?t?@?C???ɂāC?I?[?o?[???b?v???????`?q????symbol?`???ŏo?͂????C
    maxReportedOverlapGenes = 10000,?@?@#???L?I?[?o?[???b?v???????`?q???ǂ̂??炢?\???????邩???ݒ??B?S?Ăق????̂?10,000??
    removeDuplicatesInDifferentClasses =FALSE,?@#???????W???[?????ɓ??????`?q?????????ꍇ???C???̂܂܏???????
    entrySeparator = ",",?@?@#?I?[?o?[???b?v???????`?q?Q?ɂ??āC","?ŕ????ďo?͂??Ă????B?ǂ??ł??悵
    ignoreLabels = "grey", #grey???W???[???̈??`?q?Q?ɂ??Ă͖????i???͂??Ȃ??j
    combineEnrichmentTables = FALSE) 
  
  
  countsInDataSet<- analysisresult$effectiveBackgroundSize  #???̓o?b?N?O???E???h?̈??`?q?????m??
  print(table(countsInDataSet))?@?@?@?@?@?@?@?@?@?@?@?@?@?@?@?@#???????o??
  
  
  Resulttable <- analysisresult$enrichmentTable
  Resulttable<-separate(Resulttable,overlapGenes, into=as.character(c(1:1000)), sep=",")
  Resulttable2 <- subset(analysisresult$enrichmentTable, analysisresult$enrichmentTable$pValue < 0.01) #?o?͗p?ɐ??`
  list_color=unique(Group)
  list_color_1=list_color[list_color!="grey"]
  if(length(list_color_1)!=length(unique(Resulttable2$class))){
    d=setdiff(list_color_1,unique(Resulttable2$class))
    d=as.data.frame(d)
    for (i in 1:length(d$d)){
      l=data.frame("class"=d$d[i],"rank"=1,"dataSetID"="n.d.","dataSetName"="n.d.","inGroups"="n.d.","pValue"=1,
                   "Bonferroni"="n.d.","FDR"="n.d.","nCommonGenes"=0,"fracOfEffectiveClassSize"="n.d.","expectedFracOfEffectiveClassSize"="n.d.",
                   "enrichmentRatio"="n.d.","classSize"="n.d.","effectiveClassSize"="n.d.","fracOfEffectiveSetSize"="n.d.","effectiveSetSize"="n.d.",
                   "shortDataSetName"="n.d.","overlapGenes"="n.d.")
      Resulttable2=rbind(Resulttable2,l)
    }
    write.table(Resulttable2, file = "Human_PS19_comparison_result.csv",row.names = FALSE,sep=",")#???ʂ?csv?t?@?C???ɏo??
  }else{
    write.table(Resulttable2, file = "Human_PS19_comparison_result.csv",row.names = FALSE,sep=",")#???ʂ?csv?t?@?C???ɏo??
  }
  options(warn=-1) #???̃R?[?h??warning???o???B???????????????R?[?h
  #Overlap???`?q?Q???C1???`?q???ƂɃZ???ɕ????ĕ\???????邽?߂̃R?[?h(1000?͔C?ӁB)???̍??Ƃ?warning???o?邪???Ə??͖????Ȃ?
  list <- by(Resulttable, Resulttable$class, data.frame) #?O???[?v???ƂɃt?@?C???`???????X?g??
  sapply(1:dim(list), function(x){write.csv(list[x], file=paste0("CellResult_", dimnames(list)[[1]][x], ".csv"), row.names=FALSE)})?@#?O???[?v???ƂɌ??ʂ??o??
  
  #?????????C???͌??ʃf?[?^(Resulttable)?̂????O???t?ɕK?v?Ȃ??̂??Ԃ????ʂ?
  Rank_all <- sapply(1:dim(list), function(x){list(list[[x]][1:nlist_collection,c(4,6,9)])})#?O???[?v???ƂɕK?v?ȗ????????o??
  names(Rank_all) <- sapply(1:dim(list), function(x){names(Rank_all) <- names(list[x])})#???L?ŃO???[?v???????????̂ł??????????l?[??
  for (i in 1:dim(list)) Rank_all[[i]]$nCommonGenes <- paste("(",Rank_all[[i]]$nCommonGenes,")") #?O???t?̃??x???p?̖??O???`
  for (i in 1:dim(list)) Rank_all[[i]]$pValue <- -(log10(as.numeric(Rank_all[[i]]$pValue)))?@#?O???t?̃??x???p?̖??O???`
  Rank_all <- na.omit(Rank_all) #???ʂ?TOP10?ɖ????Ȃ??ꍇ??NA????
  for (i in 1:dim(list)) Rank_all[[i]] <- transform(Rank_all[[i]], "Rename"=(paste(Rank_all[[i]]$dataSetName,Rank_all[[i]]$nCommonGenes,sep="?@")))#?O???t?̃??x???p?̖??O???`
  
  #list Rename
  names(Rank_all) <- sapply(1:dim(list), function(x){names(Rank_all) <- names(list[x])})
  
  
  ##????????ggplot2???g?????}???̎w??
  #for (i in 1:dim(list)) ggsave(file=paste0("CellResult_", names(Rank_all)[i], ".pdf"),plot = ggplot(Rank_all[[i]], aes(x=reorder(Rank_all[[i]]$Rename, Rank_all[[i]]$pValue), y=Rank_all[[i]]$pValue)) +  
  #                                geom_bar(stat="identity", width=.5,fill="black") +   
  #                                coord_flip() +                                     
  #                                xlab("Cell-type\n(#Overlap genes)") + 
  #                                ylab("-log(P-value)"))
  
  #Cell-type enrichment p-value?̃q?[?g?}?b?v???}??
  #CellTyperesult <- read.csv("WGCNA_Cellenrichment_result.csv")
  #logP <- -log10(CellTyperesult$pValue)
  #CellTypeP <- data.frame(Module = CellTyperesult$class, CellType = CellTyperesult$dataSetID, logP = logP)
  #ghm <- ggplot(CellTypeP, aes(x = CellType, y = Module, fill = logP))
  #ghm <- ghm + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  #ghm <- ghm + scale_x_discrete(limits = c("astrocyte", "endothelial", "microglia", "neuron", "oligodendrocyte"))
  #ghm <- ghm + geom_tile(aes(fill = logP))
  #ghm <- ghm + xlab("Cell Type") + ylab("Module")
  #ghm <- ghm + geom_text(aes(label = round(CellTypeP$logP, 1)), size = 2)
  #ghm <- ghm + scale_fill_gradient(low = "white", high = "red", limits = c(0, 200))
  #pdf("CellType_pvalue_heatmap.pdf", width = 5, height = 7)
  #plot(ghm)
  #dev.off()
  
}
###################?֐??????܂?###################

#???ۂɉ???
#setwd("C:/Users/mjd9761/Documents/Enrichment/celltype_enrichment") 

#WGCNAgenelist <- read.table("WGCNAgenelist.txt", header = TRUE)
load("networkConstruction_StepByStep_unsigned_0.Rdata")
lnames <- load(paste("unsigned_",deepSplit,"/networkConstruction_StepByStep_unsigned_", deepSplit, "_merged.Rdata", sep = ""))
#Probe ID??Entrez Gene ID?ɕϊ?????
datExpr <- Expdata_quontile$E
keep<-rowSums(datExpr>log2(50))>=102
datExpr<-datExpr[keep,]%>% t()
probes <- colnames(datExpr)
annot<-read.csv("../1.Data/Annotation/GeneAnnotation_GPL10558.csv")
probes2annot <- match(probes, annot$ID)
#Entrez Gene ID?iLocusLinkID?j???ǂݍ???
allLLIDs <- annot$Entrez_Gene_ID[probes2annot]
#WGCNA?ɂ??????`?q?ƃ??W???[???̕\???쐬
WGCNAgenelist <- data.frame(allLLIDs, moduleColors_unsigned)
#WGCNAgenelist<-
CelltypeEnrichment(WGCNAgenelist)








decompressData.R





setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
input_dir<-"../1.Data/GSE9285/DecompressedData/All/"
output_dir<-"../1.Data/GSE9285/DecompressedData/"
files<-list.files(input_dir,full.names = T,recursive = T)
file<-files[1]
for (file in files){
  file.copy(file,paste0(output_dir,basename(file)))
}






DrawGOEnrichment_Barplot.R





setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
library(WGCNA)
library(anRichment)
library(tidyverse)
options(stringsAsFactors = FALSE)

output_file_base_path<-"../2.Output/02_signed_2_unmerged/"

lnames <- load("../2.Output/01_Compare_deepsplit_signed/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/Expdata_quontile.Rdata")
annot <- read.csv("../1.Data/Annotation/GeneAnnotation_GPL10558.csv", header = T, stringsAsFactors = F)


go_file_path<-paste(dir_path,"GOEnrichmentTable_",signed_unsigned,"_",deepSplit,merged_unmerged,"withBP",".csv",sep="")
go_data<-read_csv(go_file_path)
#go_data$termName<-str_wrap(go_data$termName, width = 40)

go_data %>% mutate()
tmp<-go_data%>% distinct(module) 
tmp<-tmp%>% mutate(index=row.names(tmp))%>% 
  mutate(col= if_else(as.double(index) < (nrow(tmp)/2), true = 1, false = 2))

ylim<-max(-log10(go_data$BonferoniP)) +1

go_data %>% 
  filter(module %in% significant_module) %>% 
  filter(BonferoniP<=0.05) %>% 
  mutate(termName=termName %>% str_replace_all("endoplasmic reticulum","ER"),
         termName=str_wrap(termName, width =45 ),
         color=module %>% 
           str_replace_all("lightyellow","yellow") %>% 
           str_replace_all("lightcyan","cyan")) %>% 
  group_by(module )%>%
  top_n(n = 5, wt = -BonferoniP ) %>% 
  ggplot(aes(x=reorder(termName,BonferoniP),y=-log10(BonferoniP)))+
  geom_bar(aes(fill=color),stat = "identity",width = 0.6)+
  facet_grid(.~module,  scales='free_x')+
  labs(y="-log(Bonferoni p-value)",x="")+
  theme_classic()+
  theme(
    axis.text=element_text(size=8),#size=text_size),
    axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")  )+
  scale_y_continuous(expand = c(0, 0),limits =c(0,27) )+
  scale_fill_identity()+
  #coord_flip()+
  guides(fill=FALSE)
savefile<-paste0(output_file_base_path,"GOEnrichment_withBP_ver7",".png")
ggsave(file = savefile, dpi = 320, width = 330, height = 140, units = "mm",  bg = "transparent")

#####yellow turquoise####################################
go_data %>% 
  filter(module=="yellow") %>% 
  mutate(color=module) %>% 
  ggplot(aes(x=reorder(termName,enrichmentP),y=-log10(enrichmentP)))+
  geom_bar(aes(fill=color),stat = "identity",width = 0.6)+
  facet_grid(.~module,  scales='free_x')+
  labs(y="-log(p-value)",x="")+
  theme_classic()+
  theme(
    axis.text=element_text(size=8),#size=text_size),
    axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")  )+
  scale_y_continuous(expand = c(0, 0),limits =c(0,10) )+
  scale_fill_identity()+
  #coord_flip()+
  guides(fill=FALSE)
savefile<-paste0(output_file_base_path,"GOEnrichment_withBP_yellow_ver1",".png")
ggsave(file = savefile, dpi = 320, width = 50, height = 200, units = "mm",  bg = "transparent")

go_data %>% 
  filter(module=="turquoise") %>% 
  mutate(color=module) %>% 
  ggplot(aes(x=reorder(termName,enrichmentP),y=-log10(enrichmentP)))+
  geom_bar(aes(fill=color),stat = "identity",width = 0.6)+
  facet_grid(.~module,  scales='free_x')+
  labs(y="-log(p-value)",x="")+
  theme_classic()+
  theme(
    axis.text=element_text(size=14),#size=text_size),
    axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")  )+
  scale_y_continuous(expand = c(0, 0),limits =c(0,10) )+
  scale_fill_identity()+
  #coord_flip()+
  guides(fill=FALSE)

savefile<-paste0(output_file_base_path,"GOEnrichment_withBP_turquoise_ver2",".png")
ggsave(file = savefile, dpi = 320, width = 100, height = 150, units = "mm",  bg = "transparent")

go_data %>% 
  filter(module=="yellow") %>% 
  mutate(color=module %>% str_replace_all("yellow","gold")) %>% 
  ggplot(aes(x=reorder(termName,-enrichmentP),y=-log10(enrichmentP)))+
  geom_bar(aes(fill=color),stat = "identity",width = 0.6)+
  facet_grid(.~module,  scales='free_x')+
  labs(y="-log(p-value)",x="")+
  theme_classic()+
  theme(
    axis.text=element_text(size=14),#size=text_size),
    axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")  )+
  scale_y_continuous(expand = c(0, 0),limits =c(0,10) )+
  scale_fill_identity()+
  coord_flip()+
  guides(fill=FALSE)

savefile<-paste0(output_file_base_path,"GOEnrichment_withBP_yellow_ver2",".png")
ggsave(file = savefile, dpi = 320, width = 150, height = 100, units = "mm",  bg = "transparent")

go_data %>% 
  filter(module=="turquoise") %>% 
  mutate(color=module) %>% 
  ggplot(aes(x=reorder(termName,-enrichmentP),y=-log10(enrichmentP)))+
  geom_bar(aes(fill=color),stat = "identity",width = 0.6)+
  facet_grid(.~module,  scales='free_x')+
  labs(y="-log(p-value)",x="")+
  theme_classic()+
  theme(
    axis.text=element_text(size=14),#size=text_size),
    axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = "transparent"), # get rid of legend bg
    legend.box.background = element_rect(fill = "transparent")  )+
  scale_y_continuous(expand = c(0, 0),limits =c(0,10) )+
  scale_fill_identity()+
  coord_flip()+
  guides(fill=FALSE)

savefile<-paste0(output_file_base_path,"GOEnrichment_withBP_turquoise_ver3",".png")
ggsave(file = savefile, dpi = 320, width = 150, height = 100, units = "mm",  bg = "transparent")






ExploratoryDataAnalysis.R





library(tidyverse)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
sample_annotation_file_path<-"../1.Data/GSE58095_all_sample_annotations.txt"
output_file_base_path<-"../2.Output/06_ExploratoryDataAnalysis/"
#annotation<-read_tsv("../1.Data/GSE58095_all_sample_annotations.txt",skip = 4,n_max = 102)%>% 
#  mutate(`Characteristics: Group`=`Characteristics: Group` %>% replace_na("SSc")) 
annotation <- read_tsv(sample_annotation_file_path,skip = 4,n_max = 102) %>% 
  rename_all(function(x) x %>% str_replace_all(" ","_")%>% 
               str_replace_all(":","_") %>% 
               str_replace_all("\\(","") %>% 
               str_replace_all("\\)","") %>% 
               str_replace_all("=","is"))%>% 
  dplyr::select(-c(title,Characteristics__Patient,Characteristics__GSE47162_Sample_name,Characteristics__Time_point,
                   molecule,label,description,platform)) %>% 
  mutate(Characteristics__Group=Characteristics__Group %>% replace_na("SSc"))

draw_piechart<-function(columnname){
  annotation %>% 
    group_by(columnname) %>% 
    summarise(n=n())%>% 
    mutate(columnname%>%replace_na("NA"),
           columnname%>% as.character()%>% factor(levels=c("NA","1","0"))) %>% 
    arrange(-n) %>% 
    mutate(prop=n/nrow(annotation)*100,
           lab.ypos = cumsum(prop) - 0.5*prop) %>% 
    ggplot(aes(x="",y=prop))+
    geom_bar(aes(fill=`Characteristics: Gender (1=male)` ) ,stat="identity",width = 1, size = 1,  color = "white")+
    geom_text(aes(y=lab.ypos,label = n), color = "white",size=8)+
    scale_fill_manual(values = c("#777777","#f79646","#0086ab","#da6272","#9bbb59","#bfbfbf"),
                      labels=c("NA","Male","Female"))+
    guides(fill=guide_legend(reverse = TRUE))+
    coord_polar("y", start=0,direction = 1)+
    theme_void()+
    theme(legend.title = element_blank(),
          axis.text.x=element_blank(),
          axis.title=element_text(size=20),
          panel.background = element_rect(fill = "transparent",color = NA),
          panel.grid.minor = element_line(color = NA), 
          panel.grid.major = element_line(color = NA),
          plot.background = element_rect(fill = "transparent",color = NA),
          plot.title = element_text(hjust = 0.5,vjust = 0.5))+
    labs(y="Gender",x="")
}
draw_piechart(columnname=`Characteristics: Gender (1=male)`)

###################################################################
#gender

annotation %>% 
  group_by(`Characteristics: Gender (1=male)`) %>% 
  summarise(n=n())%>% 
  mutate(`Characteristics: Gender (1=male)`=`Characteristics: Gender (1=male)` %>%replace_na("NA"),
         `Characteristics: Gender (1=male)`=`Characteristics: Gender (1=male)` %>% as.character()%>% factor(levels=c("NA","1","0"))) %>% 
  arrange(-n) %>% 
  mutate(prop=n/nrow(annotation)*100,
         lab.ypos = cumsum(prop) - 0.5*prop) %>% 
  ggplot(aes(x="",y=prop))+
  geom_bar(aes(fill=`Characteristics: Gender (1=male)` ) ,stat="identity",width = 1, size = 1,  color = "white")+
  geom_text(aes(y=lab.ypos,label = n), color = "white",size=8)+
  scale_fill_manual(values = c("#777777","#f79646","#0086ab","#da6272","#9bbb59","#bfbfbf"),
                    labels=c("NA","Male","Female"))+
  guides(fill=guide_legend(reverse = TRUE))+
  coord_polar("y", start=0,direction = 1)+
  theme_void()+
  theme(legend.title = element_blank(),
        axis.text.x=element_blank(),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = "transparent",color = NA),
        panel.grid.minor = element_line(color = NA), 
        panel.grid.major = element_line(color = NA),
        plot.background = element_rect(fill = "transparent",color = NA),
        plot.title = element_text(hjust = 0.5,vjust = 0.5))+
  labs(y="Gender",x="")

savefile<-paste0(output_file_base_path,"Gender_proportion",".png")
ggsave(file = savefile, dpi = 320, width = 140, height = 140, units = "mm",  bg = "transparent")
###########################################################################################
color_value=c("#0086ab","#f79646","#777777","#0086ab","#da6272","#9bbb59","#bfbfbf")
columnname="Group"
annotation %>% 
  group_by(`Characteristics: Group`) %>% 
  summarise(n=n())%>% 
  mutate(`Characteristics: Group`=`Characteristics: Group` %>%replace_na("NA")) %>% 
  arrange(-n)->tmp
tmp%>% 
  mutate(prop=n/nrow(annotation)*100,
         lab.ypos =cumsum(prop) - 0.5*prop,
         #lab.ypos =rev(lab.ypos ),
         color_value=color_value[1:nrow(tmp)]) %>% 
  ggplot(aes(x="",y=rev(prop)))+
  geom_bar(aes(fill=color_value ) ,stat="identity",width = 1, size = 1,  color = "white")+
  geom_text(aes(y=lab.ypos,label = n), color = "white",size=8)+
  
  guides(fill=guide_legend(reverse = TRUE))+
  coord_polar("y", start=0,direction = 1)+
  theme_void()+
  theme(legend.title = element_blank(),
        axis.text.x=element_blank(),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = "transparent",color = NA),
        panel.grid.minor = element_line(color = NA), 
        panel.grid.major = element_line(color = NA),
        plot.background = element_rect(fill = "transparent",color = NA),
        plot.title = element_text(hjust = 0.5,vjust = 0.5))+
  labs(y=columnname,x="")+
  scale_fill_identity(labels=rev(tmp$`Characteristics: Group`),guide = "legend")

savefile<-paste0(output_file_base_path,columnname,"_proportion",".png")
ggsave(file = savefile, dpi = 320, width = 140, height = 140, units = "mm",  bg = "transparent")
###################################################################
color_value=c("#0086ab","#f79646","#777777","#0086ab","#da6272","#9bbb59","#bfbfbf")
annotation %>% 
  group_by(`Characteristics: Group`) %>% 
  summarise(n=n())%>% 
  mutate(`Characteristics: Group`=`Characteristics: Group` %>%replace_na("NA")) %>% 
  arrange(-n) %>% 
  mutate(rev_n=rev(n),
         prop=n/nrow(annotation)*100,
         lab.ypos =cumsum(rev(prop)) - 0.5*rev(prop),
         #lab.ypos =rev(lab.ypos ),
         color_value=color_value[1:nrow(tmp)]) ->tmp
tmp
tmp%>% 
  ggplot(aes(x="",y=prop))+
  geom_bar(aes(fill=color_value ) ,stat="identity",width = 1, size = 1,  color = "white")+
  geom_text(aes(y=lab.ypos,label = rev_n), color = "white",size=8)+
  
  guides(fill=guide_legend(reverse = TRUE))+
  coord_polar("y", start=0,direction = 1)+
  theme_void()+
  theme(legend.title = element_blank(),
        axis.text.x=element_blank(),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = "transparent",color = NA),
        panel.grid.minor = element_line(color = NA), 
        panel.grid.major = element_line(color = NA),
        plot.background = element_rect(fill = "transparent",color = NA),
        plot.title = element_text(hjust = 0.5,vjust = 0.5))+
  labs(y=columnname,x="")+
  scale_fill_identity(labels=tmp$`Characteristics: Group`,guide = "legend")
savefile<-paste0(output_file_base_path,columnname,"_proportion",".png")
ggsave(file = savefile, dpi = 320, width = 140, height = 140, units = "mm",  bg = "transparent")
##########################################################################
color_value=c("#da6272","#0086ab","#f79646","#777777","#0086ab","#9bbb59","#bfbfbf")
columnname="Characteristics__Group"
annotation %>% 
  group_by(`Characteristics: Group`) %>% 
  summarise(n=n())%>% 
  mutate(`Characteristics: Group`=`Characteristics: Group` %>%replace_na("NA")) %>% 
  arrange(-n) %>% 
  mutate(rev_n=rev(n),
         prop=n/nrow(annotation)*100,
         lab.ypos =cumsum(rev(prop)) - 0.5*rev(prop),
         #lab.ypos =rev(lab.ypos ),
         color_value=color_value[1:nrow(tmp)]) ->tmp
tmp
tmp%>% 
  ggplot(aes(x="",y=prop))+
  geom_bar(aes(fill=reorder(color_value,rev_n) ) ,stat="identity",width = 1, size = 1,  color = "white")+
  geom_text(aes(y=lab.ypos,label = str_c(rev(`Characteristics: Group`),"\n",rev_n)), color = "white",size=8)+
  
  guides(fill=guide_legend(reverse = TRUE))+
  coord_polar("y", start=0,direction = 1)+
  theme_void()+
  theme(legend.title = element_blank(),
        axis.text.x=element_blank(),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = "transparent",color = NA),
        panel.grid.minor = element_line(color = NA), 
        panel.grid.major = element_line(color = NA),
        plot.background = element_rect(fill = "transparent",color = NA),
        plot.title = element_text(hjust = 0.5,vjust = 0.5))+
  labs(y=columnname,x="")+
  scale_fill_identity(labels=tmp$`Characteristics: Group`,guide = "legend")

savefile<-paste0(output_file_base_path,columnname,"_proportion2",".png")
ggsave(file = savefile, dpi = 320, width = 140, height = 140, units = "mm",  bg = "transparent")

###################################################################
color_value=c("#da6272","#0086ab","#f79646","#777777","#0086ab","#9bbb59","#bfbfbf")
columnname="Characteristics__Group"
columnname<-enquo(columnname)
annotation %>% 
  group_by(!!columnname) %>% 
  summarise(n=n())%>% 
  mutate(target=!!columnname%>%replace_na("NA")) %>% 
  arrange(-n) ->tmp
tmp%>% 
  mutate(rev_n=rev(n),
         prop=n/nrow(annotation)*100,
         lab.ypos =cumsum(rev(prop)) - 0.5*rev(prop),
         #lab.ypos =rev(lab.ypos ),
         color_value=color_value[1:nrow(tmp)]) ->tmp
tmp
tmp%>% 
  ggplot(aes(x="",y=prop))+
  geom_bar(aes(fill=reorder(color_value,rev_n) ) ,stat="identity",width = 1, size = 1,  color = "white")+
  geom_text(aes(y=lab.ypos,label = str_c(rev(`Characteristics: Group`),"\n",rev_n)), color = "white",size=8)+
  
  guides(fill=guide_legend(reverse = TRUE))+
  coord_polar("y", start=0,direction = 1)+
  theme_void()+
  theme(legend.title = element_blank(),
        axis.text.x=element_blank(),
        axis.title=element_text(size=20),
        panel.background = element_rect(fill = "transparent",color = NA),
        panel.grid.minor = element_line(color = NA), 
        panel.grid.major = element_line(color = NA),
        plot.background = element_rect(fill = "transparent",color = NA),
        plot.title = element_text(hjust = 0.5,vjust = 0.5))+
  labs(y=columnname,x="")+
  scale_fill_identity(labels=tmp$`Characteristics: Group`,guide = "legend")

savefile<-paste0(output_file_base_path,columnname,"_proportion2",".png")
ggsave(file = savefile, dpi = 320, width = 140, height = 140, units = "mm",  bg = "transparent")

columnname="Characteristics__Group"
draw_piechart<-function(annotation,columnname){
  columnname<-enquo(columnname)
  print(columnname)
  annotation %>% 
    group_by(!!columnname) %>% 
    summarise(n=n())%>% 
    arrange(-n) ->tmp
  tmp%>% 
    mutate(rev_n=rev(n),
           prop=n/nrow(annotation)*100,
           lab.ypos =cumsum(rev(prop)) - 0.5*rev(prop),
           #lab.ypos =rev(lab.ypos ),
           color_value=color_value[1:nrow(tmp)]) ->tmp
  tmp%>% 
    ggplot(aes(x="",y=prop))+
    geom_bar(aes(fill=reorder(color_value,rev_n) ) ,stat="identity",width = 1, size = 1,  color = "white")+
    geom_text(aes(y=lab.ypos,label = str_c(rev(!!columnname),"\n",rev_n)), color = "white",size=8)+
    
    guides(fill=guide_legend(reverse = TRUE))+
    coord_polar("y", start=0,direction = 1)+
    theme_void()+
    theme(legend.title = element_blank(),
          axis.text.x=element_blank(),
          axis.title=element_text(size=20),
          panel.background = element_rect(fill = "transparent",color = NA),
          panel.grid.minor = element_line(color = NA), 
          panel.grid.major = element_line(color = NA),
          plot.background = element_rect(fill = "transparent",color = NA),
          plot.title = element_text(hjust = 0.5,vjust = 0.5))+
    labs(y=columnname,x="")+
    scale_fill_identity()->g
  print(tmp)
  return(g)
}
draw_piechart(annotation,Characteristics__Group)
draw_piechart(annotation %>% 
       mutate(Characteristics__Diffuse_vs_limited=Characteristics__Diffuse_vs_limited %>% replace_na("Healthy")),Characteristics__Diffuse_vs_limited)
savefile<-paste0(output_file_base_path,"Diffuse_vs_limited","_proportion",".png")
ggsave(file = savefile, dpi = 320, width = 140, height = 140, units = "mm",  bg = "transparent")






ExploratoryDataAnalysis2.R





library(tidyverse)
library(polycor)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
sample_annotation_file_path<-"../1.Data/GSE58095_all_sample_annotations.txt"
output_file_base_path<-"../2.Output/06_ExploratoryDataAnalysis/"
#annotation<-read_tsv("../1.Data/GSE58095_all_sample_annotations.txt",skip = 4,n_max = 102)%>% 
#  mutate(`Characteristics: Group`=`Characteristics: Group` %>% replace_na("SSc")) 
annotation <- read_tsv(sample_annotation_file_path,skip = 4,n_max = 102) %>% 
  rename_all(function(x) x %>% str_replace_all(" ","_")%>% 
               str_replace_all(":","_") %>% 
               str_replace_all("\\(","") %>% 
               str_replace_all("\\)","") %>% 
               str_replace_all("=","is"))%>% 
  dplyr::select(-c(title,Characteristics__Patient,Characteristics__GSE47162_Sample_name,Characteristics__Time_point,
                   molecule,label,description,platform)) %>% 
  mutate(Characteristics__Group=Characteristics__Group %>% replace_na("SSc"))
color_value=c("#0086ab","#f79646","#777777","#0086ab","#da6272","#9bbb59","#bfbfbf")
draw_piechart<-function(annotation,columnname,color_value){
  columnname<-enquo(columnname)
  print(columnname)
  annotation %>% 
    group_by(!!columnname) %>% 
    summarise(n=n())%>% 
    arrange(-n) ->tmp
  tmp%>% 
    mutate(rev_n=rev(n),
           prop=n/nrow(annotation)*100,
           lab.ypos =cumsum(rev(prop)) - 0.5*rev(prop),
           #lab.ypos =rev(lab.ypos ),
           color_value=color_value[1:nrow(tmp)]) ->tmp
  tmp%>% 
    ggplot(aes(x="",y=prop))+
    geom_bar(aes(fill=reorder(color_value,rev_n) ) ,stat="identity",width = 1, size = 1,  color = "white")+
    geom_text(aes(y=lab.ypos,label = str_c(rev(!!columnname),"\n",rev_n)), color = "white",size=8)+
    
    guides(fill=guide_legend(reverse = TRUE))+
    coord_polar("y", start=0,direction = 1)+
    theme_void()+
    theme(legend.title = element_blank(),
          axis.text.x=element_blank(),
          axis.title=element_text(size=20),
          panel.background = element_rect(fill = "transparent",color = NA),
          panel.grid.minor = element_line(color = NA), 
          panel.grid.major = element_line(color = NA),
          plot.background = element_rect(fill = "transparent",color = NA),
          plot.title = element_text(hjust = 0.5,vjust = 0.5))+
    labs(y=columnname,x="")+
    scale_fill_identity()->g
  print(tmp)
  return(g)
}
annotation %>% 
  mutate(Characteristics__Gender_1ismale=Characteristics__Gender_1ismale %>%as.logical() %>%  if_else("Male","Female")%>% replace_na("NA")) %>% 
  draw_piechart(Characteristics__Gender_1ismale,
                c("#da6272","#0086ab","#777777","#f79646","#0086ab","#9bbb59","#bfbfbf"))
savefile<-paste0(output_file_base_path,"Gender","_proportion",".png")
ggsave(file = savefile, dpi = 320, width = 140, height = 140, units = "mm",  bg = "transparent")
draw_piechart(annotation,Characteristics__Group,color_value)
annotation %>% 
  mutate(Characteristics__Diffuse_vs_limited=Characteristics__Diffuse_vs_limited %>% replace_na("Healthy")) %>% 
  draw_piechart(Characteristics__Diffuse_vs_limited,
                c("#da6272","#0086ab","#f79646","#777777","#0086ab","#9bbb59","#bfbfbf"))
savefile<-paste0(output_file_base_path,"Diffuse_vs_limited","_proportion",".png")
ggsave(file = savefile, dpi = 320, width = 140, height = 140, units = "mm",  bg = "transparent")

annotation %>% 
  mutate(characteristics__FVC_less_than_70=characteristics__FVC_less_than_70 %>% replace_na("NA")) %>% 
  draw_piechart(characteristics__FVC_less_than_70,
                c("#777777","#0086ab","#da6272","#f79646","#0086ab","#9bbb59","#bfbfbf"))
savefile<-paste0(output_file_base_path,"FVC","_proportion",".png")
ggsave(file = savefile, dpi = 320, width = 140, height = 140, units = "mm",  bg = "transparent")

annotation %>% 
  mutate(Characteristics__Diffuse_vs_limited=Characteristics__Diffuse_vs_limited %>% replace_na("Healthy")) %>% 
  draw_piechart(Characteristics__Group,
                c("#da6272","#0086ab","#f79646","#777777","#0086ab","#9bbb59","#bfbfbf"))
savefile<-paste0(output_file_base_path,"Proportion_","Group",".png")
ggsave(file = savefile, dpi = 320, width = 140, height = 140, units = "mm",  bg = "transparent")
####################################################################################################
#correlation
color_values=c("#0086ab","#da6272","#777777","#f79646","#0086ab","#9bbb59","#bfbfbf")

annotation %>% .$characteristics__FVC_less_than_70
  
annotation %>% 
  filter(Characteristics__Group=="SSc",!is.na(characteristics__FVC_less_than_70)) %>% 
  mutate(Characteristics__total_skin_score=replace_na(Characteristics__total_skin_score,0),
         Characteristics__Skin_score_at_biopsy_site=replace_na(Characteristics__Skin_score_at_biopsy_site,0)) %>% #.$ILD
  ggplot(aes(x=characteristics__FVC_less_than_70,y=Characteristics__total_skin_score))+
  geom_boxplot(aes(fill=characteristics__FVC_less_than_70),alpha=0.5)+
  geom_point(aes(color=characteristics__FVC_less_than_70),position =position_jitter(width=0.1))+
  theme_classic()+
  theme(legend.title = element_blank(),legend.position = 'none',
        #axis.text.x=element_blank(),
        axis.title=element_text(),
        panel.background = element_rect(fill = "transparent",color = NA),
        panel.grid.minor = element_line(color = NA), 
        panel.grid.major = element_line(color = NA),
        plot.background = element_rect(fill = "transparent",color = NA),
        plot.title = element_text(hjust = 0.5,vjust = 0.5))+
  scale_fill_manual(values=color_values)+
  scale_color_manual(values=color_values)+
  labs(x="",y="Total_skin_score" %>% str_replace_all("_"," "))

polychor_result<-polychor(annotation%>% 
                            filter(Characteristics__Group=="SSc",!is.na(characteristics__FVC_less_than_70))%>% 
                            .$characteristics__FVC_less_than_70,
                          annotation%>% 
                            filter(Characteristics__Group=="SSc",!is.na(characteristics__FVC_less_than_70)) %>% 
                            .$Characteristics__Skin_score_at_biopsy_site, std.err=TRUE, ML = TRUE)
annotation %>% 
  filter(Characteristics__Group=="SSc",!is.na(characteristics__FVC_less_than_70)) %>% 
  mutate(Characteristics__total_skin_score=replace_na(Characteristics__total_skin_score,0),
         Characteristics__Skin_score_at_biopsy_site=replace_na(Characteristics__Skin_score_at_biopsy_site,0)) %>% #.$ILD
  ggplot(aes(x=characteristics__FVC_less_than_70,y=Characteristics__Skin_score_at_biopsy_site))+
  geom_boxplot(aes(fill=characteristics__FVC_less_than_70),alpha=0.5)+
  geom_point(aes(color=characteristics__FVC_less_than_70),position =position_jitter(width=0.1))+
  theme_classic()+
  theme(legend.title = element_blank(),legend.position = 'none',
        #axis.text.x=element_blank(),
        axis.title=element_text(),
        panel.background = element_rect(fill = "transparent",color = NA),
        panel.grid.minor = element_line(color = NA), 
        panel.grid.major = element_line(color = NA),
        plot.background = element_rect(fill = "transparent",color = NA),
        plot.title = element_text(hjust = 0.5,vjust = 0.5))+
  scale_fill_manual(values=color_values)+
  scale_color_manual(values=color_values)+
  labs(x="",y="Skin_score_at_biopsy_site" %>% str_replace_all("_"," "))+
  annotate("text", label=paste("polychoric correlation:",round(polychor_result$rho,3)),3*1.1,x=1)
savefile<-paste0(output_file_base_path,"FVC-SkinScoreBiopsysite","_correlation",".png")
ggsave(file = savefile, dpi = 320, width = 140, height = 140, units = "mm",  bg = "transparent")

#correlation
FVC<-annotation%>% 
  filter(Characteristics__Group=="SSc",!is.na(characteristics__FVC_less_than_70)) %>% 
  mutate(characteristics__FVC_less_than_70=if_else(characteristics__FVC_less_than_70=="No ILD",0,1)) %>% 
  .$characteristics__FVC_less_than_70
polychor(FVC,annotation%>% 
           filter(Characteristics__Group=="SSc",!is.na(characteristics__FVC_less_than_70)) %>% .$Characteristics__total_skin_score, std.err=TRUE)
polychor(FVC,annotation%>% 
           filter(Characteristics__Group=="SSc",!is.na(characteristics__FVC_less_than_70)) %>% .$Characteristics__Skin_score_at_biopsy_site, std.err=TRUE, ML = TRUE)
polychor(annotation%>% 
           filter(Characteristics__Group=="SSc",!is.na(characteristics__FVC_less_than_70))%>% 
           .$characteristics__FVC_less_than_70,
         annotation%>% 
           filter(Characteristics__Group=="SSc",!is.na(characteristics__FVC_less_than_70)) %>% .$Characteristics__total_skin_score, std.err=TRUE)
polychor(annotation%>% 
           filter(Characteristics__Group=="SSc",!is.na(characteristics__FVC_less_than_70))%>% 
           .$characteristics__FVC_less_than_70,
         annotation%>% 
           filter(Characteristics__Group=="SSc",!is.na(characteristics__FVC_less_than_70)) %>% 
           .$Characteristics__Skin_score_at_biopsy_site, std.err=TRUE, ML = TRUE)


correlation_analysis<-function(data,x,y,color_values){
  x<-enquo(x)
  y<-enquo(y)
  polychor_result<-polychor(data[[quo_name(x)]],data[[quo_name(y)]], std.err=TRUE)
  g<-data %>% 
    mutate(x=!!x,
           y=!!y) %>% 
  ggplot(aes(x=x,y=y))+
    geom_boxplot(aes(fill=x),lwd=1)+
    geom_point(aes(color=x),position =position_jitter(width=0.1,seed=1),size=3,color="white")+
    geom_point(aes(color=x),position =position_jitter(width=0.1,seed=1),size=2)+
    theme_classic()+
    theme(legend.title = element_blank(),legend.position = 'none',
          #axis.text.x=element_blank(),
          axis.title=element_text(),
          panel.background = element_rect(fill = "transparent",color = NA),
          panel.grid.minor = element_line(color = NA), 
          panel.grid.major = element_line(color = NA),
          plot.background = element_rect(fill = "transparent",color = NA),
          plot.title = element_text(hjust = 0.5,vjust = 0.5))+
    scale_fill_manual(values=color_values)+
    scale_color_manual(values=color_values)+
    labs(x="",y=quo_name(y) %>% str_replace_all("_"," ")) + 
    annotate("text", label=paste("polychoric correlation:",round(polychor_result$rho,3)),y=max(data[[quo_name(y)]])*1.1,x=1)
  return(g)
}
annotation %>% 
  filter(Characteristics__Group=="SSc",!is.na(characteristics__FVC_less_than_70)) %>% 
correlation_analysis(characteristics__FVC_less_than_70,Characteristics__total_skin_score,color_values)
savefile<-paste0(output_file_base_path,"FVC-totalSkinScore","_correlation",".png")
ggsave(file = savefile, dpi = 320, width = 140, height = 140, units = "mm",  bg = "transparent")

annotation %>% 
  filter(Characteristics__Group=="SSc",!is.na(characteristics__FVC_less_than_70),!is.na(Characteristics__Skin_score_at_biopsy_site)) %>% 
  correlation_analysis(characteristics__FVC_less_than_70,Characteristics__Skin_score_at_biopsy_site)

savefile<-paste0(output_file_base_path,"FVC-SkinScoreBiopsysite","_correlation",".png")
ggsave(file = savefile, dpi = 320, width = 140, height = 140, units = "mm",  bg = "transparent")

annotation %>% 
  filter(Characteristics__Group=="SSc",!is.na(characteristics__FVC_less_than_70)) %>% 
  correlation_analysis(Characteristics__Diffuse_vs_limited,Characteristics__total_skin_score)
savefile<-paste0(output_file_base_path,"Diffuse-totalSkinScore","_correlation",".png")
ggsave(file = savefile, dpi = 320, width = 140, height = 140, units = "mm",  bg = "transparent")

####################################################################################################
annotation %>% 
  #mutate(Characteristics__Skin_score_at_biopsy_site=Characteristics__Skin_score_at_biopsy_site %>% replace_na(0)) %>% 
  ggplot(aes(x=Characteristics__Skin_score_at_biopsy_site))+
  geom_histogram(aes(fill=Characteristics__Skin_score_at_biopsy_site %>% as.factor()),binwidth=1)+
  theme_classic()+
  theme(legend.title = element_blank(),legend.position = 'none',
        #axis.text.x=element_blank(),
        axis.title=element_text(),
        panel.background = element_rect(fill = "transparent",color = NA),
        panel.grid.minor = element_line(color = NA), 
        panel.grid.major = element_line(color = NA),
        plot.background = element_rect(fill = "transparent",color = NA),
        plot.title = element_text(hjust = 0.5,vjust = 0.5))
savefile<-paste0(output_file_base_path,"Histgram_","SkinScoreBiopsysite",".png")
ggsave(file = savefile, dpi = 320, width = 140, height = 140, units = "mm",  bg = "transparent")

annotation %>% 
  #mutate(Characteristics__Skin_score_at_biopsy_site=Characteristics__Skin_score_at_biopsy_site %>% replace_na(0)) %>% 
  ggplot(aes(x=Characteristics__total_skin_score))+
  geom_histogram(aes(fill=Characteristics__total_skin_score %>% as.factor()),binwidth=5)+
  theme_classic()+
  theme(legend.title = element_blank(),legend.position = 'none',
        #axis.text.x=element_blank(),
        axis.title=element_text(),
        panel.background = element_rect(fill = "transparent",color = NA),
        panel.grid.minor = element_line(color = NA), 
        panel.grid.major = element_line(color = NA),
        plot.background = element_rect(fill = "transparent",color = NA),
        plot.title = element_text(hjust = 0.5,vjust = 0.5))
savefile<-paste0(output_file_base_path,"Histgram_","TotalSkinScore",".png")
ggsave(file = savefile, dpi = 320, width = 140, height = 140, units = "mm",  bg = "transparent")

color_values=c("#da6272","#0086ab","#777777","#f79646","#0086ab","#9bbb59","#bfbfbf")
annotation %>% 
  filter(!is.na(Characteristics__Age_at_biopsy_date)) %>% 
  mutate(Characteristics__Gender_1ismale=Characteristics__Gender_1ismale %>% as.factor()) %>% 
  ggplot(aes(x=Characteristics__Age_at_biopsy_date))+
  geom_histogram(aes(fill=Characteristics__Gender_1ismale),binwidth = 10,position = "stack", alpha = 1)+
  facet_wrap(.~Characteristics__Group)+
  theme_classic()+
  theme(legend.title = element_blank(),
        #legend.position = 'none',
        #axis.text.x=element_blank(),
        axis.title=element_text(),
        panel.background = element_rect(fill = "transparent",color = NA),
        panel.grid.minor = element_line(color = NA), 
        panel.grid.major = element_line(color = NA),
        plot.background = element_rect(fill = "transparent",color = NA),
        plot.title = element_text(hjust = 0.5,vjust = 0.5))+
  scale_fill_manual(values=color_values,name = "Dose", labels = c("Female", "Male"))
savefile<-paste0(output_file_base_path,"Histgram_","Age",".png")
ggsave(file = savefile, dpi = 320, width = 280, height = 140, units = "mm",  bg = "transparent")

annotation %>%
  ggplot(aes(x=Characteristics__Age_at_biopsy_date,y=Characteristics__total_skin_score))+
  geom_point(aes(color=Characteristics__Gender_1ismale %>% as.factor()))
  





gene_ in_volcanoplot.R





library(tidyverse)

setwd("~/Project/20190709_SystemicSclerosis/3.Script")

logFC_data<-read_csv("../2.Output/02_signed_2_unmerged/logFC_ModuleColor.csv")

ids=c(7422,7423)

logFC_data %>% 
  ggplot(aes(x=logFC))+geom_histogram(binwidth = 0.01,aes(fill=(Entrez_Gene_ID ==7422)))+
  geom_vline(xintercept=logFC_data[logFC_data$Entrez_Gene_ID==7422 & !is.na(logFC_data$Entrez_Gene_ID==7422),"logFC"] %>% as.numeric(),color="#00BFC4")+
  ggtitle("log FC of VEGFA (Entrez Gene ID = 7422)")

logFC_data %>% 
  ggplot(aes(logFC,-log(p.value)))+geom_point(aes(color=Entrez_Gene_ID ==7422))






make_color_cell_GO_table.R





library(tidyverse)
library(dplyr)

setwd("~/Project/20190709_SystemicSclerosis/3.Script")
signed_unsigned<-"signed"
deepSplit<-2
merged_unmerged<-""
dir_path<-paste0("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/",signed_unsigned,"_", deepSplit,"/")

GO_table<-read_csv("../2.Output/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/signed_2/GOEnrichmentTable_signed_2.csv")








make_geneid2modulecolor.R





setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
library(tidyverse)

output_file_base_path<-"../2.Output/02_signed_2_unmerged/"

lnames <- load("../2.Output/01_Compare_deepsplit_signed/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/Expdata_quontile.Rdata")
annot <- read_csv("../1.Data/Annotation/GeneAnnotation_GPL10558.csv")
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor
#?????????܂?????
deepSplit <- 2
signed_unsigned<-"signed"
merged_unmerged<-""

dir_path<-paste(output_file_base_path, signed_unsigned, "_", deepSplit,"/", sep = "")

lnames <- load(paste(dir_path,"networkConstruction_StepByStep_",signed_unsigned, "_", deepSplit, merged_unmerged, ".Rdata", sep = ""))
keep<-rowSums(Expdata_quontile$E>log2(50))>=102
Expdata<-Expdata_quontile$E[keep,]%>% t()
probes<-colnames(Expdata)
probes <- colnames(Expdata)
probes2annot <- match(probes, annot$ID)
#Entrez Gene ID?iLocusLinkID?j???ǂݍ???
allLLIDs <- annot$Entrez_Gene_ID[probes2annot]
tibble(probe=probes,
       #geneid=allLLIDs,
       module=moduleColors_signed) %>% 
  left_join(annot,by=c("probe"="ID")) %>% 
  select(probe,Entrez_Gene_ID,Symbol,module,Synonyms) %>% 
  write_csv("../2.Output/02_signed_2_unmerged/geneid2modulecolor.csv")






Networkplot_by_ModuleColor.R





library(WGCNA)
library(tidyverse)
library(ggraph)
library(tidygraph)
library(gplots)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
load("../2.Output/01_Compare_deepsplit_signed/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/Expdata_quontile.Rdata")#Expdata
#Expdata_quontile$E
keep<-rowSums(Expdata_quontile$E>log2(50))>=102
Expdata<-Expdata_quontile$E[keep,]%>% t()
probes<-colnames(Expdata)
save(probes,file = "../2.Output/03_Network_Plot/probes.Rdata")
lname<-load("../2.Output/02_signed_2_unmerged/signed_2/networkConstruction_StepByStep_signed_2.Rdata")
kwithin_df<-read_csv("../2.Output/02_signed_2_unmerged/signed_2/Probe_kWitnin_color_signed_2.csv")
#TOM <- TOMsimilarityFromExpr(Expdata, networkType = "signed", power = 20, TOMType = "signed")
#save(TOM, file = "../2.Output/02_signed_2_unmerged/TOM.RData")
# ?v?Z?ς݂̏ꍇ??TOM?f?[?^???ǂݍ???
l<-load(file = "../2.Output/02_signed_2_unmerged/TOM.RData")

gene_annotation <- read_csv("../1.Data/Annotation/GeneAnnotation_GPL10558.csv")
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor

# ME???Čv?Z
MEs <- moduleEigengenes(Expdata, moduleColors_signed)$eigengenes
row.names(MEs) <- row.names(Expdata)
module_color <- str_sub(colnames(MEs), start = 3) #?擪2????ME???폜
len <- length(module_color)

modules<-significant_module %>% head(3)
modules<-c("royalblue","lightcyan")

color_code_df<-tibble(module=significant_module,
                      color_code=col2hex(significant_module))

threshold<-0.02
modTOM<-TOM[moduleColors_signed %in% modules,moduleColors_signed %in% modules]
dimnames(modTOM)<-list(colnames(Expdata)[moduleColors_signed %in% modules],colnames(Expdata)[moduleColors_signed %in% modules])
modTOM[upper.tri(modTOM,diag = TRUE)]<-NA
t_graph<-modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>% 
  filter(TO>threshold) %>% 
  select(-TO) %>% 
  as_tbl_graph(directed=FALSE) %>% 
  activate("nodes") %>% 
  left_join(gene_annotation %>% select(ID,Symbol),by=c("name"="ID")) %>% 
  left_join(kwithin_df %>% select(ProbeID,moduleColors_signed),by=c("name"="ProbeID")) %>% 
  left_join(color_code_df,by=c("moduleColors_signed"="module"))

t_graph%>% 
  ggraph(layout = "kk")+
  geom_edge_link(aes(),alpha=0.8,colour="lightgray")+
  geom_node_point(aes(),size=10,alpha=0.8,color=t_graph %>% activate("nodes") %>% as_tibble() %>%  .$color_code)+
  geom_node_text(aes(label=Symbol),color="white",size=2)+
  theme_graph()

#?}???쐬?p
color_code_df %>% 
  mutate(x=1:nrow(color_code_df),y=1) %>% 
  ggplot(aes(x,y,color=module))+
  geom_point(size=10)+
  scale_color_manual(values=color_code_df$color_code)+
  theme_classic()



##All significant modules
target_modules<-significant_module
threshold<-0.02
modTOM<-TOM[moduleColors_signed %in% target_modules,moduleColors_signed %in% target_modules]
dimnames(modTOM)<-list(colnames(Expdata)[moduleColors_signed %in% target_modules],colnames(Expdata)[moduleColors_signed %in% target_modules])
modTOM[upper.tri(modTOM,diag = TRUE)]<-NA
t_graph<-modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>% 
  filter(TO>threshold) %>% 
  select(-TO) %>% 
  as_tbl_graph(directed=FALSE) %>% 
  activate("nodes") %>% 
  left_join(gene_annotation %>% select(ID,Symbol),by=c("name"="ID")) %>% 
  left_join(kwithin_df %>% select(ProbeID,moduleColors_signed),by=c("name"="ProbeID")) 

#unused_vars<-ls()[ls()!="t_graph"]
#for (unused_var in unused_vars){
#  rm(unused_var)
#}
#save(t_graph,file = "t_graph_all_module.Rdata")
load("t_graph_all_module.Rdata")

t_graph %>% 
  ggraph(layout = "kk")+
  geom_edge_link(aes(),alpha=0.8,colour="lightgray")+
  geom_node_point(aes(color=moduleColors_signed),size=10,alpha=0.8)+
  geom_node_text(aes(label=Symbol),color="white",size=2) 

savefile<-"../2.Output/03_Network_Plot/allModule_kk.png"
ggsave(filename = savefile,width = 1000, height = 1000, units = "mm")
savefile<-"../2.Output/03_Network_Plot/allModule_kk_hr.png"
ggsave(filename = savefile,dpi=300 ,width = 1000, height = 1000, units = "mm")






Networkplot_by_ModuleColor_edgeReduced.R





library(WGCNA)
library(tidyverse)
library(ggraph)
library(tidygraph)
library(gplots)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
load("../2.Output/01_Compare_deepsplit_signed/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/Expdata_quontile.Rdata")#Expdata
#Expdata_quontile$E
keep<-rowSums(Expdata_quontile$E>log2(50))>=102
Expdata<-Expdata_quontile$E[keep,]%>% t()
probes<-colnames(Expdata)
save(probes,file = "../2.Output/03_Network_Plot/probes.Rdata")
lname<-load("../2.Output/02_signed_2_unmerged/signed_2/networkConstruction_StepByStep_signed_2.Rdata")
kwithin_df<-read_csv("../2.Output/02_signed_2_unmerged/signed_2/Probe_kWitnin_color_signed_2.csv")
#TOM <- TOMsimilarityFromExpr(Expdata, networkType = "signed", power = 20, TOMType = "signed")
#save(TOM, file = "../2.Output/02_signed_2_unmerged/TOM.RData")
# ?v?Z?ς݂̏ꍇ??TOM?f?[?^???ǂݍ???
l<-load(file = "../2.Output/02_signed_2_unmerged/TOM.RData")

gene_annotation <- read_csv("../1.Data/Annotation/GeneAnnotation_GPL10558.csv")
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor

# ME???Čv?Z
MEs <- moduleEigengenes(Expdata, moduleColors_signed)$eigengenes
row.names(MEs) <- row.names(Expdata)
module_color <- str_sub(colnames(MEs), start = 3) #?擪2????ME???폜
len <- length(module_color)

modules<-significant_module %>% head(3)
modules<-c("royalblue","lightcyan")
modules<-c("darkgrey","tan","plum1","turquoise","royalblue","lightcyan")
color_code_df<-tibble(module=significant_module,
                      color_code=col2hex(significant_module))

threshold<-0.02
modTOM<-TOM[moduleColors_signed %in% modules,moduleColors_signed %in% modules]
dimnames(modTOM)<-list(colnames(Expdata)[moduleColors_signed %in% modules],colnames(Expdata)[moduleColors_signed %in% modules])
modTOM[upper.tri(modTOM,diag = TRUE)]<-NA
t_graph<-modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>% 
  filter(TO>threshold) %>% 
  select(-TO) %>% 
  as_tbl_graph(directed=FALSE) %>% 
  activate("nodes") %>% 
  left_join(gene_annotation %>% select(ID,Symbol),by=c("name"="ID")) %>% 
  left_join(kwithin_df %>% select(ProbeID,moduleColors_signed),by=c("name"="ProbeID")) %>% 
  left_join(color_code_df,by=c("moduleColors_signed"="module"))

label2color<-t_graph%>% 
  activate("nodes") %>% 
  as_tibble() %>% 
  rownames_to_column("label") %>% 
  mutate(label=label %>% as.integer()) %>% select(c(label,moduleColors_signed))
t_graph%>% 
  activate("edges") %>% 
  left_join(label2color,by=c("from"="label")) %>% 
  rename(color_from=moduleColors_signed) %>% 
  left_join(label2color,by=c("to"="label")) %>% 
  rename(color_to=moduleColors_signed) %>% 
  mutate(issamecolor=ifelse(color_to==color_from,color_from,NA)) %>% 
  as_tibble() %>% 
  group_by(issamecolor)%>%
  group_modify(~ head(.x, 1000L)) ->a


t_graph%>% 
  ggraph(layout = "kk")+
  geom_edge_link(aes(),alpha=0.8,colour="lightgray")+
  geom_node_point(aes(),size=10,alpha=0.8,color=t_graph %>% activate("nodes") %>% as_tibble() %>%  .$color_code)+
  geom_node_text(aes(label=Symbol),color="white",size=2)+
  theme_graph()

#?}???쐬?p
color_code_df %>% 
  mutate(x=1:nrow(color_code_df),y=1) %>% 
  ggplot(aes(x,y,color=module))+
  geom_point(size=10)+
  scale_color_manual(values=color_code_df$color_code)+
  theme_classic()



##All significant modules
target_modules<-significant_module
threshold<-0.02
modTOM<-TOM[moduleColors_signed %in% target_modules,moduleColors_signed %in% target_modules]
dimnames(modTOM)<-list(colnames(Expdata)[moduleColors_signed %in% target_modules],colnames(Expdata)[moduleColors_signed %in% target_modules])
modTOM[upper.tri(modTOM,diag = TRUE)]<-NA
t_graph<-modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>% 
  filter(TO>threshold) %>% 
  select(-TO) %>% 
  as_tbl_graph(directed=FALSE) %>% 
  activate("nodes") %>% 
  left_join(gene_annotation %>% select(ID,Symbol),by=c("name"="ID")) %>% 
  left_join(kwithin_df %>% select(ProbeID,moduleColors_signed),by=c("name"="ProbeID")) 

#unused_vars<-ls()[ls()!="t_graph"]
#for (unused_var in unused_vars){
#  rm(unused_var)
#}
#save(t_graph,file = "t_graph_all_module.Rdata")
load("t_graph_all_module.Rdata")

t_graph %>% 
  ggraph(layout = "kk")+
  geom_edge_link(aes(),alpha=0.8,colour="lightgray")+
  geom_node_point(aes(color=moduleColors_signed),size=10,alpha=0.8)+
  geom_node_text(aes(label=Symbol),color="white",size=2) 

savefile<-"../2.Output/03_Network_Plot/allModule_kk.png"
ggsave(filename = savefile,width = 1000, height = 1000, units = "mm")
savefile<-"../2.Output/03_Network_Plot/allModule_kk_hr.png"
ggsave(filename = savefile,dpi=300 ,width = 1000, height = 1000, units = "mm")






Networkplot_by_ModuleColor_edgeReduced2.R





library(WGCNA)
library(tidyverse)
library(ggraph)
library(tidygraph)
library(gplots)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
load("../2.Output/01_Compare_deepsplit_signed/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/Expdata_quontile.Rdata")#Expdata
#Expdata_quontile$E
keep<-rowSums(Expdata_quontile$E>log2(50))>=102
Expdata<-Expdata_quontile$E[keep,]%>% t()
probes<-colnames(Expdata)
save(probes,file = "../2.Output/03_Network_Plot/probes.Rdata")
lname<-load("../2.Output/02_signed_2_unmerged/signed_2/networkConstruction_StepByStep_signed_2.Rdata")
kwithin_df<-read_csv("../2.Output/02_signed_2_unmerged/signed_2/Probe_kWitnin_color_signed_2.csv")
#TOM <- TOMsimilarityFromExpr(Expdata, networkType = "signed", power = 20, TOMType = "signed")
#save(TOM, file = "../2.Output/02_signed_2_unmerged/TOM.RData")
# ?v?Z?ς݂̏ꍇ??TOM?f?[?^???ǂݍ???
l<-load(file = "../2.Output/02_signed_2_unmerged/TOM.RData")

gene_annotation <- read_csv("../1.Data/Annotation/GeneAnnotation_GPL10558.csv")
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor

# ME???Čv?Z
#MEs <- moduleEigengenes(Expdata, moduleColors_signed)$eigengenes
#row.names(MEs) <- row.names(Expdata)
#module_color <- str_sub(colnames(MEs), start = 3) #?擪2????ME???폜
#len <- length(module_color)

modules<-significant_module %>% head(3)
modules<-c("royalblue","lightcyan")
modules<-c("darkgrey","tan","plum1","turquoise","royalblue","lightcyan")
modules<-significant_module
modules<-c("darkgrey","darkmagenta","darkred","lightcyan","lightyellow","midnightblue","plum1","royalblue")
modules<-c("greenyellow","orange","tan","violet","yellow","turquoise","white")
#,"turquoise","white"


color_code_df<-tibble(module=significant_module,
                      color_code=col2hex(significant_module %>% str_replace_all("white","antiquewhite")%>% str_replace_all("lightyellow","lightyellow2")))

probe2color<-tibble(ID=colnames(Expdata),
                    moduleColor=moduleColors_signed)
threshold<-0.02
modTOM<-TOM[moduleColors_signed %in% modules,moduleColors_signed %in% modules]
dimnames(modTOM)<-list(colnames(Expdata)[moduleColors_signed %in% modules],colnames(Expdata)[moduleColors_signed %in% modules])
modTOM[upper.tri(modTOM,diag = TRUE)]<-NA

tmp<-modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>% 
  filter(TO>threshold) %>% 
  arrange(-TO) %>% 
  left_join(probe2color,by=c("from"="ID")) %>% 
  rename(color_from=moduleColor) %>% 
  left_join(probe2color,by=c("to"="ID")) %>% 
  rename(color_to=moduleColor) %>% 
  mutate(issamecolor=ifelse(color_to==color_from,color_from,NA))
not_same_color<-tmp %>% 
  filter(is.na(issamecolor)) %>% 
  mutate(issamecolor=paste(!!!rlang::syms(c("color_from", "color_to")), sep="-"))
not_same_color<-not_same_color %>% 
  mutate(issamecolor=lapply(not_same_color$issamecolor,function(x){x %>% str_split("-") %>% .[[1]] %>% sort()%>% paste(collapse = "-")}) %>% as.character())%>% 
  group_by(issamecolor)%>%
  group_modify(~ head(.x, 1000))

t_graph<-tmp%>% 
  filter(!is.na(issamecolor)) %>% 
  group_by(issamecolor)%>%
  group_modify(~ head(.x, 3000)) %>% 
  bind_rows(not_same_color) %>% 
  as_tbl_graph(directed=FALSE) %>% 
  activate("nodes") %>% 
  left_join(gene_annotation %>% select(ID,Symbol),by=c("name"="ID")) %>% 
  left_join(kwithin_df %>% select(ProbeID,moduleColors_signed),by=c("name"="ProbeID")) %>% 
  left_join(color_code_df,by=c("moduleColors_signed"="module"))

t_graph%>% 
  ggraph(layout = "kk")+
  geom_edge_link(aes(),alpha=0.8,colour="lightgray")+
  geom_node_point(aes(),size=5,alpha=0.8,color=t_graph %>% activate("nodes") %>% as_tibble() %>%  .$color_code)+
  #geom_node_text(aes(label=Symbol),color="white",size=2)+
  theme_graph()

savefile<-"../2.Output/03_Network_Plot/significant_Modules_kk_EdgeReduced_ver1.png"
#ggsave(filename = savefile,width = 330, height = 160, units = "mm")
ggsave(filename = savefile,width = 500, height = 500, units = "mm")

savefile<-"../2.Output/03_Network_Plot/All_Modules_kk_EdgeReduced_ver4_hr.png"
ggsave(filename = savefile,dpi=240, width = 600, height = 600, units = "mm")


#?}???쐬?p
color_code_df %>% 
  mutate(x=1,y=1:nrow(color_code_df)) %>% 
  ggplot(aes(x,y,color=module),alpha=0.8,scale=2)+
  geom_point(size=10)+
  scale_color_manual(values=color_code_df$color_code)+
  theme_classic()
savefile<-"../2.Output/03_Network_Plot/legend.png"
ggsave(filename = savefile,width = 60, height = 180, units = "mm")


#edge?̐??̊m?F
tmp%>% filter(!is.na(issamecolor)) %>% 
  group_by(issamecolor) %>% count %>% 
  ggplot(aes(x=reorder(issamecolor,n)))+
  geom_bar(aes(y=n),stat = "identity")+coord_flip()+coord_cartesian(ylim = c(0,12000))

#threshold?̒l?̌???
modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>%
  ggplot(aes(x=TO))+
  geom_histogram(binwidth=0.001)#+
  scale_x_continuous(breaks=seq(0,0.1,0.01))+
  coord_cartesian(xlim=c(0,0.1))
#threshold?̒l?̌???(log scale)
modTOM %>% as.data.frame() %>% 
    rownames_to_column("from") %>% 
    gather(key="to",value="TO",-from) %>%
    ggplot(aes(x=log10(TO)))+
    geom_histogram(aes(fill=TO>0.02),binwidth=0.01)+
    scale_x_continuous(breaks=seq(-10,-0.1,1))+
    theme_classic()
savefile<-"../2.Output/03_Network_Plot/threshold_adjustment.png"
ggsave(filename = savefile,width = 160, height = 160, units = "mm")







Networkplot_by_ModuleColor_PlusPPI.R





library(WGCNA)
library(tidyverse)
library(ggraph)
library(tidygraph)
library(gplots)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
load("../2.Output/01_Compare_deepsplit_signed/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/Expdata_quontile.Rdata")#Expdata
#Expdata_quontile$E
keep<-rowSums(Expdata_quontile$E>log2(50))>=102
Expdata<-Expdata_quontile$E[keep,]%>% t()
probes<-colnames(Expdata)
save(probes,file = "../2.Output/03_Network_Plot/probes.Rdata")
lname<-load("../2.Output/02_signed_2_unmerged/signed_2/networkConstruction_StepByStep_signed_2.Rdata")
kwithin_df<-read_csv("../2.Output/02_signed_2_unmerged/signed_2/Probe_kWitnin_color_signed_2.csv")
#TOM <- TOMsimilarityFromExpr(Expdata, networkType = "signed", power = 20, TOMType = "signed")
#save(TOM, file = "../2.Output/02_signed_2_unmerged/TOM.RData")
# ?v?Z?ς݂̏ꍇ??TOM?f?[?^???ǂݍ???
l<-load(file = "../2.Output/02_signed_2_unmerged/TOM.RData")

gene_annotation <- read_csv("../1.Data/Annotation/GeneAnnotation_GPL10558.csv")
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor

# ME???Čv?Z
#MEs <- moduleEigengenes(Expdata, moduleColors_signed)$eigengenes
#row.names(MEs) <- row.names(Expdata)
#module_color <- str_sub(colnames(MEs), start = 3) #?擪2????ME???폜
#len <- length(module_color)

modules<-significant_module %>% head(3)
modules<-c("royalblue","lightcyan")
modules<-c("darkgrey","tan","plum1","turquoise","royalblue","lightcyan")
modules<-significant_module
modules<-c("darkgrey","darkmagenta","darkred","lightcyan","lightyellow","midnightblue","plum1","royalblue")
modules<-c("greenyellow","orange","tan","violet","yellow","turquoise","white")
#,"turquoise","white"


color_code_df<-tibble(module=significant_module,
                      color_code=col2hex(significant_module %>% str_replace_all("white","antiquewhite")%>% str_replace_all("lightyellow","lightyellow2")))

probe2color<-tibble(ID=colnames(Expdata),
                    moduleColor=moduleColors_signed)
threshold<-0.02
modTOM<-TOM[moduleColors_signed %in% modules,moduleColors_signed %in% modules]
dimnames(modTOM)<-list(colnames(Expdata)[moduleColors_signed %in% modules],colnames(Expdata)[moduleColors_signed %in% modules])
modTOM[upper.tri(modTOM,diag = TRUE)]<-NA

tmp<-modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>% 
  filter(TO>threshold) %>% 
  arrange(-TO) %>% 
  left_join(probe2color,by=c("from"="ID")) %>% 
  rename(color_from=moduleColor) %>% 
  left_join(probe2color,by=c("to"="ID")) %>% 
  rename(color_to=moduleColor) %>% 
  mutate(issamecolor=ifelse(color_to==color_from,color_from,NA))
not_same_color<-tmp %>% 
  filter(is.na(issamecolor)) %>% 
  mutate(issamecolor=paste(!!!rlang::syms(c("color_from", "color_to")), sep="-"))
not_same_color<-not_same_color %>% 
  mutate(issamecolor=lapply(not_same_color$issamecolor,function(x){x %>% str_split("-") %>% .[[1]] %>% sort()%>% paste(collapse = "-")}) %>% as.character())%>% 
  group_by(issamecolor)%>%
  group_modify(~ head(.x, 1000))

t_graph<-tmp%>% 
  filter(!is.na(issamecolor)) %>% 
  group_by(issamecolor)%>%
  group_modify(~ head(.x, 3000)) %>% 
  bind_rows(not_same_color) %>% 
  as_tbl_graph(directed=FALSE) %>% 
  activate("nodes") %>% 
  left_join(gene_annotation %>% select(ID,Symbol),by=c("name"="ID")) %>% 
  left_join(kwithin_df %>% select(ProbeID,moduleColors_signed),by=c("name"="ProbeID")) %>% 
  left_join(color_code_df,by=c("moduleColors_signed"="module"))

tibble(Symbol=t_graph %>% activate("nodes") %>% as_tibble() %>% .$Symbol) %>% 
  write_csv("../2.Output/03_Network_Plot/PPI/small_nodes.csv")

string_mapping_path<-"../2.Output/03_Network_Plot/PPI/small_nodes_string_mapping.tsv"
string_interaction_path<-"../2.Output/03_Network_Plot/PPI/small_nodes_string_interactions.tsv"
Node2Symbol<-t_graph %>% as_tibble() %>% rownames_to_column("nodeID") %>% 
  select(c(nodeID,Symbol)) %>% 
  left_join(read_tsv(string_mapping_path),by=c("Symbol"="queryItem")) %>% 
  select(c(nodeID,preferredName))

protein_edges<-read_tsv(string_interaction_path) %>% 
  select(c(`#node1`,node2)) %>% 
  left_join(Node2Symbol,by=c(`#node1`="preferredName")) %>% 
  rename(from="nodeID") %>% 
  left_join(Node2Symbol,by=c("node2"="preferredName")) %>% 
  rename(to="nodeID") %>% 
  select(c(from,to)) %>% 
  mutate(from=as.integer(from),to=as.integer(to),alpha=1) #%>% head(100)

layout=create_layout(t_graph,layout = "kk")
#alpha_value<-t_graph2 %>% activate("edges") %>% as_tibble() %>% .$alpha %>% as.numeric()
#alpha_value<-rep(c(rep(0,6239),rep(1,745)),100)
alpha_value<-c(rep(rep(0,nrow(t_graph %>% activate("edges") %>% as.data.frame())),100),rep(rep(0.8,nrow(protein_edges)),100))
t_graph2<-t_graph %>%
  activate("edges") %>% 
  mutate(alpha=0) %>% 
  select(c(from,to,alpha)) %>% 
  bind_edges(protein_edges) 

t_graph2%>% 
  ggraph(layout = 'manual',node.positions = layout)+
  geom_edge_link(edge_alpha=alpha_value,colour="lightgray")+
  geom_node_point(aes(),size=10,alpha=0.8,color=t_graph2 %>% activate("nodes") %>% as_tibble() %>%  .$color_code)+
  geom_node_text(aes(label=Symbol),color="black",size=2)+
  theme_graph()



t_graph%>% 
  ggraph(layout = "kk")+
  geom_edge_link(aes(),alpha=0.8,colour="lightgray")+
  geom_node_point(aes(),size=10,alpha=0.8,color=t_graph %>% activate("nodes") %>% as_tibble() %>%  .$color_code)+
  geom_node_text(aes(label=Symbol),color="white",size=2)+
  theme_graph()

savefile<-"../2.Output/03_Network_Plot/LatterHalf_Modules_kk_EdgeReduced_ver7.png"
ggsave(filename = savefile,width = 330, height = 160, units = "mm")
ggsave(filename = savefile,width = 500, height = 500, units = "mm")

savefile<-"../2.Output/03_Network_Plot/All_Modules_kk_EdgeReduced_ver4_hr.png"
ggsave(filename = savefile,dpi=240, width = 600, height = 600, units = "mm")


#?}???쐬?p
color_code_df %>% 
  mutate(x=1,y=1:nrow(color_code_df)) %>% 
  ggplot(aes(x,y,color=module),alpha=0.8,scale=2)+
  geom_point(size=10)+
  scale_color_manual(values=color_code_df$color_code)+
  theme_classic()
savefile<-"../2.Output/03_Network_Plot/legend.png"
ggsave(filename = savefile,width = 60, height = 180, units = "mm")


#edge?̐??̊m?F
tmp%>% filter(!is.na(issamecolor)) %>% 
  group_by(issamecolor) %>% count %>% 
  ggplot(aes(x=reorder(issamecolor,n)))+
  geom_bar(aes(y=n),stat = "identity")+coord_flip()+coord_cartesian(ylim = c(0,12000))

#threshold?̒l?̌???
modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>%
  ggplot(aes(x=TO))+
  geom_histogram(binwidth=0.001)#+
scale_x_continuous(breaks=seq(0,0.1,0.01))+
  coord_cartesian(xlim=c(0,0.1))
#threshold?̒l?̌???(log scale)
modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>%
  ggplot(aes(x=log10(TO)))+
  geom_histogram(aes(fill=TO>0.02),binwidth=0.01)+
  scale_x_continuous(breaks=seq(-10,-0.1,1))+
  theme_classic()
savefile<-"../2.Output/03_Network_Plot/threshold_adjustment.png"
ggsave(filename = savefile,width = 160, height = 160, units = "mm")







Networkplot_by_ModuleColor_PlusPPI2.R





library(WGCNA)
library(tidyverse)
library(ggraph)
library(tidygraph)
library(gplots)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
load("../2.Output/01_Compare_deepsplit_signed/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/Expdata_quontile.Rdata")#Expdata
#Expdata_quontile$E
keep<-rowSums(Expdata_quontile$E>log2(50))>=102
Expdata<-Expdata_quontile$E[keep,]%>% t()
probes<-colnames(Expdata)
save(probes,file = "../2.Output/03_Network_Plot/probes.Rdata")
lname<-load("../2.Output/02_signed_2_unmerged/signed_2/networkConstruction_StepByStep_signed_2.Rdata")
kwithin_df<-read_csv("../2.Output/02_signed_2_unmerged/signed_2/Probe_kWitnin_color_signed_2.csv")
#TOM <- TOMsimilarityFromExpr(Expdata, networkType = "signed", power = 20, TOMType = "signed")
#save(TOM, file = "../2.Output/02_signed_2_unmerged/TOM.RData")
# ?v?Z?ς݂̏ꍇ??TOM?f?[?^???ǂݍ???
l<-load(file = "../2.Output/02_signed_2_unmerged/TOM.RData")

gene_annotation <- read_csv("../1.Data/Annotation/GeneAnnotation_GPL10558.csv")
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor

#PPI information STRING ?֘A?t?@?C???p?X
target<-"LH"
base<-"../2.Output/03_Network_Plot/PPI/"
node_path<-paste0(base,target,"_nodes.csv")
string_mapping_path<-paste0(base,target,"_nodes_string_mapping.tsv")
string_interaction_path<-paste0(base,target,"_nodes_string_interactions.tsv")


#modules<-significant_module %>% head(3)
#modules<-c("royalblue","lightcyan")
#modules<-c("darkgrey","tan","plum1","turquoise","royalblue","lightcyan")
significant_modules<-significant_module
FH_modules<-c("darkgrey","darkmagenta","darkred","lightcyan","lightyellow","midnightblue","plum1","royalblue")#FirstHalf
LH_modules<-c("greenyellow","orange","tan","violet","yellow","turquoise","white")#LatterHalf

modules<-get(paste0(target,"_modules"))
#,"turquoise","white"


color_code_df<-tibble(module=significant_module,
                      color_code=col2hex(significant_module %>% str_replace_all("white","antiquewhite")%>% str_replace_all("lightyellow","lightyellow2")))

probe2color<-tibble(ID=colnames(Expdata),
                    moduleColor=moduleColors_signed)
threshold<-0.02
modTOM<-TOM[moduleColors_signed %in% modules,moduleColors_signed %in% modules]
dimnames(modTOM)<-list(colnames(Expdata)[moduleColors_signed %in% modules],colnames(Expdata)[moduleColors_signed %in% modules])
modTOM[upper.tri(modTOM,diag = TRUE)]<-NA

tmp<-modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>% 
  filter(TO>threshold) %>% 
  arrange(-TO) %>% 
  left_join(probe2color,by=c("from"="ID")) %>% 
  rename(color_from=moduleColor) %>% 
  left_join(probe2color,by=c("to"="ID")) %>% 
  rename(color_to=moduleColor) %>% 
  mutate(issamecolor=ifelse(color_to==color_from,color_from,NA))
not_same_color<-tmp %>% 
  filter(is.na(issamecolor)) %>% 
  mutate(issamecolor=paste(!!!rlang::syms(c("color_from", "color_to")), sep="-"))
not_same_color<-not_same_color %>% 
  mutate(issamecolor=lapply(not_same_color$issamecolor,function(x){x %>% str_split("-") %>% .[[1]] %>% sort()%>% paste(collapse = "-")}) %>% as.character())%>% 
  group_by(issamecolor)%>%
  group_modify(~ head(.x, 1000))

t_graph<-tmp%>% 
  filter(!is.na(issamecolor)) %>% 
  group_by(issamecolor)%>%
  group_modify(~ head(.x, 3000)) %>% 
  bind_rows(not_same_color) %>% 
  as_tbl_graph(directed=FALSE) %>% 
  activate("nodes") %>% 
  left_join(gene_annotation %>% select(ID,Symbol),by=c("name"="ID")) %>% 
  left_join(kwithin_df %>% select(ProbeID,moduleColors_signed),by=c("name"="ProbeID")) %>% 
  left_join(color_code_df,by=c("moduleColors_signed"="module"))

tibble(Symbol=t_graph %>% activate("nodes") %>% as_tibble() %>% .$Symbol) %>% 
  write_csv(node_path)




Node2Symbol<-t_graph %>% as_tibble() %>% rownames_to_column("nodeID") %>% 
  select(c(nodeID,Symbol)) %>% 
  left_join(read_tsv(string_mapping_path),by=c("Symbol"="queryItem")) %>% 
  select(c(nodeID,preferredName))

protein_edges<-read_tsv(string_interaction_path) %>% 
  select(c(`#node1`,node2)) %>% 
  left_join(Node2Symbol,by=c(`#node1`="preferredName")) %>% 
  rename(from="nodeID") %>% 
  left_join(Node2Symbol,by=c("node2"="preferredName")) %>% 
  rename(to="nodeID") %>% 
  select(c(from,to)) %>% 
  mutate(from=as.integer(from),to=as.integer(to),alpha=1) #%>% head(100)

layout=create_layout(t_graph,layout = "kk")
#alpha_value<-t_graph2 %>% activate("edges") %>% as_tibble() %>% .$alpha %>% as.numeric()
#alpha_value<-rep(c(rep(0,6239),rep(1,745)),100)
#OnlyPPI
alpha_value<-c(rep(rep(0,nrow(t_graph %>% activate("edges") %>% as.data.frame())),100),rep(rep(0.8,nrow(protein_edges)),100))
#OnlyTOM
alpha_value<-c(rep(rep(0.8,nrow(t_graph %>% activate("edges") %>% as.data.frame())),100),rep(rep(0,nrow(protein_edges)),100))

t_graph2<-t_graph %>%
  activate("edges") %>% 
  mutate(alpha=0) %>% 
  select(c(from,to,alpha)) %>% 
  bind_edges(protein_edges) 

t_graph2%>% 
  ggraph(layout = 'manual',node.positions = layout)+
  #geom_edge_link(edge_alpha=alpha_value,colour="lightgray")+
  geom_edge_link(aes(edge_alpha=alpha),colour="lightgray")+
  geom_node_point(aes(),size=10,alpha=0.8,color=t_graph2 %>% activate("nodes") %>% as_tibble() %>%  .$color_code)+
  geom_node_text(aes(label=Symbol),color="white",size=2)+
  theme_graph()
savefile<-paste0("../2.Output/03_Network_Plot/",target,"_Modules_kk_PlusPPI_TOM_ver2.png")
ggsave(filename = savefile,width = 330, height = 160, units = "mm")
ggsave(filename = savefile,width = 500, height = 500, units = "mm")


t_graph%>% 
  ggraph(layout = 'manual',node.positions = layout)+
  #geom_edge_link(edge_alpha=alpha_value,colour="lightgray")+
  geom_edge_link(aes(),colour="lightgray")+
  geom_node_point(aes(),size=10,alpha=0.8,color=t_graph %>% activate("nodes") %>% as_tibble() %>%  .$color_code)+
  geom_node_text(aes(label=Symbol),color="white",size=2)+
  theme_graph()
savefile<-paste0("../2.Output/03_Network_Plot/",target,"_Modules_kk_PlusPPI_TOM_ver1.png")
ggsave(filename = savefile,width = 330, height = 160, units = "mm")
ggsave(filename = savefile,width = 500, height = 500, units = "mm")



#####################################################################################
#?}???쐬?p
color_code_df %>% 
  mutate(x=1,y=1:nrow(color_code_df)) %>% 
  ggplot(aes(x,y,color=module),alpha=0.8,scale=2)+
  geom_point(size=10)+
  scale_color_manual(values=color_code_df$color_code)+
  theme_classic()
savefile<-"../2.Output/03_Network_Plot/legend.png"
ggsave(filename = savefile,width = 60, height = 180, units = "mm")


#edge?̐??̊m?F
tmp%>% filter(!is.na(issamecolor)) %>% 
  group_by(issamecolor) %>% count %>% 
  ggplot(aes(x=reorder(issamecolor,n)))+
  geom_bar(aes(y=n),stat = "identity")+coord_flip()+coord_cartesian(ylim = c(0,12000))

#threshold?̒l?̌???
modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>%
  ggplot(aes(x=TO))+
  geom_histogram(binwidth=0.001)#+
scale_x_continuous(breaks=seq(0,0.1,0.01))+
  coord_cartesian(xlim=c(0,0.1))
#threshold?̒l?̌???(log scale)
modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>%
  ggplot(aes(x=log10(TO)))+
  geom_histogram(aes(fill=TO>0.02),binwidth=0.01)+
  scale_x_continuous(breaks=seq(-10,-0.1,1))+
  theme_classic()
savefile<-"../2.Output/03_Network_Plot/threshold_adjustment.png"
ggsave(filename = savefile,width = 160, height = 160, units = "mm")

ID2Symbol<-gene_annotation %>% select(ID,Symbol)

Node2Symbol_ppi<-read_tsv(string_mapping_path) %>% 
  select(c(queryItem,preferredName))

protein_edge_graph<-read_tsv(string_interaction_path) %>% 
  select(c(`#node1`,node2)) %>%
  left_join(read_tsv(string_mapping_path) %>% 
              select(c(queryItem,preferredName)),by=c(`#node1`="preferredName")) %>% 
  rename(from_symbol=queryItem)%>%
  left_join(read_tsv(string_mapping_path) %>% 
              select(c(queryItem,preferredName)),by=c("node2"="preferredName")) %>% 
  rename(to_symbol=queryItem) %>% 
  left_join(gene_annotation %>% select(ID,Symbol),by=c("from_symbol"="Symbol")) %>% 
  rename(from=ID)%>% 
  left_join(gene_annotation %>% select(ID,Symbol),by=c("to_symbol"="Symbol")) %>% 
  rename(to=ID) %>% 
  select(c(from,to,from_symbol,to_symbol)) %>% 
  distinct(from,to) %>% 
  as_tbl_graph() %>% 
  activate("nodes") %>% 
  left_join(gene_annotation %>% select(ID,Symbol),by=c("name"="ID")) %>% 
  left_join(probe2color,by=c("name"="ID")) %>% 
  left_join(color_code_df,by=c("moduleColor"="module"))

protein_edge_graph %>% ggraph(layout = 'kk')+
  geom_edge_link(aes(),colour="lightgray")+
  geom_node_point(aes(),size=10,alpha=0.8,color=protein_edge_graph %>% activate("nodes") %>% as_tibble() %>%  .$color_code)+
  geom_node_text(aes(label=Symbol),color="white",size=2)+
  theme_graph()
savefile<-paste0("../2.Output/03_Network_Plot/",target,"_Modules_kk_onlyPPI_ver1.png")
ggsave(filename = savefile,width = 330, height = 160, units = "mm")
ggsave(filename = savefile,width = 500, height = 500, units = "mm")
  





Networkplot_InterModule.R





library(tidyverse)
library(WGCNA)
library(ggraph)
library(tidygraph)
library(gplots)
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/04_Immune_Disease/004_SSc/20190716_SSc/3.Script")
load("../2.Output/01_Compare_deepsplit_signed/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/Expdata_quontile.Rdata")#Expdata
#Expdata_quontile$E
keep<-rowSums(Expdata_quontile$E>log2(50))>=102
Expdata<-Expdata_quontile$E[keep,]%>% t()
probes<-colnames(Expdata)
save(probes,file = "../2.Output/03_Network_Plot/probes.Rdata")
lname<-load("../2.Output/02_signed_2_unmerged/signed_2/networkConstruction_StepByStep_signed_2.Rdata")
kwithin_df<-read_csv("../2.Output/02_signed_2_unmerged/signed_2/Probe_kWitnin_color_signed_2.csv")
#TOM <- TOMsimilarityFromExpr(Expdata, networkType = "signed", power = 20, TOMType = "signed")
#save(TOM, file = "../2.Output/02_signed_2_unmerged/TOM.RData")
# ?v?Z?ς݂̏ꍇ??TOM?f?[?^???ǂݍ???
l<-load(file = "../2.Output/02_signed_2_unmerged/TOM.RData")

gene_annotation <- read_csv("../1.Data/Annotation/GeneAnnotation_GPL10558.csv")
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor

#PPI information STRING ?֘A?t?@?C???p?X
target<-"significant"
base<-"../2.Output/03_Network_Plot/PPI/"
node_path<-paste0(base,target,"_nodes.csv")
string_mapping_path<-paste0(base,target,"_nodes_string_mapping.tsv")
string_interaction_path<-paste0(base,target,"_nodes_string_interactions.tsv")


#modules<-significant_module %>% head(3)
#modules<-c("royalblue","lightcyan")
#modules<-c("darkgrey","tan","plum1","turquoise","royalblue","lightcyan")
significant_modules<-significant_module
FH_modules<-c("darkgrey","darkmagenta","darkred","lightcyan","lightyellow","midnightblue","plum1","royalblue")#FirstHalf
LH_modules<-c("greenyellow","orange","tan","violet","yellow","turquoise","white")#LatterHalf

modules<-get(paste0(target,"_modules"))
#,"turquoise","white"


color_code_df<-tibble(module=significant_module,
                      color_code=col2hex(significant_module %>% str_replace_all("white","antiquewhite")%>% str_replace_all("lightyellow","lightyellow2")))

probe2color<-tibble(ID=colnames(Expdata),
                    moduleColor=moduleColors_signed)
threshold<-0.02
modTOM<-TOM[moduleColors_signed %in% modules,moduleColors_signed %in% modules]
dimnames(modTOM)<-list(colnames(Expdata)[moduleColors_signed %in% modules],colnames(Expdata)[moduleColors_signed %in% modules])
modTOM[upper.tri(modTOM,diag = TRUE)]<-NA

tmp<-modTOM %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="TO",-from) %>% 
  filter(TO>threshold) %>% 
  arrange(-TO) %>% 
  left_join(probe2color,by=c("from"="ID")) %>% 
  rename(color_from=moduleColor) %>% 
  left_join(probe2color,by=c("to"="ID")) %>% 
  rename(color_to=moduleColor) %>% 
  mutate(issamecolor=ifelse(color_to==color_from,color_from,NA))
not_same_color<-tmp %>% 
  filter(is.na(issamecolor)) %>% 
  mutate(issamecolor=paste(!!!rlang::syms(c("color_from", "color_to")), sep="-"))
not_same_color<-not_same_color %>% 
  mutate(issamecolor=lapply(not_same_color$issamecolor,function(x){x %>% str_split("-") %>% .[[1]] %>% sort()%>% paste(collapse = "-")}) %>% as.character())%>% 
  rename(color_combination="issamecolor") 
not_same_color%>% 
  group_by(color_combination) %>% 
  count() %>% 
  arrange(-n) %>% 
  ggplot(aes(x=reorder(color_combination,n),fill=color_from))+
  geom_hline(yintercept = 100)+
  geom_bar(aes(y=n),stat = "identity")+coord_flip()+
  scale_y_log10()+
  scale_color_identity()+
  theme_classic()

##color_combination Network####################################################
color_combination<-not_same_color%>% 
  group_by(color_combination) %>% 
  count() %>% 
  arrange(-n) %>%
  mutate(splits=strsplit(color_combination,"-")) %>% 
  mutate(from=splits[[1]][1],
         to = splits[[1]][2]) %>% 
  mutate(logn=log1p(n)) #%>% 
  #filter(n>1000)

cc_graph<-color_combination %>% 
  select(c(from,to,n,logn)) %>% 
  as_tibble() %>% 
  #arrange(from) %>% 
  as_tbl_graph(directed=FALSE) %>% 
  activate("edges") %>% 
  filter(n>100) 
  
cc_graph%>% 
  ggraph(layout="kk")+
  geom_edge_link0(aes(edge_width=n),color="#0086ab",show.legend=FALSE)+
  geom_node_point(aes(color=name%>% str_replace_all("white","antiquewhite")),size=20)+
  geom_node_text(aes(label=name),color="#404040")+
  scale_color_identity()+
  theme_graph()+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))
savefile<-"../2.Output/03_Network_Plot/InterModule_ver4.png"
ggsave(filename = savefile,width = 190, height = 190, units = "mm",bg="transparent")

probe2color %>% 
  group_by(moduleColor) %>% 
  count()

cc_graph_FH<-color_combination %>% 
  select(c(from,to,n,logn)) %>% 
  filter(from %in% FH_modules | to %in%FH_modules)%>% 
  as_tbl_graph(directed=FALSE) %>% 
  activate("nodes") %>% 
  left_join(color_code_df,by=c("name"="module")) %>% 
  activate("nodes") %>% 
  left_join(probe2color %>% 
              group_by(moduleColor) %>% 
              count(),by=c("name"="moduleColor")) %>% 
  activate("edges") %>% 
  filter(n>100)


cc_graph_FH%>% 
  ggraph(layout="kk")+
  geom_edge_link0(aes(edge_width=n),color="#0086ab",show.legend=FALSE)+
  geom_node_point(size=20,color=cc_graph_FH %>% activate("nodes") %>% as_tibble() %>%  .$color_code) +
  geom_node_text(aes(label=name),color="#404040")+
  theme_graph()
savefile<-"../2.Output/03_Network_Plot/FH_InterModule_ver2.png"
ggsave(filename = savefile,width = 190, height = 190, units = "mm")


cc_graph_LH<-color_combination %>% 
  select(c(from,to,n,logn)) %>% 
  filter(from %in% LH_modules | to %in% LH_modules)%>% 
  as_tbl_graph(directed=FALSE) %>% 
  activate("edges") %>% 
  filter(n>2000) %>% 
  activate("nodes") %>% 
  left_join(color_code_df,by=c("name"="module")) %>% 
  activate("nodes") %>% 
  left_join(probe2color %>% 
              group_by(moduleColor) %>% 
              count(),by=c("name"="moduleColor"))


cc_graph_LH%>% 
  ggraph(layout="kk")+
  geom_edge_link(aes(edge_width=n),color="#0086ab",show.legend=FALSE)+
  geom_node_point(size=20,color=cc_graph_LH %>% activate("nodes") %>% as_tibble() %>%  .$color_code) +
  geom_node_text(aes(label=name),color="#404040")+
  scale_edge_width(c(0.1,1))+
  theme_graph()
savefile<-"../2.Output/03_Network_Plot/LH_InterModule_ver3.png"
ggsave(filename = savefile,width = 190, height = 190, units = "mm")


#########color_combination Network with Protein Protein Interaction############################################

t_graph<-tmp%>% 
  filter(!is.na(issamecolor)) %>% 
  group_by(issamecolor)%>%
  group_modify(~ head(.x, 3000)) %>% 
  bind_rows(not_same_color) %>% 
  as_tbl_graph(directed=FALSE) %>% 
  activate("nodes") %>% 
  left_join(gene_annotation %>% select(ID,Symbol),by=c("name"="ID")) %>% 
  left_join(kwithin_df %>% select(ProbeID,moduleColors_signed),by=c("name"="ProbeID")) %>% 
  left_join(color_code_df,by=c("moduleColors_signed"="module"))

tibble(Symbol=t_graph %>% activate("nodes") %>% as_tibble() %>% .$Symbol) %>% 
  write_csv(node_path)


Node2Symbol<-t_graph %>% as_tibble() %>% rownames_to_column("nodeID") %>% 
  select(c(nodeID,Symbol)) %>% 
  left_join(read_tsv(string_mapping_path),by=c("Symbol"="queryItem")) %>% 
  select(c(nodeID,preferredName))

protein_edges<-read_tsv(string_interaction_path) %>% 
  select(c(`#node1`,node2)) %>% 
  left_join(Node2Symbol,by=c(`#node1`="preferredName")) %>% 
  rename(from="nodeID") %>% 
  left_join(Node2Symbol,by=c("node2"="preferredName")) %>% 
  rename(to="nodeID") %>% 
  select(c(from,to)) %>% 
  mutate(from=as.integer(from),to=as.integer(to),alpha=1) #%>% head(100)

probe2color

nodeID2geneID<-  t_graph %>% activate("nodes") %>% 
  as_tibble() %>% 
  rownames_to_column("nodeID") %>% 
  select(c(nodeID,name)) %>% 
  mutate(nodeID=as.integer(nodeID)) %>% 
  rename(ID="name")

color_combination_PPI<-protein_edges %>% 
  select(-alpha) %>% 
  left_join(nodeID2geneID,by=c("from"="nodeID")) %>% 
  rename(from_ID="ID")%>% 
  left_join(nodeID2geneID,by=c("to"="nodeID")) %>% 
  rename(to_ID="ID")%>% 
  left_join(probe2color,by=c("from_ID"="ID")) %>% 
  rename(from_moduleColor="moduleColor")%>% 
  left_join(probe2color,by=c("to_ID"="ID")) %>% 
  rename(to_moduleColor="moduleColor") %>% 
  transmute(color_combination=str_c(from_moduleColor,"-",to_moduleColor)) 
color_combination_PPI<-color_combination_PPI%>% mutate(color_combination=lapply(color_combination_PPI$color_combination,function(x){x %>% str_split("-") %>% .[[1]] %>% sort()%>% paste(collapse = "-")}) %>% as.character()) %>% 
  group_by(color_combination) %>% 
  count() %>% 
  mutate(splits=strsplit(color_combination,"-")) %>% 
  rowwise()%>% 
  mutate(from=splits[1],to = splits[2]) %>% 
  select(c(from,to,n))
#%>% 
  #arrange(from) %>% 
  #as_tbl_graph(directed=FALSE) %>% 
  #activate("nodes") %>% 
  #left_join(color_code_df,by=c("name"="module")) %>%
  #activate("edges") %>% 
  #filter(n>100)

cc_graph_layout<-create_layout(cc_graph,layout = "kk")

cc_graph_node2corlor<-cc_graph %>% 
  activate("nodes") %>% 
  as_tibble() %>% 
  select(name) %>% 
  rownames_to_column("nodeID")
color_combination_PPI_edges<-color_combination_PPI %>% 
  left_join(cc_graph_node2corlor,by=c("from"="name")) %>% 
  rename(from_nodeID="nodeID")%>% 
  left_join(cc_graph_node2corlor,by=c("to"="name")) %>% 
  rename(to_nodeID="nodeID") %>% 
  select(from_nodeID, to_nodeID, n) %>% 
  rename(from="from_nodeID",to="to_nodeID") %>% 
  mutate(from=as.integer(from),
         to=as.integer(to),
         alpha_value=0)

color_combination_PPI_graph<-cc_graph %>% 
  activate("edges") %>% 
  mutate(alpha_value=1,
         n=0) %>% 
  bind_edges(color_combination_PPI_edges) 

color_combination_PPI_graph_layout <- ggraph::create_layout(graph = cc_graph,
                               layout = "manual", node.positions = cc_graph_layout[,1:2])
color_combination_PPI_graph%>%
  filter(n>100) %>% 
  ggraph(layout=cc_graph_layout[,1:2])+
  geom_edge_link0(aes(edge_width=n),color="#0086ab",show.legend=FALSE)+
  geom_node_point(aes(color=name),size=20) +
  geom_node_text(aes(label=name),color="#404040")+
  scale_edge_width(c(0.1,1))+
  scale_color_identity()+
  theme_graph()+
  theme(panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent"))
savefile<-"../2.Output/03_Network_Plot/InterModule/InterModule_PPI_ver2.png"
ggsave(filename = savefile,width = 190, height = 190, units = "mm",bg="transparent")

color_combination_PPI_graph %>% 
  activate("edges") %>% 
  as_tibble() %>% 
  ggplot(aes(x=reorder(paste(as.character(from), as.character(to)),n)))+
  geom_bar(aes(y=n),stat = "identity")+
  coord_flip()

######################################################################################################

#edge?̐??̊m?F
tmp%>% filter(!is.na(issamecolor)) %>% 
  group_by(issamecolor) %>% count() %>% 
  ggplot(aes(x=reorder(issamecolor,n)))+
  geom_bar(aes(y=n),stat = "identity")+coord_flip()+coord_cartesian(ylim = c(0,12000))

#?}???쐬?p
color_code_df %>% 
  mutate(x=1,y=1) %>% 
  ggplot(aes(x,y,color=module %>% str_replace_all("white","antiquewhite")),alpha=0.8,scale=2, show.legend = TRUE)+
  geom_point(size=10)+
  scale_color_identity(guide = guide_legend(ncol=2))+
  theme_void()+
  theme(legend.position="right",
        panel.background = element_blank(),#fill = "transparent"), # bg of the panel
        panel.border = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        panel.grid.major = element_blank(), # get rid of major grid
        panel.grid.minor = element_blank(), # get rid of minor grid
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_blank()#fill = "transparent")
        )+
  labs(color="")+
  xlim(-5,0)

savefile<-"../2.Output/03_Network_Plot/legend_ver3.png"
ggsave(filename = savefile,width = 60, height = 180, units = "mm",dpi=500,bg="transparent")


color_combination %>% 
  select(c(color_combination,from,to,n)) %>% 
  mutate(n=n/2) %>% 
  select(-c(to)) %>%
  rename(color=from) %>% 
  bind_rows(color_combination %>% 
              select(c(color_combination,from,to,n)) %>% 
              mutate(n=n/2) %>% 
              select(-c(from)) %>% 
              rename(color=to)) %>% 
  ggplot(aes(x=reorder(color_combination,n),fill=color))+
  geom_hline(yintercept = 100)+
  geom_hline(yintercept = 1000)+
  geom_bar(aes(y=sqrt(2*n)),stat = "identity")+
  coord_flip()+
  scale_y_log10()+
  scale_fill_identity()+
  theme_classic()+ 
  labs(x="n")


color_combination_PPI2<-protein_edges %>% 
  select(-alpha) %>% 
  left_join(nodeID2geneID,by=c("from"="nodeID")) %>% 
  rename(from_ID="ID")%>% 
  left_join(nodeID2geneID,by=c("to"="nodeID")) %>% 
  rename(to_ID="ID")%>% 
  left_join(probe2color,by=c("from_ID"="ID")) %>% 
  rename(from_moduleColor="moduleColor")%>% 
  left_join(probe2color,by=c("to_ID"="ID")) %>% 
  rename(to_moduleColor="moduleColor") %>% 
  transmute(color_combination=str_c(from_moduleColor,"-",to_moduleColor)) 
color_combination_PPI2<-color_combination_PPI2%>% mutate(color_combination=lapply(color_combination_PPI2$color_combination,function(x){x %>% str_split("-") %>% .[[1]] %>% sort()%>% paste(collapse = "-")}) %>% as.character()) %>% 
  group_by(color_combination) %>% 
  count() %>% 
  mutate(splits=strsplit(color_combination,"-")) %>% 
  rowwise()%>% 
  mutate(from=splits[1],to = splits[2])
color_combination_PPI2 %>% 
  select(c(color_combination,from,to,n)) %>% 
  mutate(n=n/2) %>% 
  select(-c(to)) %>%
  rename(color=from) %>% 
  bind_rows(color_combination_PPI2 %>% 
              select(c(color_combination,from,to,n)) %>% 
              mutate(n=n/2) %>% 
              select(-c(from)) %>% 
              rename(color=to))%>% 
  mutate(logn=log10(n)) %>% 
  arrange(-n) %>% 
  filter(n>30) %>% 
  ggplot(aes(x=reorder(color_combination,n),fill=color))+
  geom_hline(yintercept = 100)+
  geom_hline(yintercept = 300)+
  geom_bar(aes(y=sqrt(2*n)),stat = "identity")+
  coord_flip()+
  scale_y_log10()+
  scale_fill_identity()+
  theme_classic()+ 
  labs(x="n")






networkplot_test.R





#reference: https://www.slideshare.net/kashitan/tidygraphggraph

library(tidyverse)
library(ggraph)
library(tidygraph)


setwd("~/Project/20190709_SystemicSclerosis/3.Script")

lname<-load("../2.Output/02_signed_2_unmerged/signed_2/networkConstruction_StepByStep_signed_2.Rdata")
kwithin_df<-read_csv("../2.Output/02_signed_2_unmerged/signed_2/Probe_kWitnin_color_signed_2.csv")

tom_lname<-load("../2.Output/01_Compare_deepsplit_signed/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/dissTOM_signed.Rdata")
#dissTOM_signed[upper.tri(dissTOM_signed,diag = TRUE)]<-NA
colnames(dissTOM_signed)<-kwithin_df$ProbeID
len<-12875

d<-dissTOM_signed[1:len,1:len] %>% 
  as.dist() %>%
  as.matrix() 
d[upper.tri(d,diag = TRUE)]<-NA
d<-d%>%
  as.data.frame() %>% 
  mutate(from=kwithin_df$ProbeID[1:len]) %>% 
  .[moduleColors_signed!="grey",c(moduleColors_signed!="grey",TRUE)] %>% 
  gather(key=to, value=distance, -from) %>% 
  filter(!is.na(distance) & distance <= 0.900)

#d$distance=1-d$distance
head(d)
g<-as_tbl_graph(d ,directed=FALSE)
g %>% 
  ggraph(layout = "kk")+
  geom_edge_link(aes(),alpha=0.8,colour="lightgray")+
  #scale_edge_width(range = c(0.1,1))+
  geom_node_point(aes())


d %>% ggplot(aes(distance))+
  geom_histogram(binwidth = 0.001)+
  coord_cartesian(xlim = c(0.90,1.00),ylim = c(0,100000))

kwithin_df %>% 
  ggplot(aes(x=kWithin))+
  geom_histogram(binwidth = 2)



#?K???ȃ}?g???b?N?X?̍쐬
m<-matrix(rbinom(25,1,0.5),5)
diag(m)<-0
#???O?p?s???̍폜
m[upper.tri(m,diag = TRUE)]<-NA
#?s???????̒??`
colnames(m)<-c(1:ncol(m))
rownames(m)<-c(1:nrow(m))


#?????????G?b?W?̃??X?g?̍쐬
edges<-m %>% 
  as.data.frame() %>% 
  rownames_to_column(var="from")%>% 
  gather(key=to, value=edge, -from) %>% 
  filter(!is.na(edge) & edge == 1)

edges





networkplot_test2.R





#reference: https://www.slideshare.net/kashitan/tidygraphggraph

library(tidyverse)
library(ggraph)
library(tidygraph)


setwd("~/Project/20190709_SystemicSclerosis/3.Script")

lname<-load("../2.Output/02_signed_2_unmerged/signed_2/networkConstruction_StepByStep_signed_2.Rdata")
kwithin_df<-read_csv("../2.Output/02_signed_2_unmerged/signed_2/Probe_kWitnin_color_signed_2.csv")

tom_lname<-load("../2.Output/01_Compare_deepsplit_signed/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/dissTOM_signed.Rdata")
gene_annotation <- read_csv("../1.Data/Annotation/GeneAnnotation_GPL10558.csv")
significant_module<-read_csv("../")
#dissTOM_signed[upper.tri(dissTOM_signed,diag = TRUE)]<-NA
colnames(dissTOM_signed)<-kwithin_df$ProbeID
len<-12875

d<-dissTOM_signed[1:len,1:len] %>% 
  as.dist() %>%
  as.matrix() 
d[upper.tri(d,diag = TRUE)]<-NA
#chosen_color<-c("plum1","red","darkorange","lightyellow")
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(qvalue<0.01) %>% 
  .$ModuleColor
  
d2<-d%>%
  as.data.frame() %>% 
  mutate(from=kwithin_df$ProbeID[1:len]) %>% 
  .[moduleColors_signed %in% significant_module,c(moduleColors_signed %in% significant_module,TRUE)] %>% 
  #.[moduleColors_signed!="grey",c(moduleColors_signed !="grey",TRUE)] %>%
  gather(key=to, value=distance, -from) %>% 
  arrange(distance) #%>% 

  #filter(distance<=0.9)

savefile<-"../2.Output/03_Network_Plot/test/allModule_graphopt.png"
d2 %>% head(nrow(d2)/100) %>% 
  as_tbl_graph(directed=FALSE)%>% 
  activate("nodes") %>% 
  left_join(gene_annotation %>% select(ID,Symbol),by=c("name"="ID")) %>% 
  left_join(kwithin_df %>% select(ProbeID,moduleColors_signed),by=c("name"="ProbeID")) %>% 
  ggraph(layout = "graphopt")+
  geom_edge_link(aes(),alpha=0.8,colour="lightgray")+
  #scale_edge_width(range = c(0.1,1))+
  geom_node_point(aes(color=moduleColors_signed),size=10,alpha=0.8)+
  geom_node_text(aes(label=Symbol),color="white",size=2) 
ggsave(filename = savefile, width = 1000, height = 1000, units = "mm")


d %>% ggplot(aes(distance))+
  geom_histogram(binwidth = 0.001)+
  coord_cartesian(xlim = c(0.90,1.00),ylim = c(0,100000))

kwithin_df %>% 
  ggplot(aes(x=kWithin))+
  geom_histogram(binwidth = 2)





networkplot_test3.R





library(WGCNA)
library(tidyverse)
library(ggraph)
library(tidygraph)
setwd("~/Project/20190709_SystemicSclerosis/3.Script")
load("../1.Data/Expdata.Rdata")#Expdata
#Expdata$E
keep<-rowSums(Expdata$E>log2(50))>=102
Expdata<-Expdata$E[keep,]%>% t()
lname<-load("../2.Output/02_signed_2_unmerged/signed_2/networkConstruction_StepByStep_signed_2.Rdata")
kwithin_df<-read_csv("../2.Output/02_signed_2_unmerged/signed_2/Probe_kWitnin_color_signed_2.csv")

TOM <- TOMsimilarityFromExpr(Expdata, networkType = "signed", power = 20, TOMType = "signed")
save(TOM, file = "../1.Data/TOM.RData")
# ?v?Z?ς݂̏ꍇ??TOM?f?[?^???ǂݍ???
load(file = "../1.Data/TOM.RData")

# ME???Čv?Z
MEs <- moduleEigengenes(Expdata, moduleColors_signed)$eigengenes
row.names(MEs) <- row.names(Expdata)
module_color <- str_sub(colnames(MEs), start = 3) #?擪2????ME???폜
len <- length(module_color)


#tom_lname<-load("../2.Output/01_Compare_deepsplit_signed/Signed_Unsigned_Deepsplit/Normalize_by_Quontile/PAM_False/dissTOM_signed.Rdata")
#TOM<-1-dissTOM_signed
gene_annotation <- read_csv("../1.Data/Annotation/GeneAnnotation_GPL10558.csv")
significant_module<-read_csv("../2.Output/02_signed_2_unmerged/modulecolor_pvalue.csv") %>% 
  filter(pvalue<=0.01) %>% 
  .$ModuleColor

modules<-significant_module %>% head(3)

# ???W???[???̃v???[?u???ƍ?
probes <- colnames(Expdata)
inModule <- is.finite(match(moduleColors_signed, modules))
modProbes <- probes[inModule]
modGenes <- gene_annotation$Symbol[match(modProbes, gene_annotation$ID)]

# ?Ή?????TOM???I??
modTOM <- TOM[inModule, inModule]
dimnames(modTOM) <- list(modProbes, modProbes)

# Cytoscape?p?Ƀl?b?g???[?N??edge??node???o??
cyt <- exportNetworkToCytoscape(modTOM,
                                edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse = "-"), ".txt", sep = ""),
                                nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse = "-"), ".txt", sep = ""),
                                weighted = TRUE,
                                threshold = 0.02,
                                nodeNames = modProbes,
                                altNodeNames = modGenes,
                                nodeAttr = moduleColors_signed[inModule])


modTOM[upper.tri(modTOM,diag = TRUE)]<-NA
(modTOM>0.02) %>% as.data.frame() %>% 
  rownames_to_column("from") %>% 
  gather(key="to",value="isEdge",-from) %>% 
  filter(isEdge) %>% 
  select(-isEdge) %>% 
  as_tbl_graph(directed=FALSE)
  

as_tbl_graph(directed=FALSE)%>% 
  
  activate("nodes") %>% 
  left_join(gene_annotation %>% select(ID,Symbol),by=c("name"="ID")) %>% 
  left_join(kwithin_df %>% select(ProbeID,moduleColors_signed),by=c("name"="ProbeID")) %>% 
  ggraph(layout = "graphopt")+
  geom_edge_link(aes(),alpha=0.8,colour="lightgray")+
  #scale_edge_width(range = c(0.1,1))+
  geom_node_point(aes(color=moduleColors_signed),size=10,alpha=0.8)+
  geom_node_text(aes(label=Symbol),color="white",size=2) 






SLE_celltype_enrichment_WGCNAgenelist.R





rm(list = ls(all.names = TRUE))
library(WGCNA)
library(anRichment)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
options(stringsAsFactors = FALSE)

##############???͂Ɏg?p???????`?q?R???N?V???????쐬????
#?܂??֋X???̋?Collection??????
collection <- newCollection()

#Collection?̂??ƂɂȂ??e???`?q???X?g?i?????C?w?b?_?[?????ň??`?q??symbol?@entrez?̏ꍇ?͕ϊ??s?v?Ȃ̂ŉ??L?X?N???v?g?ύX)
#?t?H???_????Genelist???ꊇ?œǂݍ??܂??邽?߂̐ݒ?
#genelist?̃t?@?C?????????̂܂?enrichment???͎??̈??`?q?O???[?v???ɔ??f???????̂ŁC?K?v?ł????΂??ꂼ?ꃊ?l?[??
setwd("C:/Users/ywt2100/Desktop/P_0.1/") #?t?@?C???i?[?p?X???w??
files <- list.files(path = "C:/Users/ywt2100/Desktop/P_0.1/", full.names = T)
files_names <- list.files("C:/Users/ywt2100/Desktop/P_0.1/", full.names = F, pattern="csv") 
files_names <-  gsub(".csv", "", files_names) 
files_names
cell_names <- str_replace(files_names, pattern="_expression_gene", replacement="")

#genelist?̗R???ƂȂ?Organism???w?肷??
organism <- "Human"

#Genelist??Geneset?ɂ??邽?߂̊֐????ǂݍ??܂???
Togeneset <- function (X, Y){
  colnames(X) <- "Gene" #???????????I??Gene?ɕύX
  Symbol = X$Gene  #?֋X???̃??l?[??
  #????entrez?`???̈??`?q???????????ꍇ?́C???L??Entrez = X$Gene?ɂ??āC??2?s?̃X?N???v?g?폜
  Entrez.0 = convert2entrez(organism = organism, symbol = Symbol) 
  Entrez = unique(Entrez.0)#???`?q??(Symbol)??entrez?`???ɕϊ?????
  print(table(is.finite(Entrez)))  #entrez?`???ɕϊ??ł??????`?q???̊m?F
  newGeneSet(
    geneEntrez = Entrez,
    geneEvidence = "IEP",
    geneSource = "",
    ID = Y,
    name = Y,
    description = Y,
    source = "",
    organism = organism,
    internalClassification = Y,
    groups = Y,
    lastModified = "")
  
}

#?t?H???_????Genelist???S??Togeneset?֐???geneset?ɕϊ????Ă??????L??Collection?ɓ?????
for(i in 1:length(files)){
  genes <- read.csv(files[i], fileEncoding="UTF-8")#?????????t?@?C?????ǂݍ???
  genes2<- genes$Gene.Name
  assign(cell_names[i], genes2)
  genes <- unique(genes)#?d???L???ꍇ????????
  genes_set <- Togeneset(genes, cell_names[i])#genelist??geneset??
  assign(paste(cell_names[i]), genes_set)
  collection <- addToCollection(collection, genes_set)#Geneset?????ꂽ??collection?ɓ?????
  #?????L??"Human_CNS_5cells_mouse_collection"?͎????̂??肽??collection?????ݒ?
  rm(genes, genes_set)
}

#doublecollecion <- newCollection()
#doublecollecion <- addToCollection(doublecollecion, Human_CNS_5cells_collection,Human_INF_36cells_collection_N)

##?K?v?ł????Ώ??L?ō쐬?????R???N?V???????ۑ??B
#setwd("C:/Users/mjd9761/Documents/Enrichment/celltype_enrichment/Genecollectiondata")
#save(collection, file="Human_CNS_5cells_collection.Rdata")
#???x?R???N?V???????????΁C?????????̓??[?h???Ă??ł??g????
#load("Human_Immune_Cell_Collection.RData")

######?I?v?V????
#???????R???N?V?????????????????????^?R???N?V???????쐬???āC?ꊇ?ŃG?????b?`?????g???͂??\?B
#CNS?̈??`?q???X?g?Ɖ??ǍזE?̈??`?q???X?g???ʁX?ō쐬???????ǁC?????ɃG?????b?`?????g???͂??????Ƃ??ȂǁB
#?ȉ??????̗?
#CNS_INF_collecion <- newCollection()
#CNS_INF_collecion <- addToCollection(doublecollecion, Human_CNS_5cells_collection,Human_INF_36cells_collection_N)


#enrichment???͂Ɏg?p????collection?Ɋ܂܂??????`?q???X?g?????o?^????
#???Lenrichment?֐??ɂ????āCTOP?????܂ł?p-value???ʂ??Z?o???邩?ɉe???B?C?ӂŕύX??
nlist_collection <- length(collection)

#?זEenrichment???͂Ńq?[?g?}?b?v???쐬?????ۂɍזE?̕??я????ݒ?
cellorder <- c("Colony.Forming.Unit.Granulocyte", "Colony.Forming.Unit.Megakaryocytic", "Colony.Forming.Unit.Monocyte", "Common.myeloid.progenitor",
               "Basophil", "Eosinophil", "Monocyte", "Neutrophil", "Neutrophilic.Metamyelocyte", "mDC", "pDC", "Granulocyte.monocyte.progenitor",
               "Mature.NKcell.CD56negaCD16negaCD3nega", "Mature.NKcell.CD56negaCD16posiCD3nega", "Mature.NKcell.CD56posiCD16posiCD3nega",
               "Pro.Bcell", "Early.Bcell", "Naive.Bcell", "Mature.Bcell", "Mature.Bcell.class.able.to.switch", "Mature.Bcell.class.switched",
               "Naive.CD4posi.Tcell", "Naive.CD8posi.Tcell", "CD4posi.Tcm", "CD4posi.Tem", "CD8posi.Tcm", "CD8posi.Tem", "CD8posi.Tem_RA", "NKT",
               "Erythroid.CD34negaCD71lowGlyAposi", "Erythroid.CD34negaCD71negaGlyAposi", "Erythroid.CD34negaCD71posiGlyAnega", 
               "Erythroid.CD34negaCD71posiGlyAposi", "Erythroid.CD34posiCD71posiGlyAnega", "Megakaryocyte.erythroid.progenitor",
               "Megakaryocyte", "HSC.CD133posiCD34dim", "HSC.CD38negaCD34posi")

#CelltypeEnrichment?֐????ǂݍ??܂???
#CelltypeEnrichment?֐????ǂݍ??܂???
CelltypeEnrichment <- function (X){
  names(X) <- c("Gene", "Group") #???????????I??Gene??Group?ɕύX
  symbol = X$Gene  #?֋X???̃??l?[??
  Group = X$Group?@#?֋X???̃??l?[??
  print(table(Group))  #???`?q???X?g?̐??̊m?F?p?B
  entrez = convert2entrez(organism = "human", symbol = symbol)?@#???`?q??(Symbol)??entrez?`???ɕϊ?????
  na.omit(entrez)
  print(table(is.finite(entrez)))  #entrez?`???ɕϊ??ł??????`?q???̊m?F
  #????????anRichment?p?b?P?[?W?̒???enrichmentAnalysis?֐????g?p
  analysisresult = enrichmentAnalysis(?@?@
    classLabels = Group, identifiers = entrez,?@#classLabels???????̂??̂??????̉??͑Ώۈ??`?q?Q?Ƃ??ĔF???B
    refCollection = collection,  #reference?Ƃ??Đݒ肷???f?[?^?Z?b?g?B
    useBackground = "given",?@?@#???̓o?b?N?O???E???h?̐ݒ??i?d?v?j?B?????̓C???v?b?g?????S???`?q?̂????ǂݍ??߂????̂??g?p?B
    threshold = 1,?@?@#?G?????b?`?????g???͌??ʂ??Ƃ??ďo?͂?????csv?t?@?C???ł́Cp?l?̂??????l?B?????͍L???Ƃ??Ă????B
    thresholdType = "Bonferroni",?@#???Lp?l?ɕt?????āCBonferroni?␳????p?l??臒l?Ƃ???
    getOverlapEntrez = FALSE,?@?@#?o??csv?t?@?C???ɂāC?I?[?o?[???b?v???????`?q????entrez?`???ŏo?͂????C
    getOverlapSymbols = TRUE,    #?o??csv?t?@?C???ɂāC?I?[?o?[???b?v???????`?q????symbol?`???ŏo?͂????C
    maxReportedOverlapGenes = 10000,?@?@#???L?I?[?o?[???b?v???????`?q???ǂ̂??炢?\???????邩???ݒ??B?S?Ăق????̂?10,000??
    removeDuplicatesInDifferentClasses =FALSE,?@#???????W???[?????ɓ??????`?q?????????ꍇ???C???̂܂܏???????
    entrySeparator = ",",?@?@#?I?[?o?[???b?v???????`?q?Q?ɂ??āC","?ŕ????ďo?͂??Ă????B?ǂ??ł??悵
    ignoreLabels = "grey", #grey???W???[???̈??`?q?Q?ɂ??Ă͖????i???͂??Ȃ??j
    combineEnrichmentTables = FALSE) 
  
  
  countsInDataSet<- analysisresult$effectiveBackgroundSize  #???̓o?b?N?O???E???h?̈??`?q?????m??
  print(table(countsInDataSet))?@?@?@?@?@?@?@?@?@?@?@?@?@?@?@?@#???????o??
  
  
  Resulttable <- analysisresult$enrichmentTable
  Resulttable<-separate(Resulttable,overlapGenes, into=as.character(c(1:1000)), sep=",")
  Resulttable2 <- subset(analysisresult$enrichmentTable, analysisresult$enrichmentTable$pValue < 0.05) #?o?͗p?ɐ??`
  options(warn=-1) #???̃R?[?h??warning???o???B???????????????R?[?h
  #Overlap???`?q?Q???C1???`?q???ƂɃZ???ɕ????ĕ\???????邽?߂̃R?[?h(1000?͔C?ӁB)???̍??Ƃ?warning???o?邪???Ə??͖????Ȃ?
  write.table(Resulttable2, file = "WGCNA_enrichment_result.csv",row.names = FALSE,sep=",") #???ʂ?csv?t?@?C???ɏo??
  list <- by(Resulttable, Resulttable$class, data.frame) #?O???[?v???ƂɃt?@?C???`???????X?g??
  sapply(1:dim(list), function(x){write.csv(list[x], file=paste0("Result_", dimnames(list)[[1]][x], ".csv"), row.names=FALSE)})?@#?O???[?v???ƂɌ??ʂ??o??
  
  #?????????C???͌??ʃf?[?^(Resulttable)?̂????O???t?ɕK?v?Ȃ??̂??Ԃ????ʂ?
  Rank_all <- sapply(1:dim(list), function(x){list(list[[x]][1:nlist_collection,c(4,6,9)])})#?O???[?v???ƂɕK?v?ȗ????????o??
  names(Rank_all) <- sapply(1:dim(list), function(x){names(Rank_all) <- names(list[x])})#???L?ŃO???[?v???????????̂ł??????????l?[??
  for (i in 1:dim(list)) Rank_all[[i]]$nCommonGenes <- paste("(",Rank_all[[i]]$nCommonGenes,")") #?O???t?̃??x???p?̖??O???`
  for (i in 1:dim(list)) Rank_all[[i]]$pValue <- -(log10(as.numeric(Rank_all[[i]]$pValue)))?@#?O???t?̃??x???p?̖??O???`
  Rank_all <- na.omit(Rank_all) #???ʂ?TOP10?ɖ????Ȃ??ꍇ??NA????
  for (i in 1:dim(list)) Rank_all[[i]] <- transform(Rank_all[[i]], "Rename"=(paste(Rank_all[[i]]$dataSetName,Rank_all[[i]]$nCommonGenes,sep="?@")))#?O???t?̃??x???p?̖??O???`
  
  #list Rename
  names(Rank_all) <- sapply(1:dim(list), function(x){names(Rank_all) <- names(list[x])})
  
  
  ##????????ggplot2???g?????}???̎w??
  for (i in 1:dim(list)) ggsave(file=paste0("Result_", names(Rank_all)[i], ".pdf"),plot = ggplot(Rank_all[[i]], aes(x=reorder(Rank_all[[i]]$Rename, Rank_all[[i]]$pValue), y=Rank_all[[i]]$pValue)) +  
                                  geom_bar(stat="identity", width=.5,fill="black") +   
                                  coord_flip() +                                     
                                  xlab("Cell-type\n(#Overlap genes)") + 
                                  ylab("-log(P-value)"))
  
  #Cell-type enrichment p-value?̃q?[?g?}?b?v???}??
  CellTyperesult <- read.csv("WGCNA_enrichment_result.csv")
  logP <- -log10(CellTyperesult$pValue)
  CellTypeP <- data.frame(Module = CellTyperesult$class, CellType = CellTyperesult$dataSetID, logP = logP)
  ghm <- ggplot(CellTypeP, aes(x = CellType, y = Module, fill = logP))
  ghm <- ghm + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  ghm <- ghm + scale_x_discrete(limits = cellorder)
  ghm <- ghm + geom_tile(aes(fill = logP))
  ghm <- ghm + xlab("Cell Type") + ylab("Module")
  ghm <- ghm + geom_text(aes(label = round(CellTypeP$logP, 0)), size = 1)
  ghm <- ghm + scale_fill_gradient(low = "white", high = "red", limits = c(0, 120))
  pdf("CellType_pvalue_heatmap.pdf", width = 7, height = 7)
  plot(ghm)
  dev.off()
  
}

###################?֐??????܂?###################

#???ۂɉ???
setwd("C:/Users/ywt2100/Desktop/human")

#WGCNAgenelist <- read.table("WGCNA_gene_module.txt", header = TRUE)
#WGCNAgenelist2 <- WGCNAgenelist
#colnames(WGCNAgenelist2)[1] <- "Symbol"
WGCNAgenelist<-read.csv("WGCNAgenelist3.csv",header=TRUE)
CelltypeEnrichment(WGCNAgenelist)






summarize_BaseSpace.R





library(tidyverse)

basespace_dir_path<-"C:/Users/ywt2100/Documents/Project/20190709_SystemicSclerosis/2.Output/02_signed_2_unmerged/BaseSpace/"
files<-list.files(path =basespace_dir_path,full.names = T)
file<-files[[1]]

go_result<-c()
cp_result<-c()
module_color<-c()

for (file in files){
bs_data<-read_csv(file,skip = 4)
bs_go<-bs_data %>% 
  filter(source=="GO" & `p-value`<=0.05) %>% 
  head(10) %>% 
  transmute(bg.dir=paste(!!!rlang::syms(c("Biogroup name","direction")),sep=": ")) %>% 
  .$bg.dir %>% 
  str_c(collapse="\n")
go_result<-append(go_result,bs_go)

bs_cp<-bs_data %>% 
  filter(source=="Broad MSigDB - Canonical Pathways" & `p-value`<=0.05) %>% 
  head(10) %>% 
  transmute(bg.dir=paste(!!!rlang::syms(c("Biogroup name","direction")),sep=": ")) %>% 
  .$bg.dir %>% 
  str_c(collapse="\n")
cp_result<-append(cp_result,bs_cp)

color<-basename(file) %>% 
  str_sub(4,-5)
module_color<-append(module_color,color)
}


bs_summary_df<-data.frame(module=module_color,
                          GO_BaseSpace=go_result,
                          Cannonical_Pathway_BaseSpace=cp_result)
bs_summary_df %>% 
  write_csv("../2.Output/02_signed_2_unmerged/BaseSpace_Summary.csv")

read_csv("../2.Output/02_signed_2_unmerged/ModuleFeature_Summary.csv") %>% 
  left_join(bs_summary_df,by = "module") %>% 
  write_csv("../2.Output/02_signed_2_unmerged/ModuleFeature_Summary+BS.csv")