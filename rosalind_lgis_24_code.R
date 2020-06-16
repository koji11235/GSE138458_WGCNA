1_file_rename.R





#?A???C?f?[?^?t?@?C????rename
#???f?[?^?̃t?@?C?????x?N?g?????쐬
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data/GSE89632_RAW")
filename <- list.files()
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
write.table(filename, file = "filename.txt", sep = "", eol = "\n", quote = FALSE, row.names = FALSE)

#filename.txt??excel?ŕҏW?C?ύX???̃t?@?C???????ׂɋL?ڂ??āCrename.txt?Ƃ??ĕۑ????Ă???
#rename?p?t?@?C?????ǂݍ???
rename <- read.table("rename.txt", sep = "\t", header = T)
#?ύX???̃t?@?C?????x?N?g????paste?Ŋg???q.idat???ǉ?
rename_vector <- paste(as.vector(rename$rename), ".idat", sep = "")

#???f?[?^?̃t?@?C?????x?N?g??filename???C?ύX???̃t?@?C?????x?N?g??rename_vector?ɒu??
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data/GSE89632_RAW")
file.rename(filename, rename_vector)






2_read_normalize.R





#Illumina?}?C?N???A???C?f?[?^?i.idat?j?̓ǂݍ??݁CBG?␳?C???K??
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data/GSE89632_RAW")
#limmma?p?b?P?[?W?̓ǂݍ???
library(limma)
library(illuminaio)

#GSE89632?̓o?^?????ɂ????ƃv???b?g?t?H?[????Illumina HumanHT-12 WG-DASL V4.0 R2 expression beadchip
#?Ή?????HumanHT-12 v4.0 Whole-Genome Manifest File?i.bgx?j???p????.imat???f?[?^?????荞?????Ƃ????????܂??????Ȃ?
#HumanHT-12 v4.0 Whole-Genome Manifest File???ǂݍ??݁C?o?^?v???[?u?????m?F
#bgx_WGDASL <- readBGX("HumanHT-12_V4_0_R2_15002873_B_WGDASL.bgx")
#p <- length(bgx_WGDASL$probes$Probe_Id)        #29378
#c <- length(bgx_WGDASL$controls$ILMN_3166789)  #412
#p+c                                            #29790
#?C???~?i?J?^???O?ɋL?ڂ̃q?gWG-DASL?A?b?Z?C?̃v???[?u????29285

#.imat???f?[?^???ǂݍ??݁C?f?[?^?????m?F
#HC_1 <- readIDAT("F0_HC_1.idat")
#length(HC_1$Quants$MeanBinData)?@              #48324

#???f?[?^?̃f?[?^???ƃv???b?g?t?H?[???ɂ????v???[?u?????????Ȃ?
#?????ɁCHumanHT-12 v4.0 Manifest File (.bgx?j???ǂݍ????ŁC?o?^?v???[?u?????m?F????
#bgx <- readBGX("HumanHT-12_V4_0_R2_15002873_B.bgx")
#p <- length(bgx$probes$Probe_Id)               #47323
#c <- length(bgx$controls$Probe_Id)             #887
#p+c                                            #48210
#?C???~?i?J?^???O?ɋL?ڂ?Human HT-12 v4.0?̃v???[?u????47231

#HumanHT-12 v4.0 Manifest File?̓o?^?v???[?u???͐??f?[?^?̃f?[?^???Ƃقړ???
#???ׂĂ݂??ƁCWG-DASL?͏??ʃT???v????FPE??RNA?????T???v???ɂ????????`?q???̓A?b?Z?C?@
#?v?͖?29000???ނ̓]?ʎY???ɂ???PCR???????Ă???HumanHT-12 v4 BeadChip?i??47000?̓]?ʎY?????J?o?[?j?Ƀn?C?u???_?C?Y?????Ă??邾??
#?????????āCHumanHT-12 v4 BeadChip??Manifest File?i.bgx?j???p???ăf?[?^???ǂݍ???

#?^?[?Q?b?g?t?@?C???̐ݒ?
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data/")
targets <- readTargets("Targets.txt")


#read.idat()???p????idat?t?@?C?????ǂݍ???
idatfiles <- dir(path = "//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data/GSE89632_RAW",
                 pattern = ".idat", full.names = TRUE)
bgxfile <- dir("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data/GSE89632_RAW",
               pattern = ".bgx", full.names = TRUE)

x <- read.idat(idatfiles, bgxfile)
dim(x) #48210?s

#log2_E???Z?o
log2_E <- log2(x$E)

#detectionPValues()???p????negative control?ɑ΂???detection p-value???Z?o
x$other$Detection <- detectionPValues(x)

#neqc()???p????negative control?v???[?u???p????background?␳?Cquantile???K???Clog2?ϊ??Ccontrol?v???[?u?̏??????s??
y <- neqc(x)

#HumanHT-12 v4.0 WG-DASL???Ή??????v???[?u????29378?????CHumanHT-12 v4 BeadChip?̃v???[?u????47323
#?܂???18000?v???[?u?̃f?[?^?͑S???Ӗ??̂Ȃ??l?K?e?B?u?f?[?^?ł????C?????????܂߂?quantile???K???͕s?K?؂????????Ȃ?
#x????WG-DASL?Ɋ܂܂??Ȃ??v???[?u???t?B???^?????O???Ă???BG?␳?C???K?????s??
#WG-DASL?̃v???[?u????GPL14951-11332???ǂݍ???
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
WG_platform <- read.table("GPL14951.txt", sep = "\t", header = T)

#HumanHT-12 v4 BeadChip?̃v???[?u???????ǂݍ???
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data/")
GE_platform <- read.table("HumanHT-12_V4_0_R2.txt", sep = "\t", header = T)

#WG-DASL?̃v???[?u??????HumanHT-12 v4 BeadChip?̃v???[?u??????merge
WG_GE_platform <- merge(WG_platform, GE_platform, by.x = 2, by.y = 2, all = T)
#?????o??
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
write.csv(WG_GE_platform, file = "WG_GE_platform.csv")

#x?̒??Ɋ܂܂???WG-DASL?̃v???[?u????x?̒??Ɋ܂܂????R???g???[???v???[?u?i???`?q?v???[?u?łȂ????́j??TRUE?Ƃ???flag_WG???ݒ?
flag_WG <- (x$genes$Array_Address_Id %in% WG_platform$Array_Address_Id) |x$genes$Status != "regular"

summary(x$genes$Array_Address_Id %in% WG_platform$Array_Address_Id) #TRUE=29358
summary(x$genes$Status != "regular") #TRUE=887
summary(flag_WG) #TRUE=30207 #38?̃R???g???[???v???[?u???d??

#x??flag_WG?Ńt?B???^?????O
xfilt_WG <- x[flag_WG,]

#detectionPValues()???p????negative control?ɑ΂???detection p-value???Z?o
xfilt_WG$other$Detection <- detectionPValues(xfilt_WG)

#neqc()???p????negative control?v???[?u???p????background?␳?Cquantile???K???Clog2?ϊ??Ccontrol?v???[?u?̏??????s??
y_WG <- neqc(xfilt_WG)

#detection p-value???p?????t?B???^?????O
#detection p-value < 0.01??Present?C0.01 <= detection p-value < 0.05??Marginal?Cdetection p-value >= 0.05??Absent?Ƃ???
P <- y_WG$other$Detection < 0.01
P_M <- y_WG$other$Detection < 0.05
A <- y_WG$other$Detection >= 0.05

#???ׂẴT???v????P?̃v???[?u?ɑ΂???flag_P
flag_P <- rowSums(P) == 63?@#TRUE=1?~63??
#???ׂẴT???v????P????M?̃v???[?u?ɑ΂???flag_P_M
flag_P_M <- rowSums(P_M) == 63?@#TRUE=1?~63??
#1?ł??T???v????A?̃v???[?u?ɑ΂???flag_A
flag_A <- rowSums(A) > 0?@#TRUE=1??1???ł?????

#?e?????Ńt?B???^?????O
y_WG_filt_P <- y_WG[flag_P,]
y_WG_filt_P_M <- y_WG[flag_P_M,]
y_WG_filt_A <- y_WG[flag_A,]

dim(y)              #47323?s
dim(y_WG)           #29320?s #WG?Ńt?B???^?????O?????v???[?u??30207 - ?R???g???[???v???[?u??887 = 29320
dim(y_WG_filt_P)    #1573?s
dim(y_WG_filt_P_M)  #15989?s
dim(y_WG_filt_A)?@?@#13331?s

#???ׂẴT???v????P?̏????ł?1573?v???[?u?????c?炸??????????????
#???ׂẴT???v????P????M?̏????ł?15989?v???[?u?ł????C???͂ɓK???ȃv???[?u??????
# -> ???ׂẴT???v????P????M?idetection p-value < 0.05?j??15989?v???[?u???̗p???ĉ??͂??邱?Ƃɂ???
#Illumina?}?C?N???A???C?̓??v???P?[?g?v???[?u???Ȃ??̂?avereps?K?v?Ȃ?

#?␳?E???K???O???̃V?O?i???l?ɂ???box plot???\??
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("box_plot.pdf", paper = "a4r", width = 9.5, height = 7, pointsize = 12)
boxplot(log2_E~col(x$E), main="log2 Raw data", las=2)
boxplot(y$E~col(y$E), main="Normalized data_all", las=2)
boxplot(y_WG$E~col(y_WG$E), main="Normalized data_WG", las=2)
boxplot(y_WG_filt_P_M$E~col(y_WG_filt_P_M$E), main="Normalized data_WG_P_M", las=2)
dev.off()

#?V?O?i???l?̕??z???v???b?g
pdf("plotdensities.pdf", paper = "a4r", width = 9.5, height = 7, pointsize = 14)
par(mfrow=c(1,2))
plotDensities(x, main="Raw data", legend=FALSE)
plotDensities(y, main="Normalized data_all", legend=FALSE)
plotDensities(y_WG, main="Normalized data_WG", legend=FALSE)
plotDensities(y_WG_filt_P_M, main="Normalized data_WG_P_M", legend=FALSE)
dev.off()

#?␳?E???K???E?t?B???^?????O???̃f?[?^?}?g???b?N?X??normalized_data_WG_P_M?Ƃ???
normalized_data_WG_P_M <- y_WG_filt_P_M$E

#normalized_data?̗????ύX
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
rename <- read.table("rename.txt", sep = "\t", header = T)
colnames(normalized_data_WG_P_M) <- as.vector(rename$rename)

#?????o??
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
write.csv(normalized_data_WG_P_M, file = "normalized_data_WG_P_M.csv")

#WG?t?B???^?????O?̕␳?E???K?????f?[?^???????o??
normalized_data_WG <- y_WG$E
colnames(normalized_data_WG) <- as.vector(rename$rename)
write.csv(normalized_data_WG, file = "normalized_data_WG.csv")

#?t?B???^?????O?Ȃ??̕␳?E???K?????f?[?^???????o???i????HumanHT-12 WG-DASL V4.0 R2?ɓ??ڂ????Ă????v???[?u?̂ݒ??o???Ă݂??j
normalized_data_all <- y$E
colnames(normalized_data_all) <- as.vector(rename$rename)
write.csv(normalized_data_all, file = "normalized_data_all.csv")





2_read_normalize_usual.R





#Illumina?}?C?N???A???C?f?[?^?i.idat?j?̓ǂݍ??݁CBG?␳?C???K??
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data/GSE89632_RAW")
#limmma?p?b?P?[?W?̓ǂݍ???
library(limma)
library(illuminaio)

#GSE89632?̓o?^?????ɂ????ƃv???b?g?t?H?[????Illumina HumanHT-12 WG-DASL V4.0 R2 expression beadchip
#?Ή?????HumanHT-12 v4.0 Whole-Genome Manifest File?i.bgx?j???p????.imat???f?[?^?????荞?????Ƃ????????܂??????Ȃ?
#HumanHT-12 v4.0 Whole-Genome Manifest File???ǂݍ??݁C?o?^?v???[?u?????m?F
#bgx_WGDASL <- readBGX("HumanHT-12_V4_0_R2_15002873_B_WGDASL.bgx")
#p <- length(bgx_WGDASL$probes$Probe_Id)        #29378
#c <- length(bgx_WGDASL$controls$ILMN_3166789)  #412
#p+c                                            #29790
#?C???~?i?J?^???O?ɋL?ڂ̃q?gWG-DASL?A?b?Z?C?̃v???[?u????29285

#.imat???f?[?^???ǂݍ??݁C?f?[?^?????m?F
#HC_1 <- readIDAT("F0_HC_1.idat")
#length(HC_1$Quants$MeanBinData)?@              #48324

#???f?[?^?̃f?[?^???ƃv???b?g?t?H?[???ɂ????v???[?u?????????Ȃ?
#?????ɁCHumanHT-12 v4.0 Manifest File (.bgx?j???ǂݍ????ŁC?o?^?v???[?u?????m?F????
#bgx <- readBGX("HumanHT-12_V4_0_R2_15002873_B.bgx")
#p <- length(bgx$probes$Probe_Id)               #47323
#c <- length(bgx$controls$Probe_Id)             #887
#p+c                                            #48210
#?C???~?i?J?^???O?ɋL?ڂ?Human HT-12 v4.0?̃v???[?u????47231

#HumanHT-12 v4.0 Manifest File?̓o?^?v???[?u???͐??f?[?^?̃f?[?^???Ƃقړ???
#???ׂĂ݂??ƁCWG-DASL?͏??ʃT???v????FPE??RNA?????T???v???ɂ????????`?q???̓A?b?Z?C?@
#?v?͖?29000???ނ̓]?ʎY???ɂ???PCR???????Ă???HumanHT-12 v4 BeadChip?i??47000?̓]?ʎY?????J?o?[?j?Ƀn?C?u???_?C?Y?????Ă??邾??
#?????????āCHumanHT-12 v4 BeadChip??Manifest File?i.bgx?j???p???ăf?[?^???ǂݍ???

#?^?[?Q?b?g?t?@?C???̐ݒ?
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data/")
targets <- readTargets("Targets.txt")


#read.idat()???p????idat?t?@?C?????ǂݍ???
idatfiles <- dir(path = "//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data/GSE89632_RAW",
                 pattern = ".idat", full.names = TRUE)
bgxfile <- dir("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data/GSE89632_RAW",
               pattern = ".bgx", full.names = TRUE)

x <- read.idat(idatfiles, bgxfile)
dim(x) #48210?s

#log2_E???Z?o
log2_E <- log2(x$E)

#detectionPValues()???p????negative control?ɑ΂???detection p-value???Z?o
x$other$Detection <- detectionPValues(x)

#backgroundCorrect()???p????background?␳
y <- backgroundCorrect(x, method = "normexp", offset = 50)
#normalizeBetweenArrays()???p????quantile???K??
y <- normalizeBetweenArrays(y, method = "quantile")

#detection p-value???p?????t?B???^?????O
#detection p-value < 0.01??Present?C0.01 <= detection p-value < 0.05??Marginal?Cdetection p-value >= 0.05??Absent?Ƃ???
P <- y$other$Detection < 0.01
P_M <- y$other$Detection < 0.05
A <- y$other$Detection >= 0.05

#???ׂẴT???v????P?̃v???[?u?ɑ΂???flag_P
flag_P <- rowSums(P) == 63?@#TRUE=1?~63??
#???ׂẴT???v????P????M?̃v???[?u?ɑ΂???flag_P_M
flag_P_M <- rowSums(P_M) == 63?@#TRUE=1?~63??
#1?ł??T???v????A?̃v???[?u?ɑ΂???flag_A
flag_A <- rowSums(A) > 0?@#TRUE=1??1???ł?????

#?e?????Ńt?B???^?????O
y_filt_P <- y[flag_P,]
y_filt_P_M <- y[flag_P_M,]
y_filt_A <- y[flag_A,]

dim(y)           #48210?s
dim(y_filt_P)    #1597?s
dim(y_filt_P_M)  #16290?s
dim(y_filt_A)?@?@#31920?s

#???ׂẴT???v????P?̏????ł?1597?v???[?u?????c?炸??????????????
#???ׂẴT???v????P????M?̏????ł?16290?v???[?u?ł????C???͂ɓK???ȃv???[?u??????
# -> ???ׂẴT???v????P????M?idetection p-value < 0.05?j??16290?v???[?u???̗p???ĉ??͂??邱?Ƃɂ???

#HumanHT-12 v4.0 WG-DASL???Ή??????v???[?u????29378?????CHumanHT-12 v4 BeadChip?̃v???[?u????47323
#?܂???18000?v???[?u?̃f?[?^?͑S???Ӗ??̂Ȃ??l?K?e?B?u?f?[?^?ł????C?????????܂߂?quantile???K???͕s?K?؂????????Ȃ?
#x????WG-DASL?Ɋ܂܂??Ȃ??v???[?u???t?B???^?????O???Ă???BG?␳?C???K?????s??
#WG-DASL?̃v???[?u????GPL14951-11332???ǂݍ???
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
WG_platform <- read.table("GPL14951.csv", sep = ",", header = T)

#x?̒??Ɋ܂܂???WG-DASL?̃v???[?u??TRUE?Ƃ???flag_WG???ݒ?
flag_WG <- x$genes$Array_Address_Id %in% WG_platform$Array_Address_Id
summary(flag_WG) #TRUE=29358

#x??flag_WG?Ńt?B???^?????O
xfilt_WG <- x[flag_WG,]

#detectionPValues()???p????negative control?ɑ΂???detection p-value???Z?o
xfilt_WG$other$Detection <- detectionPValues(xfilt_WG)
#HumanHT-12 v4.0 WG-DASL?̃v???[?u?????ł?detection p-value???Z?o?????̂ɕK?v??control???????????Ȃ????߂??Z?o?ł???

#backgroundCorrect()???p????background?␳
y_WG <- backgroundCorrect(xfilt_WG, method = "normexp", offset = 50)
#normalizeBetweenArrays()???p????quantile???K??
y_WG <- normalizeBetweenArrays(y_WG, method = "quantile")

#???v???P?[?g?v???[?u??avereps
y.ave <- avereps(y, ID = y$genes$Array_Address_Id)
dim(y)     #48210?s
dim(y.ave) #48107?s
y_filt_P_M.ave <- avereps(y_filt_P_M, ID = y_filt_P_M$genes$Array_Address_Id)
dim(y_filt_P_M)?@?@#16290?s
dim(yfilt_P_M.ave) #16281?s
y_WG.ave <- avereps(y_WG, ID = y_WG$genes$Array_Address_Id)
dim(y_WG)?@?@ #29358?s
dim(y_WG.ave) #29320?s

#?␳?E???K???O???̃V?O?i???l?ɂ???box plot???\??
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("usual_box_plot.pdf", paper = "a4r", width = 9.5, height = 7, pointsize = 12)
boxplot(log2_E~col(x$E), main="log2 Raw data", las=2)
boxplot(y.ave$E~col(y.ave$E), main="Normalized data_all", las=2)
boxplot(y_filt_P_M.ave$E~col(y_filt_P_M.ave$E), main="Normalized data_P_M", las=2)
boxplot(y_WG.ave$E~col(y_WG.ave$E), main="Normalized data_WG", las=2)
dev.off()

#?V?O?i???l?̕??z???v???b?g
pdf("usual_plotdensities.pdf", paper = "a4r", width = 9.5, height = 7, pointsize = 14)
par(mfrow=c(1,2))
plotDensities(x, main="Raw data", legend=FALSE)
plotDensities(y.ave, main="Normalized data_all", legend=FALSE)
plotDensities(y_filt_P_M.ave, main="Normalized data_P_M", legend=FALSE)
plotDensities(y_WG.ave, main="Normalized data_WG", legend=FALSE)
dev.off()

#?????ύX?p?x?N?g??
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
rename <- read.table("rename.txt", sep = "\t", header = T)

#?t?B???^?????O?Ȃ??̕␳?E???K?????f?[?^???????o??
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
normalized_data_all <- y.ave$E
colnames(normalized_data_all) <- as.vector(rename$rename)
write.csv(normalized_data_all, file = "usual_normalized_data_all.csv")

#P_M?t?B???^?????O?̕␳?E???K?????f?[?^???????o??
normalized_data_P_M <- y_filt_P_M.ave$E
colnames(normalized_data_P_M) <- as.vector(rename$rename)
write.csv(normalized_data_P_M, file = "usual_normalized_data_P_M.csv")

#WG?t?B???^?????O?̕␳?E???K?????f?[?^???????o??
normalized_data_WG <- y_WG.ave$E
colnames(normalized_data_WG) <- as.vector(rename$rename)
write.csv(normalized_data_WG, file = "usual_normalized_data_WG.csv")





3_QC_WG.R





#Illumina???K???f?[?^??QC
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
#limmma?p?b?P?[?W?̓ǂݍ???
library(limma)
library(stringr)
#???K?????f?[?^?̓ǂݍ???
#Array_Address_Id?̗????ł????悤??row.names=1?͎g???Ȃ?
normalized_data <- read.table("normalized_data_WG.csv", header = T, sep = ",")
#1???ڂ̗?????Array_Address_Id?ɕύX
names(normalized_data)[1] <- "Array_Address_Id"

#MDS?v???b?g
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("WG_MDS_plot.pdf")
plotMDS(normalized_data[,-1], labels = c(1:63), col = c(rep("green3", 24), rep("blue", 20), rep("red", 19)))
cols <- c("green3", "blue", "red")
legend_labels <- c("HC", "SS", "NASH")
legend("bottomleft", legend = legend_labels, text.col = cols)
dev.off()

#?听??????PCA
#PCA???͗p?ɁC1???ڂ?Array_Adress_Id?????O???āC?s?F?T???v?????C???FArray_Adress_Id?ɂȂ??悤?]?u
data_mat <- t(normalized_data[,-1])

#?s???ύX
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
rename <- read.table("rename.txt", sep = "\t", header = T)
rownames(data_mat) <- as.vector(rename$rename)

#?m?F
class(data_mat)
dim(data_mat)
data_mat[,1:5]

#PCA?????s
pca_res <- prcomp(data_mat, scale. = T)
summary(pca_res)
pca_res$x

#PC1??PC2?̃X?R?A?v???b?g
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("WG_PCA_score_plot.pdf")
plot(pca_res$x[,1], pca_res$x[,2],
     type = "n", xlab = "PC1", ylab = "PC2", main = "PCA score plot")
text(pca_res$x[,1], pca_res$x[,2], labels = c(1:63), col = c(rep("green3", 24), rep("blue", 20), rep("red", 19)))
cols <- c("green3", "blue", "red")
legend_labels <- c("HC", "SS", "NASH")
legend("bottomleft", legend = legend_labels, text.col = cols)
dev.off()

#normalized_data?ɂ??????T???v???Ԃ̋????s?????v?Z
library("colorhcplot")
d_sample <- dist(t(normalized_data[,-1]))
class(d_sample)
#?T???v???̃N???X?^?????O
c_ward_sample <- hclust(d_sample, "ward.D2")
c_ave_sample <- hclust(d_sample, "average")
#?T???v???̃f???h???O???????\??
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("WG_dendrogram_sample.pdf", paper = "a4r", width = 9.5, height = 7, pointsize = 14)
fac <- as.factor(c(rep("HC", 24), rep("SS", 20), rep("NASH", 19)))
fac <- factor(fac, levels = c("HC", "SS", "NASH"))
colorhcplot(hc = c_ward_sample, fac = fac, col = c("green3", "blue", "red"), hang = -1, lab.cex = 0.8, main = "Cluster Dendrogram (Ward)")
colorhcplot(hc = c_ave_sample, fac = fac, col = c("green3", "blue", "red"), hang = -1, lab.cex = 0.8, main = "Cluster Dendrogram (Average)")
dev.off()






3_QC_WG_P_M.R





#Illumina???K???f?[?^??QC
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
#limmma?p?b?P?[?W?̓ǂݍ???
library(limma)
library(stringr)
#???K?????f?[?^?̓ǂݍ???
#Array_Address_Id?̗????ł????悤??row.names=1?͎g???Ȃ?
normalized_data <- read.table("normalized_data_WG_P_M.csv", header = T, sep = ",")
#1???ڂ̗?????Array_Address_Id?ɕύX
names(normalized_data)[1] <- "Array_Address_Id"

#MDS?v???b?g
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("WG_P_M_MDS_plot.pdf")
plotMDS(normalized_data[,-1], labels = c(1:63), col = c(rep("green3", 24), rep("blue", 20), rep("red", 19)))
cols <- c("green3", "blue", "red")
legend_labels <- c("HC", "SS", "NASH")
legend("bottomleft", legend = legend_labels, text.col = cols)
dev.off()

#?听??????PCA
#PCA???͗p?ɁC1???ڂ?Array_Adress_Id?????O???āC?s?F?T???v?????C???FArray_Adress_Id?ɂȂ??悤?]?u
data_mat <- t(normalized_data[,-1])

#?s???ύX
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
rename <- read.table("rename.txt", sep = "\t", header = T)
rownames(data_mat) <- as.vector(rename$rename)

#?m?F
class(data_mat)
dim(data_mat)
data_mat[,1:5]

#PCA?????s
pca_res <- prcomp(data_mat, scale. = T)
summary(pca_res)
pca_res$x

#PC1??PC2?̃X?R?A?v???b?g
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("WG_P_M_PCA_score_plot.pdf")
plot(pca_res$x[,1], pca_res$x[,2],
     type = "n", xlab = "PC1", ylab = "PC2", main = "PCA score plot")
text(pca_res$x[,1], pca_res$x[,2], labels = c(1:63), col = c(rep("green3", 24), rep("blue", 20), rep("red", 19)))
cols <- c("green3", "blue", "red")
legend_labels <- c("HC", "SS", "NASH")
legend("bottomleft", legend = legend_labels, text.col = cols)
dev.off()

#normalized_data?ɂ??????T???v???Ԃ̋????s?????v?Z
library("colorhcplot")
d_sample <- dist(t(normalized_data[,-1]))
class(d_sample)
#?T???v???̃N???X?^?????O
c_ward_sample <- hclust(d_sample, "ward.D2")
c_ave_sample <- hclust(d_sample, "average")
#?T???v???̃f???h???O???????\??
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("WG_P_M_dendrogram_sample.pdf", paper = "a4r", width = 9.5, height = 7, pointsize = 14)
fac <- as.factor(c(rep("HC", 24), rep("SS", 20), rep("NASH", 19)))
fac <- factor(fac, levels = c("HC", "SS", "NASH"))
colorhcplot(hc = c_ward_sample, fac = fac, col = c("green3", "blue", "red"), hang = -1, lab.cex = 0.8, main = "Cluster Dendrogram (Ward)")
colorhcplot(hc = c_ave_sample, fac = fac, col = c("green3", "blue", "red"), hang = -1, lab.cex = 0.8, main = "Cluster Dendrogram (Average)")
dev.off()






3_QC_WG_usual.R





#Illumina???K???f?[?^_usual??QC
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
#limmma?p?b?P?[?W?̓ǂݍ???
library(limma)
library(stringr)
#???K?????f?[?^?̓ǂݍ???
#Array_Address_Id?̗????ł????悤??row.names=1?͎g???Ȃ?
normalized_data <- read.table("usual_normalized_data_WG.csv", header = T, sep = ",")
#1???ڂ̗?????Array_Address_Id?ɕύX
names(normalized_data)[1] <- "Array_Address_Id"

#MDS?v???b?g
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("usual_WG_MDS_plot.pdf")
plotMDS(normalized_data[,-1], labels = c(1:63), col = c(rep("green3", 24), rep("blue", 20), rep("red", 19)))
cols <- c("green3", "blue", "red")
legend_labels <- c("HC", "SS", "NASH")
legend("bottomleft", legend = legend_labels, text.col = cols)
dev.off()

#?听??????PCA
#PCA???͗p?ɁC1???ڂ?Array_Adress_Id?????O???āC?s?F?T???v?????C???FArray_Adress_Id?ɂȂ??悤?]?u
data_mat <- t(normalized_data[,-1])

#?s???ύX
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
rename <- read.table("rename.txt", sep = "\t", header = T)
rownames(data_mat) <- as.vector(rename$rename)

#?m?F
class(data_mat)
dim(data_mat)
data_mat[,1:5]

#PCA?????s
pca_res <- prcomp(data_mat, scale. = T)
summary(pca_res)
pca_res$x

#PC1??PC2?̃X?R?A?v???b?g
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("usual_WG_PCA_score_plot.pdf")
plot(pca_res$x[,1], pca_res$x[,2],
     type = "n", xlab = "PC1", ylab = "PC2", main = "PCA score plot")
text(pca_res$x[,1], pca_res$x[,2], labels = c(1:63), col = c(rep("green3", 24), rep("blue", 20), rep("red", 19)))
cols <- c("green3", "blue", "red")
legend_labels <- c("HC", "SS", "NASH")
legend("bottomleft", legend = legend_labels, text.col = cols)
dev.off()

#normalized_data?ɂ??????T???v???Ԃ̋????s?????v?Z
library("colorhcplot")
d_sample <- dist(t(normalized_data[,-1]))
class(d_sample)
#?T???v???̃N???X?^?????O
c_ward_sample <- hclust(d_sample, "ward.D2")
c_ave_sample <- hclust(d_sample, "average")
#?T???v???̃f???h???O???????\??
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("usual_WG_dendrogram_sample.pdf", paper = "a4r", width = 9.5, height = 7, pointsize = 14)
fac <- as.factor(c(rep("HC", 24), rep("SS", 20), rep("NASH", 19)))
fac <- factor(fac, levels = c("HC", "SS", "NASH"))
colorhcplot(hc = c_ward_sample, fac = fac, col = c("green3", "blue", "red"), hang = -1, lab.cex = 0.8, main = "Cluster Dendrogram (Ward)")
colorhcplot(hc = c_ave_sample, fac = fac, col = c("green3", "blue", "red"), hang = -1, lab.cex = 0.8, main = "Cluster Dendrogram (Average)")
dev.off()






4_DEG_selection_usual.R





#DEG?̒??o

#limma?p?b?P?[?W?̓ǂݍ???
library(limma)
library(qvalue)
#???K?????f?[?^?̓ǂݍ???
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
normalized_data <- read.table("usual_normalized_data_WG.csv", header = T, sep = ",", row.names = 1, stringsAsFactors = F)

#?^?[?Q?b?g?t?@?C???̓ǂݍ???
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
targets <- readTargets("Targets.txt")

#?f?U?C???s?????쐬
#factor()?ŃT???v?????Q???ŃO???[?v??
f <- factor(targets$Condition, levels = unique(targets$Condition))
#?f?U?C???s????design?Ƃ???
design <- model.matrix(~0+f)
#?f?U?C???s???̗??????Q???ɕύX
colnames(design) <- levels(f)
#???`???f?????K?p
fit <- lmFit(normalized_data, design)

#HC??SS?̌???
#HC??SS?̔??r?s?????쐬
contrasts_S_H <- makeContrasts(SS - HC, levels = design)
#???r?s???????f???f?[?^fit_S_H?ɓK?p
fit_S_H <- contrasts.fit(fit, contrasts_S_H)
#?o???x?C?Y?@?Ō???
fit2_S_H <- eBayes(fit_S_H)
#qvalue?̎Z?o
for_qval_S_H <- fit2_S_H$F.p.value
qval_S_H <- qvalue(p = for_qval_S_H)

#HC??NASH?̌???
#HC??NASH?̔??r?s?????쐬
contrasts_N_H <- makeContrasts(NASH - HC, levels = design)
#???r?s???????f???f?[?^fit_N_H?ɓK?p
fit_N_H <- contrasts.fit(fit, contrasts_N_H)
#?o???x?C?Y?@?Ō???
fit2_N_H <- eBayes(fit_N_H)
#qvalue?̎Z?o
for_qval_N_H <- fit2_N_H$F.p.value
qval_N_H <- qvalue(p = for_qval_N_H)

#SS??NASH?̌???
#SS??NASH?̔??r?s?????쐬
contrasts_N_S <- makeContrasts(NASH - SS, levels = design)
#???r?s???????f???f?[?^fit_N_S?ɓK?p
fit_N_S <- contrasts.fit(fit, contrasts_N_S)
#?o???x?C?Y?@?Ō???
fit2_N_S <- eBayes(fit_N_S)
#qvalue?̎Z?o
for_qval_N_S <- fit2_N_S$F.p.value
qval_N_S <- qvalue(p = for_qval_N_S)

#MA plot
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("usual_WG_MAplot.pdf")
plotMA(fit2_S_H, ylim=c(-2,2))
plotMA(fit2_N_H, ylim=c(-2,2))
plotMA(fit2_N_S, ylim=c(-2,2))
dev.off()

#logFC?̎Z?o
HC_mean <- apply(normalized_data[,1:24], 1, mean)
SS_mean <- apply(normalized_data[,25:44], 1, mean)
NASH_mean <- apply(normalized_data[,45:63], 1, mean)
logFC_S_H <- SS_mean-HC_mean
logFC_N_H <- NASH_mean-HC_mean
logFC_N_S <- NASH_mean-SS_mean

#Array_Address_Id, logFC, q-value??matrix???쐬
S_H <- cbind(logFC_S_H, qval_S_H$qvalues)
N_H <- cbind(logFC_N_H, qval_N_H$qvalues)
N_S <- cbind(logFC_N_S, qval_N_S$qvalues)
S_H_N_H_N_S <- cbind(S_H, N_H, N_S)
colnames(S_H_N_H_N_S) <- c("S_H_logFC", "S_H_q-value",
                           "N_H_logFC", "N_H_q-value",
                           "N_S_logFC", "N_S_q-value")
#?????o??
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
write.csv(S_H_N_H_N_S, "usual_WG_logFC_qvalue.csv")

#DEG???̊m?F
SS_HC_up <- sum(S_H_N_H_N_S[,1]>=0.5 & S_H_N_H_N_S[,2]<0.05) #0
SS_HC_down <- sum(S_H_N_H_N_S[,1]<=(-0.5) & S_H_N_H_N_S[,2]<0.05) #0
NASH_HC_up <- sum(S_H_N_H_N_S[,3]>=0.5 & S_H_N_H_N_S[,4]<0.05) #0
NASH_HC_down <- sum(S_H_N_H_N_S[,3]<=(-0.5) & S_H_N_H_N_S[,4]<0.05) #0
NASH_SS_up <- sum(S_H_N_H_N_S[,5]>=0.5 & S_H_N_H_N_S[,6]<0.05) #0
NASH_SS_down <- sum(S_H_N_H_N_S[,5]<=(-0.5) & S_H_N_H_N_S[,6]<0.05) #0

#Volcano plot
q_S_H <- -log10(S_H_N_H_N_S[,2])
q_N_H <- -log10(S_H_N_H_N_S[,4])
q_N_S <- -log10(S_H_N_H_N_S[,6])
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("usual_WG_volcanoplot.pdf")
plot(S_H_N_H_N_S[,1], q_S_H,
     xlab = "Log2 fold change", ylab = "-Log10 q-value", main = "SS/HC",
     xlim = c(-2,2), ylim = c(0,3))
abline(v=0.5, lty=2)
abline(v=-0.5, lty=2)
abline(h=-log10(0.05), lty=2)

plot(S_H_N_H_N_S[,3], q_N_H,
     xlab = "Log2 fold change", ylab = "-Log10 q-value", main = "NASH/HC",
     xlim = c(-2,2), ylim = c(0,3))
abline(v=0.5, lty=2)
abline(v=-0.5, lty=2)
abline(h=-log10(0.05), lty=2)

plot(S_H_N_H_N_S[,5], q_N_S,
     xlab = "Log2 fold change", ylab = "-Log10 q-value", main = "NASH/SS",
     xlim = c(-2,2), ylim = c(0,3))
abline(v=0.5, lty=2)
abline(v=-0.5, lty=2)
abline(h=-log10(0.05), lty=2)
dev.off()

#?????̐??K?????@?ɂ???fold change?͑傫???Ȃ??????C?????????Q?ԍ??͕ς????Ȃ??̂ŁCp?l?͑傫???܂?...





4_DEG_selection_WG.R





#DEG?̒??o

#limma?p?b?P?[?W?̓ǂݍ???
library(limma)
library(qvalue)
#???K?????f?[?^?̓ǂݍ???
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
normalized_data <- read.table("normalized_data_WG.csv", header = T, sep = ",", row.names = 1, stringsAsFactors = F)

#?^?[?Q?b?g?t?@?C???̓ǂݍ???
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
targets <- readTargets("Targets.txt")

#?f?U?C???s?????쐬
#factor()?ŃT???v?????Q???ŃO???[?v??
f <- factor(targets$Condition, levels = unique(targets$Condition))
#?f?U?C???s????design?Ƃ???
design <- model.matrix(~0+f)
#?f?U?C???s???̗??????Q???ɕύX
colnames(design) <- levels(f)
#???`???f?????K?p
fit <- lmFit(normalized_data, design)

#HC??SS?̌???
#HC??SS?̔??r?s?????쐬
contrasts_S_H <- makeContrasts(SS - HC, levels = design)
#???r?s???????f???f?[?^fit_S_H?ɓK?p
fit_S_H <- contrasts.fit(fit, contrasts_S_H)
#?o???x?C?Y?@?Ō???
fit2_S_H <- eBayes(fit_S_H)
#qvalue?̎Z?o
for_qval_S_H <- fit2_S_H$F.p.value
qval_S_H <- qvalue(p = for_qval_S_H)

#HC??NASH?̌???
#HC??NASH?̔??r?s?????쐬
contrasts_N_H <- makeContrasts(NASH - HC, levels = design)
#???r?s???????f???f?[?^fit_N_H?ɓK?p
fit_N_H <- contrasts.fit(fit, contrasts_N_H)
#?o???x?C?Y?@?Ō???
fit2_N_H <- eBayes(fit_N_H)
#qvalue?̎Z?o
for_qval_N_H <- fit2_N_H$F.p.value
qval_N_H <- qvalue(p = for_qval_N_H)

#SS??NASH?̌???
#SS??NASH?̔??r?s?????쐬
contrasts_N_S <- makeContrasts(NASH - SS, levels = design)
#???r?s???????f???f?[?^fit_N_S?ɓK?p
fit_N_S <- contrasts.fit(fit, contrasts_N_S)
#?o???x?C?Y?@?Ō???
fit2_N_S <- eBayes(fit_N_S)
#qvalue?̎Z?o
for_qval_N_S <- fit2_N_S$F.p.value
qval_N_S <- qvalue(p = for_qval_N_S)

#MA plot
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("WG_MAplot.pdf")
plotMA(fit2_S_H, ylim=c(-2,2))
plotMA(fit2_N_H, ylim=c(-2,2))
plotMA(fit2_N_S, ylim=c(-2,2))
dev.off()

#logFC?̎Z?o
HC_mean <- apply(normalized_data[,1:24], 1, mean)
SS_mean <- apply(normalized_data[,25:44], 1, mean)
NASH_mean <- apply(normalized_data[,45:63], 1, mean)
logFC_S_H <- SS_mean-HC_mean
logFC_N_H <- NASH_mean-HC_mean
logFC_N_S <- NASH_mean-SS_mean

#Array_Address_Id, logFC, q-value??matrix???쐬
S_H <- cbind(logFC_S_H, qval_S_H$qvalues)
N_H <- cbind(logFC_N_H, qval_N_H$qvalues)
N_S <- cbind(logFC_N_S, qval_N_S$qvalues)
S_H_N_H_N_S <- cbind(S_H, N_H, N_S)
colnames(S_H_N_H_N_S) <- c("S_H_logFC", "S_H_q-value",
                           "N_H_logFC", "N_H_q-value",
                           "N_S_logFC", "N_S_q-value")
#?????o??
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
write.csv(S_H_N_H_N_S, "WG_logFC_qvalue.csv")

#DEG???̊m?F
SS_HC_up <- sum(S_H_N_H_N_S[,1]>=0.5 & S_H_N_H_N_S[,2]<0.05) #0
SS_HC_down <- sum(S_H_N_H_N_S[,1]<=(-0.5) & S_H_N_H_N_S[,2]<0.05) #0
NASH_HC_up <- sum(S_H_N_H_N_S[,3]>=0.5 & S_H_N_H_N_S[,4]<0.05) #0
NASH_HC_down <- sum(S_H_N_H_N_S[,3]<=(-0.5) & S_H_N_H_N_S[,4]<0.05) #0
NASH_SS_up <- sum(S_H_N_H_N_S[,5]>=0.5 & S_H_N_H_N_S[,6]<0.05) #0
NASH_SS_down <- sum(S_H_N_H_N_S[,5]<=(-0.5) & S_H_N_H_N_S[,6]<0.05) #0

#Volcano plot
q_S_H <- -log10(S_H_N_H_N_S[,2])
q_N_H <- -log10(S_H_N_H_N_S[,4])
q_N_S <- -log10(S_H_N_H_N_S[,6])
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("WG_volcanoplot.pdf")
plot(S_H_N_H_N_S[,1], q_S_H,
     xlab = "Log2 fold change", ylab = "-Log10 q-value", main = "SS/HC",
     xlim = c(-2,2), ylim = c(0,3))
abline(v=0.5, lty=2)
abline(v=-0.5, lty=2)
abline(h=-log10(0.05), lty=2)

plot(S_H_N_H_N_S[,3], q_N_H,
     xlab = "Log2 fold change", ylab = "-Log10 q-value", main = "NASH/HC",
     xlim = c(-2,2), ylim = c(0,3))
abline(v=0.5, lty=2)
abline(v=-0.5, lty=2)
abline(h=-log10(0.05), lty=2)

plot(S_H_N_H_N_S[,5], q_N_S,
     xlab = "Log2 fold change", ylab = "-Log10 q-value", main = "NASH/SS",
     xlim = c(-2,2), ylim = c(0,3))
abline(v=0.5, lty=2)
abline(v=-0.5, lty=2)
abline(h=-log10(0.05), lty=2)
dev.off()





4_DEG_selection_WG_P_M.R





#DEG?̒??o

#limma?p?b?P?[?W?̓ǂݍ???
library(limma)
library(qvalue)
#???K?????f?[?^?̓ǂݍ???
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
normalized_data <- read.table("normalized_data_WG_P_M.csv", header = T, sep = ",", row.names = 1, stringsAsFactors = F)

#?^?[?Q?b?g?t?@?C???̓ǂݍ???
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
targets <- readTargets("Targets.txt")

#?f?U?C???s?????쐬
#factor()?ŃT???v?????Q???ŃO???[?v??
f <- factor(targets$Condition, levels = unique(targets$Condition))
#?f?U?C???s????design?Ƃ???
design <- model.matrix(~0+f)
#?f?U?C???s???̗??????Q???ɕύX
colnames(design) <- levels(f)
#???`???f?????K?p
fit <- lmFit(normalized_data, design)

#HC??SS?̌???
#HC??SS?̔??r?s?????쐬
contrasts_S_H <- makeContrasts(SS - HC, levels = design)
#???r?s???????f???f?[?^fit_S_H?ɓK?p
fit_S_H <- contrasts.fit(fit, contrasts_S_H)
#?o???x?C?Y?@?Ō???
fit2_S_H <- eBayes(fit_S_H)
#qvalue?̎Z?o
for_qval_S_H <- fit2_S_H$F.p.value
qval_S_H <- qvalue(p = for_qval_S_H)

#HC??NASH?̌???
#HC??NASH?̔??r?s?????쐬
contrasts_N_H <- makeContrasts(NASH - HC, levels = design)
#???r?s???????f???f?[?^fit_N_H?ɓK?p
fit_N_H <- contrasts.fit(fit, contrasts_N_H)
#?o???x?C?Y?@?Ō???
fit2_N_H <- eBayes(fit_N_H)
#qvalue?̎Z?o
for_qval_N_H <- fit2_N_H$F.p.value
qval_N_H <- qvalue(p = for_qval_N_H)

#SS??NASH?̌???
#SS??NASH?̔??r?s?????쐬
contrasts_N_S <- makeContrasts(NASH - SS, levels = design)
#???r?s???????f???f?[?^fit_N_S?ɓK?p
fit_N_S <- contrasts.fit(fit, contrasts_N_S)
#?o???x?C?Y?@?Ō???
fit2_N_S <- eBayes(fit_N_S)
#qvalue?̎Z?o
for_qval_N_S <- fit2_N_S$F.p.value
qval_N_S <- qvalue(p = for_qval_N_S)

#MA plot
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("WG_P_M_MAplot.pdf")
plotMA(fit2_S_H, ylim=c(-2,2))
plotMA(fit2_N_H, ylim=c(-2,2))
plotMA(fit2_N_S, ylim=c(-2,2))
dev.off()

#logFC?̎Z?o
HC_mean <- apply(normalized_data[,1:24], 1, mean)
SS_mean <- apply(normalized_data[,25:44], 1, mean)
NASH_mean <- apply(normalized_data[,45:63], 1, mean)
logFC_S_H <- SS_mean-HC_mean
logFC_N_H <- NASH_mean-HC_mean
logFC_N_S <- NASH_mean-SS_mean

#Array_Address_Id, logFC, q-value??matrix???쐬
S_H <- cbind(logFC_S_H, qval_S_H$qvalues)
N_H <- cbind(logFC_N_H, qval_N_H$qvalues)
N_S <- cbind(logFC_N_S, qval_N_S$qvalues)
S_H_N_H_N_S <- cbind(S_H, N_H, N_S)
colnames(S_H_N_H_N_S) <- c("S_H_logFC", "S_H_q-value",
                           "N_H_logFC", "N_H_q-value",
                           "N_S_logFC", "N_S_q-value")
#?????o??
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
write.csv(S_H_N_H_N_S, "WG_P_M_logFC_qvalue.csv")

#DEG???̊m?F
SS_HC_up <- sum(S_H_N_H_N_S[,1]>=0.5 & S_H_N_H_N_S[,2]<0.05) #0
SS_HC_down <- sum(S_H_N_H_N_S[,1]<=(-0.5) & S_H_N_H_N_S[,2]<0.05) #0
NASH_HC_up <- sum(S_H_N_H_N_S[,3]>=0.5 & S_H_N_H_N_S[,4]<0.05) #0
NASH_HC_down <- sum(S_H_N_H_N_S[,3]<=(-0.5) & S_H_N_H_N_S[,4]<0.05) #0
NASH_SS_up <- sum(S_H_N_H_N_S[,5]>=0.5 & S_H_N_H_N_S[,6]<0.05) #0
NASH_SS_down <- sum(S_H_N_H_N_S[,5]<=(-0.5) & S_H_N_H_N_S[,6]<0.05) #0

#Volcano plot
q_S_H <- -log10(S_H_N_H_N_S[,2])
q_N_H <- -log10(S_H_N_H_N_S[,4])
q_N_S <- -log10(S_H_N_H_N_S[,6])
setwd("//techdoc.jti.co.jp/DavWWWRoot/org/bpr/Shared Documents/????14G/Project/Metabolic diseases/NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("WG_P_M_volcanoplot.pdf")
plot(S_H_N_H_N_S[,1], q_S_H,
     xlab = "Log2 fold change", ylab = "-Log10 q-value", main = "SS/HC",
     xlim = c(-2,2), ylim = c(0,3))
abline(v=0.5, lty=2)
abline(v=-0.5, lty=2)
abline(h=-log10(0.05), lty=2)

plot(S_H_N_H_N_S[,3], q_N_H,
     xlab = "Log2 fold change", ylab = "-Log10 q-value", main = "NASH/HC",
     xlim = c(-2,2), ylim = c(0,3))
abline(v=0.5, lty=2)
abline(v=-0.5, lty=2)
abline(h=-log10(0.05), lty=2)

plot(S_H_N_H_N_S[,5], q_N_S,
     xlab = "Log2 fold change", ylab = "-Log10 q-value", main = "NASH/SS",
     xlim = c(-2,2), ylim = c(0,3))
abline(v=0.5, lty=2)
abline(v=-0.5, lty=2)
abline(h=-log10(0.05), lty=2)
dev.off()





5_series_matrix_QC.R





#GSE89632 series matrix??QC
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
#limmma?p?b?P?[?W?̓ǂݍ???
library(limma)
library(stringr)
#series matrix?̓ǂݍ???
#Probe_ID?̗????ł????悤??row.names=1?͎g???Ȃ?
normalized_data <- read.table("GSE89632_series_matrix.txt", row.names = 1, header = T, sep = "\t")

#box plot
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("series_matrix_box_plot.pdf", paper = "a4", width = 9.5, height = 7, pointsize = 12)
boxplot(normalized_data, main="Series matrix", las=2)
dev.off()

#?V?O?i???l?̕??z???v???b?g
pdf("series_matrix_plotdensities.pdf", width = 9.5, height = 7, pointsize = 14)
plotDensities(normalized_data, main="Series matrix", legend=FALSE)
dev.off()

#MDS?v???b?g
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("series_matrix_MDS_plot.pdf")
plotMDS(normalized_data, labels = c(1:63), col = c(rep("green3", 24), rep("blue", 20), rep("red", 19)))
cols <- c("green3", "blue", "red")
legend_labels <- c("HC", "SS", "NASH")
legend("bottomleft", legend = legend_labels, text.col = cols)
dev.off()

#?听??????PCA
#PCA???͗p?ɁC?s?F?T???v?????C???FArray_Adress_Id?ɂȂ??悤?]?u
data_mat <- t(normalized_data)

#?s???ύX
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
rename <- read.table("rename.txt", sep = "\t", header = T)
rownames(data_mat) <- as.vector(rename$rename)

#?m?F
class(data_mat)
dim(data_mat)
data_mat[,1:5]

#PCA?????s
pca_res <- prcomp(data_mat, scale. = T)
summary(pca_res)
pca_res$x

#PC1??PC2?̃X?R?A?v???b?g
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("series_matrix_PCA_score_plot.pdf")
plot(pca_res$x[,1], pca_res$x[,2],
     type = "n", xlab = "PC1 (23.23%)", ylab = "PC2(12.97%)", main = "PCA score plot")
text(pca_res$x[,1], pca_res$x[,2], labels = c(1:63), col = c(rep("green3", 24), rep("blue", 20), rep("red", 19)))
cols <- c("green3", "blue", "red")
legend_labels <- c("HC", "SS", "NASH")
legend("bottomleft", legend = legend_labels, text.col = cols)
dev.off()

#normalized_data?ɂ??????T???v???Ԃ̋????s?????v?Z
library("colorhcplot")
d_sample <- dist(t(normalized_data))
class(d_sample)
#?T???v???̃N???X?^?????O
c_ward_sample <- hclust(d_sample, "ward.D2")
c_ave_sample <- hclust(d_sample, "average")
#?T???v???̃f???h???O???????\??
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("series_matrix_dendrogram_sample.pdf", paper = "a4r", width = 9.5, height = 7, pointsize = 14)
fac <- as.factor(c(rep("HC", 24), rep("SS", 20), rep("NASH", 19)))
fac <- factor(fac, levels = c("HC", "SS", "NASH"))
colorhcplot(hc = c_ward_sample, fac = fac, col = c("green3", "blue", "red"), hang = -1, lab.cex = 0.8, main = "Cluster Dendrogram (Ward)")
colorhcplot(hc = c_ave_sample, fac = fac, col = c("green3", "blue", "red"), hang = -1, lab.cex = 0.8, main = "Cluster Dendrogram (Average)")
dev.off()






6_series_matrix_DEG_selection.R





#series matrix?ɂ?????DEG?̒??o

#limma?p?b?P?[?W?̓ǂݍ???
library(limma)
library(qvalue)
#???K?????f?[?^?̓ǂݍ???
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
normalized_data <- read.table("GSE89632_series_matrix.txt", header = T, sep = "\t", row.names = 1, stringsAsFactors = F)

#?^?[?Q?b?g?t?@?C???̓ǂݍ???
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
targets <- readTargets("Targets.txt")

#?f?U?C???s?????쐬
#factor()?ŃT???v?????Q???ŃO???[?v??
f <- factor(targets$Condition, levels = unique(targets$Condition))
#?f?U?C???s????design?Ƃ???
design <- model.matrix(~0+f)
#?f?U?C???s???̗??????Q???ɕύX
colnames(design) <- levels(f)
#???`???f?????K?p
fit <- lmFit(normalized_data, design)

#HC??SS?̌???
#HC??SS?̔??r?s?????쐬
contrasts_S_H <- makeContrasts(SS - HC, levels = design)
#???r?s???????f???f?[?^fit_S_H?ɓK?p
fit_S_H <- contrasts.fit(fit, contrasts_S_H)
#?o???x?C?Y?@?Ō???
fit2_S_H <- eBayes(fit_S_H)
#qvalue?̎Z?o
for_qval_S_H <- fit2_S_H$F.p.value
qval_S_H <- qvalue(p = for_qval_S_H)

#HC??NASH?̌???
#HC??NASH?̔??r?s?????쐬
contrasts_N_H <- makeContrasts(NASH - HC, levels = design)
#???r?s???????f???f?[?^fit_N_H?ɓK?p
fit_N_H <- contrasts.fit(fit, contrasts_N_H)
#?o???x?C?Y?@?Ō???
fit2_N_H <- eBayes(fit_N_H)
#qvalue?̎Z?o
for_qval_N_H <- fit2_N_H$F.p.value
qval_N_H <- qvalue(p = for_qval_N_H)

#SS??NASH?̌???
#SS??NASH?̔??r?s?????쐬
contrasts_N_S <- makeContrasts(NASH - SS, levels = design)
#???r?s???????f???f?[?^fit_N_S?ɓK?p
fit_N_S <- contrasts.fit(fit, contrasts_N_S)
#?o???x?C?Y?@?Ō???
fit2_N_S <- eBayes(fit_N_S)
#qvalue?̎Z?o
for_qval_N_S <- fit2_N_S$F.p.value
qval_N_S <- qvalue(p = for_qval_N_S)

#MA plot
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("series_matrix_MAplot.pdf")
plotMA(fit2_S_H, ylim=c(-5,5))
plotMA(fit2_N_H, ylim=c(-5,5))
plotMA(fit2_N_S, ylim=c(-5,5))
dev.off()

#logFC?̎Z?o
HC_mean <- apply(normalized_data[,1:24], 1, mean)
SS_mean <- apply(normalized_data[,25:44], 1, mean)
NASH_mean <- apply(normalized_data[,45:63], 1, mean)
logFC_S_H <- SS_mean-HC_mean
logFC_N_H <- NASH_mean-HC_mean
logFC_N_S <- NASH_mean-SS_mean

#Probe_ID, logFC, q-value??matrix???쐬
S_H <- cbind(logFC_S_H, qval_S_H$qvalues)
N_H <- cbind(logFC_N_H, qval_N_H$qvalues)
N_S <- cbind(logFC_N_S, qval_N_S$qvalues)
S_H_N_H_N_S <- cbind(S_H, N_H, N_S)
colnames(S_H_N_H_N_S) <- c("S_H_logFC", "S_H_q-value",
                           "N_H_logFC", "N_H_q-value",
                           "N_S_logFC", "N_S_q-value")
#?????o??
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
write.csv(S_H_N_H_N_S, "series_matrix_logFC_qvalue.csv")

#DEG???̊m?F
SS_HC_up <- sum(S_H_N_H_N_S[,1]>=1 & S_H_N_H_N_S[,2]<0.05) #207
SS_HC_down <- sum(S_H_N_H_N_S[,1]<=(-1) & S_H_N_H_N_S[,2]<0.05) #412
NASH_HC_up <- sum(S_H_N_H_N_S[,3]>=1 & S_H_N_H_N_S[,4]<0.05) #271
NASH_HC_down <- sum(S_H_N_H_N_S[,3]<=(-1) & S_H_N_H_N_S[,4]<0.05) #369
NASH_SS_up <- sum(S_H_N_H_N_S[,5]>=1 & S_H_N_H_N_S[,6]<0.05) #15
NASH_SS_down <- sum(S_H_N_H_N_S[,5]<=(-1) & S_H_N_H_N_S[,6]<0.05) #2

#Volcano plot
q_S_H <- -log10(S_H_N_H_N_S[,2])
q_N_H <- -log10(S_H_N_H_N_S[,4])
q_N_S <- -log10(S_H_N_H_N_S[,6])
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
pdf("series_matrix_volcanoplot.pdf")
plot(S_H_N_H_N_S[,1], q_S_H,
     xlab = "Log2 fold change", ylab = "-Log10 q-value", main = "SS/HC",
     xlim = c(-5,5), ylim = c(0,20))
abline(v=1, lty=2)
abline(v=-1, lty=2)
abline(h=-log10(0.05), lty=2)

plot(S_H_N_H_N_S[,3], q_N_H,
     xlab = "Log2 fold change", ylab = "-Log10 q-value", main = "NASH/HC",
     xlim = c(-5,5), ylim = c(0,20))
abline(v=1, lty=2)
abline(v=-1, lty=2)
abline(h=-log10(0.05), lty=2)

plot(S_H_N_H_N_S[,5], q_N_S,
     xlab = "Log2 fold change", ylab = "-Log10 q-value", main = "NASH/SS",
     xlim = c(-5,5), ylim = c(0,20))
abline(v=1, lty=2)
abline(v=-1, lty=2)
abline(h=-log10(0.05), lty=2)
dev.off()

#?A?m?e?[?V?????t?@?C?????ǂݍ???
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
annot <- read.csv("GeneAnnotation_GPL14951.csv", header = T, stringsAsFactors = F)
setwd("E:/IPA_gene_annotation")
IPAannot <- read.csv("GPL14951_IPAannot.csv", header = T, stringsAsFactors = F)
IPAannot <- IPAannot[,-c(2:4)]

#normalized_data??logFC_qvalue??????
normalized_data_logFC_qvalue <- cbind(normalized_data, S_H_N_H_N_S)
normalized_data_logFC_qvalue <- cbind(rownames(normalized_data), normalized_data_logFC_qvalue)
names(normalized_data_logFC_qvalue)[1] <- "ID"
rownames(normalized_data_logFC_qvalue) <- NULL

#Probe ID??key?ɃA?m?e?[?V?????t?@?C????merge
normalized_data_logFC_qvalue_annot <- merge(annot, normalized_data_logFC_qvalue, by.x = 1, by.y = 1, all.y = T)
normalized_data_logFC_qvalue_annot <- merge(normalized_data_logFC_qvalue_annot, IPAannot, by.x = 1, by.y = 1, all.x = T)

#???`?q???????A?m?e?[?V?????????t?@?C?????????o??
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
write.csv(normalized_data_logFC_qvalue_annot, "normalized_data_logFC_qvalue_annot.csv", row.names=FALSE)






7_series_matrix_pre_WGCNA.R





#series matrix???p????WGCNA
#WGCNA?Ɏg?p???锭???f?[?^?y??trait?f?[?^???쐬????
library(WGCNA)
options(stringsAsFactors = FALSE)

#???K?????f?[?^?̓ǂݍ???
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
normalized_data0 <- read.table("GSE89632_series_matrix.txt", header = T, sep = "\t", row.names = 1, stringsAsFactors = F)

#???͂Ɏg?p?????v???[?u?̑I??
#WGCNA?ɂ́C10000???`?q???x???K???Ă???
#?????ʂ??Ⴂ?v???[?u(<2^9)?????؂??ɂ?????
cutoff_filter <- rowSums(normalized_data0 < log2(512)) < 1
normalized_data <-normalized_data0[cutoff_filter,]
dim(normalized_data0)?@#29377?s
dim(normalized_data)?@ #12481?s

#WGCNA???͗p?ɁC?s???T???v???C????probe??matrix?ɓ]?u
datExpr0 <- t(normalized_data)

#?????l?̃`?F?b?N
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
if(!gsg$allOK)
{
  if(sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if(sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

#???炩?ȊO???l???Ȃ????`?F?b?N
sampleTree <- hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
pdf(file = "Plots_sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()

#HC_4, HC_16, SS_15???O???l?Ɣ??f
#PCA score plot??MDS plot?Ō??Ă??O???Ă???
#height=80?ŃJ?b?g?I?t
sizeGrWindow(12,9)
pdf(file = "Plots_sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = 80, col = "red");
dev.off()
clust <- cutreeStatic(sampleTree, cutHeight = 80, minSize = 10)
table(clust)
#clust 1???c???CdatExpr??WGCNA?̑ΏۂƂ???
keepSamples <- (clust==1)
datExpr <- datExpr0[keepSamples,]
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

#WGCNA???s???f?[?^???????o??
normalized_data_preWGCNA <- t(datExpr)
write.csv(normalized_data_preWGCNA, "normalized_data_preWGCNA.csv", quote = F, row.names = T)

#?Տ?trait?f?[?^?̓ǂݍ???
traitData <- read.csv("GSE89632_trait_data.csv")
dim(traitData)
names(traitData)
#?K?v?̂Ȃ?"Number","GEO_accession"?f?[?^???폜
allTraits <- traitData[,-c(2,3)]
dim(allTraits)
names(allTraits)

#WGCNA?Ŏg?p???锭???f?[?^?ƃT???v??????????
samples <- rownames(datExpr)
traitRows <- match(samples, allTraits$Sample_ID)
datTraits <- allTraits[traitRows, -1]
rownames(datTraits) <- allTraits[traitRows, 1];
collectGarbage();

#?T???v???f???h???O?????ƑΉ??????Տ?trait?̊֘A?????o??
#?T???v?????ăN???X?^?????O
sampleTree2 <- hclust(dist(datExpr), method = "average")
#trait???F?ŕ\???F??_???`??_???C?D?͌???
traitColors <- numbers2colors(datTraits, signed = FALSE);
#?T???v???f???h???O??????trait?̃q?[?g?}?b?v???\??
pdf("Sample_dendrogram_trait_heatmap.pdf", width = 12, height = 9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

#WGCNA?Ɏg?p???锭???f?[?^datExpr?y??trait?f?[?^datTraits???ۑ?
save(datExpr, datTraits, file = "dataInput.Rdata")





8_series_matrix_pre_WGCNA_QC_DEGselection.R





#WGCNA?Ɏg?p?????f?[?^??QC
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
#?p?b?P?[?W?̓ǂݍ???
library(limma)
library(stringr)
library(colorhcplot)
library(qvalue)
#WGCNA?Ŏg?p???鐳?K?????f?[?^?̓ǂݍ???
#Probe_ID?̗????ł????悤??row.names=1?͎g???Ȃ?
normalized_data <- read.csv("normalized_data_preWGCNA.csv", row.names = 1, header = T)

#box plot
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/WGCNA")
pdf("Sample_boxplot.pdf", width = 9.5, height = 7, pointsize = 14)
par(mar = c(6, 3, 3, 1))
boxplot(normalized_data, main="Box plot", las=2, cex.axis = 0.8)
dev.off()

#?V?O?i???l?̕??z???v???b?g
pdf("Sample_Plotdensities.pdf", width = 9.5, height = 7, pointsize = 14)
plotDensities(normalized_data, main="Plot densities", legend=FALSE)
dev.off()

#MDS?v???b?g
pdf("Sample_MDSplot.pdf")
plotMDS(normalized_data, labels = c(1:60), col = c(rep("green3", 22), rep("blue", 19), rep("red", 19)))
cols <- c("green3", "blue", "red")
legend_labels <- c("HC", "SS", "NASH")
legend("bottomleft", legend = legend_labels, text.col = cols)
dev.off()

#?听??????PCA
#PCA???͗p?ɁC?s?F?T???v?????C???FProbe_ID?ɂȂ??悤?]?u
data_mat <- t(normalized_data)

#?m?F
class(data_mat)
dim(data_mat)
data_mat[,1:5]

#PCA?????s
pca_res <- prcomp(data_mat, scale. = T)
summary(pca_res)
pca_res$x

#PC1??PC2?̃X?R?A?v???b?g
pdf("Sample_PCA_score_plot.pdf")
plot(pca_res$x[,1], pca_res$x[,2],
     type = "n", xlab = "PC1 (20.46%)", ylab = "PC2(13.03%)", main = "PCA score plot")
text(pca_res$x[,1], pca_res$x[,2], labels = c(1:60), col = c(rep("green3", 22), rep("blue", 19), rep("red", 19)))
cols <- c("green3", "blue", "red")
legend_labels <- c("HC", "SS", "NASH")
legend("bottomleft", legend = legend_labels, text.col = cols)
dev.off()

#normalized_data?ɂ??????T???v???Ԃ̋????s?????v?Z
d_sample <- dist(t(normalized_data))
class(d_sample)
#?T???v???̃N???X?^?????O
c_ward_sample <- hclust(d_sample, "ward.D2")
c_ave_sample <- hclust(d_sample, "average")
#?T???v???̃f???h???O???????\??
pdf("Sample_dendrogram.pdf", width = 9.5, height = 7, pointsize = 14)
fac <- as.factor(c(rep("HC", 22), rep("SS", 19), rep("NASH", 19)))
fac <- factor(fac, levels = c("HC", "SS", "NASH"))
colorhcplot(hc = c_ward_sample, fac = fac, col = c("green3", "blue", "red"), hang = -1, lab.cex = 0.8, main = "Cluster Dendrogram (Ward)")
colorhcplot(hc = c_ave_sample, fac = fac, col = c("green3", "blue", "red"), hang = -1, lab.cex = 0.8, main = "Cluster Dendrogram (Average)")
dev.off()


#WGCNA?Ɏg?p?????f?[?^??diagnosis?ɑ΂???DEG?I??
#?^?[?Q?b?g?t?@?C???̓ǂݍ???
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
targets <- readTargets("Targets_WGCNA.txt")

#?f?U?C???s?????쐬
#factor()?ŃT???v?????Q???ŃO???[?v??
f <- factor(targets$Condition, levels = unique(targets$Condition))
#?f?U?C???s????design?Ƃ???
design <- model.matrix(~0+f)
#?f?U?C???s???̗??????Q???ɕύX
colnames(design) <- levels(f)
#???`???f?????K?p
fit <- lmFit(normalized_data, design)

#HC??SS?̌???
#HC??SS?̔??r?s?????쐬
contrasts_S_H <- makeContrasts(SS - HC, levels = design)
#???r?s???????f???f?[?^fit_S_H?ɓK?p
fit_S_H <- contrasts.fit(fit, contrasts_S_H)
#?o???x?C?Y?@?Ō???
fit2_S_H <- eBayes(fit_S_H)
#qvalue?̎Z?o
for_qval_S_H <- fit2_S_H$F.p.value
qval_S_H <- qvalue(p = for_qval_S_H)

#HC??NASH?̌???
#HC??NASH?̔??r?s?????쐬
contrasts_N_H <- makeContrasts(NASH - HC, levels = design)
#???r?s???????f???f?[?^fit_N_H?ɓK?p
fit_N_H <- contrasts.fit(fit, contrasts_N_H)
#?o???x?C?Y?@?Ō???
fit2_N_H <- eBayes(fit_N_H)
#qvalue?̎Z?o
for_qval_N_H <- fit2_N_H$F.p.value
qval_N_H <- qvalue(p = for_qval_N_H)

#SS??NASH?̌???
#SS??NASH?̔??r?s?????쐬
contrasts_N_S <- makeContrasts(NASH - SS, levels = design)
#???r?s???????f???f?[?^fit_N_S?ɓK?p
fit_N_S <- contrasts.fit(fit, contrasts_N_S)
#?o???x?C?Y?@?Ō???
fit2_N_S <- eBayes(fit_N_S)
#qvalue?̎Z?o
for_qval_N_S <- fit2_N_S$F.p.value
qval_N_S <- qvalue(p = for_qval_N_S)

#MA plot
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/WGCNA")
pdf("Diagnosis_MAplot.pdf")
plotMA(fit2_S_H, ylim=c(-5,5))
plotMA(fit2_N_H, ylim=c(-5,5))
plotMA(fit2_N_S, ylim=c(-5,5))
dev.off()

#logFC?̎Z?o
HC_mean <- apply(normalized_data[,1:22], 1, mean)
SS_mean <- apply(normalized_data[,23:41], 1, mean)
NASH_mean <- apply(normalized_data[,42:60], 1, mean)
logFC_S_H <- SS_mean-HC_mean
logFC_N_H <- NASH_mean-HC_mean
logFC_N_S <- NASH_mean-SS_mean

#Probe_ID, logFC, q-value??matrix???쐬
S_H <- cbind(logFC_S_H, qval_S_H$qvalues)
N_H <- cbind(logFC_N_H, qval_N_H$qvalues)
N_S <- cbind(logFC_N_S, qval_N_S$qvalues)
S_H_N_H_N_S <- cbind(S_H, N_H, N_S)
colnames(S_H_N_H_N_S) <- c("S_H_logFC", "S_H_q-value",
                           "N_H_logFC", "N_H_q-value",
                           "N_S_logFC", "N_S_q-value")
#?????o??
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/WGCNA")
write.csv(S_H_N_H_N_S, "Diagnosis_logFC_qvalue.csv")

#DEG???̊m?F
SS_HC_up <- sum(S_H_N_H_N_S[,1]>=1 & S_H_N_H_N_S[,2]<0.05) #65
SS_HC_down <- sum(S_H_N_H_N_S[,1]<=(-1) & S_H_N_H_N_S[,2]<0.05) #195
NASH_HC_up <- sum(S_H_N_H_N_S[,3]>=1 & S_H_N_H_N_S[,4]<0.05) #104
NASH_HC_down <- sum(S_H_N_H_N_S[,3]<=(-1) & S_H_N_H_N_S[,4]<0.05) #190
NASH_SS_up <- sum(S_H_N_H_N_S[,5]>=1 & S_H_N_H_N_S[,6]<0.05) #3
NASH_SS_down <- sum(S_H_N_H_N_S[,5]<=(-1) & S_H_N_H_N_S[,6]<0.05) #1

#Volcano plot
q_S_H <- -log10(S_H_N_H_N_S[,2])
q_N_H <- -log10(S_H_N_H_N_S[,4])
q_N_S <- -log10(S_H_N_H_N_S[,6])
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/WGCNA")
pdf("Diagnosis_volcanoplot.pdf")
plot(S_H_N_H_N_S[,1], q_S_H,
     xlab = "Log2 fold change", ylab = "-Log10 q-value", main = "SS/HC",
     xlim = c(-5,5), ylim = c(0,20))
abline(v=1, lty=2)
abline(v=-1, lty=2)
abline(h=-log10(0.05), lty=2)

plot(S_H_N_H_N_S[,3], q_N_H,
     xlab = "Log2 fold change", ylab = "-Log10 q-value", main = "NASH/HC",
     xlim = c(-5,5), ylim = c(0,20))
abline(v=1, lty=2)
abline(v=-1, lty=2)
abline(h=-log10(0.05), lty=2)

plot(S_H_N_H_N_S[,5], q_N_S,
     xlab = "Log2 fold change", ylab = "-Log10 q-value", main = "NASH/SS",
     xlim = c(-5,5), ylim = c(0,20))
abline(v=1, lty=2)
abline(v=-1, lty=2)
abline(h=-log10(0.05), lty=2)
dev.off()

#?A?m?e?[?V?????t?@?C?????ǂݍ???
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
annot <- read.csv("GeneAnnotation_GPL14951.csv", header = T, stringsAsFactors = F)
setwd("E:/IPA_gene_annotation")
IPAannot <- read.csv("GPL14951_IPAannot.csv", header = T, stringsAsFactors = F)
IPAannot <- IPAannot[,-c(2:4)]

#normalized_data??logFC_qvalue??????
normalized_data_logFC_qvalue <- cbind(normalized_data, S_H_N_H_N_S)
normalized_data_logFC_qvalue <- cbind(rownames(normalized_data), normalized_data_logFC_qvalue)
names(normalized_data_logFC_qvalue)[1] <- "ID"
rownames(normalized_data_logFC_qvalue) <- NULL

#Probe ID??key?ɃA?m?e?[?V?????t?@?C????merge
normalized_data_logFC_qvalue_annot <- merge(annot, normalized_data_logFC_qvalue, by.x = 1, by.y = 1, all.y = T)
normalized_data_logFC_qvalue_annot <- merge(normalized_data_logFC_qvalue_annot, IPAannot, by.x = 1, by.y = 1, all.x = T)

#???`?q???????A?m?e?[?V?????????t?@?C?????????o??
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/WGCNA")
write.csv(normalized_data_logFC_qvalue_annot, "Diagnosis_normalized_data_logFC_qvalue_annot.csv", row.names=FALSE)





10_series_matrix_module_trait_relationship.R





#series matrix???p????WGCNA
#???W???[????trait?̊֘A??#
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
library(WGCNA)
options(stringsAsFactors = FALSE)
lnames <- load("dataInput.Rdata")
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
annot <- read.csv("GeneAnnotation_GPL14951.csv", header = T, stringsAsFactors = F)

#deepSplit??0?`4?????ꂼ?????͂??C???L?̃R?[?h???J???Ԃ????s

#########????????#########
for(i in 0:4){
  deepSplit <- i #???̐??l??0?`4?ɕς???

##unsigned_merge?Ȃ?
#?t?H???_?ړ?
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("networkConstruction_StepByStep_unsigned_", deepSplit, ".Rdata", sep = ""))

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
pdf(paste("Module_trait_relationships_unsigned_", deepSplit, ".pdf",sep = ""), width = 10, height = 18)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships (unsigned, deepSplit ", deepSplit, ")", sep = ""))
dev.off()

#???W???[????trait?̑??֌W???y??p-value???\?ɂ??ĕۑ?
colnames(textMatrix) <- names(datTraits)
rownames(textMatrix) <- names(MEs)
write.csv(textMatrix, paste("Module_trait_correlation_unsigned_", deepSplit, ".csv", sep = ""))

#???ڂ???trait?i???̏ꍇ??Fibrosis?j?Ɋւ???Gene significance (GS)??Module membership (MM)???Z?o????
fibrosis <- as.data.frame(datTraits$Fibrosis)
names(fibrosis) <-"fibrosis"
modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

geneTraitSignificance <- as.data.frame(cor(datExpr, fibrosis, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(fibrosis), sep = "")
names(GSPvalue) <- paste("p.GS.", names(fibrosis), sep = "")

#Probe ID??gene name?ɕϊ?
probes <- colnames(datExpr)
probes2annot <- match(probes, annot$ID)
#?A?m?e?[?V?????????Ă??Ȃ??v???[?u????0?ł??邩?m?F
sum(is.na(probes2annot)) 

#?A?m?e?[?V?????????y??fibrosis?Ɋւ???GS?̕\???쐬
geneInfo0 <- data.frame(Probe_ID = probes,
                        geneSymbol = annot$Symbol[probes2annot],
                        Entrez_gene_ID = annot$Entrez_Gene_ID[probes2annot],
                        moduleColor = moduleColors_unsigned,
                        geneTraitSignificance,
                        GSPvalue)
#fibrosis??p?l?ŕ??בւ?
modOrder <- order(-abs(cor(MEs, fibrosis, use = "p")))
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
setwd("E:/IPA_gene_annotation")
IPAannot <- read.csv("GPL14951_IPAannot.csv", header = T, stringsAsFactors = F)
IPAannot <- IPAannot[,-c(2:4)]
geneInfo0 <- merge(geneInfo0, IPAannot, by.x = 1, by.y = 1, all.x = T)

#module color?y??GS?ŕ??בւ?
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.fibrosis))
geneInfo <- geneInfo0[geneOrder, ]

#?A?m?e?[?V?????????y??fibrosis?Ɋւ???GS/MM?̕\???????o??
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, sep = ""))
write.csv(geneInfo, paste("geneInfo_unsigned_", deepSplit, ".csv",sep = "" ), row.names = FALSE)


##unsigned_merge????
#?t?H???_?ړ?
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("networkConstruction_StepByStep_unsigned_", deepSplit, "_merged.Rdata", sep = ""))

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
pdf(paste("Module_trait_relationships_unsigned_", deepSplit, "_merged.pdf",sep = ""), width = 10, height = 18)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships (unsigned, deepSplit ", deepSplit, ", merged)", sep = ""))
dev.off()

#???W???[????trait?̑??֌W???y??p-value???\?ɂ??ĕۑ?
colnames(textMatrix) <- names(datTraits)
rownames(textMatrix) <- names(MEs)
write.csv(textMatrix, paste("Module_trait_correlation_unsigned_", deepSplit, "_merged.csv", sep = ""))

#???ڂ???trait?i???̏ꍇ??Fibrosis?j?Ɋւ???Gene significance (GS)??Module membership (MM)???Z?o????
fibrosis <- as.data.frame(datTraits$Fibrosis)
names(fibrosis) <-"fibrosis"
modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

geneTraitSignificance <- as.data.frame(cor(datExpr, fibrosis, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(fibrosis), sep = "")
names(GSPvalue) <- paste("p.GS.", names(fibrosis), sep = "")

#Probe ID??gene name?ɕϊ?
probes <- colnames(datExpr)
probes2annot <- match(probes, annot$ID)
#?A?m?e?[?V?????????Ă??Ȃ??v???[?u????0?ł??邩?m?F
sum(is.na(probes2annot)) 

#?A?m?e?[?V?????????y??fibrosis?Ɋւ???GS?̕\???쐬
geneInfo0 <- data.frame(Probe_ID = probes,
                        geneSymbol = annot$Symbol[probes2annot],
                        Entrez_gene_ID = annot$Entrez_Gene_ID[probes2annot],
                        moduleColor = moduleColors_unsigned,
                        geneTraitSignificance,
                        GSPvalue)
#fibrosis??p?l?ŕ??בւ?
modOrder <- order(-abs(cor(MEs, fibrosis, use = "p")))
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
setwd("E:/IPA_gene_annotation")
IPAannot <- read.csv("GPL14951_IPAannot.csv", header = T, stringsAsFactors = F)
IPAannot <- IPAannot[,-c(2:4)]
geneInfo0 <- merge(geneInfo0, IPAannot, by.x = 1, by.y = 1, all.x = T)

#module color?y??GS?ŕ??בւ?
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.fibrosis))
geneInfo <- geneInfo0[geneOrder, ]
#?A?m?e?[?V?????????y??fibrosis?Ɋւ???GS/MM?̕\???????o??
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, sep = ""))
write.csv(geneInfo, paste("geneInfo_unsigned_", deepSplit, "_merged.csv",sep = "" ), row.names = FALSE)

##signed_merge?Ȃ?
#?t?H???_?ړ?
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("networkConstruction_StepByStep_signed_", deepSplit, ".Rdata", sep = ""))

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
pdf(paste("Module_trait_relationships_signed_", deepSplit, ".pdf",sep = ""), width = 10, height = 18)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships (signed, deepSplit ", deepSplit, ")", sep = ""))
dev.off()

#???W???[????trait?̑??֌W???y??p-value???\?ɂ??ĕۑ?
colnames(textMatrix) <- names(datTraits)
rownames(textMatrix) <- names(MEs)
write.csv(textMatrix, paste("Module_trait_correlation_signed_", deepSplit, ".csv", sep = ""))

#???ڂ???trait?i???̏ꍇ??Fibrosis?j?Ɋւ???Gene significance (GS)??Module membership (MM)???Z?o????
fibrosis <- as.data.frame(datTraits$Fibrosis)
names(fibrosis) <-"fibrosis"
modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

geneTraitSignificance <- as.data.frame(cor(datExpr, fibrosis, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(fibrosis), sep = "")
names(GSPvalue) <- paste("p.GS.", names(fibrosis), sep = "")

#Probe ID??gene name?ɕϊ?
probes <- colnames(datExpr)
probes2annot <- match(probes, annot$ID)
#?A?m?e?[?V?????????Ă??Ȃ??v???[?u????0?ł??邩?m?F
sum(is.na(probes2annot)) 

#?A?m?e?[?V?????????y??fibrosis?Ɋւ???GS?̕\???쐬
geneInfo0 <- data.frame(Probe_ID = probes,
                        geneSymbol = annot$Symbol[probes2annot],
                        Entrez_gene_ID = annot$Entrez_Gene_ID[probes2annot],
                        moduleColor = moduleColors_signed,
                        geneTraitSignificance,
                        GSPvalue)
#fibrosis??p?l?ŕ??בւ?
modOrder <- order(-abs(cor(MEs, fibrosis, use = "p")))
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
setwd("E:/IPA_gene_annotation")
IPAannot <- read.csv("GPL14951_IPAannot.csv", header = T, stringsAsFactors = F)
IPAannot <- IPAannot[,-c(2:4)]
geneInfo0 <- merge(geneInfo0, IPAannot, by.x = 1, by.y = 1, all.x = T)

#module color?y??GS?ŕ??בւ?
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.fibrosis))
geneInfo <- geneInfo0[geneOrder, ]
#?A?m?e?[?V?????????y??fibrosis?Ɋւ???GS/MM?̕\???????o??
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
write.csv(geneInfo, paste("geneInfo_signed_", deepSplit, ".csv",sep = "" ), row.names = FALSE)


##signed_merge????
#?t?H???_?ړ?
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("networkConstruction_StepByStep_signed_", deepSplit, "_merged.Rdata", sep = ""))

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
pdf(paste("Module_trait_relationships_signed_", deepSplit, "_merged.pdf",sep = ""), width = 10, height = 18)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships (signed, deepSplit ", deepSplit, ", merged)", sep = ""))
dev.off()

#???W???[????trait?̑??֌W???y??p-value???\?ɂ??ĕۑ?
colnames(textMatrix) <- names(datTraits)
rownames(textMatrix) <- names(MEs)
write.csv(textMatrix, paste("Module_trait_correlation_signed_", deepSplit, "_merged.csv", sep = ""))

#???ڂ???trait?i???̏ꍇ??Fibrosis?j?Ɋւ???Gene significance (GS)??Module membership (MM)???Z?o????
fibrosis <- as.data.frame(datTraits$Fibrosis)
names(fibrosis) <-"fibrosis"
modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

geneTraitSignificance <- as.data.frame(cor(datExpr, fibrosis, use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(fibrosis), sep = "")
names(GSPvalue) <- paste("p.GS.", names(fibrosis), sep = "")

#Probe ID??gene name?ɕϊ?
probes <- colnames(datExpr)
probes2annot <- match(probes, annot$ID)
#?A?m?e?[?V?????????Ă??Ȃ??v???[?u????0?ł??邩?m?F
sum(is.na(probes2annot)) 

#?A?m?e?[?V?????????y??fibrosis?Ɋւ???GS?̕\???쐬
geneInfo0 <- data.frame(Probe_ID = probes,
                        geneSymbol = annot$Symbol[probes2annot],
                        Entrez_gene_ID = annot$Entrez_Gene_ID[probes2annot],
                        moduleColor = moduleColors_signed,
                        geneTraitSignificance,
                        GSPvalue)
#fibrosis??p?l?ŕ??בւ?
modOrder <- order(-abs(cor(MEs, fibrosis, use = "p")))
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
setwd("E:/IPA_gene_annotation")
IPAannot <- read.csv("GPL14951_IPAannot.csv", header = T, stringsAsFactors = F)
IPAannot <- IPAannot[,-c(2:4)]
geneInfo0 <- merge(geneInfo0, IPAannot, by.x = 1, by.y = 1, all.x = T)

#module color?y??GS?ŕ??בւ?
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.fibrosis))
geneInfo <- geneInfo0[geneOrder, ]
#?A?m?e?[?V?????????y??fibrosis?Ɋւ???GS/MM?̕\???????o??
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
write.csv(geneInfo, paste("geneInfo_signed_", deepSplit, "_merged.csv",sep = "" ), row.names = FALSE)
}
#########?????܂?#########


#???ڂ???trait?ƍł????ւ??郂?W???[???i???̏ꍇcyan?j?ɂ????āC??GS??MM?̈??`?q(probe)?????肷??
trait <- "fibrosis"
module <- "cyan"
column <- match(module, modNames)
moduleGenes <- moduleColors_unsigned == module

#???ڃ??W???[???ɂ?????GS-MM?v???b?g
pdf(paste("GS_MM_unsigned_", deepSplit, "_", trait, "_", module, ".pdf",sep = ""), width = 7, height = 7)
par(mfrow = c(1, 1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Gene significance for", trait),
                       main = paste("Module membership vs. gene significance", "(unsigned_", deepSplit, ")\n", sep = ""),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()





11_series_matrix_module_GOEnrichment.R





#series matrix???p????WGCNA
#???W???[????GO enrichment????#
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
library(WGCNA)
options(stringsAsFactors = FALSE)
lnames <- load("dataInput.Rdata")
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
annot <- read.csv("GeneAnnotation_GPL14951.csv", header = T, stringsAsFactors = F)

#deepSplit??0?`4?????ꂼ?????͂??Cfor()?ɂ??艺?L?̃R?[?h???J???Ԃ????s
for(i in 0:4){
deepSplit <- i #???̐??l??0?`4?ɕς???

##unsigned_merge?Ȃ?
#?t?H???_?ړ?
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("networkConstruction_StepByStep_unsigned_", deepSplit, ".Rdata", sep = ""))
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
write.csv(tab, paste("GOEnrichmentTable_unsigned_", deepSplit, ".csv",sep = "" ), quote = TRUE, row.names = FALSE)


##unsigned_merge????
#?t?H???_?ړ?
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("networkConstruction_StepByStep_unsigned_", deepSplit, "_merged.Rdata", sep = ""))
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
write.csv(tab, paste("GOEnrichmentTable_unsigned_", deepSplit, "_merged.csv",sep = "" ), quote = TRUE, row.names = FALSE)


##signed_merge?Ȃ?
#?t?H???_?ړ?
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("networkConstruction_StepByStep_signed_", deepSplit, ".Rdata", sep = ""))
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
write.csv(tab, paste("GOEnrichmentTable_signed_", deepSplit, ".csv",sep = "" ), quote = TRUE, row.names = FALSE)


##signed_merge????
#?t?H???_?ړ?
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("networkConstruction_StepByStep_signed_", deepSplit, "_merged.Rdata", sep = ""))
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
write.csv(tab, paste("GOEnrichmentTable_signed_", deepSplit, "_merged.csv",sep = "" ), quote = TRUE, row.names = FALSE)
}





12_series_matrix_module_CellTypeEnrichment.R





#series matrix???p????WGCNA
#???W???[????Cell-Type enrichment????#
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
library(WGCNA)
library(anRichment)
library(tidyr)
library(ggplot2)
options(stringsAsFactors = FALSE)
lnames <- load("dataInput.Rdata")
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
annot <- read.csv("GeneAnnotation_GPL14951.csv", header = T, stringsAsFactors = F)


#?q?g?̍זE?ʈ??`?q???X?g?̍쐬
#???肽??Collection?????w?肷??
collection <- newCollection()
#Organism???w?肷??
organism <- "Human"

#?t?H???_????Genelist???ꊇ?œǂݍ??܂??邽?߂̐ݒ?
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data/Human_Liver_CellType")
#?t?@?C???i?[?p?X???w??
files <- list.files(path = "//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data/Human_Liver_CellType",full.names = T)
files_names <- list.files( , full.names = F, pattern="csv") 
files_names <- gsub(".csv", "", files_names) 

#Genelist??Geneset?ɂ??邽?߂̊֐????ǂݍ??܂???
Togeneset <- function (X, Y){
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
?convert2entrez

#?t?H???_????Genelist???ꊇ??geneset?ɕϊ????Ă??????L??Collection?ɓ?????
for(i in 1:length(files)){
  genes <- read.csv(files[i], fileEncoding="UTF-8")#?????????t?@?C?????ǂݍ???
  genes2<- genes$Symbol
  assign(files_names[i], genes2)
  genes <- unique(genes)#?d???L???ꍇ????????
  genes_set <- Togeneset(genes, files_names[i])#genelist??geneset??
  assign(paste(files_names[i]), genes_set)
  collection <- addToCollection(collection, genes_set)#Geneset?????ꂽ??collection?ɓ?????
  rm(genes, genes_set)
}

#???`?q?Z?b?g?ۑ?
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
save(collection, file="Human_Liver_CellType_collection.Rdata")
load("Human_Liver_CellType_collection.Rdata")


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
    refCollection = collection,  #reference?Ƃ??Đݒ肷???f?[?^?Z?b?g?B
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


#???ۂɉ???

#deepSplit??0?`4?????ꂼ?????͂??Cfor()?ɂ??艺?L?̃R?[?h???J???Ԃ????s
for(i in 0:4){
  deepSplit <- i #???̐??l??0?`4?ɕς???
  
  ##unsigned_merge?Ȃ?
  #?t?H???_?ړ?
  setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, sep = ""))
  #?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
  lnames <- load(paste("networkConstruction_StepByStep_unsigned_", deepSplit, ".Rdata", sep = ""))
  #Probe ID??Entrez Gene ID?ɕϊ?????
  probes <- colnames(datExpr)
  probes2annot <- match(probes, annot$ID)
  #Entrez Gene ID?iLocusLinkID?j???ǂݍ???
  allLLIDs <- annot$Entrez_Gene_ID[probes2annot]
  #WGCNA?ɂ??????`?q?ƃ??W???[???̕\???쐬
  WGCNAgenelist <- data.frame(allLLIDs, moduleColors_unsigned)
  #?t?H???_?쐬?ƈړ?
  ifelse(!dir.exists(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, "/CellType", sep = "")), 
         dir.create(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, "/CellType", sep = "")), FALSE)
  setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, "/CellType", sep = ""))
  #???W???[????Cell-Type enrichment????
  CelltypeEnrichment(WGCNAgenelist)


  ##unsigned_merge????
  #?t?H???_?ړ?
  setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, sep = ""))
  #?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
  lnames <- load(paste("networkConstruction_StepByStep_unsigned_", deepSplit, "_merged.Rdata", sep = ""))
  #Probe ID??Entrez Gene ID?ɕϊ?????
  probes <- colnames(datExpr)
  probes2annot <- match(probes, annot$ID)
  #Entrez Gene ID?iLocusLinkID?j???ǂݍ???
  allLLIDs <- annot$Entrez_Gene_ID[probes2annot]
  #WGCNA?ɂ??????`?q?ƃ??W???[???̕\???쐬
  WGCNAgenelist <- data.frame(allLLIDs, moduleColors_unsigned)
  #?t?H???_?쐬?ƈړ?
  ifelse(!dir.exists(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, "/CellType_merged", sep = "")), 
         dir.create(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, "/CellType_merged", sep = "")), FALSE)
  setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/unsigned_", deepSplit, "/CellType_merged", sep = ""))
  #???W???[????Cell-Type enrichment????
  CelltypeEnrichment(WGCNAgenelist)


  ##signed_merge?Ȃ?
  #?t?H???_?ړ?
  setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
  #?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
  lnames <- load(paste("networkConstruction_StepByStep_signed_", deepSplit, ".Rdata", sep = ""))
  #Probe ID??Entrez Gene ID?ɕϊ?????
  probes <- colnames(datExpr)
  probes2annot <- match(probes, annot$ID)
  #Entrez Gene ID?iLocusLinkID?j???ǂݍ???
  allLLIDs <- annot$Entrez_Gene_ID[probes2annot]
  #WGCNA?ɂ??????`?q?ƃ??W???[???̕\???쐬
  WGCNAgenelist <- data.frame(allLLIDs, moduleColors_signed)
  #?t?H???_?쐬?ƈړ?
  ifelse(!dir.exists(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, "/CellType", sep = "")), 
         dir.create(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, "/CellType", sep = "")), FALSE)
  setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, "/CellType", sep = ""))
  #???W???[????Cell-Type enrichment????
  CelltypeEnrichment(WGCNAgenelist)
  
  
  ##signed_merge????
  #?t?H???_?ړ?
  setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
  #?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
  lnames <- load(paste("networkConstruction_StepByStep_signed_", deepSplit, "_merged.Rdata", sep = ""))
  #Probe ID??Entrez Gene ID?ɕϊ?????
  probes <- colnames(datExpr)
  probes2annot <- match(probes, annot$ID)
  #Entrez Gene ID?iLocusLinkID?j???ǂݍ???
  allLLIDs <- annot$Entrez_Gene_ID[probes2annot]
  #WGCNA?ɂ??????`?q?ƃ??W???[???̕\???쐬
  WGCNAgenelist <- data.frame(allLLIDs, moduleColors_signed)
  #?t?H???_?쐬?ƈړ?
  ifelse(!dir.exists(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, "/CellType_merged", sep = "")), 
         dir.create(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, "/CellType_merged", sep = "")), FALSE)
  setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, "/CellType_merged", sep = ""))
  #???W???[????Cell-Type enrichment????
  CelltypeEnrichment(WGCNAgenelist)
}





13_series_matrix_signed4merge_table.R





#series matrix???p????WGCNA
#signed_deepSplit4_merged?????̃f?[?^????#
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
library(WGCNA)
library(anRichment)
library(tidyr)
library(ggplot2)
library(stringr)
options(stringsAsFactors = FALSE)
lnames <- load("dataInput.Rdata")
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
annot <- read.csv("GeneAnnotation_GPL14951.csv", header = T, stringsAsFactors = F)
#diagnosis??logFC, qvalue?f?[?^???ǂݍ???
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/WGCNA")
logFC_qvalue <- read.csv("Diagnosis_logFC_qvalue.csv", header = T, stringsAsFactors = F, row.names = 1)

deepSplit <- 4
  
##signed_merge????
#?t?H???_?ړ?
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("networkConstruction_StepByStep_signed_", deepSplit, "_merged.Rdata", sep = ""))
kWithin <- read.csv("Probe_kWitnin_color_signed_4_merged.csv", header = T, stringsAsFactors = F)

#Probe ID??Entrez Gene ID?ɕϊ?????
probes <- colnames(datExpr)
probes2annot <- match(probes, annot$ID)
#Probe ID, gene symbol, gene ID, moduleColor, , kWithin, ?????l, logFC, qvalue?̕\???쐬
geneID_color_kWithin_normalizedData_logFC_qvalue <- data.frame(Probe_ID = probes,
                                                  geneSymbol = annot$Symbol[probes2annot],
                                                  Entrez_gene_ID = annot$Entrez_Gene_ID[probes2annot],
                                                  moduleColor = moduleColors_signed,
                                                  kWithin = kWithin$kWithin,
                                                  t(datExpr),
                                                  logFC_qvalue)

#IPA annotation???????ǉ?
setwd("E:/IPA_gene_annotation")
IPAannot <- read.csv("GPL14951_IPAannot.csv", header = T, stringsAsFactors = F)
IPAannot <- IPAannot[,-c(2:4)]
geneID_color_kWithin_normalizedData_logFC_qvalue <- merge(geneID_color_kWithin_normalizedData_logFC_qvalue, IPAannot, by.x = 1, by.y = 1, all.x = T)

#?A?m?e?[?V?????????CmoduleColor?C?????l?CkWithin?̕\???????o??
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
write.csv(geneID_color_kWithin_normalizedData_logFC_qvalue, "geneID_color_kWithin_normalizedData_logFC_qvalue.csv", row.names = FALSE)

#BaseSpace, IPA, MetaCore?p?Ƀ??W???[??????GeneID, kWithin?̕\???쐬
module_color <- str_sub(colnames(mergedMEs_signed), start = 3) #?擪2????ME???폜
len <- length(module_color)
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, "/BaseSpace/import", sep = ""))

for(i in 1:len){
  module <- module_color[i]
  
  geneID_kWithin_symbol_color <- subset(geneID_color_kWithin_normalizedData_logFC_qvalue, geneID_color_kWithin_normalizedData_logFC_qvalue$moduleColor == module)
  geneID_kWithin_symbol_color <- data.frame(geneID_kWithin_symbol_color$Entrez_gene_ID,
                                            geneID_kWithin_symbol_color$kWithin,
                                            geneID_kWithin_symbol_color$geneSymbol,
                                            geneID_kWithin_symbol_color$moduleColor)
  colnames(geneID_kWithin_symbol_color) <- c("Gene Name", "kWithin", "Gene Symbol", "Module Color")
  write.table(geneID_kWithin_symbol_color, paste(module, ".txt", sep=""), sep="\t", append=F, quote=F, row.names = FALSE)
}





14_series_matrix_signed4merge_module_trait_relationship_CellTypeEnrichment.R





#series matrix???p????WGCNA
#signed_deepSplit4_merged?f?[?^
#???W???[????trait?̊֘A??#
#?J?e?S???[?f?[?^??Spearman?C???l?f?[?^??Peason?̑??֌W???ŎZ?o???Ȃ???
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
library(WGCNA)
options(stringsAsFactors = FALSE)
lnames <- load("dataInput.Rdata")
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
annot <- read.csv("GeneAnnotation_GPL14951.csv", header = T, stringsAsFactors = F)
deepSplit <- 4

##signed_merge????
#?t?H???_?ړ?
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("networkConstruction_StepByStep_signed_", deepSplit, "_merged.Rdata", sep = ""))

#???`?q?i?v???[?u?j???ƃT???v?????????`
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

#ME???Čv?Z
MEs <- moduleEigengenes(datExpr, moduleColors_signed)$eigengenes

#Spearman???ւ?p?l???Z?o
moduleTraitCor_Spearman <- cor(MEs, datTraits[,c(1:6, 8, 22)], use = "p", method = "spearman")
moduleTraitPvalue_Spearman <- corPvalueStudent(moduleTraitCor_Spearman, nSamples)

#Pearson???ւ?p?l???Z?o
moduleTraitCor_Pearson <- cor(MEs, datTraits[,c(7, 9:21)], use = "p", method = "pearson")
moduleTraitPvalue_Pearson <- corPvalueStudent(moduleTraitCor_Pearson, nSamples)

#???֌W????p?l???܂Ƃ߂?
moduleTraitCor <- cbind(moduleTraitCor_Spearman, moduleTraitCor_Pearson)
moduleTraitPvalue <- cbind(moduleTraitPvalue_Spearman, moduleTraitPvalue_Pearson)
datTraits <- datTraits[,c(1:6, 8, 22, 7, 9:21)]

#???W???[????trait?̊֘A?????q?[?g?}?b?v?Ŏ???
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
pdf(paste("Module_trait_relationships_signed_", deepSplit, "_merged_SP.pdf",sep = ""), width = 10, height = 16)
par(mar = c(8, 9, 3, 3))
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
dev.off()

#???W???[????trait?̑??֌W???y??p-value???\?ɂ??ĕۑ?
colnames(textMatrix) <- names(datTraits)
rownames(textMatrix) <- names(MEs)
write.csv(textMatrix, paste("Module_trait_correlation_signed_", deepSplit, "_merged_SP.csv", sep = ""))

#?זE?G?????b?`?????g???͌??ʂ??C???l?̃q?[?g?}?b?v?Ŏ???
library(reshape2)
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, "/CellType_merged", sep = ""))
CellTyperesult <- read.csv("WGCNA_enrichment_result.csv")
logP <- -log10(CellTyperesult$pValue)
CellTypeP <- data.frame(Module = CellTyperesult$class, CellType = CellTyperesult$dataSetID, logP = logP)
#wide-format?̃f?[?^?t???[???ɕϊ?
CellTypeP <- dcast(CellTypeP, Module ~ CellType, value.var = "logP")
#?זE?̏??Ԃ????בւ?
CellTypeP <- data.frame(Module = CellTypeP$Module,
                        Hep3 = CellTypeP$Hep3_cluster5, Hep5 = CellTypeP$Hep5_cluster14, Hep4 = CellTypeP$Hep4_cluster6, 
                        Hep6 = CellTypeP$Hep6_cluster15, Hep1 = CellTypeP$Hep1_cluster1, Hep2 = CellTypeP$Hep2_cluster3, 
                        PeriportalLSEC = CellTypeP$PeriportalLSEC_cluster11, CentralVenousLSEC = CellTypeP$CentralVenousLSEC_cluster12, PortalEndothelialCell = CellTypeP$PortalEndothelialCell_cluster13,
                        Cholangiocyte = CellTypeP$Cholangiocyte_cluster17, HepaticStellateCell = CellTypeP$HepaticStellateCell_cluster20,
                        InflammatoryMacs = CellTypeP$InflammatoryMacs_cluster4, NonInflammatoryMacs = CellTypeP$NonInflammatoryMacs_cluster10,
                        CD3abTcell = CellTypeP$CD3abTcell_cluster2, gdTcell1 = CellTypeP$gdTcell1_cluster9, gdTcell2 = CellTypeP$gdTcell2_cluster18, NKlikeCell = CellTypeP$NKlikeCell_cluster8,
                        MatureBcell = CellTypeP$MatureBcell_cluster16, PlasmaCell = CellTypeP$PlasmaCell_cluster7, ErythroidCell = CellTypeP$ErythroidCell_cluster19)
#?זE?G?????b?`?????g???͂Ɉ??????????Ȃ????????W???[?????ǉ?
library(stringr)
color <- colnames(MEs)
color <- str_sub(color, start = 3) #?擪2????ME???폜
dif <- setdiff(color, CellTypeP$Module)  #color?ɂ̂݊܂܂??Ă??镶?????????o
dif <- data.frame(dif, matrix(NA, nrow = 11, ncol = 20)) #color?ɂ̂݊܂܂??Ă??镶??????module?Ƃ??????f?[?^?t???[?????쐬
colnames(dif) <- colnames(CellTypeP)  #?זE?????????ɕt?^
CellTypeP <- rbind(CellTypeP, dif)
rownames(CellTypeP) <- CellTypeP$Module  #color (module)???s???ɕt?^
CellTypeP <- CellTypeP[order(CellTypeP$Module),]  #color???A???t?@?x?b?g???ɕ??בւ?
CellTypeP <- CellTypeP[, -1]  #color (Module)?̗????폜
CellTypeP <- as.matrix(CellTypeP)

#?זE?G?????b?`?????g???͂̃q?[?g?}?b?v
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
textMatrix <- paste(signif(CellTypeP, 2))
dim(textMatrix) <- dim(CellTypeP)
pdf(paste("CellType_pvalue_heatmap_signed_", deepSplit, "_merged.pdf",sep = ""), width = 10, height = 16)
par(mar = c(8, 9, 3, 3))
labeledHeatmap(Matrix = CellTypeP,
               xLabels = colnames(CellTypeP),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(100)[50:100],
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               main = paste("Cell-type enrichment (signed, deepSplit ", deepSplit, ", merged)", sep = ""))
dev.off()




#???ڂ???trait?i???̏ꍇ??Fibrosis?j?Ɋւ???Gene significance (GS)??Module membership (MM)???Z?o????
fibrosis <- as.data.frame(datTraits$Fibrosis)
names(fibrosis) <-"fibrosis"
modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

geneTraitSignificance <- as.data.frame(cor(datExpr, fibrosis, use = "p", method = "spearman"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(fibrosis), sep = "")
names(GSPvalue) <- paste("p.GS.", names(fibrosis), sep = "")

#Probe ID??gene name?ɕϊ?
probes <- colnames(datExpr)
probes2annot <- match(probes, annot$ID)
#?A?m?e?[?V?????????Ă??Ȃ??v???[?u????0?ł??邩?m?F
sum(is.na(probes2annot)) 

#?A?m?e?[?V?????????y??fibrosis?Ɋւ???GS?̕\???쐬
geneInfo0 <- data.frame(Probe_ID = probes,
                        geneSymbol = annot$Symbol[probes2annot],
                        Entrez_gene_ID = annot$Entrez_Gene_ID[probes2annot],
                        moduleColor = moduleColors_signed,
                        geneTraitSignificance,
                        GSPvalue)
#fibrosis??p?l?ŕ??בւ?
modOrder <- order(-abs(cor(MEs, fibrosis, use = "p", method = "spearman")))
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
setwd("E:/IPA_gene_annotation")
IPAannot <- read.csv("GPL14951_IPAannot.csv", header = T, stringsAsFactors = F)
IPAannot <- IPAannot[,-c(2:4)]
geneInfo0 <- merge(geneInfo0, IPAannot, by.x = 1, by.y = 1, all.x = T)
#module color?y??GS?ŕ??בւ?
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.fibrosis))
geneInfo <- geneInfo0[geneOrder, ]
#?A?m?e?[?V?????????y??fibrosis?Ɋւ???GS/MM?̕\???????o??
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
write.csv(geneInfo, paste("geneInfo_fibrosis_signed_", deepSplit, "_merged.csv",sep = "" ), row.names = FALSE)


#???ڂ???trait?ifibrosis?j?ƍł????ւ???"???זE"???W???[???itan?j?ɂ????āC??GS??MM?̈??`?q(probe)?????肷??
trait <- "fibrosis"
module <- "tan"
column <- match(module, modNames)
moduleGenes <- moduleColors_signed == module

#???ڃ??W???[???ɂ?????GS-MM?v???b?g
pdf(paste("GS_MM_signed_", deepSplit, "_", trait, "_", module, "_merged.pdf",sep = ""), width = 7, height = 7)
par(mfrow = c(1, 1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = paste("Gene significance for", trait),
                       main = paste("Module membership vs. gene significance", "(signed_", deepSplit, "_merged)\n", sep = ""),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module, pch = 16)
dev.off()

#???ڂ???trait?ifibrosis?j?ƍł????ւ???"?̍זE?H"???W???[???ithistle1?j?ɂ????āC??GS??MM?̈??`?q(probe)?????肷??
trait <- "fibrosis"
module <- "thistle1"
column <- match(module, modNames)
moduleGenes <- moduleColors_signed == module

#???ڃ??W???[???ɂ?????GS-MM?v???b?g
pdf(paste("GS_MM_signed_", deepSplit, "_", trait, "_", module, "_merged.pdf",sep = ""), width = 7, height = 7)
par(mfrow = c(1, 1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", trait),
                   main = paste("Module membership vs. gene significance", "(signed_", deepSplit, "_merged)\n", sep = ""),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module, pch = 16)
dev.off()



#???ڂ???trait?i???̏ꍇ??Diagnosis?j?Ɋւ???Gene significance (GS)??Module membership (MM)???Z?o????
diagnosis <- as.data.frame(datTraits$Diagnosis)
names(diagnosis) <-"diagnosis"
modNames <- substring(names(MEs), 3)

geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", modNames, sep = "")
names(MMPvalue) <- paste("p.MM", modNames, sep = "")

geneTraitSignificance <- as.data.frame(cor(datExpr, diagnosis, use = "p", method = "spearman"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(diagnosis), sep = "")
names(GSPvalue) <- paste("p.GS.", names(diagnosis), sep = "")

#Probe ID??gene name?ɕϊ?
probes <- colnames(datExpr)
probes2annot <- match(probes, annot$ID)
#?A?m?e?[?V?????????Ă??Ȃ??v???[?u????0?ł??邩?m?F
sum(is.na(probes2annot)) 

#?A?m?e?[?V?????????y??diagnosis?Ɋւ???GS?̕\???쐬
geneInfo0 <- data.frame(Probe_ID = probes,
                        geneSymbol = annot$Symbol[probes2annot],
                        Entrez_gene_ID = annot$Entrez_Gene_ID[probes2annot],
                        moduleColor = moduleColors_signed,
                        geneTraitSignificance,
                        GSPvalue)
#diagnosis??p?l?ŕ??בւ?
modOrder <- order(-abs(cor(MEs, diagnosis, use = "p", method = "spearman")))
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
setwd("E:/IPA_gene_annotation")
IPAannot <- read.csv("GPL14951_IPAannot.csv", header = T, stringsAsFactors = F)
IPAannot <- IPAannot[,-c(2:4)]
geneInfo0 <- merge(geneInfo0, IPAannot, by.x = 1, by.y = 1, all.x = T)
#module color?y??GS?ŕ??בւ?
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.diagnosis))
geneInfo <- geneInfo0[geneOrder, ]
#?A?m?e?[?V?????????y??fibrosis?Ɋւ???GS/MM?̕\???????o??
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
write.csv(geneInfo, paste("geneInfo_diagnosis_signed_", deepSplit, "_merged.csv",sep = "" ), row.names = FALSE)

#???ڂ???trait?idiagnosis?j?ƍł????ւ???"????T?זE2"???W???[???ilightpink3?j?ɂ????āC??GS??MM?̈??`?q(probe)?????肷??
trait <- "diagnosis"
module <- "lightpink3"
column <- match(module, modNames)
moduleGenes <- moduleColors_signed == module

#???ڃ??W???[???ɂ?????GS-MM?v???b?g
pdf(paste("GS_MM_signed_", deepSplit, "_", trait, "_", module, "_merged.pdf",sep = ""), width = 7, height = 7)
par(mfrow = c(1, 1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", trait),
                   main = paste("Module membership vs. gene significance", "(signed_", deepSplit, "_merged)\n", sep = ""),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module, pch = 16)
dev.off()

#???ڂ???trait?idiagnosis?j?ƍł??t???ւ???"?̍זE"???W???[???ifirebrick3?j?ɂ????āC??GS??MM?̈??`?q(probe)?????肷??
trait <- "diagnosis"
module <- "firebrick3"
column <- match(module, modNames)
moduleGenes <- moduleColors_signed == module

#???ڃ??W???[???ɂ?????GS-MM?v???b?g
pdf(paste("GS_MM_signed_", deepSplit, "_", trait, "_", module, "_merged.pdf",sep = ""), width = 7, height = 7)
par(mfrow = c(1, 1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", trait),
                   main = paste("Module membership vs. gene significance", "(signed_", deepSplit, "_merged)\n", sep = ""),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module, pch = 16)
dev.off()

#???ڂ???trait?idiagnosis?j?Ǝ??ɋt???ւ???"?̍זE?H"???W???[???ithistle1?j?ɂ????āC??GS??MM?̈??`?q(probe)?????肷??
trait <- "diagnosis"
module <- "thistle1"
column <- match(module, modNames)
moduleGenes <- moduleColors_signed == module

#???ڃ??W???[???ɂ?????GS-MM?v???b?g
pdf(paste("GS_MM_signed_", deepSplit, "_", trait, "_", module, "_merged.pdf",sep = ""), width = 7, height = 7)
par(mfrow = c(1, 1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", trait),
                   main = paste("Module membership vs. gene significance", "(signed_", deepSplit, "_merged)\n", sep = ""),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module, pch = 16)
dev.off()


#???ڂ???trait?idiagnosis?j?ƍł????ւ???"?̍זE"???W???[???imidnightblue?j?ɂ????āC??GS??MM?̈??`?q(probe)?????肷??
trait <- "diagnosis"
module <- "midnightblue"
column <- match(module, modNames)
moduleGenes <- moduleColors_signed == module

#???ڃ??W???[???ɂ?????GS-MM?v???b?g
pdf(paste("GS_MM_signed_", deepSplit, "_", trait, "_", module, "_merged.pdf",sep = ""), width = 7, height = 7)
par(mfrow = c(1, 1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", trait),
                   main = paste("Module membership vs. gene significance", "(signed_", deepSplit, "_merged)\n", sep = ""),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module, pch = 16)
dev.off()

#???ڂ???trait?idiagnosis?j?ƍł????ւ??郂?W???[???iblue?j?ɂ????āC??GS??MM?̈??`?q(probe)?????肷??
trait <- "diagnosis"
module <- "blue"
column <- match(module, modNames)
moduleGenes <- moduleColors_signed == module

#???ڃ??W???[???ɂ?????GS-MM?v???b?g
pdf(paste("GS_MM_signed_", deepSplit, "_", trait, "_", module, "_merged.pdf",sep = ""), width = 7, height = 7)
par(mfrow = c(1, 1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", trait),
                   main = paste("Module membership vs. gene significance", "(signed_", deepSplit, "_merged)\n", sep = ""),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module, pch = 16)
dev.off()

#???ڂ???trait?idiagnosis?j?Ƌt???ւ????̍זE???W???[???idarkseagreen3?j?ɂ????āC??GS??MM?̈??`?q(probe)?????肷??
trait <- "diagnosis"
module <- "darkseagreen3"
column <- match(module, modNames)
moduleGenes <- moduleColors_signed == module

#???ڃ??W???[???ɂ?????GS-MM?v???b?g
pdf(paste("GS_MM_signed_", deepSplit, "_", trait, "_", module, "_merged.pdf",sep = ""), width = 7, height = 7)
par(mfrow = c(1, 1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", trait),
                   main = paste("Module membership vs. gene significance", "(signed_", deepSplit, "_merged)\n", sep = ""),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module, pch = 16)
dev.off()

#???ڂ???trait?idiagnosis?j?Ƒ??ւ???T?זE???W???[???iCoral3?j?ɂ????āC??GS??MM?̈??`?q(probe)?????肷??
trait <- "diagnosis"
module <- "coral3"
column <- match(module, modNames)
moduleGenes <- moduleColors_signed == module

#???ڃ??W???[???ɂ?????GS-MM?v???b?g
pdf(paste("GS_MM_signed_", deepSplit, "_", trait, "_", module, "_merged.pdf",sep = ""), width = 7, height = 7)
par(mfrow = c(1, 1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", trait),
                   main = paste("Module membership vs. gene significance", "(signed_", deepSplit, "_merged)\n", sep = ""),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module, pch = 16)
dev.off()

#???ڂ???trait?idiagnosis?j?Ƃ??܂葊?ւ??Ȃ??}?N???t?@?[?W???W???[???isteelblue?j?ɂ????āC??GS??MM?̈??`?q(probe)?????肷??
trait <- "diagnosis"
module <- "steelblue"
column <- match(module, modNames)
moduleGenes <- moduleColors_signed == module

#???ڃ??W???[???ɂ?????GS-MM?v???b?g
pdf(paste("GS_MM_signed_", deepSplit, "_", trait, "_", module, "_merged.pdf",sep = ""), width = 7, height = 7)
par(mfrow = c(1, 1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", trait),
                   main = paste("Module membership vs. gene significance", "(signed_", deepSplit, "_merged)\n", sep = ""),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module, pch = 16)
dev.off()





15_series_matrix_signed4merge_eigengene_barplot_dendrogram_heatmap.R





#series matrix???p????WGCNA
#signed_deepSplit4_merged?f?[?^
#???W???[??eigengene?̉???#
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
library(WGCNA)
library(ggplot2)
library(reshape2)
library(stringr)
options(stringsAsFactors = FALSE)
lnames <- load("dataInput.Rdata")
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
annot <- read.csv("GeneAnnotation_GPL14951.csv", header = T, stringsAsFactors = F)
deepSplit <- 4

##signed_merge????
#?t?H???_?ړ?
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("networkConstruction_StepByStep_signed_", deepSplit, "_merged.Rdata", sep = ""))

#???`?q?i?v???[?u?j???ƃT???v?????????`
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

# ME???Čv?Z
MEs <- moduleEigengenes(datExpr, moduleColors_signed)$eigengenes
row.names(MEs) <- row.names(datExpr)
module_color <- str_sub(colnames(MEs), start = 3) #?擪2????ME???폜
len <- length(module_color)

# ?t?H???_?쐬?ƈړ?
ifelse(!dir.exists(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, "/heatmap", sep = "")), 
       dir.create(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, "/heatmap", sep = "")), FALSE)
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, "/heatmap", sep = ""))

# ?e???W???[????heatmap??eigengene???`??
for(i in 1:len){
  which.module <- module_color[i]
  sizeGrWindow(8, 7)
  pdf(paste(which.module, "_heatmap_eigengene.pdf", sep = ""), width = 8, height = 7)
  ME <- MEs[, paste("ME", which.module, sep = "")]
  par(mfrow = c(2,1), mar = c(0.3, 5.5, 3, 2))
  plotMat(t(scale(datExpr[, moduleColors_signed == which.module])),
          nrgcols = 30, rlabels = F, rcols = which.module,
          main = which.module, cex.main = 2)
  par(mar = c(5, 4.2, 0, 0.7))
  barplot(ME, col = which.module, main = "", cex.main = 2,
          ylab = "eigengene expression", xlab = "array sample")
  dev.off()
}


#ME???Čv?Z
MEs <- moduleEigengenes(datExpr, moduleColors_signed)$eigengenes
Diagnosis <- factor(datTraits$Diagnosis, levels = c(0,1,2))
Fibrosis <- factor(datTraits$Fibrosis, levels = c("NA", "0", "1", "2", "3", "4"))
MEs <- data.frame(Sample_ID = rownames(datExpr), MEs, Diagnosis, Fibrosis)
Sample_ID <- row.names(datExpr)

#ME??long?^?f?[?^?ɕϊ?
MEs <- melt(MEs, id.vars = c("Sample_ID", "Diagnosis", "Fibrosis"), variable.name = "Module", value.name = "Eigengene")
color <- str_sub(MEs$Module, start = 3) #?擪2????ME???폜
color_palette <- unique(color)

#?e???ҁC?e???W???[???ɂ?????ME??diagnosis???ɖ_?O???t?Ő}??
a <- ggplot(MEs, aes(x=Sample_ID, y=Eigengene)) +
  theme_light() +  
  geom_bar(stat = "identity", fill = color, colour = "grey", size = 0.1) +
  facet_wrap(~Module, ncol = 4, nrow = 10, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 4)) +
  theme(strip.text = element_text(size = rel(0.7), colour = "black")) +
  scale_x_discrete(limits = Sample_ID)
pdf("ME_diagnosis_barplot.pdf", width = 10, height = 16)
plot(a)
dev.off()

#?e???ҁC?e???W???[???ɂ?????ME??diagnosis???ɔ??Ђ??}?Ő}??
b <- ggplot(MEs, aes(x = Diagnosis, y = Eigengene, fill = Module)) +
  theme_light() +  
  geom_boxplot() +
  facet_wrap(~Module, ncol = 4, nrow = 10, scales = "free_y") +
  scale_x_discrete(breaks = c(0, 1, 2), labels = c("HC", "SS", "NASH")) +
  scale_fill_manual(values = color_palette) +
  guides(fill = FALSE)
pdf("ME_diagnosis_boxplot.pdf", width = 10, height = 16)
plot(b)
dev.off()

#Sample_ID??fibrosis???ɕ??בւ?
Sample_ID_fibrosis <- c("F0_HC_1",    "F0_HC_2",    "F0_HC_3",    "F0_HC_5",    "F0_HC_6",    "F0_HC_7",    "F0_HC_8",    "F0_HC_9",    "F0_HC_10",
                        "F0_SS_1",    "F0_SS_2",    "F0_SS_3",    "F0_SS_4",    "F0_SS_5",    "F0_SS_6",    "F0_SS_7",    "F0_SS_8",    "F0_SS_9",    
                        "F0_SS_10",   "F0_SS_11",   "F0_SS_12",   "F0_SS_13",   "F0_SS_14",   "F0_SS_16",   "F0_SS_17",
                        "F0_NASH_1",  "F0_NASH_2",  "F0_NASH_3", "F0_NASH_4",
                        "F1_HC_11",   "F1_HC_12",   "F1_HC_13",   "F1_HC_14",   "F1_HC_15",
                        "F1_SS_18",   "F1_SS_19",   "F1_SS_20",
                        "F1_NASH_5",  "F1_NASH_6",  "F1_NASH_7",  "F1_NASH_8",  "F1_NASH_9",
                        "F2_NASH_10", "F2_NASH_11",
                        "F3_NASH_12", "F3_NASH_13", "F3_NASH_14", "F3_NASH_15", "F4_NASH_16", "F4_NASH_17", "F4_NASH_18", "F4_NASH_19",
                        "FN_HC_17",   "FN_HC_18",   "FN_HC_19",   "FN_HC_20",   "FN_HC_21",   "FN_HC_22",   "FN_HC_23",   "FN_HC_24")

#?e???ҁC?e???W???[???ɂ?????ME??fibrosis???ɖ_?O???t?Ő}??
a <- ggplot(MEs, aes(x=Sample_ID, y=Eigengene)) +
  theme_light() +  
  geom_bar(stat = "identity", fill = color, colour = "grey", size = 0.1) +
  facet_wrap(~Module, ncol = 4, nrow = 10, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5, size = 4)) +
  theme(strip.text = element_text(size = rel(0.7), colour = "black")) +
  scale_x_discrete(limits = Sample_ID_fibrosis)
pdf("ME_fibrosis_barplot.pdf", width = 10, height = 16)
plot(a)
dev.off()

#fibrosis??NA?̍s???폜
MEs <- subset(MEs, MEs$Fibrosis != "NA")

#?e???ҁC?e???W???[???ɂ?????ME??fibrosis???ɔ??Ђ??}?Ő}??
b <- ggplot(MEs, aes(x = Fibrosis, y = Eigengene, fill = Module)) +
  theme_light() +  
  geom_boxplot() +
  facet_wrap(~Module, ncol = 4, nrow = 10, scales = "free_y") +
  scale_fill_manual(values = color_palette) +
  guides(fill = FALSE)
pdf("ME_fibrosis_boxplot.pdf", width = 10, height = 16)
plot(b)
dev.off()


#ME???Čv?Z
MEs <- moduleEigengenes(datExpr, moduleColors_signed)$eigengenes

#trait???????o??
fibrosis <- as.data.frame(datTraits$Fibrosis)
names(fibrosis) <- "fibrosis"
diagnosis <- as.data.frame(datTraits$Diagnosis)
names(diagnosis) <- "diagnosis"
NAS <- as.data.frame(datTraits$NAS)
names(NAS) <- "NAS"
BMI <- as.data.frame(datTraits$BMI)
names(BMI) <- "BMI"

#ME??trait?Ԃ̊֌W???f???h???O?????Ő}??
MET <- orderMEs(cbind(MEs, fibrosis))
pdf("ME_fibrosis_dendrogram.pdf", width = 10, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
dev.off()

MET <- orderMEs(cbind(MEs, diagnosis))
pdf("ME_diagnosis_dendrogram.pdf", width = 10, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
dev.off()

MET <- orderMEs(cbind(MEs, NAS))
pdf("ME_NAS_dendrogram.pdf", width = 10, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
dev.off()

MET <- orderMEs(cbind(MEs, BMI))
pdf("ME_BMI_dendrogram.pdf", width = 10, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
dev.off()





16_series_matrix_signed4merge_networkdata_cytoscape.R





#series matrix???p????WGCNA
#signed_deepSplit4_merged?f?[?^
#???W???[???̃l?b?g???[?N?f?[?^??Cytoscape??#
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis")
library(WGCNA)
library(stringr)
options(stringsAsFactors = FALSE)
lnames <- load("dataInput.Rdata")
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/data")
annot <- read.csv("GeneAnnotation_GPL14951.csv", header = T, stringsAsFactors = F)
deepSplit <- 4

##signed_merge????
#?t?H???_?ړ?
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, sep = ""))
#?ۑ??????l?b?g???[?N?f?[?^???ǂݍ???
lnames <- load(paste("networkConstruction_StepByStep_signed_", deepSplit, "_merged.Rdata", sep = ""))

#topological overlap???Čv?Z
TOM <- TOMsimilarityFromExpr(datExpr, networkType = "signed", power = 20, TOMType = "signed")

#???Ȃ݂ɁC???L?̃X?N???v?g?ł?????TOM??????????
#adjacency?i?אڍs???l?j??topological overlap?iTOM?j?ɕϊ?
#softPower <- 20
#adjacency <- adjacency(datExpr, power = softPower, type = "signed")
#TOM <- TOMsimilarity(adjacency, TOMType = "signed")

#ME???Čv?Z
MEs <- moduleEigengenes(datExpr, moduleColors_signed)$eigengenes
row.names(MEs) <- row.names(datExpr)
module_color <- str_sub(colnames(MEs), start = 3) #?擪2????ME???폜
len <- length(module_color)

#?t?H???_?쐬?ƈړ?
ifelse(!dir.exists(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, "/Cytoscape", sep = "")), 
       dir.create(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, "/Cytoscape", sep = "")), FALSE)
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, "/Cytoscape", sep = ""))

#?e???W???[?????Ƀl?b?g???[?N?f?[?^???Z?o

for(i in 1:len){
#???W???[???̑I??
modules <- module_color[i]

#???W???[???̃v???[?u???ƍ?
probes <- colnames(datExpr)
inModule <- is.finite(match(moduleColors_signed, modules))
modProbes <- probes[inModule]
modGenes <- annot$Symbol[match(modProbes, annot$ID)]

#?Ή?????TOM???I??
modTOM <- TOM[inModule, inModule]
dimnames(modTOM) <- list(modProbes, modProbes)

#Cytoscape?p?Ƀl?b?g???[?N??edge??node???o??
cyt <- exportNetworkToCytoscape(modTOM,
                                edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse = "-"), ".txt", sep = ""),
                                nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse = "-"), ".txt", sep = ""),
                                weighted = TRUE,
                                threshold = 0.02,
                                nodeNames = modProbes,
                                altNodeNames = modGenes,
                                nodeAttr = moduleColors_signed[inModule])
}


#???ڂ??郂?W???[???̃l?b?g???[?N?f?[?^???Z?o???悤?Ƃ??????C?f?[?^?T?C?Y???傫??????
#???W???[???̑I??
#modules <- c("tan",                                                       #HSC???W???[???Cfibrosis???W???[??
#             "lightpink3",                                                #????Tcell2???W???[???Cdiagnosis???W???[??
#             "blue",                                                      #diagnosis???W???[??
#             "darkseagreen3", "firebrick3",                               #hepatocyte???W???[??, diagnosis???W???[??
#             "thistle1","cyan", "darkviolet", "midnightblue", "thistle2", #hepatocyte???W???[??
#             "coral3",                                                    #Tcell???W???[??
#             "steelblue"                                                  #M?Ӄ??W???[??
#             )


#?e???W???[???ɂ?????20 top hub genes?????p?????悤?ɂ??Ă݂悤
#WGCNA?`???[?g???A??I-6?ł?top hub genes?̑I????softConnectivity()?ŎZ?o????intramodular connectivity???g?p???Ă???
#III-7?ł?intramodularConnectivity()?ŎZ?o??????kWithin??intramodular connectivity?ł????Ə??????Ă???
#softConnectivity()?ŎZ?o??????IMConn??kWithin???v???b?g?????ƁC?????͑??ւ??????̂́C?ꕔ???ւ????O???????̂?????
#?????Cmodule eigengene??gene expression?̑??ւŎZ?o??????module membership (MM)?͊e???W???[???ɂ?????kWithin?Ƒ??ւ??邱?Ƃ??m?F???Ă???
#?????????āC20 top hub genes?̑I???ɂ?kWithin???g?p????

#IMConn_all <- softConnectivity(datExpr, type = "signed", power = 20)
#write.csv(IMConn_all, "IMConn.csv")

#kWithin?f?[?^???ǂݍ???
setwd("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_4")
kWithin <- read.csv("Probe_kWitnin_color_signed_4_merged.csv", header = T, stringsAsFactors = F)

#ME???Čv?Z
MEs <- moduleEigengenes(datExpr, moduleColors_signed)$eigengenes
row.names(MEs) <- row.names(datExpr)
module_color <- str_sub(colnames(MEs), start = 3) #?擪2????ME???폜
len <- length(module_color)
all_topProbes <- NULL

#?e???W???[??????20 top hub genes?????o
for(i in 1:len){
  #???W???[???̑I??
  modules <- module_color[i]
  
  #???W???[???̃v???[?u???I??
  probes <- colnames(datExpr)
  inModule <- is.finite(match(moduleColors_signed, modules))
  modProbes <- probes[inModule]
  
  #20 top hub genes???I??
  nTop <- 20
  IMConn <- kWithin$kWithin[inModule]
  top <- (rank(-IMConn) <= nTop)
  topProbes <- modProbes[top]
  
  #?e???W???[????20 top hub genes?̃??X?g??????
  all_topProbes <- c(all_topProbes, topProbes)
}  

#???ׂẴ??W???[????20 top hub genes
intopProbes <- is.finite(match(probes, all_topProbes))?@?@?@   #datExpr?̕??я???20 top hub genes??TRUE?C?????ȊO??FALSE?Ƃ????x?N?g??
summary(intopProbes)
all_topProbes <- probes[intopProbes]?@?@?@?@?@?@?@?@?@ ?@?@?@?@#20 top hub genes??datExpr?̕??я??ɔ????o??
all_topGenes <- annot$Symbol[match(all_topProbes, annot$ID)]   #Probe ID??gene symbol?ɕϊ?

#?Ή?????TOM???I??
topTOM <- TOM[intopProbes, intopProbes]                        #20 top hub genes??TOM??datExpr?̕??я??ɔ????o??
dimnames(topTOM) <- list(all_topProbes, all_topProbes)
dim(topTOM)
topTOM[1:10, 1:10]

#?t?H???_?쐬?ƈړ?
ifelse(!dir.exists(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, "/Cytoscape", sep = "")), 
       dir.create(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, "/Cytoscape", sep = "")), FALSE)
setwd(paste("//cserv48.pharmacpri.jti.co.jp/smb1/01_DiseaseArea/01_Metabolic_Diseases/003_NASH/190201_NAFLD_NASH patient WGCNA_GSE89632/analysis/signed_", deepSplit, "/Cytoscape", sep = ""))

#Cytoscape?p?Ƀl?b?g???[?N??edge??node???o??
cyt <- exportNetworkToCytoscape(topTOM,
                                edgeFile = paste("CytoscapeInput-edges-top20.txt", sep = ""),
                                nodeFile = paste("CytoscapeInput-nodes-top20.txt", sep = ""),
                                weighted = TRUE,
                                threshold = 0.02,
                                nodeNames = all_topProbes,
                                altNodeNames = all_topGenes,
                                nodeAttr = moduleColors_signed[intopProbes])