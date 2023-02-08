#F0 HC RNA-Seq Dataset
#05_Gene Annotation
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-10, updated later for a few figures for the paper.


#Adding annotation:

F0_lcpm_noLowHits_DF$ENSEMBL<-row.names(F0_lcpm_noLowHits)


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Rn.eg.db")


library(org.Rn.eg.db)


#Outputting gene symbols for all ENSEMBL genes:

#grabbing all of the ensembl ids in the annotation package:
uniKeys<-keys(org.Rn.eg.db, keytype="ENSEMBL")
str(uniKeys)
#chr [1:20998] "ENSRNOG00000017701" "ENSRNOG00000028896" "ENSRNOG00000032908" "ENSRNOG00000009845" "ENSRNOG00000016924" ...
#This is a vector with all ensembl genes in the org.Rn.eg.db


##################Do not use quotations from the chat in Zoom--Do not copy and paste quotations from Zoom to R

uniKeys[1]
# [1] "ENSRNOG00000017701"

sum(uniKeys == "ENSRNOG00000052237")

sum(uniKeys == "ENSRNOG00000017701")

###########

#We would like the symbol for all of those ensembl ids
cols<-c("SYMBOL", "CHRLOC", "CHRLOCEND") 

#grabbing symbol for every ensembl id in the package:
EnsemblVsGeneSymbol<-AnnotationDbi::select(org.Rn.eg.db, keys=uniKeys, columns=cols, keytype="ENSEMBL")  
dim(EnsemblVsGeneSymbol)
#[1] 22855     5


head(EnsemblVsGeneSymbol)
# ENSEMBL SYMBOL     CHRLOC CHRLOCCHR  CHRLOCEND
# 1 ENSRNOG00000017701   Asip  150492009         3  150579870
# 2 ENSRNOG00000028896    A2m  154423167         4  154473038
# 3 ENSRNOG00000028896    A2m  154309425         4  154359138
# 4 ENSRNOG00000032908 Acaa1a  128027879         8  128036471
# 5 ENSRNOG00000009845  Acadm -260124417         2 -260148589
# 6 ENSRNOG00000016924   Acly  -88392247        10  -88442845

str(EnsemblVsGeneSymbol)
# 'data.frame':	22855 obs. of  5 variables:
#   $ ENSEMBL  : chr  "ENSRNOG00000017701" "ENSRNOG00000028896" "ENSRNOG00000028896" "ENSRNOG00000032908" ...
# $ SYMBOL   : chr  "Asip" "A2m" "A2m" "Acaa1a" ...
# $ CHRLOC   : int  150492009 154423167 154309425 128027879 -260124417 -88392247 -49836064 -49836801 -99120378 -49836685 ...
# $ CHRLOCCHR: chr  "3" "4" "4" "8" ...
# $ CHRLOCEND: int  150579870 154473038 154359138 128036471 -260148589 -88442845 -49851714 -49851648 -99121144 -49851660 ...

write.csv(EnsemblVsGeneSymbol, "EnsemblVsGeneSymbol_wCHRLOC.csv")


# Annotating our data:

install.packages('plyr')
library(plyr)


str(F0_lcpm_noLowHits_DF)
#'data.frame':	13786 obs. of  24 variables:

F0_Annotated_lcpm_NoLowHits<-join(F0_lcpm_noLowHits_DF, EnsemblVsGeneSymbol, by='ENSEMBL', type='left', match='all')
str(F0_Annotated_lcpm_NoLowHits)
# 
# 
#'data.frame':	13996 obs. of  25 variables:
# $ SL483888: num  6.69 9.12 3.69 4.56 4.34 ...
# $ SL483889: num  6.55 8.83 3.48 4.63 4.05 ...
# $ SL483890: num  6.72 8.54 3.24 4.59 3.94 ...
# $ SL483891: num  6.99 8.71 3.07 4.27 3.84 ...
# $ SL483892: num  7.03 8.64 3.22 4.32 3.7 ...
# $ SL483893: num  6.87 8.52 3.25 4.66 3.98 ...
# $ SL483894: num  6.68 9.36 3.65 4.37 4.3 ...
# $ SL483895: num  6.41 9.59 4.11 4.28 4.45 ...
# $ SL483896: num  6.31 9.48 3.9 3.6 4.42 ...
# $ SL483898: num  7.25 8.31 2.94 4.44 3.73 ...
# $ SL483899: num  6.47 9.55 4.05 4.23 4.33 ...
# $ SL483900: num  6.21 9.49 4.27 3.59 4.54 ...
# $ SL483901: num  6.76 9.92 4.4 5.14 4.68 ...
# $ SL483902: num  6.51 9.48 3.76 3.84 4.32 ...
# $ SL483903: num  5.72 10.48 4.72 4.35 4.88 ...
# $ SL483904: num  6.62 9.5 4.05 4.63 4.43 ...
# $ SL483905: num  6.05 9.97 4.05 3.82 4.48 ...
# $ SL483906: num  6.88 9.06 3.84 5.25 4.06 ...
# $ SL483907: num  6.71 9.67 4.09 4.71 4.56 ...
# $ SL483908: num  5.81 10.06 4.43 3.88 4.74 ...
# $ SL483909: num  7.03 8.75 3.67 4.04 3.91 ...
# $ SL483910: num  5.93 9.85 4.32 3.94 4.77 ...
# $ SL483911: num  6.58 9.63 4.29 4.82 4.6 ...
# $ ENSEMBL : chr  "ENSRNOG00000014303" "ENSRNOG00000014330" "ENSRNOG00000049505" "ENSRNOG00000014916" ...
# $ SYMBOL  : chr  "Lrp11" "Pcmt1" "Nup43" "Lats1" ...

#IMPORTANT: Some of the ensembl genes are mapping to more than one gene symbol
#This means that adding annotation made it so that some rows of data are present twice (or more) in the dataset
#For running analyses or fGSEA, we would need to do something about that.
#For basic plotting, not so much.

write.csv(F0_Annotated_lcpm_NoLowHits, "F0_Annotated_lcpm_NoLowHits.csv")

