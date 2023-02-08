#F2 HC RNA-Seq Dataset
#07_GeneAnnotation
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-04, updated later for a few figures for the paper.


#Adding annotation:

F2_lcpm_noLowHits_DF$ENSEMBL<-row.names(F2_lcpm_noLowHits)

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

#We would like the symbol for all of those ensembl ids
cols<-c("SYMBOL") 

#grabbing symbol for every ensembl id in the package:
EnsemblVsGeneSymbol<-AnnotationDbi::select(org.Rn.eg.db, keys=uniKeys, columns=cols, keytype="ENSEMBL")  
dim(EnsemblVsGeneSymbol)
#[1] 21556     2
head(EnsemblVsGeneSymbol)
#   ENSEMBL             SYMBOL
# 1 ENSRNOG00000017701   Asip
# 2 ENSRNOG00000028896    A2m
# 3 ENSRNOG00000032908 Acaa1a
# 4 ENSRNOG00000009845  Acadm
# 5 ENSRNOG00000016924   Acly
# 6 ENSRNOG00000005260   Acp1

str(EnsemblVsGeneSymbol)

# 'data.frame':	21556 obs. of  2 variables:
# $ ENSEMBL: chr  "ENSRNOG00000017701" "ENSRNOG00000028896" "ENSRNOG00000032908" "ENSRNOG00000009845" ...
# $ SYMBOL : chr  "Asip" "A2m" "Acaa1a" "Acadm" ...

write.csv(EnsemblVsGeneSymbol, "EnsemblVsGeneSymbol.csv" )

#Then annotate it (building on our previous annotation code):

install.packages('plyr')
library(plyr)


str(F2_lcpm_noLowHits_DF)
#'data.frame':	14061 obs. of  248 variables:

F2_Annotated_lcpm_NoLowHits<-join(F2_lcpm_noLowHits_DF, EnsemblVsGeneSymbol, by='ENSEMBL', type='left', match='all')
str(F2_Annotated_lcpm_NoLowHits)
# 
# data.frame':	14263 obs. of  249 variables:
# $ SL469967: num  8.21 7.51 2.48 5.36 2.84 ...
# $ SL469968: num  8.17 7.67 2.32 5.45 3.01 ...
# $ SL469969: num  8.26 7.45 2.38 5.47 2.95 ...
# $ SL469970: num  8.05 7.57 2.43 5.57 3.41 ...
# $ SL469971: num  8.28 7.55 2.15 5.48 3.01 .....


#IMPORTANT: Some of the ensembl genes are mapping to more than one gene symbol
#This means that adding annotation made it so that some rows of data are present twice (or more) in the dataset
#For running analyses or fGSEA, we would need to do something about that.
#For basic plotting, not so much.

write.csv(F2_Annotated_lcpm_NoLowHits, "F2_Annotated_lcpm_NoLowHits.csv")


########################