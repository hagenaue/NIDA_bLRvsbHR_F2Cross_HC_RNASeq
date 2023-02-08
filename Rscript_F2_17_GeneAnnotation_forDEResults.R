#F2 HC RNA-Seq Dataset
#17_Adding gene annotation to the differential expression results
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-04, updated later for a few figures for the paper.



#annotating the Limma results:
library(plyr)

F2_Limma_Results_EnsemblVsGeneSymbol<-join(F2_Limma_Results, EnsemblVsGeneSymbol, by="ENSEMBL", type="left", match="all")
str(F2_Limma_Results_EnsemblVsGeneSymbol)
#''data.frame':	14024 obs. of  25 variables:
#'
write.csv(F2_Limma_Results_EnsemblVsGeneSymbol, "F2_Limma_Results_EnsemblVsGeneSymbol.csv")

F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped<-F2_Limma_Results_EnsemblVsGeneSymbol[duplicated(F2_Limma_Results_EnsemblVsGeneSymbol$ENSEMBL)==FALSE, ]
str(F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped)
#'data.frame':	13814 obs. of  25 variables:
#'There we go.

write.csv(F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped, "F2_Limma_Results_TotalLocoScore_M4_EnsemblVsGeneSymbol_NoMultimapped.csv")


#######################2/28/2022  Elaine  Needed to annotate:  F2_Limma_Results_TotalLocoScore_SimplerModel_20220110.txt

#annotating the Limma results:
library(plyr)


F2_Limma_Results_TotalLocoScore_M9<-read.delim("F2_Limma_Results_TotalLocoScore_SimplerModel_20220110.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F2_Limma_Results_TotalLocoScore_M9)
colnames(F2_Limma_Results_TotalLocoScore_M9)[1]<-"ENSEMBL"
# 'data.frame':	14056 obs. of  40 variables:
#   $ ENSEMBL                        : chr  "ENSRNOG00000014303" "ENSRNOG00000014330" "ENSRNOG00000049505" "ENSRNOG00000014916" ...
# $ AveExpr                        : num  7.76 8.34 3.07 5.14 3.81 ...


F2_Limma_Results_EnsemblVsGeneSymbol_Loco_M9<-join(F2_Limma_Results_TotalLocoScore_M9, EnsemblVsGeneSymbol, by="ENSEMBL", type="left", match="all")
str(F2_Limma_Results_EnsemblVsGeneSymbol_Loco_M9)
#' data.frame':	14258 obs. of  41 variables:
#'  $ ENSEMBL                        : chr  "ENSRNOG00000014303" "ENSRNOG00000014330" "ENSRNOG00000049505" "ENSRNOG00000014916" ...
#'  $ AveExpr                        : num  7.76 8.34 3.07 5.14 3.81 ...
#'  $ Coef..Intercept.               : num  6.02 11.01 5.87 3.61 6.18 ...

write.csv(F2_Limma_Results_EnsemblVsGeneSymbol_Loco_M9, "F2_Limma_Results_EnsemblVsGeneSymbol_TotalLocoScore_M9.csv")

F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_Loco_M9<-F2_Limma_Results_EnsemblVsGeneSymbol_Loco_M9[duplicated(F2_Limma_Results_EnsemblVsGeneSymbol_Loco_M9$ENSEMBL)==FALSE, ]
str(F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_Loco_M9)
# data.frame':	14056 obs. of  41 variables:
#  $ ENSEMBL                        : chr  "ENSRNOG00000014303" "ENSRNOG00000014330" "ENSRNOG00000049505" "ENSRNOG00000014916" ...
#  $ AveExpr                        : num  7.76 8.34 3.07 5.14 3.81 ...
#  $ Coef..Intercept.               : num  6.02 11.01 5.87 3.61 6.18 ...

write.csv(F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_Loco_M9, "F2_Limma_Results_TotalLocoScore_M9_EnsemblVsGeneSymbol_NoMultimapped.csv")


