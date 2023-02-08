#F0 HC RNA-Seq Dataset
#11_Annotating DE Results
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-10, updated later for a few figures for the paper.


#############################

#annotating the Limma results:


F0_Limma_Results_EnsemblVsGeneSymbol<-join(F0_Limma_Results, EnsemblVsGeneSymbol, by="ENSEMBL", type="left", match="all")
str(F0_Limma_Results_EnsemblVsGeneSymbol)
#''data.frame':	14024 obs. of  25 variables:
#'
write.csv(F0_Limma_Results_EnsemblVsGeneSymbol, "F0_Limma_Results_EnsemblVsGeneSymbol.csv")

F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped<-F0_Limma_Results_EnsemblVsGeneSymbol[duplicated(F0_Limma_Results_EnsemblVsGeneSymbol$ENSEMBL)==FALSE, ]
str(F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped)
#'data.frame':	13814 obs. of  25 variables:
#'There we go.

write.csv(F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped, "F0_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped.csv")



F0_Limma_Results_wSexInteractions<-read.delim("F0_Limma_Results_LineageSex_Interaction_RibosomePerc.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F0_Limma_Results_wSexInteractions)

#'data.frame':	13814 obs. of  28 variables:

colnames(F0_Limma_Results_wSexInteractions)[1]<-"ENSEMBL"

F0_Limma_Results_wSexInteractions_EnsemblVsGeneSymbol<-join(F0_Limma_Results_wSexInteractions, EnsemblVsGeneSymbol, by="ENSEMBL", type="left", match="all")
str(F0_Limma_Results_wSexInteractions_EnsemblVsGeneSymbol)
#''data.frame':	14024 obs. of  29 variables:
#'
write.csv(F0_Limma_Results_wSexInteractions_EnsemblVsGeneSymbol, "F0_Limma_Results_wSexInteractions_EnsemblVsGeneSymbol.csv")

F0_Limma_Results_wSexInteractions_EnsemblVsGeneSymbol_NoMultimapped<-F0_Limma_Results_wSexInteractions_EnsemblVsGeneSymbol[duplicated(F0_Limma_Results_wSexInteractions_EnsemblVsGeneSymbol$ENSEMBL)==FALSE, ]
str(F0_Limma_Results_wSexInteractions_EnsemblVsGeneSymbol_NoMultimapped)
#'data.frame':	13814 obs. of  29 variables:
#'There we go.

write.csv(F0_Limma_Results_wSexInteractions_EnsemblVsGeneSymbol_NoMultimapped, "F0_Limma_Results_wSexInteractions_EnsemblVsGeneSymbol_NoMultimapped.csv")



#########################################