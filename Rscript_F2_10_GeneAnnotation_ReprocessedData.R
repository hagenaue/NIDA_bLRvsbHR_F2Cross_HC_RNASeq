#F2 HC RNA-Seq Dataset
#10_Adding gene annotation to reprocessed RNA-Seq data
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-04, updated later for a few figures for the paper.


##Adding annotation to our newly reprocessed data (without the mislabeled samples:)

F2_lcpm_noLowHits_DF$ENSEMBL<-row.names(F2_lcpm_noLowHits)

# (building on our previous annotation code):

str(F2_lcpm_noLowHits_DF)
#'data.frame':	14056 obs. of  246 variables:

F2_Annotated_lcpm_NoLowHits<-join(F2_lcpm_noLowHits_DF, EnsemblVsGeneSymbol, by='ENSEMBL', type='left', match='all')
str(F2_Annotated_lcpm_NoLowHits)
# 
# 'data.frame':	14258 obs. of  247 variables:
# $ SL469967: num  8.24 7.54 2.51 5.38 2.87 ...
# $ SL469968: num  8.19 7.69 2.34 5.47 3.03 ...
# $ SL469969: num  8.26 7.46 2.39 5.48 2.96 ...
# $ SL469970: num  8.08 7.6 2.45 5.6 3.43 ...
# $ SL469971: num  8.29 7.56 2.16 5.49 3.03 ...

#IMPORTANT: Some of the ensembl genes are mapping to more than one gene symbol
#This means that adding annotation made it so that some rows of data are present twice (or more) in the dataset
#For running analyses or fGSEA, we would need to do something about that.
#For basic plotting, not so much.

write.csv(F2_Annotated_lcpm_NoLowHits, "F2_Annotated_lcpm_NoLowHits.csv")

