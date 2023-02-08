#F2 HC RNA-Seq Dataset
#09_Renormalization and filtering of the RNA-Seq data following sample QC
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-04, updated later for a few figures for the paper.

#Reprocessing the RNA-Seq data again now that the mislabeled samples have been removed
#Note - we are overwriting previous objects/files:

dge<-DGEList(counts=F2_RawHitCount_ExpressionData)
F2_RawHitCount_dge_TMM<-calcNormFactors(dge, method="TMM")

str(F2_RawHitCount_dge_TMM)
# Formal class 'DGEList' [package "edgeR"] with 1 slot
# ..@ .Data:List of 2
# .. ..$ : int [1:32883, 1:245] 0 0 0 0 0 0 2 0 0 0 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:32883] "ENSRNOG00000046319" "ENSRNOG00000047964" "ENSRNOG00000050370" "ENSRNOG00000032365" ...
# .. .. .. ..$ : chr [1:245] "SL469967" "SL469968" "SL469969" "SL469970" ...
# .. ..$ :'data.frame':	245 obs. of  3 variables:
#   .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1 1 1 1 1 1 1 ...
# .. .. ..$ lib.size    : num [1:245] 21826046 19217774 20269455 26635890 21240639 ...
# .. .. ..$ norm.factors: num [1:245] 0.951 0.952 0.964 1.016 0.969 ...
# ..$ names: chr [1:2] "counts" "samples"

F2_cpm<-cpm(F2_RawHitCount_dge_TMM)
str(F2_cpm)
#num [1:32883, 1:245] 0 0 0 0 0 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:32883] "ENSRNOG00000046319" "ENSRNOG00000047964" "ENSRNOG00000050370" "ENSRNOG00000032365" ...
# ..$ : chr [1:245] "SL469967" "SL469968" "SL469969" "SL469970" ...

#Log2 transformed - an approximate way of getting rid of heteroskedasticity in the data and the traditional way of displaying RNA-Seq data in graphs
F2_lcpm<-cpm(F2_RawHitCount_dge_TMM, log=TRUE)

write.csv(F2_cpm,"F2_cpm.csv")
write.csv(F2_lcpm, "F2_lcpm.csv")

#Filtering out rows of data that are 0 cpm for at least a quarter of the animals (e.g. all male bHRs and no one else)
#This is Megan erring on the side of including too much instead of including too little

sum(rowSums(F2_cpm>1)>=1)
#[1] 15510

#For the F0s we filtered to make sure that genes were expressed in at least 1/4 of the subjects because we had 4 subgroups of interest (bHR/bLR vs. M/F)
#For consistency sake, let's use the same filter here (note: you could make other arguments for thresholds - so if some low level expressed genes of interest aren't showing up in the data, we can go back and change this)
247/4
#[1] 61.75

keep.exprs<-rowSums(F2_cpm>1)>=62
F2_RawHitCount_ExpressionData_noLowHits<-F2_RawHitCount_ExpressionData[keep.exprs,]
str(F2_RawHitCount_ExpressionData_noLowHits)
#  int [1:14056, 1:245] 6484 3991 120 894 155 344 1235 1167 186 2232 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:14056] "ENSRNOG00000014303" "ENSRNOG00000014330" "ENSRNOG00000049505" "ENSRNOG00000014916" ...
# ..$ : chr [1:245] "SL469967" "SL469968" "SL469969" "SL469970" ...

###################

#TMM normalization:

dge_noLowHits<-DGEList(counts=F2_RawHitCount_ExpressionData_noLowHits)
F2_RawHitCount_dge_noLowHits_TMM<-calcNormFactors(dge_noLowHits, method="TMM")

#create a matrix of log2 CPM (counts per million) values:

F2_cpm_noLowHits<-cpm(F2_RawHitCount_dge_noLowHits_TMM)

F2_lcpm_noLowHits<-cpm(F2_RawHitCount_dge_noLowHits_TMM, log=TRUE)


str(F2_lcpm_noLowHits)
# num [1:14056, 1:245] 8.24 7.54 2.51 5.38 2.87 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:14056] "ENSRNOG00000014303" "ENSRNOG00000014330" "ENSRNOG00000049505" "ENSRNOG00000014916" ...
# ..$ : chr [1:245] "SL469967" "SL469968" "SL469969" "SL469970" ...

#Making the log2 cpm values into a data.frame so we can join them with annotation:
F2_lcpm_noLowHits_DF<-as.data.frame(F2_lcpm_noLowHits, header=TRUE, stringAsFactors=FALSE)
#Make sure this step included the ENSEMBL ids as a column. If they are still row.names, we will need to make a column in the data.frame named ENSEMBL for those row.names
str(F2_lcpm_noLowHits_DF)

#'data.frame':	14056 obs. of  245 variables:
#' $ SL469967: num  8.24 7.54 2.51 5.38 2.87 ...
#' $ SL469968: num  8.19 7.69 2.34 5.47 3.03 ...
#' $ SL469969: num  8.26 7.46 2.39 5.48 2.96 ...
#' $ SL469970: num  8.08 7.6 2.45 5.6 3.43 ...
#' $ SL469971: num  8.29 7.56 2.16 5.49 3.03 ...

