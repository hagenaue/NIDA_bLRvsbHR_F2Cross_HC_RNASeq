#F0 HC RNA-Seq Dataset
#04_TMM Normalization
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-10, updated later for a few figures for the paper.

png("F0_boxplot_RawExpressionData_BySubject.png")
boxplot(F0_RawHitCount_ExpressionData, las=2)
title(main="Raw Reads Data by Subject", ylab="Counts")
dev.off()

#Adding in TMM normalization because PC1 (from Ahub) is strongly related to relative mRNA production (RNA concentration, mRNA percentage, %rRNA) 
#TMM normalization scales to corrrect for estimated relative RNA production levels at the point of calculating cpm

dge<-DGEList(counts=F0_RawHitCount_ExpressionData)
F0_RawHitCount_dge_TMM<-calcNormFactors(dge, method="TMM")

str(F0_RawHitCount_dge_TMM)
#Formal class 'DGEList' [package "edgeR"] with 1 slot
# ..@ .Data:List of 2
# .. ..$ : int [1:32883, 1:23] 0 0 0 0 0 0 0 0 0 0 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:32883] "ENSRNOG00000046319" "ENSRNOG00000047964" "ENSRNOG00000050370" "ENSRNOG00000032365" ...
# .. .. .. ..$ : chr [1:23] "SL483888" "SL483889" "SL483890" "SL483891" ...
# .. ..$ :'data.frame':	23 obs. of  3 variables:
#   .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1 1 1 1 1 1 1 ...
# .. .. ..$ lib.size    : num [1:23] 16824899 18993612 25568698 23034700 24765660 ...
# .. .. ..$ norm.factors: num [1:23] 1 1.07 1.13 1.28 1.3 ...
# ..$ names: chr [1:2] "counts" "samples"


F0_cpm<-cpm(F0_RawHitCount_dge_TMM)
str(F0_cpm)
#num [1:32883, 1:23] 0 0 0 0 0 0 0 0 0 0 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:32883] "ENSRNOG00000046319" "ENSRNOG00000047964" "ENSRNOG00000050370" "ENSRNOG00000032365" ...
# ..$ : chr [1:23] "SL483888" "SL483889" "SL483890" "SL483891" ...

#Log2 transformed - an approximate way of getting rid of heteroskedasticity in the data and the traditional way of displaying RNA-Seq data in graphs
F0_lcpm<-cpm(F0_RawHitCount_dge_TMM, log=TRUE)

write.csv(F0_cpm,"F0_cpm.csv")
write.csv(F0_lcpm, "F0_lcpm.csv")

#################

#Filtering out rows of data that are 0 cpm for at least a quarter of the animals (e.g. all male bHRs and no one else)
#This is Megan erring on the side of including too much instead of including too little

sum(rowSums(F0_cpm>1)>=1)
#[1] 14802

keep.exprs<-rowSums(F0_cpm>1)>=6
F0_RawHitCount_ExpressionData_noLowHits<-F0_RawHitCount_ExpressionData[keep.exprs,]
str(F0_RawHitCount_ExpressionData_noLowHits)
# int [1:13786, 1:23] 1758 9524 219 401 344 753 419 174 12 994 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:13786] "ENSRNOG00000014303" "ENSRNOG00000014330" "ENSRNOG00000049505" "ENSRNOG00000014916" ...
# ..$ : chr [1:23] "SL483888" "SL483889" "SL483890" "SL483891" ...

#TMM normalization:

dge_noLowHits<-DGEList(counts=F0_RawHitCount_ExpressionData_noLowHits)
F0_RawHitCount_dge_noLowHits_TMM<-calcNormFactors(dge_noLowHits, method="TMM")


#create a matrix of log2 CPM (counts per million) values:

F0_cpm_noLowHits<-cpm(F0_RawHitCount_dge_noLowHits_TMM)

F0_lcpm_noLowHits<-cpm(F0_RawHitCount_dge_noLowHits_TMM, log=TRUE)


str(F0_lcpm_noLowHits)
# num [1:13786, 1:23] 6.69 9.12 3.69 4.56 4.34 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:13786] "ENSRNOG00000014303" "ENSRNOG00000014330" "ENSRNOG00000049505" "ENSRNOG00000014916" ...
# ..$ : chr [1:23] "SL483888" "SL483889" "SL483890" "SL483891" ...


#Making the log2 cpm values into a data.frame so we can join them with annotation:
F0_lcpm_noLowHits_DF<-as.data.frame(F0_lcpm_noLowHits, header=TRUE, stringAsFactors=FALSE)
#Make sure this step included the ENSEMBL ids as a column. If they are still row.names, we will need to make a column in the data.frame named ENSEMBL for those row.names
str(F0_lcpm_noLowHits_DF)
# data.frame':	13786 obs. of  23 variables:
#  $ SL483888: num  6.69 9.12 3.69 4.56 4.34 ...
#  $ SL483889: num  6.55 8.83 3.48 4.63 4.05 ...
#  $ SL483890: num  6.72 8.54 3.24 4.59 3.94 ...
#  $ SL483891: num  6.99 8.71 3.07 4.27 3.84 ...
#  $ SL483892: num  7.03 8.64 3.22 4.32 3.7 ...
#  $ SL483893: num  6.87 8.52 3.25 4.66 3.98 ...
#  $ SL483894: num  6.68 9.36 3.65 4.37 4.3 ...
#  $ SL483895: num  6.41 9.59 4.11 4.28 4.45 ...
#  $ SL483896: num  6.31 9.48 3.9 3.6 4.42 ...
#  $ SL483898: num  7.25 8.31 2.94 4.44 3.73 ...
#  $ SL483899: num  6.47 9.55 4.05 4.23 4.33 ...
#  $ SL483900: num  6.21 9.49 4.27 3.59 4.54 .
