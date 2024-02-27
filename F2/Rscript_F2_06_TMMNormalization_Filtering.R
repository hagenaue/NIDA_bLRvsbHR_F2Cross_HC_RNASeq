#F2 HC RNA-Seq Dataset
#06_TMM Normalization
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-04, updated later for a few figures for the paper.


#Moving forward with processing the RNA-Seq data:

#Adding in TMM normalization because we found that our top PCs were highly correlated with overall RNA levels (in this dataset and F0s)
#TMM normalization scales to correct for estimated relative RNA production levels at the point of calculating cpm

dge<-DGEList(counts=F2_RawHitCount_ExpressionData)
F2_RawHitCount_dge_TMM<-calcNormFactors(dge, method="TMM")

str(F2_RawHitCount_dge_TMM)
#Formal class 'DGEList' [package "edgeR"] with 1 slot
# ..@ .Data:List of 2
# .. ..$ : int [1:32883, 1:247] 0 0 0 0 0 0 2 0 0 0 ...
# .. .. ..- attr(*, "dimnames")=List of 2
# .. .. .. ..$ : chr [1:32883] "ENSRNOG00000046319" "ENSRNOG00000047964" "ENSRNOG00000050370" "ENSRNOG00000032365" ...
# .. .. .. ..$ : chr [1:247] "SL469967" "SL469968" "SL469969" "SL469970" ...
# .. ..$ :'data.frame':	247 obs. of  3 variables:
#   .. .. ..$ group       : Factor w/ 1 level "1": 1 1 1 1 1 1 1 1 1 1 ...
# .. .. ..$ lib.size    : num [1:247] 21826046 19217774 20269455 26635890 21240639 ...
# .. .. ..$ norm.factors: num [1:247] 0.952 0.953 0.965 1.017 0.97 ...
# ..$ names: chr [1:2] "counts" "samples"

F2_cpm<-cpm(F2_RawHitCount_dge_TMM)
str(F2_cpm)
# num [1:32883, 1:247] 0 0 0 0 0 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:32883] "ENSRNOG00000046319" "ENSRNOG00000047964" "ENSRNOG00000050370" "ENSRNOG00000032365" ...
# ..$ : chr [1:247] "SL469967" "SL469968" "SL469969" "SL469970" ...

#Log2 transformed - an approximate way of getting rid of heteroskedasticity in the data and the traditional way of displaying RNA-Seq data in graphs
F2_lcpm<-cpm(F2_RawHitCount_dge_TMM, log=TRUE)

write.csv(F2_cpm,"F2_cpm.csv")
write.csv(F2_lcpm, "F2_lcpm.csv")

##############

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
# int [1:14061, 1:247] 6484 3991 120 894 155 344 1235 1167 186 2232 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:14061] "ENSRNOG00000014303" "ENSRNOG00000014330" "ENSRNOG00000049505" "ENSRNOG00000014916" ...
# ..$ : chr [1:247] "SL469967" "SL469968" "SL469969" "SL469970" ...

#TMM normalization:

dge_noLowHits<-DGEList(counts=F2_RawHitCount_ExpressionData_noLowHits)
F2_RawHitCount_dge_noLowHits_TMM<-calcNormFactors(dge_noLowHits, method="TMM")

#create a matrix of log2 CPM (counts per million) values:

F2_cpm_noLowHits<-cpm(F2_RawHitCount_dge_noLowHits_TMM)

F2_lcpm_noLowHits<-cpm(F2_RawHitCount_dge_noLowHits_TMM, log=TRUE)

str(F2_lcpm_noLowHits)
# num [1:14061, 1:247] 8.21 7.51 2.48 5.36 2.84 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:14061] "ENSRNOG00000014303" "ENSRNOG00000014330" "ENSRNOG00000049505" "ENSRNOG00000014916" ...
# ..$ : chr [1:247] "SL469967" "SL469968" "SL469969" "SL469970" ...


#Making the log2 cpm values into a data.frame so we can join them with annotation:
F2_lcpm_noLowHits_DF<-as.data.frame(F2_lcpm_noLowHits, header=TRUE, stringAsFactors=FALSE)
#Make sure this step included the ENSEMBL ids as a column. If they are still row.names, we will need to make a column in the data.frame named ENSEMBL for those row.names
str(F2_lcpm_noLowHits_DF)

# data.frame':	14061 obs. of  247 variables:
#  $ SL469967: num  8.21 7.51 2.48 5.36 2.84 ...
#  $ SL469968: num  8.17 7.67 2.32 5.45 3.01 ...
#  $ SL469969: num  8.26 7.45 2.38 5.47 2.95 ...
#  $ SL469970: num  8.05 7.57 2.43 5.57 3.41 ...
#  $ SL469971: num  8.28 7.55 2.15 5.48 3.01 ...
