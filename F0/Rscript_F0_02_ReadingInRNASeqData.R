#F0 HC RNA-Seq Dataset
#02_Reading In RNA-Seq Data
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-10, updated later for a few figures for the paper.



#Reading in the RNA-Seq data:

#First: What RNA-Seq data do we have?

#Files:

#Elaine_es103rn_FCgene_raw
##This one appears to be the raw counts without any normalization or filtering


#Elaine_es103rn_FCgene_4R
##This one appears to be the raw counts without any normalization or filtering, with additional annotation columns removed so that it reads into R as a data matrix



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

library(limma)



list.files()

F0_RawHitCount<-read.delim("Elaine_F0_es103rn_FCgene_4R.txt", sep="\t", header=TRUE, row.names=1)

str(F0_RawHitCount)
# 'data.frame':	32883 obs. of  24 variables:
# $ SL483888: int  0 0 0 0 0 0 0 0 0 0 ...
# $ SL483889: int  0 0 0 0 1 0 4 0 0 0 ...
# $ SL483890: int  0 0 0 0 0 0 1 0 0 0 ...
# $ SL483891: int  0 0 0 0 0 0 1 0 0 0 ...
# $ SL483892: int  0 0 0 0 1 0 3 0 0 0 ...

#Rows=genes
#Columns=subjects
#values=# of reads per gene (transcript) per subject

#Double-check that the subjects are in the same order as your metadata


#Removing the data for the subject with low RIN ("SL483897"):

str(F0_RawHitCount[,colnames(F0_RawHitCount)!= "SL483897"])
# 'data.frame':	32883 obs. of  23 variables:
#   $ SL483888: int  0 0 0 0 0 0 0 0 0 0 ...
# $ SL483889: int  0 0 0 0 1 0 4 0 0 0 ...
# $ SL483890: int  0 0 0 0 0 0 1 0 0 0 ...
# $ SL483891: int  0 0 0 0 0 0 1 0 0 0 ...
# $ SL483892: int  0 0 0 0 1 0 3 0 0 0 ...
# $ SL483893: int  0 0 0 0 0 0 1 0 0 0 ...
# $ SL483894: int  0 0 0 0 3 0 1 0 0 0 ...
# $ SL483895: int  0 0 0 0 1 0 3 0 0 3 ...
# $ SL483896: int  0 0 0 0 0 0 1 0 0 0 ...
# $ SL483898: int  0 0 0 0 1 0 3 0 0 0 ...
# $ SL483899: int  0 0 0 0 0 0 0 0 0 0 ...
# $ SL483900: int  0 0 0 0 0 0 4 0 0 0 ...
# $ SL483901: int  0 0 0 0 0 0 0 0 0 0 ...
# $ SL483902: int  0 0 0 0 0 0 0 0 0 0 ...
# $ SL483903: int  1 0 0 0 1 0 1 0 0 1 ...
# $ SL483904: int  0 0 0 0 1 0 2 0 0 0 ...
# $ SL483905: int  0 0 0 0 0 0 1 0 0 0 ...
# $ SL483906: int  0 0 0 0 0 0 1 0 0 1 ...
# $ SL483907: int  0 0 0 0 0 0 1 0 0 0 ...
# $ SL483908: int  0 0 0 0 0 0 0 0 0 0 ...
# $ SL483909: int  0 0 0 0 0 0 4 0 0 1 ...
# $ SL483910: int  0 0 0 0 0 0 1 0 0 0 ...
# $ SL483911: int  0 0 0 0 0 0 3 0 0 0 ...


F0_RawHitCount<-F0_RawHitCount[,colnames(F0_RawHitCount)!= "SL483897"]


str(F0_PCA_wMetaData)

#'data.frame':	23 obs. of  25 variables:


#Let's check visually to see if things are in the same order:
cbind(colnames(F0_RawHitCount), F0_PCA_wMetaData$row_name)
#Looks good.
# [,1]       [,2]      
# [1,] "SL483888" "sl483888"
# [2,] "SL483889" "sl483889"
# [3,] "SL483890" "sl483890"
# [4,] "SL483891" "sl483891"
# [5,] "SL483892" "sl483892"
# [6,] "SL483893" "sl483893"
# [7,] "SL483894" "sl483894"
# [8,] "SL483895" "sl483895"
# [9,] "SL483896" "sl483896"
# [10,] "SL483898" "sl483898"
# [11,] "SL483899" "sl483899"
# [12,] "SL483900" "sl483900"
# [13,] "SL483901" "sl483901"
# [14,] "SL483902" "sl483902"
# [15,] "SL483903" "sl483903"
# [16,] "SL483904" "sl483904"
# [17,] "SL483905" "sl483905"
# [18,] "SL483906" "sl483906"
# [19,] "SL483907" "sl483907"
# [20,] "SL483908" "sl483908"
# [21,] "SL483909" "sl483909"
# [22,] "SL483910" "sl483910"
# [23,] "SL483911" "sl483911"


##################

#Starting to work with the expression data:

F0_RawHitCount_ExpressionData<-as.matrix(F0_RawHitCount)
str(F0_RawHitCount_ExpressionData)
# int [1:32883, 1:23] 0 0 0 0 0 0 0 0 0 0 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:32883] "ENSRNOG00000046319" "ENSRNOG00000047964" "ENSRNOG00000050370" "ENSRNOG00000032365" ...
# ..$ : chr [1:23] "SL483888" "SL483889" "SL483890" "SL483891" ...


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

library(edgeR)


F0_LibrarySize<-apply(F0_RawHitCount_ExpressionData, 2, sum)

#This hypothetically is a graph of total reads in our data vs. total reads (following qc) in the discovery life sciences output
#but it doesn't exactly match - are we missing a qc step?
#did they take out everything that wasn't mRNA?
#... and maybe that is why things like mRNA % correlate so strongly with our top principal components of variation?


pdf("F0_LibrarySize_Vs_PFReads.pdf", height=5, width=5)
plot(F0_LibrarySize~F0_PCA_wMetaData$PF_Reads, main="F0", xlab="PF.Reads")
BestFitLine<-lm(F0_LibrarySize~F0_PCA_wMetaData$PF_Reads)
abline(BestFitLine, col=2, lwd=3)
dev.off()
#LibrarySize in our data vs. "P.F.Reads" from the output from the Discovery Life Sciences
#Doesn't quite match; does library size correlate with our PCs?

summary.lm(BestFitLine)

# Call:
#   lm(formula = F0_LibrarySize ~ F0_PCA_wMetaData$PF_Reads)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -2366209  -380991   241491   566206  1139329 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)               3.962e+05  1.710e+06   0.232    0.819    
# F0_PCA_wMetaData$PF_Reads 6.931e-01  5.340e-02  12.978 1.69e-11 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 918800 on 21 degrees of freedom
# Multiple R-squared:  0.8891,	Adjusted R-squared:  0.8839 
# F-statistic: 168.4 on 1 and 21 DF,  p-value: 1.695e-11


F0_PCA_wMetaData$LibrarySize<-F0_LibrarySize


pdf("Histogram_F0_LibrarySize.pdf", height=5, width=5)
hist(F0_LibrarySize)
dev.off()

summary(F0_LibrarySize)
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#16824899 20332912 22044668 22454305 24276095 28443049 

#Median Library Size = close to what was intended!  Whoot!
#22044668 
#... and not a huge amount of variability around that.

