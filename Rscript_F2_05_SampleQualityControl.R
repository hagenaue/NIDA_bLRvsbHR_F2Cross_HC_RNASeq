#F2 HC RNA-Seq Dataset
#05_Sample Quality Control
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-04, updated later for a few figures for the paper.


#Removing samples that failed QC:

#Removing the one sample that has crazy high library size aw well as 2 samples with low RINs compared to all of the other samples and then progressing forward again

#Low RIN: sl470590, sl470899

#High library size/PF.read: sl470850

#Note that I'm overwriting our original objects for metadata and gene expression which is lazy but means we don't have to change all of our downstream code.

str(F2_RawHitCount)
#'data.frame':	32883 obs. of  250 variables:
#'
str(F2_RawHitCount_ExpressionData)
# int [1:32883, 1:250] 0 0 0 0 0 0 2 0 0 0 ...

str(F2_PCA_wMetaData)
#'data.frame':	250 obs. of  65 variables:
#'
cbind(F2_PCA_wMetaData$row_name, colnames(F2_RawHitCount_ExpressionData))


F2_RawHitCount_2<-F2_RawHitCount[,(colnames(F2_RawHitCount)=="SL470590"|colnames(F2_RawHitCount)=="SL470899"|colnames(F2_RawHitCount)=="SL470850")==FALSE]

str(F2_RawHitCount_2)
# 'data.frame':	32883 obs. of  247 variables:
#   $ SL469967: int  0 0 0 0 0 0 2 0 0 0 ...
# $ SL469968: int  0 0 0 0 0 0 4 0 0 0 ...
# $ SL469969: int  0 0 0 0 0 0 2 0 0 0 ...
# $ SL469970: int  0 0 0 0 0 0 5 0 0 0 ...
# $ SL469971: int  0 0 0 0 0 0 0 0 0 0 ...

F2_RawHitCount_ExpressionData_2<-F2_RawHitCount_ExpressionData[,F2_LibrarySize<109000000 & F2_PCA_wMetaData$RIN>8 ]

F2_PCA_wMetaData_2<-F2_PCA_wMetaData[F2_PCA_wMetaData$LibrarySize<109000000 & F2_PCA_wMetaData$RIN>8,]

str(F2_RawHitCount_2)

#data.frame':	32883 obs. of  247 variables


str(F2_RawHitCount_ExpressionData_2)

# int [1:32883, 1:247] 0 0 0 0 0 0 2 0 0 0 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:32883] "ENSRNOG00000046319" "ENSRNOG00000047964" "ENSRNOG00000050370" "ENSRNOG00000032365" ...
# ..$ : chr [1:247] "SL469967" "SL469968" "SL469969" "SL469970" ...
# 

str(F2_PCA_wMetaData_2)

#'data.frame':	247 obs. of  65 variables

F2_RawHitCount<-F2_RawHitCount_2


F2_RawHitCount_ExpressionData<-F2_RawHitCount_ExpressionData_2

F2_PCA_wMetaData<-F2_PCA_wMetaData_2


#Now down to 247 rats after removing 1 rat with extreme library size and 2 rats with bad RIN.

cbind(F2_PCA_wMetaData$row_name, colnames(F2_RawHitCount_ExpressionData),colnames(F2_RawHitCount_2))
# 
# [,1]       [,2]       [,3]      
# [1,] "sl469967" "SL469967" "SL469967"
# [2,] "sl469968" "SL469968" "SL469968"
# [3,] "sl469969" "SL469969" "SL469969"
# [4,] "sl469970" "SL469970" "SL469970"
# [5,] "sl469971" "SL469971" "SL469971"
# [6,] "sl469972" "SL469972" "SL469972"


###########################


#This hypothetically is a graph of total reads in our data vs. total reads (following qc) in the discovery life sciences output


pdf("F2_LibrarySize_Vs_PFReads.pdf", height=5, width=5)
plot(LibrarySize~PF.Reads, data=F2_PCA_wMetaData, main="F2", xlab="PF.Reads")
BestFitLine<-lm(LibrarySize~PF.Reads, data=F2_PCA_wMetaData)
abline(BestFitLine, col=2, lwd=3)
dev.off()
#LibrarySize in our data vs. "P.F.Reads" from the output from the Discovery Life Sciences
#Does library size correlate with our PCs?

summary.lm(BestFitLine)
#Call:
# lm(formula = LibrarySize ~ PF.Reads, data = F2_PCA_wMetaData)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1011314  -270248    14330   282465   800436 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 2.442e+05  2.041e+05   1.197    0.233    
# PF.Reads    7.148e-01  6.654e-03 107.416   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 380900 on 245 degrees of freedom
# Multiple R-squared:  0.9792,	Adjusted R-squared:  0.9791 
# F-statistic: 1.154e+04 on 1 and 245 DF,  p-value: < 2.2e-16

