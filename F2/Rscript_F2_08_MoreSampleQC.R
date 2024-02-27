#F2 HC RNA-Seq Dataset
#08_More Sample Quality Control
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-04, updated later for a few figures for the paper.


#checking for sample mislabeling by comparing labeled sex to sex chromosome gene expression

colnames(F2_Annotated_lcpm_NoLowHits)

png("Eif2s3y_vs_Sex_customCDFplus2.png")
boxplot(as.numeric(F2_Annotated_lcpm_NoLowHits[which(F2_Annotated_lcpm_NoLowHits[,249]=="Eif2s3y"),c(1:247)])~F2_PCA_wMetaData$Sex, col=2)
dev.off()

png("Ddx3_vs_Sex_customCDFplus2.png")
boxplot(as.numeric(F2_Annotated_lcpm_NoLowHits[which(F2_Annotated_lcpm_NoLowHits[,249]=="Ddx3"),c(1:247)])~F2_PCA_wMetaData$Sex, col=2)
dev.off()

png("Eif2s3yvsDdx3_SexCheck.png")
plot(as.numeric(F2_Annotated_lcpm_NoLowHits[which(F2_Annotated_lcpm_NoLowHits[,249]=="Eif2s3y"),c(1:247)])~as.numeric(F2_Annotated_lcpm_NoLowHits[which(F2_Annotated_lcpm_NoLowHits[,249]=="Ddx3"),c(1:247)]), col=as.numeric(as.factor(F2_PCA_wMetaData$Sex)))
dev.off()

#Two samples have incongruous sex label vs. sex chromosome gene expression - we will consider these samples mislabeled:

colnames(F2_Annotated_lcpm_NoLowHits[,c(1:247)])[((as.numeric(F2_Annotated_lcpm_NoLowHits[which(F2_Annotated_lcpm_NoLowHits[,249]=="Eif2s3y"), c(1:247)])>4) & F2_PCA_wMetaData$Sex=="Female")]

#[1] "SL470859" "SL470871"


################

#Removing mislabeled samples: 
#Note- we are overwriting previous objects.

#ENSRNOG00000060048 is the ensembl gene name for Eif2s3y

str(F2_RawHitCount)
# 'data.frame':	32883 obs. of  247 variables:
#   $ SL469967: int  0 0 0 0 0 0 2 0 0 0 ...
# $ SL469968: int  0 0 0 0 0 0 4 0 0 0 ...
# $ SL469969: int  0 0 0 0 0 0 2 0 0 0 ...
# $ SL469970: int  0 0 0 0 0 0 5 0 0 0 ...
# $ SL469971: int  0 0 0 0 0 0 0 0 0 0 ...
# $ SL469972: int  0 0 0 0 2 0 1 0 0 0 ...

F2_RawHitCount_2<-F2_RawHitCount[,(colnames(F2_RawHitCount)=="SL470859"|colnames(F2_RawHitCount)=="SL470871")==FALSE]

str(F2_RawHitCount_2)
#'data.frame':	32883 obs. of  245 variables:

str(F2_RawHitCount_ExpressionData)
# int [1:32883, 1:247] 0 0 0 0 0 0 2 0 0 0 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:32883] "ENSRNOG00000046319" "ENSRNOG00000047964" "ENSRNOG00000050370" "ENSRNOG00000032365" ...
# ..$ : chr [1:247] "SL469967" "SL469968" "SL469969" "SL469970" ...

colnames(F2_RawHitCount_ExpressionData_2)

F2_RawHitCount_ExpressionData_2<-F2_RawHitCount_ExpressionData[,(colnames(F2_RawHitCount_ExpressionData)=="SL470859"|colnames(F2_RawHitCount_ExpressionData)=="SL470871")==FALSE]

str(F2_RawHitCount_ExpressionData_2)
#int [1:32883, 1:245] 0 0 0 0 0 0 2 0 0 0 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:32883] "ENSRNOG00000046319" "ENSRNOG00000047964" "ENSRNOG00000050370" "ENSRNOG00000032365" ...
# ..$ : chr [1:245] "SL469967" "SL469968" "SL469969" "SL469970" ...

F2_PCA_wMetaData_2<-F2_PCA_wMetaData[(colnames(F2_RawHitCount)=="SL470859"|colnames(F2_RawHitCount)=="SL470871")==FALSE,]

str(F2_PCA_wMetaData_2)

#'data.frame':	245 obs. of  65 variables

cbind(F2_PCA_wMetaData_2$row_name,colnames(F2_RawHitCount_2),colnames(F2_RawHitCount_ExpressionData_2))
# [,1]       [,2]       [,3]      
# [1,] "sl469967" "SL469967" "SL469967"
# [2,] "sl469968" "SL469968" "SL469968"
# [3,] "sl469969" "SL469969" "SL469969"
# [4,] "sl469970" "SL469970" "SL469970"
# [5,] "sl469971" "SL469971" "SL469971"
# [6,] "sl469972" "SL469972" "SL469972"
# [7,] "sl469973" "SL469973" "SL469973"
# [8,] "sl469974" "SL469974" "SL469974"
# [9,] "sl469975" "SL469975" "SL469975"
# [10,] "sl469976" "SL469976" "SL469976"

F2_RawHitCount<-F2_RawHitCount_2

F2_RawHitCount_ExpressionData<-F2_RawHitCount_ExpressionData_2

F2_PCA_wMetaData<-F2_PCA_wMetaData_2

#***************************************

str(F2_PCA_wMetaData)
#'data.frame':	245 obs. of  85 variables:

colnames(F2_PCA_wMetaData)

#Continuous variables: columns 8-13,15,18,22-28,30-38,52,54,56,60,61,65,66-68, 74-80

colnames(F2_PCA_wMetaData[,c(8:13,15,18,22:28,30:38,52,54,56,60,61,65,66:68,74:80)])

write.csv(cor(as.matrix(F2_PCA_wMetaData[,c(8:13,15,18,22:28,30:38,52,54,56,60,61,65,66:68,74:78,80)]), use="pairwise.complete.obs"), "F2_CorrelationMarix_ContinuousMetaData_AfterQC_N245.csv")

