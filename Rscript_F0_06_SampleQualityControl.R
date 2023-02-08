#F0 HC RNA-Seq Dataset
#06_Sample Quality Control
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-10, updated later for a few figures for the paper.

#Quality control: Checking for sample mislabeling by comparing labeled sex to sex chromosome gene expression


png("Eif2s3y_vs_Sex.png")
boxplot(as.numeric(F0_Annotated_lcpm_NoLowHits[which(F0_Annotated_lcpm_NoLowHits[,25]=="Eif2s3y"),c(1:23)])~F0_PCA_wMetaData$Sex, col=2)
dev.off()


png("Ddx3_vs_Sex.png")
boxplot(as.numeric(F0_Annotated_lcpm_NoLowHits[which(F0_Annotated_lcpm_NoLowHits[,25]=="Ddx3"),c(1:23)])~F0_PCA_wMetaData$Sex, col=2)
dev.off()


png("Eif2s3yvsDdx3_SexCheck.png")
plot(as.numeric(F0_Annotated_lcpm_NoLowHits[which(F0_Annotated_lcpm_NoLowHits[,25]=="Eif2s3y"),c(1:23)])~as.numeric(F0_Annotated_lcpm_NoLowHits[which(F0_Annotated_lcpm_NoLowHits[,25]=="Ddx3"),c(1:23)]), col=as.numeric(as.factor(F0_PCA_wMetaData$Sex)))
dev.off()


#Everything looks good.

#############################


#Running PCA on the QC-ed dataset:

pcaNormFilterednoOutliers<-prcomp(t(F0_lcpm_noLowHits))

tmp<-pcaNormFilterednoOutliers$x[,1:10]

write.table(tmp, "PCA_1_10.txt", sep="\t")

PCeigenvectors<-pcaNormFilterednoOutliers$rotation[ ,c(1:10)]
#PCeigenvectors2<-cbind(PCeigenvectors, F2_Annotated_lcpm_NoLowHits$SYMBOL)
#write.csv(PCeigenvectors2, "PCeigenvectors.csv")

PC1noOutliers<-pcaNormFilterednoOutliers$x[,1]
PC2noOutliers<-pcaNormFilterednoOutliers$x[,2]
PC3noOutliers<-pcaNormFilterednoOutliers$x[,3]
PC4noOutliers<-pcaNormFilterednoOutliers$x[,4]
PC5noOutliers<-pcaNormFilterednoOutliers$x[,5]
PC6noOutliers<-pcaNormFilterednoOutliers$x[,6]


str(PCeigenvectors)
# num [1:13786, 1:10] 0.00503 -0.00967 -0.00825 -0.0014 -0.00596 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:13786] "ENSRNOG00000014303" "ENSRNOG00000014330" "ENSRNOG00000049505" "ENSRNOG00000014916" ...
# ..$ : chr [1:10] "PC1" "PC2" "PC3" "PC4" ...

str(F0_Annotated_lcpm_NoLowHits)
# 'data.frame':	13996 obs. of  25 variables:
#   $ SL483888: num  6.69 9.12 3.69 4.56 4.34 ...
# $ SL483889: num  6.55 8.83 3.48 4.63 4.05 ...
# $ SL483890: num  6.72 8.54 3.24 4.59 3.94 ...
# $ SL483891: num  6.99 8.71 3.07 4.27 3.84 ...
# $ SL483892: num  7.03 8.64 3.22 4.32 3.7 ...
# $ SL483893: num  6.87 8.52 3.25 4.66 3.98 ...
# $ SL483894: num  6.68 9.36 3.65 4.37 4.3 ...
# $ SL483895: num  6.41 9.59 4.11 4.28 4.45 ...

#Output a scree plot for the PCA (no outliers):

png("10 PCA Scree Plot1.png")
plot(summary(pcaNormFilterednoOutliers)$importance[2,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off() 
png("10 PCA Scree Plot2.png")
plot(summary(pcaNormFilterednoOutliers)$importance[3,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()

######################2/24/2022 Elaine create higher quality plots

pdf("PCA_Scree_Plot1.pdf", height=5, width=5)
plot(summary(pcaNormFilterednoOutliers)$importance[2,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC Number", ylab="Proportion of Variance Explained", col=2)
dev.off() 

pdf("PCA_Scree_Plot2.pdf", height=5, width=5)
plot(summary(pcaNormFilterednoOutliers)$importance[3,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC Number", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()

pdf("PCA_Scree_Plot1_NoMainTitle.pdf", height=5, width=5)
plot(summary(pcaNormFilterednoOutliers)$importance[2,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), xlab="PC Number", ylab="Proportion of Variance Explained", col=2)
dev.off() 

pdf("PCA_Scree_Plot2_NoMainTitle.pdf", height=5, width=5)
plot(summary(pcaNormFilterednoOutliers)$importance[3,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), xlab="PC Number", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()



#Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
png("10 PC1 vs PC2.png")
plot(PC1noOutliers~PC2noOutliers, main="Principal Components Analysis of Normalized Filtered Data No Outliers")
dev.off() 

#No obvious outliers.

colnames(F0_PCA_wMetaData)

#Comparing the PCA output from Ahub to our new PCA output from QCed data:

cor(cbind(F0_PCA_wMetaData[,c(6:11)], PC1noOutliers, PC2noOutliers, PC3noOutliers, PC4noOutliers))
#                       PC1         PC2          PC3         PC4          PC5          PC6 PC1noOutliers PC2noOutliers PC3noOutliers PC4noOutliers
# PC1            1.00000000  0.04317879  0.018880124 -0.05649199  0.012799512  0.026962115  9.924250e-01 -1.179532e-01 -1.922090e-02  1.067576e-02
# PC2            0.04317879  1.00000000 -0.113457409  0.33948055 -0.076916835 -0.162024979  1.611138e-01  9.847524e-01  2.894287e-02  1.454407e-02
# PC3            0.01888012 -0.11345741  1.000000000  0.14843944 -0.033632241 -0.070846171  6.995486e-03 -1.428549e-01  5.503332e-01  6.336471e-01
# PC4           -0.05649199  0.33948055  0.148439436  1.00000000  0.100632400  0.211981718 -2.511331e-02  3.476500e-01 -1.743064e-01  6.290287e-01
# PC5            0.01279951 -0.07691683 -0.033632241  0.10063240  1.000000000 -0.048029151  3.819212e-03 -1.060935e-01  2.578526e-01  2.443028e-01
# PC6            0.02696211 -0.16202498 -0.070846171  0.21198172 -0.048029151  1.000000000  2.709621e-03 -1.819803e-01 -4.585954e-01  8.297605e-02
# PC1noOutliers  0.99242499  0.16111384  0.006995486 -0.02511331  0.003819212  0.002709621  1.000000e+00 -7.321781e-15  3.851359e-16 -5.378324e-16
# PC2noOutliers -0.11795325  0.98475244 -0.142854924  0.34765004 -0.106093518 -0.181980271 -7.321781e-15  1.000000e+00  1.318319e-14 -1.082194e-15
# PC3noOutliers -0.01922090  0.02894287  0.550333173 -0.17430638  0.257852586 -0.458595393  3.851359e-16  1.318319e-14  1.000000e+00  7.445356e-17
# PC4noOutliers  0.01067576  0.01454407  0.633647097  0.62902868  0.244302804  0.082976054 -5.378324e-16 -1.082194e-15  7.445356e-17  1.000000e+00

# Overwriting Ahubs PCA output because they were created before outliers were removed, and then we can more easily reuse (and interpret) code:

F0_PCA_wMetaData$PC1<-PC1noOutliers
F0_PCA_wMetaData$PC2<-PC2noOutliers
F0_PCA_wMetaData$PC3<-PC3noOutliers
F0_PCA_wMetaData$PC4<-PC4noOutliers
F0_PCA_wMetaData$PC5<-PC5noOutliers
F0_PCA_wMetaData$PC6<-PC6noOutliers

