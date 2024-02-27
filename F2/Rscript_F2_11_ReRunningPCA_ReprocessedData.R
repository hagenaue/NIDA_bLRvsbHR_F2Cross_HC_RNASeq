#F2 HC RNA-Seq Dataset
#11_Rerunning PCA on the reprocessed RNA-Seq data
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-04, updated later for a few figures for the paper.



#Running PCA on the QCed dataset:

pcaNormFilterednoOutliers<-prcomp(t(F2_lcpm_noLowHits))

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

# num [1:14061, 1:10] -0.00687 0.01348 0.0133 -0.00332 0.01027 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:14061] "ENSRNOG00000014303" "ENSRNOG00000014330" "ENSRNOG00000049505" "ENSRNOG00000014916" ...
# ..$ : chr [1:10] "PC1" "PC2" "PC3" "PC4" ...

str(F2_Annotated_lcpm_NoLowHits)
# data.frame':	14263 obs. of  249 variables:
#  $ SL469967: num  8.21 7.51 2.48 5.36 2.84 ...
#  $ SL469968: num  8.17 7.67 2.32 5.45 3.01 ...
#  $ SL469969: num  8.26 7.45 2.38 5.47 2.95 ...
#  $ SL469970: num  8.05 7.57 2.43 5.57 3.41 ...
#  $ SL469971: num  8.28 7.55 2.15 5.48 3.01 ...


#Output a scree plot for the PCA (no outliers):

png("10 PCA Scree Plot1.png")
plot(summary(pcaNormFilterednoOutliers)$importance[2,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off() 
png("10 PCA Scree Plot2.png")
plot(summary(pcaNormFilterednoOutliers)$importance[3,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()

###########2/24/2022 Elaine create higher quality plots 

pdf("PCA_Scree_Plot1.pdf", height=5, width=5)
plot(summary(pcaNormFilterednoOutliers)$importance[2,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC Number", ylab="Proportion of Variance Explained", col=2)
dev.off() 
pdf("PCA_Scree_Plot2.pdf", height=5, width=5)
plot(summary(pcaNormFilterednoOutliers)$importance[3,]~(c(1:length(summary(pcaNormFilterednoOutliers)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC Number", ylab="Cumulative Proportion of Variance Explained", col=3)
dev.off()


#####without main title
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


cor(cbind(F2_PCA_wMetaData[,c(2:7)], PC1noOutliers, PC2noOutliers, PC3noOutliers, PC4noOutliers))

#                         PC1          PC2           PC3          PC4          PC5           PC6 PC1noOutliers PC2noOutliers PC3noOutliers PC4noOutliers
# PC1            1.0000000000 -0.010327878 -0.0080127844  0.001271717 -0.011155607 -0.0006045998  9.984637e-01 -2.905259e-02  2.313957e-02  2.677659e-02
# PC2           -0.0103278776  1.000000000 -0.0069825235  0.008437913 -0.005805816 -0.0030803985  3.181878e-02  6.956421e-01 -4.805491e-01 -4.139370e-01
# PC3           -0.0080127844 -0.006982523  1.0000000000 -0.009479881 -0.021772185 -0.0038525463  1.770131e-04 -6.236059e-01 -5.071039e-01 -5.242836e-01
# PC4            0.0012717166  0.008437913 -0.0094798806  1.000000000 -0.020254474  0.0016666571 -5.910179e-03  8.534612e-02  3.507907e-02 -1.972726e-01
# PC5           -0.0111556066 -0.005805816 -0.0217721847 -0.020254474  1.000000000 -0.0026090005  5.482013e-03  1.529256e-01  1.591935e-01 -8.979801e-02
# PC6           -0.0006045998 -0.003080398 -0.0038525463  0.001666657 -0.002609000  1.0000000000  6.607469e-03  2.213845e-01  8.722473e-02 -1.471032e-01
# PC1noOutliers  0.9984636978  0.031818782  0.0001770131 -0.005910179  0.005482013  0.0066074693  1.000000e+00 -1.894461e-14  8.837664e-16  2.075643e-15
# PC2noOutliers -0.0290525928  0.695642072 -0.6236058592  0.085346119  0.152925572  0.2213844774 -1.894461e-14  1.000000e+00  3.046451e-16 -5.624643e-16
# PC3noOutliers  0.0231395687 -0.480549055 -0.5071038797  0.035079067  0.159193545  0.0872247310  8.837664e-16  3.046451e-16  1.000000e+00  7.188240e-16
# PC4noOutliers  0.0267765925 -0.413936987 -0.5242835835 -0.197272631 -0.089798009 -0.1471032278  2.075643e-15 -5.624643e-16  7.188240e-16  1.000000e+00

#will overwrite Fan's PC variables because they were created before outliers were removed

F2_PCA_wMetaData$PC1<-PC1noOutliers
F2_PCA_wMetaData$PC2<-PC2noOutliers
F2_PCA_wMetaData$PC3<-PC3noOutliers
F2_PCA_wMetaData$PC4<-PC4noOutliers
F2_PCA_wMetaData$PC5<-PC5noOutliers
F2_PCA_wMetaData$PC6<-PC6noOutliers

F2_PCA_wMetaData$PC7<-pcaNormFilterednoOutliers$x[,7]
F2_PCA_wMetaData$PC8<-pcaNormFilterednoOutliers$x[,8]
F2_PCA_wMetaData$PC9<-pcaNormFilterednoOutliers$x[,9]
F2_PCA_wMetaData$PC10<-pcaNormFilterednoOutliers$x[,10]


SubjectPCA<-cbind(PC1noOutliers, PC2noOutliers, PC3noOutliers, PC4noOutliers)

