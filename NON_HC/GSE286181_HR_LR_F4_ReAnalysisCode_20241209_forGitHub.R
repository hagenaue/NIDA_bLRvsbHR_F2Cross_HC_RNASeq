# Microarray study on adult male bHR/bLR from the F4 generation
# Megan Hagenauer
# Dec 9, 2024
#####################
#Notes:
# This code was adapted from "MBNI_AffymetrixRae230_F4" (Isabelle Birt)
# That code was used to re-analyze the hippocampal microarray data for the bHR/bLR hippocampal meta-analysis
# I'm returning now to also re-analyze the hypothalamus and cortex samples
########################

# This file includes code for re-analyzing the dataset of 6 HR and LR adult rats from a study performed by Dr. John Stead. The data is in CEL format and must be normalized, checked for quality, and annotated before being analyzed. A custom cdf is used for annotation.

################

# Before Running, be sure that the appropriate affymetrix packages are installed.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("affy")

library(affy)
library(org.Rn.eg.db)
library(plyr)
library(car)
library(limma)

#Custom chip definition files (.cdfs) downloaded from:

#http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/CDF_download.asp

# Version 25 (Data Sources)
# Released on Jan 5, 2021
#http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/25.0.0/entrezg.asp

#Make sure packages are installed

install.packages("rae230arnentrezg.db_25.0.0.tar.gz.crdownload", repos=NULL, type="source")
install.packages("rae230arnentrezgcdf_25.0.0.tar.gz.crdownload", repos=NULL, type="source")
install.packages("rae230arnentrezgprobe_25.0.0.tar.gz.crdownload", repos=NULL, type="source")

library(rae230arnentrezg.db)
library(rae230arnentrezgcdf)
library(rae230arnentrezgprobe)

library("affy")

#Reading in the .CEL files containing gene expression data

setwd("~/University of Michigan Dropbox/Megan Hagenauer/Wakil_Shared_Animal/bHR_bLR_Microarray_John Stead/Raw array data")

#Read .CEL files
data2<-ReadAffy(cdfname = "rae230arnentrezgcdf", celfile.path=("~/University of Michigan Dropbox/Megan Hagenauer/Wakil_Shared_Animal/bHR_bLR_Microarray_John Stead/Raw array data"))

eset2 <- rma(data2)

setwd("~/University of Michigan Dropbox/Megan Hagenauer/Wakil_Shared_Animal/bHR_bLR_Microarray_John Stead/Reanalysis_20241209")

#output customCDF
write.exprs(eset2,file="HR_LR_F4_data_customCDF.txt")

#Read in data with custom CDF so that it is in easy format
RMAExpression_customCDF<-read.delim("HR_LR_F4_data_customCDF.txt", sep="\t")

#Creating Annotation File
# The following code maps entrez ID to their respective gene symbols in order to assign gene symbols to the data
str(RMAExpression_customCDF)
#'data.frame':	10069 obs. of  37 variables:

head(RMAExpression_customCDF)

#Added "_at" at end of annotation needs to be removed
RMAExpression_EntrezID<-sub("_at", "", RMAExpression_customCDF[,1])
head(RMAExpression_EntrezID)

RMAExpression_customCDFAnnotation<-data.frame(RMAExpression_customCDF[,1], RMAExpression_EntrezID, stringsAsFactors = F )

colnames(RMAExpression_customCDFAnnotation)<-c("ProbesetID", "EntrezGeneID")
head(RMAExpression_customCDFAnnotation)

#Getting Gene Symbols:
x <- org.Rn.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

#Checking structure
xx[1]

#Assigning GeneSymbols to EntrezID
GeneSymbol<-unlist(xx, use.names=FALSE)
EntrezGeneID<-rep(names(xx), lengths(xx))

table(lengths(xx))

#There's a 1:1 mapping of gene symbol and entrez ID

EntrezVsGeneSymbol<-data.frame(EntrezGeneID, GeneSymbol, stringsAsFactors=F)

#Joining Gene symbol annotation with ProbesetID
RMAExpression_customCDFAnnotation2<-join(RMAExpression_customCDFAnnotation, EntrezVsGeneSymbol, by="EntrezGeneID", type="left")

head(RMAExpression_customCDFAnnotation2)

sum(is.na(RMAExpression_customCDFAnnotation2[,3])==F) 
#Number of probes with corresponding gene symbol name
#[1] 9910

dim(RMAExpression_customCDFAnnotation2)
#[1] 10069     3

#So almost all of the results have associated gene symbols now
#Output custom annotation file
write.csv(RMAExpression_customCDFAnnotation2,  "RMAExpression_customCDFAnnotation2.csv")


## Quality Control
# This code chunk checks the quality of the data and checks for outliers

#making numeric data matrix
tempExpressionMatrix<-as.matrix(RMAExpression_customCDF[,-1])
row.names(tempExpressionMatrix)<-RMAExpression_customCDF[,1]

SampleSample_CorMatrix<-cor(tempExpressionMatrix)

pdf("Heatmap_CorMatrix_bHR_bLR_F4.pdf", height=10, width=10)
heatmap(cor(SampleSample_CorMatrix), cex.axis=0.5)
dev.off()

pdf("Boxplots_SampleSample_CorMatrix.pdf", height=5, width=20)
boxplot(SampleSample_CorMatrix, las=2, cex.axis=0.25)
dev.off()

#Clear outliers:
#X13F4.H07M1.1.HPC.HR.1.CEL
#X22F4.L08M2.0.HPC.LR.5.CEL

median(SampleSample_CorMatrix)
#[1] 0.9524552

quantile(SampleSample_CorMatrix, probs=0.10)
#      10% 
#0.9040401

quantile(SampleSample_CorMatrix, probs=0.25)
#      25% 
# 0.9433006 

summary(SampleSample_CorMatrix)

# X13F4.H07M1.1.HPC.HR.1.CEL X14F4.L07M3.1.HPC.LR.1.CEL
# Min.   :0.6434             Min.   :0.6980            
# 1st Qu.:0.6693             1st Qu.:0.9474            
# Median :0.6885             Median :0.9785            
# Mean   :0.6980             Mean   :0.9620            
# 3rd Qu.:0.7020             3rd Qu.:0.9825            
# Max.   :1.0000             Max.   :1.0000  

# X21F4.H01M2.0.HPC.HR.5.CEL X22F4.L08M2.0.HPC.LR.5.CEL
# Min.   :0.7017             Min.   :0.8494            
# 1st Qu.:0.9533             1st Qu.:0.8670            
# Median :0.9807             Median :0.8911            
# Mean   :0.9649             Mean   :0.8874            
# 3rd Qu.:0.9848             3rd Qu.:0.9013            
# Max.   :1.0000             Max.   :1.0000 


#Output of sample-sample correlation matrix
write.csv(SampleSample_CorMatrix, "SamplevsSample_CorMatrix.csv")

## Generating boxplot to observe subject expression distribution

pdf("Boxplot_Expression_HR_LR_F4.pdf", height=5, width=20)
boxplot(data.frame(tempExpressionMatrix), cex=0.25, las=3, par(cex.axis=0.75), main="Boxplot of log signal values per sample (1 box=all filtered probes)", xlab="Sample ID", ylab="Log2 Signal")
dev.off()

#Subjects 13F4-H07M1.1-HPC-HR-1 and 22F4-L08M2.0-HPC-LR-5 are outliers identified by John Stead, also appear in the cor matrix and expression boxplots
#Need to remove both outliers and rename column "X" to ProbesetID

colnames(tempExpressionMatrix)

ExpressionMatrix_QC<-tempExpressionMatrix[,c(-13, -22)]
  
colnames(RMAExpression_customCDF)

Affy_F4_Data<-RMAExpression_customCDF[,c(-14,-23)]
colnames(Affy_F4_Data)[1]<-"ProbesetID"   #Renaming column X as probeset ID to be easily joined with annotation

str(Affy_F4_Data)


#This function breaks apart the colnames into sample metadata variables:
#To make that into an easier-to-use data.frame
SampleMetaData<-do.call(rbind.data.frame, strsplit(colnames(ExpressionMatrix_QC), "\\."))

str(SampleMetaData)
head(SampleMetaData)

#I'm going to rename these:
colnames(SampleMetaData)<-c("Sample", "Subject", "Unk", "Region", "Phenotype", "Unk2", "FileType")

str(SampleMetaData)


##########################

##PCA

#Run PCA on the transposed gene expression matrix
#In this case, ExpressionData_Subset_noBad_Filtered is the object containing the gene expression matrix
pca_output<-prcomp(t(ExpressionMatrix_QC), scale=TRUE)

#Output a scree plot for the PCA results:

#This plot illustrates the proportion of variance explained by each principal component (PC):
png("PCA Scree Plot1.png")
plot(summary(pca_output)$importance[2,]~(c(1:length(summary(pca_output)$importance[2,]))), main="Variance Explained by Each Principal Component", xlab="PC #", ylab="Proportion of Variance Explained", col=2)
dev.off()

#Extract out the first four principal components and make them named R objects for graphing purposes:
PC1<-pca_output$x[,1]
PC2<-pca_output$x[,2]
PC3<-pca_output$x[,3]
PC4<-pca_output$x[,4]
PC5<-pca_output$x[,5]
PC6<-pca_output$x[,6]
PC7<-pca_output$x[,7]

#Extract the eigenvectors for PC1-4 and make them a named R object
PCeigenvectors<-pca_output$rotation[ ,c(1:4)]

#Add the gene names (rowData) from the gene expression matrix to the object with the eigenvectors
PCeigenvectors2<-cbind(PCeigenvectors, row.names(ExpressionMatrix_QC))

#Write this information out as a .csv file:
write.csv(PCeigenvectors2, "PCeigenvectors.csv")


#Output a scatterplot illustrating the relationship between Principal components 1 & 2 (PC1 & PC2):
pdf("PC1vsPC2_ByRegion.pdf", height=5, width=5)
plot(PC1~PC2, main="Principal Components Analysis", col=as.factor(SampleMetaData$Region), pch=16)
dev.off()
#Very nice separation by region

pdf("PC1vsPC2_ByPhenotype.pdf", height=5, width=5)
plot(PC1~PC2, main="Principal Components Analysis", pch=16, col=as.factor(SampleMetaData$Phenotype))
dev.off()
#Not much separation by phenotype

pdf("PC2vsPC3_ByRegion.pdf", height=5, width=5)
plot(PC2~PC3, main="Principal Components Analysis", col=as.factor(SampleMetaData$Region), pch=16)
dev.off()
#Very pretty

pdf("PC2vsPC3_ByPhenotype.pdf", height=5, width=5)
plot(PC2~PC3, main="Principal Components Analysis", col=as.factor(SampleMetaData$Phenotype), pch=16)
dev.off()
#Nope, not much there for phenotype either.

pdf("PC3vsPC4_ByRegion.pdf", height=5, width=5)
plot(PC3~PC4, main="Principal Components Analysis", col=as.factor(SampleMetaData$Region), pch=16)
dev.off()

pdf("PC3vsPC4_ByPhenotype.pdf", height=5, width=5)
plot(PC3~PC4, main="Principal Components Analysis", col=as.factor(SampleMetaData$Phenotype), pch=16)
dev.off()

#I wonder if there isn't a phenotype effect showing up because the effect is so different in different regions or just because this is an F4 dataset

pdf("PC4vsPC5_ByRegion.pdf", height=5, width=5)
plot(PC4~PC5, main="Principal Components Analysis", col=as.factor(SampleMetaData$Region), pch=16)
dev.off()

pdf("PC4vsPC5_ByPhenotype.pdf", height=5, width=5)
plot(PC4~PC5, main="Principal Components Analysis", col=as.factor(SampleMetaData$Phenotype), pch=16)
dev.off()
#nope

pdf("PC5vsPC6_ByRegion.pdf", height=5, width=5)
plot(PC5~PC6, main="Principal Components Analysis", col=as.factor(SampleMetaData$Region), pch=16)
dev.off()

pdf("PC5vsPC6_ByPhenotype.pdf", height=5, width=5)
plot(PC5~PC6, main="Principal Components Analysis", col=as.factor(SampleMetaData$Phenotype), pch=16)
dev.off()
#PC6 might finally include phenotype, but it's not beautiful

boxplot(PC6~SampleMetaData$Phenotype)
#Nah, not super convincing. 

t.test(PC6~SampleMetaData$Phenotype)
# Welch Two Sample t-test
# 
# data:  PC6 by SampleMetaData$Phenotype
# t = -1.7768, df = 25.183, p-value = 0.08768
# alternative hypothesis: true difference in means between group HR and group LR is not equal to 0
# 95 percent confidence interval:
#   -22.290137   1.638656
# sample estimates:
#   mean in group HR mean in group LR 
# -5.16287          5.16287

pdf("PC6vsPC7_ByRegion.pdf", height=5, width=5)
plot(PC6~PC7, main="Principal Components Analysis", col=as.factor(SampleMetaData$Region), pch=16)
dev.off()

pdf("PC6vsPC7_ByPhenotype.pdf", height=5, width=5)
plot(PC6~PC7, main="Principal Components Analysis", col=as.factor(SampleMetaData$Phenotype), pch=16)
dev.off()

boxplot(PC7~SampleMetaData$Phenotype)
#Maybe?

t.test(PC7~SampleMetaData$Phenotype)

# Welch Two Sample t-test
# 
# data:  PC7 by SampleMetaData$Phenotype
# t = -1.473, df = 25.688, p-value = 0.1529
# alternative hypothesis: true difference in means between group HR and group LR is not equal to 0
# 95 percent confidence interval:
#   -19.136853   3.164528
# sample estimates:
#   mean in group HR mean in group LR 
# -3.993081         3.993081 

#Nope.

##########################

#Alright, I'm just going to move on to preparing the differential expression output for inclusion in a meta-analysis.

#Adding the additional annotation:

#Reading in annotation file
Affy_F4_Annotation<-read.csv("RMAExpression_customCDFAnnotation2.csv", header=TRUE, stringsAsFactors = FALSE)

head(Affy_F4_Annotation)

#removing unecessary first column
Affy_F4_Annotation<-Affy_F4_Annotation[,-1]

#Joining Gene Symbol annotation with expression data

Affy_F4_Data<-join(Affy_F4_Annotation, Affy_F4_Data, by="ProbesetID", type="inner")

dim(Affy_F4_Data)
#[1] 10069    37

write.csv(Affy_F4_Data, "Affy_F4_Data.csv")

###############

SampleMetaData$Phenotype_Factor<-as.factor(SampleMetaData$Phenotype)

levels(SampleMetaData$Phenotype_Factor)
#[1] "HR" "LR"
#HR is reference like usual

#Subsetting the data by brain region

colnames(ExpressionMatrix_QC)

CTXMatrix_QC<-ExpressionMatrix_QC[,c(1:12)]
HPCMatrix_QC<-ExpressionMatrix_QC[,c(13:22)]
HYPMatrix_QC<-ExpressionMatrix_QC[,c(23:34)]

SampleMetaData$Region

CTXMetaData<-SampleMetaData[SampleMetaData$Region=="CTX",]
levels(CTXMetaData$Phenotype_Factor)
#[1] "HR" "LR"

HPCMetaData<-SampleMetaData[SampleMetaData$Region=="HPC",]
levels(HPCMetaData$Phenotype_Factor)
#[1] "HR" "LR"

HYPMetaData<-SampleMetaData[SampleMetaData$Region=="HYP",]
levels(HYPMetaData$Phenotype_Factor)
#[1] "HR" "LR"

str(CTXMetaData)

#Creating statistical model matrices:
CTXdesign <- model.matrix(~Phenotype_Factor, data=CTXMetaData)
HPCdesign <- model.matrix(~Phenotype_Factor, data=HPCMetaData)
HYPdesign <- model.matrix(~Phenotype_Factor, data=HYPMetaData)

#Applying the differential expression model to all genes using limma:
CTXfit <- lmFit(CTXMatrix_QC, CTXdesign)
HPCfit <- lmFit(HPCMatrix_QC, HPCdesign)
HYPfit <- lmFit(HYPMatrix_QC, HYPdesign)

#Adding an eBayes correction to help reduce the influence of outliers/small sample size on estimates
CTXefit <- eBayes(CTXfit)
HPCefit <- eBayes(HPCfit)
HYPefit <- eBayes(HYPfit)

#Applying a false discovery rate (FDR) correction and writing out the results as a .txt file:
write.fit(CTXefit, adjust="BH", file="CTXLimma_results.txt")
write.fit(HPCefit, adjust="BH", file="HPCLimma_results.txt")
write.fit(HYPefit, adjust="BH", file="HYPLimma_results.txt")

#Writing out the gene annotation for the results as a .csv file:
write.csv(Affy_F4_Data[,c(1:3)], "Annotation_CTXLimmaResults.csv")
write.csv(Affy_F4_Data[,c(1:3)], "Annotation_HPCLimmaResults.csv")
write.csv(Affy_F4_Data[,c(1:3)], "Annotation_HYPLimmaResults.csv")

CTXLimma<-read.delim("CTXLimma_results.txt", sep="\t", header=TRUE)
HPCLimma<-read.delim("HPCLimma_results.txt", sep="\t", header=TRUE)
HYPLimma<-read.delim("HYPLimma_results.txt", sep="\t", header=TRUE)

head(CTXLimma)
head(HPCLimma)
head(HYPLimma)

sum(CTXLimma$P.value.adj.Phenotype_FactorLR<0.10)
#[1] 0
sum(HPCLimma$P.value.adj.Phenotype_FactorLR<0.10)
#[1] 0
sum(HYPLimma$P.value.adj.Phenotype_FactorLR<0.10)
#[1] 0

hist(CTXLimma$P.value.Phenotype_FactorLR, breaks=100)
hist(HPCLimma$P.value.Phenotype_FactorLR, breaks=100)
hist(HYPLimma$P.value.Phenotype_FactorLR, breaks=100)
#Nada - flat or even inverted

##############################

#Extracting results for a meta-analysis (maybe??)

head(CTXLimma)
head(Affy_F4_Data[,c(1:3)]) 

CTXLimma_DE<-data.frame(Affy_F4_Data[,c(1:3)], CTXLimma)
HPCLimma_DE<-data.frame(Affy_F4_Data[,c(1:3)], HPCLimma)
HYPLimma_DE<-data.frame(Affy_F4_Data[,c(1:3)], HYPLimma)

write.csv(CTXLimma_DE, "CTXLimma_DE.csv")
write.csv(HPCLimma_DE, "HPCLimma_DE.csv")
write.csv(HYPLimma_DE, "HYPLimma_DE.csv")


##############

setwd("~/University of Michigan Dropbox/Megan Hagenauer/Wakil_Shared_Animal/bHR_bLR_Microarray_John Stead/Reanalysis_20241209/CTX")

#Renaming things to recycle Brain Data Alchemy Project code:

DE_Results<-CTXLimma_DE

colnames(DE_Results)
# [1] "ProbesetID"                     "EntrezGeneID"                  
# [3] "GeneSymbol"                     "X"                             
# [5] "AveExpr"                        "Coef..Intercept."              
# [7] "Coef.Phenotype_FactorLR"        "t..Intercept."                 
# [9] "t.Phenotype_FactorLR"           "P.value..Intercept."           
# [11] "P.value.Phenotype_FactorLR"     "P.value.adj..Intercept."       
# [13] "P.value.adj.Phenotype_FactorLR" "F"                             
# [15] "F.p.value"     


colnames(DE_Results)[3]<-"GeneSymbol"
colnames(DE_Results)[2]<-"NCBIid"

#Using the Brain Data Alchemy Function (v.2024):

FilteringDEResults_GoodAnnotation<-function(DE_Results){
  
  print("# of rows in results")
  print(nrow(DE_Results))
  
  print("# of rows with missing NCBI annotation:")
  print(sum(DE_Results$NCBIid==""|DE_Results$NCBIid=="null"))
  
  print("# of rows with NA NCBI annotation:")
  print(sum(is.na(DE_Results$NCBIid)))
  
  print("# of rows with missing Gene Symbol annotation:")
  print(sum(DE_Results$GeneSymbol==""|DE_Results$GeneSymbol=="null"))
  
  print("# of rows mapped to multiple NCBI_IDs:")
  print(length(grep('\\|', DE_Results$NCBIid)))
  
  print("# of rows mapped to multiple Gene Symbols:")
  print(length(grep('\\|', DE_Results$GeneSymbol)))
  
  #I only want the subset of data which contains rows that do not contain an NCBI EntrezID of ""
  DE_Results_NoNA<-DE_Results[(DE_Results$NCBIid==""|DE_Results$NCBIid=="null")==FALSE & is.na(DE_Results$NCBIid)==FALSE,]
  
  #I also only want the subset of data that is annotated with a single gene (not ambiguously mapped to more than one gene)
  if(length(grep('\\|', DE_Results_NoNA$NCBIid))==0){
    DE_Results_GoodAnnotation<<-DE_Results_NoNA
  }else{
    #I only want rows annotated with a single Gene Symbol (no pipe):
    DE_Results_GoodAnnotation<<-DE_Results_NoNA[-(grep('\\|', DE_Results_NoNA$NCBIid)),]
  }
  #I used a double arrow in that conditional to place DE_Results_GoodAnnotation back out in the environment outside the function 
  
  print("# of rows with good annotation")
  print(nrow(DE_Results_GoodAnnotation))
  
  #For record keeping (sometimes useful for troubleshooting later)
  write.csv(DE_Results_GoodAnnotation, "DE_Results_GoodAnnotation.csv")
  
  rm(DE_Results_NoNA, DE_Results)
  
  print("Outputted object: DE_Results_GoodAnnotation")
}

FilteringDEResults_GoodAnnotation(DE_Results)

# [1] "# of rows in results"
# [1] 10069
# [1] "# of rows with missing NCBI annotation:"
# [1] 0
# [1] "# of rows with NA NCBI annotation:"
# [1] 0
# [1] "# of rows with missing Gene Symbol annotation:"
# [1] NA
# [1] "# of rows mapped to multiple NCBI_IDs:"
# [1] 0
# [1] "# of rows mapped to multiple Gene Symbols:"
# [1] 0
# [1] "# of rows with good annotation"
# [1] 10069
# [1] "Outputted object: DE_Results_GoodAnnotation"

colnames(DE_Results)

#We need to extract the Log2FC ("Coef") and T-statistic ("t.") columns for the statistical contrasts relevant to our meta-analysis and place them into their own matrix

FoldChanges<-cbind(DE_Results_GoodAnnotation$Coef.Phenotype_FactorLR)

Tstats<-cbind(DE_Results_GoodAnnotation$t.Phenotype_FactorLR)

#Making the row names for the Log2FC and Tstat matrices the Entrez ID gene annotation:
row.names(FoldChanges)<-DE_Results_GoodAnnotation$NCBIid

row.names(Tstats)<-DE_Results_GoodAnnotation$NCBIid

#Let's rename our columns to something nicer describing the effect of interest:
#Note - we later discovered that this name needs to include the dataset identifier (GSEID#) for later joining and plotting purposes

ComparisonsOfInterest<-c("F4_CTX_bLR_vs_bHR")

colnames(FoldChanges)<-ComparisonsOfInterest
colnames(Tstats)<-ComparisonsOfInterest

ExtractingDEResults<-function(GSE_ID, FoldChanges, Tstats){
  
  #We calculate the standard error by dividing the log2FC by the tstat
  StandardErrors<-FoldChanges/Tstats
  str(StandardErrors)
  
  #For running our meta-analysis, we are actually going to need the sampling variance instead of the standard error
  #The sampling variance is just the standard error squared.
  
  SamplingVars<-(StandardErrors)^2
  str(SamplingVars)
  
  TempMasterResults<-list(Log2FC=FoldChanges, Tstat=Tstats, SE=StandardErrors, SV=SamplingVars)
  
  assign(paste("DEResults", GSE_ID, sep="_"), TempMasterResults, envir = as.environment(1))
  
  print(paste("Output: Named DEResults", GSE_ID, sep="_"))
  
  rm(TempMasterResults, SamplingVars, StandardErrors, FoldChanges, Tstats)
  
}

ExtractingDEResults("F4_CTX", FoldChanges, Tstats)

table(table(names(DEResults_GSE86893[[1]][,1])))
# 1 
# 17887 

# num [1:10069, 1] 0.059 0.0507 0.0724 0.0679 0.0476 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:10069] "100036582" "100036765" "100125364" "100125370" ...
# ..$ : chr "bLR_vs_bHR"
# num [1:10069, 1] 0.00349 0.00257 0.00525 0.00461 0.00227 ...
# - attr(*, "dimnames")=List of 2
# ..$ : chr [1:10069] "100036582" "100036765" "100125364" "100125370" ...
# ..$ : chr "bLR_vs_bHR"
# [1] "Output: Named DEResults_F4_Cortex"

###################

DE_Results<-HPCLimma_DE

colnames(DE_Results)

colnames(DE_Results)[3]<-"GeneSymbol"
colnames(DE_Results)[2]<-"NCBIid"

FilteringDEResults_GoodAnnotation(DE_Results)

# [1] "# of rows in results"
# [1] 10069
# [1] "# of rows with missing NCBI annotation:"
# [1] 0
# [1] "# of rows with NA NCBI annotation:"
# [1] 0
# [1] "# of rows with missing Gene Symbol annotation:"
# [1] NA
# [1] "# of rows mapped to multiple NCBI_IDs:"
# [1] 0
# [1] "# of rows mapped to multiple Gene Symbols:"
# [1] 0
# [1] "# of rows with good annotation"
# [1] 10069
# [1] "Outputted object: DE_Results_GoodAnnotation"

FoldChanges<-cbind(DE_Results_GoodAnnotation$Coef.Phenotype_FactorLR)

Tstats<-cbind(DE_Results_GoodAnnotation$t.Phenotype_FactorLR)

#Making the row names for the Log2FC and Tstat matrices the Entrez ID gene annotation:
row.names(FoldChanges)<-DE_Results_GoodAnnotation$NCBIid

row.names(Tstats)<-DE_Results_GoodAnnotation$NCBIid

#Let's rename our columns to something nicer describing the effect of interest:
#Note - we later discovered that this name needs to include the dataset identifier (GSEID#) for later joining and plotting purposes

ComparisonsOfInterest<-c("F4_HPC_bLR_vs_bHR")

colnames(FoldChanges)<-ComparisonsOfInterest
colnames(Tstats)<-ComparisonsOfInterest

ExtractingDEResults("F4_HPC", FoldChanges, Tstats)

###################

DE_Results<-HYPLimma_DE

colnames(DE_Results)

colnames(DE_Results)[3]<-"GeneSymbol"
colnames(DE_Results)[2]<-"NCBIid"

FilteringDEResults_GoodAnnotation(DE_Results)

# [1] "# of rows in results"
# [1] 10069
# [1] "# of rows with missing NCBI annotation:"
# [1] 0
# [1] "# of rows with NA NCBI annotation:"
# [1] 0
# [1] "# of rows with missing Gene Symbol annotation:"
# [1] NA
# [1] "# of rows mapped to multiple NCBI_IDs:"
# [1] 0
# [1] "# of rows mapped to multiple Gene Symbols:"
# [1] 0
# [1] "# of rows with good annotation"
# [1] 10069
# [1] "Outputted object: DE_Results_GoodAnnotation"

FoldChanges<-cbind(DE_Results_GoodAnnotation$Coef.Phenotype_FactorLR)

Tstats<-cbind(DE_Results_GoodAnnotation$t.Phenotype_FactorLR)

#Making the row names for the Log2FC and Tstat matrices the Entrez ID gene annotation:
row.names(FoldChanges)<-DE_Results_GoodAnnotation$NCBIid

row.names(Tstats)<-DE_Results_GoodAnnotation$NCBIid

#Let's rename our columns to something nicer describing the effect of interest:
#Note - we later discovered that this name needs to include the dataset identifier (GSEID#) for later joining and plotting purposes

ComparisonsOfInterest<-c("F4_HYP_bLR_vs_bHR")

colnames(FoldChanges)<-ComparisonsOfInterest
colnames(Tstats)<-ComparisonsOfInterest

ExtractingDEResults("F4_HYP", FoldChanges, Tstats)

save.image("~/University of Michigan Dropbox/Megan Hagenauer/Wakil_Shared_Animal/bHR_bLR_Microarray_John Stead/Reanalysis_20241209/CTX/F4_bHRbLR_ReanalysisCode_20241209.RData")

sessionInfo()
# R version 4.4.2 (2024-10-31)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.1.1
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/Detroit
# tzcode source: internal
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices utils     datasets 
# [7] methods   base     
# 
# other attached packages:
#   [1] limma_3.62.1                 rae230arnentrezgprobe_25.0.0
# [3] rae230arnentrezgcdf_25.0.0   rae230arnentrezg.db_25.0.0  
# [5] plyr_1.8.9                   org.Rn.eg.db_3.20.0         
# [7] AnnotationDbi_1.68.0         IRanges_2.40.0              
# [9] S4Vectors_0.44.0             affy_1.84.0                 
# [11] Biobase_2.66.0               BiocGenerics_0.52.0         
# [13] BiocManager_1.30.25         
# 
# loaded via a namespace (and not attached):
#   [1] bit_4.5.0.1             preprocessCore_1.68.0  
# [3] jsonlite_1.8.9          compiler_4.4.2         
# [5] crayon_1.5.3            Rcpp_1.0.13-1          
# [7] blob_1.2.4              Biostrings_2.74.0      
# [9] png_0.1-8               statmod_1.5.0          
# [11] fastmap_1.2.0           R6_2.5.1               
# [13] XVector_0.46.0          GenomeInfoDb_1.42.1    
# [15] GenomeInfoDbData_1.2.13 DBI_1.2.3              
# [17] affyio_1.76.0           rlang_1.1.4            
# [19] KEGGREST_1.46.0         cachem_1.1.0           
# [21] bit64_4.5.2             RSQLite_2.3.9          
# [23] memoise_2.0.1           cli_3.6.3              
# [25] zlibbioc_1.52.0         rstudioapi_0.17.1      
# [27] vctrs_0.6.5             httr_1.4.7             
# [29] tools_4.4.2             pkgconfig_2.0.3        
# [31] UCSC.utils_1.2.0   


