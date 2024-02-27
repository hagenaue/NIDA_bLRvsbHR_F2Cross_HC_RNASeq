#F0 HC RNA-Seq Dataset
#08_Graphing the data for previously-identified HR/LR differentially expressed genes
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-10, updated later for a few figures for the paper.



#To make basic graphs of top bHR/bLR genes based on the previous meta-analysis:

#This section still needs to be updated following the removal of the low RIN subject.*********
#Also - do we want to update this to make it instead focus on the top genes from the RNA-Seq Meta-Analysis?


setwd("//n05-corea-cifs.umhs.med.umich.edu/Home1/Users/hebda/My Documents/F0_TMMnorm_20211119/Boxplots")


#Double check that the order of subjects in F0_PCA_wMetaData is the same as in F0_Annotated_lcpm_NoLowHits

cbind(F0_PCA_wMetaData$row_name, colnames(F0_Annotated_lcpm_NoLowHits)[c(1:24)])
#same order - yeah!
# 
# [,1]       [,2]      
# [1,] "sl483888" "SL483888"
# [2,] "sl483889" "SL483889"
# [3,] "sl483890" "SL483890"
# [4,] "sl483891" "SL483891"
# [5,] "sl483892" "SL483892"
# [6,] "sl483893" "SL483893"
# [7,] "sl483894" "SL483894"
# [8,] "sl483895" "SL483895"

#Try making a single graph first to figure out the parameters, e.g.,

#To make it easier to graph, I'm going to make a temporary data frame that contains the data for the gene of interest and the metadata for each subject:
Temp_Data<-data.frame(y=as.matrix(F0_Annotated_lcpm_NoLowHits[which(F0_Annotated_lcpm_NoLowHits$SYMBOL=="C1qa"), c(1:24)])[1,], F0_PCA_wMetaData) 
#The column numbers (currently labeled c(2:25)) should reflect the columns that now contain the gene expression data for each subject (not annotation)
#For some reason, R was treating our data for C1qa as a list and I had to force it into vector format.
#Right now, it is only grabbing the data for the first Ensembl gene with the gene symbol C1qa, if there are more than one.
#The function "which" makes it only grab rows of data with the gene symbol C1qa and not also rows with NA gene symbols

str(Temp_Data)
# 
# 'data.frame':	24 obs. of  19 variables:
#   $ y               : num  6.88 6.44 5.52 5.91 6.52 ...
# $ row_name        : chr  "sl483888" "sl483889" "sl483890" "sl483891" ...
# $ PC1             : num  -0.0118 0.1416 0.2475 0.3507 0.3147 ...
# $ PC2             : num  -0.1134 -0.1208 -0.1083 0.0685 0.0616 ...
# $ PC3             : num  0.1096 0.0784 0.3276 0.2202 -0.0863 ...
# $ PC4             : num  0.0578 0.1502 0.0723 -0.1756 -0.4392 ...
# $ PC5             : num  0.217 -0.1856 -0.3643 0.0563 -0.0302 ...
# $ PC6             : num  -0.076 0.3161 -0.0648 -0.2446 0.0626 ...
# $ Conc            : num  37 51 88.6 81.9 59.8 ...
# $ PF.Reads        : int  24601113 27048765 35810148 31060987 34021183 27895202 27671721 31535526 33899311 35707349 ...
# $ RibosomePerc    : num  0.0243 0.0239 0.0199 0.0123 0.0136 ...
# $ mRNAperc        : num  0.829 0.811 0.809 0.85 0.844 ...
# $ Sex             : chr  "Female" "Male" "Male" "Female" ...
# $ Lineage         : chr  "bLR" "bLR" "bHR" "bHR" ...
# $ Family          : chr  "L01" "L01" "H02" "H02" ...
# $ Age.days.       : int  189 162 163 190 192 165 163 190 165 187 ...
# $ Sex_Factor      : Factor w/ 2 levels "Male","Female": 2 1 1 2 2 1 1 2 1 2 ...
# $ Lineage_AsFactor: Factor w/ 2 levels "bHR","bLR": 2 2 1 1 2 2 1 1 1 1 ...
# $ LibrarySize     : int  16824899 18993612 25568698 23034700 24765660 19970873 19762442 22573091 24594779 23259327 ...

boxplot(y~Sex_Factor*Lineage_AsFactor, data=Temp_Data)
stripchart(y~Sex_Factor*Lineage_AsFactor, data=Temp_Data, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=1, col='black')

#Working on a basic level

#If that works, you can loop it over the top meta-analysis genes using code like this:

setwd("//n05-corea-cifs.umhs.med.umich.edu/Home1/Users/hebda/My Documents/F0_TMMnorm_20211119")

TopMetaAnalysisGenes_Results<-read.csv("AdultMeta_CohDandPval_GeneDatesFixed_FDR10_ForElaine.csv", header=TRUE, stringsAsFactors=FALSE)

str(TopMetaAnalysisGenes_Results)
# 'data.frame':	191 obs. of  28 variables:
# $ GeneSymbol       : chr  "Tmem144" "Asb15" "Kif15" "Pkhd1l1" ...
# $ rawpval          : num  3.04e-08 4.05e-07 4.27e-07 6.13e-07 7.87e-07 1.04e-06 2.03e-06 2.12e-06 2.86e-06 4.27e-06 ...
# $ BH               : num  0.000495 0.002317 0.002317 0.002491 0.002562 ...
# $ estimate_bLRvsbHR: num  3.57 -2.82 2.21 3.17 -2.78 ...
# $ SE               : num  0.644 0.556 0.437 0.636 0.563 ...
# $ CI_lb            : num  -4.83 1.73 -3.07 -4.42 1.68 ...
# $ CI_ub            : num  -2.3 3.91 -1.35 -1.92 3.89 ...
# $ datasets         : int  3 3 4 3 3 4 2 3 4 5 ...

#If this runs properly, you'll end up with something like 100+ graphs outputted, so you may want to set the working directory to a special folder before running this.

setwd("//n05-corea-cifs.umhs.med.umich.edu/Home1/Users/hebda/My Documents/F0_TMMnorm_20220110/Boxplots")


colnames(F0_Annotated_lcpm_NoLowHits)


for(i in c(1:nrow(TopMetaAnalysisGenes_Results))){
  
  if(length(which(F0_Annotated_lcpm_NoLowHits$SYMBOL==TopMetaAnalysisGenes_Results$GeneSymbol[i]))>0){
    
    Temp_Data<-data.frame(y=as.matrix(F0_Annotated_lcpm_NoLowHits[which(F0_Annotated_lcpm_NoLowHits$SYMBOL==TopMetaAnalysisGenes_Results$GeneSymbol[i]), c(1:23)])[1,], F0_PCA_wMetaData) 
    
    pdf(paste("F0_",TopMetaAnalysisGenes_Results$GeneSymbol[i], "_vs_LineageSex.pdf", sep=""), width=6.5, height=6)
    boxplot(y~Sex_Factor*Lineage_AsFactor, data=Temp_Data, ylab=paste(TopMetaAnalysisGenes_Results$GeneSymbol[i], ": Log2 CPM", sep=""), xlab="", col=c("green4","green2","firebrick4", "firebrick2"), cex.lab=1.5, cex.axis=1.2)
    stripchart(y~Sex_Factor*Lineage_AsFactor, data=Temp_Data, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=2, col='black')
    dev.off()
    
    rm(Temp_Data)
    
  }else{}
  
}

#Before trying the loop, I output one iteration by first running:
i<-1
#Then running the code inside the loop


#Notes from glancing at figures:
#1) a certain percentage of the top meta-analysis genes don't show bHR/bLR differences in males or females (50%)
#2) many of the genes that do replicate bHR/bLR differences in males *do not show* that relationship in females
#... so I think there might be a genuine sex difference story here... maybe?

#######################################1/31/2022

setwd("//n05-corea-cifs.umhs.med.umich.edu/Home1/Users/hebda/My Documents/F0_TMMnorm_20220110/Boxplots")

colnames(F0_Annotated_lcpm_NoLowHits)

head (F0_Annotated_lcpm_NoLowHits)



#TopHRLRF2Genes_ForGraphing<-read.csv("TopHRLRF2Genes_ForGraphing.csv", header=TRUE, stringsAsFactors = FALSE)

TopHRLRF2Genes_ForGraphing<-read.csv("Top_HRLR_F2_Genes_forGraphing.csv", header=TRUE, stringsAsFactors = FALSE)

colnames(TopHRLRF2Genes_ForGraphing)[1]<-"ENSEMBL"
colnames(TopHRLRF2Genes_ForGraphing)[2]<-"Symbol"
# 
# 
# he2d(TopHRLRF2Genes_ForGraphing)


for(i in c(1:nrow(TopHRLRF2Genes_ForGraphing))){
  
  if(length(which(F0_Annotated_lcpm_NoLowHits$ENSEMBL==TopHRLRF2Genes_ForGraphing$ENSEMBL[i]))>0){
    
    Temp_Data<-data.frame(y=as.matrix(F0_Annotated_lcpm_NoLowHits[which(F0_Annotated_lcpm_NoLowHits$ENSEMBL==TopHRLRF2Genes_ForGraphing$ENSEMBL[i]), c(1:23)])[1,], F0_PCA_wMetaData) 
    
    pdf(paste("F0_",TopHRLRF2Genes_ForGraphing$Symbol[i], "_", TopHRLRF2Genes_ForGraphing$ENSEMBL[i], "_vs_LineageSex.pdf", sep=""), width=6.5, height=6)
    boxplot(y~Sex_Factor*Lineage_AsFactor, data=Temp_Data, ylab=paste(TopHRLRF2Genes_ForGraphing$Symbol[i], ": Log2 CPM", sep=""), xlab="", col=c("green4","green2","firebrick4", "firebrick2"), cex.lab=1.5, cex.axis=1.2)
    stripchart(y~Sex_Factor*Lineage_AsFactor, data=Temp_Data, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=2, col='black')
    dev.off()
    
    print(TopHRLRF2Genes_ForGraphing$Symbol[i])
    print(Temp_Data$Symbol)
    
    print(summary.lm(lm(y~Sex_Factor*Lineage_AsFactor, data=Temp_Data)))
    
    rm(Temp_Data)
    
  }else{}
  
}

##################################################

