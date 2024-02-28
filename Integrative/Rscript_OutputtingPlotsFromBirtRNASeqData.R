#I updated this to grab example plots to compare to our NIDA U01 RNA-Seq study
#Megan Hagenauer
#last updated October 2023

library(car)
library(plyr)


setwd("~/Phenotype Project/ibirt")

######################

#Reading in Peter's data:

# PeterHR_LRdata <- read.csv("file:///C:/Users/Izzy/Documents/Phenotype Project/PeterBlandino_BasalData/Original Files/LowerCase_Akil5_geneexp_filteredandnormalized_log2fpm.csv")
# colnames(PeterHR_LRdata)

#Megan's code for reading in Peter's data:
#setwd("~/Documents/Microarray Gen/HRLR/HRLR_Studies/NIH_basalHRLR_RNAseq")
setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/HRLR_Studies/NIH_basalHRLR_RNAseq")

PeterHR_LRdata <- read.csv("LowerCase_Akil5_geneexp_filteredandnormalized_log2fpm.csv", header=T, stringsAsFactors = F)
colnames(PeterHR_LRdata)
# [1] "Gene"   "X10HR"  "X22HR"  "X26HR"  "X2HR"   "X34HR"  "X38HR"  "X149IR" "X158IR" "X213IR" "X248IR" "X261IR" "X293IR" "X16LR"  "X20LR" [16] "X32LR"  "X44LR"  "X48LR"  "X8LR"  
head(PeterHR_LRdata)
# Gene     X10HR      X22HR     X26HR      X2HR      X34HR     X38HR    X149IR    X158IR     X213IR    X248IR    X261IR    X293IR
# 1     A1bg -3.316754 -1.8308858 -6.796308 -1.877379 -2.3757866 -2.811070 -4.192394 -2.010690 -0.9914813 -3.955577 -3.383915 -1.454892
# 2      A2m  4.438876  4.6586755  4.076367  4.447728  4.6967480  4.708129  4.064050  4.814896  4.9379770  4.054531  4.887548  4.317394
# 3  A3galt2  1.314221  0.9163481  1.092436  1.502763  1.0660509  1.453633  1.193037  1.627911  1.6632514  1.716848  1.347027  1.294457
# 4 Aa926063 -2.579788 -2.5227635 -1.510905 -3.157487 -0.7908241 -3.659067 -2.607432 -2.531522 -0.9914813 -1.885188 -2.466377 -3.142948
# 5     Aaas  4.233566  4.1288402  4.180257  4.147597  4.1404387  4.234538  4.303461  4.102895  4.2221812  4.193154  4.030183  4.118625
# 6     Aacs  4.551555  4.7901195  5.218759  4.437104  4.9717609  4.885511  5.229250  4.815666  4.7927113  4.476547  4.881310  5.100226
## No need to check data quality as this data has already been processed
#Moving onto normalization - keeping IRs for behavior analyses

######
## Reading in Peter's sample descriptions/behavior scores

# sampleInfo <- read.csv("file:///C:/Users/Izzy/Documents/Phenotype Project/PeterBlandino_BasalData/Original Files/PeterHRLR_Akil5_sample_description.csv")
# colnames(sampleInfo)

#Megan's code for reading in behavioral data - note, I made sure that the behavioral data was in the same order as the gene expression data.
#setwd("~/Documents/Microarray Gen/HRLR/HRLR_Studies/NIH_basalHRLR_RNAseq/akil-5")
setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/HRLR_Studies/NIH_basalHRLR_RNAseq/akil-5")

sampleInfo <- read.csv("Akil5_sample_description.csv", header=T, stringsAsFactors = F)
colnames(sampleInfo)
# [1] "Rat_ID"             "FamilyID"           "Generation"         "Lineage"            "RatID"              "TotalLocoScore"    
# [7] "Boli"               "DistanceTraveled"   "TimeImmobile"       "PercentTimeOpenArm" "Sex"    

#################

#Reading in Cigdem's data:

#setwd("~/Documents/Microarray Gen/HRLR/HRLR_Studies/Cigdem_HRLR_RNAseq")
setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/HRLR_Studies/Cigdem_HRLR_RNAseq")

Cigdem_HRLRdata<-read.csv("CigdemNew_Vehicle.csv", header=T, stringsAsFactors = F)
colnames(Cigdem_HRLRdata)
#[1] "X"         "HR_CONT_1" "HR_CONT_2" "HR_CONT_3" "HR_CONT_4" "HR_CONT_5" "LR_CONT_1" "LR_CONT_2" "LR_CONT_3" "LR_CONT_4" "LR_CONT_5"
head(Cigdem_HRLRdata)
# X HR_CONT_1 HR_CONT_2 HR_CONT_3 HR_CONT_4 HR_CONT_5 LR_CONT_1 LR_CONT_2 LR_CONT_3 LR_CONT_4 LR_CONT_5
# 1 ENSRNOG00000029897  1.984111  1.611828  1.911635  1.879735  1.978242  1.592064  1.872764  2.052356  2.135435  2.247886
# 2 ENSRNOG00000014303  8.080741  8.114191  8.115799  8.109010  8.074836  8.069525  8.042232  8.146873  8.127431  8.047361
# 3 ENSRNOG00000014330  6.626785  6.711768  6.653282  6.678040  6.657021  6.686488  6.656599  6.742801  6.714538  6.740110
# 4 ENSRNOG00000049505  2.637977  2.645080  2.541030  2.665524  2.627412  2.746949  2.674064  2.697691  2.726000  2.758554
# 5 ENSRNOG00000014916  5.624220  5.667279  5.690249  5.689228  5.776825  5.719697  5.760098  5.726541  5.717119  5.823067
# 6 ENSRNOG00000014996  3.044383  3.196790  3.199910  3.232999  3.202155  3.065995  3.249105  3.086303  2.972160  3.228257

Cigdem_Annotation<-read.csv("GeneAnnotation_forNewData.csv", header=T, stringsAsFactors = F)
colnames(Cigdem_Annotation)
#[1] "ENSEMBLID"  "GeneSymbol" "GeneType" 
head(Cigdem_Annotation)
# ENSEMBLID     GeneSymbol       GeneType
# 1 ENSRNOG00000000007           Gad1 protein_coding
# 2 ENSRNOG00000000008           Alx4 protein_coding
# 3 ENSRNOG00000000010          Cbln1 protein_coding
# 4 ENSRNOG00000000012          Tcf15 protein_coding

Cigdem_sampleInfo <- read.csv("sampledescription_justVeh.csv", header=T, stringsAsFactors = F)
colnames(Cigdem_sampleInfo)
#[1] "Sample_ID"   "Core_ID"     "Sample_Desc" "Rat_Type"    "Treatment"   "SI_Score"    "Group"  
head(Cigdem_sampleInfo)
# Sample_ID Core_ID Sample_Desc Rat_Type Treatment SI_Score   Group
# 1        53   53816   HR_CONT_1       HR      CONT    14.00 HR_CONT
# 2        54   53817   HR_CONT_2       HR      CONT    16.67 HR_CONT
# 3       114   53825   HR_CONT_3       HR      CONT    17.67 HR_CONT
# 4       115   53826   HR_CONT_4       HR      CONT    24.33 HR_CONT
# 5       132   53829   HR_CONT_5       HR      CONT    24.00 HR_CONT
# 6        11   53807   LR_CONT_1       LR      CONT    12.33 LR_CONT
#Note: I already made sure these were in the same order as the gene expression data.

##Making some bHR/bLR plots that match the ones in our NIDA U01 


#setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/Matching_F43_F37_Plots")
setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/Matching_F43_F37_Plots")

GenesToPlot<-c("Ghdc", "Ifit1", "Sp3", "Plekhb1", "Asb15", "Fzd6", "Pex11a", "Fanci", "Vps9d1", "Idh1", "Rarres2", "Lsr", "Ist1", "Unc45a", "Gimap5", "Wdr93", "Spg7", "Afg3l1", "Ucp2", "Mfge8", "Ttc30a1", "C2cd3", "Lipt2", "Cndp1", "Zfp551", "Tmco5a","Maff", "Ankrd54", "Oard1", "Loxhd1", "Sun2", "Fcrl2", "Ndufaf6", "Tmem144", "Hba-a3", "Maff", "Mcee", "Pdxp", "Cmc1")

for(i in c(1:length(GenesToPlot))){
  
  if(length(which(Cigdem_Annotation$GeneSymbol==GenesToPlot[i]))>0){
    
    Temp_Data<-data.frame(y=as.matrix(Cigdem_HRLRdata[which(Cigdem_Annotation$GeneSymbol==GenesToPlot[i]), c(2:11)])[1,], Cigdem_sampleInfo) 
    
    pdf(paste("F43_",GenesToPlot[i], "_vs_Lineage.pdf", sep=""), width=4.5, height=6)
    boxplot(y~Group, data=Temp_Data, ylab=paste(GenesToPlot[i], ": Log2 FPM", sep=""), xlab="", col=c("green4","firebrick4"), cex.lab=1.5, cex.axis=1.2)
    stripchart(y~Group, data=Temp_Data, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=2, col='black')
    dev.off()
    
    rm(Temp_Data)
    
  }else{}
  
  if(length(which(PeterHR_LRdata$Gene==GenesToPlot[i]))>0){
    
    Temp_Data<-data.frame(y=as.matrix(PeterHR_LRdata[which(PeterHR_LRdata$Gene==GenesToPlot[i]), c(2:19)])[1,], sampleInfo) 
    
    pdf(paste("F37_",GenesToPlot[i], "_vs_Lineage.pdf", sep=""), width=5.0, height=6)
    boxplot(y~Lineage, data=Temp_Data, ylab=paste(GenesToPlot[i], ": Log2 FPM", sep=""), xlab="", col=c("green4","goldenrod2","firebrick4"), cex.lab=1.5, cex.axis=1.2)
    stripchart(y~Lineage, data=Temp_Data, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=2, col='black')
    dev.off()
    
    rm(Temp_Data)
    
  }else{}
  
}


####NOTE - COME BACK TO THIS - A WHOLE BUNCH OF THE F43 FIGS ARE OUTPUTTING INCORRECTLY
#... but maybe we won't use them in the paper anyway - easier to just talk about the F0s
