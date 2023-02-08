#F2 HC RNA-Seq Dataset
#16_VolcanoPlots
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-04, updated later for a few figures for the paper.


#Volcano plot code - hasn't been updated yet to illustrate the results for the F2s

getwd()


F2_Limma_Results<-read.delim("F2_Limma_Results_TotalLocsSCore_M4_SexCovariates.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F2_Limma_Results)
#'data.frame':	14056 obs. of  74 variables:
#'
colnames(F2_Limma_Results)[1]<-"ENSEMBL"

colnames(F2_Limma_Results)

tiff("VolcanoPlot_bLRbHR_noInteraction.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F2_Limma_Results, plot(Coef.Lineage_AsFactorbLR, -log10(P.value.Lineage_AsFactorbLR), pch=19, main="Effect of bLR vs. bHR", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-4.5,4.5), ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F2_Limma_Results, abs(Coef.Lineage_AsFactorbLR)>1), points(Coef.Lineage_AsFactorbLR, -log10(P.value.Lineage_AsFactorbLR), pch=19, col="red", cex=0.6))
with(subset(F2_Limma_Results, P.value.adj.Lineage_AsFactorbLR<0.05), points(Coef.Lineage_AsFactorbLR, -log10(P.value.Lineage_AsFactorbLR), pch=19, col="green", cex=0.6))
with(subset(F2_Limma_Results, abs(Coef.Lineage_AsFactorbLR)>1 & P.value.adj.Lineage_AsFactorbLR<0.05), points(Coef.Lineage_AsFactorbLR, -log10(P.value.Lineage_AsFactorbLR), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.05", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()


##########################2/4/2022 Elaine  --this worked but looked funny for all except ST/GT because there are no DE genes

#Volcano plot code - hasn't been updated yet to illustrate the results for the F2s

getwd()


F2_Limma_Results<-read.delim("F2_Limma_Results_TotalLocoScore_SimplerModel_20220110.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F2_Limma_Results)
#'data.frame':	14056 obs. of  40 variables:
#'
colnames(F2_Limma_Results)[1]<-"ENSEMBL"

colnames(F2_Limma_Results)

tiff("VolcanoPlot_F2_TotalLocoScore.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F2_Limma_Results, plot(Coef.Total_LocoScore, -log10(P.value.Total_LocoScore), pch=19, main="LocoScore", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-4.5,4.5), ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F2_Limma_Results, abs(Coef.Total_LocoScore)>1), points(Coef.Total_LocoScore, -log10(P.value.Total_LocoScore), pch=19, col="red", cex=0.6))
with(subset(F2_Limma_Results, P.value.adj.Total_LocoScore<0.05), points(Coef.Total_LocoScore, -log10(P.value.Total_LocoScore), pch=19, col="green", cex=0.6))
with(subset(F2_Limma_Results, abs(Coef.Total_LocoScore)>1 & P.value.adj.Total_LocoScore<0.05), points(Coef.Total_LocoScore, -log10(P.value.Total_LocoScore), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.05", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()




F2_Limma_Results<-read.delim("F2_Limma_Results_EPM_PercOA_SimplerModel_20220121.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F2_Limma_Results)
#'data.frame':	14056 obs. of  40 variables:
#'
colnames(F2_Limma_Results)[1]<-"ENSEMBL"

colnames(F2_Limma_Results)

tiff("VolcanoPlot_F2_EPM_PercOA.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F2_Limma_Results, plot(Coef.EPM_Percent_Time_Open_Arm, -log10(P.value.EPM_Percent_Time_Open_Arm), pch=19, main="EPM %Time in Open Arms", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-4.5,4.5), ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F2_Limma_Results, abs(Coef.EPM_Percent_Time_Open_Arm)>1), points(Coef.EPM_Percent_Time_Open_Arm, -log10(P.value.EPM_Percent_Time_Open_Arm), pch=19, col="red", cex=0.6))
with(subset(F2_Limma_Results, P.value.adj.EPM_Percent_Time_Open_Arm<0.05), points(Coef.EPM_Percent_Time_Open_Arm, -log10(P.value.EPM_Percent_Time_Open_Arm), pch=19, col="green", cex=0.6))
with(subset(F2_Limma_Results, abs(Coef.EPM_Percent_Time_Open_Arm)>1 & P.value.adj.EPM_Percent_Time_Open_Arm<0.05), points(Coef.EPM_Percent_Time_Open_Arm, -log10(P.value.EPM_Percent_Time_Open_Arm), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.05", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()




F2_Limma_Results<-read.delim("F2_Limma_Results_EPM_Time_Immobile_SimplerModel_20220121.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F2_Limma_Results)
#'data.frame':	14056 obs. of  40 variables:
#'
colnames(F2_Limma_Results)[1]<-"ENSEMBL"

colnames(F2_Limma_Results)

tiff("VolcanoPlot_F2_Time_Immobile.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F2_Limma_Results, plot(Coef.EPM_Time_Immobile, -log10(P.value.EPM_Time_Immobile), pch=19, main="EPM Time Immobile", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-4.5,4.5), ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F2_Limma_Results, abs(Coef.EPM_Time_Immobile)>1), points(Coef.EPM_Time_Immobile, -log10(P.value.EPM_Time_Immobile), pch=19, col="red", cex=0.6))
with(subset(F2_Limma_Results, P.value.adj.EPM_Time_Immobile<0.05), points(Coef.EPM_Time_Immobile, -log10(P.value.EPM_Time_Immobile), pch=19, col="green", cex=0.6))
with(subset(F2_Limma_Results, abs(Coef.EPM_Time_Immobile)>1 & P.value.adj.EPM_Time_Immobile<0.05), points(Coef.EPM_Time_Immobile, -log10(P.value.EPM_Time_Immobile), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.05", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()



F2_Limma_Results<-read.delim("F2_Limma_Results_EPM_DistanceTraveled_SimplerModel_20220121.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F2_Limma_Results)
#'data.frame':	14056 obs. of  40 variables:
#'
colnames(F2_Limma_Results)[1]<-"ENSEMBL"

colnames(F2_Limma_Results)

tiff("VolcanoPlot_F2_EPM_DistanceTraveled.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F2_Limma_Results, plot(Coef.EPM_DistanceTraveled, -log10(P.value.EPM_DistanceTraveled), pch=19, main="EPM Distance Traveled", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-4.5,4.5), ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F2_Limma_Results, abs(Coef.EPM_DistanceTraveled)>1), points(Coef.EPM_DistanceTraveled, -log10(P.value.EPM_DistanceTraveled), pch=19, col="red", cex=0.6))
with(subset(F2_Limma_Results, P.value.adj.EPM_DistanceTraveled<0.05), points(Coef.EPM_DistanceTraveled, -log10(P.value.EPM_DistanceTraveled), pch=19, col="green", cex=0.6))
with(subset(F2_Limma_Results, abs(Coef.EPM_DistanceTraveled)>1 & P.value.adj.EPM_DistanceTraveled<0.05), points(Coef.EPM_DistanceTraveled, -log10(P.value.EPM_DistanceTraveled), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.05", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()


#######This gives a volcano plot of ST vs No ST/GT experience
F2_Limma_Results<-read.delim("F2_Limma_Results_LearningClassification_SimplerModel_20220121.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F2_Limma_Results)
#'data.frame':	14056 obs. of  44 variables:
#'
colnames(F2_Limma_Results)[1]<-"ENSEMBL"

colnames(F2_Limma_Results)

tiff("VolcanoPlot_F2_LearningClassification_AsFactor_GivesSTvsNoSTGTexper.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F2_Limma_Results, plot(Coef.LearningClassification_AsFactor, -log10(P.value.LearningClassification_AsFactor), pch=19, main="ST/GT", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-4.5,4.5), ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F2_Limma_Results, abs(Coef.LearningClassification_AsFactor)>1), points(Coef.LearningClassification_AsFactor, -log10(P.value.LearningClassification_AsFactor), pch=19, col="red", cex=0.6))
with(subset(F2_Limma_Results, P.value.adj.LearningClassification_AsFactor<0.05), points(Coef.LearningClassification_AsFactor, -log10(P.value.LearningClassification_AsFactor), pch=19, col="green", cex=0.6))
with(subset(F2_Limma_Results, abs(Coef.LearningClassification_AsFactor)>1 & P.value.adj.LearningClassification_AsFactor<0.05), points(Coef.LearningClassification_AsFactor, -log10(P.value.LearningClassification_AsFactor), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.05", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()

str(F2_PCA_wMetaData)



##########################2/10/2022 Elaine  --changed the xlim below because the coefficients are small because the variables are continuous and the units are tiny


getwd()


F2_Limma_Results<-read.delim("F2_Limma_Results_TotalLocoScore_SimplerModel_20220110.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F2_Limma_Results)
#'data.frame':	14056 obs. of  40 variables:
#'
colnames(F2_Limma_Results)[1]<-"ENSEMBL"

colnames(F2_Limma_Results)

tiff("VolcanoPlot_F2_TotalLocoScore_changeXLIM.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F2_Limma_Results, plot(Coef.Total_LocoScore, -log10(P.value.Total_LocoScore), pch=19, main="LocoScore", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-0.1,0.1), ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F2_Limma_Results, abs(Coef.Total_LocoScore)>1), points(Coef.Total_LocoScore, -log10(P.value.Total_LocoScore), pch=19, col="red", cex=0.6))
with(subset(F2_Limma_Results, P.value.adj.Total_LocoScore<0.05), points(Coef.Total_LocoScore, -log10(P.value.Total_LocoScore), pch=19, col="green", cex=0.6))
with(subset(F2_Limma_Results, abs(Coef.Total_LocoScore)>1 & P.value.adj.Total_LocoScore<0.05), points(Coef.Total_LocoScore, -log10(P.value.Total_LocoScore), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.05", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()




F2_Limma_Results<-read.delim("F2_Limma_Results_EPM_PercOA_SimplerModel_20220121.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F2_Limma_Results)
#'data.frame':	14056 obs. of  40 variables:
#'
colnames(F2_Limma_Results)[1]<-"ENSEMBL"

colnames(F2_Limma_Results)

tiff("VolcanoPlot_F2_EPM_PercOA_ChangeXLIM.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F2_Limma_Results, plot(Coef.EPM_Percent_Time_Open_Arm, -log10(P.value.EPM_Percent_Time_Open_Arm), pch=19, main="EPM %Time in Open Arms", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-0.1,0.1), ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F2_Limma_Results, abs(Coef.EPM_Percent_Time_Open_Arm)>1), points(Coef.EPM_Percent_Time_Open_Arm, -log10(P.value.EPM_Percent_Time_Open_Arm), pch=19, col="red", cex=0.6))
with(subset(F2_Limma_Results, P.value.adj.EPM_Percent_Time_Open_Arm<0.05), points(Coef.EPM_Percent_Time_Open_Arm, -log10(P.value.EPM_Percent_Time_Open_Arm), pch=19, col="green", cex=0.6))
with(subset(F2_Limma_Results, abs(Coef.EPM_Percent_Time_Open_Arm)>1 & P.value.adj.EPM_Percent_Time_Open_Arm<0.05), points(Coef.EPM_Percent_Time_Open_Arm, -log10(P.value.EPM_Percent_Time_Open_Arm), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.05", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()




F2_Limma_Results<-read.delim("F2_Limma_Results_EPM_Time_Immobile_SimplerModel_20220121.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F2_Limma_Results)
#'data.frame':	14056 obs. of  40 variables:
#'
colnames(F2_Limma_Results)[1]<-"ENSEMBL"

colnames(F2_Limma_Results)

tiff("VolcanoPlot_F2_Time_Immobile_ChangeXLIM.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F2_Limma_Results, plot(Coef.EPM_Time_Immobile, -log10(P.value.EPM_Time_Immobile), pch=19, main="EPM Time Immobile", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-0.1,0.1), ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F2_Limma_Results, abs(Coef.EPM_Time_Immobile)>1), points(Coef.EPM_Time_Immobile, -log10(P.value.EPM_Time_Immobile), pch=19, col="red", cex=0.6))
with(subset(F2_Limma_Results, P.value.adj.EPM_Time_Immobile<0.05), points(Coef.EPM_Time_Immobile, -log10(P.value.EPM_Time_Immobile), pch=19, col="green", cex=0.6))
with(subset(F2_Limma_Results, abs(Coef.EPM_Time_Immobile)>1 & P.value.adj.EPM_Time_Immobile<0.05), points(Coef.EPM_Time_Immobile, -log10(P.value.EPM_Time_Immobile), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.05", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()



F2_Limma_Results<-read.delim("F2_Limma_Results_EPM_DistanceTraveled_SimplerModel_20220121.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F2_Limma_Results)
#'data.frame':	14056 obs. of  40 variables:
#'
colnames(F2_Limma_Results)[1]<-"ENSEMBL"

colnames(F2_Limma_Results)

tiff("VolcanoPlot_F2_EPM_DistanceTraveled_ChangeXLIM.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F2_Limma_Results, plot(Coef.EPM_DistanceTraveled, -log10(P.value.EPM_DistanceTraveled), pch=19, main="EPM Distance Traveled", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-0.1,0.1), ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F2_Limma_Results, abs(Coef.EPM_DistanceTraveled)>1), points(Coef.EPM_DistanceTraveled, -log10(P.value.EPM_DistanceTraveled), pch=19, col="red", cex=0.6))
with(subset(F2_Limma_Results, P.value.adj.EPM_DistanceTraveled<0.05), points(Coef.EPM_DistanceTraveled, -log10(P.value.EPM_DistanceTraveled), pch=19, col="green", cex=0.6))
with(subset(F2_Limma_Results, abs(Coef.EPM_DistanceTraveled)>1 & P.value.adj.EPM_DistanceTraveled<0.05), points(Coef.EPM_DistanceTraveled, -log10(P.value.EPM_DistanceTraveled), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.05", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()



F2_Limma_Results<-read.delim("F2_Limma_Results_LearningClassification_SimplerModel_20220121.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F2_Limma_Results)
#'data.frame':	14056 obs. of  44 variables:
#'
colnames(F2_Limma_Results)[1]<-"ENSEMBL"

colnames(F2_Limma_Results)

tiff("VolcanoPlot_F2_LearningClassification_AsFactor_GivesSTvsNoSTGTexper_ChangeXLIM.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F2_Limma_Results, plot(Coef.LearningClassification_AsFactor, -log10(P.value.LearningClassification_AsFactor), pch=19, main="ST/GT", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-1.5,1.5), ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F2_Limma_Results, abs(Coef.LearningClassification_AsFactor)>1), points(Coef.LearningClassification_AsFactor, -log10(P.value.LearningClassification_AsFactor), pch=19, col="red", cex=0.6))
with(subset(F2_Limma_Results, P.value.adj.LearningClassification_AsFactor<0.05), points(Coef.LearningClassification_AsFactor, -log10(P.value.LearningClassification_AsFactor), pch=19, col="green", cex=0.6))
with(subset(F2_Limma_Results, abs(Coef.LearningClassification_AsFactor)>1 & P.value.adj.LearningClassification_AsFactor<0.05), points(Coef.LearningClassification_AsFactor, -log10(P.value.LearningClassification_AsFactor), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.05", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()

##############################


##########################2/10/2022 Elaine  --changed the xlim below because the coefficients are small because the variables are continuous and the units are tiny---made XLIM even smaller


getwd()


F2_Limma_Results<-read.delim("F2_Limma_Results_TotalLocoScore_SimplerModel_20220110.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F2_Limma_Results)
#'data.frame':	14056 obs. of  40 variables:
#'
colnames(F2_Limma_Results)[1]<-"ENSEMBL"

colnames(F2_Limma_Results)

tiff("VolcanoPlot_F2_TotalLocoScore_changeXLIM2.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F2_Limma_Results, plot(Coef.Total_LocoScore, -log10(P.value.Total_LocoScore), pch=19, main="LocoScore", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-0.002,0.002), ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F2_Limma_Results, abs(Coef.Total_LocoScore)>1), points(Coef.Total_LocoScore, -log10(P.value.Total_LocoScore), pch=19, col="red", cex=0.6))
with(subset(F2_Limma_Results, P.value.adj.Total_LocoScore<0.05), points(Coef.Total_LocoScore, -log10(P.value.Total_LocoScore), pch=19, col="green", cex=0.6))
with(subset(F2_Limma_Results, abs(Coef.Total_LocoScore)>1 & P.value.adj.Total_LocoScore<0.05), points(Coef.Total_LocoScore, -log10(P.value.Total_LocoScore), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.05", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()




F2_Limma_Results<-read.delim("F2_Limma_Results_EPM_PercOA_SimplerModel_20220121.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F2_Limma_Results)
#'data.frame':	14056 obs. of  40 variables:
#'
colnames(F2_Limma_Results)[1]<-"ENSEMBL"

colnames(F2_Limma_Results)

tiff("VolcanoPlot_F2_EPM_PercOA_ChangeXLIM2.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F2_Limma_Results, plot(Coef.EPM_Percent_Time_Open_Arm, -log10(P.value.EPM_Percent_Time_Open_Arm), pch=19, main="EPM %Time in Open Arms", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-0.05,0.05), ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F2_Limma_Results, abs(Coef.EPM_Percent_Time_Open_Arm)>1), points(Coef.EPM_Percent_Time_Open_Arm, -log10(P.value.EPM_Percent_Time_Open_Arm), pch=19, col="red", cex=0.6))
with(subset(F2_Limma_Results, P.value.adj.EPM_Percent_Time_Open_Arm<0.05), points(Coef.EPM_Percent_Time_Open_Arm, -log10(P.value.EPM_Percent_Time_Open_Arm), pch=19, col="green", cex=0.6))
with(subset(F2_Limma_Results, abs(Coef.EPM_Percent_Time_Open_Arm)>1 & P.value.adj.EPM_Percent_Time_Open_Arm<0.05), points(Coef.EPM_Percent_Time_Open_Arm, -log10(P.value.EPM_Percent_Time_Open_Arm), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.05", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()




F2_Limma_Results<-read.delim("F2_Limma_Results_EPM_Time_Immobile_SimplerModel_20220121.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F2_Limma_Results)
#'data.frame':	14056 obs. of  40 variables:
#'
colnames(F2_Limma_Results)[1]<-"ENSEMBL"

colnames(F2_Limma_Results)

tiff("VolcanoPlot_F2_Time_Immobile_ChangeXLIM2.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F2_Limma_Results, plot(Coef.EPM_Time_Immobile, -log10(P.value.EPM_Time_Immobile), pch=19, main="EPM Time Immobile", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-0.02,0.02), ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F2_Limma_Results, abs(Coef.EPM_Time_Immobile)>1), points(Coef.EPM_Time_Immobile, -log10(P.value.EPM_Time_Immobile), pch=19, col="red", cex=0.6))
with(subset(F2_Limma_Results, P.value.adj.EPM_Time_Immobile<0.05), points(Coef.EPM_Time_Immobile, -log10(P.value.EPM_Time_Immobile), pch=19, col="green", cex=0.6))
with(subset(F2_Limma_Results, abs(Coef.EPM_Time_Immobile)>1 & P.value.adj.EPM_Time_Immobile<0.05), points(Coef.EPM_Time_Immobile, -log10(P.value.EPM_Time_Immobile), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.05", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()



F2_Limma_Results<-read.delim("F2_Limma_Results_EPM_DistanceTraveled_SimplerModel_20220121.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F2_Limma_Results)
#'data.frame':	14056 obs. of  40 variables:
#'
colnames(F2_Limma_Results)[1]<-"ENSEMBL"

colnames(F2_Limma_Results)

tiff("VolcanoPlot_F2_EPM_DistanceTraveled_ChangeXLIM2.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F2_Limma_Results, plot(Coef.EPM_DistanceTraveled, -log10(P.value.EPM_DistanceTraveled), pch=19, main="EPM Distance Traveled", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-0.002,0.002), ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F2_Limma_Results, abs(Coef.EPM_DistanceTraveled)>1), points(Coef.EPM_DistanceTraveled, -log10(P.value.EPM_DistanceTraveled), pch=19, col="red", cex=0.6))
with(subset(F2_Limma_Results, P.value.adj.EPM_DistanceTraveled<0.05), points(Coef.EPM_DistanceTraveled, -log10(P.value.EPM_DistanceTraveled), pch=19, col="green", cex=0.6))
with(subset(F2_Limma_Results, abs(Coef.EPM_DistanceTraveled)>1 & P.value.adj.EPM_DistanceTraveled<0.05), points(Coef.EPM_DistanceTraveled, -log10(P.value.EPM_DistanceTraveled), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.05", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()



##############################2/11/2022 Re-created ST/GT volcano plot to reflect ST vs GT not ST vs No ST/GT experience



F2_Limma_Results<-read.delim("F2_Limma_Results_LearningClassification_SimplerModel_20220121.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F2_Limma_Results)
#'data.frame':	14056 obs. of  44 variables:
#'
colnames(F2_Limma_Results)[1]<-"ENSEMBL"

colnames(F2_Limma_Results)

tiff("VolcanoPlot_F2_LearningClassification_AsFactor_STvsGT_ChangeXLIM.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F2_Limma_Results, plot(Coef.LearningClassification_AsFactorGT, -log10(P.value.LearningClassification_AsFactorGT), pch=19, main="ST vs GT", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-1,1), ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F2_Limma_Results, abs(Coef.LearningClassification_AsFactorGT)>1), points(Coef.LearningClassification_AsFactorGT, -log10(P.value.LearningClassification_AsFactorGT), pch=19, col="red", cex=0.6))
with(subset(F2_Limma_Results, P.value.adj.LearningClassification_AsFactorGT<0.05), points(Coef.LearningClassification_AsFactorGT, -log10(P.value.LearningClassification_AsFactorGT), pch=19, col="green", cex=0.6))
with(subset(F2_Limma_Results, abs(Coef.LearningClassification_AsFactorGT)>1 & P.value.adj.LearningClassification_AsFactorGT<0.05), points(Coef.LearningClassification_AsFactorGT, -log10(P.value.LearningClassification_AsFactorGT), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.05", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()


