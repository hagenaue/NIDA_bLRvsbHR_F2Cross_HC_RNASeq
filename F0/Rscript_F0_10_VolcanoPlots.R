#F0 HC RNA-Seq Dataset
#10_Volcano Plots
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-10, updated later for a few figures for the paper.



#Creating volcano plots:

getwd()


F0_Limma_Results<-read.delim("F0_Limma_Results_Sex_Lineage_rRNA_Intergenic.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F0_Limma_Results)
#'data.frame':	13814 obs. of  24 variables:
#'
colnames(F0_Limma_Results)[1]<-"ENSEMBL"

colnames(F0_Limma_Results)

tiff("VolcanoPlot_bLRbHR_noInteraction.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F0_Limma_Results, plot(Coef.Lineage_AsFactorbLR, -log10(P.value.Lineage_AsFactorbLR), pch=19, main="Effect of bLR vs. bHR", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-4.5,4.5), ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F0_Limma_Results, abs(Coef.Lineage_AsFactorbLR)>1), points(Coef.Lineage_AsFactorbLR, -log10(P.value.Lineage_AsFactorbLR), pch=19, col="red", cex=0.6))
with(subset(F0_Limma_Results, P.value.adj.Lineage_AsFactorbLR<0.05), points(Coef.Lineage_AsFactorbLR, -log10(P.value.Lineage_AsFactorbLR), pch=19, col="green", cex=0.6))
with(subset(F0_Limma_Results, abs(Coef.Lineage_AsFactorbLR)>1 & P.value.adj.Lineage_AsFactorbLR<0.05), points(Coef.Lineage_AsFactorbLR, -log10(P.value.Lineage_AsFactorbLR), pch=19, col="yellow", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.05", "both"), col=c("red", "green", "yellow"), pch=19, cex=1)
dev.off()

##########################2/11/2022 Change "yellow" to another color that can be seen better on the screen

F0_Limma_Results<-read.delim("F0_Limma_Results_Sex_Lineage_rRNA_Intergenic.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F0_Limma_Results)
#'data.frame':	13814 obs. of  24 variables:
#'
colnames(F0_Limma_Results)[1]<-"ENSEMBL"

colnames(F0_Limma_Results)

tiff("VolcanoPlot_bLRbHR_noInteraction_wGold.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F0_Limma_Results, plot(Coef.Lineage_AsFactorbLR, -log10(P.value.Lineage_AsFactorbLR), pch=19, main="Effect of bLR vs. bHR", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-4.5,4.5), ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F0_Limma_Results, abs(Coef.Lineage_AsFactorbLR)>1), points(Coef.Lineage_AsFactorbLR, -log10(P.value.Lineage_AsFactorbLR), pch=19, col="red", cex=0.6))
with(subset(F0_Limma_Results, P.value.adj.Lineage_AsFactorbLR<0.05), points(Coef.Lineage_AsFactorbLR, -log10(P.value.Lineage_AsFactorbLR), pch=19, col="green", cex=0.6))
with(subset(F0_Limma_Results, abs(Coef.Lineage_AsFactorbLR)>1 & P.value.adj.Lineage_AsFactorbLR<0.05), points(Coef.Lineage_AsFactorbLR, -log10(P.value.Lineage_AsFactorbLR), pch=19, col="gold", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.05", "both"), col=c("red", "green", "gold"), pch=19, cex=1)
dev.off()

##########################2/28/2022 Create volcano plot for sex

F0_Limma_Results<-read.delim("F0_Limma_Results_Sex_Lineage_rRNA_Intergenic.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F0_Limma_Results)
#'data.frame':	13814 obs. of  24 variables:
#'
colnames(F0_Limma_Results)[1]<-"ENSEMBL"

colnames(F0_Limma_Results)

tiff("VolcanoPlot_Sex_noInteraction_wGold.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F0_Limma_Results, plot(Coef.Sex_FactorFemale, -log10(P.value.Sex_FactorFemale), pch=19, main="Effect of Female vs. Male", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-4.5,4.5), ylim=c(0,8), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F0_Limma_Results, abs(Coef.Sex_FactorFemale)>1), points(Coef.Sex_FactorFemale, -log10(P.value.Sex_FactorFemale), pch=19, col="red", cex=0.6))
with(subset(F0_Limma_Results, P.value.adj.Sex_FactorFemale<0.05), points(Coef.Sex_FactorFemale, -log10(P.value.Sex_FactorFemale), pch=19, col="green", cex=0.6))
with(subset(F0_Limma_Results, abs(Coef.Sex_FactorFemale)>1 & P.value.adj.Sex_FactorFemale<0.05), points(Coef.Sex_FactorFemale, -log10(P.value.Sex_FactorFemale), pch=19, col="gold", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.05", "both"), col=c("red", "green", "gold"), pch=19, cex=1)
dev.off()

##########################3/28/2022 Create volcano plot for sex--adjusting x axis because of few genes that are GREATLY DE between males and females

F0_Limma_Results<-read.delim("F0_Limma_Results_Sex_Lineage_rRNA_Intergenic.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F0_Limma_Results)

######Did I not change the following line above on 2/28/2022?  Now there are 13,786 obs instead of 13,814#########################
#'data.frame':	13786 obs. of  24 variables:
#'
colnames(F0_Limma_Results)[1]<-"ENSEMBL"

colnames(F0_Limma_Results)

tiff("VolcanoPlot_Sex_noInteraction_wGold_AdjAxes.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F0_Limma_Results, plot(Coef.Sex_FactorFemale, -log10(P.value.Sex_FactorFemale), pch=19, main="Effect of Female vs. Male", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-15,15), ylim=c(0,30), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F0_Limma_Results, abs(Coef.Sex_FactorFemale)>1), points(Coef.Sex_FactorFemale, -log10(P.value.Sex_FactorFemale), pch=19, col="red", cex=0.6))
with(subset(F0_Limma_Results, P.value.adj.Sex_FactorFemale<0.05), points(Coef.Sex_FactorFemale, -log10(P.value.Sex_FactorFemale), pch=19, col="green", cex=0.6))
with(subset(F0_Limma_Results, abs(Coef.Sex_FactorFemale)>1 & P.value.adj.Sex_FactorFemale<0.05), points(Coef.Sex_FactorFemale, -log10(P.value.Sex_FactorFemale), pch=19, col="gold", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.05", "both"), col=c("red", "green", "gold"), pch=19, cex=1)
dev.off()



##########################4/14/2022 Change FDR to <0.1##########################Change "yellow" to another color that can be seen better on the screen#################

F0_Limma_Results<-read.delim("F0_Limma_Results_Sex_Lineage_rRNA_Intergenic.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F0_Limma_Results)
#'data.frame':	13786 obs. of  24 variables:
#'
colnames(F0_Limma_Results)[1]<-"ENSEMBL"

colnames(F0_Limma_Results)

tiff("VolcanoPlot_bLRbHR_noInteraction_wGold_FDR_0point1.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F0_Limma_Results, plot(Coef.Lineage_AsFactorbLR, -log10(P.value.Lineage_AsFactorbLR), pch=19, main="Effect of bLR vs. bHR", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-6,6), ylim=c(0,12), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F0_Limma_Results, abs(Coef.Lineage_AsFactorbLR)>1), points(Coef.Lineage_AsFactorbLR, -log10(P.value.Lineage_AsFactorbLR), pch=19, col="red", cex=0.6))
with(subset(F0_Limma_Results, P.value.adj.Lineage_AsFactorbLR<0.1), points(Coef.Lineage_AsFactorbLR, -log10(P.value.Lineage_AsFactorbLR), pch=19, col="green", cex=0.6))
with(subset(F0_Limma_Results, abs(Coef.Lineage_AsFactorbLR)>1 & P.value.adj.Lineage_AsFactorbLR<0.1), points(Coef.Lineage_AsFactorbLR, -log10(P.value.Lineage_AsFactorbLR), pch=19, col="gold", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.1", "both"), col=c("red", "green", "gold"), pch=19, cex=1)
dev.off()


F0_Limma_Results<-read.delim("F0_Limma_Results_Sex_Lineage_rRNA_Intergenic.txt", sep="\t", header=TRUE, stringsAsFactors = FALSE)
str(F0_Limma_Results)
#data.frame':	13786 obs. of  24 variables:
#'
colnames(F0_Limma_Results)[1]<-"ENSEMBL"

colnames(F0_Limma_Results)

tiff("VolcanoPlot_Sex_noInteraction_wGold_FDR_0point1.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F0_Limma_Results, plot(Coef.Sex_FactorFemale, -log10(P.value.Sex_FactorFemale), pch=19, main="Effect of Female vs. Male", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-2.5,2.5), ylim=c(0,8), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F0_Limma_Results, abs(Coef.Sex_FactorFemale)>1), points(Coef.Sex_FactorFemale, -log10(P.value.Sex_FactorFemale), pch=19, col="red", cex=0.6))
with(subset(F0_Limma_Results, P.value.adj.Sex_FactorFemale<0.1), points(Coef.Sex_FactorFemale, -log10(P.value.Sex_FactorFemale), pch=19, col="green", cex=0.6))
with(subset(F0_Limma_Results, abs(Coef.Sex_FactorFemale)>1 & P.value.adj.Sex_FactorFemale<0.1), points(Coef.Sex_FactorFemale, -log10(P.value.Sex_FactorFemale), pch=19, col="gold", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.1", "both"), col=c("red", "green", "gold"), pch=19, cex=1)
dev.off()


tiff("VolcanoPlot_Sex_noInteraction_wGold_FDR_0point1_axesMatch_HRLR.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F0_Limma_Results, plot(Coef.Sex_FactorFemale, -log10(P.value.Sex_FactorFemale), pch=19, main="Effect of Female vs. Male", cex.lab=1.8, cex.main=2, cex=0.6, xlim=c(-6,6), ylim=c(0,12), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F0_Limma_Results, abs(Coef.Sex_FactorFemale)>1), points(Coef.Sex_FactorFemale, -log10(P.value.Sex_FactorFemale), pch=19, col="red", cex=0.6))
with(subset(F0_Limma_Results, P.value.adj.Sex_FactorFemale<0.1), points(Coef.Sex_FactorFemale, -log10(P.value.Sex_FactorFemale), pch=19, col="green", cex=0.6))
with(subset(F0_Limma_Results, abs(Coef.Sex_FactorFemale)>1 & P.value.adj.Sex_FactorFemale<0.1), points(Coef.Sex_FactorFemale, -log10(P.value.Sex_FactorFemale), pch=19, col="gold", cex=0.6))
#legend(-1.5, 7.6, legend=c("estimate > 1", "FDR<0.1", "both"), col=c("red", "green", "gold"), pch=19, cex=1)
dev.off()






