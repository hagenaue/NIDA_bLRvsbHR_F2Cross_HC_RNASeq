#Comparing RNA-Seq Results across datasets: F0 vs. F2 vs. bHR/bLR Meta-analysis
#03_Running a candidate gene analysis for the F2 data using bHR/bLR differentially expressed genes
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-27, updated later for a few figures for the paper.


library(multtest)


tempPvalAdj<-mt.rawp2adjp(F0_Meta_F2_4Models_LineageGenes$M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore, proc="BH")
F0_Meta_F2_4Models_LineageGenes$M9_FDR.Total_LocoScore_OnlyHRLRGenes<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
rm(tempPvalAdj)

tempPvalAdj<-mt.rawp2adjp(F0_Meta_F2_4Models_LineageGenes$P.value.EPM_Percent_Time_Open_Arm, proc="BH")
F0_Meta_F2_4Models_LineageGenes$FDR.EPM_Percent_Time_Open_Arm_OnlyHRLRGenes<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
rm(tempPvalAdj)

tempPvalAdj<-mt.rawp2adjp(F0_Meta_F2_4Models_LineageGenes$P.value.EPM_DistanceTraveled, proc="BH")
F0_Meta_F2_4Models_LineageGenes$FDR.EPM_DistanceTraveled_OnlyHRLRGenes<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
rm(tempPvalAdj)

tempPvalAdj<-mt.rawp2adjp(F0_Meta_F2_4Models_LineageGenes$P.value.LearningClassification_AsFactorGT, proc="BH")
F0_Meta_F2_4Models_LineageGenes$FDR.LearningClassification_AsFactorGT_OnlyHRLRGenes<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
rm(tempPvalAdj)

tempPvalAdj<-mt.rawp2adjp(F0_Meta_F2_4Models_LineageGenes$P.value.LearningClassification_AsFactorIN, proc="BH")
F0_Meta_F2_4Models_LineageGenes$FDR.LearningClassification_AsFactorIN_OnlyHRLRGenes<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
rm(tempPvalAdj)

tempPvalAdj<-mt.rawp2adjp(F0_Meta_F2_4Models_LineageGenes$P.value.EPM_DistanceTraveled_noOutlier, proc="BH")
F0_Meta_F2_4Models_LineageGenes$FDR.EPM_DistanceTraveled_noOutlier_OnlyHRLRGenes<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
rm(tempPvalAdj)

tempPvalAdj<-mt.rawp2adjp(F0_Meta_F2_4Models_LineageGenes$P.value.EPM_Time_Immobile, proc="BH")
F0_Meta_F2_4Models_LineageGenes$FDR.P.value.EPM_Time_Immobile_OnlyHRLRGenes<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
rm(tempPvalAdj)

tempPvalAdj<-mt.rawp2adjp(F0_Meta_F2_4Models_LineageGenes$P.value.PCA_Index_Days6and.7, proc="BH")
F0_Meta_F2_4Models_LineageGenes$FDR.PCA_Index_Days6and.7_OnlyHRLRGenes<-tempPvalAdj$adjp[order(tempPvalAdj$index),2]
rm(tempPvalAdj)

#write.csv(F0_Meta_F2_4Models_LineageGenes, "F0_Meta_F2_4Models_LineageGenes.csv")

write.csv(F0_Meta_F2_4Models_LineageGenes, "F0_Meta_F2_4Models_LineageGenes_Updated.csv")

#Only locoscore and distance traveled have significant DE genes - let's pull out their info so that we can easily pull up their plots:
F0_Meta_F2_4Models_LineageGenes$ENSEMBL[which(F0_Meta_F2_4Models_LineageGenes$M9_FDR.Total_LocoScore_OnlyHRLRGenes<0.10|F0_Meta_F2_4Models_LineageGenes$FDR.EPM_DistanceTraveled_noOutlier_OnlyHRLRGenes<0.10)]
# [1] "ENSRNOG00000052237" "ENSRNOG00000017854" "ENSRNOG00000062164" "ENSRNOG00000038330" "ENSRNOG00000008416" "ENSRNOG00000052519" "ENSRNOG00000005975"
# [8] "ENSRNOG00000004660" "ENSRNOG00000015020" "ENSRNOG00000002256" "ENSRNOG00000002810" "ENSRNOG00000017577" "ENSRNOG00000017820" "ENSRNOG00000039969"
# [15] "ENSRNOG00000015150" "ENSRNOG00000028904" "ENSRNOG00000026994"

#Note - changing the Distance Traveled to the version without an outlier didn't change the final gene list, even though one more of those genes was now significant for Distance Traveled as well as locomotor activity.

F0_Meta_F2_DEGenes_FDR10_ENSEMBL<-F0_Meta_F2_4Models_LineageGenes$ENSEMBL[which(F0_Meta_F2_4Models_LineageGenes$M9_FDR.Total_LocoScore_OnlyHRLRGenes<0.10|F0_Meta_F2_4Models_LineageGenes$FDR.EPM_DistanceTraveled_noOutlier_OnlyHRLRGenes<0.10)]

F0_Meta_F2_4Models_LineageGenes$SYMBOL[which(F0_Meta_F2_4Models_LineageGenes$M9_FDR.Total_LocoScore_OnlyHRLRGenes<0.10|F0_Meta_F2_4Models_LineageGenes$FDR.EPM_DistanceTraveled_noOutlier_OnlyHRLRGenes<0.10)]
# [1] NA           "Ucp2"       NA           "RGD1359508" "Gimap5"     "Ttc30a1"    "Rpl30"      "Fzd6"       "Idh1"       "Art3"       "Gfpt2"     
# [12] "Bphl"       "Nqo2"       "Dsc2"       "Spg7"       "Vps9d1"     "Afg3l1"   

F0_Meta_F2_4Models_LineageGenes$Symbol[which(F0_Meta_F2_4Models_LineageGenes$M9_FDR.Total_LocoScore_OnlyHRLRGenes<0.10|F0_Meta_F2_4Models_LineageGenes$FDR.EPM_DistanceTraveled_noOutlier_OnlyHRLRGenes<0.10)]
#Same annotation

F0_Meta_F2_DEGenes_FDR10_SYMBOLorENSEMBL<-F0_Meta_F2_4Models_LineageGenes$SYMBOL[which(F0_Meta_F2_4Models_LineageGenes$M9_FDR.Total_LocoScore_OnlyHRLRGenes<0.10|F0_Meta_F2_4Models_LineageGenes$FDR.EPM_DistanceTraveled_noOutlier_OnlyHRLRGenes<0.10)]


F0_Meta_F2_DEGenes_FDR10_SYMBOLorENSEMBL[is.na(F0_Meta_F2_DEGenes_FDR10_SYMBOLorENSEMBL)]<-F0_Meta_F2_DEGenes_FDR10_ENSEMBL[is.na(F0_Meta_F2_DEGenes_FDR10_SYMBOLorENSEMBL)]

F0_Meta_F2_DEGenes_FDR10_SYMBOLorENSEMBL
# [1] "ENSRNOG00000052237" "Ucp2"               "ENSRNOG00000062164" "RGD1359508"         "Gimap5"             "Ttc30a1"            "Rpl30"             
# [8] "Fzd6"               "Idh1"               "Art3"               "Gfpt2"              "Bphl"               "Nqo2"               "Dsc2"              
# [15] "Spg7"               "Vps9d1"             "Afg3l1" 

#We have more up-to-date annotation than that for those two ENSEMBL genes missing symbols:

F0_Meta_F2_DEGenes_FDR10_SYMBOLorENSEMBL2<-F0_Meta_F2_DEGenes_FDR10_SYMBOLorENSEMBL

F0_Meta_F2_DEGenes_FDR10_SYMBOLorENSEMBL2[F0_Meta_F2_DEGenes_FDR10_SYMBOLorENSEMBL=="ENSRNOG00000052237"]<-"AABR07071904.1" 

F0_Meta_F2_DEGenes_FDR10_SYMBOLorENSEMBL2[F0_Meta_F2_DEGenes_FDR10_SYMBOLorENSEMBL=="ENSRNOG00000062164"]<-"Rn60_1_2216.1"

F0_Meta_F2_DEGenes_FDR10_SYMBOLorENSEMBL2
# [1] "AABR07071904.1" "Ucp2"           "Rn60_1_2216.1"  "RGD1359508"     "Gimap5"         "Ttc30a1"        "Rpl30"          "Fzd6"          
# [9] "Idh1"           "Art3"           "Gfpt2"          "Bphl"           "Nqo2"           "Dsc2"           "Spg7"           "Vps9d1"        
# [17] "Afg3l1" 

summary(F0_Meta_F2_4Models_LineageGenes$M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's 
# -0.001120 -0.000033 -0.000002 -0.000004  0.000027  0.000728        18 

F0_Meta_F2_4Models_LineageGenes$M9_FDR.Total_LocoScore_OnlyHRLRGenes

tiff("VolcanoPlot_F2_LocoScore_bHRbLRGenes.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F0_Meta_F2_4Models_LineageGenes, plot(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore, -log10(M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore), pch=19, main="Effect of LocoScore", cex.lab=1.8, cex.main=2, cex=0.6, ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F0_Meta_F2_4Models_LineageGenes, M9_FDR.Total_LocoScore_OnlyHRLRGenes<0.10), points(M9_CovSexIntergenicRrnaTechnical_Coef.Total_LocoScore, -log10(M9_CovSexIntergenicRrnaTechnical_P.value.Total_LocoScore), pch=19, col="green", cex=0.6))
dev.off()


# tiff("VolcanoPlot_F2_EPM_Distance_bHRbLRGenes.tiff", width = 5, height = 5,
#      units = 'in', res = 300, compression = "lzw")
# par(mai=c(1.02, 1,0.9,0.40))
# with(F0_Meta_F2_4Models_LineageGenes, plot(Coef.EPM_DistanceTraveled, -log10(P.value.EPM_DistanceTraveled), pch=19, main="Effect of EPM Distance", cex.lab=1.8, cex.main=2, cex=0.6, ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
# with(subset(F0_Meta_F2_4Models_LineageGenes, FDR.EPM_DistanceTraveled_OnlyHRLRGenes<0.10), points(Coef.EPM_DistanceTraveled, -log10(P.value.EPM_DistanceTraveled), pch=19, col="green", cex=0.6))
# dev.off()

tiff("VolcanoPlot_F2_EPM_Distance_NoOutliers_bHRbLRGenes.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F0_Meta_F2_4Models_LineageGenes, plot(Coef.EPM_DistanceTraveled_noOutlier, -log10(P.value.EPM_DistanceTraveled_noOutlier), pch=19, main="Effect of EPM Distance", cex.lab=1.8, cex.main=2, cex=0.6, ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F0_Meta_F2_4Models_LineageGenes, FDR.EPM_DistanceTraveled_noOutlier_OnlyHRLRGenes<0.10), points(Coef.EPM_DistanceTraveled_noOutlier, -log10(P.value.EPM_DistanceTraveled_noOutlier), pch=19, col="green", cex=0.6))
dev.off()


F0_Meta_F2_4Models_LineageGenes$FDR.EPM_Percent_Time_Open_Arm_OnlyHRLRGenes

tiff("VolcanoPlot_F2_EPM_OpenArms_bHRbLRGenes.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F0_Meta_F2_4Models_LineageGenes, plot(Coef.EPM_Percent_Time_Open_Arm, -log10(P.value.EPM_Percent_Time_Open_Arm), pch=19, main="Effect of EPM %Time Open", cex.lab=1.8, cex.main=2, cex=0.6, ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F0_Meta_F2_4Models_LineageGenes, FDR.EPM_Percent_Time_Open_Arm_OnlyHRLRGenes<0.10), points(Coef.EPM_Percent_Time_Open_Arm, -log10(P.value.EPM_Percent_Time_Open_Arm), pch=19, col="green", cex=0.6))
dev.off()

F0_Meta_F2_4Models_LineageGenes$FDR.LearningClassification_AsFactorGT

tiff("VolcanoPlot_F2_GTvsST_bHRbLRGenes.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F0_Meta_F2_4Models_LineageGenes, plot(Coef.LearningClassification_AsFactorGT, -log10(P.value.LearningClassification_AsFactorGT), pch=19, main="Effect of GT vs. ST", cex.lab=1.8, cex.main=2, cex=0.6, ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F0_Meta_F2_4Models_LineageGenes, FDR.LearningClassification_AsFactorGT<0.10), points(Coef.LearningClassification_AsFactorGT, -log10(P.value.LearningClassification_AsFactorGT), pch=19, col="green", cex=0.6))
dev.off()

tiff("VolcanoPlot_F2_INvsST_bHRbLRGenes.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F0_Meta_F2_4Models_LineageGenes, plot(Coef.LearningClassification_AsFactorIN, -log10(P.value.LearningClassification_AsFactorIN), pch=19, main="Effect of IN vs. ST", cex.lab=1.8, cex.main=2, cex=0.6, ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F0_Meta_F2_4Models_LineageGenes, FDR.LearningClassification_AsFactorIN<0.10), points(Coef.LearningClassification_AsFactorIN, -log10(P.value.LearningClassification_AsFactorIN), pch=19, col="green", cex=0.6))
dev.off()

tiff("VolcanoPlot_F2_PavCADay67_bHRbLRGenes.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F0_Meta_F2_4Models_LineageGenes, plot(Coef.PCA_Index_Days6and.7, -log10(P.value.PCA_Index_Days6and.7), pch=19, main="Effect of PavCA Index (Days 6 & 7)", cex.lab=1.8, cex.main=2, cex=0.6, ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F0_Meta_F2_4Models_LineageGenes, FDR.PCA_Index_Days6and.7_OnlyHRLRGenes<0.10), points(Coef.PCA_Index_Days6and.7, -log10(P.value.PCA_Index_Days6and.7), pch=19, col="green", cex=0.6))
dev.off()

tiff("VolcanoPlot_F2_EPM_Time_Immobile_bHRbLRGenes.tiff", width = 5, height = 5,
     units = 'in', res = 300, compression = "lzw")
par(mai=c(1.02, 1,0.9,0.40))
with(F0_Meta_F2_4Models_LineageGenes, plot(Coef.EPM_Time_Immobile, -log10(P.value.EPM_Time_Immobile), pch=19, main="Effect of EPM %Time Immobile", cex.lab=1.8, cex.main=2, cex=0.6, ylim=c(0,7), xlab="Coefficient (log2 fold change)", ylab="-log10(p-value)"))
with(subset(F0_Meta_F2_4Models_LineageGenes, FDR.P.value.EPM_Time_Immobile_OnlyHRLRGenes<0.10), points(Coef.EPM_Time_Immobile, -log10(P.value.EPM_Time_Immobile), pch=19, col="green", cex=0.6))
dev.off()
