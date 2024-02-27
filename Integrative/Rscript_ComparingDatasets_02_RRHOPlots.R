#Comparing RNA-Seq Results across datasets: F0 vs. F2 vs. bHR/bLR Meta-analysis
#02_RRHO plots
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-27, updated later for a few figures for the paper.


library(RRHO)
#For making individual RRHO plots/stats:

F0_Meta_F2_4Models$PvalNegLog10_Lineage_AsFactorbLR<-(-log10(F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR))

#Double-Checking:
plot(F0_Meta_F2_4Models$PvalNegLog10_Lineage_AsFactorbLR~F0_Meta_F2_4Models$P.value.Lineage_AsFactorbLR)

F0_Meta_F2_4Models$Coef_Direction_Lineage_AsFactorbLR[F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR<0]<-(-1)
F0_Meta_F2_4Models$Coef_Direction_Lineage_AsFactorbLR[F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR>0]<-(1)

#Double-Checking:

plot(F0_Meta_F2_4Models$Coef_Direction_Lineage_AsFactorbLR~F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR)

F0_Meta_F2_4Models$PvalNegLog10_Meta_Lineage_bLR<-(-log10(F0_Meta_F2_4Models$pval))

#Double-Checking:
plot(F0_Meta_F2_4Models$PvalNegLog10_Lineage_AsFactorbLR~F0_Meta_F2_4Models$t.Lineage_AsFactorbLR)
#Perfectly monotonic (in each direction) - I could just use t-stats to create my ranks and save myself a lot of trouble...

TempDF<-F0_Meta_F2_4Models[is.na(F0_Meta_F2_4Models$t.Lineage_AsFactorbLR)==FALSE & is.na(F0_Meta_F2_4Models$estimate)==FALSE, ]

str(TempDF)
#'data.frame':	11175 obs. of  228 variables:

list1<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=TempDF$t.Lineage_AsFactorbLR)

list2<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=(TempDF$estimate/TempDF$SE))


setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_F0vsMetaAnalysis_TwoSided")

RRHO(list1, list2, labels=c("F0: bLR vs. bHR", "Meta-Analysis: bLR vs. bHR"), plots=TRUE, alternative="two.sided", outputdir="~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_F0vsMetaAnalysis_TwoSided", BY=TRUE, log10.ind=TRUE)

# $n.items
# [1] 11175
# 
# $stepsize
# [1] 106
# 
# $log10.ind
# [1] TRUE

#I want to make sure that I am interpreting the color properly, so I'm going to make an artificial negative correlation with the same magnitude:


list1<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=TempDF$t.Lineage_AsFactorbLR)

list2<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=(TempDF$estimate/TempDF$SE)*-1)

setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_test")

RRHO(list1, list2, labels=c("F0: bLR vs. bHR", "Inverted Meta-Analysis: bLR vs. bHR"), plots=TRUE, alternative="two.sided", outputdir="~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_test", BY=TRUE, log10.ind=TRUE)
#interesting - yeah, the colors and p-values invert for a negative correlation. I'm not sure what to make of that.

# RRHO(list1, list2, labels=c("F0: bLR vs. bHR", "Inverted Meta-Analysis: bLR vs. bHR"), plots=TRUE, alternative="enrichment", outputdir="~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_test")
#This doesn't work for negative correlations at all

TempDF<-F0_Meta_F2_4Models[is.na(F0_Meta_F2_4Models$t.Lineage_AsFactorbLR)==FALSE & is.na(F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore)==FALSE, ]

str(TempDF)
#'data.frame':	13339 obs. of  228 variables:

list1<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=TempDF$t.Lineage_AsFactorbLR)

list2<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=TempDF$M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore)

setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_F0vsF2LocoScore_TwoSided")

RRHO(list1, list2, labels=c("F0: bLR vs. bHR", "F2: Total LocoScore"), plots=TRUE, alternative="two.sided", outputdir="~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_F0vsF2LocoScore_TwoSided",  BY=TRUE, log10.ind=TRUE)


TempDF<-F0_Meta_F2_4Models[is.na(F0_Meta_F2_4Models$estimate)==FALSE & is.na(F0_Meta_F2_4Models$M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore)==FALSE, ]

str(TempDF)
#'data.frame':	11055 obs. of  228 variables:

list1<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=(TempDF$estimate/TempDF$SE))

list2<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=TempDF$M9_CovSexIntergenicRrnaTechnical_t.Total_LocoScore)

setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_MetaAnalysisVsF2LocoScore_TwoSided")

RRHO(list1, list2, labels=c("Meta-Analysis: bLR vs. bHR", "F2: Total LocoScore"), plots=TRUE, alternative="two.sided", outputdir="~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_MetaAnalysisVsF2LocoScore_TwoSided", BY=TRUE, log10.ind=TRUE)
#ok


TempDF<-F0_Meta_F2_4Models[is.na(F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR)==FALSE & is.na(F0_Meta_F2_4Models$t.EPM_DistanceTraveled_noOutlier)==FALSE, ]

str(TempDF)
#'data.frame':	13339 obs. of  228 variables:

list1<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=TempDF$t.Lineage_AsFactorbLR)

list2<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=TempDF$t.EPM_DistanceTraveled_noOutlier)

setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_F0vsF2Distance_TwoSided")

RRHO(list1, list2, labels=c("F0: bLR vs. bHR", "F2: EPM Distance Traveled"), plots=TRUE, alternative="two.sided", outputdir="~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_F0vsF2Distance_TwoSided", BY=TRUE, log10.ind=TRUE)



TempDF<-F0_Meta_F2_4Models[is.na(F0_Meta_F2_4Models$estimate)==FALSE & is.na(F0_Meta_F2_4Models$t.EPM_DistanceTraveled_noOutlier)==FALSE, ]

str(TempDF)
#'data.frame':	11055 obs. of  228 variables:

list1<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=(TempDF$estimate/TempDF$SE))

list2<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=TempDF$t.EPM_DistanceTraveled_noOutlier)

setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_MetaAnalysisVsF2EPMDistance_TwoSided")

RRHO(list1, list2, labels=c("Meta-Analysis: bLR vs. bHR", "F2: EPM Distance Traveled"), plots=TRUE, alternative="two.sided", outputdir="~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_MetaAnalysisVsF2EPMDistance_TwoSided", BY=TRUE, log10.ind=TRUE)


TempDF<-F0_Meta_F2_4Models[is.na(F0_Meta_F2_4Models$Coef.Lineage_AsFactorbLR)==FALSE & is.na(F0_Meta_F2_4Models$t.EPM_Percent_Time_Open_Arm)==FALSE, ]

str(TempDF)
#'data.frame':	13339 obs. of  228 variables:

list1<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=TempDF$t.Lineage_AsFactorbLR)

list2<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=TempDF$t.EPM_Percent_Time_Open_Arm)

setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_F0vsF2EPMTimeOpenArms_TwoSided")

RRHO(list1, list2, labels=c("F0: bLR vs. bHR", "F2: EPM Percent Time Open Arms"), plots=TRUE, alternative="two.sided", outputdir="~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_F0vsF2EPMTimeOpenArms_TwoSided", BY=TRUE, log10.ind=TRUE)


TempDF<-F0_Meta_F2_4Models[is.na(F0_Meta_F2_4Models$estimate)==FALSE & is.na(F0_Meta_F2_4Models$t.EPM_Percent_Time_Open_Arm)==FALSE, ]

str(TempDF)
#'data.frame':	11055 obs. of  228 variables:

list1<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=(TempDF$estimate/TempDF$SE))

list2<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=TempDF$t.EPM_Percent_Time_Open_Arm)

setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_MetaAnalysisVsEPMTimeOpenArms")

RRHO(list1, list2, labels=c("Meta-Analysis: bLR vs. bHR", "F2: EPM Percent Time Open Arms"), plots=TRUE, alternative="two.sided", outputdir="~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_MetaAnalysisVsEPMTimeOpenArms", BY=TRUE, log10.ind=TRUE)



TempDF<-F0_Meta_F2_4Models[is.na(F0_Meta_F2_4Models$t.Lineage_AsFactorbLR)==FALSE & is.na(F0_Meta_F2_4Models$t.PCA_Index_Days6and.7)==FALSE, ]

str(TempDF)
#'data.frame':	11055 obs. of  228 variables:

list1<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=TempDF$t.Lineage_AsFactorbLR)

list2<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=TempDF$t.PCA_Index_Days6and.7)


setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_F0vsF2PavCA")

RRHO(list1, list2, labels=c("F0: bLR vs. bHR", "F2: PavCA Index"), plots=TRUE, alternative="two.sided", outputdir="~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_F0vsF2PavCA", BY=TRUE, log10.ind=TRUE)


TempDF<-F0_Meta_F2_4Models[is.na(F0_Meta_F2_4Models$estimate)==FALSE & is.na(F0_Meta_F2_4Models$t.PCA_Index_Days6and.7)==FALSE, ]

str(TempDF)
#'data.frame':	11055 obs. of  228 variables:

list1<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=(TempDF$estimate/TempDF$SE))

list2<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=TempDF$t.PCA_Index_Days6and.7)

setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_MetaAnalysisVsF2PavCA")

RRHO(list1, list2, labels=c("MetaAnalysis: bLR vs. bHR", "F2: PavCA Index"), plots=TRUE, alternative="two.sided", outputdir="~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_MetaAnalysisVsF2PavCA", BY=TRUE, log10.ind=TRUE)



TempDF<-F0_Meta_F2_4Models[is.na(F0_Meta_F2_4Models$t.Lineage_AsFactorbLR)==FALSE & is.na(F0_Meta_F2_4Models$t.EPM_Time_Immobile)==FALSE, ]

str(TempDF)
#'data.frame':	13339 obs. of  228 variables:

list1<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=TempDF$t.Lineage_AsFactorbLR)

list2<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=TempDF$t.EPM_Time_Immobile)


setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_F0vsF2EPMTimeImmobile")

RRHO(list1, list2, labels=c("F0: bLR vs. bHR", "F2: EPM Time Immobile"), plots=TRUE, alternative="two.sided", outputdir="~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_F0vsF2EPMTimeImmobile", BY=TRUE, log10.ind=TRUE)


TempDF<-F0_Meta_F2_4Models[is.na(F0_Meta_F2_4Models$estimate)==FALSE & is.na(F0_Meta_F2_4Models$t.EPM_Time_Immobile)==FALSE, ]

str(TempDF)
#'data.frame':	11055 obs. of  228 variables:

list1<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=(TempDF$estimate/TempDF$SE))

list2<-data.frame(ENSEMBL=TempDF$ENSEMBL, Metric=TempDF$t.EPM_Time_Immobile)

setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_MetaAnalysisVsF2TimeImmobile")

RRHO(list1, list2, labels=c("Meta-Analysis: bLR vs. bHR", "F2: EPM Time Immobile"), plots=TRUE, alternative="two.sided", outputdir="~/Documents/Microarray Gen/HRLR/NIDA_U01/RRHO_MetaAnalysisVsF2TimeImmobile", BY=TRUE, log10.ind=TRUE)

