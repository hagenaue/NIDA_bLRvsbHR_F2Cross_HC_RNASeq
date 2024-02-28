#NIDA U01 RNA-Seq paper:
#Adding some stats for the comparison with the Birt et al. database of Hippocampal DE Genes from other genetic internalizing-like vs. externalizing-like models
#Megan Hagenauer, 2023-08-09

#I read the F0/F2 results back in using the final supplementary table:

setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/DEResults_F0F2_SummaryTables")

F0_Meta_F2_results_wLineageGenes_forSuppl<-read.csv("F0_Meta_F2_results_wLineageGenes_ForSuppl.csv", header=TRUE, stringsAsFactors = FALSE)
str(F0_Meta_F2_results_wLineageGenes_forSuppl)
#'data.frame':	13786 obs. of  48 variables:

F0_Meta_F2_LineageGenes<-F0_Meta_F2_results_wLineageGenes_forSuppl[is.na(F0_Meta_F2_results_wLineageGenes_forSuppl$M9_FDR.Total_LocoScore_OnlyHRLRGenes)==FALSE,]


OtherRatModels<-read.csv("InternalizingGenesDatabase_fromBirt2021.csv", header=TRUE, stringsAsFactors = FALSE)

str(OtherRatModels)

# 'data.frame':	3219 obs. of  2 variables:
# $ GeneSymbol_Rat: chr  "Abhd10" "Acbd4" "Acss2" "Adam19" ...
# $ DatasetName   : chr  "Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC" "Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC" "Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC" "Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC"

length(unique(OtherRatModels$GeneSymbol_Rat))
#[1] 2580

#Oh wait - the essential comparison here is actually smaller, because we don't compare with the bHR/bLR genes

length(unique(OtherRatModels$GeneSymbol_Rat[OtherRatModels$DatasetName!="Birt_Hagenauer_2020_BredHighRespondersToNovelty_Upregulated_HC"& OtherRatModels$DatasetName!="Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC"& OtherRatModels$DatasetName!="Replicated_InternalizingBehavior_DifferentiallyExpressed_HC"]))
#[1] 2447
#I'll need to fix that in the manuscript***

ActuallyOtherRatModels<-OtherRatModels[OtherRatModels$DatasetName!="Birt_Hagenauer_2020_BredHighRespondersToNovelty_Upregulated_HC"& OtherRatModels$DatasetName!="Birt_Hagenauer_2020_BredLowRespondersToNovelty_Upregulated_HC" & OtherRatModels$DatasetName!="Replicated_InternalizingBehavior_DifferentiallyExpressed_HC",]
dim(ActuallyOtherRatModels)
#[1] 2732    2

sum(unique(ActuallyOtherRatModels$GeneSymbol_Rat)%in%F0_Meta_F2_results_wLineageGenes_forSuppl$SYMBOL.1)
#[1] 1964

table(ActuallyOtherRatModels$DatasetName)
# Andrus_2012_WistarMoreImmobileRats_Downregulated_HC     Andrus_2012_WistarMoreImmobileRats_Upregulated_HC 
# 217                                                   242 
# Blaveri_2010_FlindersSensitiveLine_Downregulated_HC     Blaveri_2010_FlindersSensitiveLine_Upregulated_HC 
# 613                                                   295 
# DiazMoran_2013_NIHHighAnxietyRats_Downregulated_HC      DiazMoran_2013_NIHHighAnxietyRats_Upregulated_HC 
# 89                                                   120 
# Garafola_Hen_2014_CongenitallyLearnedHelpless_HC                Meckes_2018_WKYvsF344_Downregulated_HC 
# 13                                                   454 
# Meckes_2018_WKYvsF344_Upregulated_HC Raghavan_2017_WistarMoreImmobileRats_Downregulated_HC 
# 389                                                    21 
# Raghavan_2017_WistarMoreImmobileRats_Upregulated_HC   Sabariego_2013_RomanLowAvoidanceRats_Upregulated_HC 
# 25                                                    10 
# Wilhelm_2013_FlindersSensitiveLine_Downregulated_HC     Wilhelm_2013_FlindersSensitiveLine_Upregulated_HC 
# 123                                                    85 
# Zhang_2005_SyracuseLowAvoidanceRats_Downregulated_HC    Zhang_2005_SyracuseLowAvoidanceRats_Upregulated_HC 
# 15                                                    21 

#Note: These may not all be unique gene SYMBOL.1s - double check?

for(i in c(1:length(names(table(ActuallyOtherRatModels$DatasetName))))){
  temp<-ActuallyOtherRatModels[ActuallyOtherRatModels$DatasetName==names(table(ActuallyOtherRatModels$DatasetName))[i],]
  print(names(table(ActuallyOtherRatModels$DatasetName))[i])
  print(length(unique(temp$GeneSymbol_Rat)))
}

# [1] "Andrus_2012_WistarMoreImmobileRats_Downregulated_HC"
# [1] 217
# [1] "Andrus_2012_WistarMoreImmobileRats_Upregulated_HC"
# [1] 242
# [1] "Blaveri_2010_FlindersSensitiveLine_Downregulated_HC"
# [1] 613
# [1] "Blaveri_2010_FlindersSensitiveLine_Upregulated_HC"
# [1] 295
# [1] "DiazMoran_2013_NIHHighAnxietyRats_Downregulated_HC"
# [1] 89
# [1] "DiazMoran_2013_NIHHighAnxietyRats_Upregulated_HC"
# [1] 120
# [1] "Garafola_Hen_2014_CongenitallyLearnedHelpless_HC"
# [1] 13
# [1] "Meckes_2018_WKYvsF344_Downregulated_HC"
# [1] 454
# [1] "Meckes_2018_WKYvsF344_Upregulated_HC"
# [1] 389
# [1] "Raghavan_2017_WistarMoreImmobileRats_Downregulated_HC"
# [1] 21
# [1] "Raghavan_2017_WistarMoreImmobileRats_Upregulated_HC"
# [1] 25
# [1] "Sabariego_2013_RomanLowAvoidanceRats_Upregulated_HC"
# [1] 10
# [1] "Wilhelm_2013_FlindersSensitiveLine_Downregulated_HC"
# [1] 123
# [1] "Wilhelm_2013_FlindersSensitiveLine_Upregulated_HC"
# [1] 85
# [1] "Zhang_2005_SyracuseLowAvoidanceRats_Downregulated_HC"
# [1] 15
# [1] "Zhang_2005_SyracuseLowAvoidanceRats_Upregulated_HC"
# [1] 21

#Note: Many of these gene symbols may not be in our dataset - double check?

for(i in c(1:length(names(table(ActuallyOtherRatModels$DatasetName))))){
  temp<-ActuallyOtherRatModels[ActuallyOtherRatModels$DatasetName==names(table(ActuallyOtherRatModels$DatasetName))[i],]
  print(names(table(ActuallyOtherRatModels$DatasetName))[i])
  print(sum(unique(temp$GeneSymbol_Rat)%in%F0_Meta_F2_results_wLineageGenes_forSuppl$SYMBOL.1))
}

# [1] "Andrus_2012_WistarMoreImmobileRats_Downregulated_HC"
# [1] 168
# [1] "Andrus_2012_WistarMoreImmobileRats_Upregulated_HC"
# [1] 161
# [1] "Blaveri_2010_FlindersSensitiveLine_Downregulated_HC"
# [1] 571
# [1] "Blaveri_2010_FlindersSensitiveLine_Upregulated_HC"
# [1] 259
# [1] "DiazMoran_2013_NIHHighAnxietyRats_Downregulated_HC"
# [1] 34
# [1] "DiazMoran_2013_NIHHighAnxietyRats_Upregulated_HC"
# [1] 44
# [1] "Garafola_Hen_2014_CongenitallyLearnedHelpless_HC"
# [1] 7
# [1] "Meckes_2018_WKYvsF344_Downregulated_HC"
# [1] 413
# [1] "Meckes_2018_WKYvsF344_Upregulated_HC"
# [1] 342
# [1] "Raghavan_2017_WistarMoreImmobileRats_Downregulated_HC"
# [1] 13
# [1] "Raghavan_2017_WistarMoreImmobileRats_Upregulated_HC"
# [1] 17
# [1] "Sabariego_2013_RomanLowAvoidanceRats_Upregulated_HC"
# [1] 6
# [1] "Wilhelm_2013_FlindersSensitiveLine_Downregulated_HC"
# [1] 95
# [1] "Wilhelm_2013_FlindersSensitiveLine_Upregulated_HC"
# [1] 70
# [1] "Zhang_2005_SyracuseLowAvoidanceRats_Downregulated_HC"
# [1] 13
# [1] "Zhang_2005_SyracuseLowAvoidanceRats_Upregulated_HC"
# [1] 12


#Alright - so the database already reflects unique gene SYMBOL.1s for each dataset.
dim(F0_Meta_F2_LineageGenes)#this was 1063, then dropped to 1045. Why?
#[1] 1045   48

dim(F0_Meta_F2_LineageGenes[is.na(F0_Meta_F2_LineageGenes$P.value.Total_LocoScore)==FALSE,])
#[1] 1045   48 
#Ah, that's right - this version is limited to the genes that are also represented in the F2 dataset.

#Double-checking: how many rows with both gene SYMBOL.1s and F2 results?
length(F0_Meta_F2_LineageGenes$SYMBOL.1[is.na(F0_Meta_F2_LineageGenes$SYMBOL.1)==FALSE & is.na(F0_Meta_F2_LineageGenes$P.value.Total_LocoScore)==FALSE])
#[1] 1045
#Are these gene SYMBOL.1s unique?
length(unique(F0_Meta_F2_LineageGenes$SYMBOL.1[is.na(F0_Meta_F2_LineageGenes$P.value.Total_LocoScore)==FALSE]))
#[1] 1045

sum(F0_Meta_F2_LineageGenes$SYMBOL.1%in%ActuallyOtherRatModels$GeneSymbol_Rat)
#[1] 247

sum(F0_Meta_F2_LineageGenes$SYMBOL.1[F0_Meta_F2_LineageGenes$P.value.Total_LocoScore<0.05]%in%ActuallyOtherRatModels$GeneSymbol_Rat)
#[1] 40

#Which of these represent different directions?
names(table(ActuallyOtherRatModels$DatasetName))
# [1] "Andrus_2012_WistarMoreImmobileRats_Downregulated_HC"   "Andrus_2012_WistarMoreImmobileRats_Upregulated_HC"    
# [3] "Blaveri_2010_FlindersSensitiveLine_Downregulated_HC"   "Blaveri_2010_FlindersSensitiveLine_Upregulated_HC"    
# [5] "DiazMoran_2013_NIHHighAnxietyRats_Downregulated_HC"    "DiazMoran_2013_NIHHighAnxietyRats_Upregulated_HC"     
# [7] "Garafola_Hen_2014_CongenitallyLearnedHelpless_HC"      "Meckes_2018_WKYvsF344_Downregulated_HC"               
# [9] "Meckes_2018_WKYvsF344_Upregulated_HC"                  "Raghavan_2017_WistarMoreImmobileRats_Downregulated_HC"
# [11] "Raghavan_2017_WistarMoreImmobileRats_Upregulated_HC"   "Sabariego_2013_RomanLowAvoidanceRats_Upregulated_HC"  
# [13] "Wilhelm_2013_FlindersSensitiveLine_Downregulated_HC"   "Wilhelm_2013_FlindersSensitiveLine_Upregulated_HC"    
# [15] "Zhang_2005_SyracuseLowAvoidanceRats_Downregulated_HC"  "Zhang_2005_SyracuseLowAvoidanceRats_Upregulated_HC" 

DatasetsDownregulated<-c(1,3,5,8,10,13,15)
DatasetsUpregulated<-c(2,4,6,9,11,12,14,16)

OtherRatModels_Down<-ActuallyOtherRatModels[ActuallyOtherRatModels$DatasetName==names(table(ActuallyOtherRatModels$DatasetName))[1],]

for(i in DatasetsDownregulated[-1]){
  temp<-ActuallyOtherRatModels[ActuallyOtherRatModels$DatasetName==names(table(ActuallyOtherRatModels$DatasetName))[i],]
  print(names(table(ActuallyOtherRatModels$DatasetName))[i])
  OtherRatModels_Down<<-rbind.data.frame(OtherRatModels_Down, temp)
}
# [1] "Blaveri_2010_FlindersSensitiveLine_Downregulated_HC"
# [1] "DiazMoran_2013_NIHHighAnxietyRats_Downregulated_HC"
# [1] "Meckes_2018_WKYvsF344_Downregulated_HC"
# [1] "Raghavan_2017_WistarMoreImmobileRats_Downregulated_HC"
# [1] "Wilhelm_2013_FlindersSensitiveLine_Downregulated_HC"
# [1] "Zhang_2005_SyracuseLowAvoidanceRats_Downregulated_HC"

table(OtherRatModels_Down$DatasetName)
# Andrus_2012_WistarMoreImmobileRats_Downregulated_HC   Blaveri_2010_FlindersSensitiveLine_Downregulated_HC 
# 217                                                   613 
# DiazMoran_2013_NIHHighAnxietyRats_Downregulated_HC    Meckes_2018_WKYvsF344_Downregulated_HC 
# 89                                                    454 
# Raghavan_2017_WistarMoreImmobileRats_Downregulated_HC   Wilhelm_2013_FlindersSensitiveLine_Downregulated_HC 
# 21                                                      123 
# Zhang_2005_SyracuseLowAvoidanceRats_Downregulated_HC 
# 15 

length(unique(OtherRatModels_Down$GeneSymbol_Rat))
#[1] 1437

OtherRatModels_Up<-ActuallyOtherRatModels[ActuallyOtherRatModels$DatasetName==names(table(ActuallyOtherRatModels$DatasetName))[2],]

for(i in DatasetsUpregulated[-1]){
  temp<-ActuallyOtherRatModels[ActuallyOtherRatModels$DatasetName==names(table(ActuallyOtherRatModels$DatasetName))[i],]
  print(names(table(ActuallyOtherRatModels$DatasetName))[i])
  OtherRatModels_Up<<-rbind.data.frame(OtherRatModels_Up, temp)
}

# [1] "Blaveri_2010_FlindersSensitiveLine_Upregulated_HC"
# [1] "DiazMoran_2013_NIHHighAnxietyRats_Upregulated_HC"
# [1] "Meckes_2018_WKYvsF344_Upregulated_HC"
# [1] "Raghavan_2017_WistarMoreImmobileRats_Upregulated_HC"
# [1] "Sabariego_2013_RomanLowAvoidanceRats_Upregulated_HC"
# [1] "Wilhelm_2013_FlindersSensitiveLine_Upregulated_HC"
# [1] "Zhang_2005_SyracuseLowAvoidanceRats_Upregulated_HC"

table(OtherRatModels_Up$DatasetName)

# Andrus_2012_WistarMoreImmobileRats_Upregulated_HC   Blaveri_2010_FlindersSensitiveLine_Upregulated_HC 
# 242                                                 295 
# DiazMoran_2013_NIHHighAnxietyRats_Upregulated_HC    Meckes_2018_WKYvsF344_Upregulated_HC 
# 120                                                 389 
# Raghavan_2017_WistarMoreImmobileRats_Upregulated_HC Sabariego_2013_RomanLowAvoidanceRats_Upregulated_HC 
# 25                                                  10 
# Wilhelm_2013_FlindersSensitiveLine_Upregulated_HC  Zhang_2005_SyracuseLowAvoidanceRats_Upregulated_HC 
# 85                                                  21 

length(unique(OtherRatModels_Up$GeneSymbol_Rat))
#[1] 1121


bLRGenes<-((F0_Meta_F2_LineageGenes$F0_Coef.Lineage_AsFactorbLR>0 & F0_Meta_F2_LineageGenes$F0_P.value.adj.Lineage_AsFactorbLR<0.10)|
             (F0_Meta_F2_LineageGenes$MetaAnalysis_estimate>0 & F0_Meta_F2_LineageGenes$MetaAnalysis_BH<0.10)|
             (F0_Meta_F2_LineageGenes$F0_P.value.Lineage_AsFactorbLR<0.05 & F0_Meta_F2_LineageGenes$MetaAnalysis_rawp<0.05 & F0_Meta_F2_LineageGenes$F0_Coef.Lineage_AsFactorbLR>0 & F0_Meta_F2_LineageGenes$MetaAnalysis_estimate>0))
sum(bLRGenes, na.rm=TRUE)
#[1] 608

length(unique(F0_Meta_F2_LineageGenes$ENSEMBL[which(bLRGenes==TRUE)]))
#[1] 608
length(unique(F0_Meta_F2_LineageGenes$SYMBOL[which(bLRGenes==TRUE)]))
#[1] 594

bHRGenes<-((F0_Meta_F2_LineageGenes$F0_Coef.Lineage_AsFactorbLR<0 & F0_Meta_F2_LineageGenes$F0_P.value.adj.Lineage_AsFactorbLR<0.10)
                    |(F0_Meta_F2_LineageGenes$MetaAnalysis_estimate<0 & F0_Meta_F2_LineageGenes$MetaAnalysis_BH<0.10)
                    |(F0_Meta_F2_LineageGenes$F0_P.value.Lineage_AsFactorbLR<0.05 & F0_Meta_F2_LineageGenes$MetaAnalysis_rawp<0.05 & F0_Meta_F2_LineageGenes$F0_Coef.Lineage_AsFactorbLR<0 & F0_Meta_F2_LineageGenes$MetaAnalysis_estimate<0))
sum(bHRGenes, na.rm=TRUE)
#[1] 437
length(unique(F0_Meta_F2_LineageGenes$ENSEMBL[which(bHRGenes==TRUE)]))
#[1] 437
length(unique(F0_Meta_F2_LineageGenes$SYMBOL[which(bHRGenes==TRUE)]))
#[1] 428

#Sanity check:
608+437
#[1] 1045

setwd("/Users/hagenaue/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/F2_eQTLs/Input")

F0_Meta_F2_DEResults<-read.csv("TableS3_F0_Meta_F2_results_wLineageGenes_ForSuppl_20230419.csv", header=TRUE, stringsAsFactors = FALSE)
str(F0_Meta_F2_DEResults)
#'data.frame':	13788 obs. of  42 variables:

bLRLikeBehaviorGenes<-((F0_Meta_F2_DEResults$F2_Pvalue_LocoScore<0.05 & F0_Meta_F2_DEResults$F2_Log2FC_LocoScore<0)
                       |(F0_Meta_F2_DEResults$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults$F2_Log2FC_EPM_Time_Immobile>0)
                       |(F0_Meta_F2_DEResults$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults$F2_Log2FC_EPM_DistanceTraveled<0)
                       |(F0_Meta_F2_DEResults$F2_Pvalue_EPM_Percent_Time_Open_Arms<0.05 & F0_Meta_F2_DEResults$F2_Log2FC_EPM_Percent_Time_Open_Arms<0)
                       |(F0_Meta_F2_DEResults$F2_Pvalue_PCA_Index <0.05 & F0_Meta_F2_DEResults$F2_Log2FC_PCA_Index<0))
sum(bLRLikeBehaviorGenes, na.rm=TRUE)
#[1] 924
length(unique(F0_Meta_F2_DEResults$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bLRLikeBehaviorGenes==TRUE)]))
#[1] 924
length(unique(F0_Meta_F2_DEResults$GENENAME..Rnor6.Ensembl.v.88.[which(bLRLikeBehaviorGenes==TRUE)]))
#[1] 924

bHRLikeBehaviorGenes<-((F0_Meta_F2_DEResults$F2_Pvalue_LocoScore<0.05 & F0_Meta_F2_DEResults$F2_Log2FC_LocoScore>0)
                       |(F0_Meta_F2_DEResults$F2_Pvalue_EPM_Time_Immobile<0.05 & F0_Meta_F2_DEResults$F2_Log2FC_EPM_Time_Immobile<0)
                       |(F0_Meta_F2_DEResults$F2_Pvalue_EPM_DistanceTraveled<0.05 & F0_Meta_F2_DEResults$F2_Log2FC_EPM_DistanceTraveled>0)
                       |(F0_Meta_F2_DEResults$F2_Pvalue_EPM_Percent_Time_Open_Arms<0.05 & F0_Meta_F2_DEResults$F2_Log2FC_EPM_Percent_Time_Open_Arms>0)
                       |(F0_Meta_F2_DEResults$F2_Pvalue_PCA_Index<0.05 & F0_Meta_F2_DEResults$F2_Log2FC_PCA_Index>0))
sum(bHRLikeBehaviorGenes, na.rm=TRUE) 
#[1] 1075

length(unique(F0_Meta_F2_DEResults$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bHRLikeBehaviorGenes==TRUE)]))
#[1] 1075
length(unique(F0_Meta_F2_DEResults$GENENAME..Rnor6.Ensembl.v.88.[which(bHRLikeBehaviorGenes==TRUE)]))
#[1] 1075

dim(F0_Meta_F2_DEResults[is.na(F0_Meta_F2_DEResults$F2_Pvalue_LocoScore)==FALSE,])
#[1] 13339    42

#############

#Just another sanity check, since I'm getting some discrepancies between workspaces/analyses:

#Using the full results instead:
bLRGenesFullResults<-((F0_Meta_F2_DEResults$F0_Log2FC_bLRvsbHR>0 & F0_Meta_F2_DEResults$F0_FDR_bLRvsbHR<0.10)|
             (F0_Meta_F2_DEResults$MetaAnalysis_estimatedD_bLRvsbHR>0 & F0_Meta_F2_DEResults$MetaAnalysis_FDR_bLRvsbHR<0.10)|
             (F0_Meta_F2_DEResults$F0_Pvalue_bLRvsbHR<0.05 & F0_Meta_F2_DEResults$MetaAnalysis_Pvalue_bLRvsbHR<0.05 & F0_Meta_F2_DEResults$F0_Log2FC_bLRvsbHR>0 & F0_Meta_F2_DEResults$MetaAnalysis_estimatedD_bLRvsbHR>0))
sum(bLRGenesFullResults, na.rm=TRUE)
#[1] 619

length(unique(F0_Meta_F2_DEResults$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bLRGenesFullResults==TRUE)]))
#[1] 619
length(unique(F0_Meta_F2_DEResults$GENENAME..Rnor6.Ensembl.v.88.[which(bLRGenesFullResults==TRUE)]))
#[1] 619

bHRGenesFullResults<-((F0_Meta_F2_DEResults$F0_Log2FC_bLRvsbHR<0 & F0_Meta_F2_DEResults$F0_FDR_bLRvsbHR<0.10)|
                        (F0_Meta_F2_DEResults$MetaAnalysis_estimatedD_bLRvsbHR<0 & F0_Meta_F2_DEResults$MetaAnalysis_FDR_bLRvsbHR<0.10)|
                        (F0_Meta_F2_DEResults$F0_Pvalue_bLRvsbHR<0.05 & F0_Meta_F2_DEResults$MetaAnalysis_Pvalue_bLRvsbHR<0.05 & F0_Meta_F2_DEResults$F0_Log2FC_bLRvsbHR<0 & F0_Meta_F2_DEResults$MetaAnalysis_estimatedD_bLRvsbHR<0))
sum(bHRGenesFullResults, na.rm=TRUE)
#[1] 444

length(unique(F0_Meta_F2_DEResults$ENSEMBL.ID..Rnor6.Ensembl.v.103.[which(bHRGenesFullResults==TRUE)]))
#[1] 444
length(unique(F0_Meta_F2_DEResults$GENENAME..Rnor6.Ensembl.v.88.[which(bHRGenesFullResults==TRUE)]))
#[1] 444

619+444
#[1] 1063 
#Ah - so the difference is because this is the full # of HR/LR genes before considering what is represented in the F2 dataset.
#The lineage genes version of the database must have already been pruned for genes represented in the F2 dataset.

#############

#When I looked at the overlap by hand before:
#I looked at the unique genes in each set.
#I used all F2 behaviors in our DE results 
#I also screened for similar direction of effect with Lineage too.
#If the F0 and MetaAnalysis conflicted, I prioritized the effect that was significant (FDR<0.10).

#I double-checked it visually (since it was a small set) - harder to do with coding
#but if I want to do a Venn diagram style analysis, I will need to do it with coding.
#Honestly, the direction of the averaged t-statistic might have been an easier way to go about this.

colnames(F0_Meta_F2_LineageGenes)
#Criteria using all F2 variables:

bLRLikeCriteria<-(((F0_Meta_F2_LineageGenes$P.value.Total_LocoScore<0.05 & F0_Meta_F2_LineageGenes$Coef.Total_LocoScore<0)
                   |(F0_Meta_F2_LineageGenes$P.value.EPM_Time_Immobile<0.05 & F0_Meta_F2_LineageGenes$Coef.EPM_Time_Immobile>0)
                   |(F0_Meta_F2_LineageGenes$P.value.EPM_DistanceTraveled_noOutlier<0.05 & F0_Meta_F2_LineageGenes$Coef.EPM_DistanceTraveled_noOutlier<0)
                   |(F0_Meta_F2_LineageGenes$P.value.EPM_Percent_Time_Open_Arm<0.05 & F0_Meta_F2_LineageGenes$Coef.EPM_Percent_Time_Open_Arm<0)
                   |(F0_Meta_F2_LineageGenes$P.value.PCA_Index_Days6and.7<0.05 & F0_Meta_F2_LineageGenes$Coef.PCA_Index_Days6and.7<0))
                  &((F0_Meta_F2_LineageGenes$F0_Coef.Lineage_AsFactorbLR>0 & F0_Meta_F2_LineageGenes$F0_P.value.adj.Lineage_AsFactorbLR<0.10)
                    |(F0_Meta_F2_LineageGenes$MetaAnalysis_estimate>0 & F0_Meta_F2_LineageGenes$MetaAnalysis_BH<0.10)
                    |(F0_Meta_F2_LineageGenes$F0_P.value.Lineage_AsFactorbLR<0.05 & F0_Meta_F2_LineageGenes$MetaAnalysis_rawp<0.05 & F0_Meta_F2_LineageGenes$F0_Coef.Lineage_AsFactorbLR>0 & F0_Meta_F2_LineageGenes$MetaAnalysis_estimate>0)))
sum(bLRLikeCriteria, na.rm=TRUE)
#[1] 111

#Just double-checking some coding:
length(bLRLikeCriteria[bLRLikeCriteria==TRUE])
#[1] 113 #with the NAs
length(bLRLikeCriteria[which(bLRLikeCriteria==TRUE)])
#[1] 111 #without the NAs
length(unique(F0_Meta_F2_LineageGenes$ENSEMBL[which(bLRLikeCriteria==TRUE)]))
#[1] 111
length(unique(F0_Meta_F2_LineageGenes$SYMBOL[which(bLRLikeCriteria==TRUE)]))
#[1] 107 - because the Birt et al. database only includes gene symbols, realistically those are the only genes that could overlap. So we'll be underestimate enrichment a little by using the full #.

bHRLikeCriteria<-(((F0_Meta_F2_LineageGenes$P.value.Total_LocoScore<0.05 & F0_Meta_F2_LineageGenes$Coef.Total_LocoScore>0)
                   |(F0_Meta_F2_LineageGenes$P.value.EPM_Time_Immobile<0.05 & F0_Meta_F2_LineageGenes$Coef.EPM_Time_Immobile<0)
                   |(F0_Meta_F2_LineageGenes$P.value.EPM_DistanceTraveled_noOutlier<0.05 & F0_Meta_F2_LineageGenes$Coef.EPM_DistanceTraveled_noOutlier>0)
                   |(F0_Meta_F2_LineageGenes$P.value.EPM_Percent_Time_Open_Arm<0.05 & F0_Meta_F2_LineageGenes$Coef.EPM_Percent_Time_Open_Arm>0)
                   |(F0_Meta_F2_LineageGenes$P.value.PCA_Index_Days6and.7<0.05 & F0_Meta_F2_LineageGenes$Coef.PCA_Index_Days6and.7>0))
                  &((F0_Meta_F2_LineageGenes$F0_Coef.Lineage_AsFactorbLR<0 & F0_Meta_F2_LineageGenes$F0_P.value.adj.Lineage_AsFactorbLR<0.10)
                    |(F0_Meta_F2_LineageGenes$MetaAnalysis_estimate<0 & F0_Meta_F2_LineageGenes$MetaAnalysis_BH<0.10)
                    |(F0_Meta_F2_LineageGenes$F0_P.value.Lineage_AsFactorbLR<0.05 & F0_Meta_F2_LineageGenes$MetaAnalysis_rawp<0.05 & F0_Meta_F2_LineageGenes$F0_Coef.Lineage_AsFactorbLR<0 & F0_Meta_F2_LineageGenes$MetaAnalysis_estimate<0)))
sum(bHRLikeCriteria, na.rm=TRUE)
#[1] 81
length(unique(F0_Meta_F2_LineageGenes$ENSEMBL[which(bHRLikeCriteria==TRUE)]))
#[1] 81
length(unique(F0_Meta_F2_LineageGenes$SYMBOL[which(bHRLikeCriteria==TRUE)]))
#[1] 79 - because the Birt et al. database only includes gene symbols, realistically those are the only genes that could overlap. So we'll be underestimate enrichment a little by using the full #.


#############################

#bLR genes: 619 (608 of which are in the F2 dataset)
#bLR-like behavior genes: 924
#Intersection: 111
#Total: 13788 w/ F0 results
#Total w/ both F0 and F2 results: 13339

#bLR genes - intersection:
#619-111
#[1] 508

#bLR genes (if only considering the ones in the F2 dataset) - intersection:
608-111
#[1] 497

#bLR-like behavior genes - intersection:
924-111
#[1] 813

#Total-other cells
#13788-111-508-813
#[1] 12356

#CrossTable_bLRGenesVsbLRLike<-rbind(c(111, 508), c(813, 12356))
#fisher.test(CrossTable_bLRGenesVsbLRLike)
# Fisher's Exact Test for Count Data
# 
# data:  CrossTable_bLRGenesVsbLRLike
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.646897 4.138193
# sample estimates:
# odds ratio 
#    3.32036 

#library(VennDiagram)
#grid.newpage() 
#draw.pairwise.venn(619,924,111,col=c("red", "red"), fill=c("pink", "pink"))


#Total-other cells - if only considering genes in F2 dataset (probably more appropriate):
13339-111-497-813
#[1] 11918

CrossTable_bLRGenesVsbLRLike<-rbind(c(111, 497), c(813, 11918))
fisher.test(CrossTable_bLRGenesVsbLRLike)
#Fisher's Exact Test for Count Data
# 
# data:  CrossTable_bLRGenesVsbLRLike
# p-value<2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.608713 4.081286
# sample estimates:
# odds ratio 
#   3.273543

library(VennDiagram)
grid.newpage() 
draw.pairwise.venn(608,924,111,col=c("red", "red"), fill=c("pink", "pink"))

############################
#bHR genes: 444 (437 of which are in the F2 dataset)
#bHR-like behavior genes: 1075
#Intersection: 81
#Total: 13788 w/ F0 results
#Total w/ both F0 and F2 results: 13339

#bHR genes - intersection
#444-81
#[1] 363

#bHR genes - intersection - when only considering genes in F2 dataset:
437-81
#[1] 356

#bHR-like behavior genes - intersection:
1075-81
#[1] 994

#Total-other cells - when only considering genes in the F2 dataset (probably most appropriate):
13339-81-356-994
#11908

CrossTable_bHRGenesVsbHRLike<-rbind(c(81, 356), c(994, 11908))
fisher.test(CrossTable_bHRGenesVsbHRLike)
# Fisher's Exact Test for Count Data
# 
# data:  CrossTable_bHRGenesVsbHRLike
# p-value = 7.141e-13
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.095674 3.510720
# sample estimates:
# odds ratio 
#   2.725419 

library(VennDiagram)
grid.newpage() 
draw.pairwise.venn(437,1075,81,col=c("darkgreen", "darkgreen"), fill=c("lightgreen", "lightgreen"))

#Total-other cells - using the full F0 dataset
# 13788-81-363-994
#12350

#CrossTable_bHRGenesVsbHRLike<-rbind(c(81, 363), c(994, 12350))
#fisher.test(CrossTable_bHRGenesVsbHRLike)
# Fisher's Exact Test for Count Data
# 
# data:  CrossTable_bHRGenesVsbHRLike
# p-value = 3.797e-13
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  2.132598 3.569164
# sample estimates:
# odds ratio 
#   2.772032

# library(VennDiagram)
# grid.newpage() 
# draw.pairwise.venn(444,1075,81,col=c("darkgreen", "darkgreen"), fill=c("lightgreen", "lightgreen"))


###########

#Full number of genes with gene symbols in our DE results (& that also have F2 results):
length(unique(F0_Meta_F2_results_wLineageGenes_forSuppl$SYMBOL.1[is.na(F0_Meta_F2_results_wLineageGenes_forSuppl$P.value.Total_LocoScore)==FALSE]))
#[1] 13301

length(unique(OtherRatModels_Up$GeneSymbol_Rat))
#[1] 1121

#Although I guess the true comparison point would be how many of these genes were in the full dataset...
sum(unique(OtherRatModels_Up$GeneSymbol_Rat)%in%F0_Meta_F2_results_wLineageGenes_forSuppl$SYMBOL.1)
#[1] 849

#Number of lineage genes that are up with bLRs and have bLR-like F2 like behaviors (p<0.05):
length(unique(F0_Meta_F2_LineageGenes$SYMBOL.1[which(bLRLikeCriteria==TRUE)]))
#[1] 111

#Up in bLRs, bLR-like F2 like behaviors (p<0.05), & up with internalizing in other models:
sum(unique(F0_Meta_F2_LineageGenes$SYMBOL.1[which(bLRLikeCriteria==TRUE)])%in%OtherRatModels_Up$GeneSymbol_Rat)
#[1] 16
unique(OtherRatModels_Up$GeneSymbol_Rat[OtherRatModels_Up$GeneSymbol_Rat%in%F0_Meta_F2_LineageGenes$SYMBOL.1[which(bLRLikeCriteria==TRUE)]])
# [1] "Ilf3"    "Serinc2" "Tmie"    "Ghdc"    "Tmem144" "B3gnt9"  "Ist1"    "Nudt4"   "Rarres2" "Bphl"    "Art3"    "Dclk2"   "Sestd1"  "Tmprss5" "Zscan18"
# [16] "Ddc" 

16/111
#[1]  0.1441441
#14%

#"Rarres2", "Serinc2",  "Ilf3" - not in our original list because they were on both upregulated and down-regulated lists

grid.newpage() 
draw.pairwise.venn(849,111, 16,col=c("red", "red"), fill=c("pink", "pink"))
#category=c("Up in Other Rat Models w/ internalizing-like behavior", "Up in bLRs and w/ bLR-like F2 behavior")

#Up in other rat models but not Up in bLRs and w/ bLR-like F2 behavior
849-16
#[1] 833

#Up in bLRs and w/ bLR-like F2 behavior but not Up in other models
111-16
#[1] 95

#Not in either
13301-833-95-16
#[1] 12357
CrossTable_OtherModelsVsbLRLike<-rbind(c(16, 833), c(95, 12357))
fisher.test(CrossTable_OtherModelsVsbLRLike)
# Fisher's Exact Test for Count Data
# 
# data:  CrossTable_OtherModelsVsbLRLike
# p-value = 0.002415
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.366029 4.293179
# sample estimates:
# odds ratio 
#   2.498152 

#Enriched for overlap, but it's a little wussy statistically

length(unique(OtherRatModels_Down$GeneSymbol_Rat))
#[1] 1437

#Although I guess the true comparison point would be how many of these genes were in the full dataset...
sum(unique(OtherRatModels_Down$GeneSymbol_Rat)%in%F0_Meta_F2_results_wLineageGenes_forSuppl$SYMBOL.1)
#[1] 1222

#Number of lineage genes that are up with bHRs and have bHR-like F2 like behaviors (p<0.05):
length(unique(F0_Meta_F2_LineageGenes$SYMBOL.1[which(bHRLikeCriteria==TRUE)]))
#[1] 81

#Up in bHRs, bHR-like F2 like behaviors (p<0.05), & down with internalizing in other models:
sum(unique(F0_Meta_F2_LineageGenes$SYMBOL.1[which(bHRLikeCriteria==TRUE)])%in%OtherRatModels_Down$GeneSymbol_Rat)
#[1] 14

14/81
#[1] 0.1728395
#17%

F0_Meta_F2_LineageGenes$SYMBOL.1[bHRLikeCriteria&(F0_Meta_F2_LineageGenes$SYMBOL.1%in%OtherRatModels_Down$GeneSymbol_Rat)]
#  [1] "Mfge8"      "Ucp2"       "RGD1359508" "Fcrl2"      "Dbt"        "Cav1"       "Ptpro"      "Pramef8"    "Glra2"     
# [10] "Sun2"       "Mxd4"       "Nqo2"       "Hook3"      "Sst"        NA    

#"Sst"   "Ptpro" 
#Not on our previous list  because they were on both upregulated and down-regulated lists
#Add filter to remove this.
# NA  
#Not a gene SYMBOL.1... - but not included in the "14" count.
#I don't see Slx4ip
#It turns out that was because Slx4ip was a typo in the original analysis

grid.newpage() 
draw.pairwise.venn(1222,81, 14,col=c("darkgreen", "darkgreen"), fill=c("lightgreen", "lightgreen"))
#category=c("Up in Other Rat Models w/ less internalizing-like behavior", "Up in bHRs and w/ bHR-like F2 behavior")

#Down in other rat models but not Up in bHRs and w/ bHR-like F2 behavior
1222-14
#[1] 1208

#Up in bHRs and w/ bHR-like F2 behavior but not Up in other models
81-14
#[1] 67

#Not in either
13301-1222-81-14
#[1] 11984
CrossTable_OtherModels_bHRLike<-rbind(c(14, 1208), c(67, 11984))
fisher.test(CrossTable_OtherModels_bHRLike)
# Fisher's Exact Test for Count Data
# 
# data:  CrossTable_OtherModels_bHRLike
# p-value = 0.0189
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.072491 3.738565
# sample estimates:
# odds ratio 
#     2.0728 

#Sig but kind-of wussy.
#Of course, the full pool of genes is overestimated, since many genes in our study weren't measured on these earlier platforms.

#Overlap with a few studies is particularly high:

#Meckes et al.

#Up in bLR/bLR-like: 111
#Total up in Meckes (subset in our data): 342
#Overlap up: 10 (I just counted these from our tables)

#Up in Meckes but not Up in bLRs and w/ bLR-like F2 behavior
342-10
#[1] 332

#Up in bLRs and w/ bLR-like F2 behavior but not Up in Meckes
111-10
#[1] 101

#Not in either
13301-332-101-10
#[1] 12858
CrossTable_MeckesVsbLRLike<-rbind(c(10, 332), c(101, 12858))
fisher.test(CrossTable_MeckesVsbLRLike)
# Fisher's Exact Test for Count Data
# 
# data:  CrossTable_MeckesVsbLRLike
# p-value = 0.0005721
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.768050 7.432642
# sample estimates:
# odds ratio 
#   3.833957 

grid.newpage() 
draw.pairwise.venn(342,111, 10,col=c("red", "red"), fill=c("pink", "pink"))
#category=c("Up in Meckes", "Up in bLRs and w/ bLR-like F2 behavior")



#save.image("~/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/DEResults_F0F2_SummaryTables/CompareF0F2Meta_Workspace_20220807_wGenetics.RData")
save.image("~/Dropbox (University of Michigan)/LaptopBackup_20221123/Microarray Gen/HRLR/NIDA_U01/DEResults_F0F2_SummaryTables/CompareF0F2Meta_Workspace_20240119_wGenetics.RData")
