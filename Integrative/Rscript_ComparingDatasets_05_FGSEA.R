#Comparing RNA-Seq Results across datasets: F0 vs. F2 vs. bHR/bLR Meta-analysis
#05_Running gene set enrichment analysis (fGSEA) using a combined version of our differential expression results from the F0 and F2 datasets
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-08, updated later for a few figures for the paper.



#############

#ok, there seems to be functional patterns, but they're not showing up using traditional tools that compare a single gene set to background using super general gene set databases (GO, reactome)
#Is there some way that we could make a continuous variable that encompasses all of our different sources of evidence and then I could use my own custom gene set database and fGSEA?
#We would need to normalize in some way, since the Log2FC aren't comparable - t-stats?

setwd("~/Documents/Microarray Gen/HRLR/NIDA_U01")

F0_Meta_F2_results<-read.csv("F0_Meta_F2_results_ForSupplTable.csv", header=TRUE, stringsAsFactors = FALSE)
str(F0_Meta_F2_results)
#Ah - no t-stat for the Meta-analysis
#The estimate divided by SE should get at that.

F0_Meta_F2_results$MetaAnalysisTstat<-F0_Meta_F2_results$MetaAnalysis_estimate/F0_Meta_F2_results$MetaAnalysis_SE

#Is there a way to combine scores so that the bHR/bLR values don't completely swamp the F2 values?
#Both averaging and PCA are likely to make it so that the bHR/bLR signature drives things...maybe?

F0_Meta_F2_results_TstatMatrix<-cbind(F0_Meta_F2_results$F0_t.Lineage_AsFactorbLR, F0_Meta_F2_results$MetaAnalysisTstat, F0_Meta_F2_results$t.EPM_Time_Immobile, F0_Meta_F2_results$t.EPM_Percent_Time_Open_Arm, F0_Meta_F2_results$t.EPM_DistanceTraveled_noOutlier, F0_Meta_F2_results$t.Total_LocoScore, F0_Meta_F2_results$t.PCA_Index_Days6and.7)

str(F0_Meta_F2_results_TstatMatrix)
row.names(F0_Meta_F2_results_TstatMatrix)<-F0_Meta_F2_results$ENSEMBL
colnames(F0_Meta_F2_results_TstatMatrix)<-c("F0_bLRvsbHR_Tstat", "Meta_bLRvsbHR_Tstat", "F2_TimeImmobile_Tstat", "F2_PercentTimeOpenArm_Tstat", "F2_EPMDistance_Tstat", "F2_LocoScore_Tstat", "F2_PCA_Tstat")

heatmap(cor(F0_Meta_F2_results_TstatMatrix, use="pairwise.complete.obs"))

hist(F0_Meta_F2_results_TstatMatrix[,1])
#Mostly between-6 to 6, but there are a few as low as -10 or as high as 25
hist(F0_Meta_F2_results_TstatMatrix[,2])
#Mostly between-5 to 5
hist(F0_Meta_F2_results_TstatMatrix[,3])
#Mostly between-3 to 3
hist(F0_Meta_F2_results_TstatMatrix[,4])
hist(F0_Meta_F2_results_TstatMatrix[,5])
hist(F0_Meta_F2_results_TstatMatrix[,6])
hist(F0_Meta_F2_results_TstatMatrix[,7])
#ditto for all F2

F0_Meta_F2_results_TstatMatrix_DirectionConsistent<-cbind(F0_Meta_F2_results_TstatMatrix[,c(1:3)], -F0_Meta_F2_results_TstatMatrix[,c(4:7)])

heatmap(cor(F0_Meta_F2_results_TstatMatrix_DirectionConsistent, use="pairwise.complete.obs"))

#Let me make sure some of these aren't lacking F2 data
sum(is.na(F0_Meta_F2_results_TstatMatrix_DirectionConsistent[,3]))
#[1] 479

F0_Meta_F2_results_TstatMatrix_DirectionConsistent_wF2<-F0_Meta_F2_results_TstatMatrix_DirectionConsistent[is.na(F0_Meta_F2_results_TstatMatrix_DirectionConsistent[,3])==FALSE,]

colnames(F0_Meta_F2_results)
F0_Meta_F2_results_TstatMatrix_wF2_Annotation<-F0_Meta_F2_results[is.na(F0_Meta_F2_results_TstatMatrix_DirectionConsistent[,3])==FALSE,c(2, 8, 17, 38:46)]

F0_Meta_F2_results_AverageTstat<-apply(F0_Meta_F2_results_TstatMatrix_DirectionConsistent_wF2, 1, function(y) mean(y, na.rm=TRUE))

hist(F0_Meta_F2_results_AverageTstat)
#Mostly falls between-2 and 2

sum(is.na(F0_Meta_F2_results_AverageTstat))
#[1] 0
length(F0_Meta_F2_results_AverageTstat)
#[1] 13669

F0_Meta_F2_results_TstatMatrix_wF2_Annotation$Symbol[order(F0_Meta_F2_results_AverageTstat)][c(1:40)]
# [1] NA             "Mfge8"        "Spg7"         "Pkib"         "RGD1359508"   "Oard1"        "Sun2"        
# [8] "Afg3l1"       "Slx4ip"       "Hmg20a"       "Aldh1a1"      "Tes"          "Dbt"          "Ucp2"        
# [15] "Asb15"        "Tmco5a"       "Cav1"         "Ttc30a1"      "Ttc30a1"      "Ier2"         "C2cd3"       
# [22] "Ptpro"        "Gfpt2"        NA             NA             NA             "Fcrl2"        "LOC100911320"
# [29] "Dpysl3"       "Rfk"          NA             "Maff"         NA             "Pxmp4"        "Glra2"       
# [36] "Slc3a2"       "Siglec5"      NA             "Cyct"         "Cyct"        

F0_Meta_F2_results_TstatMatrix_wF2_Annotation$Symbol[order(F0_Meta_F2_results_AverageTstat)][c(13629:13669)]
# [1] "Prss55"       "Fcgr3a"       "Cd200r1"      "Ube2l6"       "Acss2"        "C1qc"         NA            
# [8] "Unc45a"       "LOC108348293" "Rarres2"      NA             "Fanci"        "Herc6"        "Tmprss5"     
# [15] "Gss"          "Idh1"         NA             "C1qtnf1"      "Nudt4"        "C1qa"         "Commd9"      
# [22] "LOC689130"    NA             "C1qb"         "Art3"         "Pex11a"       "Vps9d1"       "Tmem144"     
# [29] "Ddc"          "LOC690414"    NA             "Tmem176a"     "Plekhb1"      "Wdr93"        "Bloc1s6"     
# [36] "LOC100361008" NA             "Fzd6"         "Ghdc"         NA             NA        


F0_Meta_F2_results_TstatMatrix_DirectionConsistent_wF2[which(F0_Meta_F2_results_TstatMatrix_wF2_Annotation$Symbol=="Pkib"),]

F0_Meta_F2_results_TstatMatrix_DirectionConsistent_wF2[which(F0_Meta_F2_results_TstatMatrix_wF2_Annotation$Symbol=="Oard1"),]

F0_Meta_F2_results_TstatMatrix_DirectionConsistent_wF2[which(F0_Meta_F2_results_TstatMatrix_wF2_Annotation$Symbol=="Hmg20a"),]

F0_Meta_F2_results_TstatMatrix_DirectionConsistent_wF2[which(F0_Meta_F2_results_TstatMatrix_wF2_Annotation$Symbol=="Prss55"),]

F0_Meta_F2_results_TstatMatrix_DirectionConsistent_wF2[which(F0_Meta_F2_results_TstatMatrix_wF2_Annotation$Symbol=="Fcgr3a"),]

F0_Meta_F2_results_TstatMatrix_DirectionConsistent_wF2[which(F0_Meta_F2_results_TstatMatrix_wF2_Annotation$Symbol=="Cd200r1"),]

#seems reasonably decent - top genes seem to have effects across datasets, although many wouldn't qualify as significant.

write.csv(data.frame(F0_Meta_F2_results_TstatMatrix_wF2_Annotation, F0_Meta_F2_results_AverageTstat, F0_Meta_F2_results_TstatMatrix_DirectionConsistent_wF2), "F0_Meta_F2_results_TstatMatrix_wAnnotation.csv")

###############

#fGSEA 

library(fgsea)

sum(duplicated(F0_Meta_F2_results_TstatMatrix_wF2_Annotation$ENSEMBL))
#[1] 330
sum(duplicated(F0_Meta_F2_results_TstatMatrix_wF2_Annotation$Symbol))
#[1] 995
sum(is.na(F0_Meta_F2_results_TstatMatrix_wF2_Annotation$Symbol))
#[1] 658


library(dplyr)

Temp<-data.frame(ENSEMBL=F0_Meta_F2_results_TstatMatrix_wF2_Annotation$ENSEMBL, Symbol=F0_Meta_F2_results_TstatMatrix_wF2_Annotation$Symbol, Tstat=F0_Meta_F2_results_AverageTstat, stringsAsFactors = FALSE)
Temp_NoNA<-Temp[is.na(Temp$Symbol)==FALSE,]

sum(duplicated(Temp_NoNA$ENSEMBL))
#[1] 316
sum(duplicated(Temp_NoNA$Symbol))
#[1] 338
#Pretty close.
sum(duplicated(Temp_NoNA$Tstat))
#[1] 316

338-316
F0_Meta_F2_results_AverageTstat_forGSEA<-distinct(Temp_NoNA, ENSEMBL,.keep_all=TRUE)
str(F0_Meta_F2_results_AverageTstat_forGSEA)
# 'data.frame':	12695 obs. of  3 variables:
# $ ENSEMBL: chr  "ENSRNOG00000014303" "ENSRNOG00000014330" "ENSRNOG00000049505" "ENSRNOG00000014916" ...
# $ Symbol : chr  "Lrp11" "Pcmt1" "Nup43" "Lats1" ...
# $ Tstat  : num  0.293 0.474 0.603 0.487 -0.17 ...

F0_Meta_F2_results_AverageTstat_forGSEA[duplicated(F0_Meta_F2_results_AverageTstat_forGSEA$Symbol),]
#Several ENSEMBL genes map to the same symbol

F0_Meta_F2_results_AverageTstat_forGSEA_NoDups<-tapply(X=F0_Meta_F2_results_AverageTstat_forGSEA$Tstat, INDEX=F0_Meta_F2_results_AverageTstat_forGSEA$Symbol, FUN=mean)

str(F0_Meta_F2_results_AverageTstat_forGSEA_NoDups)
# num [1:12673(1d)] 0.838 -1.463 0.481 1.267 -0.578 ...
# - attr(*, "dimnames")=List of 1
# ..$ : chr [1:12673] "A3galt2" "Aaas" "Aacs" "Aadat" ...

names(F0_Meta_F2_results_AverageTstat_forGSEA_NoDups)<-names(table(F0_Meta_F2_results_AverageTstat_forGSEA$Symbol))

F0_Meta_F2_results_AverageTstat_forGSEA_Ranked<-F0_Meta_F2_results_AverageTstat_forGSEA_NoDups[order(F0_Meta_F2_results_AverageTstat_forGSEA_NoDups)]

head(F0_Meta_F2_results_AverageTstat_forGSEA_Ranked)
# Mfge8       Spg7       Pkib RGD1359508      Oard1       Sun2 
# -2.487051  -2.477041  -2.474786  -2.382322  -2.376376  -2.318062 

write.csv(F0_Meta_F2_results_AverageTstat_forGSEA_Ranked, "F0_Meta_F2_results_AverageTstat_forGSEA_Ranked.csv")

###########

GMT_ForRats<-gmtPathways("c5withBrainCellTypesFunctionGemmaGeneWeaver_RatOrtholog_HC_ForElaine.txt")
str(GMT_ForRats)
# 
# List of 15866
# $ GOBP_MITOCHONDRIAL_GENOME_MAINTENANCE                                                                                                                                                                                                                                                                                                       : chr [1:2082] "Akt3" "Ppargc1a" "Polg2" "Parp1" ...
# $ GOBP_REPRODUCTION                                                                                                                                                                                                                                                                                                                           : chr [1:2082] "Ada" "Gnpda1" "Zglp1" "Syce1l" ...
# $ GOBP_SINGLE_STRAND_BREAK_REPAIR                                                                                                                                                                                                                                                                                                             : chr [1:2082] "Ercc8" "Parp1" "Aplf" "Ercc6" ...

#############

temp1<-fgsea(GMT_ForRats, F0_Meta_F2_results_AverageTstat_forGSEA_Ranked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "fGSEA_F0_Meta_F2_results_AverageTstat_GMTc5withBrainCellTypesFunctionGemmaGeneWeaver_RatOrtholog_HC_F.csv")

##############

F0_Meta_F2_results_AverageTstat_ABS_forGSEA_Ranked<-abs(F0_Meta_F2_results_AverageTstat_forGSEA_NoDups)[order(abs(F0_Meta_F2_results_AverageTstat_forGSEA_NoDups))]

head(F0_Meta_F2_results_AverageTstat_ABS_forGSEA_Ranked)
# Ccdc66         Vrk1       Unc119          Ipp       Nt5dc2          Ca4 
# 6.849115e-06 9.037175e-05 9.502621e-05 1.832394e-04 3.205448e-04 4.728053e-04 
tail(F0_Meta_F2_results_AverageTstat_ABS_forGSEA_Ranked)
# Plekhb1        Wdr93      Bloc1s6 LOC100361008         Fzd6         Ghdc 
# 2.562954     2.628368     2.672913     2.704944     2.809315     2.923333 

write.csv(F0_Meta_F2_results_AverageTstat_ABS_forGSEA_Ranked, "F0_Meta_F2_results_AverageTstat_ABS_forGSEA_Ranked.csv")

rm(temp1)

temp1<-fgsea(GMT_ForRats, F0_Meta_F2_results_AverageTstat_ABS_forGSEA_Ranked, nperm=10000, minSize = 10, maxSize = 1000)
str(temp1)

temp1$leadingEdge<-vapply(temp1$leadingEdge, paste, collapse= ",", character(1L))

write.csv(temp1, "fGSEA_F0_Meta_F2_results_AverageTstat_ABS_GMTc5withBrainCellTypesFunctionGemmaGeneWeaver_RatOrtholog_HC_F.csv")

