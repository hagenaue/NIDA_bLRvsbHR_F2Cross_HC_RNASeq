#F2 HC RNA-Seq Dataset
#18_Comparing our results to meta-analysis HR/LR results (maybe old) and outputting graphs for top HRLR DE genes
#Megan Hagenauer and Elaine Hebda-Bauer
#2022-01-04, updated later for a few figures for the paper.


#Comparison with meta-analysis output - I'm not sure if this section has been updated.
#Also: I ran other comparisons in a separate code document.

MetaAnalysisOutput_JustLateGen<-read.csv("metaOutput_JustLateGen_FDR.csv", header=TRUE, stringsAsFactors = FALSE)
str(MetaAnalysisOutput_JustLateGen)

# 'data.frame':	13319 obs. of  11 variables:
#   $ X         : int  19 20 21 22 23 24 25 26 27 28 ...
# $ GeneSymbol: chr  "A2m" "A3galt2" "Aaas" "Aacs" ...
# $ estimate  : num  0.482 -0.486 0.817 -1.048 -0.436 ...
# $ SE        : num  0.438 0.441 0.443 0.467 0.457 ...
# $ pval      : num  0.2716 0.2713 0.0649 0.0249 0.3399 ...
# $ CI_lb     : num  -0.3773 -1.3509 -0.0504 -1.9638 -1.3312 ...
# $ CI_ub     : num  1.341 0.38 1.685 -0.132 0.459 ...
# $ datasets  : int  2 2 2 2 2 2 2 2 2 2 ...
# $ rawp      : num  0.2716 0.2713 0.0649 0.0249 0.3399 ...
# $ BH        : num  0.599 0.599 0.312 0.197 0.659 ...
# $ BY        : num  1 1 1 1 1 1 1 1 1 1 ..


colnames(MetaAnalysisOutput_JustLateGen)[2]<-"SYMBOL"


F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen<-join(F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped, MetaAnalysisOutput_JustLateGen, by="SYMBOL", type="left")

str(F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen)

write.csv(F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen, "F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen.csv")



#Ah - the estimates are reversed (bLR = reference group) - I just double-checked C1qa and Bmp4 in the output

F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$estimate<-(F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$estimate)*-1

F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$CI_lb<-(F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$CI_lb)*-1

F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$CI_ub<-(F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$CI_ub)*-1


write.csv(F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen, "F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen.csv")

pdf("U01_F2_Locomotor_RNASeq_Vs_LateGeneMetaAnalysis.pdf", width=6, height=6)
plot(Coef.Total_LocoScore~estimate, data=F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen, ylab="U01 RNA-Seq: F2 Total LocoScore Log2FC", xlab="Meta-Analysis of F37 & F43 RNA-Seq: bLR vs. bHR Estimated Effect Size")
BestFitLine<-lm(Coef.Total_LocoScore~estimate, data=F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen)
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(BestFitLine)

# Call:
#   lm(formula = Coef.Total_LocoScore ~ estimate, data = F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -1.132e-03 -2.453e-05  3.740e-06  2.942e-05  8.022e-04 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  2.885e-07  5.810e-07   0.497     0.62    
# estimate    -9.868e-06  8.013e-07 -12.315   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 6.237e-05 on 11551 degrees of freedom
# (2508 observations deleted due to missingness)
# Multiple R-squared:  0.01296,	Adjusted R-squared:  0.01287 
# F-statistic: 151.7 on 1 and 11551 DF,  p-value: < 2.2e-16


pdf("U01_F2_Locomotor_RNASeq_Vs_LateGeneMetaAnalysis_bLRvsbHR_Tstat.pdf", width=6, height=6)
plot(Coef.Total_LocoScore~estimate, data=F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen, ylab="U01 RNA-Seq: F2 Total LocoScore Tstat", xlab="Meta-Analysis of F37 & F43 RNA-Seq: bLR vs. bHR Estimated Effect Size")
BestFitLine<-lm(Coef.Total_LocoScore~estimate, data=F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen)
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(BestFitLine)
# 
# Call:
#   lm(formula = t.Total_LocoScore ~ estimate, data = F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.8610 -0.7066  0.0354  0.7106  4.1148 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.081825   0.009434   8.673   <2e-16 ***
#   estimate    -0.194374   0.013011 -14.939   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 1.013 on 11551 degrees of freedom
# (2508 observations deleted due to missingness)
# Multiple R-squared:  0.01895,	Adjusted R-squared:  0.01887 
# F-statistic: 223.2 on 1 and 11551 DF,  p-value: < 2.2e-16


pdf("U01_F2_Locomotor_RNASeq__Vs_LateGeneMetaAnalysis_bLRvsbHR_FDR10Meta.pdf", width=6, height=6)
plot(Coef.Total_LocoScore~estimate, data=F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen[F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH<0.1,], ylab="U01 RNA-Seq: F2 Total LocoScore Log2FC", xlab="Meta-Analysis of F37 & F43 RNA-Seq: bLR vs. bHR Estimated Effect Size")
BestFitLine<-lm(Coef.Total_LocoScore~estimate, data=F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen[F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH<0.1,])
abline(BestFitLine, col=2, lwd=3)
dev.off()

summary.lm(BestFitLine)

# Call:
#   lm(formula = Coef.Total_LocoScore ~ estimate, data = F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen[F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH < 
#                                                                                                                                               0.1, ])
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -4.172e-04 -2.794e-05  1.400e-06  3.071e-05  7.765e-04 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.477e-06  2.628e-06   0.562    0.574    
# estimate    -1.361e-05  1.421e-06  -9.583   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 7.287e-05 on 784 degrees of freedom
# (2508 observations deleted due to missingness)
# Multiple R-squared:  0.1049,	Adjusted R-squared:  0.1037 
# F-statistic: 91.84 on 1 and 784 DF,  p-value: < 2.2e-16


pdf("U01_F2_Locomotor_RNASeq_Vs_LateGeneMetaAnalysis_bLRvsbHR_Tstat_FDR10Meta.pdf", width=6, height=6)
plot(t.Total_LocoScore~estimate, data=F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen[F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH<0.1,], ylab="U01 RNA-Seq: F2 Total LocoScore Tstat", xlab="Meta-Analysis of F37 & F43 RNA-Seq: bLR vs. bHR Estimated Effect Size")
BestFitLine<-lm(t.Total_LocoScore~estimate, data=F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen[F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH<0.1,])
abline(BestFitLine, col=2, lwd=3)
abline(a=0, b=0)
dev.off()

summary.lm(BestFitLine)

# Call:
#   lm(formula = t.Total_LocoScore ~ estimate, data = F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen[F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH < 
#                                                                                                                                            0.1, ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.4679 -0.7423  0.0501  0.7629  3.6020 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.03080    0.03921   0.786    0.432    
# estimate    -0.23927    0.02120 -11.287   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 1.087 on 784 degrees of freedom
# (2508 observations deleted due to missingness)
# Multiple R-squared:  0.1398,	Adjusted R-squared:  0.1387 
# F-statistic: 127.4 on 1 and 784 DF,  p-value: < 2.2e-16


pdf("U01_F2_Locomotor_RNASeq_Vs_LateGeneMetaAnalysis_bLRvsbHR_Tstat_FDR05Meta.pdf", width=6, height=6)
plot(t.Total_LocoScore~estimate, data=F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen[F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH<0.05,], ylab="U01 RNA-Seq: F2 Total LocoScore Tstat", xlab="Meta-Analysis of F37 & F43 RNA-Seq: bLR vs. bHR Estimated Effect Size")
BestFitLine<-lm(t.Total_LocoScore~estimate, data=F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen[F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH<0.05,])
abline(BestFitLine, col=2, lwd=3)
abline(a=0, b=0)
dev.off()

summary.lm(BestFitLine)

# Call:
#   lm(formula = t.Total_LocoScore ~ estimate, data = F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen[F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH < 
#                                                                                                                                            0.05, ])
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -3.2444 -0.7795  0.0652  0.7691  3.1871 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.02053    0.05241   0.392    0.696    
# estimate    -0.22402    0.02454  -9.129   <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Residual standard error: 1.086 on 440 degrees of freedom
# (2508 observations deleted due to missingness)
# Multiple R-squared:  0.1592,	Adjusted R-squared:  0.1573 
# F-statistic: 83.33 on 1 and 440 DF,  p-value: < 2.2e-16


#Write out genes that previously were associated with bHR/bLR in the late generation meta-analysis and are also now showing nominal associations with locomotor score:

GenesWithFDR05inLateGenMeta_AndF2NominallySigLocomotorEffects<-F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen[is.na(F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH)==FALSE & F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH<0.05 & F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$P.value.Total_LocoScore<0.05,]

str(GenesWithFDR05inLateGenMeta_AndF2NominallySigLocomotorEffects)

#'data.frame':	42 obs. of  59 variables:
#'
sum(is.na(F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH)==FALSE & F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH<0.05)
#[1] 442

GenesWithFDR05inLateGenMeta_AndF2NominallySigLocomotorEffects$SYMBOL
# [1] "Mfge8"    "Pex11a"   "Unc45a"   "C2cd3"    "Ucp2"     "Plekhb1"  "Aldh1a1"  "Lipa"     "Tes"      "Cav1"     "Tmem176a" "Cycs"     "Tex261"   "Txnrd3"   "Fkbp4"    "Sp3"      "Ttc30a1"  "Slx4ip"   "Tek"     
# [20] "Serinc2"  "C1qc"     "Pramef8"  "Akap8l"   "Cyp4f4"   "Nudt4"    "Sun2"     "Acad11"   "Gls"      "Ddc"      "Cyb5r1"   "Tmem101"  "Fn3krp"   "Prss55"   "Nqo2"     "Tmem161a" "Mau2"     "Camk4"    "Rpl17"   
# [39] "Spg7"     "Vps9d1"   "Eif4ebp2" "Ift81"  

42/442
#[1] 0.09502262

#So almost twice as many nominally significant associations with locomotor score than we would expect by random chance.


sum(is.na(F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH)==FALSE & F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH<0.05 & F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$P.value.Total_LocoScore<0.005)
#[1] 12

12/442
#[1] 0.02714932

0.005*442
#[1] 2.21

12/2.21
#5.429864

F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$SYMBOL[is.na(F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH)==FALSE & F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$BH<0.05 & F2_Limma_Results_EnsemblVsGeneSymbol_NoMultimapped_vs_MetaAnalysisOutput_JustLateGen$P.value.Total_LocoScore<0.005]
# [1] "Mfge8"    "Pex11a"   "C2cd3"    "Ucp2"     "Plekhb1"  "Tmem176a" "Ttc30a1"  "Serinc2"  "Nudt4"    "Nqo2"     "Spg7"     "Vps9d1"  

#Includes genes we already identified as overlapping with locomotor QTLs in the past, including Mfge8 (QTL: RNor5 1:19-142MB), C2cd3 & Ucp2 (QTL: RNor5: 1:156-171), Ttc30a1 (QTL: RNor 3: 43-65)

#########################

#To make basic graphs of top bHR/bLR genes based on the previous meta-analysis:

setwd("//n05-corea-cifs.umhs.med.umich.edu/Home1/Users/hebda/My Documents/F2_20211207/HRLR_Scatterplots")

#Double check that the order of subjects in F2_PCA_wMetaData is the same as in F2_Annotated_lcpm_NoLowHits

cbind(F2_PCA_wMetaData$row_name, colnames(F2_Annotated_lcpm_NoLowHits)[c(1:245)])
#same order - yeah!
# [,1]       [,2]      
# [1,] "sl469967" "SL469967"
# [2,] "sl469968" "SL469968"
# [3,] "sl469969" "SL469969"
# [4,] "sl469970" "SL469970"
# [5,] "sl469971" "SL469971"
# [6,] "sl469972" "SL469972"
# [7,] "sl469973" "SL469973"
# [8,] "sl469974" "SL469974"
# [9,] "sl469975" "SL469975"
# [10,] "sl469976" "SL469976"

# looping graphing over the top meta-analysis genes:

setwd("//n05-corea-cifs.umhs.med.umich.edu/Home1/Users/hebda/My Documents/F2_20211207")

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

setwd("//n05-corea-cifs.umhs.med.umich.edu/Home1/Users/hebda/My Documents/F2_20211207/HRLR_Scatterplots")

str(F2_Annotated_lcpm_NoLowHits)
colnames(F2_Annotated_lcpm_NoLowHits)

for(i in c(1:nrow(TopMetaAnalysisGenes_Results))){
  
  if(length(which(F2_Annotated_lcpm_NoLowHits$SYMBOL==TopMetaAnalysisGenes_Results$GeneSymbol[i]))>0){
    
    Temp_Data<-data.frame(y=as.matrix(F2_Annotated_lcpm_NoLowHits[which(F2_Annotated_lcpm_NoLowHits$SYMBOL==TopMetaAnalysisGenes_Results$GeneSymbol[i]), c(1:245)])[1,], F2_PCA_wMetaData) 
    
    pdf(paste("F2_",TopMetaAnalysisGenes_Results$GeneSymbol[i], "_vs_LocoScore_bySex.pdf", sep=""), width=6.5, height=6)
    plot(y~Total_LocoScore, data=Temp_Data, col=Sex_AsFactor)
    BestFitLine<-lm(y~Total_LocoScore, data=Temp_Data) 
    abline(BestFitLine, lwd=3)
    BestFitLineM<-lm(y~Total_LocoScore, data=Temp_Data[Temp_Data$Sex_AsFactor=="Male",]) 
    abline(BestFitLineM)
    BestFitLineF<-lm(y~Total_LocoScore, data=Temp_Data[Temp_Data$Sex_AsFactor=="Female",]) 
    abline(BestFitLineF, col=2)
    dev.off()
    
    rm(Temp_Data)
    
  }else{}
  
}

#Before trying the loop, I output one iteration by first running:
i<-1
#Then running the code inside the loop


#Notes from glancing at figures:
#1) the data seems super noisy - not many genes showing much relationship with Locoscore - this seems unlikely
#2) Some massive bimodality happening in the expression data for some genes (e.g., Asb15, chd1l, slc39a12, Bmp4) - dissector?

#Clearly we should see what we can do about noise.

########################################1/31/2022


#TopHRLRF2Genes_ForGraphing<-read.csv("TopHRLRF2Genes_ForGraphing.csv", header=TRUE, stringsAsFactors = FALSE)

TopHRLRF2Genes_ForGraphing<-read.csv("Top_HRLR_F2_Genes_forGraphing.csv", header=TRUE, stringsAsFactors = FALSE)

colnames(TopHRLRF2Genes_ForGraphing)[1]<-"ENSEMBL"
colnames(TopHRLRF2Genes_ForGraphing)[2]<-"Symbol"


setwd("//n05-corea-cifs.umhs.med.umich.edu/Home1/Users/hebda/My Documents/F2_TMMnorm_20220104updated/Plots_20220131_Loco")
list.files()

#double check colnames of (F2_Annotated_lcpm_NoLowHits)

for(i in c(1:nrow(TopHRLRF2Genes_ForGraphing))){
  
  if(length(which(F2_Annotated_lcpm_NoLowHits$ENSEMBL==TopHRLRF2Genes_ForGraphing$ENSEMBL[i]))>0){
    
    Temp_Data<-data.frame(y=as.matrix(F2_Annotated_lcpm_NoLowHits[which(F2_Annotated_lcpm_NoLowHits$ENSEMBL==TopHRLRF2Genes_ForGraphing$ENSEMBL[i]), c(1:245)])[1,], F2_PCA_wMetaData) 
    
    pdf(paste("F2_",TopHRLRF2Genes_ForGraphing$Symbol[i], "_", TopHRLRF2Genes_ForGraphing$ENSEMBL[i], "_vs_LocoScore_bySex.pdf", sep=""), width=6.5, height=6)
    plot(y~Total_LocoScore, data=Temp_Data, ylab=paste(TopHRLRF2Genes_ForGraphing$Symbol[i], ": Log2 CPM", sep=""), col=Sex_AsFactor)
    BestFitLine<-lm(y~Total_LocoScore, data=Temp_Data) 
    abline(BestFitLine, lwd=3)
    BestFitLineM<-lm(y~Total_LocoScore, data=Temp_Data[Temp_Data$Sex_AsFactor=="Male",]) 
    abline(BestFitLineM)
    BestFitLineF<-lm(y~Total_LocoScore, data=Temp_Data[Temp_Data$Sex_AsFactor=="Female",]) 
    abline(BestFitLineF, col=2)
    dev.off()
    
    rm(Temp_Data)
    
  }else{}
  
}


setwd("//n05-corea-cifs.umhs.med.umich.edu/Home1/Users/hebda/My Documents/F2_TMMnorm_20220104updated/Plots_20220131_PercTimeOA")


# TopHRLRF2Genes_ForGraphing<-read.csv("TopHRLRF2Genes_ForGraphing.csv", header=TRUE, stringsAsFactors = FALSE)
# 
# 
# colnames(TopHRLRF2Genes_ForGraphing)[1]<-"ENSEMBL"
# 
# colnames(F2_PCA_wMetaData)

#double check colnames of (F2_Annotated_lcpm_NoLowHits)

for(i in c(1:nrow(TopHRLRF2Genes_ForGraphing))){
  
  if(length(which(F2_Annotated_lcpm_NoLowHits$ENSEMBL==TopHRLRF2Genes_ForGraphing$ENSEMBL[i]))>0){
    
    Temp_Data<-data.frame(y=as.matrix(F2_Annotated_lcpm_NoLowHits[which(F2_Annotated_lcpm_NoLowHits$ENSEMBL==TopHRLRF2Genes_ForGraphing$ENSEMBL[i]), c(1:245)])[1,], F2_PCA_wMetaData) 
    
    pdf(paste("F2_",TopHRLRF2Genes_ForGraphing$Symbol[i], "_", TopHRLRF2Genes_ForGraphing$ENSEMBL[i], "_vs_EPM_Percent_Time_Open_Arm_bySex.pdf", sep=""), width=6.5, height=6)
    plot(y~EPM_Percent_Time_Open_Arm, data=Temp_Data, ylab=paste(TopHRLRF2Genes_ForGraphing$Symbol[i], ": Log2 CPM", sep=""), col=Sex_AsFactor)
    BestFitLine<-lm(y~EPM_Percent_Time_Open_Arm, data=Temp_Data) 
    abline(BestFitLine, lwd=3)
    BestFitLineM<-lm(y~EPM_Percent_Time_Open_Arm, data=Temp_Data[Temp_Data$Sex_AsFactor=="Male",]) 
    abline(BestFitLineM)
    BestFitLineF<-lm(y~EPM_Percent_Time_Open_Arm, data=Temp_Data[Temp_Data$Sex_AsFactor=="Female",]) 
    abline(BestFitLineF, col=2)
    dev.off()
    
    rm(Temp_Data)
    
  }else{}
  
}


setwd("//n05-corea-cifs.umhs.med.umich.edu/Home1/Users/hebda/My Documents/F2_TMMnorm_20220104updated/Plots_20220509_TimeImmobile")


# TopHRLRF2Genes_ForGraphing<-read.csv("TopHRLRF2Genes_ForGraphing.csv", header=TRUE, stringsAsFactors = FALSE)
# 
# 
# colnames(TopHRLRF2Genes_ForGraphing)[1]<-"ENSEMBL"
# 
# colnames(F2_PCA_wMetaData)

#double check colnames of (F2_Annotated_lcpm_NoLowHits)

for(i in c(1:nrow(TopHRLRF2Genes_ForGraphing))){
  
  if(length(which(F2_Annotated_lcpm_NoLowHits$ENSEMBL==TopHRLRF2Genes_ForGraphing$ENSEMBL[i]))>0){
    
    Temp_Data<-data.frame(y=as.matrix(F2_Annotated_lcpm_NoLowHits[which(F2_Annotated_lcpm_NoLowHits$ENSEMBL==TopHRLRF2Genes_ForGraphing$ENSEMBL[i]), c(1:245)])[1,], F2_PCA_wMetaData) 
    
    pdf(paste("F2_",TopHRLRF2Genes_ForGraphing$Symbol[i], "_", TopHRLRF2Genes_ForGraphing$ENSEMBL[i], "_vs_EPM_Percent_Time_Immobile_bySex.pdf", sep=""), width=6.5, height=6)
    plot(y~EPM_Time_Immobile , data=Temp_Data, ylab=paste(TopHRLRF2Genes_ForGraphing$Symbol[i], ": Log2 CPM", sep=""), col=Sex_AsFactor)
    BestFitLine<-lm(y~EPM_Time_Immobile, data=Temp_Data) 
    abline(BestFitLine, lwd=3)
    BestFitLineM<-lm(y~EPM_Time_Immobile, data=Temp_Data[Temp_Data$Sex_AsFactor=="Male",]) 
    abline(BestFitLineM)
    BestFitLineF<-lm(y~EPM_Time_Immobile, data=Temp_Data[Temp_Data$Sex_AsFactor=="Female",]) 
    abline(BestFitLineF, col=2)
    dev.off()
    
    rm(Temp_Data)
    
  }else{}
  
}



#######################################

setwd("//n05-corea-cifs.umhs.med.umich.edu/Home1/Users/hebda/My Documents/F2_TMMnorm_20220104updated/Plots_20220131_EPMdistanceTraveled")


# TopHRLRF2Genes_ForGraphing<-read.csv("TopHRLRF2Genes_ForGraphing.csv", header=TRUE, stringsAsFactors = FALSE)
# 
# 
# colnames(TopHRLRF2Genes_ForGraphing)[1]<-"ENSEMBL"
# 
# colnames(F2_PCA_wMetaData)

#double check colnames of (F2_Annotated_lcpm_NoLowHits)

for(i in c(1:nrow(TopHRLRF2Genes_ForGraphing))){
  
  if(length(which(F2_Annotated_lcpm_NoLowHits$ENSEMBL==TopHRLRF2Genes_ForGraphing$ENSEMBL[i]))>0){
    
    Temp_Data<-data.frame(y=as.matrix(F2_Annotated_lcpm_NoLowHits[which(F2_Annotated_lcpm_NoLowHits$ENSEMBL==TopHRLRF2Genes_ForGraphing$ENSEMBL[i]), c(1:245)])[1,], F2_PCA_wMetaData) 
    
    pdf(paste("F2_",TopHRLRF2Genes_ForGraphing$Symbol[i], "_", TopHRLRF2Genes_ForGraphing$ENSEMBL[i], "_vs_EPM_DistanceTraveled_bySex.pdf", sep=""), width=6.5, height=6)
    plot(y~EPM_DistanceTraveled, data=Temp_Data, ylab=paste(TopHRLRF2Genes_ForGraphing$Symbol[i], ": Log2 CPM", sep=""), col=Sex_AsFactor)
    BestFitLine<-lm(y~EPM_DistanceTraveled, data=Temp_Data) 
    abline(BestFitLine, lwd=3)
    BestFitLineM<-lm(y~EPM_DistanceTraveled, data=Temp_Data[Temp_Data$Sex_AsFactor=="Male",]) 
    abline(BestFitLineM)
    BestFitLineF<-lm(y~EPM_DistanceTraveled, data=Temp_Data[Temp_Data$Sex_AsFactor=="Female",]) 
    abline(BestFitLineF, col=2)
    dev.off()
    
    rm(Temp_Data)
    
  }else{}
  
}


setwd("//n05-corea-cifs.umhs.med.umich.edu/Home1/Users/hebda/My Documents/F2_TMMnorm_20220104updated/Plots_20220509_EPMdistanceTraveled_noOutlier")


# TopHRLRF2Genes_ForGraphing<-read.csv("TopHRLRF2Genes_ForGraphing.csv", header=TRUE, stringsAsFactors = FALSE)
# 
# 
# colnames(TopHRLRF2Genes_ForGraphing)[1]<-"ENSEMBL"
# 
# colnames(F2_PCA_wMetaData)

#double check colnames of (F2_Annotated_lcpm_NoLowHits)

F2_PCA_wMetaData$EPM_DistanceTraveled>5000
F2_PCA_wMetaData$EPM_DistanceTraveled[18]
#[1] 5208.82 - this is the outlier

for(i in c(1:nrow(TopHRLRF2Genes_ForGraphing))){
  
  if(length(which(F2_Annotated_lcpm_NoLowHits$ENSEMBL==TopHRLRF2Genes_ForGraphing$ENSEMBL[i]))>0){
    
    Temp_Data<-data.frame(y=as.matrix(F2_Annotated_lcpm_NoLowHits[which(F2_Annotated_lcpm_NoLowHits$ENSEMBL==TopHRLRF2Genes_ForGraphing$ENSEMBL[i]), c(1:17, 19:245)])[1,], F2_PCA_wMetaData[c(1:17, 19:245),]) 
    
    pdf(paste("F2_",TopHRLRF2Genes_ForGraphing$Symbol[i], "_", TopHRLRF2Genes_ForGraphing$ENSEMBL[i], "_vs_EPM_DistanceTraveled_bySex.pdf", sep=""), width=6.5, height=6)
    plot(y~EPM_DistanceTraveled, data=Temp_Data, ylab=paste(TopHRLRF2Genes_ForGraphing$Symbol[i], ": Log2 CPM", sep=""), col=Sex_AsFactor)
    BestFitLine<-lm(y~EPM_DistanceTraveled, data=Temp_Data) 
    abline(BestFitLine, lwd=3)
    BestFitLineM<-lm(y~EPM_DistanceTraveled, data=Temp_Data[Temp_Data$Sex_AsFactor=="Male",]) 
    abline(BestFitLineM)
    BestFitLineF<-lm(y~EPM_DistanceTraveled, data=Temp_Data[Temp_Data$Sex_AsFactor=="Female",]) 
    abline(BestFitLineF, col=2)
    dev.off()
    
    rm(Temp_Data)
    
  }else{}
  
}


####################################################


setwd("//n05-corea-cifs.umhs.med.umich.edu/Home1/Users/hebda/My Documents/F2_TMMnorm_20220104updated/Plots_20220509_PavCADay67")


# TopHRLRF2Genes_ForGraphing<-read.csv("TopHRLRF2Genes_ForGraphing.csv", header=TRUE, stringsAsFactors = FALSE)
# 
# 
# colnames(TopHRLRF2Genes_ForGraphing)[1]<-"ENSEMBL"
# 
# colnames(F2_PCA_wMetaData)

#double check colnames of (F2_Annotated_lcpm_NoLowHits)

for(i in c(1:nrow(TopHRLRF2Genes_ForGraphing))){
  
  if(length(which(F2_Annotated_lcpm_NoLowHits$ENSEMBL==TopHRLRF2Genes_ForGraphing$ENSEMBL[i]))>0){
    
    Temp_Data<-data.frame(y=as.matrix(F2_Annotated_lcpm_NoLowHits[which(F2_Annotated_lcpm_NoLowHits$ENSEMBL==TopHRLRF2Genes_ForGraphing$ENSEMBL[i]), c(1:245)])[1,], F2_PCA_wMetaData) 
    
    pdf(paste("F2_",TopHRLRF2Genes_ForGraphing$Symbol[i], "_", TopHRLRF2Genes_ForGraphing$ENSEMBL[i], "_vs_PavCA_Day67_bySex.pdf", sep=""), width=6.5, height=6)
    plot(y~PCA_Index_Days6and.7, data=Temp_Data, ylab=paste(TopHRLRF2Genes_ForGraphing$Symbol[i], ": Log2 CPM", sep=""), col=Sex_AsFactor)
    BestFitLine<-lm(y~PCA_Index_Days6and.7, data=Temp_Data) 
    abline(BestFitLine, lwd=3)
    BestFitLineM<-lm(y~PCA_Index_Days6and.7, data=Temp_Data[Temp_Data$Sex_AsFactor=="Male",]) 
    abline(BestFitLineM)
    BestFitLineF<-lm(y~PCA_Index_Days6and.7, data=Temp_Data[Temp_Data$Sex_AsFactor=="Female",]) 
    abline(BestFitLineF, col=2)
    dev.off()
    
    rm(Temp_Data)
    
  }else{}
  
}


setwd("//n05-corea-cifs.umhs.med.umich.edu/Home1/Users/hebda/My Documents/F2_TMMnorm_20220104updated/Plots_20220131_STGT")

# 
# TopHRLRF2Genes_ForGraphing<-read.csv("TopHRLRF2Genes_ForGraphing.csv", header=TRUE, stringsAsFactors = FALSE)
# 
# 
# colnames(TopHRLRF2Genes_ForGraphing)[1]<-"ENSEMBL"
# 
# colnames(F2_PCA_wMetaData)

#double check colnames of (F2_Annotated_lcpm_NoLowHits)

for(i in c(1:nrow(TopHRLRF2Genes_ForGraphing))){
  
  if(length(which(F2_Annotated_lcpm_NoLowHits$ENSEMBL==TopHRLRF2Genes_ForGraphing$ENSEMBL[i]))>0){
    
    Temp_Data<-data.frame(y=as.matrix(F2_Annotated_lcpm_NoLowHits[which(F2_Annotated_lcpm_NoLowHits$ENSEMBL==TopHRLRF2Genes_ForGraphing$ENSEMBL[i]), c(1:245)])[1,], F2_PCA_wMetaData) 
    
    pdf(paste("F2_",TopHRLRF2Genes_ForGraphing$Symbol[i], "_", TopHRLRF2Genes_ForGraphing$ENSEMBL[i], "_vs_LearningClassification_bySex.pdf", sep=""), width=12, height=6)
    boxplot(y~Sex_AsFactor*LearningClassification_AsFactor, data=Temp_Data, ylab=paste(TopHRLRF2Genes_ForGraphing$Symbol[i], ": Log2 CPM", sep=""), xlab="", col=c("green4","green2", "grey48", "grey","firebrick4", "firebrick2", "gold3", "gold"), cex.lab=1.5, cex.axis=1.2)
    stripchart(y~Sex_AsFactor*LearningClassification_AsFactor, data=Temp_Data, vertical=TRUE, method="jitter", add=TRUE, pch=18, cex=2, col='black')
    
    dev.off()
    
    rm(Temp_Data)
    
  }else{}
  
}

#########################################2/11/2022
